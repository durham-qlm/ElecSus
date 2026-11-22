"""
fit_data's own arguments.

- fit_data(**kw) accepted and ignored any keyword, so datatype='S1' fitted
  the S0 model to S1 data (T = 159 °C for data at 40 °C), and the documented
  verbose did nothing.
- An unknown fit_algorithm ran random restart: 'differential evolution' (with
  a space) dropped the DE bounds and ran RR.
- The Jones vector [Ex, Ey, 0] that calculate() takes raised TypeError, and a
  fit-form phase outside [0, pi) was clamped by lmfit, so right-circular
  light given as [a, [b, -pi/2]] was fitted as linear light.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import _fit_form, calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs.param_validation import ParameterRangeError
from elecsus.libs import RRFittingRoutine as RR

X = np.linspace(-4000, 4000, 120)
TRUE = dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, rb85frac=72.17, Bfield=500.)
START = dict(TRUE, T=35.)
RIGHT = np.array([1, -1j, 0]) / np.sqrt(2)
XPOL = [1., [0., 0.]]          # the fit form, which the old code took too


def _spectrum(E, output='S0', **extra):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.real(calculate(X, E, dict(TRUE, **extra), outputs=[output])[0])


# --- keyword arguments ---------------------------------------------------------

def test_misspelt_keyword_raises():
    """The reported failure: datatype='S1' was ignored and S0 fitted."""
    with pytest.raises(TypeError, match='datatype'):
        fit_data([X, _spectrum([1, 0, 0], 'S1')], dict(START), {'T': True},
                 E_in=XPOL, datatype='S1')


def test_verbose_prints(capsys):
    fit_data([X, _spectrum([1, 0, 0])], dict(START), {'T': True}, E_in=XPOL,
             verbose=True)
    out = capsys.readouterr().out
    assert 'Starting parameter dictionary' in out and 'Fit Statistics' in out


def test_quiet_by_default(capsys):
    fit_data([X, _spectrum([1, 0, 0])], dict(START), {'T': True}, E_in=XPOL)
    out = capsys.readouterr().out
    assert 'Starting parameter dictionary' not in out and 'Fit Statistics' not in out


def test_strict_and_verbose_are_keyword_only():
    with pytest.raises(TypeError):
        fit_data([X, _spectrum([1, 0, 0])], dict(START), {'T': True}, XPOL, None,
                 'S0', 'ML', True)


# --- fit_algorithm ---------------------------------------------------------------

@pytest.mark.parametrize('name, hint', [
    ('differential evolution', 'differential_evolution'), ('de', 'DE'),
    ('lm', 'LM'), ('Nelder', None), (None, None), ('', None)])
def test_unknown_algorithm_raises(name, hint):
    """Each of these used to run random restart."""
    match = r"unknown fit_algorithm" + (fr".*did you mean '{hint}'" if hint else '')
    with pytest.raises(ValueError, match=match):
        fit_data([X, _spectrum([1, 0, 0])], dict(START), {'T': True}, E_in=XPOL,
                 fit_algorithm=name)


class _Stop(Exception):
    pass


@pytest.mark.parametrize('name, routine', [
    ('ML', 'ML'), ('LM', 'ML'), ('leastsq', 'ML'),
    ('Marquardt-Levenberg', 'ML'), ('DE', 'differential_evolution'),
    ('differential_evolution', 'differential_evolution'), ('RR', 'RR')])
def test_documented_names_run_their_routine(monkeypatch, name, routine):
    ran = []

    def stand_in(label):
        """Records the routine (ML_fit's method, if not the default) and stops."""
        def stop(*args, method=None, **kw):
            ran.append(method or label)
            raise _Stop
        return stop

    monkeypatch.setattr(ML, 'ML_fit', stand_in('ML'))
    monkeypatch.setattr(RR, 'RR_fit', stand_in('RR'))
    with pytest.raises(_Stop):
        fit_data([X, _spectrum([1, 0, 0])], dict(START), {'T': True}, E_in=XPOL,
                 fit_algorithm=name)
    assert ran == [routine]


@pytest.mark.parametrize('name', ['SA', 'sa', 'simulated annealing', 'simulated_annealing'])
def test_simulated_annealing_was_removed(name):
    """SA crashed on Linux and macOS and was slower than DE."""
    with pytest.raises(ValueError, match=r"simulated annealing has been removed.*'DE'"):
        fit_data([X, _spectrum([1, 0, 0])], dict(START), {'T': True}, E_in=XPOL,
                 fit_algorithm=name)


# --- E_in ------------------------------------------------------------------------

def _jones(E):
    """E as a Jones vector, from either form."""
    if len(E) == 2:
        return np.array([E[0], E[1][0] * np.exp(1j * E[1][1]), 0])
    return np.asarray(E, dtype=complex)


_STATES = [
    [1, 0, 0], [0, 1, 0], [1, 1, 0], [1, -1, 0], [1, 1j, 0], [1, -1j, 0],
    [np.cos(-0.5), np.sin(-0.5), 0], 1j * np.array([1, 2 - 1j, 0]),
    [-1, 0.3j, 0], [0, -1j, 0],
    # fit forms with the phase outside [0, pi), and negative amplitudes
    [0.6, [0.8, -np.pi / 2]], [0.6, [0.8, 3 * np.pi / 2]], [0.6, [0.8, np.pi]],
    [0.6, [0.8, 2 * np.pi]], [0.6, [-0.8, -3 * np.pi / 2]], [-0.6, [0.8, 5.]],
    [0.6, [0.8, -7.]]]


@pytest.mark.parametrize('E', _STATES)
def test_fit_form_is_the_same_polarisation_state(E):
    """The fit form, rebuilt as a Jones vector, gives the same outputs as E,
    with its phase inside the range the fit allows."""
    E_x, (E_y, phase) = _fit_form(E)
    assert 0 <= phase <= ML.E_PHASE_MAX
    back = _jones([E_x, [E_y, phase]])
    for output in ('S0', 'S1', 'S2', 'S3'):
        np.testing.assert_allclose(_spectrum(back, output),
                                   _spectrum(_jones(E), output), rtol=0, atol=1e-12)


@pytest.mark.parametrize('E', [
    [1, [1, np.pi - 5e-5]], [1, np.exp(1j * (np.pi - 5e-5)), 0], [1, [1, -1e-6]]])
def test_phase_just_below_pi_is_clamped_to_the_bound(E):
    """The default range leaves out (pi - 1e-4, pi): there the phase is
    clamped, as lmfit would, which moves the state by under 1e-4 rad."""
    E_x, (E_y, phase) = _fit_form(E)
    assert 0 <= phase <= ML.E_PHASE_MAX
    np.testing.assert_allclose(_spectrum(_jones([E_x, [E_y, phase]]), 'S3'),
                               _spectrum(_jones(E), 'S3'), rtol=0, atol=1e-3)


def test_fit_form_in_range_is_unchanged():
    """Including forms the Jones round trip would rewrite."""
    assert _fit_form([1, [-1, 0.5]]) == [1., [-1., 0.5]]
    assert _fit_form([-1, [1, 0.5]]) == [-1., [1., 0.5]]
    assert _fit_form([0.7, [0.7, 0.]]) == [0.7, [0.7, 0.]]


@pytest.mark.parametrize('E, expected', [
    ([0, 1j, 0], [0., [1., 0.]]), ([0, -1j, 0], [0., [1., 0.]]),
    ([1j, 0, 0], [1., [0., 0.]]), ([0, [1, -np.pi / 2]], [0., [1., 0.]])])
def test_zero_amplitude_gives_zero_phase(E, expected):
    """With one amplitude zero the relative phase means nothing. Taking it
    from the other component's global phase made a fit of the zero amplitude
    explore only ellipses, never the linear states between."""
    assert _fit_form(E) == expected


def test_fit_data_takes_a_jones_vector():
    """It raised TypeError: 'int' object is not subscriptable."""
    best, rms, _ = fit_data([X, _spectrum(RIGHT)], dict(START), {'T': True},
                            E_in=RIGHT)
    assert rms < 1e-8 and best['T'] == pytest.approx(40., abs=1e-6)


def test_fit_form_phase_below_zero_is_not_clamped():
    """lmfit clamped -pi/2 to 0: right-circular light became linear."""
    linear = _spectrum(np.array([1, 1, 0]) / np.sqrt(2))
    assert np.max(abs(linear - _spectrum(RIGHT))) > 0.05   # the test can tell
    best, rms, _ = fit_data([X, _spectrum(RIGHT)], dict(START), {'T': True},
                            E_in=[2 ** -0.5, [2 ** -0.5, -np.pi / 2]])
    assert rms < 1e-8 and best['T'] == pytest.approx(40., abs=1e-6)


@pytest.mark.parametrize('E', [
    [1, 0, 0.5], [1, 2], 'x', [1, [1]], [[1, 0, 0]], [1, 0, np.nan],
    [1, [0, np.inf]], [np.complex128(1j), [1, 0]], [1, np.array([1j, 0])],
    [1j, [1, 0]]])
def test_malformed_E_in_raises(E):
    with pytest.raises(ValueError, match=r"E_in must be the fit form"):
        fit_data([X, _spectrum(RIGHT)], dict(START), {'T': True}, E_in=E)


# --- the caller's own E_phase bounds, and the ES3 path ------------------------------

NEAR_RIGHT = np.array([1, np.exp(-1j * (np.pi / 2 - 0.2)), 0]) / np.sqrt(2)


@pytest.mark.parametrize('E', [NEAR_RIGHT, [2 ** -0.5, [2 ** -0.5, -np.pi / 2 + 0.2]]])
def test_phase_stays_inside_the_callers_E_phase_bounds(E):
    """ML_fit replaces its default [0, pi) with the caller's bounds. Moving the
    phase into [0, pi) regardless put it outside [-pi, 0], where lmfit
    clamped it to the opposite handedness (found in review).

    The start is right-handed but 0.2 rad from circular. At exactly +-pi/2
    the Faraday spectrum is even in the phase, so its gradient is zero and
    leastsq can stop at once: for one handedness or the other, depending on
    the sense of B (the 4.0.0b1 geometry fix swapped which)."""
    best, rms, _ = fit_data([X, _spectrum(RIGHT)], dict(START),
                            {'T': True, 'E_phase': True}, E_in=E,
                            p_dict_bounds={'E_phase': [-np.pi, 0.]})
    # the opposite handedness leaves rms ~0.3
    assert rms < 1e-8 and best['T'] == pytest.approx(40., abs=1e-6)
    assert best['E_phase'] == pytest.approx(-np.pi / 2, abs=1e-4)


def test_bounds_that_cannot_hold_the_state_raise():
    with pytest.raises(ValueError, match=r"p_dict_bounds\['E_phase'\].*cannot hold"):
        fit_data([X, _spectrum(RIGHT)], dict(START), {'T': True, 'E_phase': True},
                 E_in=RIGHT, p_dict_bounds={'E_phase': [0.2, 0.5]})


def test_es3_field_in_p_dict_is_converted():
    """E_in=None takes E_x/E_y/E_phase from p_dict: right-circular light there
    was clamped to linear, as in the fit form."""
    p = dict(START, E_x=2 ** -0.5, E_y=2 ** -0.5, E_phase=-np.pi / 2)
    best, rms, _ = fit_data([X, _spectrum(RIGHT)], p, {'T': True})
    assert rms < 1e-8 and best['T'] == pytest.approx(40., abs=1e-6)


def test_bad_es3_field_in_p_dict_is_named():
    with pytest.raises(ParameterRangeError, match="'E_x'"):
        fit_data([X, _spectrum(RIGHT)], dict(START, E_x=np.nan, E_y=0., E_phase=0.),
                 {'T': True})
