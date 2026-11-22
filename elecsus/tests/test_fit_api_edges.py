"""
Fitting and API edges.

- _fit_form crashed (OverflowError) or hung on infinite or very wide
  E_phase bounds, refused an x-polarised start under phase bounds that
  exclude 0, raised UnboundLocalError for a scalar E_in, accepted a zero
  field, and flipped E_x past the caller's own bounds on E_x.
- fit_data let E_x and E_y both vary, though only their ratio is
  observable; the fit returned arbitrary amplitudes with ~1e6 error bars.
- a p_dict with E_x but not E_y/E_phase raised a bare KeyError.
- the docstrings said data was Nx2; an Nx2 array fitted two points.
- fit_data raised KeyError: 'Elem' for a p_dict without Elem or Dline.
- calculate() ignored BufferPressure without BufferGas in silence.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import _fit_form, calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs.param_validation import FIT_PARAMETERS

X = np.linspace(-4000, 4000, 120)
TRUE = dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, Bfield=500.)
START = dict(TRUE, T=35.)
RIGHT = np.array([1, -1j, 0]) / np.sqrt(2)
# right-handed elliptical light, as the fit form with E_x and E_y both > 0
ELLIPSE = [0.6, [0.8, -np.pi / 2]]


def _spectrum(E, output='S0', p=TRUE):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.real(calculate(X, E, p, outputs=[output])[0])


def _jones(E):
    if len(E) == 2:
        return np.array([E[0], E[1][0] * np.exp(1j * E[1][1]), 0])
    return np.asarray(E, dtype=complex)


def _same_state(E, F, atol=1e-12):
    """E and F are the same polarisation state, with the same intensity:
    |Ex|^2, |Ey|^2 and Ex conj(Ey) do not depend on the global phase."""
    a, b = _jones(E), _jones(F)
    np.testing.assert_allclose(
        [abs(a[0]) ** 2, abs(a[1]) ** 2, a[0] * np.conj(a[1])],
        [abs(b[0]) ** 2, abs(b[1]) ** 2, b[0] * np.conj(b[1])], rtol=0, atol=atol)


# --- _fit_form ------------------------------------------------------------------

_STATES = [[1, 1j, 0], [1, -1j, 0], [0.6, 0.8 * np.exp(2.5j), 0], [1, [1, -0.5]],
           [-0.6, [0.8, 5.]], ELLIPSE, [0.6, [-0.8, 7.]]]
_WIDE = [(-np.inf, np.inf), (0., np.inf), (-np.inf, 0.), (-1e15, 1e15), (-1e6, -1e6 + 4.)]


@pytest.mark.parametrize('phase_range', _WIDE)
@pytest.mark.parametrize('E', _STATES)
def test_infinite_and_wide_phase_bounds(E, phase_range):
    """These raised OverflowError (infinite) or hung (1e15: the old code walked
    every multiple of pi in the range)."""
    E_x, (E_y, phase) = _fit_form(E, phase_range)
    assert phase_range[0] <= phase <= phase_range[1]
    # a phase near -1e6 rad carries ~1e-10 rad of float rounding
    _same_state(E, [E_x, [E_y, phase]], atol=1e-12 * max(1., abs(phase)))


def test_wide_bounds_keep_the_phase_nearest_its_start():
    """The multiple of pi nearest 0: a phase already in range is not moved."""
    assert _fit_form([1, 1j, 0], (-1e15, 1e15)) == [1., [1., np.pi / 2]]
    assert _fit_form([1, [1, -0.5]], (-np.inf, np.inf)) == [1., [1., -0.5]]


@pytest.mark.parametrize('E', [[1, [0, 0]], [1, 0, 0], [0, 1j, 0], [0, [1, 3.]]])
def test_zero_amplitude_takes_the_nearest_phase_in_range(E):
    """With one amplitude zero every phase is the same state. Bounds that
    exclude 0 raised 'cannot hold'."""
    E_x, (E_y, phase) = _fit_form(E, (0.5, 1.0))
    assert phase == 0.5
    _same_state(E, [E_x, [E_y, phase]])


def test_in_range_fit_form_is_returned_unchanged_against_the_callers_range():
    assert _fit_form([1, [1, -0.5]], (-1., 0.)) == [1., [1., -0.5]]
    assert _fit_form([-1, [2, -0.5]], (-1., 0.), (-2., 0.), (0., 3.)) == [-1., [2., -0.5]]
    # in the phase range but not E_x's: the sign moves to E_y, which can take it
    assert _fit_form([-1, [2, 0.]], (-np.pi, np.pi), (0., 2.)) == [1., [-2., 0.]]


@pytest.mark.parametrize('E, phase_range, expected', [
    # within 1e-4 rad of the caller's range: clamped, as lmfit would
    ([1, [1, 1.00005]], (0., 1.), [1., [1., 1.]]),
    ([1, [1, -0.00005]], (0., 1.), [1., [1., 0.]]),
    ([1, [1, 1.00005 + np.pi]], (0., 1.), [-1., [1., 1.]]),
])
def test_phase_just_outside_the_callers_range_is_clamped(E, phase_range, expected):
    out = _fit_form(E, phase_range)
    assert out[0] == expected[0] and out[1][0] == expected[1][0]
    assert out[1][1] == pytest.approx(expected[1][1], abs=1e-12)


@pytest.mark.parametrize('phi', [1.0002, -0.0002, 1.5, np.pi - 0.2])
def test_phase_further_outside_the_callers_range_raises(phi):
    with pytest.raises(ValueError, match=r"p_dict_bounds\['E_phase'\] = \[0, 1\] cannot hold"):
        _fit_form([1, [1, phi]], (0., 1.))


def test_the_flip_respects_the_callers_E_x_bounds():
    """E_x cannot go negative, so E_y changes sign instead."""
    E_x, (E_y, phase) = _fit_form(ELLIPSE, None, (0., 1.))
    assert (E_x, E_y) == (0.6, -0.8) and 0 <= phase <= ML.E_PHASE_MAX
    _same_state(ELLIPSE, [E_x, [E_y, phase]])
    # and E_x where E_y cannot
    E_x, (E_y, phase) = _fit_form(ELLIPSE, None, None, (0., 1.))
    assert (E_x, E_y) == (-0.6, 0.8)
    _same_state(ELLIPSE, [E_x, [E_y, phase]])


def test_bounds_that_hold_no_form_of_the_state_raise():
    """With E_x, E_y >= 0 and the phase in [0, pi) the fit form holds no
    right-handed light at all."""
    with pytest.raises(ValueError, match=r"p_dict_bounds\['E_x'\] = \[0, 1\], "
                                         r"p_dict_bounds\['E_y'\] = \[0, 1\]"):
        _fit_form(ELLIPSE, None, (0., 1.), (0., 1.))
    with pytest.raises(ValueError, match=r"p_dict_bounds\['E_x'\] = \[0.7, 1\]"):
        _fit_form([0.6, [0.8, 0.]], None, (0.7, 1.))


@pytest.mark.parametrize('E', [1.0, 0, np.float64(2.), None])
def test_scalar_E_in_raises_the_documented_error(E):
    """It raised UnboundLocalError."""
    with pytest.raises(ValueError, match=r"E_in must be the fit form"):
        _fit_form(E)


@pytest.mark.parametrize('E', [[0, 0, 0], [0, [0, 0]], [0., [0., 7.]], np.zeros(3)])
def test_zero_field_raises(E):
    """It was accepted, and the fit then failed inside lmfit."""
    with pytest.raises(ValueError, match=r"zero field"):
        _fit_form(E)


def test_fit_with_infinite_phase_bounds():
    """The first case: OverflowError from fit_data. (The start is 0.2
    rad from circular: at exactly circular the gradient in E_phase is zero,
    and leastsq can stop at once; see test_fit_arguments.)"""
    start = np.array([1, np.exp(-1j * (np.pi / 2 - 0.2)), 0]) / np.sqrt(2)
    best, rms, _ = fit_data([X, _spectrum(RIGHT)], dict(START),
                            {'T': True, 'E_phase': True}, E_in=start,
                            p_dict_bounds={'E_phase': [-np.inf, np.inf]})
    assert best['T'] == pytest.approx(40., abs=1e-6) and rms < 1e-8
    assert best['E_phase'] == pytest.approx(-np.pi / 2, abs=1e-4)


def test_fit_of_E_x_bounded_to_positive_values():
    """The reported fit: the flipped start E_x = -0.6 was clamped to the bound
    0, and the fit stayed at T = 35 with rms 0.28."""
    data = [X, _spectrum(_jones(ELLIPSE))]
    best, rms, _ = fit_data(data, dict(START), {'T': True, 'E_x': True},
                            E_in=ELLIPSE, p_dict_bounds={'E_x': [0., 1.]})
    assert best['T'] == pytest.approx(40., abs=1e-6) and rms < 1e-8
    assert best['E_x'] >= 0


def test_differential_evolution_keeps_the_callers_E_x_bounds(monkeypatch):
    """DE takes the caller's E_x bounds too, at its start and at the ML
    round-off that follows it. A real DE run takes over a minute, so the DE
    stage is stubbed: it records its start, and returns the optimum in the
    fit form with E_x > 0 and the phase below 0; the round-off is real."""
    real_ml_fit = ML.ML_fit
    starts = []

    def de_stub_then_real(data, E_in, p_dict, p_dict_bools, method='leastsq', **kw):
        starts.append(E_in)
        if method == 'differential_evolution':
            found = {k: p_dict[k] for k in FIT_PARAMETERS if k in p_dict}
            return dict(found, T=40., E_x=0.6, E_y=0.8, E_phase=-np.pi / 2), None
        return real_ml_fit(data, E_in, p_dict, p_dict_bools, method=method, **kw)

    monkeypatch.setattr(ML, 'ML_fit', de_stub_then_real)
    best, rms, _ = fit_data([X, _spectrum(_jones(ELLIPSE))], dict(START),
                            {'T': True, 'E_x': True}, E_in=ELLIPSE,
                            p_dict_bounds={'T': [30., 50.], 'E_x': [0., 1.]},
                            fit_algorithm='DE')
    assert [E[0] for E in starts] == [0.6, 0.6]      # the DE start was -0.6
    assert best['T'] == pytest.approx(40., abs=1e-6) and rms < 1e-8
    assert best['E_x'] >= 0


def test_bounds_that_hold_no_form_raise_from_fit_data():
    with pytest.raises(ValueError, match=r"cannot be written as"):
        fit_data([X, _spectrum(RIGHT)], dict(START), {'T': True, 'E_x': True, 'E_y': False},
                 E_in=RIGHT, p_dict_bounds={'E_x': [0., 1.], 'E_y': [0., 1.]})


# --- E_x and E_y both varied --------------------------------------------------

def test_E_x_and_E_y_both_varied_raise():
    with pytest.raises(ValueError, match=r"both 'E_x' and 'E_y'.*only their ratio"):
        fit_data([X, _spectrum([0.6, 0.8, 0], 'S1')], dict(TRUE),
                 {'E_x': True, 'E_y': True}, E_in=[0.5, [0.9, 0.]], data_type='S1')


@pytest.mark.parametrize('key', ['E_x', 'E_y'])
def test_one_amplitude_varied_fits_the_ratio(key):
    """The supported way: vary one amplitude (as the GUI does)."""
    E_in = [0.6, [0.9, 0.]] if key == 'E_x' else [0.5, [0.8, 0.]]
    best, rms, res = fit_data([X, _spectrum([0.6, 0.8, 0], 'S1')], dict(TRUE),
                              {key: True}, E_in=E_in, data_type='S1')
    assert best['E_y'] / best['E_x'] == pytest.approx(0.8 / 0.6, rel=1e-6)
    assert rms < 1e-8


# --- a partial set of ES3's field keys ----------------------------------------

@pytest.mark.parametrize('keys, missing', [
    (dict(E_x=1.), "'E_y', 'E_phase'"), (dict(E_x=1., E_y=0.), "'E_phase'"),
    (dict(E_x=1., E_phase=0.), "'E_y'"), (dict(E_y=1.), "'E_x', 'E_phase'"),
    (dict(E_y=1., E_phase=0.5), "'E_x'"), (dict(E_phase=0.5), "'E_x', 'E_y'")])
def test_partial_es3_field_raises_naming_the_missing_keys(keys, missing):
    """E_x alone raised KeyError: 'E_y'; E_y/E_phase alone were ignored."""
    with pytest.raises(ValueError, match=fr"but not {missing}: ES3's input field takes all three"):
        fit_data([X, _spectrum([1, 0, 0])], dict(START, **keys), {'T': True})


def test_full_es3_field_still_works():
    best, rms, _ = fit_data([X, _spectrum([1, 0, 0])],
                            dict(START, E_x=1., E_y=0., E_phase=0.), {'T': True})
    assert best['T'] == pytest.approx(40., abs=1e-6) and rms < 1e-8


# --- the shape of data ---------------------------------------------------------

def test_Nx2_data_raises():
    """It was fitted as two data points: T = 4545 for data at 50."""
    y = _spectrum([1, 0, 0])
    with pytest.raises(ValueError, match=r"shape \(120, 2\), which is Nx2: pass its transpose"):
        fit_data(np.column_stack([X, y]), dict(START), {'T': True})


@pytest.mark.parametrize('data', [
    [X[:-1], _spectrum([1, 0, 0])], [X], [X, X, X], [np.zeros((2, 3)), np.zeros((2, 3))],
    _spectrum([1, 0, 0]), 'data'])
def test_malformed_data_raises(data):
    with pytest.raises(ValueError, match=r"data must be \[x, y\]"):
        fit_data(data, dict(START), {'T': True})


@pytest.mark.parametrize('wrap', [
    lambda x, y: [x, y], lambda x, y: (x, y), lambda x, y: np.array([x, y]),
    lambda x, y: [list(x), list(y)], lambda x, y: (tuple(x), tuple(y))])
def test_2xN_data_fits_in_any_container(wrap):
    best, rms, res = fit_data(wrap(X, _spectrum([1, 0, 0])), dict(START), {'T': True})
    assert res.ndata == len(X)
    assert best['T'] == pytest.approx(40., abs=1e-6) and rms < 1e-8


# --- Elem and Dline default to Rb D2 --------------------------------------------

@pytest.mark.parametrize('drop', [('Elem',), ('Dline',), ('Elem', 'Dline')])
def test_fit_without_Elem_or_Dline_is_Rb_D2(drop):
    """calculate() defaults to Rb D2; fit_data raised KeyError: 'Elem'."""
    data = [X, _spectrum([1, 0, 0])]
    given, _, _ = fit_data(data, dict(START), {'T': True})
    start = {key: value for key, value in START.items() if key not in drop}
    best, rms, _ = fit_data(data, start, {'T': True})
    assert best['Elem'] == 'Rb' and best['Dline'] == 'D2'
    assert best['T'] == pytest.approx(given['T'], abs=1e-9) and rms < 1e-8


# --- buffer-gas settings without BufferGas ---------------------------------------

P_BUF = dict(Elem='Rb', Dline='D2', T=60.)
X_BUF = np.linspace(-3000, 3000, 7)


def _caught(p):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        calculate(X_BUF, [1, 0, 0], p, outputs=['S0'])
    return [str(w.message) for w in caught if 'without BufferGas' in str(w.message)]


@pytest.mark.parametrize('extra, named', [
    (dict(BufferPressure=60.), 'BufferPressure = 60.0 Torr has'),
    (dict(BufferTemp_fill=300.), 'BufferTemp_fill = 300.0 K has'),
    (dict(BufferPressure=60., BufferTemp_fill=300.), 'BufferPressure = 60.0 Torr and '
                                                     'BufferTemp_fill = 300.0 K have'),
    (dict(BufferPressure=60., BufferGas=None), 'BufferPressure = 60.0 Torr has')])
def test_buffer_settings_without_a_gas_warn(extra, named):
    """The reported case: no warning, and the spectrum without the gas."""
    caught = _caught(dict(P_BUF, **extra))
    assert len(caught) == 1 and named in caught[0]


@pytest.mark.parametrize('extra', [
    dict(), dict(BufferPressure=0.), dict(BufferTemp_fill=293.15),
    dict(BufferGas='Ne', BufferPressure=60., BufferTemp_fill=300.)])
def test_no_warning_where_nothing_is_ignored(extra):
    assert _caught(dict(P_BUF, **extra)) == []


def test_fits_without_a_buffer_gas_do_not_warn():
    """The fit model always passes BufferPressure = 0 and the default fill
    temperature, with no gas."""
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        fit_data([X, _spectrum([1, 0, 0])], dict(START), {'T': True})
