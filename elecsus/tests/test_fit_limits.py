"""
Fits and the physical range of their parameters.

A fit stopped at the first trial value outside the range where the model
is defined. For a differential-evolution fit of GammaBuf from 0 the GUI
searched [-100, 100] MHz, and the fit failed with scipy's "map-like callable"
RuntimeError, which named neither the parameter nor the cause. Bounds are now
narrowed to param_validation.FIT_LIMITS (with a warning), the GUI builds its
search windows inside them (test_gui_model_controls.py), and a value that
still gets through is reported by name.

And fit_data left T, lcell and Bfield to lmfit when p_dict lacked them, and
lmfit starts a missing parameter at -inf, so a p_dict that calculate()
accepts failed with "'Bfield' must be a finite real number, got -inf".

A real differential-evolution fit takes ~17 s, so the tests that need one
either stub that stage or stop at its first trial value.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs.param_validation import (
    FIT_LIMITS, ParameterRangeError, check_ranges, clip_bounds_to_limits)

X = np.linspace(-4000, 4000, 120)
TRUE = dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, rb85frac=72.17,
            GammaBuf=20., Bfield=0.)


def _spectrum(p):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.real(calculate(X, [1, 0, 0], p, outputs=['S0'])[0])


# --- the limits table --------------------------------------------------------

def _alone(key, value):
    """A dict setting only `key`. The other potassium fraction is set to 0,
    so that the joint check (K40frac + K41frac <= 100) stays out of the way."""
    other = {'K40frac': 'K41frac', 'K41frac': 'K40frac'}.get(key)
    return {key: value, **({other: 0.} if other else {})}


@pytest.mark.parametrize('key', sorted(FIT_LIMITS))
def test_fit_limits_agree_with_check_ranges(key):
    """Each limit is the edge check_ranges enforces: the limit itself passes
    and the next float beyond it fails."""
    lo, hi = FIT_LIMITS[key]
    for limit, outward in ((lo, -np.inf), (hi, np.inf)):
        if not np.isfinite(limit):
            continue
        check_ranges(_alone(key, limit))
        with pytest.raises(ParameterRangeError, match=key):
            check_ranges(_alone(key, float(np.nextafter(limit, outward))))


def test_range_errors_are_still_value_errors():
    with pytest.raises(ValueError):
        check_ranges({'GammaBuf': -1.})


# --- bounds ------------------------------------------------------------------

def test_bounds_past_a_limit_are_narrowed_with_a_warning():
    bounds = {'GammaBuf': [-100., 100.], 'T': [-300., 50.],
              'rb85frac': [90., 110.], 'Bfield': [-500., 500.]}
    with pytest.warns(UserWarning) as record:
        clipped = clip_bounds_to_limits(bounds)
    warned = {key for key in bounds for w in record if f"['{key}']" in str(w.message)}
    assert warned == {'GammaBuf', 'T', 'rb85frac'}
    assert clipped['GammaBuf'] == [0., 100.]
    assert clipped['T'] == [FIT_LIMITS['T'][0], 50.]
    assert clipped['rb85frac'] == [90., 100.]
    assert clipped['Bfield'] == [-500., 500.]       # no physical limit
    assert bounds['GammaBuf'] == [-100., 100.]      # the caller's dict is untouched


def test_bounds_inside_the_limits_pass_silently():
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert clip_bounds_to_limits({'GammaBuf': [0., 50.]}) == {'GammaBuf': [0., 50.]}


def test_bounds_wholly_outside_the_limits_raise():
    with pytest.raises(ParameterRangeError, match='GammaBuf'):
        clip_bounds_to_limits({'GammaBuf': [-10., -5.]})


def test_fit_data_hands_the_narrowed_bounds_to_differential_evolution(monkeypatch):
    """The reported failure: DE of GammaBuf from 0 with [-100, 100] MHz."""
    seen = {}
    real_ml_fit = ML.ML_fit

    def de_stub_then_real(data, E_in, p_dict, p_dict_bools, method='leastsq', **kw):
        if method == 'differential_evolution':
            seen['bounds'] = kw['p_dict_bounds']
            return dict(p_dict, GammaBuf=20.), None
        return real_ml_fit(data, E_in, p_dict, p_dict_bools, method=method, **kw)

    monkeypatch.setattr(ML, 'ML_fit', de_stub_then_real)
    with pytest.warns(UserWarning, match='GammaBuf'):
        best, rms, _ = fit_data(
            [X, _spectrum(TRUE)], dict(TRUE, GammaBuf=0.), {'GammaBuf': True},
            E_in=[1, [0, 0]], p_dict_bounds={'GammaBuf': [-100., 100.]},
            fit_algorithm='differential_evolution')
    assert seen['bounds']['GammaBuf'] == [0., 100.]
    assert best['GammaBuf'] == pytest.approx(20., abs=1e-3)


# --- a value that gets through is named --------------------------------------

def test_out_of_range_trial_is_reported_by_name():
    """ML_fit called directly, past fit_data's narrowing: scipy used to turn
    the error into a RuntimeError about its map-like callable."""
    p = dict(TRUE, GammaBuf=0., Btheta=0., Bphi=0., shift=0.,
             E_x=1., E_y=0., E_phase=0.)
    with pytest.raises(ParameterRangeError,
                       match=r"(?s)GammaBuf.*negative.*p_dict_bounds"):
        ML.ML_fit([X, _spectrum(TRUE)], [1, [0, 0]], p, {'GammaBuf': True},
                  p_dict_bounds={'GammaBuf': [-100., 100.]},
                  method='differential_evolution')


def test_unphysical_starting_value_is_rejected_before_fitting():
    with pytest.raises(ParameterRangeError, match='GammaBuf'):
        fit_data([X, _spectrum(TRUE)], dict(TRUE, GammaBuf=-5.),
                 {'GammaBuf': True}, E_in=[1, [0, 0]])


# --- parameters missing from p_dict -------------------------------------------

def test_fit_accepts_a_p_dict_without_bfield_or_lcell():
    """calculate() fills Bfield = 0 and lcell = 75 mm by default; so must the
    fit. Before, the missing values started at -inf."""
    p = dict(Elem='Rb', Dline='D2', T=35., rb85frac=72.17)
    y = _spectrum(dict(p, T=40.))
    best, rms, _ = fit_data([X, y], p, {'T': True}, E_in=[1, [0, 0]])
    assert best['T'] == pytest.approx(40., abs=1e-3)
    assert best['Bfield'] == 0.
    assert best['lcell'] == 75e-3
