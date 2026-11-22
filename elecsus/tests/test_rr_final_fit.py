"""
Random restart's last fit runs the model the caller asked for.

RR_fit ends with an ML fit from the best start, for the lmfit result. Its
p_dict was the start's, updated with all of lmfit's best_values, which also
holds the fit function's independent variables at their signature defaults:
Elem='Rb', Dline='D2', BufferGas='', and Constrain as a float. So the last fit
was a Rb D2 fit without the buffer gas (Cs D1 at 50 C came back as T = 4.0),
and a p_dict that set Constrain stopped with a ParameterRangeError about
'Constrain' = 1.0. fit_data then wrote the caller's Elem, Dline and BufferGas
into the result, so it named the model asked for, with values from another.

The DE round-off had been patched for this key by key; both now take
the fit parameters only (MLFittingRoutine.fitted_values).

The searches run for real, in this process and with three starts.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs import RRFittingRoutine as RR

from .test_random_restart import _SerialPool

CASES = {
    'Cs D1': (dict(Elem='Cs', Dline='D1', T=50., lcell=75e-3), 'T', 45., 5.),
    'Rb D2 buffer gas': (dict(Elem='Rb', Dline='D2', T=60., lcell=75e-3,
                              BufferGas='N2', BufferPressure=10.), 'BufferPressure', 8., 3.),
    'Constrain True': (dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, Constrain=True),
                       'T', 35., 5.),
    'Constrain False': (dict(Elem='K', Dline='D2', T=80., lcell=20e-3, Constrain=False,
                             DoppTemp=120.), 'T', 70., 5.),
}


def _data(p):
    x = np.linspace(-8000, 8000, 200)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        y = calculate(x, [1, 0, 0], p, outputs=['S0'])[0].real
    return [x, y + np.random.default_rng(0).normal(0, 0.002, len(x))]


def _fit(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return fit_data(*args, **kwargs)


@pytest.fixture
def three_serial_starts(monkeypatch):
    monkeypatch.setattr(RR, 'Pool', _SerialPool)
    rr_fit = RR.RR_fit
    monkeypatch.setattr(RR, 'RR_fit', lambda *a, **k: rr_fit(*a, **dict(k, no_evals=3)))
    np.random.seed(0)


@pytest.mark.parametrize('case', CASES)
def test_rr_agrees_with_ml(three_serial_starts, case):
    true, key, start, width = CASES[case]
    data = _data(true)
    ml, ml_rms, ml_result = _fit(data, dict(true, **{key: start}), {key: True})
    rr, rr_rms, rr_result = _fit(data, dict(true, **{key: start}), {key: True},
                                 p_dict_bounds={key: width}, fit_algorithm='RR')
    stderr = rr_result.params[key].stderr
    assert stderr is not None and stderr > 0   # it was None: P had no effect
    assert rr[key] == pytest.approx(true[key], abs=3 * stderr)
    assert rr[key] == pytest.approx(ml[key], rel=1e-5)
    assert rr_rms == pytest.approx(ml_rms, rel=1e-3)
    np.testing.assert_allclose(rr_result.best_fit, ml_result.best_fit, atol=1e-5)
    for name in ('Elem', 'Dline', 'BufferGas'):
        assert rr.get(name) == true.get(name, rr.get(name))


def test_fitted_values_holds_the_fit_parameters_only():
    best_values = {'T': 50., 'BufferPressure': 9.9, 'E_x': 1., 'Elem': 'Rb', 'Dline': 'D2',
                   'Constrain': 1.0, 'BufferGas': '', 'output': 'S0'}
    assert ML.fitted_values(best_values) == {'T': 50., 'BufferPressure': 9.9, 'E_x': 1.}


def test_de_round_off_keeps_the_model():
    """The DE round-off, rebuilt the same way: Cs D1, a buffer gas, and
    Constrain=False with its own Doppler temperature."""
    true = dict(Elem='Cs', Dline='D1', T=55., lcell=20e-3, BufferGas='Ne',
                BufferPressure=20., Constrain=False, DoppTemp=80.)
    data = _data(true)
    np.random.seed(0)
    de, rms, result = _fit(data, dict(true, T=45., BufferPressure=10.),
                           {'T': True, 'BufferPressure': True},
                           p_dict_bounds={'T': [30, 80], 'BufferPressure': [1, 50]},
                           fit_algorithm='DE')
    for key in ('T', 'BufferPressure'):
        assert de[key] == pytest.approx(true[key], abs=3 * result.params[key].stderr)
    assert rms < 0.0025
    assert de['Constrain'] is False and de['BufferGas'] == 'Ne' and de['Elem'] == 'Cs'
