"""
fit_data's best parameters go straight back into calculate().

lmfit's best_values also holds the fit model's independent variable 'output'
(the data type), which is no p_dict key, so calculate(X, E, best) raised
"unknown parameter key: 'output'", and notebook 06 failed where it plots the
best fit. fit_data now returns p_dict keys only, whichever routine ran.
After differential evolution it also reported Constrain as lmfit's 1.0 when
p_dict did not set it, which calculate() rejects as a switch.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs import RRFittingRoutine as RR
from elecsus.libs.param_validation import FIT_ONLY_KEYS, KNOWN_KEYS

X = np.linspace(-4000, 4000, 120)
TRUE = dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, rb85frac=72.17, Bfield=0.)


def _data():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return [X, np.real(calculate(X, [1, 0, 0], TRUE, outputs=['S0'])[0])]


def _one_ml_fit(data, E_in, p_dict, p_dict_bools, **kw):
    """Stands in for a search routine: one ML fit, which is how each ends."""
    return ML.ML_fit(data, E_in, p_dict, p_dict_bools)


@pytest.mark.parametrize('algo', ['ML', 'DE', 'RR'])
def test_best_parameters_go_back_into_calculate(monkeypatch, algo):
    """DE and RR end with an ML fit too; their searches are stubbed (a real
    one takes 15 s or more)."""
    real_ml_fit = ML.ML_fit

    def de_stub_then_real(data, E_in, p_dict, p_dict_bools, method='leastsq', **kw):
        if method == 'differential_evolution':
            return real_ml_fit(data, E_in, p_dict, p_dict_bools)
        return real_ml_fit(data, E_in, p_dict, p_dict_bools, method=method, **kw)

    monkeypatch.setattr(ML, 'ML_fit', de_stub_then_real)
    monkeypatch.setattr(RR, 'RR_fit', _one_ml_fit)
    bounds = {'T': [30., 50.]} if algo == 'DE' else None
    best, rms, result = fit_data(_data(), dict(TRUE, T=35.), {'T': True},
                                 p_dict_bounds=bounds, fit_algorithm=algo)
    assert 'output' not in best and set(best) <= KNOWN_KEYS
    assert best['Constrain'] is True      # DE reported 1.0, which calculate rejects
    with warnings.catch_warnings():
        # E_x/E_y/E_phase, the fitted input field, warn: calculate takes E_in
        warnings.filterwarnings('ignore', message=".*set the input polarisation")
        [S0] = calculate(X, [1, 0, 0], best, outputs=['S0'])
    np.testing.assert_allclose(S0, result.best_fit, rtol=0, atol=1e-12)


def test_fitted_field_is_still_reported():
    best, rms, _ = fit_data(_data(), dict(TRUE, T=35.), {'T': True})
    assert set(FIT_ONLY_KEYS) <= set(best)
    assert best['T'] == pytest.approx(40., abs=1e-4)
