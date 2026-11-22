"""
fit_data checks its keys as calculate() does.

calculate() rejects a p_dict key it does not read. fit_data handed the dict to
lmfit, which ignores unknown names, so a mistyped key left its parameter at
the default and the fit converged on the wrong model in silence: T = 70.7 °C
for data at 40 °C with 'Bfeild' for 'Bfield'. Related gaps:

- BoltzmannFactor never reached the fit model, which applied the Boltzmann
  factor whatever p_dict said;
- a mistyped key in p_dict_bools raised a bare KeyError from lmfit;
- a mistyped key in p_dict_bounds was ignored, leaving the fit unbounded;
- the ES3-only keys theta0/Pol warned in calculate() but not in fits.
"""
import warnings

import lmfit
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs.param_validation import FIT_PARAMETERS, KNOWN_KEYS

X = np.linspace(-4000, 4000, 120)
TRUE = dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, rb85frac=72.17, Bfield=0.)
START = dict(TRUE, T=35.)


def _data(**extra):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return [X, np.real(calculate(X, [1, 0, 0], dict(TRUE, **extra),
                                     outputs=['S0'])[0])]


def _fit(p_dict, p_dict_bools=None, data=None, **kw):
    return fit_data(data or _data(), p_dict, p_dict_bools or {'T': True},
                    E_in=[1, [0, 0]], **kw)


def test_fit_parameters_are_lmfits():
    """FIT_PARAMETERS is what p_dict_bools and p_dict_bounds may name; it
    must match the parameters lmfit makes from the fit function."""
    assert set(FIT_PARAMETERS) == set(lmfit.Model(ML.fit_function).param_names)
    assert set(FIT_PARAMETERS) <= KNOWN_KEYS


# --- p_dict ------------------------------------------------------------------

def test_mistyped_key_is_refused():
    """The reported failure: the fit used Bfield = 0 and converged."""
    p = {k: v for k, v in START.items() if k != 'Bfield'}
    with pytest.raises(ValueError, match=r"'Bfeild' \(did you mean 'Bfield'\?\)"):
        _fit(dict(p, Bfeild=500.))


def test_strict_false_carries_an_extra_key_with_a_warning():
    with pytest.warns(UserWarning, match='run_label'):
        best, rms, _ = _fit(dict(START, run_label='cell 3'), strict=False)
    assert best['T'] == pytest.approx(40., abs=1e-4)


def test_strict_false_keeps_an_extra_key_from_the_routines(monkeypatch):
    """The routine gets the p_dict without the extra key. Simulated
    annealing, removed since, passed it on to get_spectra, which would
    have raised on it."""
    from elecsus.libs import RRFittingRoutine as RR
    seen = {}

    def stand_in(data, E_in, p_dict, p_dict_bools, **kw):
        seen['keys'] = set(p_dict)
        return ML.ML_fit(data, E_in, p_dict, p_dict_bools)

    monkeypatch.setattr(RR, 'RR_fit', stand_in)
    with pytest.warns(UserWarning, match='run_label'):
        _fit(dict(START, run_label='cell 3'), fit_algorithm='RR', strict=False)
    assert 'run_label' not in seen['keys'] and 'T' in seen['keys']


def test_deprecated_key_warns_in_fits_too():
    with pytest.warns(UserWarning, match='theta0'):
        _fit(dict(START, theta0=45.))


# --- BoltzmannFactor -----------------------------------------------------------

def _boltzmann_off_model(best):
    return _data(T=best['T'], BoltzmannFactor=False)[1]


def test_BoltzmannFactor_reaches_the_fit():
    data = _data(BoltzmannFactor=False)
    best, rms, result = _fit(dict(START, BoltzmannFactor=False), data=data)
    np.testing.assert_allclose(result.best_fit, _boltzmann_off_model(best),
                               rtol=0, atol=1e-10)
    assert best['BoltzmannFactor'] is False
    assert best['T'] == pytest.approx(40., abs=1e-4)


def test_differential_evolution_round_off_keeps_BoltzmannFactor(monkeypatch):
    """DE's best values hold only the fit parameters, and fit_data rebuilds
    the p_dict for the ML round-off from them. A real DE run takes ~15 s, so
    the DE stage is stubbed to return exactly that; the round-off is real."""
    real_ml_fit = ML.ML_fit

    def de_stub_then_real(data, E_in, p_dict, p_dict_bools, method='leastsq', **kw):
        if method == 'differential_evolution':
            return {k: p_dict[k] for k in FIT_PARAMETERS if k in p_dict}, None
        return real_ml_fit(data, E_in, p_dict, p_dict_bools, method=method, **kw)

    monkeypatch.setattr(ML, 'ML_fit', de_stub_then_real)
    best, rms, result = _fit(dict(START, BoltzmannFactor=False),
                             data=_data(BoltzmannFactor=False),
                             p_dict_bounds={'T': [30., 50.]},
                             fit_algorithm='differential_evolution')
    np.testing.assert_allclose(result.best_fit, _boltzmann_off_model(best),
                               rtol=0, atol=1e-10)


# --- p_dict_bools and p_dict_bounds -------------------------------------------

@pytest.mark.parametrize('bools, bad', [
    ({'Temp': True}, 'Temp'),          # was a bare KeyError from lmfit
    ({'T': True, 'Elem': False}, 'Elem'),
])
def test_p_dict_bools_names_only_fit_parameters(bools, bad):
    with pytest.raises(ValueError, match=fr"p_dict_bools: '{bad}'.* is not a "
                                         r"fit parameter"):
        _fit(dict(START), bools)


@pytest.mark.parametrize('algo, bounds', [
    ('ML', {'Temp': [30., 50.]}),
    ('differential_evolution', {'T': [30., 50.], 'Temp': [30., 50.]}),
    ('RR', {'Temp': 5.}),
])
def test_p_dict_bounds_names_only_fit_parameters(algo, bounds):
    """A misnamed bound used to be ignored: the fit ran unbounded, or with
    the default search width."""
    with pytest.raises(ValueError, match=r"p_dict_bounds: 'Temp'.* is not a "
                                         r"fit parameter"):
        _fit(dict(START), p_dict_bounds=bounds, fit_algorithm=algo)


def test_strict_false_drops_a_misnamed_fit_key_with_a_warning():
    with pytest.warns(UserWarning, match="p_dict_bools: 'Tmep'"):
        best, rms, _ = _fit(dict(START), {'T': True, 'Tmep': True}, strict=False)
    assert best['T'] == pytest.approx(40., abs=1e-4)
