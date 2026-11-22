"""
The DE round-off and RR's starts use their own input field.

ML_fit starts E_x, E_y and E_phase from its E_in argument, not from p_dict.

- After differential evolution, fit_data ran the ML round-off from DE's best
  values but with the original E_in, so DE's result for the polarisation was
  thrown away and ML restarted it from the user's first guess.
- RR_fit scattered E_x/E_y/E_phase in each start's p_dict when they are
  fitted, then passed every start (and the final fit) the original E_in, so
  the scatter was undone.

A real DE search takes ~15 s and RR runs its starts in a process pool, so the
searches are stubbed: DE "finds" the optimum, and a stand-in pool takes RR's
start points or returns made-up results.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs import RRFittingRoutine as RR
from elecsus.libs.param_validation import FIT_PARAMETERS

X = np.linspace(-4000, 4000, 120)
TRUE = dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, rb85frac=72.17, Bfield=500.)
LEFT = np.array([1, 1j, 0]) / np.sqrt(2)     # fit form [a, [a, pi/2]]
RIGHT = np.array([1, -1j, 0]) / np.sqrt(2)   # [a, [a, -pi/2]] within [-pi, 0]


def _spectrum(E):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.real(calculate(X, E, TRUE, outputs=['S0'])[0])


# --- DE --------------------------------------------------------------------------

@pytest.mark.parametrize('E, phase, bounds', [
    (LEFT, np.pi / 2, [0., 3.]),
    (RIGHT, -np.pi / 2, [-np.pi, 0.]),     # inside the caller's own range
])
def test_de_round_off_starts_from_des_field(monkeypatch, E, phase, bounds):
    real_ml_fit = ML.ML_fit

    def de_finds_the_optimum(data, E_in, p_dict, p_dict_bools, method='leastsq', **kw):
        if method == 'differential_evolution':
            found = {k: p_dict[k] for k in FIT_PARAMETERS if k in p_dict}
            return dict(found, T=40., E_phase=phase), None
        return real_ml_fit(data, E_in, p_dict, p_dict_bools, method=method, **kw)

    monkeypatch.setattr(ML, 'ML_fit', de_finds_the_optimum)
    best, rms, result = fit_data(
        [X, _spectrum(E)], dict(TRUE, T=35.), {'T': True, 'E_phase': True},
        E_in=[2 ** -0.5, [2 ** -0.5, 0.3 * np.sign(phase)]],
        p_dict_bounds={'T': [30., 50.], 'E_phase': bounds}, fit_algorithm='DE')
    assert result.init_values['E_phase'] == pytest.approx(phase)   # was 0.3
    assert rms < 1e-8 and best['E_phase'] == pytest.approx(phase, abs=1e-6)


# --- RR --------------------------------------------------------------------------

class _Stop(Exception):
    """Raised by the capturing pool once it holds the start points."""


class _CapturingPool:
    """Stands in for multiprocessing.Pool: keeps the start points' arguments
    and stops RR_fit there."""
    args = []

    def __init__(self, *args, **kwargs):
        pass

    def map_async(self, func, args_list):
        type(self).args = list(args_list)
        raise _Stop


def test_rr_starts_carry_their_own_field(monkeypatch):
    monkeypatch.setattr(RR, 'Pool', _CapturingPool)
    np.random.seed(0)
    with pytest.raises(_Stop):
        RR.RR_fit([X, _spectrum(LEFT)], [0.7, [0.7, 1.0]], dict(TRUE),
                  {'E_phase': True}, p_dict_bounds={'E_phase': 0.3}, no_evals=50)
    fields = [args[1] for args in _CapturingPool.args]
    phases = np.array([field[1][1] for field in fields])
    assert phases == pytest.approx([args[2]['E_phase'] for args in _CapturingPool.args])
    assert np.ptp(phases) > 0.3            # all 1.0 before: the scatter was undone


class _ResultsPool:
    """Stands in for multiprocessing.Pool: returns made-up downhill results,
    the best of which (lowest reduced chi-squared) has E_phase = 1.3."""

    def __init__(self, *args, **kwargs):
        pass

    def map_async(self, func, args_list):
        results = [(10. + k, {'E_phase': 0.5}, None) for k in range(len(args_list))]
        results[3] = (1., {'E_phase': 1.3}, None)

        class _Async:
            def get(self):
                return results
        return _Async()

    def close(self):
        pass

    def join(self):
        pass


def test_rr_final_fit_starts_from_the_best_start(monkeypatch):
    monkeypatch.setattr(RR, 'Pool', _ResultsPool)
    seen = {}

    def final_fit(data, E_in, p_dict, p_dict_bools, data_type='S0', **kw):
        seen['E_in'] = E_in
        return dict(p_dict), None

    monkeypatch.setattr(ML, 'ML_fit', final_fit)
    RR.RR_fit([X, _spectrum(LEFT)], [0.7, [0.7, 1.0]], dict(TRUE),
              {'E_phase': True}, no_evals=8)
    assert seen['E_in'][1][1] == pytest.approx(1.3)      # was the original 1.0
