"""
What p_dict_bounds means in each fitting algorithm.

ML and differential evolution take a [min, max] pair for each parameter.
Random restart (RR) takes the width of the random search around each starting
value, as its docstring and ES3's fitting demo have it
(test_fitting._compare_fit_methods_demo: {'T': [30]}). So did simulated
annealing, until it was removed.

- fit_data used to pass no widths to RR, which ignored them in
  silence.
- Then every entry was narrowed to the physical limits as a pair, so the demo's
  form failed with "not enough values to unpack".
- RR itself raised KeyError when only some of the fit parameters had a
  width, and TypeError for the one-element form.

Each form is now checked where it is used, and the other form is refused by
name. RR evaluates its start points in a multiprocessing pool, so the tests
take the start points from a stand-in pool rather than running the search.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs import RRFittingRoutine as RR
from elecsus.libs.param_validation import clip_bounds_to_limits, search_widths

X = np.linspace(-4000, 4000, 120)
TRUE = dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, rb85frac=72.17, Bfield=0.)
START = dict(TRUE, T=35.)


def _data():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return [X, np.real(calculate(X, [1, 0, 0], TRUE, outputs=['S0'])[0])]


class _Stop(Exception):
    """Raised by the stand-in pool once it holds the start points."""


class _CapturingPool:
    """Stands in for multiprocessing.Pool: keeps each start point's p_dict
    and stops the routine there."""
    starts = []

    def __init__(self, *args, **kwargs):
        pass

    def map_async(self, func, args_list):
        type(self).starts = [args[2] for args in args_list]
        raise _Stop


def _starts(monkeypatch, p_dict_bools, p_dict_bounds):
    """The start points that RR draws."""
    monkeypatch.setattr(RR, 'Pool', _CapturingPool)
    np.random.seed(0)
    with pytest.raises(_Stop):
        RR.RR_fit(_data(), [1, [0, 0]], dict(START), p_dict_bools,
                  p_dict_bounds=p_dict_bounds, no_evals=200)
    return _CapturingPool.starts


# --- random restart: search widths -----------------------------------------------

def test_start_points_lie_within_the_given_widths(monkeypatch):
    """A width for T only, in the demo's one-element form: T is scattered
    over it, and Bfield over its default width (it used to raise TypeError,
    and KeyError for Bfield)."""
    starts = _starts(monkeypatch, {'T': True, 'Bfield': True}, {'T': [0.5]})
    T = np.array([p['T'] for p in starts])
    B = np.array([p['Bfield'] for p in starts])
    assert np.all(abs(T - 35.) <= 0.5) and np.ptp(T) > 0.5
    width_B = RR.p_dict_bounds_default['Bfield']
    assert np.all(abs(B) <= width_B) and np.ptp(B) > width_B


def test_fit_data_passes_the_widths_on(monkeypatch):
    """The reported failure: the demo's {'T': [30]} raised "not enough values
    to unpack", and before that never reached the routine."""
    seen = {}

    def stand_in(data, E_in, p_dict, p_dict_bools, p_dict_bounds=None, **kw):
        seen['bounds'] = p_dict_bounds
        return ML.ML_fit(data, E_in, p_dict, p_dict_bools)

    monkeypatch.setattr(RR, 'RR_fit', stand_in)
    best, rms, _ = fit_data(_data(), dict(START), {'T': True}, E_in=[1, [0, 0]],
                            p_dict_bounds={'T': [30]}, fit_algorithm='RR')
    assert seen['bounds'] == {'T': [30]}
    assert best['T'] == pytest.approx(40., abs=1e-3)


def test_random_search_refuses_a_min_max_pair():
    """A pair used to be ignored by RR in silence: the default width ran."""
    with pytest.raises(ValueError, match=r"p_dict_bounds\['T'\] = \[70, 130\] "
                                         r"is not a search width"):
        fit_data(_data(), dict(START), {'T': True}, E_in=[1, [0, 0]],
                 p_dict_bounds={'T': [70, 130]}, fit_algorithm='RR')


@pytest.mark.parametrize('entry, width', [
    (30, 30.), (30., 30.), ([30], 30.), ((30,), 30.), (np.array([30.]), 30.),
    (np.float32(30), 30.), (0, 0.)])
def test_width_forms(entry, width):
    assert search_widths({'T': entry}) == {'T': width}


@pytest.mark.parametrize('entry', [
    [70, 130], -1., np.inf, np.nan, 'abc', True, [], [[30]], None])
def test_malformed_widths_raise(entry):
    with pytest.raises(ValueError, match=r"p_dict_bounds\['T'\] = .* is not a "
                                         r"search width"):
        search_widths({'T': entry})


# --- ML and differential evolution: [min, max] pairs -------------------------

@pytest.mark.parametrize('algo', ['ML', 'differential_evolution'])
def test_ranged_fits_refuse_a_width(algo):
    with pytest.raises(ValueError, match=r"p_dict_bounds\['T'\] = \[30\] is not "
                                         r"a \[min, max\] pair"):
        fit_data(_data(), dict(START), {'T': True}, E_in=[1, [0, 0]],
                 p_dict_bounds={'T': [30]}, fit_algorithm=algo)


@pytest.mark.parametrize('key', ['T', 'Bfield'])
@pytest.mark.parametrize('entry', [[130., 70.], [70., 70.]])
def test_pair_must_increase(key, entry):
    """A reversed pair was reported as outside the physical range (T), or
    passed to lmfit, which swaps it in silence (Bfield)."""
    with pytest.raises(ValueError, match=fr"p_dict_bounds\['{key}'\].*minimum "
                                         r"must be below the maximum"):
        clip_bounds_to_limits({key: entry})


@pytest.mark.parametrize('entry', [
    [70, 130], (70., 130.), np.array([70, 130]), [np.float32(70), np.int64(130)]])
def test_pair_forms(entry):
    assert clip_bounds_to_limits({'Bfield': entry}) == {'Bfield': [70., 130.]}


def test_unbounded_pair_is_allowed():
    assert clip_bounds_to_limits({'Bfield': [-np.inf, np.inf]}) == {
        'Bfield': [-np.inf, np.inf]}


@pytest.mark.parametrize('entry', [
    30, [30], [1, 2, 3], [], ['70', 130], [True, 130], [np.nan, 130], None,
    'abc', [[70, 130]]])
def test_malformed_pairs_raise(entry):
    with pytest.raises(ValueError, match=r"p_dict_bounds\['Bfield'\] = .* is not "
                                         r"a \[min, max\] pair"):
        clip_bounds_to_limits({'Bfield': entry})
