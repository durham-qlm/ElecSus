"""
Random restart keeps to the physical range, and survives a fit that leaves it
.

RR_fit scatters its starts uniformly by +- a width around the starting value,
and runs a downhill (Marquardt-Levenberg) fit from each, in a process pool.

- Nothing kept the starts inside the physical range: from GammaBuf = 0, half
  the starts had GammaBuf below 0, as did lcell below its 1 mm width and a
  pure-isotope cell's fraction. The first fit from such a start raised.
- A fit from a valid start can still step outside the range: leastsq tried
  GammaBuf = -474 MHz in ES3's fitting demo (test_fitting).
- Either error stopped the whole search, since Pool.map_async raises the
  first worker exception once every start has run.

Each start is now drawn from the part of start +- width inside the range
(param_validation.FIT_LIMITS), and a start whose fit leaves it is dropped
with a warning. The pool is a stand-in here: one that stops RR_fit once it
holds the starts, or one that runs them in this process.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs import RRFittingRoutine as RR
from elecsus.libs.param_validation import FIT_LIMITS, ParameterRangeError, check_ranges

X = np.linspace(-4000, 4000, 120)
TRUE = dict(Elem='Rb', Dline='D2', T=40., lcell=75e-3, rb85frac=72.17, Bfield=0.)
# every parameter of the fit function, as fit_data fills them in
START = dict(TRUE, Btheta=0., Bphi=0., GammaBuf=0., shift=0.)


def _data(**extra):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return [X, np.real(calculate(X, [1, 0, 0], dict(TRUE, **extra), outputs=['S0'])[0])]


class _Stop(Exception):
    """Raised by the capturing pool once it holds the starts."""


class _CapturingPool:
    """Stands in for multiprocessing.Pool: keeps each start's p_dict and
    stops RR_fit there."""
    starts = []

    def __init__(self, *args, **kwargs):
        pass

    def map_async(self, func, args_list):
        type(self).starts = [args[2] for args in args_list]
        raise _Stop


class _SerialPool:
    """Stands in for multiprocessing.Pool: runs the starts one after another,
    in this process."""

    def __init__(self, *args, **kwargs):
        pass

    def map_async(self, func, args_list):
        results = [func(args) for args in args_list]

        class _Async:
            def get(self):
                return results
        return _Async()

    def close(self):
        pass

    def join(self):
        pass


def _starts(monkeypatch, p_dict, p_dict_bools, p_dict_bounds=None, n=400):
    monkeypatch.setattr(RR, 'Pool', _CapturingPool)
    np.random.seed(0)
    with pytest.raises(_Stop):
        RR.RR_fit(_data(), [1, [0, 0]], p_dict, p_dict_bools,
                  p_dict_bounds=p_dict_bounds, no_evals=n)
    return _CapturingPool.starts


# --- the starts ------------------------------------------------------------------

@pytest.mark.parametrize('key, value, low, high', [
    ('GammaBuf', 0., 0., 20.),         # the reported failure
    ('lcell', 0.5e-3, FIT_LIMITS['lcell'][0], 1.5e-3),
    ('rb85frac', 99.5, 98.5, 100.),
    ('rb85frac', 100., 99., 100.),     # a pure 85Rb cell
    ('rb85frac', 0., 0., 1.),        # a pure 87Rb cell
    ('BufferPressure', 5., 0., 15.),   # had no default width: KeyError
])
def test_starts_lie_inside_the_physical_range(monkeypatch, key, value, low, high):
    """Each drawn with RR's default width, which reaches past the range.
    The starts fill the part inside it, rather than piling up on the limit."""
    starts = np.array([p[key] for p in _starts(monkeypatch, dict(START, **{key: value}),
                                               {key: True})])
    assert np.all((starts >= low) & (starts <= high))
    assert np.ptp(starts) > 0.95 * (high - low)
    assert len(set(starts)) == len(starts)
    check_ranges({key: starts.min()})
    check_ranges({key: starts.max()})


def test_starts_away_from_the_limits_are_drawn_as_before(monkeypatch):
    """Uniform on [start - width, start + width], from the same random
    numbers, so a seeded search starts where it did."""
    widths = {'T': 30., 'Bfield': 200.}
    starts = _starts(monkeypatch, dict(START), {'T': True, 'Bfield': True}, widths, n=50)
    np.random.seed(0)
    for key in ('T', 'Bfield'):                  # the order RR_fit draws them in
        before = [START[key] + np.random.uniform(-1, 1) * widths[key] for _ in starts]
        assert [p[key] for p in starts] == pytest.approx(before, rel=1e-12, abs=1e-12)


def test_a_start_further_outside_the_range_than_its_width_raises(monkeypatch):
    """fit_data refuses such a start before it reaches RR_fit; called
    directly, RR_fit would have drawn from the reversed interval."""
    monkeypatch.setattr(RR, 'Pool', _CapturingPool)
    with pytest.raises(ParameterRangeError, match=r"'GammaBuf' = -50.0 lies further than "
                                                  r"its search width 20.0"):
        RR.RR_fit(_data(), [1, [0, 0]], dict(START, GammaBuf=-50.), {'GammaBuf': True},
                  no_evals=4)


# --- fits that leave the range -----------------------------------------------------

class _Fit:
    def __init__(self, redchi):
        self.redchi = redchi


def test_a_start_whose_fit_leaves_the_range_is_dropped(monkeypatch):
    """ML_fit is a stand-in: from a start with GammaBuf above 10 MHz the fit
    steps below 0 and raises as the real one does; from the rest it
    converges, with a cost that grows away from 5 MHz."""
    monkeypatch.setattr(RR, 'Pool', _SerialPool)
    seen = []

    def downhill(data, E_in, p_dict, p_dict_bools, data_type='S0', **kw):
        seen.append(p_dict['GammaBuf'])
        if p_dict['GammaBuf'] > 10.:
            try:
                check_ranges({'GammaBuf': -474.})
            except ParameterRangeError as cause:
                raise ParameterRangeError(f"The fit (leastsq) stopped at a trial value "
                                          f"outside the physical range: {cause}") from cause
        return dict(p_dict), _Fit(abs(p_dict['GammaBuf'] - 5.))

    monkeypatch.setattr(ML, 'ML_fit', downhill)
    np.random.seed(0)
    with pytest.warns(UserWarning, match=r"of the 16 random-restart starts failed and "
                                         r"were dropped") as record:
        best, _ = RR.RR_fit(_data(), [1, [0, 0]], dict(START, GammaBuf=10.),
                            {'GammaBuf': True}, no_evals=16)
    starts, final = np.array(seen[:16]), seen[16]
    failed = int(np.sum(starts > 10.))
    assert 0 < failed < 16
    [warning] = record
    assert str(warning.message).startswith(f"{failed} of the 16")
    assert "The first failure: 'GammaBuf' = -474.0 must not be negative" in str(warning.message)
    assert warning.filename == __file__              # the caller's line
    kept = starts[starts <= 10.]
    assert final == best['GammaBuf'] == kept[np.argmin(abs(kept - 5.))]


def test_every_start_failing_raises(monkeypatch):
    """A real search: the data's lines are narrower than the model can make
    them at T = 40 C (they come from a Doppler temperature of -150 C), so
    every downhill fit takes GammaBuf below 0 to narrow them."""
    monkeypatch.setattr(RR, 'Pool', _SerialPool)
    np.random.seed(0)
    with pytest.raises(ParameterRangeError, match=r"Every one of the 4 random-restart "
                                                  r"starts failed.*The first failure: "
                                                  r"'GammaBuf' = -[0-9.]+ must not be "
                                                  r"negative"):
        RR.RR_fit(_data(Constrain=False, DoppTemp=-150.), [1, [0, 0]], dict(START),
                  {'GammaBuf': True}, no_evals=4)


def test_a_failed_start_is_reported_not_raised():
    """evaluate() runs in the worker processes; it used to raise."""
    data = _data(Constrain=False, DoppTemp=-150.)
    p = dict(START, GammaBuf=5., E_x=1., E_y=0., E_phase=0.)
    cost, best, error = RR.evaluate((data, [1, [0, 0]], p, {'GammaBuf': True}, 'S0'))
    assert cost == np.inf and best is p
    assert error.startswith("'GammaBuf' = -") and 'must not be negative' in error
