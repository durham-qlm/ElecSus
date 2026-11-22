"""
Small API fixes.

* A single detuning is a one-point axis. It used to crash: calc_chi took
  len() of it before get_spectra's scalar branch ran, and that branch missed
  numpy scalars anyway.
* calculate_cascade takes strict, as calculate does.
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, calculate_cascade

P = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=100.)
X = np.linspace(-1e3, 1e3, 5)


@pytest.mark.parametrize('detuning', [5.0, 5, np.float64(5.0)])
def test_a_single_detuning_is_a_one_point_axis(detuning):
    [s] = calculate(detuning, [1, 0, 0], P, outputs=['S0'])
    [reference] = calculate([5.0], [1, 0, 0], P, outputs=['S0'])
    assert np.shape(s) == (1,)
    np.testing.assert_array_equal(s, reference)


def test_a_single_detuning_in_a_cascade():
    (s1,), (s2,) = calculate_cascade(5.0, [1, 0, 0], P, P, outputs=['S0'])
    (r1,), (r2,) = calculate_cascade([5.0], [1, 0, 0], P, P, outputs=['S0'])
    np.testing.assert_array_equal([s1, s2], [r1, r2])


def test_cascade_strict_false_downgrades_an_unknown_key():
    mistyped = dict(P, n_windows=1.51)
    with pytest.raises(ValueError, match='n_windows'):
        calculate_cascade(X, [1, 0, 0], P, mistyped, outputs=['S0'])
    with pytest.warns(UserWarning, match='n_windows') as record:
        calculate_cascade(X, [1, 0, 0], P, mistyped, outputs=['S0'],
                          strict=False)
    assert {w.filename for w in record} == {__file__}


def test_the_default_outputs_are_the_documented_eleven():
    """[ES4] (an untested path found by the coverage report) With no
    outputs argument, get_spectra returns the eleven outputs its docstring
    lists, in that order; the GUI relies on it."""
    X = np.linspace(-3000, 3000, 31)
    p = dict(Elem='Rb', Dline='D2', T=60., Bfield=200., Btheta=0.7, Bphi=0.2)
    names = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45',
             'alphaPlus', 'alphaMinus', 'alphaZ']
    default = calculate(X, [1, 0.5j, 0], p)
    named = calculate(X, [1, 0.5j, 0], p, outputs=names)
    assert len(default) == len(names)
    for d, n in zip(default, named):
        np.testing.assert_array_equal(d, np.asarray(n).real)
