"""'theta0' and 'Pol' are ES3-only polarisation keys with no effect
in ElecSus 4, where the input polarisation is the explicit E_in vector.

A p_dict ported from ES3 that still sets them must warn (so the polarisation
angle is not silently dropped), and must not change the computed spectrum.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate

_X = np.linspace(-4e3, 4e3, 20)
_BASE = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=100.)


def test_theta0_warns():
    with pytest.warns(UserWarning, match="theta0"):
        calculate(_X, [1, 0, 0], dict(_BASE, theta0=45.), outputs=['S0'])


def test_Pol_warns():
    with pytest.warns(UserWarning, match="Pol"):
        calculate(_X, [1, 0, 0], dict(_BASE, Pol=50.), outputs=['S0'])


def test_no_warning_when_absent(recwarn):
    calculate(_X, [1, 0, 0], dict(_BASE), outputs=['S0'])
    assert not any("theta0" in str(w.message) or "'Pol'" in str(w.message)
                   for w in recwarn)


def test_theta0_has_no_effect():
    """The spectrum must be unchanged with and without theta0 set: guards the
    removal contract against an accidental reimplementation.

    Compared to within a few ulp rather than exactly.  Two *identical* calls to
    calculate() are not bit-reproducible across call sites -- the same spectrum
    computed twice in one process can differ in the last bit (~2.2e-16 on an S0
    of order 1), because numpy dispatches different SIMD kernels depending on
    the alignment of the freshly allocated intermediates.  assert_array_equal
    here therefore failed for ~40% of runs for reasons having nothing to do
    with theta0; a 1e-12 relative tolerance is still ~4 orders of magnitude
    tighter than any effect a real polarisation-angle key could have.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        [s_with] = calculate(_X, [1, 0, 0], dict(_BASE, theta0=90.), outputs=['S0'])
    [s_without] = calculate(_X, [1, 0, 0], dict(_BASE), outputs=['S0'])
    np.testing.assert_allclose(s_with, s_without, rtol=1e-12, atol=0)
