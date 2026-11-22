"""
The default propagation is the non-orthogonal solver.

ElecSus 3 projected the input field onto the vapour's two modes with the
conjugate transpose, which is right only when the modes are orthogonal: B
along the beam (Faraday) or across it (Voigt). In an oblique field the modes
are not orthogonal, and the projection is wrong by up to ~0.5 in S0 (Briscoe
et al.'s Regime I cell at Btheta = 80 deg). Briscoe et al.'s linear solve,
use_non_orthogonal=True, holds in any geometry, and is the default since
4.0.0b1; use_non_orthogonal=False is ES3's method, kept for compatibility.

The expected values come from a first-principles propagation (the transverse
wave equation for the lab-frame dielectric tensor,
test_geometry_conventions._reference), never from ElecSus's own output.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate
from elecsus.tests.test_geometry_conventions import _reference

X = np.linspace(-6000, 6000, 121)
# Briscoe et al. 2025, Regime I, without the buffer-gas broadening
CELL = dict(Elem='Rb', Dline='D2', T=120., lcell=2e-3, Bfield=239., rb85frac=72.17)
E_INS = ([1., 0., 0.], [np.cos(np.radians(50)), np.sin(np.radians(50)), 0.],
         [0.6, 0.8 * np.exp(0.7j), 0.])
OUTPUTS = ['S0', 'S1', 'S2', 'S3']


def _calc(E_in, p):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return [np.real(np.asarray(o)) for o in calculate(X, E_in, p, outputs=OUTPUTS)]


def _B_hat(t, f):
    return np.array([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)])


@pytest.mark.parametrize('theta_deg, phi_deg', [(30, 0), (45, 20), (80, 0), (135, -60)])
def test_the_default_is_the_first_principles_propagation(theta_deg, phi_deg):
    """No propagation key given. Before 4.0.0b1 this was ES3's method, off by
    up to 0.48 in S0 here."""
    t, f = np.radians(theta_deg), np.radians(phi_deg)
    for E_in in E_INS:
        ref = _reference(X, CELL, np.asarray(E_in, dtype=complex), _B_hat(t, f))
        code = _calc(E_in, dict(CELL, Btheta=t, Bphi=f))
        ok = ref[0] > 1e-3                 # S1..S3 are 0/0 where nothing is transmitted
        for r, c in zip(ref, code):
            np.testing.assert_allclose(c[ok], r[ok], rtol=0, atol=1e-5)


def test_es3_propagation_is_still_there_and_is_wrong_in_oblique_fields():
    """use_non_orthogonal=False still gives ES3's method (the ES3 golden
    tests pin it), and its error is what the manual and README quote."""
    t = np.radians(80.)
    ref, = _reference(X, CELL, np.array([1, 0, 0], dtype=complex), _B_hat(t, 0.))[:1]
    S0 = _calc([1., 0., 0.], dict(CELL, Btheta=t, use_non_orthogonal=False))[0]
    assert 0.4 < np.max(abs(S0 - ref)) < 0.55


@pytest.mark.parametrize('theta_deg, x_error, worst_error', [(30, 0.03, 0.07), (45, 0.07, 0.13),
                                                             (80, 0.48, 0.53)])
def test_the_quoted_errors_of_es3s_method(theta_deg, x_error, worst_error):
    """The README and CHANGELOG quote the error of ES3's method in this cell
    (no buffer gas), to two places: for x-polarised light, and for the worst
    linear polarisation (searched in 5 degree steps)."""
    X7 = np.linspace(-7000, 7000, 1401)
    p = dict(CELL, Btheta=np.radians(theta_deg), Bphi=0.)

    def err(a):
        E_in = [np.cos(a), np.sin(a), 0.]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            new, = calculate(X7, E_in, dict(p, use_non_orthogonal=True), outputs=['S0'])
            old, = calculate(X7, E_in, dict(p, use_non_orthogonal=False), outputs=['S0'])
        return np.max(abs(np.real(new) - np.real(old)))

    assert round(err(0.), 2) == x_error
    assert round(max(err(a) for a in np.radians(np.arange(0, 180, 5))), 2) == worst_error


@pytest.mark.parametrize('theta', [0., np.pi / 2, np.pi], ids=['faraday', 'voigt', 'anti-faraday'])
@pytest.mark.parametrize('phi', [0., 0.7])
def test_faraday_and_voigt_results_do_not_change(theta, phi):
    """There the modes are orthogonal, so the two methods agree. They agree
    to round-off now that the default's field leaves the cell transverse
: the |E_z|^2 it used to keep made Voigt S0 differ by ~2e-6."""
    for E_in in E_INS:
        p = dict(CELL, Btheta=theta, Bphi=phi)
        new = _calc(E_in, dict(p, use_non_orthogonal=True))
        old = _calc(E_in, dict(p, use_non_orthogonal=False))
        ok = old[0] > 1e-3
        np.testing.assert_allclose(new[0], old[0], rtol=0, atol=1e-12)
        for a, b in zip(new[1:], old[1:]):
            np.testing.assert_allclose(a[ok], b[ok], rtol=0, atol=1e-12)
