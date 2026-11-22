"""
The dielectric solver's shortcuts.

solve_diel has closed forms for B along the beam (Faraday) and across it
(Voigt), and a general oblique solution for every other angle. Each shortcut
must be the limit of the oblique solution, so that nothing jumps at its edge:

- below 0.01 G the Faraday form used to take every angle but Voigt,
  modelling the whole field along +z or -z. The Faraday rotation of an
  oblique field came out 1/|cos Btheta| too large, and the spectra jumped at
  0.01 G. Now the geometry alone picks the form, at any field.
- the Voigt form's mode-1 vector had the wrong sign of its longitudinal
  part, so it did not solve the wave equation and flipped against the
  oblique solution at the edge of the shortcut. (Now the field leaving
  the cell is transverse, so this no longer shows in E_out.)

The expected values come from the wave equation and from a first-principles
propagation (test_geometry_conventions._reference), not from ElecSus's output.
"""
import warnings

import numpy as np
import pytest

from elecsus.libs import spectra as SP
from elecsus.libs import solve_dielectric as SD
from elecsus.tests.test_geometry_conventions import _reference

PATHS = [True, False]   # use_non_orthogonal


def _spectra(X, E_in, p, outputs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return [np.real(np.asarray(o))
                for o in SP.get_spectra(X, np.asarray(E_in, dtype=complex), p, outputs=outputs)]


# --- small fields at oblique angles ------------------------------------------

X_HOT = np.linspace(-4000, 4000, 201)
HOT = dict(Elem='Rb', Dline='D2', T=140., lcell=75e-3, Bphi=0.)
E_30 = [np.cos(np.pi / 6), np.sin(np.pi / 6), 0.]


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('theta_deg', [30, 60, 85])
def test_spectra_are_continuous_in_bfield_across_0p01_gauss(theta_deg, nonortho):
    """Before the fix, psi jumped by 1e-3 to 7e-3 rad here, and S0 changed too."""
    p = dict(HOT, Btheta=np.radians(theta_deg), use_non_orthogonal=nonortho)
    below = _spectra(X_HOT, E_30, dict(p, Bfield=0.0099999), ['S0', 'psi'])
    above = _spectra(X_HOT, E_30, dict(p, Bfield=0.0100001), ['S0', 'psi'])
    ok = below[0] > 1e-3                 # psi is undefined where nothing is transmitted
    assert np.max(abs(below[0] - above[0])) < 1e-8
    assert np.max(abs(below[1] - above[1])[ok]) < 1e-6


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('theta_deg', [60, 85])
def test_a_small_oblique_field_matches_a_first_principles_propagation(theta_deg, nonortho):
    """0.005 G: only B cos(Btheta) makes circular birefringence. Before the fix the
    code returned the Btheta = 0 result, off by 1e-3 to 1e-2 here."""
    t = np.radians(theta_deg)
    p = dict(HOT, Bfield=0.005, rb85frac=72.17)
    ref = _reference(X_HOT, p, np.asarray(E_30, dtype=complex),
                     np.array([np.sin(t), 0., np.cos(t)]))
    along_z = _reference(X_HOT, p, np.asarray(E_30, dtype=complex), np.array([0., 0., 1.]))
    code = _spectra(X_HOT, E_30, dict(p, Btheta=t, use_non_orthogonal=nonortho),
                    ['S0', 'S1', 'S2', 'S3'])
    ok = ref[0] > 1e-3                   # S1..S3 are 0/0 where nothing is transmitted
    # the old answer, B along z, is far outside the tolerance
    assert np.max(abs(along_z[2] - ref[2])[ok]) > 1e-3
    # <= 1e-9 for the non-orthogonal path; the orthogonal one assumes the
    # modes orthogonal, which at such fields they nearly are (~1e-6)
    tol = 1e-7 if nonortho else 1e-5
    for r, c in zip(ref, code):
        np.testing.assert_allclose(c[ok], r[ok], rtol=0, atol=tol)


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('theta_deg', [30, 85, 150])
def test_a_zero_field_at_an_oblique_angle(theta_deg, nonortho):
    """Bfield = 0 goes through the oblique solver too, which calc_chi feeds
    with its stand-in for zero, 1e-4 G: the limit of small fields, with
    nothing special at 0. The reference makes the same stand-in."""
    t = np.radians(theta_deg)
    p = dict(HOT, Bfield=0., rb85frac=72.17)
    ref = _reference(X_HOT, p, np.asarray(E_30, dtype=complex),
                     np.array([np.sin(t), 0., np.cos(t)]))
    code = _spectra(X_HOT, E_30, dict(p, Btheta=t, use_non_orthogonal=nonortho),
                    ['S0', 'S1', 'S2', 'S3'])
    ok = ref[0] > 1e-3
    for r, c in zip(ref, code):
        np.testing.assert_allclose(c[ok], r[ok], rtol=0, atol=1e-7)


# --- the Voigt form's mode vector ---------------------------------------------

X_V = np.linspace(-3000, 3000, 31)
P_V = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=500., Bphi=0.)


def _modes(theta, force_numeric=False):
    """The solver's two modes at Btheta = theta, and the lab-frame dielectric
    tensor built from first principles: eps = 1 + R chi R^T, R turning z to
    B = (sin theta, 0, cos theta), with calc_chi's sigma+ susceptibility on
    the (x + iy)/sqrt2 component about B (as test_geometry_conventions checks
    by the Zeeman shift)."""
    p = dict(P_V, Btheta=theta)
    cp, cm, cz = SP.calc_chi(X_V, p)
    RM, n1, n2 = SD.solve_diel(cp, cm, cz, theta, p['Bfield'], force_numeric=force_numeric)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    e_p = np.array([1, 1j, 0]) / np.sqrt(2)
    e_m = np.array([1, -1j, 0]) / np.sqrt(2)
    z = np.array([0, 0, 1])
    eps = np.array([np.eye(3) + R @ (a * np.outer(e_p, e_p.conj()) + b * np.outer(e_m, e_m.conj())
                                     + d * np.outer(z, z)) @ R.T
                    for a, b, d in zip(cp, cm, cz)])
    return RM, (n1, n2), eps


@pytest.mark.parametrize('theta', [np.pi / 2, 3 * np.pi / 2, -np.pi / 2, np.radians(70)])
def test_the_modes_solve_the_wave_equation(theta):
    """[eps - n^2 (1 - zz)] u = 0 for light along z. Before the fix the Voigt
    form's mode 1 left a residual of order 1 at pi/2 and 3 pi/2. (70 deg goes
    through the oblique solver, as a control.)"""
    RM, ns, eps = _modes(theta)
    transverse = np.diag([1., 1., 0.])
    for mode in (0, 1):
        for i in range(len(X_V)):
            u = RM[mode, :, i]
            residual = (eps[i] - ns[mode][i] ** 2 * transverse) @ u
            assert np.linalg.norm(residual) <= 1e-9 * np.linalg.norm(u), (mode, i)


@pytest.mark.parametrize('theta0', [np.pi / 2, 3 * np.pi / 2])
def test_the_voigt_mode_is_the_limit_of_the_oblique_one(theta0):
    """Mode 1's longitudinal part relative to its transverse part, u_z/u_y,
    on both sides of the shortcut (|cos Btheta| < 1e-4) and in it. Before
    the fix it changed sign at the edge."""
    RM_v, ns_v, _ = _modes(theta0)
    ratio_v = RM_v[0, 2] / RM_v[0, 1]
    for d in (-3e-4, 3e-4):
        RM_o, ns, _ = _modes(theta0 + d)
        # at each detuning, the oblique mode with the Voigt mode 1's index (the
        # other is the pi mode, along x); the solver's order can change
        mode = np.where(abs(ns[0] - ns_v[0]) < abs(ns[1] - ns_v[0]), 0, 1)
        idx = np.arange(len(X_V))
        ratio_o = RM_o[mode, 2, idx] / RM_o[mode, 1, idx]
        np.testing.assert_allclose(ratio_o, ratio_v, rtol=2e-2, atol=0)
