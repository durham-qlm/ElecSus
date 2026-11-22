"""
ElecSus computes the geometry its documents describe.

The documents put B along (sin t cos p, sin t sin p, cos t) for Btheta = t,
Bphi = p, with fields varying as exp(i(kz - wt)). ElecSus 3, and the ElecSus 4
alpha pre-releases, computed the mirror image (y -> -y) of that problem, through three
sign choices:

- the field was turned by +Bphi to bring B into the x-z plane, which put B at
  azimuth -Bphi;
- the [1,+i]/sqrt2 field was given the sigma- susceptibility, which it has
  only when B points along -z;
- the dielectric solver oriented B at (-sin t, 0, cos t).

Its S3 also had the opposite sign to Born & Wolf's, which the manual states.
The tests here pin the documented conventions with checks that need no
reference to the old code: the pi eigenmode of Voigt geometry, the Zeeman
shift of the lines [1,+i] light drives, a first-principles propagation of
the documented geometry , and
Briscoe et al.'s Regime I data at input angles where the mirror matters.
"""
import os
import warnings

import numpy as np
import pytest

from elecsus.libs import spectra as SP
from elecsus.libs import AtomConstants as AC

K_RB_D2 = AC.transitions['RbD2'].wavevectorMagnitude


def _spectra(X, E_in, p, outputs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return SP.get_spectra(X, np.asarray(E_in, dtype=complex), p, outputs=outputs)


# --- Bphi turns B from x towards y ---------------------------------------------

@pytest.mark.parametrize('bphi_deg', [30., -50., 120.])
def test_light_along_B_in_voigt_geometry_is_the_pi_mode(bphi_deg):
    """With B transverse, light polarised along B is a pure pi eigenmode: it
    transmits as exp(-2 Im(n_pi) k L) and keeps its polarisation. This needs
    no convention. In ES3 that was light at -Bphi."""
    X = np.linspace(-6000, 6000, 241)
    p = dict(Elem='Rb', Dline='D2', T=120., lcell=2e-3, Bfield=239., rb85frac=72.17,
             Btheta=np.pi / 2, Bphi=np.radians(bphi_deg))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        n_pi = np.sqrt(1 + SP.calc_chi(X, p)[2])
    pi_mode = np.exp(-2 * n_pi.imag * K_RB_D2 * p['lcell'])
    a = np.radians(bphi_deg)
    S0, psi = _spectra(X, [np.cos(a), np.sin(a), 0], p, ['S0', 'psi'])
    np.testing.assert_allclose(S0, pi_mode, rtol=0, atol=1e-12)
    np.testing.assert_allclose(np.cos(2 * (psi - a)), 1, rtol=0, atol=1e-12)
    [S0_mirror] = _spectra(X, [np.cos(a), -np.sin(a), 0], p, ['S0'])
    assert np.max(abs(S0_mirror - pi_mode)) > 0.1


# --- [1,+i] drives sigma+ when B points along +z -------------------------------

@pytest.mark.parametrize('btheta, sign', [(0., 1), (np.pi, -1)])
def test_left_circular_light_drives_sigma_plus_along_B(btheta, sign):
    """At 1 T the sigma+ lines of Rb-87 D2 sit ~16 GHz above the sigma- ones.
    With fields varying as exp(i(kz - wt)), [1,+i] light turns positively
    about +z, so it drives sigma+ when B is along +z, and sigma- when B is
    along -z. In ES3 it was the other way round."""
    X = np.linspace(-40000, 40000, 4001)
    p = dict(Elem='Rb', Dline='D2', Bfield=10000., T=50., lcell=2e-3, Btheta=btheta,
             rb85frac=0.)
    for E, expected in (([1, 1j, 0], sign), ([1, -1j, 0], -sign)):
        [S0] = _spectra(X, np.array(E) / np.sqrt(2), p, ['S0'])
        centroid = (X * (1 - S0)).sum() / (1 - S0).sum()
        assert np.sign(centroid) == expected and abs(centroid) > 10000


def test_chiplus_is_the_blue_shifted_manifold():
    X = np.linspace(-40000, 40000, 4001)
    p = dict(Elem='Rb', Dline='D2', Bfield=10000., T=50., rb85frac=0.)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        chi_plus, chi_minus, _ = SP.calc_chi(X, p)
    centroid = lambda chi: (X * chi.imag).sum() / chi.imag.sum()
    assert centroid(chi_plus) > 10000 > -10000 > centroid(chi_minus)


# --- Stokes S3 and Il/Ir follow Born & Wolf --------------------------------------

@pytest.mark.parametrize('E, S3, Il', [([1, 1j, 0], -1., 1.), ([1, -1j, 0], 1., 0.)])
def test_S3_is_born_and_wolf(E, S3, Il):
    """[1,+i] turns anticlockwise for an observer facing the source: left
    circular in Born & Wolf, where S3 > 0 means right circular. ES3 had this
    sign; the 4.0.0 alphas reversed it."""
    p = dict(Elem='Rb', Dline='D2', Bfield=0., T=20., lcell=1e-9)
    out = _spectra(np.array([-1e5, 1e5]), np.array(E) / np.sqrt(2), p, ['S3', 'Il', 'Ir'])
    np.testing.assert_allclose(out[0], S3, atol=1e-6)
    np.testing.assert_allclose(out[1], Il, atol=1e-6)
    np.testing.assert_allclose(out[2], 1 - Il, atol=1e-6)


# --- a first-principles propagation of the documented geometry ------------------

def _reference(X, p, E_in, B_hat):
    """Homogeneous medium, light along +z: the transverse field equation from
    eps = 1 + R (chi+ e+ e+^dag + chi- e- e-^dag + chi_pi z z^T) R^T, where R
    takes z to B and e+ = -(x + iy)/sqrt2 drives delta m = +1. Which of
    calc_chi's circular manifolds is sigma+ is decided here by its Zeeman
    shift (sigma+ is the blue-shifted one), not by its position."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        chi_a, chi_b, chi_pi = SP.calc_chi(X, p)
        # the order is by position, so a field large enough to show the shift
        # tells which position holds sigma+
        probe_a, probe_b, _ = SP.calc_chi(X, dict(p, Bfield=max(p['Bfield'], 100.)))
    centroid = lambda chi: (X * chi.imag).sum() / chi.imag.sum()
    chi_plus, chi_minus = ((chi_a, chi_b) if centroid(probe_a) > centroid(probe_b)
                           else (chi_b, chi_a))
    z = np.array([0., 0., 1.])
    v, c = np.cross(z, B_hat), z @ B_hat
    if np.linalg.norm(v) < 1e-12:
        R = np.eye(3) if c > 0 else np.diag([1., -1., -1.])
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx * ((1 - c) / np.linalg.norm(v) ** 2)
    e_plus = -np.array([1, 1j, 0]) / np.sqrt(2)
    e_minus = np.array([1, -1j, 0]) / np.sqrt(2)
    E = np.zeros((2, len(X)), complex)
    for i in range(len(X)):
        chi = (chi_plus[i] * np.outer(e_plus, e_plus.conj())
               + chi_minus[i] * np.outer(e_minus, e_minus.conj())
               + chi_pi[i] * np.outer(z, z))
        eps = np.eye(3) + R @ chi @ R.T
        eff = eps[:2, :2] - np.outer(eps[:2, 2], eps[2, :2]) / eps[2, 2]
        w, V = np.linalg.eig(eff)
        E[:, i] = V @ (np.exp(1j * np.sqrt(w) * K_RB_D2 * p['lcell'])
                       * np.linalg.solve(V, np.asarray(E_in)[:2]))
    I = abs(E[0]) ** 2 + abs(E[1]) ** 2
    return (I, (abs(E[0]) ** 2 - abs(E[1]) ** 2) / I, 2 * np.real(E[0] * E[1].conj()) / I,
            -2 * np.imag(E[0].conj() * E[1]) / I)       # Born & Wolf's S3


@pytest.mark.parametrize('theta_deg, phi_deg', [
    (0, 0), (180, 0), (30, 0), (60, 30), (80, -50), (90, 30), (120, 75), (150, -120)])
def test_matches_a_first_principles_propagation(theta_deg, phi_deg):
    """ES3 matched this only with B mirrored."""
    X = np.linspace(-6000, 6000, 121)
    p = dict(Elem='Rb', Dline='D2', T=120., lcell=2e-3, Bfield=239., rb85frac=72.17)
    t, f = np.radians(theta_deg), np.radians(phi_deg)
    B_hat = np.array([np.sin(t) * np.cos(f), np.sin(t) * np.sin(f), np.cos(t)])
    for E_in in ([np.cos(np.pi / 6), np.sin(np.pi / 6), 0], [0.6, 0.8 * np.exp(0.7j), 0]):
        ref = _reference(X, p, np.asarray(E_in, dtype=complex), B_hat)
        code = _spectra(X, E_in, dict(p, Btheta=t, Bphi=f, use_non_orthogonal=True),
                        ['S0', 'S1', 'S2', 'S3'])
        ok = ref[0] > 1e-3                 # S1..S3 are 0/0 where nothing is transmitted
        for r, c in zip(ref, code):
            np.testing.assert_allclose(c[ok], r[ok], rtol=0, atol=1e-5)


@pytest.mark.parametrize('theta', [0.3, np.pi - 0.3])
def test_near_zero_field_takes_the_nearer_axis(theta):
    """Below 0.01 G the solver took B along z. It took -z at every oblique
    angle, which, with B modelled reversed, meant +z; the geometry fix made it take
    whichever of +z and -z is nearer. The wrong one gets the (tiny) rotation
    backwards. The oblique solver now takes such fields
    (test_solver_shortcuts has steeper angles, where B along z is wrong too)."""
    X = np.linspace(-6000, 6000, 121)
    p = dict(Elem='Rb', Dline='D2', T=120., lcell=5e-3, Bfield=0.005, rb85frac=72.17)
    E_in = np.array([1, 0, 0], dtype=complex)
    ref = _reference(X, p, E_in, np.array([np.sin(theta), 0., np.cos(theta)]))
    S2, = _spectra(X, E_in, dict(p, Btheta=theta), ['S2'])
    ok = ref[0] > 1e-3
    assert np.max(abs(ref[2][ok])) > 4e-4              # a rotation to get right
    np.testing.assert_allclose(S2[ok], ref[2][ok], rtol=0, atol=1e-4)


# --- an experiment where the mirror matters --------------------------------------

BRISCOE_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'paper', 'data',
                           'briscoe2025', 'Regime I')


def _briscoe_rms(fn, E_in):
    d = np.loadtxt(os.path.join(BRISCOE_DIR, fn), delimiter=',')
    x, y = d[:, 0], d[:, 1]
    keep = (x >= -6.25) & (x <= 7.75)
    x, y = x[keep], y[keep]
    edges = np.linspace(x.min(), x.max(), 401)
    idx = np.clip(np.digitize(x, edges) - 1, 0, 399)
    count = np.bincount(idx, minlength=400)
    ok = count > 0
    xb = np.bincount(idx, x, minlength=400)[ok] / count[ok]
    yb = np.bincount(idx, y, minlength=400)[ok] / count[ok]
    p = dict(Elem='Rb', Dline='D2', T=119.69, lcell=2e-3, rb85frac=72.17,
             Btheta=np.radians(80.), Bphi=0., Bfield=239., GammaBuf=7.9,
             use_non_orthogonal=True)
    [S0] = _spectra(xb * 1e3, E_in, p, ['S0'])
    return 100 * np.sqrt(np.mean((yb - S0) ** 2))


@pytest.mark.skipif(not os.path.isdir(BRISCOE_DIR), reason='needs the repository checkout')
@pytest.mark.parametrize('fn, theta_E', [
    ('tau_thetaB80_thetaE50_1.csv', 50.), ('tau_thetaB80_thetaEminus50_2.csv', -50.),
    ('tau_thetaB80_thetaE20_2.csv', 20.)])
def test_briscoe_regime_one_data_take_the_natural_input_angle(fn, theta_E):
    """Briscoe et al.'s transmission at theta_B = 80 deg is reproduced to
    ~1 % RMS with E_in = [cos theta_E, sin theta_E], theta_E measured from x
    towards y. The mirror input [cos, -sin] gives 7-8 % here; it is the one
    that 4.0.0a1, which had the mirror image and the non-orthogonal modes,
    needed. (ES3, with the mirror image and the orthogonal projection, gives
    0.9-1.2 % with the natural input and 7-8 % with the mirror one.)"""
    t = np.radians(theta_E)
    assert _briscoe_rms(fn, [np.cos(t), np.sin(t), 0]) < 1.5
    assert _briscoe_rms(fn, [np.cos(t), -np.sin(t), 0]) > 4.
