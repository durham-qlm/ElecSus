"""
Symmetries and bounds any solution must obey (a convention-free physics test).

None of these needs a sign convention, a reference spectrum or ElecSus's own
output:

- **B -> -B is a mirror image.** A reflection sigma in a plane containing the
  wavevector k and B keeps k and the polar vector E's in-plane parts, but
  reverses the axial vector B (an axial vector's in-plane components flip
  under a reflection). Maxwell's equations and the atoms are
  reflection-invariant, so E_out(-B, sigma E_in) = sigma E_out(B, E_in)
  exactly, for either propagation model, with or without windows, and for a
  reversed pass.
- **B = 0 is isotropic.** With no field there is no preferred direction: the
  three susceptibilities are equal, every input polarisation is transmitted
  alike and keeps its polarisation, and the field angles do nothing.
- **Passivity.** Unexcited atoms only absorb, so S0 <= 1.
- **A pure state stays pure.** A coherent field remains fully polarised:
  S1^2 + S2^2 + S3^2 = 1 for the fractional Stokes parameters, wherever any
  light is left. Outside the cell the field is a plane wave along z, which
  has no z component.

The reversed pass's non-reciprocity (Faraday rotation doubles on the way
back) and Lorentz reciprocity, J_back(B) = J_forward(-B)^T, are pinned in
test_reversed_pass.py.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate
from elecsus.libs.spectra import calc_chi

X = np.linspace(-8e3, 8e3, 401)


def _random_cases(n, seed):
    """Oblique fields, both propagation models, windows, both directions,
    elements and lines, temperatures from transparent to opaque."""
    rng = np.random.default_rng(seed)
    lines = [('Rb', 'D2'), ('Rb', 'D1'), ('Cs', 'D2'), ('K', 'D1'), ('Na', 'D2')]
    cases = []
    for i in range(n):
        elem, dline = lines[i % len(lines)]
        nonortho = bool(i % 2)
        E = rng.normal(size=3) + 1j * rng.normal(size=3)
        E[2] = 0
        p = dict(Elem=elem, Dline=dline,
                 T=float(rng.uniform(20, 120 if elem != 'Na' else 200)),
                 lcell=float(rng.choice([1e-3, 5e-3, 75e-3])),
                 Bfield=float(rng.choice([0.5, 50., 500., 5000.])),
                 Btheta=float(rng.uniform(0, np.pi)),
                 Bphi=float(rng.uniform(-np.pi, np.pi)),
                 use_non_orthogonal=nonortho,
                 n_window=1.5 if (nonortho and i % 3 == 0) else None,
                 wavevector_dirn=-1 if i % 4 == 1 else 1)
        cases.append(pytest.param(p, E / np.linalg.norm(E), id=f'case{i}'))
    return cases


def _calc(E_in, p, outputs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return [np.asarray(o) for o in calculate(X, E_in, p, outputs=outputs)]


def _mirror(phi):
    """Reflection in the plane containing z and the azimuth phi."""
    c, s = np.cos(phi), np.sin(phi)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return R @ np.diag([1., -1., 1.]) @ R.T


# --- B -> -B ----------------------------------------------------------------------

@pytest.mark.parametrize('p, E_in', _random_cases(24, seed=11))
def test_reversing_B_is_the_mirror_image(p, E_in):
    sigma = _mirror(p['Bphi'])
    reversed_B = dict(p, Btheta=np.pi - p['Btheta'], Bphi=p['Bphi'] + np.pi)
    [E_out] = _calc(E_in, p, ['E_out'])
    [E_mirror] = _calc(sigma @ E_in, reversed_B, ['E_out'])
    np.testing.assert_allclose(E_mirror, sigma @ E_out, rtol=0, atol=1e-9)


@pytest.mark.parametrize('nonortho', [False, True])
@pytest.mark.parametrize('theta', [0., np.pi])
def test_reversing_a_longitudinal_B_reverses_the_rotation(theta, nonortho):
    """The mirror image's best-known case: x-polarised light along B. S0 and S1
    are unchanged, S2 and S3 change sign, so the Faraday rotation is reversed."""
    p = dict(Elem='Rb', Dline='D2', T=80., lcell=5e-3, Bfield=300., Btheta=theta,
             use_non_orthogonal=nonortho)
    outputs = ['S0', 'S1', 'S2', 'S3']
    fwd = _calc([1, 0, 0], p, outputs)
    rev = _calc([1, 0, 0], dict(p, Btheta=np.pi - theta), outputs)
    for a, b, sign in zip(fwd, rev, (1, 1, -1, -1)):
        np.testing.assert_allclose(b, sign * a, rtol=0, atol=1e-12)
    assert np.max(np.abs(fwd[2])) > 0.1           # there is a rotation to reverse


# --- B = 0 ------------------------------------------------------------------------

@pytest.mark.parametrize('elem, dline', [('Rb', 'D2'), ('Cs', 'D1'), ('K', 'D2'),
                                         ('Na', 'D1')])
def test_zero_field_susceptibilities_are_equal(elem, dline):
    """To 1e-4 of the peak. At B = 0 each hyperfine level is degenerate, the
    eigenvectors within it are an arbitrary choice, and FreqStren's cut of
    transitions with |<e|q,g>|^2 < 0.0005 (test_sum_rule.CUT) then drops
    different strength from each manifold: 3e-5 of the peak for Rb D2. With
    the cut removed the three agree to 4e-7."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        chi_plus, chi_minus, chi_pi = calc_chi(X, dict(Elem=elem, Dline=dline, T=60.))
    scale = np.max(np.abs(chi_pi))
    np.testing.assert_allclose(chi_plus, chi_pi, rtol=0, atol=1e-4 * scale)
    np.testing.assert_allclose(chi_minus, chi_pi, rtol=0, atol=1e-4 * scale)


INPUTS = [[1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 1j, 0], [1, -1j, 0], [0.6, 0.8 * np.exp(0.7j), 0]]


@pytest.mark.parametrize('nonortho', [False, True])
@pytest.mark.parametrize('theta, phi', [(0., 0.), (np.pi / 2, 0.), (1.0, 0.4), (2.5, -2.0)])
def test_zero_field_is_isotropic(theta, phi, nonortho):
    """Every input is transmitted as exp(-2 Im(n) k L) of the common index,
    whatever its polarisation and whatever the (meaningless) field angles, and
    leaves with the polarisation it had. S0 to 1e-5: the manifolds differ by
    the cut's 3e-5 (see above), 5e-7 in S0 here."""
    p = dict(Elem='Rb', Dline='D2', T=70., lcell=5e-3, Bfield=0., Btheta=theta, Bphi=phi,
             use_non_orthogonal=nonortho)
    [S0_ref] = _calc([1, 0, 0], dict(p, Btheta=0., Bphi=0.), ['S0'])
    assert S0_ref.min() < 0.5                     # the vapour absorbs here
    for E in INPUTS:
        E = np.asarray(E, dtype=complex) / np.linalg.norm(E)
        S0, E_out = _calc(E, p, ['S0', 'E_out'])
        np.testing.assert_allclose(S0, S0_ref, rtol=0, atol=1e-5)
        # same polarisation: E_out is E_in times one complex number
        overlap = np.abs(E.conj() @ E_out) ** 2 / S0
        np.testing.assert_allclose(overlap, 1, rtol=0, atol=1e-9)


# --- passivity and purity --------------------------------------------------------

CASES = _random_cases(40, seed=5)


@pytest.mark.parametrize('p, E_in', CASES)
def test_transmission_never_exceeds_one(p, E_in):
    outputs = ['S0', 'S0_P'] if p['use_non_orthogonal'] else ['S0']
    for s in _calc(E_in, p, outputs):
        assert np.all(np.isfinite(s))
        assert s.max() <= 1 + 1e-12


# Fixed oblique cases for the non-orthogonal path, on top of the random ones.
# It used to return the field inside the vapour at z = L, with its
# longitudinal component, rather than the transmitted (transverse) field, and
# S0 counted |E_z|^2: here |E_z| was ~5e-5 of the input and the
# purity was off by ~1e-7. In Faraday and Voigt geometry the modes are
# transverse, and the path always passed.
NON_ORTHOGONAL_CASES = [
    pytest.param(dict(Elem='Rb', Dline='D2', T=100., lcell=5e-3, Bfield=B, Btheta=theta,
                      Bphi=0.7, use_non_orthogonal=True, n_window=n_window,
                      wavevector_dirn=dirn),
                 np.array([0.6, 0.8j, 0]),
                 id=f'nonortho-theta{theta}-B{B:g}-window{n_window}-dirn{dirn}')
    for theta in (0.5, 2.2) for B in (500., 5000.) for n_window in (None, 1.5)
    for dirn in (1, -1)]


@pytest.mark.parametrize('p, E_in', CASES + NON_ORTHOGONAL_CASES)
def test_a_pure_state_stays_fully_polarised(p, E_in):
    S0, S1, S2, S3, E_out = _calc(E_in, p, ['S0', 'S1', 'S2', 'S3', 'E_out'])
    lit = S0 > 1e-30            # below this the Stokes ratios are 0/0 round-off
    assert lit.any()
    np.testing.assert_allclose((S1 ** 2 + S2 ** 2 + S3 ** 2)[lit], 1, rtol=0, atol=1e-10)


@pytest.mark.parametrize('p, E_in', CASES + NON_ORTHOGONAL_CASES)
def test_the_transmitted_field_is_transverse(p, E_in):
    [E_out] = _calc(E_in, p, ['E_out'])
    np.testing.assert_allclose(E_out[2], 0, rtol=0, atol=1e-12)
