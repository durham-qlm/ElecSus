"""
Field angle Btheta is periodic, and B along -x is a valid Voigt geometry
.

solve_diel took its Voigt shortcut only for Btheta = pi/2 (mod 2 pi). B along
-x (Btheta = 3pi/2 or -pi/2) fell through to the oblique solver, whose Ez = 1
eigenvector normalisation is singular when B is perpendicular to k, and S0
came out at ~1e10 with Stokes parameters wrong by O(1). No test, and none of
the ES3 golden cases, had Btheta above pi/2.

Reversing B in Voigt geometry swaps the sigma+ and sigma- susceptibilities
but leaves both refractive indices and the transverse mode directions
unchanged, so every Stokes spectrum must be identical. The same holds for
any angle that differs by a whole turn.
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate

X = np.linspace(-1e4, 1e4, 801)
P = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, rb85frac=72.17,
         Bfield=1000.)
STOKES = ['S0', 'S1', 'S2', 'S3']
E_IN = [0.6, 0.8 * np.exp(0.4j), 0]     # elliptical: excites both modes
PATHS = [False, True]                   # use_non_orthogonal


def _stokes(theta, nonortho, **extra):
    p = dict(P, Btheta=theta, use_non_orthogonal=nonortho, **extra)
    return np.array([np.asarray(s).real
                     for s in calculate(X, E_IN, p, outputs=STOKES)])


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('theta', [3 * np.pi / 2, -np.pi / 2, 5 * np.pi / 2,
                                   np.pi / 2 - 2 * np.pi])
def test_voigt_same_for_either_direction_of_B(theta, nonortho):
    np.testing.assert_allclose(_stokes(theta, nonortho),
                               _stokes(np.pi / 2, nonortho),
                               rtol=0, atol=1e-12)


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('theta, same_as', [
    (2 * np.pi, 0.), (-2 * np.pi, 0.), (-np.pi, np.pi), (3 * np.pi, np.pi),
])
def test_faraday_periodic(theta, same_as, nonortho):
    np.testing.assert_allclose(_stokes(theta, nonortho),
                               _stokes(same_as, nonortho),
                               rtol=0, atol=1e-12)


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('centre', [np.pi / 2, 3 * np.pi / 2])
@pytest.mark.parametrize('offset', [-1e-3, -2e-4, -1e-6, -1e-12, 1e-12,
                                    1e-9, 1e-4, 2e-4, 1e-3])
def test_no_gain_near_voigt(centre, offset, nonortho):
    """Either side of the shortcut tolerance (1e-4 rad), both directions of B.
    The physical transmission is S0 on the orthogonal path and S0_P on the
    non-orthogonal path (where |E|^2-based S0 is not bounded by 1)."""
    key = 'S0_P' if nonortho else 'S0'
    p = dict(P, Btheta=centre + offset, use_non_orthogonal=nonortho)
    transmission = np.asarray(calculate(X, E_IN, p, outputs=[key])[0]).real
    assert np.all(np.isfinite(transmission))
    assert transmission.max() <= 1 + 1e-9
