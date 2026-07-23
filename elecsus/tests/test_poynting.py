"""
Poynting-flux transmission S0_P for the non-orthogonal path (issue #28).

Background: |E|^2 is not proportional to transported power inside an
anisotropic medium, and near mode degeneracies the biorthogonal amplitudes
grow large, so S0 = |E_out|^2/|E_in|^2 is not guaranteed <= 1. (The
catastrophic S0 ~ 41 reported in issue #28 predates the 2026-06-28
eigenvector fix and no longer reproduces; parameter scans below pin that.)
S0_P is the time-averaged Poynting-vector ratio <S_z(L)>/<S_z(0)> of the
two-mode field, which passivity bounds by 1 exactly.
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate

X = np.linspace(-4e3, 4e3, 500)

# Briscoe 2025 Regime I: near-Voigt, strong mode non-orthogonality
P_REGIME1 = dict(Elem='Rb', Dline='D2', T=120., lcell=2e-3, Bfield=240.,
                 Btheta=np.radians(84.), Bphi=0.,
                 use_non_orthogonal=True, n_window=None)


def test_regime1_poynting_bounded_and_consistent():
    [S0, S0_P] = calculate(X, [1, 0, 0], P_REGIME1, outputs=['S0', 'S0_P'])
    S0, S0_P = np.asarray(S0).real, np.asarray(S0_P).real

    assert np.all(np.isfinite(S0_P))
    assert np.all(S0_P <= 1 + 1e-9), 'passivity bound violated'
    assert np.all(S0_P >= -1e-12)
    # where the field picture is physical the two definitions agree closely
    assert np.max(np.abs(S0 - S0_P)) < 0.05


@pytest.mark.parametrize('B, theta_deg', [
    (0.5, 84.), (1., 45.), (2., 89.5), (5., 84.), (240., 84.),
])
def test_small_field_oblique_passivity(B, theta_deg):
    """The issue-#28 danger zone (small B + oblique angle): S0_P must respect
    passivity, and post-eigenvector-fix S0 itself must stay physical."""
    p = dict(Elem='Rb', Dline='D2', T=100., lcell=75e-3, Bfield=B,
             Btheta=np.radians(theta_deg), Bphi=0.,
             use_non_orthogonal=True, n_window=None)
    [S0, S0_P] = calculate(X, [1, 0, 0], p, outputs=['S0', 'S0_P'])
    S0, S0_P = np.asarray(S0).real, np.asarray(S0_P).real

    assert np.all(np.isfinite(S0_P))
    assert np.all(S0_P <= 1 + 1e-9)
    # regression pin: the catastrophic |E|^2 blow-up stays fixed
    assert S0.max() < 1.001, f'S0 exceeded 1: {S0.max():.4f} (issue #28 regression)'


def test_faraday_poynting_equals_intensity():
    """In Faraday geometry the modes are orthogonal circular polarisations and
    n ~ 1, so the flux and |E|^2 definitions must agree."""
    p = dict(Elem='Rb', Dline='D2', T=60., lcell=75e-3, Bfield=100.,
             Btheta=0., Bphi=0., use_non_orthogonal=True, n_window=None)
    [S0, S0_P] = calculate(X, [1, 0, 0], p, outputs=['S0', 'S0_P'])
    np.testing.assert_allclose(np.asarray(S0).real, np.asarray(S0_P).real,
                               atol=1e-3)


def test_poynting_requires_non_orthogonal_path():
    p = dict(Elem='Rb', Dline='D2', T=60., lcell=75e-3, Bfield=100.,
             use_non_orthogonal=False)
    with pytest.raises(ValueError, match='use_non_orthogonal'):
        calculate(X, [1, 0, 0], p, outputs=['S0_P'])
