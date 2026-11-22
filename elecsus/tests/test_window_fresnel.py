"""
Window Fresnel factors in the non-orthogonal propagation path.

Pins the numerical claims of the window Fresnel derivation:
far from resonance (n_j -> 1) the four-interface chain
air -> glass -> vapour -> glass -> air multiplies the transmission by exactly
|t1 t2 t3 t4|^2, and n_window=None reproduces the bare (window-free) Briscoe
propagation.
"""
import numpy as np

from elecsus.elecsus_methods import calculate

# far wing: ~80 GHz detuned, vapour index -> 1 to ~1e-9
X_FAR = np.linspace(80e3, 90e3, 5)
BASE = dict(Elem='Rb', Dline='D2', T=20., lcell=75e-3, Bfield=100.,
            Btheta=np.pi / 4, Bphi=0., use_non_orthogonal=True)

NW = 1.51


def _S0(**kw):
    [S0] = calculate(X_FAR, [1, 0, 0], dict(BASE, **kw), outputs=['S0'])
    return np.asarray(S0).squeeze().real


def test_far_wing_window_factor_matches_fresnel_product():
    """S0(with windows) / S0(bare) equals |t1 t2 t3 t4|^2 with n_j = 1."""
    ratio = _S0(n_window=NW) / _S0(n_window=None)
    t_product = (2 / (1 + NW)) * (2 * NW / (NW + 1.0)) * \
                (2 * 1.0 / (1.0 + NW)) * (2 * NW / (NW + 1))
    np.testing.assert_allclose(ratio, abs(t_product) ** 2, rtol=1e-6)
    # headline number for borosilicate: ~15.5% reflection loss
    np.testing.assert_allclose(ratio, 0.844808, atol=1e-4)


def test_no_window_is_lossless_far_from_resonance():
    np.testing.assert_allclose(_S0(n_window=None), 1.0, atol=1e-6)


def test_index_matched_windows_equal_none_far_from_resonance():
    """n_window=1 removes the outer factors exactly and the inner factors
    tend to 1 as n_j -> 1; far from resonance it must agree with None.
    (Near resonance they differ at O(chi) - index-matched glass still
    presents a glass-vapour step; see the derivation note.)"""
    np.testing.assert_allclose(_S0(n_window=1.0), _S0(n_window=None), atol=1e-9)


def test_default_is_windowless():
    """[ES4 MOD] n_window defaults to None, not 1.51.

    With a 1.51 default, flipping use_non_orthogonal -- a switch about mode
    decomposition -- also switched Fresnel window losses on, moving absolute S0
    by 15.5% far from resonance, because the orthogonal path ignores n_window
    altogether. The two are independent and are now decoupled: windows are
    opt-in, so the solvers agree on absolute transmission and
    use_non_orthogonal changes only what its name says.
    """
    np.testing.assert_allclose(_S0(), _S0(n_window=None), atol=1e-12)


def test_both_solvers_agree_far_from_resonance_by_default():
    """The orthogonal and non-orthogonal paths must give the same absolute
    transmission where the mode decomposition cannot matter (n_j -> 1)."""
    [orthogonal] = calculate(X_FAR, [1, 0, 0],
                             dict(BASE, use_non_orthogonal=False), outputs=['S0'])
    np.testing.assert_allclose(_S0(), np.asarray(orthogonal).squeeze().real,
                               rtol=1e-9)
