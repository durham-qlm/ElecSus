"""
Theory reproduction of Briscoe et al., Opt. Express 33, 32099-32114 (2025).

Reproduces Figures 3 and 5 of that paper using our implementation:
  Fig 3 — Regime I:  Rb D2, T=120°C, B=240G,  θ_B=84°,  θ_E=0°  (horizontal linear)
  Fig 5 — Regime II: Rb D2, T=80°C,  B=2500G, θ_B=130°, θ_E=90° (vertical linear)

No experimental data required.  Both figures are pure theory comparisons between
ElecSus 3 (orthogonal propagation, old) and ElecSus 4 (non-orthogonal, new).

Notes:
  - Briscoe's formulation does not include Fresnel transmission at cell windows
    (their Eq. 9 is pure propagation with no interface corrections).  We set
    n_window=None here to reproduce their figures exactly.  Fresnel corrections
    are an additional feature of our code enabled by supplying a numeric n_window.
  - Btheta is in radians throughout ElecSus.
  - The eigenvectors in the oblique branch are computed by solving the 2×2
    characteristic matrix system (rows 0 and 1 of the DielMat wave equation)
    for each eigenvalue.  This replaces the Webber Appendix-B analytic formulas
    which gave incorrect eigenvectors at near-Voigt angles (large residuals).

Run:
    python -m pytest elecsus/tests/briscoe_2025_theory_repro.py -v   (CI checks)
    python elecsus/tests/briscoe_2025_theory_repro.py                 (produces plots)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from elecsus.elecsus_methods import calculate
from elecsus.libs.spectra import calc_chi
from elecsus.libs.solve_dielectric import solve_diel


# ── Helpers ───────────────────────────────────────────────────────────────────

def mode_overlap(X_MHz, p_dict):
    """
    Compute |<m1|m2>| at each detuning — the degree of non-orthogonality.

    Returns an (N,) real array in [0, 1].  A value of 0 means fully orthogonal
    (Faraday or Voigt geometry); ~1 means nearly degenerate modes.
    """
    Chi = calc_chi(X_MHz, p_dict)
    ChiPlus, ChiMinus, ChiZ = Chi
    Btheta = p_dict.get('Btheta', 0.0)
    Bfield = p_dict.get('Bfield', 0.0)

    RM_ary, _, _ = solve_diel(ChiPlus, ChiMinus, ChiZ, Btheta, Bfield)

    # RM_ary[i, :, k] is the i-th eigenvector at frequency index k.
    # Normalised overlap: |<m̂1|m̂2>| where m̂i = mi / |mi|
    inner = np.abs(np.einsum('jk,jk->k', RM_ary[0].conj(), RM_ary[1]))
    norm1 = np.linalg.norm(RM_ary[0], axis=0)
    norm2 = np.linalg.norm(RM_ary[1], axis=0)
    return inner / (norm1 * norm2)


def compute_spectra(X_MHz, E_in, p_dict):
    """Return (tau_ES3, tau_ES4, overlap) arrays."""
    p_es3 = dict(p_dict, use_non_orthogonal=False)
    # n_window=None: skip Fresnel interface corrections to match Briscoe et al. Eq. (9)
    p_es4 = dict(p_dict, use_non_orthogonal=True, n_window=None)

    tau_es3 = calculate(X_MHz, E_in, p_es3, outputs=['S0'])[0].real
    tau_es4 = calculate(X_MHz, E_in, p_es4, outputs=['S0'])[0].real
    overlap  = mode_overlap(X_MHz, p_dict)

    return tau_es3, tau_es4, overlap


def make_figure(X_GHz, E_in, p_dict, title='', outfile=None):
    """
    Two-panel plot reproducing the structure of Briscoe et al. Figs 3 / 5:
      top    — transmission τ, ElecSus 4 (blue solid) vs ElecSus 3 (black dashed)
      bottom — mode overlap |<m1|m2>| (red)
    """
    X_MHz = X_GHz * 1e3
    tau_es3, tau_es4, overlap = compute_spectra(X_MHz, E_in, p_dict)

    fig, (ax_t, ax_o) = plt.subplots(
        2, 1, figsize=(8, 5), sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})

    ax_t.plot(X_GHz, tau_es4, color='#1f77b4', lw=1.5,
              label='ElecSus 4 (non-orthogonal)')
    ax_t.plot(X_GHz, tau_es3, 'k--', lw=1.0, dashes=(4, 2),
              label='ElecSus 3 (orthogonal)')
    ax_t.set_ylabel('Transmission (τ)')
    ax_t.set_ylim(-0.05, 1.1)
    ax_t.legend(fontsize=9, loc='upper right')
    if title:
        ax_t.set_title(title, fontsize=9)

    ax_o.plot(X_GHz, overlap, color='#d62728', lw=1.2)
    ax_o.set_ylabel('|⟨m₁|m₂⟩|', fontsize=9)
    ax_o.set_xlabel('Linear Detuning (GHz)')
    ax_o.set_ylim(0, 1.05)
    ax_o.axhline(0, color='k', lw=0.5, ls=':')

    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f"Saved {outfile}")
    return fig, tau_es3, tau_es4, overlap


# ── Figure 3 — Regime I ───────────────────────────────────────────────────────
# Rb D2, L=2mm, T=120°C, B=240G, θ_B=84°, E_in horizontal (θ_E=0°)
P_REG1 = {
    'Elem':   'Rb',
    'Dline':  'D2',
    'T':      120.0,
    'lcell':  2e-3,
    'Bfield': 240.0,
    'Btheta': np.deg2rad(84.0),
    'Bphi':   0.0,
}
E_REG1 = np.array([1, 0, 0])   # θ_E = 0° horizontal linear
X_REG1 = np.linspace(-7, 7, 2000)  # GHz


# ── Figure 5 — Regime II ──────────────────────────────────────────────────────
# Rb D2, L=2mm, T=80°C, B=2500G, θ_B=130°, E_in vertical (θ_E=90°)
P_REG2 = {
    'Elem':   'Rb',
    'Dline':  'D2',
    'T':      80.0,
    'lcell':  2e-3,
    'Bfield': 2500.0,
    'Btheta': np.deg2rad(130.0),
    'Bphi':   0.0,
}
E_REG2 = np.array([0, 1, 0])   # θ_E = 90° vertical linear
X_REG2 = np.linspace(-15, 15, 3000)  # GHz


# ── Pytest-compatible tests (CI) ──────────────────────────────────────────────

def test_regime1_es3_es4_differ_significantly():
    """
    Regime I: for near-Voigt (θ_B=84°) Rb D2, ES3 and ES4 give different
    spectra because ES4 uses the correct biorthogonal mode decomposition while
    ES3 uses conjugate-transpose (which is wrong for non-orthogonal modes).

    With correct eigenvectors the maximum difference is ~0.5 in transmission,
    and both models give physically sensible results (0 ≤ τ ≤ 1).
    """
    X_MHz = X_REG1 * 1e3
    tau_es3, tau_es4, _ = compute_spectra(X_MHz, E_REG1, P_REG1)

    max_diff = np.max(np.abs(tau_es3 - tau_es4))
    assert max_diff > 0.1, \
        f"ES3 and ES4 should differ meaningfully at oblique angles; got {max_diff:.4f}"

    # ES4 should give physically reasonable transmission values
    assert tau_es4.max() <= 1.05, \
        f"ES4 tau should not significantly exceed 1.0; got {tau_es4.max():.4f}"
    assert tau_es4.min() >= -0.05, \
        f"ES4 tau should not be significantly negative; got {tau_es4.min():.4f}"


def test_regime1_mode_overlap_peak():
    """Mode overlap |<m1|m2>| should peak above 0.5 near the artefact frequency."""
    X_MHz = X_REG1 * 1e3
    overlap = mode_overlap(X_MHz, P_REG1)
    assert overlap.max() > 0.5, \
        f"Expected peak mode overlap > 0.5, got {overlap.max():.3f}"


def test_regime2_large_disagreement():
    """
    Regime II: ES3 and ES4 should disagree significantly over a wide frequency
    range (|tau_ES3 - tau_ES4| > 0.2 for many points).
    """
    X_MHz = X_REG2 * 1e3
    tau_es3, tau_es4, _ = compute_spectra(X_MHz, E_REG2, P_REG2)
    n_disagreeing = np.sum(np.abs(tau_es3 - tau_es4) > 0.2)
    assert n_disagreeing > 50, \
        f"Expected >50 points with >20% disagreement in Regime II, got {n_disagreeing}"


def test_voigt_models_agree():
    """At θ_B = 90° (Voigt), modes are orthogonal — ES3 and ES4 must agree."""
    p_voigt = dict(P_REG1, Btheta=np.pi / 2, Bfield=200.0)
    X_MHz = np.linspace(-5, 5, 500) * 1e3
    tau_es3, tau_es4, overlap = compute_spectra(X_MHz, E_REG1, p_voigt)

    assert overlap.max() < 0.01, \
        f"Voigt modes should be orthogonal (overlap < 0.01), got {overlap.max():.4f}"
    np.testing.assert_allclose(
        tau_es3, tau_es4, atol=1e-3,
        err_msg="ES3 and ES4 must agree in the Voigt geometry")


def test_faraday_models_agree():
    """At θ_B = 0° (Faraday), modes are circular and orthogonal — ES3 and ES4 agree."""
    p_faraday = dict(P_REG1, Btheta=0.0, Bfield=200.0)
    X_MHz = np.linspace(-5, 5, 500) * 1e3
    tau_es3, tau_es4, overlap = compute_spectra(X_MHz, E_REG1, p_faraday)

    assert overlap.max() < 0.01, \
        f"Faraday modes should be orthogonal (overlap < 0.01), got {overlap.max():.4f}"
    np.testing.assert_allclose(
        tau_es3, tau_es4, atol=1e-3,
        err_msg="ES3 and ES4 must agree in the Faraday geometry")


# ── Standalone plot generation ────────────────────────────────────────────────

if __name__ == '__main__':
    out_dir = os.path.dirname(__file__)

    print("Computing Regime I (Fig 3 reproduction)...")
    make_figure(
        X_REG1, E_REG1, P_REG1,
        title=(r'Regime I: Rb D$_2$, $T$=120°C, $B$=240 G, '
               r'$\theta_B$=84°, $\theta_E$=0°'
               '\n(reproduces Briscoe et al. 2025, Fig. 3)'),
        outfile=os.path.join(out_dir, 'briscoe_fig3_repro.pdf'))

    print("Computing Regime II (Fig 5 reproduction)...")
    make_figure(
        X_REG2, E_REG2, P_REG2,
        title=(r'Regime II: Rb D$_2$, $T$=80°C, $B$=2500 G, '
               r'$\theta_B$=130°, $\theta_E$=90°'
               '\n(reproduces Briscoe et al. 2025, Fig. 5)'),
        outfile=os.path.join(out_dir, 'briscoe_fig5_repro.pdf'))

    plt.show()
