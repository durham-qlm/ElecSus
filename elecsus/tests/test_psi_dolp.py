"""
Tests for the Stokes-based generalised rotation angle 'psi' and degree of
linear polarisation 'DoLP' outputs (added because the dispersion-only 'phi'
output keeps reporting a rotation angle in regimes where circular dichroism
has made the light circular and no linear rotation is measurable).

Pinned properties:

1. Weak-dichroism agreement: psi == phi (mod pi) wherever DoLP ~ 1, since
   with negligible differential absorption the Stokes azimuth IS the
   dispersive Faraday angle.

2. Coherence identity: the propagated fields are fully polarised, so
   DoLP^2 + S3^2 = 1 at machine precision, everywhere, in any geometry.

3. Dichroic collapse (the motivating physics): at high optical depth in the
   HPB (hyperfine Paschen-Back) regime DoLP -> 0
   at the absorption clusters while |phi| grows without bound -- psi/DoLP
   flag exactly where phi stops being measurable.

4. psi is the principal-value azimuth in (-pi/2, pi/2].
"""
import numpy as np

from elecsus.elecsus_methods import calculate

FARADAY = dict(Btheta=0.0, Bphi=0.0)


def _wrap_to_half_pi(angle):
    """Wrap an angle into (-pi/2, pi/2] (mod-pi principal value)."""
    return np.mod(angle + np.pi / 2, np.pi) - np.pi / 2


def test_psi_equals_phi_mod_pi_at_weak_dichroism():
    X = np.linspace(-15000, 15000, 1500)
    p = dict(Elem='Rb', Dline='D2', T=60., lcell=75e-3, Bfield=300., **FARADAY)
    phi, psi, DoLP = calculate(X, [1, 0, 0], p, outputs=['phi', 'psi', 'DoLP'])
    phi = np.real(phi)

    mask = DoLP > 0.995
    assert mask.sum() > 500, 'not enough weak-dichroism points to test'
    diff = np.abs(psi[mask] - _wrap_to_half_pi(phi[mask]))
    diff = np.minimum(diff, np.pi - diff)  # allow the +-pi/2 branch point
    assert diff.max() < 1e-8


def test_dolp_s3_coherence_identity():
    X = np.linspace(-10000, 10000, 800)
    # oblique geometry too -- psi/DoLP are field-derived, valid anywhere
    for geom in (FARADAY, dict(Btheta=np.pi / 3, Bphi=0.0)):
        p = dict(Elem='Rb', Dline='D2', T=100., lcell=75e-3, Bfield=1000., **geom)
        S3, DoLP = calculate(X, [1, 0, 0], p, outputs=['S3', 'DoLP'])
        np.testing.assert_allclose(DoLP**2 + np.real(S3)**2, 1.0, atol=1e-12)


def test_dolp_collapses_where_phi_keeps_growing():
    """HPB regime (Rb-87, 1.6 T, sigma-minus cluster): at high T the light is
    fully circularised at the absorption clusters (DoLP -> 0) while the
    dispersion-only phi reports ever-larger 'rotation'."""
    X = np.linspace(-28000, -16000, 1500)
    p_lo = dict(Elem='Rb', Dline='D2', T=92.1, lcell=1e-3, Bfield=16000.,
                rb85frac=0., **FARADAY)
    p_hi = dict(p_lo, T=132.0)

    phi_lo, DoLP_lo = calculate(X, [1, 0, 0], p_lo, outputs=['phi', 'DoLP'])
    phi_hi, DoLP_hi = calculate(X, [1, 0, 0], p_hi, outputs=['phi', 'DoLP'])

    # moderate T: rotation is genuinely measurable everywhere
    assert DoLP_lo.min() > 0.5
    # high T: dichroism destroys the linear component at the clusters...
    assert DoLP_hi.min() < 0.01
    # ...while phi grows regardless -- the unmeasurable regime phi can't see
    assert np.abs(np.real(phi_hi)).max() > 2.0
    assert np.abs(np.real(phi_hi)).max() > np.abs(np.real(phi_lo)).max()


def test_psi_principal_value_range():
    X = np.linspace(-28000, -16000, 1000)
    p = dict(Elem='Rb', Dline='D2', T=125., lcell=1e-3, Bfield=16000.,
             rb85frac=0., **FARADAY)
    [psi] = calculate(X, [1, 0, 0], p, outputs=['psi'])
    assert psi.min() > -np.pi / 2 - 1e-12
    assert psi.max() <= np.pi / 2 + 1e-12
