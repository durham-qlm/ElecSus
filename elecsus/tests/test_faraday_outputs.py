"""
nPlus, nMinus, phi output keys (requested by users for the GUI: refractive
indices and Faraday rotation angle).

Background: these circular-basis quantities were already computed inside
get_spectra() (used internally to derive TransLeft/TransRight and
alphaPlus/alphaMinus) but were not in the returned output dict, so
outputs=['nPlus'] raised KeyError. They are exposed here with the same
caveat as the pre-existing TransLeft/alphaPlus outputs: physically
meaningful only in Faraday geometry (Btheta=0); the code computes them
unconditionally for any Btheta, an ES3-inherited behaviour kept for
backwards compatibility (Sec. 4.3 of the manual).
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate
from elecsus.libs.spectra import calc_chi
import elecsus.libs.AtomConstants as AC

X = np.linspace(-4e3, 4e3, 200)

P_FARADAY = dict(Elem='Rb', Dline='D2', T=60., lcell=75e-3, Bfield=300.,
                 Btheta=0., Bphi=0.)


def test_nplus_nminus_match_susceptibility_definition():
    """nPlus/nMinus must equal sqrt(1 + chiPlus/Minus) at every detuning,
    for any geometry (this is a definitional identity, not Faraday-specific)."""
    nPlus, nMinus, ChiPlus, ChiMinus = calculate(
        X, [1, 0, 0], P_FARADAY, outputs=['nPlus', 'nMinus', 'ChiPlus', 'ChiMinus'])
    np.testing.assert_allclose(nPlus, np.sqrt(1.0 + ChiPlus))
    np.testing.assert_allclose(nMinus, np.sqrt(1.0 + ChiMinus))


def test_phi_matches_faraday_rotation_formula():
    """phi is the standard Faraday-rotation half-difference of the real
    refractive indices, scaled by wavenumber*lcell."""
    nPlus, nMinus, phi = calculate(
        X, [1, 0, 0], P_FARADAY, outputs=['nPlus', 'nMinus', 'phi'])
    wavenumber = AC.transitions['RbD2'].wavevectorMagnitude
    lcell = P_FARADAY['lcell']
    expected = 0.5 * wavenumber * lcell * (nMinus.real - nPlus.real)
    np.testing.assert_allclose(np.asarray(phi).real, expected)


def test_outputs_available_at_oblique_geometry_too():
    """The keys must not raise even off-Faraday (they're computed
    unconditionally, per the documented caveat) -- this pins that the
    lookup itself is geometry-independent, not that the physics is valid."""
    p = dict(P_FARADAY, Btheta=np.radians(45.))
    nPlus, nMinus, phi, phiP, phiM = calculate(
        X, [1, 0, 0], p, outputs=['nPlus', 'nMinus', 'phi', 'phiPlus', 'phiMinus'])
    assert np.all(np.isfinite(nPlus)) and np.all(np.isfinite(nMinus))
    for a in (phi, phiP, phiM):
        assert np.all(np.isfinite(np.asarray(a).real))


# ── Convention pins: the Plus/Minus suffix on Chi/n is the sigma label of the
#    transitions, and with B along +z the [1,+i]/sqrt2 field drives sigma+
#    (in ES3 and the 4.0.0 alphas ChiPlus was the sigma- manifold and B was
#    modelled reversed). Lock it so a future refactor can't silently flip it. ──────────

def test_chiplus_is_sigma_plus_manifold():
    """High-field Faraday: the sigma+ (dm=+1) manifold sits at higher frequency
    (more positive detuning) than sigma-. ChiPlus must be the one that shifts
    positive."""
    Xw = np.linspace(-15e3, 15e3, 4000)
    p = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=6000.,
             Btheta=0., Bphi=0., rb85frac=0.)  # pure Rb87
    ChiPlus, ChiMinus = calculate(Xw, [1, 0, 0], p, outputs=['ChiPlus', 'ChiMinus'])
    peak_plus = Xw[np.argmax(ChiPlus.imag)]
    peak_minus = Xw[np.argmax(ChiMinus.imag)]
    assert peak_plus > peak_minus, 'ChiPlus (sigma+) must sit at higher detuning'


def test_nplus_governs_1plusi_circular_light():
    """In Faraday geometry the circular states are the eigenmodes. nPlus (via
    TransLeft = exp(-2 Im(nPlus) k L)) must govern the [1,+i]/sqrt2 mode."""
    import elecsus.libs.spectra as sp
    p = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=6000.,
             Btheta=0., Bphi=0., rb85frac=0.)
    nPlus, nMinus = calculate(X, [1, 0, 0], p, outputs=['nPlus', 'nMinus'])
    kL = sp.AC.transitions['RbD2'].wavevectorMagnitude * p['lcell']
    trans_nplus = np.exp(-2 * nPlus.imag * kL)
    [S0_Rmode] = calculate(X, np.array([1, 1j, 0]) / np.sqrt(2), p, outputs=['S0'])
    np.testing.assert_allclose(np.asarray(S0_Rmode).real, trans_nplus, atol=1e-9)


@pytest.mark.parametrize('name', ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy',
                                  'I_P45', 'I_M45', 'Il', 'Ir', 'psi', 'DoLP'])
def test_intensity_outputs_are_real(name):
    """[ES4 MOD] Intensities and Stokes parameters are real by construction
    (E . E*), but were returned as complex128 with a ~1e-17 imaginary residue.
    A complex dtype pushes a ComplexWarning into user plotting/fitting code and
    makes np.max compare lexicographically rather than by magnitude, so the
    real part is now taken where each observable is formed."""
    p = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=100.,
             Btheta=np.pi / 4, Bphi=0.)
    [out] = calculate(X, [1, 0, 0], p, outputs=[name])
    assert np.asarray(out).dtype == np.float64, f'{name} is {np.asarray(out).dtype}'


@pytest.mark.parametrize('wavevector_dirn', [1, -1])
@pytest.mark.parametrize('Btheta', [0., np.pi])
def test_phi_is_the_measured_rotation_for_b_along_plus_and_minus_z(Btheta, wavevector_dirn):
    """[ES4 FIX] In the line wings (weak dichroism) phi, the
    rotation from x towards y, must equal the ellipse azimuth psi, which the
    field propagation gives. With B along -z, [1,+i] light is sigma-, so the
    rotation is -(phiMinus - phiPlus)/2; phi kept the +z sign and came out
    opposite to psi. A reversed pass returns the forward pass's lab-frame
    field, so its phi is the same."""
    Xw = np.array([-8000., -6000., 6000., 8000.])
    p = dict(Elem='Rb', Dline='D2', T=60., lcell=75e-3, Bfield=200.,
             Btheta=Btheta, wavevector_dirn=wavevector_dirn)
    phi, psi, phiP, phiM = calculate(Xw, [1, 0, 0], p,
                                     outputs=['phi', 'psi', 'phiPlus', 'phiMinus'])
    assert np.min(np.abs(psi)) > 0.04                     # a rotation to get right
    np.testing.assert_allclose(phi, psi, rtol=0, atol=1e-9)
    np.testing.assert_allclose(phi, np.sign(np.cos(Btheta)) * (phiM - phiP) / 2,
                               rtol=0, atol=1e-12)
