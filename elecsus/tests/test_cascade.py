"""
Pytest validation of calculate_cascade() (converted from
dev-notes/validate_cascade.py, which also generates a figure-based report;
this file is the authoritative CI check - issue #26).

Three properties are pinned:

1. Mathematical identity: two cells of length L/2 with identical parameters
   must reproduce a single cell of length L exactly, in any geometry.

2. Output combination rules:
     intensity outputs (S0, Ix, Iy, I_P45, I_M45):  X_total = S0_cell1 * X_cell2
     fractional Stokes (S1, S2, S3):                X_total = X_cell2

3. Polarisation coupling: the cascade propagates the full E-field between
   cells, so a birefringent cell 1 followed by a dichroic cell 2 must differ
   measurably from the naive product of independent S0 spectra.

NOTE: the core API expects Btheta in RADIANS. The original dev-notes script
used Btheta=90 (radians!) for its "Voigt" case - actually an oblique
geometry of 90 mod 2pi ~ 2.04 rad. Here Voigt is np.pi/2.
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, calculate_cascade

X = np.linspace(-6000, 6000, 600)  # MHz

FARADAY = dict(Btheta=0.0, Bphi=0.0)
VOIGT = dict(Btheta=np.pi / 2, Bphi=0.0)

INTENSITY_OUTPUTS = ['S0', 'Ix', 'Iy', 'I_P45', 'I_M45']
STOKES_OUTPUTS = ['S1', 'S2', 'S3']


@pytest.mark.parametrize('geom, Bfield', [
    pytest.param({}, 0, id='zero-field'),
    pytest.param(FARADAY, 200, id='faraday-200G'),
    pytest.param(VOIGT, 200, id='voigt-200G'),
], )
def test_cascade_identity(geom, Bfield):
    """Two half-length cells must reproduce one full-length cell exactly."""
    L = 75e-3
    p_full = dict(Elem='Rb', Dline='D2', T=20.0, lcell=L, Bfield=Bfield, **geom)
    p_half = dict(Elem='Rb', Dline='D2', T=20.0, lcell=L / 2, Bfield=Bfield, **geom)

    [S0_single] = calculate(X, [1, 0, 0], p_full, outputs=['S0'])
    (S0_1,), (S0_2,) = calculate_cascade(X, [1, 0, 0], p_half, p_half, outputs=['S0'])
    residual = np.max(np.abs((S0_1 * S0_2).real - S0_single.real))
    assert residual < 1e-10


def test_cascade_combination_rules():
    """Cell-2 outputs combine with cell-1 outputs by the documented rules."""
    L = 75e-3
    p_full = dict(Elem='Rb', Dline='D2', T=25.0, lcell=L, Bfield=150, **FARADAY)
    p_half = dict(Elem='Rb', Dline='D2', T=25.0, lcell=L / 2, Bfield=150, **FARADAY)

    outs = INTENSITY_OUTPUTS + STOKES_OUTPUTS
    ref = dict(zip(outs, calculate(X, [1, 0, 0], p_full, outputs=outs)))
    cell1_raw, cell2_raw = calculate_cascade(X, [1, 0, 0], p_half, p_half, outputs=outs)
    c1 = dict(zip(outs, cell1_raw))
    c2 = dict(zip(outs, cell2_raw))

    for key in INTENSITY_OUTPUTS:
        residual = np.max(np.abs((c1['S0'] * c2[key]).real - ref[key].real))
        assert residual < 1e-10, f'{key}: intensity combination rule violated'

    for key in STOKES_OUTPUTS:
        residual = np.max(np.abs(c2[key].real - ref[key].real))
        assert residual < 1e-10, f'{key}: fractional-Stokes combination rule violated'


def test_cascade_polarisation_coupling():
    """A Voigt (birefringent) cell followed by a Faraday (dichroic) cell must
    differ from the naive product of independently-computed S0 spectra.

    The input must NOT be an eigenmode of cell 1: x-polarised light in Voigt
    geometry (B along x) drives only the pi transition and exits with its
    polarisation unchanged, making cascade == naive product identically.
    Light at 45 degrees to B drives both linear eigenmodes, so cell 1's
    birefringence produces elliptical light and the coupling is real."""
    p_voigt = dict(Elem='Rb', Dline='D2', T=60.0, lcell=50e-3, Bfield=300, **VOIGT)
    p_faraday = dict(Elem='Rb', Dline='D2', T=40.0, lcell=75e-3, Bfield=100, **FARADAY)
    E_in = [1 / np.sqrt(2), 1 / np.sqrt(2), 0]

    (S0_1t,), (S0_2t,) = calculate_cascade(X, E_in, p_voigt, p_faraday, outputs=['S0'])
    cascade = (S0_1t * S0_2t).real

    [S0_1i] = calculate(X, E_in, p_voigt, outputs=['S0'])
    [S0_2i] = calculate(X, E_in, p_faraday, outputs=['S0'])
    naive = (S0_1i * S0_2i).real

    assert np.max(np.abs(cascade - naive)) > 1e-3
    # and both stay physical
    for arr in (cascade, naive):
        assert np.all(np.isfinite(arr))
        assert arr.min() > -1e-9 and arr.max() < 1.0 + 1e-9


# ─────────────────────────────────────────────────────────────────────────────
# Intermediate polariser (MOF topology): cell 1 → linear polariser → cell 2
# ─────────────────────────────────────────────────────────────────────────────

def test_polariser_malus_law_off_resonance():
    """Far off resonance the vapour does nothing, so the polariser is the only
    active element and the total transmission must follow Malus's law."""
    X_far = np.linspace(80e3, 90e3, 20)  # MHz, ~80 GHz detuned
    p = dict(Elem='Rb', Dline='D2', T=20.0, lcell=10e-3, Bfield=0)

    for angle in (0.0, np.pi / 6, np.pi / 4, np.pi / 3):
        (_,), (S0_tot,) = calculate_cascade(X_far, [1, 0, 0], p, p,
                                            outputs=['S0'], polariser_angle=angle)
        np.testing.assert_allclose(S0_tot, np.cos(angle) ** 2, atol=1e-3,
                                   err_msg=f'Malus law violated at angle {angle}')


def test_crossed_polariser_no_rotation_blocks_everything():
    """B = 0: no optical rotation in cell 1, so a crossed polariser transmits
    nothing and the instrument total must be ~0 at ALL detunings - finitely,
    with no NaN from the 0/0 wings (the reason totals are normalised to the
    original input)."""
    p = dict(Elem='Rb', Dline='D2', T=20.0, lcell=75e-3, Bfield=0)
    (_,), (S0_tot,) = calculate_cascade(X, [1, 0, 0], p, p,
                                        outputs=['S0'], polariser_angle=np.pi / 2)
    assert np.all(np.isfinite(S0_tot))
    assert np.max(np.abs(S0_tot)) < 1e-12


def test_mof_crossed_polariser_cascade():
    """Same-species crossed-polariser MOF pair (laser input): finite wing
    transmission where cell 1 rotates the polarisation, fully blocked far
    wings, and no NaNs anywhere. Faraday geometry, as in a standard MOF."""
    p1 = dict(Elem='Rb', Dline='D2', T=60.0, lcell=75e-3, Bfield=200, **FARADAY)
    p2 = dict(Elem='Rb', Dline='D2', T=80.0, lcell=75e-3, Bfield=50, **FARADAY)

    outs1, outs2 = calculate_cascade(X, [1, 0, 0], p1, p2,
                                     outputs=['S0', 'Ix', 'Iy'],
                                     polariser_angle=np.pi / 2)
    S0_1 = outs1[0]
    S0_tot, Ix_tot, Iy_tot = outs2

    for arr in (S0_tot, Ix_tot, Iy_tot):
        assert np.all(np.isfinite(arr)), 'NaN/inf in crossed-polariser totals'

    # physical bounds: 0 <= totals <= 1, and components sum to the total
    assert S0_tot.min() > -1e-12 and S0_tot.max() < 1.0 + 1e-9
    np.testing.assert_allclose(Ix_tot + Iy_tot, S0_tot, atol=1e-12)

    # the MOF transmits near the line where cell 1 rotates the polarisation
    assert S0_tot.max() > 0.01, 'no MOF transmission peak found'

    # far from resonance (~80 GHz) there is no rotation left and the crossed
    # polarisers must block everything. (At only +-6 GHz residual Faraday
    # rotation still leaks a few percent - that is real MOF wing behaviour,
    # so the blocking check belongs at genuinely large detuning.)
    X_far = np.linspace(80e3, 85e3, 10)
    (_,), (S0_far,) = calculate_cascade(X_far, [1, 0, 0], p1, p2,
                                        outputs=['S0'], polariser_angle=np.pi / 2)
    assert np.all(np.isfinite(S0_far))
    assert np.max(S0_far) < 1e-4, 'far wings should be blocked by crossed polarisers'

    # totals can never exceed what cell 1 alone transmits
    assert np.all(S0_tot <= S0_1.real + 1e-9)


def test_polariser_at_zero_equals_open_path_for_x_input():
    """A polariser at 0 degrees is transparent to purely-x-polarised light: in
    Faraday-free (B=0) geometry the totals must equal the no-polariser
    cascade combination exactly."""
    p1 = dict(Elem='Rb', Dline='D2', T=30.0, lcell=25e-3, Bfield=0)
    p2 = dict(Elem='Rb', Dline='D2', T=50.0, lcell=25e-3, Bfield=0)

    (S0_1,), (S0_2,) = calculate_cascade(X, [1, 0, 0], p1, p2, outputs=['S0'])
    reference_total = (S0_1 * S0_2).real

    (_,), (S0_tot,) = calculate_cascade(X, [1, 0, 0], p1, p2,
                                        outputs=['S0'], polariser_angle=0.0)
    np.testing.assert_allclose(S0_tot, reference_total, atol=1e-12)


def test_polariser_mode_rejects_medium_outputs():
    """Medium properties are per-cell quantities with no instrument-total
    meaning; requesting them with a polariser must raise, not mislead."""
    p = dict(Elem='Rb', Dline='D2', T=20.0, lcell=10e-3, Bfield=0)
    with pytest.raises(ValueError, match='polariser'):
        calculate_cascade(X, [1, 0, 0], p, p,
                          outputs=['S0', 'alphaPlus'], polariser_angle=np.pi / 2)


# ─────────────────────────────────────────────────────────────────────────────
# reset_polarisation (wing-selector topology): cell 1 -> [polariser] -> fixed
# circular QWP -> cell 2. Models real MOF instruments such as VAMOS
# (Oliviero 1998 PhD thesis, Fig. 3.3) where a quarter-wave plate after the
# crossed-polariser filter always prepares the same circular state for the
# second (wing-selector) cell, regardless of wavelength.
# ─────────────────────────────────────────────────────────────────────────────

def test_reset_polarisation_rejects_bad_value():
    p = dict(Elem='Rb', Dline='D2', T=20.0, lcell=10e-3, Bfield=0)
    with pytest.raises(ValueError, match='reset_polarisation'):
        calculate_cascade(X, [1, 0, 0], p, p, outputs=['S0'],
                          reset_polarisation='bogus')


def test_reset_polarisation_left_right_differ_with_real_signal():
    """With an actual Zeeman-split MOF signal, left_cp and right_cp couple to
    opposite wings and must give visibly different totals (this is the whole
    point of a wing selector). Without a B-field there is no asymmetry to
    select between, so this needs a real field."""
    p1 = dict(Elem='Rb', Dline='D2', T=90.0, lcell=75e-3, Bfield=500, **FARADAY)
    p2 = dict(Elem='Rb', Dline='D2', T=90.0, lcell=75e-3, Bfield=2000, **FARADAY)

    (_,), (S0_left,) = calculate_cascade(X, [1, 0, 0], p1, p2, outputs=['S0'],
                                         polariser_angle=np.pi / 2,
                                         reset_polarisation='left_cp')
    (_,), (S0_right,) = calculate_cascade(X, [1, 0, 0], p1, p2, outputs=['S0'],
                                          polariser_angle=np.pi / 2,
                                          reset_polarisation='right_cp')
    assert np.max(np.abs(S0_left - S0_right)) > 1e-3


def test_reset_polarisation_preserves_circularity_through_faraday_cell():
    """A pure circular state is an eigenmode of a Faraday-geometry cell, so
    if reset_polarisation genuinely prepares fixed circular light, the
    instrument-total S3 (circular Stokes parameter) coming out of a Faraday
    cell 2 must be exactly +-1 wherever transmission is non-negligible,
    regardless of whatever polarisation state cell 1 produced."""
    p1 = dict(Elem='Rb', Dline='D2', T=60.0, lcell=75e-3, Bfield=300, **VOIGT)
    p2 = dict(Elem='Rb', Dline='D2', T=90.0, lcell=75e-3, Bfield=1000, **FARADAY)
    E_in = [1 / np.sqrt(2), 1 / np.sqrt(2), 0]  # not an eigenmode of cell 1

    (_, _), (S0_tot, S3_tot) = calculate_cascade(
        X, E_in, p1, p2, outputs=['S0', 'S3'],
        polariser_angle=np.pi / 3, reset_polarisation='right_cp')

    mask = S0_tot > 1e-6 * S0_tot.max()
    assert mask.sum() > 10, 'not enough non-negligible transmission to test'
    np.testing.assert_allclose(S3_tot[mask], -1.0, atol=1e-6)


def test_reset_polarisation_matches_independent_single_cell_product():
    """The whole point of reset_polarisation: a linear polariser followed by
    a QWP prepares a FIXED circular state for cell 2, independent of
    wavelength. This means the instrument total must equal the simple
    product of two INDEPENDENT single-cell calculations -- cell 1 read
    through a crossed linear analyser (Iy, x-polarised input), and cell 2
    driven directly with the corresponding fixed circular input -- exactly
    the manual two-run-and-multiply approach this kwarg replaces with a
    single call. This is the crossed-polariser-Filtro + circular-input-WS
    topology of a real magneto-optical filter instrument, e.g. VAMOS
    (Oliviero 1998 PhD thesis, Fig. 3.3: Trasmissione FILTRO / WS / VAMOS)."""
    p1 = dict(Elem='Na', Dline='D2', T=250.0, lcell=100e-3, Bfield=1500, **FARADAY)
    p2 = dict(Elem='Na', Dline='D2', T=250.0, lcell=100e-3, Bfield=4000, **FARADAY)

    (_,), (S0_tot,) = calculate_cascade(X, [1, 0, 0], p1, p2, outputs=['S0'],
                                        polariser_angle=np.pi / 2,
                                        reset_polarisation='left_cp')

    [T_filtro] = calculate(X, [1, 0, 0], p1, outputs=['Iy'])
    E_in_left_cp = [1 / np.sqrt(2), 1j / np.sqrt(2), 0]
    [T_ws] = calculate(X, E_in_left_cp, p2, outputs=['S0'])

    manual_total = T_filtro.real * T_ws.real
    np.testing.assert_allclose(S0_tot, manual_total, atol=1e-10)
