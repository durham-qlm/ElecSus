"""
Pytest validation of calculate_cascade() (converted from a script that also
generated a figure-based report; this file is the authoritative CI check).

Three properties are pinned:

1. Mathematical identity: two cells of length L/2 with identical parameters
   must reproduce a single cell of length L exactly, in any geometry.

2. Output combination rules:
     intensity outputs (S0, Ix, Iy, I_P45, I_M45):  X_total = S0_cell1 * X_cell2
     fractional Stokes (S1, S2, S3):                X_total = X_cell2

3. Polarisation coupling: the cascade propagates the full E-field between
   cells, so a birefringent cell 1 followed by a dichroic cell 2 must differ
   measurably from the naive product of independent S0 spectra.

NOTE: the core API expects Btheta in RADIANS. The original script
used Btheta=90 (radians!) for its "Voigt" case - actually an oblique
geometry of 90 mod 2pi ~ 2.04 rad. Here Voigt is np.pi/2.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, calculate_cascade
from elecsus.libs import spectra

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
    # right_cp is [1,-i]/sqrt2: S3 = +1 in Born & Wolf's sign. This asserted -1,
    # the sign the instrument totals kept from the 4.0.0 alphas.
    np.testing.assert_allclose(S3_tot[mask], 1.0, atol=1e-6)


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


@pytest.mark.parametrize('reset, S3', [('left_cp', -1.0), ('right_cp', 1.0)])
def test_instrument_total_S3_follows_born_and_wolf(reset, S3):
    """[ES4] Far from resonance both cells are transparent, so the
    instrument-total S3 is that of the prepared state: left_cp is [1, +i]/sqrt2,
    left-circular in Born & Wolf (Sec. 1.4.2) for fields varying as
    exp(i(kz - wt)), with S3 = (Ir - Il)/S0 = -1; right_cp gives +1. The manual
    (ch. 5) and output_list state this convention; calculate() follows it.
    _totals_from_field kept the alphas' sign, corrected elsewhere. (The
    polariser at 0 passes the x input whole; reset_polarisation needs one.)"""
    X_far = np.array([-30e3, 30e3])
    p = dict(Elem='Rb', Dline='D2', T=20., lcell=1e-3)
    (_,), (S3_tot,) = calculate_cascade(X_far, [1, 0, 0], p, p, outputs=['S3'],
                                        polariser_angle=0., reset_polarisation=reset)
    np.testing.assert_allclose(S3_tot, S3, rtol=0, atol=1e-6)


def test_instrument_totals_are_calculate_s_outputs_of_the_same_field():
    """[ES4] With the polariser at 0 and an x input far from
    resonance, the light reaching cell 2 is cell 1's output, so the instrument
    totals must be calculate()'s outputs for cell 2 (normalised to the same
    input), every one of them, S3 included, in an oblique field that makes all
    of them nonzero."""
    X = np.linspace(-5e3, 5e3, 61)
    p1 = dict(Elem='Rb', Dline='D2', T=20., lcell=1e-3)
    p2 = dict(Elem='Rb', Dline='D2', T=80., lcell=5e-3, Bfield=300.,
              Btheta=np.radians(50.), Bphi=0.3)
    names = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45']
    E_in = [1., 0., 0.]
    (E1,), _ = calculate_cascade(X, E_in, p1, p2, outputs=['E_out'],
                                  polariser_angle=0.)
    _, totals = calculate_cascade(X, E_in, p1, p2, outputs=names, polariser_angle=0.)
    E_mid = np.array(E1)
    E_mid[1] = 0.        # what the polariser at 0 passes (E_y is 0 here anyway)
    direct = calculate(X, E_mid, p2, outputs=names)
    S0_1 = np.real((E_mid * E_mid.conj()).sum(axis=0))
    for name, t, d in zip(names, totals, direct):
        expected = d if name in ('S1', 'S2', 'S3') else S0_1 * d
        np.testing.assert_allclose(t, expected, rtol=0, atol=1e-12, err_msg=name)
    assert np.max(np.abs(totals[3])) > 0.05          # a circular part to get right


def test_reset_polarisation_without_a_polariser_is_refused():
    """[ES4] Only a polariser and
    a quarter-wave plate turn every state into one circular state; without the
    polariser the reset kept the full intensity, which no passive optic does."""
    p = dict(Elem='Rb', Dline='D2', T=20., lcell=1e-3)
    for reset in ('left_cp', 'right_cp'):
        with pytest.raises(ValueError, match='needs polariser_angle'):
            calculate_cascade(X, [1, 0, 0], p, p, outputs=['S0'],
                              reset_polarisation=reset)


@pytest.mark.parametrize('alpha', [0.0, np.pi / 8, np.pi / 4, 2 * np.pi / 3])
def test_instrument_totals_of_a_linear_polariser(alpha):
    """[ES4] (S1, S2, I_P45 and I_M45 of the instrument-total path were
    untested.) Far from resonance the cells are transparent, so after x
    light meets a polariser at alpha: Malus, S0 = cos^2 alpha; the state is
    linear at alpha, S1 = cos 2 alpha, S2 = sin 2 alpha, S3 = 0; and a +-45 deg
    analyser passes cos^2 alpha cos^2(alpha -+ pi/4)."""
    X_far = np.array([-30e3, 30e3])
    p = dict(Elem='Rb', Dline='D2', T=20., lcell=1e-3)
    names = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45']
    _, totals = calculate_cascade(X_far, [1, 0, 0], p, p, outputs=names,
                                  polariser_angle=alpha)
    c2 = np.cos(alpha) ** 2
    expected = [c2, np.cos(2 * alpha), np.sin(2 * alpha), 0., c2 * np.cos(alpha) ** 2,
                c2 * np.sin(alpha) ** 2, c2 * np.cos(alpha - np.pi / 4) ** 2,
                c2 * np.cos(alpha + np.pi / 4) ** 2]
    for name, value, want in zip(names, totals, expected):
        np.testing.assert_allclose(value, want, rtol=0, atol=1e-6, err_msg=name)


# --- No light into cell 2 ----------------------------------------------------
# Behind a crossed polariser, and behind an opaque cell 1, cell 2 gets no
# light at line centre (|E|^2 = 0, or a denormal). Its per-stage outputs,
# normalised by that intensity, were 0/0, and numpy warned "overflow / invalid
# value encountered in divide" on every such run, which the GUI showed on its
# cascade tab. There is no warning now. Without a polariser, those per-stage
# outputs were NaN, and so was the documented total S0_cell1 * X_cell2; they
# are 0 there now, as no light leaves an opaque filter. No other number
# changes.

MOF = dict(Elem='Rb', Dline='D2', T=120., lcell=75e-3, Bfield=100.)
X_MOF = np.linspace(-12e3, 12e3, 800)


@pytest.mark.parametrize('kw', [
    pytest.param(dict(polariser_angle=np.pi / 2), id='crossed-polariser'),
    pytest.param(dict(polariser_angle=np.pi / 2, reset_polarisation='left_cp'),
                 id='polariser-and-qwp'),
    pytest.param({}, id='no-polariser')])
def test_no_light_into_cell_2_raises_no_runtime_warning(kw):
    [E1] = calculate(X_MOF, [1, 0, 0], MOF, outputs=['E_out'])
    E1 = np.asarray(E1)
    into_cell_2 = abs(E1[1]) ** 2 if kw else (abs(E1) ** 2).sum(axis=0)
    assert np.any(into_cell_2 < np.finfo(float).tiny)     # the case is there
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        calculate_cascade(X_MOF, [1, 0, 0], MOF, MOF,
                          outputs=['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy'], **kw)


def test_crossed_polariser_totals_are_the_plain_ratios():
    """The instrument totals are the field's intensities over the input's
    (here 1), and S1-S3 those over S0 wherever S0 is above calculate_cascade's
    floor of 1e-15; below it they are not divided, so stay below the floor.
    Both kinds of point occur."""
    names = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45', 'E_out']
    _, out = calculate_cascade(X_MOF, [1, 0, 0], MOF, MOF, outputs=names,
                               polariser_angle=np.pi / 2)
    got = dict(zip(names, out))
    Ex, Ey = np.asarray(got['E_out'])[:2]
    Ix, Iy = abs(Ex) ** 2, abs(Ey) ** 2
    I_P45, I_M45 = abs(Ex + Ey) ** 2 / 2, abs(Ex - Ey) ** 2 / 2
    Il, Ir = abs(Ex - 1j * Ey) ** 2 / 2, abs(Ex + 1j * Ey) ** 2 / 2
    S0 = Ix + Iy
    for name, want in (('S0', S0), ('Ix', Ix), ('Iy', Iy), ('I_P45', I_P45),
                       ('I_M45', I_M45)):
        # atol: a denormal (~1e-314) carries few significant bits
        np.testing.assert_allclose(got[name], want, rtol=1e-12, atol=1e-300,
                                   err_msg=name)
    lit = S0 > 1e-15
    assert lit.any() and (~lit).any()
    for name, diff in (('S1', Ix - Iy), ('S2', I_P45 - I_M45), ('S3', Ir - Il)):
        np.testing.assert_allclose(got[name][lit], diff[lit] / S0[lit], rtol=0,
                                   atol=1e-12, err_msg=name)
        assert np.all(abs(got[name][~lit]) <= 1e-15), name


def test_division_by_the_input_intensity_is_unchanged_where_light_enters():
    """spectra._over_input(I, I_in) is I / I_in to the last bit wherever I_in
    is a normal float, and NaN where it is 0 or denormal, where the plain
    division gave 0/0 or, inside numpy's complex division, 0 * inf."""
    rng = np.random.default_rng(3)
    n = 600
    I = (rng.normal(size=n) + 1j * rng.normal(size=n)) * 10.0 ** rng.uniform(-300, 3, n)
    I_in = 10.0 ** rng.uniform(-300, 3, n)
    I[:40], I_in[:40] = 0, 0                 # no light in, and none out
    I[40:60], I_in[40:60] = 0, 5e-320        # a denormal input intensity
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        got = spectra._over_input(I, I_in)
        got_scalar = spectra._over_input(I[60:], 2.5)
    with np.errstate(all='ignore'):
        plain = I / I_in
    normal = I_in >= np.finfo(float).tiny
    assert got[normal].tobytes() == plain[normal].tobytes()
    assert np.all(np.isnan(got[~normal])) and not np.any(np.isfinite(plain[~normal]))
    assert got_scalar.tobytes() == (I[60:] / 2.5).tobytes()


_OPEN_OUTPUTS = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45', 'Il', 'Ir',
                 'psi', 'DoLP', 'S0_P', 'alphaPlus', 'ChiPlus', 'E_out']
_INTENSITIES = ('S0', 'Ix', 'Iy', 'I_P45', 'I_M45', 'Il', 'Ir', 'S0_P')
_FRACTIONAL = ('S1', 'S2', 'S3', 'psi', 'DoLP')


def test_open_cascade_total_is_zero_where_cell_1_passes_no_light():
    """No polariser: where cell 1 is opaque (its output intensity 0 or a
    denormal), the documented totals S0_cell1 * X_cell2 and X_cell2 are 0,
    not NaN, and no RuntimeWarning is raised. Everywhere else every output is
    the plain per-stage value, to the last bit, and the medium properties are
    untouched."""
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        out_1, out_2 = calculate_cascade(X_MOF, [1, 0, 0], MOF, MOF,
                                         outputs=_OPEN_OUTPUTS)
    cell_1, cell_2 = dict(zip(_OPEN_OUTPUTS, out_1)), dict(zip(_OPEN_OUTPUTS, out_2))
    E1 = np.asarray(cell_1['E_out'])
    dark = (abs(E1) ** 2).sum(axis=0) < np.finfo(float).tiny
    assert 0 < dark.sum() < len(X_MOF)                    # the case is there

    for name in _INTENSITIES:
        total = cell_1['S0'] * cell_2[name]
        assert np.all(np.isfinite(total)), name
        assert np.all(total[dark] == 0), name
    for name in _FRACTIONAL:
        assert np.all(np.isfinite(cell_2[name])), name
        assert np.all(cell_2[name][dark] == 0), name

    # elsewhere: cell 2's outputs for the same input field, from get_spectra
    with np.errstate(all='ignore'):
        plain = dict(zip(_OPEN_OUTPUTS, spectra.get_spectra(
            X_MOF, E1, MOF, outputs=_OPEN_OUTPUTS)))
    for name in _INTENSITIES + _FRACTIONAL:
        got, want = np.asarray(cell_2[name]), np.asarray(plain[name])
        assert got[~dark].tobytes() == want[~dark].tobytes(), name
    for name in ('alphaPlus', 'ChiPlus', 'E_out'):
        assert np.asarray(cell_2[name]).tobytes() == np.asarray(plain[name]).tobytes(), name


def test_poynting_transmission_of_no_light_is_nan_without_a_warning():
    """S0_P for an input field that is 0 at some detunings: NaN there, as the
    plain division gave, but without its RuntimeWarning; unchanged elsewhere."""
    E_in = np.zeros((3, len(X_MOF)), dtype=complex)
    E_in[0] = 1.
    E_in[0, 100:110] = 0.
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        [S0_P] = calculate(X_MOF, E_in, MOF, outputs=['S0_P'])
    [S0_P_x] = calculate(X_MOF, [1, 0, 0], MOF, outputs=['S0_P'])
    assert np.all(np.isnan(S0_P[100:110]))
    lit = np.ones(len(X_MOF), bool)
    lit[100:110] = False
    assert S0_P[lit].tobytes() == np.asarray(S0_P_x)[lit].tobytes()
