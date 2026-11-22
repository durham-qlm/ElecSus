"""
Cascade-tab callback tests for the Dash GUI (same pattern as
test_gui_angles.py: the registered callbacks are plain functions, so the
GUI code path is exercised headlessly, minus only the HTTP layer).

The GUI must be a thin wrapper over the core calculate_cascade() - these
tests pin GUI output to the core bit-for-bit in both polariser modes.
"""
import numpy as np
import pytest

from elecsus import elecsus_dash_gui as gui
from elecsus.elecsus_methods import calculate, calculate_cascade

X = np.linspace(-6e3, 6e3, 400)

P1 = dict(Elem='Rb', Dline='D2', T=60., lcell=0.075, Bfield=200., GammaBuf=0.,
          shift=0., Btheta=0., Bphi=0., Constrain=True, rb85frac=72.17)
P2 = dict(P1, T=80., Bfield=50.)


def _run_gui_cascade(use_pol, outputs, reset_pol='none'):
    """Invoke the registered cascade callback with GUI-equivalent arguments."""
    return gui.compute_cascade(
        1,
        'Rb', 'D2', 60., 0.075, 200., 0., 0., 0., 0., 0., 'Linear (x)',
        'None', 0., 293.15,
        72.17, 0.01, 6.73,
        80., 0.075, 50., 0., 0., 0., 0., 'None', 0., 293.15,
        use_pol, 90., reset_pol,
        -6., 6., 400,
        outputs,
    )[:3]


def test_cell2_species_synced_to_cell1():
    assert gui.sync_cell2_species('K', 'D1') == ('K', 'D1')
    assert gui.sync_cell2_species('Cs', 'D2') == ('Cs', 'D2')


def test_gui_mof_cascade_matches_core():
    """Crossed-polariser (MOF) mode: GUI totals == core instrument totals."""
    fig, store, err = _run_gui_cascade(['on'], ['S0', 'Iy'])
    assert err == ''

    tot_S0 = np.array(store['total']['S0'])
    tot_Iy = np.array(store['total']['Iy'])
    assert np.all(np.isfinite(tot_S0)) and np.all(np.isfinite(tot_Iy))

    (_, _), (S0_core, Iy_core) = calculate_cascade(
        X, [1, 0, 0], P1, P2, outputs=['S0', 'Iy'], polariser_angle=np.pi / 2)
    np.testing.assert_allclose(tot_S0, S0_core, atol=1e-12, rtol=0)
    np.testing.assert_allclose(tot_Iy, Iy_core, atol=1e-12, rtol=0)

    # MOF physics: finite transmission near the line for these parameters
    assert tot_S0.max() > 0.5


def test_gui_open_cascade_matches_combination_rule():
    """No polariser: GUI totals == S0_cell1 * S0_cell2 (documented rule)."""
    fig, store, err = _run_gui_cascade([], ['S0'])
    assert err == ''

    tot = np.array(store['total']['S0'])
    (S0_1,), (S0_2,) = calculate_cascade(X, [1, 0, 0], P1, P2, outputs=['S0'])
    np.testing.assert_allclose(tot, (S0_1 * S0_2).real, atol=1e-12, rtol=0)


def test_gui_cascade_fractional_stokes_pass_through():
    """No polariser: fractional Stokes totals are the cell-2 values unchanged."""
    fig, store, err = _run_gui_cascade([], ['S0', 'S1'])
    assert err == ''

    (_, _), (_, S1_2) = calculate_cascade(X, [1, 0, 0], P1, P2, outputs=['S0', 'S1'])
    np.testing.assert_allclose(np.array(store['total']['S1']), S1_2.real,
                               atol=1e-12, rtol=0)


def test_gui_qwp_cascade_matches_core():
    """Wing-selector (QWP) mode: GUI totals == core instrument totals with
    reset_polarisation set (the VAMOS-style topology, Oliviero 1998)."""
    fig, store, err = _run_gui_cascade(['on'], ['S0'], reset_pol='left_cp')
    assert err == ''

    tot_S0 = np.array(store['total']['S0'])
    assert np.all(np.isfinite(tot_S0))

    (_,), (S0_core,) = calculate_cascade(
        X, [1, 0, 0], P1, P2, outputs=['S0'],
        polariser_angle=np.pi / 2, reset_polarisation='left_cp')
    np.testing.assert_allclose(tot_S0, S0_core, atol=1e-12, rtol=0)


def test_gui_qwp_without_polariser_is_refused():
    """[ES4] The core refuses reset_polarisation without a
    polariser; the GUI says so in its own terms instead of computing."""
    fig, store, err = _run_gui_cascade([], ['S0'], reset_pol='right_cp')
    assert 'needs the polariser' in err and 'Polariser between cells' in err


def test_gui_open_cascade_total_is_finite_behind_an_opaque_cell():
    """No polariser, two hot cells (120 C, 75 mm, 100 G): cell 1 passes no
    light at line centre. The total there is 0, not NaN, and the callback
    reports no warning."""
    args = (1,
            'Rb', 'D2', 120., 0.075, 100., 0., 0., 0., 0., 0., 'Linear (x)',
            'None', 0., 293.15,
            72.17, 0.01, 6.73,
            120., 0.075, 100., 0., 0., 0., 0., 'None', 0., 293.15,
            [], 90., 'none',
            -15., 15., 800,
            ['S0', 'S1', 'Ix'])
    fig, store, err, messages = gui.compute_cascade(*args)[:4]
    assert err == ''
    assert not messages
    cell_1 = np.array(store['cell1']['S0'])
    for name in ('S0', 'S1', 'Ix'):
        total = np.array(store['total'][name], dtype=float)
        assert np.all(np.isfinite(total)), name
    assert np.any(cell_1 < 1e-300)                        # the case is there
    total_S0 = np.array(store['total']['S0'], dtype=float)
    assert np.all(total_S0[cell_1 == 0] == 0)


# --- The lower panel: cell 2 on its own ------------------------------------------
# It shows the light reaching cell 2 and cell 2's transmission of it; their
# product is the total S0, in every optic arrangement.

def _alone(store):
    into = np.array(store['cell2']['into'], dtype=float)
    T_2 = np.array([np.nan if v is None else v for v in store['cell2']['transmission']])
    return into, T_2


@pytest.mark.parametrize('use_pol, reset_pol, state', [
    (['on'], 'left_cp', [1 / np.sqrt(2), 1j / np.sqrt(2), 0.]),
    (['on'], 'right_cp', [1 / np.sqrt(2), -1j / np.sqrt(2), 0.]),
    (['on'], 'none', [0., 1., 0.]),          # linear at the polariser's 90 deg
])
def test_with_the_polariser_cell_2_alone_is_its_spectrum_for_the_fixed_state(
        use_pol, reset_pol, state):
    """Light reaches cell 2 in one state at every detuning, so the lower panel
    is cell 2's own spectrum for that input -- the core's, called by name --
    and the light reaching it is cell 1's I_y (crossed polariser)."""
    fig, store, err = _run_gui_cascade(use_pol, ['S0'], reset_pol)
    assert err == ''
    into, T_2 = _alone(store)
    [want] = calculate(X, state, P2, outputs=['S0'])
    np.testing.assert_allclose(T_2, np.real(want), rtol=1e-12, atol=1e-15)
    (Iy_1,), _ = calculate_cascade(X, [1, 0, 0], P1, P2, outputs=['Iy'])
    np.testing.assert_allclose(into, np.real(Iy_1), rtol=1e-9, atol=1e-15)
    np.testing.assert_allclose(into * T_2, store['total']['S0'], rtol=1e-9, atol=1e-15)
    # its own panel, below the first
    names = [t.name for t in fig.data]
    assert names[-2].startswith('Light reaching cell 2') and names[-1].startswith('Cell 2 alone')
    assert fig.data[-1].yaxis != fig.data[0].yaxis


def test_without_the_polariser_cell_2_alone_is_its_per_stage_transmission():
    fig, store, err = _run_gui_cascade([], ['S0'])
    assert err == ''
    into, T_2 = _alone(store)
    (S0_1,), (S0_2,) = calculate_cascade(X, [1, 0, 0], P1, P2, outputs=['S0'])
    np.testing.assert_allclose(into, np.real(S0_1), rtol=1e-12)
    np.testing.assert_allclose(T_2, np.real(S0_2), rtol=1e-12)
    np.testing.assert_allclose(into * T_2, store['total']['S0'], rtol=1e-12)


def test_cell_2_alone_has_a_gap_where_no_light_reaches_it():
    """No polariser, cell 1 opaque at line centre: cell 2's transmission of
    no light is undefined, a gap (None in the store, NaN in the plot and the
    export), and finite everywhere else."""
    fig, store, err = gui.compute_cascade(
        1, 'Rb', 'D2', 120., 0.075, 100., 0., 0., 0., 0., 0., 'Linear (x)',
        'None', 0., 293.15, 72.17, 0.01, 6.73,
        120., 0.075, 100., 0., 0., 0., 0., 'None', 0., 293.15,
        [], 90., 'none', -15., 15., 800, ['S0'])[:3]
    assert err == ''
    into, T_2 = _alone(store)
    dark = into < np.finfo(float).tiny
    assert dark.any() and not dark.all()
    assert np.all(np.isnan(T_2[dark])) and np.all(np.isfinite(T_2[~dark]))
    assert np.all(np.isnan(np.asarray(fig.data[-1].y, dtype=float)[dark]))
    rows = gui.export_cascade(1, store)['content'].splitlines()[1:]
    assert sum(r.endswith(',nan') for r in rows) == dark.sum()


def test_the_wing_selector_absorbs_the_sigma_plus_wing_for_left_circular_light():
    """VAMOS-like K D1 wing selector at 2700 G: with the QWP set to left
    circular, cell 2 alone absorbs at +5.0 GHz (sigma+, B along +z) and passes
    -5.0 GHz; right circular the other way round."""
    def ws(reset):
        _, store, err = gui.compute_cascade(
            1, 'K', 'D1', 120., 0.025, 1400., 0., 0., 0., 0., 0., 'Linear (x)',
            'None', 0., 293.15, 72.17, 0.01, 6.73,
            120., 0.025, 2700., 0., 0., 0., 0., 'None', 0., 293.15,
            ['on'], 90., reset, -5.04, 5.04, 11, ['S0'])[:3]
        assert err == ''
        return _alone(store)[1]
    red, blue = ws('left_cp')[[0, -1]]          # -5.04 and +5.04 GHz
    assert blue < 1e-3 and red > 0.99
    red, blue = ws('right_cp')[[0, -1]]
    assert red < 1e-3 and blue > 0.99
