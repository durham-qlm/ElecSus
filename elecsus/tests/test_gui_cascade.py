"""
Cascade-tab callback tests for the Dash GUI (same pattern as
test_gui_angles.py: the registered callbacks are plain functions, so the
GUI code path is exercised headlessly, minus only the HTTP layer).

The GUI must be a thin wrapper over the core calculate_cascade() - these
tests pin GUI output to the core bit-for-bit in both polariser modes.
"""
import numpy as np

from elecsus import elecsus_dash_gui as gui
from elecsus.elecsus_methods import calculate_cascade

X = np.linspace(-6e3, 6e3, 400)

P1 = dict(Elem='Rb', Dline='D2', T=60., lcell=0.075, Bfield=200., GammaBuf=0.,
          shift=0., Btheta=0., Bphi=0., theta0=0., Constrain=True, rb85frac=72.17)
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
    )


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


def test_gui_qwp_without_polariser_matches_core():
    """QWP alone (no explicit intermediate polariser toggled) also routes to
    the core's reset_polarisation path and matches bit-for-bit."""
    fig, store, err = _run_gui_cascade([], ['S0'], reset_pol='right_cp')
    assert err == ''

    tot_S0 = np.array(store['total']['S0'])
    (_,), (S0_core,) = calculate_cascade(
        X, [1, 0, 0], P1, P2, outputs=['S0'],
        polariser_angle=None, reset_polarisation='right_cp')
    np.testing.assert_allclose(tot_S0, S0_core, atol=1e-12, rtol=0)
