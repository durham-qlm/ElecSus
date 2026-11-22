"""Tests for the GUI's v4 model controls and fit-panel additions.

Propagation (non-orthogonal solver, cell windows, S0_P), the independent
Doppler temperature, the bundled sample spectrum, and the extra fit
parameters (Doppler T, linear polarisation angle). As in
test_gui_callbacks.py the callbacks are called directly: the point is that
the GUI passes exactly what the core API would be given, and reports rather
than raises when the core refuses.
"""
import base64
import warnings

import numpy as np
import pytest

from elecsus import elecsus_dash_gui as gui
from elecsus.elecsus_methods import calculate, calculate_cascade

ON = ['on']

SINGLE = dict(
    elem='Rb', dline='D2', T=120., lcell=2e-3, B=239., gbuf=7.9, shift=0.,
    btheta=80., bphi=0., theta0=0., pol='Linear (y)',
    bufgas='None', bufP=0., bufTf=293.15,
    rb85frac=72.17, K40frac=0.01, K41frac=6.73,
    det_start=-6., det_stop=8., det_N=300,
)


def _compute_single(outputs=('S0',), nonortho=gui.NONORTHO_DEFAULT, windows=None, n_window=None,
                    doppT=None, **overrides):
    a = dict(SINGLE, **overrides)
    return gui.compute_single(
        1, a['elem'], a['dline'], a['T'], a['lcell'], a['B'], a['gbuf'],
        a['shift'], a['btheta'], a['bphi'], a['theta0'], a['pol'],
        a['bufgas'], a['bufP'], a['bufTf'],
        a['rb85frac'], a['K40frac'], a['K41frac'],
        a['det_start'], a['det_stop'], a['det_N'], list(outputs),
        nonortho, windows, n_window, doppT)[:3]


def _core(X_ghz, E_in, outputs=('S0',), **extra):
    p = dict(Elem='Rb', Dline='D2', T=120., lcell=2e-3, Bfield=239.,
             GammaBuf=7.9, shift=0., Btheta=np.radians(80.), Bphi=0.,
             rb85frac=72.17, **extra)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return calculate(np.asarray(X_ghz) * 1e3, E_in, p, outputs=list(outputs))


# --- the p_dict the controls produce --------------------------------------

def _p(**kw):
    return gui._p('Rb', 'D2', 60, 0.075, 0, 0, 0, 0, 0, 0, 'None', 0, 293.15, **kw)


def test_defaults_are_the_non_orthogonal_solver_without_windows():
    """As in the core API since 4.0.0b1. The choice is always sent, so the GUI
    does not depend on the core's default."""
    p = _p()
    assert p['use_non_orthogonal'] is True and 'n_window' not in p
    assert p['Constrain'] is True and 'DoppTemp' not in p
    assert _p(nonortho=[])['use_non_orthogonal'] is False


def test_the_layout_starts_with_the_solver_on_and_the_windows_free():
    def find(node, id_):
        if getattr(node, 'id', None) == id_:
            return node
        children = getattr(node, 'children', None)
        for child in (children if isinstance(children, (list, tuple)) else [children]):
            if child is not None and not isinstance(child, (str, int, float)):
                found = find(child, id_)
                if found is not None:
                    return found
    for prefix in ('s', 'c1', 'c2'):
        assert find(gui.app.layout, f'{prefix}-nonortho').value == ['on']
        assert find(gui.app.layout, f'{prefix}-windows').options[0]['disabled'] is False
    S0_P = [o for o in find(gui.app.layout, 's-outputs').options if o['value'] == 'S0_P']
    assert S0_P[0]['disabled'] is False


def test_windows_need_the_non_orthogonal_solver():
    """The orthogonal path ignores n_window, so the GUI never sends it alone."""
    assert 'n_window' not in _p(nonortho=[], windows=ON, n_window=1.51)
    assert _p(nonortho=ON) == _p()
    assert _p(nonortho=ON, windows=ON, n_window=1.6)['n_window'] == 1.6
    # a blank index with the box ticked falls back to borosilicate
    assert _p(nonortho=ON, windows=ON, n_window=None)['n_window'] == gui.N_WINDOW_DEFAULT


def test_doppler_temperature_unties_the_doppler_width():
    p = _p(doppT=80.)
    assert p['Constrain'] is False and p['DoppTemp'] == 80.
    assert _p(doppT='')['Constrain'] is True


# --- single-cell compute ---------------------------------------------------

@pytest.mark.parametrize('windows', [None, ON], ids=['no-windows', 'windows'])
def test_non_orthogonal_compute_matches_the_core(windows):
    fig, store, err = _compute_single(nonortho=ON, windows=windows, n_window=1.51)
    assert err == ''
    extra = dict(use_non_orthogonal=True)
    if windows:
        extra['n_window'] = 1.51
    [expected] = _core(store['X'], gui._pol_efield('Linear (y)', 0.), **extra)
    np.testing.assert_allclose(store['results']['S0'], np.real(expected), rtol=1e-12)


def test_the_switch_changes_the_spectrum():
    """Regime I is where the two propagation models differ (by ~0.47)."""
    _, ortho, _ = _compute_single(nonortho=[])
    _, nonortho, _ = _compute_single(nonortho=ON)
    diff = np.abs(np.subtract(ortho['results']['S0'], nonortho['results']['S0']))
    assert diff.max() > 0.3


def test_S0_P_needs_non_orthogonal_and_is_reported_not_raised():
    _, _, err = _compute_single(outputs=('S0_P',), nonortho=[])
    assert 'S0_P' in err
    _, store, err = _compute_single(outputs=('S0', 'S0_P'), nonortho=ON,
                                    windows=ON, n_window=1.51)
    assert err == ''
    assert np.max(store['results']['S0_P']) <= 1 + 1e-9    # passive vapour column
    # Far wing (+80 to +90 GHz, as in the window review): S0 carries the four
    # window factors, |t1 t2 t3 t4|^2 = 0.844808 at n = 1.51; S0_P does not.
    _, far, err = _compute_single(outputs=('S0', 'S0_P'), nonortho=ON,
                                  windows=ON, n_window=1.51,
                                  det_start=80., det_stop=90., det_N=11)
    assert err == ''
    np.testing.assert_allclose(far['results']['S0'], 0.844808, rtol=1e-4)
    np.testing.assert_allclose(far['results']['S0_P'], 1.0, rtol=1e-4)


def test_doppler_temperature_compute_matches_the_core():
    _, store, err = _compute_single(T=100., doppT=140., btheta=0.)
    assert err == ''
    p = dict(Elem='Rb', Dline='D2', T=100., DoppTemp=140., Constrain=False,
             lcell=2e-3, Bfield=239., GammaBuf=7.9, shift=0., Btheta=0., Bphi=0.,
             rb85frac=72.17)
    [expected] = calculate(np.array(store['X']) * 1e3,
                           gui._pol_efield('Linear (y)', 0.), p, outputs=['S0'])
    np.testing.assert_allclose(store['results']['S0'], np.real(expected), rtol=1e-12)


def test_plot_title_names_the_model():
    """It names the choices that are not the default."""
    fig, _, _ = _compute_single(nonortho=ON, windows=ON, n_window=1.51, doppT=130.)
    title = fig.layout.title.text
    assert 'windows n=1.51' in title and 'orthogonal' not in title
    assert 'Doppler T=130' in title
    fig, _, _ = _compute_single(nonortho=[])
    assert 'orthogonal (ES3)' in fig.layout.title.text
    fig, _, _ = _compute_single()
    assert 'orthogonal' not in fig.layout.title.text


# --- control enabling ------------------------------------------------------

def test_window_controls_follow_the_solver_switch():
    opts, nwin_off, _ = gui._prop_control_state([], [])
    assert opts[0]['disabled'] is True and nwin_off is True
    opts, nwin_off, _ = gui._prop_control_state(ON, [])
    assert opts[0]['disabled'] is False and nwin_off is True
    opts, nwin_off, style = gui._prop_control_state(ON, ON)
    assert nwin_off is False and style == gui.INPUT_STYLE
    # ticked windows stay disabled once the solver is switched off
    _, nwin_off, _ = gui._prop_control_state([], ON)
    assert nwin_off is True


def test_S0_P_is_greyed_out_and_unticked_without_the_solver():
    opts, value = gui.nonortho_output_options([], ['S0', 'S0_P'])
    assert value == ['S0']
    assert {o['value']: o['disabled'] for o in opts}['S0_P'] is True
    opts, value = gui.nonortho_output_options(ON, ['S0', 'S0_P'])
    assert value == ['S0', 'S0_P']
    assert {o['value']: o['disabled'] for o in opts}['S0_P'] is False
    # unticking the only selection falls back to S0, not to nothing
    assert gui.nonortho_output_options([], ['S0_P'])[1] == ['S0']


# --- cascade ---------------------------------------------------------------

def test_cascade_cells_take_their_own_propagation_model():
    X_ghz = np.linspace(-6, 8, 200)
    fig, store, err = gui.compute_cascade(
        1,
        'Rb', 'D2', 120., 2e-3, 239., 7.9, 0., 80., 0., 0., 'Linear (y)',
        'None', 0., 293.15, 72.17, 0.01, 6.73,
        100., 2e-3, 500., 0., 0., 45., 0., 'None', 0., 293.15,
        ['on'], 90., 'none',
        -6., 8., 200, ['S0'],
        ON, ON, 1.51, None,          # cell 1: non-orthogonal with windows
        [], [], None, 130.)[:3]      # cell 2: orthogonal, Doppler T freed
    assert err == ''
    base = dict(Elem='Rb', Dline='D2', rb85frac=72.17, shift=0., Bphi=0.)
    p1 = dict(base, T=120., lcell=2e-3, Bfield=239., GammaBuf=7.9,
              Btheta=np.radians(80.), use_non_orthogonal=True, n_window=1.51)
    p2 = dict(base, T=100., DoppTemp=130., Constrain=False, lcell=2e-3,
              Bfield=500., GammaBuf=0., Btheta=np.radians(45.), use_non_orthogonal=False)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _, (S0_tot,) = calculate_cascade(X_ghz * 1e3, gui._pol_efield('Linear (y)', 0.),
                                         p1, p2, outputs=['S0'],
                                         polariser_angle=np.pi / 2)
    np.testing.assert_allclose(store['total']['S0'], np.real(S0_tot).squeeze(),
                               rtol=1e-10, atol=1e-14)
    title = fig.layout.title.text
    assert 'cell 1 (windows n=1.51)' in title
    assert 'cell 2 (Doppler T=130°C · orthogonal (ES3))' in title


# --- fit panel: data in ----------------------------------------------------

def _upload(text):
    return 'data:text/csv;base64,' + base64.b64encode(text.encode()).decode()


def test_upload_reads_MHz_from_the_header():
    x, y = gui._parse_fit_csv(_upload('Detuning (MHz),T\n-1500,0.9\n2500,0.8\n'))
    np.testing.assert_allclose(x, [-1.5, 2.5])
    x, _ = gui._parse_fit_csv(_upload('Detuning (GHz),T\n-1.5,0.9\n2.5,0.8\n'))
    np.testing.assert_allclose(x, [-1.5, 2.5])
    x, _ = gui._parse_fit_csv(_upload('-1.5,0.9\n2.5,0.8\n'))   # no header: GHz
    np.testing.assert_allclose(x, [-1.5, 2.5])


def test_load_sample_sets_the_cell_it_was_taken_with():
    store, status, elem, dline, lcell, B = gui.load_fit_sample(1)[:6]
    X = np.array(store['X'])
    assert 3 < X.max() < 10 and -10 < X.min() < -3       # GHz, not MHz
    assert len(X) == len(store['Y']) == 3000
    assert (elem, dline, lcell, B) == ('Rb', 'D2', 0.075, 0.0)
    assert 'Sample' in status


def test_sample_then_fit_recovers_the_temperature():
    """The demo path: load the sample, run the default fit (T only)."""
    store = gui.load_fit_sample(1)[0]
    out = _fit(store, ['T'], T=60., lcell=0.075, B=0., btheta=0., gbuf=0.,
               pol='Linear (x)')
    assert out[2] == ''
    assert 37.0 < out[3] < 43.0                     # written-back T
    assert 'RMS residual' in str(out[1])


# --- fit panel: parameters -------------------------------------------------

FIT = {k: v for k, v in SINGLE.items() if not k.startswith('det_')}


def _fit(store, vary, nonortho=gui.NONORTHO_DEFAULT, windows=None, n_window=None, doppT=None,
         data_type='S0', algo='ML', **overrides):
    a = dict(FIT, **overrides)
    return gui.run_fit_single(
        1, store, data_type, algo, vary, None, None,
        a['elem'], a['dline'], a['T'], a['lcell'], a['B'], a['gbuf'], a['shift'],
        a['btheta'], a['bphi'], a['theta0'], a['pol'],
        a['bufgas'], a['bufP'], a['bufTf'], a['rb85frac'], a['K40frac'], a['K41frac'],
        nonortho, windows, n_window, doppT)


def _synthetic(E_in, X_ghz=np.linspace(-6, 8, 250), **p):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        [y] = calculate(X_ghz * 1e3, E_in, p, outputs=['S0'])
    return {'X': X_ghz.tolist(), 'Y': np.real(y).tolist()}


def test_vary_options_grey_out_what_cannot_be_fitted():
    def disabled(**kw):
        return {o['value'] for o in gui._fit_vary_options(**kw) if o['disabled']}
    assert disabled(bufgas='None') == {'DoppTemp', 'BufferPressure'}
    assert disabled(bufgas='None', doppT=80.) == {'BufferPressure'}
    assert 'rb85frac' in disabled(bufgas='None', elem='K')
    assert 'pol_angle' in disabled(bufgas='None', pol='Right CP')
    assert disabled(bufgas='Ne', doppT=80., bufP=5.) == {'GammaBuf'}
    # the model is on only for a gas with data for the line, at P > 0
    assert disabled(bufgas='Ne', doppT=80., bufP=0.) == {'BufferPressure'}
    assert disabled(bufgas='Ne', doppT=80., elem='Na', dline='D2', bufP=5.,
                    ) == {'BufferPressure', 'rb85frac'}


def test_a_disabled_choice_still_ticked_is_dropped_not_fatal():
    store = _synthetic(np.array([0., 1., 0.]), Elem='Rb', Dline='D2', T=120.,
                       lcell=2e-3, Bfield=239., Btheta=np.radians(80.),
                       GammaBuf=7.9)
    out = _fit(store, ['T', 'DoppTemp'], T=118.)        # Doppler T is blank
    assert out[2] == ''
    assert 'Doppler' not in str(out[1])


def test_differential_evolution_searches_only_physical_values(monkeypatch):
    """For a GammaBuf guess of 0 the search window was [-100, 100]
    MHz, and the fit stopped at the first negative width it tried. Every
    window must lie inside param_validation.FIT_LIMITS. (A real DE fit takes
    ~17 s, so fit_data is stubbed to capture what the GUI would hand it.)"""
    from types import SimpleNamespace
    from elecsus.libs.param_validation import FIT_LIMITS
    seen = {}

    def capture(data, p0, vary, p_dict_bounds=None, **kw):
        seen['bounds'] = p_dict_bounds
        return dict(p0), 0.0, SimpleNamespace(best_fit=np.zeros(len(data[0])))

    monkeypatch.setattr(gui, 'fit_data', capture)
    store = _synthetic(np.array([1., 0., 0.]), Elem='Rb', Dline='D2', T=40.,
                       lcell=75e-3, Bfield=0., Btheta=0., GammaBuf=20.)
    out = _fit(store, ['GammaBuf', 'rb85frac', 'T'], algo='differential_evolution',
               T=40., lcell=75e-3, B=0., btheta=0., gbuf=0., rb85frac=100.,
               pol='Linear (x)')
    assert out[2] == ''
    bounds = seen['bounds']
    assert bounds['GammaBuf'] == [0., 100.]          # was [-100, 100]
    assert bounds['rb85frac'] == [50., 100.]         # was [50, 150] before clipping
    for key, (lo, hi) in bounds.items():
        if key in FIT_LIMITS:
            assert FIT_LIMITS[key][0] <= lo < hi <= FIT_LIMITS[key][1], key


def test_fit_through_the_gui_uses_the_non_orthogonal_model():
    """Depends on fit_data honouring use_non_orthogonal (test_fit_propagation)."""
    store = _synthetic(np.array([0., 1., 0.]), Elem='Rb', Dline='D2', T=119.69,
                       lcell=2e-3, Bfield=239., Btheta=np.radians(80.),
                       GammaBuf=7.9, use_non_orthogonal=True)
    out = _fit(store, ['T'], nonortho=ON, T=115.)
    assert out[2] == ''
    assert out[3] == pytest.approx(119.69, abs=0.05)
    assert 'non-orthogonal' in str(out[1])


def test_fit_recovers_and_writes_back_the_doppler_temperature():
    store = _synthetic(np.array([1., 0., 0.]), Elem='Rb', Dline='D2', T=60.,
                       DoppTemp=110., Constrain=False, lcell=0.075, Bfield=0.,
                       X_ghz=np.linspace(-5, 5, 200))
    out = _fit(store, ['DoppTemp'], doppT=80., T=60., lcell=0.075, B=0.,
               btheta=0., gbuf=0., pol='Linear (x)')
    assert out[2] == ''
    assert out[11] == pytest.approx(110., abs=0.5)      # written back


def test_fit_recovers_and_writes_back_the_buffer_gas_pressure():
    """With a gas selected the pressure can vary, and so can the shift, as an
    offset of the frequency axis on top of the model's pressure shift."""
    store = _synthetic(np.array([1., 0., 0.]), Elem='Rb', Dline='D2', T=60.,
                       lcell=0.075, Bfield=0., BufferGas='N2', BufferPressure=10.,
                       shift=30.)
    out = _fit(store, ['BufferPressure', 'shift'], T=60., lcell=0.075, B=0.,
               btheta=0., gbuf=0., pol='Linear (x)', bufgas='N2', bufP=6.)
    assert out[2] == ''
    assert out[15] == pytest.approx(10., abs=0.05)      # written back
    assert out[6] == pytest.approx(30., abs=1.)
    assert 'Buffer P' in str(out[1])


def test_differential_evolution_searches_positive_buffer_pressures(monkeypatch):
    """The search window around the pressure guess stops at 0 Torr, as for
    the other parameters in test_differential_evolution_searches_only_physical_values."""
    from types import SimpleNamespace
    seen = {}

    def capture(data, p0, vary, p_dict_bounds=None, **kw):
        seen['bounds'] = p_dict_bounds
        return dict(p0), 0.0, SimpleNamespace(best_fit=np.zeros(len(data[0])))

    monkeypatch.setattr(gui, 'fit_data', capture)
    store = _synthetic(np.array([1., 0., 0.]), Elem='Rb', Dline='D2', T=60.,
                       lcell=0.075, Bfield=0., BufferGas='N2', BufferPressure=10.)
    out = _fit(store, ['BufferPressure'], algo='differential_evolution', T=60.,
               lcell=0.075, B=0., btheta=0., gbuf=0., pol='Linear (x)',
               bufgas='N2', bufP=6.)
    assert out[2] == ''
    assert seen['bounds'] == {'BufferPressure': [0., 16.]}


@pytest.mark.parametrize('true_deg,start_deg', [(30., 15.), (-35., -20.)])
def test_fit_recovers_a_linear_polarisation_angle(true_deg, start_deg):
    """Voigt geometry with B rotated in the transverse plane, so transmission
    depends on the sign of the angle as well as its size.

    Transmission is symmetric about B's azimuth, so the start and the true
    angle must lie on the same side of it: B is at 340 (-20) degrees, the
    GUI's range being 0-360. (The test used Bphi = 20 degrees when the code
    put B at -Bphi.)"""
    p = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=1500.,
             Btheta=np.pi / 2, Bphi=np.radians(340.), GammaBuf=0.)
    store = _synthetic(gui._pol_efield('Linear (x)', true_deg), **p)
    out = _fit(store, ['pol_angle'], T=80., lcell=2e-3, B=1500., btheta=90.,
               bphi=340., gbuf=0., pol='Linear (x)', theta0=start_deg)
    assert out[2] == ''
    assert out[13] == 'Linear (x)'
    assert out[14] == pytest.approx(true_deg, abs=0.2)
    assert 'Pol. θ₀ [deg]' in str(out[1])


def test_right_circular_fit_uses_right_circular_light():
    """fit_data bounds the E_y phase to [0, pi); the old conversion passed -pi/2
    for right circular light, lmfit clamped it to 0, and the fit used linear
    45-degree light. Faraday geometry, where the two differ strongly."""
    p = dict(Elem='Rb', Dline='D2', T=70., lcell=0.075, Bfield=600., Btheta=0.,
             GammaBuf=0.)
    E_rcp = gui._pol_efield('Right CP', 0.)
    store = _synthetic(E_rcp, X_ghz=np.linspace(-6, 6, 200), **p)
    out = _fit(store, ['T'], T=66., lcell=0.075, B=600., btheta=0., gbuf=0.,
               pol='Right CP')
    assert out[2] == ''
    assert out[3] == pytest.approx(70., abs=0.05)


@pytest.mark.parametrize('pol,theta0', [
    ('Linear (x)', 0.), ('Linear (x)', -30.), ('Linear (y)', 0.),
    ('Linear 45°', 0.), ('Left CP', 0.), ('Right CP', 0.)])
def test_fit_field_is_the_same_polarisation_state(pol, theta0):
    """_pol_fit_setup's [E_x, [E_y, phase]] must describe the GUI's field up
    to a global phase, with the phase inside fit_data's [0, pi) bound."""
    E = np.array(gui._pol_efield(pol, theta0), dtype=complex)[:2]
    (E_x, (E_y, phase)), _ = gui._pol_fit_setup(gui._pol_efield(pol, theta0))
    assert 0.0 <= phase < np.pi
    F = np.array([E_x, E_y * np.exp(1j * phase)])
    overlap = abs(np.vdot(E, F)) / (np.linalg.norm(E) * np.linalg.norm(F))
    assert overlap == pytest.approx(1.0, abs=1e-12)
