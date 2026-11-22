"""Callback tests for the Dash GUI.

test_gui_angles.py and test_gui_cascade.py cover the helpers and the cascade
compute path. This module covers the rest of the registered callbacks -- the
single-cell compute, the CSV exports, the upload parser and the small
enable/disable and visibility callbacks -- which between them were most of the
GUI's uncovered lines.

Dash callbacks are ordinary module-level functions, so they are called directly
here: no browser, no HTTP layer, no Selenium. What is being pinned is that the
GUI stays a thin wrapper over the core API and reports failures instead of
dying, since a user meets a traceback in the browser as an empty panel.
"""
import base64

import numpy as np
import pytest

from elecsus import elecsus_dash_gui as gui
from elecsus.elecsus_methods import calculate

# The single-cell compute callback takes its inputs positionally in the order
# declared in its @app.callback State list.
SINGLE_ARGS = dict(
    elem='Rb', dline='D2', T=60., lcell=0.075, B=200., gbuf=0., shift=0.,
    btheta=0., bphi=0., theta0=0., pol='Linear (x)',
    bufgas='None', bufP=0., bufTf=293.15,
    rb85frac=72.17, K40frac=0.01, K41frac=6.73,
    det_start=-6., det_stop=6., det_N=400,
)


def _compute_single(outputs=('S0',), **overrides):
    a = dict(SINGLE_ARGS, **overrides)
    return gui.compute_single(
        1, a['elem'], a['dline'], a['T'], a['lcell'], a['B'], a['gbuf'],
        a['shift'], a['btheta'], a['bphi'], a['theta0'], a['pol'],
        a['bufgas'], a['bufP'], a['bufTf'],
        a['rb85frac'], a['K40frac'], a['K41frac'],
        a['det_start'], a['det_stop'], a['det_N'], list(outputs))[:3]


# --- single-cell compute ------------------------------------------------

def test_compute_single_matches_the_core_api():
    """The GUI must be a thin wrapper: same numbers as calculate()."""
    fig, store, err = _compute_single()
    assert err == ''
    X = np.array(store['X']) * 1e3          # the store keeps GHz
    p = dict(Elem='Rb', Dline='D2', T=60., lcell=0.075, Bfield=200.,
             GammaBuf=0., shift=0., Btheta=0., Bphi=0.,
             rb85frac=72.17, K40frac=0.01, K41frac=6.73)
    [expected] = calculate(X, gui._pol_efield('Linear (x)', 0.), p, outputs=['S0'])
    np.testing.assert_allclose(store['results']['S0'], np.real(expected), rtol=1e-12)
    assert fig.data, 'no traces on the returned figure'


def test_compute_single_returns_several_outputs():
    _, store, err = _compute_single(outputs=('S0', 'S1', 'Ix'))
    assert err == ''
    assert set(store['results']) == {'S0', 'S1', 'Ix'}
    assert all(len(v) == SINGLE_ARGS['det_N'] for v in store['results'].values())


def test_compute_single_defaults_to_S0_when_no_output_is_ticked():
    _, store, err = _compute_single(outputs=[])
    assert err == '' and list(store['results']) == ['S0']


def test_compute_single_reports_a_bad_parameter_instead_of_crashing():
    """A negative cell length now raises in the core (parameter validation).
    The GUI must surface that as an error string, not propagate the exception:
    in the browser an uncaught callback exception is an empty panel."""
    fig, store, err = _compute_single(lcell=-0.075)
    assert 'lcell' in err
    assert fig is not None


def test_compute_single_reports_an_unsupported_element():
    _, _, err = _compute_single(elem='Xe')
    assert 'Xe' in err or 'atomic data' in err


# --- CSV exports --------------------------------------------------------

def test_export_single_round_trips_the_store():
    _, store, _ = _compute_single(outputs=('S0', 'S1'))
    out = gui.export_single(1, store)
    lines = out['content'].splitlines()
    assert out['filename'].endswith('.csv')
    assert lines[0] == 'Detuning_GHz,S0,S1'
    assert len(lines) == 1 + len(store['X'])
    first = lines[1].split(',')
    assert float(first[0]) == pytest.approx(store['X'][0])
    assert float(first[1]) == pytest.approx(store['results']['S0'][0], abs=1e-8)


def test_export_single_without_a_computation_does_nothing():
    assert gui.export_single(1, None) is gui.no_update


def test_export_fit_includes_the_residual_column():
    store = {'X': [-1., 0., 1.], 'data': [0.9, 0.5, 0.9],
             'fit': [0.89, 0.52, 0.91], 'type': 'S1'}
    out = gui.export_fit(1, store)
    lines = out['content'].splitlines()
    assert lines[0] == 'Detuning_GHz,S1_data,S1_fit,residual'
    residual = float(lines[1].split(',')[3])
    assert residual == pytest.approx(0.9 - 0.89, abs=1e-8)


def test_export_cascade_labels_cell1_and_total_columns():
    store = {'X': [-1., 0., 1.],
             'cell1': {'S0': [0.9, 0.5, 0.9]},
             'total': {'S0': [0.8, 0.3, 0.8]}}
    out = gui.export_cascade(1, store)
    lines = out['content'].splitlines()
    assert lines[0] == 'Detuning_GHz,Cell1_S0,Total_S0'
    assert len(lines) == 4


@pytest.mark.parametrize('exporter', ['export_single', 'export_fit', 'export_cascade'])
def test_exports_are_inert_without_data(exporter):
    assert getattr(gui, exporter)(1, None) is gui.no_update


# --- upload parsing -----------------------------------------------------

def _data_uri(text):
    return 'data:text/csv;base64,' + base64.b64encode(text.encode()).decode()


def test_parse_fit_upload_accepts_a_plain_two_column_csv():
    store, status = gui.parse_fit_upload(_data_uri('-1.0,0.9\n0.0,0.5\n1.0,0.9\n'),
                                         'scan.csv')
    assert store['X'] == [-1.0, 0.0, 1.0]
    assert store['Y'] == [0.9, 0.5, 0.9]
    assert '3 points' in status


def test_parse_fit_upload_accepts_a_header_row():
    store, status = gui.parse_fit_upload(
        _data_uri('detuning,signal\n-1.0,0.9\n0.0,0.5\n'), 'with_header.csv')
    assert store['X'] == [-1.0, 0.0]
    assert 'points loaded' in status


def test_parse_fit_upload_reports_an_unparseable_file():
    store, status = gui.parse_fit_upload(_data_uri('not a csv at all\n'), 'junk.txt')
    assert store is gui.no_update
    assert 'Failed to parse junk.txt' in status


def test_parse_fit_upload_rejects_a_single_column():
    store, status = gui.parse_fit_upload(_data_uri('1.0\n2.0\n3.0\n'), 'one_col.csv')
    assert store is gui.no_update
    assert 'two columns' in status


# --- small state callbacks ---------------------------------------------

@pytest.mark.parametrize('callback', ['grey_gbuf_single', 'grey_gbuf_c1', 'grey_gbuf_c2'])
@pytest.mark.parametrize('bufgas, expected_disabled',
                         [('None', False), (None, False), ('N2', True), ('Ne', True)])
def test_gammabuf_is_greyed_out_when_a_buffer_gas_is_selected(callback, bufgas,
                                                              expected_disabled):
    """GammaBuf is overridden by the buffer-gas model, so the GUI must not
    invite the user to set both (see the manual's precedence section)."""
    disabled, style = getattr(gui, callback)(bufgas, 'Rb', 'D2', 5.)
    assert disabled is expected_disabled
    assert isinstance(style, dict)


@pytest.mark.parametrize('callback', ['grey_gbuf_single', 'grey_gbuf_c1', 'grey_gbuf_c2'])
@pytest.mark.parametrize('bufgas, elem, dline, bufP', [
    ('N2', 'Rb', 'D2', 0.),       # a gas but no pressure: the core uses GammaBuf
    ('N2', 'Rb', 'D2', None),     # pressure box empty
    ('N2', 'Na', 'D2', 5.),       # no data for Na: the core uses GammaBuf
    ('Ne', 'Na', 'D1', 5.),
])
def test_gammabuf_stays_open_when_the_model_is_not_used(callback, bufgas, elem, dline, bufP):
    """The box was greyed out for any gas, though the core uses the
    static GammaBuf at P = 0 and for a gas it has no data for."""
    disabled, style = getattr(gui, callback)(bufgas, elem, dline, bufP)
    assert disabled is False
    assert style == gui.INPUT_STYLE


@pytest.mark.parametrize('callback', ['iso_visibility_s', 'iso_visibility_c1'])
@pytest.mark.parametrize('elem, rb_shown, k_shown', [
    ('Rb', True, False), ('K', False, True),
    ('Cs', False, False), ('Na', False, False),
])
def test_isotope_rows_follow_the_element(callback, elem, rb_shown, k_shown):
    rb_style, k_style = getattr(gui, callback)(elem)
    assert (rb_style['display'] == 'block') is rb_shown
    assert (k_style['display'] == 'block') is k_shown


def test_fit_vary_options_track_the_buffer_gas_selection():
    """With a buffer gas active the model's broadening replaces GammaBuf, and
    fit_data() raises if GammaBuf is also marked to vary; without one the
    buffer pressure has no effect, and fit_data() raises if it is marked. The
    GUI keeps both in the list but disables them, so the user sees the option
    and cannot pick a combination the core rejects. The shift stays: the
    model's pressure shift is added to it."""
    def by_key(bufgas):
        return {o['value']: o for o in
                gui.update_fit_vary_options(bufgas, 'Rb', None, 'Linear (x)', 'D2', 5.)}

    without, with_gas = by_key('None'), by_key('N2')
    assert 'T' in without and 'T' in with_gas
    assert without['GammaBuf']['disabled'] is False
    assert with_gas['GammaBuf']['disabled'] is True
    assert without['BufferPressure']['disabled'] is True
    assert with_gas['BufferPressure']['disabled'] is False
    assert without['shift']['disabled'] is False
    assert with_gas['shift']['disabled'] is False
    assert with_gas['T']['disabled'] is False


def test_show_error_opens_the_toast_only_when_there_is_a_message():
    assert gui.show_error('boom') == ('boom', True)
    assert gui.show_error('') == ('', False)
    assert gui.show_error(None) == ('', False)
