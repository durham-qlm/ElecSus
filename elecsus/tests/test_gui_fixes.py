"""GUI fixes made before the beta: a failed run keeps nothing to export; a fit
that did not run is not reported as one; the core's warnings reach the page
(GammaBuf is greyed out only when the model replaces it); the pinned signs and
labels, and the values that reach the export; a fit range with one end; CSV
precision; the upload parser; Load sample spectrum; and the polariser angle's
label.

The wiring tests are in test_gui_wiring.py.
"""
import base64
import io
import contextlib
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from elecsus import elecsus_dash_gui as gui
from elecsus.elecsus_methods import calculate
from elecsus.tests.test_gui_wiring import post_callback, SINGLE


# compute_single's positional arguments, as in test_gui_inputs
ARGS = dict(
    elem='Rb', dline='D2', T=60., lcell=0.075, B=200., gbuf=0., shift=0.,
    btheta=0., bphi=0., theta0=0., pol='Linear (x)',
    bufgas='None', bufP=0., bufTf=293.15,
    rb85frac=72.17, K40frac=0.01, K41frac=6.73,
    det_start=-3., det_stop=3., det_N=40,
)


def _compute(outputs=('S0',), **over):
    a = dict(ARGS, **over)
    return gui.compute_single(1, *[a[k] for k in ARGS], list(outputs))


CASCADE = dict(
    e1='Rb', d1='D2', T1=60., l1=0.075, B1=200., g1=0., sh1=0., bt1=0., bp1=0.,
    th1=0., pol1='Linear (x)', bg1='None', bP1=0., bTf1=293.15,
    rb85=72.17, k40=0.01, k41=6.73,
    T2=80., l2=0.075, B2=50., g2=0., sh2=0., bt2=0., bp2=0.,
    bg2='None', bP2=0., bTf2=293.15,
    use_pol=['on'], pol_angle=90., reset_pol='none',
    det_start=-3., det_stop=3., det_N=40, outputs=['S0'])


def _cascade(**over):
    return gui.compute_cascade(1, *dict(CASCADE, **over).values())


SAMPLE = gui.load_fit_sample(1)[0]


def _fit(store=SAMPLE, vary=('T',), fit_min=None, fit_max=None, algo='ML', **over):
    a = dict(elem='Rb', dline='D2', T=40., lcell=0.075, B=0., gbuf=0., shift=0.,
             btheta=0., bphi=0., theta0=0., pol='Linear (x)', bufgas='None', bufP=0.,
             bufTf=293.15, rb85frac=72.17, K40frac=0.01, K41frac=6.73)
    a.update(over)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('ignore')
        return gui.run_fit_single(1, store, 'S0', algo, list(vary), fit_min, fit_max,
                                  *[a[k] for k in a])


def _text(component):
    return str(component)


# --- a failed run leaves nothing to export -------------------------------

def test_a_failed_compute_empties_the_store_so_export_writes_nothing():
    _, store, err, _ = _compute(lcell=-1.)
    assert err and store is None
    assert gui.export_single(1, store) is gui.no_update


def test_a_failed_cascade_empties_the_store_so_export_writes_nothing():
    _, store, err, _ = _cascade(l2=-1.)
    assert err and store is None
    assert gui.export_cascade(1, store) is gui.no_update


def test_a_failed_fit_empties_the_result_store_so_export_writes_nothing():
    good = _fit()
    assert good[2] == '' and good[10]
    out = _fit(T=None)
    assert out[2] and out[10] is None
    assert gui.export_fit(1, out[10]) is gui.no_update
    out = _fit(store=None)                        # no data at all
    assert out[10] is None


# --- a fit that did not run says so ----------------------------------------

def test_a_fit_that_cannot_move_from_its_start_says_so():
    """At 400 C the model is 0 everywhere: lmfit stops after 2 evaluations with
    success=True and no error bar."""
    out = _fit(T=400.)
    assert out[2] == ''
    assert 'did not move' in _text(out[1])
    assert out[3] == 400.


def test_an_undetermined_parameter_is_flagged_and_not_written_back():
    """Bphi does nothing in zero field: lmfit reports a value of 1e8 degrees."""
    out = _fit(vary=('T', 'Bphi'))
    assert out[2] == ''
    assert 'Bφ' in _text(out[1]) and 'not determined' in _text(out[1])
    assert out[9] is gui.no_update                     # the Bφ box keeps what was typed
    assert out[3] == pytest.approx(40.1, abs=0.2)      # T, which was determined, is written


def test_a_missing_error_estimate_is_flagged():
    out = _fit(vary=('T', 'Btheta'))
    assert 'no error estimate' in _text(out[1])


def test_a_good_fit_has_no_warning_lines():
    out = _fit()
    assert out[2] == '' and '⚠' not in _text(out[1])
    assert out[3] == pytest.approx(40.1, abs=0.2)


def test_a_fit_with_nothing_ticked_is_refused():
    out = _fit(vary=())
    assert 'Nothing to vary' in out[2]
    assert out[10] is None


def test_a_fit_with_only_greyed_out_choices_ticked_says_which():
    out = _fit(vary=('DoppTemp',))                    # Doppler T is blank
    assert 'Nothing to vary' in out[2] and 'Doppler T' in out[2]


def test_the_stderr_of_the_field_angles_is_shown_in_degrees(monkeypatch):
    """The fit works in radians; the box and the table are in degrees."""
    def stub(data, p0, vary, **kw):
        params = {k: SimpleNamespace(value=v, init_value=v - 0.1, stderr=0.01)
                  for k, v in (('Btheta', 0.5), ('Bphi', 0.25))}
        best = dict(p0, Btheta=0.5, Bphi=0.25)
        return best, 0.01, SimpleNamespace(best_fit=np.zeros(len(data[0])), params=params,
                                           success=True, message='', errorbars=True)
    monkeypatch.setattr(gui, 'fit_data', stub)
    out = _fit(vary=('Btheta', 'Bphi'))
    text = _text(out[1])
    assert '± 0.573' in text                           # 0.01 rad
    assert out[8] == pytest.approx(np.degrees(0.5))
    assert out[9] == pytest.approx(np.degrees(0.25))


# --- the core's warnings reach the page ------------------------------------

def test_a_warning_from_the_core_is_returned_for_the_page():
    with warnings.catch_warnings(record=True) as leaked:
        warnings.simplefilter('always')
        fig, store, err, notes = _compute(bufgas='N2', bufP=0., gbuf=20.)
    assert err == '' and store is not None
    assert len(notes) == 1 and "BufferGas='N2'" in notes[0] and 'BufferPressure' in notes[0]
    assert not leaked, 'the warning also went to the server'


def test_the_no_data_fallback_warning_is_returned():
    _, _, err, notes = _compute(elem='Na', bufgas='N2', bufP=5., gbuf=20., T=100.)
    assert err == ''
    assert any('No buffer gas data' in m for m in notes)


def test_a_clean_compute_returns_no_warnings():
    assert _compute()[3] == []


def test_warnings_that_come_before_an_error_are_still_shown():
    _, store, err, notes = _compute(bufgas='N2', bufP=0., gbuf=20., lcell=-1.)
    assert store is None and err
    assert isinstance(notes, list)


def test_the_cascade_and_the_fit_return_the_cores_warnings_too():
    notes = _cascade(bg1='N2', bP1=0., g1=20.)[3]
    assert any("BufferGas='N2'" in m for m in notes)
    notes = _fit(bufgas='N2', bufP=0., gbuf=20.)[-1]
    assert any("BufferGas='N2'" in m for m in notes)
    assert _fit()[-1] == []


def test_the_warning_toast_follows_the_warning_store():
    children, is_open = gui.show_warnings(['one', 'two'])
    assert is_open and 'one' in _text(children) and 'two' in _text(children)
    assert gui.show_warnings([]) == ('', False)
    assert gui.show_warnings(None) == ('', False)


def test_the_page_receives_the_warning_through_dash():
    values = dict(SINGLE, **{'s-bufgas': 'N2', 's-bufP': 0., 's-gbuf': 20.})
    got = post_callback('s-plot', [('s-compute', 'n_clicks')], values)
    assert any("BufferGas='N2'" in m for m in got[('warning-store', 'data')])
    got = post_callback('s-plot', [('s-compute', 'n_clicks')], SINGLE)
    assert got[('warning-store', 'data')] == []


def test_a_blank_gammabuf_is_refused_when_the_model_is_not_used():
    """Blank means 'the gas supplies it' only while the model is on."""
    _, _, err, _ = _compute(bufgas='N2', bufP=0., gbuf=None)
    assert 'GammaBuf' in err
    _, _, err, _ = _compute(elem='Na', bufgas='N2', bufP=5., gbuf=None)
    assert 'GammaBuf' in err
    _, _, err, _ = _compute(bufgas='N2', bufP=5., gbuf=None)
    assert err == ''


def test_the_gammabuf_boxes_follow_gas_pressure_and_element_through_dash():
    spec = {k: s for k, s in gui.app.callback_map.items() if k.startswith('..s-gbuf.')}
    [(key, s)] = spec.items()
    assert {i['id'] for i in s['inputs']} == {'s-bufgas', 's-elem', 's-dline', 's-bufP'}


# --- the signs and the labels the GUI pins -----------------------------------

def test_the_circular_presets_are_the_documented_states():
    np.testing.assert_allclose(gui._pol_efield('Left CP', 0.), np.array([1, 1j, 0]) / np.sqrt(2))
    np.testing.assert_allclose(gui._pol_efield('Right CP', 0.), np.array([1, -1j, 0]) / np.sqrt(2))


def test_the_sigma_labels_are_the_physical_ones():
    """TransLeft is the [1, +i] transmission, sigma+ for B along +z."""
    assert 'σ⁺' in gui.OUTPUT_LABELS['TransLeft'] and 'σ⁻' in gui.OUTPUT_LABELS['TransRight']
    assert 'σ⁺' in gui.OUTPUT_LABELS['alphaPlus'] and 'σ⁻' in gui.OUTPUT_LABELS['alphaMinus']
    assert 'σ⁺' in gui.OUTPUT_LABELS['nPlus'] and 'σ⁻' in gui.OUTPUT_LABELS['nMinus']
    # and the physics behind the label: in B > 0 along +z, [1, +i] light is
    # absorbed on the blue-shifted lines (sigma+), the high-frequency side
    X = np.linspace(-6e3, 6e3, 400)
    p = dict(Elem='Rb', Dline='D2', T=60., lcell=0.075, Bfield=600., Btheta=0.)
    [left] = calculate(X, gui._pol_efield('Left CP', 0.), p, outputs=['S0'])
    [right] = calculate(X, gui._pol_efield('Right CP', 0.), p, outputs=['S0'])
    assert np.mean(np.real(left)[X > 0]) < np.mean(np.real(right)[X > 0])


def test_the_isotope_fractions_of_potassium_reach_the_core():
    _, a, _, _ = _compute(elem='K', dline='D1', T=150., K41frac=6.73)
    _, b, _, _ = _compute(elem='K', dline='D1', T=150., K41frac=20.)
    assert not np.allclose(a['results']['S0'], b['results']['S0'], atol=1e-4)
    _, c, _, _ = _compute(elem='K', dline='D1', T=150., K40frac=0.5)
    assert not np.allclose(a['results']['S0'], c['results']['S0'], atol=1e-6)


def test_the_cascade_export_holds_each_columns_values():
    _, store, err, _ = _cascade(outputs=['S0', 'S3'])
    assert err == ''
    lines = gui.export_cascade(1, store)['content'].splitlines()
    # the last two columns are the lower panel: cell 2 on its own
    assert lines[0] == ('Detuning_GHz,Cell1_S0,Cell1_S3,Total_S0,Total_S3,'
                        'Into_Cell2,Cell2_alone_S0')
    rows = np.array([[float(v) for v in l.split(',')] for l in lines[1:]])
    np.testing.assert_allclose(rows[:, 0], store['X'], atol=1e-6)
    np.testing.assert_allclose(rows[:, 1], store['cell1']['S0'], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(rows[:, 2], store['cell1']['S3'], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(rows[:, 3], store['total']['S0'], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(rows[:, 4], store['total']['S3'], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(rows[:, 5], store['cell2']['into'], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(rows[:, 6], store['cell2']['transmission'],
                               rtol=1e-8, atol=1e-12)
    assert not np.allclose(rows[:, 1], rows[:, 3])          # cell 1 differs from the total


# --- a fit range with one end ----------------------------------------------

def _n_fitted(fit_min, fit_max):
    return len(_fit(fit_min=fit_min, fit_max=fit_max)[10]['X'])


def test_each_end_of_the_fit_range_works_alone():
    X = np.array(SAMPLE['X'])
    assert _n_fitted(None, None) == len(X) == 3000
    assert _n_fitted(-1.0, None) == np.sum(X >= -1.0)
    assert _n_fitted(None, 1.0) == np.sum(X <= 1.0)
    assert _n_fitted(-1.0, 1.0) == np.sum((X >= -1.0) & (X <= 1.0))
    assert 0 < _n_fitted(-1.0, None) < 3000 and 0 < _n_fitted(None, 1.0) < 3000


def test_the_results_say_how_many_points_were_fitted():
    assert 'all 3000 points fitted' in _text(_fit()[1])
    assert 'of 3000 points fitted' in _text(_fit(fit_min=-1.0)[1])


def test_an_empty_fit_range_is_refused():
    assert 'Not enough data points' in _fit(fit_min=50.)[2]


# --- the CSV keeps small values ---------------------------------------------

def test_the_exports_keep_values_below_5e_9():
    store = {'X': [0.0, 0.5, 1.0], 'results': {'S0': [1e-12, 9.5e-9, 0.9971138]}}
    rows = gui.export_single(1, store)['content'].splitlines()[1:]
    got = [float(r.split(',')[1]) for r in rows]
    np.testing.assert_allclose(got, store['results']['S0'], rtol=1e-6)

    fit = {'X': [0.0, 0.5], 'data': [3e-11, 0.5], 'fit': [1e-11, 0.5], 'type': 'S0'}
    rows = gui.export_fit(1, fit)['content'].splitlines()[1:]
    cols = np.array([[float(v) for v in r.split(',')] for r in rows])
    np.testing.assert_allclose(cols[:, 1:], [[3e-11, 1e-11, 2e-11], [0.5, 0.5, 0.0]],
                               rtol=1e-6, atol=1e-15)

    cas = {'X': [0.0], 'cell1': {'S0': [2e-13]}, 'total': {'S0': [4e-15]}}
    row = gui.export_cascade(1, cas)['content'].splitlines()[1]
    np.testing.assert_allclose([float(v) for v in row.split(',')[1:]], [2e-13, 4e-15], rtol=1e-6)


def test_a_hot_spectrum_exports_with_its_line_centre():
    _, store, err, _ = _compute(T=150., det_start=-1.5, det_stop=1.5, det_N=11)
    assert err == ''
    got = [float(l.split(',')[1]) for l in gui.export_single(1, store)['content'].splitlines()[1:]]
    exact = np.array(store['results']['S0'])
    assert 0 < exact.max() < 1e-5 and exact.max() > 0
    np.testing.assert_allclose(got, exact, rtol=1e-6, atol=0)


# --- the upload parser ------------------------------------------------------

def _uri(text):
    return 'data:text/csv;base64,' + base64.b64encode(text.encode()).decode()


@pytest.mark.parametrize('name, text, x', [
    ('whitespace and a header', 'Detuning (GHz) T\n-1.0 0.9\n0.0 0.5\n1.0 0.9\n', [-1, 0, 1]),
    ('tabs and a header (MHz)', 'Detuning (MHz)\tT\n-1000\t0.9\n0\t0.5\n1000\t0.9\n', [-1, 0, 1]),
    ('two header rows', 'run 12\nDetuning (MHz),T\n-1000,0.9\n0,0.5\n1000,0.9\n', [-1, 0, 1]),
    ('GHz, a MHz in another column', 'Detuning (GHz),T (10 MHz BW)\n-1.0,0.9\n0.0,0.5\n1.0,0.9\n',
     [-1, 0, 1]),
    ('whitespace, MHz, longer names', 'Detuning (MHz) Transmission (a.u.)\n-1000 0.9\n1000 0.8\n',
     [-1, 1]),
    ('comment header', '# Detuning (MHz), T\n-1000,0.9\n1000,0.8\n', [-1, 1]),
    ('byte-order mark and CRLF', '﻿Detuning (MHz),T\r\n-1000,0.9\r\n1000,0.8\r\n', [-1, 1]),
    ('no header: GHz', '-1.5,0.9\n2.5,0.8\n', [-1.5, 2.5]),
    ('a third column is ignored', 'Detuning (GHz),T,err\n-1,0.9,0.01\n1,0.8,0.01\n', [-1, 1]),
])
def test_the_upload_parser_accepts(name, text, x):
    X, Y = gui._parse_fit_csv(_uri(text))
    np.testing.assert_allclose(X, x)
    assert len(Y) == len(x)


def test_the_status_line_names_the_unit_and_the_range():
    _, status = gui.parse_fit_upload(_uri('Detuning (MHz),T\n-1500,0.9\n2500,0.8\n'), 's.csv')
    assert 'points loaded' in status and 'MHz' in status and '-1.5 to 2.5 GHz' in status
    _, status = gui.parse_fit_upload(_uri('-1.5,0.9\n2.5,0.8\n'), 's.csv')
    assert 'GHz' in status and 'default' in status


@pytest.mark.parametrize('text, message', [
    ('not,a,csv\n', 'Could not parse'),
    ('1.0\n2.0\n', 'two columns'),
    ('1,2\nfoo\n3,4\n', 'line 2'),
    ('1,2\n3,4,5\n', 'same number of columns'),
])
def test_the_upload_parser_still_refuses_what_is_not_two_columns(text, message):
    with pytest.raises(ValueError, match=message):
        gui._parse_fit_csv(_uri(text))


# --- Load sample spectrum -----------------------------------------------------

def _walk(c):
    yield c
    kids = getattr(c, 'children', None)
    for k in (kids if isinstance(kids, (list, tuple)) else [kids]):
        if hasattr(k, 'to_plotly_json'):
            yield from _walk(k)


def test_load_sample_resets_the_other_boxes_to_the_fresh_page():
    defaults = {c.id: c.value for c in _walk(gui.app.layout)
                if getattr(c, 'id', None) and hasattr(c, 'value') and c.value is not None
                or getattr(c, 'id', None) in ('s-doppT', 's-fit-min', 's-fit-max')}
    for id_, value in gui.SAMPLE_PANEL:
        assert defaults.get(id_) == value, id_


def test_run_fit_after_load_sample_gives_the_documented_result_on_a_used_panel():
    """A panel left with a GammaBuf, a shift, a gas and angles from earlier."""
    dirty = dict(SINGLE, **{
        's-T': 40., 's-gbuf': 30., 's-shift': 60., 's-btheta': 20., 's-bphi': 40.,
        's-bufgas': 'Ne', 's-bufP': 5., 's-doppT': 70., 's-pol': 'Left CP',
        's-theta0': 33., 's-nonortho': [], 's-windows': ['on'], 's-rb85frac': 50.,
        's-fit-datatype': 'S1', 's-fit-vary': ['T', 'shift'], 's-fit-min': 2., 's-fit-max': 3.,
        's-fit-algo': 'ML', 's-fit-data-store': None})
    loaded = post_callback('s-fit-data-store', [('s-fit-sample', 'n_clicks')], {})
    panel = dict(dirty)
    for (id_, prop), value in loaded.items():
        if prop == 'value':
            panel[id_] = value
    panel['s-fit-data-store'] = loaded[('s-fit-data-store', 'data')]
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('ignore')
        got = post_callback('s-fit-plot', [('s-fit-run', 'n_clicks')], panel)
    assert got[('error-store', 'data')] == ''
    assert got[('s-T', 'value')] == pytest.approx(40.10, abs=0.05)
    assert 'RMS residual: 0.0085' in _text(got[('s-fit-results', 'children')])
    assert panel['s-gbuf'] == 0 and panel['s-bufgas'] == 'None' and panel['s-pol'] == 'Linear (x)'


# --- the polariser angle is measured from x --------------------------------------

def test_the_polariser_hint_does_not_claim_it_is_always_crossed():
    [field] = [c for c in _walk(gui.app.layout)
               if getattr(c, 'id', None) == 'cas-pol-angle']
    text = ' '.join(_text(c) for c in _walk(gui.app.layout)
                    if 'crossed' in _text(getattr(c, 'children', '')))
    assert 'from x' in text and 'x-polarised' in text
    assert '90 = crossed (MOF)' not in text and field.value == 90.0


def test_ninety_degrees_is_crossed_only_for_x_polarised_light():
    """The behaviour the label describes (the angle stays absolute)."""
    def after_both(pol):
        _, st, err, _ = _cascade(pol1=pol, T1=60., T2=60., B1=0., B2=0.,
                                 det_start=-3., det_stop=3., det_N=11)
        assert err == ''
        return np.array(st['total']['S0'])[5]
    assert after_both('Linear (x)') < 1e-6          # crossed
    assert after_both('Linear (y)') > 0.5           # the polariser passes y-light


# --- the header gives the version, and says a pre-release is not final

@pytest.mark.parametrize('version, note', [
    ('4.0.0b1', 'v4.0.0b1 · beta, not the final release'),
    ('4.0.0a3', 'v4.0.0a3 · alpha, not the final release'),
    ('4.0.0rc2', 'v4.0.0rc2 · release candidate, not the final release'),
    ('4.1.0.dev3', 'v4.1.0.dev3 · development version, not the final release'),
    ('4.0.0', 'v4.0.0'), ('4.0.1.post1', 'v4.0.1.post1'),
    ('0+unknown', ''),                  # a source tree that was never installed
])
def test_the_release_note_follows_the_version(version, note):
    assert gui._release_note(version) == note


def test_the_header_shows_the_installed_version():
    [line] = [c for c in _walk(gui.app.layout) if getattr(c, 'id', None) == 'release-note']
    assert line.children == gui._release_note(gui.__version__)
