"""The GUI's wiring: the order of each callback's Input/State list.

Every other GUI test calls a callback as a plain function, in the order the
*function* declares its arguments. The order the *decorator* declares its
State list in, which is what Dash feeds the function, was never run. Swap two
`State(...)` lines, or add an input in the wrong place, and the page breaks
while every such test passes.

These tests post to Dash's own endpoint (`/_dash-update-component`), as the
browser does. Each box gets a value of its own, so a transposed pair changes
the result; the expected numbers come from the core API called by name.

`post_callback` insists that the test names every State the callback
declares, so a new input cannot be added without a value for it here.
"""
import base64
import warnings

import numpy as np
import pytest

from elecsus import elecsus_dash_gui as gui
from elecsus.elecsus_methods import calculate, calculate_cascade


def _spec(output_id, inputs):
    """The callback with this first output (a duplicate output's id carries a
    hash after the property) and these inputs."""
    ids = {i for i, _ in inputs}
    for key, spec in gui.app.callback_map.items():
        first = key.split('...')[0].lstrip('.')
        if first.split('.')[0] == output_id and ids == {i['id'] for i in spec['inputs']}:
            return key, spec
    raise KeyError(output_id)


def post_callback(output_id, inputs, values):
    """Run the callback whose first output is `output_id` through Dash's HTTP
    endpoint. `inputs` is [(id, property)]; `values` maps each State's id to
    its value. Returns {(id, property): value} for the callback's outputs."""
    key, spec = _spec(output_id, inputs)
    missing = [s['id'] for s in spec['state'] if s['id'] not in values]
    assert not missing, f'give these State boxes a value: {missing}'
    body = {
        'output': key,
        'outputs': [{'id': o['id'], 'property': o['property']} for o in spec['output']]
                   if isinstance(spec['output'][0], dict) else
                   [{'id': o.component_id, 'property': o.component_property}
                    for o in spec['output']],
        'inputs': [{'id': i, 'property': p, 'value': 1} for i, p in inputs],
        'state': [{'id': s['id'], 'property': s['property'], 'value': values[s['id']]}
                  for s in spec['state']],
        'changedPropIds': [f'{i}.{p}' for i, p in inputs],
    }
    r = gui.server.test_client().post('/_dash-update-component', json=body)
    assert r.status_code == 200, r.data
    resp = r.get_json()['response']
    return {(i, p): v for i, props in resp.items() for p, v in props.items()}


# One value in every single-cell box, none of them a default
SINGLE = {
    's-elem': 'Rb', 's-dline': 'D2', 's-T': 61., 's-lcell': 0.033, 's-B': 123.,
    's-gbuf': 7., 's-shift': 11., 's-btheta': 30., 's-bphi': 20., 's-theta0': 25.,
    's-pol': 'Linear (x)', 's-bufgas': 'None', 's-bufP': 3., 's-bufTf': 301.,
    's-rb85frac': 60., 's-K40frac': 0.5, 's-K41frac': 9.,
    's-det-start': -4., 's-det-stop': 3., 's-det-N': 50, 's-outputs': ['S0'],
    's-nonortho': ['on'], 's-windows': [], 's-nwin': 1.6, 's-doppT': None,
    's-xunit': 'GHz',
}


def _core_single(over=None, outputs=('S0',), E_in=None, **p):
    base = dict(Elem='Rb', Dline='D2', T=61., lcell=0.033, Bfield=123., GammaBuf=7.,
                shift=11., Btheta=np.radians(30.), Bphi=np.radians(20.), rb85frac=60.)
    X = np.linspace(-4e3, 3e3, 50)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = calculate(X, E_in if E_in is not None else gui._pol_efield('Linear (x)', 25.),
                        dict(base, **p), outputs=list(outputs))
    return [np.real(r) for r in res]


def test_compute_button_uses_the_values_in_the_boxes():
    got = post_callback('s-plot', [('s-compute', 'n_clicks')], SINGLE)
    [want] = _core_single()
    np.testing.assert_allclose(got[('s-store', 'data')]['results']['S0'], want, rtol=1e-9)
    assert got[('error-store', 'data')] == ''


def test_compute_button_sends_the_gas_the_isotopes_and_the_polariser():
    gas = dict(SINGLE, **{'s-bufgas': 'N2', 's-bufP': 8., 's-bufTf': 301., 's-pol': 'Left CP'})
    got = post_callback('s-plot', [('s-compute', 'n_clicks')], gas)
    [want] = _core_single(E_in=[1 / np.sqrt(2), 1j / np.sqrt(2), 0.], BufferGas='N2',
                          BufferPressure=8., BufferTemp_fill=301.)
    np.testing.assert_allclose(got[('s-store', 'data')]['results']['S0'], want, rtol=1e-9)

    k = dict(SINGLE, **{'s-elem': 'K', 's-dline': 'D1'})
    got = post_callback('s-plot', [('s-compute', 'n_clicks')], k)
    [want] = _core_single(Elem='K', Dline='D1', K40frac=0.5, K41frac=9.)
    np.testing.assert_allclose(got[('s-store', 'data')]['results']['S0'], want, rtol=1e-9)


# Cell 1 and cell 2 differ in every box
CASCADE = {
    'c1-elem': 'Rb', 'c1-dline': 'D2', 'c1-T': 70., 'c1-lcell': 0.02, 'c1-B': 90.,
    'c1-gbuf': 4., 'c1-shift': 3., 'c1-btheta': 10., 'c1-bphi': 15., 'c1-theta0': 0.,
    'c1-pol': 'Linear (x)', 'c1-bufgas': 'None', 'c1-bufP': 0., 'c1-bufTf': 293.15,
    'c1-rb85frac': 72.17, 'c1-K40frac': 0.01, 'c1-K41frac': 6.73,
    'c2-T': 95., 'c2-lcell': 0.045, 'c2-B': 210., 'c2-gbuf': 9., 'c2-shift': -5.,
    'c2-btheta': 0., 'c2-bphi': 0., 'c2-bufgas': 'None', 'c2-bufP': 0.,
    'c2-bufTf': 293.15,
    'cas-use-pol': ['on'], 'cas-pol-angle': 70., 'cas-reset-pol': 'none',
    'cas-det-start': -3., 'cas-det-stop': 4., 'cas-det-N': 40, 'cas-outputs': ['S0', 'S3'],
    'c1-nonortho': ['on'], 'c1-windows': [], 'c1-nwin': 1.51, 'c1-doppT': None,
    'c2-nonortho': ['on'], 'c2-windows': [], 'c2-nwin': 1.51, 'c2-doppT': None,
    'cas-xunit': 'GHz',
}


def test_cascade_button_uses_the_values_in_the_boxes():
    got = post_callback('cas-plot', [('cas-compute', 'n_clicks')], CASCADE)
    p1 = dict(Elem='Rb', Dline='D2', T=70., lcell=0.02, Bfield=90., GammaBuf=4., shift=3.,
              Btheta=np.radians(10.), Bphi=np.radians(15.), rb85frac=72.17)
    p2 = dict(Elem='Rb', Dline='D2', T=95., lcell=0.045, Bfield=210., GammaBuf=9., shift=-5.,
              Btheta=0., Bphi=0., rb85frac=72.17)
    X = np.linspace(-3e3, 4e3, 40)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r1, r2 = calculate_cascade(X, [1, 0, 0], p1, p2, outputs=['S0', 'S3'],
                                   polariser_angle=np.radians(70.))
    store = got[('cas-store', 'data')]
    np.testing.assert_allclose(store['total']['S0'], np.real(np.squeeze(r2[0])), rtol=1e-9)
    np.testing.assert_allclose(store['total']['S3'], np.real(np.squeeze(r2[1])), rtol=1e-9)
    np.testing.assert_allclose(store['cell1']['S0'], np.real(np.squeeze(r1[0])), rtol=1e-9)
    assert got[('error-store', 'data')] == ''



def _trace_x(trace):
    """A trace's x from the JSON the page receives: plotly sends numbers as
    base64 typed arrays."""
    x = trace['x']
    if isinstance(x, dict):
        return np.frombuffer(base64.b64decode(x['bdata']), dtype=x['dtype'])
    return np.asarray(x, dtype=float)


def test_the_plot_against_selector_reaches_both_compute_buttons():
    """The wavelength axis: both Compute buttons read "Plot against", and the
    stores keep GHz and the line for redrawing."""
    got = post_callback('s-plot', [('s-compute', 'n_clicks')], dict(SINGLE, **{'s-xunit': 'mA'}))
    store = got[('s-store', 'data')]
    want = gui._x_on(store['X'], 'mA', 'Rb', 'D2')
    np.testing.assert_allclose(_trace_x(got[('s-plot', 'figure')]['data'][0]), want)
    assert (store['elem'], store['dline'], store['X'][0]) == ('Rb', 'D2', -4.)

    got = post_callback('cas-plot', [('cas-compute', 'n_clicks')],
                        dict(CASCADE, **{'cas-xunit': 'mA'}))
    store = got[('cas-store', 'data')]
    want = gui._x_on(store['X'], 'mA', 'Rb', 'D2')
    for trace in got[('cas-plot', 'figure')]['data']:
        np.testing.assert_allclose(_trace_x(trace), want)


# The fit panel: its own boxes, and the single-cell boxes it reads
FIT_T_TRUE = 52.
FIT = {**SINGLE, 's-fit-datatype': 'S0', 's-fit-algo': 'ML', 's-fit-vary': ['T'],
       's-fit-min': None, 's-fit-max': None, 's-fit-data-store': None}


def _fit_data(**p):
    base = dict(Elem='Rb', Dline='D2', lcell=0.033, Bfield=123., GammaBuf=7.,
                shift=11., Btheta=np.radians(30.), Bphi=np.radians(20.), rb85frac=60.)
    X = np.linspace(-4e3, 3e3, 160)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        [S0] = calculate(X, gui._pol_efield('Linear (x)', 25.), dict(base, **p), outputs=['S0'])
    return {'X': list(X / 1e3), 'Y': list(np.real(S0))}


def test_fit_button_uses_the_values_in_the_boxes():
    """The data are the model at T = 52 with every other box as typed; a fit of
    T recovers 52 only if each State reaches the model in the right place."""
    values = dict(FIT, **{'s-fit-data-store': _fit_data(T=FIT_T_TRUE), 's-T': 45.})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        got = post_callback('s-fit-plot', [('s-fit-run', 'n_clicks')], values)
    assert got[('error-store', 'data')] == ''
    assert got[('s-T', 'value')] == pytest.approx(FIT_T_TRUE, abs=0.01)
    # the other boxes are written back as typed
    assert got[('s-B', 'value')] == pytest.approx(123.)
    assert got[('s-gbuf', 'value')] == pytest.approx(7.)
    assert got[('s-shift', 'value')] == pytest.approx(11.)
    assert got[('s-lcell', 'value')] == pytest.approx(0.033)
    assert got[('s-btheta', 'value')] == pytest.approx(30.)
    assert got[('s-bphi', 'value')] == pytest.approx(20.)


def test_fit_button_honours_the_fit_range_boxes():
    """Fit min and Fit max reach the right variables (not swapped)."""
    values = dict(FIT, **{'s-fit-data-store': _fit_data(T=FIT_T_TRUE),
                          's-fit-min': -1.5, 's-fit-max': 1.0})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        got = post_callback('s-fit-plot', [('s-fit-run', 'n_clicks')], values)
    store = got[('s-fit-result-store', 'data')]
    assert min(store['X']) >= -1.5 and max(store['X']) <= 1.0
    assert 0 < len(store['X']) < 160


def test_fit_button_sends_the_gas_and_the_isotopes():
    gas = dict(BufferGas='N2', BufferPressure=8., BufferTemp_fill=301.)
    values = dict(FIT, **{'s-fit-data-store': _fit_data(T=FIT_T_TRUE, **gas), 's-T': 47.,
                          's-bufgas': 'N2', 's-bufP': 8., 's-bufTf': 301.})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        got = post_callback('s-fit-plot', [('s-fit-run', 'n_clicks')], values)
    assert got[('error-store', 'data')] == ''
    assert got[('s-T', 'value')] == pytest.approx(FIT_T_TRUE, abs=0.01)

    k = dict(Elem='K', Dline='D1', K40frac=0.5, K41frac=9.)
    values = dict(FIT, **{'s-fit-data-store': _fit_data(T=FIT_T_TRUE, **k), 's-T': 47.,
                          's-elem': 'K', 's-dline': 'D1'})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        got = post_callback('s-fit-plot', [('s-fit-run', 'n_clicks')], values)
    assert got[('error-store', 'data')] == ''
    assert got[('s-T', 'value')] == pytest.approx(FIT_T_TRUE, abs=0.01)


def test_cascade_button_sends_each_cells_gas_and_the_isotopes():
    values = dict(CASCADE, **{
        'c1-elem': 'K', 'c1-dline': 'D1', 'c1-T': 150., 'c2-T': 160.,
        'c1-K40frac': 0.5, 'c1-K41frac': 9.,
        'c1-bufgas': 'N2', 'c1-bufP': 6., 'c1-bufTf': 300.,
        'c2-bufgas': 'He', 'c2-bufP': 12., 'c2-bufTf': 310.,
        'cas-use-pol': [], 'cas-outputs': ['S0']})
    got = post_callback('cas-plot', [('cas-compute', 'n_clicks')], values)
    iso = dict(Elem='K', Dline='D1', K40frac=0.5, K41frac=9.)
    p1 = dict(iso, T=150., lcell=0.02, Bfield=90., GammaBuf=4., shift=3.,
              Btheta=np.radians(10.), Bphi=np.radians(15.),
              BufferGas='N2', BufferPressure=6., BufferTemp_fill=300.)
    p2 = dict(iso, T=160., lcell=0.045, Bfield=210., GammaBuf=9., shift=-5.,
              Btheta=0., Bphi=0., BufferGas='He', BufferPressure=12., BufferTemp_fill=310.)
    X = np.linspace(-3e3, 4e3, 40)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r1, r2 = calculate_cascade(X, [1, 0, 0], p1, p2, outputs=['S0'])
    store = got[('cas-store', 'data')]
    assert got[('error-store', 'data')] == ''
    S0_1 = np.real(np.squeeze(r1[0]))
    np.testing.assert_allclose(store['cell1']['S0'], S0_1, rtol=1e-9)
    np.testing.assert_allclose(store['total']['S0'], S0_1 * np.real(np.squeeze(r2[0])), rtol=1e-9)
