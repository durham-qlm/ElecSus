"""The GUI's Energy Levels tab: the levels of one isotope in a field, the lines
between them and the spectrum, from elecsus.energy_levels (tested in
test_energy_levels.py)."""
import warnings

import numpy as np
import pytest

from elecsus import elecsus_dash_gui as gui
from elecsus import energy_levels as EL
from elecsus.tests.test_gui_wiring import post_callback

# the boxes as the page starts, but for the field
BOXES = dict(iso='Rb87', dline='D1', B=6000., geometry='faraday', min_strength=0.01,
             labels='mJmI', T=75., lcell=0.002, det_start=None, det_stop=None, natural=[])


def _compute(**over):
    a = dict(BOXES, **over)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return gui.compute_levels(1, *a.values())


def test_the_tab_draws_the_lines_of_the_core():
    fig, table, store, err = _compute()
    assert err == ''
    want = EL.transitions('Rb87', 'D1', 6000., kinds=('sigma+', 'sigma-'), min_strength=0.01)
    assert len(store['rows']) == len(want) == 14
    np.testing.assert_allclose([r['detuning'] for r in store['rows']],
                               sorted(line['detuning'] for line in want))
    # an arrow (a shaft and a head) for each line, and its details on hover
    shafts = [s for s in fig.layout.shapes if s.yref == 'paper']
    heads = [a for a in fig.layout.annotations if a.showarrow]
    [hover] = [t for t in fig.data if t.name == 'lines']
    assert len(shafts) == len(heads) == len(hover.x) == 14
    assert {t.split()[0] for t in hover.text} == {'σ+', 'σ−'}
    assert '14 lines' in table[0].children
    assert 'B = 6000 G (0.60 T)' in fig.layout.title.text


def test_voigt_adds_the_pi_lines_and_the_labels_follow_the_choice():
    fig, table, store, err = _compute(iso='Cs', dline='D2', B=1000., geometry='voigt',
                                      labels='FmF')
    assert err == ''
    assert {r['kind'] for r in store['rows']} == {'σ+', 'σ−', 'π'}
    texts = [a.text for a in fig.layout.annotations]
    assert '|4, 4⟩' in texts and '|F, mF⟩ at B = 0' in texts


def test_the_level_labels_stay_inside_their_panels():
    fig, _, _, _ = _compute(iso='Cs', dline='D2', B=1000., geometry='voigt')
    for a in fig.layout.annotations:
        if a.xref == 'paper' and a.yref.startswith('y') and 'domain' not in a.yref:
            r0, r1 = fig.layout[a.yref.replace('y', 'yaxis')].range
            assert r0 <= a.y <= r1, a.text


def test_the_detuning_range_and_natural_abundance():
    fig, _, _, err = _compute(det_start=-20., det_stop=20.)
    assert err == ''
    spectrum = fig.data[0]
    assert spectrum.x[0] == -20. and spectrum.x[-1] == 20.
    natural, _, _, _ = _compute(det_start=-20., det_stop=20., natural=['on'])
    assert 'natural abundance' in natural.layout.title.text
    assert not np.allclose(natural.data[0].y, spectrum.y)     # 85Rb's lines appear


@pytest.mark.parametrize('over, message', [
    (dict(B=None), 'B is empty'), (dict(B=-5.), '0 G or more'),
    (dict(min_strength=2.), '0 to 1'), (dict(det_start=-5.), 'both ends'),
    (dict(det_start=5., det_stop=-5.), 'below Stop'), (dict(min_strength=None), 'Min. strength'),
    (dict(T=None), 'T is empty'),
])
def test_bad_boxes_are_named_and_nothing_is_kept(over, message):
    fig, table, store, err = _compute(**over)
    assert message in err and store is None and table == []


def test_the_export_writes_one_row_per_line():
    _, _, store, _ = _compute()
    out = gui.export_levels(1, store)
    rows = out['content'].splitlines()
    assert rows[0].startswith('Kind,Detuning_GHz,Relative_strength,Strength,Ground_mJmI')
    assert len(rows) == 15 and out['filename'] == 'elecsus_levels_Rb87_D1_6000G.csv'
    assert rows[1].startswith('σ−,-14.037028,1.000000,')
    assert gui.export_levels(1, None) is gui.no_update


def test_the_tab_button_shows_the_tab():
    import dash
    from contextvars import copy_context
    from dash._callback_context import context_value
    from dash._utils import AttributeDict

    def run():
        context_value.set(AttributeDict(triggered_inputs=[{'prop_id': 'tab-levels-btn.n_clicks',
                                                           'value': 1}]))
        return gui.switch_tab(0, 0, 1)
    styles = copy_context().run(run)
    assert [s['display'] for s in styles[:3]] == ['none', 'none', 'block']
    assert styles[5]['borderBottom'] != 'none' and styles[3]['borderBottom'] == 'none'


def test_the_compute_button_uses_the_values_in_the_boxes():
    """Through Dash's endpoint: each box reaches its argument."""
    values = {'lv-iso': 'K39', 'lv-dline': 'D2', 'lv-B': 1234., 'lv-geometry': 'voigt',
              'lv-min': 0.05, 'lv-labels': 'FmF', 'lv-T': 61., 'lv-lcell': 0.004,
              'lv-det-start': -3., 'lv-det-stop': 4., 'lv-natural': []}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = post_callback('lv-plot', [('lv-compute', 'n_clicks')], values)
    assert out[('error-store', 'data')] == ''
    store = out[('lv-store', 'data')]
    assert (store['iso'], store['dline'], store['B']) == ('K39', 'D2', 1234.)
    want = EL.transitions('K39', 'D2', 1234., min_strength=0.05)
    assert len(store['rows']) == len(want)
    title = out[('lv-plot', 'figure')]['layout']['title']['text']
    assert 'T = 61°C, l = 0.004 m' in title
    x = out[('lv-plot', 'figure')]['data'][0]['x']
    if isinstance(x, dict):              # plotly 6+ sends arrays as binary
        x = np.frombuffer(__import__('base64').b64decode(x['bdata']), dtype=x['dtype'])
    assert x[0] == -3. and x[-1] == 4.


@pytest.mark.parametrize('over', [dict(), dict(iso='Cs', dline='D2', B=1000., geometry='voigt')])
def test_each_arrow_runs_from_its_ground_level_to_its_excited_level(over):
    """The shaft from the ground level to the excited level (in paper
    coordinates, from the axis ranges and domains), and the head's tip on the
    excited level itself (a triangle marker, centred on the level, put its
    tip above it)."""
    fig, _, store, err = _compute(**over)
    assert err == ''
    a = dict(BOXES, **over)
    levels = EL.energy_levels(a['iso'], a['dline'], a['B'])
    lines = EL.transitions(a['iso'], a['dline'], a['B'], kinds=EL.default_kinds(a['geometry']),
                           min_strength=a['min_strength'])

    def paper(yaxis, y):
        lo, hi = fig.layout[yaxis].range
        d0, d1 = fig.layout[yaxis].domain
        return d0 + (y - lo) / (hi - lo) * (d1 - d0)

    shafts = sorted((s for s in fig.layout.shapes if s.yref == 'paper'), key=lambda s: s.x0)
    heads = sorted((h for h in fig.layout.annotations if h.showarrow), key=lambda h: h.x)
    lines = sorted(lines, key=lambda l: l['detuning'])
    assert len(shafts) == len(heads) == len(lines)
    exc_axis = fig.get_subplot(2, 2).yaxis.plotly_name
    gnd_axis = fig.get_subplot(3, 2).yaxis.plotly_name
    for shaft, head, line in zip(shafts, heads, lines):
        e = levels['excited'][line['excited']]['energy']
        g = levels['ground'][line['ground']]['energy']
        assert shaft.x0 == shaft.x1 == pytest.approx(line['detuning'])
        assert shaft.y0 == pytest.approx(paper(gnd_axis, g))
        assert shaft.y1 == pytest.approx(paper(exc_axis, e))
        assert head.yref == exc_axis.replace('axis', '') and head.y == pytest.approx(e)
        assert head.x == pytest.approx(line['detuning']) and head.ay > 0 and head.standoff == 0


def test_a_detuning_range_is_honoured_and_lines_outside_it_are_not_drawn():
    """Start and Stop set every panel's detuning axis (the column's axes are
    matched, and a range on one of them was overridden by the others'
    autorange), and a line outside the range has no arrow or mark: drawn,
    its shaft crossed the gaps between the panels. The table still lists it."""
    fig, table, store, err = _compute(det_start=-5., det_stop=5.)
    assert err == ''
    for r in (1, 2, 3):
        assert tuple(fig.layout[fig.get_subplot(r, 2).xaxis.plotly_name].range) == (-5., 5.)
    inside = [r for r in store['rows'] if -5 <= r['detuning'] <= 5]
    shafts = [s for s in fig.layout.shapes if s.yref == 'paper']
    heads = [a for a in fig.layout.annotations if a.showarrow]
    marks = [s for s in fig.layout.shapes if s.type == 'line' and s.yref != 'paper'
             and s.xref == fig.get_subplot(1, 2).xaxis.plotly_name.replace('axis', '')]
    assert len(shafts) == len(heads) == len(marks) == len(inside) == 1
    assert len(store['rows']) == 14
    assert '13 outside the detuning range' in table[0].children


def test_zero_field_draws_a_field_axis_that_can_be_drawn():
    fig, _, store, err = _compute(B=0.)
    assert err == '' and store['rows']
    lo, hi = fig.layout[fig.get_subplot(3, 1).xaxis.plotly_name].range
    assert lo == 0 and hi > 0


def test_the_stick_spectrum_shows_every_line_in_the_range():
    fig, _, _, err = _compute(iso='Na', B=2440., T=176.)
    assert err == ''
    every = EL.transitions('Na', 'D1', 2440., kinds=('sigma+', 'sigma-'))
    stems = [t for t in fig.data if t.mode == 'lines' and t.name in ('σ+', 'σ−')]
    tips = [t for t in fig.data if t.mode == 'markers' and t.name.endswith('strength')]
    xs = sorted(x for t in tips for x in t.x)
    np.testing.assert_allclose(xs, sorted(l['detuning'] for l in every))
    assert sum(len(t.x) for t in stems) == 3 * len(every)       # x, x, None per stick
    assert {t.name for t in stems} == {'σ+', 'σ−'} and all(t.showlegend for t in stems)
    yaxis = fig.layout[fig.get_subplot(4, 2).yaxis.plotly_name]
    assert yaxis.type == 'log' and yaxis.range[0] == pytest.approx(-4.)
    # and only those inside a range that is set
    fig, _, _, _ = _compute(det_start=-5., det_stop=5.)
    tips = [x for t in fig.data if t.mode == 'markers' and t.name.endswith('strength')
            for x in t.x]
    assert tips and all(-5 <= x <= 5 for x in tips)
