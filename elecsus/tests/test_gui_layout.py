"""GUI presentation: per-unit plot panels and offline styling.

Outputs of different units share no axis: an absorption coefficient in m^-1
or an angle in radians would flatten a 0-1 transmission. And the page must
render the same with no network, so the theme and font come from
elecsus/assets rather than a CDN.
"""
import os
import re

import numpy as np

from elecsus import elecsus_dash_gui as gui


def _compute(outputs):
    return gui.compute_single(
        1, 'Rb', 'D2', 60., 0.075, 300., 0., 0., 0., 0., 0., 'Linear (x)',
        'None', 0., 293.15, 72.17, 0.01, 6.73, -6., 6., 200, list(outputs))[:3]


def _axis_titles(fig):
    return {name: fig.layout[name].title.text
            for name in fig.layout if name.startswith('yaxis')}


def test_outputs_of_different_units_get_separate_panels():
    fig, _, err = _compute(['S0', 'S1', 'phi', 'alphaPlus', 'nMinus'])
    assert err == ''
    titles = _axis_titles(fig)
    assert sorted(titles.values()) == sorted(
        ['Signal', 'Rotation (rad)', 'Absorption coeff. (m⁻¹)', 'Refractive index'])
    axis_of = {t.name: t.yaxis for t in fig.data}
    # dimensionless outputs share the first panel ...
    assert axis_of['S₀ (Transmission)'] == axis_of['S₁'] == 'y'
    # ... and each other unit has its own, in OUTPUT_GROUPS order
    assert titles['yaxis' + axis_of[gui.OUTPUT_LABELS['phi']][1:]] == 'Rotation (rad)'
    assert titles['yaxis' + axis_of[gui.OUTPUT_LABELS['alphaPlus']][1:]] == \
        'Absorption coeff. (m⁻¹)'
    # one shared detuning axis, labelled once at the bottom
    assert sum(1 for name in fig.layout if name.startswith('xaxis')
               and fig.layout[name].title.text == 'Detuning (GHz)') == 1


def test_a_single_unit_keeps_a_single_panel():
    fig, _, err = _compute(['S0', 'S1', 'DoLP'])
    assert err == ''
    assert _axis_titles(fig) == {'yaxis': 'Signal'}


def test_every_output_has_a_label_and_a_panel():
    for key in gui.OUTPUTS:
        assert key in gui.OUTPUT_LABELS
        assert 0 <= gui._group_of(key) < len(gui.OUTPUT_GROUPS)


def _components(node):
    yield node
    children = getattr(node, 'children', None)
    for child in (children if isinstance(children, (list, tuple)) else [children]):
        if child is not None and hasattr(child, 'to_plotly_json'):
            yield from _components(child)


def _number_inputs():
    return [c for c in _components(gui.app.layout)
            if type(c).__name__ == 'Input' and getattr(c, 'type', None) == 'number']


def test_number_inputs_accept_any_value():
    """A number input with a step, min or max marks values outside them
    invalid, and Dash then sends None, which the callbacks used to replace with
    a default: B = 239 G (step 10) silently computed at B = 0, and Btheta = 200
    (max 180) at 0. Every numeric input takes step='any' and no
    min or max; the callbacks and the core's validation refuse bad values."""
    numbers = _number_inputs()
    assert len(numbers) > 30
    stepped = [c.id for c in numbers if getattr(c, 'step', None) != 'any']
    assert not stepped, f'inputs with a fixed step: {stepped}'
    limited = [c.id for c in numbers
               if getattr(c, 'min', None) is not None or getattr(c, 'max', None) is not None]
    assert not limited, f'inputs with a min or max: {limited}'


ASSETS = os.path.join(os.path.dirname(gui.__file__), 'assets')


def test_styles_come_from_local_assets_not_a_cdn():
    assert not gui.app.config.external_stylesheets
    css = [f for f in os.listdir(ASSETS) if f.endswith('.css')]
    assert 'bootstrap-flatly.min.css' in css and 'inter.css' in css
    for name in css:
        text = open(os.path.join(ASSETS, name), encoding='utf-8').read()
        assert '@import' not in text, f'{name} imports another stylesheet'
        for ref in re.findall(r'url\(\s*[\'"]?([^)\'"]+)', text):
            if ref.startswith('data:'):
                continue
            assert '://' not in ref, f'{name} fetches {ref}'
            assert os.path.exists(os.path.join(ASSETS, ref)), f'{name}: missing {ref}'


def test_bundled_fonts_and_theme_carry_their_licences():
    files = set(os.listdir(ASSETS))
    assert {'LICENSE-Bootswatch.txt', 'LICENSE-Inter-OFL.txt'} <= files
    assert {f'inter-latin-{w}-normal.woff2' for w in (400, 500, 700)} <= files
