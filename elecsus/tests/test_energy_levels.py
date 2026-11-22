"""Energy levels, transitions and level diagrams (elecsus.energy_levels).

The lines must be the ones the spectra are made of: the same detunings and
strengths as spectra.FreqStren, for every isotope and line. The labels are
checked against the figure of Rb87 D1 at 0.6 T in a former student's thesis
(made by hand from ES3's eigenvectors), the zero-field levels against the
measured hyperfine splittings, and the strengths against the sum rule.
"""
import contextlib
import io

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pytest

from elecsus import energy_levels as EL
from elecsus.libs import EigenSystem as ES
from elecsus.libs.spectra import FreqStren

HAND = {'sigma+': 'Right', 'sigma-': 'Left', 'pi': 'Z'}


@pytest.mark.parametrize('isotope', sorted(EL.ISOTOPES))
@pytest.mark.parametrize('Dline', ['D1', 'D2'])
@pytest.mark.parametrize('B', [0., 300., 6000.])
def test_the_lines_are_those_of_the_spectra(isotope, Dline, B):
    lines = EL.transitions(isotope, Dline, B)
    system = ES.Hamiltonian(isotope, Dline, 1.0, B)
    for kind, hand in HAND.items():
        f, s, n = FreqStren(system.groundManifold, system.excitedManifold, system.ds,
                            system.dp, Dline, hand, BoltzmannFactor=False)
        got = [(line['detuning'], line['strength']) for line in lines if line['kind'] == kind]
        # at zero field the levels are degenerate and the split of a line's
        # strength among a degenerate set depends on the basis: compare the
        # strength at each frequency. FreqStren leaves out the pieces below
        # cleb^2 = 5e-4 in its own basis, so the totals differ by those.
        np.testing.assert_allclose(_by_frequency(got), _by_frequency(zip(f[:n] / 1e3, s[:n])),
                                   atol=1e-9 if B else 1e-4)


def _by_frequency(pairs):
    total = {}
    for f, s in pairs:
        key = round(f, 6)
        total[key] = total.get(key, 0.) + s
    return sorted(total.items())


def test_the_labels_are_those_of_the_thesis_figure():
    """Rb87 D1 at 0.6 T: the ground states' |mJ, mI> components, top to
    bottom, as the figure gives them (the figure drops the signs)."""
    ground = [level['label'] for level in EL.energy_levels('Rb87', 'D1', 6000.)['ground'][::-1]]
    assert ground == ['|1/2, 3/2⟩',
                      '0.99 |1/2, 1/2⟩ + 0.14 |-1/2, 3/2⟩',
                      '0.98 |1/2, -1/2⟩ + 0.19 |-1/2, 1/2⟩',
                      '0.98 |1/2, -3/2⟩ + 0.21 |-1/2, -1/2⟩',
                      '|-1/2, -3/2⟩',
                      '0.98 |-1/2, -1/2⟩ − 0.21 |1/2, -3/2⟩',
                      '0.98 |-1/2, 1/2⟩ − 0.19 |1/2, -1/2⟩',
                      '0.99 |-1/2, 3/2⟩ − 0.14 |1/2, 1/2⟩']
    excited = [level['label'] for level in EL.energy_levels('Rb87', 'D1', 6000.)['excited'][::-1]]
    assert excited == ['|1/2, 3/2⟩', '|1/2, 1/2⟩', '|1/2, -1/2⟩', '|1/2, -3/2⟩',
                       '|-1/2, -3/2⟩', '|-1/2, -1/2⟩', '|-1/2, 1/2⟩', '|-1/2, 3/2⟩']


def test_the_thesis_lines_sit_where_the_figure_puts_them():
    """The eight strong lines of the figure: sigma- from -14 to -9 GHz, sigma+
    from +8.4 to +14.6 GHz."""
    lines = EL.transitions('Rb87', 'D1', 6000., kinds=('sigma+', 'sigma-'), min_strength=0.5)
    minus = sorted(line['detuning'] for line in lines if line['kind'] == 'sigma-')
    plus = sorted(line['detuning'] for line in lines if line['kind'] == 'sigma+')
    np.testing.assert_allclose(minus, [-14.04, -12.56, -10.88, -8.91], atol=0.01)
    np.testing.assert_allclose(plus, [8.41, 10.96, 12.90, 14.55], atol=0.01)


@pytest.mark.parametrize('isotope, F, splitting', [
    ('Rb87', (1, 2), 6.834682611), ('Rb85', (2, 3), 3.035732439),
    ('Cs', (3, 4), 9.192631770), ('K39', (1, 2), 0.461719720), ('Na', (1, 2), 1.771626129),
    ('K40', (9 / 2, 7 / 2), 1.285790708),     # inverted: the larger F lies lower
])
def test_zero_field_levels_are_the_hyperfine_levels(isotope, F, splitting):
    ground = EL.energy_levels(isotope, 'D1', 0.)['ground']
    lower, upper = F
    for f in F:
        assert sum(np.isclose(level['F'], f) for level in ground) == 2 * f + 1
    E = {f: np.mean([level['energy'] for level in ground if np.isclose(level['F'], f)])
         for f in F}
    # ElecSus's hyperfine constants give the measured splittings to 3 kHz
    assert E[upper] - E[lower] == pytest.approx(splitting, abs=3e-6)
    # the energies are from the centre of gravity of the manifold
    assert sum(level['energy'] for level in ground) == pytest.approx(0., abs=1e-9)


@pytest.mark.parametrize('isotope', ['Rb85', 'Cs', 'K40'])
@pytest.mark.parametrize('B', [0., 450., 6000.])
def test_the_strengths_obey_the_sum_rule(isotope, B):
    """Over both D lines and all three kinds, each ground state's lines add up
    to 1/3 of the full strength per polarisation: in all, 1."""
    total = sum(line['strength'] for Dline in ('D1', 'D2')
                for line in EL.transitions(isotope, Dline, B))
    assert total == pytest.approx(1., abs=5e-3)


def test_breit_rabi_ends_on_the_levels_and_starts_at_zero_field():
    curves = EL.breit_rabi('Rb87', 'D2', 2000., n=81)
    levels = EL.energy_levels('Rb87', 'D2', 2000.)
    zero = EL.energy_levels('Rb87', 'D2', 0.)
    for which in ('ground', 'excited'):
        np.testing.assert_allclose(curves[which][-1],
                                   [level['energy'] for level in levels[which]], atol=1e-9)
        np.testing.assert_allclose(sorted(curves[which][0]),
                                   sorted(level['energy'] for level in zero[which]), atol=1e-6)
        # each column is one state followed in field: no jumps between steps
        assert np.max(np.abs(np.diff(curves[which], axis=0))) < 0.1


def test_mF_is_kept_and_F_is_the_zero_field_one():
    for level in EL.energy_levels('Cs', 'D2', 3000.)['excited']:
        assert level['F'] in (2, 3, 4, 5) and abs(level['mF']) <= level['F']


def test_the_jmj_basis_is_orthonormal():
    for L, I in [(0, 1.5), (1, 1.5), (1, 3.5), (1, 4.)]:
        C, _ = EL._jmj_basis(L, I)
        np.testing.assert_allclose(C.T @ C, np.identity(len(C)), atol=1e-12)


def test_the_spectrum_dips_at_the_strong_lines():
    lines = EL.transitions('Rb87', 'D1', 6000., kinds=EL.default_kinds('faraday'),
                           min_strength=0.5)
    X = np.linspace(-17, 17, 3401)
    [S0] = EL.spectrum('Rb87', 'D1', 6000., X, T=75., lcell=2e-3)
    for line in lines:
        near = np.abs(X - line['detuning']) < 0.6
        assert abs(X[near][np.argmin(S0[near])] - line['detuning']) < 0.05


def test_a_weak_line_threshold_and_the_kinds_filter():
    sigma = ('sigma+', 'sigma-')
    every = EL.transitions('Rb87', 'D1', 6000., kinds=sigma)
    strong = EL.transitions('Rb87', 'D1', 6000., kinds=sigma, min_strength=0.5)
    assert len(strong) == 8 < len(every)
    assert {line['kind'] for line in EL.transitions('Rb87', 'D2', 100., kinds=['pi'])} == {'pi'}
    assert max(line['relative'] for line in every) == 1.


@pytest.mark.parametrize('bad, message', [
    (dict(isotope='Rb86'), 'isotope'), (dict(Dline='D3'), 'Dline'),
    (dict(kinds=['sigma']), 'kinds'),
])
def test_bad_arguments_are_refused(bad, message):
    args = dict(isotope='Rb87', Dline='D1', Bfield=100.)
    args.update(bad)
    with pytest.raises(ValueError, match=message):
        EL.transitions(**args)


@pytest.mark.parametrize('isotope, Dline, geometry, labels', [
    ('Rb87', 'D1', 'faraday', 'mJmI'), ('Cs', 'D2', 'voigt', 'FmF'),
    ('K40', 'D1', 'faraday', 'FmF'), ('Na', 'D2', 'voigt', 'mJmI'),
])
def test_the_figure_is_drawn(isotope, Dline, geometry, labels):
    import matplotlib.pyplot as plt
    fig = EL.level_diagram(isotope, Dline, 1500., geometry=geometry, labels=labels,
                           n_points=400)
    assert len(fig.axes) == 6            # with the stick spectrum
    plt.close(fig)
    fig = EL.level_diagram(isotope, Dline, 1500., geometry=geometry, sticks=False, n_points=400)
    assert len(fig.axes) == 5
    plt.close(fig)


def test_the_script_prints_the_lines_and_saves_the_figure(tmp_path):
    out = tmp_path / 'fig.png'
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        EL.main(['Rb87', 'D1', '6000', '--out', str(out)])
    text = buffer.getvalue()
    assert 'Rb87 D1, B = 6000 G: 14 lines' in text and out.stat().st_size > 10000
    assert 'σ+    +8.414     1.000  |-1/2, -3/2⟩  ->  |1/2, -3/2⟩' in text


@pytest.mark.parametrize('isotope, Dline, B, geometry', [
    ('Rb87', 'D1', 6000., 'faraday'), ('Cs', 'D2', 1000., 'voigt'), ('K39', 'D2', 300., 'voigt'),
])
def test_each_arrow_runs_from_its_ground_level_to_its_excited_level(isotope, Dline, B, geometry):
    """The arrow's tail on the ground level and its tip on the excited level,
    at the line's detuning. (A stroked head overshot by half its line width,
    up to a level spacing at 0.6 T.)"""
    import matplotlib.pyplot as plt
    fig = EL.level_diagram(isotope, Dline, B, geometry=geometry, n_points=200)
    fig.canvas.draw()
    ax_e, ax_g = fig.axes[2], fig.axes[4]
    levels = EL.energy_levels(isotope, Dline, B)
    arrows = [a for a in fig.artists if hasattr(a, 'line')]
    assert arrows
    for arrow in arrows:
        line = arrow.line
        paths, _ = arrow._get_path_in_displaycoord()
        v = np.concatenate([p.vertices for p in (paths if isinstance(paths, list) else [paths])])
        tip = ax_e.transData.inverted().transform(v[np.argmax(v[:, 1])])
        tail = ax_g.transData.inverted().transform(v[np.argmin(v[:, 1])])
        assert tip[1] == pytest.approx(levels['excited'][line['excited']]['energy'], abs=1e-9)
        assert tail[1] == pytest.approx(levels['ground'][line['ground']]['energy'], abs=1e-9)
        assert tip[0] == pytest.approx(line['detuning'], abs=1e-9)
    plt.close(fig)


def test_the_figure_leaves_out_lines_outside_its_detuning_range():
    import matplotlib.pyplot as plt
    fig = EL.level_diagram('Rb87', 'D1', 6000., det_range=(-5., 5.), n_points=200)
    drawn = [a.line['detuning'] for a in fig.artists if hasattr(a, 'line')]
    assert drawn and all(-5 <= d <= 5 for d in drawn)
    plt.close(fig)


def test_zero_field_and_no_lines():
    import warnings
    import matplotlib.pyplot as plt
    with warnings.catch_warnings():
        warnings.simplefilter('error', UserWarning)     # 'identical low and high xlims'
        fig = EL.level_diagram('Rb87', 'D1', 0., n_points=200)
    assert fig.axes[1].get_xlim()[1] > 0
    plt.close(fig)
    with pytest.raises(ValueError, match='min_strength'):
        EL.level_diagram('Rb87', 'D1', 6000., min_strength=1.5)


def _sticks(ax):
    """(x, top) of each stick drawn on ``ax``."""
    from matplotlib.collections import LineCollection
    out = []
    for c in ax.collections:
        if isinstance(c, LineCollection):
            out += [(seg[0][0], max(seg[0][1], seg[1][1])) for seg in c.get_segments()]
    return sorted(out)


def test_the_stick_spectrum_shows_every_line_on_a_log_scale():
    """Under the levels: every line in the range, however weak (the arrows
    keep min_strength), at its detuning and up to its relative strength, as in
    the student's stick plots. Na D1 at 2440 G has 4 + 4 strong lines and 12
    weak ones down to 1.2e-3."""
    import matplotlib.pyplot as plt
    fig = EL.level_diagram('Na', 'D1', 2440., n_points=200)
    ax = fig.axes[5]
    assert ax.get_yscale() == 'log' and ax.get_ylim()[0] == pytest.approx(1e-4)
    every = EL.transitions('Na', 'D1', 2440., kinds=('sigma+', 'sigma-'))
    assert len(every) == 20
    np.testing.assert_allclose(_sticks(ax), sorted((l['detuning'], l['relative']) for l in every))
    arrows = [a for a in fig.artists if hasattr(a, 'line')]
    assert len(arrows) == 12                      # min_strength 0.01 drops the 1e-3 lines
    plt.close(fig)


def test_the_stick_spectrum_alone_and_in_a_range():
    import matplotlib.pyplot as plt
    ax = EL.stick_spectrum('Rb87', 'D1', 6000., det_range=(-5., 5.))
    assert all(-5 <= x <= 5 for x, _ in _sticks(ax)) and _sticks(ax)
    assert ax.get_xlim() == (-5., 5.)
    plt.close(ax.figure)
