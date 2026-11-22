"""
Energy levels, transitions and level diagrams: what a spectrum is made of.

For one isotope and D line in a magnetic field, this module gives

* ``energy_levels``: the eigenstates of the ground and excited Hamiltonians
  that the spectra are computed from (``libs.EigenSystem``), each with its
  energy, its components in the |mJ, mI> basis and the |F, mF> state it joins
  at zero field;
* ``transitions``: the sigma+, sigma- and pi lines between them, with their
  detunings and strengths, as ``spectra.FreqStren`` computes them;
* ``breit_rabi``: the levels against field, from zero to the field asked for;
* ``level_diagram``: a figure of all three with the transmission spectrum
  above: the levels against field on the left, the levels at the field across
  the detuning axis, an arrow for each transition at its detuning, and a line
  for each on the spectrum; below, the stick spectrum (``stick_spectrum``):
  every line's relative strength on a log scale.

The field defines the quantisation axis, so the levels and lines are the same
for any field direction; only which lines a beam drives depends on the
geometry (sigma+ and sigma- along the field, pi as well across it).

Run as a script for a figure and the table of lines::

    python -m elecsus.energy_levels Rb87 D1 6000 --T 75 --out rb87_d1.pdf
"""
import argparse
from fractions import Fraction
from math import factorial, sqrt

import numpy as np

from .libs import AtomConstants as AC
from .libs import EigenSystem as ES
from .libs.ang_mom import jx, jy, jz

# isotope -> (element, the isotope fractions of a cell holding only it)
ISOTOPES = {
    'Na': ('Na', {}),
    'K39': ('K', {'K40frac': 0., 'K41frac': 0.}),
    'K40': ('K', {'K40frac': 100., 'K41frac': 0.}),
    'K41': ('K', {'K40frac': 0., 'K41frac': 100.}),
    'Rb85': ('Rb', {'rb85frac': 100.}),
    'Rb87': ('Rb', {'rb85frac': 0.}),
    'Cs': ('Cs', {}),
}
KINDS = ('sigma+', 'sigma-', 'pi')
_DELTA_ML = {'sigma+': 1, 'pi': 0, 'sigma-': -1}
_SYMBOL = {'sigma+': 'σ+', 'sigma-': 'σ−', 'pi': 'π'}
_TEX = {'sigma+': r'$\sigma_+$', 'sigma-': r'$\sigma_-$', 'pi': r'$\pi$'}
_STATE = {'D1': ('S1/2', 'P1/2'), 'D2': ('S1/2', 'P3/2')}

# A line weaker than this (cleb^2) is left out, as spectra.FreqStren does.
_CLEB2_MIN = 0.0005
# A component of a state smaller than this is left out of its label.
_LABEL_MIN = 0.1


def _atom(isotope):
    if isotope not in ISOTOPES:
        raise ValueError(f'isotope must be one of {", ".join(ISOTOPES)}, not {isotope!r}')
    return getattr(AC, isotope)


def _check_line(Dline):
    if Dline not in ('D1', 'D2'):
        raise ValueError(f"Dline must be 'D1' or 'D2', not {Dline!r}")


def _ms(j):
    """m = j, j - 1, ..., -j: the order of the ang_mom matrices."""
    return [j - k for k in range(int(round(2 * j)) + 1)]


def _basis(L, I):
    """(mL, mS, mI) of each basis state, in the order sz_lsi builds them
    (kron of L, S and I, each m from +j down to -j)."""
    return [(mL, mS, mI) for mL in _ms(L) for mS in _ms(0.5) for mI in _ms(I)]


def _cg(j1, m1, j2, m2, j, m):
    """The Clebsch-Gordan coefficient <j1 m1; j2 m2 | j m> (Condon-Shortley)."""
    if abs(m1 + m2 - m) > 1e-9 or not abs(j1 - j2) - 1e-9 <= j <= j1 + j2 + 1e-9:
        return 0.
    if abs(m1) > j1 or abs(m2) > j2 or abs(m) > j:
        return 0.
    r = lambda x: int(round(x))
    pre = sqrt((2 * j + 1) * factorial(r(j + j1 - j2)) * factorial(r(j - j1 + j2))
               * factorial(r(j1 + j2 - j)) / factorial(r(j1 + j2 + j + 1)))
    pre *= sqrt(factorial(r(j + m)) * factorial(r(j - m)) * factorial(r(j1 - m1))
                * factorial(r(j1 + m1)) * factorial(r(j2 - m2)) * factorial(r(j2 + m2)))
    total = 0.
    for k in range(0, r(j1 + j2 - j) + 1):
        terms = (k, r(j1 + j2 - j - k), r(j1 - m1 - k), r(j2 + m2 - k),
                 r(j - j2 + m1 + k), r(j - j1 - m2 + k))
        if min(terms) < 0:
            continue
        total += (-1) ** k / np.prod([factorial(t) for t in terms])
    return pre * total


def _jmj_basis(L, I):
    """Columns: the states |J, mJ, mI> in the product basis of ``_basis``;
    and their (J, mJ, mI)."""
    basis = _basis(L, I)
    cols, labels = [], []
    for J in ([0.5] if L == 0 else [L - 0.5, L + 0.5]):
        for mJ in _ms(J):
            for mI in _ms(I):
                col = [(_cg(L, mL, 0.5, mS, J, mJ) if mi == mI else 0.)
                       for mL, mS, mi in basis]
                cols.append(col)
                labels.append((J, mJ, mI))
    return np.array(cols).T, labels


def _f_squared(L, I):
    """F^2 = (L + S + I)^2 in the product basis."""
    def op(make, j, which):
        dims = [int(round(2 * L + 1)), 2, int(round(2 * I + 1))]
        parts = [np.identity(d) for d in dims]
        parts[which] = make(j)
        return np.kron(np.kron(parts[0], parts[1]), parts[2])
    total = 0
    for make in (jx, jy, jz):
        F = sum(op(make, j, k) for k, j in enumerate((L, 0.5, I)))
        total = total + F @ F
    return np.real(total)


def _hamiltonians(isotope, Dline, Bfield):
    """The ground and excited Hamiltonian matrices (MHz), rebuilt from the
    eigensystem that ``spectra.calc_chi`` uses."""
    system = ES.Hamiltonian(isotope, Dline, 1.0, float(Bfield))
    out = []
    for manifold in (system.groundManifold, system.excitedManifold):
        rows = np.array(manifold).real
        E, V = rows[:, 0], rows[:, 1:].T
        out.append((V * E) @ V.T)
    return out


def _eigen(H, mF):
    """Eigenstates of H block by block in mF (the field is along z, so mF is a
    good quantum number at any field). Returns energies, vectors (columns),
    mF and the rank of each state among those of its mF, all sorted by energy."""
    E, V, M, K = [], [], [], []
    for m in sorted(set(mF)):
        idx = np.flatnonzero(np.isclose(mF, m))
        e, v = np.linalg.eigh(H[np.ix_(idx, idx)])
        full = np.zeros((len(H), len(idx)))
        full[idx] = v
        E += list(e)
        V += list(full.T)
        M += [m] * len(idx)
        K += list(range(len(idx)))
    order = np.argsort(E, kind='stable')
    return (np.array(E)[order], np.array(V)[order].T,
            np.array(M)[order], np.array(K)[order])


def _frac(m):
    return str(Fraction(m).limit_denominator(2))


def _ket(*ms, tex=False):
    if tex:   # {-}: a sign, not mathtext's spaced binary minus
        inner = ', '.join(_frac(m).replace('-', '{-}') for m in ms)
        return rf'$|{inner}\rangle$'
    return f"|{', '.join(_frac(m) for m in ms)}⟩"


class _Manifold:
    """One manifold (ground, or the excited states of one D line) at one field."""

    def __init__(self, isotope, Dline, which, Bfield, H=None):
        atom = _atom(isotope)
        self.L = 0 if which == 'ground' else 1
        self.I = atom.I
        basis = _basis(self.L, self.I)
        self.mL = np.array([b[0] for b in basis])
        mF = np.array([sum(b) for b in basis])
        if H is None:
            H = _hamiltonians(isotope, Dline, Bfield)[0 if which == 'ground' else 1]
        E, V, M, K = _eigen(H, mF)
        if which == 'excited':
            # the J = 1/2 states are the lowest (the fine structure is THz)
            n = len(basis) // 3
            keep = slice(0, n) if Dline == 'D1' else slice(n, None)
            E, V, M = E[keep], V[:, keep], M[keep]
            # rank among the line's own states of the same mF (E is sorted)
            K = np.array([np.sum(np.isclose(M[:i], M[i])) for i in range(len(M))])
        self.E, self.V, self.mF, self.rank = E, V, M, K
        self.J = 0.5 if (which == 'ground' or Dline == 'D1') else 1.5


def _zero_field_F(isotope, Dline, which):
    """F of each (mF, rank) state at zero field."""
    m = _Manifold(isotope, Dline, which, 0.)
    F2 = _f_squared(m.L, m.I)
    out = {}
    for k in range(len(m.E)):
        v = m.V[:, k]
        f2 = v @ F2 @ v
        out[(round(2 * m.mF[k]), m.rank[k])] = (-1 + sqrt(1 + 4 * f2)) / 2
    return out, float(np.mean(m.E))


def _components(manifold, k):
    """[(amplitude, mJ, mI)] of state k, largest first, without components
    below _LABEL_MIN; the overall sign makes the largest one positive."""
    C, labels = _jmj_basis(manifold.L, manifold.I)
    amps = C.T @ manifold.V[:, k]
    order = np.argsort(-np.abs(amps))
    sign = np.sign(amps[order[0]]) or 1.
    return [(float(sign * amps[i]), labels[i][1], labels[i][2]) for i in order
            if abs(amps[i]) >= _LABEL_MIN and np.isclose(labels[i][0], manifold.J)]


def _label(level, basis='mJmI', tex=False):
    if basis == 'FmF':
        return _ket(level['F'], level['mF'], tex=tex)
    comps = level['components']
    if len(comps) == 1 and abs(comps[0][0]) > 0.995:
        return _ket(comps[0][1], comps[0][2], tex=tex)
    parts = []
    for i, (a, mJ, mI) in enumerate(comps):
        sign = ('−' if a < 0 else '') if i == 0 else (' − ' if a < 0 else ' + ')
        ket = _ket(mJ, mI, tex=tex)
        parts.append(f'{sign}{abs(a):.2f} {ket}')
    return ''.join(parts)


def energy_levels(isotope, Dline, Bfield):
    """The ground and excited levels of ``isotope`` (Na, K39, K40, K41, Rb85,
    Rb87 or Cs) on ``Dline`` in a field of ``Bfield`` gauss.

    Returns ``{'ground': [...], 'excited': [...]}``, each level a dict, lowest
    energy first, with
      ``energy``     GHz, from the centre of gravity of its manifold;
      ``absolute``   MHz, as the spectra use it (a line's detuning is the
                     excited state's minus the ground state's);
      ``F``, ``mF``  the state it joins at zero field (mF is exact at any field);
      ``components`` [(amplitude, mJ, mI)], largest first, those of 0.1 and above;
      ``label``      the components as text, e.g. '0.99 |1/2, 1/2⟩ + 0.14 |-1/2, 3/2⟩'.
    """
    _check_line(Dline)
    _atom(isotope)
    H = _hamiltonians(isotope, Dline, Bfield)
    out = {}
    for which, Hm in (('ground', H[0]), ('excited', H[1])):
        m = _Manifold(isotope, Dline, which, Bfield, H=Hm)
        F_of, centre = _zero_field_F(isotope, Dline, which)
        levels = []
        for k in range(len(m.E)):
            level = {'absolute': float(m.E[k]), 'energy': float(m.E[k] - centre) / 1e3,
                     'mF': float(m.mF[k]),
                     'F': round(2 * F_of[(round(2 * m.mF[k]), m.rank[k])]) / 2,
                     'components': _components(m, k), 'index': k}
            level['label'] = _label(level)
            levels.append(level)
        out[which] = levels
    out.update(isotope=isotope, Dline=Dline, Bfield=float(Bfield))
    return out


def transitions(isotope, Dline, Bfield, kinds=KINDS, min_strength=0.):
    """The lines of ``isotope`` on ``Dline`` at ``Bfield`` gauss, strongest
    first, as spectra.FreqStren computes them (without the Boltzmann factor,
    which is below 1e-2 at any field up to 1 T at room temperature).

    Each line is a dict: ``kind`` ('sigma+', 'sigma-' or 'pi', for
    Delta mL = +1, -1, 0), ``detuning`` (GHz, as on a spectrum's axis),
    ``strength`` (FreqStren's, cleb^2 / 3 / ground-state count),
    ``relative`` (to the strongest line of the kinds asked for), and ``ground`` and
    ``excited`` (indices into ``energy_levels``' lists). Lines with
    ``relative`` below ``min_strength`` are left out.
    """
    _check_line(Dline)
    _atom(isotope)
    bad = set(kinds) - set(KINDS)
    if bad:
        raise ValueError(f'kinds must be among {KINDS}, not {sorted(bad)}')
    H = _hamiltonians(isotope, Dline, Bfield)
    g = _Manifold(isotope, Dline, 'ground', Bfield, H=H[0])
    e = _Manifold(isotope, Dline, 'excited', Bfield, H=H[1])
    ds = len(g.E)
    lines = []
    for kind in KINDS:
        block = e.V[e.mL == _DELTA_ML[kind]]
        cleb = g.V.T @ block                      # (ground, excited)
        for gi, ei in zip(*np.nonzero(cleb ** 2 > _CLEB2_MIN)):
            lines.append({'kind': kind, 'ground': int(gi), 'excited': int(ei),
                          'detuning': float(e.E[ei] - g.E[gi]) / 1e3,
                          'strength': float(cleb[gi, ei] ** 2 / 3 / ds)})
    lines = [line for line in lines if line['kind'] in kinds]
    top = max((line['strength'] for line in lines), default=1.)
    for line in lines:
        line['relative'] = line['strength'] / top
    lines = [line for line in lines if line['relative'] >= min_strength]
    return sorted(lines, key=lambda line: -line['strength'])


def breit_rabi(isotope, Dline, Bfield, n=121):
    """The levels from zero field to ``Bfield`` gauss: ``{'B': (n,), 'ground':
    (n, levels), 'excited': (n, levels)}``, energies in GHz from each
    manifold's centre of gravity. Column k follows level k of
    ``energy_levels(isotope, Dline, Bfield)`` back to zero field (the field
    keeps mF, and levels of the same mF never cross)."""
    _check_line(Dline)
    B = np.linspace(0., float(Bfield), int(n))
    out = {'B': B}
    for which in ('ground', 'excited'):
        final = _Manifold(isotope, Dline, which, Bfield)
        keys = [(round(2 * m), r) for m, r in zip(final.mF, final.rank)]
        centre = _zero_field_F(isotope, Dline, which)[1]
        rows = []
        for b in B:
            m = _Manifold(isotope, Dline, which, b)
            at = {(round(2 * mf), r): E for mf, r, E in zip(m.mF, m.rank, m.E)}
            rows.append([(at[key] - centre) / 1e3 for key in keys])
        out[which] = np.array(rows)
    return out


def spectrum(isotope, Dline, Bfield, X_GHz, T=75., lcell=2e-3, geometry='faraday',
             natural=False, outputs=('S0',)):
    """The transmission (or other ``outputs``) of a cell of ``isotope`` (or of
    the element in natural abundance) on the detunings ``X_GHz``. 'faraday':
    the field along the beam, light linear along x (it drives sigma+ and
    sigma-); 'voigt': the field along x, light linear at 45 deg to it (pi and
    sigma)."""
    from .elecsus_methods import calculate
    elem, fractions = ISOTOPES[isotope]
    p = dict(Elem=elem, Dline=Dline, T=float(T), lcell=float(lcell), Bfield=float(Bfield))
    if not natural:
        p.update(fractions)
    if geometry == 'faraday':
        E_in = [1., 0., 0.]
    elif geometry == 'voigt':
        p.update(Btheta=np.pi / 2, Bphi=0.)
        E_in = [1 / sqrt(2), 1 / sqrt(2), 0.]
    else:
        raise ValueError(f"geometry must be 'faraday' or 'voigt', not {geometry!r}")
    return [np.real(s) for s in calculate(np.asarray(X_GHz) * 1e3, E_in, p, outputs=list(outputs))]


def default_kinds(geometry):
    return ('sigma+', 'sigma-') if geometry == 'faraday' else KINDS


def in_range(lines, lo, hi):
    """The lines whose detuning is inside [lo, hi] (GHz): the ones a figure
    with that detuning axis can show."""
    return [line for line in lines if lo <= line['detuning'] <= hi]


def field_axis_end(Bfield):
    """The right end of the levels-against-field axis: the field, or 1 G at
    zero field (an axis from 0 to 0 cannot be drawn)."""
    return max(float(Bfield), 1.)


def detuning_range(lines, margin=0.15):
    """A symmetric detuning range (GHz) that holds every line, with a margin."""
    reach = max(abs(line['detuning']) for line in lines)
    reach = reach + max(1.5, margin * reach)
    return -reach, reach


def _spread(y, gap, lo=-np.inf, hi=np.inf):
    """Positions near ``y`` (sorted ascending), at least ``gap`` apart and
    within [lo, hi]: labels that would overlap are grouped and each group is
    centred on its levels."""
    y = np.asarray(y, float)
    groups = [[i] for i in range(len(y))]

    def place(group):
        centre = np.mean(y[group])
        return centre + gap * (np.arange(len(group)) - (len(group) - 1) / 2)

    merged = True
    while merged:
        merged = False
        for k in range(len(groups) - 1):
            if place(groups[k])[-1] + gap > place(groups[k + 1])[0] + 1e-12:
                groups[k:k + 2] = [groups[k] + groups[k + 1]]
                merged = True
                break
    out = np.concatenate([place(g) for g in groups]) if groups else y
    # then inside [lo, hi]: up from the bottom, down from the top
    for i in range(len(out)):
        out[i] = max(out[i], lo + i * gap)
    for i in range(len(out) - 1, -1, -1):
        out[i] = min(out[i], hi - (len(out) - 1 - i) * gap)
    return out


def stick_floor(lines):
    """The bottom of a stick spectrum's log axis: the decade below its
    weakest line, and no higher than 1e-4."""
    weakest = min(line['relative'] for line in lines)
    return min(1e-4, 10 ** np.floor(np.log10(weakest)))


def _draw_sticks(ax, lines, kinds, colour, tex=True):
    """Draw ``lines`` as sticks from the axis floor up to their relative
    strength, log scale, weakest first."""
    floor = stick_floor(lines)
    for kind in kinds:
        mine = sorted((l for l in lines if l['kind'] == kind), key=lambda l: l['relative'])
        if not mine:
            continue
        x = [l['detuning'] for l in mine]
        y = [l['relative'] for l in mine]
        ax.vlines(x, floor, y, color=colour[kind], lw=1.6,
                  label=_TEX[kind] if tex else _SYMBOL[kind])
        ax.plot(x, y, 'o', ms=3, color=colour[kind])
    ax.set_yscale('log')
    ax.set_ylim(floor, 2.)
    ax.set_ylabel('Relative strength')
    ax.set_xlabel('Detuning (GHz)')


def stick_spectrum(isotope, Dline, Bfield, geometry='faraday', kinds=None, det_range=None,
                   ax=None):
    """Every line's strength, relative to the strongest of the kinds shown, as
    a stick at its detuning on a log scale (matplotlib): down to the cut-off
    the spectra use (cleb^2 = 5e-4), so the weak lines that the field mixes
    in can be read off. Draws into ``ax`` if given; returns the Axes."""
    import matplotlib.pyplot as plt
    from .libs.durhamcolours import d_blue, d_purple, d_red

    kinds = default_kinds(geometry) if kinds is None else tuple(kinds)
    lines = transitions(isotope, Dline, Bfield, kinds=kinds)
    lo, hi = det_range if det_range is not None else detuning_range(lines)
    lines = in_range(lines, lo, hi)
    if ax is None:
        ax = plt.figure(figsize=(8, 3.5)).add_subplot()
        ax.set_title(f'{isotope} {Dline}, B = {Bfield:g} G')
    _draw_sticks(ax, lines, kinds, {'sigma+': d_blue, 'sigma-': d_purple, 'pi': d_red})
    ax.set_xlim(lo, hi)
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    return ax


def level_diagram(isotope='Rb87', Dline='D1', Bfield=6000., T=75., lcell=2e-3,
                  geometry='faraday', kinds=None, labels='mJmI', min_strength=0.01,
                  det_range=None, natural=False, n_points=2000, ax_title=True, sticks=True):
    """A figure of the levels, the transitions between them and the spectrum
    they make (matplotlib). The arrows and the marks on the spectrum are drawn
    darker for stronger lines; lines weaker than ``min_strength`` of the
    strongest are left out. ``labels`` 'mJmI' names each level by its
    components in |mJ, mI>; 'FmF' by the |F, mF> state it joins at zero
    field. ``sticks`` adds the stick spectrum underneath: every line, however
    weak, on a log scale (``stick_spectrum``). Returns the Figure."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import ConnectionPatch
    from .libs.durhamcolours import d_blue, d_purple, d_red, d_black, d_grey

    kinds = default_kinds(geometry) if kinds is None else tuple(kinds)
    colour = {'sigma+': d_blue, 'sigma-': d_purple, 'pi': d_red}
    levels = energy_levels(isotope, Dline, Bfield)
    lines = transitions(isotope, Dline, Bfield, kinds=kinds, min_strength=min_strength)
    if not lines:
        raise ValueError(f'no line is as strong as min_strength = {min_strength} of the '
                         'strongest: lower it (1 keeps only the strongest)')
    curves = breit_rabi(isotope, Dline, Bfield)
    lo, hi = det_range if det_range is not None else detuning_range(lines)
    every = in_range(transitions(isotope, Dline, Bfield, kinds=kinds), lo, hi)
    lines = in_range(lines, lo, hi)       # the others would be drawn outside the axes
    X = np.linspace(lo, hi, int(n_points))
    [S0] = spectrum(isotope, Dline, Bfield, X, T=T, lcell=lcell, geometry=geometry,
                    natural=natural)

    fig = plt.figure(figsize=(11.5, 9.6 if sticks else 7.5))
    # the stick spectrum gets a block of its own, below the field axis's label
    outer = fig.add_gridspec(2 if sticks else 1, 1, height_ratios=[3.8, 1.1] if sticks else [1],
                             hspace=0.2, left=0.08, right=0.70,
                             top=0.88 if sticks else 0.86, bottom=0.06 if sticks else 0.08)
    grid = outer[0].subgridspec(3, 2, width_ratios=[1, 5.5], height_ratios=[1, 1.4, 1.4],
                                wspace=0.0, hspace=0.12)
    ax_s = fig.add_subplot(grid[0, 1])
    ax_be, ax_e = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    ax_bg, ax_g = fig.add_subplot(grid[2, 0], sharex=ax_be), fig.add_subplot(grid[2, 1], sharex=ax_e)
    ax_e.sharex(ax_s)
    ax_e.sharey(ax_be)
    ax_g.sharey(ax_bg)
    ax_k = (fig.add_subplot(outer[1].subgridspec(1, 2, width_ratios=[1, 5.5], wspace=0.0)[0, 1],
                            sharex=ax_e) if sticks else None)

    ax_s.plot(X, S0, color=d_black, lw=1)
    ax_s.set_xlim(lo, hi)
    ax_s.set_ylabel('Transmission')
    ax_s.xaxis.set_label_position('top')
    ax_s.xaxis.tick_top()
    ax_s.set_xlabel('Detuning (GHz)')
    for kind, x in (('sigma-', 0.06), ('sigma+', 0.94), ('pi', 0.5)):
        if kind in kinds:
            ax_s.text(x, 0.12 if kind != 'pi' else 0.88, _TEX[kind], transform=ax_s.transAxes,
                      color=colour[kind], ha='center', fontsize=13)

    B = curves['B']
    for which, ax_b, ax in (('excited', ax_be, ax_e), ('ground', ax_bg, ax_g)):
        ax_b.plot(B, curves[which], color=d_black, lw=0.9)
        ax_b.set_xlim(0, field_axis_end(Bfield))
        for level in levels[which]:
            ax.axhline(level['energy'], color=d_black, lw=0.9)
        ax.tick_params(left=False, labelleft=False)
        ax.spines['left'].set_visible(False)
        ax_b.spines['right'].set_visible(False)
        e = [level['energy'] for level in levels[which]]
        all_e = np.concatenate([curves[which].ravel(), e])
        pad = 0.08 * (all_e.max() - all_e.min() or 1.)
        ax_b.set_ylim(all_e.min() - pad, all_e.max() + pad)
    state = _STATE[Dline]
    n = {'Na': 3, 'K39': 4, 'K40': 4, 'K41': 4, 'Rb85': 5, 'Rb87': 5, 'Cs': 6}[isotope]
    ax_be.set_ylabel(f'{n}{state[1][0]}$_{{{state[1][1:]}}}$ energy/h (GHz)')
    ax_bg.set_ylabel(f'{n}{state[0][0]}$_{{{state[0][1:]}}}$ energy/h (GHz)')
    ax_bg.set_xlabel('Magnetic field (G)')
    ax_be.tick_params(labelbottom=False, bottom=False)
    ax_e.tick_params(labelbottom=False, bottom=False)
    ax_g.tick_params(labelbottom=False, bottom=False)
    ax_e.spines['bottom'].set_visible(False)
    ax_g.spines['top'].set_visible(False)
    ax_be.spines['bottom'].set_visible(False)
    ax_bg.spines['top'].set_visible(False)
    if sticks:
        _draw_sticks(ax_k, every, kinds, colour)
        ax_k.set_xlim(lo, hi)
        ax_k.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=11)

    for line in sorted(lines, key=lambda l: l['relative']):
        g = levels['ground'][line['ground']]['energy']
        e = levels['excited'][line['excited']]['energy']
        alpha = 0.3 + 0.7 * line['relative']
        c = colour[line['kind']]
        ax_s.axvline(line['detuning'], color=c, lw=1, alpha=alpha)
        # a filled arrow with no outline and no shrink, so that its tip is
        # exactly on the excited level and its tail on the ground level (a
        # stroked head overshoots by about half its line width)
        style = (f"simple,head_length=0.9,head_width={0.55 + 0.25 * line['relative']:.2f},"
                 f"tail_width={0.12 + 0.2 * line['relative']:.2f}")
        arrow = ConnectionPatch((line['detuning'], g), (line['detuning'], e), 'data', 'data',
                                axesA=ax_g, axesB=ax_e, arrowstyle=style, mutation_scale=10,
                                shrinkA=0, shrinkB=0, lw=0, facecolor=c, edgecolor='none',
                                alpha=alpha)
        arrow.line = line
        fig.add_artist(arrow)

    # labels to the right of the levels, spread so that they do not overlap
    fig.canvas.draw()
    for which, ax in (('excited', ax_e), ('ground', ax_g)):
        ordered = sorted(levels[which], key=lambda level: level['energy'])
        ylo, yhi = ax.get_ylim()
        per_point = (yhi - ylo) / (ax.get_window_extent().height * 72 / fig.dpi)
        # 7.5 pt where there is room; smaller, down to 4 pt, where there is not
        gap = min(1.25 * 7.5 * per_point, 0.98 * (yhi - ylo) / len(ordered))
        fs = max(4., gap / (1.25 * per_point))
        y = _spread([level['energy'] for level in ordered], gap, ylo + gap / 2, yhi - gap / 2)
        for level, yl in zip(ordered, y):
            ax.annotate(_label(level, labels, tex=True), xy=(1.0, level['energy']),
                        xycoords=('axes fraction', 'data'), xytext=(1.02, yl),
                        textcoords=('axes fraction', 'data'), fontsize=fs, va='center',
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle='-', lw=0.4, color=d_grey,
                                        shrinkA=0, shrinkB=0))
    head = r'$|m_J, m_I\rangle$' if labels == 'mJmI' else r'$|F, m_F\rangle$ at $B = 0$'
    ax_s.text(1.02, 0.0, head, transform=ax_s.transAxes, fontsize=9, va='bottom')
    if ax_title:
        wavelength = getattr(AC, {'Na': 'Na', 'K39': 'K', 'K40': 'K', 'K41': 'K', 'Rb85': 'Rb',
                                  'Rb87': 'Rb', 'Cs': 'Cs'}[isotope] + Dline + 'Transition').wavelength
        mass = ''.join(ch for ch in isotope if ch.isdigit())
        sym = ''.join(ch for ch in isotope if ch.isalpha())
        name = f'$^{{{mass}}}${sym}' if mass else sym
        cell = 'natural abundance' if natural else 'isotopically pure'
        fig.suptitle(f'{name} {Dline} transitions, {wavelength * 1e9:.0f} nm: B = {Bfield:g} G '
                     f'({Bfield / 1e4:.2f} T); spectrum at {T:g} °C, {lcell * 1e3:g} mm, {cell}',
                     x=0.08, y=0.97, ha='left', fontsize=12)
    return fig


def table(lines, levels, basis='mJmI'):
    """The lines as rows of text: kind, detuning, relative strength, from, to."""
    rows = []
    for line in sorted(lines, key=lambda line: line['detuning']):
        g = _label(levels['ground'][line['ground']], basis)
        e = _label(levels['excited'][line['excited']], basis)
        rows.append((_SYMBOL[line['kind']], f"{line['detuning']:+.3f}",
                     f"{line['relative']:.3f}", g, e))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='python -m elecsus.energy_levels',
        description='Energy levels, transitions and spectrum of one isotope in a field.')
    parser.add_argument('isotope', choices=list(ISOTOPES))
    parser.add_argument('Dline', choices=['D1', 'D2'])
    parser.add_argument('Bfield', type=float, help='field (G)')
    parser.add_argument('--T', type=float, default=75., help='cell temperature (C), default 75')
    parser.add_argument('--lcell', type=float, default=2., help='cell length (mm), default 2')
    parser.add_argument('--geometry', choices=['faraday', 'voigt'], default='faraday')
    parser.add_argument('--labels', choices=['mJmI', 'FmF'], default='mJmI')
    parser.add_argument('--min-strength', type=float, default=0.01,
                        help='leave out lines weaker than this fraction of the strongest')
    parser.add_argument('--natural', action='store_true',
                        help='spectrum of the element in natural abundance')
    parser.add_argument('--no-sticks', action='store_true',
                        help='leave out the stick spectrum under the levels')
    parser.add_argument('--out', help='save the figure here (.pdf, .png, ...)')
    a = parser.parse_args(argv)
    kinds = default_kinds(a.geometry)
    levels = energy_levels(a.isotope, a.Dline, a.Bfield)
    lines = transitions(a.isotope, a.Dline, a.Bfield, kinds=kinds, min_strength=a.min_strength)
    print(f'{a.isotope} {a.Dline}, B = {a.Bfield:g} G: {len(lines)} lines')
    print(f"{'':3}{'GHz':>9}  {'strength':>8}  from -> to")
    for kind, det, rel, g, e in table(lines, levels, a.labels):
        print(f'{kind:3}{det:>9}  {rel:>8}  {g}  ->  {e}')
    if a.out:
        import matplotlib
        matplotlib.use('Agg')
        fig = level_diagram(a.isotope, a.Dline, a.Bfield, T=a.T, lcell=a.lcell / 1e3,
                            geometry=a.geometry, labels=a.labels, min_strength=a.min_strength,
                            natural=a.natural, sticks=not a.no_sticks)
        fig.savefig(a.out, dpi=200)
        print(f'wrote {a.out}')


if __name__ == '__main__':
    main()
