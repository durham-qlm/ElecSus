# Copyright 2026 S. A. Wrathmall.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
ElecSus Dash GUI  [ES4 NEW]
===========================

Cross-platform web-based GUI for ElecSus, new in ES4.
Replaces dependency on wxPython with Plotly Dash (runs in any browser).

Features:
  - Single Cell tab: forward theory calculation with all parameters, and
    fitting to uploaded data (or the bundled sample spectrum)
  - Cascaded Two Cells tab: sequential propagation through two cells
  - Non-orthogonal propagation (Briscoe et al. 2025) and the optional
    four-interface cell-window Fresnel model
  - Buffer gas parameters exposed (BufferGas, BufferPressure, BufferTemp_fill)
  - Independent Doppler temperature for dual-temperature cells
  - Interactive Plotly plots (one panel per unit), Durham University light theme
  - CSV export

Usage:
    python -m elecsus.elecsus_dash_gui
    # Then open http://127.0.0.1:8050 in a browser

Requirements:
    pip install -e ".[gui]"
"""

import base64
import io
import os
import re
import threading
import warnings

import numpy as np
import dash
from dash import dcc, html, Input, Output, State, no_update
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from elecsus import __version__
from elecsus import energy_levels as EL
from elecsus.elecsus_methods import calculate, fit_data, _fit_form
from elecsus.libs.param_validation import FIT_LIMITS
from elecsus.libs.BufferGasData import BUFFER_GAS_DATA
from elecsus.libs import AtomConstants, spectra

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
ELEMENTS   = ['Rb', 'Cs', 'K', 'Na']
DLINES     = ['D1', 'D2']
# Output order: Stokes/intensities, then the circular-basis quantities paired
# sigma+ then sigma-. [ES4 FIX] The API keys now carry the physical
# sigma too: 'TransLeft' and 'nPlus' are sigma+ (the [1,+i] mode for B along
# +z). In ES3 and the 4.0.0 alphas the 'Plus'/'Left' keys held sigma-, and the
# labels crossed over.
OUTPUTS    = ['S0', 'S0_P', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45',
              'TransLeft', 'TransRight', 'alphaPlus', 'alphaMinus', 'alphaZ',
              'nPlus', 'nMinus', 'phi', 'psi', 'DoLP']
OUTPUT_LABELS = {
    'S0': 'S₀ (Transmission)', 'S1': 'S₁', 'S2': 'S₂', 'S3': 'S₃',
    # S0_P: Poynting-flux transmission of the bare vapour column. It needs the
    # non-orthogonal solver and, unlike S0, excludes the window factors -- with
    # windows on the two differ by ~15.5% even far from resonance.
    'S0_P': 'S₀ Poynting (vapour only, no windows)',
    'Ix': 'Iₓ', 'Iy': 'Iᵧ', 'I_P45': 'I₊₄₅', 'I_M45': 'I₋₄₅',
    # Labelled by physical sigma: 'TransLeft' (=exp(-2 Im(nPlus) k L)) is the
    # sigma+ transmission, 'TransRight' the sigma-.
    'TransLeft': 'Trans σ⁺', 'TransRight': 'Trans σ⁻',
    'alphaPlus': 'α σ⁺', 'alphaMinus': 'α σ⁻', 'alphaZ': 'α π',
    # n(sigma±): physically meaningful only for Btheta=0 (Faraday); computed
    # unconditionally. Key 'nPlus' is the sigma+ index, 'nMinus' the sigma-.
    'nPlus': 'n(σ⁺) — Faraday only', 'nMinus': 'n(σ⁻) — Faraday only',
    'phi': 'Faraday rotation φ — dispersion only',
    # psi/DoLP: Stokes-based (dichroism-aware) generalisation of phi.
    # psi = 0.5*atan2(S2,S1); only meaningful where DoLP is appreciable.
    'psi': 'Rotation ψ (Stokes, incl. dichroism)',
    'DoLP': 'Degree of linear pol.',
}
# Outputs the core computes only with use_non_orthogonal=True.
NONORTHO_ONLY_OUTPUTS = {'S0_P'}
# One plot panel per unit, so an absorption coefficient in m^-1 or an angle in
# radians does not flatten a 0-1 transmission on a shared axis. Anything not
# listed here is dimensionless and goes in the first panel.
OUTPUT_GROUPS = [
    ('Signal', None),
    ('Rotation (rad)', {'phi', 'psi'}),
    ('Absorption coeff. (m⁻¹)', {'alphaPlus', 'alphaMinus', 'alphaZ'}),
    ('Refractive index', {'nMinus', 'nPlus'}),
    # [ES4 NEW] the cascade tab's lower panel: cell 2 on its own
    ('Into and through cell 2', {'cell2_into', 'cell2_T'}),
]
# [ES4 NEW] x axes the plots can be drawn against. Astronomers give filter
# profiles in wavelength: the offset from line centre in mÅ, or the air
# wavelength in Å. The detuning boxes stay in GHz (see _x_on). The offset in
# air and in vacuum differ by 0.03 %, so the mÅ axis does not name either.
X_UNITS = {'GHz': 'Detuning (GHz)', 'mA': 'Δλ from line centre (mÅ)',
           'A': 'Wavelength in air (Å)'}
# The air-wavelength axis ('A') is hidden for now; the
# code behind it stays, and uncommenting its option brings it back.
X_UNIT_OPTIONS = [{'label': 'Detuning [GHz]', 'value': 'GHz'},
                  {'label': 'Δλ [mÅ]', 'value': 'mA'},
                  # {'label': 'λ air [Å]', 'value': 'A'},
                  ]
X_HOVER_FORMAT = {'GHz': '', 'mA': '.2f', 'A': '.5f'}
# plotly's own ticks on the Å axis read "7,698.6" and "7699"
X_TICK_FORMAT = {'GHz': '', 'mA': '', 'A': '.3~f'}
# Default window index offered when the window model is switched on
# (borosilicate glass, as in the window derivation in the manual).
N_WINDOW_DEFAULT = 1.51
# [ES4 MOD] The non-orthogonal solver is on by default, as in the
# core API: the orthogonal (ES3) path is wrong in oblique fields. This is the
# switch's initial checklist value, and what the callbacks assume when they
# are called without it.
NONORTHO_DEFAULT = ('on',)
BUFFER_GAS_OPTIONS = ['None', 'Ne', 'Ar', 'Kr', 'Xe', 'N2', 'He']
POL_OPTIONS = ['Linear (x)', 'Linear (y)', 'Left CP', 'Right CP', 'Linear 45°']

# [ES4 NEW] Fitting panel — data-type dropdown reuses OUTPUTS above.
# The RR fitting algorithm uses multiprocessing.Pool internally, which is fragile
# when invoked from inside a running Dash server process, so the GUI only exposes
# the two algorithms that run in-process (ML directly, DE via scipy through lmfit).
FIT_ALGO_OPTIONS = [
    {'label': 'Marquardt-Levenberg', 'value': 'ML'},
    {'label': 'Differential Evolution', 'value': 'differential_evolution'},
]
# (p_dict key, display label) for parameters the GUI allows the user to vary.
# The buffer-gas pressure can vary while a gas is selected; the fill
# temperature cannot, as the model depends on the ratio of the two (see the
# fit_data() docstring).
# 'pol_angle' is not a p_dict key: the GUI fits the linear input angle through
# the E_x/E_y field amplitudes (see _pol_fit_setup).
VARY_PARAM_LABELS = [
    ('T', 'T'), ('Bfield', 'B'), ('GammaBuf', 'GammaBuf'), ('shift', 'Shift'),
    ('BufferPressure', 'Buffer P'),
    ('lcell', 'l_cell'), ('Btheta', 'Bθ'), ('Bphi', 'Bφ'),
    ('DoppTemp', 'Doppler T'), ('rb85frac', '⁸⁵Rb %'), ('pol_angle', 'Pol. θ₀'),
]
LINEAR_POLS = ('Linear (x)', 'Linear (y)', 'Linear 45°')

# Colours — Durham University light theme
BG        = '#ffffff'           # page / plot background
PANEL     = '#f6f2f5'           # sidebar panels (light Palatinate tint)
BORDER    = '#d4c6d3'           # panel borders
ACCENT    = '#00638A'           # d_blue — secondary buttons
HIGHLIGHT = '#7E317B'           # d_purple (Palatinate) — primary accent
TEXT      = '#231F20'           # d_black — body text
SUBTEXT   = '#68625F'           # mid-grey — secondary / annotation text
PLOT_BG   = '#ffffff'           # plot area background
WARN_COLOUR = '#AA2B4A'         # d_red — notes under a fit result
LEVELS_HEIGHT = 1150            # px: the energy-level figure is drawn to a fixed height

# Durham colour cycle — matches plot_style.py (d_yellow omitted; too pale on white)
PALETTE = [
    '#7E317B',  # d_purple (Palatinate)
    '#00638A',  # d_blue   (Durham Blue)
    '#AA2B4A',  # d_red
    '#A277A6',  # d_lightpurple
    '#00848E',  # d_midblue (teal)
    '#B44876',  # d_pink
    '#231F20',  # d_black
    '#716A64',  # mid-grey
]

PLOT_LAYOUT = dict(
    paper_bgcolor=PLOT_BG, plot_bgcolor=PLOT_BG,
    font=dict(color=TEXT, family='Times New Roman, Times, serif', size=13),
    xaxis=dict(gridcolor='#ece8ec', zerolinecolor='#b0a8b0', linecolor='#b0a8b0',
               showgrid=False),
    yaxis=dict(gridcolor='#ece8ec', zerolinecolor='#b0a8b0', linecolor='#b0a8b0',
               showgrid=False),
    legend=dict(bgcolor='rgba(255,255,255,0.9)', bordercolor=BORDER, borderwidth=1),
    margin=dict(l=60, r=20, t=50, b=60),
    hovermode='x unified',
)

DROPDOWN_STYLE = {
    'backgroundColor': '#ffffff', 'color': '#000',
    'border': f'1px solid {BORDER}', 'borderRadius': '4px', 'fontSize': '0.85rem',
}
INPUT_STYLE = {
    'width': '100%', 'backgroundColor': '#ffffff', 'color': TEXT,
    'border': f'1px solid {BORDER}', 'borderRadius': '4px',
    'padding': '5px 8px', 'fontSize': '0.84rem',
}
INPUT_DISABLED_STYLE = {
    **INPUT_STYLE,
    'backgroundColor': '#ede8ec', 'color': '#aaa',
    'cursor': 'not-allowed',
}
BTN = {
    'backgroundColor': HIGHLIGHT, 'color': '#fff', 'border': 'none',
    'borderRadius': '6px', 'padding': '8px 20px', 'fontWeight': '700',
    'cursor': 'pointer', 'fontSize': '0.9rem', 'letterSpacing': '0.04em',
}
BTN_SEC = {**BTN, 'backgroundColor': ACCENT}


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def lbl(text, sub=None):
    children = [html.Span(text, style={'color': TEXT, 'fontSize':'0.82rem','fontWeight':'500'})]
    if sub:
        children.append(html.Span(f' {sub}', style={'color':SUBTEXT,'fontSize':'0.73rem'}))
    return html.Label(children, style={'marginBottom':'2px','display':'block'})

def ni(id_, val):
    # No step, min or max: the browser marks a value outside any of them
    # invalid, and Dash then sends None -- so B = 239 G (step 10) computed at
    # B = 0, T = 119.69 at 20 C, and Btheta = 200 at 0 (min/max).
    # step='any' keeps every value valid. Limits are enforced by the callbacks
    # (see _num and _X) and the core's param_validation, which say so.
    return dcc.Input(id=id_, type='number', value=val,
                     step='any', debounce=True, style=INPUT_STYLE)

def row(*kids, gap='8px'):
    # flex-end: when one field's label wraps (narrow cascade columns), the
    # inputs still line up along the bottom of the row
    return html.Div(kids, style={'display':'flex','gap':gap,'alignItems':'flex-end','marginBottom':'6px'})

def field(label_text, ctrl, sub=None, flex='1'):
    return html.Div([lbl(label_text, sub), ctrl], style={'flex':flex,'minWidth':'70px'})

def sec(text):
    return html.Div([
        html.Span(text, style={'color':HIGHLIGHT,'fontWeight':'700','fontSize':'0.72rem',
                               'letterSpacing':'0.1em','textTransform':'uppercase'}),
        html.Hr(style={'borderColor':BORDER,'margin':'3px 0 8px 0'}),
    ])


def empty_fig(msg='Set parameters and click Compute'):
    fig = go.Figure()
    fig.update_layout(**PLOT_LAYOUT,
        annotations=[dict(text=msg, xref='paper', yref='paper',
                          x=0.5, y=0.5, showarrow=False,
                          font=dict(size=16, color=SUBTEXT))])
    return fig


def _num(value, label, cell=None):
    """A number from a required field. A field that is empty, or holds text
    that is not a number, reaches the callbacks as None; 0 is a value like any
    other (it used to be read as blank). Raises ValueError naming
    the field, which the callbacks show in the error toast rather than
    computing with a default the user did not choose."""
    if value is None or isinstance(value, (bool, str)):
        where = f'{cell}: ' if cell else ''
        raise ValueError(f'{where}{label} is empty or not a number: enter a value.')
    return float(value)


def _require_isotopes(elem, rb85frac, K40frac, K41frac):
    """The isotope fields of the element's own row, which are shown only for
    Rb and K; for the other elements they are not used and may be anything."""
    if elem == 'Rb':
        _num(rb85frac, '⁸⁵Rb (%)')
    elif elem == 'K':
        _num(K40frac, '⁴⁰K (%)')
        _num(K41frac, '⁴¹K (%)')


def _pol_efield(pol, theta0_deg):
    # theta0 is the rotation of the two linear presets only
    t = np.radians(_num(theta0_deg, 'θ₀ [deg]')) if pol in ('Linear (x)', 'Linear (y)') else 0.
    m = {
        'Linear (x)':  [np.cos(t), np.sin(t), 0.],
        'Linear (y)':  [-np.sin(t), np.cos(t), 0.],
        'Left CP':     [1/np.sqrt(2), 1j/np.sqrt(2), 0.],
        'Right CP':    [1/np.sqrt(2), -1j/np.sqrt(2), 0.],
        'Linear 45°':  [1/np.sqrt(2), 1/np.sqrt(2), 0.],
    }
    return m.get(pol, [1., 0., 0.])


def _on(checklist_value):
    """A single-option dcc.Checklist is ticked when its value contains 'on'."""
    return 'on' in (checklist_value or [])


def _buffer_model_on(bufgas, elem, dline, bufP):
    """Whether the core's buffer-gas model is used for these settings, as in
    spectra.calc_chi: a gas with data for this element and line, and a
    pressure above 0. Only then does the model replace GammaBuf; otherwise
    the static GammaBuf is what the core uses."""
    if bufgas in (None, 'None') or (elem, dline, bufgas) not in BUFFER_GAS_DATA:
        return False
    return isinstance(bufP, (int, float)) and not isinstance(bufP, bool) and bufP > 0


# The core's warnings, which the page would otherwise never show.
# warnings.catch_warnings swaps process-wide state, so callbacks running in
# two server threads at once could lose each other's warnings, or leave the
# filters changed; the lock lets one at a time use it.
_WARNINGS_LOCK = threading.Lock()
_NOT_FOR_USERS = (DeprecationWarning, PendingDeprecationWarning, ImportWarning,
                  ResourceWarning)


class _Collected:
    """`with _collect_warnings() as notices:` records the warnings raised in
    the block; `notices.messages` then lists what the user should see, once
    each, in the order raised. Nothing reaches the server's stderr."""

    def __enter__(self):
        _WARNINGS_LOCK.acquire()
        self._ctx = warnings.catch_warnings(record=True)
        self._caught = self._ctx.__enter__()
        warnings.simplefilter('always')
        return self

    def __exit__(self, *exc):
        try:
            self._ctx.__exit__(*exc)
        finally:
            _WARNINGS_LOCK.release()
        return False

    @property
    def messages(self):
        out = []
        for w in self._caught:
            msg = str(w.message)
            if not issubclass(w.category, _NOT_FOR_USERS) and msg not in out:
                out.append(msg)
        return out


def _collect_warnings():
    return _Collected()


def _p(elem, dline, T, lcell, B, gbuf, shift, btheta, bphi, theta0,
       bufgas, bufP, bufTf,
       rb85frac=None, K40frac=None, K41frac=None,
       nonortho=NONORTHO_DEFAULT, windows=None, n_window=None, doppT=None, cell=None):
    # Every required field must hold a number (_num); `cell` names the cell in
    # the cascade's error messages. Blank is allowed where it means something:
    # Doppler T (tied to T), the window index (borosilicate), GammaBuf while a
    # buffer-gas model supplies the broadening (a gas with data for the line and
    # P > 0), and the isotope fields not shown.
    gas = bufgas and bufgas != 'None'
    # GUI angle inputs are in degrees (matching the ES3 wx GUI and the other
    # angle controls in the panel); the core p_dict convention is radians.
    d = dict(Elem=elem, Dline=dline,
             T=_num(T, 'T [°C]', cell), lcell=_num(lcell, 'l [m]', cell),
             Bfield=_num(B, 'B [G]', cell),
             GammaBuf=(0. if gbuf is None and _buffer_model_on(bufgas, elem, dline, bufP)
                       else _num(gbuf, 'GammaBuf [MHz]', cell)),
             shift=_num(shift, 'Shift [MHz]', cell),
             Btheta=np.radians(_num(btheta, 'Bθ [deg]', cell)),
             Bphi=np.radians(_num(bphi, 'Bφ [deg]', cell)),
             Constrain=True)
    # Blank Doppler T keeps the Doppler width tied to T (the core's Constrain);
    # a value frees it, for cells whose body and stem are at different
    # temperatures (T then sets only the number density).
    if doppT not in (None, ''):
        d.update(DoppTemp=float(doppT), Constrain=False)
    if gas:
        d.update(BufferGas=bufgas,
                 BufferPressure=_num(bufP, 'P [Torr]', cell),
                 BufferTemp_fill=_num(bufTf, 'Fill T [K]', cell))
    if elem == 'Rb' and rb85frac is not None:
        d['rb85frac'] = float(rb85frac)
    if elem == 'K':
        if K40frac is not None:
            d['K40frac'] = float(K40frac)
        if K41frac is not None:
            d['K41frac'] = float(K41frac)
    # The solver choice is always sent, so the GUI does not depend on the
    # core's default. The window model lives only in the non-orthogonal
    # propagator: the orthogonal path ignores n_window, so it is sent only
    # when both are on.
    d['use_non_orthogonal'] = _on(nonortho)
    if _on(nonortho):
        if _on(windows):
            d['n_window'] = float(n_window if n_window not in (None, '')
                                  else N_WINDOW_DEFAULT)
    return d


def _model_tag(p):
    """Plot-title suffix naming the non-default model choices, so an exported
    figure says how it was computed."""
    tags = []
    if not p.get('Constrain', True):
        tags.append(f"Doppler T={p['DoppTemp']:g}°C")
    if not p.get('use_non_orthogonal', True):
        tags.append('orthogonal (ES3)')
    elif p.get('n_window') is not None:
        tags.append(f"windows n={p['n_window']:g}")
    return ''.join(f' · {t}' for t in tags)


N_POINTS_MIN, N_POINTS_MAX = 10, 5000


def _X(start, stop, N):
    """Detuning axis in MHz from the panel's GHz start and stop and point count."""
    lo, hi = _num(start, 'Start [GHz]'), _num(stop, 'Stop [GHz]')
    n = _num(N, 'N pts')
    if n != int(n) or not N_POINTS_MIN <= n <= N_POINTS_MAX:
        raise ValueError(f'N pts must be a whole number from {N_POINTS_MIN} to '
                         f'{N_POINTS_MAX}, not {n:g}.')
    return np.linspace(lo * 1e3, hi * 1e3, int(n))


# ─────────────────────────────────────────────────────────────────────────────
# Wavelength axes  [ES4 NEW]
# ─────────────────────────────────────────────────────────────────────────────
_C = 299792458.0    # speed of light, m/s


def _air_index(lambda_vac):
    """Refractive index of standard air (dry, 15 C, 101.325 kPa) at a vacuum
    wavelength in metres: the IAU standard, Edlen (1966) as given by Morton
    (2000, ApJS 130, 403). Astronomers quote air wavelengths above 2000 Å."""
    s2 = (1e-6 / lambda_vac) ** 2          # vacuum wavenumber squared, um^-2
    return 1. + 8.34213e-5 + 2.406030e-2 / (130. - s2) + 1.5997e-4 / (38.9 - s2)


def _line_centre(elem, dline):
    """(frequency in Hz, air wavelength in Å) of zero detuning: the line centre
    of the natural-abundance mixture, the transition's v0 in AtomConstants."""
    nu0 = getattr(AtomConstants, f'{elem}{dline}Transition').v0
    return nu0, 1e10 * _C / (nu0 * _air_index(_C / nu0))


def _x_on(X_GHz, unit, elem, dline):
    """Detunings X_GHz on the x axis `unit` (a key of X_UNITS): unchanged for
    'GHz', else the air wavelength in Å ('A') or its offset from line centre
    in mÅ ('mA'). The wavelength is lambda0 nu0 / (nu0 + detuning), with the
    air index at line centre (it changes by ~1e-10 across a scan). A higher
    frequency is a shorter wavelength: positive detuning is negative Δλ."""
    X = np.asarray(X_GHz, dtype=float)
    if unit not in ('mA', 'A'):
        return X
    nu0, lam0 = _line_centre(elem, dline)
    dlam = -lam0 * (X * 1e9) / (nu0 + X * 1e9)         # Å
    return dlam * 1e3 if unit == 'mA' else lam0 + dlam


def _conversion_text(elem, dline, start, stop):
    """The line under a Detuning Axis panel: mÅ per GHz, and the GHz range in
    mÅ. (The line centre's air wavelength is left out while the air axis is
    hidden; see X_UNIT_OPTIONS.)"""
    try:
        nu0, lam0 = _line_centre(elem, dline)
    except (AttributeError, TypeError):
        return ''
    per_GHz = -lam0 * 1e12 / nu0                        # mÅ per GHz at line centre
    # with the air axis: f'{elem} {dline}: λ₀ = {lam0:.3f} Å in air; ...'
    text = (f'{elem} {dline}: +1 GHz ≙ {per_GHz:.2f} mÅ, '
            f'+100 mÅ ≙ {100 / per_GHz:.3f} GHz.')
    try:
        lo, hi = _x_on([_num(start, 'Start'), _num(stop, 'Stop')], 'mA', elem, dline)
        text += f' This range is {lo:+.1f} to {hi:+.1f} mÅ.'
    except (ValueError, TypeError):
        pass
    return text


def _set_xunit(fig, X_GHz, unit, elem, dline):
    """`fig` (a Figure or its dict) with every trace drawn against the x axis
    `unit`. All the traces of one figure share the detuning grid X_GHz. A
    figure the page sends back may carry layout that plotly.js added; anything
    plotly.py does not know is dropped rather than refused."""
    fig = go.Figure(fig, skip_invalid=True)
    fig.update_traces(x=_x_on(X_GHz, unit, elem, dline))
    title = X_UNITS.get(unit, X_UNITS['GHz'])
    fig.for_each_xaxis(lambda axis: axis.update(title_text=title)
                       if axis.title.text else None)
    fig.update_xaxes(autorange=True, hoverformat=X_HOVER_FORMAT.get(unit, ''),
                     tickformat=X_TICK_FORMAT.get(unit, ''))
    return fig


def _redraw(fig, store, unit):
    """A figure the x-axis selector redraws: the last result's, on the new
    axis. Its store holds the detuning grid and the line; without a result
    the figure is left as it is."""
    if not fig or not store or 'elem' not in store:
        return no_update
    return _set_xunit(fig, store['X'], unit, store['elem'], store['dline'])


def _group_of(key):
    """Index into OUTPUT_GROUPS of the plot panel an output key belongs in."""
    for i, (_title, keys) in enumerate(OUTPUT_GROUPS):
        if keys and key in keys:
            return i
    return 0


def _fig(X_GHz, traces, keys=None):
    """Plot traces (name, y, dotted, opacity) against detuning.

    With `keys` (the output key of each trace), outputs of different units go
    in separate panels stacked on a shared detuning axis, in OUTPUT_GROUPS
    order; without, everything shares one 'Signal' axis."""
    groups = [_group_of(k) for k in keys] if keys else [0] * len(traces)
    panels = sorted(set(groups))
    row_of = {g: r + 1 for r, g in enumerate(panels)}
    fig = make_subplots(rows=len(panels), cols=1, shared_xaxes=True,
                        vertical_spacing=0.07 if len(panels) > 1 else 0)
    for i, ((name, y, dash_, alpha), g) in enumerate(zip(traces, groups)):
        fig.add_trace(go.Scatter(
            x=X_GHz, y=np.real(y),
            name=name, mode='lines', opacity=alpha,
            line=dict(color=PALETTE[i % len(PALETTE)], width=2,
                      dash='dot' if dash_ else 'solid'),
        ), row=row_of[g], col=1)
    layout = PLOT_LAYOUT.copy()
    axis_style = {k: v for k, v in layout.pop('xaxis').items()}
    layout.pop('yaxis')
    fig.update_layout(**layout)
    fig.update_xaxes(**axis_style)
    fig.update_yaxes(**axis_style)
    fig.update_xaxes(title_text='Detuning (GHz)', row=len(panels), col=1)
    for g in panels:
        fig.update_yaxes(title_text=OUTPUT_GROUPS[g][0], row=row_of[g], col=1)
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Fitting helpers  [ES4 NEW]
# ─────────────────────────────────────────────────────────────────────────────
def _fit_vary_disabled(bufgas, elem='Rb', doppT=None, pol='Linear (x)',
                       dline='D2', bufP=0.):
    """Fit parameters that cannot be varied in the current configuration:
    GammaBuf while the buffer-gas model gives the broadening, and the buffer
    pressure while it does not (fit_data raises otherwise: the model is on for
    a gas with data for this element and line, and P > 0); Doppler T unless it
    is decoupled from T (with Constrain on, the core ignores DoppTemp); the
    85Rb fraction for other elements; and the polarisation angle for
    circular input. The shift stays: the model's pressure shift is added to
    it."""
    disabled = set()
    if _buffer_model_on(bufgas, elem, dline, bufP):
        disabled.add('GammaBuf')
    else:
        disabled.add('BufferPressure')
    if doppT in (None, ''):
        disabled.add('DoppTemp')
    if elem != 'Rb':
        disabled.add('rb85frac')
    if pol not in LINEAR_POLS:
        disabled.add('pol_angle')
    return disabled


def _fit_vary_options(bufgas, elem='Rb', doppT=None, pol='Linear (x)',
                      dline='D2', bufP=0.):
    """Greyed-out rather than hidden, so the user sees why an option is off."""
    disabled = _fit_vary_disabled(bufgas, elem, doppT, pol, dline, bufP)
    return [{'label': f' {label}', 'value': key, 'disabled': key in disabled}
            for key, label in VARY_PARAM_LABELS]


def _fit_scales(p0):
    """A size for each fit parameter, in the units the fit uses: the half-width
    of the differential-evolution window when the guess is near zero, and the
    error above which the data are taken not to fix the parameter (a fitted
    value with an error estimate larger than this has not been determined)."""
    return {'T': 20.0, 'Bfield': 200.0, 'GammaBuf': 100.0, 'shift': 100.0,
            'BufferPressure': 10.0,
            'lcell': 0.5 * float(p0.get('lcell', 0.075)),
            'Btheta': 0.5, 'Bphi': 0.5, 'DoppTemp': 20.0, 'rb85frac': 10.0}


def _pol_fit_setup(E_in):
    """E_in in fit_data's [E_x, [E_y, E_phase]] form, and the amplitude to vary
    to fit a linear input angle.

    fit_data takes a real E_x, and E_y as a magnitude and a relative phase
    that the ML fit bounds to [0, pi). A relative phase in [pi, 2pi) -- right
    circular light, or a linear angle below zero -- is carried by the sign of
    E_x instead: (-a, b e^{i(phi-pi)}) = -(a, b e^{i phi}) is the same
    polarisation state. (Passing the raw phase let lmfit clamp it to the
    bound, so a right-circular fit used linear light.) That conversion is
    elecsus_methods._fit_form, which fit_data now applies itself.

    S0 and every output are normalised to the input intensity, so only the
    ratio E_y/E_x matters: to fit a linear angle one amplitude is fixed and
    the other varied. The larger component is the one fixed, so the ratio
    stays well conditioned wherever the angle starts."""
    E_fit = _fit_form(E_in)
    a, b = E_fit[0], E_fit[1][0]
    vary_key = 'E_y' if abs(a) >= b else 'E_x'
    return E_fit, vary_key


def _pol_angle_deg(E_x, E_y):
    """Linear polarisation angle from the x axis, folded into (-90, 90]."""
    ang = np.degrees(np.arctan2(E_y, E_x))
    ang = (ang + 90.0) % 180.0 - 90.0
    return 90.0 if ang == -90.0 else ang


def _fit_fig(X_GHz, Y_data, Y_fit, data_label, rms):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=X_GHz, y=np.real(Y_data), name='Data', mode='markers',
                              marker=dict(size=4, color=SUBTEXT, opacity=0.6)))
    fig.add_trace(go.Scatter(x=X_GHz, y=np.real(Y_fit), name=f'Fit ({data_label})',
                              mode='lines', line=dict(color=HIGHLIGHT, width=2)))
    layout = PLOT_LAYOUT.copy()
    layout['xaxis'] = {**layout['xaxis'], 'title': 'Detuning (GHz)'}
    layout['yaxis'] = {**layout['yaxis'], 'title': 'Signal'}
    fig.update_layout(**layout, title=dict(
        text=f'Fit result — RMS = {rms:.5g}', font=dict(color=TEXT, size=13)))
    return fig


def _split_row(line):
    """The fields of one line: comma-separated if it has a comma, else
    separated by tabs or spaces."""
    if ',' in line:
        return [f.strip() for f in line.split(',')]
    return line.split()


def _unit_scale(header, ncols):
    """(scale to GHz, how it was decided) from the header line.

    Only the first column's label is read, so a signal column that mentions
    'MHz' (a bandwidth, say) does not rescale the detuning. In a
    comma- or tab-separated header that label is the text before the first
    delimiter; in a whitespace-separated one it is all but the last ncols - 1
    words."""
    label = ''
    if header:
        if ',' in header:
            label = header.split(',', 1)[0]
        elif '\t' in header:
            label = header.split('\t', 1)[0]
        else:
            words = header.split()
            label = ' '.join(words[:max(len(words) - (ncols - 1), 1)])
    label = label.lower()
    if re.search(r'mhz', label) and not re.search(r'ghz', label):
        return 1e-3, 'MHz, from the header'
    if re.search(r'ghz', label):
        return 1.0, 'GHz, from the header'
    return 1.0, 'GHz, the default: the header names no unit'


def _read_fit_csv(contents):
    """Parse an uploaded 2-column file (detuning, signal) into
    (X_ghz, Y, how), `how` saying which unit the detuning was read in.

    Lines before the first row of numbers are the header (and may be any
    number of them, as may `#` comment lines); the delimiter is a comma, tabs
    or spaces. Detuning is taken as GHz unless the first column's header names
    MHz, as the bundled sample (elecsus/sample_data) and the scope2elecsus
    files do. Columns after the second are ignored."""
    _, content_string = contents.split(',', 1)
    text = base64.b64decode(content_string).decode('utf-8', errors='replace').lstrip('\ufeff')
    rows, header, comment = [], None, None
    for number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            if line.startswith('#'):
                raise ValueError
            rows.append([float(f) for f in _split_row(line)])
        except ValueError:
            if rows and not line.startswith('#'):
                raise ValueError(f'Could not parse file as a 2-column CSV (detuning, '
                                 f'signal): line {number} is not numbers: {line[:30]!r}')
            if not rows:
                if line.startswith('#'):
                    comment = line.lstrip('#').strip()
                else:
                    header = line
    if not rows:
        raise ValueError('Could not parse file as a 2-column CSV (detuning, signal)')
    if len({len(r) for r in rows}) != 1:
        raise ValueError('Could not parse file as a 2-column CSV (detuning, signal): '
                         'the rows do not all have the same number of columns')
    arr = np.array(rows)
    if arr.shape[1] < 2:
        raise ValueError('CSV must have two columns: detuning (GHz), signal')
    scale, how = _unit_scale(header or comment, arr.shape[1])
    return arr[:, 0] * scale, arr[:, 1], how


def _parse_fit_csv(contents):
    """(X_ghz, Y) from an uploaded file; see _read_fit_csv."""
    X_ghz, Y, _ = _read_fit_csv(contents)
    return X_ghz, Y


# ─────────────────────────────────────────────────────────────────────────────
# Reusable param panel  (prefix determines all ids)
# ─────────────────────────────────────────────────────────────────────────────
def _prop_control_state(nonortho, windows):
    """(windows checklist options, n_window disabled, n_window style)."""
    windows_label = html.Span(' Cell windows (Fresnel)',
                              style={'color': TEXT if _on(nonortho) else SUBTEXT,
                                     'fontSize': '0.82rem'})
    options = [{'label': windows_label, 'value': 'on', 'disabled': not _on(nonortho)}]
    nwin_off = not (_on(nonortho) and _on(windows))
    return options, nwin_off, INPUT_DISABLED_STYLE if nwin_off else INPUT_STYLE


def prop_section(prefix):
    """Propagation model: the non-orthogonal solver and the window model.

    Defaults match the core API: the solver on, the windows off. The window
    controls take effect only with the non-orthogonal solver, since the
    orthogonal path ignores n_window; they are disabled while it is unticked
    (see _prop_control_state)."""
    note = {'color': SUBTEXT, 'fontSize': '0.72rem', 'margin': '0 0 6px 0'}
    return [
        sec('Propagation'),
        dcc.Checklist(
            id=f'{prefix}-nonortho',
            options=[{'label': html.Span(' Non-orthogonal modes (Briscoe 2025)',
                                         style={'color': TEXT, 'fontSize': '0.82rem'}),
                      'value': 'on'}],
            value=list(NONORTHO_DEFAULT), inputStyle={'accentColor': HIGHLIGHT}),
        html.P('On by default. Off is ElecSus 3\'s orthogonal method: exact '
               'only for B along or across the beam, and wrong in oblique '
               'fields (by up to ~0.5 in S₀).', style=note),
        dcc.Checklist(
            id=f'{prefix}-windows',
            options=_prop_control_state(NONORTHO_DEFAULT, [])[0],
            value=[], inputStyle={'accentColor': HIGHLIGHT}),
        row(field('Window n', dcc.Input(
                id=f'{prefix}-nwin', type='number', value=N_WINDOW_DEFAULT,
                step='any', debounce=True, disabled=True,
                style=INPUT_DISABLED_STYLE),
            sub='1.51 borosilicate')),
        html.P('Four interfaces, normal incidence, single pass: '
               'the window etalon is not modelled.', style=note),
    ]


def param_panel(prefix, show_pol_selector=True):
    return html.Div([
        sec('Element & Line'),
        row(
            field('Element', dcc.Dropdown(ELEMENTS, value='Rb', id=f'{prefix}-elem',
                             clearable=False, style=DROPDOWN_STYLE)),
            field('D-line',  dcc.Dropdown(DLINES,   value='D2', id=f'{prefix}-dline',
                             clearable=False, style=DROPDOWN_STYLE)),
        ),
        sec('Cell Parameters'),
        row(
            field('T [°C]',     ni(f'{prefix}-T',    20.0)),
            field('l [m]',      ni(f'{prefix}-lcell', 0.075)),
        ),
        row(
            field('B [G]',      ni(f'{prefix}-B',    0.0)),
            field('GammaBuf [MHz]', ni(f'{prefix}-gbuf', 0.0)),
        ),
        row(
            field('Shift [MHz]', ni(f'{prefix}-shift', 0.0)),
            field('Doppler T [°C]', ni(f'{prefix}-doppT', None),
                  sub='blank = T'),
        ),
        sec('Magnetic Field Angles'),
        row(
            field('Bθ [deg]', ni(f'{prefix}-btheta', 0.0),
                  sub='0=Faraday  90=Voigt'),
            field('Bφ [deg]', ni(f'{prefix}-bphi',   0.0)),
        ),
        *prop_section(prefix),
        *([
            sec('Polarisation'),
            row(
                field('Pol.', dcc.Dropdown(POL_OPTIONS, value='Linear (x)',
                              id=f'{prefix}-pol', clearable=False, style=DROPDOWN_STYLE)),
                field('θ₀ [deg]', ni(f'{prefix}-theta0', 0.0),
                      sub='linear angle'),
            ),
        ] if show_pol_selector else [
            html.Div(id=f'{prefix}-pol', style={'display':'none'}),
            html.Div(id=f'{prefix}-theta0', style={'display':'none'}),
        ]),
        sec('Buffer Gas'),
        row(
            field('Gas', dcc.Dropdown(BUFFER_GAS_OPTIONS, value='None',
                         id=f'{prefix}-bufgas', clearable=False, style=DROPDOWN_STYLE)),
            field('P [Torr]', ni(f'{prefix}-bufP', 0.0)),
        ),
        row(
            field('Fill T [K]', ni(f'{prefix}-bufTf', 293.15),
                  sub='sealed cell'),
        ),
        # Isotope abundance controls — shown only for elements with non-trivial mixtures
        html.Div(id=f'{prefix}-iso-rb-row', style={'display':'none'}, children=[
            sec('Isotope Abundance'),
            row(field('⁸⁵Rb (%)', ni(f'{prefix}-rb85frac', 72.17),
                      sub='natural 72.17')),
        ]),
        html.Div(id=f'{prefix}-iso-k-row', style={'display':'none'}, children=[
            sec('Isotope Abundances'),
            row(
                field('⁴⁰K (%)', ni(f'{prefix}-K40frac', 0.01),
                      sub='natural 0.01'),
                field('⁴¹K (%)', ni(f'{prefix}-K41frac', 6.73),
                      sub='natural 6.73'),
            ),
        ]),
    ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
              'border':f'1px solid {BORDER}','marginBottom':'10px'})


def det_panel(prefix):
    return html.Div([
        sec('Detuning Axis'),
        row(
            field('Start [GHz]', ni(f'{prefix}-det-start', -5.0)),
            field('Stop [GHz]',  ni(f'{prefix}-det-stop',   5.0)),
            field('N pts',       ni(f'{prefix}-det-N', 500)),
        ),
        # [ES4 NEW] plot against wavelength; the boxes above stay in GHz
        row(field('Plot against', dcc.Dropdown(X_UNIT_OPTIONS, value='GHz',
                                               id=f'{prefix}-xunit', clearable=False,
                                               style=DROPDOWN_STYLE))),
        html.Div(id=f'{prefix}-det-conv',
                 style={'color': SUBTEXT, 'fontSize': '0.72rem', 'marginTop': '2px'}),
    ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
              'border':f'1px solid {BORDER}','marginBottom':'10px'})


def _output_options(opts, nonortho_on):
    """Checklist options; outputs that need the non-orthogonal solver are
    greyed out while it is off (the core raises for them otherwise)."""
    return [{'label': html.Span(OUTPUT_LABELS[o],
                                style={'color': TEXT, 'fontSize': '0.82rem',
                                       'marginLeft': '4px'}),
             'value': o,
             'disabled': (o in NONORTHO_ONLY_OUTPUTS) and not nonortho_on}
            for o in opts]


def out_panel(prefix, allowed=None):
    opts = allowed if allowed is not None else OUTPUTS
    return html.Div([
        sec('Outputs to Plot'),
        dcc.Checklist(
            id=f'{prefix}-outputs',
            options=_output_options(opts, nonortho_on=_on(NONORTHO_DEFAULT)),
            value=['S0'], inline=False,
            style={'columnCount':2, 'padding':'4px 0'},
            inputStyle={'accentColor': HIGHLIGHT},
        ),
    ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
              'border':f'1px solid {BORDER}','marginBottom':'10px'})


# Cascade instrument totals are computed from the final E-field, so only
# field-derived quantities are meaningful there (no alphas/TransLeft/Right,
# which are per-cell medium properties).
CASCADE_OUTPUTS = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45']


def fit_panel(prefix):
    """[ES4 NEW] Fit-to-data panel. Entirely additive: uses only new component
    ids (prefix-fit-*), so it does not touch the existing forward-calculation
    param/det/out panels or the compute callback for this prefix."""
    return html.Div([
        sec('Fit to Experimental Data'),
        dcc.Upload(
            id=f'{prefix}-fit-upload',
            children=html.Div(['Drag & drop or ', html.A('select CSV')],
                              style={'fontSize':'0.8rem','color':SUBTEXT}),
            style={'width':'100%','padding':'10px','borderWidth':'1px','borderStyle':'dashed',
                   'borderColor':BORDER,'borderRadius':'6px','textAlign':'center',
                   'backgroundColor':'#fff','marginBottom':'6px','boxSizing':'border-box'},
            multiple=False,
        ),
        html.Button('Load sample spectrum', id=f'{prefix}-fit-sample', n_clicks=0,
                    title='Rb D2, 75 mm natural-abundance cell, zero field, ~40 °C '
                          '(Pizzey et al., New J. Phys. 24, 125001, 2022)',
                    style={**BTN_SEC, 'width': '100%', 'padding': '5px 10px',
                           'fontSize': '0.8rem', 'marginBottom': '6px'}),
        html.Div(id=f'{prefix}-fit-upload-status',
                 style={'fontSize':'0.75rem','color':SUBTEXT,'marginBottom':'8px'}),
        row(
            field('Data type', dcc.Dropdown(OUTPUTS, value='S0', id=f'{prefix}-fit-datatype',
                                            clearable=False, style=DROPDOWN_STYLE)),
            field('Algorithm', dcc.Dropdown(FIT_ALGO_OPTIONS, value='ML',
                                            id=f'{prefix}-fit-algo', clearable=False,
                                            style=DROPDOWN_STYLE)),
        ),
        row(
            field('Fit min [GHz]', ni(f'{prefix}-fit-min', None),
                  sub='blank = full range'),
            field('Fit max [GHz]', ni(f'{prefix}-fit-max', None)),
        ),
        lbl('Vary in fit', sub='greyed out: not usable with these settings'),
        dcc.Checklist(
            id=f'{prefix}-fit-vary',
            options=_fit_vary_options('None'),
            value=['T'], inline=True,
            style={'marginBottom':'4px'},
            inputStyle={'accentColor': HIGHLIGHT, 'marginRight':'3px'},
            labelStyle={'marginRight':'10px','fontSize':'0.8rem','color':TEXT},
        ),
        html.Div(style={'height':'6px'}),
        html.Button('Run Fit', id=f'{prefix}-fit-run', n_clicks=0, style={**BTN,'width':'100%'}),
        dcc.Loading(type='circle', color=HIGHLIGHT,
                    children=html.Div(id=f'{prefix}-fit-results', style={'marginTop':'8px'})),
        html.Div(style={'height':'6px'}),
        html.Button('Export fit CSV', id=f'{prefix}-fit-export', n_clicks=0,
                    style={**BTN_SEC,'width':'100%'}),
        dcc.Download(id=f'{prefix}-fit-download'),
        dcc.Store(id=f'{prefix}-fit-data-store'),
        dcc.Store(id=f'{prefix}-fit-result-store'),
    ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
              'border':f'1px solid {BORDER}','marginBottom':'10px'})


# ─────────────────────────────────────────────────────────────────────────────
# App layout — all components exist in DOM at all times (tabs hide via CSS)
# ─────────────────────────────────────────────────────────────────────────────
# The Bootstrap theme (Flatly, as dbc.themes.FLATLY) and the Inter font are
# served from elecsus/assets/, which Dash adds to the page automatically, so
# the GUI looks the same on a machine with no network.
app = dash.Dash(
    __name__,
    suppress_callback_exceptions=False,
    title='ElecSus GUI',
)

# Underlying Flask WSGI instance, exposed for production servers (gunicorn):
#   gunicorn 'elecsus.elecsus_dash_gui:server'
# The __main__ block below still uses Dash's built-in dev server for local runs.
server = app.server

# ── Header ────────────────────────────────────────────────────────────────────
_STAGES = {'a': 'alpha', 'b': 'beta', 'rc': 'release candidate',
           '.dev': 'development version'}


def _release_note(version):
    """The line under the title: the version, and for a pre-release (4.0.0b1,
    4.0.0rc2, 4.0.0.dev3) that it is not the final release. Nothing for a
    source tree that was never installed, which has no version."""
    if 'unknown' in version:
        return ''
    stage = re.match(r'\d+(?:\.\d+)*(a|b|rc|\.dev)\d', version)
    if stage:
        return f'v{version} · {_STAGES[stage[1]]}, not the final release'
    return f'v{version}'


header = html.Div([
    html.Div([
        html.H1('ElecSus', style={'color':HIGHLIGHT,'fontWeight':'700','margin':'0',
                                   'fontSize':'1.7rem','letterSpacing':'-0.02em'}),
        html.P('Atomic Spectroscopy · Durham QLM',
               style={'color':SUBTEXT,'margin':'0','fontSize':'0.82rem'}),
        html.P(_release_note(__version__), id='release-note',
               style={'color':WARN_COLOUR,'margin':'0','fontSize':'0.78rem'}),
    ]),
], style={'padding':'12px 24px','backgroundColor':PANEL,
          'borderBottom':f'2px solid {HIGHLIGHT}'})

# ── Tab bar ───────────────────────────────────────────────────────────────────
tab_bar = html.Div([
    html.Button('Single Cell',         id='tab-single-btn', n_clicks=0,
                style={**BTN_SEC,'marginRight':'4px','borderRadius':'6px 6px 0 0'}),
    html.Button('Cascaded Two Cells',  id='tab-cascade-btn', n_clicks=0,
                style={**BTN_SEC,'marginRight':'4px','borderRadius':'6px 6px 0 0'}),
    html.Button('Energy Levels',       id='tab-levels-btn', n_clicks=0,
                style={**BTN_SEC,'borderRadius':'6px 6px 0 0'}),
], style={'padding':'10px 16px 0 16px','backgroundColor':BG,
          'borderBottom':f'1px solid {BORDER}'})

# ── Single-cell tab ───────────────────────────────────────────────────────────
single_tab = html.Div([
    html.Div([
        # Left panel
        html.Div([
            det_panel('s'), out_panel('s'),
            param_panel('s'),
            html.Button('Compute', id='s-compute', n_clicks=0, style={**BTN,'width':'100%'}),
            html.Div(style={'height':'6px'}),
            html.Button('Export CSV', id='s-export', n_clicks=0,
                        style={**BTN_SEC,'width':'100%'}),
            dcc.Download(id='s-download'),
            dcc.Store(id='s-store'),
            fit_panel('s'),
        ], style={'width':'310px','minWidth':'280px','padding':'14px',
                  'overflowY':'auto','maxHeight':'calc(100vh - 110px)'}),
        # Plot
        html.Div([
            dcc.Loading(type='circle', color=HIGHLIGHT, children=[
                dcc.Graph(id='s-plot', figure=empty_fig(),
                          style={'height':'calc(100vh - 140px)'},
                          config={'toImageButtonOptions':{'format':'svg','filename':'elecsus_single'}}),
            ]),
            html.Hr(style={'borderColor':BORDER,'margin':'16px 0'}),
            html.H4('Fit Result', style={'color':HIGHLIGHT,'margin':'0 0 8px 0','fontSize':'0.95rem'}),
            dcc.Loading(type='circle', color=HIGHLIGHT, children=[
                dcc.Graph(id='s-fit-plot', figure=empty_fig('Upload data and click Run Fit'),
                          style={'height':'calc(100vh - 140px)'},
                          config={'toImageButtonOptions':{'format':'svg','filename':'elecsus_fit'}}),
            ]),
        ], style={'flex':'1','padding':'12px'}),
    ], style={'display':'flex'}),
], id='single-tab-content')

# ── Cascade tab ───────────────────────────────────────────────────────────────
cascade_tab = html.Div([
    html.Div([
        # Cell 1 panel
        html.Div([
            html.H4('Cell 1', style={'color':PALETTE[0],'margin':'0 0 8px 0','fontSize':'0.95rem'}),
            param_panel('c1'),
        ], style={'width':'265px','minWidth':'250px','padding':'14px',
                  'overflowY':'auto','maxHeight':'calc(100vh - 110px)'}),
        # Cell 2 + controls panel
        html.Div([
            html.H4('Cell 2', style={'color':PALETTE[1],'margin':'0 0 8px 0','fontSize':'0.95rem'}),
            # Cell 2 params (no separate polarisation selector — inherits from Cell 1 output)
            html.Div([
                sec('Element & Line'),
                # Same-species cascaded filters: element/D-line/isotopes are
                # locked to cell 1 (synced by callback, shown greyed out).
                row(
                    field('Element', dcc.Dropdown(ELEMENTS, value='Rb', id='c2-elem',
                                     clearable=False, disabled=True, style=DROPDOWN_STYLE)),
                    field('D-line',  dcc.Dropdown(DLINES,   value='D2', id='c2-dline',
                                     clearable=False, disabled=True, style=DROPDOWN_STYLE)),
                ),
                html.Div('matched to Cell 1 (same-species filters)',
                         style={'color':SUBTEXT,'fontSize':'0.72rem','marginTop':'-4px',
                                'marginBottom':'6px'}),
                sec('Cell Parameters'),
                row(field('T [°C]', ni('c2-T', 20.0)),
                    field('l [m]',  ni('c2-lcell', 0.075))),
                row(field('B [G]',  ni('c2-B', 0.0)),
                    field('GammaBuf [MHz]', ni('c2-gbuf', 0.0))),
                row(field('Shift [MHz]', ni('c2-shift', 0.0)),
                    field('Doppler T [°C]', ni('c2-doppT', None),
                          sub='blank = T')),
                sec('Magnetic Field Angles'),
                row(field('Bθ [deg]', ni('c2-btheta', 0.0),
                          sub='0=Faraday  90=Voigt'),
                    field('Bφ [deg]', ni('c2-bphi',   0.0))),
                *prop_section('c2'),
                sec('Buffer Gas'),
                row(field('Gas', dcc.Dropdown(BUFFER_GAS_OPTIONS, value='None', id='c2-bufgas',
                                 clearable=False, style=DROPDOWN_STYLE)),
                    field('P [Torr]', ni('c2-bufP', 0.0))),
                row(field('Fill T [K]', ni('c2-bufTf', 293.15), sub='sealed cell')),
            ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
                      'border':f'1px solid {BORDER}','marginBottom':'10px'}),

            # Intermediate polariser — default ON at 90° from the x axis, which
            # is crossed for an x-polarised input (the default; the standard MOF
            # filter-section topology) and only then.
            html.Div([
                sec('Intermediate Optics'),
                dcc.Checklist(id='cas-use-pol',
                    options=[{'label': html.Span(' Polariser between cells',
                                                  style={'color':TEXT,'fontSize':'0.83rem'}),
                              'value':'on'}],
                    value=['on'], inputStyle={'accentColor':HIGHLIGHT}),
                html.Div(style={'height':'4px'}),
                row(field('Angle from x [deg]', ni('cas-pol-angle', 90.0),
                          sub='90 = crossed for x-polarised input')),
                html.Div(style={'height':'6px'}),
                # QWP after the polariser: prepares a FIXED circular state for
                # cell 2 regardless of wavelength -- the wing-selector topology
                # of instruments like VAMOS (Oliviero 1998), where the second
                # cell must always see the same helicity to select one Zeeman
                # wing. 'None' keeps today's behaviour (no QWP).
                field('QWP after polariser', dcc.Dropdown(
                    [{'label': 'None', 'value': 'none'},
                     {'label': 'Left CP',  'value': 'left_cp'},
                     {'label': 'Right CP', 'value': 'right_cp'}],
                    value='none', id='cas-reset-pol', clearable=False,
                    style=DROPDOWN_STYLE),
                    sub='needs the polariser; fixed circular into Cell 2'),
            ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
                      'border':f'1px solid {BORDER}','marginBottom':'10px'}),

            det_panel('cas'), out_panel('cas', allowed=CASCADE_OUTPUTS),

            html.Div([
                html.Button('Compute Cascade', id='cas-compute', n_clicks=0,
                            style={**BTN,'marginRight':'6px'}),
                html.Button('Export CSV',      id='cas-export',  n_clicks=0,
                            style=BTN_SEC),
            ], style={'display':'flex','marginTop':'4px'}),
            dcc.Download(id='cas-download'),
            dcc.Store(id='cas-store'),
        ], style={'width':'278px','minWidth':'260px','padding':'14px',
                  'overflowY':'auto','maxHeight':'calc(100vh - 110px)'}),
        # Plot
        html.Div([
            dcc.Loading(type='circle', color=HIGHLIGHT, children=[
                dcc.Graph(id='cas-plot', figure=empty_fig(),
                          style={'height':'calc(100vh - 140px)'},
                          config={'toImageButtonOptions':{'format':'svg','filename':'elecsus_cascade'}}),
            ]),
        ], style={'flex':'1','padding':'12px'}),
    ], style={'display':'flex'}),
], id='cascade-tab-content', style={'display':'none'})

# ── Energy levels tab  [ES4 NEW] ─────────────────────────────────────────────
# The levels of one isotope in a field, the lines between them and the spectrum
# they make (elecsus.energy_levels): the figure students used to draw by hand.
ISOTOPE_OPTIONS = list(EL.ISOTOPES)
GEOMETRY_OPTIONS = [{'label': 'Faraday (B ∥ k): σ+, σ−', 'value': 'faraday'},
                    {'label': 'Voigt (B ⊥ k): σ+, σ−, π', 'value': 'voigt'}]
LABEL_OPTIONS = [{'label': '|mJ, mI⟩ components', 'value': 'mJmI'},
                 {'label': '|F, mF⟩ at zero field', 'value': 'FmF'}]
_LEVELS_NOTE = {'color': SUBTEXT, 'fontSize': '0.72rem', 'margin': '0 0 6px 0'}

levels_tab = html.Div([
    html.Div([
        # Left panel
        html.Div([
            html.Div([
                sec('Atom'),
                row(
                    field('Isotope', dcc.Dropdown(ISOTOPE_OPTIONS, value='Rb87', id='lv-iso',
                                                  clearable=False, style=DROPDOWN_STYLE)),
                    field('D-line', dcc.Dropdown(DLINES, value='D1', id='lv-dline',
                                                 clearable=False, style=DROPDOWN_STYLE)),
                ),
                sec('Field'),
                row(field('B [G]', ni('lv-B', 6000.0))),
                sec('Lines'),
                row(field('Geometry', dcc.Dropdown(GEOMETRY_OPTIONS, value='faraday',
                                                   id='lv-geometry', clearable=False,
                                                   style=DROPDOWN_STYLE))),
                row(field('Min. strength', ni('lv-min', 0.01), sub='of the strongest')),
                html.P('Lines weaker than this are left out. Arrows and marks are darker '
                       'for stronger lines.', style=_LEVELS_NOTE),
                row(field('Label levels by', dcc.Dropdown(LABEL_OPTIONS, value='mJmI',
                                                          id='lv-labels', clearable=False,
                                                          style=DROPDOWN_STYLE))),
                sec('Spectrum'),
                row(
                    field('T [°C]', ni('lv-T', 75.0)),
                    field('l [m]', ni('lv-lcell', 0.002)),
                ),
                row(
                    field('Start [GHz]', ni('lv-det-start', None)),
                    field('Stop [GHz]', ni('lv-det-stop', None)),
                ),
                html.P('Blank: a range that holds every line.', style=_LEVELS_NOTE),
                dcc.Checklist(
                    id='lv-natural',
                    options=[{'label': html.Span(' Natural abundance (arrows: this isotope only)',
                                                 style={'color': TEXT, 'fontSize': '0.82rem'}),
                              'value': 'on'}],
                    value=[], inputStyle={'accentColor': HIGHLIGHT}),
            ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
                      'border':f'1px solid {BORDER}','marginBottom':'10px'}),
            html.Button('Compute', id='lv-compute', n_clicks=0, style={**BTN,'width':'100%'}),
            html.Div(style={'height':'6px'}),
            html.Button('Export lines CSV', id='lv-export', n_clicks=0,
                        style={**BTN_SEC,'width':'100%'}),
            dcc.Download(id='lv-download'),
            dcc.Store(id='lv-store'),
        ], style={'width':'290px','minWidth':'260px','padding':'14px',
                  'overflowY':'auto','maxHeight':'calc(100vh - 110px)'}),
        # Plot and table of lines
        html.Div([
            dcc.Loading(type='circle', color=HIGHLIGHT, children=[
                dcc.Graph(id='lv-plot', figure=empty_fig(), style={'height': f'{LEVELS_HEIGHT}px'},
                          config={'toImageButtonOptions':{'format':'svg',
                                                          'filename':'elecsus_levels'}}),
            ]),
            html.Div(id='lv-table', style={'padding':'0 12px'}),
        ], style={'flex':'1','padding':'12px','overflowX':'auto'}),
    ], style={'display':'flex'}),
], id='levels-tab-content', style={'display':'none'})

# ── Error and warning toasts ──────────────────────────────────────────────────
# An error stops the calculation; a warning is the core saying that the result
# is not quite what was asked for. Both sit in one fixed column.
error_toast = dbc.Toast(id='error-toast', header='Calculation Error',
                        is_open=False, dismissable=True, duration=7000, icon='danger',
                        style={'maxWidth':'420px'})
warning_toast = dbc.Toast(id='warning-toast', header='Warning',
                          is_open=False, dismissable=True, duration=15000, icon='warning',
                          style={'maxWidth':'420px'})
toasts = html.Div([warning_toast, error_toast],
                  style={'position':'fixed','bottom':'20px','right':'20px','zIndex':9999,
                         'display':'flex','flexDirection':'column','gap':'8px'})

app.layout = html.Div([
    header, tab_bar,
    html.Div([single_tab, cascade_tab, levels_tab],
             style={'backgroundColor':BG,'minHeight':'calc(100vh - 110px)'}),
    toasts,
    dcc.Store(id='error-store', data=''),
    dcc.Store(id='warning-store', data=[]),
], style={'backgroundColor':BG,'minHeight':'100vh',
          'fontFamily':'Inter,Arial,sans-serif','margin':'0'})


# ─────────────────────────────────────────────────────────────────────────────
# Tab switching
# ─────────────────────────────────────────────────────────────────────────────
TABS = ('single', 'cascade', 'levels')


@app.callback(
    *[Output(f'{t}-tab-content', 'style') for t in TABS],
    *[Output(f'tab-{t}-btn', 'style') for t in TABS],
    *[Input(f'tab-{t}-btn', 'n_clicks') for t in TABS],
)
def switch_tab(*clicks):
    trigger = dash.callback_context.triggered_id
    active = next((t for t in TABS if trigger == f'tab-{t}-btn'), 'single')
    shape = {'borderRadius':'6px 6px 0 0'}
    gap = {'marginRight':'4px'}
    contents = [{'display':'block' if t == active else 'none'} for t in TABS]
    buttons = [{**(BTN if t == active else BTN_SEC), **shape,
                **(gap if t != TABS[-1] else {}),
                'borderBottom': f'2px solid {HIGHLIGHT}' if t == active else 'none'}
               for t in TABS]
    return (*contents, *buttons)


# ─────────────────────────────────────────────────────────────────────────────
# Single-cell compute
# ─────────────────────────────────────────────────────────────────────────────
@app.callback(
    Output('s-plot',  'figure'),
    Output('s-store', 'data'),
    Output('error-store', 'data'),
    # [ES4] appended last, so the indices above stay put
    Output('warning-store', 'data'),
    Input('s-compute', 'n_clicks'),
    State('s-elem','value'),  State('s-dline','value'),
    State('s-T','value'),     State('s-lcell','value'),
    State('s-B','value'),     State('s-gbuf','value'),
    State('s-shift','value'),
    State('s-btheta','value'),State('s-bphi','value'),
    State('s-theta0','value'),State('s-pol','value'),
    State('s-bufgas','value'),State('s-bufP','value'),State('s-bufTf','value'),
    State('s-rb85frac','value'),State('s-K40frac','value'),State('s-K41frac','value'),
    State('s-det-start','value'),State('s-det-stop','value'),State('s-det-N','value'),
    State('s-outputs','value'),
    # [ES4] appended last, with defaults below, so positional callers keep working
    State('s-nonortho','value'), State('s-windows','value'), State('s-nwin','value'),
    State('s-doppT','value'),
    State('s-xunit','value'),
    prevent_initial_call=True,
)
def compute_single(n, elem, dline, T, lcell, B, gbuf, shift,
                   btheta, bphi, theta0, pol,
                   bufgas, bufP, bufTf,
                   rb85frac, K40frac, K41frac,
                   det_start, det_stop, det_N, outputs,
                   nonortho=NONORTHO_DEFAULT, windows=None, n_window=None, doppT=None,
                   xunit='GHz'):
    # On a failure the data store is emptied (None), not left as it was: Export
    # CSV would otherwise write the previous run's numbers under a page that
    # says the calculation failed.
    with _collect_warnings() as caught:
        try:
            _require_isotopes(elem, rb85frac, K40frac, K41frac)
            p = _p(elem, dline, T, lcell, B, gbuf, shift,
                   btheta, bphi, theta0, bufgas, bufP, bufTf,
                   rb85frac=rb85frac, K40frac=K40frac, K41frac=K41frac,
                   nonortho=nonortho, windows=windows, n_window=n_window, doppT=doppT)
            E_in = _pol_efield(pol, theta0)
            X = _X(det_start, det_stop, det_N)
            outs = outputs or ['S0']
            res = calculate(X, E_in, p, outputs=outs)
            traces = [(OUTPUT_LABELS.get(k,k), np.real(v), False, 1.0)
                      for k, v in zip(outs, res)]
            fig = _set_xunit(_fig(X/1e3, traces, keys=outs), X/1e3, xunit, elem, dline)
            fig.update_layout(title=dict(
                # the values used, which are the typed ones (_num refuses a blank)
                text=f"{elem} {dline} · T={p['T']:g}°C · B={p['Bfield']:g} G · "
                     f"l={p['lcell']:g} m" + _model_tag(p),
                font=dict(color=TEXT, size=13)))
            # elem and dline let the x-axis selector redraw it (_redraw)
            store = {'X': list(X/1e3), 'results': {k: list(np.real(v)) for k,v in zip(outs,res)},
                     'elem': elem, 'dline': dline}
            result = (fig, store, '')
        except Exception as e:
            result = (empty_fig('Calculation failed'), None, str(e))
    return (*result, caught.messages)


# ─────────────────────────────────────────────────────────────────────────────
# Fit to experimental data  [ES4 NEW]
# All new component ids (s-fit-*) — does not touch compute_single or its
# Output/State signature, so the existing forward-calculation path is untouched.
# ─────────────────────────────────────────────────────────────────────────────
@app.callback(
    Output('s-fit-data-store', 'data'),
    Output('s-fit-upload-status', 'children'),
    Input('s-fit-upload', 'contents'),
    State('s-fit-upload', 'filename'),
    prevent_initial_call=True,
)
def parse_fit_upload(contents, filename):
    if not contents:
        return no_update, no_update
    try:
        X_ghz, Y, how = _read_fit_csv(contents)
        store = {'X': X_ghz.tolist(), 'Y': Y.tolist()}
        # the unit and range, so a detuning read in the wrong unit shows up here
        return store, (f'{filename}: {len(X_ghz)} points loaded, {X_ghz.min():.4g} to '
                       f'{X_ghz.max():.4g} GHz (detuning read as {how})')
    except Exception as e:
        return no_update, f'Failed to parse {filename}: {e}'


@app.callback(
    Output('s-fit-vary', 'options'),
    Input('s-bufgas', 'value'),
    Input('s-elem', 'value'),
    Input('s-doppT', 'value'),
    Input('s-pol', 'value'),
    Input('s-dline', 'value'),
    Input('s-bufP', 'value'),
)
def update_fit_vary_options(bufgas, elem='Rb', doppT=None, pol='Linear (x)',
                            dline='D2', bufP=0.):
    return _fit_vary_options(bufgas, elem, doppT, pol, dline, bufP)


# The bundled sample spectrum (see elecsus/sample_data/README.md): Rb D2 through
# a 75 mm natural-abundance cell in zero field at ~40 C, detuning in MHz.
SAMPLE_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'sample_data', 'rb_d2_doppler_40C.csv')


# What 'Load sample spectrum' puts in the other boxes: the conditions the sample
# was taken in (a natural-abundance Rb D2 cell in zero field, no buffer gas, the
# light x-polarised), so that Run Fit gives the documented result whatever the
# panel held before. T is left alone: it is what the fit finds.
SAMPLE_PANEL = [
    # (component id, value)
    ('s-gbuf', 0.0), ('s-shift', 0.0), ('s-btheta', 0.0), ('s-bphi', 0.0),
    ('s-doppT', None), ('s-pol', 'Linear (x)'), ('s-theta0', 0.0),
    ('s-bufgas', 'None'), ('s-bufP', 0.0), ('s-bufTf', 293.15), ('s-rb85frac', 72.17),
    ('s-nonortho', list(NONORTHO_DEFAULT)), ('s-windows', []),
    ('s-fit-datatype', 'S0'), ('s-fit-min', None), ('s-fit-max', None),
    ('s-fit-vary', ['T']),
]


@app.callback(
    Output('s-fit-data-store', 'data', allow_duplicate=True),
    Output('s-fit-upload-status', 'children', allow_duplicate=True),
    Output('s-elem', 'value'), Output('s-dline', 'value'),
    Output('s-lcell', 'value', allow_duplicate=True),
    Output('s-B', 'value', allow_duplicate=True),
    *[Output(id_, 'value', allow_duplicate=id_ in ('s-gbuf', 's-shift', 's-btheta', 's-bphi',
                                                   's-doppT', 's-pol', 's-theta0', 's-bufP',
                                                   's-rb85frac'))
      for id_, _ in SAMPLE_PANEL],
    Input('s-fit-sample', 'n_clicks'),
    prevent_initial_call=True,
)
def load_fit_sample(n):
    """Load the sample spectrum and set the panel to the conditions it was taken
    in, so that 'Load sample spectrum' then 'Run Fit' works with no other input
    (T, the quantity the fit finds, is left as it is)."""
    try:
        with open(SAMPLE_CSV, 'rb') as fh:
            contents = 'data:text/csv;base64,' + base64.b64encode(fh.read()).decode()
        X_ghz, Y = _parse_fit_csv(contents)
    except Exception as e:
        return (no_update, f'Could not load the sample spectrum: {e}',
                *(no_update,) * (4 + len(SAMPLE_PANEL)))
    return ({'X': X_ghz.tolist(), 'Y': Y.tolist()},
            f'Sample: Rb D2, 75 mm cell, B = 0, ~40 °C ({len(X_ghz)} points). '
            'The other boxes are set to match: no buffer gas, no shift or '
            'broadening, B along the beam, x-polarised light.',
            'Rb', 'D2', 0.075, 0.0, *[v for _, v in SAMPLE_PANEL])


def _fit_warnings(result, vary_keys, labels, scales):
    """What the user needs to know about a fit that did not do what a fit
    does, and the fitted parameters whose values the data did not determine
.

    lmfit reports success for a fit that never left its start (the model was
    flat there) and for a parameter whose error bar is 10^17 times its size.
    `vary_keys` are the lmfit parameter names that varied, `labels` their names
    on the page, `scales` the error above which a value is taken not to have
    been determined. Returns (lines to show, set of undetermined keys)."""
    lines, undetermined = [], set()
    message = str(getattr(result, 'message', '') or '').strip()
    if not getattr(result, 'success', True):
        lines.append(f'The optimiser did not report success: {message or "no reason given"}')
    params = getattr(result, 'params', None)
    if params is None:
        return lines, undetermined
    vary_keys = [k for k in vary_keys if k in params]
    moved = [k for k in vary_keys
             if not np.isclose(params[k].value, params[k].init_value, rtol=1e-9, atol=1e-12)]
    if vary_keys and not moved:
        lines.append('The fit did not move from the starting values: the model may be '
                     'flat there (a cell too hot to transmit, say), so the result is '
                     'the starting model, not a fit. Try other starting values.'
                     + (f' ({message})' if message else ''))
    for k in vary_keys:
        err = params[k].stderr
        if err is None or not np.isfinite(err):
            if moved:
                lines.append(f'{labels[k]}: no error estimate.')
        elif k in scales and err > scales[k]:
            undetermined.add(k)
            lines.append(f'{labels[k]}: not determined by these data (error estimate '
                         f'{err:.3g}); its box keeps the value you typed.')
    return lines, undetermined


@app.callback(
    Output('s-fit-plot', 'figure'),
    Output('s-fit-results', 'children'),
    Output('error-store', 'data', allow_duplicate=True),
    Output('s-T', 'value'), Output('s-B', 'value'),
    Output('s-gbuf', 'value'), Output('s-shift', 'value'),
    Output('s-lcell', 'value'), Output('s-btheta', 'value'), Output('s-bphi', 'value'),
    Output('s-fit-result-store', 'data'),
    # [ES4] appended last so the indices above stay put
    Output('s-doppT', 'value'), Output('s-rb85frac', 'value'),
    Output('s-pol', 'value'), Output('s-theta0', 'value'),
    Output('s-bufP', 'value'),
    Output('warning-store', 'data', allow_duplicate=True),
    Input('s-fit-run', 'n_clicks'),
    State('s-fit-data-store', 'data'),
    State('s-fit-datatype', 'value'), State('s-fit-algo', 'value'),
    State('s-fit-vary', 'value'),
    State('s-fit-min', 'value'), State('s-fit-max', 'value'),
    State('s-elem','value'),  State('s-dline','value'),
    State('s-T','value'),     State('s-lcell','value'),
    State('s-B','value'),     State('s-gbuf','value'),
    State('s-shift','value'),
    State('s-btheta','value'),State('s-bphi','value'),
    State('s-theta0','value'),State('s-pol','value'),
    State('s-bufgas','value'),State('s-bufP','value'),State('s-bufTf','value'),
    State('s-rb85frac','value'),State('s-K40frac','value'),State('s-K41frac','value'),
    # [ES4] appended last, with defaults below, so positional callers keep working
    State('s-nonortho','value'), State('s-windows','value'), State('s-nwin','value'),
    State('s-doppT','value'),
    State('s-xunit','value'),
    prevent_initial_call=True,
)
def run_fit_single(n, data, data_type, algo, vary_keys, fit_min, fit_max,
                   elem, dline, T, lcell, B, gbuf, shift, btheta, bphi, theta0, pol,
                   bufgas, bufP, bufTf, rb85frac, K40frac, K41frac,
                   nonortho=NONORTHO_DEFAULT, windows=None, n_window=None, doppT=None,
                   xunit='GHz'):
    no_change = (no_update,) * 7
    no_extra = (no_update,) * 5
    if not data:
        return (empty_fig('Upload data first'), '',
                'No experimental data uploaded yet.', *no_change, None, *no_extra, [])
    with _collect_warnings() as caught:
        try:
            result = _run_fit(data, data_type, algo, vary_keys, fit_min, fit_max,
                              elem, dline, T, lcell, B, gbuf, shift, btheta, bphi,
                              theta0, pol, bufgas, bufP, bufTf, rb85frac, K40frac,
                              K41frac, nonortho, windows, n_window, doppT)
            fit_store = result[10]
            result = (_set_xunit(result[0], fit_store['X'], xunit, elem, dline),
                      *result[1:])
        except Exception as e:
            # the result store is emptied, so Export fit CSV cannot write the
            # previous fit under a page that says this one failed
            result = (empty_fig('Fit failed'), '', str(e), *no_change, None, *no_extra)
    return (*result, caught.messages)


def _run_fit(data, data_type, algo, vary_keys, fit_min, fit_max,
             elem, dline, T, lcell, B, gbuf, shift, btheta, bphi, theta0, pol,
             bufgas, bufP, bufTf, rb85frac, K40frac, K41frac,
             nonortho, windows, n_window, doppT):
    """run_fit_single's work: the callback's 16 outputs before the warnings."""
    no_change = (no_update,) * 7
    X_ghz = np.array(data['X'], dtype=float)
    Y     = np.array(data['Y'], dtype=float)
    n_all = len(X_ghz)
    # a blank end of the range is open
    lo = -np.inf if fit_min is None else float(fit_min)
    hi = np.inf if fit_max is None else float(fit_max)
    mask = (X_ghz >= lo) & (X_ghz <= hi)
    X_ghz, Y = X_ghz[mask], Y[mask]
    if len(X_ghz) < 2:
        raise ValueError('Not enough data points in the selected fit range')
    X_mhz = X_ghz * 1e3

    _require_isotopes(elem, rb85frac, K40frac, K41frac)
    p0 = _p(elem, dline, T, lcell, B, gbuf, shift, btheta, bphi, theta0,
           bufgas, bufP, bufTf, rb85frac=rb85frac, K40frac=K40frac, K41frac=K41frac,
           nonortho=nonortho, windows=windows, n_window=n_window, doppT=doppT)
    E_fit, pol_key = _pol_fit_setup(_pol_efield(pol, theta0))

    # A choice the current setup cannot use is dropped even if it is still
    # ticked from before: the checklist only greys such options out.
    disabled = _fit_vary_disabled(bufgas, elem, doppT, pol, dline, bufP)
    chosen = [k for k in (vary_keys or []) if k not in disabled]
    if not chosen:
        dropped = [label for key, label in VARY_PARAM_LABELS if key in (vary_keys or [])]
        raise ValueError(
            ('Nothing to vary: ' + ', '.join(dropped) + ' cannot be varied with these '
             'settings (greyed out). Tick another parameter under "Vary in fit".')
            if dropped else
            'Nothing to vary: tick at least one parameter under "Vary in fit".')
    vary = {key: (key in chosen) for key, _ in VARY_PARAM_LABELS
            if key != 'pol_angle'}
    if 'pol_angle' in chosen:
        vary[pol_key] = True

    # differential_evolution requires finite bounds on every varying parameter,
    # and the GUI deliberately doesn't collect explicit bounds - build a search
    # window around each initial guess instead: ±max(50% of |value|, floor),
    # where the floor keeps the window usable when the guess is near zero.
    scales = _fit_scales(p0)
    p_dict_bounds = None
    if algo != 'ML':
        p_dict_bounds = {}
        for key, on in vary.items():
            if not on:
                continue
            if key in ('E_x', 'E_y'):
                # the other amplitude is the larger and held fixed: +-4x it
                # spans linear angles to ~76 deg either side of that axis
                fixed = abs(E_fit[0]) if key == 'E_y' else E_fit[1][0]
                p_dict_bounds[key] = [-4.0 * fixed, 4.0 * fixed]
                continue
            v0 = float(p0.get(key, 0.0))
            half = max(abs(v0) * 0.5, scales.get(key, 1.0))
            lo_b, hi_b = v0 - half, v0 + half
            # keep the window inside the physical range:
            # GammaBuf from 0 searched [-100, 100] MHz, and the first
            # negative width tried stopped the fit
            if key in FIT_LIMITS:
                lo_b = max(lo_b, FIT_LIMITS[key][0])
                hi_b = min(hi_b, FIT_LIMITS[key][1])
            p_dict_bounds[key] = [lo_b, hi_b]

    best, rms, result = fit_data([X_mhz, Y], p0, vary, E_in=E_fit,
                                 p_dict_bounds=p_dict_bounds,
                                 data_type=data_type, fit_algorithm=algo)

    fig = _fit_fig(X_ghz, Y, result.best_fit, OUTPUT_LABELS.get(data_type, data_type), rms)

    def stderr_of(key):
        param = result.params.get(key) if hasattr(result, 'params') else None
        return getattr(param, 'stderr', None) if param is not None else None

    pol_angle = None
    if vary.get(pol_key) and 'pol_angle' in chosen:
        E_x, E_y = float(best['E_x']), float(best['E_y'])
        pol_angle = _pol_angle_deg(E_x, E_y)

    # the lmfit parameter behind each row, and the row's name on the page
    row_params = {k: (pol_key if k == 'pol_angle' else k) for k in chosen}
    row_labels = {k: dict(VARY_PARAM_LABELS)[k] for k in chosen}
    notes, undetermined_params = _fit_warnings(
        result, list(row_params.values()),
        {row_params[k]: row_labels[k] for k in chosen},
        {row_params[k]: scales[k] for k in chosen if k in scales})
    undetermined = {k for k in chosen if row_params[k] in undetermined_params}

    rows = []
    for key, label in VARY_PARAM_LABELS:
        if key not in chosen:
            continue
        if key == 'pol_angle':
            val = pol_angle
            # d(atan2(E_y, E_x)) = (E_x dE_y - E_y dE_x) / (E_x^2 + E_y^2)
            s = stderr_of(pol_key)
            lever = (E_x if pol_key == 'E_y' else -E_y) / (E_x**2 + E_y**2)
            stderr = abs(np.degrees(lever * s)) if s else None
            label = f'{label} [deg]'
        else:
            val, stderr = best.get(key), stderr_of(key)
            if key in ('Btheta', 'Bphi'):
                # fit runs in radians (core convention); GUI displays degrees
                val = np.degrees(val)
                stderr = np.degrees(stderr) if stderr else stderr
                label = f'{label} [deg]'
        rows.append(html.Tr([
            html.Td(label, style={'paddingRight':'10px'}),
            html.Td(f'{val:.6g}'),
            html.Td(f'± {stderr:.3g}' if stderr else '—',
                    style={'color':SUBTEXT,'paddingLeft':'8px'}),
        ]))
    used = (f'{len(X_ghz)} of {n_all} points fitted' if len(X_ghz) < n_all
            else f'all {n_all} points fitted')
    results_children = html.Div([
        html.P(f'RMS residual: {rms:.6g}',
               style={'fontWeight':'700','color':TEXT,'fontSize':'0.85rem'}),
        *[html.P(f'⚠ {line}', style={'color': WARN_COLOUR, 'fontSize': '0.8rem',
                                     'margin': '2px 0'}) for line in notes],
        html.Table([html.Tbody(rows)], style={'fontSize':'0.8rem','color':TEXT}),
        html.P(f'{used}, {X_ghz.min():.4g} to {X_ghz.max():.4g} GHz. Model: '
               + (_model_tag(p0).lstrip(' ·') or 'non-orthogonal'),
               style={'color':SUBTEXT,'fontSize':'0.75rem','marginTop':'4px'}),
    ])

    fit_store = {'X': X_ghz.tolist(), 'data': Y.tolist(),
                 'fit': list(np.real(np.asarray(result.best_fit))),
                 'type': data_type, 'elem': elem, 'dline': dline}

    def written(key, value):
        """The fitted value for a box, unless the data left it undetermined:
        then the box keeps what was typed."""
        return no_update if key in undetermined else value

    return (fig, results_children, '',
            written('T', best.get('T', T)), written('Bfield', best.get('Bfield', B)),
            written('GammaBuf', best.get('GammaBuf', gbuf)),
            written('shift', best.get('shift', shift)),
            written('lcell', best.get('lcell', lcell)),
            # best-fit angles are radians (core convention); inputs are degrees
            written('Btheta', np.degrees(best['Btheta']) if 'Btheta' in best else btheta),
            written('Bphi', np.degrees(best['Bphi']) if 'Bphi' in best else bphi),
            fit_store,
            written('DoppTemp', best['DoppTemp']) if 'DoppTemp' in chosen else no_update,
            written('rb85frac', best['rb85frac']) if 'rb85frac' in chosen else no_update,
            # a fitted angle is written back as Linear (x) rotated by it
            'Linear (x)' if pol_angle is not None and 'pol_angle' not in undetermined
            else no_update,
            pol_angle if pol_angle is not None and 'pol_angle' not in undetermined
            else no_update,
            written('BufferPressure', best['BufferPressure'])
            if 'BufferPressure' in chosen else no_update)


# ─────────────────────────────────────────────────────────────────────────────
# Cascade compute
# ─────────────────────────────────────────────────────────────────────────────
def _cell2_alone(X, E_in, E_1, S0_1, S0_2, p2, pol_rad, reset):
    """The cascade's lower panel [ES4 NEW]: the light reaching cell 2, as a
    fraction of the input, and cell 2's transmission of it (out / in). Their
    product is the total S0.

    With the polariser (and the QWP), light reaches cell 2 in the same state at
    every detuning -- linear at the polariser angle, or circular -- so cell 2's
    transmission is its spectrum for that input, defined even where no light
    arrives: for a wing selector, the WS curve. The light reaching it is cell
    1's field E_1 through the polariser, |E_1 . p|^2, which the QWP keeps (see
    calculate_cascade). Without the polariser, cell 2 takes whatever cell 1
    sends on: the light reaching it is cell 1's S0_1, and its transmission is
    the core's per-stage S0_2, which is undefined where cell 1 passes no light
    (the core's own test, spectra._lit): NaN there, a gap in the plot."""
    if pol_rad is None:
        lit = spectra._lit(np.real((E_1 * E_1.conjugate()).sum(axis=0)))
        return S0_1, np.where(lit, S0_2, np.nan)
    c, s = np.cos(pol_rad), np.sin(pol_rad)
    into = np.abs(E_1[0] * c + E_1[1] * s) ** 2 / np.real(np.vdot(E_in, E_in))
    if reset is None:
        state = [c, s, 0.]
    else:
        state = [1 / np.sqrt(2), (1j if reset == 'left_cp' else -1j) / np.sqrt(2), 0.]
    [T_2] = calculate(X, state, p2, outputs=['S0'])
    return into, np.asarray(T_2).squeeze().real


# Cell 2 species locked to cell 1 (same-species cascaded filters)
@app.callback(
    Output('c2-elem',  'value'),
    Output('c2-dline', 'value'),
    Input('c1-elem',  'value'),
    Input('c1-dline', 'value'),
)
def sync_cell2_species(elem, dline):
    return elem, dline


@app.callback(
    Output('cas-plot',  'figure'),
    Output('cas-store', 'data'),
    Output('error-store', 'data', allow_duplicate=True),
    Output('warning-store', 'data', allow_duplicate=True),
    Input('cas-compute', 'n_clicks'),
    # Cell 1
    State('c1-elem','value'), State('c1-dline','value'),
    State('c1-T','value'),    State('c1-lcell','value'),
    State('c1-B','value'),    State('c1-gbuf','value'),  State('c1-shift','value'),
    State('c1-btheta','value'),State('c1-bphi','value'), State('c1-theta0','value'),
    State('c1-pol','value'),
    State('c1-bufgas','value'),State('c1-bufP','value'), State('c1-bufTf','value'),
    State('c1-rb85frac','value'),State('c1-K40frac','value'),State('c1-K41frac','value'),
    # Cell 2 (element/D-line/isotopes inherited from cell 1)
    State('c2-T','value'),    State('c2-lcell','value'),
    State('c2-B','value'),    State('c2-gbuf','value'),  State('c2-shift','value'),
    State('c2-btheta','value'),State('c2-bphi','value'),
    State('c2-bufgas','value'),State('c2-bufP','value'), State('c2-bufTf','value'),
    # Polariser + QWP + detuning + outputs
    State('cas-use-pol','value'), State('cas-pol-angle','value'),
    State('cas-reset-pol','value'),
    State('cas-det-start','value'),State('cas-det-stop','value'),State('cas-det-N','value'),
    State('cas-outputs','value'),
    # [ES4] per-cell propagation model and Doppler T; appended last, with
    # defaults below, so positional callers keep working
    State('c1-nonortho','value'), State('c1-windows','value'), State('c1-nwin','value'),
    State('c1-doppT','value'),
    State('c2-nonortho','value'), State('c2-windows','value'), State('c2-nwin','value'),
    State('c2-doppT','value'),
    State('cas-xunit','value'),
    prevent_initial_call=True,
)
def compute_cascade(n,
                    e1, d1, T1, l1, B1, g1, sh1, bt1, bp1, th1, pol1, bg1, bP1, bTf1,
                    rb85, k40, k41,
                    T2, l2, B2, g2, sh2, bt2, bp2, bg2, bP2, bTf2,
                    use_pol, pol_angle, reset_pol,
                    det_start, det_stop, det_N, outputs,
                    no1=None, win1=None, nw1=None, dT1=None,
                    no2=None, win2=None, nw2=None, dT2=None, xunit='GHz'):
    """Cascaded same-species filters via the core calculate_cascade().

    With the intermediate polariser on (and the QWP after it, if chosen),
    cell-2 outputs from the core are instrument totals normalised to the
    laser input. The QWP needs the polariser. With neither,
    totals are formed by the documented combination rules (intensities x
    S0_cell1; fractional Stokes pass through). The QWP option models a
    wing-selector topology (e.g. VAMOS, Oliviero 1998): it prepares a FIXED
    circular state for cell 2 regardless of wavelength, so cell 2 always
    absorbs the same Zeeman component.

    [ES4 NEW] A lower panel shows cell 2 on its own (_cell2_alone): the light
    reaching it and its transmission of that light, whose product is the
    total S0."""
    with _collect_warnings() as caught:
        try:
            from elecsus.elecsus_methods import calculate_cascade

            # Both cells share element/D-line/isotope fractions (same vapour species)
            _require_isotopes(e1, rb85, k40, k41)
            p1 = _p(e1, d1, T1, l1, B1, g1, sh1, bt1, bp1, th1, bg1, bP1, bTf1,
                    rb85frac=rb85, K40frac=k40, K41frac=k41,
                    nonortho=no1, windows=win1, n_window=nw1, doppT=dT1,
                    cell='Cell 1')
            p2 = _p(e1, d1, T2, l2, B2, g2, sh2, bt2, bp2, 0., bg2, bP2, bTf2,
                    rb85frac=rb85, K40frac=k40, K41frac=k41,
                    nonortho=no2, windows=win2, n_window=nw2, doppT=dT2,
                    cell='Cell 2')

            E_in = np.array(_pol_efield(pol1, th1), dtype='complex')
            X    = _X(det_start, det_stop, det_N)
            outs = [k for k in (outputs or ['S0']) if k in CASCADE_OUTPUTS] or ['S0']

            pol_on = 'on' in (use_pol or [])
            pol_deg = _num(pol_angle, 'Polariser angle [deg]') if pol_on else None
            pol_rad = np.radians(pol_deg) if pol_on else None
            reset_val = reset_pol if reset_pol in ('left_cp', 'right_cp') else None
            # [ES4 FIX] the core refuses a QWP without the polariser
            if reset_val is not None and not pol_on:
                raise ValueError('The QWP needs the polariser between the cells in '
                                 'front of it: tick "Polariser between cells", or set '
                                 'the QWP to None.')

            # The lower panel needs cell 1's field and, without the polariser,
            # cell 2's per-stage S0: asked for here, plotted only if ticked
            asked = outs + [k for k in ('S0', 'E_out') if k not in outs]
            res1, res2 = calculate_cascade(X, E_in, p1, p2, outputs=asked,
                                           polariser_angle=pol_rad,
                                           reset_polarisation=reset_val)
            got1 = dict(zip(asked, res1))
            r1 = {k: np.asarray(got1[k]).squeeze().real for k in outs}
            S0_1 = np.asarray(got1['S0']).squeeze().real
            r2 = {k: np.asarray(v).squeeze().real for k, v in zip(asked, res2)
                  if k != 'E_out'}

            if pol_on:
                # core returns instrument totals directly with an intermediate optic
                total = {k: r2[k] for k in outs}
            else:
                # per-stage outputs: combine by the documented rules
                total = {k: (r2[k] if k in ('S1', 'S2', 'S3') else S0_1 * r2[k])
                         for k in outs}
            into, T_2 = _cell2_alone(X, E_in, np.asarray(got1['E_out']), S0_1,
                                     r2['S0'], p2, pol_rad, reset_val)

            traces, keys = [], []
            for k in outs:
                traces.append((f'After Cell 1 — {OUTPUT_LABELS.get(k,k)}', r1[k], True,  0.7))
                traces.append((f'Total after both filters — {OUTPUT_LABELS.get(k,k)}',
                               total[k], False, 1.0))
                keys += [k, k]
            reaches = (' (after the polariser and QWP)' if reset_val is not None else
                       ' (after the polariser)' if pol_on else ' (cell 1 output)')
            state = ({'left_cp': 'left-circular', 'right_cp': 'right-circular'}
                     .get(reset_val, f'linear at {pol_deg:g}°') if pol_on
                     else 'what cell 1 sends')
            traces += [(f'Light reaching cell 2{reaches}', into, True, 1.0),
                       (f'Cell 2 alone: out ÷ in, for {state}', T_2, False, 1.0)]
            keys += ['cell2_into', 'cell2_T']

            fig = _set_xunit(_fig(X/1e3, traces, keys=keys), X/1e3, xunit, e1, d1)
            pol_txt = f' · polariser at {pol_deg:g}°' if pol_on else ' · no polariser'
            qwp_txt = f' · QWP → {reset_val}' if reset_val is not None else ''
            fig.update_layout(title=dict(
                text=f'Cascaded filters: {e1} {d1} → {e1} {d1}{pol_txt}{qwp_txt}'
                     + ''.join(f' · cell {i} ({tag[3:]})'
                               for i, tag in ((1, _model_tag(p1)), (2, _model_tag(p2)))
                               if tag),
                font=dict(color=TEXT, size=13)))

            store = {'X': list(X/1e3),
                     'cell1': {k: list(v) for k, v in r1.items()},
                     'total': {k: list(v) for k, v in total.items()},
                     # NaN (no light reaching cell 2) as None, which JSON has
                     'cell2': {'into': list(into),
                               'transmission': [None if np.isnan(v) else v for v in T_2]},
                     'elem': e1, 'dline': d1}
            result = (fig, store, '')
        except Exception as e:
            result = (empty_fig('Calculation failed'), None, str(e))
    return (*result, caught.messages)


# ─────────────────────────────────────────────────────────────────────────────
# Error and warning toasts (one callback each, watching the shared stores)
# ─────────────────────────────────────────────────────────────────────────────
@app.callback(
    Output('error-toast', 'children'),
    Output('error-toast', 'is_open'),
    Input('error-store', 'data'),
)
def show_error(msg):
    if msg:
        return msg, True
    return '', False


@app.callback(
    Output('warning-toast', 'children'),
    Output('warning-toast', 'is_open'),
    Input('warning-store', 'data'),
)
def show_warnings(messages):
    """The core's warnings from the last Compute, Run Fit or Compute Cascade."""
    if messages:
        return [html.Div(m, style={'marginBottom': '4px'}) for m in messages], True
    return '', False


# ─────────────────────────────────────────────────────────────────────────────
# CSV exports
# ─────────────────────────────────────────────────────────────────────────────
def _csv_value(v):
    """A data value for the CSV, to 10 significant digits whatever its size: a
    fixed 8 decimals wrote every transmission below 5e-9 as 0."""
    return f'{v:.10g}'


@app.callback(Output('s-download','data'), Input('s-export','n_clicks'),
              State('s-store','data'), prevent_initial_call=True)
def export_single(n, store):
    if not store: return no_update
    X = np.array(store['X'])
    keys = list(store['results'].keys())
    lines = ['Detuning_GHz,' + ','.join(keys)]
    for i,x in enumerate(X):
        vals = [_csv_value(np.array(store["results"][k])[i]) for k in keys]
        lines.append(f'{x:.6f},' + ','.join(vals))
    return dict(content='\n'.join(lines), filename='elecsus_single.csv')


@app.callback(Output('s-fit-download','data'), Input('s-fit-export','n_clicks'),
              State('s-fit-result-store','data'), prevent_initial_call=True)
def export_fit(n, store):
    if not store: return no_update
    X = np.array(store['X']); data = np.array(store['data']); fit = np.array(store['fit'])
    label = store.get('type', 'S0')
    lines = [f'Detuning_GHz,{label}_data,{label}_fit,residual']
    for i, x in enumerate(X):
        lines.append(f'{x:.6f},' + ','.join(
            _csv_value(v) for v in (data[i], fit[i], data[i] - fit[i])))
    return dict(content='\n'.join(lines), filename='elecsus_fit.csv')


@app.callback(Output('cas-download','data'), Input('cas-export','n_clicks'),
              State('cas-store','data'), prevent_initial_call=True)
def export_cascade(n, store):
    if not store: return no_update
    X = np.array(store['X'])
    keys = list(store['total'].keys())
    hdr = 'Detuning_GHz,' + ','.join([f'Cell1_{k}' for k in keys]) + ',' + \
          ','.join([f'Total_{k}' for k in keys])
    # [ES4 NEW] the lower panel, cell 2 on its own, in the last two columns
    alone = store.get('cell2')
    if alone:
        hdr += ',Into_Cell2,Cell2_alone_S0'
    lines = [hdr]
    for i,x in enumerate(X):
        c1 = [_csv_value(np.array(store["cell1"].get(k,[0]*len(X)))[i]) for k in keys]
        c2 = [_csv_value(np.array(store["total"][k])[i]) for k in keys]
        if alone:
            c2 += [_csv_value(alone['into'][i]),
                   _csv_value(np.nan if alone['transmission'][i] is None
                              else alone['transmission'][i])]
        lines.append(f'{x:.6f},' + ','.join(c1) + ',' + ','.join(c2))
    return dict(content='\n'.join(lines), filename='elecsus_cascade.csv')


# ─────────────────────────────────────────────────────────────────────────────
# Wavelength axes  [ES4 NEW] -- the "Plot against" selector redraws the last
# result without computing again, and each Detuning Axis panel shows the
# conversion between GHz and mÅ for its line
# ─────────────────────────────────────────────────────────────────────────────
@app.callback(
    Output('s-plot', 'figure', allow_duplicate=True),
    Output('s-fit-plot', 'figure', allow_duplicate=True),
    Input('s-xunit', 'value'),
    State('s-plot', 'figure'), State('s-store', 'data'),
    State('s-fit-plot', 'figure'), State('s-fit-result-store', 'data'),
    prevent_initial_call=True,
)
def redraw_single(unit, fig, store, fit_fig, fit_store):
    return _redraw(fig, store, unit), _redraw(fit_fig, fit_store, unit)


@app.callback(
    Output('cas-plot', 'figure', allow_duplicate=True),
    Input('cas-xunit', 'value'),
    State('cas-plot', 'figure'), State('cas-store', 'data'),
    prevent_initial_call=True,
)
def redraw_cascade(unit, fig, store):
    return _redraw(fig, store, unit)


@app.callback(
    Output('s-det-conv', 'children'),
    Input('s-elem', 'value'), Input('s-dline', 'value'),
    Input('s-det-start', 'value'), Input('s-det-stop', 'value'),
)
def conversion_single(elem, dline, start, stop):
    return _conversion_text(elem, dline, start, stop)


@app.callback(
    Output('cas-det-conv', 'children'),
    Input('c1-elem', 'value'), Input('c1-dline', 'value'),
    Input('cas-det-start', 'value'), Input('cas-det-stop', 'value'),
)
def conversion_cascade(elem, dline, start, stop):
    return _conversion_text(elem, dline, start, stop)


# ─────────────────────────────────────────────────────────────────────────────
# GammaBuf greying — disabled + styled while the buffer-gas model supplies the
# broadening, which replaces GammaBuf. That is when a gas with data for the
# element and line is selected and P > 0; with P = 0, or a gas the model has no
# data for (Na), the core uses the static GammaBuf, so the box stays open
#. (The shift input always stays: the model's pressure shift is
# added to it.)
# ─────────────────────────────────────────────────────────────────────────────
def _gbuf_state(bufgas, elem, dline, bufP):
    active = _buffer_model_on(bufgas, elem, dline, bufP)
    style  = INPUT_DISABLED_STYLE if active else INPUT_STYLE
    return active, style

@app.callback(
    Output('s-gbuf', 'disabled'),
    Output('s-gbuf', 'style'),
    Input('s-bufgas', 'value'),
    Input('s-elem', 'value'), Input('s-dline', 'value'), Input('s-bufP', 'value'),
)
def grey_gbuf_single(bufgas, elem, dline, bufP):
    return _gbuf_state(bufgas, elem, dline, bufP)

@app.callback(
    Output('c1-gbuf', 'disabled'),
    Output('c1-gbuf', 'style'),
    Input('c1-bufgas', 'value'),
    Input('c1-elem', 'value'), Input('c1-dline', 'value'), Input('c1-bufP', 'value'),
)
def grey_gbuf_c1(bufgas, elem, dline, bufP):
    return _gbuf_state(bufgas, elem, dline, bufP)

@app.callback(
    Output('c2-gbuf', 'disabled'),
    Output('c2-gbuf', 'style'),
    Input('c2-bufgas', 'value'),
    Input('c2-elem', 'value'), Input('c2-dline', 'value'), Input('c2-bufP', 'value'),
)
def grey_gbuf_c2(bufgas, elem, dline, bufP):
    return _gbuf_state(bufgas, elem, dline, bufP)


# ─────────────────────────────────────────────────────────────────────────────
# Propagation controls — the window model needs the non-orthogonal solver
# (the orthogonal path ignores n_window), and the index input needs the window
# model, so each is enabled only when what it depends on is ticked (the
# state itself is _prop_control_state, above prop_section, which also sets
# the layout's initial state).
# ─────────────────────────────────────────────────────────────────────────────
def _register_prop_controls(prefix):
    @app.callback(
        Output(f'{prefix}-windows', 'options'),
        Output(f'{prefix}-nwin', 'disabled'),
        Output(f'{prefix}-nwin', 'style'),
        Input(f'{prefix}-nonortho', 'value'),
        Input(f'{prefix}-windows', 'value'),
    )
    def prop_controls(nonortho, windows):
        return _prop_control_state(nonortho, windows)
    return prop_controls


prop_controls_s, prop_controls_c1, prop_controls_c2 = (
    _register_prop_controls(p) for p in ('s', 'c1', 'c2'))


@app.callback(
    Output('s-outputs', 'options'),
    Output('s-outputs', 'value'),
    Input('s-nonortho', 'value'),
    State('s-outputs', 'value'),
)
def nonortho_output_options(nonortho, selected):
    """Grey out S0_P while the non-orthogonal solver is off, and untick it, so
    a compute never asks the core for an output it would refuse."""
    on = _on(nonortho)
    kept = [o for o in (selected or [])
            if on or o not in NONORTHO_ONLY_OUTPUTS] or ['S0']
    return _output_options(OUTPUTS, nonortho_on=on), kept


# ─────────────────────────────────────────────────────────────────────────────
# Isotope row visibility — Rb shows rb85frac; K shows K40frac + K41frac;
# Cs, Na have no mixed-isotope parameters so both rows stay hidden.
# ─────────────────────────────────────────────────────────────────────────────
def _iso_styles(elem):
    show = {'display': 'block'}
    hide = {'display': 'none'}
    return (show if elem == 'Rb' else hide,
            show if elem == 'K'  else hide)

@app.callback(
    Output('s-iso-rb-row', 'style'),
    Output('s-iso-k-row',  'style'),
    Input('s-elem', 'value'),
)
def iso_visibility_s(elem):
    return _iso_styles(elem)

@app.callback(
    Output('c1-iso-rb-row', 'style'),
    Output('c1-iso-k-row',  'style'),
    Input('c1-elem', 'value'),
)
def iso_visibility_c1(elem):
    return _iso_styles(elem)


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────
def run(debug=False, port=8050, host='127.0.0.1'):
    print(f'\n  ElecSus Dash GUI  —  http://{host}:{port}\n')
    app.run(debug=debug, port=port, host=host)


# ─────────────────────────────────────────────────────────────────────────────
# Energy levels  [ES4 NEW]
# The levels of one isotope in a field, the lines between them and the
# spectrum, drawn with Plotly from elecsus.energy_levels (which draws the same
# figure with matplotlib for a paper or a thesis).
# ─────────────────────────────────────────────────────────────────────────────
KIND_COLOUR = {'sigma+': ACCENT, 'sigma-': HIGHLIGHT, 'pi': WARN_COLOUR}
_LEVELS_MARGIN = dict(l=70, t=100, b=50)


def _rgba(hex_colour, alpha):
    h = hex_colour.lstrip('#')
    return f'rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{alpha:.2f})'


def _axis_refs(fig, r, c):
    """('x2', 'y3', 'xaxis2', 'yaxis3') for the subplot at row r, column c."""
    sub = fig.get_subplot(r, c)
    xn, yn = sub.xaxis.plotly_name, sub.yaxis.plotly_name
    return xn.replace('axis', ''), yn.replace('axis', ''), xn, yn


def _levels_ranges(curves, levels):
    out = {}
    for which in ('ground', 'excited'):
        e = np.concatenate([curves[which].ravel(), [l['energy'] for l in levels[which]]])
        pad = 0.08 * ((e.max() - e.min()) or 1.)
        out[which] = (float(e.min() - pad), float(e.max() + pad))
    return out


def _levels_fig(iso, dline, B, T, lcell, geometry, labels, min_strength, natural, det_range):
    """The figure, and the lines and levels it shows."""
    kinds = EL.default_kinds(geometry)
    levels = EL.energy_levels(iso, dline, B)
    lines = EL.transitions(iso, dline, B, kinds=kinds, min_strength=min_strength)
    if not lines:
        raise ValueError('No line is that strong: lower Min. strength (1 keeps only the '
                         'strongest).')
    curves = EL.breit_rabi(iso, dline, B, n=81)
    lo, hi = det_range if det_range is not None else EL.detuning_range(lines)
    all_lines, lines = lines, EL.in_range(lines, lo, hi)   # the figure shows these
    sticks = EL.in_range(EL.transitions(iso, dline, B, kinds=kinds), lo, hi)
    X = np.linspace(lo, hi, 1500)
    [S0] = EL.spectrum(iso, dline, B, X, T=T, lcell=lcell, geometry=geometry, natural=natural)

    # rows: the spectrum, the excited and ground levels, the stick spectrum
    fig = make_subplots(rows=4, cols=2, column_widths=[0.16, 0.84],
                        row_heights=[0.21, 0.29, 0.29, 0.21], shared_xaxes=True,
                        shared_yaxes=True, horizontal_spacing=0.0, vertical_spacing=0.045,
                        specs=[[None, {}], [{}, {}], [{}, {}], [None, {}]])
    xs, ys, _, ysn = _axis_refs(fig, 1, 2)
    row_of = {'excited': 2, 'ground': 3}
    ranges = _levels_ranges(curves, levels)

    # the spectrum, with a mark at each line
    fig.add_trace(go.Scatter(x=X, y=S0, mode='lines', line=dict(color=TEXT, width=1.4),
                             name='Transmission', showlegend=False,
                             hovertemplate='%{x:.3f} GHz: %{y:.4f}<extra></extra>'), 1, 2)
    for line in sorted(lines, key=lambda l: l['relative']):
        fig.add_vline(x=line['detuning'], row=1, col=2, layer='below',
                      line=dict(color=_rgba(KIND_COLOUR[line['kind']],
                                            0.3 + 0.7 * line['relative']), width=1.2))
    for kind, x in (('sigma-', 0.04), ('sigma+', 0.96), ('pi', 0.5)):
        if kind in kinds:
            fig.add_annotation(text=EL._SYMBOL[kind], xref=f'{xs} domain', yref=f'{ys} domain',
                               x=x, y=0.12 if kind != 'pi' else 0.92, showarrow=False,
                               font=dict(size=18, color=KIND_COLOUR[kind]))

    # the levels against field, and across the detuning axis at the field
    for which, r in row_of.items():
        for k, level in enumerate(levels[which]):
            fig.add_trace(go.Scatter(x=curves['B'], y=curves[which][:, k], mode='lines',
                                     line=dict(color=TEXT, width=1), showlegend=False,
                                     hoverinfo='skip'), r, 1)
            text = f"{level['label']}<br>|F, mF⟩ = {EL._ket(level['F'], level['mF'])} at B = 0" \
                   f"<br>{level['energy']:+.4f} GHz"
            fig.add_trace(go.Scatter(x=[lo, hi], y=[level['energy']] * 2, mode='lines',
                                     line=dict(color=TEXT, width=1.2), showlegend=False,
                                     hovertemplate=text + '<extra></extra>'), r, 2)
        fig.update_yaxes(range=ranges[which], row=r, col=1)
        fig.update_yaxes(range=ranges[which], row=r, col=2)

    # an arrow for each line: the shaft a shape across the two panels (in paper
    # coordinates, from the axis ranges and domains) ending on the excited
    # level, and the head an annotation whose tip is on that level (a marker
    # is centred on its point, so a triangle's tip would sit above the level).
    # An invisible marker there carries the line's details on hover.
    def paper(which, y):
        dom = fig.layout[_axis_refs(fig, row_of[which], 2)[3]].domain
        r0, r1 = ranges[which]
        return dom[0] + (y - r0) / (r1 - r0) * (dom[1] - dom[0])

    xl, yl, _, _ = _axis_refs(fig, 2, 2)
    hover = dict(x=[], y=[], text=[])
    for line in sorted(lines, key=lambda l: l['relative']):
        g = levels['ground'][line['ground']]
        e = levels['excited'][line['excited']]
        a = 0.3 + 0.7 * line['relative']
        colour = _rgba(KIND_COLOUR[line['kind']], a)
        width = 1 + 3 * line['relative']
        fig.add_shape(type='line', xref=xl, yref='paper', x0=line['detuning'],
                      x1=line['detuning'], y0=paper('ground', g['energy']),
                      y1=paper('excited', e['energy']), line=dict(color=colour, width=width))
        fig.add_annotation(x=line['detuning'], y=e['energy'], xref=xl, yref=yl, text='',
                           showarrow=True, ax=0, ay=14, arrowhead=2, arrowsize=1.2,
                           arrowwidth=width, arrowcolor=colour, standoff=0, startstandoff=0)
        hover['x'].append(line['detuning'])
        hover['y'].append(e['energy'])
        hover['text'].append(f"{EL._SYMBOL[line['kind']]}  {line['detuning']:+.3f} GHz, "
                             f"strength {line['relative']:.3f}<br>{g['label']}<br>→ {e['label']}")
    fig.add_trace(go.Scatter(x=hover['x'], y=hover['y'], mode='markers', name='lines',
                             marker=dict(size=12, color='rgba(0,0,0,0)'),
                             text=hover['text'], hovertemplate='%{text}<extra></extra>',
                             showlegend=False), 2, 2)

    # the stick spectrum: every line in the range, however weak, log scale
    floor = EL.stick_floor(sticks)
    for kind in kinds:
        mine = sorted((l for l in sticks if l['kind'] == kind), key=lambda l: l['relative'])
        if not mine:
            continue
        x = [v for l in mine for v in (l['detuning'], l['detuning'], None)]
        y = [v for l in mine for v in (floor, l['relative'], None)]
        fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name=EL._SYMBOL[kind],
                                 line=dict(color=KIND_COLOUR[kind], width=2),
                                 hoverinfo='skip', legendgroup=kind, showlegend=True), 4, 2)
        fig.add_trace(go.Scatter(
            x=[l['detuning'] for l in mine], y=[l['relative'] for l in mine], mode='markers',
            name=f'{EL._SYMBOL[kind]} strength', legendgroup=kind, showlegend=False,
            marker=dict(size=6, color=KIND_COLOUR[kind]),
            text=[f"{EL._SYMBOL[kind]}  {l['detuning']:+.3f} GHz, strength {l['relative']:.2e}"
                  f"<br>{levels['ground'][l['ground']]['label']}"
                  f"<br>→ {levels['excited'][l['excited']]['label']}" for l in mine],
            hovertemplate='%{text}<extra></extra>'), 4, 2)
    fig.update_yaxes(type='log', range=[np.log10(floor), np.log10(2.)], dtick=1,
                     exponentformat='power', title_text='Relative strength', row=4, col=2)

    # labels to the right of the levels, spread so that they do not overlap
    plot_h = LEVELS_HEIGHT - _LEVELS_MARGIN['t'] - _LEVELS_MARGIN['b']
    longest = 0
    for which, r in row_of.items():
        _, _, _, yname = _axis_refs(fig, r, 2)
        dom = fig.layout[yname].domain
        r0, r1 = ranges[which]
        px_per_unit = (dom[1] - dom[0]) * plot_h / (r1 - r0)
        ordered = sorted(levels[which], key=lambda l: l['energy'])
        gap_px = min(15., 0.98 * (r1 - r0) * px_per_unit / len(ordered))
        size = max(7., min(12., 0.85 * gap_px))
        gap = gap_px / px_per_unit
        ys_ = EL._spread([l['energy'] for l in ordered], gap, r0 + gap / 2, r1 - gap / 2)
        yref = yname.replace('axis', '')
        for level, y in zip(ordered, ys_):
            text = EL._label(level, labels)
            longest = max(longest, len(text) * size)
            fig.add_shape(type='line', xref='paper', yref=yref, x0=1.0, x1=1.012,
                          y0=level['energy'], y1=y, line=dict(color='#b0a8b0', width=0.7))
            fig.add_annotation(text=text, xref='paper', yref=yref, x=1.014, y=y,
                               xanchor='left', showarrow=False, font=dict(size=size))
    head = '|mJ, mI⟩' if labels == 'mJmI' else '|F, mF⟩ at B = 0'
    fig.add_annotation(text=head, xref='paper', yref=f'{ys} domain', x=1.014, y=0.0,
                       xanchor='left', yanchor='bottom', showarrow=False, font=dict(size=13))

    n = {'Na': 3, 'K39': 4, 'K40': 4, 'K41': 4, 'Rb85': 5, 'Rb87': 5, 'Cs': 6}[iso]
    ground, excited = EL._STATE[dline]
    sub = lambda term: f'{n}{term[0]}<sub>{term[1:]}</sub>'
    fig.update_yaxes(title_text=f'{sub(excited)} energy/h (GHz)', row=2, col=1)
    fig.update_yaxes(title_text=f'{sub(ground)} energy/h (GHz)', row=3, col=1)
    fig.update_yaxes(title_text='Transmission', row=1, col=2)
    fig.update_yaxes(showticklabels=False, ticks='', row=2, col=2)
    fig.update_yaxes(showticklabels=False, ticks='', row=3, col=2)
    # every panel of the column: the axes are matched, and a range set on one
    # of them only was overridden by the others' autorange
    for r in (1, 2, 3, 4):
        fig.update_xaxes(range=[lo, hi], row=r, col=2)
    fig.update_xaxes(side='top', title_text='Detuning (GHz)', showticklabels=True, row=1, col=2)
    for r in (2, 3):
        fig.update_xaxes(showticklabels=False, ticks='', showline=False, row=r, col=2)
    fig.update_xaxes(showticklabels=True, title_text='Detuning (GHz)', row=4, col=2)
    fig.update_xaxes(zeroline=False)
    fig.update_yaxes(zeroline=False)
    fig.update_xaxes(range=[0, EL.field_axis_end(B)], row=2, col=1)
    fig.update_xaxes(range=[0, EL.field_axis_end(B)], title_text='B (G)', row=3, col=1)
    mass = ''.join(ch for ch in iso if ch.isdigit())
    name = f"<sup>{mass}</sup>{''.join(ch for ch in iso if ch.isalpha())}" if mass else iso
    cell = 'natural abundance' if natural else 'isotopically pure'
    # the legend names the sticks' colours, beside the stick spectrum
    legend = dict(x=1.01, xanchor='left', yanchor='top', borderwidth=0, font=dict(size=15),
                  y=fig.layout[_axis_refs(fig, 4, 2)[3]].domain[1])
    fig.update_layout(**{**PLOT_LAYOUT, 'hovermode': 'closest', 'legend': legend,
                         'margin': dict(**_LEVELS_MARGIN, r=int(min(460, 40 + 0.55 * longest)))},
                      height=LEVELS_HEIGHT, showlegend=True,
                      title_y=0.985, title_yanchor='top',
                      title=dict(text=f'{name} {dline} · B = {B:g} G ({B / 1e4:.2f} T) · '
                                      f'spectrum: T = {T:g}°C, l = {lcell:g} m, {cell}',
                                 font=dict(color=TEXT, size=13)))
    return fig, all_lines, levels


def _levels_rows(lines, levels):
    """One dict per line, by detuning, for the table and the CSV."""
    rows = []
    for line in sorted(lines, key=lambda l: l['detuning']):
        g = levels['ground'][line['ground']]
        e = levels['excited'][line['excited']]
        rows.append({'kind': EL._SYMBOL[line['kind']], 'detuning': line['detuning'],
                     'relative': line['relative'], 'strength': line['strength'],
                     'ground': g['label'], 'ground_FmF': EL._ket(g['F'], g['mF']),
                     'excited': e['label'], 'excited_FmF': EL._ket(e['F'], e['mF'])})
    return rows


def _levels_table(rows, labels, outside=0):
    cell = {'padding': '2px 10px', 'fontSize': '0.8rem', 'borderBottom': f'1px solid {BORDER}',
            'textAlign': 'left', 'whiteSpace': 'nowrap'}
    g, e = ('ground', 'excited') if labels == 'mJmI' else ('ground_FmF', 'excited_FmF')
    head = ['', 'Detuning [GHz]', 'Strength', 'From (ground)', 'To (excited)']
    body = [html.Tr([html.Td(r['kind'], style=cell), html.Td(f"{r['detuning']:+.4f}", style=cell),
                     html.Td(f"{r['relative']:.3f}", style=cell), html.Td(r[g], style=cell),
                     html.Td(r[e], style=cell)]) for r in rows]
    title = f'{len(rows)} lines' + (f' ({outside} outside the detuning range, so not drawn)'
                                    if outside else '')
    return [html.H4(title, style={'color': HIGHLIGHT, 'fontSize': '0.95rem', 'margin': '8px 0'}),
            html.Table([html.Thead(html.Tr([html.Th(h, style={**cell, 'color': HIGHLIGHT})
                                            for h in head])), html.Tbody(body)],
                       style={'borderCollapse': 'collapse', 'color': TEXT})]


@app.callback(
    Output('lv-plot', 'figure'),
    Output('lv-table', 'children'),
    Output('lv-store', 'data'),
    Output('error-store', 'data', allow_duplicate=True),
    Input('lv-compute', 'n_clicks'),
    State('lv-iso', 'value'), State('lv-dline', 'value'), State('lv-B', 'value'),
    State('lv-geometry', 'value'), State('lv-min', 'value'), State('lv-labels', 'value'),
    State('lv-T', 'value'), State('lv-lcell', 'value'),
    State('lv-det-start', 'value'), State('lv-det-stop', 'value'),
    State('lv-natural', 'value'),
    prevent_initial_call=True,
)
def compute_levels(n, iso, dline, B, geometry, min_strength, labels, T, lcell,
                   det_start, det_stop, natural):
    try:
        B = _num(B, 'B')
        if B < 0:
            raise ValueError('B must be 0 G or more: the levels do not depend on its sign.')
        min_strength = _num(min_strength, 'Min. strength')
        if not 0 <= min_strength <= 1:
            raise ValueError('Min. strength is a fraction of the strongest line: 0 to 1.')
        if (det_start is None) != (det_stop is None):
            raise ValueError('Give both ends of the detuning range, or leave both blank.')
        det_range = None
        if det_start is not None:
            det_range = (_num(det_start, 'Start'), _num(det_stop, 'Stop'))
            if det_range[0] >= det_range[1]:
                raise ValueError('Start must be below Stop.')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fig, lines, levels = _levels_fig(iso, dline, B, _num(T, 'T'), _num(lcell, 'l'),
                                             geometry, labels, min_strength, _on(natural),
                                             det_range)
        rows = _levels_rows(lines, levels)
        outside = 0 if det_range is None else \
            sum(not det_range[0] <= r['detuning'] <= det_range[1] for r in rows)
        store = {'iso': iso, 'dline': dline, 'B': B, 'rows': rows}
        return fig, _levels_table(rows, labels, outside), store, ''
    except Exception as e:
        return empty_fig('Calculation failed'), [], None, str(e)


@app.callback(Output('lv-download', 'data'), Input('lv-export', 'n_clicks'),
              State('lv-store', 'data'), prevent_initial_call=True)
def export_levels(n, store):
    if not store:
        return no_update
    out = ['Kind,Detuning_GHz,Relative_strength,Strength,Ground_mJmI,Ground_FmF_at_B0,'
           'Excited_mJmI,Excited_FmF_at_B0']
    for r in store['rows']:
        out.append(','.join([r['kind'], f"{r['detuning']:.6f}", f"{r['relative']:.6f}",
                             _csv_value(r['strength'])]
                            + [f'"{r[k]}"' for k in ('ground', 'ground_FmF', 'excited',
                                                     'excited_FmF')]))
    name = f"elecsus_levels_{store['iso']}_{store['dline']}_{store['B']:g}G.csv"
    return dict(content='\n'.join(out), filename=name)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='ElecSus Dash GUI')
    parser.add_argument('--port',  type=int, default=8050, help='TCP port (default 8050)')
    parser.add_argument('--host',  default='127.0.0.1',   help='Bind address (default 127.0.0.1)')
    parser.add_argument('--debug', action='store_true',   help='Enable Dash debug/reloader mode')
    args = parser.parse_args()
    run(debug=args.debug, port=args.port, host=args.host)
