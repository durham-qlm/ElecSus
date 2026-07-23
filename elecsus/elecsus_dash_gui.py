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
  - Single Cell tab: forward theory calculation with all parameters
  - Cascaded Two Cells tab: sequential propagation through two cells
  - Buffer gas parameters exposed (BufferGas, BufferPressure, BufferTemp_fill)
  - Interactive Plotly plots, Durham University light theme
  - CSV export

Usage:
    python -m elecsus.elecsus_dash_gui
    # Then open http://127.0.0.1:8050 in a browser

Requirements:
    pip install dash dash-bootstrap-components plotly numpy
"""

import base64
import io

import numpy as np
import dash
from dash import dcc, html, Input, Output, State, no_update, MATCH, ALL
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from elecsus.elecsus_methods import calculate, fit_data

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
ELEMENTS   = ['Rb', 'Cs', 'K', 'Na']
DLINES     = ['D1', 'D2']
# Output order: Stokes/intensities, then the circular-basis quantities paired
# sigma+ then sigma-. GUI labels show the PHYSICAL sigma; the API keys keep the
# historical (inverted) circular-basis names -- e.g. the key 'TransRight' is the
# sigma+ transmission, 'nMinus' the sigma+ index (see spectra.output_list()).
OUTPUTS    = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45',
              'TransRight', 'TransLeft', 'alphaPlus', 'alphaMinus', 'alphaZ',
              'nMinus', 'nPlus', 'phi', 'psi', 'DoLP']
OUTPUT_LABELS = {
    'S0': 'S₀ (Transmission)', 'S1': 'S₁', 'S2': 'S₂', 'S3': 'S₃',
    'Ix': 'Iₓ', 'Iy': 'Iᵧ', 'I_P45': 'I₊₄₅', 'I_M45': 'I₋₄₅',
    # Labelled by physical sigma. Key 'TransRight' (=exp(-2 Im(nMinus) k L)) is the
    # sigma+ transmission; 'TransLeft' (nPlus) the sigma-. alphaPlus/alphaMinus
    # already carry the physical sigma (they cross over internally).
    'TransRight': 'Trans σ⁺', 'TransLeft': 'Trans σ⁻',
    'alphaPlus': 'α σ⁺', 'alphaMinus': 'α σ⁻', 'alphaZ': 'α π',
    # n(sigma±): physically meaningful only for Btheta=0 (Faraday); computed
    # unconditionally. Key 'nMinus' is the sigma+ index, 'nPlus' the sigma-.
    'nMinus': 'n(σ⁺) — Faraday only', 'nPlus': 'n(σ⁻) — Faraday only',
    'phi': 'Faraday rotation φ — dispersion only',
    # psi/DoLP: Stokes-based (dichroism-aware) generalisation of phi.
    # psi = 0.5*atan2(S2,S1); only meaningful where DoLP is appreciable.
    'psi': 'Rotation ψ (Stokes, incl. dichroism)',
    'DoLP': 'Degree of linear pol.',
}
BUFFER_GAS_OPTIONS = ['None', 'Ne', 'Ar', 'Kr', 'Xe', 'N2', 'He']
POL_OPTIONS = ['Linear (x)', 'Linear (y)', 'Left CP', 'Right CP', 'Linear 45°']

# [ES4 NEW] Fitting panel — data-type dropdown reuses OUTPUTS above.
# RR/SA fitting algorithms use multiprocessing.Pool internally, which is fragile
# when invoked from inside a running Dash server process, so the GUI only exposes
# the two algorithms that run in-process (ML directly, DE via scipy through lmfit).
FIT_ALGO_OPTIONS = [
    {'label': 'Marquardt-Levenberg', 'value': 'ML'},
    {'label': 'Differential Evolution', 'value': 'differential_evolution'},
]
# (p_dict key, display label) for parameters the GUI allows the user to vary.
# BufferGas/BufferPressure/BufferTemp_fill are deliberately absent - they are
# fixed cell-fabrication properties, not fit parameters (see fit_data() docstring).
VARY_PARAM_LABELS = [
    ('T', 'T'), ('Bfield', 'B'), ('GammaBuf', 'GammaBuf'), ('shift', 'Shift'),
    ('lcell', 'l_cell'), ('Btheta', 'Bθ'), ('Bphi', 'Bφ'),
]

# Colours — Durham University light theme
BG        = '#ffffff'           # page / plot background
PANEL     = '#f6f2f5'           # sidebar panels (light Palatinate tint)
BORDER    = '#d4c6d3'           # panel borders
ACCENT    = '#00638A'           # d_blue — secondary buttons
HIGHLIGHT = '#7E317B'           # d_purple (Palatinate) — primary accent
TEXT      = '#231F20'           # d_black — body text
SUBTEXT   = '#68625F'           # mid-grey — secondary / annotation text
PLOT_BG   = '#ffffff'           # plot area background

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

def ni(id_, val, mn=None, mx=None, step=None):
    return dcc.Input(id=id_, type='number', value=val, min=mn, max=mx,
                     step=step, debounce=True, style=INPUT_STYLE)

def row(*kids, gap='8px'):
    return html.Div(kids, style={'display':'flex','gap':gap,'alignItems':'flex-start','marginBottom':'6px'})

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


def _pol_efield(pol, theta0_deg):
    t = np.radians(float(theta0_deg or 0))
    m = {
        'Linear (x)':  [np.cos(t), np.sin(t), 0.],
        'Linear (y)':  [-np.sin(t), np.cos(t), 0.],
        'Left CP':     [1/np.sqrt(2), 1j/np.sqrt(2), 0.],
        'Right CP':    [1/np.sqrt(2), -1j/np.sqrt(2), 0.],
        'Linear 45°':  [1/np.sqrt(2), 1/np.sqrt(2), 0.],
    }
    return m.get(pol, [1., 0., 0.])


def _p(elem, dline, T, lcell, B, gbuf, shift, btheta, bphi, theta0,
       bufgas, bufP, bufTf,
       rb85frac=None, K40frac=None, K41frac=None):
    # GUI angle inputs are in degrees (matching the ES3 wx GUI and the other
    # angle controls in the panel); the core p_dict convention is radians.
    d = dict(Elem=elem, Dline=dline,
             T=float(T or 20), lcell=float(lcell or 0.075),
             Bfield=float(B or 0), GammaBuf=float(gbuf or 0),
             shift=float(shift or 0),
             Btheta=np.radians(float(btheta or 0)), Bphi=np.radians(float(bphi or 0)),
             theta0=float(theta0 or 0), Constrain=True)
    if bufgas and bufgas != 'None':
        d.update(BufferGas=bufgas,
                 BufferPressure=float(bufP or 0),
                 BufferTemp_fill=float(bufTf or 293.15))
    if elem == 'Rb' and rb85frac is not None:
        d['rb85frac'] = float(rb85frac)
    if elem == 'K':
        if K40frac is not None:
            d['K40frac'] = float(K40frac)
        if K41frac is not None:
            d['K41frac'] = float(K41frac)
    return d


def _X(start, stop, N):
    return np.linspace(float(start or -5)*1e3, float(stop or 5)*1e3, int(N or 500))


def _fig(X_GHz, traces):
    fig = go.Figure()
    for i, (name, y, dash_, alpha) in enumerate(traces):
        fig.add_trace(go.Scatter(
            x=X_GHz, y=np.real(y),
            name=name, mode='lines', opacity=alpha,
            line=dict(color=PALETTE[i % len(PALETTE)], width=2,
                      dash='dot' if dash_ else 'solid'),
        ))
    layout = PLOT_LAYOUT.copy()
    layout['xaxis'] = {**layout['xaxis'], 'title': 'Detuning (GHz)'}
    layout['yaxis'] = {**layout['yaxis'], 'title': 'Signal'}
    fig.update_layout(**layout)
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Fitting helpers  [ES4 NEW]
# ─────────────────────────────────────────────────────────────────────────────
def _fit_vary_options(bufgas):
    """Grey out GammaBuf/shift as vary-able whenever BufferGas is active - the
    dynamic buffer gas model silently overrides them (see fit_data() docstring),
    so letting a user tick them here would produce a confusing ValueError."""
    disable_buffer = bufgas not in (None, 'None')
    return [{'label': f' {label}', 'value': key,
             'disabled': disable_buffer and key in ('GammaBuf', 'shift')}
            for key, label in VARY_PARAM_LABELS]


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


def _parse_fit_csv(contents):
    """Parse an uploaded 2-column CSV (detuning GHz, signal) into (X_ghz, Y) arrays."""
    _, content_string = contents.split(',', 1)
    text = base64.b64decode(content_string).decode('utf-8', errors='replace')
    for kwargs in ({'delimiter': ','}, {'delimiter': ',', 'skiprows': 1}, {}):
        try:
            arr = np.loadtxt(io.StringIO(text), **kwargs)
            break
        except ValueError:
            continue
    else:
        raise ValueError('Could not parse file as a 2-column CSV (detuning, signal)')
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError('CSV must have two columns: detuning (GHz), signal')
    return arr[:, 0], arr[:, 1]


# ─────────────────────────────────────────────────────────────────────────────
# Reusable param panel  (prefix determines all ids)
# ─────────────────────────────────────────────────────────────────────────────
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
            field('T [°C]',     ni(f'{prefix}-T',    20.0,  step=1)),
            field('l [m]',      ni(f'{prefix}-lcell', 0.075, step=0.001)),
        ),
        row(
            field('B [G]',      ni(f'{prefix}-B',    0.0,   step=10)),
            field('GammaBuf [MHz]', ni(f'{prefix}-gbuf', 0.0, step=0.5)),
        ),
        row(field('Shift [MHz]', ni(f'{prefix}-shift', 0.0, step=1))),
        sec('Magnetic Field Angles'),
        row(
            field('Bθ [deg]', ni(f'{prefix}-btheta', 0.0, mn=0, mx=180, step=1),
                  sub='0=Faraday  90=Voigt'),
            field('Bφ [deg]', ni(f'{prefix}-bphi',   0.0, mn=0, mx=360, step=1)),
        ),
        *([
            sec('Polarisation'),
            row(
                field('Pol.', dcc.Dropdown(POL_OPTIONS, value='Linear (x)',
                              id=f'{prefix}-pol', clearable=False, style=DROPDOWN_STYLE)),
                field('θ₀ [deg]', ni(f'{prefix}-theta0', 0.0, step=1),
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
            field('P [Torr]', ni(f'{prefix}-bufP', 0.0, mn=0, step=1)),
        ),
        row(
            field('Fill T [K]', ni(f'{prefix}-bufTf', 293.15, mn=0, step=1),
                  sub='sealed cell'),
        ),
        # Isotope abundance controls — shown only for elements with non-trivial mixtures
        html.Div(id=f'{prefix}-iso-rb-row', style={'display':'none'}, children=[
            sec('Isotope Abundance'),
            row(field('⁸⁵Rb (%)', ni(f'{prefix}-rb85frac', 72.17, mn=0, mx=100, step=0.01),
                      sub='natural 72.17')),
        ]),
        html.Div(id=f'{prefix}-iso-k-row', style={'display':'none'}, children=[
            sec('Isotope Abundances'),
            row(
                field('⁴⁰K (%)', ni(f'{prefix}-K40frac', 0.01, mn=0, mx=100, step=0.001),
                      sub='natural 0.01'),
                field('⁴¹K (%)', ni(f'{prefix}-K41frac', 6.73, mn=0, mx=100, step=0.01),
                      sub='natural 6.73'),
            ),
        ]),
    ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
              'border':f'1px solid {BORDER}','marginBottom':'10px'})


def det_panel(prefix):
    return html.Div([
        sec('Detuning Axis'),
        row(
            field('Start [GHz]', ni(f'{prefix}-det-start', -5.0, step=0.5)),
            field('Stop [GHz]',  ni(f'{prefix}-det-stop',   5.0, step=0.5)),
            field('N pts',       ni(f'{prefix}-det-N', 500, mn=10, mx=5000, step=50)),
        ),
    ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
              'border':f'1px solid {BORDER}','marginBottom':'10px'})


def out_panel(prefix, allowed=None):
    opts = allowed if allowed is not None else OUTPUTS
    return html.Div([
        sec('Outputs to Plot'),
        dcc.Checklist(
            id=f'{prefix}-outputs',
            options=[{'label': html.Span(OUTPUT_LABELS[o],
                                         style={'color':TEXT,'fontSize':'0.82rem','marginLeft':'4px'}),
                      'value': o} for o in opts],
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
            field('Fit min [GHz]', ni(f'{prefix}-fit-min', None, step=0.5),
                  sub='blank = full range'),
            field('Fit max [GHz]', ni(f'{prefix}-fit-max', None, step=0.5)),
        ),
        lbl('Vary in fit'),
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
app = dash.Dash(
    __name__,
    external_stylesheets=[
        dbc.themes.FLATLY,
        'https://fonts.googleapis.com/css2?family=Inter:wght@400;500;700&display=swap',
    ],
    suppress_callback_exceptions=False,
    title='ElecSus GUI',
)

# ── Header ────────────────────────────────────────────────────────────────────
header = html.Div([
    html.Div([
        html.H1('ElecSus', style={'color':HIGHLIGHT,'fontWeight':'700','margin':'0',
                                   'fontSize':'1.7rem','letterSpacing':'-0.02em'}),
        html.P('Atomic Spectroscopy · Durham QLM',
               style={'color':SUBTEXT,'margin':'0','fontSize':'0.82rem'}),
    ]),
], style={'padding':'12px 24px','backgroundColor':PANEL,
          'borderBottom':f'2px solid {HIGHLIGHT}'})

# ── Tab bar ───────────────────────────────────────────────────────────────────
tab_bar = html.Div([
    html.Button('Single Cell',         id='tab-single-btn', n_clicks=0,
                style={**BTN_SEC,'marginRight':'4px','borderRadius':'6px 6px 0 0'}),
    html.Button('Cascaded Two Cells',  id='tab-cascade-btn', n_clicks=0,
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
                row(field('T [°C]', ni('c2-T', 20.0, step=1)),
                    field('l [m]',  ni('c2-lcell', 0.075, step=0.001))),
                row(field('B [G]',  ni('c2-B', 0.0, step=10)),
                    field('GammaBuf [MHz]', ni('c2-gbuf', 0.0, step=0.5))),
                row(field('Shift [MHz]', ni('c2-shift', 0.0, step=1))),
                sec('Magnetic Field Angles'),
                row(field('Bθ [deg]', ni('c2-btheta', 0.0, mn=0, mx=180, step=1),
                          sub='0=Faraday  90=Voigt'),
                    field('Bφ [deg]', ni('c2-bphi',   0.0, mn=0, mx=360, step=1))),
                sec('Buffer Gas'),
                row(field('Gas', dcc.Dropdown(BUFFER_GAS_OPTIONS, value='None', id='c2-bufgas',
                                 clearable=False, style=DROPDOWN_STYLE)),
                    field('P [Torr]', ni('c2-bufP', 0.0, mn=0, step=1))),
                row(field('Fill T [K]', ni('c2-bufTf', 293.15, mn=0, step=1), sub='sealed cell')),
            ], style={'padding':'12px','backgroundColor':PANEL,'borderRadius':'8px',
                      'border':f'1px solid {BORDER}','marginBottom':'10px'}),

            # Intermediate polariser — default ON at 90° (crossed w.r.t. an
            # x-polarised input), the standard MOF filter-section topology.
            html.Div([
                sec('Intermediate Optics'),
                dcc.Checklist(id='cas-use-pol',
                    options=[{'label': html.Span(' Polariser between cells',
                                                  style={'color':TEXT,'fontSize':'0.83rem'}),
                              'value':'on'}],
                    value=['on'], inputStyle={'accentColor':HIGHLIGHT}),
                html.Div(style={'height':'4px'}),
                row(field('Angle [deg]', ni('cas-pol-angle', 90.0, step=1),
                          sub='90 = crossed (MOF)')),
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
                    sub='wing-selector: fixed circular into Cell 2'),
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

# ── Error toast ───────────────────────────────────────────────────────────────
error_toast = dbc.Toast(id='error-toast', header='Calculation Error',
                        is_open=False, dismissable=True, duration=7000, icon='danger',
                        style={'position':'fixed','bottom':'20px','right':'20px',
                               'zIndex':9999,'maxWidth':'420px'})

app.layout = html.Div([
    header, tab_bar,
    html.Div([single_tab, cascade_tab],
             style={'backgroundColor':BG,'minHeight':'calc(100vh - 110px)'}),
    error_toast,
    dcc.Store(id='error-store', data=''),
], style={'backgroundColor':BG,'minHeight':'100vh',
          'fontFamily':'Inter,Arial,sans-serif','margin':'0'})


# ─────────────────────────────────────────────────────────────────────────────
# Tab switching
# ─────────────────────────────────────────────────────────────────────────────
@app.callback(
    Output('single-tab-content',  'style'),
    Output('cascade-tab-content', 'style'),
    Output('tab-single-btn',  'style'),
    Output('tab-cascade-btn', 'style'),
    Input('tab-single-btn',  'n_clicks'),
    Input('tab-cascade-btn', 'n_clicks'),
)
def switch_tab(ns, nc):
    single_active  = (dash.callback_context.triggered_id != 'tab-cascade-btn')
    active_style   = {**BTN, 'marginRight':'4px','borderRadius':'6px 6px 0 0',
                      'borderBottom':f'2px solid {HIGHLIGHT}'}
    active_cas     = {**BTN, 'borderRadius':'6px 6px 0 0',
                      'borderBottom':f'2px solid {HIGHLIGHT}'}
    inactive_style = {**BTN_SEC,'marginRight':'4px','borderRadius':'6px 6px 0 0',
                      'borderBottom':'none'}
    inactive_cas   = {**BTN_SEC,'borderRadius':'6px 6px 0 0','borderBottom':'none'}
    if single_active:
        return ({'display':'block'}, {'display':'none'},
                active_style, inactive_cas)
    return ({'display':'none'}, {'display':'block'},
            inactive_style, active_cas)


# ─────────────────────────────────────────────────────────────────────────────
# Single-cell compute
# ─────────────────────────────────────────────────────────────────────────────
@app.callback(
    Output('s-plot',  'figure'),
    Output('s-store', 'data'),
    Output('error-store', 'data'),
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
    prevent_initial_call=True,
)
def compute_single(n, elem, dline, T, lcell, B, gbuf, shift,
                   btheta, bphi, theta0, pol,
                   bufgas, bufP, bufTf,
                   rb85frac, K40frac, K41frac,
                   det_start, det_stop, det_N, outputs):
    try:
        p = _p(elem, dline, T, lcell, B, gbuf, shift,
               btheta, bphi, theta0, bufgas, bufP, bufTf,
               rb85frac=rb85frac, K40frac=K40frac, K41frac=K41frac)
        E_in = _pol_efield(pol, theta0)
        X = _X(det_start, det_stop, det_N)
        outs = outputs or ['S0']
        res = calculate(X, E_in, p, outputs=outs)
        traces = [(OUTPUT_LABELS.get(k,k), np.real(v), False, 1.0)
                  for k, v in zip(outs, res)]
        fig = _fig(X/1e3, traces)
        fig.update_layout(title=dict(
            text=f'{elem} {dline} · T={T}°C · B={B} G · l={lcell} m',
            font=dict(color=TEXT, size=13)))
        store = {'X': list(X/1e3), 'results': {k: list(np.real(v)) for k,v in zip(outs,res)}}
        return fig, store, ''
    except Exception as e:
        return empty_fig('Calculation failed'), no_update, str(e)


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
        X_ghz, Y = _parse_fit_csv(contents)
        store = {'X': X_ghz.tolist(), 'Y': Y.tolist()}
        return store, f'{filename}: {len(X_ghz)} points loaded'
    except Exception as e:
        return no_update, f'Failed to parse {filename}: {e}'


@app.callback(
    Output('s-fit-vary', 'options'),
    Input('s-bufgas', 'value'),
)
def update_fit_vary_options(bufgas):
    return _fit_vary_options(bufgas)


@app.callback(
    Output('s-fit-plot', 'figure'),
    Output('s-fit-results', 'children'),
    Output('error-store', 'data', allow_duplicate=True),
    Output('s-T', 'value'), Output('s-B', 'value'),
    Output('s-gbuf', 'value'), Output('s-shift', 'value'),
    Output('s-lcell', 'value'), Output('s-btheta', 'value'), Output('s-bphi', 'value'),
    Output('s-fit-result-store', 'data'),
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
    prevent_initial_call=True,
)
def run_fit_single(n, data, data_type, algo, vary_keys, fit_min, fit_max,
                   elem, dline, T, lcell, B, gbuf, shift, btheta, bphi, theta0, pol,
                   bufgas, bufP, bufTf, rb85frac, K40frac, K41frac):
    no_change = (no_update,) * 7
    if not data:
        return (empty_fig('Upload data first'), '',
                'No experimental data uploaded yet.', *no_change, no_update)
    try:
        X_ghz = np.array(data['X'], dtype=float)
        Y     = np.array(data['Y'], dtype=float)
        if fit_min is not None and fit_max is not None:
            mask = (X_ghz >= float(fit_min)) & (X_ghz <= float(fit_max))
            X_ghz, Y = X_ghz[mask], Y[mask]
        if len(X_ghz) < 2:
            raise ValueError('Not enough data points in the selected fit range')
        X_mhz = X_ghz * 1e3

        p0 = _p(elem, dline, T, lcell, B, gbuf, shift, btheta, bphi, theta0,
               bufgas, bufP, bufTf, rb85frac=rb85frac, K40frac=K40frac, K41frac=K41frac)
        E_in = np.array(_pol_efield(pol, theta0), dtype='complex')
        E_in_angle = [E_in[0].real, [abs(E_in[1]), np.angle(E_in[1])]]

        vary = {key: (key in (vary_keys or [])) for key, _ in VARY_PARAM_LABELS}

        # differential_evolution requires finite bounds on every varying parameter,
        # and the GUI deliberately doesn't collect explicit bounds - build a search
        # window around each initial guess instead: ±max(50% of |value|, floor),
        # where the floor keeps the window usable when the guess is near zero.
        p_dict_bounds = None
        if algo != 'ML':
            _floor = {'T': 20.0, 'Bfield': 200.0, 'GammaBuf': 100.0, 'shift': 100.0,
                      'lcell': 0.5 * float(p0.get('lcell', 0.075)),
                      'Btheta': 0.5, 'Bphi': 0.5}
            p_dict_bounds = {}
            for key, _label in VARY_PARAM_LABELS:
                if vary.get(key):
                    v0 = float(p0.get(key, 0.0))
                    half = max(abs(v0) * 0.5, _floor.get(key, 1.0))
                    p_dict_bounds[key] = [v0 - half, v0 + half]

        best, rms, result = fit_data([X_mhz, Y], p0, vary, E_in=E_in_angle,
                                     p_dict_bounds=p_dict_bounds,
                                     data_type=data_type, fit_algorithm=algo)

        fig = _fit_fig(X_ghz, Y, result.best_fit, OUTPUT_LABELS.get(data_type, data_type), rms)

        rows = []
        for key, label in VARY_PARAM_LABELS:
            if vary.get(key):
                val = best.get(key)
                param = result.params.get(key) if hasattr(result, 'params') else None
                stderr = getattr(param, 'stderr', None) if param is not None else None
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
        results_children = html.Div([
            html.P(f'RMS residual: {rms:.6g}',
                   style={'fontWeight':'700','color':TEXT,'fontSize':'0.85rem'}),
            html.Table([html.Tbody(rows)], style={'fontSize':'0.8rem','color':TEXT}),
        ])

        fit_store = {'X': X_ghz.tolist(), 'data': Y.tolist(),
                     'fit': list(np.real(np.asarray(result.best_fit))),
                     'type': data_type}

        return (fig, results_children, '',
                best.get('T', T), best.get('Bfield', B),
                best.get('GammaBuf', gbuf), best.get('shift', shift),
                best.get('lcell', lcell),
                # best-fit angles are radians (core convention); inputs are degrees
                np.degrees(best['Btheta']) if 'Btheta' in best else btheta,
                np.degrees(best['Bphi'])   if 'Bphi'   in best else bphi,
                fit_store)
    except Exception as e:
        return empty_fig('Fit failed'), '', str(e), *no_change, no_update


# ─────────────────────────────────────────────────────────────────────────────
# Cascade compute
# ─────────────────────────────────────────────────────────────────────────────
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
    prevent_initial_call=True,
)
def compute_cascade(n,
                    e1, d1, T1, l1, B1, g1, sh1, bt1, bp1, th1, pol1, bg1, bP1, bTf1,
                    rb85, k40, k41,
                    T2, l2, B2, g2, sh2, bt2, bp2, bg2, bP2, bTf2,
                    use_pol, pol_angle, reset_pol,
                    det_start, det_stop, det_N, outputs):
    """Cascaded same-species filters via the core calculate_cascade().

    With the intermediate polariser and/or QWP on, cell-2 outputs from the
    core are instrument totals normalised to the laser input. With neither,
    totals are formed by the documented combination rules (intensities x
    S0_cell1; fractional Stokes pass through). The QWP option models a
    wing-selector topology (e.g. VAMOS, Oliviero 1998): it prepares a FIXED
    circular state for cell 2 regardless of wavelength, so cell 2 always
    absorbs the same Zeeman component."""
    try:
        from elecsus.elecsus_methods import calculate_cascade

        # Both cells share element/D-line/isotope fractions (same vapour species)
        p1 = _p(e1, d1, T1, l1, B1, g1, sh1, bt1, bp1, th1, bg1, bP1, bTf1,
                rb85frac=rb85, K40frac=k40, K41frac=k41)
        p2 = _p(e1, d1, T2, l2, B2, g2, sh2, bt2, bp2, 0., bg2, bP2, bTf2,
                rb85frac=rb85, K40frac=k40, K41frac=k41)

        E_in = np.array(_pol_efield(pol1, th1), dtype='complex')
        X    = _X(det_start, det_stop, det_N)
        outs = [k for k in (outputs or ['S0']) if k in CASCADE_OUTPUTS] or ['S0']

        pol_on = 'on' in (use_pol or [])
        pol_rad = np.radians(float(pol_angle or 0)) if pol_on else None
        reset_val = reset_pol if reset_pol in ('left_cp', 'right_cp') else None

        res1, res2 = calculate_cascade(X, E_in, p1, p2, outputs=outs,
                                       polariser_angle=pol_rad,
                                       reset_polarisation=reset_val)
        r1 = {k: np.asarray(v).squeeze().real for k, v in zip(outs, res1)}

        if pol_on or reset_val is not None:
            # core returns instrument totals directly with an intermediate optic
            total = {k: np.asarray(v).squeeze().real for k, v in zip(outs, res2)}
        else:
            # per-stage outputs: combine by the documented rules
            r2 = {k: np.asarray(v).squeeze().real for k, v in zip(outs, res2)}
            S0_1 = r1.get('S0')
            if S0_1 is None:
                [S0_1] = calculate(X, E_in, p1, outputs=['S0'])
                S0_1 = np.asarray(S0_1).squeeze().real
            total = {k: (r2[k] if k in ('S1', 'S2', 'S3') else S0_1 * r2[k])
                     for k in outs}

        traces = []
        for k in outs:
            traces.append((f'After Cell 1 — {OUTPUT_LABELS.get(k,k)}', r1[k], True,  0.7))
            traces.append((f'Total after both filters — {OUTPUT_LABELS.get(k,k)}',
                           total[k], False, 1.0))

        fig = _fig(X/1e3, traces)
        pol_txt = f' · polariser at {float(pol_angle or 0):g}°' if pol_on else ' · no polariser'
        qwp_txt = f' · QWP → {reset_val}' if reset_val is not None else ''
        fig.update_layout(title=dict(
            text=f'Cascaded filters: {e1} {d1} → {e1} {d1}{pol_txt}{qwp_txt}',
            font=dict(color=TEXT, size=13)))

        store = {'X': list(X/1e3),
                 'cell1': {k: list(v) for k, v in r1.items()},
                 'total': {k: list(v) for k, v in total.items()}}
        return fig, store, ''
    except Exception as e:
        return empty_fig('Calculation failed'), no_update, str(e)


# ─────────────────────────────────────────────────────────────────────────────
# Error toast (single callback watching the shared error-store)
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


# ─────────────────────────────────────────────────────────────────────────────
# CSV exports
# ─────────────────────────────────────────────────────────────────────────────
@app.callback(Output('s-download','data'), Input('s-export','n_clicks'),
              State('s-store','data'), prevent_initial_call=True)
def export_single(n, store):
    if not store: return no_update
    X = np.array(store['X'])
    keys = list(store['results'].keys())
    lines = ['Detuning_GHz,' + ','.join(keys)]
    for i,x in enumerate(X):
        vals = [f'{np.array(store["results"][k])[i]:.8f}' for k in keys]
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
        lines.append(f'{x:.6f},{data[i]:.8f},{fit[i]:.8f},{data[i]-fit[i]:.8f}')
    return dict(content='\n'.join(lines), filename='elecsus_fit.csv')


@app.callback(Output('cas-download','data'), Input('cas-export','n_clicks'),
              State('cas-store','data'), prevent_initial_call=True)
def export_cascade(n, store):
    if not store: return no_update
    X = np.array(store['X'])
    keys = list(store['total'].keys())
    hdr = 'Detuning_GHz,' + ','.join([f'Cell1_{k}' for k in keys]) + ',' + \
          ','.join([f'Total_{k}' for k in keys])
    lines = [hdr]
    for i,x in enumerate(X):
        c1 = [f'{np.array(store["cell1"].get(k,[0]*len(X)))[i]:.8f}' for k in keys]
        c2 = [f'{np.array(store["total"][k])[i]:.8f}' for k in keys]
        lines.append(f'{x:.6f},' + ','.join(c1) + ',' + ','.join(c2))
    return dict(content='\n'.join(lines), filename='elecsus_cascade.csv')


# ─────────────────────────────────────────────────────────────────────────────
# GammaBuf greying — disabled + styled when a BufferGas species is selected,
# because the dynamic Pitz model silently overwrites GammaBuf and shift.
# ─────────────────────────────────────────────────────────────────────────────
def _gbuf_state(bufgas):
    active = bufgas not in (None, 'None')
    style  = INPUT_DISABLED_STYLE if active else INPUT_STYLE
    return active, style

@app.callback(
    Output('s-gbuf', 'disabled'),
    Output('s-gbuf', 'style'),
    Input('s-bufgas', 'value'),
)
def grey_gbuf_single(bufgas):
    disabled, style = _gbuf_state(bufgas)
    return disabled, style

@app.callback(
    Output('c1-gbuf', 'disabled'),
    Output('c1-gbuf', 'style'),
    Input('c1-bufgas', 'value'),
)
def grey_gbuf_c1(bufgas):
    disabled, style = _gbuf_state(bufgas)
    return disabled, style

@app.callback(
    Output('c2-gbuf', 'disabled'),
    Output('c2-gbuf', 'style'),
    Input('c2-bufgas', 'value'),
)
def grey_gbuf_c2(bufgas):
    disabled, style = _gbuf_state(bufgas)
    return disabled, style


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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='ElecSus Dash GUI')
    parser.add_argument('--port',  type=int, default=8050, help='TCP port (default 8050)')
    parser.add_argument('--host',  default='127.0.0.1',   help='Bind address (default 127.0.0.1)')
    parser.add_argument('--debug', action='store_true',   help='Enable Dash debug/reloader mode')
    args = parser.parse_args()
    run(debug=args.debug, port=args.port, host=args.host)
