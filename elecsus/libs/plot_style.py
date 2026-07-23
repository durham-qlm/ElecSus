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
ElecSus house plot style.

Applies Times New Roman typography and the Durham University colour palette to
all subsequent matplotlib figures.  Call ``apply_style()`` once at the top of
any script or notebook before creating axes.

Usage
-----
from elecsus.libs.plot_style import apply_style, DURHAM_CYCLE
from elecsus.libs.plot_style import d_purple, d_blue, d_red  # individual colours

apply_style()               # sets rcParams globally
fig, ax = plt.subplots()    # colour cycle picks Durham colours automatically
"""

import matplotlib.pyplot as plt
from cycler import cycler

from elecsus.libs.durhamcolours import (
    d_purple, d_blue, d_red, d_yellow, d_grey, d_black,
    d_lightpurple, d_midblue, d_pink, d_olive, d_lightblue, d_lightgrey, d_midgrey,
)

# Default colour cycle: purple → blue → red → yellow → grey → black
DURHAM_CYCLE = [d_purple, d_blue, d_red, d_yellow, d_grey, d_black]

_BASE = {
    # Typography
    'font.family':       'Times New Roman',
    'mathtext.fontset':  'custom',
    'mathtext.rm':       'Times New Roman',
    'mathtext.it':       'Times New Roman:italic',
    'mathtext.bf':       'Times New Roman:bold',
    # Sizes (overrideable per-script with plt.rcParams.update after apply_style())
    'font.size':         14,
    'axes.labelsize':    14,
    'axes.titlesize':    12,
    'xtick.labelsize':   14,
    'ytick.labelsize':   14,
    'legend.fontsize':   12,
    # Lines
    'lines.linewidth':   1.5,
    'lines.markersize':  5,
    # Axes
    'axes.grid':         False,
    # Output quality
    'figure.dpi':        150,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'savefig.facecolor': 'white',
    # Font embedding in PDF/PS (required for journal submission)
    'pdf.fonttype':      42,
    'ps.fonttype':       42,
}

# Smaller preset for multi-panel paper figures at ~6–7 inch column width
_PAPER = dict(_BASE, **{
    'font.size':       11,
    'axes.labelsize':  11,
    'axes.titlesize':  10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize':  9,
})


def apply_style(size='default'):
    """Apply the ElecSus/Durham house style to all subsequent matplotlib figures.

    Parameters
    ----------
    size : {'default', 'paper'}
        'default' — 14 pt, suitable for notebooks and standalone figures.
        'paper'   — 11 pt, suited for ~6–7 inch two-column paper figures.
    """
    params = _PAPER if size == 'paper' else _BASE
    plt.rcParams.update(params)
    plt.rcParams['axes.prop_cycle'] = cycler(color=DURHAM_CYCLE)
