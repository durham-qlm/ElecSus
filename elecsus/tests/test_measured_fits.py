"""
Fits to measured spectra, against the published values.

test_measured_data.py checks the forward model against these spectra with the
published parameters. Here the parameters are *fitted*, by fit_data and by the
GUI's fit callback, from starting values well away from the answer, and the
fitted values are compared with what the papers report. The data and their
setup are those of test_measured_data.py (CC BY 4.0, Durham Research Data
Repository); they live in the development repository only, so these tests
skip in an installed wheel.

- Briscoe et al., Opt. Express 33, 32099 (2025). Rb D2 through a 2 mm cell in
  an oblique field.
  - Regime I (239 G, theta_B = 80 deg, 119.69 C, 7.9 MHz): T, B, theta_B and
    GammaBuf fitted to each of the four polarisations, with Briscoe's input
    angle; and once with the input angle fitted as well.
  - Regime II (about 2580 G, theta_B = 99.7-129.7 deg, about 80 C): T, B and
    theta_B fitted to all eight spectra. For linear input along x or y,
    theta_B and 180 deg - theta_B give the same spectrum exactly, so the fit
    may land on either; angles are compared folded into [0, 90] deg. The
    y-polarised spectra fix theta_B only loosely (4 deg changes S0 by about
    0.01 there, against 0.01-0.04 for x), so they are held to 6 deg.
- Alqarni et al., J. Phys. B 58, 135201 (2025), Fig. 3. K D1 through a 75 mm
  cell, body 370 K and stem 350 K, with 0, 57.4 and 96 Torr of Ne filled at
  293 K.
  - The stem temperature, fitted to the spectrum without gas.
  - The Ne pressure, fitted with T, against the fill pressure.
  - Their own fits: ES4 reproduces their published fit curves (fitting T,
    GammaBuf and shift to the curves), and fitting their data with the same
    free parameters does at least as well as their fits did. Their Table 1
    values (417 MHz / -87 MHz and 700 MHz / -144 MHz) come from a separate
    temperature series (their Fig. 4); the Fig. 3 spectrum at 96 Torr sits
    about 50 MHz from it, so the shift is not compared with the table.
"""

import contextlib
import io
import math
import warnings
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import binned_statistic

from elecsus import elecsus_dash_gui as gui
from elecsus.elecsus_methods import fit_data

DATA = Path(__file__).resolve().parents[2] / 'paper' / 'data'
BRISCOE = DATA / 'briscoe2025'
ALQARNI = DATA / 'alqarni2025' / 'figure 3'

needs_briscoe = pytest.mark.skipif(
    not BRISCOE.is_dir(), reason="Briscoe et al.'s data are not shipped")
needs_alqarni = pytest.mark.skipif(
    not ALQARNI.is_dir(), reason="Alqarni et al.'s data are not shipped")


def _fit(x_ghz, y, p, vary, E_in, **kw):
    """fit_data, quietly: (best values, RMS in %)."""
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('ignore')
        best, rms, _ = fit_data((np.asarray(x_ghz) * 1e3, np.asarray(y)), p,
                                {k: True for k in vary}, E_in=E_in, data_type='S0', **kw)
    return best, 100 * rms


def _binned(path, window, bins=400):
    """Briscoe's spectra as in his figure scripts: clipped and binned."""
    d = np.loadtxt(path, delimiter=',')
    x, y = d[:, 0], d[:, 1]
    keep = (x >= window[0]) & (x <= window[1])
    x, y = x[keep], y[keep]
    return (binned_statistic(x, x, statistic='mean', bins=bins)[0],
            binned_statistic(x, y, statistic='mean', bins=bins)[0])


def _linear(theta_deg):
    t = math.radians(theta_deg)
    return [math.cos(t), math.sin(t), 0.0]


def _folded(theta_rad):
    """theta_B folded into [0, 90] deg: theta, 180 - theta, -theta and
    theta + 360 give the same spectrum for linear input along x or y."""
    return math.degrees(math.acos(abs(math.cos(theta_rad))))


# --- Briscoe et al. 2025, Regime I --------------------------------------------

REGIME_ONE = [('tau_thetaB80_thetaE90_1.csv', 90.9), ('tau_thetaB80_thetaE20_2.csv', 19.8),
              ('tau_thetaB80_thetaE50_1.csv', 49.1), ('tau_thetaB80_thetaEminus50_2.csv', -48.1)]
REGIME_ONE_START = dict(Elem='Rb', Dline='D2', T=110., lcell=2e-3, rb85frac=72.17,
                        Btheta=math.radians(70.), Bfield=210., Bphi=0., shift=0., GammaBuf=20.)


@needs_briscoe
@pytest.mark.parametrize('fname, theta_E', REGIME_ONE)
def test_regime_one_fit_recovers_briscoes_values(fname, theta_E):
    """From 110 C, 210 G, 70 deg and 20 MHz, each polarisation fits to
    Briscoe's 119.69 C, 239 G, 80 deg and 7.9 MHz (found: 119.0-120.1 C,
    239-244 G, 79.5-79.9 deg, 6.1-8.8 MHz, at 0.42-0.50 % RMS)."""
    x, y = _binned(BRISCOE / 'Regime I' / fname, (-6.25, 7.75))
    best, rms = _fit(x, y, REGIME_ONE_START, ['T', 'Bfield', 'Btheta', 'GammaBuf'],
                     _linear(theta_E))
    assert rms < 0.6
    assert abs(best['T'] - 119.69) < 1.0
    assert abs(best['Bfield'] - 239.0) < 6.0
    assert abs(math.degrees(best['Btheta']) - 80.0) < 1.0
    assert abs(best['GammaBuf'] - 7.9) < 2.5


@needs_briscoe
def test_regime_one_fit_recovers_the_input_angle_too():
    """With the input angle free as well, started 9 deg away (found: 48.6 deg
    against Briscoe's 49.1)."""
    x, y = _binned(BRISCOE / 'Regime I' / 'tau_thetaB80_thetaE50_1.csv', (-6.25, 7.75))
    t = math.radians(40.)
    best, rms = _fit(x, y, REGIME_ONE_START, ['T', 'Bfield', 'Btheta', 'GammaBuf', 'E_y'],
                     [math.cos(t), [math.sin(t), 0.]])
    assert rms < 0.6
    assert abs(math.degrees(math.atan2(best['E_y'], best['E_x'])) - 49.1) < 1.5
    assert abs(best['T'] - 119.69) < 1.0


# --- Briscoe et al. 2025, Regime II -------------------------------------------

# (theta_B label, Briscoe's fitted theta_B and B, the repeats of his Fig. 6)
REGIME_TWO = [(100, 99.7, 2578.2, 5, 3), (110, 109.7, 2581.3, 1, 3),
              (120, 119.7, 2584.9, 5, 3), (130, 129.7, 2590.2, 3, 1)]
# (theta_E label, Briscoe's fitted T and theta_E)
REGIME_TWO_E = {90: (79.52, 90.8), 0: (80.0, 0.8)}


@needs_briscoe
@pytest.mark.parametrize('theta_E', [0, 90])
@pytest.mark.parametrize('label, theta_B, B, repeat_90, repeat_0', REGIME_TWO)
def test_regime_two_fit_recovers_briscoes_values(label, theta_B, B, repeat_90, repeat_0, theta_E):
    """From 75 C, 2400 G and theta_B + 5 deg (GammaBuf at Briscoe's 7.9 MHz):
    T within 0.5 C and B within 0.3 % for all eight spectra (found: 0.25 C and
    0.16 %); theta_B within 1 deg for x-polarised light (found: 0.5) and 6 deg
    for y (found: 2-5)."""
    repeat = repeat_90 if theta_E == 90 else repeat_0
    T, E_angle = REGIME_TWO_E[theta_E]
    x, y = _binned(BRISCOE / 'Regime II' / f'tau_thetaB{label}_thetaE{theta_E}_{repeat}.csv',
                   (-12.5, 12.5))
    start = dict(Elem='Rb', Dline='D2', T=75., lcell=2e-3, rb85frac=72.17,
                 Btheta=math.radians(label + 5.), Bfield=2400., Bphi=0., shift=0., GammaBuf=7.9)
    best, rms = _fit(x, y, start, ['T', 'Bfield', 'Btheta'], _linear(E_angle))
    assert rms < 0.5
    assert abs(best['T'] - T) < 0.5
    assert abs(best['Bfield'] - B) < 0.003 * B
    off = abs(_folded(best['Btheta']) - _folded(math.radians(theta_B)))
    assert off < (1.0 if theta_E == 0 else 6.0)


# --- Alqarni et al. 2025, Fig. 3 ----------------------------------------------

BODY, STEM = 370.0 - 273.15, 350.0 - 273.15
ALQ = dict(Elem='K', Dline='D1', lcell=75e-3, DoppTemp=BODY, Constrain=False, Bfield=0.)


def _alqarni(fname, column=2):
    """(detuning in GHz, column) of a Fig. 3 file; their header row, where
    there is one, is skipped."""
    rows = []
    with open(ALQARNI / fname) as f:
        for line in f:
            parts = line.strip().split(',')
            try:
                rows.append((float(parts[1]), float(parts[column])))
            except (ValueError, IndexError):
                continue
    d = np.array(rows)
    return d[:, 0], d[:, 1]


@needs_alqarni
def test_alqarni_stem_temperature_from_the_spectrum_without_gas():
    """From 60 C, with the body at 370 K: the stem's 76.85 C (found: 76.0 C)."""
    x, y = _alqarni('No buffer gas.csv')
    best, rms = _fit(x, y, dict(ALQ, T=60.), ['T'], [1, 0, 0])
    assert rms < 0.6
    assert abs(best['T'] - STEM) < 1.5


@needs_alqarni
@pytest.mark.parametrize('fname, fill', [('60 Torr.csv', 57.4), ('100 torr.csv', 96.0)])
def test_alqarni_ne_pressure_against_the_fill(fname, fill):
    """The buffer-gas model, with the pressure and T fitted from 30 Torr and
    the stem temperature: within 8 % of the fill (found: 59.8 and 100.7 Torr,
    +4 % and +5 %), T within 2.5 C of the stem (found: 75.0 and 76.0 C)."""
    x, y = _alqarni(fname)
    p = dict(ALQ, T=STEM, BufferGas='Ne', BufferPressure=30., BufferTemp_fill=293.)
    best, rms = _fit(x, y, p, ['BufferPressure', 'T'], [1, 0, 0])
    assert rms < 1.6
    assert abs(best['BufferPressure'] - fill) < 0.08 * fill
    assert abs(best['T'] - STEM) < 2.5


@needs_alqarni
@pytest.mark.parametrize('fname, fit_name, Gamma, shift', [
    ('60 Torr.csv', '60 Torr fit.csv', 417., -87.),
    ('100 torr.csv', '100 torr fit.csv', 700., -144.),
])
def test_es4_fits_alqarnis_data_as_well_as_their_own_fits(fname, fit_name, Gamma, shift):
    """ES4 reproduces Alqarni et al.'s fit curves (T, GammaBuf and shift free:
    0.10 and 0.03 % RMS), and fitting their data with the same free
    parameters is at least as good as their fits (0.81 against 0.81 %, and
    0.96 against 1.02 %)."""
    x, y = _alqarni(fname)
    xf, yf = _alqarni(fit_name)
    p = dict(ALQ, T=STEM, GammaBuf=Gamma, shift=shift)
    _, rms_curve = _fit(xf, yf, p, ['GammaBuf', 'shift', 'T'], [1, 0, 0])
    assert rms_curve < 0.2
    theirs = 100 * np.sqrt(np.mean((np.interp(x, xf, yf) - y) ** 2))
    _, rms_data = _fit(x, y, p, ['GammaBuf', 'shift', 'T'], [1, 0, 0])
    assert rms_data <= theirs + 0.02


# --- the same through the GUI's fit button ------------------------------------

def _gui_fit(store, vary, **cell):
    a = dict(elem='Rb', dline='D2', T=110., lcell=2e-3, B=210., gbuf=20., shift=0.,
             btheta=70., bphi=0., theta0=0., pol='Linear (x)', bufgas='None', bufP=0.,
             bufTf=293.15, rb85frac=72.17, K40frac=0.01, K41frac=6.73)
    doppT = cell.pop('doppT', None)
    a.update(cell)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('ignore')
        return gui.run_fit_single(1, store, 'S0', 'ML', vary, None, None, *a.values(),
                                  doppT=doppT)


@needs_briscoe
def test_gui_fit_of_regime_one_writes_back_briscoes_values():
    """The fit button on Briscoe's theta_E = 90.9 deg spectrum: the boxes get
    T, B, theta_B (in degrees) and GammaBuf near his values."""
    x, y = _binned(BRISCOE / 'Regime I' / 'tau_thetaB80_thetaE90_1.csv', (-6.25, 7.75))
    out = _gui_fit({'X': x.tolist(), 'Y': y.tolist()}, ['T', 'Bfield', 'Btheta', 'GammaBuf'],
                   theta0=90.9)
    assert out[2] == ''
    T, B, gbuf, btheta = out[3], out[4], out[5], out[8]
    assert abs(T - 119.69) < 1.0 and abs(B - 239.0) < 6.0
    assert abs(btheta - 80.0) < 1.0 and abs(gbuf - 7.9) < 2.5


@needs_alqarni
def test_gui_fit_of_the_ne_pressure_writes_back_the_fill():
    """The fit button on Alqarni's 57.4 Torr spectrum, with Ne selected, the
    Doppler T box at the body temperature and Buffer P and T ticked."""
    x, y = _alqarni('60 Torr.csv')
    out = _gui_fit({'X': x.tolist(), 'Y': y.tolist()}, ['BufferPressure', 'T'],
                   elem='K', dline='D1', T=STEM, lcell=75e-3, B=0., gbuf=0., btheta=0.,
                   bufgas='Ne', bufP=30., bufTf=293., doppT=BODY)
    assert out[2] == ''
    assert abs(out[15] - 57.4) < 0.08 * 57.4
