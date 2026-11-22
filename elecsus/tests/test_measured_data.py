"""
ElecSus against measured spectra, with no parameter fitted.

Two published datasets (both CC BY 4.0, Durham Research Data Repository):

- Briscoe et al., Opt. Express 33, 32099 (2025), doi:10.15128/r24m90dv58f:
  Rb D2 in Regime I (2 mm natural-abundance cell, Btheta = 80 deg, B = 239 G,
  T = 119.69 C, a residual Lorentzian broadening of 7.9 MHz), at four linear
  input polarisations. The input is [cos theta_E, +sin theta_E, 0] at
  Briscoe's fitted angle theta_E, measured from x towards y. ES4's default
  propagation (non-orthogonal modes, no windows) must reproduce each spectrum
  to 1 % RMS (it gives 0.47-0.54 %); ES3's orthogonal projection
  (use_non_orthogonal=False) must miss each by more than 3 % (3.6-8.1 %).
- Alqarni et al., J. Phys. B 58, 135201 (2025), doi:10.15128/r2kw52j810w,
  Fig. 3: K D1 through a 75 mm cell, body at 370 K (DoppTemp, with
  Constrain=False) and stem at 350 K (T), with 0, 57.4 and 96 Torr of Ne
  filled at 293 K. The broadening and shift are predicted by the buffer-gas
  model, evaluated at the body temperature. The spectra must be reproduced to
  0.8 / 2.1 / 1.8 % RMS (0.67 / 1.88 / 1.56 %), and the predicted broadening
  and shift must lie within the uncertainties of the values Alqarni et al.
  fitted to each spectrum (their Table 1).

The data are kept in the development repository only, two levels above this
package, and are not part of the release: there, and in an installed wheel,
these tests skip. The setup follows the scripts that make the paper's
figures from the same data (Briscoe's spectra binned to 400 points, as there).
"""

import functools
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import binned_statistic

from elecsus.elecsus_methods import calculate
from elecsus.libs.BufferGasData import get_buffer_broadening_and_shift

DATA = Path(__file__).resolve().parents[2] / 'paper' / 'data'
BRISCOE = DATA / 'briscoe2025' / 'Regime I'
ALQARNI = DATA / 'alqarni2025' / 'figure 3'


def _rms_percent(y, model):
    return float(np.sqrt(np.nanmean((100.0 * (y - model)) ** 2)))


# --- Briscoe et al. 2025, Regime I ------------------------------------------

BRISCOE_P = dict(Elem='Rb', Dline='D2', T=119.69, lcell=2e-3, rb85frac=72.17,
                 Btheta=np.radians(80.0), Bfield=239.0, Bphi=0.0, shift=0,
                 GammaBuf=7.9)
# (file: the repeat Briscoe et al. chose for their Fig. 4, their fitted theta_E)
BRISCOE_SPECTRA = [('tau_thetaB80_thetaE90_1.csv', 90.9),
                   ('tau_thetaB80_thetaE20_2.csv', 19.8),
                   ('tau_thetaB80_thetaE50_1.csv', 49.1),
                   ('tau_thetaB80_thetaEminus50_2.csv', -48.1)]
BRISCOE_WINDOW = (-6.25, 7.75)      # GHz, Briscoe's plot window
BRISCOE_BINS = 400


@functools.lru_cache(maxsize=None)
def _briscoe_binned(fname):
    d = np.loadtxt(BRISCOE / fname, delimiter=',')
    x, y = d[:, 0], d[:, 1]
    keep = (x >= BRISCOE_WINDOW[0]) & (x <= BRISCOE_WINDOW[1])
    x, y = x[keep], y[keep]
    xb = binned_statistic(x, x, statistic='mean', bins=BRISCOE_BINS)[0]
    yb = binned_statistic(x, y, statistic='mean', bins=BRISCOE_BINS)[0]
    return xb, yb


def _briscoe_rms(fname, theta_E, **model):
    xb, yb = _briscoe_binned(fname)
    t = np.radians(theta_E)
    E_in = np.array([np.cos(t), np.sin(t), 0.0])
    [S0] = calculate(xb * 1e3, E_in, dict(BRISCOE_P, **model), outputs=['S0'])
    return _rms_percent(yb, np.real(S0))


needs_briscoe = pytest.mark.skipif(
    not BRISCOE.is_dir(), reason="Briscoe et al.'s data are not shipped")


@needs_briscoe
@pytest.mark.parametrize('fname, theta_E', BRISCOE_SPECTRA)
def test_briscoe_regime_one_es4_default(fname, theta_E):
    """ES4's default propagation reproduces each polarisation to 1 % RMS."""
    assert _briscoe_rms(fname, theta_E) <= 1.0


@needs_briscoe
@pytest.mark.parametrize('fname, theta_E', BRISCOE_SPECTRA)
def test_briscoe_regime_one_es3_method_misses(fname, theta_E):
    """ES3's orthogonal projection misses each polarisation by over 3 % RMS:
    the data tell the two methods apart."""
    assert _briscoe_rms(fname, theta_E, use_non_orthogonal=False) > 3.0


# --- Alqarni et al. 2025, Fig. 3 --------------------------------------------

TC = 370.0 - 273.15     # cell body, deg C: Doppler width and buffer-gas model
TS = 350.0 - 273.15     # stem, deg C: number density
T_FILL = 293.0          # K
# (file, Ne fill pressure in Torr, RMS bound in %, Table 1 of Alqarni et al.:
#  (Gamma, its uncertainty), (shift, its uncertainty), in MHz)
ALQARNI_SPECTRA = [
    ('No buffer gas.csv', 0.0, 0.8, None, None),
    ('60 Torr.csv', 57.4, 2.1, (417.0, 3.0), (-87.0, 2.0)),
    ('100 torr.csv', 96.0, 1.8, (700.0, 3.0), (-144.0, 2.0)),
]

needs_alqarni = pytest.mark.skipif(
    not ALQARNI.is_dir(), reason="Alqarni et al.'s data are not shipped")


def _alqarni_p(P):
    p = dict(Elem='K', Dline='D1', lcell=75e-3, T=TS, DoppTemp=TC,
             Constrain=False, Bfield=0.0)
    if P > 0:
        p.update(BufferGas='Ne', BufferPressure=P, BufferTemp_fill=T_FILL)
    return p


@needs_alqarni
@pytest.mark.parametrize('fname, P, bound, gamma, shift', ALQARNI_SPECTRA)
def test_alqarni_fig3_spectra(fname, P, bound, gamma, shift):
    """The forward calculation reproduces each spectrum within its bound."""
    d = np.loadtxt(ALQARNI / fname, delimiter=',', skiprows=1, usecols=(1, 2))
    x, y = d[:, 0], d[:, 1]
    [S0] = calculate(x * 1e3, [1, 0, 0], _alqarni_p(P), outputs=['S0'])
    assert _rms_percent(y, np.real(S0)) <= bound


@needs_alqarni
@pytest.mark.parametrize('fname, P, bound, gamma, shift',
                         [s for s in ALQARNI_SPECTRA if s[1] > 0])
def test_alqarni_predicted_broadening_and_shift(fname, P, bound, gamma, shift):
    """The model's broadening and shift at the body temperature lie within
    Alqarni et al.'s fitted values and their uncertainties."""
    g, s = get_buffer_broadening_and_shift('K', 'D1', 'Ne', P, T_FILL,
                                           TC + 273.15)
    assert abs(g - gamma[0]) <= gamma[1]
    assert abs(s - shift[0]) <= shift[1]
