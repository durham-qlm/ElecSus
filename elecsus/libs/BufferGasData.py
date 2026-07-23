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

# elecsus/libs/BufferGasData.py
# [ES4 NEW] — entire file; ES3 had no buffer gas database or dynamic pressure model.
#
# Buffer gas broadening and shift reference data.
#
# All entries are stored at T_ref = 328 K (the primary measurement temperature
# of Pitz et al. 2014).  Values originally measured at other temperatures have
# been converted using the hard-sphere T^(1/2) scaling:
#
#   gamma(T_ref) = gamma(T_meas) * (T_meas / T_ref)^0.5
#
# Source: Pitz, G. A. et al., J. Quant. Spectrosc. Radiat. Transfer 140, 18-29
#         (2014).  DOI 10.1016/j.jqsrt.2014.01.024
#         Tables 2-9 (K D1/D2 measured directly at 328 K; Rb D1 rare-gas gamma
#         converted from 394 K via the gamma_T=328K column in Table 2; Rb D2
#         and Cs data converted from their respective measurement temperatures).
#         K D1 Ne additionally validated in Alqarni et al., J. Phys. B 58
#         (2025).
#
# Dictionary structure:
#   (Element, D-line, BufferGas): {
#       'gamma_0' : broadening rate (MHz/Torr) at T_ref,
#       'shift_0' : frequency shift rate (MHz/Torr) at T_ref,  [+ = blue shift]
#       'n_gamma' : broadening temperature exponent (see model below),
#       'n_shift' : shift temperature exponent,
#       'T_ref'   : reference temperature (K)
#   }
#
# Model for a sealed cell (fixed number of buffer gas atoms):
#
#   P(T) = P_fill * (T / T_fill)                       [ideal gas law]
#
#   gamma_rate(T) = gamma_0 * (T_ref / T)^(1 - n)     [collision cross-section
#                                                        scales as T^(-0.5)
#                                                        in hard-sphere limit]
#
#   Gamma(T) = gamma_rate(T) * P(T)
#            = gamma_0 * P_fill * (T_ref / T_fill) * (T / T_ref)^n
#
# with n = n_gamma = 0.5 (hard sphere).  Same formula applies for shift.
#
# In spectra.py:
#   P_eq   = P_fill * (T_ref / T_fill)   [effective pressure at T_ref]
#   Gamma  = gamma_0 * P_eq * (T / T_ref)^n_gamma
#   Shift  = shift_0 * P_eq * (T / T_ref)^n_shift
#
# Na is NOT in the database — no consistent modern dataset exists at present.
# An unrecognised (element, dline, gas) key returns (None, None) from
# get_buffer_broadening_and_shift(), which is treated as "model not available"
# in spectra.py (a warning is issued and the static GammaBuf/shift is used).
#
# Gas name keys (use these strings in the parameter dictionary):
#   noble gases  : 'He', 'Ne', 'Ar', 'Kr', 'Xe'
#   molecular    : 'N2', 'H2', 'CH4', 'C2H6', 'C3H8', 'C4H10'

BUFFER_GAS_DATA = {

    # =========================================================================
    # Potassium D1  (4^2P_{1/2} <- 4^2S_{1/2}, 770 nm)
    # All measured at 328 K — Pitz 2014 Tables 2-5 (this work, citing Pitz 2012)
    # =========================================================================
    ('K', 'D1', 'He'):    {'gamma_0': 13.08, 'shift_0':  1.63, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'Ne'):    {'gamma_0':  6.14, 'shift_0': -1.27, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'Ar'):    {'gamma_0': 19.45, 'shift_0': -6.44, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'Kr'):    {'gamma_0': 16.64, 'shift_0': -5.24, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'Xe'):    {'gamma_0': 20.02, 'shift_0': -6.54, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'H2'):    {'gamma_0': 22.15, 'shift_0': -5.34, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'N2'):    {'gamma_0': 17.78, 'shift_0': -6.80, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'CH4'):   {'gamma_0': 29.35, 'shift_0': -7.41, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'C2H6'):  {'gamma_0': 26.63, 'shift_0': -8.32, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'C3H8'):  {'gamma_0': 27.27, 'shift_0': -8.59, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'C4H10'): {'gamma_0': 27.85, 'shift_0': -8.80, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},

    # =========================================================================
    # Potassium D2  (4^2P_{3/2} <- 4^2S_{1/2}, 766 nm)
    # All measured at 328 K — Pitz 2014 Tables 6-9 (this work)
    # =========================================================================
    ('K', 'D2', 'He'):    {'gamma_0': 19.84, 'shift_0':  0.52, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'Ne'):    {'gamma_0':  8.88, 'shift_0': -2.06, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'Ar'):    {'gamma_0': 18.65, 'shift_0': -5.52, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'Kr'):    {'gamma_0': 19.17, 'shift_0': -5.42, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'Xe'):    {'gamma_0': 22.19, 'shift_0': -7.01, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'N2'):    {'gamma_0': 18.98, 'shift_0': -5.66, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'CH4'):   {'gamma_0': 27.78, 'shift_0': -8.38, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'C2H6'):  {'gamma_0': 27.60, 'shift_0': -8.04, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'C3H8'):  {'gamma_0': 27.70, 'shift_0': -9.22, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D2', 'C4H10'): {'gamma_0': 33.48, 'shift_0': -9.37, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},

    # =========================================================================
    # Rubidium D1  (5^2P_{1/2} <- 5^2S_{1/2}, 795 nm)
    # Rare-gas gamma: Pitz 2014 Table 2, gamma_T=328K column (converted from
    #   394 K Rotondaro & Perram [6] via (394/328)^0.5 = 1.096).
    # Rare-gas shift: converted from Rotondaro & Perram [6] 394 K values
    #   except He where Pitz 2014 "this work" gives 5.80 directly at 328 K.
    # Molecular (CH4, C2H6, C3H8, C4H10): Pitz 2014 "this work" at 328 K
    #   (Tables 3, 5).
    # H2, N2 shifts: converted from Rotondaro & Perram 394 K values.
    # =========================================================================
    ('Rb', 'D1', 'He'):    {'gamma_0': 20.71, 'shift_0':  5.80, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'Ne'):    {'gamma_0': 10.75, 'shift_0': -0.99, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'Ar'):    {'gamma_0': 19.84, 'shift_0': -7.42, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'Kr'):    {'gamma_0': 18.74, 'shift_0': -5.61, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'Xe'):    {'gamma_0': 20.39, 'shift_0': -6.40, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'H2'):    {'gamma_0': 18.96, 'shift_0': -2.38, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'N2'):    {'gamma_0': 17.86, 'shift_0': -8.12, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'CH4'):   {'gamma_0': 32.78, 'shift_0': -6.96, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'C2H6'):  {'gamma_0': 30.49, 'shift_0': -7.88, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'C3H8'):  {'gamma_0': 33.05, 'shift_0': -8.61, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'C4H10'): {'gamma_0': 29.61, 'shift_0': -9.43, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},

    # =========================================================================
    # Rubidium D2  (5^2P_{3/2} <- 5^2S_{1/2}, 780 nm)
    # Rare-gas gamma: Pitz 2014 Table 6 gamma_T=328K column.
    # Rare-gas shift: converted from Rotondaro & Perram [6] 394 K values
    #   via (394/328)^0.5 = 1.096.
    # Molecular: mix of [6] (394 K → converted) and [5] Zameroski (314 K → converted).
    # NOTE — Rb D2 Kr shift: Rotondaro & Perram [6] reports +5.50 MHz/Torr at
    #   394 K (anomalously positive; likely a sign error in [6]).  The value
    #   below is from McCartan & Lwin [17] at 400-420 K, converted to 328 K.
    #   Treat with caution until independently verified.
    # =========================================================================
    ('Rb', 'D2', 'He'):    {'gamma_0': 21.92, 'shift_0':  0.40, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'Ne'):    {'gamma_0': 10.38, 'shift_0': -2.67, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'Ar'):    {'gamma_0': 19.40, 'shift_0': -6.31, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'Kr'):    {'gamma_0': 18.85, 'shift_0': -4.86, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'Xe'):    {'gamma_0': 19.51, 'shift_0': -6.78, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'N2'):    {'gamma_0': 20.06, 'shift_0': -6.34, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'CH4'):   {'gamma_0': 28.72, 'shift_0': -7.67, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'C2H6'):  {'gamma_0': 27.49, 'shift_0': -8.61, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'C3H8'):  {'gamma_0': 29.84, 'shift_0': -9.49, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'C4H10'): {'gamma_0': 30.62, 'shift_0': -9.78, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},

    # =========================================================================
    # Caesium D1  (6^2P_{1/2} <- 6^2S_{1/2}, 894 nm)
    # From Pitz 2009 [15] and Pitz 2014 "this work", reported in Tables 2-5.
    # He gamma: "this work" 328 K directly.  He shift: converted from 323 K.
    # Noble gas data mostly from [15] at 313 K → converted to 328 K.
    # Molecular data from [15] at 313-333 K → converted.
    # Note: Cs D1/H2 shift is +1.11 MHz/Torr (blue shift — correct per [15]).
    # =========================================================================
    ('Cs', 'D1', 'He'):   {'gamma_0': 20.80, 'shift_0':  4.21, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'Ne'):   {'gamma_0': 10.60, 'shift_0': -1.56, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'Ar'):   {'gamma_0': 17.89, 'shift_0': -6.32, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'Kr'):   {'gamma_0': 17.41, 'shift_0': -5.33, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'Xe'):   {'gamma_0': 19.74, 'shift_0': -6.28, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'H2'):   {'gamma_0': 20.81, 'shift_0':  1.11, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'N2'):   {'gamma_0': 16.23, 'shift_0': -7.65, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'CH4'):  {'gamma_0': 29.22, 'shift_0': -9.35, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},

    # =========================================================================
    # Caesium D2  (6^2P_{3/2} <- 6^2S_{1/2}, 852 nm)
    # From Pitz 2009 [16], reported in Pitz 2014 Tables 6-9.
    # All originally at 313 K, converted to 328 K via (313/328)^0.5 = 0.9770.
    # =========================================================================
    ('Cs', 'D2', 'He'):   {'gamma_0': 21.83, 'shift_0':  0.67, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D2', 'Ne'):   {'gamma_0':  9.58, 'shift_0': -2.52, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D2', 'Ar'):   {'gamma_0': 16.09, 'shift_0': -6.04, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D2', 'Kr'):   {'gamma_0': 15.18, 'shift_0': -5.95, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D2', 'Xe'):   {'gamma_0': 17.98, 'shift_0': -6.59, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D2', 'N2'):   {'gamma_0': 18.74, 'shift_0': -6.06, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D2', 'CH4'):  {'gamma_0': 25.24, 'shift_0': -8.66, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D2', 'C2H6'): {'gamma_0': 25.54, 'shift_0': -9.16, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
}


def get_buffer_broadening_and_shift(element, dline, buffer_gas, P_fill_torr,
                                    T_fill_K, T_K):
    """
    Calculate the temperature-dependent Lorentzian broadening and shift
    imparted by a buffer gas in a sealed cell.

    Parameters
    ----------
    element : str
        'Rb', 'K', or 'Cs'.  Na is not supported.
    dline : str
        'D1' or 'D2'.
    buffer_gas : str
        Gas species key, e.g. 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'N2', 'CH4'.
    P_fill_torr : float
        Fill pressure of the buffer gas (Torr) at temperature T_fill_K.
    T_fill_K : float
        Temperature at which the cell was filled (K).
    T_K : float
        Current operating temperature of the vapour (K).

    Returns
    -------
    (gamma_MHz, shift_MHz) : (float, float) or (None, None)
        Lorentzian FWHM broadening (MHz) and line-centre shift (MHz).
        Returns (None, None) if the (element, dline, buffer_gas) combination
        is not in the database.
    """
    key = (element, dline, buffer_gas)
    if key not in BUFFER_GAS_DATA:
        return None, None

    d = BUFFER_GAS_DATA[key]
    T_ref = d['T_ref']

    # Effective pressure at T_ref (number-density conserved in sealed cell)
    P_eq = P_fill_torr * (T_ref / T_fill_K)

    gamma = d['gamma_0'] * P_eq * (T_K / T_ref) ** d['n_gamma']
    shift = d['shift_0'] * P_eq * (T_K / T_ref) ** d['n_shift']

    return gamma, shift
