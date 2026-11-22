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
# All entries are stored at T_ref = 328 K (the measurement temperature of
# the potassium data of Pitz et al. 2012 and 2014). A value measured at another
# temperature T_meas is converted with the hard-sphere T^(1/2) scaling of the
# collision rate per Torr (the collision velocity grows as T^(1/2), and the
# gas density at a given pressure falls as 1/T):
#
#   gamma(T_ref) = gamma(T_meas) * (T_meas / T_ref)^0.5     (likewise the shift)
#
# and rounded to 0.01 MHz/Torr. The conversion is the one the papers use for
# their own comparisons (e.g. Pitz 2014, Eq. 6; Pitz 2009, Eq. 6).
#
# Sources (every number below was checked against the printed tables):
#   [P12] G. A. Pitz et al., J. Quant. Spectrosc. Radiat. Transfer 113, 387
#         (2012). K D1, Tables 2 and 3, measured at 328 K.
#   [P14] G. A. Pitz et al., J. Quant. Spectrosc. Radiat. Transfer 140, 18
#         (2014). K D2 (Tables 6-9, measured at 328 K) and the Rb D1 "this work"
#         rows (He, CH4, C2H6, C3H8, C4H10; Tables 2-5, 328 K). Its Table 4
#         prints the K D1/Kr shift as -5.24; [P12] has -5.42 (Table 3 and the
#         abstract), which is used.
#   [RP]  M. D. Rotondaro and G. P. Perram, J. Quant. Spectrosc. Radiat.
#         Transfer 57, 497 (1997). Rb D1 and D2, rare gases, H2, N2, CH4, at
#         394 K (Tables 1-4).
#   [Z]   N. D. Zameroski et al., J. Quant. Spectrosc. Radiat. Transfer 112, 59
#         (2011). Rb D2, He and the four hydrocarbons, at 314 K (Tables 1-2).
#   [P09] G. A. Pitz, D. E. Wertepny and G. P. Perram, Phys. Rev. A 80, 062718
#         (2009). Cs D1, Tables II and III, measured at 313-333 K.
#   [P10] G. A. Pitz, C. D. Fox and G. P. Perram, Phys. Rev. A 82, 042502
#         (2010). Cs D2, Tables II and III, measured at 313 K.
#   (Lwin & McCartan, J. Phys. B 11, 3841 (1978), measured K at 400-420 K; it
#   is not a source for any entry here: the K data are those of [P12], [P14].)
#
# Choosing between sources. Where more than one measurement of an entry exists
# in these papers, the one made nearest to T_ref is used, because the
# conversion is then the smallest and its assumed exponent matters least; a tie
# goes to the smaller quoted uncertainty. Both the broadening and the shift of
# an entry come from the same row. So: Rb D1/He is the 328 K measurement of
# [P14] (not [RP] at 394 K); Rb D2/He and Rb D2/CH4 are [Z] at 314 K (not [RP]
# at 394 K), which also puts CH4 with the other three hydrocarbons; Cs D1/N2
# is the 323 K row of [P09] (318, 323 and 333 K were measured; 323 and 333 K are
# equally near, and the 323 K slope is four times better determined).
# elecsus/tests/test_buffer_gas_data.py repeats every entry from the printed
# values (elecsus/tests/literature.py) and applies this rule.
#
# K D1/Ne is also validated against Alqarni et al., J. Phys. B 58, 135201
# (2025) (elecsus/tests/test_buffer_gas.py).
#
# Known limits of the data:
#   * The exponent 0.5 is checked for the broadening only (Pitz 2009, Fig. 6:
#     the Cs-N2 broadening cross-section is the same at 318, 323 and 333 K to
#     a few per cent). For the shift the papers say it may differ (Pitz 2009,
#     sec. IV; Pitz 2014, sec. 4): Rb D1/He is 5.80 measured at 328 K [P14]
#     but 4.71 x 1.096 = 5.16 from the 394 K value of [RP]: 12 % apart. The
#     shift scaled away from the temperatures measured (313-394 K) is the
#     least certain output of the model.
#   * Rb D2/Kr shift: [RP] print -5.50 (negative, as every other Rb D2 shift
#     in their table). [P14]'s Table 8 prints "5.50", with the minus sign
#     lost; the entry here is -5.50 x 1.096 = -6.03. There was never a sign
#     error to correct.
#   * Cs D1/N2: the abstract of [P09] gives 16.64 and -7.76, which are no row
#     of its Tables II and III (318 K: 15.82, -7.69; 323 K: 16.36, -7.71; 333 K:
#     15.66, -7.41); the table is used.
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
#   gamma_rate(T) = gamma_0 * (T_ref / T)^(1 - n)     [rate per Torr; it falls
#                                                        as T^(-0.5) when the
#                                                        cross-section does not
#                                                        depend on T (hard
#                                                        spheres)]
#
#   Gamma(T) = gamma_rate(T) * P(T)
#            = gamma_0 * P_fill * (T_ref / T_fill) * (T / T_ref)^n
#
# with n = n_gamma = 0.5 (hard sphere).  Same formula applies for shift (where
# 0.5 is an assumption, see "Known limits" above).
#
# In spectra.py:
#   P_eq   = P_fill * (T_ref / T_fill)   [effective pressure at T_ref]
#   Gamma  = gamma_0 * P_eq * (T / T_ref)^n_gamma
#   Shift  = shift_0 * P_eq * (T / T_ref)^n_shift
# with T the temperature of the gas in the beam, in kelvin: p_dict['T'], or
# p_dict['DoppTemp'] when Constrain=False gives the cell body its own
# temperature. Gamma replaces GammaBuf; Shift is added to
# p_dict['shift'].
#
# The result depends on P_fill and T_fill only through P_fill / T_fill (the
# density of the gas), so a fit can find one of them, not both: fit_data
# varies BufferPressure and refuses BufferTemp_fill.
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
    # All measured at 328 K: Pitz 2012 [P12], Tables 2 and 3. (Pitz 2014
    # repeats them in its Tables 2-5, but prints the Kr shift as -5.24.)
    # =========================================================================
    ('K', 'D1', 'He'):    {'gamma_0': 13.08, 'shift_0':  1.63, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'Ne'):    {'gamma_0':  6.14, 'shift_0': -1.27, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'Ar'):    {'gamma_0': 19.45, 'shift_0': -6.44, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'Kr'):    {'gamma_0': 16.64, 'shift_0': -5.42, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'Xe'):    {'gamma_0': 20.02, 'shift_0': -6.54, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'H2'):    {'gamma_0': 22.15, 'shift_0': -5.34, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'N2'):    {'gamma_0': 17.78, 'shift_0': -6.80, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'CH4'):   {'gamma_0': 29.35, 'shift_0': -7.41, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'C2H6'):  {'gamma_0': 26.63, 'shift_0': -8.32, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'C3H8'):  {'gamma_0': 27.27, 'shift_0': -8.59, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('K', 'D1', 'C4H10'): {'gamma_0': 27.85, 'shift_0': -8.80, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},

    # =========================================================================
    # Potassium D2  (4^2P_{3/2} <- 4^2S_{1/2}, 766 nm)
    # All measured at 328 K: Pitz 2014 [P14], Tables 6-9.
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
    # He, CH4, C2H6, C3H8, C4H10: Pitz 2014 [P14] "this work", measured at
    #   328 K, no conversion (Tables 2-5).
    # Ne, Ar, Kr, Xe, H2, N2: Rotondaro & Perram [RP] at 394 K (Tables 1, 2),
    #   times (394/328)^0.5 = 1.096.
    # =========================================================================
    ('Rb', 'D1', 'He'):    {'gamma_0': 20.80, 'shift_0':  5.80, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D1', 'Ne'):    {'gamma_0': 10.78, 'shift_0': -0.99, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
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
    # He, CH4, C2H6, C3H8, C4H10: Zameroski [Z] at 314 K (Tables 1, 2), times
    #   (314/328)^0.5 = 0.978. [RP] also measured He and CH4, at 394 K; the
    #   314 K row is nearer to T_ref and puts CH4 with the other hydrocarbons
    #   (the two differ by 10 % for the He width, 5 % for the CH4 width and
    #   7 % for the CH4 shift).
    # Ne, Ar, Kr, Xe, N2: [RP] at 394 K (Tables 3, 4), times
    #   (394/328)^0.5 = 1.096.
    # Kr shift: [RP]'s -5.50, a negative number as printed (see the header).
    # =========================================================================
    ('Rb', 'D2', 'He'):    {'gamma_0': 19.86, 'shift_0':  0.38, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'Ne'):    {'gamma_0': 10.38, 'shift_0': -2.67, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'Ar'):    {'gamma_0': 19.40, 'shift_0': -6.31, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'Kr'):    {'gamma_0': 18.85, 'shift_0': -6.03, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'Xe'):    {'gamma_0': 19.51, 'shift_0': -6.78, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'N2'):    {'gamma_0': 20.06, 'shift_0': -6.35, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'CH4'):   {'gamma_0': 27.40, 'shift_0': -8.22, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'C2H6'):  {'gamma_0': 27.49, 'shift_0': -8.61, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'C3H8'):  {'gamma_0': 29.84, 'shift_0': -9.49, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Rb', 'D2', 'C4H10'): {'gamma_0': 30.62, 'shift_0': -9.78, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},

    # =========================================================================
    # Caesium D1  (6^2P_{1/2} <- 6^2S_{1/2}, 894 nm)
    # Pitz 2009 [P09], Tables II and III, converted to 328 K with
    #   (T_meas/328)^0.5: He 323 K, Ne/Ar/Kr/Xe 313 K, H2 328 K (no
    #   conversion), N2 323 K (the middle of three rows, see the header),
    #   CH4 333 K.
    # Cs D1/H2 shift is +1.11 MHz/Torr (a blue shift, as in [P09]).
    # =========================================================================
    ('Cs', 'D1', 'He'):   {'gamma_0': 23.95, 'shift_0':  4.21, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'Ne'):   {'gamma_0': 10.60, 'shift_0': -1.56, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'Ar'):   {'gamma_0': 17.89, 'shift_0': -6.32, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'Kr'):   {'gamma_0': 17.41, 'shift_0': -5.33, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'Xe'):   {'gamma_0': 19.28, 'shift_0': -6.28, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'H2'):   {'gamma_0': 20.81, 'shift_0':  1.11, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'N2'):   {'gamma_0': 16.23, 'shift_0': -7.65, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
    ('Cs', 'D1', 'CH4'):  {'gamma_0': 29.22, 'shift_0': -9.35, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},

    # =========================================================================
    # Caesium D2  (6^2P_{3/2} <- 6^2S_{1/2}, 852 nm)
    # Pitz 2010 [P10], Tables II and III. All measured at 313 K, converted to
    #   328 K with (313/328)^0.5 = 0.9769.
    # He = 4He (gamma 20.59 -> 20.11); do NOT use the 3He value 22.35.
    # =========================================================================
    ('Cs', 'D2', 'He'):   {'gamma_0': 20.11, 'shift_0':  0.67, 'n_gamma': 0.5, 'n_shift': 0.5, 'T_ref': 328.0},
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
        Temperature (K) at which to evaluate them. spectra.calc_chi passes
        the temperature of the gas the light passes through, where the
        collisions happen: T, or DoppTemp when Constrain=False gives the
        cell body its own temperature.

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
