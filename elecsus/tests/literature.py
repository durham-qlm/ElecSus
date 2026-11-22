"""
Literature values for the convention-free physics tests.

Every expected value in test_sum_rule.py, test_kramers_kronig.py,
test_hyperfine_levels.py and test_symmetries.py comes from here or from
CODATA (scipy.constants), never from ElecSus's own constants or output. The
tests import nothing from elecsus.libs.AtomConstants for that reason: a typo
there must show up as a disagreement, not be copied into the expectation.

Sources
-------
[Steck]  D. A. Steck, "Rubidium 85 D Line Data", "Rubidium 87 D Line Data",
         "Cesium D Line Data", "Sodium D Line Data", http://steck.us/alkalidata
         (revisions 2.x).
[Tiecke] T. G. Tiecke, "Properties of Potassium", University of Amsterdam
         (2010, revised 2019).
[DN06]   D. Das and V. Natarajan, Eur. Phys. J. D 37, 313 (2006): Rb 5P1/2.
[P15]    B. M. Patterson et al., Phys. Rev. A 91, 012506 (2015): Cs 6P3/2
         lifetime.

Where measurements disagree by more than the tests' tolerance needs, the
spread is noted next to the value and the tests are written to tolerate it.
"""
from scipy.constants import physical_constants

# The free-electron g-factor, positive (CODATA). ElecSus's Hamiltonian uses
# it for the electron spin (FundamentalConstants.gs), so the Breit-Rabi
# formula with it is the model's own first-principles prediction.
G_S = -physical_constants['electron g factor'][0]

# Bohr magneton in MHz per gauss.
MU_B_MHZ_PER_G = physical_constants['Bohr magneton in Hz/T'][0] * 1e-10


# --- isotopes: nuclear spin I, ground-state A(S1/2) [MHz], measured gJ(S1/2),
# gI (sign convention of Steck/Tiecke: H_Z = muB (gJ Jz + gI Iz) B) ----------
ISOTOPES = {
    'Rb85': dict(I=2.5, A_S=1011.910813, gJ=2.00233113, gI=-0.00029364000),
    'Rb87': dict(I=1.5, A_S=3417.341305452145, gJ=2.00233113, gI=-0.0009951414),
    'Cs':   dict(I=3.5, A_S=2298.1579425, gJ=2.00254032, gI=-0.00039885395),
    'Na':   dict(I=1.5, A_S=885.8130644, gJ=2.0022960, gI=-0.00080461080),
    'K39':  dict(I=1.5, A_S=230.8598601, gJ=2.00229421, gI=-0.00014193489),
    'K40':  dict(I=4.0, A_S=-285.7308, gJ=2.00229421, gI=0.000176490),
    'K41':  dict(I=1.5, A_S=127.0069352, gJ=2.00229421, gI=-0.00007790600),
}

# --- excited-state hyperfine constants [MHz]: A(P1/2); A and B of P3/2 ------
# [Steck] for Rb, Cs, Na; [Tiecke] for K.
# Spreads between published values (these set the tests' 1 % tolerance):
#   Rb87 A(P1/2): 406.119(7) [DN06] to 408.328(15) [Steck], 0.5 %;
#   Rb85 A(P1/2): 120.527(56) [Steck] to 120.645(5) [DN06], 0.1 %;
#   Rb85 P3/2: A = 25.0020(99), B = 25.790(93) [Steck]; later saturated-
#     absorption values (A = 25.038, B = 26.011) are 0.15 % and 0.9 % higher.
EXCITED = {
    'Rb85': dict(A_P12=120.527, A_P32=25.0020, B_P32=25.790),
    'Rb87': dict(A_P12=408.328, A_P32=84.7185, B_P32=12.4965),
    'Cs':   dict(A_P12=291.9201, A_P32=50.28827, B_P32=-0.4934),
    'Na':   dict(A_P12=94.44, A_P32=18.534, B_P32=2.724),
    'K39':  dict(A_P12=27.775, A_P32=6.093, B_P32=2.786),
    'K40':  dict(A_P12=-34.523, A_P32=-7.585, B_P32=-3.445),
    'K41':  dict(A_P12=15.245, A_P32=3.363, B_P32=3.351),
}

# --- D lines: vacuum wavelength [m] and excited-state lifetime [s] ----------
# [Steck] for Rb, Cs D1 and Na; [P15] for Cs D2; [Tiecke] for K. Measured
# lifetimes of one level differ by up to ~0.5 %, which with the 1 % tolerance
# of the sum-rule tests leaves room for the choice of measurement.
LINES = {
    ('Rb', 'D1'): dict(wavelength=794.978851156e-9, tau=27.679e-9),
    ('Rb', 'D2'): dict(wavelength=780.241209686e-9, tau=26.2348e-9),
    ('Cs', 'D1'): dict(wavelength=894.59295986e-9, tau=34.791e-9),
    ('Cs', 'D2'): dict(wavelength=852.34727582e-9, tau=30.462e-9),
    ('Na', 'D1'): dict(wavelength=589.7558147e-9, tau=16.299e-9),
    ('Na', 'D2'): dict(wavelength=589.1583264e-9, tau=16.2492e-9),
    ('K', 'D1'):  dict(wavelength=770.108385049e-9, tau=26.72e-9),
    ('K', 'D2'):  dict(wavelength=766.700921822e-9, tau=26.37e-9),
}

# Statistical-weight ratio g(J')/g(J) of the upper and lower levels.
G_RATIO = {'D1': 1.0, 'D2': 2.0}

# --- saturated vapour pressure at 25 C [torr] ---------------------------------
# [Steck]: Rb 3.92(20)e-7 torr; Cs 1.488e-6 torr (quoted to ~5 %, the
# accuracy of the vapour-pressure models these come from).
VAPOUR_PRESSURE_25C_TORR = {'Rb': 3.92e-7, 'Cs': 1.488e-6}


def oscillator_strength(elem, dline):
    """Absorption oscillator strength from the lifetime, the standard relation
    A = (2 pi e^2 / (eps0 m_e c lambda^2)) (g/g') f (e.g. Steck, Sec. 2.3;
    Hilborn, Am. J. Phys. 50, 982 (1982))."""
    from scipy.constants import e, epsilon_0, m_e, c, pi
    line = LINES[(elem, dline)]
    A = 1 / line['tau']
    return (G_RATIO[dline] * A * epsilon_0 * m_e * c * line['wavelength'] ** 2
            / (2 * pi * e ** 2))


# --- buffer-gas broadening and shift rates, as printed -----------
# Every measurement of an entry of the database that the papers below give,
# transcribed from their tables: (T_meas [K], broadening gamma [MHz/Torr],
# shift delta [MHz/Torr], source). Used by test_buffer_gas_data.py. Not copied
# from BufferGasData.py: a typo there must show up as a disagreement.
#
# [P12] Pitz, Sandoval, Zameroski, Klennert, Hostutler, JQSRT 113, 387 (2012):
#       K D1, Tables 2 and 3 (and the abstract), 328 K.
# [P14] Pitz, Sandoval, Tafoya, Klennert, Hostutler, JQSRT 140, 18 (2014):
#       K D2, Tables 6-9; Rb D1 "this work", Tables 2-5; all 328 K. (Its Table 4
#       prints the K D1/Kr shift as -5.24 and its Table 8 the Rb D2/Kr shift
#       of [RP] as 5.50: neither is used.)
# [RP]  Rotondaro and Perram, JQSRT 57, 497 (1997): Rb D1 (Tables 1, 2) and
#       D2 (Tables 3, 4), 394 K. Table 3 prints the Kr shift as -5.50.
# [Z]   Zameroski, Hager, Rudolph, Erickson, Hostutler, JQSRT 112, 59 (2011):
#       Rb D2, Tables 1 and 2, 314 K.
# [P09] Pitz, Wertepny, Perram, PRA 80, 062718 (2009): Cs D1, Tables II, III.
#       (Its abstract gives N2 as 16.64 and -7.76, which are no row of them.)
# [P10] Pitz, Fox, Perram, PRA 82, 042502 (2010): Cs D2, Tables II, III, 313 K.
BUFFER_GAS_T_REF = 328.0


def _rows(elem, dline, T, ref, table):
    return {(elem, dline, gas): [(T, g, d, ref)] for gas, g, d in table}


BUFFER_GAS_MEASUREMENTS = {}
BUFFER_GAS_MEASUREMENTS.update(_rows('K', 'D1', 328., 'P12', [
    ('He', 13.08, 1.63), ('Ne', 6.14, -1.27), ('Ar', 19.45, -6.44),
    ('Kr', 16.64, -5.42), ('Xe', 20.02, -6.54), ('H2', 22.15, -5.34),
    ('N2', 17.78, -6.80), ('CH4', 29.35, -7.41), ('C2H6', 26.63, -8.32),
    ('C3H8', 27.27, -8.59), ('C4H10', 27.85, -8.80)]))
BUFFER_GAS_MEASUREMENTS.update(_rows('K', 'D2', 328., 'P14', [
    ('He', 19.84, 0.52), ('Ne', 8.88, -2.06), ('Ar', 18.65, -5.52),
    ('Kr', 19.17, -5.42), ('Xe', 22.19, -7.01), ('N2', 18.98, -5.66),
    ('CH4', 27.78, -8.38), ('C2H6', 27.60, -8.04), ('C3H8', 27.70, -9.22),
    ('C4H10', 33.48, -9.37)]))
BUFFER_GAS_MEASUREMENTS.update(_rows('Rb', 'D1', 394., 'RP', [
    ('He', 18.9, 4.71), ('Ne', 9.84, -0.90), ('Ar', 18.1, -6.77),
    ('Kr', 17.1, -5.12), ('Xe', 18.6, -5.84), ('H2', 17.3, -2.17),
    ('N2', 16.3, -7.41), ('CH4', 29.1, -7.92)]))
for _gas, _g, _d in [('He', 20.80, 5.80), ('CH4', 32.78, -6.96),
                     ('C2H6', 30.49, -7.88), ('C3H8', 33.05, -8.61),
                     ('C4H10', 29.61, -9.43)]:
    BUFFER_GAS_MEASUREMENTS.setdefault(('Rb', 'D1', _gas), []).append(
        (328., _g, _d, 'P14'))
BUFFER_GAS_MEASUREMENTS.update(_rows('Rb', 'D2', 394., 'RP', [
    ('He', 20.0, 0.37), ('Ne', 9.47, -2.44), ('Ar', 17.7, -5.76),
    ('Kr', 17.2, -5.50), ('Xe', 17.8, -6.19), ('N2', 18.3, -5.79),
    ('CH4', 26.2, -7.00)]))
for _gas, _g, _d in [('He', 20.3, 0.39), ('CH4', 28.0, -8.4),
                     ('C2H6', 28.1, -8.8), ('C3H8', 30.5, -9.7),
                     ('C4H10', 31.3, -10.0)]:
    BUFFER_GAS_MEASUREMENTS.setdefault(('Rb', 'D2', _gas), []).append(
        (314., _g, _d, 'Z'))
BUFFER_GAS_MEASUREMENTS.update({
    ('Cs', 'D1', 'He'): [(323., 24.13, 4.24, 'P09')],
    ('Cs', 'D1', 'Ne'): [(313., 10.85, -1.60, 'P09')],
    ('Cs', 'D1', 'Ar'): [(313., 18.31, -6.47, 'P09')],
    ('Cs', 'D1', 'Kr'): [(313., 17.82, -5.46, 'P09')],
    ('Cs', 'D1', 'Xe'): [(313., 19.74, -6.43, 'P09')],
    ('Cs', 'D1', 'H2'): [(328., 20.81, 1.11, 'P09')],
    # Three temperatures; 323 K and 333 K are equally far from 328 K. The
    # 323 K slope errors are 0.02 and 0.01; the 333 K ones 0.08 and 0.01.
    ('Cs', 'D1', 'N2'): [(318., 15.82, -7.69, 'P09'),
                         (323., 16.36, -7.71, 'P09'),
                         (333., 15.66, -7.41, 'P09')],
    ('Cs', 'D1', 'CH4'): [(333., 29.00, -9.28, 'P09')],
})
BUFFER_GAS_MEASUREMENTS.update(_rows('Cs', 'D2', 313., 'P10', [
    ('He', 20.59, 0.69), ('Ne', 9.81, -2.58), ('Ar', 16.47, -6.18),
    ('Kr', 15.54, -6.09), ('Xe', 18.41, -6.75), ('N2', 19.18, -6.20),
    ('CH4', 25.84, -8.86), ('C2H6', 26.14, -9.38)]))

# The quoted slope error of each row where two rows are equally near T_ref.
BUFFER_GAS_TIE_BREAK = {('Cs', 'D1', 'N2'): 323.}
