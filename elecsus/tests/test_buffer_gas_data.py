"""
The buffer-gas database against the printed tables of its sources.

The expected numbers are in literature.BUFFER_GAS_MEASUREMENTS, transcribed
from the papers. Each entry must be the measurement made nearest to T_ref
(328 K), converted with (T_meas / T_ref)^0.5 and rounded to 0.01 MHz/Torr; for
an entry with several rows the broadening and the shift come from the same
row. See the header of BufferGasData.py for the rule, and for the sources.

Before this check, eight entries differed: K D1/Kr (the shift,
-5.24 for -5.42, a misprint of Pitz 2014), Cs D1/Xe (not converted), Cs D2/He
(rounding), Rb D1/He and Ne, Rb D2/He, N2 and CH4.
"""
import numpy as np
import pytest

from elecsus.libs.BufferGasData import BUFFER_GAS_DATA, get_buffer_broadening_and_shift
from elecsus.tests import literature as lit

T_REF = lit.BUFFER_GAS_T_REF


def _chosen_row(key):
    """The measurement the rule picks: nearest to T_ref, then the stated tie-break."""
    rows = lit.BUFFER_GAS_MEASUREMENTS[key]
    nearest = min(abs(r[0] - T_REF) for r in rows)
    tied = [r for r in rows if abs(r[0] - T_REF) == nearest]
    if len(tied) > 1:
        tied = [r for r in tied if r[0] == lit.BUFFER_GAS_TIE_BREAK[key]]
    assert len(tied) == 1, key
    return tied[0]


def test_every_entry_has_a_printed_source_and_the_reverse():
    assert set(BUFFER_GAS_DATA) == set(lit.BUFFER_GAS_MEASUREMENTS)
    assert len(BUFFER_GAS_DATA) == 58


@pytest.mark.parametrize('key', sorted(BUFFER_GAS_DATA), ids=lambda k: '/'.join(k))
def test_entry_is_the_nearest_measurement_converted_to_328_K(key):
    T, gamma, shift, _ = _chosen_row(key)
    factor = (T / T_REF) ** 0.5
    entry = BUFFER_GAS_DATA[key]
    assert entry['T_ref'] == T_REF
    assert entry['n_gamma'] == 0.5 and entry['n_shift'] == 0.5
    # rounded to 0.01 MHz/Torr, so within half a unit of the last digit
    assert entry['gamma_0'] == pytest.approx(gamma * factor, abs=0.0051)
    assert entry['shift_0'] == pytest.approx(shift * factor, abs=0.0051)


def test_tie_breaks_are_only_used_where_two_rows_are_equally_near():
    for key, rows in lit.BUFFER_GAS_MEASUREMENTS.items():
        nearest = min(abs(r[0] - T_REF) for r in rows)
        tied = [r for r in rows if abs(r[0] - T_REF) == nearest]
        assert (len(tied) > 1) == (key in lit.BUFFER_GAS_TIE_BREAK), key


# --- the individual findings, by name -----------------------------------------

def test_K_D1_Kr_shift_is_minus_5_42_not_the_transposed_5_24():
    # Pitz 2012 Table 3 and abstract; Pitz 2014 Table 4 misprints it.
    assert BUFFER_GAS_DATA[('K', 'D1', 'Kr')]['shift_0'] == -5.42


def test_Cs_D1_Xe_broadening_is_converted_from_313_K_like_its_shift():
    d = BUFFER_GAS_DATA[('Cs', 'D1', 'Xe')]
    assert d['gamma_0'] == pytest.approx(19.74 * (313 / 328) ** 0.5, abs=0.0051)
    assert d['shift_0'] == pytest.approx(-6.43 * (313 / 328) ** 0.5, abs=0.0051)


def test_Cs_D1_N2_comes_from_one_row_of_the_table():
    """Both numbers from the 323 K row (16.36, -7.71): the abstract's 16.64 and
    -7.76 are no row, and mixing it with the table gave two different factors."""
    d = BUFFER_GAS_DATA[('Cs', 'D1', 'N2')]
    f = (323 / 328) ** 0.5
    assert d['gamma_0'] / 16.36 == pytest.approx(f, abs=5e-4)
    assert d['shift_0'] / -7.71 == pytest.approx(f, abs=5e-4)


def test_Rb_D2_Kr_shift_is_negative_as_Rotondaro_and_Perram_print_it():
    """R&P's Table 3 has -5.50; there was never a sign to correct."""
    d = BUFFER_GAS_DATA[('Rb', 'D2', 'Kr')]
    assert d['shift_0'] == pytest.approx(-5.50 * (394 / 328) ** 0.5, abs=0.0051)


def test_Rb_D1_He_takes_both_numbers_from_the_328_K_measurement():
    d = BUFFER_GAS_DATA[('Rb', 'D1', 'He')]
    assert (d['gamma_0'], d['shift_0']) == (20.80, 5.80)


def test_Rb_D2_helium_and_methane_follow_Zameroski_like_the_other_hydrocarbons():
    f = (314 / 328) ** 0.5
    he = BUFFER_GAS_DATA[('Rb', 'D2', 'He')]
    ch4 = BUFFER_GAS_DATA[('Rb', 'D2', 'CH4')]
    assert he['gamma_0'] == pytest.approx(20.3 * f, abs=0.0051)
    assert ch4['gamma_0'] == pytest.approx(28.0 * f, abs=0.0051)
    assert ch4['shift_0'] == pytest.approx(-8.4 * f, abs=0.0051)


# --- convention-free checks on the whole table ---------------------------------

@pytest.mark.parametrize('elem', ['K', 'Rb', 'Cs'])
@pytest.mark.parametrize('dline', ['D1', 'D2'])
def test_heavy_rare_gases_and_hydrocarbons_shift_to_the_red(elem, dline):
    """Ar, Kr, Xe and every hydrocarbon red-shift every alkali D line in every
    one of these papers; a positive entry is a lost or flipped sign."""
    for (e, l, gas), d in BUFFER_GAS_DATA.items():
        if (e, l) == (elem, dline) and gas in ('Ar', 'Kr', 'Xe', 'CH4', 'C2H6', 'C3H8', 'C4H10'):
            assert d['shift_0'] < 0, (e, l, gas)


def test_helium_is_the_only_rare_gas_with_a_blue_shift():
    for (e, l, gas), d in BUFFER_GAS_DATA.items():
        if gas == 'He':
            assert d['shift_0'] > 0, (e, l)
        if gas in ('Ne', 'Ar', 'Kr', 'Xe'):
            assert d['shift_0'] < 0, (e, l, gas)


def test_broadening_is_positive_and_of_the_published_size():
    """5 to 35 MHz/Torr for every entry (the papers' tables run 6 to 34)."""
    for key, d in BUFFER_GAS_DATA.items():
        assert 5 < d['gamma_0'] < 35, key


def test_the_rate_at_the_reference_temperature_is_the_stored_one():
    for key, d in BUFFER_GAS_DATA.items():
        gamma, shift = get_buffer_broadening_and_shift(*key, 1.0, T_REF, T_REF)
        assert np.isclose(gamma, d['gamma_0']) and np.isclose(shift, d['shift_0'])
