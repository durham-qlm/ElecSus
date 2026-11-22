"""
Test script for assessing the new BufferGas functionality in ElecSus.
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate
from elecsus.libs.BufferGasData import get_buffer_broadening_and_shift

def test_buffer_gas_data_extraction():
    # Test K D1 Ne — the one set of confirmed values (Pitz 2012, validated Alqarni 2025).
    # At T = T_fill = T_ref, P_eq = P_fill exactly, so Gamma = gamma_0 * P_fill.
    elem = 'K'
    dline = 'D1'
    bg = 'Ne'

    gamma_0 = 6.14   # MHz/Torr  (Pitz et al. 2012, Table 1)
    shift_0 = -1.27  # MHz/Torr
    T_ref   = 328.0  # K
    P_fill  = 100.0  # Torr
    T_fill  = T_ref  # fill at T_ref => P_eq = P_fill exactly

    gamma, shift = get_buffer_broadening_and_shift(elem, dline, bg, P_fill, T_fill, T_ref)
    assert np.isclose(gamma, gamma_0 * P_fill)
    assert np.isclose(shift, shift_0 * P_fill)

    # Temperature scaling: n_gamma = 1 - eta = 0.5
    T_hot = T_ref * 1.5
    gamma_hot, _ = get_buffer_broadening_and_shift(elem, dline, bg, P_fill, T_fill, T_hot)
    expected_gamma = gamma_0 * P_fill * (1.5 ** 0.5)
    assert np.isclose(gamma_hot, expected_gamma)

def test_buffer_gas_alqarni_fitted_values():
    """The dynamic K D1/Ne model must reproduce the broadening and shift that
    Alqarni et al. (J. Phys. B 58, 135201 (2025), Fig. 3 / Table 1) obtained by
    fitting each measured spectrum independently, to within experiment.

    Conditions: 75 mm natural-abundance K cell at cell-body temperature
    Tc = 370 K; neon fills with measured pressures 57.4 Torr and 96 Torr at a
    293 K room-temperature fill. Fitted values (Table 1):
        57.4 Torr:  Gamma = 417 +/- 3 MHz,  shift = -87 +/- 2 MHz
        96   Torr:  Gamma = 700 +/- 3 MHz,  shift = -144 +/- 2 MHz
    """
    T_fill_K = 293.0
    Tc_K     = 370.0
    cases = [
        # P_fill (Torr), Gamma_fit, shift_fit
        (57.4, 417.0, -87.0),
        (96.0, 700.0, -144.0),
    ]
    for P_fill, gamma_fit, shift_fit in cases:
        gamma, shift = get_buffer_broadening_and_shift('K', 'D1', 'Ne',
                                                       P_fill, T_fill_K, Tc_K)
        # Agreement with the independently fitted experimental values, within
        # the model-vs-fit gap plus the quoted +/-3 MHz experimental spread.
        assert abs(gamma - gamma_fit) < 5.0, (P_fill, gamma, gamma_fit)
        assert abs(shift - shift_fit) < 3.0, (P_fill, shift, shift_fit)


def test_buffer_gas_integration():
    """Integration test: BufferGas path in calculate() must match explicit GammaBuf/shift."""
    detuning = np.linspace(-10000, 10000, 1000)
    E_in = np.array([1, 0, 0])

    # Use K D1 Ne — the only gas with confirmed (non-placeholder) constants.
    # Fill at room temperature (293.15 K), operate at 60 C (333.15 K).
    T_op_C = 60.0
    T_fill_K = 293.15
    P_fill = 60.0  # Torr

    p_buffer = {
        'Elem': 'K', 'Dline': 'D1', 'T': T_op_C, 'lcell': 75e-3, 'Bfield': 0.0,
        'BufferGas': 'Ne', 'BufferPressure': P_fill, 'BufferTemp_fill': T_fill_K
    }

    # Compute what Gamma and shift the model predicts at T_op
    gamma_p, shift_p = get_buffer_broadening_and_shift('K', 'D1', 'Ne', P_fill, T_fill_K, T_op_C + 273.15)

    p_equivalent = {
        'Elem': 'K', 'Dline': 'D1', 'T': T_op_C, 'lcell': 75e-3, 'Bfield': 0.0,
        'GammaBuf': gamma_p, 'shift': shift_p
    }

    spec_buffer = calculate(detuning, E_in, p_buffer, outputs=['S0'])[0]
    spec_equiv  = calculate(detuning, E_in, p_equivalent, outputs=['S0'])[0]

    np.testing.assert_allclose(spec_buffer, spec_equiv, rtol=1e-5, atol=1e-8)

if __name__ == '__main__':
    pytest.main([__file__])
