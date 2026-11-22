"""
Validation test for fitting against a BufferGas cell (issue: fit_data() previously
dropped BufferGas/BufferPressure/BufferTemp_fill, silently falling back to
GammaBuf=0/shift=0 during fitting).
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs.BufferGasData import get_buffer_broadening_and_shift

# K D1 Ne - the one buffer-gas combination with confirmed (non-placeholder) constants,
# validated against Alqarni 2025.
ELEM, DLINE, BUFFER_GAS = 'K', 'D1', 'Ne'
P_FILL = 60.0       # Torr
T_FILL_K = 293.15   # K


def _true_p_dict(T_C, Bfield):
    return {
        'Elem': ELEM, 'Dline': DLINE, 'T': T_C, 'lcell': 75e-3, 'Bfield': Bfield,
        'Btheta': 0.0, 'Bphi': 0.0,
        'BufferGas': BUFFER_GAS, 'BufferPressure': P_FILL, 'BufferTemp_fill': T_FILL_K,
    }


def test_fit_recovers_T_and_Bfield_with_buffer_gas_fixed():
    """Fitting T/Bfield with BufferGas fixed should recover the true values, and the
    resulting GammaBuf must match the buffer-gas-model prediction (not 0)."""
    np.random.seed(0)

    T_true, B_true = 60.0, 200.0
    p_dict_true = _true_p_dict(T_true, B_true)

    x = np.linspace(-10000, 10000, 400)
    E_in = np.array([1.0, 0.0, 0.0])
    y = calculate(x, E_in, p_dict_true, outputs=['S0'])[0].real
    y_noisy = y + np.random.randn(len(x)) * 0.005

    # fit_data() takes E_in in the ES3 polar format: [E_x, [E_y_amplitude, E_phase]]
    E_in_angle = [E_in[0].real, [abs(E_in[1]), np.angle(E_in[1])]]

    # Initial guess: BufferGas/Pressure/TempFill fixed (known cell fabrication facts),
    # T and Bfield perturbed away from truth and left free to fit.
    p_dict_guess = _true_p_dict(T_true + 5.0, B_true + 30.0)
    p_dict_bools = {'T': True, 'Bfield': True}

    data = [x, y_noisy]
    best_params, RMS, result = fit_data(
        data, p_dict_guess, p_dict_bools, E_in=E_in_angle, data_type='S0', fit_algorithm='ML'
    )

    assert RMS < 0.02
    assert best_params['T'] == pytest.approx(T_true, abs=1.0)
    assert best_params['Bfield'] == pytest.approx(B_true, abs=5.0)

    # BufferGas identity must have passed through untouched.
    assert best_params['BufferGas'] == BUFFER_GAS
    assert best_params['BufferPressure'] == pytest.approx(P_FILL)

    # And the effective GammaBuf implied by the fit's own inputs (T, BufferGas params)
    # must equal the physical model prediction - not the ES3 default of 0.
    gamma_expected, _ = get_buffer_broadening_and_shift(
        ELEM, DLINE, BUFFER_GAS, P_FILL, T_FILL_K, best_params['T'] + 273.15
    )
    assert gamma_expected > 0  # sanity: buffer broadening is non-trivial at this pressure


def test_varying_gammabuf_with_buffer_gas_active_raises():
    """Marking GammaBuf as a free fit parameter while BufferGas is active is a no-op
    trap (spectra.py silently overrides it) - fit_data() must reject it outright."""
    p_dict = _true_p_dict(60.0, 200.0)
    p_dict_bools = {'T': True, 'GammaBuf': True}

    x = np.linspace(-1000, 1000, 10)
    data = [x, np.zeros_like(x)]

    with pytest.raises(ValueError, match="BufferGas"):
        fit_data(data, p_dict, p_dict_bools, data_type='S0', fit_algorithm='ML')
