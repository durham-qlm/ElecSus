"""
The buffer-gas model in calculations and fits.

- The model is evaluated at the temperature of the gas in the beam: T, or
  DoppTemp when Constrain=False gives the cell body its own temperature.
- Its broadening replaces GammaBuf; its pressure shift is now added to shift,
  which it used to replace. shift is then a frequency offset on top of the
  model, which a fit can vary.
- A fit can vary BufferPressure, to find the amount of gas in the cell. It
  cannot vary BufferTemp_fill: the model depends on BufferPressure /
  BufferTemp_fill only, so the two are degenerate. Nor can it vary
  BufferPressure while the model is off, where it has no effect, or GammaBuf
  while it is on. fit_data refuses each of these by name.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import _check_buffer_gas_fit, calculate, fit_data
from elecsus.libs import RRFittingRoutine as RR
from elecsus.libs.BufferGasData import get_buffer_broadening_and_shift
from elecsus.libs.param_validation import FIT_ONLY_KEYS, FIT_PARAMETERS

X = np.linspace(-6000, 8000, 200)
E_IN = [1, 0, 0]
CELL = dict(Elem='Rb', Dline='D2', T=60., lcell=75e-3, Bfield=0.)
GAS = dict(BufferGas='N2', BufferPressure=10., BufferTemp_fill=293.15)
TRUE = dict(CELL, **GAS)


def _S0(X, p_dict):
    return np.real(calculate(X, E_IN, p_dict, outputs=['S0'])[0])


def _data(**extra):
    rng = np.random.default_rng(1)
    return [X, _S0(X, dict(TRUE, **extra)) + rng.normal(0, 0.002, len(X))]


def _fit(data, start, bools):
    with warnings.catch_warnings():
        warnings.simplefilter('error')   # the model is on: nothing to warn about
        return fit_data(data, dict(TRUE, **start), bools)


# --- calculate --------------------------------------------------------------------

def _static(cell, T_eval):
    """The model's broadening and shift at T_eval (C), given as GammaBuf and shift."""
    gamma, shift = get_buffer_broadening_and_shift(
        'Rb', 'D2', 'N2', GAS['BufferPressure'], GAS['BufferTemp_fill'], T_eval + 273.15)
    return _S0(X, dict(cell, GammaBuf=gamma, shift=shift))


def test_the_model_is_evaluated_at_DoppTemp_in_a_two_temperature_cell():
    """An 80 C stem and a 100 C body (Constrain=False). The
    collisions happen in the beam, at the body's temperature. The model used to
    take T, the stem's."""
    cell = dict(CELL, T=80., DoppTemp=100., Constrain=False)
    model = _S0(X, dict(cell, **GAS))
    assert np.array_equal(model, _static(cell, 100.))
    assert np.max(np.abs(model - _static(cell, 80.))) > 1e-3


def test_with_constrain_the_model_is_evaluated_at_T():
    """One temperature, as before: a DoppTemp is ignored with Constrain=True."""
    cell = dict(CELL, T=80., DoppTemp=100.)
    assert np.array_equal(_S0(X, dict(cell, **GAS)), _static(cell, 80.))


def test_a_gas_without_data_falls_back_to_the_static_model():
    """(an untested path found by the coverage report) Na has no
    buffer-gas data: calc_chi warns and uses GammaBuf and shift as given."""
    X = np.linspace(-3000, 3000, 61)
    static = dict(Elem='Na', Dline='D2', T=150., lcell=5e-3, GammaBuf=20., shift=5.)
    with pytest.warns(UserWarning, match="No buffer gas data for"):
        [with_gas] = calculate(X, [1, 0, 0], dict(static, BufferGas='N2', BufferPressure=10.),
                               outputs=['S0'])
    [without] = calculate(X, [1, 0, 0], static, outputs=['S0'])
    np.testing.assert_array_equal(with_gas, without)


def test_shift_is_added_to_the_pressure_shift():
    """It was replaced by it, with a warning. Now the spectrum moves by shift,
    as it does without the model, and nothing warns."""
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        shifted = _S0(X, dict(TRUE, shift=50.))
    assert np.array_equal(shifted, _S0(X - 50., TRUE))


def test_gammabuf_is_still_replaced_by_the_model():
    with pytest.warns(UserWarning, match='GammaBuf=30.0 MHz is overridden'):
        with_gammabuf = _S0(X, dict(TRUE, GammaBuf=30.))
    assert np.array_equal(with_gammabuf, _S0(X, TRUE))


# --- fits that work ------------------------------------------------------------------

def test_a_fit_finds_the_fill_pressure():
    best, rms, result = _fit(_data(), {'BufferPressure': 6.}, {'BufferPressure': True})
    assert best['BufferPressure'] == pytest.approx(10., abs=0.1)
    assert 0 < result.params['BufferPressure'].stderr < 0.1
    assert rms < 0.003


def test_a_fit_finds_the_pressure_with_a_frequency_offset():
    """The width fixes the pressure, so a free shift takes only what the
    pressure shift does not account for: here a 30 MHz offset of the axis."""
    best, _, _ = _fit(_data(shift=30.), {'BufferPressure': 6.},
                      {'BufferPressure': True, 'shift': True})
    assert best['BufferPressure'] == pytest.approx(10., abs=0.1)
    assert best['shift'] == pytest.approx(30., abs=1.5)


def test_a_fit_finds_the_pressure_and_the_temperature():
    best, _, _ = _fit(_data(), {'BufferPressure': 6., 'T': 55.},
                      {'BufferPressure': True, 'T': True})
    assert best['BufferPressure'] == pytest.approx(10., abs=0.15)
    assert best['T'] == pytest.approx(60., abs=0.2)


def test_the_fitted_pressure_goes_back_into_calculate():
    """As every fit result does: with the gas, the fill
    temperature and the fitted pressure."""
    best, _, result = _fit(_data(), {'BufferPressure': 6.}, {'BufferPressure': True})
    assert (best['BufferGas'], best['BufferTemp_fill']) == ('N2', 293.15)
    p_dict = {key: value for key, value in best.items() if key not in FIT_ONLY_KEYS}
    assert np.allclose(_S0(X, p_dict), result.best_fit, rtol=0, atol=1e-12)


# --- fits that are refused -------------------------------------------------------

def _without(key):
    return {k: v for k, v in TRUE.items() if k != key}


@pytest.mark.parametrize('p_dict, bools, match', [
    # degenerate with the pressure, with or without the model on
    (TRUE, {'BufferTemp_fill': True}, 'degenerate with BufferPressure'),
    (TRUE, {'BufferPressure': True, 'BufferTemp_fill': True}, 'degenerate'),
    # no effect: the model is off
    (_without('BufferGas'), {'BufferPressure': True}, 'BufferGas=None'),
    (dict(TRUE, BufferPressure=0.), {'BufferPressure': True}, 'above 0 Torr'),
    (dict(TRUE, Elem='Na', T=150.), {'BufferPressure': True}, 'for Na D2'),
    # replaced by the model's broadening
    (TRUE, {'GammaBuf': True}, 'replaces GammaBuf'),
    (TRUE, {'T': True, 'GammaBuf': True}, 'replaces GammaBuf'),
])
def test_a_parameter_the_model_ignores_is_refused(p_dict, bools, match):
    with pytest.raises(ValueError, match=match):
        fit_data([X, np.ones_like(X)], p_dict, bools)


@pytest.mark.parametrize('p_dict, bools', [
    (TRUE, {'shift': True}),
    (TRUE, {'BufferPressure': True, 'shift': True, 'T': True}),
    (_without('BufferGas'), {'GammaBuf': True}),
    # the model is off: GammaBuf is used, as the warning says
    (dict(TRUE, BufferPressure=0.), {'GammaBuf': True}),
    (dict(TRUE, Elem='Na'), {'GammaBuf': True}),
])
def test_a_parameter_the_model_uses_is_allowed(p_dict, bools):
    _check_buffer_gas_fit(p_dict, bools)


def test_random_restart_has_a_search_width_for_the_pressure():
    """RR read each varied parameter's width from its defaults, which had
    none for BufferPressure: KeyError. BufferTemp_fill cannot vary."""
    assert set(FIT_PARAMETERS) - {'BufferTemp_fill'} <= set(RR.p_dict_bounds_default)
