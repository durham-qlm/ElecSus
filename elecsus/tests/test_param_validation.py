"""[ES4] Parameter-dictionary validation (libs/param_validation.py).

Before this, `p_dict` keys the code did not recognise were ignored in silence.
Five plausible typos -- including `n_windows` for `n_window` and
`use_nonorthogonal` for `use_non_orthogonal` -- each produced a spectrum
bit-identical to one computed with the ES4 window-Fresnel and non-orthogonal
corrections switched off, with no warning of any kind.  That is the one failure
mode of the package that yields confidently wrong numbers rather than an error,
so it is pinned here.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate
from elecsus.libs import spectra
from elecsus.libs.param_validation import (KNOWN_KEYS, DEPRECATED_KEYS,
                                           check_keys, check_ranges)

_X = np.linspace(-4e3, 4e3, 20)
_BASE = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=100.)


def _calc(**extra):
    return calculate(_X, [1, 0, 0], dict(_BASE, **extra), outputs=['S0'])


# --- unknown keys -------------------------------------------------------

@pytest.mark.parametrize('typo, intended', [
    ('n_windows',         'n_window'),
    ('use_nonorthogonal', 'use_non_orthogonal'),
    ('Btheta_deg',        'Btheta'),
    ('lcel',              'lcell'),
    ('Bfeild',            'Bfield'),
])
def test_mistyped_key_raises(typo, intended):
    """Each of these silently disabled the parameter it was meant to set."""
    with pytest.raises(ValueError, match=typo):
        _calc(**{typo: 1.0})


def test_error_names_the_intended_key():
    """The message has to be actionable: name the key that was meant."""
    with pytest.raises(ValueError, match="did you mean 'n_window'"):
        _calc(n_windows=1.51)


def test_unknown_key_without_a_near_match_still_raises():
    with pytest.raises(ValueError, match='zzz_not_a_parameter'):
        _calc(zzz_not_a_parameter=1.0)


def test_strict_false_downgrades_to_a_warning():
    """The escape hatch, for workflows that carry extra keys in the dict."""
    with pytest.warns(UserWarning, match='n_windows'):
        [s] = calculate(_X, [1, 0, 0], dict(_BASE, n_windows=1.51),
                        outputs=['S0'], strict=False)
    [reference] = _calc()
    np.testing.assert_allclose(s, reference, rtol=1e-12, atol=0)


def test_deprecated_keys_are_still_accepted():
    """'theta0'/'Pol' warn (test_deprecated_params.py) but must not raise."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _calc(theta0=45.)
        _calc(Pol=50.)


def test_known_keys_covers_the_defaults():
    """A key added to p_dict_defaults must be added to KNOWN_KEYS too."""
    assert set(spectra.p_dict_defaults) <= KNOWN_KEYS


def test_known_and_deprecated_do_not_overlap():
    assert not (KNOWN_KEYS & DEPRECATED_KEYS)


def test_every_known_key_is_accepted():
    """Guards against a typo in KNOWN_KEYS itself. The deprecated keys are
    accepted with a warning, which comes from check_keys
    so that fit_data gives it too."""
    with pytest.warns(UserWarning, match='theta0'):
        check_keys({key: None for key in KNOWN_KEYS | DEPRECATED_KEYS})


# --- physical ranges ----------------------------------------------------

def test_negative_cell_length_raises():
    """Returned S0 > 1 -- gain from a passive absorber -- before this check."""
    with pytest.raises(ValueError, match='lcell'):
        _calc(lcell=-2e-3)


def test_zero_cell_length_raises():
    with pytest.raises(ValueError, match='lcell'):
        _calc(lcell=0.0)


@pytest.mark.parametrize('value', [-1.0, 150.0])
def test_isotope_fraction_outside_0_100_raises(value):
    with pytest.raises(ValueError, match='rb85frac'):
        _calc(rb85frac=value)


def test_temperature_below_absolute_zero_raises():
    """Previously a bare LinAlgError naming neither parameter nor problem."""
    with pytest.raises(ValueError, match='absolute zero'):
        _calc(T=-300.)


def test_nan_temperature_raises():
    with pytest.raises(ValueError, match='finite'):
        _calc(T=float('nan'))


def test_window_index_below_one_raises():
    with pytest.raises(ValueError, match='n_window'):
        _calc(n_window=0.0)


def test_window_index_none_is_allowed():
    """None is the documented way to ask for a windowless cell."""
    _calc(use_non_orthogonal=True, n_window=None)


def test_unsupported_element_names_what_is_supported():
    with pytest.raises(ValueError, match='RbD2'):
        _calc(Elem='Xe')


def test_negative_buffer_pressure_raises():
    with pytest.raises(ValueError, match='BufferPressure'):
        check_ranges(dict(_BASE, BufferGas='N2', BufferPressure=-1.0))


def test_buffer_fill_temperature_must_be_kelvin():
    with pytest.raises(ValueError, match='kelvin'):
        check_ranges(dict(_BASE, BufferTemp_fill=-20.))


# --- values that used to pass, then crash or be ignored -------

def test_ideal_atom_is_rejected_by_name():
    """AtomConstants describes an ideal atom, but calc_chi has no model for
    it: Elem='Ideal' passed validation, then raised UnboundLocalError."""
    with pytest.raises(ValueError, match="Elem='Ideal'") as err:
        _calc(Elem='Ideal', Dline='D1')
    assert 'IdealD1' not in str(err.value)


@pytest.mark.parametrize('fractions', [
    dict(K40frac=80., K41frac=50.),
    dict(K41frac=99.995),               # with the default K40frac of 0.01
])
def test_potassium_fractions_may_not_sum_past_100(fractions):
    """Each within 0-100, but the 39K fraction would be negative."""
    with pytest.raises(ValueError, match='39K'):
        _calc(Elem='K', Dline='D2', **fractions)


def test_potassium_fractions_summing_to_100_are_allowed():
    _calc(Elem='K', Dline='D2', K40frac=50., K41frac=50.)


@pytest.mark.parametrize('key, value', [
    ('use_non_orthogonal', 'False'),   # a non-empty string: switched the solver ON
    ('Constrain', 'no'),
    ('Constrain', None),
    ('BoltzmannFactor', 3),
])
def test_switches_must_be_true_or_false(key, value):
    with pytest.raises(ValueError, match=key):
        _calc(**{key: value})


@pytest.mark.parametrize('value', [True, False, 1, 0, np.True_])
def test_switches_accept_booleans_and_one_or_zero(value):
    _calc(use_non_orthogonal=value)


def test_fit_only_polarisation_keys_warn_in_calculate():
    """calculate() takes the field from E_in: these keys gave x-polarised
    light whatever they said."""
    with pytest.warns(UserWarning, match='E_in'):
        _calc(E_x=0., E_y=1., E_phase=0.)


def test_fits_do_not_trigger_the_polarisation_key_warning():
    from elecsus.elecsus_methods import fit_data
    [y] = _calc()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit_data([_X, np.real(y)], dict(_BASE, T=78.), {'T': True},
                 E_in=[1, [0, 0]])
    assert not [w for w in caught if 'fit_data only' in str(w.message)]


def test_window_index_without_the_solver_warns():
    """Ignored in silence before: only the non-orthogonal solver models
    windows."""
    with pytest.warns(UserWarning, match='use_non_orthogonal=True'):
        _calc(use_non_orthogonal=False, n_window=1.51)


def test_window_index_with_the_solver_is_silent():
    """Also with the solver left at its default, which is on."""
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        _calc(use_non_orthogonal=True, n_window=1.51)
        _calc(n_window=1.51)


# --- where warnings point -------------------------------------

@pytest.mark.parametrize('extra', [
    dict(n_windows=1.51),                            # unknown key, strict=False
    dict(BufferGas='N2', BufferPressure=0.),         # from inside calc_chi
    dict(E_x=1.),                                    # fit-only key
])
def test_warnings_point_at_the_callers_line(extra):
    """They named library frames (spectra.get_spectra) through fixed
    stacklevels; now they name the line that called calculate()."""
    with pytest.warns(UserWarning) as record:
        calculate(_X, [1, 0, 0], dict(_BASE, **extra), outputs=['S0'],
                  strict=False)
    assert {w.filename for w in record} == {__file__}


# --- the checks must not disturb a valid calculation --------------------

def test_valid_dict_is_unchanged_by_validation():
    p_dict = dict(_BASE, Btheta=np.pi / 4, use_non_orthogonal=True)
    before = dict(p_dict)
    calculate(_X, [1, 0, 0], p_dict, outputs=['S0'])
    assert p_dict == before
