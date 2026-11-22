# Copyright 2014-2026 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams, I. G. Hughes and S. A. Wrathmall.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#	 http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
[ES4 NEW] Validation of the parameter dictionary.

Up to and including ES3, `p_dict` was consumed by direct indexing and `.get()`
with defaults, and unknown keys were never checked.  A mistyped key was
therefore silently ignored: `n_windows` instead of `n_window`, or
`use_nonorthogonal` instead of `use_non_orthogonal`, gave a spectrum
bit-identical to one computed with the ES4 window-Fresnel and non-orthogonal
corrections switched off, with no warning of any kind.  In a fitting context
that produces converged, publishable, wrong parameters.

This module supplies the two checks that close that route:

* `check_keys`   -- every key must be one the code actually reads.
* `check_ranges` -- every value must lie in the range where the model is
                    defined (a negative `lcell` otherwise returns a spectrum
                    with S0 > 1, i.e. gain from an absorbing vapour).

Both are applied by `spectra.get_spectra`, which is the single entry point for
`elecsus_methods.calculate`, `calculate_cascade` and all three fitting
backends.  Pass `strict=False` to `get_spectra`/`calculate` to downgrade the
unknown-key error to a warning if a workflow genuinely needs to carry extra
keys through the dictionary; the range checks always apply.
"""

import difflib
import numbers
import os
import sys
import warnings

import numpy as np

from . import AtomConstants as AC

_PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_TESTS_DIR = os.path.join(_PACKAGE_DIR, 'tests')


def _in_library(filename):
    """True for a frame in ElecSus's own code; its tests count as user code."""
    path = os.path.abspath(filename)
    return (path.startswith(_PACKAGE_DIR + os.sep)
            and not path.startswith(_TESTS_DIR + os.sep))


def warn(message, category=UserWarning):
    """`warnings.warn`, attributed to the first caller outside ElecSus.

    A fixed stacklevel points wherever the call chain happens to be that deep:
    check_keys' stacklevel=3 named spectra.get_spectra rather than the user's
    line. Walking out of the package finds the user's line from
    any entry point (calculate, calculate_cascade, fit_data, ...).
    """
    frame, level = sys._getframe(1), 2
    while frame is not None and _in_library(frame.f_code.co_filename):
        frame, level = frame.f_back, level + 1
    warnings.warn(message, category, stacklevel=level)

__all__ = ['KNOWN_KEYS', 'DEPRECATED_KEYS', 'SUPPORTED_ELEMENTS',
           'FIT_ONLY_KEYS', 'FIT_PARAMETERS', 'FIT_LIMITS',
           'ParameterRangeError', 'validate_p_dict', 'check_keys',
           'check_fit_keys', 'check_ranges', 'check_consistency',
           'clip_bounds_to_limits', 'search_widths', 'warn', 'BUFFER_TEMP_FILL_DEFAULT']


class ParameterRangeError(ValueError):
    """A parameter value outside the range where the model is defined.

    Raised by `check_ranges`. It is a ValueError, so existing handlers still
    catch it; the fitting code tells it apart from other ValueErrors so that a
    trial value an optimiser steps to is reported as such.
    """

#: Every key the package reads out of a parameter dictionary.  Keep in step
#: with `spectra.calc_chi`, `spectra.get_Efield`, `spectra.get_spectra` and the
#: `p_dict.setdefault` block in `elecsus_methods.fit_data`; the test suite
#: (`test_param_validation.py`) asserts that this set covers
#: `spectra.p_dict_defaults`.
KNOWN_KEYS = frozenset({
    # atom / line
    'Elem', 'Dline',
    # cell and vapour
    'T', 'lcell', 'DoppTemp', 'Constrain', 'BoltzmannFactor',
    'rb85frac', 'K40frac', 'K41frac',
    # magnetic field
    'Bfield', 'Btheta', 'Bphi',
    # line broadening / shift
    'GammaBuf', 'shift',
    # dynamic buffer-gas model (Pitz 2014)
    'BufferGas', 'BufferPressure', 'BufferTemp_fill',
    # [ES4] non-orthogonal propagation and cell windows (Briscoe et al. 2025)
    'use_non_orthogonal', 'n_window', 'wavevector_dirn',
    # input polarisation, as used by the fitting routines
    'E_x', 'E_y', 'E_phase',
})

#: ES3-only keys that are accepted (with a warning from `get_spectra`) but have
#: no effect in ES4, where the input polarisation is the explicit E_in vector.
DEPRECATED_KEYS = frozenset({'theta0', 'Pol'})

#: The elements `spectra.calc_chi` models. `AtomConstants.transitions` also
#: describes an ideal atom, which calc_chi does not.
SUPPORTED_ELEMENTS = ('Rb', 'Cs', 'K', 'Na')

#: Keys only `fit_data` reads: it builds the input field from them. `calculate`
#: takes the field from its E_in argument instead.
FIT_ONLY_KEYS = ('E_x', 'E_y', 'E_phase')

#: The parameters a fit can vary: the numeric arguments of
#: `MLFittingRoutine.fit_function`, which lmfit turns into its Parameters.
#: `p_dict_bools` and `p_dict_bounds` may name only these;
#: `test_fit_keys.py` checks the list against lmfit's.
FIT_PARAMETERS = ('E_x', 'E_y', 'E_phase', 'T', 'lcell', 'Bfield', 'Btheta',
                  'Bphi', 'GammaBuf', 'shift', 'DoppTemp', 'rb85frac',
                  'K40frac', 'K41frac', 'BufferPressure', 'BufferTemp_fill')

_ABOVE_ABSOLUTE_ZERO = float(np.nextafter(-273.15, np.inf))   # degrees C
_ABOVE_ZERO = float(np.nextafter(0.0, np.inf))

#: The values a fit may try for each numeric parameter that has a physical
#: limit, as (lowest, highest). These are the limits `check_ranges` enforces;
#: where it requires a strict inequality (lcell > 0), the limit here is the
#: nearest float inside. `fit_data` and the GUI keep fit bounds within them
#:; `test_fit_limits.py` checks they agree with `check_ranges`.
#: One limit is joint and so not in the table: K40frac + K41frac <= 100. A fit
#: of both that crosses it stops with ParameterRangeError, naming it.
FIT_LIMITS = {
    'T': (_ABOVE_ABSOLUTE_ZERO, np.inf),
    'DoppTemp': (_ABOVE_ABSOLUTE_ZERO, np.inf),
    'lcell': (_ABOVE_ZERO, np.inf),
    'GammaBuf': (0.0, np.inf),
    'rb85frac': (0.0, 100.0),
    'K40frac': (0.0, 100.0),
    'K41frac': (0.0, 100.0),
    'BufferPressure': (0.0, np.inf),
    'BufferTemp_fill': (_ABOVE_ZERO, np.inf),
}


def _is_finite_number(value):
    """True for anything that behaves as a finite real scalar."""
    if isinstance(value, bool):
        return False
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def _described(keys, candidates):
    """The keys, quoted, each with the closest of `candidates` as a hint."""
    described = []
    for key in keys:
        close = difflib.get_close_matches(str(key), candidates, n=1, cutoff=0.6)
        described.append(f"'{key}'" + (f" (did you mean '{close[0]}'?)" if close else ''))
    return ', '.join(described)


def check_keys(p_dict, strict=True):
    """Raise on keys the code does not read; warn on deprecated ones.

    A mistyped key is otherwise ignored in silence, which is the one failure
    mode of this package that produces confidently wrong numbers rather than an
    error.  With `strict=False` an unknown key downgrades to a UserWarning.
    """
    # [ES4] 'theta0'/'Pol' are ES3-only polarisation keys with no effect here:
    # the input polarisation is the explicit E_in vector.  Warn
    # rather than silently drop a ported p_dict's polarisation angle. (This
    # was in get_spectra, which fits bypass with a rebuilt dict.)
    if set(p_dict) & DEPRECATED_KEYS:
        warn("ElecSus 4 ignores the 'theta0'/'Pol' p_dict keys: set the input "
             "polarisation through the explicit E_in vector (and fit it via the "
             "'E_x'/'E_y'/'E_phase' parameters). These keys have no effect and "
             "will be removed in a future release.")

    unknown = sorted(set(p_dict) - KNOWN_KEYS - DEPRECATED_KEYS)
    if not unknown:
        return

    message = (
        'unknown parameter ' + ('keys' if len(unknown) > 1 else 'key') + ': '
        + _described(unknown, sorted(KNOWN_KEYS)) + '. '
        'ElecSus ignores keys it does not recognise, so a mistyped key silently '
        'leaves that parameter at its default -- see spectra.p_dict_defaults for '
        'the accepted names. Pass strict=False to downgrade this to a warning.'
    )
    if strict:
        raise ValueError(message)
    warn(message)


def check_fit_keys(entries, name, strict=True):
    """Return `entries` (`p_dict_bools` or `p_dict_bounds`) if every key is a
    fit parameter.

    A key that is not one used to raise a bare KeyError from lmfit
    (`p_dict_bools`) or to be ignored in silence (`p_dict_bounds`), so a
    mistyped name left its parameter fixed, or unbounded. It now raises
    ValueError naming the key; with `strict=False` it warns instead, and the
    copy returned leaves the key out.
    """
    unknown = sorted(set(entries) - set(FIT_PARAMETERS), key=str)
    if not unknown:
        return entries
    message = (
        f"{name}: " + _described(unknown, FIT_PARAMETERS)
        + (' are not fit parameters' if len(unknown) > 1 else ' is not a fit parameter')
        + '. The fit parameters are ' + ', '.join(FIT_PARAMETERS) + '. '
        'Pass strict=False to downgrade this to a warning.'
    )
    if strict:
        raise ValueError(message)
    warn(message)
    return {key: value for key, value in entries.items() if key not in unknown}


def check_ranges(p_dict):
    """Raise on values outside the range where the model is defined.

    These are hard physical limits, not fitting bounds: a negative cell length
    or an out-of-range isotope fraction returns a spectrum that looks plausible
    but violates energy conservation (S0 > 1 from a passive absorber), and a
    temperature below absolute zero reaches the solver as a bare
    `LinAlgError: Singular matrix` naming neither the parameter nor the problem.
    """
    def get(key):
        return p_dict[key] if key in p_dict else None

    # -- element and line ------------------------------------------------
    # AtomConstants also describes an ideal atom ('IdealD1'), but calc_chi has
    # no susceptibility model for it: Elem='Ideal' passed here and then raised
    # UnboundLocalError inside calc_chi.
    elem, dline = get('Elem'), get('Dline')
    if elem is not None or dline is not None:
        elem = elem if elem is not None else 'Rb'
        dline = dline if dline is not None else 'D2'
        if elem not in SUPPORTED_ELEMENTS or elem + str(dline) not in AC.transitions:
            supported = sorted(t for t in AC.transitions
                               if t[:-2] in SUPPORTED_ELEMENTS)
            raise ParameterRangeError(
                f"no atomic data for Elem='{elem}', Dline='{dline}'. "
                f"Supported combinations are: {', '.join(supported)}.")

    # -- finite real scalars ---------------------------------------------
    # 'n_window' is an exception: None is the documented way to ask for a
    # windowless cell (see spectra.propagate_non_orthogonal). So is
    # 'DoppTemp': None, its default, means it follows T.
    for key in ('T', 'DoppTemp', 'lcell', 'Bfield', 'Btheta', 'Bphi',
                'GammaBuf', 'shift', 'rb85frac', 'K40frac', 'K41frac',
                'n_window', 'BufferPressure', 'BufferTemp_fill',
                'E_x', 'E_y', 'E_phase'):
        if key not in p_dict:
            continue
        if key in ('n_window', 'DoppTemp') and p_dict[key] is None:
            continue
        if not _is_finite_number(p_dict[key]):
            raise ParameterRangeError(
                f"'{key}' must be a finite real number, got {p_dict[key]!r}.")

    # -- temperatures (degrees Celsius) ----------------------------------
    for key in ('T', 'DoppTemp'):
        value = get(key)
        if value is not None and float(value) <= -273.15:
            raise ParameterRangeError(
                f"'{key}' = {value} is at or below absolute zero. ElecSus takes "
                "cell temperatures in degrees Celsius (not kelvin).")

    # -- cell length ------------------------------------------------------
    lcell = get('lcell')
    if lcell is not None and float(lcell) <= 0:
        raise ParameterRangeError(
            f"'lcell' = {lcell} must be positive. ElecSus takes the cell length "
            "in metres (75e-3 for a 75 mm cell); a negative length returns a "
            "spectrum with transmission above unity.")

    # -- isotope fractions (per cent) -------------------------------------
    for key in ('rb85frac', 'K40frac', 'K41frac'):
        value = get(key)
        if value is not None and not 0.0 <= float(value) <= 100.0:
            raise ParameterRangeError(
                f"'{key}' = {value} is outside 0-100. Isotope abundances are "
                "given as a percentage (rb85frac=72.17 for natural rubidium).")

    # The 39K fraction is what is left, 100 - K40frac - K41frac: each could
    # be within 0-100 while their sum passed 100. A fraction not
    # given takes its default, as in calc_chi.
    if get('K40frac') is not None or get('K41frac') is not None:
        from .spectra import p_dict_defaults   # here: spectra imports this module
        k40 = float(p_dict.get('K40frac', p_dict_defaults['K40frac']))
        k41 = float(p_dict.get('K41frac', p_dict_defaults['K41frac']))
        if k40 + k41 > 100.0:
            raise ParameterRangeError(
                f"'K40frac' + 'K41frac' = {k40:g} + {k41:g} exceeds 100: the "
                "39K fraction, 100 - K40frac - K41frac, would be negative.")

    # -- switches -----------------------------------------------------------
    # A string is truthy, so use_non_orthogonal='False' switched the solver
    # on. True/False, or 1/0, only; None is not "off".
    for key in ('use_non_orthogonal', 'Constrain', 'BoltzmannFactor'):
        if key not in p_dict:
            continue
        value = p_dict[key]
        if not (isinstance(value, (bool, np.bool_))
                or (isinstance(value, numbers.Integral) and value in (0, 1))):
            raise ParameterRangeError(
                f"'{key}' = {value!r} must be True or False.")

    # -- non-negative quantities ------------------------------------------
    for key, unit in (('GammaBuf', 'MHz'), ('BufferPressure', 'Torr')):
        value = get(key)
        if value is not None and float(value) < 0:
            raise ParameterRangeError(f"'{key}' = {value} must not be negative ({unit}).")

    fill_temp = get('BufferTemp_fill')
    if fill_temp is not None and float(fill_temp) <= 0:
        raise ParameterRangeError(
            f"'BufferTemp_fill' = {fill_temp} must be positive: unlike 'T', the "
            "buffer-gas fill temperature is in kelvin (293.15 by default).")

    # -- propagation direction --------------------------------------------
    # Anything but +1 or -1 used to scale the wavenumber: 2 doubled the
    # optical depth in silence.
    # A string such as '-1' is rejected too: the propagation code compares
    # the value with -1, so it would quietly run a forward pass.
    dirn = get('wavevector_dirn')
    if dirn is not None and (isinstance(dirn, (bool, np.bool_))
                             or not isinstance(dirn, numbers.Real)
                             or dirn not in (1, -1)):
        raise ParameterRangeError(
            f"'wavevector_dirn' = {dirn!r} must be 1 (light travelling along "
            "+z, the default) or -1 (the return pass of a double-pass "
            "geometry, with B given in the same lab frame).")

    # -- window refractive index ------------------------------------------
    n_window = get('n_window')  # None means "no windows", and is allowed
    if n_window is not None and float(n_window) < 1.0:
        raise ParameterRangeError(
            f"'n_window' = {n_window} must be at least 1: it is the refractive "
            "index of the cell windows (1.51 for borosilicate). Set "
            "n_window=None to model a cell without windows.")


def _bound_numbers(entry):
    """The numbers in one `p_dict_bounds` entry, as floats: a single number,
    or a list, tuple or array of them. None if the entry holds anything else
    (a string, a bool, NaN, a nested list)."""
    if isinstance(entry, np.ndarray):
        entry = entry.tolist()
    values = []
    for item in (entry if isinstance(entry, (list, tuple)) else [entry]):
        if isinstance(item, (bool, np.bool_)) or not isinstance(item, numbers.Real):
            return None
        if np.isnan(float(item)):
            return None
        values.append(float(item))
    return values


def clip_bounds_to_limits(p_dict_bounds):
    """Return a copy of fit bounds as [min, max] pairs narrowed to `FIT_LIMITS`.

    A fit stops at the first trial value outside the range where the model is
    defined: bounds of [-100, 100] MHz on GammaBuf made differential
    evolution try negative widths and fail. Bounds that reach
    past a physical limit are narrowed to it, with a warning; bounds lying
    wholly outside it raise.

    This is the form ML and differential-evolution fits take: each entry must
    be a [min, max] pair with min < max, or ValueError names the parameter
. A single search width, the form for random restart
    (`search_widths`), used to fail here with "not enough values to unpack",
    and a reversed pair was reported as out of range.
    """
    clipped = {}
    for key, entry in p_dict_bounds.items():
        pair = _bound_numbers(entry)
        if pair is None or len(pair) != 2:
            raise ValueError(
                f"p_dict_bounds['{key}'] = {entry!r} is not a [min, max] pair. "
                "ML and differential-evolution fits take the range each "
                "parameter may take, e.g. 'T': [70, 130]; a single width such "
                "as 'T': 30 is the form for fit_algorithm 'RR'.")
        lo, hi = pair
        if not lo < hi:
            raise ValueError(
                f"p_dict_bounds['{key}'] = {entry!r}: the minimum must be "
                "below the maximum.")
        if key in FIT_LIMITS:
            low_limit, high_limit = FIT_LIMITS[key]
            new_lo, new_hi = max(lo, low_limit), min(hi, high_limit)
            if new_lo > new_hi:
                raise ParameterRangeError(
                    f"p_dict_bounds['{key}'] = [{lo}, {hi}] lies outside the "
                    f"physical range of '{key}', [{low_limit}, {high_limit}].")
            if (new_lo, new_hi) != (lo, hi):
                warn(f"p_dict_bounds['{key}'] = [{lo}, {hi}] reaches outside the "
                     f"physical range of '{key}'; fitting within "
                     f"[{new_lo:.10g}, {new_hi:.10g}] instead.")
            lo, hi = new_lo, new_hi
        clipped[key] = [lo, hi]
    return clipped


def search_widths(p_dict_bounds):
    """Return random-search widths as floats.

    For random restart, `p_dict_bounds[key]` is the width of the random
    search around the starting value of `key`, not a hard limit: one number
    of at least 0, or a list holding one ([30], as ES3's fitting demo
    writes it). A [min, max] pair, the form ML and
    differential evolution take, raises ValueError naming the parameter
    rather than being misread.
    """
    widths = {}
    for key, entry in p_dict_bounds.items():
        values = _bound_numbers(entry)
        if values is None or len(values) != 1 or not 0 <= values[0] < np.inf:
            raise ValueError(
                f"p_dict_bounds['{key}'] = {entry!r} is not a search width. "
                "Random-restart fits take the width of the random search "
                "around each starting value: one number of at least 0, e.g. "
                "'T': 30. A [min, max] pair is the form for ML and "
                "differential evolution.")
        widths[key] = values[0]
    return widths


#: The fill temperature the buffer-gas model takes when BufferTemp_fill is
#: not given (K): spectra.calc_chi's, and the fit model's.
BUFFER_TEMP_FILL_DEFAULT = 293.15

#: Warned when Constrain=False comes without DoppTemp, which then has nothing
#: to decouple; `fit_data` warns the same.
NO_DOPPTEMP = ("Constrain=False has no effect without DoppTemp: the Doppler width "
               "and the buffer-gas model follow T. Set DoppTemp to the temperature "
               "of the vapour in the beam (the cell body).")


def check_consistency(p_dict):
    """Warn about settings that are valid but have no effect.

    These cases used to pass in silence, leaving the user believing the
    setting had been applied.
    """
    given = [key for key in FIT_ONLY_KEYS if key in p_dict]
    if given:
        warn(f"{', '.join(repr(k) for k in given)} set the input polarisation "
             "for fit_data only. calculate() takes the input field from its "
             "E_in argument, so these keys have no effect here.")
    if (p_dict.get('n_window') is not None
            and not p_dict.get('use_non_orthogonal', True)):
        warn(f"n_window = {p_dict['n_window']} has no effect: the cell windows "
             "are modelled only by the non-orthogonal solver. Set "
             "use_non_orthogonal=True to include them.")
    if not p_dict.get('Constrain', True) and p_dict.get('DoppTemp') is None:
        warn(NO_DOPPTEMP)
    # The buffer-gas model runs only with BufferGas set; fit_data
    # already refused to vary BufferPressure without it.
    if p_dict.get('BufferGas') is None:
        given = []
        if p_dict.get('BufferPressure', 0.0) > 0:
            given.append(f"BufferPressure = {p_dict['BufferPressure']} Torr")
        if p_dict.get('BufferTemp_fill', BUFFER_TEMP_FILL_DEFAULT) != BUFFER_TEMP_FILL_DEFAULT:
            given.append(f"BufferTemp_fill = {p_dict['BufferTemp_fill']} K")
        if given:
            warn(f"{' and '.join(given)} {'has' if len(given) == 1 else 'have'} no "
                 "effect without BufferGas: the spectrum has no buffer-gas broadening "
                 "or shift. Set 'BufferGas' to the gas the cell was filled with.")


def validate_p_dict(p_dict, strict=True):
    """Run all three checks. Called by `spectra.get_spectra`."""
    check_keys(p_dict, strict=strict)
    check_ranges(p_dict)
    check_consistency(p_dict)
