# Copyright 2015-2026 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams, I. G. Hughes and S. A. Wrathmall.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This is the highest-level interface to calculate atomic spectra. It provides two methods,
to either calculate spectra or fit experimental data to the model.

Note: This is incompatible with the old way of running elecsus, using the runcard.py module

Usage, and examples, are detailed in each method separately

Last updated 2018-07-04 MAZ
"""
# py 2.7 compatibility

import difflib

import numpy as np
from numpy import arange, zeros, array, sqrt

# import elecsus modules
from .libs import spectra
from .libs.BufferGasData import BUFFER_GAS_DATA
from .libs.param_validation import (DEPRECATED_KEYS, KNOWN_KEYS, NO_DOPPTEMP, check_fit_keys,
                                    check_keys, check_ranges, clip_bounds_to_limits, warn)
from .libs import MLFittingRoutine as ML
from .libs import RRFittingRoutine as RR


def calculate(detuning_range, E_in=[1, 0, 0], p_dict={}, outputs=None, strict=True):
    """
    Alias for the get_spectra() method in libs.spectra.

    Inputs:
        detuning_range [ numpy 1D array ]
            The independent variable and defines the detuning points over which to calculate. Values in MHz

        E_in [ numpy 1/2D array ]
            Defines the input electric field vector in the xyz basis.
            The z-axis is always the direction of propagation (independent of the magnetic field axis), and therefore
            the electric field should be a plane wave in the x,y plane. The array passed to this method
            should be in one of two formats:
                (1) A 1D array of (Ex,Ey,Ez) which is the input electric field for all detuning values;
                or
                (2) A 2D array with dimensions (3,len(detuning_range)) - i.e. each detuning has a different electric
                field associated with it - which will happen on propagation through a birefringent/dichroic medium

        p_dict [ dictionary ]
            Dictionary containing all parameters (the order of parameters is therefore not important)
                Dictionary keys:

                Key				DataType	Unit		Description
                ---				---------	----		-----------
                Elem	   			str			--			The chosen alkali element.
                Dline	  			str			--			Specifies which D-line transition to calculate for (D1 / D2)

                # Magnetic field parameters
                Bfield	 			float			Gauss	Magnitude of the applied magnetic field
                Btheta			float			radians	Angle of B from +z, the forward pass's direction
                Bphi				float			radians	Azimuth of B, from x towards y: B is along
                                                                (sin Btheta cos Bphi, sin Btheta sin Bphi, cos Btheta)
                                                                [ES4 FIX] ES3 computed the mirror image:
                                                                its result for (Btheta, Bphi) is ES4's for
                                                                (pi - Btheta, -Bphi)

                # Experimental parameters
                T		 			float			Celsius	Temperature used to calculate atomic number density; with
                                                                Constrain=True (the default), also the temperature of the
                                                                vapour in the beam (see DoppTemp)
                GammaBuf   	float			MHz		Extra lorentzian broadening (usually from buffer gas
                                                                but can be any extra homogeneous broadening).
                                                                The buffer-gas model's broadening replaces it.
                BufferGas   		str			--		Optional string (e.g. 'He', 'N2') to dynamically model GammaBuf/shift
                BufferPressure   	float			Torr		Fill pressure of BufferGas to act on
                BufferTemp_fill   	float			K		Temperature cell was filled with BufferGas (default 293.15)
                shift	  			float			MHz		A global frequency shift of the atomic resonance frequencies.
                                                                [ES4 MOD] The buffer-gas model's pressure shift is added
                                                                to it (it used to replace it).
                DoppTemp   	float			Celsius	Temperature of the vapour in the beam (the cell body, when a
                                                                colder stem sets T): the Doppler width, and the buffer-gas
                                                                model's broadening and shift. Used only with
                                                                Constrain=False; not given, it is T
                Constrain  		bool			--			If True, overides the DoppTemp value and sets it to T
                lcell	  			float			m			length of the vapour cell

                # Elemental abundancies, where applicable
                rb85frac   		float			%			percentage of rubidium-85 atoms
                K40frac			float			%			percentage of potassium-40 atoms
                K41frac			float			%			percentage of potassium-41 atoms


                NOTE: If keys are missing from p_dict, default values contained in p_dict_defaults will be loaded.

        outputs [ list of strings ]
            Keyword argument that defines the quantities that are returned.
            If not specified, defaults to None, in which case a default set of outputs is returned, which are:
                S0, S1, S2, S3, Ix, Iy, I_P45, I_M45, alphaPlus, alphaMinus, alphaZ

    Returns:
        A list of output arrays as defined by the 'outputs' keyword argument.


    Example usage:
        To calculate the room temperature absorption of a 75 mm long Cs reference cell in an applied magnetic field of
        100 G aligned along the direction of propagation (Faraday geometry), between -10 and +10 GHz, with an input
        electric field aligned along the x-axis:

        > detuning_range = np.linspace(-10,10,1000)*1e3 # GHz to MHz conversion
        > E_in = np.array([1,0,0])
        > p_dict = {'Elem':'Cs', 'Dline':'D2', 'Bfield':100, 'T':21, 'lcell':75e-3}
        > [Transmission] = calculate(detuning_range,E_in,p_dict,outputs=['S0'])

        More examples available in the /tests/ directory
    """

    return spectra.get_spectra(detuning_range, E_in, p_dict, outputs, strict=strict)


# [ES4 BEGIN] calculate_cascade — two cascaded vapour cells (not in ES3)
def calculate_cascade(detuning_range, E_in=[1, 0, 0], p_dict_1={}, p_dict_2={},
                      outputs=None, polariser_angle=None, reset_polarisation=None,
                      strict=True):
    """
    Calculate spectra through two cascaded vapour cells of the same species.

    The output electric field from cell 1 is used as the input to cell 2,
    preserving both the polarisation state and the intensity attenuation.
    Optionally a linear polariser is inserted between the cells — the
    magneto-optical filter (MOF) topology, in which the first cell sits
    between crossed polarisers and only light whose polarisation the vapour
    has rotated reaches the second (wing-selector) cell. A quarter-wave
    plate immediately after that polariser can also be modelled via
    reset_polarisation, for wing-selector designs where cell 2 must always
    receive a fixed circular polarisation (e.g. VAMOS, Oliviero 1998).

    Inputs:
        detuning_range  numpy 1D array    Detuning points in MHz
        E_in            array             Input electric field [Ex, Ey, Ez]
        p_dict_1        dict              Parameters for cell 1 (see calculate())
        p_dict_2        dict              Parameters for cell 2 (see calculate())
        outputs         list of str       Quantities to return per cell (see calculate()).
                                          If None, returns the default set:
                                          S0, S1, S2, S3, Ix, Iy, I_P45, I_M45,
                                          alphaPlus, alphaMinus, alphaZ
        polariser_angle float or None     If not None, a linear polariser at this
                                          angle (RADIANS, from the x-axis) is applied
                                          to the cell-1 output field before cell 2.
                                          90° crossed w.r.t. x-polarised input is
                                          np.pi/2.
        reset_polarisation str or None   'left_cp', 'right_cp' or None (default).
                                          Needs polariser_angle. If not None, the
                                          field reaching cell 2 (after the
                                          polariser_angle projection) is replaced
                                          with a FIXED left ([1,+i]/sqrt2) or right
                                          ([1,-i]/sqrt2) circular state carrying the
                                          same transverse (x,y) intensity, and the
                                          z-component is left untouched. This
                                          models an ideal quarter-wave plate at 45°
                                          to the linear polariser, which converts
                                          ANY input polarisation state to the SAME
                                          fixed circular output regardless of
                                          wavelength — the standard "polariser +
                                          QWP" circular-polariser combination of a
                                          wing-selector's P2+QWP stage.
                                          [ES4 FIX] Without
                                          polariser_angle it raises ValueError: no
                                          passive optic maps every state onto one
                                          circular state and keeps its intensity.
        strict          bool              As in calculate(): False downgrades an
                                          unknown p_dict key (in either cell) from
                                          an error to a warning.

    Returns:
        (cell1_outputs, cell2_outputs) -- two lists matching the outputs argument.

        cell1_outputs: results for cell 1 alone, normalised to E_in.

        cell2_outputs, polariser_angle=None and reset_polarisation=None (default):
            results for cell 2, normalised to E_out from cell 1
            (per-stage convention; combination rules below).

        cell2_outputs, polariser_angle set (with or without reset_polarisation):
            INSTRUMENT TOTALS — cell-2 outputs computed from the final field
            but normalised to the ORIGINAL E_in, so they already include the
            cell-1 attenuation and any polariser/QWP projection loss. No
            combination rule is needed (and per-stage normalisation would be
            0/0 in the crossed-polariser wings, where almost no light enters
            cell 2). Only field-derived outputs (S0, S1, S2, S3, Ix, Iy,
            I_P45, I_M45) are available in this mode; medium properties
            (alphaPlus/alphaMinus/alphaZ, n's, chi's, E_out) are not.

        Combination rules for the cascaded total:

          Intensity outputs (S0, Ix, Iy, I_P45, I_M45, Il, Ir):
            X_total = S0_cell1 * X_cell2

          Fractional Stokes parameters (S1, S2, S3 — already normalised by
          S0 inside get_spectra):
            X_total = X_cell2   (no prefactor needed)

          Medium properties (alphaPlus, alphaMinus, alphaZ):
            These are intrinsic to each cell; no simple combination rule
            applies.

          Where cell 1 lets no light through (its output intensity below the
          smallest normal float, ~2e-308, as at an opaque line centre), cell 2
          has no input to normalise by. Its per-stage field outputs (the
          intensities, S1-S3, psi, DoLP and S0_P) are 0 there, so the rules
          give a total of 0, as for an opaque filter. [ES4 FIX] They were NaN.

    Example usage:
        import numpy as np
        from elecsus import elecsus_methods as em

        detuning = np.linspace(-10000, 10000, 1000)
        p1 = {'Elem': 'Rb', 'Dline': 'D2', 'T': 20, 'lcell': 75e-3}
        p2 = {'Elem': 'Rb', 'Dline': 'D2', 'T': 80, 'lcell': 50e-3}

        (S0_1,), (S0_2,) = em.calculate_cascade(detuning, [1, 0, 0], p1, p2, outputs=['S0'])
        S0_combined = S0_1 * S0_2

        # Wing-selector topology (VAMOS-style): cell 1 between crossed linear
        # polarisers, then a QWP preparing fixed left-circular light for cell 2.
        (S0_1,), (S0_2,) = em.calculate_cascade(
            detuning, [1, 0, 0], p1, p2, outputs=['S0'],
            polariser_angle=np.pi / 2, reset_polarisation='left_cp')
    """
    _default_outputs = [
        'S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45',
        'alphaPlus', 'alphaMinus', 'alphaZ',
    ]
    _field_derived = list(_CASCADE_FIELD_OUTPUTS) + ['E_out']

    if reset_polarisation not in (None, 'left_cp', 'right_cp'):
        raise ValueError(
            f"reset_polarisation must be None, 'left_cp' or 'right_cp', got "
            f"{reset_polarisation!r}.")
    # [ES4 FIX] Only a polariser followed by a quarter-wave plate
    # maps every input onto one circular state. Without the polariser the reset
    # kept the full transverse intensity, which no passive optic can do.
    if reset_polarisation is not None and polariser_angle is None:
        raise ValueError(
            f"reset_polarisation={reset_polarisation!r} needs polariser_angle: it "
            "models a quarter-wave plate at 45 deg to a linear polariser between "
            "the cells, which passes the polariser's share of the light. Without "
            "the polariser no passive optic turns every polarisation into the "
            "same circular state. Set polariser_angle (radians from x), or leave "
            "reset_polarisation as None.")
    intermediate_optic = polariser_angle is not None   # a QWP needs one (above)

    if not intermediate_optic:
        outputs_list = list(outputs) if outputs is not None else _default_outputs
    else:
        outputs_list = list(outputs) if outputs is not None else _field_derived[:-1]
        bad = [k for k in outputs_list if k not in _field_derived]
        if bad:
            raise ValueError(
                f'outputs {bad} are not available with an intermediate polariser/QWP: '
                'cell-2 quantities are computed from the final field (instrument '
                f'totals), so only {_field_derived} are supported.')

    # Always retrieve E_out from cell 1 to feed into cell 2.
    # Append it to the request only if the user didn't ask for it explicitly.
    user_wants_E_out = 'E_out' in outputs_list
    outputs_1_internal = outputs_list if user_wants_E_out else outputs_list + ['E_out']

    raw_1 = spectra.get_spectra(detuning_range, E_in, p_dict_1,
                                outputs=outputs_1_internal, strict=strict)
    E_out_1 = raw_1[-1]
    cell1_results = list(raw_1) if user_wants_E_out else list(raw_1[:-1])

    if not intermediate_optic:
        cell2_results = list(spectra.get_spectra(detuning_range, E_out_1, p_dict_2,
                                                 outputs=outputs_list, strict=strict))
        # [ES4 FIX] Where cell 1 lets no light through, cell 2's per-stage field
        # outputs were 0/0 = NaN, so the documented total S0_cell1 * X_cell2 was
        # NaN where no light leaves the cascade. They are 0 there, and the
        # fractional Stokes parameters 0, as in the instrument totals below.
        E_1 = array(E_out_1)
        dark = ~spectra._lit(np.real((E_1 * E_1.conjugate()).sum(axis=0)))
        if dark.any():
            cell2_results = [np.where(dark, 0., value) if key in _PER_STAGE_FIELD_OUTPUTS
                             else value
                             for key, value in zip(outputs_list, cell2_results)]
        return cell1_results, cell2_results

    # ── Intermediate polariser / QWP path (MOF topology) ──────────────────────
    # Project the cell-1 output field through a linear polariser at
    # polariser_angle (radians from x). If reset_polarisation is also
    # set, replace the field with a fixed left/right circular state of the
    # same transverse intensity — see reset_polarisation above.
    # Then propagate through cell 2 and compute the requested quantities from
    # the final field, normalised to the ORIGINAL input intensity. Per-stage
    # normalisation is deliberately not used here: with crossed polarisers
    # essentially no light enters cell 2 in the far wings and stage-normalised
    # outputs would be 0/0 there.
    E_out_1 = array(E_out_1)
    c, s = np.cos(polariser_angle), np.sin(polariser_angle)
    E_mid = np.empty_like(E_out_1)
    E_mid[0] = E_out_1[0] * c * c + E_out_1[1] * c * s
    E_mid[1] = E_out_1[0] * c * s + E_out_1[1] * s * s
    E_mid[2] = E_out_1[2]

    if reset_polarisation is not None:
        I_xy = np.real(E_mid[0] * E_mid[0].conjugate() + E_mid[1] * E_mid[1].conjugate())
        amp = np.sqrt(I_xy / 2.)
        sign = 1j if reset_polarisation == 'left_cp' else -1j
        E_reset = np.empty_like(E_mid)
        E_reset[0] = amp
        E_reset[1] = amp * sign
        E_reset[2] = E_out_1[2]
        E_mid = E_reset

    [E_out_2] = spectra.get_spectra(detuning_range, E_mid, p_dict_2, outputs=['E_out'],
                                    strict=strict)
    E_out_2 = array(E_out_2)

    # original input intensity (scalar 3-vector or per-detuning (3, N) array)
    E_in_arr = array(E_in, dtype=complex)
    if E_in_arr.ndim == 1:
        I_in = float(np.real((E_in_arr * E_in_arr.conjugate()).sum()))
    else:
        I_in = np.real((E_in_arr * E_in_arr.conjugate()).sum(axis=0))

    cell2_results = [_totals_from_field(E_out_2, I_in, key) for key in outputs_list]
    return cell1_results, cell2_results


def _totals_from_field(E_field, I_in, key):
    """Compute one field-derived output from a (3, N) E-field, normalised to
    the original input intensity I_in (instrument-total convention).
    Fractional Stokes are normalised by the total S0 with a safe divide, so
    the crossed-polariser wings (S0 -> 0) return at most 1e-15, not NaN.

    [ES4 FIX] Computed by spectra.field_outputs, as get_spectra's
    are. This function had its own copy, which kept the 4.0.0 alphas' sign of S3:
    left-circular light gave S3 = +1 here and -1 from calculate()."""
    if key == 'E_out':
        return E_field
    if key not in _CASCADE_FIELD_OUTPUTS:
        raise ValueError(f'unknown field-derived output {key!r}')
    return spectra.field_outputs(E_field, I_in, S0_floor=1e-15)[key]


#: The outputs calculate_cascade returns as instrument totals.
_CASCADE_FIELD_OUTPUTS = ('S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45')

#: Cell 2's per-stage outputs that are normalised by the light it receives,
#: and so are 0 where it receives none (calculate_cascade without a polariser).
_PER_STAGE_FIELD_OUTPUTS = _CASCADE_FIELD_OUTPUTS + ('Il', 'Ir', 'psi', 'DoLP', 'S0_P')
# [ES4 END] calculate_cascade


#: The fit_algorithm names fit_data accepts, and the routine each one runs.
_FIT_ALGORITHMS = {'ML': 'ML', 'Marquardt-Levenberg': 'ML', 'LM': 'ML', 'leastsq': 'ML',
                   'DE': 'DE', 'differential_evolution': 'DE', 'RR': 'RR'}


def _fit_routine(fit_algorithm):
    """The routine ('ML', 'DE' or 'RR') that `fit_algorithm` names.

    [ES4 FIX] Any other name used to run random restart, so
    'differential evolution' (with a space) ran RR and dropped the DE bounds.
    It now raises, suggesting the closest name.

    [ES4] Simulated annealing ('SA') was removed: it had
    crashed on Linux and macOS since ES3 3.0.0, and differential evolution
    does the same job faster. Asking for it raises, pointing to DE.
    """
    if isinstance(fit_algorithm, str) and fit_algorithm in _FIT_ALGORITHMS:
        return _FIT_ALGORITHMS[fit_algorithm]
    if (isinstance(fit_algorithm, str) and fit_algorithm.strip().lower().replace(' ', '_')
            in ('sa', 'simulated_annealing')):
        raise ValueError(
            f"fit_algorithm={fit_algorithm!r}: simulated annealing has been removed from "
            "ElecSus 4. For a global fit use 'DE' (differential evolution), with "
            "[min, max] bounds in p_dict_bounds.")
    names = list(_FIT_ALGORITHMS)
    lower = [name.lower() for name in names]
    close = difflib.get_close_matches(str(fit_algorithm).lower(), lower, n=1, cutoff=0.6)
    hint = f" (did you mean '{names[lower.index(close[0])]}'?)" if close else ''
    raise ValueError(f"unknown fit_algorithm {fit_algorithm!r}{hint}. Use one of "
                     + ', '.join(repr(name) for name in names) + '.')


def _phase_shifts(phi, lo, hi):
    """The multiples k of pi that bring phi - k pi into [lo, hi]: the one
    nearest 0, then the nearest of the other parity, if the range holds one.
    If none does, the one that lands nearest the range, if that is within
    1e-4 rad (fit_data then clamps the phase, as lmfit would); else []."""
    # floats: an infinite bound gives an infinite end, which the clamp handles
    k_min, k_max = np.ceil((phi - hi) / np.pi), np.floor((phi - lo) / np.pi)
    if k_min <= k_max:
        k = min(max(0.0, k_min), k_max)
        other = [j for j in (k - 1, k + 1) if k_min <= j <= k_max]
        return [int(k)] + ([int(min(other, key=lambda j: (abs(j), j)))] if other else [])
    below = int(np.floor((phi - hi) / np.pi))   # lands just above hi
    above = below + 1                            # lands just below lo
    k = below if (phi - below * np.pi) - hi <= lo - (phi - above * np.pi) else above
    if min(abs(phi - k * np.pi - hi), abs(phi - k * np.pi - lo)) > 1e-4 + 1e-12:
        return []
    return [k]


def _fit_form(E_in, phase_range=None, x_range=None, y_range=None):
    """E_in in fit_data's form [E_x, [E_y, E_phase]], with E_phase in the
    range the fit allows, and E_x and E_y in theirs.

    [ES4 FIX] fit_data takes the input field as a real E_x, and E_y
    as an amplitude and a phase relative to E_x. The ML fit bounds the phase
    to `phase_range`: [0, ML.E_PHASE_MAX] unless the caller bounds E_phase in
    p_dict_bounds. Two inputs went wrong:

    - the Jones vector [Ex, Ey, 0] that calculate() takes raised TypeError;
    - a phase outside that range, such as right-circular light, was clamped
      to the bound by lmfit, so the fit modelled a different polarisation.

    A phase is moved into range by a multiple of pi, with the sign of E_x
    flipped for an odd multiple: (-a, b e^{i(phi-pi)}) = -(a, b e^{i phi}) is
    the same polarisation state. Within 1e-4 rad of the range (the gap below
    pi that the default range leaves out) it is clamped, as lmfit would. A
    fit-form E_in already in range is returned unchanged. When either
    amplitude is zero the relative phase means nothing and is set to 0 (or
    the nearest end of the range), so that a fit of the other amplitude
    explores linear states. (Moved here from the GUI's _pol_fit_setup, which
    uses it.)

    [ES4 FIX] `x_range` and `y_range` are the caller's bounds on
    E_x and E_y (None: unbounded), which lmfit would clamp the start to. The
    flip respects them: where E_x cannot change sign, E_y does, since
    (a, -b e^{i(phi-pi)}) is the same state too. A state that no form can
    hold inside all three ranges raises ValueError naming the bounds; so do
    a zero field and an E_in that is neither form. Infinite or very wide
    phase bounds work: the multiple of pi is computed, not searched for.
    """
    lo, hi = (0.0, ML.E_PHASE_MAX) if phase_range is None else phase_range
    x_lo, x_hi = (-np.inf, np.inf) if x_range is None else x_range
    y_lo, y_hi = (-np.inf, np.inf) if y_range is None else y_range
    form = ("E_in must be the fit form [E_x, [E_y, E_phase]] or a Jones vector "
            "[Ex, Ey, 0] as calculate() takes")
    fit_form = False
    try:
        fit_form = (len(E_in) == 2 and np.ndim(E_in[0]) == 0
                    and np.ndim(E_in[1]) == 1 and len(E_in[1]) == 2)
        if fit_form:
            raw = [E_in[0], E_in[1][0], E_in[1][1]]
            # float() would drop a numpy complex part with only a warning
            if any(np.iscomplexobj(c) and np.imag(c) != 0 for c in raw):
                raise ValueError
            values = [float(np.real(c)) for c in raw]
        else:
            values = np.asarray(E_in, dtype=complex)
            if values.shape != (3,):
                raise ValueError
    except (TypeError, ValueError):
        raise ValueError(f"{form}, with real fit-form components; got {E_in!r}."
                         if fit_form else f"{form}; got {E_in!r}.") from None
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{form}, with finite components; got {E_in!r}.")
    if not np.any(values[:2] if fit_form else values):
        raise ValueError(f"E_in is a zero field: there is no light to fit; got {E_in!r}.")
    if fit_form:
        a, b, phi = values
        if lo <= phi <= hi and x_lo <= a <= x_hi and y_lo <= b <= y_hi:
            return [a, [b, phi]]
    else:
        E = values
        if abs(E[2]) > 1e-9 * max(np.linalg.norm(E[:2]), 1e-300):
            raise ValueError(f"{form}: the input field of a beam along z has no "
                             f"z component; got {E_in!r}.")
        a, b = float(abs(E[0])), float(abs(E[1]))
        phi = float(np.angle(E[1]) - np.angle(E[0]))

    if a == 0 or b == 0:
        # no relative phase, and the sign of the other amplitude is a global
        # phase: any sign will do, and no shift is needed
        shifts = [0]
        phi = min(max(0.0, lo), hi)
        signs = {0: ((1, 1), (-1, 1), (1, -1), (-1, -1))}
    else:
        shifts = _phase_shifts(phi, lo, hi)
        if not shifts:
            raise ValueError(
                f"E_in has relative phase {phi:.6g} rad (mod pi), which "
                f"p_dict_bounds['E_phase'] = [{lo:.6g}, {hi:.6g}] cannot hold: "
                "widen the bounds, or start from a state inside them.")
        # an odd multiple of pi flips one sign (E_x first), an even one none or both
        signs = {0: ((1, 1), (-1, -1)), 1: ((-1, 1), (1, -1))}
    for k in shifts:
        for s_a, s_b in signs[k % 2]:
            if x_lo <= s_a * a <= x_hi and y_lo <= s_b * b <= y_hi:
                return [float(s_a * a + 0.0), [float(s_b * b + 0.0),
                                               float(min(max(phi - k * np.pi, lo), hi))]]
    bounds = ', '.join(f"p_dict_bounds['{key}'] = [{r[0]:.6g}, {r[1]:.6g}]"
                       for key, r in (('E_x', x_range), ('E_y', y_range),
                                      ('E_phase', phase_range)) if r is not None)
    if phase_range is None:
        bounds += f" and the fit's E_phase range [0, {hi:.6g}]"
    raise ValueError(
        f"E_in = {E_in!r} cannot be written as [E_x, [E_y, E_phase]] inside "
        f"{bounds}: changing the sign of E_x or E_y (with E_phase moved by pi) "
        "gives the same state, but none of its forms fits. Widen the bounds, or "
        "start from a state inside them.")


def _check_buffer_gas_fit(p_dict, p_dict_bools):
    """Raise ValueError if p_dict_bools varies a buffer-gas parameter that
    the data cannot fix, having no effect on the model or being degenerate
    with another: the fit would report a value and an error bar for it that
    the data never constrained.

    The model is on when BufferGas names a gas with data for the element and
    line and BufferPressure > 0 (spectra.calc_chi). Its broadening then
    replaces GammaBuf, and its pressure shift is added to shift, which can
    still vary. It depends on BufferPressure / BufferTemp_fill only (the
    density of the gas), so a fit can find the pressure but not both."""
    gas = p_dict.get('BufferGas')
    elem = p_dict.get('Elem', spectra.p_dict_defaults['Elem'])
    dline = p_dict.get('Dline', spectra.p_dict_defaults['Dline'])
    pressure = p_dict.get('BufferPressure', 0.0)
    model_on = (gas is not None and pressure > 0
                and (elem, dline, gas) in BUFFER_GAS_DATA)

    if p_dict_bools.get('BufferTemp_fill'):
        raise ValueError(
            "'BufferTemp_fill' cannot be fitted: the buffer-gas model depends only "
            "on BufferPressure / BufferTemp_fill (the density of the gas), so it is "
            "degenerate with BufferPressure. Fit 'BufferPressure', with "
            "'BufferTemp_fill' fixed at the temperature the cell was filled at.")
    if p_dict_bools.get('BufferPressure') and not model_on:
        raise ValueError(
            "'BufferPressure' can be fitted only with the buffer-gas model on: "
            "'BufferGas' set to a gas with data for this element and line, and a "
            f"starting 'BufferPressure' above 0 Torr (got BufferGas={gas!r}, "
            f"BufferPressure={pressure} for {elem} {dline}).")
    if p_dict_bools.get('GammaBuf') and model_on:
        raise ValueError(
            f"Cannot vary 'GammaBuf' while the buffer-gas model is on (BufferGas={gas!r}, "
            f"BufferPressure={pressure}): its broadening replaces GammaBuf. Vary "
            "'BufferPressure' to fit the broadening, or leave out 'BufferGas' to fit "
            "GammaBuf instead.")


def _fit_data_arrays(data):
    """data as [x, y], two 1-D float arrays of equal length, or ValueError.

    [ES4 FIX] The docstrings said data was "an Nx2 iterable";
    the routines read data[0] as x and data[1] as y, so an Nx2 array was
    fitted as two data points, with no error. (A 2x2 array is taken as
    [x, y], as before.)"""
    form = ("data must be [x, y]: the detunings and the measured values, two "
            "1-D sequences of equal length (2xN)")
    try:
        shape = np.shape(data)
    except ValueError:      # ragged: x and y of different lengths
        shape = None
    if shape is not None and len(shape) == 2 and shape[1] == 2 and shape[0] != 2:
        raise ValueError(f"{form}; got an array of shape {shape}, which is Nx2: "
                         "pass its transpose (data.T).")
    try:
        x, y = (np.asarray(column) for column in data)
    except (TypeError, ValueError):
        raise ValueError(f"{form}; got {type(data).__name__} of shape {shape}.") from None
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError(f"{form}; got x of shape {x.shape} and y of shape {y.shape}.")
    return [x, y]


def fit_data(data, p_dict, p_dict_bools, E_in=None, p_dict_bounds=None, data_type='S0', fit_algorithm='ML',
             *, strict=True, verbose=False):
    """
    Method to compare and fit experimental data to ElecSus.

    *** Example use cases can be found in /tests/test_fitting.py

    Arguments:
        data:					[x_array, y_array]: the detunings (MHz) and the measured values, two 1-D
                        sequences of equal length (2xN). [ES4 FIX] This said "an Nx2
                        iterable", which the routines read as two data points; anything but two 1-D
                        sequences of equal length now raises ValueError (for an Nx2 array: pass data.T).

        p_dict:					dictionary containing all the calculation (initial) parameters
        p_dict_bools:		dictionary with the same keys as p_dict, with Boolean values representing each parameter
                            that is to be varied in the fitting

    Options:
        E_in:			the initial electric field input: [E_x, [E_y, E_phase]], a real E_x and E_y as an
                        amplitude with a phase relative to E_x. [ES4] The Jones vector [Ex, Ey, 0] that
                        calculate() takes works too. A phase outside the range the fit allows is moved into
                        it by a multiple of pi, with the sign of E_x (or, where p_dict_bounds keeps E_x from
                        changing sign, of E_y) carrying the rest: the same polarisation state. That range is
                        [0, pi), or p_dict_bounds['E_phase'] for 'ML' and 'DE' when 'E_phase' is a key of
                        p_dict_bools; bounds on E_x and E_y apply the same way. A state that no form can
                        hold inside those bounds raises ValueError, as does a zero field.
                        Default: x-polarised, or E_x/E_y/E_phase from p_dict if it has them (ES3); [ES4]
                        p_dict must then have all three, or ValueError names the missing ones.
                        [ES4] Every output is normalised by the input intensity, so only the ratio E_y/E_x
                        and the phase E_phase are observable: p_dict_bools may vary E_x or E_y, not both
                        (ValueError). To fit a linear angle, vary the smaller amplitude, as the
                        GUI does; to fit an ellipse, vary one amplitude and E_phase.
        p_dict_bounds:	dictionary keyed by fit parameter. What each value means depends on fit_algorithm:
                        'ML' and 'DE': a [min, max] pair, the range the parameter can take.
                            Optional, except for when using 'differential_evolution' fitting method,
                            when bounds must be provided on fit parameters.
                            [ES4] Bounds reaching past a parameter's physical range (GammaBuf below 0,
                            an isotope fraction outside 0-100, ...: see param_validation.FIT_LIMITS)
                            are narrowed to it, with a warning.
                        'RR': the width of the random search around the starting value, one
                            number such as {'T': 30} (not a hard limit). Parameters without one keep
                            the routine's default width.
                        [ES4] An entry in the wrong form raises ValueError, naming the parameter.
        strict:         [ES4] p_dict keys are checked as calculate() checks them, and the keys of p_dict_bools
                        and p_dict_bounds must be fit parameters (param_validation.FIT_PARAMETERS). An
                        unknown key raises ValueError naming it; strict=False downgrades that to a warning
                        (the key is then ignored), as in calculate(). The range checks always apply.
        data_type:			Data type to fit experimental data to. Can be one of:
                                        'S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', ...
        verbose:				Boolean - more print statements provided as the program progresses.
                        [ES4] A real argument now; it was a module-level flag.

        fit_algorithm:	One of the following:
                            'ML',  standard Marquardt-Levenberg fitting ('LM', 'leastsq' and
                                   'Marquardt-Levenberg' are the same)
                            'RR', Random-restart. [ES4] Its starts are drawn inside the physical
                                   range, and a start whose downhill fit leaves that range is
                                   dropped with a warning (ParameterRangeError if every start fails).
                            'DE', Differential Evolution (or 'differential_evolution')
                        [ES4] Any other name raises ValueError; it used to run random restart.
                        [ES4] 'SA' (simulated annealing) was removed: use 'DE' for a global fit.

                        In principle the other methods supported by lmfit are possible, but these have
                        not so far been implemented here

        [ES4] strict and verbose are keyword-only, and a misspelt argument raises TypeError: fit_data
        used to accept and ignore any keyword (datatype='S1' fitted S0).

    Note on buffer gas fitting [ES4 NEW]:
        With BufferGas set and BufferPressure > 0, the buffer-gas model gives the
        Lorentzian broadening, which replaces GammaBuf, and a pressure shift, which
        is added to shift; both are evaluated at the temperature of the vapour in the
        beam: T, or DoppTemp when Constrain=False.
        In p_dict_bools:
          - BufferPressure may vary: the fit then finds the amount of gas in the cell.
            It needs BufferGas set to a gas with data for the element and line, and
            a starting BufferPressure above 0 Torr; otherwise ValueError.
          - shift may vary, as a frequency offset on top of the pressure shift.
          - GammaBuf may not (ValueError): the model replaces it.
          - BufferTemp_fill may not (ValueError): the model depends on
            BufferPressure / BufferTemp_fill only, so the two are degenerate. Keep it
            at the temperature the cell was filled at.

    Returns:
        best_params, RMS, result:
            best_params: dict of best-fit parameters with the fixed ones. [ES4] It holds p_dict keys
                only, so it can go straight back into calculate(); its E_x/E_y/E_phase,
                the fitted input field, make calculate() warn, as it takes the field from E_in.
            RMS: root-mean-square residual of the fit.
            result: the lmfit ModelResult (best_fit, fit_report(), parameter uncertainties).
    """

    # [ES4 FIX] An unknown algorithm name raises here: it ran RR.
    routine = _fit_routine(fit_algorithm)
    data = _fit_data_arrays(data)

    # [ES4 FIX] calculate() rejects a key it does not read, but
    # fit_data handed the dict to lmfit, which ignores unknown names: a
    # mistyped key ('Bfeild') left its parameter at the default, and the fit
    # converged on the wrong model in silence. Check the keys as calculate()
    # does, and the keys of p_dict_bools and p_dict_bounds against the fit
    # parameters (a bare KeyError from lmfit, and ignored, before).
    check_keys(p_dict, strict=strict)
    # With strict=False an unknown key has only warned: leave it out, as the
    # fit model does, so that no routine passes it on.
    p_dict = {key: value for key, value in p_dict.items()
              if key in KNOWN_KEYS or key in DEPRECATED_KEYS}
    p_dict_bools = check_fit_keys(p_dict_bools, 'p_dict_bools', strict=strict)
    if p_dict_bounds is not None:
        p_dict_bounds = check_fit_keys(p_dict_bounds, 'p_dict_bounds', strict=strict)

    # [ES4 FIX] The outputs are normalised by the input intensity,
    # so scaling E_x and E_y together changes nothing: a fit of both had a flat
    # direction, and returned arbitrary amplitudes with ~1e6 error bars.
    if p_dict_bools.get('E_x') and p_dict_bools.get('E_y'):
        raise ValueError(
            "p_dict_bools varies both 'E_x' and 'E_y', but only their ratio is "
            "observable: every output is normalised by the input intensity. Vary "
            "one of them (with 'E_phase' for an elliptical state), as the GUI does.")

    if verbose:
        print('Starting parameter dictionary:\n', p_dict)

    # Work on a copy so we don't mutate the caller's dict
    p_dict = dict(p_dict)

    # [ES4 FIX] Reject an unphysical starting value here, by name,
    # rather than from inside the optimiser. (Before E_in is built from
    # p_dict's E_x/E_y/E_phase, so that a bad one is named.)
    check_ranges(p_dict)

    # [ES4 FIX] ES3's field keys come as a set: E_x alone raised a
    # bare KeyError: 'E_y', and E_y/E_phase without E_x were ignored in silence.
    given = [key for key in ('E_x', 'E_y', 'E_phase') if key in p_dict]
    if 0 < len(given) < 3:
        missing = [key for key in ('E_x', 'E_y', 'E_phase') if key not in p_dict]
        raise ValueError(
            f"p_dict has {', '.join(repr(k) for k in given)} but not "
            f"{', '.join(repr(k) for k in missing)}: ES3's input field takes all three "
            "(E_x, and E_y with its phase E_phase relative to E_x). Give all three, "
            "or none and pass the field as E_in.")

    # [ES4 NEW] Which buffer-gas parameters can vary (after check_ranges, so
    # that BufferPressure is a number here).
    _check_buffer_gas_fit(p_dict, p_dict_bools)

    # [ES4 FIX] p_dict_bounds holds [min, max] pairs for ML
    # and DE, which are kept inside the range where the model is defined
    # (bounds past a physical limit are narrowed to it, with a warning), so no
    # trial value can stop the fit. For RR it holds search widths, which RR
    # checks; fit_data used never to pass them on.
    if p_dict_bounds is not None and routine in ('ML', 'DE'):
        p_dict_bounds = clip_bounds_to_limits(p_dict_bounds)

    if E_in is None:
        if 'E_x' in p_dict:
            # ES3 back-compat: electric field components embedded in parameter dict
            E_in = [p_dict['E_x'], [p_dict['E_y'], p_dict['E_phase']]]
        else:
            # Default: x-polarised input
            E_in = [1.0, [0.0, 0.0]]
    # [ES4 FIX] also takes calculate()'s Jones vector, and keeps
    # the phase inside the range the fit will allow: ML_fit's [0, E_PHASE_MAX],
    # or the caller's own bounds on E_phase, which ML_fit applies instead.
    # And E_x and E_y inside the caller's bounds on them: ML_fit
    # sets the bounds of every key of p_dict_bools, varied or not, and lmfit
    # clamps a start outside them.
    def caller_range(key):
        if (p_dict_bounds is not None and routine in ('ML', 'DE')
                and key in p_dict_bounds and key in p_dict_bools):
            return tuple(p_dict_bounds[key])
        return None
    field_ranges = [caller_range(key) for key in ('E_phase', 'E_x', 'E_y')]
    E_in = _fit_form(E_in, *field_ranges)

    # Fitting routines use lmfit, which requires all positional parameters of
    # fit_function to be in p_dict (or lmfit initialises them to -inf).
    p_dict.setdefault('E_x',     E_in[0])
    p_dict.setdefault('E_y',     E_in[1][0])
    p_dict.setdefault('E_phase', E_in[1][1])
    # [ES4 FIX] T, lcell and Bfield used to be missing here, so a
    # p_dict that calculate() accepts (one without Bfield, say) started the fit
    # at Bfield = -inf. Fill every one from the defaults calculate() uses.
    # Elem and Dline too: ML_fit read them, and raised KeyError.
    for key in ('Elem', 'Dline', 'T', 'lcell', 'Bfield', 'Btheta', 'Bphi', 'GammaBuf', 'shift'):
        p_dict.setdefault(key, spectra.p_dict_defaults[key])
    # [ES4 FIX] A DoppTemp not given (or None) is T, as in
    # calculate(). The fit model always passes one, at its start value, which
    # was 20 C when p_dict had none: with Constrain=False the Doppler width sat
    # at 20 C whatever T was. A fit of DoppTemp starts at T; otherwise nothing
    # decouples the two temperatures, so the fit keeps them together, as
    # calculate() does, and does not report the model's unused DoppTemp.
    report_dopptemp = True
    if p_dict.get('DoppTemp') is None:
        p_dict.pop('DoppTemp', None)
        if p_dict_bools.get('DoppTemp'):
            p_dict['DoppTemp'] = p_dict['T']
        else:
            report_dopptemp = False
            if not p_dict.get('Constrain', True):
                warn(NO_DOPPTEMP)
                p_dict['Constrain'] = True

    # Call different fitting routines
    if routine == 'ML':
        print('\nPerforming Marquardt-Levenberg fitting routine.')
        optParams, result = ML.ML_fit(data, E_in, p_dict, p_dict_bools, p_dict_bounds=p_dict_bounds,
                                      data_type=data_type)
        print('ML Fit completed')
    elif routine == 'DE':
        print('\nPerforming Differential Evolution fitting routine.')
        # Run with differential evolution
        optParams_DE, result = ML.ML_fit(data, E_in, p_dict, p_dict_bools, p_dict_bounds=p_dict_bounds,
                                         data_type=data_type, method='differential_evolution')
        # Then to get errors on parameters, run ML fit with optimised parameters
        print('DE fitting finished - rounding off with ML fit with optimised parameters...')
        # [ES4 FIX] The ML round-off takes the fit parameters DE
        # found and everything else from p_dict: lmfit's best_values reports
        # the element, line, buffer gas and Constrain at the fit function's
        # defaults, and holds none of the model settings. (They were patched
        # back one by one; random restart, which was not, fitted Rb D2.)
        optParams_DE = {**p_dict, **ML.fitted_values(optParams_DE)}
        # then do ML fit on the end to get error bars ...
        # [ES4 FIX] ... starting from DE's input field as well:
        # ML_fit starts E_x/E_y/E_phase from its E_in argument, and the
        # original E_in threw away what DE found for the polarisation.
        E_in_DE = _fit_form(ML.E_in_from(optParams_DE), *field_ranges)
        optParams, result = ML.ML_fit(data, E_in_DE, optParams_DE, p_dict_bools, p_dict_bounds=p_dict_bounds,
                                      data_type=data_type)
    else:  # routine == 'RR'
        print('\nPerforming fitting by Random-Restart hill climbing method.')
        # The more parameters to fit, the more evaluations we need to do.
        nparameters = 0
        for key in p_dict_bools:
            if p_dict_bools[key]: nparameters += 1

        nevaluations = no_evals = 2 ** (4 + nparameters)  # integer
        optParams, result = RR.RR_fit(data, E_in, p_dict, p_dict_bools, p_dict_bounds=p_dict_bounds,
                                      no_evals=nevaluations, data_type=data_type)

    # Add fixed parameters back to the optParams dictionary.
    # [ES4 FIX] with the default where p_dict has none, as the fit
    # used: lmfit reports a Constrain the DE round-off was given as 1.0, which
    # calculate() rejects (a switch must be True or False).
    for key in ['Elem', 'Dline', 'Constrain']:
        optParams[key] = p_dict.get(key, spectra.p_dict_defaults[key])
    # [ES4 NEW] BufferGas is an independent variable in the underlying lmfit Model, so
    # result.best_values reports the function's default ('') rather than the value
    # actually used during fitting - patch it back from the input p_dict, same as above.
    optParams['BufferGas'] = p_dict.get('BufferGas')
    # [ES4 FIX] report the model settings the fit actually used
    for key in ML.SETTING_KEYS:
        if key in p_dict:
            optParams[key] = p_dict[key]
    # [ES4 FIX] lmfit's best_values also holds the fit model's
    # independent variable 'output' (the data type), which is no p_dict key,
    # so calculate(X, E, optParams) raised and notebook 06 failed. Return
    # p_dict keys only, so the result goes straight back into calculate().
    optParams = {key: value for key, value in optParams.items() if key in KNOWN_KEYS}
    if not report_dopptemp:
        optParams.pop('DoppTemp', None)    # it follows T

    if verbose: print(result.fit_report())

    Spec = result.best_fit
    ydata = data[1]
    RMS = sqrt(((ydata - Spec) ** 2).sum() / float(len(ydata)))

    return optParams, RMS, result
