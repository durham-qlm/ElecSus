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

import os
import sys
import inspect
import warnings
from shutil import copyfile
import numpy as np
from numpy import arange, zeros, array, sqrt

# import elecsus modules
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"libs")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
from libs import spectra
from libs import MLFittingRoutine as ML
from libs import SAFittingRoutine as SA
from libs import RRFittingRoutine as RR

# Stop warnings about casting complex numbers
warnings.simplefilter("ignore")

verbose = False


def calculate(detuning_range, E_in=[1, 0, 0], p_dict={}, outputs=None):
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
                Btheta			float			radians	Angle B-field makes with the y-z plane
                Bphi				float			radians	Angle B-field makes with the x-z plane

                # Experimental parameters
                T		 			float			Celsius	Temperature used to calculate atomic number density
                GammaBuf   	float			MHz		Extra lorentzian broadening (usually from buffer gas
                                                                but can be any extra homogeneous broadening)
                BufferGas   		str			--		Optional string (e.g. 'He', 'N2') to dynamically model GammaBuf/shift
                BufferPressure   	float			Torr		Fill pressure of BufferGas to act on
                BufferTemp_fill   	float			K		Temperature cell was filled with BufferGas (default 293.15)
                shift	  			float			MHz		A global frequency shift of the atomic resonance frequencies
                DoppTemp   	float			Celsius	Temperature linked to the Doppler width (used for
                                                                independent Doppler width and number density)
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

    return spectra.get_spectra(detuning_range, E_in, p_dict, outputs)


# [ES4 BEGIN] calculate_cascade — two cascaded vapour cells (not in ES3)
def calculate_cascade(detuning_range, E_in=[1, 0, 0], p_dict_1={}, p_dict_2={},
                      outputs=None, polariser_angle=None, reset_polarisation=None):
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
                                          If not None, the field reaching cell 2
                                          (after the polariser_angle projection, if
                                          also set) is replaced with a FIXED
                                          left/right circular state carrying the
                                          same transverse (x,y) intensity, and the
                                          z-component is left untouched. This
                                          models an ideal quarter-wave plate at 45°
                                          to a preceding linear polariser, which
                                          converts ANY input polarisation state to
                                          the SAME fixed circular output regardless
                                          of wavelength — the standard "polariser +
                                          QWP" circular-polariser combination.
                                          Combine with polariser_angle to match a
                                          real wing-selector's P2+QWP stage; using
                                          it alone is a simplified/idealised model
                                          of any element with the same net effect
                                          acting directly on cell 1's raw output.

    Returns:
        (cell1_outputs, cell2_outputs) -- two lists matching the outputs argument.

        cell1_outputs: results for cell 1 alone, normalised to E_in.

        cell2_outputs, polariser_angle=None and reset_polarisation=None (default):
            results for cell 2, normalised to E_out from cell 1
            (per-stage convention; combination rules below).

        cell2_outputs, polariser_angle and/or reset_polarisation set:
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
    _field_derived = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45', 'E_out']

    if reset_polarisation not in (None, 'left_cp', 'right_cp'):
        raise ValueError(
            f"reset_polarisation must be None, 'left_cp' or 'right_cp', got "
            f"{reset_polarisation!r}.")
    intermediate_optic = polariser_angle is not None or reset_polarisation is not None

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

    raw_1 = spectra.get_spectra(detuning_range, E_in, p_dict_1, outputs=outputs_1_internal)
    E_out_1 = raw_1[-1]
    cell1_results = list(raw_1) if user_wants_E_out else list(raw_1[:-1])

    if not intermediate_optic:
        cell2_results = list(spectra.get_spectra(detuning_range, E_out_1, p_dict_2,
                                                 outputs=outputs_list))
        return cell1_results, cell2_results

    # ── Intermediate polariser / QWP path (MOF topology) ──────────────────────
    # Project the cell-1 output field through a linear polariser at
    # polariser_angle (radians from x) if set. If reset_polarisation is also
    # (or instead) set, replace the field with a fixed left/right circular
    # state of the same transverse intensity — see reset_polarisation above.
    # Then propagate through cell 2 and compute the requested quantities from
    # the final field, normalised to the ORIGINAL input intensity. Per-stage
    # normalisation is deliberately not used here: with crossed polarisers
    # essentially no light enters cell 2 in the far wings and stage-normalised
    # outputs would be 0/0 there.
    E_out_1 = array(E_out_1)
    if polariser_angle is not None:
        c, s = np.cos(polariser_angle), np.sin(polariser_angle)
        E_mid = np.empty_like(E_out_1)
        E_mid[0] = E_out_1[0] * c * c + E_out_1[1] * c * s
        E_mid[1] = E_out_1[0] * c * s + E_out_1[1] * s * s
        E_mid[2] = E_out_1[2]
    else:
        E_mid = E_out_1

    if reset_polarisation is not None:
        I_xy = np.real(E_mid[0] * E_mid[0].conjugate() + E_mid[1] * E_mid[1].conjugate())
        amp = np.sqrt(I_xy / 2.)
        sign = 1j if reset_polarisation == 'left_cp' else -1j
        E_reset = np.empty_like(E_mid)
        E_reset[0] = amp
        E_reset[1] = amp * sign
        E_reset[2] = E_out_1[2]
        E_mid = E_reset

    [E_out_2] = spectra.get_spectra(detuning_range, E_mid, p_dict_2, outputs=['E_out'])
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
    the crossed-polariser wings (S0 -> 0) return 0 rather than NaN."""
    from elecsus.libs import JonesMatrices as JM
    from elecsus.libs import BasisChanger as BC

    if key == 'E_out':
        return E_field

    S0 = np.real((E_field * E_field.conjugate()).sum(axis=0)) / I_in
    if key == 'S0':
        return S0

    Ex = np.array(JM.HorizPol_xy * E_field[:2])
    Ey = np.array(JM.VertPol_xy * E_field[:2])
    Ix = np.real((Ex * Ex.conjugate()).sum(axis=0)) / I_in
    Iy = np.real((Ey * Ey.conjugate()).sum(axis=0)) / I_in
    if key == 'Ix':
        return Ix
    if key == 'Iy':
        return Iy

    EP45 = np.array(JM.LPol_P45_xy * E_field[:2])
    EM45 = np.array(JM.LPol_M45_xy * E_field[:2])
    IP45 = np.real((EP45 * EP45.conjugate()).sum(axis=0)) / I_in
    IM45 = np.real((EM45 * EM45.conjugate()).sum(axis=0)) / I_in
    if key == 'I_P45':
        return IP45
    if key == 'I_M45':
        return IM45

    Elrz = BC.xyz_to_lrz(E_field)
    El = np.array(JM.CPol_L_lr * Elrz[:2])
    Er = np.array(JM.CPol_R_lr * Elrz[:2])
    Il = np.real((El * El.conjugate()).sum(axis=0)) / I_in
    Ir = np.real((Er * Er.conjugate()).sum(axis=0)) / I_in

    safe = np.where(S0 > 1e-15, S0, 1.)
    if key == 'S1':
        return (Ix - Iy) / safe
    if key == 'S2':
        return (IP45 - IM45) / safe
    if key == 'S3':
        return -1 * (Ir - Il) / safe
    raise ValueError(f'unknown field-derived output {key!r}')
# [ES4 END] calculate_cascade


def fit_data(data, p_dict, p_dict_bools, E_in=None, p_dict_bounds=None, data_type='S0', fit_algorithm='ML', **kw):
    """
    Method to compare and fit experimental data to ElecSus.

    *** Example use cases can be found in /tests/test_fitting.py

    Arguments:
        data:					an Nx2 iterable for the x and y data to be fitted

        p_dict:					dictionary containing all the calculation (initial) parameters
        p_dict_bools:		dictionary with the same keys as p_dict, with Boolean values representing each parameter
                            that is to be varied in the fitting

    Options:
        E_in:			the initial electric field input. See docstring for the spectra.py module for details.
        p_dict_bounds:	dictionary with the same keys as p_dict, with values that are pairs of min/max values that
                        each parameter can take.
                        Optional, except for when using 'differential_evolution' fitting method,
                        when bounds must be provided on fit parameters
        data_type:			Data type to fit experimental data to. Can be one of:
                                        'S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', ...
        verbose:				Boolean - more print statements provided as the program progresses

        fit_algorithm:	One of the following:
                            'ML',  standard Marquardt-Levenberg fitting
                            'RR', Random-restart
                            'SA', Simulated Annealing
                            'DE', Differential Evolution

                        In principle the other methods supported by lmfit are possible, but these have
                        not so far been implemented here

    Note on buffer gas fitting [ES4 NEW]:
        BufferGas, BufferPressure and BufferTemp_fill are fixed cell-fabrication
        properties (known from how the cell was filled) - pass them in p_dict as usual,
        but do not mark them True in p_dict_bools. When BufferGas is set (with
        BufferPressure > 0), GammaBuf/shift are derived from the Pitz buffer-gas model
        as a function of T, so marking GammaBuf/shift as varying in p_dict_bools
        alongside an active BufferGas model raises ValueError (see spectra.py's
        override precedence).
    """

    # The more parameters to fit, the more evaluations we need to do.
    nparameters = 0
    for key in p_dict_bools:
        if p_dict_bools[key]:
            nparameters += 1

    if verbose:
        print('Starting parameter dictionary:\n', p_dict)

    # Work on a copy so we don't mutate the caller's dict
    p_dict = dict(p_dict)

    # [ES4 NEW] BufferGas model derives GammaBuf/shift from T (see spectra.py); if a
    # caller also asks to vary GammaBuf/shift directly, those fit params would silently
    # have no effect on the model output, wasting iterations and misleading error bars.
    if p_dict.get('BufferGas') is not None and p_dict.get('BufferPressure', 0.0) > 0:
        if p_dict_bools.get('GammaBuf') or p_dict_bools.get('shift'):
            raise ValueError(
                "Cannot vary 'GammaBuf'/'shift' while BufferGas is active (BufferGas="
                f"{p_dict.get('BufferGas')!r}, BufferPressure={p_dict.get('BufferPressure')}). "
                "The buffer-gas model derives GammaBuf/shift from T, so they are not "
                "independent fit parameters in this mode - remove them from p_dict_bools."
            )

    if E_in is None:
        if 'E_x' in p_dict:
            # ES3 back-compat: electric field components embedded in parameter dict
            E_in = [p_dict['E_x'], [p_dict['E_y'], p_dict['E_phase']]]
        else:
            # Default: x-polarised input
            E_in = [1.0, [0.0, 0.0]]

    # Fitting routines use lmfit, which requires all positional parameters of
    # fit_function to be in p_dict (or lmfit initialises them to -inf).
    p_dict.setdefault('E_x',     E_in[0])
    p_dict.setdefault('E_y',     E_in[1][0] if hasattr(E_in[1], '__len__') else 0.0)
    p_dict.setdefault('E_phase', E_in[1][1] if hasattr(E_in[1], '__len__') else 0.0)
    p_dict.setdefault('Btheta',   0.0)
    p_dict.setdefault('Bphi',     0.0)
    p_dict.setdefault('GammaBuf', 0.0)
    p_dict.setdefault('shift',    0.0)

    # Call different fitting routines
    if fit_algorithm in ('ML', 'Marquardt-Levenberg', 'LM', 'leastsq'):
        print('\nPerfoming Marquardt-Levenberg fitting routine.')
        optParams, result = ML.ML_fit(data, E_in, p_dict, p_dict_bools, p_dict_bounds=p_dict_bounds,
                                      data_type=data_type)
        print('ML Fit completed')
    elif fit_algorithm == 'SA':
        print('\nPerforming fitting by simulated annealing.')
        nevaluations = 2 ** (8 + 2 * nparameters)
        optParams, result = SA.SA_fit(data, E_in, p_dict, p_dict_bools, data_type=data_type, no_evals=nevaluations)
    elif fit_algorithm in ('DE', 'differential_evolution'):
        print('\nPerfoming Differential Evolution fitting routine.')
        # Run with differential evolution
        optParams_DE, result = ML.ML_fit(data, E_in, p_dict, p_dict_bools, p_dict_bounds=p_dict_bounds,
                                         data_type=data_type, method='differential_evolution')
        # Then to get errors on parameters, run ML fit with optimised parameters
        print('DE fitting finished - rounding off with ML fit with optimised parameters...')
        try:
            for key in ['Elem', 'Dline', 'Constrain']:
                optParams_DE[key] = p_dict[key]
        except KeyError:
            optParams_DE[key] = spectra.p_dict_defaults[key]
        # [ES4 NEW] BufferGas has no entry in p_dict_defaults (it's legitimately absent
        # when no buffer gas is used), so patch it separately with a plain default of None.
        optParams_DE['BufferGas'] = p_dict.get('BufferGas')
        # then do ML fit on the end to get error bars ...
        optParams, result = ML.ML_fit(data, E_in, optParams_DE, p_dict_bools, p_dict_bounds=p_dict_bounds,
                                      data_type=data_type)
    else:
        print('\nPerforming fitting by Random-Restart hill climbing method.')
        # The more parameters to fit, the more evaluations we need to do.
        nparameters = 0
        for key in p_dict_bools:
            if p_dict_bools[key]: nparameters += 1

        nevaluations = no_evals = 2 ** (4 + nparameters)  # integer
        optParams, result = RR.RR_fit(data, E_in, p_dict, p_dict_bools,
                                      no_evals=nevaluations, data_type=data_type)

    # Add fixed parameters back to the optParams dictionary if they exist
    try:
        for key in ['Elem', 'Dline', 'Constrain']:
            optParams[key] = p_dict[key]
    except KeyError:
        pass
    # [ES4 NEW] BufferGas is an independent variable in the underlying lmfit Model, so
    # result.best_values reports the function's default ('') rather than the value
    # actually used during fitting - patch it back from the input p_dict, same as above.
    optParams['BufferGas'] = p_dict.get('BufferGas')

    if verbose: print(result.fit_report())

    Spec = result.best_fit
    ydata = data[1]
    RMS = sqrt(((ydata - Spec) ** 2).sum() / float(len(ydata)))

    return optParams, RMS, result
