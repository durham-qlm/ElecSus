# Copyright 2015 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams and I. G. Hughes.

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

Last updated 2018-02-19 JK
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)

import os
import sys
import inspect
import warnings
from shutil import copyfile
from numpy import arange, zeros, array, sqrt

# import elecsus modules
from .libs import spectra
from .libs import MLFittingRoutine as ML
from .libs import SAFittingRoutine as SA
from .libs import RRFittingRoutine as RR

# if os.name == 'posix':
# 	from time import time as timing #Timing for linux or apple
# else:
# 	from time import clock as timing #Timing for windows

# import time

# Stop warnings about casting complex numbers
warnings.simplefilter("ignore")

verbose = True


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
                Btheta			float			degrees	Angle B-field makes with the y-z plane
                Bphi				float			degrees	Angle B-field makes with the x-z plane

                # Experimental parameters
                T		 			float			Celsius	Temperature used to calculate atomic number density
                GammaBuf   	float			MHz		Extra lorentzian broadening (usually from buffer gas
                                                                but can be any extra homogeneous broadening)
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


def fit_data(data, p_dict, p_dict_bools, E_in=None, p_dict_bounds=None, data_type='S0', fit_algorithm='ML', **kw):
    """
    Method to compare and fit experimental data to ElecSus.

    *** Example use cases can be found in /tests/fitting_tests.py

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
    """

    # The more parameters to fit, the more evaluations we need to do.
    nparameters = 0
    for key in p_dict_bools:
        if p_dict_bools[key]:
            nparameters += 1

    if verbose:
        print('Starting parameter dictionary:\n', p_dict)

    if E_in is None:
        try:
            E_in = [p_dict['E_x'], [p_dict['E_y'], p_dict['E_phase']]]
        except:
            # E_in not in p_dict or specified otherwise...
            raise

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

    if verbose: print(result.fit_report())

    Spec = result.best_fit
    ydata = data[1]
    RMS = sqrt(((ydata - Spec) ** 2).sum() / float(len(ydata)))

    return optParams, RMS, result
