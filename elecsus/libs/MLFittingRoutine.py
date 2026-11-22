# Copyright 2014 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
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
Marquardt-Levenberg fit module

Overhauled October 2016 to use the lmfit module 
(https://lmfit.github.io/lmfit-py/index.html; pip install lmfit; also in the enthought package manager index)
 instead of curve_fit, which makes it much more elegant and easy
 to fit only the selected parameters, and include bounds on parameters
 
Author: JK

differential evolution needs lmfit version >= 0.9.3

Last updated 2018-07-04 MAZ
"""
# py 2.7 compatibility



import numpy as np
import lmfit as lm

from .spectra import get_spectra
from .param_validation import FIT_PARAMETERS, ParameterRangeError

# [ES4 FIX] p_dict keys that choose how light propagates through the cell rather
# than describing the vapour. None of them is a fit parameter, but every fit
# evaluation must pass them to get_spectra.
PROPAGATION_KEYS = ('use_non_orthogonal', 'n_window', 'wavevector_dirn')

# [ES4 FIX] The model settings that reach fit_function by closure,
# because lmfit cannot carry them: the propagation settings, and whether the
# ground-state populations carry the Boltzmann factor, which every fit used
# to apply whatever p_dict said.
SETTING_KEYS = PROPAGATION_KEYS + ('BoltzmannFactor',)

# The largest relative phase of E_y the fit allows: E_phase is bounded to
# [0, E_PHASE_MAX] (elecsus_methods._fit_form keeps a start value inside).
E_PHASE_MAX = np.pi - 1e-4


def fitted_values(best_values):
	"""The fit parameters in lmfit's best_values, without the rest.

	best_values also holds the fit function's independent variables, with
	their signature defaults rather than the values the fit used: Elem='Rb',
	Dline='D2', BufferGas='', and Constrain as a float. A p_dict updated with
	all of it fitted a different model: random restart's last fit
	ran Rb D2 without the buffer gas, whatever the caller asked for."""
	return {key: value for key, value in best_values.items() if key in FIT_PARAMETERS}

def E_in_from(p_dict):
	"""The input field that p_dict's E_x, E_y and E_phase describe, in ML_fit's
	E_in form. ML_fit starts those three parameters from its E_in argument, so
	a fit that should start from a p_dict's field must pass this."""
	return [p_dict['E_x'], [p_dict['E_y'], p_dict['E_phase']]]

def fit_function(x,E_x,E_y,E_phase,T,lcell,Bfield,Btheta,Bphi,GammaBuf,shift,
							DoppTemp=20,rb85frac=72.17,K40frac=0.01,K41frac=6.73,
							Elem='Rb',Dline='D2',Constrain=True,
							BufferGas='',BufferPressure=0.0,BufferTemp_fill=293.15,
							output='S0', verbose=False, settings=None):
	"""
	Fit function that is used by lmfit. Essentially just a wrapper for the get_spectra method,
	but the arguments are specified explicitly rather than being inside a parameter dictionary.

	Also allows for the polarisation components to be fitted

	# [ES4 NEW] BufferGas / BufferPressure / BufferTemp_fill are forwarded to
	# get_spectra so the dynamic buffer-gas model (see BufferGasData.py) can
	# compute the broadening and pressure shift from T during fitting.
	# BufferPressure can be a fit parameter; BufferTemp_fill cannot, as the
	# model depends on their ratio only (fit_data checks both).
	# BufferGas defaults to '' rather than None: lmfit's make_params() reinterprets any
	# independent variable whose *signature* default is None/bool as a fittable numeric
	# Parameter if it's also passed as a kwarg, which crashes on a non-numeric gas name
	# like 'Ne'. An empty-string default keeps BufferGas classified as a pure
	# independent variable; it's normalised back to None below.

	# [ES4 FIX] settings: dict of the SETTING_KEYS settings (non-orthogonal
	# solver, window index, beam direction, Boltzmann factor), merged into the
	# p_dict below. They are bool/None-valued, so they cannot go through lmfit
	# at all (see above); ML_fit passes them by closure. Before this, fit_data
	# silently fitted with the orthogonal model whatever use_non_orthogonal
	# said, and with the Boltzmann factor whatever BoltzmannFactor said.

	For explanation of parameters, see the documentation for get_spectra (in spectra.py)
	"""

	if verbose:
		print(('Parameters: ', Bfield, T, lcell, E_x, E_y, E_phase, Btheta, Bphi, GammaBuf, shift, DoppTemp, rb85frac))

	# Ex / Ey separated to allow for fitting polarisation
	E_in = np.array([E_x,E_y*np.exp(1.j*E_phase),0.])

	#print E_in

	#reconstruct parameter dictionary from arguments
	p_dict = {'Elem':Elem,'Dline':Dline,'T':T,'lcell':lcell,
			'Bfield':Bfield,'Btheta':Btheta,'Bphi':Bphi,'GammaBuf':GammaBuf,
			'shift':shift,'DoppTemp':DoppTemp,'Constrain':Constrain,
			'rb85frac':rb85frac,'K40frac':K40frac,'K41frac':K41frac,
			'BufferGas':(BufferGas or None),'BufferPressure':BufferPressure,
			'BufferTemp_fill':BufferTemp_fill}
	if settings:
		p_dict.update(settings)

	#for key in p_dict.keys():
	#	print key, p_dict[key]
	outputs = [output]

	y_out = get_spectra(x,E_in,p_dict,outputs)[0].real

	return y_out

def ML_fit(data,E_in,p_dict,p_dict_bools,data_type='S0',p_dict_bounds=None,method='leastsq',verbose=False):
	"""
	Main fitting method.
	
	*** Example use cases can be found in the /tests/test_fitting.py file ***

	data:					[x_array, y_array]: two 1-D sequences of equal length (2xN), the x and y data
							to be fitted. [ES4 FIX] This said "an Nx2 iterable": data[0] is x, data[1] is y.
	E_in:					the initial electric field input. See docstring for the spectra.py module for details.

	p_dict:					dictionary containing all the calculation (initial) parameters
	p_dict_bools:		dictionary with the same keys as p_dict, with Boolean values representing each parameter that is to be varied in the fitting
	p_dict_bounds:	dictionary with the same keys as p_dict, with values that are pairs of min/max values that each parameter can take.
									Optional, except for when using 'differential_evolution' fitting method, when bounds must be provided on fit parameters
	
	data_type:			Data type to fit experimental data to. Can be one of:
									'S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', ...
	verbose:				Boolean - more print statements provided as the program progresses
	method:				passed to lmfit, the fitting algorithm that is used. Can be anything that is supported by lmfit, which is currently (as of 2017-June-01):	
									'leastsq': Levenberg-Marquardt (default)
									'least_squares': Least-Squares minimization, using Trust Region Reflective method by default
									'differential_evolution': differential evolution
									'brute': brute force method
									'nelder': Nelder-Mead
									'lbfgsb': L-BFGS-B
									'powell': Powell
									'cg': Conjugate-Gradient
									'newton': Newton-Congugate-Gradient
									'cobyla': Cobyla
									'tnc': Truncate Newton
									'trust-ncg': Trust Newton-Congugate-Gradient
									'dogleg': Dogleg
									'slsqp': Sequential Linear Squares Programming
										see https://lmfit.github.io/lmfit-py/fitting.html for more information

	"""
	
	x = np.array(data[0])
	y = np.array(data[1])
	
	# Non-numeric arguments to pass to fitting function
	kwargz = {'Elem':p_dict['Elem'],'Dline':p_dict['Dline']}
	if 'Constrain' in list(p_dict.keys()):
		kwargz['Constrain'] = p_dict['Constrain']
	else:
		kwargz['Constrain'] = True
	# [ES4 NEW] BufferGas is a fixed, non-numeric cell-fabrication property (which gas
	# was used to fill the cell) - never a fit parameter, so it goes through kwargz
	# alongside Elem/Dline/Constrain rather than through lmfit's Parameters.
	kwargz['BufferGas'] = p_dict.get('BufferGas')
	kwargz['output'] = data_type

	# [ES4 MOD] lmfit ≥1.3 compatibility fix.  lmfit ≥1.3 treats non-numeric defaults as
	# independent variables and also owns a 'verbose' kwarg on Model.fit() itself, so
	# 'verbose' can't be passed via **kwargz without being consumed by Model.fit() rather
	# than forwarded to fit_function.  Wrap fit_function to drop verbose from its signature
	# so lmfit never sees it.
	# The model settings are captured by closure for the same reason.
	settings = {key: p_dict[key] for key in SETTING_KEYS if key in p_dict}
	def _fit_fn(x,E_x,E_y,E_phase,T,lcell,Bfield,Btheta,Bphi,GammaBuf,shift,
				DoppTemp=20,rb85frac=72.17,K40frac=0.01,K41frac=6.73,
				Elem='Rb',Dline='D2',Constrain=True,
				BufferGas='',BufferPressure=0.0,BufferTemp_fill=293.15,output='S0'):
		return fit_function(x,E_x,E_y,E_phase,T,lcell,Bfield,Btheta,Bphi,GammaBuf,shift,
							DoppTemp,rb85frac,K40frac,K41frac,Elem,Dline,Constrain,
							BufferGas,BufferPressure,BufferTemp_fill,output,False,
							settings)

	model = lm.Model(_fit_fn)
	params = model.make_params(**p_dict)
	
	params['E_x'].value = E_in[0]
	params['E_y'].value = E_in[1][0]
	params['E_phase'].value = E_in[1][1]
	params['E_phase'].min = 0
	params['E_phase'].max = E_PHASE_MAX
	
	# Turn off all parameters varying by default, unless specified in p_dict_bools
	allkeys = params.valuesdict()
	for key in allkeys:
		params[key].vary = False
	
	# Turn on fitting parameters as specified in p_dict_bools
	for key in p_dict_bools:
		params[key].vary = p_dict_bools[key]
		
		if p_dict_bounds is not None:
			if key in p_dict_bounds:
				params[key].min = p_dict_bounds[key][0]
				params[key].max = p_dict_bounds[key][1]
	
	if verbose: print(params)

	try:
		result = model.fit(y, x=x, params=params, method=method, **kwargz)
	except (ParameterRangeError, RuntimeError) as err:
		# [ES4 FIX] The optimiser tried a value outside the range
		# where the model is defined. scipy's differential evolution re-raises
		# any ValueError from the objective as a RuntimeError about its
		# "map-like callable", which named neither the parameter nor the cause.
		cause = err if isinstance(err, ParameterRangeError) else err.__cause__
		if not isinstance(cause, ParameterRangeError):
			raise
		raise ParameterRangeError(
			f"The fit ({method}) stopped at a trial value outside the physical "
			f"range: {cause} Keep the fit inside that range with p_dict_bounds."
			) from cause

	return result.best_values, result
