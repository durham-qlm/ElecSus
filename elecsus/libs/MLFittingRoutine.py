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
from __future__ import (division, print_function, absolute_import)



import matplotlib.pyplot as plt

import numpy as np
import lmfit as lm

from spectra import get_spectra

import time

def fit_function(x,E_x,E_y,E_phase,T,lcell,Bfield,Btheta,Bphi,GammaBuf,shift,
							DoppTemp=20,rb85frac=72.17,K40frac=0.01,K41frac=6.73,
							Elem='Rb',Dline='D2',Constrain=True,bwf_precision='high',
							laserPower=0,laserWaist=5e-3,output='S0', verbose=False):
	"""
	Fit function that is used by lmfit. Essentially just a wrapper for the get_spectra method,
	but the arguments are specified explicitly rather than being inside a parameter dictionary.

	Also allows for the polarisation components to be fitted

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
			'laserPower': laserPower, 'laserWaist': laserWaist, 'bwf_precision': bwf_precision,
			'rb85frac':rb85frac,'K40frac':K40frac,'K41frac':K41frac}

	#for key in p_dict.keys():
	#	print key, p_dict[key]
	outputs = [output]

	y_out = get_spectra(x,E_in,p_dict,outputs)[0].real

	return y_out

def ML_fit(data,E_in,p_dict,p_dict_bools,data_type='S0',p_dict_bounds=None,method='leastsq',verbose=False):
	"""
	Main fitting method.

	*** Example use cases can be found in the /tests/fitting_tests.py file ***

	data:					an Nx2 iterable for the x and y data to be fitted
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
	kwargz['output'] = data_type
	if 'bwf_precision' in list(p_dict.keys()):
		kwargz['bwf_precision'] = p_dict['bwf_precision']

	model = lm.Model(fit_function)
	params = model.make_params(**p_dict)

	params['E_x'].value = E_in[0]
	params['E_y'].value = E_in[1][0]
	params['E_phase'].value = E_in[1][1]
	params['E_phase'].min = 0
	params['E_phase'].max = np.pi-1e-4
	params['laserPower'].min = 0
	params['laserWaist'].min = 1e-12

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

	result = model.fit(y, x=x, params=params, method=method, **kwargz)

	return result.best_values, result