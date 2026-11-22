# Copyright 2014-2019 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
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

Random restart fitting routine. 

Fit by taking a random sample around parameters and then
fit using Marquardt-Levenberg.

Complete rebuild of the original RR fitting module now using lmfit

Author: JK
Last updated 2018-02-21 MAZ
"""
# py 2.7 compatibility


import numpy as np
import matplotlib.pyplot as plt
import warnings
import sys
import copy

from multiprocessing import Pool

from . import MLFittingRoutine as ML
import lmfit as lm

from .spectra import get_spectra
from .param_validation import FIT_LIMITS, ParameterRangeError, search_widths, warn

p_dict_bounds_default = {'lcell':1e-3,'Bfield':100., 'T':20., 
							'GammaBuf':20., 'shift':100.,
							# Polarisation of light
							'E_x':0.05, 'E_y':0.05, 'E_phase':0.01, 
							# B-field angle w.r.t. light k-vector
							'Btheta':10*3.14/180, 'Bphi':10*3.14/180,
							'DoppTemp':20.,
							'rb85frac':1, 'K40frac':1, 'K41frac':1,
							# [ES4 NEW] Torr. fit_data refuses to vary
							# BufferTemp_fill, so it has no width.
							'BufferPressure':10.,
							}
	
def evaluate(args):
	data = args[0]
	E_in = args[1] 
	p_dict = args[2]
	p_dict_bools = args[3]
	data_type = args[4]
	try:
		best_params, result = ML.ML_fit(data, E_in, p_dict, p_dict_bools, data_type)
	except ParameterRangeError as error:
		# [ES4 FIX] A downhill fit that steps outside the physical
		# range (leastsq tried GammaBuf = -474 MHz in ES3's fitting demo) is a
		# failed start. Its error stopped the whole search: Pool.map_async
		# raises the first worker exception once every start has run.
		return np.inf, p_dict, str(error.__cause__ or error)
	#print 'Eval_ML COmplete'

	# returns reduced chi-squared value, best fit parameters, and
	# the error that stopped a failed start, or None
	return result.redchi, best_params, None

def RR_fit(data,E_in,p_dict,p_dict_bools,p_dict_bounds=None,no_evals=None,data_type='S0',verbose=False):
	"""
	Random restart fitting method.

	data:					[x_array, y_array]: two 1-D sequences of equal length (2xN), the x and y data
							to be fitted. [ES4 FIX] This said "an Nx2 iterable": data[0] is x, data[1] is y.
	E_in:					the initial electric field input. See docstring for the spectra.py module for details.

	no_evals:			The number of randomly-selected start points for downhill fitting. Defaults to 2**(3+2*nFitParams) where nFitParams is
									the number of varying fit parameters
	
	p_dict:					dictionary containing all the calculation (initial) parameters
	p_dict_bools:		dictionary with the same keys as p_dict, with Boolean values representing each parameter that is to be varied in the fitting
	p_dict_bounds:	dictionary keyed by fit parameter, with values that are the width of the random search around each starting
									value: one number, e.g. {'T': 30}. Parameters without one keep their width in p_dict_bounds_default.
									NOTE: this works differently to p_dict_bounds in the ML fitting methods, where each value is a [min, max] pair.
									In RR fitting, the widths select the range in parameter space that is randomly explored as starting parameters
									for a downhill fit, rather than being strict bounds on the fit parameters.
									[ES4] Each start is drawn uniformly from the part of [start - width, start + width] that lies inside
									the physical range (param_validation.FIT_LIMITS).

	data_type:			Data type to fit experimental data to. Can be one of:
									'S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', ...
	verbose:				Boolean - more print statements provided as the program progresses

	[ES4] The downhill fits are not bounded. A start whose fit steps outside the physical range is dropped, with a
	warning; if every start fails that way, ParameterRangeError is raised.

	"""

	# [ES4 FIX] Widths given for only some of the fit parameters
	# used to raise KeyError for the rest, and ES3's one-element form [30] a
	# TypeError. Each width is now checked, and the rest keep their defaults.
	p_dict_bounds = {**p_dict_bounds_default, **search_widths(p_dict_bounds or {})}

	print('Starting Random Restart Fitting Routine')
	x = np.array(data[0])
	y = np.array(data[1])

	p_dict['E_x'] = E_in[0]
	p_dict['E_y'] = E_in[1][0]
	p_dict['E_phase'] = E_in[1][1]

	
	# count number of fit parameters
	nFitParams = 0
	for key in p_dict_bools:
		if p_dict_bools[key]: nFitParams += 1
	
	# default number of iterations based on number of fit parameters
	if no_evals == None:
		no_evals = nFitParams**2 + 5 # 2**(3+2*nFitParams)
		
	# Create random array of starting parameters based on parameter ranges given in p_dict range dictionary
	# Scattered uniformly over the parameter space
	
	#clone the parameter dictionary
	p_dict_list = []
	for i in range(no_evals):
		p_dict_list.append(copy.deepcopy(p_dict))
		
	# [ES4 FIX] Each start is drawn from the part of start +- width
	# that lies inside the physical range (param_validation.FIT_LIMITS). A start
	# outside it, such as GammaBuf below 0 from a start of 0, made its downhill
	# fit raise. Where the range does not cut the interval the draw is the same
	# as before: uniform on [start - width, start + width].
	for key in p_dict_bools:
		if p_dict_bools[key]==True:
			start_vals = p_dict[key]
			low_limit, high_limit = FIT_LIMITS.get(key, (-np.inf, np.inf))
			low = max(start_vals - p_dict_bounds[key], low_limit)
			high = min(start_vals + p_dict_bounds[key], high_limit)
			if low > high:
				raise ParameterRangeError(
					f"'{key}' = {start_vals} lies further than its search width "
					f"{p_dict_bounds[key]} outside its physical range "
					f"[{low_limit}, {high_limit}].")
			for i in range(len(p_dict_list)):
				p_dict_list[i][key] = np.random.uniform(low, high)
				
	if verbose:
		print('List of initial parameter dictionaries:')
		for pd in p_dict_list:
			print(pd)
		#print p_dict_list
		print('\n\n')
	
	#Do parallel ML fitting by utilising multiple cores
	po = Pool() # Pool() uses all cores, Pool(3) uses 3 cores for example.
	
	# [ES4 FIX] each start with its own input field: ML_fit starts
	# E_x/E_y/E_phase from its E_in argument, so passing the original E_in to
	# every start undid their scatter when the polarisation is fitted.
	args_list = [(data, ML.E_in_from(p_dict_list[k]), p_dict_list[k], p_dict_bools, data_type)
				 for k in range(no_evals)]
	
	Res = po.map_async(evaluate,args_list)
	result = Res.get()
	po.close()
	po.join()
	
	if verbose: print('RR calculation complete')

	#Find best fit
	# [ES4 FIX] The starts whose downhill fit left the physical
	# range are dropped, with a warning; if none is left, say so.
	failures = [error for cost, best, error in result if error is not None]
	result = [(cost, best) for cost, best, error in result if error is None]
	if not result:
		raise ParameterRangeError(
			f"Every one of the {len(failures)} random-restart starts failed: each "
			"downhill fit stepped outside the physical range. Start closer to the "
			"data, narrow the search widths in p_dict_bounds, or fit with 'ML' or "
			"'DE' and [min, max] bounds, which keep the fit inside the range. The "
			f"first failure: {failures[0]}")
	if failures:
		warn(f"{len(failures)} of the {len(failures) + len(result)} random-restart "
			 "starts failed and were dropped: their downhill fits stepped outside "
			 f"the physical range. The first failure: {failures[0]}")
	lineMin = np.argmin([cost for cost, best in result]) ## pick the fit with the lowest cost value

	best_values = result[lineMin][1] # best parameter dictionary
	if verbose:
		print('\n\n\n')
		print(best_values)

	# [ES4 FIX] the fit parameters only: best_values reports the
	# element, line, buffer gas and Constrain at the fit function's defaults
	p_dict_best = copy.deepcopy(p_dict)
	p_dict_best.update(ML.fitted_values(best_values))
	
	# Finally run the ML fitting one more time, using the best parameters 
	# (so we get the final_result object, which cannot be pickled and therefore isn't supported in multiprocessing)
	# [ES4 FIX] from the best start's input field, not the original
	best_values, final_result = ML.ML_fit(data, ML.E_in_from(p_dict_best), p_dict_best, p_dict_bools, data_type)
	
	# return best fit parameters, and the lmfit result object
	return best_values, final_result

def test_fit():
	p_dict = {'Elem':'Rb','Dline':'D2','T':80.,'lcell':2e-3,'Bfield':600.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':0.,'shift':0.}
	
	# only need to specify parameters that are varied
	p_dict_bools = {'T':True,'Bfield':True,'E_x':True}
	p_dict_bounds = {'T':10,'Bfield':100,'E_x':0.01}
	
	E_in = np.array([0.7,0.7,0])
	E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]
	
	print(E_in_angle)
	
	x = np.linspace(-10000,10000,100)
	[y] = get_spectra(x,E_in,p_dict,outputs=['S1']) + np.random.randn(len(x))*0.015
		
	data = [x,y.real]
	
	best_params, result = RR_fit(data, E_in_angle, p_dict, p_dict_bools, p_dict_bounds, no_evals = 8, data_type='S1')
	report = result.fit_report()
	fit = result.best_fit
	
	print(report)
	plt.plot(x,y,'ko')
	plt.plot(x,fit,'r-',lw=2)
	plt.show()
	
if __name__ == '__main__':
	test_fit()
	
