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

Random restart fitting routine. 

Fit by taking a random sample around parameters and then
fit using Marquardt-Levenberg.

Complete rebuild of the original RR fitting module now using lmfit

Author: JK
2016-10-17
"""

import numpy as np
import matplotlib.pyplot as plt
import warnings
import sys
import copy

import psutil
from multiprocessing import Pool

import MLFittingRoutine as ML
import lmfit as lm

from spectra import get_spectra

p_dict_bounds_default = {'lcell':1e-3,'Bfield':100., 'T':20., 
							'GammaBuf':20., 'shift':100.,
							# Polarisation of light
							'theta0':10., 'E_x':0.05, 'E_y':0.05, 'E_phase':0.01, 
							# B-field angle w.r.t. light k-vector
							'Btheta':10*3.14/180, 'Bphi':10*3.14/180,
							'DoppTemp':20.,
							'rb85frac':1, 'K40frac':1, 'K41frac':1,
							}
	
def evaluate(args):
	data = args[0]
	E_in = args[1] 
	p_dict = args[2]
	p_dict_bools = args[3]
	data_type = args[4]
	best_params, result = ML.ML_fit(data, E_in, p_dict, p_dict_bools, data_type)
	#print 'Eval_ML COmplete'

	# returns reduced chi-squared value and best fit parameters
	return result.redchi, best_params #, result

def RR_fit(data,E_in,p_dict,p_dict_bools,p_dict_bounds=None,no_evals=None,data_type='S0',verbose=False):
	"""
	Random restart fitting method.

	data:					an Nx2 iterable for the x and y data to be fitted
	E_in:					the initial electric field input. See docstring for the spectra.py module for details.

	no_evals:			The number of randomly-selected start points for downhill fitting. Defaults to 2**(3+2*nFitParams) where nFitParams is
									the number of varying fit parameters
	
	p_dict:					dictionary containing all the calculation (initial) parameters
	p_dict_bools:		dictionary with the same keys as p_dict, with Boolean values representing each parameter that is to be varied in the fitting
	p_dict_bounds:	dictionary with the same keys as p_dict, with values that are pairs of min/max values that each parameter can take.
									NOTE: this works slightly differently to p_dict_bounds in the other fitting methods. In RR fitting, the bounds
									select the range in parameter space that is randomly explored as starting parameters for a downhill fit, rather than being
									strict bounds on the fit parameters.
	
	data_type:			Data type to fit experimental data to. Can be one of:
									'S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', ...
	verbose:				Boolean - more print statements provided as the program progresses

	"""

	if p_dict_bounds is None:
		p_dict_bounds = p_dict_bounds_default

	print 'Starting Random Restart Fitting Routine'
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
		no_evals = 2**(3+2*nFitParams)
		
	# Create random array of starting parameters based on parameter ranges given in p_dict range dictionary
	# Scattered uniformly over the parameter space
	
	#clone the parameter dictionary
	p_dict_list = []
	for i in range(no_evals):
		p_dict_list.append(copy.deepcopy(p_dict))
		
	for key in p_dict_bools:
		if p_dict_bools[key]==True:
			start_vals = p_dict[key]
			#print start_vals
			for i in range(len(p_dict_list)):
				p_dict_list[i][key] = start_vals + np.random.uniform(-1,1) * p_dict_bounds[key]
				
	if verbose:
		print 'List of initial parameter dictionaries:'
		for pd in p_dict_list:
			print pd
		#print p_dict_list
		print '\n\n'
	
	#Do parallel ML fitting by utilising multiple cores
	po = Pool() # Pool() uses all cores, Pool(3) uses 3 cores for example.
	
	## use lower process priority so computer is still responsive while calculating!!
	parent = psutil.Process()
	parent.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
	for child in parent.children():
		child.nice(psutil.IDLE_PRIORITY_CLASS)
	
	args_list = [(data, E_in, p_dict_list[k], p_dict_bools, data_type) for k in range(no_evals)]
	
	Res = po.map_async(evaluate,args_list)
	result = Res.get()
	po.close()
	po.join()
	
	if verbose: print 'RR calculation complete'

	#Find best fit
	result = np.array(result)
	#print result
	#result = result.astype(np.float64)
	lineMin = np.argmin(result, axis=0)[0] ## pick the fit with the lowest cost value

	best_values = result[lineMin][1] # best parameter dictionary
	if verbose:
		print '\n\n\n'
		print best_values

	p_dict_best = copy.deepcopy(p_dict)
	p_dict_best.update(best_values)
	
	# Finally run the ML fitting one more time, using the best parameters 
	# (so we get the final_result object, which cannot be pickled and therefore isn't supported in multiprocessing)
	best_values, final_result = ML.ML_fit(data, E_in, p_dict_best, p_dict_bools, data_type)
	
	# return best fit parameters, and the lmfit result object
	return best_values, final_result

def test_fit():
	p_dict = {'elem':'Rb','Dline':'D2','T':80.,'lcell':2e-3,'Bfield':600.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':0.,'shift':0.}
	
	# only need to specify parameters that are varied
	p_dict_bools = {'T':True,'Bfield':True,'E_x':True}
	p_dict_bounds = {'T':10,'Bfield':100,'E_x':0.01}
	
	E_in = np.array([0.7,0.7,0])
	E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]
	
	print E_in_angle
	
	x = np.linspace(-10000,10000,100)
	[y] = get_spectra(x,E_in,p_dict,outputs=['S1']) + np.random.randn(len(x))*0.015
		
	data = [x,y.real]
	
	best_params, result = RR_fit(data, E_in_angle, p_dict, p_dict_bools, p_dict_bounds, no_evals = 8, data_type='S1')
	report = result.fit_report()
	fit = result.best_fit
	
	print report
	plt.plot(x,y,'ko')
	plt.plot(x,fit,'r-',lw=2)
	plt.show()
	
if __name__ == '__main__':
	test_fit()
	
