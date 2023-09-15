# Copyright 2014 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams and I. G. Hughes.

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

Simulated annealing fitting routine. 

Complete rebuild of the original SA fitting module written by Mark A. Zentile, now using lmfit

Last updated 2018-07-12 MAZ
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)

import numpy as np
import matplotlib.pyplot as plt
import warnings
import sys
import copy

import pickle as pickle

import psutil
from multiprocessing import Pool

from . import MLFittingRoutine as ML
import lmfit as lm

from .spectra import get_spectra

p_dict_bounds_default = {'lcell':1e-3,'Bfield':100., 'T':20., 
							'GammaBuf':20., 'shift':100.,
							# Polarisation of light
							'theta0':10., 'E_x':0.05, 'E_y':0.05, 'E_phase':0.01, 
							# B-field angle w.r.t. light k-vector
							'Btheta':10*3.14/180, 'Bphi':10*3.14/180,
							'DoppTemp':20.,
							'rb85frac':1, 'K40frac':1, 'K41frac':1,
							}

def chisq(yo,ye,err=None):
	""" 
	Evaluate the chi-squared value of a given data/theory combination
	
	Inputs:
		yo : observed data
		ye : expected (theory) data
		err : Optional, error bars - array of length(x)
	
	Returns:
		float of chi-squared value
	"""
	if err is None:
		err = np.ones_like(yo)
	return (((yo-ye)/err)**2).sum()

def evaluate(args):
	""" Evaluate chi-squared value for a given set of parameters """
	
	warnings.simplefilter("ignore")
	data = args[0]
	
	p_dict = args[2]
	E_in = np.array([p_dict['E_x'],p_dict['E_y']*np.exp(1.j*p_dict['E_phase']),0.])
	p_dict_bools = args[3]
	data_type = args[4]
	theory_vals = get_spectra(data[0],E_in,p_dict,outputs=[data_type])[0].real
	chisq_val = chisq(data[1],theory_vals)

	return chisq_val, p_dict

def SA_fit(data,E_in,p_dict,p_dict_bools,p_dict_bounds=None,no_evals=None,data_type='S0',verbose=False):
	"""
	Simulated annealing fitting method.
	
	Before simulated annealing starts, the parameter space is randomly sampled to find good starting conditions for the SA fit.

	data:					an Nx2 iterable for the x and y data to be fitted
	E_in:					the initial electric field input. See docstring for the spectra.py module for details.

	no_evals:			The number of randomly-selected start points for downhill fitting. Defaults to 2**(3+2*nFitParams) where nFitParams is
									the number of varying fit parameters
	
	p_dict:					dictionary containing all the calculation (initial) parameters
	p_dict_bools:		dictionary with the same keys as p_dict, with Boolean values representing each parameter that is to be varied in the fitting
	p_dict_bounds:	dictionary with the same keys as p_dict, with values that represent the deviation each parameter can take in the initial parameter search
									NOTE: this works slightly differently to p_dict_bounds in the ML fitting methods. In RR and SA fitting, the bounds
									select the range in parameter space that is randomly explored to find good starting parameters for the 
									SA routine, rather than being strict bounds on the fit parameters.
	
	data_type:			Data type to fit experimental data to. Can be one of:
									'S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', ...
	verbose:				Boolean - more print statements provided as the program progresses
	
	"""
	
	if p_dict_bounds is None:
		p_dict_bounds = p_dict_bounds_default
	
	print('Starting Simulated Annealing Fitting Routine')
	x = np.array(data[0])
	y = np.array(data[1])

	p_dict['E_x'] = E_in[0]
	p_dict['E_y'] = E_in[1][0]
	p_dict['E_phase'] = E_in[1][1]
	E_in_vector = np.array([p_dict['E_x'],p_dict['E_y']*np.exp(1.j*p_dict['E_phase']),0.])
	
	# count number of fit parameters
	nFitParams = 0
	for key in p_dict_bools:
		if p_dict_bools[key]: nFitParams += 1
	
	# default number of iterations based on number of fit parameters
	if no_evals == None:
		no_evals = 2**(8+2*nFitParams)
		
	# Create random array of starting parameters based on parameter ranges given in p_dict_bounds dictionary
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
		print('List of initial parameter dictionaries:')
		for pd in p_dict_list:
			print(pd)
		#print p_dict_list
		print('\n\n')
	
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
		
	result = np.array(result)
	lineMin = np.argmin(result[:,0]) # position of lowest chi-squared value from initial guesses
	SA_params = result[lineMin][1] # parameter dictionary associated with lowest chi-squared
	best_chi_sq = result[lineMin][0] # lowest chi-squared value
	current_chi_sq = best_chi_sq
	current_params = copy.deepcopy(SA_params)
	
	print('\n\n Initial parameter space evaluations completed...')
	
	# Measure of how much the chi-squared values change over the parameter space
	spread = np.std(result[:,0])
	
	# Initial 'temperature' for the annealing
	T = 1500.
	# Cold temperature extreme where the algorithm will stop searching
	minimum_T = 50.
	
	# sets the effective range of a given temperature
	range_scaling = 3
	k = abs(spread)/T / range_scaling

	# Cooling rate - how much the temperature decreases on each run.
	cooling_rate = 7e-6

	# Number of times a jump uphill is rejected before switching to ML fit.
	uphill_escape_chance = 0
	uphill_escape_threshold = 300
	
	# Number of iterations to investigate a plateau in parameter space (i.e. chi-squared does not change)
	plateau_escape_chance = 0
	plateau_escape_threshold = 300
	
	trial_params = copy.deepcopy(current_params)
	best_params = copy.deepcopy(current_params)
	best_chi_sq = current_chi_sq
	#print trial_params
	
	## testing - log algorithm path
	chi_sq_log = [current_chi_sq]
	T_log = [trial_params['T']]
	B_log = [trial_params['Bfield']]

	hot_iterations = 250
	iteration = 0
	
	while T > minimum_T: 
		# Generate a new set of parameters that are close to the last evaluated parameters
		for key in p_dict_bools:
			if p_dict_bools[key]==True: # (extra check)
				## parameters varied over 10% of range specified before...	
				trial_params[key] = copy.deepcopy(current_params[key]) + np.random.normal(0,p_dict_bounds[key]) # * p_dict_bounds[key]	#* 0.1
		
		print('Before calculation:')
		print(('Trial: ', trial_params['T'], trial_params['Bfield']))
		print(('Current: ', current_params['T'], current_params['Bfield']))

		trial_theory = get_spectra(x,E_in_vector,trial_params,outputs=[data_type])[0].real
		trial_chi_sq = chisq(data[1],trial_theory)
		
		# log best chi_squared values and parameters
		if trial_chi_sq < best_chi_sq:
			best_chi_sq = trial_chi_sq
			best_params = copy.deepcopy(trial_params)
		
		# Calculate energy difference
		delta_chi_sq = trial_chi_sq - current_chi_sq

		# and convert to probability of acceptance
		if delta_chi_sq < 1e-5:
			# if there's no change in chi-squared (i.e. parameter space is locally flat), accept with small (5%) probability
			print('WARNING - no change in chi-squared - probably due to plateau in parameter space!')
			prob = 0.05
			plateau_escape_chance += 1
		else:
			prob = np.exp(-delta_chi_sq/(k*T))
		
		if verbose:
			print(('\n\tBest chi-squared so far:', best_chi_sq))
			print(('\tBest Parameters so far (T, B): ', best_params['T'], best_params['Bfield']))

			print(('\tCurrent Parameters (T, B): ', current_params['T'], current_params['Bfield']))
			print(('\tTrial Parameters (T, B): ', trial_params['T'], trial_params['Bfield']))
			print(('\n\tCurrent chi-squared:', current_chi_sq))
			print(('\tTrial chi-squared:', trial_chi_sq))
			print(('\tChange in chi-squared from previous:', delta_chi_sq))
			print(('\tTemperature parameter: ', T))
			print(('\tProbability that new parameters will be accepted (>1 == 1):', prob))
		
		if (delta_chi_sq < 0) or (prob > np.random.random()):
			# accept downhill movement, or uphill movement with probability prob - update chi_squared and parameters
			current_chi_sq = trial_chi_sq
			current_params = copy.deepcopy(trial_params)
			
			chi_sq_log.append(trial_chi_sq) ## keep log of chi-squared values (on succesful iterations only)
			T_log.append(trial_params['T'])
			B_log.append(trial_params['Bfield'])

			print('\t...Values accepted. Current parameters updated.')
			print('\n')
			
			# reset escape chance
			uphill_escape_chance = 0
		else:
			print(('\t...Values rejected. Escape threshold:', uphill_escape_chance, ' / ', uphill_escape_threshold))
			print('\n')
			uphill_escape_chance += 1
			
		
		# Cool system
		# Hold T constant for first N iterations
		if iteration > hot_iterations:
			T = T/(1 + (cooling_rate*T)) #Lundy's Method (http://link.springer.com/article/10.1007/BF01582166)
		
		iteration += 1
		
		# Exit annealing loop if conditions are correct
		if (uphill_escape_chance > uphill_escape_threshold):
			print(('Simulated annealing completed ( No improvement found after {:d} iterations)'.format(uphill_escape_threshold)))
			print('Switching to downhill fit using best found parameters...\n')
			break

		
		if (T < minimum_T): #No jumps up hill for a while, or minimum temperature reached
			print('Simulated annealing completed ( Temperature reached minimum threshold )')
			print('Switching to downhill fit using best found parameters...\n')
			break
		
		if (plateau_escape_chance > plateau_escape_threshold):
			print('!!!\n\tCAUTION :: Annealing has not converged.')
			print('\tAnnealing algorithm found plateau in parameter space')
			print('\tSwitching to downhill fit using best found parameters...\n!!!')
			break

	#### Marquardt-Levenberg fit #####

	try:
		print(('Downhill fit with initial parameters: (T,B) ', best_params['T'], best_params['Bfield']))
		ML_best_params, result = ML.ML_fit(data, E_in, best_params, p_dict_bools, data_type)
		MLchi_sq = result.chisqr
		if MLchi_sq <= best_chi_sq:
			best_params = ML_best_params
			final_result = result
			success = True
			print('Downhill fit converged successfully.\n')
		else:
			print('Downhill fit did not find further improvement. Continuing with simulated annealing result.\n')
			success = False
			final_result = 1			
	except:
		print('Downhill fit failed to converge. Continuing with simulated annealing result.\n')
		success = False
		final_result = 1
		
	return best_params, final_result #, chi_sq_log, T_log, B_log, success
	
def test_fit(calc=False):
	p_dict = {'elem':'Rb','Dline':'D2','T':80.,'lcell':2e-3,'Bfield':600.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':0.,'shift':0.}
	
	# only need to specify parameters that are varied
	p_dict_bools = {'T':True,'Bfield':True}
	p_dict_bounds = {'T':20,'Bfield':200}
	
	E_in = np.array([0.7,0.7,0])
	E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]
	
	print(E_in_angle)
	
	x = np.linspace(-10000,10000,200)
	if calc:
		[y] = get_spectra(x,E_in,p_dict,outputs=['S1']) + np.random.randn(len(x))*0.01
		pickle.dump([x,y],open('pickle_xydata.pkl','wb'))
	else:
		x,y = pickle.load(open('pickle_xydata.pkl','rb'))

	data = [x,y.real]

	# Map chi-squared surface:
	Tvals = np.linspace(60,100,200)
	Bvals = np.linspace(400,800,220)
	
	T2D, B2D = np.meshgrid(Tvals,Bvals)
	
	CSQ_map = np.zeros((len(Tvals),len(Bvals)))
	
	if calc:
		for i,TT in enumerate(Tvals):
			print((i, len(Tvals)))
			for j, BB in enumerate(Bvals):
				p_dict['T'] = TT
				p_dict['Bfield'] = BB
				[ye] = get_spectra(x,E_in,p_dict,outputs=['S1'])
				CSQ_map[i,j] = chisq(y,ye)
		pickle.dump([T2D,B2D,CSQ_map],open('pickle_CSQmap.pkl','wb'))
	else:
		T2D, B2D, CSQ_map = pickle.load(open('pickle_CSQmap.pkl','rb'))
	
	fig_map = plt.figure()
	ax = fig_map.add_subplot(111)
	ax.imshow(CSQ_map.T,origin='lower',aspect='auto',
					extent=(Tvals[0],Tvals[-1],Bvals[0],Bvals[-1]),
					cmap=plt.cm.jet,alpha=0.7)
	ax.contour(T2D, B2D, CSQ_map.T,7,lw=2,color='k')

	#plt.show()
	
	## Do SA fitting
	best_params, result, chi_sq_log, T_log, B_log = SA_fit(data, E_in_angle, p_dict, p_dict_bools,
																					p_dict_bounds, no_evals = 4, data_type='S1')
	report = result.fit_report()
	fit = result.best_fit
	
	fig_data = plt.figure()
	print(report)
	plt.plot(x,y,'ko')
	plt.plot(x,fit,'r-',lw=2)
	
	# Chi-squared log with iteration number
	fig_chisqlog = plt.figure()
	plt.plot(chi_sq_log)
	plt.ylabel('Chi_squared value')
	plt.xlabel('Iteration')
	
	
	#plt.show()
	
	ax.plot(T_log, B_log, 'ko', ms=1)
	
	plt.show()

def test_fit2(calc=False):
	""" Alternate data set to test """
	
	p_dict = {'elem':'Rb','Dline':'D2','T':80.,'lcell':2e-3,'Bfield':250.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':0.,'shift':0.}
	
	# only need to specify parameters that are varied
	p_dict_bools = {'T':True,'Bfield':True}
	p_dict_bounds = {'T':15,'Bfield':100}
	
	E_in = np.array([1.0,0.0,0.0])
	E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]
	
	print(E_in_angle)
	
	x = np.linspace(-7000,8000,200)
	if calc:
		[y] = get_spectra(x,E_in,p_dict,outputs=['S0']) + np.random.randn(len(x))*0.02
		pickle.dump([x,y],open('pickle_xydata.pkl','wb'))
	else:
		x,y = pickle.load(open('pickle_xydata.pkl','rb'))

	data = [x,y.real]

	# Map chi-squared surface:
	Tvals = np.linspace(40,120,300)
	Bvals = np.linspace(0,500,250)
	
	T2D, B2D = np.meshgrid(Tvals,Bvals)
	CSQ_map = np.zeros((len(Tvals),len(Bvals)))
	
	if calc:
		for i,TT in enumerate(Tvals):
			print((i, len(Tvals)))
			for j, BB in enumerate(Bvals):
				p_dict['T'] = TT
				p_dict['Bfield'] = BB
				[ye] = get_spectra(x,E_in,p_dict,outputs=['S0'])
				CSQ_map[i,j] = chisq(y,ye)
		pickle.dump([T2D,B2D,CSQ_map],open('pickle_CSQmap.pkl','wb'))
	else:
		T2D, B2D, CSQ_map = pickle.load(open('pickle_CSQmap.pkl','rb'))
	
	fig_map = plt.figure()
	ax = plt.subplot2grid((1,14),(0,0),colspan=13)
	ax_cb = plt.subplot2grid((1,14),(0,13),colspan=1)
	im = ax.imshow(CSQ_map.T,origin='lower',aspect='auto',
					extent=(Tvals[0],Tvals[-1],Bvals[0],Bvals[-1]),
					cmap=plt.cm.jet,alpha=0.7)
	ax.contour(T2D, B2D, CSQ_map.T,7,lw=2,color='k')
	
	cb = fig_map.colorbar(im,cax=ax_cb)
	
	#plt.show()
	
	## Do SA fitting
	best_params, result, chi_sq_log, T_log, B_log = SA_fit(data, E_in_angle, p_dict, p_dict_bools,
																					p_dict_bounds, no_evals = 2, data_type='S0')
	report = result.fit_report()
	fit = result.best_fit
	
	fig_data = plt.figure()
	print(report)
	plt.plot(x,y,'ko')
	plt.plot(x,fit,'r-',lw=2)
	
	# Chi-squared log with iteration number
	fig_chisqlog = plt.figure()
	plt.plot(chi_sq_log)
	plt.ylabel('Chi_squared value')
	plt.xlabel('Iteration')
	
	
	#plt.show()
	
	ax.plot(T_log, B_log, 'ko', ms=1)
	ax.plot(T_log[0], B_log[0], 'rs', ms=6)
	ax.plot(T_log[-1], B_log[-1], 'bs', ms=6)
	
	plt.show()

def test_fit3(calc=False):
	""" Alternate data set to test """
	
	p_dict = {'elem':'Rb','Dline':'D2','T':155.,'lcell':5e-3,'Bfield':1100.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':0.,'shift':0.}
	
	# only need to specify parameters that are varied
	p_dict_bools = {'T':True,'Bfield':True}
	p_dict_bounds = {'T':25,'Bfield':300}
	
	E_in = np.array([1.0,0.0,0.0])
	E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]
	
	print(E_in_angle)
	
	x = np.linspace(-13000,-6000,150)
	if calc:
		[y] = get_spectra(x,E_in,p_dict,outputs=['S1']) + np.random.randn(len(x))*0.02
		pickle.dump([x,y],open('pickle_xydata.pkl','wb'))
	else:
		x,y = pickle.load(open('pickle_xydata.pkl','rb'))

	data = [x,y.real]

	# Map chi-squared surface:
	#Tvals = np.linspace(p_dict['T']-p_dict_bounds['T'],p_dict['T']+p_dict_bounds['T'],150)
	#Bvals = np.linspace(p_dict['Bfield']-p_dict_bounds['Bfield'],p_dict['Bfield']+p_dict_bounds['Bfield'],125)
	
	Tvals = np.linspace(80,250,450)
	Bvals = np.linspace(0,2500,300)
	
	T2D, B2D = np.meshgrid(Tvals,Bvals)
	CSQ_map = np.zeros((len(Tvals),len(Bvals)))
	
	if calc:
		for i,TT in enumerate(Tvals):
			print((i, len(Tvals)))
			for j, BB in enumerate(Bvals):
				p_dict['T'] = TT
				p_dict['Bfield'] = BB
				[ye] = get_spectra(x,E_in,p_dict,outputs=['S1'])
				CSQ_map[i,j] = chisq(y,ye)
		pickle.dump([T2D,B2D,CSQ_map],open('pickle_CSQmap.pkl','wb'))
	else:
		T2D, B2D, CSQ_map = pickle.load(open('pickle_CSQmap.pkl','rb'))
	
	print(T2D)
	print(B2D)
	fig_map = plt.figure()
	ax = plt.subplot2grid((1,14),(0,0),colspan=13)
	ax_cb = plt.subplot2grid((1,14),(0,13),colspan=1)
	im = ax.imshow(CSQ_map.T,origin='lower',aspect='auto',
					extent=(T2D[0][0],T2D[-1][-1],B2D[0][0],B2D[-1][-1]),
					cmap=plt.cm.jet,alpha=0.7)
	ax.contour(T2D, B2D, CSQ_map.T,9,lw=2,color='k')
	
	cb = fig_map.colorbar(im,cax=ax_cb)
	
	#plt.show()
	
	## Do SA fitting
	best_params, result, chi_sq_log, T_log, B_log, success = SA_fit(data, E_in_angle, p_dict, p_dict_bools,
																					p_dict_bounds, no_evals = 32, data_type='S1')
	if success:
		report = result.fit_report()
		fit = result.best_fit
		print(report)
	
	fig_data = plt.figure()
	plt.plot(x,y,'ko')
	if success:
		plt.plot(x,fit,'r-',lw=2)
	
	# Chi-squared log with iteration number
	fig_chisqlog = plt.figure()
	plt.plot(chi_sq_log)
	plt.ylabel('Chi_squared value')
	plt.xlabel('Iteration')
	
	
	#plt.show()
	
	ax.plot(T_log, B_log, 'ko', ms=1)
	ax.plot(T_log, B_log, 'k--', lw=1.5)
	ax.plot(T_log[0], B_log[0], 'rs', ms=6)
	ax.plot(T_log[-1], B_log[-1], 'bs', ms=6)
	
	plt.show()

def test_fit4(calc=False):
	""" Alternate data set to test """
	
	#actual params 110 / 7000
	
	p_dict = {'elem':'Rb','Dline':'D2','T':100.,'lcell':5e-3,'Bfield':6100.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':0.,'shift':0.,'rb85frac':1.0}
	
	# only need to specify parameters that are varied
	p_dict_bools = {'T':True,'Bfield':True}
	p_dict_bounds = {'T':5,'Bfield':400}
	
	E_in = np.array([1.0,0.0,0.0])
	E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]
	
	print(E_in_angle)
	
	x = np.linspace(-21000,-13000,150)
	if calc:
		[y] = get_spectra(x,E_in,p_dict,outputs=['S0']) + np.random.randn(len(x))*0.02
		pickle.dump([x,y],open('pickle_xydata.pkl','wb'))
	else:
		x,y = pickle.load(open('pickle_xydata.pkl','rb'))

	data = [x,y.real]

	# Map chi-squared surface:
	#Tvals = np.linspace(p_dict['T']-p_dict_bounds['T'],p_dict['T']+p_dict_bounds['T'],150)
	#Bvals = np.linspace(p_dict['Bfield']-p_dict_bounds['Bfield'],p_dict['Bfield']+p_dict_bounds['Bfield'],125)
	
	Tvals = np.linspace(70,150,350)
	Bvals = np.linspace(5000,10000,300)
	
	T2D, B2D = np.meshgrid(Tvals,Bvals)
	CSQ_map = np.zeros((len(Tvals),len(Bvals)))
	
	if calc:
		for i,TT in enumerate(Tvals):
			print((i, len(Tvals)))
			for j, BB in enumerate(Bvals):
				p_dict['T'] = TT
				p_dict['Bfield'] = BB
				[ye] = get_spectra(x,E_in,p_dict,outputs=['S0'])
				CSQ_map[i,j] = chisq(y,ye)
		pickle.dump([T2D,B2D,CSQ_map],open('pickle_CSQmap.pkl','wb'))
	else:
		T2D, B2D, CSQ_map = pickle.load(open('pickle_CSQmap.pkl','rb'))
	
	print(T2D)
	print(B2D)
	fig_map = plt.figure()
	ax = plt.subplot2grid((1,14),(0,0),colspan=13)
	ax_cb = plt.subplot2grid((1,14),(0,13),colspan=1)
	im = ax.imshow(CSQ_map.T,origin='lower',aspect='auto',
					extent=(T2D[0][0],T2D[-1][-1],B2D[0][0],B2D[-1][-1]),
					cmap=plt.cm.jet,alpha=0.7)
	ax.contour(T2D, B2D, CSQ_map.T,9,lw=2,color='k')
	
	cb = fig_map.colorbar(im,cax=ax_cb)
	
	#plt.show()
	
	## Do SA fitting
	best_params, result, chi_sq_log, T_log, B_log, success = SA_fit(data, E_in_angle, p_dict, p_dict_bools,
																					p_dict_bounds, no_evals = 32, data_type='S0')
		
	fig_data = plt.figure()
	plt.plot(x,y,'ko')
	if success:
		report = result.fit_report()
		fit = result.best_fit
		print(report)
	else:
		print(best_params)
		fit = get_spectra(x,E_in,best_params,outputs=['S0'])[0].real

	plt.plot(x,fit,'r-',lw=2)
	
	# Chi-squared log with iteration number
	fig_chisqlog = plt.figure()
	plt.plot(chi_sq_log)
	plt.ylabel('Chi_squared value')
	plt.xlabel('Iteration')
	
	
	#plt.show()
	
	ax.plot(T_log, B_log, 'ko', ms=1)
	ax.plot(T_log, B_log, 'k--', lw=1.5)
	ax.plot(T_log[0], B_log[0], 'rs', ms=6)
	ax.plot(T_log[-1], B_log[-1], 'bs', ms=6)
	
	plt.show()
	
if __name__ == '__main__':
	test_fit4(False)
	
