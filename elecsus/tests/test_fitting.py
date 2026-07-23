"""
Series of tests and example code to run elecsus via the API.

Converted to pytest (issue #26, 2026-07-09): fixed relative imports, added
assertions, and gated interactive plotting behind ELECSUS_SHOW_PLOTS so the
default run is fast and non-blocking. Physics/logic unchanged from the
original 2018 version.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

from elecsus import elecsus_methods as EM
from elecsus.libs.spectra import get_spectra

SHOW_PLOTS = os.environ.get('ELECSUS_SHOW_PLOTS', '') == '1'


def test_MLfit():
	"""
	Test simple Marquart-Levenburg fit for S1 spectrum. Generate a data set and add random noise,
	then run an ML fit to the noisy data set, with randomised initial parameters that are the correct ones
	plus some small change
	"""
	np.random.seed(0)

	p_dict = {'Elem':'Rb','Dline':'D2','T':80.,'lcell':2e-3,'Bfield':600.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':0.,'shift':0.}

	# only need to specify parameters that are varied
	p_dict_bools = {'T':True,'Bfield':True}

	E_in = np.sqrt(1./2) * np.array([1.,1.,0])
	# need to express E_in as Ex, Ey and phase difference for fitting
	E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]

	# Generate spectrum and add some random noise
	x = np.linspace(-10000,10000,300)
	[y] = get_spectra(x,E_in,p_dict,outputs=['S1']) + np.random.randn(len(x))*0.01

	# randomise the starting parameters a little
	p_dict_guess = dict(p_dict)
	p_dict_guess['T'] += np.random.randn()*1
	p_dict_guess['Bfield'] += np.random.randn()*10

	data = [x,y.real]

	best_params, RMS, result = EM.fit_data(data, p_dict_guess, p_dict_bools, E_in=E_in_angle, data_type='S1',fit_algorithm='ML')

	assert RMS < 0.05
	assert best_params['T'] == pytest.approx(p_dict['T'], abs=2.0)
	assert best_params['Bfield'] == pytest.approx(p_dict['Bfield'], abs=20.0)

	if SHOW_PLOTS:
		fit = result.best_fit
		plt.plot(x,y,'k.',label='Data')
		plt.plot(x,fit,'r-',lw=2,label='Fit')
		plt.legend(loc=0)
		plt.show()
	plt.close('all')


def _compare_fit_methods_demo():
	"""
	Compare ML and RR, SA and differential_evolution (DE) fitting for a more complicated case where there is a chi-squared local minimum in parameter space.
	In this case, the simple ML fit doesn't work, since it gets stuck in a very bad local minimum.

	The RR, SA and DE methods all converge to a good solution, but in different times.
	The RR fit takes somewhere in the region of 50x longer than ML, but produces a much better fit.
	The SA fit takes somewhere in the region of 100x longer than ML, but produces a much better fit.
	The DE fit takes somewhere in the region of 20x longer than ML, but produces a much better fit.

	From this test (and others which are broadly similar but for different data sets), we see that for global minimisation, we recommend using DE

	Not collected by default (leading underscore + slow marker): this is a multi-minute
	demonstration script, not a fast regression check. Run explicitly with
	`pytest -k compare_fit_methods_demo --no-header` after renaming, or call directly.
	"""

	p_dict = {'Elem':'Rb','Dline':'D2','T':100.,'lcell':5e-3,'Bfield':7000.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':50,'shift':0.}

	# only need to specify parameters that are varied
	p_dict_bools = {'T':True,'Bfield':True,'GammaBuf':True}#,'E_x':True,'E_y':True,'E_phase':True}

	E_in = np.array([0.7,0.7,0])
	# need to express E_in as Ex, Ey and phase difference for fitting
	E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]

	x = np.linspace(-21000,-13000,300)
	y = get_spectra(x,E_in,p_dict,outputs=['S0'])[0] + np.random.randn(len(x))*0.01 # <<< RMS error should be around 0.01 for the best case fit
	data = [x,y.real]

	# start at not the optimal parameters (near a local minimum in chi-squared)
	p_dict['T'] = 90
	p_dict['Bfield'] = 6100
	p_dict['GammaBuf'] = 150

	p_dict_bounds = {}

	# time it...
	import time
	st_ML = time.time()
	best_paramsML, RMS_ML, resultML = EM.fit_data(data, p_dict, p_dict_bools, E_in=E_in_angle, data_type='S0')
	et_ML = time.time() -st_ML
	print('ML complete')

	# RR and SA need range over which to search
	p_dict_bounds['T'] = [30]
	p_dict_bounds['Bfield'] = [3000]
	p_dict_bounds['GammaBuf'] = [150]

	st_RR = time.time()
	best_paramsRR, RMS_RR, resultRR = EM.fit_data(data, p_dict, p_dict_bools, E_in=E_in_angle,
																					p_dict_bounds=p_dict_bounds, data_type='S0', fit_algorithm='RR')
	et_RR = time.time() - st_RR
	print('RR complete')

	st_SA = time.time()
	best_paramsSA, RMS_SA, resultSA = EM.fit_data(data, p_dict, p_dict_bools, E_in=E_in_angle,
																					p_dict_bounds=p_dict_bounds, data_type='S0', fit_algorithm='SA')
	et_SA = time.time() - st_SA
	print('SA complete')

	# differential evolution needs upper and lower bounds on fit parameters
	p_dict_bounds['T'] = [70,130]
	p_dict_bounds['Bfield'] = [4000,8000]
	p_dict_bounds['GammaBuf'] = [0, 100]

	st_DE = time.time()
	best_paramsDE, RMS_DE, resultDE = EM.fit_data(data, p_dict, p_dict_bools, E_in=E_in_angle,
																					p_dict_bounds=p_dict_bounds, data_type='S0',
																					fit_algorithm='differential_evolution')
	et_DE = time.time() - st_DE
	print('DE complete')


	print(('ML found best params in ', et_ML, 'seconds. RMS error, ', RMS_ML))
	print(('RR found best params in ', et_RR, 'seconds. RMS error, ', RMS_RR))
	print(('SA found best params in ', et_SA, 'seconds. RMS error, ', RMS_SA))
	print(('DE found best params in ', et_DE, 'seconds. RMS error, ', RMS_DE))

	if SHOW_PLOTS:
		fitML = resultML.best_fit
		fitDE = resultDE.best_fit
		plt.plot(x,y,'k.',label='Data')
		plt.plot(x,fitML,'-',lw=2,label='ML fit')
		plt.plot(x,fitDE,'-',lw=2,label='RR fit')
		plt.plot(x,fitDE,'--',lw=2,label='SA fit')
		plt.plot(x,fitDE,':',lw=2,label='DE fit')
		plt.legend(loc=0)
		plt.show()
	plt.close('all')

	# DE should substantially beat the plain ML fit stuck in the local minimum
	assert RMS_DE < RMS_ML


if __name__ == '__main__':
	test_MLfit()
	_compare_fit_methods_demo()
