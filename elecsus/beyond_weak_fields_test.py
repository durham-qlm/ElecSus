from elecsus_methods import spectra, fit_data
import datetime
import matplotlib.pyplot as plt
import numpy as np


#################################################################################################################
# Define some necessary variables
#################################################################################################################
x = np.linspace(-4000, 6000, 1000)
E_in = np.sqrt(1/2) * np.array([1,1,0])
p_dict =     {'Elem': 'Rb', 'Dline': 'D1', 'T': 40, 'DoppTemp': -80, 'Constrain': True, 'lcell': 75e-3, 'Bfield':0.,'Btheta':0.,
		'Bphi':0.,'GammaBuf':0.,'shift':0.}
p_dict_bwf = {**p_dict, 'laserPower': 1e-15, 'laserWaist': 5e-3}

#################################################################################################################
# Show that BWF correctly uses same p_dict interface as ElecSus
#################################################################################################################
y = spectra.get_spectra(x, E_in=E_in, p_dict=p_dict, outputs=['S0'])[0]
y_bwf = spectra.get_spectra(x, E_in=E_in, p_dict=p_dict_bwf, outputs=['S0'])[0]

plt.figure(tight_layout=True)
plt.plot(x, y, label='ElecSus')
plt.plot(x, y_bwf, '--', label='BWF')
plt.plot(x, 1 + (y - y_bwf) / y, label='relative error')
plt.legend()

#################################################################################################################
# Now fit to show that ElecSus has problems with light > saturation intensity
#################################################################################################################
p_dict_bwf_fit = {**p_dict_bwf}
p_dict_bwf_fit['laserPower'] = 4.5e-2
p_dict_bwf_fit['shift'] = 500
y_noised = spectra.get_spectra(x, E_in=E_in, p_dict=p_dict_bwf_fit, outputs=['S0'])[0] + np.random.randn(len(x))*0.002

# need to express E_in as Ex, Ey and phase difference for fitting
E_in_angle = [E_in[0].real,[abs(E_in[1]),np.angle(E_in[1])]]

# Normal ElecSus
p_dict_bools = {}
p_dict_bools['T'] = True
p_dict_bools['shift'] = True
p_dict_guess = p_dict
p_dict_guess['T'] += np.random.randn()*1

t0 = datetime.datetime.now()
best_params, _, result = fit_data([x, y_noised], p_dict_guess, p_dict_bools, E_in=E_in_angle, data_type='S0',fit_algorithm='ML')
print(datetime.datetime.now() - t0)

# With BWF extension
p_dict_bools['laserPower'] = True
p_dict_guess = p_dict_bwf_fit
p_dict_guess['T'] += np.random.randn()*1
p_dict_guess['laserPower'] *= np.random.randn()*1

t0 = datetime.datetime.now()
best_params_bwf, _, result_bwf = fit_data([x, y_noised], p_dict_guess, p_dict_bools, E_in=E_in_angle, data_type='S0',fit_algorithm='ML')
print(datetime.datetime.now() - t0)

plt.figure(tight_layout=True)
plt.plot(x, y_noised, label='noised generated spectrum')
plt.plot(x, result.best_fit, '--', label='best fit ElecSus')
plt.plot(x, result_bwf.best_fit, '--', label='best fit BWF')
plt.xlabel('Detuning [MHz]')
plt.ylabel('Transmission [a.u.]')
plt.legend()
plt.show()