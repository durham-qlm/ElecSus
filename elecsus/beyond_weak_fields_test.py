from elecsus_methods import spectra
import matplotlib.pyplot as plt
import numpy as np


detuning = np.linspace(-4000, 6000, 1000)
E_in = [1,0,0]
p_dict =     {'Elem': 'Rb', 'Dline': 'D2', 'T': 40, 'DoppTemp': -80, 'Constrain': True, 'lcell': 75e-3, 'shift':0}
p_dict_bwf = {**p_dict, 'laserPower': 1e-15, 'laserWaist': 5e-3, 'bwf_precision': 'low'}
transmission = spectra.get_spectra(detuning, E_in=E_in, p_dict=p_dict, outputs=['S0'])[0]
transmission_bwf = spectra.get_spectra(detuning, E_in=E_in, p_dict=p_dict_bwf, outputs=['S0'])[0]

plt.figure(tight_layout=True)
plt.plot(detuning, transmission, label='ElecSus')
plt.plot(detuning, transmission_bwf, '--', label='BWF')
plt.plot(detuning, 1 + (transmission - transmission_bwf)/transmission_bwf)
plt.legend()
plt.show()
