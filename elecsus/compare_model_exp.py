from elecsus_methods import spectra, fit_data
from libs.beyond_weak_fields import atomic_system as ats
import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
import seaborn as sns

cmap1 = sns.color_palette('mako_r', as_cmap=True)
cmap2 = sns.color_palette('flare', as_cmap=True)
groundState = ats.state(5, 0, 1/2)     # 5S1/2
excitedState_D2 = ats.state(5, 1, 3/2)  # 5P1/2


def fit(x, x0, a, b):
	return a * voigt_profile(x+x0, b, 6.0) + 1


def get_exp_peaks(set):
	files = glob.glob(f'{set}/*.npz')
	powers = np.empty((len(files)))
	for i, f in enumerate(files):
		powers[i] = float(f.split('/')[-1].split('W')[0][:-1])
		if f.split('/')[-1].split('W')[0][-1] == 'u':
			powers[i] *= 1e-6
		else:
			powers[i] *= 1e-3
	powers *= 1e6 # uW
	
	D = float(f.split('/')[2].split('mm')[0]) / 10 # cm
	TF1 = np.empty_like(powers)
	TF2 = np.empty_like(powers)

	for i, file in enumerate(files):
		with np.load(file) as f:
			x = f['x']
			y = f['y']
		peak = np.argmin(y)
		x = x - x[peak] - 2500
		poptF1, _, = curve_fit(fit, x[x>0], y[x>0], p0=(-4000, -350, 278))
		poptF2, _, = curve_fit(fit, x[x<0], y[x<0], p0=(+2500, -350, 278))
		TF1[i] = np.min(fit(x, *poptF1))
		TF2[i] = np.min(fit(x, *poptF2))

	alphaF1 = -np.log(TF1) / 2e-3
	alphaF2 = -np.log(TF2) / 2e-3
	alphaF1 /= alphaF1.max()
	alphaF2 /= alphaF2.max()

	idx = np.argsort(powers)
	Isat = 1.669e3 # uW / cm2
	P2I = lambda P: 2 * P / np.pi / (D/2)**2
	I = P2I(powers) / Isat
	return I[idx], alphaF1[idx], alphaF2[idx]


def get_sim_peaks(D):
	xF1 = np.linspace(4000, 4400, 100) * 1e6
	xF2 = np.linspace(-2600, -2200, 100) * 1e6
	Isat_factor = np.geomspace(0.01, 1000, 40, endpoint=True)
	F1 = np.empty(Isat_factor.size)
	F2 = np.empty(Isat_factor.size)
	Isat = 16.69 # W/m2
	Psat = np.pi * Isat * D**2 / 8
	p_dict = {'T': 54., 'laserWaist': D, 'rb85frac': 0.7}
	Rb87_D2 = ats.atomicSystem('Rb87', [groundState, excitedState_D2], p_dict=p_dict)
	for j, factor in enumerate(Isat_factor):
		print(f'{j/Isat_factor.size * 100:.0f} %')
		a1 = -Rb87_D2.optical_depth([ats.beam(w=xF1, P=factor * Psat, D=D)], doppler=True)
		a2 = -Rb87_D2.optical_depth([ats.beam(w=xF2, P=factor * Psat, D=D)], doppler=True)
		F1[j] = np.max(a1)
		F2[j] = np.max(a2)
	F1 /= F1.max()
	F2 /= F2.max()
	return Isat_factor, F1, F2

	
# import sys
# p_dict = {'T': 55., 'laserWaist': 0.5e-3, 'rb85frac': 0.7}
# Rb87_D2 = ats.atomicSystem('Rb87', [groundState, excitedState_D2], p_dict=p_dict)
# Rb87_D2.optical_depth([ats.beam(w=0, P=5e-3, D=0.5e-3)])
# sys.exit()

basefolder = 'durham_analysis/peak_measurements'
dataset05 = f'{basefolder}/0.5mm_iris_8.5V'
dataset26 = f'{basefolder}/2.6mm_iris_8.5V'

I05e, a1_05e, a2_05e = get_exp_peaks(dataset05)
I26e, a1_26e, a2_26e = get_exp_peaks(dataset26)
I05m, a1_05m, a2_05m = get_sim_peaks(0.5e-3)
I26m, a1_26m, a2_26m = get_sim_peaks(2.6e-3)

plt.figure(tight_layout=True)
plt.plot(I05m, a1_05m, c=cmap1(0.2))
plt.plot(I05m, a2_05m, c=cmap2(0.2))
plt.plot(I26m, a1_26m, c=cmap1(0.5))
plt.plot(I26m, a2_26m, c=cmap2(0.5))
plt.plot(I05e, a1_05e, 'x', c=cmap1(0.2))
plt.plot(I05e, a2_05e, 'x', c=cmap2(0.2))
plt.plot(I26e, a1_26e, 'x', c=cmap1(0.5))
plt.plot(I26e, a2_26e, 'x', c=cmap2(0.5))
plt.xlabel('I/Isat')
plt.ylabel(r'α/α0')
plt.xscale('log')
plt.savefig(f'{basefolder}/model_vs_measurement_peaks.png', dpi=200)
plt.show()
