# Copyright 2017 J. Keaveney

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
A series of examples, calculating various spectra for different parameter regimes.

Converted to pytest (issue #26, 2026-07-09): gated interactive plotting behind
ELECSUS_SHOW_PLOTS and added physical sanity assertions. Physics/logic
otherwise unchanged from the original 2018 version.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from elecsus.elecsus_methods import calculate as get_spectra
from elecsus.libs.durhamcolours import *

SHOW_PLOTS = os.environ.get('ELECSUS_SHOW_PLOTS', '') == '1'


def test_Siddons():
	"""
	The Classic:
	Reproduce the theory plots for the Siddons 2008 paper (Fig 9)
	DOI: 10.1088/0953-4075/41/15/155004
	"""
	np.random.seed(0)

	d = np.arange(-4000,5000,15) # MHz
	#Voigt
	p_dict = {'Bfield':0,'rb85frac':72.17,'Btheta':0,'lcell':75e-3,'T':16.5,'Dline':'D2','Elem':'Rb'}

	[S0_16] = get_spectra(d,[1,0,0],p_dict,outputs=['S0'])
	p_dict['T'] = 25
	[S0_25] = get_spectra(d,[1,0,0],p_dict,outputs=['S0'])
	p_dict['T'] = 36.6
	[S0_36] = get_spectra(d,[1,0,0],p_dict,outputs=['S0'])

	# transmission must stay within physical bounds, and heating should deepen
	# the absorption feature (more atoms -> lower minimum transmission)
	for S0 in (S0_16, S0_25, S0_36):
		assert np.all(np.isfinite(S0))
		assert S0.real.min() >= -1e-6
		assert S0.real.max() <= 1.0 + 1e-6
	assert S0_36.real.min() < S0_25.real.min() < S0_16.real.min()

	S0_16_noisy = S0_16 + np.random.randn(len(d))*0.0025
	S0_25_noisy = S0_25 + np.random.randn(len(d))*0.0025
	S0_36_noisy = S0_36 + np.random.randn(len(d))*0.0025

	if SHOW_PLOTS:
		fig = plt.figure("Faraday comparison")
		yy = 5
		xx = 1
		ax1 = plt.subplot2grid((yy,xx), (0,0), rowspan=yy-1)
		axR = plt.subplot2grid((yy,xx), (yy-1,0), sharex=ax1)
		plt.setp(ax1.get_xticklabels(), visible=False)
		ax1.plot(d/1e3, S0_16_noisy, 'r-')
		ax1.plot(d/1e3, S0_16, 'k--')
		ax1.plot(d/1e3, S0_25_noisy, 'r-')
		ax1.plot(d/1e3, S0_25, 'k--')
		ax1.plot(d/1e3, S0_36_noisy, 'r-')
		ax1.plot(d/1e3, S0_36, 'k--')
		axR.plot(d/1e3, 100*(S0_16_noisy - S0_16), 'k-')
		axR.set_xlabel('Detuning (GHz)')
		ax1.set_ylabel('Transmission')
		axR.set_ylabel(r'Difference ($\times 100$)')
		ax1.set_xlim(-4,5)
		ax1.set_ylim(0,1.03)
		plt.show()
	plt.close('all')


def test_Weller():
	"""
	Reproduce Weller's data from the 2011 paper (Fig 5)
	DOI: 10.1088/0953-4075/44/19/195006

	Shows self-broadening in action...
	"""
	np.random.seed(0)

	d = np.arange(-8000,-3000,15) # MHz
	#Voigt
	p_dict = {'Bfield':0,'rb85frac':72.17,'Btheta':0,'lcell':75e-3,'T':80,'Dline':'D2','Elem':'Rb'}

	Ts = [70,80,90,100,110,120]

	if SHOW_PLOTS:
		fig = plt.figure("Faraday comparison")
		ax1 = fig.add_subplot(111)

	for T in Ts:
		p_dict['T'] = T
		[S0] = get_spectra(d,[1,0,0],p_dict,outputs=['S0'])
		assert np.all(np.isfinite(S0))
		assert -1e-6 <= S0.real.min() and S0.real.max() <= 1.0 + 1e-6

		S0_noisy = S0 + np.random.randn(len(d))*0.0025

		if SHOW_PLOTS:
			ax1.plot(d/1e3, S0_noisy, 'k-')
			ax1.plot(d/1e3, S0, 'r--')

	if SHOW_PLOTS:
		ax1.set_xlabel('Detuning (GHz)')
		ax1.set_ylabel('Transmission')
		ax1.set_xlim(-8,-3)
		ax1.set_ylim(-0.02,1.02)
		plt.show()
	plt.close('all')


def test_Zentile():
	"""
	Reproduce the sub-GHz Faraday filter from Zentile's 2015 Opt. Lett. paper (Fig. 3)
	DOI: 10.1364/OL.40.002000
	"""
	np.random.seed(0)

	d = np.arange(-8000,9000,25) # MHz
	#Voigt
	p_dict = {'Bfield':45.7,'Btheta':0,'lcell':75e-3,'T':67.8,'Dline':'D1','Elem':'Cs'}

	[Iy] = get_spectra(d,[1,0,0],p_dict,outputs=['Iy'])

	assert Iy.shape == d.shape
	assert np.all(np.isfinite(Iy))
	assert np.ptp(Iy.real) > 0.01  # a Faraday filter must actually filter something

	Iy_noisy = Iy + np.random.randn(len(d))*0.005

	if SHOW_PLOTS:
		fig = plt.figure("Faraday filtering, Cs D1")
		yy = 5
		xx = 1
		ax1 = plt.subplot2grid((yy,xx), (0,0), rowspan=yy-1)
		axR = plt.subplot2grid((yy,xx), (yy-1,0), sharex=ax1)
		plt.setp(ax1.get_xticklabels(), visible=False)
		ax1.plot(d/1e3, Iy_noisy, '-', color='k', lw=2.5)
		ax1.plot(d/1e3, Iy, 'r--', lw=2)
		axR.plot(d/1e3, 100*(Iy_noisy - Iy), '-', color='k')
		axR.set_xlabel('Detuning (GHz)')
		ax1.set_ylabel('Transmission')
		axR.set_ylabel(r'R ($\times 100$')
		ax1.set_xlim(-8,9)
		ax1.set_ylim(0,0.8)
		plt.show()
	plt.close('all')


def test_Keaveney():
	"""
	Reproduce the Faraday filter from the 2016 Rev. Sci. Inst. paper (Fig. 2)
	DOI: 10.1063/1.4963230
	"""
	np.random.seed(0)

	d = np.arange(-7000,7000,25) # MHz
	#Voigt
	p_dict = {'Bfield':248.5,'rb85frac':72.17,'Btheta':0,'lcell':5e-3,'T':90.3,'Dline':'D2','Elem':'Rb'}

	[Iy] = get_spectra(d,[1,0,0],p_dict,outputs=['Iy'])
	assert np.all(np.isfinite(Iy))

	Iy2 = Iy**2
	Iy_noisy = Iy + np.random.randn(len(d))*0.005

	p_dict['T'] = 20
	p_dict['Bfield'] = 0
	p_dict['lcell'] = 75e-3
	[S0_25] = get_spectra(d,[1,0,0],p_dict,outputs=['S0'])
	assert np.all(np.isfinite(S0_25))
	assert -1e-6 <= S0_25.real.min() and S0_25.real.max() <= 1.0 + 1e-6

	if SHOW_PLOTS:
		fig = plt.figure("Faraday filtering...")
		yy = 5
		xx = 1
		ax1 = plt.subplot2grid((yy,xx), (0,0), rowspan=yy-1)
		axR = plt.subplot2grid((yy,xx), (yy-1,0), sharex=ax1)
		plt.setp(ax1.get_xticklabels(), visible=False)
		blueish = (0.25, 0.25, 1.0)
		ax1.fill_between(d/1e3, 1, S0_25, color=blueish, alpha=0.25, lw=0)
		ax1.plot(d/1e3, Iy_noisy, '-', color=d_blue, lw=2.5)
		ax1.plot(d/1e3, Iy, 'k--', lw=2.5)
		ax1.plot(d/1e3, Iy2,':', color=d_red, lw=2)
		axR.plot(d/1e3, 100*(Iy_noisy - Iy), '-', color=d_blue)
		axR.set_xlabel('Detuning (GHz)')
		ax1.set_ylabel('Transmission')
		axR.set_ylabel(r'$100R$')
		ax1.set_xlim(-7,7)
		ax1.set_ylim(0,1.0)
		plt.show()
	plt.close('all')


def test_published_figure_anchors():
	"""
	Quantitative anchors for the published-figure configurations above,
	pinned to genuine ElecSus v3.0.8 output (elecsus/tests/golden/,
	generated from upstream commit 065fa28). The sanity checks in the tests
	above catch crashes; these catch physics regressions - a wrong linewidth
	or Zeeman shift changes these numbers immediately.
	"""
	golden = np.load(os.path.join(os.path.dirname(__file__), 'golden',
	                              'es3_reference.npz'))

	# (case, output, reduce, ES3 v3.0.8 value) - values in comments are the
	# figure each configuration reproduces
	anchors = [
		('siddons_rb_d2_16C', 'S0', 'min', 0.594422),   # Siddons 2008 Fig 9
		('siddons_rb_d2_25C', 'S0', 'min', 0.273472),
		('siddons_rb_d2_36C', 'S0', 'min', 0.015748),
		('weller_rb_d2_80C',  'S0', 'min', 0.000000),   # Weller 2011 Fig 5 (opaque)
		('zentile_cs_d1_filter',  'Iy', 'max', 0.764744),  # Zentile 2015 OL Fig 3
		('keaveney_rb_d2_filter', 'Iy', 'max', 0.918866),  # Keaveney 2018 filter
	]
	for case, output, reduce_, es3_val in anchors:
		X = golden[f'{case}__X']
		E_in = golden[f'{case}__E_in']
		p = eval(golden[f'{case}__p'].tobytes().decode())
		[arr] = get_spectra(X, E_in, p, outputs=[output])
		v4_val = getattr(np.asarray(arr).squeeze().real, reduce_)()
		g_val = getattr(golden[f'{case}__{output}'].real, reduce_)()
		# v4 must reproduce the stored ES3 spectrum feature...
		assert abs(v4_val - g_val) < 1e-6, \
			f'{case}: {reduce_}({output}) = {v4_val:.6f}, ES3 gave {g_val:.6f}'
		# ...and the stored value itself must match this file (guards against
		# a silently regenerated/corrupted reference file)
		assert abs(g_val - es3_val) < 1e-5, \
			f'{case}: golden file value {g_val:.6f} != pinned {es3_val:.6f}'


if __name__ == '__main__':
	print('Running Spectroscopy Test Cases...')
	test_Siddons()
	test_Weller()
	test_Zentile()
	test_Keaveney()
	test_published_figure_anchors()
