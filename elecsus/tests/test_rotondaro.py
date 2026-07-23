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
A series of plots to test ElecSus for the parameters specified in the Rotondaro JOSA B 2015 paper.

Converted to pytest (issue #26, 2026-07-09): gated interactive plotting behind
ELECSUS_SHOW_PLOTS, reduced point counts on the two slow oblique-field cases for
CI speed, and added physical sanity assertions (finite, correctly shaped,
bounded output). Physics/logic otherwise unchanged from the original 2018 version.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from elecsus.libs.spectra import get_spectra

SHOW_PLOTS = os.environ.get('ELECSUS_SHOW_PLOTS', '') == '1'


def _assert_sane_spectrum(arr, d, bound=1.2):
	"""Common sanity checks for a returned spectrum array."""
	arr = np.asarray(arr)
	assert arr.shape == d.shape
	assert np.all(np.isfinite(arr))
	assert np.max(np.abs(arr)) < bound


def test_faraday_transmission():
	"""
	1. Fig 3 of Generalised treatment ... Rotondaro JOSAB 2015 paper
	Normal Faraday spectrum
	"""

	d = np.arange(-10000,10000,10) # MHz
	#Voigt
	p_dict = {'Bfield':300,'rb85frac':1,'Btheta':0,'lcell':75e-3,'T':58,'Dline':'D2','Elem':'Cs'}

	TF = get_spectra(d,[1,0,0],p_dict,outputs=['Iy'])

	_assert_sane_spectrum(TF[0], d)
	# Faraday rotation signal should not be flat - the whole point of the filter
	assert np.ptp(TF[0].real) > 0.01

	if SHOW_PLOTS:
		fig = plt.figure("Faraday comparison")
		ax1 = fig.add_subplot(111)
		ax1.plot(d,TF[0],'r',lw=2,label='Faraday')
		ax1.set_xlabel('Detuning (MHz)')
		ax1.set_ylabel('Transmission')
		plt.show()
	plt.close('all')


def test_voigt_filter():
	"""
	2. Fig 4/5 of Rotondaro paper
	Voigt Filter
	"""

	d = np.linspace(-15000,15000,2500)
	## 700 G, 84 C, Cs, 75mm
	p_dict = {'Bfield':700,'Btheta':90*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	pol = 1./np.sqrt(2)*np.array([1.0,1.0,0.0])
	TVx = get_spectra(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0'])

	for arr in TVx:
		_assert_sane_spectrum(arr, d)

	# I_+45 + I_-45 should equal Ix + Iy (energy conservation between the two bases)
	np.testing.assert_allclose((TVx[0]+TVx[1]).real, (TVx[2]+TVx[3]).real, atol=1e-6)

	if SHOW_PLOTS:
		fig2 = plt.figure()
		ax1a = fig2.add_subplot(411)
		ax2a = fig2.add_subplot(412,sharex=ax1a)
		ax3a = fig2.add_subplot(413,sharex=ax1a)
		ax4a = fig2.add_subplot(414,sharex=ax1a)
		ax1a.plot(d,TVx[0],'r',lw=2,label=r'$I_{-45}$')
		ax2a.plot(d,TVx[1],'b',lw=2,label=r'$I_{+45}$')
		ax3a.plot(d,TVx[2],'r',lw=2,label=r'$I_x$')
		ax4a.plot(d,TVx[3],'b',lw=2,label=r'$I_y$')
		ax4a.plot(d,TVx[0]+TVx[1],'r:',lw=3.5,label=r'$I_{+45}+I_{-45}$')
		ax4a.plot(d,TVx[2]+TVx[3],'k:',lw=2.5,label=r'$I_x + I_y$')
		ax4a.plot(d,TVx[4],'g--',lw=1.5,label='$S_0$')
		ax4a.set_xlabel('Detuning (MHz)')
		ax1a.set_ylabel('I -45')
		ax2a.set_ylabel('I +45')
		ax3a.set_ylabel('Ix')
		ax4a.set_ylabel('Iy')
		ax4a.set_xlim(d[0],d[-1]+3000)
		ax4a.legend(loc=0)
		plt.show()
	plt.close('all')


def test_arbitrary_filter_unoptimised():
	"""
	3. Fig 7 of Rotondaro paper
	Arbitrary Filter - non-optimised

	Point count reduced from the original 300 to 60 for CI speed; the oblique
	(Btheta=87deg) calculation path is the same either way.
	"""

	d = np.linspace(-15000,15000,60)
	p_dict = {'Bfield':500,'rb85frac':1,'Btheta':87*np.pi/180,'Bphi':00*np.pi/180,'lcell':75e-3,'T':100,'Dline':'D2','Elem':'Cs'}
	pol = np.array([1.0,0.0,0.0])
	TVx = get_spectra(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0'])

	for arr in TVx:
		_assert_sane_spectrum(arr, d)

	if SHOW_PLOTS:
		fig2 = plt.figure()
		ax1a = fig2.add_subplot(411)
		ax2a = fig2.add_subplot(412,sharex=ax1a)
		ax3a = fig2.add_subplot(413,sharex=ax1a)
		ax4a = fig2.add_subplot(414,sharex=ax1a)
		ax1a.plot(d,TVx[0],'r',lw=2,label=r'$I_{-45}$')
		ax2a.plot(d,TVx[1],'b',lw=2,label=r'$I_{+45}$')
		ax3a.plot(d,TVx[2],'r',lw=2,label=r'$I_x$')
		ax4a.plot(d,TVx[3],'b',lw=2,label=r'$I_y$')
		ax4a.plot(d,TVx[0]+TVx[1],'r:',lw=3.5,label=r'$I_{+45}+I_{-45}$')
		ax4a.plot(d,TVx[2]+TVx[3],'k:',lw=2.5,label=r'$I_x + I_y$')
		ax4a.plot(d,TVx[4],'g--',lw=1.5,label='$S_0$')
		ax4a.set_xlabel('Detuning (MHz)')
		ax1a.set_ylabel('I -45')
		ax2a.set_ylabel('I +45')
		ax3a.set_ylabel('Ix')
		ax4a.set_ylabel('Iy')
		ax4a.set_xlim(d[0],d[-1]+3000)
		ax4a.legend(loc=0)
		plt.show()
	plt.close('all')


def test_arbitrary_filter_optimised():
	"""
	4. Fig 8 of Rotondaro paper
	Arbitrary Filter - optimised

	Point count reduced from the original 300 to 60 for CI speed.
	"""

	d = np.linspace(-15000,15000,60)
	p_dict = {'Bfield':1000,'rb85frac':1,'Btheta':88*np.pi/180,'Bphi':00*np.pi/180,'lcell':75e-3,'T':93,'Dline':'D2','Elem':'Cs'}
	pol = np.array([1.0,0.0,0.0])
	TVx = get_spectra(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0'])

	for arr in TVx:
		_assert_sane_spectrum(arr, d)

	if SHOW_PLOTS:
		fig2 = plt.figure()
		ax1a = fig2.add_subplot(411)
		ax2a = fig2.add_subplot(412,sharex=ax1a)
		ax3a = fig2.add_subplot(413,sharex=ax1a)
		ax4a = fig2.add_subplot(414,sharex=ax1a)
		ax1a.plot(d,TVx[0],'r',lw=2,label=r'$I_{-45}$')
		ax2a.plot(d,TVx[1],'b',lw=2,label=r'$I_{+45}$')
		ax3a.plot(d,TVx[2],'r',lw=2,label=r'$I_x$')
		ax4a.plot(d,TVx[3],'b',lw=2,label=r'$I_y$')
		ax4a.plot(d,TVx[0]+TVx[1],'r:',lw=3.5,label=r'$I_{+45}+I_{-45}$')
		ax4a.plot(d,TVx[2]+TVx[3],'k:',lw=2.5,label=r'$I_x + I_y$')
		ax4a.plot(d,TVx[4],'g--',lw=1.5,label='$S_0$')
		ax4a.set_xlabel('Detuning (MHz)')
		ax1a.set_ylabel('I -45')
		ax2a.set_ylabel('I +45')
		ax3a.set_ylabel('Ix')
		ax4a.set_ylabel('Iy')
		ax4a.set_xlim(d[0],d[-1]+3000)
		ax4a.legend(loc=0)
		plt.show()
	plt.close('all')


def test_rotondaro_filter_anchor():
	"""
	Quantitative anchor for the Rotondaro JOSA B 2015 Fig 3 Faraday-filter
	configuration, pinned to genuine ElecSus v3.0.8 output
	(elecsus/tests/golden/es3_reference.npz, upstream commit 065fa28).
	Covers open issue #27 at the level of ES3 equivalence.
	"""
	golden = np.load(os.path.join(os.path.dirname(__file__), 'golden',
	                              'es3_reference.npz'))
	case = 'rotondaro_cs_d2_faraday'
	X = golden[f'{case}__X']
	E_in = golden[f'{case}__E_in']
	p = eval(golden[f'{case}__p'].tobytes().decode())

	[Iy] = get_spectra(X, E_in, p, outputs=['Iy'])
	Iy = np.asarray(Iy).squeeze().real
	g_Iy = golden[f'{case}__Iy'].real

	# full-spectrum agreement with ES3, plus human-readable peak values
	np.testing.assert_allclose(Iy, g_Iy, atol=1e-7, rtol=0)
	assert abs(Iy.max() - 0.985282) < 1e-5          # peak filter transmission
	assert abs(X[Iy.argmax()] / 1e3 - 2.735) < 0.1  # peak position (GHz)


if __name__ == '__main__':
	test_faraday_transmission()
	test_voigt_filter()
	test_arbitrary_filter_unoptimised()
	test_arbitrary_filter_optimised()
	test_rotondaro_filter_anchor()
