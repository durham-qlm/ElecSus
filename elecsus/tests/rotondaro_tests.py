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
A series of plots to test ElecSus for the parameters specified in the Rotondaro JOSA B 2015 paper 

Last updated 2018-02-19 JK
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)


import matplotlib.pyplot as plt
import numpy as np

import time

from elecsus.libs.spectra import get_spectra


def test1():
	""" 
	1. Fig 3 of Generalised treatment ... Rotondaro JOSAB 2015 paper
	Normal Faraday spectrum
	"""
	
	d = np.arange(-10000,10000,10) # MHz
	#Voigt
	p_dict = {'Bfield':300,'rb85frac':1,'Btheta':0,'lcell':75e-3,'T':58,'Dline':'D2','Elem':'Cs'}
	
	#timing:
	st = time.clock()
	TF = get_spectra(d,[1,0,0],p_dict,outputs=['Iy'])
	et = time.clock() - st
	print(('E-field - Elapsed time (s):', et))

	'''
	#check vs old elecsus
	from elecsus_v2.libs import spectra as old_spec
	st = time.clock()
	TF_old = old_spec.get_spectra(d,p_dict,outputs=['Iy'])
	et = time.clock() - st
	print 'Old elecsus - Elapsed time (s):', et
	'''
	
	fig = plt.figure("Faraday comparison")
	ax1 = fig.add_subplot(111)
	ax1.plot(d,TF[0],'r',lw=2,label='Faraday')
	#ax1.plot(d,TF_old[0],'k--',lw=2,label='Vanilla ElecSus')
	
	#ax1.legend(loc=0)
	
	ax1.set_xlabel('Detuning (MHz)')
	ax1.set_ylabel('Transmission')
	
	plt.show()

def test2():
	""" 
	2. Fig 4/5 of Rotondaro paper
	Voigt Filter
	"""
	
	d = np.linspace(-15000,15000,2500)
	#Voigt


	## 700 G, 84 C, Cs, 75mm

	p_dict = {'Bfield':700,'Btheta':90*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	pol = 1./np.sqrt(2)*np.array([1.0,1.0,0.0])
	TVx = get_spectra(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0','Iz'])
	
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
	#	ax4a.plot(d,TVx[5],'c--',lw=2.5,label='$I_z$')
	
	
	ax4a.set_xlabel('Detuning (MHz)')
	ax1a.set_ylabel('I -45')
	ax2a.set_ylabel('I +45')
	ax3a.set_ylabel('Ix')
	ax4a.set_ylabel('Iy')
	
	ax4a.set_xlim(d[0],d[-1]+3000)
	ax4a.legend(loc=0)
	
	plt.show()

def test3():
	"""
	3. Fig 7 of Rotondaro paper
	Arbitrary Filter - non-optimised
	"""
	
	print('This takes a while to compute - be patient!')
	
	d = np.linspace(-15000,15000,300)
	p_dict = {'Bfield':500,'rb85frac':1,'Btheta':87*np.pi/180,'Bphi':00*np.pi/180,'lcell':75e-3,'T':100,'Dline':'D2','Elem':'Cs'}
	pol = np.array([1.0,0.0,0.0])
	TVx = get_spectra(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0','Iz'])
	
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
#	ax4a.plot(d,TVx[5],'c--',lw=2.5,label='$I_z$')
	
	
	ax4a.set_xlabel('Detuning (MHz)')
	ax1a.set_ylabel('I -45')
	ax2a.set_ylabel('I +45')
	ax3a.set_ylabel('Ix')
	ax4a.set_ylabel('Iy')
	
	ax4a.set_xlim(d[0],d[-1]+3000)
	ax4a.legend(loc=0)
	
	plt.show()

def test4():
	"""
	4. Fig 8 of Rotondaro paper
	Arbitrary Filter - optimised
	"""
	
	print('This takes a while to compute - be patient!')
	
	d = np.linspace(-15000,15000,300)
	#Voigt
	#p_dict = {'Bfield':700,'rb85frac':1,'Btheta':90*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	p_dict = {'Bfield':1000,'rb85frac':1,'Btheta':88*np.pi/180,'Bphi':00*np.pi/180,'lcell':75e-3,'T':93,'Dline':'D2','Elem':'Cs'}
	pol = np.array([1.0,0.0,0.0])
	TVx = get_spectra(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0','Iz'])
	
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
#	ax4a.plot(d,TVx[5],'c--',lw=2.5,label='$I_z$')
	
	
	ax4a.set_xlabel('Detuning (MHz)')
	ax1a.set_ylabel('I -45')
	ax2a.set_ylabel('I +45')
	ax3a.set_ylabel('Ix')
	ax4a.set_ylabel('Iy')
	
	ax4a.set_xlim(d[0],d[-1]+3000)
	ax4a.legend(loc=0)
	
	plt.show()

if __name__ == '__main__':
	test1()
	test2()
	test3()
	test4()