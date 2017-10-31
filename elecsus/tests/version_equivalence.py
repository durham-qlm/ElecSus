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

""" A series of examples, calculating various spectra for different parameter regimes """

import matplotlib.pyplot as plt
import numpy as np

import time
import sys

sys.path.append('../')
from elecsus_methods import calculate as get_spectra
#sys.path.append('../libs')
from libs.durhamcolours import *
#from elecsus.libs.spectra import get_spectra

sys.path.append('C:\Users\James\Documents\Programming\elecsus 2.2\elecsus\libs')
import spectra as spec_v2

def test_transmission():
	""" 
	ETC
	"""
	
	d = np.arange(-4000,5000,15) # MHz
	p_dict = {'Bfield':0,'rb85frac':72.17,'Btheta':0,'lcell':75e-3,'T':16.5,'Dline':'D2','Elem':'Rb'}
	
	# 3.0
	[S0_v3] = get_spectra(d,[1,0,0],p_dict,outputs=['S0'])
	
	A_v3 = -np.log(S0_v3)
	A_v3h = -np.log(S0_v3) / 2
	# 2.2
	[S0_v2] = spec_v2.get_spectra(d,p_dict,outputs=['S0'])
	A_v2 = -np.log(S0_v2)
	
	
	fig = plt.figure("Faraday comparison")
	
	yy = 5
	xx = 1
	ax1 = plt.subplot2grid((yy,xx), (0,0), rowspan=yy-1)
	axR = plt.subplot2grid((yy,xx), (yy-1,0), sharex=ax1)
	
	plt.setp(ax1.get_xticklabels(), visible=False)
	
	ax1.plot(d/1e3, A_v2, 'k-', label='v2.2')
	ax1.plot(d/1e3, A_v3, 'r--', label='v3')
	ax1.plot(d/1e3, A_v3h, 'b--', label='v3, scaled')
	
	ax1.legend(loc=0)
	
	axR.set_xlabel('Detuning (GHz)')
	ax1.set_ylabel('Transmission')
	axR.set_ylabel(r'Difference ($\times 100$)')
	
	ax1.set_xlim(-4,5)
	ax1.set_ylim(0,1.03)
	
	plt.show()
	


if __name__ == '__main__':
	test_transmission()