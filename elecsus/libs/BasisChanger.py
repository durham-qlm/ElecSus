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
Change electric field basis between xyz and +/-/z bases.
this module provides methods to convert in both directions

JK
"""
import numpy as np

def xyz_to_lrz(E_in):
	""" Convert from linear to circular bases """
	# create output array
	E_out = np.zeros_like(E_in,dtype='complex')
	
	# z-component doesn't change
	E_out[2] = E_in[2]
	
	## Following sign convention in 
	## 'Optically Polarised Atoms' by Auzinsh, Budker and Rochester, eq 6.32
	## OUP, 2010
	# L = 1./sqrt(2) * (x - iy)
	# R = 1./sqrt(2) * (x + iy) 
	E_out[0] = 1./np.sqrt(2) * (E_in[0] - 1.j*E_in[1])
	E_out[1] = 1./np.sqrt(2) * (E_in[0] + 1.j*E_in[1])
	
	return E_out
	
def lrz_to_xyz(E_in):
	""" Convert from circular to linear bases """

	# create output array
	E_out = np.zeros_like(E_in,dtype='complex')
	
	# z-component doesn't change
	E_out[2] = E_in[2]

	# x = 1. / sqrt(2) * [L + R]
	# y = 1.j / sqrt(2) * [L - R]
	E_out[0] = 1./np.sqrt(2) * (E_in[0] + E_in[1])
	E_out[1] = 1.j/np.sqrt(2) * (E_in[0] - E_in[1])
	
	return E_out