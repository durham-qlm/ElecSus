# Copyright 2014 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams and I. G. Hughes.

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
Functions to calculate the different number densities

All formulae take in temperature in kelvin and return the
atomic number densities.

Based on values given in the following publication:

C. B. Alcock, V. P. Itkin and M. K. Horrigan,
Can. Metall. Q. 23 (1984) 309-313.

"""

from FundamentalConstants import kB
from numpy import log10

def CalcNumberDensity(T,atom):
	""" Helper function to tidy up code in spectra_SPD.py 
		Calls one of the other functions in this module, based on atom parameter
		
		Temperature in Kelvin
		Number density returned in inverse cubic metres
	"""
	
	if atom in ['Rb85','Rb87','Rb']:
		return numDenRb(T)
	elif atom=='Cs':
		return numDenCs(T)
	elif atom in ['K','K39','K40','K41']:
		return numDenK(T)
	elif atom=='Na':
		return numDenNa(T)


def numDenRb(T):
    """Calculates the rubidium number density"""
    if T<312.46:
        p=10.0**(4.857-4215./T)
    else:
        p=10.0**(8.316-4275./T-1.3102*log10(T))
    NumberDensity=101325.0*p/(kB*T)
    return NumberDensity

def numDenK(T):
    '''Potassium number density'''
    if T<336.8:
        p=10.0**(4.961-4646.0/T)
    else:
        p=10.0**(8.233-4693.0/T-1.2403*log10(T))
    NumberDensity=101325.0*p/(kB*T)
    return NumberDensity

def numDenCs(T):
    '''Caesium number density'''
    if T<301.65:
        p=10.0**(4.711-3999./T)
    else:
        p=10.0**(8.232-4062./T-1.3359*log10(T))
    NumberDensity=101325.0*p/(kB*T)
    return NumberDensity

def numDenNa(T):
    '''Sodium number density'''
    if T<370.95:
        p=10.0**(5.298-5603./T)
    else:
        p=10.0**(8.400-5634./T-1.1748*log10(T))
    NumberDensity=101325.0*p/(kB*T)
    return NumberDensity
