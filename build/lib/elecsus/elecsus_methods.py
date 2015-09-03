# Copyright 2015 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
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

"""ElecSus: Program to calculate the electric susceptibility

This module takes user inputs from the run card (default runcard.py),
and outputs data, performs fitting and shows plots.

Usage:

.......
.........
........
.......	
	

"""

import os
import sys
import inspect
import warnings
from shutil import copyfile
from numpy import arange, zeros, array, sqrt

# import elecsus modules
from libs import spectraEfield
from libs import runcardCheck
from libs.tools import fileOutput, plotOutput, smoother, read_in_twoColumn
from libs.numberDensityEqs import *

import libs.MLFittingRoutine as ML
import libs.SAFittingRoutine as SA
import libs.RRFittingRoutine as RR

if os.name == 'posix':
	from time import time as timing #Timing for linux or apple
else:
	from time import clock as timing #Timing for windows

#Stop warnings about casting complex numbers
warnings.simplefilter("ignore")

def calculate(detuning_range, p, OutputType='All'):
	""" 
	Detuning range in GHz.
	p is a list of parameters comprising (in order):
	[Element, Dline, Bfield, T, lcell, rb85frac, DoppTemp, Theta0, Pol, shift, GammaBuf,
	ConstrainDoppT, K40frac, K41frac]
	
	Returns list of arrays (in order):
	
	S0,S1,S2,S3,Ix,Iy,nplus,nminus,phi,alphaplus,alphaminus
	
	"""
	startTime = timing()
	spec_data = spectraEfield.spectrum(detuning_range, # convert to MHz
					Elem=p[0], Dline=p[1], Bfield=p[2], T=p[3], lcell=p[4],
					rb85frac=p[5], DoppTemp=p[6], theta0=p[7], Pol=p[8], shift=p[9],
					GammaBuf=p[10], Constrain=p[11], K40frac=p[12], K41frac=p[13],
					OutputType=OutputType)
	
	print 'Time taken (Calculation only):', timing() - startTime

	return spec_data
	
def fit_data(data,parameters,paramBoolList,experimental_datatype='S0',fit_algorithm='ML',**kw):
	""" 
	Blurb about this:
	
	data is a 2-element list containing x and y data 1-d arrays.
	parameters (and parameter order) is same as for the calculate() method.
	paramBoolList is ...
	experimental_datatype is ...
	fit_algorith is ...
	keywords are
	"""
	
	
	## alter the parameter order. again.
	## Really need to rewrite the fitting modules for a more sensible parameter order
	## but it's a pain and doesn't add any functionality, so
	## it is therefore on the 'to do' list!
	
	## Order as given in the 'parameters' argument
	## Element, Dline, B, T, L, Rb85, DoppT, Theta0, Pol, shift,
	## GammaBuf, Constrain, K40, K41
	
	## Order required by Fitting routines:
	## Element, OutputType, B, T, L, Rb85%, DoppT, Theta0, Pol, Shift, 
	## GammaBuf, Constrain, Dline, Precision, K40%, K41%

	parameters_new = [None]*16
	parameters_new[0] = parameters[0]
	parameters_new[1] = experimental_datatype
	parameters_new[2:12] = parameters[2:12]
	parameters_new[12] = parameters[1]
	parameters_new[13] = 10
	parameters_new[14:] = parameters[12:]
	
	print parameters_new
	
	
	startTime = timing()
	xdata, ydata = data
	print xdata, ydata
	print type(xdata), type(ydata)
	
	
	# Call different fitting routines        
	if fit_algorithm == 'Marquardt-Levenberg':
		#print '\nPerfoming Marquardt-Levenberg fitting routine.'
		optParams, Spec = ML.MLfit(xdata,ydata,parameters_new,
												 paramBoolList,**kw)
	elif fit_algorithm == 'Simulated Annealing':
		#print '\nPerforming fitting by simulated annealing.'
		optParams, Spec = SA.SAFit(xdata,ydata,parameters_new,
												 paramBoolList,**kw)
	else:
		print '\nPerforming fitting by Random-Restart hill climbing method.'
		#The more parameters to fit, the more evaluations we need to do.
		factor = sum(paramBoolList) 
		evaluationNumber = factor**2 + 5 #integer
		optParams, Spec = RR.RRFit(xdata,ydata,parameters_new,paramBoolList,
								   evaluationNumber,**kw)
	
	
	RMS = sqrt(((ydata - Spec)**2).sum()/float(len(ydata)))
	
	# Write the fit parameters to a file
	parameterLabels = ['Magnetic field in Gauss =',
					   'Reservoir temperature in Celsius =',
					   'Cell Length in mm =',
					   'Rb85 percentage =',
					   'Doppler temperature in Celsius =',
					   'Theta0 in degrees =',
					   'Initial sigma minus polarisation percentage = ',
					   'Shift in MHz =', 
					   'Extra Lorentzian width broadening (MHz) =',
					   'Potassium-40 percentage =',
					   'Potassium-41 percentage =']
	#f_Parameters = open(os.path.join(outputDirectory,DAT+'_Parameters.txt')
	#					,'w')
	#optParams[2] = optParams[2]
	#optParams[3] = optParams[3]
	#optParams[6] = optParams[6]
	#optParams[5] = optParams[5]

	## re-order the optimal params back to the way they were given in the 'parameters' argument....
	## this is getting pretty tedious...
	optParams_out = [0]*len(parameters)
	optParams_out[0] = optParams[0]
	optParams_out[1] = optParams[12]
	optParams_out[2:12] = optParams[2:12]
	optParams_out[12:] = optParams[14:]
	
	print 'Optimum parameters found !'
	
	return optParams_out, RMS