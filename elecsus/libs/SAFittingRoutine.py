# Copyright 2014 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams and I. G. Hughes.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#	 http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Program to fit parameters to data using a simulated annealing algorithm.

Author: MAZ
Updated 2015-09-03 - added keyword arg support - Author: JK"""

import numpy as N
import random as R
import MLFittingRoutine as ML
import warnings

from multiprocessing import Pool
from spectra import spectrum

def Cost(yData,yTheory):
	'''
	Function to evaluate the goodness of the fit, larger value is worse!
	'''
	SumDiff = abs(yTheory-yData).sum()
	return SumDiff

def evaluate(arguments):
	warnings.simplefilter("ignore")
	theorySpec = spectrum(arguments[0],*arguments[2:-1])
	costValue = Cost(arguments[1],theorySpec)
	return arguments[-1], costValue

##### Begining of Code ######

def SAFit(xdata,ydata,initParams,paramBools,noEvals=1000,verbose=False,**kw):
	""" Document this method!!! """
	
	x = N.array(xdata)
	y = N.array(ydata)

	Bfield	 = initParams[2]
	NTemp	  = initParams[3]
	cellLength = initParams[4]
	Rb85	   = initParams[5]
	DopTemp	= initParams[6]
	Theta0	 = initParams[7]
	Pol		= initParams[8]
	Shift	  = initParams[9]
	Gamma	  = initParams[10]

	if paramBools[0]:
		errBfield	 = 3*Bfield/100.0 + 1.0
	else:
		errBfield	 = 0.0

	if paramBools[1]:
		errNTemp	  = 3*(NTemp+273.15)/100.0 + 0.1
	else:
		errNTemp	  = 0.0

	if paramBools[2]:
		errcellLength = 3*cellLength/100.0
	else:
		errcellLength = 0.0

	if paramBools[3]:
		errRb85	   = initParams[5]/20.0 + 0.01
	else:
		errRb85	   = 0.0

	if paramBools[4]:
		errDopTemp	= 3*(DopTemp+273.15)/100.0 + 0.1
	else:
		errDopTemp	= 0.0
		
	if paramBools[5]:
		errTheta0	 = 5
	else:
		errTheta0	 = 0.0

	if paramBools[6]:
		errPol		= 3*Pol/100.0 + 0.01
	else:
		errPol		= 0.0

	if paramBools[7]:
		errShift	  = 2*Shift/8.0 + 0.3
	else:
		errShift	  = 0.0
		
	if paramBools[8]:
		errGamma	  = 2*Gamma/10.0 + 0.5
	else:
		errGamma	  = 0.0

	print ''
	print 'Making initial evaluations for scaling parameter\n',

	#Fill arrays of random values scattered around given start
	BfieldVals		= N.zeros(noEvals)
	BfieldVals[0]	 = Bfield
	
	NTempVals		 = N.zeros(noEvals)
	NTempVals[0]	  = NTemp
	
	cellLengthVals	= N.zeros(noEvals)
	cellLengthVals[0] = cellLength 
	
	Rb85Vals		  = N.zeros(noEvals)
	Rb85Vals[0]	   = Rb85
	
	DopTempVals	   = N.zeros(noEvals)
	DopTempVals[0]	= DopTemp
	
	Theta0Vals		= N.zeros(noEvals)
	Theta0Vals[0]	 = Theta0

	PolVals		   = N.zeros(noEvals)
	PolVals[0]		= Pol

	ShiftVals		 = N.zeros(noEvals)
	ShiftVals[0]	  = Shift

	GammaVals		 = N.zeros(noEvals)
	GammaVals[0]	  = Gamma
	
	for i in range(1,noEvals):
		BfieldVals[i]	 = initParams[2] + errBfield*R.uniform(-1,1)
		NTempVals[i]	  = initParams[3] + errNTemp*R.uniform(-1,1)
		cellLengthVals[i] = initParams[4] + errcellLength*R.uniform(-1,1)
		Rb85Vals[i]	   = initParams[5] + errRb85*R.uniform(-1,1)
		DopTempVals[i]	= initParams[6] + errDopTemp*R.uniform(-1,1)
		Theta0Vals[i]	 = initParams[7] + errTheta0*R.uniform(-1,1)
		PolVals[i]		= initParams[8] + errPol*R.uniform(-1,1)
		ShiftVals[i]	  = initParams[9] + errShift*R.uniform(-1,1)
		GammaVals[i]	  = initParams[10] + errGamma*R.uniform(-1,1)
	
	po = Pool() # Pool() uses all cores, Pool(3) uses 3 cores for example.
	res = po.map_async(evaluate,((x,y,initParams[0],initParams[1],BfieldVals[k],
								  NTempVals[k],cellLengthVals[k],Rb85Vals[k],
								  DopTempVals[k],Theta0Vals[k],PolVals[k],ShiftVals[k],
								  GammaVals[k],initParams[11],initParams[12],
								  initParams[13],initParams[14],initParams[15],k) for k in xrange(noEvals)))
	
	Result = res.get()
	po.close()
	po.join()
		
	Result = N.array(Result)
	Result = Result.astype(N.float64)
	lineMin = N.argmin(Result, axis=0)[1] ## pick the fit with the lowest cost value
	indexmin = Result[lineMin][0]
	bestCostValue = Result[lineMin][1]
	
	bestBfield	 = BfieldVals[indexmin]
	bestNTemp	  = NTempVals[indexmin]
	bestcellLength = cellLengthVals[indexmin]
	bestRb85	   = Rb85Vals[indexmin]
	bestDopTemp	= DopTempVals[indexmin]
	bestTheta0	 = Theta0Vals[indexmin]
	bestPol		= PolVals[indexmin]
	bestShift	  = ShiftVals[indexmin]
	bestGamma	  = GammaVals[indexmin]
	
	spread = N.std(Result, axis=0)[1]
	
	costValue = Result[0][1] # starting point

	T = 999.0
	k = abs(spread)/T
	Beta = 0.000005 # This value should be small. Smaller gives better fit, but takes longer.
	breakThreshold = 140 # Number of times a jump uphill is rejected before switching to ML fit.
	probabilities = N.ones(breakThreshold)
	prob_j = 0
	Tthreshold = 5.0 #Cold
	while T > Tthreshold: 
		tryNTemp = NTemp + errNTemp*R.uniform(-1,1) #Vary Parameters
		trycellLength = cellLength + errcellLength*R.uniform(-1,1)
		if trycellLength <= 0:
			trycellLength = 1.0e-9 #very short
		tryRb85 = Rb85 + errRb85*R.uniform(-1,1)
		if tryRb85 > 1:
			tryRb85 = 1.0
		elif tryRb85 < 0:
			tryRb85 = 0.0
		tryBfield = Bfield + errBfield*R.uniform(-1,1)
		if tryBfield < 0:
			tryBfield = 0.0
		tryDopTemp = DopTemp + errDopTemp*R.uniform(-1,1)
		tryTheta0  = Theta0 + errTheta0*R.uniform(-1,1)
		if tryTheta0 < 0.0:
			tryTheta0 = 0.0
		elif tryTheta0 > 1.0:
			tryTheta0 = 1.0
		tryPol = Pol + errPol*R.uniform(-1,1)
		if tryPol < 0.0:
			tryPol = 0.0
		elif tryPol > 1.0:
			tryPol = 1.0
		tryShift = Shift + errShift*R.uniform(-1,1)
		tryGamma = Gamma + errGamma*R.uniform(-1,1)
		if tryGamma < 0:
			tryGamma = 0.0
		newTransTheory = spectrum(x,initParams[0],initParams[1],tryBfield,tryNTemp,
								  trycellLength,tryRb85,tryDopTemp,tryTheta0,tryPol,
								  tryShift,tryGamma,initParams[11],initParams[12],
								  initParams[13],initParams[14],initParams[15])
		newCostValue = Cost(y,newTransTheory)
		if newCostValue < bestCostValue:
			bestNTemp	  = tryNTemp
			bestRb85	   = tryRb85
			bestcellLength = trycellLength
			bestBfield	 = tryBfield
			bestDopTemp	= tryDopTemp
			bestTheta0	 = tryTheta0
			bestPol		= tryPol
			bestShift	  = tryShift
			bestGamma	  = tryGamma
			bestCostValue  = newCostValue
		DeltaCost = costValue - newCostValue
		if verbose:
			print ''
			print 'Current cost:', newCostValue
			print 'Change in cost:', -1*DeltaCost
			print 'T parameter: ', T
			print ''
			if paramBools[0]:
				print 'magnetic field		= ', tryBfield, 'Gauss'
			if paramBools[1]:
				print 'Reservoir temperature = ', tryNTemp, 'degrees Celsius'
			if paramBools[2]:
				print 'Cell length		   = ', trycellLength*1000.0, 'mm'
			if paramBools[3]:
				print 'Rubidium-85 fraction  = ', tryRb85*100.0, '%'
			if paramBools[4]:
				print 'Doppler temperature   = ', tryDopTemp, 'degrees Celsius'
			if paramBools[5]:
				print 'Polarization angle	= ', tryTheta0, 'pi radians'
			if paramBools[6]:
				print 'Sigma minus fraction  = ', tryPol*100, '%'
			if paramBools[7]:
				print 'Global Shift		  = ', tryShift, 'MHz'
			if paramBools[8]:
				print 'Buffer gas broadening = ', tryGamma, 'MHz'
		prob = N.exp(DeltaCost/(k*T))
		if verbose:
			print 'Probability that these will be accepted (>1 == 1):', prob
			print ''

		if DeltaCost >= 0:
			NTemp = tryNTemp
			Rb85 = tryRb85
			cellLength = trycellLength
			Bfield = tryBfield
			DopTemp = tryDopTemp
			Theta0 = tryTheta0
			Pol = tryPol
			Shift = tryShift
			Gamma = tryGamma
			costValue = newCostValue
			if verbose:
				print 'Values accepted.'
				print '\n\n'
		elif prob > R.random():
			NTemp = tryNTemp
			Rb85 = tryRb85
			cellLength = trycellLength
			Bfield = tryBfield
			DopTemp = tryDopTemp
			Theta0 = tryTheta0
			Pol = tryPol
			Shift  = tryShift
			Gamma  = tryGamma
			costValue = newCostValue
			probabilities[prob_j] = 1.0
			prob_j += 1
			if prob_j == breakThreshold:
				prob_j = 0
			if verbose:
				print 'Values accepted.'
				print '\n\n'
		else:
			probabilities[prob_j] = 0.0
			prob_j += 1
			if prob_j == breakThreshold:
				prob_j = 0
			if verbose:
				print 'Values rejected.'
				print '\n\n'
		T = T/(1+(Beta*T)) #Lundy's Method
		if (N.sum(probabilities) < 0.9) or (T < Tthreshold): #No jumps up hill for a while
			print '.. Simulated annealing completed ..'
			print '.... Switching to downhill fit ....\n'
			break

	bestParams = [initParams[0],initParams[1],bestBfield,bestNTemp,
				  bestcellLength,bestRb85,bestDopTemp,bestTheta0,
				  bestPol,bestShift,bestGamma,initParams[11],
				  initParams[12],initParams[13],initParams[14],
				  initParams[15]]
	bestFitParams = [bestBfield,bestNTemp,bestcellLength,bestRb85,
				  bestDopTemp,bestTheta0,bestPol,bestShift,bestGamma]

	#### Marquardt-Levenberg fit #####
	try:
		Poptimal, Spectr = ML.MLfit(xdata,ydata,bestParams,paramBools,**kw)
		newCostValue = Cost(y,Spectr)
		if newCostValue < bestCostValue:
			bestParams = Poptimal
			FinalTheory = Spectr
			print 'Downhill fit converged successfully.\n'
		else:
			print 'Downhill fit did not find further improvement. Continuing with simulated annealing result.\n'
			
	except:
		print 'Downhill fit failed to converge. Continuing with simulated annealing result.\n'

	FinalTheory = spectrum(x,*bestParams)
	return bestParams, FinalTheory