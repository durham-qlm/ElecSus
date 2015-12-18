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

Random restart fitting routine. 

Fit by taking a random sample around parameters and then
fit using Marquardt-Levenberg.

Author: MAZ
Updated 20/08/2015 to reflect changes in MLFittingRoutine (Author: JK)

"""

import numpy as N
import random as R
import warnings
import sys
from multiprocessing import Pool

from spectra import spectrum
import MLFittingRoutine as ML

def Cost(yData,yTheory):
    '''
    Function to evaluate the goodness of the fit, larger value is worse!
    '''
    #Sum of absolute value of difference
    SumDiff = N.sum(N.absolute(yTheory-yData))
    return SumDiff

def evaluate(args):
    #arguments = args
    arguments, kw = args[:-1], args[-1]
    warnings.simplefilter("ignore") #Stops the casting complex to real warning
    parametersToFit = arguments[2:-2] 
    bFitPar, spectr = ML.MLfit(arguments[0],arguments[1],
                                parametersToFit,arguments[-2],**kw)
    print 'Eval_ML COmplete'
	
    costValue = Cost(arguments[1],spectr)
    #sys.stdout.flush()
    #print arguments[-1],costValue,bFitPar[2],bFitPar[3],bFitPar[4], \
    #       bFitPar[5],bFitPar[6],bFitPar[7],bFitPar[8],bFitPar[9], \
    #       bFitPar[10]
    return arguments[-1],costValue,bFitPar[2],bFitPar[3],bFitPar[4], \
           bFitPar[5],bFitPar[6],bFitPar[7],bFitPar[8],bFitPar[9], \
           bFitPar[10]

def RRFit(xdata,ydata,initParams,paramBools,noEvals,**kw):
    """ 
	Write some stuff about this here...
	
	kwargs must be dictionary 
	"""
	
    print '\n\nStarting Random Restart Fitting Routine'
    x = N.array(xdata)
    y = N.array(ydata)
    
    #Default uncertainty on values
    if paramBools[0]:
        errBfield     = 1.5*initParams[2]/10.0 + 1.0
    else:
        errBfield     = 0.0

    if paramBools[1]:
        errNTemp      = 4*(initParams[3]+273.15)/100.0 + 0.1
    else:
        errNTemp      = 0.0

    if paramBools[2]:
        errcellLength = initParams[4]/25.0
    else:
        errcellLength = 0.0

    if paramBools[3]:
        errRb85       = initParams[5]/20.0 + 0.01
    else:
        errRb85       = 0.0

    if paramBools[4]:
        errDopTemp    = (initParams[6]+273.15)/50.0 + 0.1
    else:
        errDopTemp    = 0.0
        
    if paramBools[5]:
        errTheta0     = 0.01
    else:
        errTheta0     = 0.0

    if paramBools[6]:
        errPol        = initParams[8]/10.0 + 0.01
    else:
        errPol        = 0.0
    
    if paramBools[7]:
        errShift      = initParams[9]/4.0 + 0.5
    else:
        errShift      = 0.0
        
    if paramBools[8]:
        errGamma      = initParams[10]/5.0 + 1.0
    else:
        errGamma      = 0.0
    
    #Fill arrays of random values scattered around given start
    BfieldVals = N.zeros(noEvals)
    BfieldVals[0] = initParams[2]
    
    NTempVals = N.zeros(noEvals)
    NTempVals[0] = initParams[3]
    
    cellLengthVals = N.zeros(noEvals)
    cellLengthVals[0] = initParams[4] 
    
    Rb85Vals = N.zeros(noEvals)
    Rb85Vals[0] = initParams[5]
    
    DopTempVals = N.zeros(noEvals)
    DopTempVals[0] = initParams[6]
    
    Theta0Vals = N.zeros(noEvals)
    Theta0Vals[0] = initParams[7]

    PolVals = N.zeros(noEvals)
    PolVals[0] = initParams[8]
    
    ShiftVals = N.zeros(noEvals)
    ShiftVals[0] = initParams[9]
    
    GammaVals = N.zeros(noEvals)
    GammaVals[0] = initParams[10]       
    
    for i in xrange(1,noEvals):
        BfieldVals[i]     = initParams[2] + errBfield*R.uniform(-1,1)
        NTempVals[i]      = initParams[3] + errNTemp*R.uniform(-1,1)
        cellLengthVals[i] = initParams[4] + errcellLength*R.uniform(-1,1)
        Rb85Vals[i]       = initParams[5] + errRb85*R.uniform(-1,1)
        DopTempVals[i]    = initParams[6] + errDopTemp*R.uniform(-1,1)
        Theta0Vals[i]     = initParams[7] + errTheta0*R.uniform(-1,1)
        PolVals[i]        = initParams[8] + errPol*R.uniform(-1,1)
        ShiftVals[i]      = initParams[9] + errShift*R.uniform(-1,1)
        GammaVals[i]      = initParams[10] + errGamma*R.uniform(-1,1)
    
    #Do parallel ML fitting by utilising multiple cores
    po = Pool() # Pool() uses all cores, Pool(3) uses 3 cores for example.
    res = po.map_async(evaluate,((x,y,initParams[0],initParams[1],BfieldVals[k],
                                  NTempVals[k],cellLengthVals[k],Rb85Vals[k],
                                  DopTempVals[k],Theta0Vals[k],PolVals[k],ShiftVals[k],
                                  GammaVals[k],initParams[11],initParams[12],
                                  initParams[13],initParams[14],initParams[15],
                                  paramBools,k, kw) for k in xrange(noEvals))
								  )
    Results = res.get()
    po.close()
    po.join()
    print 'RR calculation complete'

    #Find best fit
    Results = N.array(Results)
    Results = Results.astype(N.float64)
    lineMin = N.argmin(Results, axis=0)[1] ## pick the fit with the lowest cost value
    indexmin = Results[lineMin][0]
    bestCostValue = Results[lineMin][1]

    bestFitParams = [initParams[0],initParams[1],
                     Results[lineMin][2],Results[lineMin][3],Results[lineMin][4],
                     Results[lineMin][5],Results[lineMin][6],Results[lineMin][7],
                     Results[lineMin][8],Results[lineMin][9],Results[lineMin][10],
					 initParams[11],initParams[12],initParams[13],initParams[14],initParams[15]]
    
    FinalTheory = spectrum(x,*bestFitParams)
	
    return bestFitParams, FinalTheory

