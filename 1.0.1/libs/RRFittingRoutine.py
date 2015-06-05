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

"""Random restart fitting routine. 

Fit by taking a random sample around parameters and then
fit using Marquardt-Levenberg.

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

def evaluate(arguments):
    warnings.simplefilter("ignore") #Stops the casting complex to real warning
    parametersToFit = arguments[2:-2] 
    bFitPar, spectr = ML.MLfit(arguments[0],arguments[1],
                                parametersToFit,arguments[-2])
    costValue = Cost(arguments[1],spectr)
    sys.stdout.flush()
    return arguments[-1],costValue,bFitPar[0],bFitPar[1],bFitPar[2], \
           bFitPar[3],bFitPar[4],bFitPar[5],bFitPar[6],bFitPar[7], \
           bFitPar[8]

def RRFit(xdata,ydata,initParams,paramBools,noEvals):
    x = N.array(xdata)
    y = N.array(ydata)
    
    #Default uncertainty on values
    if paramBools[0]:
        errBfield     = initParams[2]/10.0 + 1.0
    else:
        errBfield     = 0.0

    if paramBools[1]:
        errNTemp      = (initParams[3]+273.15)/50.0 + 0.1
    else:
        errNTemp      = 0.0

    if paramBools[2]:
        errcellLength = initParams[4]/25.0
    else:
        errcellLength = 0.0

    if paramBools[3]:
        errRb85       = initParams[5]/20.0 + 1.0
    else:
        errRb85       = 0.0

    if paramBools[4]:
        errDopTemp    = (initParams[6]+273.15)/50.0 + 0.1
    else:
        errDopTemp    = 0.0
        
    if paramBools[5]:
        errTheta0     = 0.2
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
        if errBfield > initParams[2]:
            BfieldVals[i] = R.uniform(0,2.0*errBfield)
        else:
            BfieldVals[i] = initParams[2] + errBfield*R.uniform(-1,1)

        if errNTemp > (initParams[3]+273.15):
            NTempVals[i]  = R.uniform(-273.15,2.0*errNTemp)
        else:
            NTempVals[i]  = initParams[3] + errNTemp*R.uniform(-1,1)

        cellLengthVals[i] = initParams[4] + errcellLength*R.uniform(-1,1)

        if errRb85 > initParams[5]:
            Rb85Vals[i]   = R.uniform(0,2.0*errRb85)
        elif (errRb85 + initParams[5]) > 100:
            Rb85Vals[i]   = R.uniform(100.0-2.0*errRb85,100.0)
        else:
            Rb85Vals[i]   = initParams[5] + errRb85*R.uniform(-1,1)

        if errDopTemp > (initParams[6]+273.15):
            DopTempVals[i]= R.uniform(-273.15,2.0*errDopTemp)
        else:
            DopTempVals[i]= initParams[6] + errDopTemp*R.uniform(-1,1)

        Theta0Vals[i]     = initParams[7] + errTheta0*R.uniform(-1,1)

        if errPol > initParams[8]:
            PolVals[i]    = R.uniform(0,2.0*errPol)
        elif (errPol + initParams[8]) > 100:
            PolVals[i]    = R.uniform(100.0-2.0*errPol,100.0)
        else:
            PolVals[i]    = initParams[8] + errPol*R.uniform(-1,1)

        ShiftVals[i]      = initParams[9] + errShift*R.uniform(-1,1)

        if errGamma > initParams[10]:
            GammaVals[i]  = R.uniform(0,2.0*errGamma)
        else:
            GammaVals[i]  = initParams[10] + errGamma*R.uniform(-1,1)
    
    #Do parallel ML fitting by utilising multiple cores
    po = Pool() # Pool() uses all cores, Pool(3) uses 3 cores for example.
    res = po.map_async(evaluate,((x,y,initParams[0],initParams[1],BfieldVals[k],
                                  NTempVals[k],cellLengthVals[k],Rb85Vals[k],
                                  DopTempVals[k],Theta0Vals[k],PolVals[k],ShiftVals[k],
                                  GammaVals[k],initParams[11],initParams[12],
                                  initParams[13],initParams[14],initParams[15],
                                  paramBools,k,) for k in xrange(noEvals)))
    costValues = res.get()
    po.close()
    po.join()

    #Find best fit
    costValues = N.array(costValues)
    costValues = costValues.astype(N.float64)
    lineMin = N.argmin(costValues, axis=0)[1]
    indexmin = costValues[lineMin][0]
    bestCostValue = costValues[lineMin][1]

    bestFitParams = [costValues[lineMin][2],costValues[lineMin][3],costValues[lineMin][4],
                     costValues[lineMin][5],costValues[lineMin][6],costValues[lineMin][7],
                     costValues[lineMin][8],costValues[lineMin][9],costValues[lineMin][10]]
    
    FinalTheory = spectrum(x,initParams[0],initParams[1],bestFitParams[0],
                           bestFitParams[1],bestFitParams[2],bestFitParams[3],
                           bestFitParams[4],bestFitParams[5],bestFitParams[6],
                           bestFitParams[7],bestFitParams[8],initParams[11],
                           initParams[12],initParams[13],
                           initParams[14],initParams[15])
    return bestFitParams, FinalTheory

