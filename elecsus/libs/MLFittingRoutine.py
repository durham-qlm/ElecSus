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
Marquardt-Levenberg fit module

Updated 20/08/2015 
	Now output bestFitParams includes all parameters that are input in initParams
"""

from numpy import array
from scipy.optimize import curve_fit

import time


def MLfit(xdata,ydata,initParams,paramBools,verbose=True,**kw):
    from spectra import spectrum
    
    x = array(xdata)
    y = array(ydata)
    
    #print 'length of init par:', len(initParams)
    #print 'par bools:', paramBools
	
    functionString = "def ForFit(x,"
    spectrumString = "spec = spectrum(x,initParams[0],initParams[1],"
    initialGuesses = []
    counter = 0
    for boolian in paramBools:
        if boolian and (counter == 0):
            spectrumString = spectrumString + "B,"
            functionString = functionString + "B,"
            initialGuesses.append(initParams[2])
        elif (not boolian) and (counter == 0):
            spectrumString = spectrumString + "initParams[2],"
        if boolian and (counter == 1):
            functionString = functionString + "T,"
            spectrumString = spectrumString + "T,"
            initialGuesses.append(initParams[3])
        elif (not boolian) and (counter == 1):
            spectrumString = spectrumString + "initParams[3],"
        if boolian and (counter == 2):
            functionString = functionString + "LenCell,"
            spectrumString = spectrumString + "LenCell,"
            initialGuesses.append(initParams[4])
        elif (not boolian) and (counter == 2):
            spectrumString = spectrumString + "initParams[4],"
        if boolian and (counter == 3):
            spectrumString = spectrumString + "Rb85fr,"
            functionString = functionString + "Rb85fr,"
            initialGuesses.append(initParams[5])
        elif (not boolian) and (counter == 3):
            spectrumString = spectrumString + "initParams[5],"
        if boolian and (counter == 4):
            functionString = functionString + "DopT,"
            spectrumString = spectrumString + "DopT,"
            initialGuesses.append(initParams[6])
        elif (not boolian) and (counter == 4):
            spectrumString = spectrumString + "initParams[6],"
        if boolian and (counter == 5):
            functionString = functionString + "Theta0,"
            spectrumString = spectrumString + "Theta0,"
            initialGuesses.append(initParams[7])
        elif (not boolian) and (counter == 5):
            spectrumString = spectrumString + "initParams[7],"
        if boolian and (counter == 6):
            functionString = functionString + "Pol,"
            spectrumString = spectrumString + "Pol,"
            initialGuesses.append(initParams[8])
        elif (not boolian) and (counter == 6):
            spectrumString = spectrumString + "initParams[8],"
        if boolian and (counter == 7):
            functionString = functionString + "Sh,"
            spectrumString = spectrumString + "Sh,"
            initialGuesses.append(initParams[9])
        elif (not boolian) and (counter == 7):
            spectrumString = spectrumString + "initParams[9],"
        if boolian and (counter == 8):
            functionString = functionString + "GamBuf,"
            spectrumString = spectrumString + "GamBuf,"
            initialGuesses.append(initParams[10])
        elif (not boolian) and (counter == 8):
            spectrumString = spectrumString + "initParams[10],"
        counter += 1
    functionString = functionString + "):"
    spectrumString = spectrumString + "initParams[11],initParams[12],initParams[13],initParams[14],initParams[15])"
    code = functionString + "\n    " + spectrumString + "\n    " + "return spec"

    exec code in locals()

    Popt, Pcov, infodict, errmsg, ier  = curve_fit(ForFit,x,y,initialGuesses,full_output=True,**kw)
    if verbose: print 'Curvefit routine finished'
	
    spectr = ForFit(x,*Popt)
    if verbose: print 'Spectrum generated'
    #print 'Len spectra: ',len(spectr)
	
    bestFitParams = list(initParams)
	
    ## error in here??
    i = 0
    j = 0
    #print len(paramBools), len(bestFitParams), len(initParams)
    for boo in paramBools:
        if boo:
            bestFitParams[j+2] = Popt[i]
            i+=1
        j+=1
    if verbose: print 'Fit params stored'

    return bestFitParams, spectr