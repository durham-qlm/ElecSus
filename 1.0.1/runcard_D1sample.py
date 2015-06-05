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

"""Run card.

Format should be:

Variable = input

Change the inputs, not the variables."""

#### Essential options ####
########################### 

#Alkali element       # Choices are 'Rb' for rubidium, 'Cs' for caesium,
element = 'Rb'        # 'K' for potassium, 'Na' for sodium

#Which transition? (D1 or D2)
Transition = 'D1' 

#Spectrum type?       # Choices: 'S0' for transmission, 'S1' for Ix - Iy,
Spectrum = 'S0'       #          'S2' similar to S1, 
                      #          'S3' difference between sigma transmissions,
                      #          'Ix' Intensity through PBS (FADOF),
                      #          'Iy' Intensity reflected at PBS,
                      #          'RI+' Sigma+ refractive index,
                      #          'RI-' Sigma- refractive index,
                      #          'GI+' Sigma+ group index,
                      #          'GI-' Sigma- group index.

#Do you want to fit to experiment or just a theoretical presiction?
FitType = 'RR'        # Choices: 'T' for just theory (later fitting options ignored)
                      #          'ML' for Marquardt-Levenberg fit (standard fitting)
                      #          'RR' for Random restart (3 or more fit parameters)
                      #          'SA' for Simulated annealing (many parameters) 

#Experimental data file name?
DataFilename = 'SampleDataRbD1' #For theory plots this is only used as a label for the data.

#### Essential cell parameters ####
###################################

#Temperature 
NdenTemp = 130.0      #in degrees Celsius.
NdenTempBool = True   #Choose to fit? (True or False) (NOT true or false)

#Cell length
CellLength = 1.0      #in millimetres
CellLenBool = False   #Choose to fit? (True or False) (NOT true or false)

##################################################
######### END of essential input section #########
######### less commonly needed options   #########
######### and parameters are given below #########
##################################################

#### Less commonly used parameters ####
#######################################

#Magnetic field
Bfield = 1000.0        #in Gauss
BfieldBool = True      #Choose to fit? (True or False) (NOT true or false)

#Initial polarization angle
Theta0 = 45.0          #in degrees. Only important for S1, S2, Ix, Iy spectrum types
Theta0Bool = False     #Choose to fit? (True or False) (NOT true or false)

#Initial polarization (Change if not using linearly polarized light)
Polar = 50.0           #As a percentage of light that drives sigma minus transitions
PolarBool = False      #Choose to fit? (True or False)

#Shift
Shift = 0.0            #Global shift from linecentre in MHz
ShiftBool = False      #Choose whether to fit or not (True or False)

#Extra lorentzian broadening (buffer gas etc.)
Gamma = 20.0           #Lorentzian width due to buffer gas, in MHz
GammaBool = True       #Choose whether to fit or not (True or False)

#Doppler temperature
ConstrainBool = True  #Constain the doppler temp to the number density temp?
DoppTemp = 132.0      #in deg C. Overided by the number density temp if constrained
DoppTempBool = False  #Choose to fit? (True or False). Choose False if constrained

## Isotopic abundances ##
#########################

#Values given for elements not being calculated are not important,
#but do not leave blank.

#Rubidium 85 percentage
Rb85 = 1.0             #in %. 72.17 is the natural abundance
Rb85Bool = False       #Choose to fit? (True or False)

#Potassium-41 and 40 percentages
K41 = 6.7302           #in %. Use K41 = 6.7302 and K40 = 0.0117 for natural abundance
K40 = 0.0117
#Not possible to fit in this version.

#### Further Options ####
#########################

#Precision to which the theoretical calulation is performed in MHz.
Precision  = 10.0      #10.0 is recommended.

#Smooth input data? (True or False)  (NOT true or false)
SmoothBool = True      #Not important for theory only.
SmoothFactor = 19      #Number of local points to average over (must be an odd number)

#Detuning range and offset (not looked at if experimental data is fitted to)
DetStart = -10.0       #Detuning start value in GHz
DetStop  = 10.0        #Detuning end value in GHz

#Show plot on screen?
PlotBool = True #Set to False if working remotely

#Save plot?
SavePlotBool = True #(True or False)  (NOT true or false)

#Plot file type? ('.pdf', '.eps', '.png' etc)
PlotFileType = '.png'  #Not important when not saving plots

#Plot hyperfine stick spectra for no magnetic field?
#Useful when comparing to sub-doppler features for calibration.
PlotHypBool = False #(True or False)
