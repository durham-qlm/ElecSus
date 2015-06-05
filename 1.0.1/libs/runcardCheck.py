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

"""Module to check user inputs from runcard"""

def value_type_check(label,obj):
    try:
        obj = float(obj)
    except:
        print "ERROR: Value for", label, "not recognised."
        print 'Please input as a number (float or integer).'
        exit(1)
    return obj

def boolean_type_check(label,obj):
    if obj in ['true','True', 'TRUE']:
        obj = True
    elif obj in ['false', 'False', 'FALSE']:
        obj = False
    elif not (type(obj) is bool):
        print 'ERROR:', label, 'option not recognised.'
        print "Please input as a python boolean or string"
        exit(1)
    return obj

def string_type_check(label,obj):
    if not (type(obj) is str):
        print 'ERROR:', label, 'choice not recognised.'
        print "Please input as a python string (in quotation marks)."
        exit(1)
    return obj

def sanity_check(Elem,NTemp,NTempBool,ConstrainBool,DopplerTemp,DopplerTempBool,
                 Bfield,BfieldBool,Rb85frac,Rb85fracBool,K40frac,K41frac,theta0,
                 theta0Bool,cellLength,cellLengthBool,Polar,PolarBool,shift,
                 shiftBool,Gamma,GammaBool,precision,Transition,SpectrumType,
                 DataFilename,SmoothBool,SmoothFactor,FitType,DetStart,
                 DetStop,PlotBool,SavePlotBool,PlotFileType,PlotHypBool):
    
    #Check if inputs are the correct object type and rectify if possible.
    Elem = string_type_check('Element',Elem)

    if not Elem in ['Rb', 'Cs', 'K', 'Na']:
        print 'ERROR: Element type not recognised.\n'
        print "Options are: 'Rb', 'Cs', 'K' or 'Na', corresponding"
        print "to rubidium, caesium or potassium\n"
        print "See manual for more details."

    if Elem == 'Rb':
        K40frac = 0.0
        K41frac = 0.0
        Rb85frac = value_type_check('rubidium-85 fraction', Rb85frac)
        Rb85fracBool = boolean_type_check('Rubidium-85 fraction fitting', Rb85fracBool)
        if (Rb85frac < 0.0) or (Rb85frac > 100.0):
            print "ERROR: Please set a rubidium-85 percentage from zero to 100."
            exit(1)
    elif (Elem == 'Cs') or (Elem == 'Na'):
        Rb85frac = 0.0
        Rb85fracBool = False
        K40frac = 0.0
        K41frac = 0.0

    elif Elem == 'K':
        Rb85frac = 0.0
        Rb85fracBool = False
        K40frac = value_type_check('potassium-40 fraction', K40frac)
        if (K40frac < 0.0) or (K40frac > 100.0):
            print "ERROR: Please set a potassium-40 percentage from zero to 100."
            exit(1)
        K41frac = value_type_check('potassium-41 fraction', K41frac)
        if (K41frac < 0.0) or (K41frac > 100.0):
            print "ERROR: Please set a potassium-41 percentage from zero to 100."
            exit(1)
        if (K40frac+K41frac) > 100.0:
            print "ERROR: the chosen potassium abundances are larger than 100 %."

    SpectrumType = string_type_check('Spectrum',SpectrumType)
    NTemp = value_type_check('number density temperature', NTemp)
    NTempBool = boolean_type_check('Number density fitting', NTempBool)
    ConstrainBool = boolean_type_check('Temperature constraint', ConstrainBool)
    if ConstrainBool:
        DopplerTemp = NTemp
        DopplerTempBool = False
    else:
        DopplerTemp = value_type_check('doppler temperature', DopplerTemp)
        DopplerTempBool = boolean_type_check('Doppler temperature fitting', DopplerTempBool)
    Bfield = value_type_check('magnetic field', Bfield)
    BfieldBool = boolean_type_check('Magnetic field fitting', BfieldBool)
    if (SpectrumType == 'S1') or (SpectrumType == 'S2') or (SpectrumType == 'Ix') or (SpectrumType == 'Iy'):
        theta0 = value_type_check('linear polarization', theta0)
        theta0Bool = boolean_type_check('Linear polarization fitting', theta0Bool)
    cellLength = value_type_check('cell length', cellLength)
    cellLengthBool = boolean_type_check('Cell length fitting', cellLengthBool)
    if Bfield != 0.0:
        Polar = value_type_check('sigma minus fraction', Polar)
        PolarBool = boolean_type_check('Sigma minus fraction fitting', PolarBool)
    else:
        Polar = 50.0
        PolarBool = False
    shift = value_type_check('shift', shift)
    shiftBool = boolean_type_check('Shift fitting', shiftBool)
    Gamma = value_type_check('extra lorentzian broadening',Gamma)
    GammaBool = boolean_type_check('Lorentzian broadening fitting',GammaBool)
    precision = value_type_check('precision', precision)
    Transition = string_type_check('Transition', Transition)
    DataFilename = string_type_check('Data file name', DataFilename)
    if '.csv' in DataFilename:
        try:
            DataFilename = DataFilename.rstrip('.csv')
        except:
            pass
    SmoothBool = boolean_type_check('smoothing',SmoothBool)
    SmoothFactor = int(value_type_check('Smoothing factor', SmoothFactor))
    FitType = string_type_check('Fitting', FitType)
    if not (FitType in ['T','SA','ML','RR']):
        print "Error: Fitting option not recognised."
        print "Options are: 'T', 'ML', 'RR', 'SA'"
        print "where 'T' is to output a theory plot, and the rest are choices for fitting techniques"
        print "See manual for details"
        exit(1)
    if FitType == 'T':
        DetStart = value_type_check('detuning start point', DetStart)
        DetStop = value_type_check('detuning end point', DetStop)
    PlotBool = boolean_type_check('Plotting', PlotBool)
    PlotHypBool = boolean_type_check('Hyperfine feature plotting',PlotHypBool)
    SavePlotBool = boolean_type_check('Save plot', SavePlotBool)
    if SavePlotBool:
        PlotFileType = string_type_check('Plot file type', PlotFileType)
        if not ('.' in PlotFileType):
            PlotFileType = '.' + PlotFileType

    #Check if numbers make sense and warn user if they look wrong.
    if (Elem == 'Rb' or Elem == 'Cs') and (NTemp < 24.85 or NTemp > 276.85):
        print "Warning: Temperature is outside quoted range of the number density formula."
        print "         The verfied range is 25 to 277 Celsius."
    if (Elem == 'K') and (NTemp < 24.85 or NTemp > 326.85):
        print "Warning: Temperature is outside quoted range of the number density formula."
        print "         The verfied range is 25 to 327 Celsius."
    if (Elem == 'Na') and (NTemp < 24.85 or NTemp > 426.85):
        print "Warning: Temperature is outside quoted range of the number density formula."
        print "         The verfied range is 25 to 427 Celsius."
    if DopplerTemp < -273.149999999:
        print "ERROR: Doppler temperature is below the range that can be handled."
        print "Please set a doppler temperature greater than -273.149999999 Celsius."
        exit(1)
    if Bfield < 0.0:
        print "ERROR: Please set a non-negative magnetic field."
        exit(1)
    if cellLength < 0.0:
        print "ERROR: Please input a non-negative cell length"
        exit(1)
    elif cellLength > 1000.0:
        print "Warning: cell length is very long, (",cellLength/1000.0," metres) is this correct?."
        contin = raw_input('Continue? (y/n): ')
        if contin in ['n','N','no','No','NO']:
            exit(0)
    if (Polar < 0.0) or (Polar > 100.0):
        print "ERROR: Please input a sigma minus percentage from zero to 100."
        exit(1)
    if (abs(shift)) > 10000:
        print 'Warning: Shift is very large, is this correct?'
    if (Gamma < 0.0):
        print "Warning: Extra lorentzian broadening inputed as a negative value."
    if Gamma > 150.0:
        print "Warning: Extra lorentzian broadening is very large, is this correct?"
    if precision < 0.1:
        print "Warning: High precision may use large amounts of memory and calculation time."
        contin = raw_input('Continue? (y/n): ')
        if contin in ['n','N','no','No','NO']:
            exit(0)
    if not (Transition in ['D1','D2']):
        print "ERROR: Transition option not recognised."
        print "Valid options are: 'D1' or 'D2'"

    if not SpectrumType in ['S0','S1','S2','S3','RI+','RI-','GI+','GI-','Ix', 'Iy']:
        print 'ERROR: Spectrum type not recognised.'
        print "Options are: 'S0', 'S1', 'S2', 'S3', 'RI+', 'RI-', 'GI+', 'GI-', 'Ix' and 'Iy'."
        print "See manual for more details."
        exit(1)

    if ((SmoothFactor % 2) == 0):
        print 'Warning: Smoothing factor is not an odd number.'
        print 'This will cause a small shift in the data.'

    if FitType == 'T':
        FitType = 'No fit, just theory'
    elif FitType == 'SA':
        FitType = 'Simulated annealing'
    elif FitType == 'ML':
        FitType = 'Marquardt-Levenberg'
    elif FitType == 'RR':
        FitType = 'Random restart'

    return Elem,NTemp,NTempBool,ConstrainBool,DopplerTemp,DopplerTempBool, \
           Bfield,BfieldBool,Rb85frac,Rb85fracBool,K40frac,K41frac,theta0,theta0Bool, \
           cellLength,cellLengthBool,Polar,PolarBool,shift,shiftBool, \
           Gamma,GammaBool,precision,Transition,SpectrumType, \
           DataFilename,SmoothBool,SmoothFactor,FitType,DetStart, \
           DetStop,PlotBool,SavePlotBool,PlotFileType,PlotHypBool