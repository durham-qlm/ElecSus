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

import spectra as S 

def stickSpectra(X,Elem,Dline):
    '''
    Function to return the hyperfine structure at zero magnetic field
    (stick spectra).
    '''
    StokesType = 'Sticks'
    Bfield = 0.0
    T = 20.0
    lcell = 75.0
    rb85frac = 0.7217
    DoppTemp = 20.0
    theta0 = 0.0
    Pol = 50.0
    shift = 0.0
    GammaBuf = 0.0
    Constrain = True
    precision = 10.0
    K40frac = 0.0117
    K41frac = 6.7302
    
    energies, strengths = S.spectrum(X,Elem,StokesType,Bfield,T,lcell,rb85frac,
                                   DoppTemp,theta0,Pol,shift,GammaBuf,Constrain,
                                   Dline,precision,K40frac,K41frac)
    return energies/1000.0, strengths