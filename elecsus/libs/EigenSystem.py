# Copyright 2014 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams and I. G. Hughes.

# Updated 2017 JK

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
Calculates the atomic Hamiltonian for a given Isotope and magnetic field

Modules called:

FundamentalConstants -- fundamental physical constants from CODATA
AtomConstants        -- All isotope and D-line specific constants
sz_lsi               -- 

Last updated 2018-02-19 JK
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)



# Calulates the ground state manifold and the excited state manifold.
# Is called by spectra.py

from scipy.linalg import eig, eigh
from numpy import pi, append, transpose, identity

from .AtomConstants import *
from .FundamentalConstants import *
from .sz_lsi import sz, lz, Iz
from .fs_hfs import Hfs,Hhfs,Bbhfs


class Hamiltonian(object):
    """Functions to create the atomic hamiltonian."""

    def __init__(self, Isotope, Trans, gL, Bfield):
        """Ground and excited state Hamiltonian for an isotope"""
        if Isotope=='Rb87':
            atom = Rb87
        elif Isotope=='Rb85':
            atom = Rb85
        elif Isotope=='Cs':
            atom = Cs
        elif Isotope=='K39':
            atom = K39
        elif Isotope=='K40':
            atom = K40
        elif Isotope=='K41':
            atom = K41
        elif Isotope=='Na':
            atom = Na
        elif Isotope=='IdealAtom':
            atom = IdealAtom
            transition = IdealD1Transition
            atom_transition = Ideal_D1

        self.atom = atom
		
        if (Trans=='D1') and (Isotope=='Rb85'):
            transition = RbD1Transition
            atom_transition = Rb85_D1
        elif (Trans=='D2') and (Isotope=='Rb85'):
            transition = RbD2Transition
            atom_transition = Rb85_D2
        elif (Trans=='D1') and (Isotope=='Rb87'):
            transition = RbD1Transition
            atom_transition = Rb87_D1
        elif (Trans=='D2') and (Isotope=='Rb87'):
            transition = RbD2Transition
            atom_transition = Rb87_D2
        elif (Trans=='D1') and (Isotope=='Cs'):
            transition = CsD1Transition
            atom_transition = Cs_D1
        elif (Trans=='D2') and (Isotope=='Cs'):
            transition = CsD2Transition
            atom_transition = Cs_D2
        elif (Trans=='D1') and (Isotope=='Na'):
            transition = NaD1Transition
            atom_transition = Na_D1
        elif (Trans=='D2') and (Isotope=='Na'):
            transition = NaD2Transition
            atom_transition = Na_D2
        elif (Trans=='D1') and (Isotope=='K39'):
            transition = KD1Transition
            atom_transition = K39_D1
        elif (Trans=='D2') and (Isotope=='K39'):
            transition = KD2Transition
            atom_transition = K39_D2
        elif (Trans=='D1') and (Isotope=='K40'):
            transition = KD1Transition
            atom_transition = K40_D1
        elif (Trans=='D2') and (Isotope=='K40'):
            transition = KD2Transition
            atom_transition = K40_D2
        elif (Trans=='D1') and (Isotope=='K41'):
            transition = KD1Transition
            atom_transition = K41_D1
        elif (Trans=='D2') and (Isotope=='K41'):
            transition = KD2Transition
            atom_transition = K41_D2
			
        if Bfield == 0.0:
            Bfield += 1e-5 # avoid degeneracy problem..?

        #Useful quantities to return
        self.ds=int((2*S+1)*(2*atom.I+1)) #Dimension of S-term matrix
        self.dp=int(3*(2*S+1)*(2*atom.I+1)) #Dimension of P-term matrix

        self.groundManifold, self.groundEnergies = self.groundStateManifold(atom.gI,atom.I,atom.As,
                                atom_transition.IsotopeShift,Bfield)
        self.excitedManifold, self.excitedEnergies = self.excitedStateManifold(gL,atom.gI,atom.I,
                                atom_transition.Ap,atom_transition.Bp,Bfield)
    
    def groundStateManifold(self,gI,I,A_hyp_coeff,IsotopeShift,Bfield):
        """Function to produce the ground state manifold"""
        ds = int((2*S+1)*(2*I+1))  # total dimension of matrix
        #print 'Matrix dim:', ds
        As = A_hyp_coeff
        # Add the S-term hyperfine interaction
        S_StateHamiltonian = As*Hhfs(0.0,S,I)+IsotopeShift*identity(ds)
        Ez = muB*Bfield*1.e-4/(hbar*2.0*pi*1.0e6)
        S_StateHamiltonian += Ez*(gs*sz(0.0,S,I)+gI*Iz(0.0,S,I)) # Add Zeeman
        EigenSystem = eigh(S_StateHamiltonian)
        EigenValues = EigenSystem[0].real
        EigenVectors = EigenSystem[1]
        stateManifold = append([EigenValues],EigenVectors,axis=0)
        sortedManifold = sorted(transpose(stateManifold),key=(lambda i:i[0]))
        return sortedManifold, EigenValues

    def excitedStateManifold(self,gL,gI,I,A_hyp_coeff,B_hyp_coeff,Bfield):
        """Function to produce the excited state manifold"""
        dp = int(3*(2*S+1)*(2*I+1))  # total dimension of matrix
        # The actual value of FS is unimportant.
        FS = self.atom.FS # Fine structure splitting
        Ap = A_hyp_coeff
        Bp = B_hyp_coeff
        # Add P-term fine and hyperfine interactions
        if Bp==0.0:
            P_StateHamiltonian=FS*Hfs(1.0,S,I)+FS*identity(dp)+Ap*Hhfs(1.0,S,I)
        if Bp!=0.0:
            P_StateHamiltonian=FS*Hfs(1.0,S,I)-(FS/2.0)*identity(dp)+Ap*Hhfs(1.0,S,I)
            P_StateHamiltonian+=Bp*Bbhfs(1.0,S,I) # add p state quadrupole
        E=muB*(Bfield*1.0e-4)/(hbar*2.0*pi*1.0e6)
        # Add magnetic interaction
        P_StateHamiltonian+=E*(gL*lz(1.0,S,I)+gs*sz(1.0,S,I)+gI*Iz(1.0,S,I))
        ep=eigh(P_StateHamiltonian)
        EigenValues=ep[0].real
        EigenVectors=ep[1]
        stateManifold=append([EigenValues],EigenVectors,axis=0)
        sortedManifold=sorted(transpose(stateManifold),key=(lambda i:i[0]))
        return sortedManifold, EigenValues
