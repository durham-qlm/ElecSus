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


#######
# modified to allow returning electric field and raw susceptibilities
# JK - 1 Nov. 2014

# modified for Am J Phys FADOF paper, to model an 'ideal' atom
# JK - 26 Nov 2014

# modified for Sr(87) 1S0 to 3P1 J=0 to 1 transition - Sr FADOF Clock

#######

"""Module containing functions to calculate the spectra

Constructs the electric susceptibility and then returns
the spectrum requested

Calls numberDensityEqs, tools, EigenSystem and AtomConstants
"""

from numpy import zeros,sqrt,pi,dot,exp,sin,cos,array,amax,arange,concatenate
from scipy.special import wofz
from scipy.interpolate import interp1d
from FundamentalConstants import *
from numberDensityEqs import *
from tools import derivative

import EigenSystem_Sr88 as ES

class Sr88Atom:
	""" Constants for a J=0 -> J=1 atom """
	I = 0
	As = 0
	gI = 0
	mass = 87*amu

class SrNarrowTransition:
    """Constants relating to the Sr88 intercombination line transition"""
    wavelength = 689.45e-9 #The weighted linecentre of the rubidium D1 line in m
    wavevectorMagnitude = 2.0*pi/wavelength #Magnitude of the wavevector
    NatGamma = 7.6e-3 # Natural linewidth in MHz
    dipoleStrength = 3.0*sqrt(e0*hbar*(2.0*NatGamma*(10.0**6))*(wavelength**3)/(8.0*pi))
	#print dipoleStrength / e / a0
    v0 = 435.13e12 #The weighted linecentre of the rubidium D1 line in Hz
    Ap = 0
    Bp = 0
    IsotopeShift = 0. # 21.734 #MHz
	
def FreqStren_noLS(groundLevels,excitedLevels,groundDim,
              excitedDim,precision,hand):
    transitionFrequency = zeros(groundDim*2*groundDim) #Initialise lists
    transitionStrength = zeros(groundDim*2*groundDim)
    transNo = 0
    if hand=='Left':
        bottom = 2*groundDim+1
        top = excitedDim+1
    elif hand=='Right':
        bottom = 1
        top = groundDim + 1
    interatorList = xrange(3)
    for gg in xrange(groundDim):
        for ee in interatorList:
            cleb = dot(groundLevels[gg][1:],excitedLevels[ee][bottom:top]).real
            cleb2 = cleb*cleb
            if cleb2 > 0.0005: #If negligible don't calculate.
                transitionFrequency[transNo] = int((-groundLevels[gg][0].real
                                              +excitedLevels[ee][0].real)/precision)
                # We choose to perform the ground manifold reduction (see
                # equation (4) in manual) here for convenience.
                transitionStrength[transNo] = 1./3*1./groundDim*cleb2        ######## the factor of 1/3 only for L=0 to 1?
                transNo += 1
    return transitionFrequency, transitionStrength, transNo

def add_voigt(d,DoppTemp,atomMass,wavenumber,gamma,voigtwidth,
        precision,ltransno,lenergy,lstrength,rtransno,
        renergy,rstrength):
    xpts = len(d)
    npts = 2*voigtwidth+1
    detune = 2.0*pi*1.0e6*(arange(npts)-voigtwidth)*precision #Angular detuning
    wavenumber =  wavenumber + detune/c #Allow the wavenumber to change
    u = sqrt(2.0*kB*DoppTemp/atomMass)
    ku = wavenumber*u
    a = gamma/ku
    b = detune/ku
    y = 1.0j*(0.5*sqrt(pi)/ku)*wofz(b+0.5j*a)
    ab = y.imag
    disp = y.real
    #interpolate lineshape functions
    f_ab = interp1d(detune,ab)
    f_disp = interp1d(detune,disp)
    #Add contributions from all transitions to user defined detuning axis
    lab = zeros(xpts)
    ldisp = zeros(xpts)
    for line in xrange(ltransno+1):
        xc = lenergy[line]
        lab += lstrength[line]*f_ab(2.0*pi*(d-xc)*precision*1.0e6)
        ldisp += lstrength[line]*f_disp(2.0*pi*(d-xc)*precision*1.0e6)
    rab = zeros(xpts)
    rdisp = zeros(xpts)
    for line in xrange(rtransno+1):
        xc = renergy[line]
        rab += rstrength[line]*f_ab(2.0*pi*(d-xc)*precision*1.0e6)
        rdisp += rstrength[line]*f_disp(2.0*pi*(d-xc)*precision*1.0e6)
    return lab, ldisp, rab, rdisp

	
##### main method #####
def spectrum(
        X,
		Bfield=0,T=20,lcell=75, #experimental conditions
		shift=0,GammaBuf=0,
		
		theta0=0,Pol=50, #default linear pol
		PolVector=None,

		precision=10
		):
				
    """Returns a spectrum as a numpy array

    Arguments:
    X          -- the detuning axis in GHz (numpy 1d array, list or float).
    
	Options:

	Elem       -- the chosen alkali element.
    Dline      -- specifies which D-line transition to calculate for
	OutputType -- specifies the type of spectrum required

	### Experimental parameters
    Bfield     -- the magnitude of the magnetic field in Gauss
    T          -- temperature (Celsius) linked to the number density
    lcell      -- length of the vapour cell in millimetres
    GammaBuf   -- Extra lorentzian broadening in MHz (usually from buffer gas 
					but can be any extra homogeneous broadening
    shift      -- a global frequency shift in MHz
	
	### Initial polarisation
    theta0     -- Linear polarisation angle (in degrees) w.r.t to the x-axis
    Pol        -- percentage of probe beam that is polarised to drive sigma minus
	PolVector  -- Polarisation in the Jones vector circular basis (array of size (len(X),2))

    DoppTemp   -- temperature (Celsius) linked to the Doppler width
    Constrain  -- if True, overides the DoppTemp value and sets it to T
	SkipVoigt  -- If True, returns the Lorentzian absorption, i.e. without Doppler broadening
					for use with nano-cells where the normal Maxwellian isn't valid
					Overrides DoppTemp and Constrain options

	### Elemental abundancies, where applicable
	Atomfrac   -- percentage of rubidium-85 atoms
    K40frac    -- percentage of potassium-40 atoms
    K41frac    -- percentage of potassium-41 atoms

    precision  -- the required precision of the calculation
    """

	#testing
    #print Bfield, T,
	
    ##convert X to array
    X = array(X)
	
    #print '.' # Print a dot to show the user that the program hasn't crashed (for fitting 
               # routines where this method is called many times)

    #Change units to more useful ones
    Pol      = Pol/100.0
    theta0   = theta0/180.0
    lcell    = lcell/1000.0 # mm to m
    X        = X*1.0e3 # GHz to MHz

    # Make sure that the fundamental precision does not affect the result at
    # the specified order. Done by setting precision to an order of 
    # magnetidue better than specified.

    precision = precision*0.0001 #Essentially defines the new frequency units.

    if Bfield==0.0:
        Bfield = 0.0001 #To avoid degeneracy problem at B = 0.

    # Rubidium energy levels
    Atom = Sr88Atom()
    Transition = SrNarrowTransition()
    Atom_ES = ES.Hamiltonian(Atom, Transition, 1.0, Bfield)

    # Rb-85 allowed transitions for light driving sigma minus
    lenergy, lstrength, ltransno = FreqStren_noLS(
                                            Atom_ES.groundManifold,
                                            Atom_ES.excitedManifold,
                                            Atom_ES.ds,Atom_ES.dp,
                                            precision,'Left')          

    # Rb-85 allowed transitions for light driving sigma plus
    renergy, rstrength, rtransno = FreqStren_noLS(
                                            Atom_ES.groundManifold,
                                            Atom_ES.excitedManifold,
                                            Atom_ES.ds,Atom_ES.dp,
                                            precision,'Right')

    AllEnergyLevels = concatenate((lenergy,renergy))
    transitionConst = SrNarrowTransition()

#Calculate Voigt
    DoppTemp = 200e-6 # Cold atomic beam - Doppler suppressed; T sets number density only

    d = (array(X)-shift)/precision #Linear detuning
    xpts = len(d)
    maxdev = amax(abs(d))

    NDensity = T*1e6 # convert number density from cm-3 to m-3

    #Calculate lorentzian broadening and shifts
    gamma0 = 2.0*pi*transitionConst.NatGamma*1.e6
    gammaself = 2.0*pi*gamma0*NDensity*sqrt(2)*\
                    (transitionConst.wavelength/(2.0*pi))**(3)
    gamma = gamma0 # + gammaself ## no self-broadening
    gamma = gamma + (2.0*pi*GammaBuf*1.e6) #Add extra lorentzian broadening
        
    maxShiftedEnergyLevel = amax(abs(AllEnergyLevels)) #integer value in MHz
    voigtwidth = int(1.1*(maxdev+maxShiftedEnergyLevel))
    wavenumber = transitionConst.wavevectorMagnitude
    dipole = transitionConst.dipoleStrength
    prefactor = 2.0*NDensity*dipole**2/(hbar*e0)

##calc absorption / dispersion
    lab, ldisp, rab, rdisp = add_voigt(d,DoppTemp,Atom.mass,wavenumber,
                                       gamma,voigtwidth,precision,ltransno,
                                       lenergy,lstrength,rtransno,renergy,
                                       rstrength)
    ChiRealLeft = prefactor*ldisp
    ChiRealRight = prefactor*rdisp
    ChiImLeft = prefactor*lab
    ChiImRight = prefactor*rab

# Reconstruct total susceptibility and index of refraction
    totalChiLeft = ChiRealLeft + 1.j*ChiImLeft
    totalChiRight = ChiRealRight + 1.j*ChiImRight
    
    total_n_Left = sqrt(1.0+totalChiLeft) #Complex refractive index left hand
    total_n_Right = sqrt(1.0+totalChiRight) #Complex refractive index right hand

# calc outputs
    TransLeft = exp(-2.0*total_n_Left.imag*wavenumber*lcell)
    TransRight = exp(-2.0*total_n_Right.imag*wavenumber*lcell)
    THETA0 = theta0*pi
    phiPlus = wavenumber*total_n_Left.real*lcell
    phiMinus = wavenumber*total_n_Right.real*lcell
    PHI = (phiMinus-phiPlus)/2.0 + THETA0
    PolFactor = 2.0*sqrt(Pol-(Pol**2))

    S0 = (Pol*TransLeft) + ((1.0-Pol)*TransRight)
    S1 = PolFactor*cos(2.0*PHI)*exp(-wavenumber*lcell*(total_n_Right.imag
                                                           +total_n_Left.imag))
    S2 = PolFactor*sin(2.0*PHI)*exp(-wavenumber*lcell*(total_n_Right.imag
                                                           +total_n_Left.imag))
    S3 = (Pol*TransLeft) - ((1.0-Pol)*TransRight)
    Ix = (S1+S0)/2.0
    Iy = (S0-S1)/2.0
    nminus = total_n_Left.real
    nplus = total_n_Right.real
	
    alphaplus = 2.0*total_n_Right.imag*wavenumber
    alphaminus = 2.0*total_n_Left.imag*wavenumber
    
    return S0,S1,S2,S3,Ix,Iy,nminus,nplus,PHI,alphaplus,alphaminus