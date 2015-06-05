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

import EigenSystem as ES
import AtomConstants as AC

def FreqStren(groundLevels,excitedLevels,groundDim,
              excitedDim,Dline,precision,hand):
    transitionFrequency = zeros(groundDim*2*groundDim) #Initialise lists
    transitionStrength = zeros(groundDim*2*groundDim)
    transNo = 0
    if hand=='Left':
        bottom = 2*groundDim+1
        top = excitedDim+1
    elif hand=='Right':
        bottom = 1
        top = groundDim + 1
    if Dline=='D1':
        interatorList = xrange(groundDim)
    elif Dline=='D2':
        interatorList = xrange(groundDim,excitedDim)
    for gg in xrange(groundDim):
        for ee in interatorList:
            cleb = dot(groundLevels[gg][1:],excitedLevels[ee][bottom:top]).real
            cleb2 = cleb*cleb
            if cleb2 > 0.0005: #If negligable don't calculate.
                transitionFrequency[transNo] = int((-groundLevels[gg][0]
                                              +excitedLevels[ee][0])/precision)
                # We choose to perform the ground manifold reduction (see
                # equation (5) in manual) here for convenience.
                transitionStrength[transNo] = 1./3*1./groundDim*cleb2
                transNo += 1
    return transitionFrequency, transitionStrength, transNo

def add_voigt(
        d,DoppTemp,atomMass,wavenumber,gamma,voigtwidth,
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
    lab=zeros(xpts)
    ldisp=zeros(xpts)
    for line in xrange(ltransno+1):
        xc=lenergy[line]
        lab += lstrength[line]*f_ab(2.0*pi*(d-xc)*precision*1.0e6)
        ldisp += lstrength[line]*f_disp(2.0*pi*(d-xc)*precision*1.0e6)
    rab=zeros(xpts)
    rdisp=zeros(xpts)
    for line in xrange(rtransno+1):
        xc=renergy[line]
        rab += rstrength[line]*f_ab(2.0*pi*(d-xc)*precision*1.0e6)
        rdisp += rstrength[line]*f_disp(2.0*pi*(d-xc)*precision*1.0e6)
    return lab, ldisp, rab, rdisp

def spectrum(
        X,Elem,StokesType,Bfield,T,lcell,rb85frac,DoppTemp,
        theta0,Pol,shift,GammaBuf,Constrain,Dline,precision,
        K40frac,K41frac):
    """Returns a spectrum as a numpy array

    Arguments:
    X          -- the detuning axis in GHz.
    Elem       -- the chosen alkali element.
    StokesType -- specifies the type of spectrum required
    Bfield     -- the magnitude of the magnetic field in Gauss
    T          -- temperature (Celsius) linked to the number density
    lcell      -- length of the vapour cell in millimetres
    rb85frac   -- percentage of rubidium-85 atoms
    DoppTemp   -- temperature (Celsius) linked to the Doppler width
    theta0     -- Linear polarisation angle (in degrees) w.r.t to the x-axis
    Pol        -- percentage of probe beam that is polarised to drive sigma minus
    shift      -- a global frequency shift in MHz
    GammaBuf   -- Extra lorentzian broadening in MHz
    Constrain  -- if True, overides the DoppTemp value and sets it to T
    Dline      -- specifies which D-line transition to calculate for
    precision  -- the required precision of the calculation
    K40frac    -- percentage of potassium-40 atoms
    K41frac    -- percentage of potassium-41 atoms
    """

    print '.', # Print a dot to show the user that the program hasn't crashed

    #Change units to more useful ones
    rb85frac = rb85frac/100.0
    K40frac  = K40frac/100.0
    K41frac  = K41frac/100.0
    Pol      = Pol/100.0
    theta0   = theta0/180.0
    lcell    = lcell/1000.0
    X        = X*1.0e3 # MHz

    # Make sure that the fundamental precision does not affect the result at
    # the specified order. Done by setting precision to an order of 
    # magnetidue better than specified.

    precision = precision*0.1 #Essentially defines the new frequency units.

    if Bfield==0.0:
        Bfield = 0.0001 #To avoid degeneracy problem at B = 0.

    # Rubidium energy levels
    if Elem=='Rb':
        rb87frac=1.0-rb85frac  # Rubidium-87 fraction
        if rb85frac!=0.0: #Save time if no rubidium-85 required
            Rb85atom = AC.Rb85
            #Hamiltonian(isotope,transition,gL,Bfield)
            Rb85_ES = ES.Hamiltonian('Rb85',Dline,1.0,Bfield)

            # Rb-85 allowed transitions for light driving sigma minus
            lenergy85, lstrength85, ltransno85 = FreqStren(
                                                    Rb85_ES.groundManifold,
                                                    Rb85_ES.excitedManifold,
                                                    Rb85_ES.ds,Rb85_ES.dp,
                                                    Dline,precision,'Left')          

            # Rb-85 allowed transitions for light driving sigma plus
            renergy85, rstrength85, rtransno85 = FreqStren(
                                                    Rb85_ES.groundManifold,
                                                    Rb85_ES.excitedManifold,
                                                    Rb85_ES.ds,Rb85_ES.dp,
                                                    Dline,precision,'Right')

        if rb87frac!=0.0:
            Rb87atom = AC.Rb87
            #Hamiltonian(isotope,transition,gL,Bfield)
            Rb87_ES = ES.Hamiltonian('Rb87',Dline,1.0,Bfield)
            # Rb-87 allowed transitions for light driving sigma minus
            lenergy87, lstrength87, ltransno87 = FreqStren(
                                                    Rb87_ES.groundManifold,
                                                    Rb87_ES.excitedManifold,
                                                    Rb87_ES.ds,Rb87_ES.dp,
                                                    Dline,precision,'Left')

            # Rb-87 allowed transitions for light driving sigma plus
            renergy87, rstrength87, rtransno87 = FreqStren(
                                                    Rb87_ES.groundManifold,
                                                    Rb87_ES.excitedManifold,
                                                    Rb87_ES.ds,Rb87_ES.dp,
                                                    Dline,precision,'Right')

        if Dline=='D1':
            transitionConst = AC.RbD1Transition
        elif Dline=='D2':
            transitionConst = AC.RbD2Transition

        if (rb85frac!=0.0) and (rb87frac!=0.0):
            AllEnergyLevels = concatenate((lenergy87,lenergy85,renergy87,
                                           renergy85))
        elif (rb85frac!=0.0) and (rb87frac==0.0):
            AllEnergyLevels = concatenate((lenergy85,renergy85))
        elif (rb85frac==0.0) and (rb87frac!=0.0):
            AllEnergyLevels = concatenate((lenergy87,renergy87))

    # Caesium energy levels
    elif Elem=='Cs':
        CsAtom = AC.Cs
        Cs_ES = ES.Hamiltonian('Cs',Dline,1.0,Bfield)

        lenergy, lstrength, ltransno = FreqStren(Cs_ES.groundManifold,
                                                 Cs_ES.excitedManifold,
                                                 Cs_ES.ds,Cs_ES.dp,Dline,
                                                 precision,'Left')
        renergy, rstrength, rtransno = FreqStren(Cs_ES.groundManifold,
                                                 Cs_ES.excitedManifold,
                                                 Cs_ES.ds,Cs_ES.dp,Dline,
                                                 precision,'Right')

        if Dline=='D1':
            transitionConst = AC.CsD1Transition
        elif Dline =='D2':
            transitionConst = AC.CsD2Transition
        AllEnergyLevels = concatenate((lenergy,renergy))

    # Sodium energy levels
    elif Elem=='Na':
        NaAtom = AC.Na
        Na_ES = ES.Hamiltonian('Na',Dline,1.0,Bfield)

        lenergy, lstrength, ltransno = FreqStren(Na_ES.groundManifold,
                                                 Na_ES.excitedManifold,
                                                 Na_ES.ds,Na_ES.dp,Dline,
                                                 precision,'Left')
        renergy, rstrength, rtransno = FreqStren(Na_ES.groundManifold,
                                                 Na_ES.excitedManifold,
                                                 Na_ES.ds,Na_ES.dp,Dline,
                                                 precision,'Right')

        if Dline=='D1':
            transitionConst = AC.NaD1Transition
        elif Dline=='D2':
            transitionConst = AC.NaD2Transition
        AllEnergyLevels = concatenate((lenergy,renergy))

    #Potassium energy levels
    elif Elem=='K':
        K39frac=1.0-K40frac-K41frac #Potassium-39 fraction
        if K39frac!=0.0:
            K39atom = AC.K39
            K39_ES = ES.Hamiltonian('K39',Dline,1.0,Bfield)
            lenergy39, lstrength39, ltransno39 = FreqStren(
                                                    K39_ES.groundManifold,
                                                    K39_ES.excitedManifold,
                                                    K39_ES.ds,K39_ES.dp,Dline,
                                                    precision,'Left')
            renergy39, rstrength39, rtransno39 = FreqStren(
                                                    K39_ES.groundManifold,
                                                    K39_ES.excitedManifold,
                                                    K39_ES.ds,K39_ES.dp,Dline,
                                                    precision,'Right')

        if K40frac!=0.0:
            K40atom = AC.K40
            K40_ES = ES.Hamiltonian('K40',Dline,1.0,Bfield)
            lenergy40, lstrength40, ltransno40 = FreqStren(
                                                    K40_ES.groundManifold,
                                                    K40_ES.excitedManifold,
                                                    K40_ES.ds,K40_ES.dp,Dline,
                                                    precision,'Left')
            renergy40, rstrength40, rtransno40 = FreqStren(
                                                    K40_ES.groundManifold,
                                                    K40_ES.excitedManifold,
                                                    K40_ES.ds,K40_ES.dp,Dline,
                                                    precision,'Right')
        if K41frac!=0.0:
            K41atom = AC.K41
            K41_ES = ES.Hamiltonian('K41',Dline,1.0,Bfield)
            lenergy41, lstrength41, ltransno41 = FreqStren(
                                                    K41_ES.groundManifold,
                                                    K41_ES.excitedManifold,
                                                    K41_ES.ds,K41_ES.dp,Dline,
                                                    precision,'Left')
            renergy41, rstrength41, rtransno41 = FreqStren(
                                                    K41_ES.groundManifold,
                                                    K41_ES.excitedManifold,
                                                    K41_ES.ds,K41_ES.dp,Dline,
                                                    precision,'Right')
        if Dline=='D1':
            transitionConst = AC.KD1Transition
        elif Dline=='D2':
            transitionConst = AC.KD2Transition
        if K39frac!=0.0:
            AllEnergyLevels = concatenate((lenergy39,renergy39))
        if K40frac!=0.0 and K39frac!=0.0:
            AllEnergyLevels = concatenate((AllEnergyLevels,lenergy40,renergy40))
        elif K40frac!=0.0 and K39frac==0.0:
            AllEnergyLevels = concatenate((lenergy40,renergy40))
        if K41frac!=0.0 and (K39frac!=0.0 or K40frac!=0.0):
            AllEnergyLevels = concatenate((AllEnergyLevels,lenergy41,renergy41))
        elif K41frac!=0.0 and (K39frac==0.0 and K40frac==0.0):
            AllEnergyLevels = concatenate((lenergy41,renergy41))

#Calculate Voigt
    if Constrain:
         DoppTemp = T #Set doppler temperature to the number density temperature
    T+=273.15
    DoppTemp+=273.15

    d = (array(X)-shift)/precision #Linear detuning
    xpts = len(d)
    maxdev = amax(abs(d))

    if Elem=='Rb':
        NDensity=numDenRb(T) #Calculate number density
    elif Elem=='Cs':
        NDensity=numDenCs(T)
    elif Elem=='K':
        NDensity=numDenK(T)
    elif Elem=='Na':
        NDensity=numDenNa(T)

    #Calculate lorentzian broadening and shifts
    gamma0 = 2.0*pi*transitionConst.NatGamma*1.e6
    if Dline=='D1': # D1 self-broadening parameter
        gammaself = 2.0*pi*gamma0*NDensity*\
                    (transitionConst.wavelength/(2.0*pi))**(3)
    elif Dline=='D2': # D2 self-broadening parameter
        gammaself = 2.0*pi*gamma0*NDensity*1.414213562373095*\
                    (transitionConst.wavelength/(2.0*pi))**(3)
    gamma = gamma0 + gammaself
    gamma = gamma + (2.0*pi*GammaBuf*1.e6) #Add extra lorentzian broadening
        
    maxShiftedEnergyLevel = amax(abs(AllEnergyLevels)) #integer value in MHz
    voigtwidth = int(1.1*(maxdev+maxShiftedEnergyLevel))
    wavenumber = transitionConst.wavevectorMagnitude
    dipole = transitionConst.dipoleStrength
    prefactor=2.0*NDensity*dipole*dipole/(hbar*e0)

    if Elem=='Rb':
        if rb85frac!=0.0:
            lab85, ldisp85, rab85, rdisp85 = add_voigt(d,DoppTemp,
                                                       Rb85atom.mass,
                                                       wavenumber,gamma,
                                                       voigtwidth,precision,
                                                       ltransno85,lenergy85,
                                                       lstrength85,rtransno85,
                                                       renergy85,rstrength85)
        if rb87frac!=0.0:
            lab87, ldisp87, rab87, rdisp87 = add_voigt(d,DoppTemp,
                                                       Rb87atom.mass,
                                                       wavenumber,gamma,
                                                       voigtwidth,precision,
                                                       ltransno87,lenergy87,
                                                       lstrength87,rtransno87,
                                                       renergy87,rstrength87)
        # Make the parts of the susceptibility
        if (rb85frac!=0.0) and (rb87frac!=0.0):
            ChiRealLeft= prefactor*(rb85frac*ldisp85+rb87frac*ldisp87)
            ChiRealRight= prefactor*(rb85frac*rdisp85+rb87frac*rdisp87)
            ChiImLeft = prefactor*(rb85frac*lab85+rb87frac*lab87)
            ChiImRight = prefactor*(rb85frac*rab85+rb87frac*rab87)
        elif (rb85frac!=0.0) and (rb87frac==0.0):
            ChiRealLeft= prefactor*rb85frac*ldisp85
            ChiRealRight= prefactor*rb85frac*rdisp85
            ChiImLeft = prefactor*rb85frac*lab85
            ChiImRight = prefactor*rb85frac*rab85
        elif (rb85frac==0.0) and (rb87frac!=0.0):
            ChiRealLeft= prefactor*rb87frac*ldisp87
            ChiRealRight= prefactor*rb87frac*rdisp87
            ChiImLeft = prefactor*rb87frac*lab87
            ChiImRight = prefactor*rb87frac*rab87
    elif Elem=='Cs':
        lab, ldisp, rab, rdisp = add_voigt(d,DoppTemp,CsAtom.mass,wavenumber,
                                           gamma,voigtwidth,precision,ltransno,
                                           lenergy,lstrength,rtransno,renergy,
                                           rstrength)
        ChiRealLeft= prefactor*ldisp
        ChiRealRight= prefactor*rdisp
        ChiImLeft = prefactor*lab
        ChiImRight = prefactor*rab
    elif Elem=='Na':
        lab, ldisp, rab, rdisp = add_voigt(d,DoppTemp,NaAtom.mass,wavenumber,
                                           gamma,voigtwidth,precision,ltransno,
                                           lenergy,lstrength,rtransno,renergy,
                                           rstrength)
        ChiRealLeft= prefactor*ldisp
        ChiRealRight= prefactor*rdisp
        ChiImLeft = prefactor*lab
        ChiImRight = prefactor*rab
    elif Elem=='K':
        if K39frac!=0.0:
            lab39, ldisp39, rab39, rdisp39 = add_voigt(d,DoppTemp,K39atom.mass,
                                                       wavenumber,gamma,
                                                       voigtwidth,precision,
                                                       ltransno39,lenergy39,
                                                       lstrength39,rtransno39,
                                                       renergy39,rstrength39)
        if K40frac!=0.0:
            lab40, ldisp40, rab40, rdisp40 = add_voigt(d,DoppTemp,K40atom.mass,
                                                       wavenumber,gamma,
                                                       voigtwidth,precision,
                                                       ltransno40,lenergy40,
                                                       lstrength40,rtransno40,
                                                       renergy40,rstrength40)
        if K41frac!=0.0:
            lab41, ldisp41, rab41, rdisp41 = add_voigt(d,DoppTemp,K41atom.mass,
                                                       wavenumber,gamma,
                                                       voigtwidth,precision,
                                                       ltransno41,lenergy41,
                                                       lstrength41,rtransno41,
                                                       renergy41,rstrength41)

        if (K39frac!=0.0) and (K40frac!=0.0) and (K41frac!=0.0):
            ChiRealLeft = prefactor*(K39frac*ldisp39+K40frac\
                                    *ldisp40+K41frac*ldisp41)
            ChiRealRight = prefactor*(K39frac*rdisp39+K40frac\
                                    *rdisp40+K41frac*rdisp41)
            ChiImLeft = prefactor*(K39frac*lab39+K40frac*lab40+K41frac*lab41)
            ChiImRight = prefactor*(K39frac*rab39+K40frac*rab40+K41frac*rab41)
        elif (K39frac!=0.0) and (K40frac!=0.0) and (K41frac==0.0):
            ChiRealLeft= prefactor*(K39frac*ldisp39+K40frac*ldisp40)
            ChiRealRight= prefactor*(K39frac*rdisp39+K40frac*rdisp40)
            ChiImLeft = prefactor*(K39frac*lab39+K40frac*lab40)
            ChiImRight = prefactor*(K39frac*rab39+K40frac*rab40)
        elif (K39frac!=0.0) and (K40frac==0.0) and (K41frac!=0.0):
            ChiRealLeft= prefactor*(K39frac*ldisp39+K41frac*ldisp41)
            ChiRealRight= prefactor*(K39frac*rdisp39+K41frac*rdisp41)
            ChiImLeft = prefactor*(K39frac*lab39+K41frac*lab41)
            ChiImRight = prefactor*(K39frac*rab39+K41frac*rab41)
        elif (K39frac==0.0) and (K40frac!=0.0) and (K41frac!=0.0):
            ChiRealLeft= prefactor*(K40frac*ldisp40+K41frac*ldisp41)
            ChiRealRight= prefactor*(K40frac*rdisp40+K41frac*rdisp41)
            ChiImLeft = prefactor*(K40frac*lab40+K41frac*lab41)
            ChiImRight = prefactor*(K40frac*rab40+K41frac*rab41)
        elif (K39frac!=0.0) and (K40frac==0.0) and (K41frac==0.0):
            ChiRealLeft= prefactor*ldisp39
            ChiRealRight= prefactor*rdisp39
            ChiImLeft = prefactor*lab39
            ChiImRight = prefactor*rab39
        elif (K39frac==0.0) and (K40frac!=0.0) and (K41frac==0.0):
            ChiRealLeft= prefactor*ldisp40
            ChiRealRight= prefactor*rdisp40
            ChiImLeft = prefactor*lab40
            ChiImRight = prefactor*rab40
        elif (K39frac==0.0) and (K40frac==0.0) and (K41frac!=0.0):
            ChiRealLeft= prefactor*ldisp41
            ChiRealRight= prefactor*rdisp41
            ChiImLeft = prefactor*lab41
            ChiImRight = prefactor*rab41

    # Reconstruct total susceptibility and index of refraction
    totalChiLeft = ChiRealLeft + 1.j*ChiImLeft
    totalChiRight = ChiRealRight + 1.j*ChiImRight
    
    total_n_Left = sqrt(1.0+totalChiLeft) #Complex refractive index left hand
    total_n_Right = sqrt(1.0+totalChiRight) #Complex refractive index right hand

    # Return user specified spectrum
    if StokesType=='S0':
        TransLeft = exp(-2.0*total_n_Left.imag*wavenumber*lcell)
        TransRight = exp(-2.0*total_n_Right.imag*wavenumber*lcell)
        return (Pol*TransLeft) + ((1.0-Pol)*TransRight)
    elif StokesType=='S1':
        THETA0 = theta0*pi
        phiPlus = wavenumber*total_n_Left.real*lcell
        phiMinus = wavenumber*total_n_Right.real*lcell
        PHI = (phiMinus-phiPlus)/2.0 + THETA0
        PolFactor = 2.0*sqrt(Pol-(Pol**2))
        return PolFactor*cos(2.0*PHI)*exp(-wavenumber*lcell*(total_n_Right.imag
                                                           +total_n_Left.imag))
    elif StokesType=='S2':
        THETA0 = theta0*pi
        phiPlus = wavenumber*total_n_Left.real*lcell
        phiMinus = wavenumber*total_n_Right.real*lcell
        PHI = (phiMinus-phiPlus)/2.0 + THETA0
        PolFactor = 2.0*sqrt(Pol-(Pol**2))
        return PolFactor*sin(2.0*PHI)*exp(-wavenumber*lcell*(total_n_Right.imag
                                                           +total_n_Left.imag))
    elif StokesType=='S3':
        TransLeft = exp(-2.0*total_n_Left.imag*wavenumber*lcell)
        TransRight = exp(-2.0*total_n_Right.imag*wavenumber*lcell)
        return (Pol*TransLeft) - ((1.0-Pol)*TransRight)
    elif StokesType=='Ix':
        TransLeft = exp(-2.0*total_n_Left.imag*wavenumber*lcell) 
        TransRight = exp(-2.0*total_n_Right.imag*wavenumber*lcell) 
        s0 = (Pol*TransLeft) + ((1.0-Pol)*TransRight)
        THETA0 = theta0*pi
        phiPlus = wavenumber*total_n_Left.real*lcell
        phiMinus = wavenumber*total_n_Right.real*lcell
        PHI = (phiMinus-phiPlus)/2.0 + THETA0
        PolFactor = 2.0*sqrt(Pol-(Pol**2))
        s1 = PolFactor*cos(2.0*PHI)*exp(-wavenumber*lcell*(total_n_Right.imag
                                                          +total_n_Left.imag))
        return (s1+s0)/2.0
    elif StokesType=='Iy':
        TransLeft = exp(-2.0*total_n_Left.imag*wavenumber*lcell)
        TransRight = exp(-2.0*total_n_Right.imag*wavenumber*lcell)
        s0 = (Pol*TransLeft) + ((1.0-Pol)*TransRight)
        THETA0 = theta0*pi
        phiPlus = wavenumber*total_n_Left.real*lcell
        phiMinus = wavenumber*total_n_Right.real*lcell
        PHI = (phiMinus-phiPlus)/2.0 + THETA0
        PolFactor = 2.0*sqrt(Pol-(Pol**2))
        s1 = PolFactor*cos(2.0*PHI)*exp(-wavenumber*lcell*(total_n_Right.imag
                                                          +total_n_Left.imag))
        return (s0-s1)/2.0
    elif StokesType=='RI-':
        return total_n_Left.real
    elif StokesType=='RI+':
        return total_n_Right.real
    elif StokesType=='GI-':
        dnWRTv = derivative(d,total_n_Left.real)
        GroupIndex = total_n_Left.real + (X + transitionConst.v0*1.0e-6)*dnWRTv
        return GroupIndex
    elif StokesType=='GI+':
        dnWRTv = derivative(d,total_n_Right.real)
        GroupIndex = total_n_Right.real + (X + transitionConst.v0*1.0e-6)*dnWRTv
        return GroupIndex
    elif StokesType=='Sticks':
        if (Elem=='Cs' or Elem=='Na'):
            return lenergy, lstrength
        elif Elem=='Rb':
            return concatenate((lenergy87,lenergy85)), concatenate(
                   (0.2783*lstrength87,0.7217*lstrength85))
        elif Elem=='K':
            return concatenate((lenergy39,lenergy40,lenergy41)),concatenate((
                   0.932581*lstrength39,0.000117*lstrength40,0.067302*lstrength41))
