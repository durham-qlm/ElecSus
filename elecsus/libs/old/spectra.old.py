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

"""Module containing functions to calculate the spectra

Constructs the electric susceptibility and then returns
the spectrum requested

Calls numberDensityEqs, tools, EigenSystem and AtomConstants

Updated 2015-12-15 JK
"""

from numpy import zeros,sqrt,pi,dot,exp,sin,cos,array,amax,arange,concatenate
import numpy as np
from scipy.special import wofz
from scipy.interpolate import interp1d
from FundamentalConstants import *
from numberDensityEqs import *
from tools import derivative

import EigenSystem as ES
import AtomConstants as AC
import RotationMatrices as RM
import BasisChanger as BC
import JonesMatrices as JM

#testing
import matplotlib.pyplot as plt

# Default values for parameters
p_dict_defaults = {	'Elem':'Rb', 'Dline':'D2', 
							'lcell':75e-3,'Bfield':0., 'T':20., 
							'GammaBuf':0., 'shift':0.,
							# Polarisation of light
							'theta0':0., 'Pol':50., 
							# B-field angle w.r.t. light k-vector
							'Btheta':0, 'Bphi':0,
							'Constrain':True, 'DoppTemp':20.,
							'rb85frac':72.17, 'K40frac':0.01, 'K41frac':6.73 }

def FreqStren(groundLevels,excitedLevels,groundDim,
			  excitedDim,Dline,hand):
	""" Calculate transition frequencies and strengths by taking dot products of relevant parts of the ground / excited state eigenvectors """
	
	transitionFrequency = zeros(groundDim*2*groundDim) #Initialise lists
	transitionStrength = zeros(groundDim*2*groundDim)
	transNo = 0
	
	# select correct columns of matrix corresponding to delta mL = -1, 0, +1
	if hand=='Right':
		bottom = 1
		top = groundDim + 1
	elif hand=='Z':
		bottom = groundDim + 1
		top = 2*groundDim + 1
	elif hand=='Left':
		bottom = 2*groundDim+1
		top = excitedDim+1

	# select correct rows of matrix corresponding to D1 / D2 lines
	if Dline=='D1':
		interatorList = xrange(groundDim)
	elif Dline=='D2':
		interatorList = xrange(groundDim,excitedDim)\
	
	# find difference in energies and do dot product between (all) ground states
	# 	and selected parts of excited state matrix
	for gg in xrange(groundDim):
		for ee in interatorList:
			cleb = dot(groundLevels[gg][1:],excitedLevels[ee][bottom:top]).real
			cleb2 = cleb*cleb
			if cleb2 > 0.0005: #If negligable don't calculate.
				transitionFrequency[transNo] = int((-groundLevels[gg][0].real
											  +excitedLevels[ee][0].real))
				# We choose to perform the ground manifold reduction (see
				# equation (4) in manual) here for convenience.
				transitionStrength[transNo] = 1./3*1./groundDim*cleb2
				transNo += 1

	#print 'No transitions (ElecSus): ',transNo
	return transitionFrequency, transitionStrength, transNo

def add_voigt(d,DoppTemp,atomMass,wavenumber,gamma,voigtwidth,
		ltransno,lenergy,lstrength,
		rtransno,renergy,rstrength,
		ztransno,zenergy,zstrength):
	xpts = len(d)
	npts = 2*voigtwidth+1
	detune = 2.0*pi*1.0e6*(arange(npts)-voigtwidth) # Angular detuning (2pi Hz)
	wavenumber =  wavenumber + detune/c # Allow the wavenumber to change (large detuning)
	u = sqrt(2.0*kB*DoppTemp/atomMass)
	ku = wavenumber*u
	
	# Fadeeva function:
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
		lab += lstrength[line]*f_ab(2.0*pi*(d-xc)*1.0e6)
		ldisp += lstrength[line]*f_disp(2.0*pi*(d-xc)*1.0e6)
	rab = zeros(xpts)
	rdisp = zeros(xpts)
	for line in xrange(rtransno+1):
		xc = renergy[line]
		rab += rstrength[line]*f_ab(2.0*pi*(d-xc)*1.0e6)
		rdisp += rstrength[line]*f_disp(2.0*pi*(d-xc)*1.0e6)
	zab = zeros(xpts)
	zdisp = zeros(xpts)
	for line in xrange(ztransno+1):
		xc = zenergy[line]
		zab += zstrength[line]*f_ab(2.0*pi*(d-xc)*1.0e6)
		zdisp += zstrength[line]*f_disp(2.0*pi*(d-xc)*1.0e6)
	return lab, ldisp, rab, rdisp, zab, zdisp

	
##### main method #####
def calc_chi(X, p_dict):			   
	"""Returns the complex susceptibility for left and right-hand circularly polarised light as a 1D array

	Arguments:
	
		X: 		Detuning axis (float, list, or numpy array) in MHz
		p_dict: 	Dictionary containing all parameters (the order of parameters is therefore not important)
				
			Dictionary keys:

			Key				DataType	Unit		Description
			---				---------	----		-----------
			Elem	   			str			--			The chosen alkali element.
			Dline	  			str			--			Specifies which D-line transition to calculate for (D1 or D2)
			
			# Experimental parameters
			Bfield	 			float			Gauss	Magnitude of the applied magnetic field
			T		 			float			Celsius	Temperature used to calculate atomic number density
			GammaBuf   	float			MHz		Extra lorentzian broadening (usually from buffer gas 
															but can be any extra homogeneous broadening)
			shift	  			float			MHz		A global frequency shift of the atomic resonance frequencies

			DoppTemp   	float			Celsius	Temperature linked to the Doppler width (used for
															independent Doppler width and number density)
			Constrain  		bool			--			If True, overides the DoppTemp value and sets it to T

			# Elemental abundancies, where applicable
			rb85frac   		float			%			percentage of rubidium-85 atoms
			K40frac			float			%			percentage of potassium-40 atoms
			K41frac			float			%			percentage of potassium-41 atoms
			
			lcell	  			float			m			length of the vapour cell
			theta0	 		float			degrees	Linear polarisation angle w.r.t. to the x-axis
			Pol				float			%			Percentage of probe beam that drives sigma minus (50% = linear polarisation)
			
			NOTE: If keys are missing from p_dict, default values contained in p_dict_defaults will be loaded.
			
			Any additional keys in the dict are ignored.

	"""
	
	# get parameters from dictionary
	if 'Elem' in p_dict.keys():
		Elem = p_dict['Elem']
	else:
		Elem = p_dict_defaults['Elem']
	if 'Dline' in p_dict.keys():
		Dline = p_dict['Dline']
	else:
		Dline = p_dict_defaults['Dline']
	if 'T' in p_dict.keys():
		T = p_dict['T']
	else:
		T = p_dict_defaults['T']
	if 'Bfield' in p_dict.keys():
		Bfield = p_dict['Bfield']
	else:
		Bfield = p_dict_defaults['Bfield']
	if 'GammaBuf' in p_dict.keys():
		GammaBuf = p_dict['GammaBuf']
	else:
		GammaBuf = p_dict_defaults['GammaBuf']
	if 'shift' in p_dict.keys():
		shift = p_dict['shift']
	else:
		shift = p_dict_defaults['shift']
	if 'Constrain' in p_dict.keys():
		Constrain = p_dict['Constrain']
	else:
		Constrain = p_dict_defaults['Constrain']
	if 'rb85frac' in p_dict.keys():
		rb85frac = p_dict['rb85frac']
	else:
		rb85frac = p_dict_defaults['rb85frac']
	if 'DoppTemp' in p_dict.keys():
		DoppTemp = p_dict['DoppTemp']
	else:
		DoppTemp = p_dict_defaults['DoppTemp']
	if 'K40frac' in p_dict.keys():
		K40frac = p_dict['K40frac']
	else:
		K40frac = p_dict_defaults['K40frac']
	if 'K41frac' in p_dict.keys():
		K41frac = p_dict['K41frac']
	else:
		K41frac = p_dict_defaults['K41frac']
		
	# convert X to array if needed (does nothing otherwise)
	X = array(X)
   
	# Change to fraction from %
	rb85frac = rb85frac/100.0
	K40frac  = K40frac/100.0
	K41frac  = K41frac/100.0

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
													Dline,'Left')		  

			# Rb-85 allowed transitions for light driving sigma plus
			renergy85, rstrength85, rtransno85 = FreqStren(
													Rb85_ES.groundManifold,
													Rb85_ES.excitedManifold,
													Rb85_ES.ds,Rb85_ES.dp,
													Dline,'Right')
			
			# Rb-85 allowed transitions for light driving sigma plus
			zenergy85, zstrength85, ztransno85 = FreqStren(
													Rb85_ES.groundManifold,
													Rb85_ES.excitedManifold,
													Rb85_ES.ds,Rb85_ES.dp,
													Dline,'Z')
			

		if rb87frac!=0.0:
			Rb87atom = AC.Rb87
			#Hamiltonian(isotope,transition,gL,Bfield)
			Rb87_ES = ES.Hamiltonian('Rb87',Dline,1.0,Bfield)
			# Rb-87 allowed transitions for light driving sigma minus
			lenergy87, lstrength87, ltransno87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Left')

			# Rb-87 allowed transitions for light driving sigma plus
			renergy87, rstrength87, rtransno87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Right')

			# Rb-87 allowed transitions for light driving sigma plus
			zenergy87, zstrength87, ztransno87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Z')	
													
		if Dline=='D1':
			transitionConst = AC.RbD1Transition
		elif Dline=='D2':
			transitionConst = AC.RbD2Transition

		if (rb85frac!=0.0) and (rb87frac!=0.0):
			AllEnergyLevels = concatenate((lenergy87,lenergy85,
														renergy87,renergy85,
														zenergy87,zenergy85))
		elif (rb85frac!=0.0) and (rb87frac==0.0):
			AllEnergyLevels = concatenate((lenergy85,renergy85,zenergy85))
		elif (rb85frac==0.0) and (rb87frac!=0.0):
			AllEnergyLevels = concatenate((lenergy87,renergy87,zenergy87))

	# Caesium energy levels
	elif Elem=='Cs':
		CsAtom = AC.Cs
		Cs_ES = ES.Hamiltonian('Cs',Dline,1.0,Bfield)

		lenergy, lstrength, ltransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Left')
		renergy, rstrength, rtransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Right')		
		zenergy, zstrength, ztransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Z')

		if Dline=='D1':
			transitionConst = AC.CsD1Transition
		elif Dline =='D2':
			transitionConst = AC.CsD2Transition
		AllEnergyLevels = concatenate((lenergy,renergy,zenergy))

	# Sodium energy levels
	elif Elem=='Na':
		NaAtom = AC.Na
		Na_ES = ES.Hamiltonian('Na',Dline,1.0,Bfield)

		lenergy, lstrength, ltransno = FreqStren(Na_ES.groundManifold,
												 Na_ES.excitedManifold,
												 Na_ES.ds,Na_ES.dp,Dline,
												 'Left')
		renergy, rstrength, rtransno = FreqStren(Na_ES.groundManifold,
												 Na_ES.excitedManifold,
												 Na_ES.ds,Na_ES.dp,Dline,
												 'Right')

		zenergy, zstrength, ztransno = FreqStren(Na_ES.groundManifold,
												 Na_ES.excitedManifold,
												 Na_ES.ds,Na_ES.dp,Dline,
												 'Z')												 
		if Dline=='D1':
			transitionConst = AC.NaD1Transition
		elif Dline=='D2':
			transitionConst = AC.NaD2Transition
		AllEnergyLevels = concatenate((lenergy,renergy))

	#Potassium energy levels <<<<< NEED TO ADD Z-COMPONENT >>>>>
	elif Elem=='K':
		K39frac=1.0-K40frac-K41frac #Potassium-39 fraction
		if K39frac!=0.0:
			K39atom = AC.K39
			K39_ES = ES.Hamiltonian('K39',Dline,1.0,Bfield)
			lenergy39, lstrength39, ltransno39 = FreqStren(
													K39_ES.groundManifold,
													K39_ES.excitedManifold,
													K39_ES.ds,K39_ES.dp,Dline,
													'Left')
			renergy39, rstrength39, rtransno39 = FreqStren(
													K39_ES.groundManifold,
													K39_ES.excitedManifold,
													K39_ES.ds,K39_ES.dp,Dline,
													'Right')

		if K40frac!=0.0:
			K40atom = AC.K40
			K40_ES = ES.Hamiltonian('K40',Dline,1.0,Bfield)
			lenergy40, lstrength40, ltransno40 = FreqStren(
													K40_ES.groundManifold,
													K40_ES.excitedManifold,
													K40_ES.ds,K40_ES.dp,Dline,
													'Left')
			renergy40, rstrength40, rtransno40 = FreqStren(
													K40_ES.groundManifold,
													K40_ES.excitedManifold,
													K40_ES.ds,K40_ES.dp,Dline,
													'Right')
		if K41frac!=0.0:
			K41atom = AC.K41
			K41_ES = ES.Hamiltonian('K41',Dline,1.0,Bfield)
			lenergy41, lstrength41, ltransno41 = FreqStren(
													K41_ES.groundManifold,
													K41_ES.excitedManifold,
													K41_ES.ds,K41_ES.dp,Dline,
													'Left')
			renergy41, rstrength41, rtransno41 = FreqStren(
													K41_ES.groundManifold,
													K41_ES.excitedManifold,
													K41_ES.ds,K41_ES.dp,Dline,
													'Right')
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
	T += 273.15
	DoppTemp += 273.15

	# For thin cells: Don't add Doppler effect, by setting DopplerTemperature to near-zero
	# can then convolve with different velocity distribution later on
		
	d = (array(X)-shift) #Linear detuning
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
	prefactor=2.0*NDensity*dipole**2/(hbar*e0)

	if Elem=='Rb':
		lab85, ldisp85, rab85, rdisp85 = 0,0,0,0
		lab87, ldisp87, rab87, rdisp87 = 0,0,0,0
		if rb85frac!=0.0:
			lab85, ldisp85, rab85, rdisp85, zab85, zdisp85 = add_voigt(d,DoppTemp,
													   Rb85atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno85,lenergy85,lstrength85,
													   rtransno85,renergy85,rstrength85,
													   ztransno85,zenergy85,zstrength85)
		if rb87frac!=0.0:
			lab87, ldisp87, rab87, rdisp87, zab87, zdisp87 = add_voigt(d,DoppTemp,
													   Rb87atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno87,lenergy87,lstrength87,
													   rtransno87,renergy87,rstrength87,
													   ztransno87,zenergy87,zstrength87)
		# Make the parts of the susceptibility
		ChiRealLeft= prefactor*(rb85frac*ldisp85+rb87frac*ldisp87)
		ChiRealRight= prefactor*(rb85frac*rdisp85+rb87frac*rdisp87)
		ChiRealZ = prefactor*(rb85frac*zdisp85 + rb87frac*zdisp87)
		ChiImLeft = prefactor*(rb85frac*lab85+rb87frac*lab87)
		ChiImRight = prefactor*(rb85frac*rab85+rb87frac*rab87)
		ChiImZ = prefactor*(rb85frac*zab85 + rb87frac*zab87)
		

	elif Elem=='Cs':
		lab, ldisp, rab, rdisp, zab, zdisp = add_voigt(d,DoppTemp,CsAtom.mass,wavenumber,
										   gamma,voigtwidth,
										   ltransno,lenergy,lstrength,
										   rtransno,renergy,rstrength,
										   ztransno,zenergy,zstrength)
		ChiRealLeft = prefactor*ldisp
		ChiRealRight = prefactor*rdisp
		ChiRealZ = prefactor*zdisp
		ChiImLeft = prefactor*lab
		ChiImRight = prefactor*rab
		ChiImZ = prefactor*zab
	elif Elem=='Na':
		lab, ldisp, rab, rdisp = add_voigt(d,DoppTemp,NaAtom.mass,wavenumber,
										   gamma,voigtwidth,
										   ltransno,lenergy,lstrength,
										   rtransno,renergy,rstrength,
										   ztransno,zenergy,zstrength)
		ChiRealLeft = prefactor*ldisp
		ChiRealRight = prefactor*rdisp
		ChiRealZ = prefactor*zdisp
		ChiImLeft = prefactor*lab
		ChiImRight = prefactor*rab
		ChiImZ = prefactor*zab
	elif Elem=='K':
		if K39frac!=0.0:
			lab39, ldisp39, rab39, rdisp39 = add_voigt(d,DoppTemp,K39atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno39,lenergy39,
													   lstrength39,rtransno39,
													   renergy39,rstrength39)
		if K40frac!=0.0:
			lab40, ldisp40, rab40, rdisp40 = add_voigt(d,DoppTemp,K40atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno40,lenergy40,
													   lstrength40,rtransno40,
													   renergy40,rstrength40)
		if K41frac!=0.0:
			lab41, ldisp41, rab41, rdisp41 = add_voigt(d,DoppTemp,K41atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
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
	totalChiZ = ChiRealZ + 1.j*ChiImZ
	
	return totalChiLeft, totalChiRight, totalChiZ
	
def get_spectra(X, p_dict, outputs=None):
	""" 
	Calls calc_chi() to get susceptibility, then processes the raw 
	susceptibility into experimentally useful quantities.
	
	Arguments:
	
		X: 		Detuning axis (float, list, or numpy array) in MHz
		p_dict: 	Dictionary containing all parameters (the order of parameters is therefore not important)
				Dictionary keys:
	
				Key				DataType	Unit		Description
				---				---------	----		-----------
				Elem	   			str			--			The chosen alkali element.
				Dline	  			str			--			Specifies which D-line transition to calculate for (D1 or D2)
				
				# Experimental parameters
				Bfield	 			float			Gauss	Magnitude of the applied magnetic field
				T		 			float			Celsius	Temperature used to calculate atomic number density
				GammaBuf   	float			MHz		Extra lorentzian broadening (usually from buffer gas 
																but can be any extra homogeneous broadening)
				shift	  			float			MHz		A global frequency shift of the atomic resonance frequencies

				DoppTemp   	float			Celsius	Temperature linked to the Doppler width (used for
																independent Doppler width and number density)
				Constrain  		bool			--			If True, overides the DoppTemp value and sets it to T

				# Elemental abundancies, where applicable
				rb85frac   		float			%			percentage of rubidium-85 atoms
				K40frac			float			%			percentage of potassium-40 atoms
				K41frac			float			%			percentage of potassium-41 atoms
				
				lcell	  			float			m			length of the vapour cell
				theta0	 		float			degrees	Linear polarisation angle w.r.t. to the x-axis
				Pol				float			%			Percentage of probe beam that drives sigma minus (50% = linear polarisation)
				
				NOTE: If keys are missing from p_dict, default values contained in p_dict_defaults will be loaded.
				
				Extra p_dict values are ignored
	
	Options:
	
		outputs: an iterable (list,tuple...) of strings that defines which spectra are returned, and in which order.
					 If outputs is None (default), a default list of spectra are returned. These are:
						S0,S1,S2,S3,Ix,Iy,nPlus-1,nMinus-1,phi,alphaPlus,alphaMinus,GIPlus,GIMinus
						
	Example Use:
	
		To calculate and quickly plot a transmission spectrum:
			
			import matplotlib.pyplot as plt
			import numpy as np
			
			X = np.arange(-10,10,0.01)*1e3
			p_dict = {'Elem':'Rb', 'Dline':'D2', 'Bfield':0., 'T':100., 'GammaBuf':0., 'shift':0.,
						  'Constrain':True,'rb85frac':1.,'lcell':2e-3 }
			Transmission = get_spectra(X, p_dict, outputs=['S0'])
			
			plt.plot(X,Transmission)
			plt.show()
		
	"""
	
	# get some parameters from p dictionary
	
	# need in try/except or equiv.
	if 'Elem' in p_dict.keys():
		Elem = p_dict['Elem']
	else:
		Elem = p_dict_defaults['Elem']
	if 'Dline' in p_dict.keys():
		Dline = p_dict['Dline']
	else:
		Dline = p_dict_defaults['Dline']
	if 'shift' in p_dict.keys():
		shift = p_dict['shift']
	else:
		shift = p_dict_defaults['shift']
	if 'lcell' in p_dict.keys():
		lcell = p_dict['lcell']
	else:
		lcell = p_dict_defaults['lcell']
	if 'theta0' in p_dict.keys():
		theta0 = p_dict['theta0']
	else:
		theta0 = p_dict_defaults['theta0']
	if 'Pol' in p_dict.keys():
		Pol = p_dict['Pol']
	else:
		Pol = p_dict_defaults['Pol']

	# get wavenumber
	exec('transition = AC.'+Elem+Dline+'Transition')

	wavenumber = transition.wavevectorMagnitude
	
	# convert units
	Pol /= 100. # convert to fraction from %
	theta0 *= pi/180. # convert to radians
	
	
	# Calculate susceptibility
	ChiLeft, ChiRight, ChiZ = calc_chi(X, p_dict)
	
	# Complex refractive index
	nLeft = sqrt(1.0+ChiLeft) #Complex refractive index left hand
	nRight = sqrt(1.0+ChiRight) #Complex refractive index right hand
	nZ = sqrt(1.0+ChiZ) # Complex index for pi transitions


	#Transmission
	TransLeft = exp(-2.0*nLeft.imag*wavenumber*lcell)
	TransRight = exp(-2.0*nRight.imag*wavenumber*lcell)
	
	# Faraday rotation angle (including incident linear polarisation angle)
	phiPlus = wavenumber*nLeft.real*lcell
	phiMinus = wavenumber*nRight.real*lcell
	phi = (phiMinus-phiPlus)/2.0 + theta0 
	
	#Stokes parameters
	PolFactor = 2.0*sqrt(Pol-(Pol**2))
	S0 = (Pol*TransLeft) + ((1.0-Pol)*TransRight)
	# alias:
	Transmission = S0
	S1 = PolFactor*cos(2.0*phi)*exp(-wavenumber*lcell*(nRight.imag
													   +nLeft.imag))
	S2 = PolFactor*sin(2.0*phi)*exp(-wavenumber*lcell*(nRight.imag
													   +nLeft.imag))
	S3 = (Pol*TransLeft) - ((1.0-Pol)*TransRight)
	
	# Intensity after PBS aligned with the x-axis
	Ix = (S1+S0)/2.0
	Iy = (S0-S1)/2.0
	
	# (Real part) refractive indices
	nMinus = nLeft.real
	nPlus = nRight.real

	# Absorption coefficients
	alphaPlus = 2.0*nRight.imag*wavenumber
	alphaMinus = 2.0*nLeft.imag*wavenumber
	
	# Group indices
	d = (array(X)-shift) #Linear detuning
	dnWRTv = derivative(d,nRight.real)
	GIPlus = nRight.real + (X + transition.v0*1.0e-6)*dnWRTv
	dnWRTv = derivative(d,nLeft.real)
	GIMinus = nLeft.real + (X + transition.v0*1.0e-6)*dnWRTv
	 
	if (outputs == None) or ('All' in outputs):
		# Default - return 'all' outputs (as used by GUI)
		return S0,S1,S2,S3,Ix,Iy,nPlus-1,nMinus-1,phi,alphaPlus,alphaMinus,GIPlus,GIMinus
	else:
		# Return the variable names mentioned in the outputs list of strings
		# the strings in outputs must exactly match the local variable names here!
		return [locals()[output_str] for output_str in outputs]

def get_Efield(X, E_in, Chi, p_dict):
	""" 
	Most general form of calculation - return the electric field vector E_out. 
	(Can use Jones matrices to calculate all other experimental spectra from this,
	which may be implemented in future.)
	
	Electric field is in the lab frame, in the X/Y/Z basis:
		- light propagation is along the Z axis, X is horizontal and Y is vertical dimension
		
	To change between x/y and L/R bases, one may use:
			E_left = 1/sqrt(2) * ( E_x + i.E_y )
			E_right = 1/sqrt(2) * ( E_x - i.E_y )
	
	Allows calculation with non-uniform B fields by slicing the cell with total length L into 
	n smaller parts with length L/n - assuming that the field is uniform over L/n,
	which can be checked by convergence testing. E_out can then be used as the new E_in 
	for each subsequent cell slice.
	
	Different to get_spectra() in that the input electric field, E_in, must be defined, 
	and the only quantity returned is the output electric field, E_out.
	
	Arguments:
	
		X:			Detuning in MHz
		E_in:		Array of [E_x, E_y, E_z], complex
						If E_x/y/z are each 1D arrays (with the same dimensions as X) 
						then the polarisation depends on detuning - used e.g. when simulating
						a non-uniform magnetic field
						If E_x/y/z are single (complex) values, then the input polarisation is
						assumed to be independent of the detuning.
		p:			Parameter dictionary - see get_spectra() docstring for details.
	
	Returns:
		
		E_out:	Detuning-dependent output Electric-field vector.
					2D-Array of [E_x, E_y, E_z] where E_x/y/z are 1D arrays with same dimensions
					as the detuning axis.
					The output polarisation always depends on the detuning, in the presence of an
					applied magnetic field.
	"""

	# if not already numpy arrays, convert
	X = array(X)
	E_in = array(E_in)
	print 'Input Efield shape: ',E_in.shape
	# check detuning axis X has the correct dimensions. If not, make it so.
	if E_in.shape != (3,len(X)):
		print 'E field not same size as detuning'
		if E_in.shape == (3,):
			E_in = np.array([np.ones(len(X))*E_in[0],np.ones(len(X))*E_in[1],np.ones(len(X))*E_in[2]])
		else:
			print 'ERROR in method get_Efield(): INPUT ELECTRIC FIELD E_in BADLY DEFINED'
			raise ValueError
	
	print 'New Efield shape: ', E_in.shape

	# fill in required dictionary keys from defaults if not given
	if 'lcell' in p_dict.keys():
		lcell = p_dict['lcell']
	else:
		lcell = p_dict_defaults['lcell']
		
	if 'Elem' in p_dict.keys():
		Elem = p_dict['Elem']
	else:
		Elem = p_dict_defaults['Elem']
	if 'Dline' in p_dict.keys():
		Dline = p_dict['Dline']
	else:
		Dline = p_dict_defaults['Dline']

	# get wavenumber
	exec('transition = AC.'+Elem+Dline+'Transition')
	wavenumber = transition.wavevectorMagnitude
	
	## get magnetic field spherical coordinates
	# defaults to 0,0 i.e. B aligned with kvector of light (Faraday)
	if 'Btheta' in p_dict.keys():
		Btheta = p_dict['Btheta']
	else:
		Btheta = p_dict_defaults['Btheta']
	if 'Bphi' in p_dict.keys():
		Bphi = p_dict['Bphi']
	else:
		Bphi = p_dict_defaults['Bphi']
		
	# Find eigen-vectors for propagation and create rotation matrix
	RotMat = SD.solve_diel(ChiLeft,ChiRight,ChiZ,theta)
	
	#
	### Rotate Electric field of light into frame along B-field (x,y,z)
	#print E_in.T
	#E_fieldframe_xyz = RM.rotate_forward(E_in.T,Bphi,Btheta).T
	#print 'Eff_xyz: ',E_fieldframe_xyz
	#print 'Efieldframe xyz shape: ', E_fieldframe_xyz.shape
	#
	## convert Efield to eigen basis (l/r/z for Faraday)
	##E_fieldframe_lrz = BC.xyz_to_lrz(E_fieldframe_xyz)
	
	# susceptibility (already calculated)
	ChiLeft, ChiRight, ChiZ = Chi

	# Complex refractive index
	nLeft = sqrt(1.0+ChiLeft) #Complex refractive index left hand
	nRight = sqrt(1.0+ChiRight) #Complex refractive index right hand
	nZ = sqrt(1.0+ChiZ) # Complex ref index - pi transitions
	
	
	## propagate electric field in field frame, in eigen basis
	## sign convention <<CHECK>> 
	
	E_Lout = E_fieldframe_lrz[0]*exp(1.j*nLeft*wavenumber*lcell)
	E_Rout = E_fieldframe_lrz[1]*exp(1.j*nRight*wavenumber*lcell)
	E_Zout = E_fieldframe_lrz[2]*exp(1.j*nZ*wavenumber*lcell)
	
	print 'N_L:\n',nLeft
	print 'N_R:\n',nRight
	print 'N_Z:\n',nZ
	
	E_out_lrz = np.array([E_Lout,E_Rout,E_Zout])
	
	## Transform electric field back to field frame, x/y/z basis
	E_out_fieldframe_xyz = BC.lrz_to_xyz(E_out_lrz)
	
	## Rotate electric field back to lab frame
	E_out_lab = RM.rotate_back(E_out_fieldframe_xyz.T,Bphi,Btheta).T
	
	## return electric field vector - can then use Jones matrices to do everything else
	return E_out_lab

def get_spectra2(X, E_in, p_dict, outputs=None):
	""" 
	Calls get_Efield() to get Electric field, then use Jones matrices 
	to calculate experimentally useful quantities.
	
	Arguments:
	
		X: 		Detuning axis (float, list, or numpy array) in MHz
		p_dict: 	Dictionary containing all parameters (the order of parameters is therefore not important)
				Dictionary keys:
	
				Key				DataType	Unit		Description
				---				---------	----		-----------
				Elem	   			str			--			The chosen alkali element.
				Dline	  			str			--			Specifies which D-line transition to calculate for (D1 or D2)
				
				# Experimental parameters
				Bfield	 			float			Gauss	Magnitude of the applied magnetic field
				T		 			float			Celsius	Temperature used to calculate atomic number density
				GammaBuf   	float			MHz		Extra lorentzian broadening (usually from buffer gas 
																but can be any extra homogeneous broadening)
				shift	  			float			MHz		A global frequency shift of the atomic resonance frequencies

				DoppTemp   	float			Celsius	Temperature linked to the Doppler width (used for
																independent Doppler width and number density)
				Constrain  		bool			--			If True, overides the DoppTemp value and sets it to T

				# Elemental abundancies, where applicable
				rb85frac   		float			%			percentage of rubidium-85 atoms
				K40frac			float			%			percentage of potassium-40 atoms
				K41frac			float			%			percentage of potassium-41 atoms
				
				lcell	  			float			m			length of the vapour cell
				theta0	 		float			degrees	Linear polarisation angle w.r.t. to the x-axis
				Pol				float			%			Percentage of probe beam that drives sigma minus (50% = linear polarisation)
				
				NOTE: If keys are missing from p_dict, default values contained in p_dict_defaults will be loaded.
				
				Extra p_dict values are ignored
	
	Options:
	
		outputs: an iterable (list,tuple...) of strings that defines which spectra are returned, and in which order.
					 If outputs is None (default), a default list of spectra are returned. These are:
						S0,S1,S2,S3,Ix,Iy,nPlus-1,nMinus-1,phi,alphaPlus,alphaMinus,GIPlus,GIMinus
						
	Example Use:
	
		To calculate and quickly plot a transmission spectrum:
			
			import matplotlib.pyplot as plt
			import numpy as np
			
			X = np.arange(-10,10,0.01)*1e3
			p_dict = {'Elem':'Rb', 'Dline':'D2', 'Bfield':0., 'T':100., 'GammaBuf':0., 'shift':0.,
						  'Constrain':True,'rb85frac':1.,'lcell':2e-3 }
			Transmission = get_spectra(X, p_dict, outputs=['S0'])
			
			plt.plot(X,Transmission)
			plt.show()
		
	"""
	
	# get some parameters from p dictionary
	
	# need in try/except or equiv.
	if 'Elem' in p_dict.keys():
		Elem = p_dict['Elem']
	else:
		Elem = p_dict_defaults['Elem']
	if 'Dline' in p_dict.keys():
		Dline = p_dict['Dline']
	else:
		Dline = p_dict_defaults['Dline']
	if 'shift' in p_dict.keys():
		shift = p_dict['shift']
	else:
		shift = p_dict_defaults['shift']
	if 'lcell' in p_dict.keys():
		lcell = p_dict['lcell']
	else:
		lcell = p_dict_defaults['lcell']
	if 'theta0' in p_dict.keys():
		theta0 = p_dict['theta0']
	else:
		theta0 = p_dict_defaults['theta0']
	if 'Pol' in p_dict.keys():
		Pol = p_dict['Pol']
	else:
		Pol = p_dict_defaults['Pol']

	# get wavenumber
	exec('transition = AC.'+Elem+Dline+'Transition')

	wavenumber = transition.wavevectorMagnitude
	
	# Calculate Susceptibility
	ChiLeft, ChiRight, ChiZ = calc_chi(X, p_dict)
	Chi = [ChiLeft, ChiRight, ChiZ]
	
	# Complex refractive index
	nLeft = sqrt(1.0+ChiLeft) #Complex refractive index left hand
	nRight = sqrt(1.0+ChiRight) #Complex refractive index right hand
	nZ = sqrt(1.0+ChiZ) # Complex index for pi transitions

	# Calculate E_field
	E_out = get_Efield(X, E_in, Chi, p_dict)
	

	## Apply Jones matrices
	
	# Transmission - total intensity - just E_out**2 / E_in**2
	E_in = np.array(E_in)
	if E_in.shape == (3,):
			E_in = np.array([np.ones(len(X))*E_in[0],np.ones(len(X))*E_in[1],np.ones(len(X))*E_in[2]])
	
	print 'Input E field shape: ',E_in.shape
	print 'Output E field shape: ',E_out.shape
	
	I_in = (E_in * E_in.conjugate()).sum(axis=0)
	
	print type(E_out)
	S0 = (E_out * E_out.conjugate()).sum(axis=0) / I_in
	Transmission = S0
	
## Doesn't make sense when B and k not aligned
	#TransLeft = exp(-2.0*nLeft.imag*wavenumber*lcell)
	#TransRight = exp(-2.0*nRight.imag*wavenumber*lcell)
	#TransZ = exp(-2.0*nZ.imag*wavenumber*lcell)
	
	# Faraday rotation angle (including incident linear polarisation angle)
	#phiPlus = wavenumber*nLeft.real*lcell
	#phiMinus = wavenumber*nRight.real*lcell
	#phi = (phiMinus-phiPlus)/2.0 + theta0 
##
	
	print ' Check for z-axis e-field == 0?:\n',E_out[2]
	
	#Stokes parameters
	Ex = np.array(JM.HorizPol_xy * E_out[:2])
	Ix =  (Ex * Ex.conjugate()).sum(axis=0) / I_in
	Ey =  np.array(JM.VertPol_xy * E_out[:2])
	#print type(Ey)
	#print Ey.shape
	#print Ey
	Iy =  (Ey * Ey.conjugate()).sum(axis=0) / I_in
	#print Iy.real
	
	# S1 = Ix - Iy
	S1 = Ix - Iy
	
	# S2 = I_plus45 - I_minus45
	E_P45 =  np.array(JM.LPol_P45_xy * E_out[:2])
	E_M45 =  np.array(JM.LPol_M45_xy * E_out[:2])
	I_P45 = (E_P45 * E_P45.conjugate()).sum(axis=0) / I_in
	I_M45 = (E_M45 * E_M45.conjugate()).sum(axis=0) / I_in
	
	S2 = I_P45 - I_M45
	
	# S3 = I- - I+
	# change to circular basis
	# LCP(RCP) light drives sigma+(-) transitions
	E_out_lrz = BC.xyz_to_lrz(E_out)
	El =  np.array(JM.CPol_L_lr * E_out_lrz[:2])
	Er =  np.array(JM.CPol_R_lr * E_out_lrz[:2])
	Il = (El * El.conjugate()).sum(axis=0) / I_in
	Ir = (Er * Er.conjugate()).sum(axis=0) / I_in
	
	S3 = Ir - Il
	
	Ix = Ix.real
	Iy = Iy.real
##
	## (Real part) refractive indices
	#nMinus = nLeft.real
	#nPlus = nRight.real

	## Absorption coefficients
	#alphaPlus = 2.0*nRight.imag*wavenumber
	#alphaMinus = 2.0*nLeft.imag*wavenumber
	#
	# Group indices
	#d = (array(X)-shift) #Linear detuning
	#dnWRTv = derivative(d,nRight.real)
	#GIPlus = nRight.real + (X + transition.v0*1.0e-6)*dnWRTv
	#dnWRTv = derivative(d,nLeft.real)
	#GIMinus = nLeft.real + (X + transition.v0*1.0e-6)*dnWRTv
	# 
	if (outputs == None) or ('All' in outputs):
		# Default - return 'all' outputs (as used by GUI)
		return S0.real,S1.real,S2.real,S3.real,Ix.real,Iy.real,Il.real,Ir.real
	else:
		# Return the variable names mentioned in the outputs list of strings
		# the strings in outputs must exactly match the local variable names here!
		return [locals()[output_str] for output_str in outputs]
		
def spectrum(X,Elem='Rb',OutputType='Ix',Bfield=0,T=10,lcell=0,rb85frac=0,DoppTemp=10,
		theta0=0,Pol=0.5,shift=0,GammaBuf=0,Constrain=1,Dline='D2',precision=10,
		K40frac=0,K41frac=0):
	""" 
	! LEGACY !
	wrapper for get_spectra() that takes in the same 
	parameters as the old spectrum() method, for backwards-compatibility
	
	For any new user, we recommend using the get_spectra() method instead,
	or using the calculate() method in elecsus_methods (which is essentially a wrapper for get_spectra()).
	
	The precision keyword is deprecated and has no effect whatsoever.
	"""
	if OutputType == 'Transmission (S0)':
		OutputType = 'S0'
	
	# create parameter dictionary
	p_dict = {'Elem':Elem, 'Dline':Dline, 'Bfield':Bfield, 'T':T, 'GammaBuf':GammaBuf, 'shift':shift,
						  'Constrain':Constrain,'rb85frac':rb85frac,'lcell':lcell*1e-3,'DoppTemp':DoppTemp,
						  'theta0':theta0, 'Pol':Pol, 'K40frac':K40frac, 'K41frac':K41frac }
	
	# convert units
	Xin = array(X)*1e3
	
	# call get_spectra
	spec = get_spectra(Xin,p_dict,outputs=[OutputType])
	if len(spec) == 1:
		return spec[0] ## Output single spectrum
	else:
		return spec ## Output == 'All'
	
def output_list():
	""" Helper method that returns a list of all possible variables that get_spectra can return """
	tstr = " \
	All possible outputs from the get_spectra method: \n\n\
	Variable Name		Description \n \
	S0						Total transmission through the cell (Ix + Iy) \n\
	S1						Stokes parameter - Ix - Iy \n\
	S2						Stokes parameter - I_45 - I_-45 \n\
	S3						Stokes parameter - I- - I+ \n\
	TransLeft			Transmission of only left-circularly polarised light \n\
	TransRight			Transmission of only right-circularly polarised light \n\
	ChiLeft				Complex susceptibility of left-circularly polarised light \n\
	ChiRight				Complex susceptibility of right-circularly polarised light \n\
	nLeft					Complex Refractive index of left-circularly polarised light \n\
	nRight				Complex Refractive index of right-circularly polarised light \n\
	phiPlus				Rotation of linear polarisation caused by sigma-plus transitions \n\
	phiMinus				Rotation of linear polarisation caused by sigma-minus transitions \n\
	phi					Total rotation of linear polarisation \n\
	Ix						Intensity of light transmitted through a linear polariser aligned with the x-axis \n\
	Iy						Intensity of light transmitted through a linear polariser aligned with the y-axis \n\
	alphaPlus			Absorption coefficient due to sigma-plus transitions \n\
	alphaMinus			Absorption coefficient due to sigma-minus transitions \n\
	GIMinus				Group index of left-circularly polarised light \n\
	GIPlus				Group index of right-circularly polarised light \n\
	"	
	print tstr

def test1():
### 1. Fig 3 of Generalised treatment ... Rotondaro JOSAB 2015 paper
	d = np.arange(-10000,10000,100)
	#Voigt
	p_dict = {'Bfield':300,'rb85frac':1,'Btheta':0,'lcell':75e-3,'T':58,'Dline':'D2','Elem':'Cs'}
	TF = get_spectra2(d,[1,0,0],p_dict,outputs=['Iy'])
	
	#check vs old elecsus
	from elecsus.libs import spectra as old_spec
	
	TF_old = old_spec.get_spectra(d,p_dict,outputs=['Iy'])
	
	index = 0 # Iy
	
	fig = plt.figure("Faraday comparison")
	ax1 = fig.add_subplot(111)
	ax1.plot(d,TF[index],'r',lw=2,label='Faraday')
	ax1.plot(d,TF_old[0],'k--',lw=2,label='Vanilla ElecSus')
	
	ax1.legend(loc=0)
	
	ax1.set_xlabel('Detuning (MHz)')
	ax1.set_ylabel('Transmission')
	
	plt.show()

def test2():
### 2. Fig 4 of General.... paper
	d = np.arange(-12000,12000,25)
	#Voigt
	#p_dict = {'Bfield':700,'rb85frac':1,'Btheta':90*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	p_dict = {'Bfield':700,'rb85frac':1,'Btheta':np.pi/2,'Bphi':0*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	pol = 1./np.sqrt(2)*np.array([0.0,1.0,0.0])
	TVx = get_spectra2(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0'])
	
	fig2 = plt.figure()
	ax1a = fig2.add_subplot(411)
	ax2a = fig2.add_subplot(412,sharex=ax1a)
	ax3a = fig2.add_subplot(413,sharex=ax1a)
	ax4a = fig2.add_subplot(414,sharex=ax1a)
	
	ax1a.plot(d,TVx[0],'r',lw=2,label='Voigt, I-45')
	ax2a.plot(d,TVx[1],'b',lw=2,label='Voigt, X-polarised')
	ax3a.plot(d,TVx[2],'r',lw=2,label='Voigt, X-polarised')
	ax4a.plot(d,TVx[3],'b',lw=2,label='Voigt, X-polarised')
	ax4a.plot(d,TVx[2]+TVx[3],'k:',lw=2.5,label='Voigt, X-polarised')
	ax4a.plot(d,TVx[4],'g--',lw=1.5,label='Voigt, X-polarised')
	
	
	ax4a.set_xlabel('Detuning (MHz)')
	ax1a.set_ylabel('I -45')
	ax2a.set_ylabel('I +45')
	ax3a.set_ylabel('Ix')
	ax4a.set_ylabel('Iy')
	
	plt.show()
	
if __name__ == '__main__':
	test2()