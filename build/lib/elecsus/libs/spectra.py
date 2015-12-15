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
from scipy.special import wofz
from scipy.interpolate import interp1d
from FundamentalConstants import *
from numberDensityEqs import *
from tools import derivative

import EigenSystem as ES
import AtomConstants as AC

# Default values for parameters
p_dict_defaults = {'Elem':'Rb', 'Dline':'D2', 'Bfield':0., 'T':20., 'GammaBuf':0., 'shift':0.,
						  'Constrain':True,'rb85frac':72.17,'lcell':75e-3,'DoppTemp':20.,
						  'theta0':0., 'Pol':50., 'K40frac':0.01, 'K41frac':6.73 }

def FreqStren(groundLevels,excitedLevels,groundDim,
			  excitedDim,Dline,hand):
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
				transitionFrequency[transNo] = int((-groundLevels[gg][0].real
											  +excitedLevels[ee][0].real))
				# We choose to perform the ground manifold reduction (see
				# equation (4) in manual) here for convenience.
				transitionStrength[transNo] = 1./3*1./groundDim*cleb2
				transNo += 1

	#print 'No transitions (ElecSus): ',transNo
	return transitionFrequency, transitionStrength, transNo

def add_voigt(d,DoppTemp,atomMass,wavenumber,gamma,voigtwidth,
		ltransno,lenergy,lstrength,rtransno,
		renergy,rstrength):
	xpts = len(d)
	npts = 2*voigtwidth+1
	detune = 2.0*pi*1.0e6*(arange(npts)-voigtwidth) #Angular detuning
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
		lab += lstrength[line]*f_ab(2.0*pi*(d-xc)*1.0e6)
		ldisp += lstrength[line]*f_disp(2.0*pi*(d-xc)*1.0e6)
	rab = zeros(xpts)
	rdisp = zeros(xpts)
	for line in xrange(rtransno+1):
		xc = renergy[line]
		rab += rstrength[line]*f_ab(2.0*pi*(d-xc)*1.0e6)
		rdisp += rstrength[line]*f_disp(2.0*pi*(d-xc)*1.0e6)
	return lab, ldisp, rab, rdisp

	
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
												 'Left')
		renergy, rstrength, rtransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Right')

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
												 'Left')
		renergy, rstrength, rtransno = FreqStren(Na_ES.groundManifold,
												 Na_ES.excitedManifold,
												 Na_ES.ds,Na_ES.dp,Dline,
												 'Right')

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
			lab85, ldisp85, rab85, rdisp85 = add_voigt(d,DoppTemp,
													   Rb85atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno85,lenergy85,
													   lstrength85,rtransno85,
													   renergy85,rstrength85)
		if rb87frac!=0.0:
			lab87, ldisp87, rab87, rdisp87 = add_voigt(d,DoppTemp,
													   Rb87atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno87,lenergy87,
													   lstrength87,rtransno87,
													   renergy87,rstrength87)
		# Make the parts of the susceptibility
		ChiRealLeft= prefactor*(rb85frac*ldisp85+rb87frac*ldisp87)
		ChiRealRight= prefactor*(rb85frac*rdisp85+rb87frac*rdisp87)
		ChiImLeft = prefactor*(rb85frac*lab85+rb87frac*lab87)
		ChiImRight = prefactor*(rb85frac*rab85+rb87frac*rab87)

	elif Elem=='Cs':
		lab, ldisp, rab, rdisp = add_voigt(d,DoppTemp,CsAtom.mass,wavenumber,
										   gamma,voigtwidth,ltransno,
										   lenergy,lstrength,rtransno,renergy,
										   rstrength)
		ChiRealLeft= prefactor*ldisp
		ChiRealRight= prefactor*rdisp
		ChiImLeft = prefactor*lab
		ChiImRight = prefactor*rab
	elif Elem=='Na':
		lab, ldisp, rab, rdisp = add_voigt(d,DoppTemp,NaAtom.mass,wavenumber,
										   gamma,voigtwidth,ltransno,
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
	
	return totalChiLeft, totalChiRight
	
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
	ChiLeft, ChiRight = calc_chi(X, p_dict)
	
	# Complex refractive index
	nLeft = sqrt(1.0+ChiLeft) #Complex refractive index left hand
	nRight = sqrt(1.0+ChiRight) #Complex refractive index right hand


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

def get_Efield(X, E_in, p_dict):
	""" 
	Most general form of calculation - return the electric field vector E_out. 
	(Can use Jones matrices to calculate all other experimental spectra from this,
	which may be implemented in future.)
	
	Electric field is in the Left/Right circular basis (for the moment - this may change
	when a vectorial B-field is introduced in a later version). 
	To change between x/y and L/R bases, one may use:
			E_left = 1/sqrt(2) * ( E_x + i.E_y )
			E_right = 1/sqrt(2) * ( E_x - i.E_y )
	Ignoring Z-component of E-field (assuming plane waves)
	
	Allows calculation with non-uniform B fields by slicing the cell with total length L into 
	n smaller parts with length L/n - assuming that the field is uniform over L/n,
	which can be checked by convergence testing. E_out can then be used as the new E_in 
	for each subsequent cell slice.
	
	Different to get_spectra() in that the input electric field, E_in, must be defined, 
	and the only quantity returned is the output electric field, E_out.
	
	Arguments:
	
		X:			Detuning in MHz
		E_in:		2D Array of [E_left, E_right] where E_left/right are each 1D arrays 
					(in general, with complex elements) with the same dimensions as X.
		p:			Parameter dictionary - see get_spectra() docstring for details.
	
	Returns:
		
		E_out:	Same dimensions as E_in, but propagated for a length L in a uniform B field
	"""

	lcell = p_dict['lcell']

	# get wavenumber
	Elem = p_dict['elem']
	Dline = p_dict['Dline']
	exec('transition = AC.'+Elem+Dline+'Transition')
	wavenumber = transition.wavevectorMagnitude

	# get susceptibility
	ChiLeft, ChiRight = calc_chi(X, p_dict)

	# Complex refractive index
	nLeft = sqrt(1.0+ChiLeft) #Complex refractive index left hand
	nRight = sqrt(1.0+ChiRight) #Complex refractive index right hand
	
	## propagate electric field
	E_L_out = E_in[0]*exp(1.j*nLeft*wavenumber*lcell)
	E_R_out = E_in[1]*exp(1.j*nRight*wavenumber*lcell)
	E_out = [E_L_out, E_R_out]
	
	## return electric field vector - can then use Jones matrices to do everything else
	return E_out
	
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
