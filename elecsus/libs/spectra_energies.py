# Copyright 2014-2016 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
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

"""
Module containing functions to calculate the spectra

Constructs the electric susceptibility and then returns
the spectrum requested

Last Updated 2016-10-12 JK
"""

from numpy import zeros,sqrt,pi,dot,exp,sin,cos,array,amax,arange,concatenate,argmin
import numpy as np
from scipy.special import wofz
from scipy.interpolate import interp1d
from FundamentalConstants import *
from numberDensityEqs import *
#from tools import derivative

import EigenSystem as ES
import AtomConstants as AC
import RotationMatrices as RM
import BasisChanger as BC
import JonesMatrices as JM
import solve_dielectric as SD

# direct testing
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
							'rb85frac':72.17, 'K40frac':0.01, 'K41frac':6.73,
							'BoltzmannFactor':True}

def FreqStren(groundLevels,excitedLevels,groundDim,
			  excitedDim,Dline,hand,BoltzmannFactor=True,T=293.16):
	""" 
	Calculate transition frequencies and strengths by taking dot
	products of relevant parts of the ground / excited state eigenvectors 
	"""
	
	transitionFrequency = zeros(groundDim*2*groundDim) #Initialise lists
	transitionStrength = zeros(groundDim*2*groundDim)
	transNo = 0
	
	groundLevelsT = []
	excitedLevelsT = []
	transitionFrequencyT = []
	
	## Boltzmann factor: -- only really needed when energy splitting of ground state is large
	if BoltzmannFactor:
		groundEnergies = array(groundLevels)[:,0].real
		lowestEnergy = groundLevels[0][argmin(groundLevels,axis=0)[0]].real
		BoltzDist = exp(-(groundEnergies-lowestEnergy)*h*1.e6/(kB*T)) #Not normalised
		#print BoltzDist
		BoltzDist = BoltzDist/BoltzDist.sum() # Normalised

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
			if cleb2 > 0.0005: #If negligible don't calculate.
				transitionFrequency[transNo] = int((-groundLevels[gg][0].real
											  +excitedLevels[ee][0].real))
				# We choose to perform the ground manifold reduction (see
				# equation (4) in manual) here for convenience.
				transitionFrequencyT.append(transitionFrequency[transNo])
				groundLevelsT.append(groundLevels[gg][0].real)
				excitedLevelsT.append(excitedLevels[ee][0].real)
				
				## Boltzmann factor:
				if BoltzmannFactor: 
					transitionStrength[transNo] = 1./3 * cleb2 * BoltzDist[gg]
				else:
					transitionStrength[transNo] = 1./3 * cleb2 * 1./groundDim
				
				transNo += 1

	#print 'No transitions (ElecSus): ',transNo
	return np.array(transitionFrequencyT), np.array(transitionStrength), transNo, np.array(groundLevelsT), np.array(excitedLevelsT)

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

def calc_chi_energies(X, p_dict,verbose=False):			   
	"""
	Computes the complex susceptibility for sigma plus/minus and pi transitions as a 1D array

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
	if 'BoltzmannFactor' in p_dict.keys():
		BoltzmannFactor =  p_dict['BoltzmannFactor']
	else:
		BoltzmannFactor =  p_dict_defaults['BoltzmannFactor']
	
	
	if verbose: print 'Temperature: ', T, '\tBfield: ', Bfield
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
													Dline,'Left',BoltzmannFactor,T+273.16)		  

			# Rb-85 allowed transitions for light driving sigma plus
			renergy85, rstrength85, rtransno85 = FreqStren(
													Rb85_ES.groundManifold,
													Rb85_ES.excitedManifold,
													Rb85_ES.ds,Rb85_ES.dp,
													Dline,'Right',BoltzmannFactor,T+273.16)
			
			# Rb-85 allowed transitions for light driving sigma plus
			zenergy85, zstrength85, ztransno85 = FreqStren(
													Rb85_ES.groundManifold,
													Rb85_ES.excitedManifold,
													Rb85_ES.ds,Rb85_ES.dp,
													Dline,'Z',BoltzmannFactor,T+273.16)
			
		
		if rb87frac!=0.0:
			Rb87atom = AC.Rb87
			#Hamiltonian(isotope,transition,gL,Bfield)
			Rb87_ES = ES.Hamiltonian('Rb87',Dline,1.0,Bfield)
			# Rb-87 allowed transitions for light driving sigma minus
			lenergy87, lstrength87, ltransno87, lgl87, lel87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Left',BoltzmannFactor,T+273.16)

			# Rb-87 allowed transitions for light driving sigma plus
			renergy87, rstrength87, rtransno87, rgl87, rel87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Right',BoltzmannFactor,T+273.16)

			# Rb-87 allowed transitions for light driving pi
			zenergy87, zstrength87, ztransno87, zgl87, zel87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Z',BoltzmannFactor,T+273.16)
													
			
													
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
			
		return lenergy87, lstrength87, ltransno87, lgl87, lel87, \
					renergy87, rstrength87, rtransno87, rgl87, rel87, \
					zenergy87, zstrength87, ztransno87, zgl87, zel87

	'''
	# Caesium energy levels
	elif Elem=='Cs':
		CsAtom = AC.Cs
		Cs_ES = ES.Hamiltonian('Cs',Dline,1.0,Bfield)

		lenergy, lstrength, ltransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Left',BoltzmannFactor,T+273.16)
		renergy, rstrength, rtransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Right',BoltzmannFactor,T+273.16)		
		zenergy, zstrength, ztransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Z',BoltzmannFactor,T+273.16)

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
												 'Left',BoltzmannFactor,T+273.16)
		renergy, rstrength, rtransno = FreqStren(Na_ES.groundManifold,
												 Na_ES.excitedManifold,
												 Na_ES.ds,Na_ES.dp,Dline,
												 'Right',BoltzmannFactor,T+273.16)
		zenergy, zstrength, ztransno = FreqStren(Na_ES.groundManifold,
												 Na_ES.excitedManifold,
												 Na_ES.ds,Na_ES.dp,Dline,
												 'Z',BoltzmannFactor,T+273.16)												 
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
													'Left',BoltzmannFactor,T+273.16)
			renergy39, rstrength39, rtransno39 = FreqStren(
													K39_ES.groundManifold,
													K39_ES.excitedManifold,
													K39_ES.ds,K39_ES.dp,Dline,
													'Right',BoltzmannFactor,T+273.16)
			zenergy39, zstrength39, ztransno39 = FreqStren(K39_ES.groundManifold,
													K39_ES.excitedManifold,
													K39_ES.ds,K39_ES.dp,Dline,
													'Z',BoltzmannFactor,T+273.16)
		if K40frac!=0.0:
			K40atom = AC.K40
			K40_ES = ES.Hamiltonian('K40',Dline,1.0,Bfield)
			
			lenergy40, lstrength40, ltransno40 = FreqStren(
													K40_ES.groundManifold,
													K40_ES.excitedManifold,
													K40_ES.ds,K40_ES.dp,Dline,
													'Left',BoltzmannFactor,T+273.16)
			renergy40, rstrength40, rtransno40 = FreqStren(
													K40_ES.groundManifold,
													K40_ES.excitedManifold,
													K40_ES.ds,K40_ES.dp,Dline,
													'Right',BoltzmannFactor,T+273.16)
			zenergy40, zstrength40, ztransno40 = FreqStren(K40_ES.groundManifold,
													K40_ES.excitedManifold,
													K40_ES.ds,K40_ES.dp,Dline,
													'Z',BoltzmannFactor,T+273.16)
		if K41frac!=0.0:
			K41atom = AC.K41
			K41_ES = ES.Hamiltonian('K41',Dline,1.0,Bfield)
			
			lenergy41, lstrength41, ltransno41 = FreqStren(
													K41_ES.groundManifold,
													K41_ES.excitedManifold,
													K41_ES.ds,K41_ES.dp,Dline,
													'Left',BoltzmannFactor,T+273.16)
			renergy41, rstrength41, rtransno41 = FreqStren(
													K41_ES.groundManifold,
													K41_ES.excitedManifold,
													K41_ES.ds,K41_ES.dp,Dline,
													'Right',BoltzmannFactor,T+273.16)
			zenergy41, zstrength41, ztransno41 = FreqStren(K41_ES.groundManifold,
													K41_ES.excitedManifold,
													K41_ES.ds,K41_ES.dp,Dline,
													'Z',BoltzmannFactor,T+273.16)
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
		lab85, ldisp85, rab85, rdisp85, zab85, zdisp85 = 0,0,0,0,0,0
		lab87, ldisp87, rab87, rdisp87, zab87, zdisp87 = 0,0,0,0,0,0
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
		lab, ldisp, rab, rdisp, zab, zdisp = 0,0,0,0,0,0
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
		lab, ldisp, rab, rdisp, zab, zdisp = 0,0,0,0,0,0
		lab, ldisp, rab, rdisp, zab, zdisp = add_voigt(d,DoppTemp,NaAtom.mass,wavenumber,
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
		lab39, ldisp39, rab39, rdisp39, zab39, zdisp39 = 0,0,0,0,0,0
		lab40, ldisp40, rab40, rdisp40, zab40, zdisp40 = 0,0,0,0,0,0
		lab41, ldisp41, rab41, rdisp41, zab41, zdisp41 = 0,0,0,0,0,0
		if K39frac!=0.0:
			lab39, ldisp39, rab39, rdisp39, zab39, zdisp39 = add_voigt(d,DoppTemp,K39atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno39,lenergy39,lstrength39,
													   rtransno39,renergy39,rstrength39,
													   ztransno39,zenergy39,zstrength39)
		if K40frac!=0.0:
			lab40, ldisp40, rab40, rdisp40, zab40, zdisp40 = add_voigt(d,DoppTemp,K40atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno40,lenergy40,lstrength40,
													   rtransno40,renergy40,rstrength40,
													   ztransno40,zenergy40,zstrength40)
		if K41frac!=0.0:
			lab41, ldisp41, rab41, rdisp41, zab41, zdisp41 = add_voigt(d,DoppTemp,K41atom.mass,
													   wavenumber,gamma,
													   voigtwidth,
													   ltransno41,lenergy41,lstrength41,
													   rtransno41,renergy41,rstrength41,
													   ztransno41,zenergy41,zstrength41)

		ChiRealLeft = prefactor*(K39frac*ldisp39+K40frac\
									*ldisp40+K41frac*ldisp41)
		ChiRealRight = prefactor*(K39frac*rdisp39+K40frac\
									*rdisp40+K41frac*rdisp41)
		ChiRealZ = prefactor*(K39frac*zdisp39+K40frac\
									*zdisp40+K41frac*zdisp41)
		ChiImLeft = prefactor*(K39frac*lab39+K40frac*lab40+K41frac*lab41)
		ChiImRight = prefactor*(K39frac*rab39+K40frac*rab40+K41frac*rab41)
		ChiImZ = prefactor*(K39frac*zab39+K40frac*zab40+K41frac*zab41)

	# Reconstruct total susceptibility and index of refraction
	totalChiPlus = ChiRealLeft + 1.j*ChiImLeft
	totalChiMinus = ChiRealRight + 1.j*ChiImRight
	totalChiZ = ChiRealZ + 1.j*ChiImZ
	
	return totalChiPlus, totalChiMinus, totalChiZ
	'''

def get_Efield(X, E_in, Chi, p_dict, verbose=False):
	""" 
	Most general form of calculation - return the electric field vector E_out. 
	Can use Jones matrices to calculate all other experimental spectra from 
	this, as in the get_spectra2() method
	
	Electric field is in the lab frame, in the X/Y/Z basis:
		- light propagation is along the Z axis, X is horizontal and Y is vertical dimension
		
	To change between x/y and L/R bases, one may use:
			E_left = 1/sqrt(2) * ( E_x - i.E_y )
			E_right = 1/sqrt(2) * ( E_x + i.E_y )
	(See BasisChanger.py module)
	
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
		Chi:		(3,len(X)) array of susceptibility values, for sigma+, sigma-, and pi transitions
		p_dict:	Parameter dictionary - see get_spectra() docstring for details.
	
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
	if verbose:
		print 'Electric field input:'
		print E_in
	
	#print 'Input Efield shape: ',E_in.shape
	# check detuning axis X has the correct dimensions. If not, make it so.
	if E_in.shape != (3,len(X)):
		#print 'E field not same size as detuning'
		if E_in.shape == (3,):
			E_in = np.array([np.ones(len(X))*E_in[0],np.ones(len(X))*E_in[1],np.ones(len(X))*E_in[2]])
		else:
			raise ValueError( 'ERROR in method get_Efield(): INPUT ELECTRIC FIELD E_in BADLY DEFINED' )
	
	#print 'New Efield shape: ', E_in.shape

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
		
	# get susceptibility (already calculated, input to this method)
	ChiPlus, ChiMinus, ChiZ = Chi

	# Rotate initial Electric field so that B field lies in x-z plane
	# (Effective polarisation rotation)
	E_xz = RM.rotate_around_z(E_in.T,Bphi)
		
	# Find eigen-vectors for propagation and create rotation matrix
	RM_ary, n1, n2 = SD.solve_diel(ChiPlus,ChiMinus,ChiZ,Btheta)

	# propagation matrix
	PropMat = np.array(
				[	[exp(1.j*n1*wavenumber*lcell),np.zeros(len(n1)),np.zeros(len(n1))],
					[np.zeros(len(n1)),exp(1.j*n2*wavenumber*lcell),np.zeros(len(n1))],
					[np.zeros(len(n1)),np.zeros(len(n1)),np.ones(len(n1))]	])
	#print 'prop matrix shape:',PropMat.T.shape
	#print 'prop mat [0]: ', PropMat.T[0]
	
	# calcualte output field - a little messy to make it work nicely with array operations
	# - need to play around with matrix dimensions a bit
	# Effectively this does this, element-wise: E_out_xz = RotMat.I * PropMat * RotMat * E_xz
	
	E_xz = np.reshape(E_xz.T, (len(X),3,1))
	
	E_out_xz = np.zeros((len(X),3,1),dtype='complex')
	E_out = np.zeros_like(E_out_xz)
	for i in range(len(X)):
		#print 'Propagation Matrix:\n',PropMat.T[i]
		#print 'Rotation matrix:\n',RM_ary.T[i]
		#inverse rotation matrix
		RMI_ary = np.matrix(RM_ary.T[i].T).I
		#print 'Inverse rotation matrix:\n',RMI_ary
		
		E_out_xz[i] = RMI_ary * np.matrix(PropMat.T[i]) * np.matrix(RM_ary.T[i].T)*np.matrix(E_xz[i]) 
		#print 'E out xz i: ',E_out_xz[i].T
		E_out[i] = RM.rotate_around_z(E_out_xz[i].T[0],-Bphi)
	
	#print 'E out [0]: ',E_out[0]
	#print 'E out shape: ',E_out.shape
	
	## return electric field vector - can then use Jones matrices to do everything else
	return E_out.T[0], np.matrix(RM_ary.T[i])

def get_spectra2(X, E_in, p_dict, outputs=None):
	""" 
	Calls get_Efield() to get Electric field, then use Jones matrices 
	to calculate experimentally useful quantities.
	
	Alias for the get_spectra2 method in libs.spectra.
	
	Inputs:
		detuning_range [ numpy 1D array ] 
			The independent variable and defines the detuning points over which to calculate. Values in MHz
		
		E_in [ numpy 1/2D array ] 
			Defines the input electric field vector in the xyz basis. The z-axis is always the direction of propagation (independent of the magnetic field axis), and therefore the electric field should be a plane wave in the x,y plane. The array passed to this method should be in one of two formats:
				(1) A 1D array of (Ex,Ey,Ez) which is the input electric field for all detuning values;
				or
				(2) A 2D array with dimensions (3,len(detuning_range)) - i.e. each detuning has a different electric field associated with it - which will happen on propagation through a birefringent/dichroic medium
		
		p_dict [ dictionary ]
			Dictionary containing all parameters (the order of parameters is therefore not important)
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
		
		outputs: an iterable (list,tuple...) of strings that defines which spectra are returned, and in which order.
			If not specified, defaults to None, in which case a default set of outputs is returned, which are:
				S0, S1, S2, S3, Ix, Iy, I_P45, I_M45, alphaPlus, alphaMinus, alphaZ
	
	Returns:
		A list of output arrays as defined by the 'outputs' keyword argument.
		
		
	Example usage:
		To calculate the room temperature absorption of a 75 mm long Cs reference cell in an applied magnetic field of 100 G aligned along the direction of propagation (Faraday geometry), between -10 and +10 GHz, with an input electric field aligned along the x-axis:
		
		detuning_range = np.linspace(-10,10,1000)*1e3 # GHz to MHz conversion
		E_in = np.array([1,0,0])
		p_dict = {'Elem':'Cs', 'Dline':'D2', 'Bfield':100, 'T':21, 'lcell':75e-3}
		
		[Transmission] = calculate(detuning_range,E_in,p_dict,outputs=['S0'])
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
	ChiPlus, ChiMinus, ChiZ = calc_chi(X, p_dict)
	Chi = [ChiPlus, ChiMinus, ChiZ]
	
	# Complex refractive index
	nPlus = sqrt(1.0+ChiPlus) #Complex refractive index driving sigma plus transitions
	nMinus = sqrt(1.0+ChiMinus) #Complex refractive index driving sigma minus transitions
	nZ = sqrt(1.0+ChiZ) # Complex index driving pi transitions

	# convert (if necessary) detuning axis X to np array
	if type(X) in (int, float, long):
		X = np.array([X])
	else:
		X = np.array(X)
	
	# Calculate E_field
	E_out, R = get_Efield(X, E_in, Chi, p_dict)
	#print 'Output E field (Z): \n', E_out[2]
	

	## Apply Jones matrices
	
	# Transmission - total intensity - just E_out**2 / E_in**2
	E_in = np.array(E_in)
	if E_in.shape == (3,):
			E_in = np.array([np.ones(len(X))*E_in[0],np.ones(len(X))*E_in[1],np.ones(len(X))*E_in[2]])
	
	# normalised by input intensity
	I_in = (E_in * E_in.conjugate()).sum(axis=0)
	
	S0 = (E_out * E_out.conjugate()).sum(axis=0) / I_in
	
	Iz = (E_out[2] * E_out[2].conjugate()).real / I_in
	
	Transmission = S0
	
	
	## Some quantities from Faraday geometry don't make sense when B and k not aligned, but leave them here for historical reasons
	TransLeft = exp(-2.0*nPlus.imag*wavenumber*lcell)
	TransRight = exp(-2.0*nMinus.imag*wavenumber*lcell)
	
	# Faraday rotation angle (including incident linear polarisation angle)
	phiPlus = wavenumber*nPlus.real*lcell
	phiMinus = wavenumber*nMinus.real*lcell
	phi = (phiMinus-phiPlus)/2.0 
	##
	
	#Stokes parameters

	#S1#
	Ex = np.array(JM.HorizPol_xy * E_out[:2])
	Ix =  (Ex * Ex.conjugate()).sum(axis=0) / I_in
	Ey =  np.array(JM.VertPol_xy * E_out[:2])
	Iy =  (Ey * Ey.conjugate()).sum(axis=0) / I_in
	
	S1 = Ix - Iy
	
	#S2#
	E_P45 =  np.array(JM.LPol_P45_xy * E_out[:2])
	E_M45 =  np.array(JM.LPol_M45_xy * E_out[:2])
	I_P45 = (E_P45 * E_P45.conjugate()).sum(axis=0) / I_in
	I_M45 = (E_M45 * E_M45.conjugate()).sum(axis=0) / I_in
	
	S2 = I_P45 - I_M45
	
	#S3#
	# change to circular basis
	E_out_lrz = BC.xyz_to_lrz(E_out)
	El =  np.array(JM.CPol_L_lr * E_out_lrz[:2])
	Er =  np.array(JM.CPol_R_lr * E_out_lrz[:2])
	Il = (El * El.conjugate()).sum(axis=0) / I_in
	Ir = (Er * Er.conjugate()).sum(axis=0) / I_in
	
	S3 = Ir - Il
	
	Ir = Ir.real
	Il = Il.real
	Ix = Ix.real
	Iy = Iy.real

	## (Real part) refractive indices
	#nMinus = nPlus.real
	#nPlus = nMinus.real

	## Absorption coefficients - again not a physically relevant quantity anymore since propagation is not as simple as k * Im(Chi) * L in a non-Faraday geometry
	alphaPlus = 2.0*nMinus.imag*wavenumber
	alphaMinus = 2.0*nPlus.imag*wavenumber
	alphaZ = 2.0*nZ.imag*wavenumber

	# Refractive/Group indices for left/right/z also no longer make any sense
	#d = (array(X)-shift) #Linear detuning
	#dnWRTv = derivative(d,nMinus.real)
	#GIPlus = nMinus.real + (X + transition.v0*1.0e-6)*dnWRTv
	#dnWRTv = derivative(d,nPlus.real)
	#GIMinus = nPlus.real + (X + transition.v0*1.0e-6)*dnWRTv
		
	if (outputs == None) or ('All' in outputs):
		# Default - return 'all' outputs (as used by GUI)
		return S0.real,S1.real,S2.real,S3.real,Ix.real,Iy.real,I_P45.real,I_M45.real,alphaPlus,alphaMinus,alphaZ
	else:
	# Return the variable names mentioned in the outputs list of strings
		# the strings in outputs must exactly match the local variable names here!
		return [locals()[output_str] for output_str in outputs]
	
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
	ChiPlus				Complex susceptibility of left-circularly polarised light \n\
	ChiMinus				Complex susceptibility of right-circularly polarised light \n\
	nPlus					Complex Refractive index of left-circularly polarised light \n\
	nMinus				Complex Refractive index of right-circularly polarised light \n\
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
### Normal Faraday spectrum
	import time
	d = np.arange(-10000,10000,10)
	#Voigt
	p_dict = {'Bfield':300,'rb85frac':1,'Btheta':0,'lcell':75e-3,'T':58,'Dline':'D2','Elem':'Cs'}
	
	#timing:
	st = time.clock()
	TF = get_spectra2(d,[1,0,0],p_dict,outputs=['Iy'])
	et = time.clock() - st
	print 'E-field - Elapsed time (s):', et

	#check vs old elecsus
	from elecsus.libs import spectra as old_spec
	
	st = time.clock()
	TF_old = old_spec.get_spectra(d,p_dict,outputs=['Iy'])
	et = time.clock() - st
	print 'Old elecsus - Elapsed time (s):', et
	
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
### 2. Fig 4/5 of General.... paper
### Voigt Filter
	d = np.linspace(-65000,65000,1500)
	#Voigt


## 700 G, 84 C, Cs, 75mm

	#p_dict = {'Bfield':700,'rb85frac':1,'Btheta':90*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	p_dict = {'Bfield':10000,'rb85frac':1,'Btheta':np.pi/2,'Bphi':90*np.pi/180,'lcell':5e-3,'T':124,'Dline':'D2','Elem':'Rb'}
	pol = 1./sqrt(2)*np.array([1.0,1.0,0.0])
	TVx = get_spectra2(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0','Iz'])
	
	fig2 = plt.figure()
	ax1a = fig2.add_subplot(411)
	ax2a = fig2.add_subplot(412,sharex=ax1a)
	ax3a = fig2.add_subplot(413,sharex=ax1a)
	ax4a = fig2.add_subplot(414,sharex=ax1a)
	
	ax1a.plot(d,TVx[0],'r',lw=2,label=r'$I_{-45}$')
	ax2a.plot(d,TVx[1],'b',lw=2,label=r'$I_{+45}$')
	ax3a.plot(d,TVx[2],'r',lw=2,label=r'$I_x$')
	ax4a.plot(d,TVx[3],'b',lw=2,label=r'$I_y$')
	ax4a.plot(d,TVx[0]+TVx[1],'r:',lw=3.5,label=r'$I_{+45}+I_{-45}$')
	ax4a.plot(d,TVx[2]+TVx[3],'k:',lw=2.5,label=r'$I_x + I_y$')
	ax4a.plot(d,TVx[4],'g--',lw=1.5,label='$S_0$')
#	ax4a.plot(d,TVx[5],'c--',lw=2.5,label='$I_z$')
	
	
	ax4a.set_xlabel('Detuning (MHz)')
	ax1a.set_ylabel('I -45')
	ax2a.set_ylabel('I +45')
	ax3a.set_ylabel('Ix')
	ax4a.set_ylabel('Iy')
	
	ax4a.set_xlim(d[0],d[-1]+3000)
	ax4a.legend(loc=0)
	
	plt.show()

def test3():
### 3. Fig 7 of General.... paper
### Arbitrary Filter - non-optimised
	d = np.linspace(-15000,15000,300)
	#Voigt
	#p_dict = {'Bfield':700,'rb85frac':1,'Btheta':90*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	p_dict = {'Bfield':500,'rb85frac':1,'Btheta':87*np.pi/180,'Bphi':00*np.pi/180,'lcell':75e-3,'T':100,'Dline':'D2','Elem':'Cs'}
	pol = 1./sqrt(2)*np.array([1.0,0.0,0.0])
	TVx = get_spectra2(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0','Iz'])
	
	fig2 = plt.figure()
	ax1a = fig2.add_subplot(411)
	ax2a = fig2.add_subplot(412,sharex=ax1a)
	ax3a = fig2.add_subplot(413,sharex=ax1a)
	ax4a = fig2.add_subplot(414,sharex=ax1a)
	
	ax1a.plot(d,TVx[0],'r',lw=2,label=r'$I_{-45}$')
	ax2a.plot(d,TVx[1],'b',lw=2,label=r'$I_{+45}$')
	ax3a.plot(d,TVx[2],'r',lw=2,label=r'$I_x$')
	ax4a.plot(d,TVx[3],'b',lw=2,label=r'$I_y$')
	ax4a.plot(d,TVx[0]+TVx[1],'r:',lw=3.5,label=r'$I_{+45}+I_{-45}$')
	ax4a.plot(d,TVx[2]+TVx[3],'k:',lw=2.5,label=r'$I_x + I_y$')
	ax4a.plot(d,TVx[4],'g--',lw=1.5,label='$S_0$')
#	ax4a.plot(d,TVx[5],'c--',lw=2.5,label='$I_z$')
	
	
	ax4a.set_xlabel('Detuning (MHz)')
	ax1a.set_ylabel('I -45')
	ax2a.set_ylabel('I +45')
	ax3a.set_ylabel('Ix')
	ax4a.set_ylabel('Iy')
	
	ax4a.set_xlim(d[0],d[-1]+3000)
	ax4a.legend(loc=0)
	
	plt.show()

def test4():
### 4. Fig 8 of General.... paper
### Arbitrary Filter - optimised
	d = np.linspace(-15000,15000,300)
	#Voigt
	#p_dict = {'Bfield':700,'rb85frac':1,'Btheta':90*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	p_dict = {'Bfield':1000,'rb85frac':1,'Btheta':88*np.pi/180,'Bphi':00*np.pi/180,'lcell':75e-3,'T':93,'Dline':'D2','Elem':'Cs'}
	pol = 1./sqrt(2)*np.array([1.0,0.0,0.0])
	TVx = get_spectra2(d,pol,p_dict,outputs=['I_M45','I_P45','Ix','Iy','S0','Iz'])
	
	fig2 = plt.figure()
	ax1a = fig2.add_subplot(411)
	ax2a = fig2.add_subplot(412,sharex=ax1a)
	ax3a = fig2.add_subplot(413,sharex=ax1a)
	ax4a = fig2.add_subplot(414,sharex=ax1a)
	
	ax1a.plot(d,TVx[0],'r',lw=2,label=r'$I_{-45}$')
	ax2a.plot(d,TVx[1],'b',lw=2,label=r'$I_{+45}$')
	ax3a.plot(d,TVx[2],'r',lw=2,label=r'$I_x$')
	ax4a.plot(d,TVx[3],'b',lw=2,label=r'$I_y$')
	ax4a.plot(d,TVx[0]+TVx[1],'r:',lw=3.5,label=r'$I_{+45}+I_{-45}$')
	ax4a.plot(d,TVx[2]+TVx[3],'k:',lw=2.5,label=r'$I_x + I_y$')
	ax4a.plot(d,TVx[4],'g--',lw=1.5,label='$S_0$')
#	ax4a.plot(d,TVx[5],'c--',lw=2.5,label='$I_z$')
	
	
	ax4a.set_xlabel('Detuning (MHz)')
	ax1a.set_ylabel('I -45')
	ax2a.set_ylabel('I +45')
	ax3a.set_ylabel('Ix')
	ax4a.set_ylabel('Iy')
	
	ax4a.set_xlim(d[0],d[-1]+3000)
	ax4a.legend(loc=0)
	
	plt.show()

def test5():
### 5. Durham Lab tests
	d = np.linspace(-15000,15000,1500)
	
	BFIELD = 500
	TEMP = 60
	LCELL = 2e-3
	
	ELEM = 'Rb'
	DLINE = 'D2'
	RB85FRAC = 99.0 # Nat abundance?
	
	# Voigt, H pol
	pol = 1./sqrt(2)*np.array([1.0,0.0,0.0])
	p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':90*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM}
	[T_H] = get_spectra2(d,pol,p_dict,outputs=['S0'])
		
	# Voigt, V-pol
	pol = 1./sqrt(2)*np.array([0.0,1.0,0.0])
	p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':90*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM}
	[T_V] = get_spectra2(d,pol,p_dict,outputs=['S0'])

	# Voigt, zero-field
	p_dict = {'Bfield':0,'rb85frac':RB85FRAC,'Btheta':90.*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM}
	[T_0] = get_spectra2(d,pol,p_dict,outputs=['S0'])
	
	## Bastard offspring of Voigt and Faraday
	#p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':45*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM}
	#[T_45Deg] = get_spectra2(d,pol,p_dict,outputs=['S0'])
	
	
	fig2 = plt.figure()
	ax1a = fig2.add_subplot(311)
	ax2a = fig2.add_subplot(312,sharex=ax1a)
	ax3a = fig2.add_subplot(313,sharex=ax1a)
	
	ax1a.plot(d,T_V,'r:',lw=1,label=r'$\theta=0$')
	ax1a.plot(d,T_H,'b:',lw=1,label=r'$\theta=0$')
	ax1a.plot(d,T_0,'k',lw=2,label=r'$\theta=0$')
	
	ax2a.plot(d,T_V,'k:',lw=1,label=r'$\theta=0$')
	ax2a.plot(d,T_0,'r:',lw=1,label=r'$\theta=0$')
	ax2a.plot(d,T_H,'b',lw=2,label=r'$\theta=0$')

	ax3a.plot(d,T_0,'k:',lw=1,label=r'$\theta=0$')
	ax3a.plot(d,T_H,'b:',lw=1,label=r'$\theta=0$')
	ax3a.plot(d,T_V,'r',lw=2,label=r'$\theta=0$')

	#	ax4a.plot(d,TVx[5],'c--',lw=2.5,label='$I_z$')
	
	
	ax3a.set_xlabel('Detuning (MHz)')
	ax1a.set_ylabel('Zero-field')
	ax2a.set_ylabel('H-polarised')
	ax3a.set_ylabel('V-polarised')
	
	ax3a.set_xlim(d[0],d[-1])
	#ax4a.legend(loc=0)
	
	plt.show()

def test6():
### 5. Boltzmann Factor tests
	d = np.linspace(-15,15,10000)*1e3
	
	BFIELD = 15500
	TEMP = 130
	LCELL = 1e-3
	
	ELEM = 'Rb'
	DLINE = 'D2'
	RB85FRAC = 1.0 # Nat abundance?
	
	# Voigt, 45 degree linear pol
	
	#Boltzmann ON
	pol = 1./sqrt(2)*np.array([1.0,1.0,0.0])
	p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':90*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM,'BoltzmannFactor':True}
	[S0_on,S2_on,S3_on,I_P45_on,I_M45_on] = get_spectra2(d,pol,p_dict,outputs=['S0','S2','S3','I_P45','I_M45'])
		
	#Boltzmann OFF
	pol = 1./sqrt(2)*np.array([1.0,1.0,0.0])
	p_dict = {'Bfield':BFIELD,'rb85frac':RB85FRAC,'Btheta':90*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':TEMP,'Dline':DLINE,'Elem':ELEM,'BoltzmannFactor':False}
	[S0_off,S2_off,S3_off,I_P45_off,I_M45_off] = get_spectra2(d,pol,p_dict,outputs=['S0','S2','S3','I_P45','I_M45'])
	
	fig2 = plt.figure()
	ax1a = fig2.add_subplot(311)
	ax2a = fig2.add_subplot(312,sharex=ax1a)
	ax3a = fig2.add_subplot(313,sharex=ax1a)
	
	ax1a.plot(d/1e3,S0_on,'r',lw=2,label=r'$\theta=0$')
	ax1a.plot(d/1e3,S0_off,'b--',lw=2,label=r'$\theta=0$')
	ax2a.plot(d/1e3,S2_on,'r',lw=2,label=r'$\theta=0$')
	ax2a.plot(d/1e3,S2_off,'b--',lw=2,label=r'$\theta=0$')
	ax3a.plot(d/1e3,S3_on,'r',lw=2,label=r'ON')
	ax3a.plot(d/1e3,S3_off,'b--',lw=2,label=r'OFF')
		
	ax3a.set_xlabel('Detuning (GHz)')
	ax1a.set_ylabel('S0')
	ax2a.set_ylabel('S2')
	ax3a.set_ylabel('S3')
	
	ax3a.set_xlim(d[0]/1e3,d[-1]/1e3)
	ax3a.legend(loc=0)
	
	plt.show()
	
if __name__ == '__main__':
	#test1()
	#test2()
	#test3()
	#test4()
	test6()