# Copyright 2014-2026 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams, I. G. Hughes and S. A. Wrathmall.

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

Last updated 2018-07-12 MAZ
"""
# py 2.7 compatibility

from numpy import zeros,sqrt,pi,dot,exp,sin,cos,array,amax,arange,concatenate,argmin,matmul
import numpy as np
from scipy.special import wofz
from scipy.interpolate import interp1d

import os
import sys
import inspect
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from FundamentalConstants import *
from numberDensityEqs import *

import EigenSystem as ES
import AtomConstants as AC
import RotationMatrices as RM
import BasisChanger as BC
import JonesMatrices as JM
import solve_dielectric as SD

# direct testing

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
							'BoltzmannFactor':True,
							# [ES4] non-orthogonal propagation switch and window refractive index (Briscoe et al. 2025)
							'use_non_orthogonal':False, 'n_window':None}

def FreqStren(groundLevels,excitedLevels,groundDim,
			  excitedDim,Dline,hand,BoltzmannFactor=True,T=293.16):
	""" 
	Calculate transition frequencies and strengths by taking dot
	products of relevant parts of the ground / excited state eigenvectors 
	"""
	
	transitionFrequency = zeros(groundDim*2*groundDim) #Initialise lists
	transitionStrength = zeros(groundDim*2*groundDim)
	transNo = 0
	
	## Boltzmann factor: -- only really needed when energy splitting of ground state is large
	if BoltzmannFactor:
		groundEnergies = array(groundLevels)[:,0].real
		lowestEnergy = groundLevels[0][argmin(groundLevels,axis=0)[0]].real
		BoltzDist = exp(-(groundEnergies-lowestEnergy)*h*1.e6/(kB*T)) #Not normalised
		#print BoltzDist
		BoltzDist = BoltzDist/BoltzDist.sum() # Normalised

	# select correct columns of matrix corresponding to delta mL = +1, 0, -1.
	# (Block order verified from the lz-operator basis: the excited manifold is
	# ordered mL = +1, 0, -1, so 'Right' picks the delta m = +1 (sigma-plus)
	# manifold and 'Left' the delta m = -1 (sigma-minus) manifold. The 'Left'/
	# 'Right' hand labels are the circular-polarisation names, NOT the sigma
	# labels -- see calc_chi and the naming caveat in output_list().)
	if hand=='Right':          # delta m = +1  ->  sigma-plus manifold
		bottom = 1
		top = groundDim + 1
	elif hand=='Z':            # delta m =  0  ->  pi
		bottom = groundDim + 1
		top = 2*groundDim + 1
	elif hand=='Left':         # delta m = -1  ->  sigma-minus manifold
		bottom = 2*groundDim+1
		top = excitedDim+1

	# select correct rows of matrix corresponding to D1 / D2 lines
	if Dline=='D1':
		interatorList = list(range(groundDim))
	elif Dline=='D2':
		interatorList = list(range(groundDim,excitedDim))\
	
	# find difference in energies and do dot product between (all) ground states
	# 	and selected parts of excited state matrix
	for gg in range(groundDim):
		for ee in interatorList:
			cleb = dot(groundLevels[gg][1:],excitedLevels[ee][bottom:top]).real
			cleb2 = cleb*cleb
			if cleb2 > 0.0005: #If negligable don't calculate.
				transitionFrequency[transNo] = int((-groundLevels[gg][0].real
											  +excitedLevels[ee][0].real))
				# We choose to perform the ground manifold reduction (see
				# equation (4) in manual) here for convenience.
				
				## Boltzmann factor:
				if BoltzmannFactor: 
					transitionStrength[transNo] = 1./3 * cleb2 * BoltzDist[gg]
				else:
					transitionStrength[transNo] = 1./3 * cleb2 * 1./groundDim
				
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
	for line in range(ltransno+1):
		xc = lenergy[line]
		lab += lstrength[line]*f_ab(2.0*pi*(d-xc)*1.0e6)
		ldisp += lstrength[line]*f_disp(2.0*pi*(d-xc)*1.0e6)
	rab = zeros(xpts)
	rdisp = zeros(xpts)
	for line in range(rtransno+1):
		xc = renergy[line]
		rab += rstrength[line]*f_ab(2.0*pi*(d-xc)*1.0e6)
		rdisp += rstrength[line]*f_disp(2.0*pi*(d-xc)*1.0e6)
	zab = zeros(xpts)
	zdisp = zeros(xpts)
	for line in range(ztransno+1):
		xc = zenergy[line]
		zab += zstrength[line]*f_ab(2.0*pi*(d-xc)*1.0e6)
		zdisp += zstrength[line]*f_disp(2.0*pi*(d-xc)*1.0e6)
	return lab, ldisp, rab, rdisp, zab, zdisp

def calc_chi(X, p_dict,verbose=False):			   
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
			GammaBuf   	float			MHz		Extra phenomenological Lorentzian broadening (MHz).
															Added on top of the natural linewidth. Use for residual
															broadening not captured by the BufferGas model (e.g.
															impurities, laser linewidth contribution).
															IGNORED (silently overridden) if BufferGas is set and
															BufferPressure > 0 — the Pitz model value is used instead.
															Set to 0 and use BufferGas parameters for known buffer gases.
			shift	  			float			MHz		Global frequency shift of all resonance lines (MHz).
															IGNORED (silently overridden) if BufferGas is set and
															BufferPressure > 0 — the Pitz pressure-shift is used instead.

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
	if 'Elem' in list(p_dict.keys()):
		Elem = p_dict['Elem']
	else:
		Elem = p_dict_defaults['Elem']
	if 'Dline' in list(p_dict.keys()):
		Dline = p_dict['Dline']
	else:
		Dline = p_dict_defaults['Dline']
	if 'T' in list(p_dict.keys()):
		T = p_dict['T']
	else:
		T = p_dict_defaults['T']
	if 'Bfield' in list(p_dict.keys()):
		Bfield = p_dict['Bfield']
	else:
		Bfield = p_dict_defaults['Bfield']
	if 'GammaBuf' in list(p_dict.keys()):
		GammaBuf = p_dict['GammaBuf']
	else:
		GammaBuf = p_dict_defaults['GammaBuf']
	if 'shift' in list(p_dict.keys()):
		shift = p_dict['shift']
	else:
		shift = p_dict_defaults['shift']

	# [ES4 NEW] Dynamic buffer gas model — overrides static GammaBuf/shift when 'BufferGas' is set.
	# Computes temperature-dependent broadening and pressure shift via the Pitz 2014 model
	# (see elecsus/libs/BufferGasData.py).  To use the static model, omit 'BufferGas'
	# from the parameter dictionary or set it to None.
	if p_dict.get('BufferGas') is not None:
		import warnings
		from elecsus.libs.BufferGasData import get_buffer_broadening_and_shift
		_bg   = p_dict['BufferGas']
		_bg_P = p_dict.get('BufferPressure', 0.0)
		_bg_Tf = p_dict.get('BufferTemp_fill', 293.15)
		if _bg_P <= 0.0:
			warnings.warn(
				f"BufferGas='{_bg}' is set but BufferPressure={_bg_P} Torr. "
				"Set 'BufferPressure' to the fill pressure in Torr.", stacklevel=2)
		else:
			_gamma_phys, _shift_phys = get_buffer_broadening_and_shift(
				Elem, Dline, _bg, _bg_P, _bg_Tf, T + 273.15)
			if _gamma_phys is not None:
				if GammaBuf != 0.0:
					warnings.warn(
						f"BufferGas model active: GammaBuf={GammaBuf} MHz is overridden "
						f"by the Pitz model value ({_gamma_phys:.2f} MHz). "
						"Remove GammaBuf from your parameter dict to suppress this warning.",
						stacklevel=2)
				if shift != 0.0:
					warnings.warn(
						f"BufferGas model active: shift={shift} MHz is overridden "
						f"by the Pitz model pressure-shift ({_shift_phys:.2f} MHz). "
						"Remove shift from your parameter dict to suppress this warning.",
						stacklevel=2)
				GammaBuf = _gamma_phys
				shift    = _shift_phys
			else:
				warnings.warn(
					f"No buffer gas data for ({Elem}, {Dline}, '{_bg}'). "
					"Falling back to static GammaBuf/shift.", stacklevel=2)

	if 'Constrain' in list(p_dict.keys()):
		Constrain = p_dict['Constrain']
	else:
		Constrain = p_dict_defaults['Constrain']
	if 'rb85frac' in list(p_dict.keys()):
		rb85frac = p_dict['rb85frac']
	else:
		rb85frac = p_dict_defaults['rb85frac']
	if 'DoppTemp' in list(p_dict.keys()):
		DoppTemp = p_dict['DoppTemp']
	else:
		DoppTemp = p_dict_defaults['DoppTemp']
	if 'K40frac' in list(p_dict.keys()):
		K40frac = p_dict['K40frac']
	else:
		K40frac = p_dict_defaults['K40frac']
	if 'K41frac' in list(p_dict.keys()):
		K41frac = p_dict['K41frac']
	else:
		K41frac = p_dict_defaults['K41frac']
	if 'BoltzmannFactor' in list(p_dict.keys()):
		BoltzmannFactor =  p_dict['BoltzmannFactor']
	else:
		BoltzmannFactor =  p_dict_defaults['BoltzmannFactor']
	
	
	if verbose: print(('Temperature: ', T, '\tBfield: ', Bfield))
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
			
			# Rb-85 allowed transitions for light driving pi
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
			lenergy87, lstrength87, ltransno87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Left',BoltzmannFactor,T+273.16)

			# Rb-87 allowed transitions for light driving sigma plus
			renergy87, rstrength87, rtransno87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Right',BoltzmannFactor,T+273.16)

			# Rb-87 allowed transitions for light driving sigma plus
			zenergy87, zstrength87, ztransno87 = FreqStren(
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
		AllEnergyLevels = concatenate((lenergy,renergy,zenergy))

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
			AllEnergyLevels = concatenate((lenergy39,renergy39,zenergy39))
		if K40frac!=0.0 and K39frac!=0.0:
			AllEnergyLevels = concatenate((AllEnergyLevels,lenergy40,renergy40,zenergy40))
		elif K40frac!=0.0 and K39frac==0.0:
			AllEnergyLevels = concatenate((lenergy40,renergy40,zenergy40))
		if K41frac!=0.0 and (K39frac!=0.0 or K40frac!=0.0):
			AllEnergyLevels = concatenate((AllEnergyLevels,lenergy41,renergy41,zenergy41))
		elif K41frac!=0.0 and (K39frac==0.0 and K40frac==0.0):
			AllEnergyLevels = concatenate((lenergy41,renergy41,zenergy41))

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

# [ES4 BEGIN] propagate_non_orthogonal — biorthogonal mode propagation (Briscoe et al. 2025)
def propagate_non_orthogonal(E_xz, RM_ary, n1, n2, wavenumber, lcell, n_window=None):
	"""
	Non-orthogonal mode propagation following Briscoe et al.,
	Opt. Express 33, 32099-32114 (2025).

	Solves E_inc = c1*u1 + c2*u2 via linear system instead of assuming
	orthogonal eigenmodes.

	Arguments:
		E_xz:       (N, 3, 1) input E-field vectors (rotated frame)
		RM_ary:     (3, 3, N) eigenvector array from solve_diel (eigenvectors as rows)
		n1, n2:     (N,) complex refractive indices for the two eigenmodes
		wavenumber: scalar, wavevector magnitude
		lcell:      scalar, cell length in metres
		n_window:   refractive index of cell windows (default None — Fresnel window correction disabled; set e.g. 1.51 for borosilicate glass).

		            When set, a 4-interface Fresnel model is applied at normal
		            incidence, assuming outer medium is air (n=1):

		              air ──[t_outer_in]──> glass ──[t_inner_in]──> vapour
		                                          <── propagate ──>
		              air <──[t_outer_out]── glass <──[t_inner_out]── vapour

		            Outer interfaces (air ↔ glass): scalar, frequency-independent.
		              t_outer_in  = 2 / (1 + n_window)
		              t_outer_out = 2·n_window / (n_window + 1)

		            Inner interfaces (glass ↔ vapour): per eigenmode,
		            frequency-dependent near resonance.
		              t_inner_in  = 2·n_window / (n_window + n_i)
		              t_inner_out = 2·n_i      / (n_i + n_window)

		            The combined window factor per mode at normal incidence is:
		              T_window,i = |t_outer_in · t_inner_in · t_inner_out · t_outer_out|²
		            which equals the physical intensity transmission (T = 4n_a n_b/(n_a+n_b)²
		            at each interface) because the input and output media are both air.

		            For n_window=1.51 far from resonance (n_i ≈ 1):
		              T_window ≈ 0.959⁴ ≈ 0.845  (~15% total reflection loss)

		            The outer surfaces contribute a constant factor that cancels in
		            spectra normalised to a reference measurement, but is required
		            for rigorous absolute transmission modelling.

		            Set to None to disable all Fresnel corrections — this reproduces
		            the bare Briscoe et al. formalism (Eq. 9) and is appropriate for
		            theory-only comparisons or when windows are not modelled.

		            Note: the legacy orthogonal propagation path (use_non_orthogonal=False)
		            contains no window correction. This inconsistency is documented in
		            dev-notes/window_fresnel_correction.md.

	Returns:
		E_out: (N, 3, 1) output E-field vectors
	"""
	N = len(n1)

	# Extract transverse (x,y) components of eigenvectors
	# RM_ary[0,:,i] = eigenvector 1 at point i
	# RM_ary[1,:,i] = eigenvector 2 at point i
	u1_xy = RM_ary[0, :2, :].T  # shape (N, 2)
	u2_xy = RM_ary[1, :2, :].T  # shape (N, 2)
	u1_z = RM_ary[0, 2, :]       # shape (N,)
	u2_z = RM_ary[1, 2, :]       # shape (N,)

	# Build 2x2 basis matrices: M[i] = [[u1_x, u2_x], [u1_y, u2_y]]
	M = np.zeros((N, 2, 2), dtype='complex')
	M[:, 0, 0] = u1_xy[:, 0]
	M[:, 1, 0] = u1_xy[:, 1]
	M[:, 0, 1] = u2_xy[:, 0]
	M[:, 1, 1] = u2_xy[:, 1]

	# Input E-field transverse components
	E_in_xy = E_xz[:, :2, 0]  # shape (N, 2)

	# Solve for mode amplitudes: M @ c = E_in_xy  (Briscoe et al. Eq. 11)
	try:
		c = np.linalg.solve(M, E_in_xy[..., np.newaxis])[..., 0]
	except np.linalg.LinAlgError:
		# Fallback for singular matrices near resonances
		c = np.zeros((N, 2), dtype='complex')
		for i in range(N):
			c[i] = np.linalg.lstsq(M[i], E_in_xy[i], rcond=None)[0]

	if n_window is not None:
		# Outer entrance: air (n=1) → glass — scalar, same for both modes
		t_outer_in = 2.0 / (1.0 + n_window)
		# Inner entrance: glass → vapour mode i — per-mode, frequency-dependent
		t_inner_in = np.column_stack([
			2.0 * n_window / (n_window + n1),
			2.0 * n_window / (n_window + n2),
		])
		c *= t_outer_in * t_inner_in   # scalar broadcasts over (N, 2) array

	# Propagate through vapour (Briscoe et al. Eq. 9: t(n_i) = exp(i n_i k_0 L))
	phase = np.column_stack([
		np.exp(1j * wavenumber * n1 * lcell),
		np.exp(1j * wavenumber * n2 * lcell),
	])
	c *= phase

	if n_window is not None:
		# Inner exit: vapour mode i → glass — per-mode, frequency-dependent
		t_inner_out = np.column_stack([
			2.0 * n1 / (n1 + n_window),
			2.0 * n2 / (n2 + n_window),
		])
		# Outer exit: glass → air — scalar, same for both modes
		t_outer_out = 2.0 * n_window / (n_window + 1.0)
		c *= t_inner_out * t_outer_out   # scalar broadcasts over (N, 2) array

	# Reconstruct output E-field
	E_out = np.zeros((N, 3, 1), dtype='complex')
	E_out[:, 0, 0] = c[:, 0] * u1_xy[:, 0] + c[:, 1] * u2_xy[:, 0]  # E_x
	E_out[:, 1, 0] = c[:, 0] * u1_xy[:, 1] + c[:, 1] * u2_xy[:, 1]  # E_y
	E_out[:, 2, 0] = c[:, 0] * u1_z + c[:, 1] * u2_z                  # E_z

	return E_out
# [ES4 END] propagate_non_orthogonal


# [ES4 BEGIN] Poynting-flux transmission for the non-orthogonal path (issue #28)
def poynting_transmission_non_orthogonal(E_xz, RM_ary, n1, n2, wavenumber, lcell):
	"""
	Physical (Poynting-flux) transmission through the vapour for the
	non-orthogonal propagation path:

	    S0_P = <S_z(L)> / <S_z(0)>,
	    <S_z(z)> ∝ Re[ Σ_jk c_j(z) c_k*(z) n_k* (m_j · m_k*) ]

	where c_j(z) = c_j exp(i n_j k0 z) and m_j are the transverse eigenvector
	components. This is the time-averaged z-component of the Poynting vector
	of the two-mode field, using H_k = (n_k/Z0) ẑ×E_k per mode.

	Rationale (issue #28): |E|² is not proportional to power flux inside an
	anisotropic medium. Near mode degeneracies the biorthogonal amplitudes
	c_j grow without bound and |E_out|²/|E_in|² ('S0') can exceed unity by
	orders of magnitude even though no energy is created. The Poynting ratio
	is the physically transported power and satisfies S0_P ≤ 1 for a passive
	medium (∂<S_z>/∂z = −absorption ≤ 0 for genuine Maxwell solutions, which
	the two-mode superposition is).

	Window Fresnel factors are NOT included: this is the flux transmission of
	the bare vapour column (compare n_window=None), so that the passivity
	bound applies exactly.
	"""
	N = len(n1)
	u1_xy = RM_ary[0, :2, :].T   # (N, 2)
	u2_xy = RM_ary[1, :2, :].T

	M = np.zeros((N, 2, 2), dtype='complex')
	M[:, 0, 0] = u1_xy[:, 0]
	M[:, 1, 0] = u1_xy[:, 1]
	M[:, 0, 1] = u2_xy[:, 0]
	M[:, 1, 1] = u2_xy[:, 1]

	E_in_xy = E_xz[:, :2, 0]
	try:
		c0 = np.linalg.solve(M, E_in_xy[..., np.newaxis])[..., 0]
	except np.linalg.LinAlgError:
		c0 = np.zeros((N, 2), dtype='complex')
		for i in range(N):
			c0[i] = np.linalg.lstsq(M[i], E_in_xy[i], rcond=None)[0]

	cL = c0 * np.column_stack([
		np.exp(1j * wavenumber * n1 * lcell),
		np.exp(1j * wavenumber * n2 * lcell),
	])

	n_arr = np.column_stack([n1, n2])                 # (N, 2)
	m = np.stack([u1_xy, u2_xy], axis=1)              # (N, mode j, xy)
	G = np.einsum('njx,nkx->njk', m, m.conj())        # m_j · m_k*   (N, 2, 2)

	def _flux(c):
		return np.real(np.einsum('nj,nk,nk,njk->n', c, c.conj(), n_arr.conj(), G))

	return _flux(cL) / _flux(c0)
# [ES4 END] poynting_transmission_non_orthogonal


def get_Efield(X, E_in, Chi, p_dict, verbose=False, return_poynting=False):
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
		print('Electric field input:')
		print(E_in)
	
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
	if 'lcell' in list(p_dict.keys()):
		lcell = p_dict['lcell']
	else:
		lcell = p_dict_defaults['lcell']
		
	if 'Elem' in list(p_dict.keys()):
		Elem = p_dict['Elem']
	else:
		Elem = p_dict_defaults['Elem']
	if 'Dline' in list(p_dict.keys()):
		Dline = p_dict['Dline']
	else:
		Dline = p_dict_defaults['Dline']
	
	## get magnetic field spherical coordinates
	# defaults to 0,0 i.e. B aligned with kvector of light (Faraday)
	if 'Bfield' in list(p_dict.keys()):
		Bfield = p_dict['Bfield']
	else:
		Bfield = p_dict_defaults['Bfield']
	if 'Btheta' in list(p_dict.keys()):
		Btheta = p_dict['Btheta']
	else:
		Btheta = p_dict_defaults['Btheta']
	if 'Bphi' in list(p_dict.keys()):
		Bphi = p_dict['Bphi']
	else:
		Bphi = p_dict_defaults['Bphi']
		
	# [ES4] read non-orthogonal propagation parameters from p_dict
	use_non_orthogonal = p_dict.get('use_non_orthogonal',
		p_dict_defaults.get('use_non_orthogonal', False))
	n_window = p_dict.get('n_window',
		p_dict_defaults.get('n_window', None))

	# direction of wavevector (for double-pass geometry simulations)
	if 'wavevector_dirn' in list(p_dict.keys()):
		wavevector_dirn = p_dict['wavevector_dirn']
	else:
		wavevector_dirn = 1
	
	if wavevector_dirn == -1:
		# light propagating backwards, but always calculate assuming light travelling in the +z direction,
		# therefore flip B-field 180 degrees around (can use either theta or phi for this)
		Btheta += np.pi
	
	# get wavenumber
	transition = AC.transitions[Elem+Dline]
	wavenumber = transition.wavevectorMagnitude * wavevector_dirn

	# get susceptibility (already calculated, input to this method)
	ChiPlus, ChiMinus, ChiZ = Chi

	# Rotate initial Electric field so that B field lies in x-z plane
	# (Effective polarisation rotation)
	E_xz = RM.rotate_around_z(E_in.T,Bphi)
		
	# Find eigen-vectors for propagation and create rotation matrix
	RM_ary, n1, n2 = SD.solve_diel(ChiPlus,ChiMinus,ChiZ,Btheta,Bfield)

	# Reshape E_xz for array operations
	E_xz = np.reshape(E_xz.T, (len(X),3,1))

	# [ES4 BEGIN] branch: non-orthogonal (Briscoe 2025) vs ES3 orthogonal propagation
	if use_non_orthogonal:
		# Non-orthogonal mode propagation (Briscoe et al. 2025)
		# Uses linear system solver instead of matrix inverse for eigenmode projection,
		# and includes Fresnel transmission at cell windows.
		E_out_xz = propagate_non_orthogonal(E_xz, RM_ary, n1, n2,
											 wavenumber, lcell, n_window)
		fast_E_out = RM.rotate_around_z2(
			np.transpose(E_out_xz,(0,2,1))[::,0], -Bphi)
		if return_poynting:
			# [ES4] Poynting-flux transmission (issue #28): physical power
			# ratio, bounded by 1 for the passive vapour, unlike |E|²-based
			# S0 near mode degeneracies.
			S0_P = poynting_transmission_non_orthogonal(
				E_xz, RM_ary, n1, n2, wavenumber, lcell)
			return fast_E_out.T[0], RM_ary.T[-1], S0_P
	else:
		if return_poynting:
			raise ValueError(
				"S0_P (Poynting-flux transmission) requires "
				"use_non_orthogonal=True: it is defined for the "
				"non-orthogonal mode decomposition, where |E|²-based S0 "
				"can exceed 1 near mode degeneracies (issue #28).")
		# Original orthogonal propagation
		# take complex conjugate of rotation matrix
		RM_ary = RM_ary.conjugate()

		# propagation matrix
		PropMat = np.array(
					[	[exp(1.j*n1*wavenumber*lcell),np.zeros(len(n1)),np.zeros(len(n1))],
						[np.zeros(len(n1)),exp(1.j*n2*wavenumber*lcell),np.zeros(len(n1))],
						[np.zeros(len(n1)),np.zeros(len(n1)),np.ones(len(n1))]	])

		ary_of_RMI_ary = np.linalg.inv(np.transpose(RM_ary.T,(0,2,1)))
		fast_E_out_xz = matmul(matmul(matmul(ary_of_RMI_ary,PropMat.T),
							np.transpose(RM_ary.T,(0,2,1))),E_xz)
		fast_E_out = RM.rotate_around_z2(
			np.transpose(fast_E_out_xz,(0,2,1))[::,0], -Bphi)
	# [ES4 END] propagation branch

	## return electric field vector - can then use Jones matrices to do everything else
	return fast_E_out.T[0], RM_ary.T[-1]  # [ES4 MOD] np.matrix deprecated; RM_ary is already ndarray

def get_spectra(X, E_in, p_dict, outputs=None):
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
				GammaBuf   	float			MHz		Extra phenomenological Lorentzian broadening (MHz).
																Added on top of the natural linewidth. Use for residual
																broadening not captured by the BufferGas model.
																IGNORED if BufferGas is set and BufferPressure > 0.
				shift	  			float			MHz		Global frequency shift of all resonance lines (MHz).
																IGNORED if BufferGas is set and BufferPressure > 0
																(the Pitz pressure-shift is used instead).
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
		
		
		More examples are available in the /tests/ folder
	"""
	
	# get some parameters from p dictionary
	
	# need in try/except or equiv.
	if 'Elem' in list(p_dict.keys()):
		Elem = p_dict['Elem']
	else:
		Elem = p_dict_defaults['Elem']
	if 'Dline' in list(p_dict.keys()):
		Dline = p_dict['Dline']
	else:
		Dline = p_dict_defaults['Dline']
	if 'shift' in list(p_dict.keys()):
		shift = p_dict['shift']
	else:
		shift = p_dict_defaults['shift']
	if 'lcell' in list(p_dict.keys()):
		lcell = p_dict['lcell']
	else:
		lcell = p_dict_defaults['lcell']
	if 'theta0' in list(p_dict.keys()):
		theta0 = p_dict['theta0']
	else:
		theta0 = p_dict_defaults['theta0']
	if 'Pol' in list(p_dict.keys()):
		Pol = p_dict['Pol']
	else:
		Pol = p_dict_defaults['Pol']

	# get wavenumber
	transition = AC.transitions[Elem+Dline]

	wavenumber = transition.wavevectorMagnitude
	
	# Calculate Susceptibility
	ChiPlus, ChiMinus, ChiZ = calc_chi(X, p_dict)
	Chi = [ChiPlus, ChiMinus, ChiZ]
	
	# Complex refractive index.
	# NB the 'Plus'/'Minus' suffix here is the circular-basis label, not the
	# sigma label: ChiPlus is the sigma-MINUS (delta m=-1) manifold and ChiMinus
	# the sigma-PLUS (delta m=+1) manifold (verified from the lz basis ordering
	# and the high-field Zeeman shift; see output_list() naming caveat). The
	# alphaPlus/alphaMinus below cross over so they carry the physical-sigma label.
	nPlus = sqrt(1.0+ChiPlus)    # index of the [1,+i]/sqrt2 mode (sigma-minus manifold)
	nMinus = sqrt(1.0+ChiMinus)  # index of the [1,-i]/sqrt2 mode (sigma-plus manifold)
	nZ = sqrt(1.0+ChiZ) # Complex index driving pi transitions

	# convert (if necessary) detuning axis X to np array
	if type(X) in (int, float, int):
		X = np.array([X])
	else:
		X = np.array(X)


	
	# Calculate E_field
	# [ES4] S0_P (Poynting-flux transmission, issue #28) is computed inside
	# get_Efield where the mode data are available; request it only if asked.
	if outputs is not None and 'S0_P' in outputs:
		E_out, R, S0_P = get_Efield(X, E_in, Chi, p_dict, return_poynting=True)
	else:
		E_out, R = get_Efield(X, E_in, Chi, p_dict)
		S0_P = None
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
	# [ES4 MOD] S1 now fractional (normalised by S0); ES3 returned Ix - Iy (absolute intensity
	# difference, units of input intensity).  Stokes parameters are dimensionless by definition
	# and must satisfy S0² ≥ S1² + S2² + S3² (Born & Wolf §1.4).
	S1 = (Ix - Iy) / S0

	#S2#
	E_P45 =  np.array(JM.LPol_P45_xy * E_out[:2])
	E_M45 =  np.array(JM.LPol_M45_xy * E_out[:2])
	I_P45 = (E_P45 * E_P45.conjugate()).sum(axis=0) / I_in
	I_M45 = (E_M45 * E_M45.conjugate()).sum(axis=0) / I_in

	# [ES4 MOD] S2 now fractional (normalised by S0); ES3 returned I_P45 - I_M45 (absolute).
	# Same correction as S1 — fractional Stokes is the standard polarimetric convention.
	S2 = (I_P45 - I_M45) / S0

	#S3#
	# change to circular basis
	E_out_lrz = BC.xyz_to_lrz(E_out)
	El =  np.array(JM.CPol_L_lr * E_out_lrz[:2])
	Er =  np.array(JM.CPol_R_lr * E_out_lrz[:2])
	Il = (El * El.conjugate()).sum(axis=0) / I_in
	Ir = (Er * Er.conjugate()).sum(axis=0) / I_in

	# [ES4 MOD] S3 now fractional and sign-corrected.  ES3: S3 = Ir - Il (absolute, positive for
	# right-circular dominant).  Born & Wolf §1.4 defines S3 = I_l - I_r, i.e. positive S3
	# means left-circularly polarised light is dominant.  ES3 had the sign inverted.
	S3 = -1 * (Ir - Il) / S0

	Ir = Ir.real
	Il = Il.real
	Ix = Ix.real
	Iy = Iy.real

	# [ES4] Generalised polarisation rotation angle (psi) and degree of linear
	# polarisation (DoLP), from the fractional Stokes parameters.
	#
	# 'phi' above is a DISPERSION-ONLY quantity: (k L / 2)(Re nMinus - Re nPlus)
	# ignores the imaginary parts entirely, so it keeps reporting a "rotation
	# angle" even when strong circular dichroism near resonance has converted
	# the light to (nearly) circular polarisation and no linear alignment
	# exists to rotate.  psi is instead the polarisation-ellipse azimuth
	# (Born & Wolf §1.4) computed from the full complex field propagation, so
	# absorption and dispersion are both included.  In pure Faraday geometry
	# the circular modes are exact eigenmodes and the azimuth is set by the
	# phase difference alone, so psi == phi (mod pi) at ALL optical depths
	# there (verified to 1e-12): dichroism alters the ellipticity, not the
	# azimuth -- what it destroys is the linear component itself, which DoLP
	# tracks.  In oblique geometry phi is not the propagation phase difference
	# at all and psi is the correct angle.  psi is returned as the principal
	# value in (-pi/2, pi/2]; rotations beyond +-pi/2 wrap (unwrap with
	# np.unwrap(2*psi)/2 if tracking a continuous scan).
	#
	# DoLP = sqrt(S1^2 + S2^2) (fractional Stokes, so already /S0) is the
	# companion validity diagnostic: psi is physically meaningful only where
	# DoLP is appreciable.  DoLP -> 0 marks exactly the dichroic regime where
	# a rotation angle becomes unmeasurable.  For the coherent (fully
	# polarised) fields propagated here, DoLP^2 + S3^2 = 1.
	psi = 0.5 * np.arctan2(S2.real, S1.real)
	DoLP = np.sqrt(S1.real**2 + S2.real**2)

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
		
	# Valid outputs
	op = {'S0':S0, 'S1':S1, 'S2':S2, 'S3':S3, 'Ix':Ix, 'Iy':Iy, 'Il':Il, 'Ir':Ir,
				'I_P45':I_P45, 'I_M45':I_M45,
				'alphaPlus':alphaPlus, 'alphaMinus':alphaMinus, 'alphaZ':alphaZ,
				'E_out':E_out,
				'Chi':Chi, 'ChiPlus':ChiPlus, 'ChiMinus':ChiMinus, 'ChiZ':ChiZ,
				# [ES4] TransLeft/TransRight added to output dict (were computed but not accessible in ES3)
				'TransLeft': TransLeft, 'TransRight': TransRight,
				# [ES4] nPlus/nMinus/phi/phiPlus/phiMinus added to output dict (were
				# computed but not accessible in ES3, so requesting them raised
				# KeyError). Like TransLeft/alphaPlus above, these are circular-basis
				# quantities only physically meaningful in Faraday geometry (Btheta=0);
				# the code computes them unconditionally. See output_list() for the
				# Plus/Minus-vs-sigma naming caveat.
				'nPlus': nPlus, 'nMinus': nMinus,
				'phi': phi, 'phiPlus': phiPlus, 'phiMinus': phiMinus,
				# [ES4] Stokes-based rotation angle + linear-polarisation degree
				# (dichroism-aware; valid in any geometry — see comment above)
				'psi': psi, 'DoLP': DoLP,
				# [ES4] Poynting-flux transmission (issue #28); requires
				# use_non_orthogonal=True and excludes window Fresnel factors
				'S0_P': S0_P
			}
	
	if (outputs == None) or ('All' in outputs):
		# Default - return 'all' outputs (as used by GUI)
		return S0.real,S1.real,S2.real,S3.real,Ix.real,Iy.real,I_P45.real,I_M45.real,alphaPlus,alphaMinus,alphaZ
	else:
	# Return the variable names mentioned in the outputs list of strings
		# the strings in outputs must exactly match the local variable names here!
		return [op[output_str] for output_str in outputs]
	
def output_list():
	""" Helper method that returns a list of all possible variables that get_spectra can return """
	# NAMING CAVEAT (verified 2026-07-11): the historical 'Plus'/'Minus' suffix on
	# Chi/n/phi refers to the CIRCULAR-BASIS component, not the sigma transition:
	# ChiPlus/nPlus are built from the sigma-MINUS (delta m = -1) manifold and govern
	# [1,+i]/sqrt2 light; ChiMinus/nMinus are the sigma-PLUS (delta m = +1) manifold
	# and govern [1,-i]/sqrt2 light. alphaPlus/alphaMinus deliberately cross over
	# (alphaPlus = 2*nMinus.imag*k) so that they DO carry the physical-sigma label.
	# See Sec. 5 of the manual and the 'Complex refractive index' block above.
	tstr = " \
	All possible outputs from the get_spectra method: \n\n\
	Variable Name		Description \n \
	S0						Total transmission through the cell (Ix + Iy) \n\
	S1						Stokes parameter - Ix - Iy \n\
	S2						Stokes parameter - I_45 - I_-45 \n\
	S3						Stokes parameter - I- - I+ \n\
	TransLeft			Single-mode transmission exp(-2 Im(nPlus) k L) of the [1,+i]/sqrt2 (sigma-minus) mode \n\
	TransRight			Single-mode transmission exp(-2 Im(nMinus) k L) of the [1,-i]/sqrt2 (sigma-plus) mode \n\
	S0_P					Poynting-flux (physical power) transmission of the bare vapour column. \n\
							Requires use_non_orthogonal=True; excludes window Fresnel factors; \n\
							bounded by 1 for the passive medium even at mode degeneracies where \n\
							|E|^2-based S0 is not proportional to transported power (issue #28) \n\
	ChiPlus				Complex susceptibility of the sigma-minus (delta m=-1) manifold [governs [1,+i]/sqrt2 light] \n\
	ChiMinus				Complex susceptibility of the sigma-plus (delta m=+1) manifold [governs [1,-i]/sqrt2 light] \n\
	nPlus					Complex refractive index sqrt(1+ChiPlus) of the sigma-minus manifold [Faraday geometry only] \n\
	nMinus				Complex refractive index sqrt(1+ChiMinus) of the sigma-plus manifold [Faraday geometry only] \n\
	phiPlus				Single-pass phase k Re(nPlus) L of the sigma-minus mode [Faraday geometry only] \n\
	phiMinus				Single-pass phase k Re(nMinus) L of the sigma-plus mode [Faraday geometry only] \n\
	phi					Net Faraday rotation of linear polarisation, (phiMinus - phiPlus)/2 [Faraday geometry only]. \n\
								DISPERSION-ONLY: ignores circular dichroism, so it is not a measurable \n\
								rotation near resonance at high optical depth — use psi (and check DoLP) there \n\
	psi					Polarisation-ellipse azimuth 0.5*arctan2(S2, S1), from the full complex field \n\
								propagation (includes dichroism). Equals phi mod pi at weak dichroism; principal \n\
								value in (-pi/2, pi/2]. Only meaningful where DoLP is appreciable \n\
	DoLP					Degree of linear polarisation sqrt(S1^2 + S2^2) (fractional, in [0, 1]). \n\
								DoLP -> 0 flags where a linear rotation angle ceases to exist (dichroism \n\
								has made the light circular); DoLP^2 + S3^2 = 1 for the coherent fields here \n\
	Ix						Intensity of light transmitted through a linear polariser aligned with the x-axis \n\
	Iy						Intensity of light transmitted through a linear polariser aligned with the y-axis \n\
	Ir						Intensity of the [1,+i]/sqrt2 (circular 'R') component\n\
	Il						Intensity of the [1,-i]/sqrt2 (circular 'L') component\n\
	alphaPlus			Absorption coefficient 2 Im(nMinus) k of the sigma-plus manifold \n\
	alphaMinus			Absorption coefficient 2 Im(nPlus) k of the sigma-minus manifold \n\
	alphaZ				Absorption coefficient 2 Im(nZ) k of the pi (delta m=0) manifold \n\
	"
	print(tstr)