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
"""
# py 2.7 compatibility

from numpy import zeros,sqrt,pi,dot,exp,sin,cos,array,amax,arange,concatenate,argmin,matmul
import numpy as np
from scipy.special import wofz
from scipy.interpolate import interp1d


from .FundamentalConstants import *
from .numberDensityEqs import *

from . import EigenSystem as ES
from . import AtomConstants as AC
from . import RotationMatrices as RM
from . import BasisChanger as BC
from . import JonesMatrices as JM
from . import solve_dielectric as SD
from .param_validation import BUFFER_TEMP_FILL_DEFAULT, validate_p_dict
from .param_validation import warn as _warn   # attributed to the user's line

# direct testing

# Default values for parameters
p_dict_defaults = {	'Elem':'Rb', 'Dline':'D2', 
							'lcell':75e-3,'Bfield':0., 'T':20., 
							'GammaBuf':0., 'shift':0.,
							# B-field angle w.r.t. the +z axis (the forward light's direction)
							'Btheta':0, 'Bphi':0,
							# [ES4 FIX] DoppTemp None, or not given, is T (it was 20 C)
							'Constrain':True, 'DoppTemp':None,
							'rb85frac':72.17, 'K40frac':0.01, 'K41frac':6.73,
							'BoltzmannFactor':True,
							# [ES4] non-orthogonal propagation switch and window refractive index (Briscoe et al. 2025).
							# n_window defaults to None (no windows): it used to default to 1.51, which
							# coupled two independent things -- switching use_non_orthogonal on then also
							# switched Fresnel window losses on, dropping absolute S0 by 15.5% far from
							# resonance, while the orthogonal path ignores n_window entirely. The two
							# solvers now agree on absolute transmission and windows are opt-in.
							# [ES4 MOD] use_non_orthogonal defaults to True. The
							# orthogonal path (ES3's, kept for compatibility as False) is exact
							# only in Faraday and Voigt geometry: in oblique fields it is off by
							# up to ~0.5 in S0 (Rb D2, 2 mm, 120 C, 239 G, Btheta = 80 deg).
							'use_non_orthogonal':True, 'n_window':None}

# [ES4 FIX] Field components in the frame of a reversed pass
# (wavevector_dirn = -1): turning by pi about x reverses y and z. Applying it
# twice is the identity, so the same factor takes the field in and back out.
_REVERSED_FRAME = np.array([[1.], [-1.], [-1.]])

# [ES4 FIX] ES3 truncated each transition frequency to whole MHz.
# ES4 keeps the exact value. Set this to True only to reproduce ES3's output
# exactly: the ES3-equivalence tests (test_es3_golden and the published-figure
# anchors) do, so that they still check everything else to ~1e-9.
_ES3_WHOLE_MHZ_LINES = False

# [ES4 FIX] ES3 added 273.16 to the temperature in C for the
# Boltzmann populations and 273.15 everywhere else. ES4 uses 273.15 for both
# (the populations move by ~1e-8). Set this to 273.16 only to reproduce ES3's
# output to ~1e-9: the ES3 golden tests do.
_KELVIN_BOLTZMANN = 273.15

def FreqStren(groundLevels,excitedLevels,groundDim,
			  excitedDim,Dline,hand,BoltzmannFactor=True,T=293.15):
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
				# [ES4 FIX] ES3 truncated this to whole MHz
				# with int(), which put each line up to 1 MHz too close to
				# the line centre and made the spectrum a staircase in B:
				# a fit's finite-difference Jacobian then missed the line
				# shift, and ML fits of Bfield stalled near their start.
				frequency = (-groundLevels[gg][0].real
							 +excitedLevels[ee][0].real)
				transitionFrequency[transNo] = (int(frequency) if _ES3_WHOLE_MHZ_LINES
												else frequency)
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

	Returns (ChiPlus, ChiMinus, ChiZ): the sigma+ (delta m = +1), sigma- (delta m = -1)
	and pi (delta m = 0) susceptibilities. [ES4 FIX] ES3, and the 4.0.0
	alphas, returned the sigma- manifold first, under the name ChiPlus.

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
			T		 			float			Celsius	Temperature used to calculate atomic number density; with
															Constrain=True (the default), also the temperature of the
															vapour in the beam (see DoppTemp)
			GammaBuf   	float			MHz		Extra Lorentzian broadening (FWHM, MHz), added to the
															natural linewidth: the static alternative to the
															buffer-gas model, or any other homogeneous broadening.
															Replaced by the model's broadening when BufferGas is set
															and BufferPressure > 0 (with a warning if it is not 0).
			shift	  			float			MHz		Global frequency shift of all resonance lines (MHz).
															[ES4 MOD] With the buffer-gas model on, its pressure
															shift is added to this one (it used to replace it).

			DoppTemp   	float			Celsius	Temperature of the vapour in the beam (the cell body, when a
															colder stem sets T): the Doppler width, and the buffer-gas
															model's broadening and shift. Used only with
															Constrain=False; not given, it is T
			Constrain  		bool			--			If True, overides the DoppTemp value and sets it to T

			# Elemental abundancies, where applicable
			rb85frac   		float			%			percentage of rubidium-85 atoms
			K40frac			float			%			percentage of potassium-40 atoms
			K41frac			float			%			percentage of potassium-41 atoms
			
			lcell	  			float			m			length of the vapour cell
			
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
	if 'Constrain' in list(p_dict.keys()):
		Constrain = p_dict['Constrain']
	else:
		Constrain = p_dict_defaults['Constrain']

	# [ES4 FIX] The temperature of the vapour in the beam, which sets
	# the Doppler width and the buffer-gas model: T, unless Constrain=False gives
	# the cell body its own DoppTemp (T then being the stem's, which sets the
	# number density). A DoppTemp not given (or None) is T: it used to be 20 C, so
	# Constrain=False without DoppTemp put the Doppler width at 20 C whatever T was.
	if Constrain or (p_dict.get('DoppTemp') is None):
		DoppTemp = T
	else:
		DoppTemp = p_dict['DoppTemp']

	# [ES4 NEW] Dynamic buffer gas model, used when 'BufferGas' is set: the
	# temperature-dependent broadening and pressure shift of the Pitz 2014 model
	# (see elecsus/libs/BufferGasData.py). Its broadening replaces GammaBuf, and
	# its shift is added to shift. To use the static model, omit 'BufferGas'
	# from the parameter dictionary or set it to None.
	# [ES4 FIX] It is evaluated at DoppTemp, the temperature of the
	# gas the light passes through, where the collisions happen. With
	# Constrain=True it is T.
	if p_dict.get('BufferGas') is not None:
		from elecsus.libs.BufferGasData import get_buffer_broadening_and_shift
		_bg   = p_dict['BufferGas']
		_bg_P = p_dict.get('BufferPressure', 0.0)
		_bg_Tf = p_dict.get('BufferTemp_fill', BUFFER_TEMP_FILL_DEFAULT)
		if _bg_P <= 0.0:
			_warn(
				f"BufferGas='{_bg}' is set but BufferPressure={_bg_P} Torr. "
				"Set 'BufferPressure' to the fill pressure in Torr.")
		else:
			_gamma_phys, _shift_phys = get_buffer_broadening_and_shift(
				Elem, Dline, _bg, _bg_P, _bg_Tf, DoppTemp + 273.15)
			if _gamma_phys is not None:
				if GammaBuf != 0.0:
					_warn(
						f"BufferGas model active: GammaBuf={GammaBuf} MHz is overridden "
						f"by the Pitz model value ({_gamma_phys:.2f} MHz). "
						"Remove GammaBuf from your parameter dict to suppress this warning.")
				GammaBuf = _gamma_phys
				# [ES4 MOD] shift is an offset of the frequency axis (the laser's
				# calibration, say), which a fit may need on top of the pressure
				# shift. It was overridden, with a warning.
				shift    = shift + _shift_phys
			else:
				_warn(
					f"No buffer gas data for ({Elem}, {Dline}, '{_bg}'). "
					"Falling back to static GammaBuf/shift.")

	if 'rb85frac' in list(p_dict.keys()):
		rb85frac = p_dict['rb85frac']
	else:
		rb85frac = p_dict_defaults['rb85frac']
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
													Dline,'Left',BoltzmannFactor,T+_KELVIN_BOLTZMANN)		  

			# Rb-85 allowed transitions for light driving sigma plus
			renergy85, rstrength85, rtransno85 = FreqStren(
													Rb85_ES.groundManifold,
													Rb85_ES.excitedManifold,
													Rb85_ES.ds,Rb85_ES.dp,
													Dline,'Right',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
			
			# Rb-85 allowed transitions for light driving pi
			zenergy85, zstrength85, ztransno85 = FreqStren(
													Rb85_ES.groundManifold,
													Rb85_ES.excitedManifold,
													Rb85_ES.ds,Rb85_ES.dp,
													Dline,'Z',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
			

		if rb87frac!=0.0:
			Rb87atom = AC.Rb87
			#Hamiltonian(isotope,transition,gL,Bfield)
			Rb87_ES = ES.Hamiltonian('Rb87',Dline,1.0,Bfield)
			# Rb-87 allowed transitions for light driving sigma minus
			lenergy87, lstrength87, ltransno87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Left',BoltzmannFactor,T+_KELVIN_BOLTZMANN)

			# Rb-87 allowed transitions for light driving sigma plus
			renergy87, rstrength87, rtransno87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Right',BoltzmannFactor,T+_KELVIN_BOLTZMANN)

			# Rb-87 allowed transitions for light driving sigma plus
			zenergy87, zstrength87, ztransno87 = FreqStren(
													Rb87_ES.groundManifold,
													Rb87_ES.excitedManifold,
													Rb87_ES.ds,Rb87_ES.dp,
													Dline,'Z',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
													
			
													
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
												 'Left',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
		renergy, rstrength, rtransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Right',BoltzmannFactor,T+_KELVIN_BOLTZMANN)		
		zenergy, zstrength, ztransno = FreqStren(Cs_ES.groundManifold,
												 Cs_ES.excitedManifold,
												 Cs_ES.ds,Cs_ES.dp,Dline,
												 'Z',BoltzmannFactor,T+_KELVIN_BOLTZMANN)

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
												 'Left',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
		renergy, rstrength, rtransno = FreqStren(Na_ES.groundManifold,
												 Na_ES.excitedManifold,
												 Na_ES.ds,Na_ES.dp,Dline,
												 'Right',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
		zenergy, zstrength, ztransno = FreqStren(Na_ES.groundManifold,
												 Na_ES.excitedManifold,
												 Na_ES.ds,Na_ES.dp,Dline,
												 'Z',BoltzmannFactor,T+_KELVIN_BOLTZMANN)												 
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
													'Left',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
			renergy39, rstrength39, rtransno39 = FreqStren(
													K39_ES.groundManifold,
													K39_ES.excitedManifold,
													K39_ES.ds,K39_ES.dp,Dline,
													'Right',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
			zenergy39, zstrength39, ztransno39 = FreqStren(K39_ES.groundManifold,
													K39_ES.excitedManifold,
													K39_ES.ds,K39_ES.dp,Dline,
													'Z',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
		if K40frac!=0.0:
			K40atom = AC.K40
			K40_ES = ES.Hamiltonian('K40',Dline,1.0,Bfield)
			
			lenergy40, lstrength40, ltransno40 = FreqStren(
													K40_ES.groundManifold,
													K40_ES.excitedManifold,
													K40_ES.ds,K40_ES.dp,Dline,
													'Left',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
			renergy40, rstrength40, rtransno40 = FreqStren(
													K40_ES.groundManifold,
													K40_ES.excitedManifold,
													K40_ES.ds,K40_ES.dp,Dline,
													'Right',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
			zenergy40, zstrength40, ztransno40 = FreqStren(K40_ES.groundManifold,
													K40_ES.excitedManifold,
													K40_ES.ds,K40_ES.dp,Dline,
													'Z',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
		if K41frac!=0.0:
			K41atom = AC.K41
			K41_ES = ES.Hamiltonian('K41',Dline,1.0,Bfield)
			
			lenergy41, lstrength41, ltransno41 = FreqStren(
													K41_ES.groundManifold,
													K41_ES.excitedManifold,
													K41_ES.ds,K41_ES.dp,Dline,
													'Left',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
			renergy41, rstrength41, rtransno41 = FreqStren(
													K41_ES.groundManifold,
													K41_ES.excitedManifold,
													K41_ES.ds,K41_ES.dp,Dline,
													'Right',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
			zenergy41, zstrength41, ztransno41 = FreqStren(K41_ES.groundManifold,
													K41_ES.excitedManifold,
													K41_ES.ds,K41_ES.dp,Dline,
													'Z',BoltzmannFactor,T+_KELVIN_BOLTZMANN)
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
	# (DoppTemp was resolved against Constrain above)
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
	# [ES4 FIX] The 'Right' transitions are the sigma+ (delta m = +1)
	# manifold, blue-shifted at high field, and the 'Left' ones sigma-. ChiPlus
	# was the 'Left' (sigma-) sum, so a [1,+i] field, which drives sigma+ about
	# B, was given the sigma- susceptibility: the code modelled B reversed.
	totalChiPlus = ChiRealRight + 1.j*ChiImRight
	totalChiMinus = ChiRealLeft + 1.j*ChiImLeft
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
		n_window:   refractive index of cell windows (None, the default, models a
		            windowless cell; 1.51 is typical for borosilicate glass).

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
		            contains no window correction; that path is kept as ES3 computed it.

	Returns:
		E_out: (N, 3, 1) output E-field vectors, transmitted into air: E_z = 0
	"""
	N = len(n1)

	# Extract transverse (x,y) components of eigenvectors
	# RM_ary[0,:,i] = eigenvector 1 at point i
	# RM_ary[1,:,i] = eigenvector 2 at point i
	u1_xy = RM_ary[0, :2, :].T  # shape (N, 2)
	u2_xy = RM_ary[1, :2, :].T  # shape (N, 2)

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
	# [ES4 FIX] The field leaving the cell is a plane wave along z
	# in air, so it has no z component: E_x and E_y are continuous across the
	# exit face (and carry the window factors above), but the modes' E_z,
	# c1 u1_z + c2 u2_z, exists only inside the anisotropic vapour, where D,
	# not E, is transverse. It used to be returned, so S0 counted |E_z|^2 and
	# exceeded Ix + Iy (by up to ~2e-6 of S0 in oblique fields) and a pure
	# state came out not quite fully polarised.
	E_out = np.zeros((N, 3, 1), dtype='complex')
	E_out[:, 0, 0] = c[:, 0] * u1_xy[:, 0] + c[:, 1] * u2_xy[:, 0]  # E_x
	E_out[:, 1, 0] = c[:, 0] * u1_xy[:, 1] + c[:, 1] * u2_xy[:, 1]  # E_y

	return E_out
# [ES4 END] propagate_non_orthogonal


# [ES4 BEGIN] Poynting-flux transmission for the non-orthogonal path
def poynting_transmission_non_orthogonal(E_xz, RM_ary, n1, n2, wavenumber, lcell):
	"""
	Physical (Poynting-flux) transmission through the vapour for the
	non-orthogonal propagation path:

	    S0_P = <S_z(L)> / <S_z(0)>,
	    <S_z(z)> ∝ Re[ Σ_jk c_j(z) c_k*(z) n_k* (m_j · m_k*) ]

	where c_j(z) = c_j exp(i n_j k0 z) and m_j are the transverse eigenvector
	components. This is the time-averaged z-component of the Poynting vector
	of the two-mode field, using H_k = (n_k/Z0) ẑ×E_k per mode.

	Rationale: |E|² is not proportional to power flux inside an
	anisotropic medium, and near mode degeneracies the biorthogonal
	amplitudes c_j grow without bound. Before the eigenvector fix (Webber) and
	the transverse E_out, 'S0' could exceed unity by orders of
	magnitude this way although no energy is created; it is now ≤ 1 in every
	scan the test suite makes (test_symmetries.py). The Poynting ratio is the
	physically transported power and satisfies S0_P ≤ 1 for a passive
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

	# [ES4 FIX] NaN where no light enters (cell 2 of a cascade behind an opaque
	# cell 1), as the plain division gave, but without its RuntimeWarning
	return _over_input(_flux(cL), _flux(c0))
# [ES4 END] poynting_transmission_non_orthogonal


def get_Efield(X, E_in, Chi, p_dict, verbose=False, return_poynting=False):
	""" 
	Most general form of calculation - return the electric field vector E_out. 
	Can use Jones matrices to calculate all other experimental spectra from 
	this, as get_spectra() does
	
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
		Chi:		(3,len(X)) array of susceptibility values, for sigma+, sigma-, and pi transitions,
						as calc_chi returns them
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
		p_dict_defaults['use_non_orthogonal'])
	n_window = p_dict.get('n_window',
		p_dict_defaults.get('n_window', None))

	# direction of wavevector (for double-pass geometry simulations)
	if 'wavevector_dirn' in list(p_dict.keys()):
		wavevector_dirn = p_dict['wavevector_dirn']
	else:
		wavevector_dirn = 1
	
	# [ES4 FIX] Light travelling along -z is solved as a forward
	# problem in the frame turned by pi about the x-axis, (x, y, z) ->
	# (x, -y, -z), where it travels along +z. B is given in the lab frame for
	# both passes, so in that frame it has polar angle pi - Btheta and azimuth
	# -Bphi, and the field's y and z components change sign on the way in and
	# again on the way out. The wavenumber stays positive. ES3 instead turned
	# B through pi and made the wavenumber negative: the Faraday rotation came
	# out right, but exp(i n k L) then grew wherever the vapour absorbs, and
	# the transmission reached ~1e62.
	reverse = (wavevector_dirn == -1)
	if reverse:
		Btheta = np.pi - Btheta
		Bphi = -Bphi
		E_in = E_in * _REVERSED_FRAME

	# get wavenumber
	transition = AC.transitions[Elem+Dline]
	wavenumber = transition.wavevectorMagnitude

	# get susceptibility (already calculated, input to this method)
	ChiPlus, ChiMinus, ChiZ = Chi

	# Rotate initial Electric field so that B field lies in x-z plane
	# (Effective polarisation rotation)
	# [ES4 FIX] B has azimuth Bphi, measured from x towards y, so
	# the frame turns by -Bphi to bring it into the x-z plane, and back by
	# +Bphi on the way out. Turning by +Bphi put B at azimuth -Bphi.
	E_xz = RM.rotate_around_z(E_in.T,-Bphi)
		
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
			np.transpose(E_out_xz,(0,2,1))[::,0], Bphi)
		if return_poynting:
			# [ES4] Poynting-flux transmission: physical power
			# ratio, bounded by 1 for the passive vapour, unlike |E|²-based
			# S0 near mode degeneracies.
			S0_P = poynting_transmission_non_orthogonal(
				E_xz, RM_ary, n1, n2, wavenumber, lcell)
	else:
		if return_poynting:
			raise ValueError(
				"S0_P (Poynting-flux transmission) requires "
				"use_non_orthogonal=True: it is defined for the "
				"non-orthogonal mode decomposition, where |E|²-based S0 "
				"can exceed 1 near mode degeneracies.")
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
			np.transpose(fast_E_out_xz,(0,2,1))[::,0], Bphi)
	# [ES4 END] propagation branch

	E_out = fast_E_out.T[0]
	if reverse:
		# [ES4 FIX] back from the reversed-pass frame to the lab
		E_out = E_out * _REVERSED_FRAME

	## return electric field vector - can then use Jones matrices to do everything else
	# [ES4 MOD] np.matrix deprecated; RM_ary is already ndarray. On a reversed
	# pass it is in the solver's frame (y and z reversed), not the lab frame.
	if return_poynting:
		return E_out, RM_ary.T[-1], S0_P
	return E_out, RM_ary.T[-1]

def _over_input(I, I_in):
	"""
	I / I_in: an intensity relative to the input intensity I_in (a scalar or
	an (N,) array).

	[ES4 FIX] Where no light enters the result is NaN, as the plain division
	gave, but without numpy's 'invalid value' and 'overflow' RuntimeWarnings.
	'No light' is I_in below the smallest normal float (~2e-308), 0 included:
	numpy divides a complex number by forming 1/I_in, which overflows there,
	so the plain division returned 0/0 or 0*inf. Behind a crossed polariser
	the second cell of calculate_cascade gets no light wherever the first is
	opaque, and the GUI showed those warnings on its cascade tab. Every other
	point is divided exactly as before, to the last bit.
	"""
	lit = _lit(I_in)
	out = np.full(np.broadcast(I, I_in).shape, np.nan, dtype=np.result_type(I, I_in))
	return np.divide(I, I_in, out=out, where=lit)


def _lit(I_in):
	"""
	True where light enters: the input intensity I_in is at least the smallest
	normal float (~2e-308). The test _over_input and calculate_cascade share.
	"""
	return np.abs(np.asarray(I_in)) >= np.finfo(float).tiny


# [ES4] The intensities and Stokes parameters of a field, shared
# by get_spectra and calculate_cascade's instrument totals. The cascade had its
# own copy, which kept the 4.0.0 alphas' S3 sign after it was corrected here.
def field_outputs(E_out, I_in, S0_floor=0.):
	"""
	Intensities and fractional Stokes parameters of a (3, N) field E_out,
	normalised to the input intensity I_in (a scalar or an (N,) array).

	S1, S2 and S3 are divided by S0 where S0 > S0_floor, and are left
	undivided elsewhere, where |S_i| <= S0 <= S0_floor: 0 to within the floor.
	S0_floor = 0 guards only a fully opaque point; calculate_cascade passes
	1e-15, so that the crossed-polariser wings return at most 1e-15 rather
	than noise.
	Where no light enters (I_in ~ 0) every output is NaN (see _over_input).
	S3 = (Ir - Il)/S0 follows Born & Wolf (see get_spectra).

	Returns a dict with S0, S1, S2, S3, Ix, Iy, I_P45, I_M45, Il and Ir, all real.
	"""
	# [ES4 MOD] Intensities and Stokes parameters are real by construction
	# (E . E*), but the sum over a complex field leaves a complex128 dtype with
	# a ~1e-17 imaginary residue.  Take the real part here, at the point where
	# each observable is formed, so that every intensity-like output is float64:
	# a complex dtype propagates a ComplexWarning into user plotting/fitting code
	# and makes np.max compare lexicographically rather than by magnitude.
	S0 = _over_input((E_out * E_out.conjugate()).sum(axis=0), I_in).real

	#S1#
	Ex = JM.HorizPol_xy @ E_out[:2]
	Ix =  _over_input((Ex * Ex.conjugate()).sum(axis=0), I_in).real
	Ey =  JM.VertPol_xy @ E_out[:2]
	Iy =  _over_input((Ey * Ey.conjugate()).sum(axis=0), I_in).real

	#S2#
	E_P45 =  JM.LPol_P45_xy @ E_out[:2]
	E_M45 =  JM.LPol_M45_xy @ E_out[:2]
	I_P45 = _over_input((E_P45 * E_P45.conjugate()).sum(axis=0), I_in).real
	I_M45 = _over_input((E_M45 * E_M45.conjugate()).sum(axis=0), I_in).real

	#S3#
	# change to circular basis
	E_out_lrz = BC.xyz_to_lrz(E_out)
	El =  JM.CPol_L_lr @ E_out_lrz[:2]
	Er =  JM.CPol_R_lr @ E_out_lrz[:2]
	Il = _over_input((El * El.conjugate()).sum(axis=0), I_in).real
	Ir = _over_input((Er * Er.conjugate()).sum(axis=0), I_in).real

	# [ES4 MOD] S1, S2 and S3 fractional (normalised by S0); ES3 returned Ix - Iy,
	# I_P45 - I_M45 and Ir - Il (absolute intensity differences, units of input
	# intensity).  Stokes parameters are dimensionless by definition and must
	# satisfy S0² ≥ S1² + S2² + S3² (Born & Wolf §1.4).  [ES4] guard only
	# S0 == 0 (fully opaque) by default; at nonzero S0 the light is still fully
	# polarised, so the true fractional Stokes ratio must be kept (do not floor S0).
	safe = np.where(S0 > S0_floor, S0, 1.0)
	S1 = (Ix - Iy) / safe
	S2 = (I_P45 - I_M45) / safe
	# [ES4 FIX] S3 in Born & Wolf's sign again, as in ES3:
	# S3 = (Ir - Il)/S0, positive when right-circular light dominates. With
	# fields varying as exp(i(kz - wt)), Il is the [1,+i]/sqrt2 component, whose
	# field turns anticlockwise for an observer facing the source: left-circular
	# in Born & Wolf (Sec. 1.4.2), where right means clockwise. The 4.0.0 alphas
	# returned -(Ir - Il)/S0, citing Born & Wolf for the opposite sign.
	# Il, Ir and S3 are about the lab +z axis, as E_in and E_out
	# are: on a reversed pass (wavevector_dirn = -1), S3 > 0 is left-circular
	# for the returning beam itself.
	S3 = (Ir - Il) / safe
	return {'S0': S0, 'S1': S1, 'S2': S2, 'S3': S3, 'Ix': Ix, 'Iy': Iy,
			'I_P45': I_P45, 'I_M45': I_M45, 'Il': Il, 'Ir': Ir}
# [ES4 END] field_outputs


def get_spectra(X, E_in, p_dict, outputs=None, strict=True):
	""" 
	Calls get_Efield() to get Electric field, then use Jones matrices 
	to calculate experimentally useful quantities.
	
	
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
				T		 			float			Celsius	Temperature used to calculate atomic number density; with
																Constrain=True (the default), also the temperature of the
																vapour in the beam (see DoppTemp)
				GammaBuf   	float			MHz		Extra Lorentzian broadening (FWHM, MHz), added to the
																natural linewidth. Replaced by the buffer-gas model's
																broadening when BufferGas is set and BufferPressure > 0.
				shift	  			float			MHz		Global frequency shift of all resonance lines (MHz).
																With the buffer-gas model on, its pressure shift is
																added to this one.
				DoppTemp   	float			Celsius	Temperature of the vapour in the beam (the cell body, when a
																colder stem sets T): the Doppler width, and the buffer-gas
																model's broadening and shift. Used only with
																Constrain=False; not given, it is T
				Constrain  		bool			--			If True, overides the DoppTemp value and sets it to T

				# Elemental abundancies, where applicable
				rb85frac   		float			%			percentage of rubidium-85 atoms
				K40frac			float			%			percentage of potassium-40 atoms
				K41frac			float			%			percentage of potassium-41 atoms
				
				lcell	  			float			m			length of the vapour cell
				
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
	
	# [ES4 NEW] Reject mistyped and unphysical parameters before anything is
	# computed.  Unknown keys are otherwise ignored in silence, which quietly
	# disables whichever feature was being asked for (see libs/param_validation.py).
	validate_p_dict(p_dict, strict=strict)

	# [ES4 FIX] A single detuning is a 1-point axis. The scalar
	# branch below used to run only after calc_chi had already taken len(X),
	# and it missed numpy scalars anyway.
	X = np.atleast_1d(np.asarray(X, dtype=float))

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
	# [ES4] 'theta0'/'Pol' (ES3-only, no effect) warn from validate_p_dict's
	# check_keys, which fit_data also runs.

	# get wavenumber
	transition = AC.transitions[Elem+Dline]

	wavenumber = transition.wavevectorMagnitude
	
	# Calculate Susceptibility
	# [ES4 FIX] ChiPlus is the sigma+ (delta m = +1) manifold, which
	# a [1,+i]/sqrt2 field drives when B is along +z, and ChiMinus sigma-. In ES3
	# and the 4.0.0 alphas they were the other way round, and the code modelled B reversed.
	ChiPlus, ChiMinus, ChiZ = calc_chi(X, p_dict)
	Chi = [ChiPlus, ChiMinus, ChiZ]
	
	# Complex refractive index.
	nPlus = sqrt(1.0+ChiPlus)    # sigma+: the [1,+i]/sqrt2 mode for B along +z
	nMinus = sqrt(1.0+ChiMinus)  # sigma-: the [1,-i]/sqrt2 mode for B along +z
	nZ = sqrt(1.0+ChiZ) # Complex index driving pi transitions


	
	# Calculate E_field
	# [ES4] S0_P (Poynting-flux transmission) is computed inside
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
	I_in = (E_in * E_in.conjugate()).sum(axis=0).real
	
	# Intensities and fractional Stokes parameters (shared with
	# calculate_cascade's instrument totals)
	fo = field_outputs(E_out, I_in)
	S0 = fo['S0']
	
	Transmission = S0
	
	
	## Some quantities from Faraday geometry don't make sense when B and k not aligned, but leave them here for historical reasons
	TransLeft = exp(-2.0*nPlus.imag*wavenumber*lcell)
	TransRight = exp(-2.0*nMinus.imag*wavenumber*lcell)
	
	# Faraday rotation angle (including incident linear polarisation angle)
	phiPlus = wavenumber*nPlus.real*lcell
	phiMinus = wavenumber*nMinus.real*lcell
	# [ES4 FIX] The rotation from x towards y is
	# (k L / 2) Re(n_[1,-i] - n_[1,+i]). With B along +z, [1,+i] is sigma+ and
	# that is (phiMinus - phiPlus)/2; with B along -z (Btheta = pi, where the
	# ES3 recipe puts ES3's Faraday geometry) [1,+i] is sigma-, and the
	# rotation is the opposite. It used to keep the +z sign there, opposite to
	# psi. phiPlus and phiMinus stay the sigma+ and sigma- phases. A reversed
	# pass needs nothing: its field, in lab components, is the forward pass's.
	if 'Btheta' in list(p_dict.keys()):
		Btheta = p_dict['Btheta']
	else:
		Btheta = p_dict_defaults['Btheta']
	phi = (1.0 if np.cos(Btheta) >= 0 else -1.0) * (phiMinus-phiPlus)/2.0
	##
	
	#Stokes parameters, computed in field_outputs above
	S1, S2, S3 = fo['S1'], fo['S2'], fo['S3']
	Ix, Iy, I_P45, I_M45, Il, Ir = (fo[k] for k in ('Ix', 'Iy', 'I_P45', 'I_M45', 'Il', 'Ir'))

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
	# there, Btheta = 0 and pi alike (verified to 1e-12 in
	# test_faraday_outputs): dichroism alters the ellipticity, not the
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
	# [ES4 FIX] nPlus is now the sigma+ index, so alphaPlus no
	# longer crosses over to nMinus: the values are unchanged
	alphaPlus = 2.0*nPlus.imag*wavenumber
	alphaMinus = 2.0*nMinus.imag*wavenumber
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
				# [ES4] Poynting-flux transmission; requires
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
	# [ES4 FIX] 'Plus'/'Minus' on Chi/n/phi/alpha is the sigma label:
	# ChiPlus/nPlus/phiPlus/alphaPlus are the sigma+ (delta m = +1) manifold, which
	# a [1,+i]/sqrt2 field drives when B is along +z. ES3, and the 4.0.0 alphas,
	# returned the sigma- quantities under the 'Plus' names (except alphaPlus).
	# Il/Ir and S3 follow Born & Wolf: Il is the [1,+i]/sqrt2 component, which is
	# left-circular for fields varying as exp(i(kz - wt)). See Sec. 5 of the manual.
	# [ES4] All are in lab components, about +z, so on a reversed
	# pass the left and right are reversed for the returning beam.
	tstr = " \
	All possible outputs from the get_spectra method: \n\n\
	Variable Name		Description \n \
	S0						Total transmission through the cell (Ix + Iy) \n\
	S1						Stokes parameter (Ix - Iy)/S0 \n\
	S2						Stokes parameter (I_45 - I_-45)/S0 \n\
	S3						Stokes parameter (Ir - Il)/S0: positive when right-circular light dominates (Born & Wolf). \n\
								Il, Ir and S3 are about the lab +z axis: on a reversed pass (wavevector_dirn=-1), \n\
								S3 > 0 is left-circular for the returning beam \n\
	TransLeft			Single-mode transmission exp(-2 Im(nPlus) k L) of the sigma+ mode [the [1,+i]/sqrt2 (left-circular) \n\
								mode for B along +z; Faraday geometry only] \n\
	TransRight			Single-mode transmission exp(-2 Im(nMinus) k L) of the sigma- mode [the [1,-i]/sqrt2 \n\
								(right-circular) mode for B along +z; Faraday geometry only] \n\
	S0_P					Poynting-flux (physical power) transmission of the bare vapour column. \n\
							Requires use_non_orthogonal=True; excludes window Fresnel factors; \n\
							bounded by 1 for the passive medium even at mode degeneracies where \n\
							|E|^2-based S0 is not proportional to transported power \n\
	ChiPlus				Complex susceptibility of the sigma+ (delta m=+1) manifold [governs [1,+i]/sqrt2 light for B along +z] \n\
	ChiMinus				Complex susceptibility of the sigma- (delta m=-1) manifold [governs [1,-i]/sqrt2 light for B along +z] \n\
	nPlus					Complex refractive index sqrt(1+ChiPlus) of the sigma+ manifold [Faraday geometry only] \n\
	nMinus				Complex refractive index sqrt(1+ChiMinus) of the sigma- manifold [Faraday geometry only] \n\
	phiPlus				Single-pass phase k Re(nPlus) L of the sigma+ mode [Faraday geometry only] \n\
	phiMinus				Single-pass phase k Re(nMinus) L of the sigma- mode [Faraday geometry only] \n\
	phi					Net Faraday rotation of linear polarisation, positive from x towards y: (phiMinus - phiPlus)/2 \n\
									for B along +z, and its negative for B along -z [Faraday geometry only]. \n\
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
	Ir						Intensity of the [1,-i]/sqrt2 component: right-circular (Born & Wolf), clockwise for an \n\
								observer facing the source [lab components: for a reversed pass, left-circular] \n\
	Il						Intensity of the [1,+i]/sqrt2 component: left-circular (Born & Wolf) [lab components, as Ir] \n\
	alphaPlus			Absorption coefficient 2 Im(nPlus) k of the sigma+ manifold \n\
	alphaMinus			Absorption coefficient 2 Im(nMinus) k of the sigma- manifold \n\
	alphaZ				Absorption coefficient 2 Im(nZ) k of the pi (delta m=0) manifold \n\
	"
	print(tstr)