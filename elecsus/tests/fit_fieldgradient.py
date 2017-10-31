""" 
Example code to fit a magnetic field gradient (figure 8 of 2017 ElecSus paper)
This script will reproduce the entirety of figure 8, including the fitting to data and figure layout

JK, 2017
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import copy
import lmfit as lm
import cPickle as pickle

# update matplotlib fonts etc
plt.rc('font',**{'family':'Serif','serif':['Times New Roman']})
params={'axes.labelsize':13,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize': 11,'mathtext.fontset':'cm','mathtext.rm':'serif'}
plt.rcParams.update(params)

## Import elecsus modules
# if elecsus is installed...
from elecsus.libs.durhamcolours import *
from elecsus.libs.spectra import get_spectra
import elecsus.libs.JonesMatrices as JM

# alternately, if elecsus is not installed, add the source directory to sys.path
#sys.path.append('C:\Path\To\ElecSus\...\elecsus\libs')
#from durhamcolours import *
#from spectra import get_spectra
#import JonesMatrices as JM


def ringMagnetBfield(z,z0,IR,OR,length,strength=1.42/2):
	""" 
	Calculate the on axis axial (z-axis) magnetic field from a ring magnet
	See Lee Weller's thesis, appendix F, for derivation.
	Calculation is based on axial field from a solid cylinder; the ring magnet is the outer body inner (hole) solid is subtracted from the outer solid
	
	All dimensions in metres
	
	Inputs:
		
		z			:	1D numpy array of points along the z-axis, at which the magnetic field is evaluated
		z0		:	positional offset of the magnet
		IR			:	inner radius of the ring magnet
		OR		:	outer radius of the ring magnet
		length	:	cylindrical length of the ring magnet
		
	Options:
	
		strength	:	Remnant magnetic field of the magnetic material, in Tesla. Defaults to the value for N52 grade.
	
	Outputs:
	
		B			:	1D numpy array of magnetic field strength, in Tesla
		
	"""
	outerB1 = strength*(z-z0+length)/np.sqrt(((z-z0+length)**2)+(OR**2))- strength*(z-z0-length)/np.sqrt(((z-z0-length)**2)+(OR**2))
	innerB1 = strength*(z-z0+length)/np.sqrt(((z-z0+length)**2)+(IR**2))- strength*(z-z0-length)/np.sqrt(((z-z0-length)**2)+(IR**2))
	B = - (outerB1-innerB1)
	return B

def tophat_profile(z,halfsep=46e-3):
	""" 
	Calculate B-field for a pair of top-hat magnets, separated by a distance 2*halfsep (m), centred around z=0 (each magnet is at +/-halfsep)
	
	Each top-hat magnet is made up from two cylindrical ring magnets, with the same inner diameter but different outer diameters
	
	The parameters given here are the relevant ones to reproduce the field profile in figure 8 of the Elecsus (2017) paper
	"""
	
	# specific top-hat cross-section
	IR = 2.1e-3
	OR1 = 14.5e-3
	OR2 = 19.5e-3
	L1 = 13.3e-3
	L2 = 9.5e-3
		
	B_R = np.zeros_like(z)
	B_L = np.zeros_like(z)
	#first top-hat magnet
	B_R += ringMagnetBfield(z,halfsep+L1/2,IR,OR1,L1/2)
	B_R += ringMagnetBfield(z,halfsep+L1+L2/2,IR,OR2,L2/2)
	#second top-hat (axis is flipped)
	B_L += ringMagnetBfield(-z,halfsep+L1+L2/2,IR,OR2,L2/2)
	B_L += ringMagnetBfield(-z,halfsep+L1/2,IR,OR1,L1/2)
	
	B_tot = B_L + B_R
	
	return B_tot
	
def test_fieldprofile():
	""" Test method to generate magnetic field profile """
	z = np.linspace(-1e-1,1e-1,500)
	
	Bz = tophat_profile(z)
	
	plt.plot(z,-Bz) # Bz is inverted to agree with fig 8 - it's just a choice of direction for the z-axis.
	plt.xlabel('axial position (mm)')
	plt.ylabel('magnetic field (T)')
	plt.show()

def fieldgrad_fitfn(x, sep, offset_adj, temperature, n_segments=25,return_S1=False, verbose=False):
	""" 
	Fit function to calculate transmission spectra in Faraday geometry, with a magnetic field gradient
	
	The magnetic field is from a pair of Nd top-hat magnets placed next to a 75 mm reference cell,
	calculated using the tophat_profile() method above
	
	The calculation is based on splitting the cell into multiple (n_segments) elements, and propagating the electric field through each element separately with a local value of the magnetic field. Convergence tests should be run to ensure that a suitable number of elements are used.
	
	Floating parameters are the separation between the magnets, their position relative to the vapour cell, and the cell temperature.
	
	sep and offset_adj are scaled to mm so that all fit parameters are around the same order of magnitude to prevent potential issues with levenburg-marquardt fitting routines
	
	Inputs:
	
		x						:	1D-array, detuning axis for the calculation, in GHz
		sep					:	Separation between the two top-hat magnets (in mm)
		offset_adj 		:	Adjust the position between the magnets and the vapour cell (in mm)
		temperature	:	cell temperature in C
	
	Options:
	
		n_segments	: 	number of segments cell is split into for calculation
		return_S1		:	If True, returns both S0 and S1 spectra (need False to run fit)
		verbose			:	If True, give more information about present cell segments
		
	Outputs:
		
		S0					:	1D numpy array, Transmission after the cell
		S1					:	Only if return_S1 is True. 1D numpy array, Faraday rotation signal after the cell
			
	"""
	
	# so we can see what the fit is doing...
	print sep, offset_adj, temperature
	
	sep *= 1e-3 				# convert to m from mm
	offset_adj *= 1e-3	# convert to m from mm
	
	# cell length
	LCELL = 75e-3
	#magnetic field angles
	BTHETA = 0
	BPHI = 0
	#element, Dline and isotopic abundance
	RB85FRAC = 72.17
	DLINE = 'D2'
	ELEM = 'Rb'

	# cell starts at z=0, ends at z=LCELL. Offset needed to adjust the position of the magnets such that
	# they are initially centred at the middle of the cell. This offset is adjustable through offset_adj to adjust the relative position of cell and magnet pair.
	z0 = sep/2
	offset = z0 - LCELL/2 - offset_adj
	
	#parameter dictionary
	p_dict = {'Bfield':0,'rb85frac':RB85FRAC,'Btheta':BTHETA*np.pi/180,'Bphi':00*np.pi/180,'lcell':LCELL,'T':temperature,'Dline':DLINE,'Elem':ELEM}
	
	# input polarisation
	pol_in = 1./np.sqrt(2) * np.array([1,1,0]) # from expt
	I_in = 1.
	
	# number of parts to split the cell into (see convergence testing method below)
	#n_segments = 25
	seg_length = LCELL/n_segments
	
	# centre position of each segment
	edges = np.linspace(0,LCELL,n_segments+1)
	seg_centres = (edges[:-1] + edges[1:])/2
	
	# calculate field at centre of each segment
	#seg_fields = ringMagnetBfield(seg_centres,z0,M_ID,M_OD,M_T)
	seg_fields = tophat_profile(seg_centres-z0+offset,z0)
	
	# Evolve electric field over each segment
	for i in range(n_segments):
		if verbose: print 'Cell segment ', i+1, '/', n_segments
		
		# average magnetic field in the segment
		p_dict['Bfield'] = seg_fields[i] * 1e4 # convert to Gauss from Tesla
		p_dict['lcell'] = seg_length
		
		[E_out] = get_spectra(x*1e3,pol_in,p_dict,outputs=['E_out'])
		
		# input field for the next iteration is the output field of this iteration
		pol_in = E_out
	
	# Calculate Stokes parameters...
	Ey =  np.array(JM.VertPol_xy * E_out[:2])
	Iy =  (Ey * Ey.conjugate()).sum(axis=0) / I_in
	Ex =  np.array(JM.HorizPol_xy * E_out[:2])
	Ix =  (Ex * Ex.conjugate()).sum(axis=0) / I_in

	S0 = Ix + Iy
	
	#Other Stokes parameters...
	S1 = Ix - Iy
	
	if return_S1:
		return S0, S1
	else:
		return S0
		
def test_fitfn():
	""" Test method for fitting function """
	
	# parameters which are close to the optimum
	sep = 92.
	temperature = 17.7 #50
	offset_adj = -0.7

	n_segments = 25
	
	detuning = np.linspace(-10,10,1000)
	S0, S1 = fieldgrad_fitfn(detuning, sep, offset_adj, temperature, n_segments,return_S1=True)

	plt.plot(detuning, S0)

	# Try to add data as well from example file...
	try:
		# assumes test data file is in same directory
		d_expt, S0_expt = np.loadtxt('./S0_Bgradient.csv',delimiter=',').T
		d_expt = d_expt[::200]
		S0_expt = S0_expt[::200]
		
		plt.plot(d_expt, S0_expt, '.')
	except:
		print 'Cannot find test data file... '
	
	plt.xlabel('Detuning (GHz)')
	plt.ylabel('Transmission, S0')
	plt.show()
	
def field_gradient_fit(fit=True):
	""" 
	Fit theory curve to experimental data using the non-uniform magnetic field model above (fieldgrad_fitfn).
	
	Uses lmfit module for fitting. Differential evolution is recommended to find global optimum parameters
	(leastsq method may work depending on choice of intiial parameters, but not guaranteed). Fitting takes a while,
	around 10-15 minutes depending on computer speed.
	
	The 'fit' keyword argument is used to either perform the fitting (True) or load in best-fit parameters from previous
	calculations (False).
	
	This method will reproduce figure 8 of the ElecSus paper completely.
	For this the file S0_Bgradient.csv from the 'ElecSusTestData\Field Gradient' GitHub repository is required to be 
	in the same directory as this file.
	
	"""
	## initial guess params for fit
	sep = 92.
	offset_adj = 0.
	temperature = 17.7

	# experimental data; load and crop
	d_expt, S0_expt = np.loadtxt('./S0_Bgradient.csv',delimiter=',').T
	d_expt = d_expt[::200]
	S0_expt = S0_expt[::200]

	#######################################
	
	if fit:
		# Fitting takes a while..!
		
		x = np.array(d_expt)
		y = np.array(S0_expt)
		
		p_dict = {'sep':sep, 'offset_adj':offset_adj, 'temperature':temperature}
		
		# Have a quick look at the guess parameter curve to see if it's close
		S0_trial = fieldgrad_fitfn(d_expt, sep, offset_adj, temperature)
		plt.plot(d_expt,S0_trial)
		plt.plot(d_expt,S0_expt,'.')
		plt.title('Curve based on initial parameters:')
		plt.show() # halts the program here until plot window is closed...
		
		model = lm.Model(fieldgrad_fitfn)
		params = model.make_params(**p_dict)
		
		# Turn off all parameters varying by default, unless specified in p_dict_bools
		allkeys = params.valuesdict()
		for key in allkeys:
			params[key].vary = False
		
		p_dict_bools = {'sep':True, 'offset_adj':True, 'temperature':True}
		p_dict_bounds = {'sep': [75., 120.], 'offset_adj':[-10,10], 'temperature':[15,25]} # sensible boundaries based on lab conditions
			
		# Turn on fitting parameters as specified in p_dict_bools, and set boundaries
		for key in p_dict_bools:
			params[key].vary = p_dict_bools[key]
			
			if p_dict_bounds is not None:
				if key in p_dict_bounds:
					params[key].min = p_dict_bounds[key][0]
					params[key].max = p_dict_bounds[key][1]
		
		
		# need to use global solver - there are local minima in parameter space
		#method = 'leastsq'
		method = 'differential_evolution'
		
		# run the fit
		result = model.fit(y, x=x, params=params, method=method)
		
		print result.fit_report()
		S0_thy = result.best_fit
		# make quick residual plot to look for any remaining structure
		result.plot_residuals()
			
		# save parameters to pickled file so we don't have to re-run the fit every time we want to look at the plot
		fit_params = result.best_values
		pickle.dump(fit_params, open('./bgradient_fitparams.pkl','wb'))
	else:
		# Load in previously calculated fit parameters
		fit_params = pickle.load(open('./bgradient_fitparams.pkl','rb'))
	
	# extract parameters from dictionary for convenience
	sep = fit_params['sep']
	offset_adj = fit_params['offset_adj']
	temperature= fit_params['temperature']

	# detuning axis for theory arrays
	d = np.linspace(-10,10,1000)
	S0, S1 = fieldgrad_fitfn(d, sep, offset_adj, temperature,return_S1=True) # evaluated for best-fit parameters

	# scale parameters back into metres, from mm
	sep *= 1e-3
	offset_adj *= 1e-3
	
	# plot magnetic field with optimised parameters
	z0 = sep/2
	LCELL = 75e-3
	offset = z0 - LCELL/2 - offset_adj
	zs = np.linspace(-60e-3,120e-3,1400)
	Bf = tophat_profile((zs-z0+offset),z0)
	
	# Calculate average magnetic field across cell
	zs_cell = np.linspace(0,LCELL,1000)
	Bf2 = tophat_profile(zs_cell-z0+offset,z0)
	B_avg = Bf2.mean() * 1e4 # in Gauss
	
	# Compare against single pass with average magnetic field - i.e. no gradient
	p_dict = {'rb85frac':72.17,'Btheta':0,'Bphi':0,'lcell':LCELL,'T':temperature,'Dline':'D2','Elem':'Rb'}
	p_dict['Bfield'] = B_avg
	print 'Average field,', B_avg
	print 'Min / Max field across cell: ', Bf2.min()*1e4, Bf2.max()*1e4
	
	# calculate specra with average field
	pol_in = 1./np.sqrt(2) * np.array([1,1,0])
	[Iy_avg, S1_avg, S0_avg] = get_spectra(d*1e3,pol_in,p_dict,outputs=['Iy','S1','S0'])
	
	# also calculate zero-field spectra for comparison (olive shading in fig 8)
	p_dict['Bfield'] = 0
	[Iy_0, S1_0, S0_0] = get_spectra(d*1e3,pol_in,p_dict,outputs=['Iy','S1','S0'])
	# arbitrary scaling of zero-field spectra to fit nicely on the same axes limits as the other data.
	S0_0 = S0_0 * 0.5 + 0.5
	
	
	
	# Set up figure panels with subplot2grid
	fig2 = plt.figure("Figure 8 of ElecSus paper",figsize=(5,6))
	yy = 3
	xx = 9
	axM = plt.subplot2grid((yy,xx),(0,1),colspan=7)
	ax = plt.subplot2grid((yy,xx),(1,0), colspan=xx)
	ax2 = plt.subplot2grid((yy,xx),(2,0),colspan=xx,sharex=ax)
	
	# Plot magnetic field profile
	axM.plot(zs*1e3,-Bf,color=d_blue)
	# format axes for this sub-plot
	axM.xaxis.set_label_position('top')
	axM.tick_params(axis='x',bottom=True,top=True,labelbottom=False,labeltop=True)
	axM.set_xlabel('Axial position, $z$ (mm)')
	axM.set_ylabel('$B_z(z)$ (T)')
	axM.set_yticks([-.4000,-.2000,0,.2000,.4000,.6000])
	axM.set_ylim(-0.1,0.45)
	axM.set_xlim(-50,120)
	# magnet length parameters (for shading on axes)
	L1 = 13.3e-3
	L2 = 9.5e-3
	# Add shading for cell and axial extent of magnets
	axM.axvspan(1e3*(-offset),1e3*(-offset-L1-L2),alpha=0.5,color=d_midblue)
	axM.axvspan(1e3*(sep-offset),1e3*(sep-offset+L1+L2),alpha=0.5,color=d_midblue)
	axM.axvspan(0*1e3,LCELL*1e3,alpha=0.35, color=d_purple)
	
	
	# Add S0 spectral data
	ax.fill_between(d,1,S0_0,label=r'Zero field',color=d_olive,alpha=0.25)
	ax.plot(d,S0_avg,label=r'Uniform, $\langle B_z(z) \rangle$',color=2.5*np.array(d_black),linestyle='dashed')
	ax.plot(d,S0,label='Gradient, $B_z(z)$',color=d_blue)
	ax.plot(d_expt,S0_expt,'.', ms=4, label='Gradient, $B_z(z)$',color=d_purple)
	
	# Add S1 curves
	ax2.plot(d,S1_avg,label=r'Uniform, ($\langle B_z(z) \rangle$)',color=2.5*np.array(d_black),linestyle='dashed')
	ax2.plot(d,S1,label='Gradient, ($B_z(z)$)',color=d_blue)
	
	# Format axes
	ax2.set_xlabel('Detuning (GHz)')
	ax.set_ylabel('$S_0$')
	ax2.set_ylabel('$S_1$')
	plt.setp(ax.get_xticklabels(),visible=False)	
	ax.set_xlim(-10,10)
	
	# Scale to full figure canvas
	plt.tight_layout()
	
	# Save figure images
	fig2.savefig('./field_gradient_fit.png')
	fig2.savefig('./field_gradient_fit.pdf')
	
	# Show figure interactively
	plt.show()

	
def gradient_convergenceTest():
	""" 
	Convergence test for gradient examples
	
	Calculate difference from 'best' case (most segments) - compare RMS between them and plot RMS as function of number of segments
	"""
	
	# maximum number of segments to break the cell into
	N_segs_max = 100
	
	# calculation parameters
	sep = 100.
	temperature = 20. #50
	offset_adj = 0.
	
	# detuning
	detuning = np.linspace(-10,10,1000)
	
	S0_arrays = []
	N_segs = range(N_segs_max,2,-1)
	rms_vals = np.zeros(len(N_segs))
	
	# Calculate spectra
	for ii,n_segments in enumerate(N_segs):
		print 'NUMBER OF SEGMENTS:', ii, n_segments
		S0_arrays.append(fieldgrad_fitfn(detuning, sep, offset_adj, temperature, n_segments))
				
	S0_maxSegs = S0_arrays[-1]
	S0_arrays = np.array(S0_arrays)
		
	for ii in range(len(N_segs)):
		diff = S0_arrays[:,ii] - S0_maxSegs
		
		rms_vals[ii] = np.sqrt(np.mean(diff**2).real)
		
	plt.plot(N_segs, 100*rms_vals)
	plt.ylabel('RMS difference (%)')
	plt.xlabel('Number of segments')
	
	plt.tight_layout()
	
	plt.savefig('./field_gradient_convergencetest.png')
	plt.savefig('./field_gradient_convergencetest.pdf')
	
	plt.show()
	
	
	
if __name__ == '__main__':
	field_gradient_fit()