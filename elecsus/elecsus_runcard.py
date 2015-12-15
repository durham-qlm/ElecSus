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

"""ElecSus: Program to calculate the electric susceptibility

This module takes user inputs from the run card (default runcard.py),
and outputs data, performs fitting and shows plots.

Usage:

	python elecsus.py
or
	python elecsus.py [run card filename]

Modules called from libs folder:

spectra		  -- calculates a spectrum
runcardCheck	 -- checks user inputs
tools			-- functions to read/write files, smooth data and plot.
numberDensityEqs -- find the number density from temperature.

Updated 2015-12-14: JK

"""
import os
import sys
import inspect
import warnings
import time

from shutil import copyfile
from numpy import arange, zeros, array, sqrt

cmd_subfolder = os.path.realpath(
					os.path.abspath(os.path.join(
					os.path.split(inspect.getfile(
					inspect.currentframe()))[0],"libs")))
if cmd_subfolder not in sys.path:
	sys.path.insert(0, cmd_subfolder) #Put libs directory in path

import spectra
import runcardCheck
from tools import fileOutput, plotOutput, smoother, read_in_twoColumn
from numberDensityEqs import *

if os.name == 'posix':
	from time import time as timing #Timing for linux or apple
else:
	from time import clock as timing #Timing for windows

esc_seq = ('\n', '\'','\"','\a','\b','\f','\n','\N','\r','\t','\u','\U','\v','\o')

def main(argv):
	""" 
	Main method. Imports specified runcard.py file from the second system argument.
	Program behaviour is then determined by contents of the runcard.py file. See
	runcard documentation for more information on program options.
	
	Usage examples:
	
		From command line / terminal:
			>> python elecsus_runcard.py path_to/runcard_file.py

		Within python interpreter:
		
			> python
			>>> import elecsus_runcard as er
			>>> er.main([None,'path_to/runcard_file.py'])
			
	!!! NOTE:: Use forward slashes to separate directories, even on Windows !!! 
	
	"""
	
	initialFiles = os.listdir('.') #List files existing before program started


	# Try to import runcard 
	if (len(argv) > 1): #If user specifies a run card name.
		runcardPath = argv[1]
		# if user has entered filename in another directory, check if any backslashes are used in string
		# try to convert string to raw string - replace backslashes (directories) with forward slashes
		# this doesn't work when there are escape sequences in the string, 
		# so warn user of errors and exit if found.
		runcardPath = runcardPath.replace("\\", "/")
		if any(es in runcardPath for es in esc_seq):
			print 'There are escape characters in the runcard file name string.\
					\nUse Forward slashes (/) to separate directories instead of back slashes.'
			sys.exit(1)
			
		# if file is in other directory, add to sys.path
		runcardDir = os.path.dirname(runcardPath)
		sys.path.append(runcardDir)
		
		# get only the file name
		runcardFilename = os.path.split(runcardPath)[1]
		RCFNstripped = runcardFilename.rstrip('.py')
		code1 = 'import ' + RCFNstripped + ' as G'
		try:
			exec code1 in globals()#Import user specified run card
		except:
			print 'Run card file not recognised.'
			sys.exit(1)
	else:
		print 'Run card file name must be specified. Exiting Python...' 
		sys.exit(1)

	print '\n'
	banner = open(os.path.join(cmd_subfolder,'NOTICE'),'r')
	for line in banner: #This loop prints the banner to the terminal
		print line,
	print '\n\n'
	banner.close()

	
	#run user inputs through the run card checker
	ELE,NDT,NDTB,CB,DT,DTB,BF,BFB,R85,R85B,K40,K41,T0,T0B,CL,CLB,POL,POLB,SHI,\
	SHIB,GAM,GAMB,PRE,TRA,SPE,DAT,SMB,SMF,FT,DStart,DStop,PB,SB,PFT =\
	runcardCheck.sanity_check(
		G.element,G.NdenTemp,G.NdenTempBool,G.ConstrainBool,G.DoppTemp,
		G.DoppTempBool,G.Bfield,G.BfieldBool,G.Rb85,G.Rb85Bool,G.K40,G.K41,
		G.Theta0,G.Theta0Bool,G.CellLength,G.CellLenBool,G.Polar,G.PolarBool,
		G.Shift,G.ShiftBool,G.Gamma,G.GammaBool,G.Precision,G.Transition,
		G.Spectrum,G.DataFilename,G.SmoothBool,G.SmoothFactor,G.FitType,
		G.DetStart,G.DetStop,G.PlotBool,G.SavePlotBool,G.PlotFileType
		)

	if TRA == 'D1': #If specified D-line is D1.
		print 'Calculating D1 spectrum...\n' #Tell user program is running.
	elif TRA == 'D2':					   
		print 'Calculating D2 spectrum...\n'

	#Stop warnings about casting complex numbers
	warnings.simplefilter("ignore")

	#Create results output directory
	tryDir = os.path.join(os.getcwd(),DAT+"_output")
	while True:
		outputDirectory = tryDir
		try:
			os.makedirs(outputDirectory)
			break
		except:
			tryDir = tryDir + "1"

	#Copy runcard to output directory for user reference
	try:
		copyfile(runcardPath,os.path.join(
									outputDirectory,runcardFilename))
	except:
		copyfile('runcard.py',os.path.join(outputDirectory,'runcard.py'))

	if FT == 'No fit, just theory':
		startTime = timing()
		GHzPrecision = PRE/1000.0
		#Detuning in GHz
		X	= arange(DStart,DStop+GHzPrecision,GHzPrecision)
		XMHz = X*1.0e3 #Detuning in MHz
		Spec = spectra.spectrum(XMHz,ELE,SPE,BF,NDT,CL,R85,DT,T0,POL,SHI,GAM,
								CB,TRA,PRE,K40,K41)
		print 'Calculation complete.'
	else: #If a fit to experimental data was requested
		try:
			xdata, ydata = read_in_twoColumn(DAT)
			file_read_success = True
		except:
			try:
				xdata, ydata = read_in_twoColumn(os.path.join(runcardDir,DAT))
				file_read_success = True
			except:
				print 'Data file not found. Check file exists in specified location with name: ', runcardFilename
				sys.exit(1)
			
		startTime = timing()
		if SMB:
			print '\nData binned for faster fitting.\n'
			#Take average of local points (low pass filter).
			xdata, ydata = smoother(xdata,ydata,SMF)
		
		X = xdata*1e-3 #Convert back to GHz for plotting
		
		#Generate parameter lists
		paramBoolList = [BFB,NDTB,CLB,R85B,DTB,T0B,POLB,SHIB,GAMB]
		parameters = [ELE,SPE,BF,NDT,CL,R85,DT,T0,POL,SHI,GAM,CB,\
					  TRA,PRE,K40,K41]
		
		# Call different fitting routines		
		if FT == 'Marquardt-Levenberg':
			import MLFittingRoutine as FR
			print '\nPerfoming Marquardt-Levenberg fitting routine.'
			optParams, Spec = FR.MLfit(xdata,ydata,parameters,
													 paramBoolList)
		elif FT == 'Simulated annealing':
			import SAFittingRoutine as FR
			print '\nPerforming fitting by simulated annealing.'
			optParams, Spec = FR.SAFit(xdata,ydata,parameters,
													 paramBoolList)
		else:
			import RRFittingRoutine as FR
			print '\nPerforming fitting by random-restart hill climbing method.'
			#The more parameters to fit, the more evaluations we need to do.
			factor = sum(paramBoolList) 
			evaluationNumber = factor**2 + 5 #integer
			optParams, Spec = FR.RRFit(xdata,ydata,parameters,paramBoolList,
									   evaluationNumber)
		
		optParams = optParams[2:] # compatibility with newer fitting routine (elecsus >V2.0.0)
		
		# Write the fit parameters to a file
		parameterLabels = ['Magnetic field in Gauss =',
						   'Reservoir temperature in Celsius =',
						   'Cell Length in mm =','Rb85 percentage =',
						   'Doppler temperature in Celsius =',
						   'Theta0 in degrees =',
						   'Initial sigma minus polarisation percentage = ',
						   'Shift in MHz =', 
						   'Extra Lorentzian width broadening (MHz) =',
						   'Potassium-40 percentage =',
						   'Potassium-41 percentage =']
		f_Parameters = open(os.path.join(outputDirectory,os.path.split(DAT)[1]+'_Parameters.txt')
							,'w')
		i = 0
		for par in parameterLabels:
			if i == 2:
				optParams[i] = optParams[i]*1000.0
				print >> f_Parameters, parameterLabels[i], optParams[i]
			elif (i==3 and ELE=='Rb') or i==6:
				optParams[i] = optParams[i]*100.0
				print >> f_Parameters, parameterLabels[i], optParams[i]
			elif i == 4 and CB:
				print >> f_Parameters, parameterLabels[i],\
						 "Reservoir temperature (constrained)"
			elif i == 5 and (SPE in ['Ix','Iy','S1','S2']):
				optParams[i] = optParams[i]*180.0
				print >> f_Parameters, parameterLabels[i], optParams[i]
			elif (i in [0,1,4,7,8]):
				print >> f_Parameters, parameterLabels[i], optParams[i]
			elif (i in [9,10]) and ELE == 'K':
				print >> f_Parameters, parameterLabels[i], parameters[i+5]
			i += 1
		numDenVal = numDenRb(optParams[1]+273.15)/1.0e18
		#Find root mean square deviation
		RMS = sqrt(((ydata - Spec)**2).sum()/float(len(ydata)))
		print >> f_Parameters, 'Number density of vapour =', numDenVal,\
				 'x 10^12 cm^-3'
		print >> f_Parameters, "\n\n"
		print >> f_Parameters, "RMS deviation of experimental data and theory =", RMS
		f_Parameters.close()
		#Print parameters to screen
		os.system('cls' if os.name == 'nt' else 'clear') #Clear the terminal
		print ''
		i=0
		for Boolian in paramBoolList:
			if Boolian:
				print parameterLabels[i], optParams[i]
			i+=1

	# Output theory to csv file
	fileDirectory = os.path.join(outputDirectory,os.path.split(DAT)[1]+'_theory.csv')
	fileOutput(fileDirectory,X,Spec)

	# Calculate and output residuals
	if FT != 'No fit, just theory':
		residuals = ydata - Spec
		fileOutput(os.path.join(outputDirectory,os.path.split(DAT)[1]+'_residuals.csv'),
				   X,residuals)
	
	print '\n\nTime taken:', timing() - startTime

	if SB:
		PlotName = os.path.join(outputDirectory,os.path.split(DAT)[1]+'_output'+PFT)
	else:
		PlotName = False

	# Plot results
	if FT == 'No fit, just theory':
		plotOutput(X,Spec,SPE,0,PBool=PB,path=PlotName)
	else:
		plotOutput(X,Spec,SPE,ydata,PB,True,residuals,PlotName)

	#Clean up pyc files
	if os.path.isfile("elecsus.pyc"):
		os.remove("elecsus.pyc")
	if os.path.isfile("runcard.pyc"):
		os.remove("runcard.pyc")
	if (len(argv) > 1):
		NewruncardFilename = runcardFilename + "c"
		BOOL1 = os.path.isfile(NewruncardFilename)
		BOOL2 = (NewruncardFilename in initialFiles)
		if BOOL1 and not BOOL2:
			os.remove(argv[1]+"c")

if __name__ =="__main__": # This line is needed for multiprocessing to work
	main(sys.argv)