==============
ElecSus v2.0.0
==============

A program to calculate the electric susceptibility of an atomic ensemble.
The program is designed to model weak probe laser spectra on the D-lines
of thermal alkali metal vapour cells. The program also includes fitting 
routines which allow experimental parameters to be extracted from 
experimental spectroscopic data.

--------------------
New in version 2.0.0
--------------------

	-	Significantly improved user-friendliness in the form of a 
		GUI to aid in calculating theory spectra and fitting 
		experimental data. 
		
		Works on Windows and Linux, tested on Windows 7, 8.1, 10
		and Ubuntu 14.04. Currently not tested on Mac.

	-	Rewritten fitting routines MLFittingRouine.py,
		RRFittingRoutine.py and SAFittingRoutine.py to support
		keyword arguments, passed to scipy.curve_fit / leastsq
		methods

	-	Added new support modules:
		- elecsus_methods.py
		- libs/data_proc.py
		
		elecsus_new.py contains two simplified methods for either
		calculating spectra or fitting data, and should be easier to 
		interface with external code for, e.g., batch processing / 
		fitting of data or generating 2D plots.

		data_proc.py contains methods for binning (reducing the 
		number of data points by local averaging) and moving-average 
		smoothing data traces

		both of these new modules are used by the GUI program
		
	-	Renamed the old elecsus.py module for added clarity
	
		elecsus.py --> elecsus_runcard.py
		
		This is the old method of calling elecsus with <runcard>.py files as system arguments.
		This way is now obsolete, being replaced by either the GUI or the methods contained in
		elecsus_methods.py. For backwards compatibility, the elecsus_runcard.py module allows
		the runcards to be used in the same way as before.
		
		The example runcards and example data have been moved to sub-directories, 
		/runcard and /sample_data, respectively.
		
	-	Added a new module, spectra_Efield.py
		
		This module allows calcualtion of spectra with the output of electric field vectors, rather
		than spectroscopic quantities. This should allow calculation of spectra in cells with non-uniform
		magnetic fields, by splitting the cell into sufficiently small parts that the field variation
		across any one part is negligible. Spectroscopic quantities can be calculated from the electric
		field vector by using Jones matrices.
		
		<< add example for magnetic field calculation? >>

-------------
Prerequisites
-------------

Must have the python programming language installed with the following 
packages:

	- Scipy version 0.12.0 or later
	- Numpy
	- Matplotlib
	- wxPython 2.8 (for GUI)


------------
Installation
------------

Unzip the elecsus.zip file << to update with pip install info >>

-----
Usage
-----

	1. For GUI operation:

	- After package installation, from the python interpreter type:
	
		from elecsus import elecsus_gui
		
		elecsus_gui.start()

	- In windows, double-click on the run_gui.bat file in the elecsus directory

	- Alternately, open a terminal or command-line window and navigate to the ElecSus directory. Type:

		python elecsus_gui.py
			
	2. For Runcard operation:

	- Open a terminal window and move to the directory where the files have been extracted to.

	- To run the program taking parameters from runcard.py type:

		python elecsus_runcard.py

	- To run using parameters from a particular runcard type

		python elecsus_runcard.py <run card file name>

	- So to run a the first D1 example, type

		python elecsus.py runcard_D1sample.py

	- To run the second example, type

		python elecsus.py runcard_D2sample.py
		

	3. For integration into external code:
	
	- The elecsus_methods.py module contains two methods, calculate() and fit_data(), 
	  which allow for easy integration into external codes. See the elecsus_methods.py source
	  for more details

------
Manual
------

For GUI documentation, see docs/ElecSus_GUI_UserGuide.pdf

For the ElecSus paper, go to http://dx.doi.org/10.1016/j.cpc.2014.11.023
and download the paper. It is published open-access and therefore freely available.

-------
License
-------

All the files distributed with this program are provided subject to the
Apache License, Version 2.0. A Copy of the license is provided.
