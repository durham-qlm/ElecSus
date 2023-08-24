> **Version Incompatibility**: There is currently a known issue with wxpython and python 3.10 that means 
> ElecSus only supports python 3.9 and below. This will be updated once these issues are resolved.


ElecSus v3.1.0(beta): Extension to Arbitrary magnetic field vectors
==============

A program to calculate the electric susceptibility of an atomic ensemble. The program is designed to model weak-probe laser spectra on the D-lines of thermal alkali metal vapour cells. The program also includes fitting routines which allow experimental parameters to be extracted from experimental spectroscopic data.

--------------------
New in version 3.0
--------------------

-	**Major overhaul** that adds much additional functionality. See the new paper (https://doi.org/10.1016/j.cpc.2017.12.001) for the full details. In brief, the additions are
    - To include calculation for an arbitrary angle between magnetic field axis and light propagation vector, 
	in part based on the publication by Rotondaro, Zhdanov and Knize [ J. Opt. Soc. Am. B 32, 12, 2507 (2015) ] and references therein.
	
    - A more general approach to light propagation, with the input polarisation more rigorously defined. 
	*This change makes this version backwards incompatible with the previous (v1.x, v2.x) versions*, since now we explicitly calculate the electric field at the exit of the medium and use Jones matrices/vectors to compute the Stokes parameters and other derived quantities. This also has the benefit of enabling the simulation of, for example: magnetic field gradients; polarising optical elements; birefringence or other optical imperfections. Previous versions are still available on the releases section of these pages, but are no longer maintained.

    - The GUI has been significantly changed to accomodate the extra magnetic field and polarisation options.

    - The old (v1.x) runcard method of calling elecsus is now dead and buried, since it's no longer compatible with the current version.

    - The fitting routines have been completely rewritten to use parameter dictionaries. We now utilise the lmfit module (https://lmfit.github.io/lmfit-py/), instead of the vanilla scipy least-square minimisation modules. Performance is broadly similar (since lmfit also runs over the top of scipy.optimize modules), but there are many advantages of this model: Firstly, all parameters can now be selected to vary during a fit or be held constant, and bounds for these parameters can be given to prevent unphysical values from being returned. In addition, the differential_evolution solver is now availablle, which is very efficient at finding the global minimum for multi-parameter fits (we leave in random_restart and simulated_annealing for now, though these might disappear in future versions as it appears differential_evolution is much better in all tested cases so far - see elecsus/tests/fitting_tests.py for examples).
		
-------------
Prerequisites
-------------

Must have the python programming language installed with the following 
packages:

- Scipy version 0.12.0 or later
- Numpy
- Matplotlib
- wxPython >= 2.8 (for GUI)
- lmfit 0.9.5 or later (https://lmfit.github.io/lmfit-py/index.html)

------------
Installation
------------

Python and required packages must be installed prior to installing ElecSus.

- Download the zip file and extract the ElecSus directory.

- For windows, there are pre-built binaries which will install ElecSus for you. Simply double click on the installer exe/msi file and follow the instructions.

- For linux-based systems, download or clone this repository and navigate to the download location in a terminal window. Install using the setup.py file by typing
	
    ```
    python setup.py install
    ```

- Note the GUI part of ElecSus is currently untested on Mac OSX!

-----
Usage
-----
1. For running 'Basic Calculations' notebook (RECOMMENDED FOR NEW USERS):
- After package installation, open 'Basic Calculations.ipynb' using Jupyter, JupyterLab etc.
- Run all cells to see all the basic calculations ElecSus can perform.
- Adjust parameters in the second cell of the notebook.
- A list of parameters and keywords can be found commented out in the first cell.

2. For GUI operation:

- After package installation, from the python interpreter type:

    ```
    from elecsus import elecsus_gui
    elecsus_gui.start()
    ```

- In windows, double-click on the run_gui.bat file in the elecsus directory

- Alternately, open a terminal or command-line window and navigate to the ElecSus directory. Type:

    ```
    python elecsus_gui.py
    ```


3. For integration into external code:

- The elecsus_methods.py module contains two methods, calculate() and fit_data(), 
  which allow for easy integration into external codes. See the elecsus_methods.py doc strings
  for more details and example usage

------
Manual
------

For GUI documentation, see the docs/ sub-folder.

For the original ElecSus paper, go to http://dx.doi.org/10.1016/j.cpc.2014.11.023
and download the paper.
For the ElecSus paper that goes with version 3, see https://doi.org/10.1016/j.cpc.2017.12.001
Both papers are published open-access and therefore freely available.

A notebook which goes processes and fits data with ElecSus can be found here:
https://github.com/durham-qlm/scope2elecsus. This is accompanied by its own published paper which is also open access and freely available:
https://doi.org/10.1088/1367-2630/ac9cfe


-------
License
-------

All the files distributed with this program are provided subject to the
Apache License, Version 2.0. A Copy of the license is provided.

-----------
Change Log
-----------
v 3.0.8
- All references to deprecated 'time' module have been removed (excluding GUI).
- New 'Notebook' folder added. It is designed to house jupyter notebooks which showcase various functions in ElecSus
- The folder currently contains 'Basic Calculations.ipynb' which calculates spectra, susceptibilities and refractive indices.

v 3.0.7
- Large speed improvement for electric field calculations.
- Fixed bug to allow data saving in python 3.
- Bug with relative paths fixed.
- Reduced number of initial points in RR fitting routine.
- Fixed compatibility issue with matplotlib versions > 3.
- Fixed the fitting test modules

V 3.0.6

- Bug in the data processing module (libs/data_proc.py) fixed
- Minor change to support a change in matplotlib v2.2
	
V 3.0.5

- Bug fix with some menu items not working properly

V 3.0.4 

- Support for python 3.x added (maintains compatibility with python 2.7)
	
V 3.0.3

- Minor fixes to GUI for file input/output not working properly and an error that stopped fitting working

V 3.0.2

- Minor fix: changed a wx.OPEN to wx.FD_OPEN that affected newer versions of wx FileDialogs not opening properly

V 3.0.1

- Minor fix: missing square-root in the solve_dielectric.py which was causing errors for the non-analytic magnetic field angles
- Updated documentation in solve_dielectric.py
- Added a new example to the /tests/ subdirectory for non-uniform magnetic field fitting (to reproduce figure 8 of the new paper)
- Minor addition: use proper fine-structure constants for all alkalis, instead of just the value for Rb

V 3.0.0

- Main program overhaul to include arbitrary angle between magnetic field axis and light propagation vector.
- Large update to the GUI to support the above change (see above)
- Initial ground state populations are now calculated via the Boltzmann factor rather than assuming equal populations
- Some housekeeping on many of the supporting files in the libs/ directory - mostly to tidy up redundant code, added more comments and examples etc.
	
V 2.2.0

- GUI version number and ElecSus version number are now the same
- Since Enthought Canopy now ships with wxPython version 3.0, GUI has been
	updated to work with this version of wxPython. All previous functionality should 
	remain the same, but a few things have changed:
    - Theory/Fit plot selections are no longer Transient Popups - they are now Dialog Boxes
	- Default plot scaling may look a bit odd when viewing on small resolution monitors -
			not sure what the real fix for this is, but as a workaround, resizing the GUI window
			should fix this
- Added ability to use experimental detuning axis on theory curve, 
	for direct comparison when outputting data (Menu > Edit > Use Experimental Detuning Axis)
- Added ability to turn off automatic scaling of axes (Menu > View > Autoscale)
- Fixed an issue where save files would not save correctly after performing a fit
- Minor fix for an issue where starting from the python interpreter would not exit properly on clicking the 'X'
- Corrected some incorrect tooltips
- Added show_versions() method to elecsus_gui.py, which shows the currently used version numbers of 
	modules important to running this program
- An alternate GUI based on the TraitsUI interface has been developed by Matthias Widman and co-authors and is available here: https://github.com/matwid/ElecSus (Note this is based on ElecSus v2.2.0)
		
V 2.1.0

- Cleaned up a lot of code in the spectra.py module. Added new methods (calc_chi, get_spectra) to spectra.py that allow users to easily create wrapper methods to return custom data types that aren't returned by default. Separated the method to calculate electric field propagation (get_Efield), for use with non-uniform B fields. Backward compatibility is preserved with the spectrum() method, which is now just a wrapper for get_spectra().
	
V 2.0.3

- Bug fix for the runcard method of using ElecSus. Now supports runcards in directories other than the local directory.
	
V 2.0.2

- Minor bug fix for the GUI - fixed an issue where the Phi plots would not be plotted
	
V 2.0.1

- Updated Eigensystem.py with correct fine structure constant for Rb.
