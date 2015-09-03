from setuptools import setup
import time


print "\n\n NOTE: wxPython 2.8 needs to be installed for the\n\
GUI part of this program to work! It is not currently possible\n\
to install this automatically through setuptools / pip / easy_install.\n\
\
This is included with Enthought Canopy. For Windows systems, this is all that's needed.\
For Linux systems, wxPython is not supported in Canopy and needs to be installed separately.\
To install wxPython, please visit the wxPython website:\n\
http://www.wxpython.org/download.php \n\n"
time.sleep(1)

setup(
	name='ElecSus',
	description='(Atomic Physics) Calculate the weak-probe electric susceptibility for alkali-metal vapours',
	author='James Keaveney, Mark Zentile et. al.',
	author_email='james.keaveney@durham.ac.uk',
	url='https://github.com/jameskeaveney/ElecSus',	
	version='2.0.0.rc1',
	packages=['elecsus', 'elecsus.libs', 'elecsus.runcard'],
	package_data={'elecsus':['images/elecsus_group.ico',
							'images/elecsus_t_group.ico',
							'images/elecsus-logo.png',
							'images/jqc-logo.png',
							'sample_data/SampleData_Ix_RbD2.csv',
							'sample_data/SampleData_S1_RbD2.csv',
							'sample_data/SampleDataRbD1.csv',
							'docs/ElecSus_GUI_UserGuide.pdf', 
							'docs/Manual.pdf']},
	license='Apache License, Version 2.0',
	long_description=open('README.txt').read(),
	install_requires=[
		#'wxPython == 2.8.11.0' ### wxPython is needed but can't be installed from PyPi.
		'numpy',
		'scipy >= 0.12.0',
		'matplotlib >= 1.3.1',
		],
	classifiers=[
		'Development Status :: 1 - Planning',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: Apache Software License',
		'Programming Language :: Python :: 2.7',
		'Topic :: Scientific/Engineering :: Physics']
	
	)