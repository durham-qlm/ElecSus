# Copyright 2014 M. A. Zentile, J. Keaveney, L. Weller, D. Whiting,
# C. S. Adams and I. G. Hughes.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
	version='2.0.1',
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
	long_description=open('README.md').read(),
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