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

setup(
    name='ElecSus',
    description='(Atomic Physics) Calculate the weak-probe electric susceptibility for alkali-metal vapours',
    author='James Keaveney et. al.',
    author_email='jkeaveney23@gmail.com',
    url='https://github.com/Quantum-Light-and-Matter/ElecSus',
    version='3.0.8',
    packages=['elecsus', 'elecsus.libs', 'elecsus.tests'],
    package_data={'elecsus': ['images/elecsus_group.ico',
                              'images/elecsus_t_group.ico',
                              'images/elecsus-logo.png',
                              'images/jqc-logo.png',
                              'docs/ElecSus_GUI_UserGuide.pdf']},
    license='Apache License, Version 2.0',
    long_description=open('README.md').read(),
    install_requires=[
        'wxPython >= 4.1.0',
        'numpy',
        'scipy >= 0.12.0',
        'matplotlib >= 1.3.1',
        'lmfit >= 0.9.5',
        'psutil',
        'cycler',
        'sympy'
    ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3.x',
        'Topic :: Scientific/Engineering :: Physics'
    ]
)
