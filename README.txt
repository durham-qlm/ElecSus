ElecSus
--------------

A program to calculate the electric susceptibility of an atomic ensemble.
The program is designed to model weak probe laser spectra on the D-lines
of thermal alkali metal vapour cells. The program also includes fitting 
routines which allow experimental parameters to be extracted from 
experimental spectroscopic data.

An academic paper detailing the physics behind this code is available from
http://dx.doi.org/10.1016/j.cpc.2014.11.023
or the pre-print is available here
http://arxiv.org/abs/1409.1873

Authors
-------

Mark A. Zentile
James Keaveney
Lee Weller
Daniel J. Whiting
Charles S. Adams
Ifan G. Hughes


Prerequisites
-------------

Must have python installed with the following 
packages:

- Scipy version 0.12.0 or later
- Numpy
- Matplotlib

Installation
------------

Unpack the archived file (.zip or .tar.gz)

Usage
-----

- Open a terminal/command line window and move to the directory where the files
  have been extracted to.

- to run the program taking parameters from runcard.py type

python elecsus.py

- To run using parameters from a particular runcard type

python elecsus.py (run card file name)

- So to run a the first D1 example, type

python elecsus.py runcard_D1sample.py

- To run the second example, type

python elecsus.py runcard_D2sample.py

- Note that ElecSus can be run from an ipython terminal using

%run elecsus.py

License
-------

All the files distributed with this program are provided subject to the
Apache License, Version 2.0. A Copy of the license is provided.
