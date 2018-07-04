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

"""return Jx, Jy, Jz matrices for an angular momentum j 

Calculates and returns the angular momentum matrices 
for x, y and z projections.

Jz called by fs_hfs.py and sz_lsi.py

Last updated 2018-07-04 MAZ
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)

from numpy import transpose,dot
import ang_mom_p

def jx(jj):
    jp=ang_mom_p.jp(jj)
    jm=transpose(jp)
    jx=0.5*(jp+jm)
    return jx

def jy(jj):
    jp=ang_mom_p.jp(jj)
    jm=transpose(jp)
    jy=0.5j*(jm-jp)
    return jy

def jz(jj):
    jp=ang_mom_p.jp(jj)
    jm=transpose(jp)
    jz=0.5*(dot(jp,jm)-dot(jm,jp))
    return jz
