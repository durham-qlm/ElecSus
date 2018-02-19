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

"""
Functions to calculate the z-projection matrices.

Calculates the full size matrices for the projection on the quantization 
axis for electron spin, orbital angular momentum, nuclear spin, and the 
coupled angular momentum F. Essentially takes results from jz and puts
them in the full Hilbert space.

Calls jz from ang_mom.

Last updated 2018-02-19 JK
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)


from numpy import identity
from scipy.linalg import kron
from .ang_mom import jz

def sz(L,S,I):
    Sz=jz(S)
    gL=int(2*L+1)
    Li=identity(gL)
    gI=int(2*I+1)
    Ii=identity(gI)
    sz=kron(kron(Li,Sz),Ii)
    return sz

def lz(L,S,I):
    gS=int(2*S+1)
    Si=identity(gS)
    Lz=jz(L)
    gI=int(2*I+1)
    Ii=identity(gI)
    lz=kron(kron(Lz,Si),Ii)
    return lz

def Iz(L,S,I):
    gS=int(2*S+1)
    gL=int(2*L+1)
    Si=identity(gS)
    Li=identity(gL)
    Iz_num=jz(I)
    Iz=kron(kron(Li,Si),Iz_num)
    return Iz

def fz(L,S,I):
    gS=int(2*S+1)
    Sz=jz(S)
    Si=identity(gS)
    gL=int(2*L+1)
    Lz=jz(L)
    Li=identity(gL)
    gJ=gL*gS
    Jz=kron(Lz,Si)+kron(Li,Sz)
    Ji=identity(gJ)
    gI=int(2*I+1)
    Iz=jz(I)
    Ii=identity(gI)
    Fz=kron(Jz,Ii)+kron(Ji,Iz)
    return Fz
