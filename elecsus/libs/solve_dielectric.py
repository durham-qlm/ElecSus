# Copyright 2017-2026 J. Keaveney and S. A. Wrathmall

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
Solve the dielectric tensor for the roots of the complex refractive index by setting determinant to 0

Use analytic solutions for the 'easy' geometries - Faraday (B-field aligned with wavevector) and Voigt (B-field orthogonal to k-vector)
Use sympy to calculate solutions for all other non-trivial geometries.
Since the solutions for the non-trivial geometries depend on the susceptibility, array operations don't completely work, so it's *much* slower to calculate

Last updated 2018-07-04 MAZ
"""
# py 2.7 compatibility


from sympy import Symbol, cos, sin, pi, simplify, eye, powsimp, powdenest, lambdify, solve, solveset
#from sympy.solvers.solveset import linsolve as solve
from sympy.matrices import det, Matrix

import numpy as np
import scipy.linalg as la
from scipy.linalg import qr
import scipy

from FundamentalConstants import e0

def solve_diel(chiL, chiR, chiZ, THETA, Bfield, verbose=False,force_numeric=False):
	''' 
	Solves the wave equation to find the two propagating normal modes of the system, 
	for a given magnetic field angle THETA. For the general case, use symbolic python to 
	solve for the roots of n-squared.
	(Escapes this slow approach for the two analytic cases for the Voigt and Faraday geometries)
	
	Returns the rotation matrix to transform the coordinate system into the normal mode basis,
	and returns the two refractive index arrays.
	
	Inputs:
	
		chiL, chiR, chiZ	:	1D lists or numpy arrays, of length N, that are the frequency-dependent electric susceptibilities
		THETA				:	Float, Magnetic field angle in radians
		Bfield				:	Float, Magnitude of applied magnetic field (skips slow approach if magnetic field is very close to zero)
		
	Options:
		
		verbose			:	Boolean to output more print statements (timing reports mostly)
		force_numeric	:	If True, forces all angles to go through the numeric approach, rather than escaping for the analytic cases (THETA=0, THETA=pi/2...)
	
	Outputs:
		RotMat	:	Rotation matrix to transform coordinate system, dimensions (3, 3, N)
		n1		:	First solution for refractive index, dimensions (N)
		n2		:	Second solution for refractive index, dimensions (N)
	
	'''


	if verbose: 
		print(('B-field angle (rad, pi rad): ',THETA, THETA/np.pi))
	
	# make chiL,R,Z arrays if not already
	chiL = np.array(chiL)
	chiR = np.array(chiR)
	chiZ = np.array(chiZ)
	
	#### Escape the slow loop for analytic (Faraday and Voigt) cases
	## For these analytic cases we can use array operations and it is therefore
	## much faster to compute
	if (abs(THETA%(2*np.pi) - np.pi/2) < 1e-4) and (not force_numeric):
		# ANALYTIC SOLNS FOR VOIGT
		if verbose: print('Voigt - analytic')
		
		# solutions for elements of the dielectric tensor:
		ex = 0.5 * (2. + chiL + chiR)
		exy = 0.5j * (chiR - chiL)
		ez = 1.0 + chiZ

		# refractive indices to propagate
		n1 = np.sqrt(ex + exy**2/ex)
		n2 = np.sqrt(ez)
		
		ev1 = [np.zeros(len(ex)),ex/exy,np.ones(len(ex))]
		ev2 = [np.ones(len(ex)),np.zeros(len(ex)),np.zeros(len(ex))]
		ev3 = [np.zeros(len(ex)),np.zeros(len(ex)),np.ones(len(ex))]
		
		RotMat = np.array([ev1,ev2,ev3])
		
		if verbose:
			print('Shortcut:')
			print((RotMat.shape))
			print((n1.shape))
			print((n2.shape))
		
	elif ((abs(THETA) < 1e-4) or ((abs(THETA - np.pi)) < 1e-4) or abs(Bfield)<1e-2)  and (not force_numeric): ## Use Faraday geometry if Bfield is very close to zero
		# ANALYTIC SOLNS FOR FARADAY
		#if verbose: 
		if verbose: print('Faraday - analytic TT')
		
		ex = 0.5 * (2. + chiL + chiR)
		exy = 0.5j * (chiR - chiL)
		e_z = 1.0 + chiZ
		
		n1 = np.sqrt(ex + 1.j*exy)
		n2 = np.sqrt(ex - 1.j*exy)

		ev1 = np.array([-1.j*np.ones(len(ex)),np.ones(len(ex)),np.zeros(len(ex))]) 
		ev2 = np.array([1.j*np.ones(len(ex)),np.ones(len(ex)),np.zeros(len(ex))])
		ev3 = [np.zeros(len(ex)),np.zeros(len(ex)),np.ones(len(ex))]

		if (abs(THETA) < 1e-4):
			RotMat = np.array([ev1,ev2,ev3])
		else:
			#if anti-aligned, swap the two eigenvectors
			RotMat = np.array([ev2,ev1,ev3])
			
		if verbose:
			print('Shortcut:')
			print((RotMat.shape))
			print((n1.shape))
			print((n2.shape))

	else:
		# [ES4 MOD] Oblique solver: replaced slow SymPy + Python for-loop with Webber's fully
		# vectorised analytic solution (A. Webber, Level 4 Project, Durham 2019, Appendix A & B;
		# benchmarked by Logue 2023 PhD thesis, Appendix B, Fig. B.7: ~128× faster).
		# ES3 method produced unphysical transmission > 100% at high fields due to numerical
		# instability. Eigenvectors further corrected 2026-06-28 (Webber Appendix-B formulas
		# gave wrong vectors near Voigt; replaced with direct characteristic-matrix solve).
		# Original SymPy code is preserved at the bottom of this block for reference.
		if verbose: print('Oblique - Webber analytic [ES4]')

		# Dielectric tensor elements (vectorised over all detuning points)
		ex  = 0.5 * (2. + chiL + chiR)
		exy = 0.5j * (chiR - chiL)
		ez  = 1.0 + chiZ

		ct  = np.cos(THETA)
		st  = np.sin(THETA)
		ct2 = ct**2
		st2 = st**2

		# Refractive indices (Webber eq. 7 / Logue eq. B.3-B.4)
		A = 2.*ex*ez + (ex**2 - ex*ez + exy**2)*st2
		B = (ex**2 - ex*ez + exy**2)**2 * st2**2 - 4.*exy**2 * ez**2 * ct2
		C = 2.*(ex*st2 + ez*ct2)

		sqrtB = np.sqrt(B.astype(complex))
		n1 = np.sqrt((A + sqrtB) / C)
		n2 = np.sqrt((A - sqrtB) / C)

		# Eigenvectors: solve the characteristic-matrix system derived from the
		# wave equation for propagation along z in a medium with B at angle THETA.
		# For each mode with eigenvalue n², the characteristic matrix is:
		#   [(ex-n²)ct,  exy,     ex·st ]
		#   [-exy·ct,    ex-n²,  -exy·st]
		#   [(n²-ez)st,  0,       ez·ct ]
		# Setting Ez=1 (Webber z=1 normalisation), rows 0 and 1 give a 2×2 system
		# for [Ex, Ey].  This is fully vectorised over all N frequency points.
		# Replaces Webber Appendix-B analytic formulas which give incorrect
		# eigenvectors near the Voigt geometry (large residuals confirmed numerically).
		def _solve_ev(n_sq):
			"""Return [Ex, Ey] for eigenvector with eigenvalue n_sq, shape (N,2)."""
			Amat = np.empty((len(n_sq), 2, 2), dtype='complex')
			Amat[:, 0, 0] = (ex - n_sq) * ct
			Amat[:, 0, 1] = exy
			Amat[:, 1, 0] = -exy * ct
			Amat[:, 1, 1] = ex - n_sq
			bvec = np.column_stack([-ex * st, exy * st])
			# numpy ≥2.0 gufunc signature (m,m),(m,n)->(m,n) requires b to carry an
			# explicit column dimension; squeeze it back to (N,2) on return.
			return np.linalg.solve(Amat, bvec[:, :, np.newaxis])[:, :, 0]  # (N, 2)

		ev_a = _solve_ev(n1**2)   # [Ex, Ey] for mode n1
		ev_b = _solve_ev(n2**2)   # [Ex, Ey] for mode n2
		e1a, e2a = ev_a[:, 0], ev_a[:, 1]
		e1b, e2b = ev_b[:, 0], ev_b[:, 1]

		# e3 = 1 for both modes (Webber's normalisation)
		zeros = np.zeros(len(chiL), dtype='complex')
		ones  = np.ones(len(chiL),  dtype='complex')

		# Assemble RotMat (3, 3, N) — no Python for-loop
		# Convention (matches Faraday/Voigt branches):
		#   RotMat[mode, component, freq_idx]
		# Row 0 = eigenvector for mode n1 (a): [e1a, e2a, 1]   (x, y, z components)
		# Row 1 = eigenvector for mode n2 (b): [e1b, e2b, 1]
		# Row 2 = z-direction placeholder:     [0,   0,   1]
		RotMat = np.zeros((3, 3, len(chiL)), dtype='complex')
		RotMat[0, 0, :] = e1a   # mode a, x-component
		RotMat[0, 1, :] = e2a   # mode a, y-component
		RotMat[0, 2, :] = ones  # mode a, z-component (e3=1)
		RotMat[1, 0, :] = e1b   # mode b, x-component
		RotMat[1, 1, :] = e2b   # mode b, y-component
		RotMat[1, 2, :] = ones  # mode b, z-component (e3=1)
		RotMat[2, 0, :] = zeros
		RotMat[2, 1, :] = zeros
		RotMat[2, 2, :] = ones

		# [ES4 REF] Original ES3 SymPy + for-loop code below — retained for regression comparison.
		# Original: if verbose: print('Non-analytic angle.. This will take a while...')  ##### THIS IS THE ONE THAT's WRONG....
		# theta = Symbol('theta',real=True); n_sq = Symbol('n_sq')
		# e_x = Symbol('e_x'); e_xy = Symbol('e_xy'); e_z = Symbol('e_z')
		# DielMat = Matrix((
		#     [(e_x - n_sq)*cos(theta), e_xy, e_x*sin(theta)],
		#     [-e_xy * cos(theta), e_x - n_sq, -e_xy*sin(theta)],
		#     [(n_sq - e_z)*sin(theta), 0, e_z*cos(theta)] ))
		# DielMat_sub = DielMat.subs(theta, pi*THETA/np.pi)
		# solns = solve(det(DielMat_sub), n_sq)
		# n1 = np.sqrt(nsub1(0.5*(2.+chiL+chiR), 0.5j*(chiR-chiL), (1.0+chiZ)))
		# n2 = np.sqrt(nsub2(0.5*(2.+chiL+chiR), 0.5j*(chiR-chiL), (1.0+chiZ)))
		# for i, (cL, cR, cZ) in enumerate(zip(chiL,chiR,chiZ)):
		#     DMaNP = Dsub1(0.5*(2.+cL+cR), 0.5j*(cR-cL), (1.0+cZ))
		#     ev1 = null(DMaNP).T
		#     DMaNP = Dsub2(0.5*(2.+cL+cR), 0.5j*(cR-cL), (1.0+cZ))
		#     ev2 = null(DMaNP).T
		#     RotMat[:,:,i] = [ev1, ev2, [0,0,1]]

	if verbose: print('SD done')
	return RotMat, n1, n2
	

def null(A,tol=1e-6):
	ee, ev = la.eig(A)
	
	#for E,V in zip(ee,ev.T):
	#	print 'Eigs:',abs(E), '\t', E#, '\t', V
	#print '\n'
	
	z = list(zip(ee,ev.T))
	zs = sorted(z, key=lambda f: abs(f[0])) # sort by absolute value of eigenvectors
	ees, evs = list(zip(*zs))
	
	#for E,V in zip(ee,ev):
	#	print abs(E), '\t', E, '::', V
	
	if abs(ees[0]<tol):
		return evs[0].T
	else:
		print('No null eigenvector found! List of eigenvalules:')
		for E,V in zip(ee,ev.T):
			print(('Eigs:',abs(E), '\t', E, '\n\t', V))
		print('\n')
		return 0
		
def test_null():
	A = np.matrix([[2,3,5],[-4,2,3],[0,0,0]])
	SymA = Matrix(A)
	
	nv = null(A)
	nvold = nullOld(A)
	
	print((nv.T))
	print((nvold.T[0]))
	print((SymA.nullspace()[0].evalf()))
	
	print((A * nv))
	
def test_solveset():
	x = Symbol('x')
	A = Matrix([[x,2,x*x],[4,5,x],[x,8,9]])
	
	solns = solve(det(A), x)
	solns_set = list(solveset(det(A), x))
	
	print(solns)
	print('\n')
	print(solns_set)
	
	print('\n\n\n')
	print((solns[0]))
	print('\n')
	print((solns_set[0]))
	
	soln_sub = solns[0].subs(x, 1)
	solnset_sub = solns_set[0].subs(x, 1)
	
	s1 = soln_sub.evalf()
	s1set = solnset_sub.evalf()
	
	s2set = solns_set[1].subs(x, 1).evalf()
	
	print(s1)
	print(s1set)
	print(s2set)

def nullOld(A, eps=1e-14):
	""" Find the null eigenvector x of matrix A, such that Ax=0"""
	# Taken with gratitude from http://stackoverflow.com/questions/5889142/python-numpy-scipy-finding-the-null-space-of-a-matrix
	u, s, vh = la.svd(A)
	null_mask = (s <= eps)
	null_space = scipy.compress(null_mask, vh, axis=0)
	return scipy.transpose(null_space)	


'''
def null(A, atol=1e-15, rtol=0):
    """Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A : ndarray
        A should be at most 2-D.  A 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    ns : ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.
    """

    A = np.atleast_2d(A)
    u, s, vh = la.svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    print nnz
    ns = vh[nnz:].conj().T
    return ns
'''
	
def main():
	""" General test method """
	from . import spectra as sp
	p_dict = {'Bfield':700,'rb85frac':1,'Btheta':88*np.pi/180,'Bphi':0*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	chiL,chiR,chiZ = sp.calc_chi(np.linspace(-3500,3500,10),p_dict)
	
	#print 'ez: ',chiZ + 1 # ez / e0
	#print 'ex: ',0.5*(2+chiL+chiR) # ex / e0
	#print 'exy: ',0.5j*(chiR-chiL) # exy / e0
	
	RotMat, n1, n2 = solve_diel(chiL,chiR,chiZ,88*np.pi/180)
	print((RotMat.shape))

def calculation_time_analysis():
	""" Test method for looking at timing performance """
	from . import spectra as sp
	p_dict = {'Bfield':700,'rb85frac':1,'Btheta':88*np.pi/180,'Bphi':0*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	chiL,chiR,chiZ = sp.calc_chi([-3500],p_dict)
	
	for angle in [0, np.pi/32, np.pi/16, np.pi/8, np.pi/4, np.pi/2]:
		print(('Angle (degrees): ',angle*180/np.pi))
		RotMat, n1, n2 = solve_diel(chiL,chiR,chiZ,angle)
		
def test_equivalence():
	""" Test numeric vs analytic solutions """
	
	from . import spectra as sp
	
	#analytic
	p_dict = {'Bfield':15000,'rb85frac':1,'Btheta':0*np.pi/180,'Bphi':0*np.pi/180,'lcell':1e-3,'T':84,'Dline':'D2','Elem':'Rb'}
	chiL1,chiR1,chiZ1 = sp.calc_chi([-18400],p_dict)
	RotMat1, n11, n21 = solve_diel(chiL1,chiR1,chiZ1,0,150,force_numeric=False)
	
	#numeric
	chiL2, chiR2, chiZ2 = chiL1, chiR1, chiZ1
	#chiL2,chiR2,chiZ2 = sp.calc_chi([-18400],p_dict)
	RotMat2, n12, n22 = solve_diel(chiL2,chiR2,chiZ2,0,150,force_numeric=True)
	
	print('RM 1')
	print(RotMat1)

	print('RM 2')
	print(RotMat2)	
	
	print('n1_1 (analytic)')
	print(n11)
	print('n1_2')
	print(n12)
	print('n2_1 (analytic)')
	print(n21)
	print('n2_2')
	print(n22)
	
	print('chi1')
	print((chiL1, chiR1, chiZ1))

	print('chi2')
	print((chiL2, chiR2, chiZ2))
	
if __name__ == '__main__':
	test_equivalence()
