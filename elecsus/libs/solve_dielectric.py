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

Use analytic solutions for the 'easy' geometries - Faraday (B-field aligned with
wavevector) and Voigt (B-field orthogonal to k-vector).

[ES4 MOD] All other geometries now go through Webber's vectorised analytic
eigen-solution below.  ES3 solved those with SymPy, one detuning point at a
time; that path was removed, and with it this module's SymPy dependency -- the
symbolic derivation it performed is kept as comments in solve_diel() so the
result stays checkable.  Removed alongside it: the development scratch that had
accumulated here (null/nullOld/test_null/test_solveset/main/
calculation_time_analysis/test_equivalence), none of which was reachable from
solve_diel().
"""

import numpy as np

def solve_diel(chiL, chiR, chiZ, THETA, Bfield, verbose=False,force_numeric=False):
	''' 
	Solves the wave equation to find the two propagating normal modes of the system, 
	for a given magnetic field angle THETA. The Faraday (B along k) and Voigt
	(B perpendicular to k) geometries have closed forms; every other angle goes
	through the vectorised oblique solution further down. ([ES4] ES3 solved the
	oblique case with SymPy, one detuning point at a time.)
	
	Returns the rotation matrix to transform the coordinate system into the normal mode basis,
	and returns the two refractive index arrays.
	
	Inputs:
	
		chiL, chiR, chiZ	:	1D lists or numpy arrays, of length N, that are the frequency-dependent electric susceptibilities
								of the sigma+, sigma- and pi manifolds, as spectra.calc_chi returns them. About B, chiL is
								the susceptibility of the [1,+i]/sqrt2 component and chiR that of [1,-i]/sqrt2.
		THETA				:	Float, Magnetic field angle in radians: B lies along (sin THETA, 0, cos THETA), in the x-z plane
		Bfield				:	Float, Magnitude of applied magnetic field. [ES4 FIX] Unused: below 0.01 G it
								used to select the Faraday solution at any angle but Voigt. Kept for callers' signatures
		
	Options:
		
		verbose			:	Boolean to output more print statements (timing reports mostly)
		force_numeric	:	If True, every angle goes through the oblique solution, even where the Faraday or Voigt form applies
	
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
	# [ES4 FIX] Both geometry tests are periodic in THETA. The
	# Voigt test used to catch THETA = pi/2 (mod 2 pi) only, so B along -x
	# (THETA = 3pi/2 or -pi/2, and a reversed pass at pi/2 in ES3's
	# convention) fell through to the oblique solver. Its Ez = 1 eigenvector
	# normalisation is singular when B is perpendicular to k, and S0 came out
	# at ~1e10. The Faraday test likewise missed THETA = 2 pi and -pi.
	cos_theta = np.cos(THETA)
	# [ES4 FIX] The oblique eigenvectors below are written for B
	# along (-sin THETA, 0, cos THETA) (Webber's characteristic matrix). B lies
	# along (+sin THETA, 0, cos THETA), so sin THETA enters with its sign
	# reversed. The two frames differ by a turn of pi about z, which reverses
	# the transverse part of each eigenvector relative to its longitudinal part
	# and leaves the indices unchanged. (The Voigt form is written for the lab
	# frame; see the note there.)
	sin_theta = -np.sin(THETA)
	if (abs(cos_theta) < 1e-4) and (not force_numeric):
		# ANALYTIC SOLNS FOR VOIGT
		if verbose: print('Voigt - analytic')

		# solutions for elements of the dielectric tensor:
		ex = 0.5 * (2. + chiL + chiR)
		exy = 0.5j * (chiR - chiL)
		ez = 1.0 + chiZ

		# refractive indices to propagate
		n1 = np.sqrt(ex + exy**2/ex)
		n2 = np.sqrt(ez)

		# [ES4 FIX] With B along -x (sin THETA < 0) chiL and chiR trade places,
		# so exy changes sign: the indices are unchanged, and the longitudinal
		# part of mode 1 changes sign relative to its transverse part.
		# [ES4 FIX] For B along +x the z row of the wave equation,
		# [eps - n^2 (1 - zz)] u = 0, reads -exy u_y + ex u_z = 0, so
		# u_y/u_z = +ex/exy: the sign of sin THETA, not of sin_theta (which is
		# -sin THETA, for the oblique solver's frame). ES3's vector had the
		# opposite sign, so E_out[2] (non-orthogonal path) was exactly reversed
		# at THETA = pi/2 and 3 pi/2 and flipped at the edge of this branch.
		ev1 = [np.zeros(len(ex)),np.sign(np.sin(THETA))*ex/exy,np.ones(len(ex))]
		ev2 = [np.ones(len(ex)),np.zeros(len(ex)),np.zeros(len(ex))]
		ev3 = [np.zeros(len(ex)),np.zeros(len(ex)),np.ones(len(ex))]
		
		RotMat = np.array([ev1,ev2,ev3])
		
		if verbose:
			print('Shortcut:')
			print((RotMat.shape))
			print((n1.shape))
			print((n2.shape))
		
	# [ES4 FIX] Below 0.01 G this branch used to take every angle
	# but Voigt, modelling the whole field along +z or -z. The circular
	# birefringence is set by B cos THETA, so an oblique field's Faraday
	# rotation came out 1/|cos THETA| too large, and the spectra jumped at
	# 0.01 G. The oblique solver is accurate at small fields (to 1e-9 against
	# a first-principles propagation down to 1e-6 G, and at Bfield = 0, which
	# calc_chi models as 1e-4 G), so the geometry alone now decides.
	elif (abs(sin_theta) < 1e-4) and (not force_numeric):
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

		# [ES4 FIX] aligned means THETA = 0 (mod 2 pi), not just THETA = 0.
		# [ES4 FIX] B along -z takes the swapped order. (Before 4.0.0b1
		# this also took a field below 0.01 G at an oblique angle, along the
		# nearer of +z and -z.)
		if cos_theta > 0:
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

		ct  = cos_theta
		st  = sin_theta
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
