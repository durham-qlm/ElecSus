from sympy import Symbol, solve, cos, sin, pi, simplify, eye, powsimp, powdenest
from sympy.matrices import det, Matrix

import numpy as np
import scipy.linalg as la
import scipy

import time

from FundamentalConstants import e0
''' Use symbolic python to solve for the roots of n-squared in the dielectric tensor matrix '''

def solve_diel(chiL, chiR, chiZ, THETA):


	print 'B-field angle (rad): ',THETA
	
	stt = time.clock()
	
	# make chiL,R,Z arrays if not already
	chiL = np.array(chiL)
	chiR = np.array(chiR)
	chiZ = np.array(chiZ)
	
	theta = Symbol('theta',real=True)
	n_sq = Symbol('n_sq')

	e_x = Symbol('e_x')
	e_xy = Symbol('e_xy')
	e_z = Symbol('e_z')

	# General form of the dielectric tensor
	DielMat = Matrix (( 	[(e_x - n_sq)*cos(theta), e_xy, e_x*sin(theta)],
								[-e_xy * cos(theta), e_x - n_sq, -e_xy*sin(theta)],
								[(n_sq - e_z)*sin(theta), 0, e_z*cos(theta)] 			))

	et1 = time.clock() - stt
	
	# Substitute in angle
	#if THETA == np.pi/2:
	#	THETA = pi/2 # sympy equivalent
		
	DielMat_sub = DielMat.subs(theta, pi*THETA/np.pi)
	
	
	et2 = time.clock() - stt
	
	#### Escape the slow loop for analytic (Faraday and Voigt) cases
	
	if THETA == np.pi:
		# ANALYTIC SOLNS FOR VOIGT
		solns = [0,0]
		solns[0] = e_x + e_xy**2/e_x
		solns[1] = e_z
	elif THETA == 0:
		# ANALYTIC SOLNS FOR FARADAY
		solns = [0,0]
		solns[0] = e_x + 1.j*e_xy
		solns[1] = e_x - 1.j*e_xy
	else:
		pass
	
	# Find solutions for complex indices for a given angle
	solns = solve(det(DielMat_sub), n_sq)
	
	# Find first eigenvector
	DielMat_sub1 = DielMat_sub.subs(n_sq, solns[0])
	#ev1 = np.zeros((len(chiL),3),dtype='complex')
	n1 = np.zeros(len(chiL),dtype='complex')
	# Find second eigenvector
	DielMat_sub2 = DielMat_sub.subs(n_sq, solns[1])
	#ev2 = np.zeros((len(chiL),3),dtype='complex')
	n2 = np.zeros(len(chiL),dtype='complex')
	
	RotMat = np.zeros((3,3,len(chiL)),dtype='complex')
	
	et3 = time.clock() - stt
	
	print 'setup time:', et1, et1
	print 'solve nsq: ', et2, et2-et1
	print 'sub in: ', et3, et3-et2
	
	# loop over all elements of chiL,R,Z ---- !!!!!!!!!! massive impact on speed !!!!!!!!!!
	for i, (cL, cR, cZ) in enumerate(zip(chiL,chiR,chiZ)):
		print 'i: ',i
		
		#time diagnostics
		st = time.clock()
		
		DielMat_sub1a = DielMat_sub1.subs(e_x, 0.5*(2.+cL+cR))
		DielMat_sub1a = DielMat_sub1a.subs(e_xy, 0.5j*(cR-cL))
		DielMat_sub1a = DielMat_sub1a.subs(e_z, (1.0+cZ))
		
		et1 = time.clock() - st
		
	
		DM = np.array(DielMat_sub1a.evalf())
		DMa = np.zeros((3,3),dtype='complex')
		for ii in range(3):
			for jj in range(3):
				DMa[ii,jj] = np.complex128(DM[ii,jj])
				
		et2 = time.clock() - st
	
		# use scipy to find eigs
		
		#print null(DMa).T[0]
		ev1 = null(DMa).T[0]
		
		et3 = time.clock() - st
		
		# sub in for ref. index
		n1soln = solns[0].subs(e_x, 0.5*(2.+cL+cR))
		n1soln = n1soln.subs(e_xy, 0.5j*(cR-cL))
		n1soln = n1soln.subs(e_z, (1.0+cZ))
		
		n1[i] = np.sqrt(np.complex128(n1soln.evalf()))
		
		et4 = time.clock() - st
		
		DielMat_sub2a = DielMat_sub2.subs(e_x, 0.5*(2.+cL+cR))
		DielMat_sub2a = DielMat_sub2a.subs(e_xy, 0.5j*(cR-cL))
		DielMat_sub2a = DielMat_sub2a.subs(e_z, (1.0+cZ))
	
		et5 = time.clock() - st
		
		DM = np.array(DielMat_sub2a.evalf())
		DMa = np.zeros((3,3),dtype='complex')
		for ii in range(3):
			for jj in range(3):
				DMa[ii,jj] = np.complex128(DM[ii,jj])
	
		et6 = time.clock() - st
		
		# use numpy to find eigs
		ev2 = null(DMa).T[0]

		et7 = time.clock() - st
		
		# sub in for ref. index
		n2soln = solns[1].subs(e_x, 0.5*(2.+cL+cR))
		n2soln = n2soln.subs(e_xy, 0.5j*(cR-cL))
		n2soln = n2soln.subs(e_z, (1.0+cZ))
		
		n2[i] = np.sqrt(np.complex128(n2soln.evalf()))

		et8 = time.clock() - st
		
		RotMat[:,:,i] = [ev1, ev2, [0,0,1]]
		
		print 'Time elapsed:'
		print '1 Sub into matrix:', et1, et1
		print '1 Eval and convert to numpy:', et2, et2-et1
		print '1 Get eig vector:', et3, et3-et2
		print '1 Get ref index array:', et4, et4-et3
		print '2 Sub into matrix:', et5, et5 - et4
		print '2 Eval and convert to numpy:', et6, et6 - et5
		print '2 Get eig vector:', et7, et7 - et6
		print '2 Get ref index array:', et8, et8 - et7
		
	'''
	###### all numpy version (analytic soln only)
	EX = e0*0.5*(2.+chiL+chiR)
	EXY = e0*0.5j*(chiR-chiL)
	EZ = e0*(1.0+chiZ)

	print 'Numpy version...........'
	# Voigt
	THETA = np.pi/2
	NSQ = EX + EXY**2/EX
	#NSQ = EZ
	# Faraday
	#THETA = 0
	#NSQ = EX + 1.j*EXY
	#NSQ = EX - 1.j*EXY
	DielMat = np.matrix ([ 	[(EX - NSQ)*np.cos(THETA), EXY, EX*np.sin(THETA)],
								[-EXY*np.cos(THETA), EX - NSQ, -EXY*np.sin(THETA)],
								[(NSQ - EZ)*np.sin(THETA), 0, EZ*np.cos(THETA)] 			])
	
	print DielMat
	
	ev = null(DielMat)
	print ev/abs(ev[2])
	print type(ev[0][0])
	print 'Null::', ev
	print EX/EXY
	scl = 1.j/abs(EX/EXY)
	print ev/scl.conjugate()
	#ee, ev = la.eig(DielMat)
	#print 'EigenVals:'
	#print ee[0], ee[1], ee[2]
	#print 'EigenVects:'
	#print ev[0]
	#print ev[1]
	#print ev[2]
	'''
	
	
	############ MAKE THIS WORK ON ARRAYS OF CHI_R/L/Z
	return RotMat, n1, n2
	

def null(A, eps=1e-15):
	# Taken from http://stackoverflow.com/questions/5889142/python-numpy-scipy-finding-the-null-space-of-a-matrix
	u, s, vh = la.svd(A)
	null_mask = (s <= eps)
	null_space = scipy.compress(null_mask, vh, axis=0)
	return scipy.transpose(null_space)	
	
def main():
	""" testing """
	import spectra as sp
	p_dict = {'Bfield':700,'rb85frac':1,'Btheta':88*np.pi/180,'Bphi':0*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	chiL,chiR,chiZ = sp.calc_chi(np.linspace(-3500,3500,10),p_dict)
	
	#print 'ez: ',chiZ + 1 # ez / e0
	#print 'ex: ',0.5*(2+chiL+chiR) # ex / e0
	#print 'exy: ',0.5j*(chiR-chiL) # exy / e0
	
	RotMat, n1, n2 = solve_diel(chiL,chiR,chiZ,88*np.pi/180)
	print RotMat.shape

def calculation_time_analysis():
	import spectra as sp
	p_dict = {'Bfield':700,'rb85frac':1,'Btheta':88*np.pi/180,'Bphi':0*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	chiL,chiR,chiZ = sp.calc_chi([-3500],p_dict)
	
	for angle in [0, np.pi/32, np.pi/16, np.pi/8, np.pi/4, np.pi/2]:
		print 'Angle (degrees): ',angle*180/np.pi
		RotMat, n1, n2 = solve_diel(chiL,chiR,chiZ,angle)