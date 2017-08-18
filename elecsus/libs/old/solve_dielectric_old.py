from sympy import Symbol, solve, cos, sin, pi, simplify, eye, powsimp, powdenest
from sympy.matrices import det, Matrix

import numpy as np
import scipy.linalg as la
import scipy

from FundamentalConstants import e0
''' Use symbolic python to solve for the roots of n-squared in the dielectric tensor matrix '''

def solve_diel(chiL, chiR, chiZ, THETA):

	
	theta = Symbol('theta',real=True)
	n_sq = Symbol('n_sq')

	e_x = Symbol('e_x')
	e_xy = Symbol('e_xy')
	e_z = Symbol('e_z')

	DielMat = Matrix (( 	[(e_x - n_sq)*cos(theta), e_xy, e_x*sin(theta)],
								[-e_xy * cos(theta), e_x - n_sq, -e_xy*sin(theta)],
								[(n_sq - e_z)*sin(theta), 0, e_z*cos(theta)] 			))

#	DielMat = Matrix (( 	[e_x - n_sq*cos(theta)**2, e_xy, n_sq*sin(theta)*cos(theta)],
	#							[-e_xy, e_x - n_sq, 0],
		#						[n_sq * sin(theta) * cos(theta), 0, e_z - n_sq*sin(theta)**2]			))
								
	

	DielMat_sub = DielMat.subs(theta, THETA)
	#DielMat_sub = DielMat_sub.subs(e_x, 0.5*(2.+chiL+chiR))
	#DielMat_sub = DielMat_sub.subs(e_xy, 0.5j*(chiR-chiL))
	#DielMat_sub = DielMat_sub.subs(e_z, 1.0+chiZ)
								
	solns = solve(det(DielMat_sub), n_sq)
	
	#print powdenest(solns[0],force=True)
	
	#solns[0] = e_x + e_xy**2/e_x
	#solns[1] = e_z
	
	#solns[0] = e_x + 1.j*e_xy
	#solns[1] = e_x - 1.j*e_xy

	DielMat_sub1 = DielMat_sub.subs(n_sq, solns[0])
	#print DielMat_sub1
	DielMat_sub1 = DielMat_sub1.subs(e_x, e0*0.5*(2.+chiL+chiR))
	DielMat_sub1 = DielMat_sub1.subs(e_xy, e0*0.5j*(chiR-chiL))
	DielMat_sub1 = DielMat_sub1.subs(e_z, e0*(1.0+chiZ))
	
	#print DielMat_sub1
	#print 'QRsolve:', DielMat_sub.QRsolve(0)
	
	#print '\nEigenvalues (1):'
	#print DielMat_sub1.eigenvals()
	
	#print 'Eigenvectors (1):'
	#print DielMat_sub1.eigenvects()[0]
	
	DielMat_sub2 = DielMat_sub.subs(n_sq, solns[1])
	#print DielMat_sub
	#print DielMat_sub2
	DielMat_sub2 = DielMat_sub2.subs(e_x, e0*0.5*(2.+chiL+chiR))
	DielMat_sub2 = DielMat_sub2.subs(e_xy, e0*0.5j*(chiR-chiL))
	DielMat_sub2 = DielMat_sub2.subs(e_z, e0*(1.0+chiZ))
	#DielMat_sub2 = DielMat_sub2.add(eye(3)*2)
	#print DielMat_sub2
	
	#print '\n\nEigenvalues (2):'
	#print DielMat_sub2.eigenvals()
	#print 'Eigenvectors (2):'
	#print DielMat_sub2.eigenvects()[0][2][0] # matrix form
	
	#print '\n'
	
	
	print 'np version::'
	
	DM = np.array(DielMat_sub2.evalf())
	print 'Matrix:'
	print DM
	print 'Array:'
	DMa = np.zeros((3,3),dtype='complex')
	for i in range(3):
		for j in range(3):
			DMa[i,j] = np.complex128(DM[i,j])
	
	print DMa
	print DMa.shape
	print type(DMa)
	print type(DMa[0][0])
	
	# use numpy to find eigs
	ev = null(DMa)
	print ev/abs(ev[2])
	print 'Null (SYMPY)::', ev
	

	print '\n\n\n'

	
	
	###### all numpy version
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

def null(A, eps=1e-15):
	#http://stackoverflow.com/questions/5889142/python-numpy-scipy-finding-the-null-space-of-a-matrix
	u, s, vh = la.svd(A)
	null_mask = (s <= eps)
	null_space = scipy.compress(null_mask, vh, axis=0)
	return scipy.transpose(null_space)	
	#print ev[2][0]**2 + ev[2][1]**2 + ev[2][2]**2
	
def main():
	import spectra as sp
	p_dict = {'Bfield':700,'rb85frac':1,'Btheta':0*np.pi/2,'Bphi':0*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	chiL,chiR,chiZ = sp.calc_chi([-3500],p_dict)
	
	print 'ez: ',chiZ + 1 # ez / e0

	print 'ex: ',0.5*(2+chiL+chiR) # ex / e0
	print 'exy: ',0.5j*(chiR-chiL) # exy / e0
	
	solve_diel(chiL[0],chiR[0],chiZ[0],0)