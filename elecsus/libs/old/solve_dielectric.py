from sympy import Symbol, solve, cos, sin, pi, simplify, eye, powsimp, powdenest
from sympy.matrices import det, Matrix

import numpy as np
import scipy.linalg as la

from FundamentalConstants import e0
''' Use symbolic python to solve for the roots of n-squared in the dielectric tensor matrix '''

def solve_diel(chiL, chiR, chiZ, theta):
	"""
	Solve dielectric tensor for arbitrary angle theta.
	
	Get eigenvectors of the dielectric tensor matrix, and use these to form a rotataion
	matrix.
	Find the two complex indices that solve the equation det(DielTensor) = 0
	
	Returns rotation matrix and two associated complex refractive indices
	
	**
	NOTE::
	Having issues with this so currently only implemented for theta=0 and theta=pi/2
	**
	"""
	
	'''
	# Sympy symbols
	theta = Symbol('theta',real=True)
	n_sq = Symbol('n_sq')
	e_x = Symbol('e_x')
	e_xy = Symbol('e_xy')
	e_z = Symbol('e_z')

	# Dielectric tensor in the rotated frame
	DielTensor = Matrix (( 	[(e_x - n_sq)*cos(theta), e_xy, e_x*sin(theta)],
								[-e_xy * cos(theta), e_x - n_sq, -e_xy*sin(theta)],
								[(n_sq - e_z)*sin(theta), 0, e_z*cos(theta)] 			))								

	# make the substitution for the angle theta
	DielMat_sub = DielMat.subs(theta, theta2)

	## also substitute in for ex, exy and ez
	#DielMat_sub = DielMat_sub.subs(e_x, 0.5*(2.+chiL+chiR))
	#DielMat_sub = DielMat_sub.subs(e_xy, 0.5j*(chiR-chiL))
	#DielMat_sub = DielMat_sub.subs(e_z, 1.0+chiZ)
								
	# Get solutions for the complex refractive indices
	solns = solve(det(DielMat_sub), n_sq)
	print 'Solutions for n_sq:'
	print solns[0]-1,'\n'
	print solns[1]-1,'\n'

	## Analytic solutions for theta = pi/2
	solns[0] = e_x + e_xy**2/e_x
	solns[1] = e_z

	## Analytic solutions for theta = 0
	#solns[0] = e_x + 1.j*e_xy
	#solns[1] = e_x - 1.j*e_xy

	# To find eigenvectors - plug solutions for n_sq back in and diagonalise
	# Get the eigenvector corresponding to the zero-energy eigenvalue
	DielMat_sub1 = DielMat_sub.subs(n_sq, solns[0])
	DielMat_sub1 = DielMat_sub1.subs(e_x, e0*0.5*(2.+chiL+chiR))
	DielMat_sub1 = DielMat_sub1.subs(e_xy, e0*0.5j*(chiR-chiL))
	DielMat_sub1 = DielMat_sub1.subs(e_z, e0*(1.0+chiZ))	
	Evect1 = DielMat_sub1.eigenvects()[0][2][0]

	DielMat_sub2 = DielMat_sub.subs(n_sq, solns[1])
	DielMat_sub2 = DielMat_sub2.subs(e_x, e0*0.5*(2.+chiL+chiR))
	DielMat_sub2 = DielMat_sub2.subs(e_xy, e0*0.5j*(chiR-chiL))
	DielMat_sub2 = DielMat_sub2.subs(e_z, e0*(1.0+chiZ))
	Evect2 = DielMat_sub2.eigenvects()[0][2][0] # matrix form

	evects = [Evect1, Evect2]



	###### all numpy version
	print 'Numpy version...........'
	NSQ = EX + EXY**2/EX
	NSQ = EZ
	#NSQ = EX + 1.j*EXY
	#NSQ = EX - 1.j*EXY
	DielMat = np.matrix ([ 	[(EX - NSQ)*np.cos(THETA), EXY, EX*np.sin(THETA)],
								[-EXY*np.cos(THETA), EX - NSQ, -EXY*np.sin(THETA)],
								[(NSQ - EZ)*np.sin(THETA), 0, EZ*np.cos(THETA)] 			])

	print DielMat
	ee, ev = la.eig(DielMat)
	print 'EigenVals:'
	print ee[0], ee[1], ee[2]
	print 'EigenVects:'
	print ev[0]
	print ev[1]
	print ev[2]
	'''
	
	solns = [0,0]
	EX = 0.5*(2.+chiL+chiR)
	EXY = 0.5j*(chiR-chiL)
	EZ = 1.0+chiZ

	if theta == np.pi/2:
		solns[0] = EX + EXY**2/EX
		solns[1] = EZ
		
		# Eigenvectors
		evects = [ [np.zeros(len(EX)), EX/EXY, np.ones(len(EX))], [np.ones(len(EX)), np.zeros(len(EX)), np.zeros(len(EX))], [np.zeros(len(EX)), np.zeros(len(EX)), np.ones(len(EX)) ]]
	elif theta == 0.0:
		solns[0] = EX - 1.j*EXY
		solns[1] = EX + 1.j*EXY
		
		# Eigenvectors
		evects = [ [1.j*np.ones(len(EX)), 1.0*np.ones(len(EX)), 0.0*np.ones(len(EX))], [-1.j*np.ones(len(EX)), 1.0*np.ones(len(EX)), 0.0*np.ones(len(EX))], [np.zeros(len(EX)), np.zeros(len(EX)), np.ones(len(EX)) ] ]
	else:
		raise ValueError(" theta is currently only supported for 0 or pi/2 ")
	
	RotMat = np.array(evects)
	#print 'RM:::', RotMat
	
	# diagnostics for general version
	print 'Rotation Matrices shape: ',RotMat.shape
	print 'Refractive indices shape:',solns[0].shape
	return RotMat, np.sqrt(solns[0]), np.sqrt(solns[1])
	
	

	
def main():
	import spectra as sp
	p_dict = {'Bfield':700,'rb85frac':1,'Btheta':np.pi/2,'Bphi':0*np.pi/180,'lcell':75e-3,'T':84,'Dline':'D2','Elem':'Cs'}
	chiL,chiR,chiZ = sp.calc_chi(np.arange(-3500,3500,100),p_dict)
	
	print 'ez: ',chiZ + 1 # ez / e0

	print 'ex: ',0.5*(2+chiL+chiR) # ex / e0
	print 'exy: ',0.5j*(chiR-chiL) # exy / e0
	
	solve_diel(chiL,chiR,chiZ,np.pi/2)