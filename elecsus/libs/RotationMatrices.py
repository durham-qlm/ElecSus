# Copyright 2017 J. Keaveney

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
Rotation matrices for rotating cartesian coordinate systems around a particular axis

Last updated 2018-02-19 JK
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)


import numpy as np
from numpy import cos, sin

# rotation matrices:

# {{1,0,0},{0,cos(-b),-sin(-b)},{0,sin(-b),cos(-b)}}

# {{cos(a), 0, sin(a)},{0,1,0},{-sin(a),0,cos(a)}}

def rotate_forward(input_vec,phi,theta,test=False):
	"""
	With respect to the lab frame (x,y,z), the magnetic field vector is
	( 	cos(phi) sin(theta)
		sin(phi)
		cos(phi) cos(theta) )
	or in python:
	BVec = np.matrix([[cos(phi)*sin(theta)],[sin(phi)],[cos(phi)*cos(theta)]])
	
	We need to rotate this so the field points along the z' direction in the new coordinates (x',y',z')
	
	To do this we need two transforms - one around the y-axis and one around the x-axis
	
	input_vec must be either 
		[x,y,z] as list or array {input_vec.shape == (3,)} 
		or 
		[[x1,y1,z1],[x2,y2,z2],...] {input_vec.shape == (n,3)}
	"""
	BVec = np.matrix([[cos(phi)*sin(theta)],[sin(phi)],[cos(phi)*cos(theta)]])

	# input_vec is given as a row vector [Ex,Ey,Ez] - translate to column vector [[Ex],[Ey],[Ez]]
	input_vec_col = np.array([input_vec]).T # note - need the extra []'s around input_vec!
	
	# rotation around the x-axis into xz plane
	R1 = np.matrix([[1,0,0],[0,cos(phi),-sin(phi)],[0,sin(phi),cos(phi)]])
	# rotation around the y-axis oriented around z
	R2 = np.matrix([[cos(-theta), 0, sin(-theta)],[0,1,0],[-sin(-theta),0,cos(-theta)]])
	
	if test:
		return R1 * R2 * np.matrix(input_vec_col), R1 * R2 * BVec
	else:
		## apply rotation matrices, and output as row vector [Ex,Ey,Ez] - need another transpose
		return np.array(R1 * R2 * np.matrix(input_vec_col)).T
		
def rotate_back(input_vec,phi,theta):
	""" Reverse the rotation that rotate_foward() method does """
	# input_vec is given as a row vector [Ex,Ey,Ez] - translate to column vector [[Ex],[Ey],[Ez]]
	input_vec_col = np.array([input_vec]).T # note - need the extra []'s around input_vec!
	
	R1 = np.matrix([[1,0,0],[0,cos(phi),-sin(phi)],[0,sin(phi),cos(phi)]])
	R2 = np.matrix([[cos(-theta), 0, sin(-theta)],[0,1,0],[-sin(-theta),0,cos(-theta)]])
	
	# apply inverse matrices in the reverse order
	return np.array(R2.I * R1.I * np.matrix(input_vec_col))

def rotate_around_z(input_vec,phi):
	"""
	Rotate 3d vector around the 3rd dimension, counter-clockwise
	"""
	#print 'input shape: ',input_vec[0].shape
	
	# input_vec is given as a row vector [Ex,Ey,Ez] - translate to column vector [[Ex],[Ey],[Ez]]
	input_vec_col = np.array([input_vec]).T # note - need the extra []'s around input_vec!
	
	# rotation around the z-axis (rotation in the xy plane)
	R1 = np.matrix([[cos(phi),-sin(phi),0],[sin(phi),cos(phi),0],[0,0,1]])
	
	#print np.dot(R1,np.matrix(input_vec_col))
	return np.array(R1 * np.matrix(input_vec_col))
		
def test_forward_rotn():	
	""" Testing ... """

# field along Z
	theta = 0
	phi = 0
	
	print('B-Field along Z')
	
	# X-axis expressed in x',y',z' coords
	X_new, test = rotate_foward([[1],[0],[0]], phi, theta, True)
	print(('This should be (0,0,1) by definition... \n',test))
	print(('X_new: \n',X_new))
	# Y-axis expressed in x',y',z' coords
	Y_new, test = rotate_foward([[0],[1],[0]], phi, theta, True)
	print(('Y_new: \n',Y_new))
	# Z-axis expressed in x',y',z' coords
	Z_new, test = rotate_foward([[0],[0],[1]], phi, theta, True)
	print(('Z_new: \n',Z_new))

# field along Y
	theta = 0
	phi = np.pi/2
	
	print('\n\nB-Field along Y')

	# X-axis expressed in x',y',z' coords
	X_new, test = rotate_foward([[1],[0],[0]], phi, theta, True)
	print(('This should be (0,0,1) by definition... \n',test))
	print(('X_new: \n',X_new))
	# Y-axis expressed in x',y',z' coords
	Y_new, test = rotate_foward([[0],[1],[0]], phi, theta, True)
	print(('Y_new: \n',Y_new))
	# Z-axis expressed in x',y',z' coords
	Z_new, test = rotate_foward([[0],[0],[1]], phi, theta, True)
	print(('Z_new: \n',Z_new))

# field along X
	theta = np.pi/2
	phi = 0
	
	print('\n\nB-Field along X')

	# X-axis expressed in x',y',z' coords
	X_new, test = rotate_foward([[1],[0],[0]], phi, theta, True)
	print(('This should be (0,0,1) by definition... \n',test))
	print(('X_new: \n',X_new))
	# Y-axis expressed in x',y',z' coords
	Y_new, test = rotate_foward([[0],[1],[0]], phi, theta, True)
	print(('Y_new: \n',Y_new))
	# Z-axis expressed in x',y',z' coords
	Z_new, test = rotate_foward([[0],[0],[1]], phi, theta, True)
	print(('Z_new: \n',Z_new))

def test_reverse_rotn():
	""" Testing reverse rotation ... """
	
	lab_frame_E = np.matrix([[0],[1],[0]])
	phi = np.random.random()*np.pi
	theta = np.random.random()*2*np.pi
	
	print('Lab frame input: ')
	print(lab_frame_E)
	
	E_new = rotate_foward(lab_frame_E, phi, theta)
	print('Field: ')
	print(E_new)
	
	E_original = rotate_back(E_new,phi,theta)
	print('Rotated back:')
	print(E_original)
