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

Last updated 2018-07-12 MAZ
"""
# py 2.7 compatibility


import numpy as np
from numpy import cos, sin

# [ES4 MOD] Converted off the deprecated np.matrix subclass (NumPy has signalled
# its removal): plain ndarrays with the @ operator and np.linalg.inv.  np.matrix
# silently collapsed a trailing length-1 axis when handed a 3-D array -- which is
# how the stacked (n,3) input path reached these functions -- so that behaviour
# is reproduced explicitly by _to_2d() below.
def _to_2d(vec):
	"""Reproduce np.matrix()'s collapse of trailing length-1 axes."""
	arr = np.asarray(vec)
	while arr.ndim > 2 and arr.shape[-1] == 1:
		arr = arr[..., 0]
	return arr


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
	BVec = np.array([[cos(phi)*sin(theta)],[sin(phi)],[cos(phi)*cos(theta)]])
	
	We need to rotate this so the field points along the z' direction in the new coordinates (x',y',z')
	
	To do this we need two transforms - one around the y-axis and one around the x-axis
	
	input_vec must be either 
		[x,y,z] as list or array {input_vec.shape == (3,)} 
		or 
		[[x1,y1,z1],[x2,y2,z2],...] {input_vec.shape == (n,3)}
	"""
	BVec = np.array([[cos(phi)*sin(theta)],[sin(phi)],[cos(phi)*cos(theta)]])

	# input_vec is given as a row vector [Ex,Ey,Ez] - translate to column vector [[Ex],[Ey],[Ez]]
	input_vec_col = np.array([input_vec]).T # note - need the extra []'s around input_vec!
	
	# rotation around the x-axis into xz plane
	R1 = np.array([[1,0,0],[0,cos(phi),-sin(phi)],[0,sin(phi),cos(phi)]])
	# rotation around the y-axis oriented around z
	R2 = np.array([[cos(-theta), 0, sin(-theta)],[0,1,0],[-sin(-theta),0,cos(-theta)]])
	
	if test:
		return R1 @ R2 @ _to_2d(input_vec_col), R1 @ R2 @ BVec
	else:
		## apply rotation matrices, and output as row vector [Ex,Ey,Ez] - need another transpose
		return (R1 @ R2 @ _to_2d(input_vec_col)).T
		
def rotate_back(input_vec,phi,theta):
	""" Reverse the rotation that rotate_forward() method does """
	# input_vec is given as a row vector [Ex,Ey,Ez] - translate to column vector [[Ex],[Ey],[Ez]]
	input_vec_col = np.array([input_vec]).T # note - need the extra []'s around input_vec!
	
	R1 = np.array([[1,0,0],[0,cos(phi),-sin(phi)],[0,sin(phi),cos(phi)]])
	R2 = np.array([[cos(-theta), 0, sin(-theta)],[0,1,0],[-sin(-theta),0,cos(-theta)]])
	
	# apply inverse matrices in the reverse order
	return np.linalg.inv(R2) @ np.linalg.inv(R1) @ _to_2d(input_vec_col)

def rotate_around_z(input_vec,phi):
	"""
	Rotate 3d vector around the 3rd dimension, counter-clockwise
	"""
	#print 'input shape: ',input_vec[0].shape
	
	# input_vec is given as a row vector [Ex,Ey,Ez] - translate to column vector [[Ex],[Ey],[Ez]]
	input_vec_col = np.array([input_vec]).T # note - need the extra []'s around input_vec!
	
	# rotation around the z-axis (rotation in the xy plane)
	R1 = np.array([[cos(phi),-sin(phi),0],[sin(phi),cos(phi),0],[0,0,1]])
	
	return R1 @ _to_2d(input_vec_col)

def rotate_around_z2(input_vec,phi):
	"""
	Rotate 3d vector around the 3rd dimension, counter-clockwise
	"""
	#This function adds to rotate_around_z the capability to act on a stacked array of vectors.
	#print 'input shape: ',input_vec[0].shape
	
	# input_vec is given as an array n of row vectors, eg. n x [Ex,Ey,Ez] - translate to array of n column vectors n x [[Ex],[Ey],[Ez]]
	input_vec_col = np.transpose(np.array([input_vec]),(1,2,0)) # note - need the extra []'s around input_vec!
	
	# rotation around the z-axis (rotation in the xy plane)
	R1 = np.array([[cos(phi),-sin(phi),0],[sin(phi),cos(phi),0],[0,0,1]])
	
	#print np.dot(R1,np.matrix(input_vec_col))
	return np.array(R1 @ input_vec_col) #Return rotated vectors as a stacked array
		
