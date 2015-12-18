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

""" Data processing tools - binning and moving average smoothing """

import numpy as np

def bin_data(x,y,blength):
	""" Takes 2 arrays x and y and bins them into groups of blength. """
	if blength % 2 == 0: 
		blength -= 1
	nobins = len(x)/blength
	xmid = (blength-1)/2
	xbinmax = nobins*blength - xmid
	a=0
	binned = np.zeros((nobins,3))
	xout,yout,yerrout = np.array([]), np.array([]), np.array([])
	for i in range(int(xmid),int(xbinmax),int(blength)):
		xmin = i-int(xmid)
		xmax = i+int(xmid)
		xout = np.append(xout,sum(x[xmin:xmax+1])/blength)
		yout = np.append(yout,sum(y[xmin:xmax+1])/blength)
		yerrout = np.append(yerrout,np.std(y[xmin:xmax+1]))
	return xout,yout,yerrout

def smooth_data(data,degree,dropVals=False):
	"""performs moving triangle smoothing with a variable degree."""
	"""note that if dropVals is False, output length will be identical
	to input length, but with copies of data at the flanking regions"""
	triangle = np.array(range(degree)+[degree]+range(degree)[::-1])+1
	smoothed = []
	for i in range(degree,len(data)-degree*2):
		point = data[i:i+len(triangle)]*triangle
		smoothed.append(sum(point)/sum(triangle))
	if dropVals: return smoothed
	smoothed = [smoothed[0]]*(degree+degree/2)+smoothed

	j = len(data)-len(smoothed)
	if j%2==1:
		for i in range(0,(j-1)/2):
			smoothed.append(data[-1-(j-1)/2+i])
			smoothed.insert(0,data[(j-1)/2-i])
		smoothed.append(data[-1])
	else:
		for i in range(0,j/2):
			smoothed.append(data[-1-i])
			smoothed.insert(0,data[i])
	#print j,len(data),len(smoothed)
	return np.array(smoothed)