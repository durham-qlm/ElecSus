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
To use:

In any other python file, or from the interpreter:

>>> import sys
>>> sys.path.append('<path to this file>')
>>> from durhamcolours import *

alternately, if ElecSus is installed:

>>> from elecsus.libs.durhamcolours import *

From any python interactive session or script,

>>> durhamcolours.help_dcols()

will print a plot with all of the colours on, named (requires matplotlib).

Last updated 2018-02-19 JK
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)

##
# Colour numbers									Colour				Code (website)
##	
d_purple=[126.0/255.0,49.0/255.0,123.0/255.0] 		# 	Palatinate Purple	255C
d_lightpurple=[216.0/255.0,172.0/255.0,244.0/255.0]	# 	Light purple		257C
d_midblue = [145./255,184.0/255.0,189.0/255.0]		# 	Mid Blue			5493C
d_olive	= [159.0/255.0,161.0/255.0,97.0/255.0] 		# 	Olive Green			5835C	
d_blue	= [0,99.0/255.0,136.0/255.0] 			# 	Blue				634C
d_red	= [170.0/255.0,43.0/255.0,74.0/255.0] 		#  	Red					201C
d_yellow= [232.0/255.0,227.0/255.0,145.0/255.0] 	#  	Yellow				459C
d_lightblue = [196./255,229.0/255.0,250.0/255.0]	#	Light Blue			290C
d_pink	= [196.0/255.0,59.0/255.0,142.0/255.0] 		#  	Pink				675C
d_grey = [150./255,142./255,133./255]			#	Warm Grey			8C
d_black = [35./255,31./255,32./255]			#	Black				BlackC
d_lightgrey = [207./255,218./255,209./255]		#	Near White/L. Grey	5655C
d_midgrey = [110./255, 100./255, 100./255]

cols = [d_purple,d_blue,d_lightpurple,d_midblue,d_olive,d_red,d_yellow,d_lightblue,d_grey,d_black,d_lightgrey,d_pink]

colname = ['d_purple', 'd_blue', 'd_lightpurple', 'd_midblue', 'd_olive', 'd_red', 'd_yellow', 'd_lightblue', 'd_grey', 'd_black', 'd_lightgrey', 'd_pink']

def help_dcols():
	""" Show a plot demonstrating all the colours """
	from pylab import figure,clf,subplots,show
	from matplotlib.patches import Rectangle
	fig, axes = subplots(4,3)
	i=0
	print(axes)
	for axA in axes:
		for ax in axA:
			patch = Rectangle([0,0],1,1,color=cols[i])
			ax.add_patch(patch)
			ax.set_xticks([])
			ax.set_yticks([])
			ax.set_xticklabels([])
			ax.set_yticklabels([])			
			if colname[i] is not 'd_black':
				ax.text(0.2,0.4,colname[i],color='k')
			else:
				ax.text(0.2,0.4,colname[i],color='w')
			i+=1
	fig.subplots_adjust(bottom=0.25,left=0.03,right=0.98,top=0.98,wspace=0.05,hspace=0.05)
	fig.text(0.08,0.18,'To use: from durhamcolours import *')
	fig.text(0.08,0.12,'then e.g. ax.plot(x,y,colour=d_olive)')
	fig.text(0.08,0.06,"Note: no parentheses i.e. d_olive not 'd_olive'")
	
	fig.savefig('durhamcolors.png')
	fig.savefig('durhamcolors.pdf')
	
	show()
	
if __name__ == '__main__':
	help_dcols()
