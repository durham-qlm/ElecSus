"""
Polarisation animation...

the animate_vectors() method creates an interactive 3D plot, visualising the resultant polarisation for a
given input of Ex, Ey and the phase difference (in radians) between them
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

#replace default matplotlib text and color sequence with durham colours
plt.rc('font',**{'family':'Serif','serif':['Times New Roman']})
params={'axes.labelsize':13,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize': 11,'mathtext.fontset':'cm','mathtext.rm':'serif'}
plt.rcParams.update(params)

def update_lines(num, Ex, Ey, z, time, curve):
	# NOTE: there is no .set_data() for 3 dim data...
	curve.set_data([Ex*np.exp(-1.j*time[num]*2*np.pi),z])
	curve.set_3d_properties(Ey*np.exp(-1.j*time[num]*2*np.pi))
	return curve

def animate_vectors(Exi,Eyi,phase):
	# Attaching 3D axis to the figure
	
	E = [Exi, Eyi*np.exp(1.j*phase)]
	
	fig = plt.figure("Polarisation animation")
	ax = p3.Axes3D(fig)

	k = 2*np.pi / 2

	# 100 resultant vectors at various z
	z_axis_curve = np.linspace(-2.5,2.5,500)

	# t = 0 data, all z
	Ex = E[0] * np.exp(1.j*(k*z_axis_curve))
	Ey = E[1] * np.exp(1.j*(k*z_axis_curve))

	nframes = 100
	time = np.linspace(0,4,nframes)


	# Creating fifty line objects.
	# NOTE: Can't pass empty arrays into 3d version of plot()
	curve = ax.plot(Ex, z_axis_curve, Ey, color='k', lw=2)[0]

	spokes = 6
	lines = [ax.plot([0,Ex[::spokes][i]],[z_axis_curve[::spokes][i],z_axis_curve[::spokes][i]],[0,Ey[::spokes][i]],color='k',alpha=0.3, lw=1)[0] for i in range(len(Ex[::spokes]))]

	x_quiver = ax.quiver3D([0],[0],[0],[1],[0],[0],length=2,arrow_length_ratio=0.05,pivot='middle',color='k',lw=2)
	y_quiver = ax.quiver3D([0],[0],[0],[0],[0],[1],length=2,arrow_length_ratio=0.05,pivot='middle',color='k',lw=2)
	#z_quiver = ax.quiver3D([0],[0],[0],[0],[1],[0],length=5,arrow_length_ratio=0.05,pivot='middle',color='k',lw=2)

	k_quiver = ax.quiver3D([0],[0],[0],[0],[1],[0],length=5,arrow_length_ratio=0.05,pivot='middle',color='r',lw=3,alpha=0.6)


	ax.text(0, 2.2, 0.15, r"$\vec{k}, \vec{z}$", (0,1,0), color='red', size=18)
	ax.text(0.9, 0, 0.1, r"$\vec{x}$", (1,0,0), color='k', size=18)
	ax.text(0.1, 0, 0.9, r"$\vec{y}$", (1,0,0), color='k', size=18)

	ax.set_xlim3d(-1,1)
	ax.set_zlim3d(-1,1)

	# Creating the Animation object
	line_ani = animation.FuncAnimation(fig, update_lines, nframes, fargs=(Ex, Ey, z_axis_curve,time, curve),
								   interval=50, blit=False)
	
	plt.show()
	
if __name__ == '__main__':
	animate_vectors(1,1.j,0)