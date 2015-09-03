# -*- coding: utf-8 -*-

import spectra
import numpy as np
import matplotlib.pyplot as plt
import time; start_time = time.time()
import renju_libraries
import csv
import time
import math
import warnings
from matplotlib import cm

warnings.simplefilter("ignore")

plt.ion()

def displayContour(param_data,T_points,T_step,B_points,B_step):
	Temp = []
	B_field = []
	Trans = []
	Uniq = []
	for line in param_data:
		Temp.append(line[1])
		B_field.append(line[0])
		Trans.append(line[4])
		Uniq.append(line[6])

	Temp = np.array(Temp)
	B_field = np.array(B_field)
	Trans = np.array(Trans)
	Uniq = np.array(Uniq)

	xy_Tran = np.zeros((T_points,B_points))
	xy_Uniq = np.ones((T_points,B_points))

	i = 0
	
	for el in Trans:
		j = int((Temp[i]-Temp.min())/T_step)
		k = int((B_field[i]-B_field.min())/B_step)
		xy_Tran[j][k] = el
		i += 1
	
	i = 0
	
	for el in Uniq:
		j = int((Temp[i]-Temp.min())/T_step)
		k = int((B_field[i]-B_field.min())/B_step)
		xy_Uniq[j][k] = el
		i += 1

	x = 25
	y = 1
	end = 2

	xy_combo = np.multiply(xy_Uniq,xy_Tran)
	
	fig = plt.figure()

	levels = [Uniq.max()*0.7,Uniq.max()*0.8,Uniq.max()*0.9,2.0]

	ax1 = plt.subplot2grid((y,x),(0,0),colspan=x-end)
	ax2 = plt.subplot2grid((y,x),(0,x-end),colspan=end-1)
	extent = (Temp.min(),Temp.max(),B_field.min(),B_field.max())

	heatmap = ax1.imshow(xy_Uniq.transpose(),cmap=plt.get_cmap('Greys'),extent=extent,aspect='auto',origin='lower')
	imCON = ax1.contour(xy_Uniq.transpose(),levels=[2],extent=extent,linewidths = 2,colors='r',linestyles=('dashed'))
	cb = fig.colorbar(heatmap,cax=ax2)
	

def findNextPeak(X, Y_elec):

	peak1_x = np.argmax(Y_elec)
	peak1_y = Y_elec.max()
	peaks_detected = 1
	h = 0.99

	if peak1_y < 0.9:
		return peak1_y

	while(peaks_detected == 1):
		old_sign = -1
		changes = 0.0
		y_scan_height = peak1_y*h
		for i in Y_elec:
			sign = int((i-y_scan_height)/abs(i-y_scan_height))
			if sign != old_sign:
				changes += 1
				old_sign = sign
		peaks_detected = int(changes/2.0)
		h -= 0.01
		if h < 0.1:
			break
	return y_scan_height
		


def write_csv(xy,filename):
	""" newer, better module for writing csv data with arbitrary
		number of columns to filename.
		takes in xy, which should be of the form [[x1,y1],[x2,y2] ...]
		this can be done by zipping arrays, e.g.
			xy = zip(x,y,z)
			where x,y and z are 1d arrays
	"""	
	
	with open(filename, 'wb') as csvfile:
		csv_writer = csv.writer(csvfile,delimiter=',')
		for xy_line in xy:
			csv_writer.writerow(xy_line)

if __name__ == '__main__':
	"""
	Created on Fri Jul 18 14:23:58 2014

	@author: R
	"""

	''' Testing RRFittingRoutine '''


	X = np.arange(-20,20,0.02) # the detuning axis in GHz.
	Elem = 'Rb' # the chosen alkali element.
	StokesType = 'Ix'   # specifies the type of spectrum required
	Bfield = 230.0# the magnitude of the magnetic field in Gauss
	T = 90.0   # temperature (Celsius) linked to the number density
	lcell = 5.0  # length of the vapour cell in millimetres
	rb85frac = 72.17   # % of rubidium-85 atoms
	DoppTemp = -250.0 # temperature (Celsius) linked to the Doppler width
	theta0 = 90.0 # Linear polarisation angle (in degrees) w.r.t to the x-axis
	Pol = 50.0  # percentage of probe beam that is polarised to drive sigma minus
	shift = 0.0   # a global frequency shift in MHz
	GammaBuf = 0.0# Extra lorentzian broadening in MHz
	Constrain = True# if True, overides the DoppTemp value and sets it to T
	Dline = 'D2'# specifies which D-line transition to calculate for
	precision = 10.0   # the required precision of the calculation
	K40frac = 6.7302# number of potassium-40 atoms divided by the total atom number
	K41frac = 0.0117# number of potassium-40 atoms divided by the total atom number

	all_params = []

	T_begin=25
	T_end=170
	T_step=1
	T_points = math.ceil((T_end - T_begin)/float(T_step))

	B_begin=2500
	B_end=5000
	B_step=10
	B_points = math.ceil((B_end - B_begin)/float(B_step))
	
	kk = 0

	Total_points = T_points*B_points

	print 'Total iterations:  ' + str(int(Total_points))
	if raw_input('Proceed? (y/n)') =='y':
#		print 'A'

		end_times = []
		iter_times = []
		start_time = time.time()
		end_times.append(start_time)
#		print 'B'
	
		for t in np.arange(T_begin,T_end,T_step):
#			print 'C'
			#print 'Temperature:  ' + str(t)
			for B in np.arange(B_begin,B_end,B_step):
				#print 'B field:  ' + str(B)
#				print 'D'
				x_data1 =  X 
				Bfield = B
				T = t
#				print 'E'
				bestFitParams = [Elem,StokesType,Bfield,T,lcell,
						rb85frac,DoppTemp,theta0,Pol,shift,
						GammaBuf,Constrain,Dline,precision,
						K40frac,K41frac]
#				print 'F'				
				paramBools = [True,True,False,False,False,False,False,False,False]
#				print 'G'
				x_indices = np.arange(0,len(x_data1),1)
#				print 'H'

				Yfit = spectra.spectrum(X,Elem,StokesType, Bfield,T,lcell,
							rb85frac,DoppTemp,theta0,Pol,shift,
							GammaBuf,Constrain,Dline,precision,
							K40frac,K41frac)
#				print 'I'				
#				plt.plot(X ,Yfit , label = str(Bfield)+'G'+str(T)+'C')
#				plt.legend(loc = 'best')
#				plt.show()
				
				FWHM = renju_libraries.findFWHM(X, Yfit)
#				print 'J'
				ENBW = renju_libraries.findENBW(X, Yfit)
#				print 'K'
				Transmission_max, idx_max = renju_libraries.findMax(X,Yfit)
				next_best = findNextPeak(X, Yfit)
				uniqueness = Transmission_max/float(next_best)
				if uniqueness<1.0:
					uniqueness = 1					
#				print "----------------------------------------------"
#				print "Bfield (Gauss):", Bfield
#				print "Temp (degrees C)", T
#				print "ENBW", ENBW
#				print "Max Transmission", Transmission_max
#				print "Max occurs at detuning (GHz) of:", X[idx_max]/1000.0
#				print "FWHM (MHz): ", FWHM/1000.0
				all_params.append([Bfield,T,(ENBW/1000.0),(FWHM/1000.0),Transmission_max,X[idx_max],uniqueness])
				#print "ENBW (Arb units): ", ENBW
#				print all_params
				
				end_times.append(time.time())
				iter_times.append(end_times[kk+1]-end_times[kk])
				kk += 1
				avg_time = np.mean(iter_times)
				time_display = int(avg_time*(Total_points-kk))
				
				mins = int(time_display/60.0)
				hrs = int(time_display/3600.0)
				mins = mins - 60*hrs
				secs = time_display - 60*mins - 3600*hrs
				
				print str(kk) + '/' + str(int(Total_points)),
				print '    Approx.  ',
				if hrs != 0:
					print str(hrs)+'hrs ',
				if mins != 0:
					print str(mins)+'m ',
				if hrs == 0:
					print str(secs) + 's ',
				print'   remaining'
	write_csv(all_params,Elem+StokesType+'dataT'+str(T_begin)+'to'+str(T_end)+'B'+str(B_begin)+'to'+str(B_end)+'.csv')
	if raw_input("Display contour? (y/n)") =='y':	
		displayContour(all_params,T_points,T_step,B_points,B_step)
	raw_input("Press any key to end")
	 
