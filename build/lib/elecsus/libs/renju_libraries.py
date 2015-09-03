# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 17:13:07 2014

@author: R
"""

import csv
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
#print '---------renju_libraries imported----------'

def getColumn(filename, column):
    results = csv.reader(open(filename), delimiter=",")
    return [result[column] for result in results]
    
def make2ArraysFromCSV(filename):
    x_data1 = np.array(getColumn(filename,0)); x_data1 = x_data1.astype(float)
    y_data1 = np.array(getColumn(filename,1)); y_data1 = y_data1.astype(float)
    return x_data1, y_data1
    
def makeCSVFrom2Arrays(x_data, y_data, filename):
    file1 = np.transpose(np.vstack((x_data, y_data)))
    np.savetxt(filename, file1, delimiter=","); #print "Saved", filename


def findFWHM(x,y):
#    return 0
    '''Finds the full width at half-maximum of an array'''
    
    idx_max = np.argmax(y) # The index of y.max()
    #print "idx_max", idx_max
    left_edge = x[idx_max]
    right_edge = x[idx_max]
    ''' The next bit of code only works for symmetrical one peaked functions'''
#    idx_half_max = np.abs(y - y.max()/2 ).argmin()
#    FWHM = 2 * (x[idx_max] - x[idx_half_max])
#    print FWHM
    ''' For local peaks '''
    ''' moves downhill leftwards'''
    for i in range(idx_max, 0, -1):
        left_edge = x[i]
        if (y[i] <= y.max()/2) :
            left_edge = x[i]
            #print "left_edge", left_edge
            break
    ''' moves downhill rightwards'''
    for i in range(idx_max, len(x), 1):
        if (y[i] <= y.max()/2) : 
            right_edge = x[i]
            #print "right_edge", right_edge
            break
    
    FWHM = right_edge - left_edge
    return FWHM
    
def findENBW(x,y):
    distance = (x[len(x)-1] - x[0])/(len(x)-1)      #    print "Space between x_values"
    # Compute the area using the composite trapezoidal rule.
    area = scipy.integrate.trapz(y, dx=distance)    #    print "area =", area
    Transmission_max = y.max()                      # print "Max Transmission:", Transmission_max
    ENBW = area/Transmission_max  
    return ENBW
    
def findMax(x,y):
    Transmission_max = y.max()
#    print "Max Transmission:", Transmission_max
    #    idx_max = np.abs(y - y.max() ).argmin()   # The index of the value nearest y.max()
    idx_max = np.argmax(y) # The index of y.max()
    #print "Max occurs at detuning (GHz) of:", x[idx_max]
    return Transmission_max, idx_max

def makeRectangle(x, height, centre, width):
    y = np.zeros(len(x))
    left_edge = centre - width/2
    right_edge = centre + width/2
    y[left_edge:right_edge] = height
    return y
    
def lorentzian(x, amplitude, centre, width, offset):
    """ Analytic Lorentzian function with amplitude 'a', center 'c', width 'w'.
    The FWHM of this fn is 2*w 
    NOT NORMALISED """
    amplitude = float(amplitude); centre = float(centre); width = float(width); offset = float(offset)
    L = amplitude* ( width**2 / ( (x-centre)**2 + width**2 )) + offset   
    return L

def plotCSV(csvFile, x_label, y_label):
    x, y = make2ArraysFromCSV(csvFile)
    plt.figure()
    plt.xlabel(x_label); plt.ylabel(y_label)
    plt.plot(x,y)
    plt.show()

if __name__ == '__main__':
    x_data = np.arange(-10000,10000,10) # the detuning axis in MHz.
    y_data1 = makeRectangle(x_data, 1, 570, 14)
    y_data2 = lorentzian(x_data, 1.0, -4300.0, 2.0, 0)
    y_data3 = lorentzian(x_data, 1.0, -4300.0, 20.0, 0)
    import matplotlib.pyplot as plt
    plt.plot(x_data, y_data1)
    plt.plot(x_data, y_data2)
    plt.plot(x_data, y_data3)
    plt.show()
