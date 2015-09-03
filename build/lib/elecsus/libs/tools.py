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

"""Module to provide functions for data manipulation and plotting"""

from numpy import zeros, average, amax, amin, array

def derivative(xs,f):
    '''Numerical derivative. Last point is simply a copy of previous'''
    flen = len(f)
    df=zeros(flen)
    for y in xrange(flen-1):
        dx = xs[y+1]-xs[y]
        df[y] = (f[y+1]-f[y])/dx
    df[-1] = df[-2] # Last point is set to the value of the one before
    return df

def smoother(xData,yData,factor):
    '''Function to smooth over the data'''

    oldlength = len(xData)
    print 'Original data length:', oldlength
    truncation = oldlength % factor
    if truncation:
        print 'Data truncated at the end by', truncation, 'points.'
        xData = xData[0:-truncation]
        yData = yData[0:-truncation]

    newlength = (oldlength-truncation)/factor
    print 'Smoothed data length:', newlength
    smoothedX = zeros(newlength)
    smoothedY = zeros(newlength)

    i = 0
    while i < newlength:
        smoothedX[i] = average(xData[(i*factor):((i+1)*factor)])
        smoothedY[i] = average(yData[(i*factor):((i+1)*factor)])
        i+=1

    return smoothedX,smoothedY

def fileOutput(filename,x,y):
    """Outputs two arrays or lists as columns in a csv file"""
    f_csv = open(filename, 'w')
    for i in xrange(len(x)):
        print >> f_csv, x[i], ',', y[i]
    f_csv.close()
    return None

def read_in_twoColumn(filename):
    """Takes data from a two column file and returns two arrays."""
    try:
        data_file = open(filename + '.csv','r')
    except:
        print 'Data file not found.'
        exit(1)
    xData = []
    yData = []
    for line in data_file:
        newline = line[0:-1].split(',')
        xData.append(float(newline[0]))
        yData.append(float(newline[1]))
    data_file.close()
    xData = array(xData)
    yData = array(yData)
    if amax(abs(xData)) > 500.0:
        print "The data file detuning values are very large (The detuning axis should be given in GHz)."
        convert = raw_input("Do you want to convert the detuning axis from MHz to GHz? (y/n): ")
        if convert in ['n','N','no','No','NO','nO']:
            print "This may use large amounts of memory and calculation time"
            contin = raw_input('Continue? (y/n): ')
            if contin in ['n','N','no','No','NO']:
                exit(0)
        elif not (convert in ['n','N','no','No','NO','nO']):
            xData = xData/1000.0
    return xData, yData

def plotOutput(X,Y,spectrumLabel,ydata,PBool=True,ydataBool=False,residuals=False,path=False):
    """Plots the program outputs for the user"""
    import matplotlib
    if PBool == False:
        matplotlib.use('Agg') #non-interactive backend
    import pylab as P
    P.ioff() #Ensure interactivity mode is off so that graph does not dissapear immediately
    fig = P.figure()
    maxYval = amax(Y)
    minYval = amin(Y)
    DynamicRange = maxYval - minYval
    if not ydataBool:
        P.plot(X,Y,'g', linewidth = 2.0)
        P.xlabel(r'Detuning (GHz)')
        P.ylabel(spectrumLabel)
        P.xlim(X[0],X[-1])
        P.ylim(minYval-0.02*DynamicRange,maxYval+0.02*DynamicRange)
    else:
        ax1 = fig.add_axes([0.15,0.30,0.75,0.65])
        ax1.plot(X,ydata,'k')
        ax1.plot(X,Y,'r--', linewidth=1.8)
        ax1.set_xlim(X[0],X[-1])
        ax1.set_ylim(minYval-0.02*DynamicRange,maxYval+0.02*DynamicRange)
        ax1.set_xticklabels([])
        P.ylabel(spectrumLabel)
        ax2 = fig.add_axes([0.15,0.10,0.75,0.15])
        ax2.plot(X,residuals*100.0,'k')
        ax2.set_xlim(X[0],X[-1])
        ax2.axhline(color='r', linestyle = '--', linewidth=1.8)
        P.xlabel(r'Detuning (GHz)')
        P.ylabel(r'Residuals $(\times 100)$')
    if path:
        P.savefig(path)
    if PBool:
        P.show()


