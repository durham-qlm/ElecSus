"""

ElecSus GUI

A GUI based on wxpython for ElecSus, intended to replace/augment the 
runcard method of calling ElecSus.

More blurb here...

Requirements:
	python, wxpython, matplotlib, numpy, scipy
	versions unknown, tested on...

LICENSE info

James Keaveney and co-authors
2011-2015
"""

#!/usr/bin/env python
import matplotlib
matplotlib.use('WxAgg')
import pylab
pylab.ioff()

from matplotlib import rc
rc('text', usetex=False)
rc('font',**{'family':'serif'})

import wx
import os
import sys
import csv

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg as Toolbar

###IMPORT ELECSUS MODULES
##
##
##


# We use sizer constructs for placement of all GUI elements. 
# This makes the whole window scalable in a predictable way.

# Set default parameters

# Button size
BtnSize = 30

# Sizer indent defaults
Lindent = 20
Rindent = 20
indent = 20

## IMPORTANT! Master list of all output types that will be referenced for dynamic plotting
OutputTypes = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix', 'Iy']

class FittingPanel(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self,parent)
		
		t = wx.StaticText(self,-1,"This is some text", (20,20))

class TheoryPanel(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self,parent)
		
		t = wx.StaticText(self,-1,"This is some more text", (20,20))
		
class ElecSus_GUI_Frame(wx.Frame):
	""" Main window """
	def __init__(self,parent,title):
		wx.Frame.__init__(self,None,title=title,size=(1200,700))
		
		#if the window is closed, exit
		self.Bind(wx.EVT_CLOSE,self.OnExit)

		self.panel = wx.Panel(self)
		
		self.panel.SetBackgroundColour(wx.Colour(240,240,240))
		
		self._init_default_values()
	
		self._init_panel()
		
		self._init_menus()
		
		self.OnCreateAxes(1)
		
		
	def _init_default_values(self):
		""" Initialise default values for various things ... """
		
		self.PlotHold = False
		
		self.plot_outputs = ['Transmission (S0)', 'Ix']
		self.plot_output_indices = [OutputTypes.index(po) for po in self.plot_outputs]
		
		self.xrange = [-10,10] # detuning range, in GHz
		self.npoints = 10000 # number of detuning points to calculate
		
		# default data for plots
		self.x_array = np.linspace(self.xrange[0],self.xrange[1],self.npoints)
		self.y_arrays = [None, None]
		self.x_expt_arrays = [None, None]
		self.y_expt_arrays = [None, None]
		
		# ...
	
	def _init_menus(self):
		""" Initialise menu bar items """
		
		# Create menuBar object
		menuBar = wx.MenuBar()
		
		# File
		fileMenu = wx.Menu()
		fM_open = fileMenu.Append(wx.ID_OPEN, "&Open\tCtrl+O", "Open file for plotting and/or fitting.")
		self.Bind(wx.EVT_MENU, self.OnFileOpen, fM_open)
		fM_about = fileMenu.Append(wx.ID_ABOUT, "&About", "About this program.")
		self.Bind(wx.EVT_MENU, self.OnAbout, fM_about)
		fM_exit = fileMenu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", "Close window and exit program.")
		self.Bind(wx.EVT_MENU, self.OnExit, fM_exit)
		#
		menuBar.Append(fileMenu, "&File")
		
		# Plot
		plotMenu = wx.Menu()
		pM_plotholdon = plotMenu.Append(wx.ID_ANY, "&Hold data on plot update", "Select whether to hold or clear current plot data on updating the figure", kind=wx.ITEM_CHECK)
		pM_plotholdon.Check(self.PlotHold)
		self.Bind(wx.EVT_MENU, self.OnPlotHold, pM_plotholdon)
		pM_clearplot = plotMenu.Append(wx.ID_ANY, "&Clear current plot data", "Clear current plot data on all axes")
		self.Bind(wx.EVT_MENU, self.OnClearPlot, pM_clearplot)
		#
		menuBar.Append(plotMenu, "&Plot")
		
		# Fit
		fitMenu = wx.Menu()
		fitM_dummy = fitMenu.Append(wx.ID_ANY, "&Dummy", "Dummy")
		fittypeSubMenu = wx.Menu()
		fittypeSubMenu.Append(wx.ID_ANY, "&Marquardt-Levenburg", "Use ML Fitting", kind=wx.ITEM_CHECK)
		fittypeSubMenu.Append(wx.ID_ANY, "&Random Restart", "Use RR Fitting", kind=wx.ITEM_CHECK)
		fittypeSubMenu.Append(wx.ID_ANY, "&Simulated Annealing", "Use SA Fitting", kind=wx.ITEM_CHECK)
		
		### BIND EVENTS TO THESE !!!!
		
		fitM_fittype = fitMenu.AppendMenu(wx.ID_ANY, "Fit &Type", fittypeSubMenu)
		
		#
		menuBar.Append(fitMenu, "F&it")
		
			
		#
		# ... other menus
		#
		
		
		# Add Menu Bar to the main panel
		self.SetMenuBar(menuBar)
		
	def _init_panel(self):
		""" Initialise panel with matplotlib window, buttons, text boxes etc. Doesn't really need to be in a separate function... """
		
	## Create plot part of the window
		self.fig = plt.figure(1,facecolor=(240./255,240./255,240./255))
		
		# display some text in the middle of the window to begin with
		#self.fig.text(0.5,0.5,"Welcome to the ElecSus GUI\n\nSelect some options to begin",ha='center',va='center')
		
		# create the wx objects to hold the figure
		self.canvas = FigureCanvasWxAgg(self.panel, wx.ID_ANY, self.fig)
		self.toolbar = Toolbar(self.canvas) #matplotlib toolbar (pan, zoom, save etc)
		
		# Create vertical sizer to hold figure and toolbar - dynamically expand with window size
		plot_sizer = wx.BoxSizer(wx.VERTICAL)
		plot_sizer.Add(self.canvas, 1, wx.LEFT|wx.RIGHT|wx.GROW,border=0)
		plot_sizer.Add(self.toolbar, 0, wx.LEFT|wx.RIGHT|wx.EXPAND,border=0)

	## Create button part of the window
		
		# elecsus and JQC logos at the top
		elecsuslogo = wx.Image('images/elecsus-logo.png',wx.BITMAP_TYPE_ANY)
		elecsuslogo.Rescale(108/2,138/2)
		#elecsus_bmp = wx.StaticBitmap(self.panel,wx.ID_ANY,wx.BitmapFromImage(elecsuslogo),size=(108/2,-1))
		jqclogo = wx.Image('images/jqc-logo.png',wx.BITMAP_TYPE_ANY)
		#jqclogo.Rescale(191*3/4,84*3/4)
		jqc_bmp = wx.StaticBitmap(self.panel,wx.ID_ANY,wx.BitmapFromImage(jqclogo),size=(191,-1))
		
		image_sizer = wx.BoxSizer(wx.HORIZONTAL)
		image_sizer.Add((10,-1),1,wx.EXPAND)
		image_sizer.Add((40,-1),0,wx.EXPAND)
		#image_sizer.Add(elecsus_bmp,0,wx.EXPAND)
		#image_sizer.Add((10,-1),0,wx.EXPAND)
		image_sizer.Add(jqc_bmp,0,wx.EXPAND)
		image_sizer.Add((40,-1),0,wx.EXPAND)
		image_sizer.Add((10,-1),1,wx.EXPAND)
		
	## Menus are tabbed into theory and fitting sections
	## use wx.Notebook for this
		TabPanel = wx.Notebook(self.panel)
		
		# add tabs
		ThyPanel = TheoryPanel(TabPanel)
		FitPanel = FittingPanel(TabPanel)
		
		TabPanel.AddPage(ThyPanel, "Theory Settings")
		TabPanel.AddPage(FitPanel, "Fit Settings")
		
		tab_sizer = wx.BoxSizer(wx.HORIZONTAL)
		#tab_sizer.Add((40,-1),0,wx.EXPAND)
		tab_sizer.Add(TabPanel, 1, wx.EXPAND|wx.LEFT|wx.RIGHT, border=10)
		#tab_sizer.Add((40,-1),0,wx.EXPAND)
		
		
	## Create sizer for right-hand side of panel
		button_sizer = wx.BoxSizer(wx.VERTICAL)
		button_sizer.Add((-1,20),0,wx.EXPAND)
		button_sizer.Add(image_sizer,0,wx.EXPAND)
		button_sizer.Add((-1,20),0,wx.EXPAND)
		button_sizer.Add(tab_sizer,1,wx.EXPAND)
		button_sizer.Add((-1,10),0,wx.EXPAND)
		
		
	## Put plot part and button part together
		main_sizer = wx.BoxSizer(wx.HORIZONTAL)
		main_sizer.Add(plot_sizer,1,wx.EXPAND)
		main_sizer.Add(wx.StaticLine(self.panel,-1,size=(1,-1),style=wx.LI_VERTICAL),0,wx.EXPAND)
		main_sizer.Add(button_sizer,0,wx.EXPAND)

		self.panel.SetSizer(main_sizer)
		self.panel.Layout()
				
		
## Actions for events
	## General Actions
	def OnExit(self,event):
		self.Destroy()
		app.ExitMainLoop()
		
	def OnAbout(self,event):
		## do stuff
		dlg = wx.MessageDialog(self, "Blurb about this program", "About", wx.OK)
		if dlg.ShowModal() == wx.ID_OK:
			dlg.Destroy()
			
	def OnFileOpen(self,event):
		self.dirname= ''
		dlg_choice = wx.SingleChoiceDialog(self,"Choose type of data to be imported","Data import",choices=OutputTypes)
		
		# wait for OK to be clicked
		if dlg_choice.ShowModal() == wx.ID_OK:
			choice = dlg_choice.GetSelection()
			# use the choice index to select which axes the data appears on
			if choice not in self.plot_output_indices:
				self.AddPlotType(choice)
				
			choice_index = self.plot_output_indices.index(choice)
		dlg_choice.Destroy()
		
		dlg_open = wx.FileDialog(self,"Choose 2-column csv file (Detuning, Transmission)",
								self.dirname,"","*.csv",wx.OPEN)
		
		# if OK button clicked, open and read file
		if dlg_open.ShowModal() == wx.ID_OK:
			self.filename = dlg_open.GetFilename()
			self.dirname = dlg_open.GetDirectory()
			#call read
			self.x_expt_arrays[choice_index],self.y_expt_arrays[choice_index] = read_CSV(os.path.join(self.dirname,self.filename),spacing=0)
			
			self.OnUpdatePlot(1)
			
		dlg_open.Destroy()	
		
	
	## Plot-specific actions
	def OnCreateAxes(self,event):
		""" Create as many sets of axes in the current figure as are needed, and label them all """
		n_axes = len(self.plot_outputs)
		
		# create bare axes
		self.fig.add_subplot(n_axes,1,1)
		i=0
		for i in range(1,len(self.plot_outputs)):
			self.fig.add_subplot(n_axes,1,i+1,sharex=self.fig.axes[0])
		
		# add info
		ylabels = self.plot_outputs
		
		for ax,ydata,ylabel in zip(self.fig.axes,self.y_arrays,ylabels):
			if ydata is not None:
				ax.plot(self.x_array,ydata)
			ax.set_ylabel(ylabel)
		
		self.fig.axes[-1].set_xlabel('Detuning (GHz)')
		self.fig.axes[-1].set_xlim(self.xrange)
		
		self.draw_fig()
	
	def OnCreateColorMapPlot(self,event):
		""" Create the plot for colormap/contour data ... """
		pass
	
	def OnAddAxes(self,event):
		""" adjust subplot layout when axes are added/removed """
		
		# use matplotlib ax.change_geometry(<new subplot params>)
		n_axes_now = len(self.fig.axes)
		n_axes_needed = len(self.plot_outputs)
		
		# adjust old subplots
		for i, ax in enumerate(self.fig.axes):
			ax.change_geometry(n_axes_needed,1,i+1)
		
		# add new subplots
		for j in range(n_axes_now,n_axes_needed):
			self.fig.add_subplot(n_axes_needed,1,j+1,sharex=self.fig.axes[0])
		
		#remove all xaxis labels, and redraw on bottom plot. add new ylabels
		for ax,ylabel in zip(self.fig.axes,self.plot_outputs):
			ax.set_xlabel('')
			ax.set_ylabel(ylabel)
			
		self.fig.axes[-1].set_xlabel('Detuning (GHz)')
		
		self.draw_fig()
		
	def OnCreateResultsPlot(self,event):
		""" Create a new plot for fit results - external to the main window? """
		pass
		
	def OnPlotHold(self,event):
		""" Toggle plot hold (keep data on updating figure) on/off """
		self.PlotHold = bool(event.Checked())
		self.fig.hold(self.PlotHold)
			
	def OnClearPlot(self,event):
		""" Get list of all axes in figure and clear them all """
		for ax in self.fig.axes:
			ax.cla()
	
	def OnPlotLegend(self,event):
		""" Toggle plot legend on/off """
		self.legendOn = bool(event.Checked())
		
	def OnUpdatePlot(self,event):
		""" Update plot with latest xy data """
		
		# For now, simple x-y plots only - will add colormaps/contour plots later
		
		# n-sets of x,y data possible (theory and experiment)
		for ax, y, xe, ye in zip(self.fig.axes, self.y_arrays, self.x_expt_arrays, self.y_expt_arrays):
			if y is not None:
				ax.plot(self.x_array,y)
			if ye is not None:
				ax.plot(xe,ye,color='k')
			
		self.draw_fig()
	
	def draw_fig(self):
		self.fig.tight_layout()
		self.canvas.draw()
		
	def AddPlotType(self,plotindex):
		# if not a currently selected output type, add it to the list of output types
		self.plot_output_indices.append(plotindex)
		self.plot_outputs.append(OutputTypes[plotindex])
		
		self.y_arrays.append(None)
		self.x_expt_arrays.append(None)
		self.y_expt_arrays.append(None)
		
		### Re-size plot axes to add another subplot
		self.OnAddAxes(1)

##
## Other global functions
##

def read_CSV(filename,spacing=0,columns=2):
	""" Read an n-column csv file """
	f=open(filename,'U')
	fid=[]
	for line in f:
		fid.append(line)
	f.close()
	# -1: ignore last (blank) line in csv file (lecroy)
	fid = fid[spacing:-1]
	inData=csv.reader(fid,delimiter=',')
	# spacing : skips lines if needed (e.g., on oscilloscope files)
	
	data=[]
	for i in range(0,columns): data.append([])
	
	for row in inData:
		for i in range(0,columns):
			data[i].append(float(row[i]))
	
	for i in range(0,columns):
		data[i] = np.array(data[i])
		
	return data		
		
app = wx.App(redirect=False)
frame = ElecSus_GUI_Frame(None,"ElecSus GUI BareBones")
#frame.Maximize()
frame.Show()
app.MainLoop()