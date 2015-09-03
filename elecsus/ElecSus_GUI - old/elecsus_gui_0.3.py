"""

ElecSus GUI

A GUI based on wxpython for ElecSus, intended to augment/replace the 
runcard method of calling ElecSus.

More blurb here...

Requirements:
	python, matplotlib, numpy, scipy
	versions unknown, 
	
	tested on:
	Windows 8.1
		python 2.7.6 - 64-bit
		wxpython 3.0.2.0
		matplotlib 1.4.3
		numpy 1.9.2
		scipy 0.15.1
	
	
	Requires
	--------
	
	wxpython version 2.8 (Note: Enthought Canopy on Windows only comes with version 2.8 as of 2015-08-09)
		In newer versions, there are bugs with either matplotlib or wxpython that cause figures 
		not to fill the entire panel.
		
	
LICENSE info

James Keaveney, Mark Zentile and co-authors
2011-2015
"""

#!/usr/bin/env python
import matplotlib
matplotlib.use('WxAgg')
import pylab
pylab.ioff()

import wx
# Ignore warnings
wx.Log.SetLogLevel(0)
import wx.lib.scrolledpanel as scrolled
from wx.lib.agw.floatspin import FloatSpin


import os
import sys
import csv
import time

import numpy as np
import matplotlib.pylab as plt

# Matplotlib/wx integration
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg as Toolbar

## fits run in a separate thread
import threading
# define a new event type which will be triggered when fitting is completed
import wx.lib.newevent
FitCompleteEvent, EVT_FIT_COMPLETE = wx.lib.newevent.NewEvent()

# import elecsus modules
import elecsus_new
from libs import spectra
from libs import NOTICE
from libs.durhamcolours import cols as durhamcols


#replace default matplotlib text and color sequence with durham colours
from matplotlib import rc
rc('text', usetex=False)
rc('font',**{'family':'serif', 'size':14})
rc('lines', linewidth=2)
rc('axes', color_cycle=durhamcols)

# Button size
BtnSize = 30

## IMPORTANT! Master list of all output types that will be referenced for dynamic plotting
OutputTypes = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'Alpha Plus', 'Alpha Minus', 'N Plus', 'N Minus', 'Phi']
OutputTypes_index = [0,1,2,3,4,4,5,5,6,6,7]
OutputPlotTypes = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix, Iy', 'Alpha Plus/Minus', 'N Plus/Minus', 'Phi']

fixed_parameterlist = ['Element','D-line','Constrain Dopp./Density Temps']
element_list = ['Na', 'K', 'Rb', 'Cs']
D_line_list = ['D1', 'D2']

fittable_parameterlist = ['Bfield','Temperature','Cell length','Shift','Additional-Broadening','Theta-0','LC Polarisation','Doppler-Temperature','Rb-85','K40','K41']
units_parameterlist = [' [G]',' [C]',' [mm]',' [MHz]',' [MHz]',' [deg]',' [%]',' [C]',' [%]',' [%]',' [%]']
defaultvals_parameterlist = [0,20,75,0,0,0,50,20,72.17,1,0]
defaultvals_increments = [1, 0.5, 1, 1, 1, 10, 5, 0.5, 5, 5, 5]
detuning_defaults = [-10, 10, 5000]
detuning_increments = [1,1,1000]

class FittingThread(threading.Thread):
	""" 
	Fitting takes a long time, so we run the fit in another thread, leaving the GUI active in the meantime.
	When the fitting completes, this thread triggers an event in the GUI that then completes the rest of the
	fitting routine - i.e. updates the plot, writes text ... etc
	"""
	def __init__(self,parent):
		threading.Thread.__init__(self)
		
		self.mainwin = parent
		#self.start() ## this calls the run() method - called from the GUI
	
	def run(self):
		""" Run fitting thread """
		mainwin = self.mainwin
		try:
			mainwin.opt_params, mainwin.rms = elecsus_new.fit_data((mainwin.x_fit_array,mainwin.y_fit_array),\
								mainwin.params,mainwin.fit_bools_ordered,\
								mainwin.fit_datatype,mainwin.fit_algorithm,\
								**mainwin.advanced_fitoptions)
								

			
		except RuntimeError: 
			## specify which kind of exception! - maybe add custom exception to cancel fitting?
			self.fitting_dlg.Destroy()
			
			raise RuntimeError
			dlg = wx.MessageDialog(mainwin,"Fitting encountered an error! Check the errors panel for details.","Fitting Error",style=wx.OK|wx.CENTRE|wx.ICON_ERROR)
			
			if dlg.ShowModal() == wx.ID_OK:
				pass

		print 'Fit completed'
		evt = FitCompleteEvent()
		wx.PostEvent(mainwin, evt)			

class ProgressThread(threading.Thread):
	""" update the progress bar continually to show user something is happening ..."""
	def __init__(self,parent):
		threading.Thread.__init__(self)
		self.mainwin = parent
		
	def run(self):
		while self.mainwin.already_fitting:
			for i in range(100):
				self.mainwin.fitting_dlg.Update(i)
				time.sleep(0.1)
		print 'Quitting progress bar update'
		
#class ProgressDialog_OnOff(wx.ProgressDialog):
#	def __init__(self,
	
class OptionsPanel(scrolled.ScrolledPanel):
	def __init__(self, parent, mainwin, paneltype):
		""" 
		Most of the options are the same for theory and fitting, but there are some differences. 
		Hence, the same panel class with an argument 'paneltype' which will set what is displayed.
		
		This ScrolledPanel class will automatically add horizontal and vertical scroll bars when required, and dynamically turns them on and off as window is resized. Neat!
		""" 
		# mainwin is the main panel so we can bind buttons to actions in the main frame """
		
		scrolled.ScrolledPanel.__init__(self, parent)
		#wx.Panel.__init__(self, parent)
		
		#self.parent = parent
		self.paneltype = paneltype
		
		# Fitting only:
		if self.paneltype == 'Fit':
			# Import data from csv
			
			# Booleans (check boxes) for fitting
			self.fit_paramlist_bools = \
				[ wx.CheckBox(self, label='') for parameter in fittable_parameterlist ]
						
			# Fitting algorithm selection
			# ties in with menu items
			fit_ML = wx.RadioButton(self,label="Marquardt-Levenberg", style=wx.RB_GROUP)
			fit_RR = wx.RadioButton(self,label="Random-Restart")
			fit_SA = wx.RadioButton(self,label="Simulated Annealing")
			
			self.Bind(wx.EVT_RADIOBUTTON, mainwin.OnFitTypeChange, fit_ML)
			self.Bind(wx.EVT_RADIOBUTTON, mainwin.OnFitTypeChange, fit_RR)
			self.Bind(wx.EVT_RADIOBUTTON, mainwin.OnFitTypeChange, fit_SA)
		
			# Run Fit Button
			RunFitButton = wx.Button(self,wx.ID_ANY, 'Run Fit', size=(120,1.5*BtnSize))
			self.Bind(wx.EVT_BUTTON, mainwin.OnFitButton, RunFitButton)
			
		# Theory only:
		elif self.paneltype == 'Theory':
			# Detuning range selection
			Detuning_Labels = [ wx.StaticText(self,wx.ID_ANY,lab) for lab in ["Start [GHz]", "Stop [GHz]", "No. of points"]]
			self.DetuningCtrl = [ FloatSpin(self,value=str(defval),increment=definc,size=(80,-1),digits=2) for defval,definc in zip(detuning_defaults,detuning_increments) ]

			# Calculate spectrum
			ComputeButton = wx.Button(self,wx.ID_ANY, 'Compute Spectrum', size=(120,1.5*BtnSize))
			self.Bind(wx.EVT_BUTTON, mainwin.OnComputeButton, ComputeButton)

			#dummy variable
			self.fit_paramlist_bools = [0]*len(fittable_parameterlist)
			
		# Common options to both experiment and theory
		ImportButton = wx.Button(self,wx.ID_OPEN,label="Import Data", size=(120,1.5*BtnSize))
		self.Bind(wx.EVT_BUTTON,mainwin.OnFileOpen,ImportButton)
			
		fixed_paramlist_labels = [ wx.StaticText(self,wx.ID_ANY,fixed_param) for fixed_param in fixed_parameterlist ]
		self.fixed_paramlist_inputs = [ \
		wx.ComboBox(self,wx.ID_ANY,choices=element_list,style=wx.CB_READONLY,size=(80,-1)), \
		wx.ComboBox(self,wx.ID_ANY,choices=D_line_list, style=wx.CB_READONLY,size=(80,-1)), \
		wx.CheckBox(self, label="")]
		self.fixed_paramlist_inputs[0].SetSelection(0)
		self.fixed_paramlist_inputs[1].SetSelection(0)
		self.fixed_paramlist_inputs[2].SetValue(True)
		
		# create list of parameters - labels and spin-control boxes
		fit_paramlist_labels = [ wx.StaticText(self,wx.ID_ANY,fit_param+unit) for fit_param,unit in zip(fittable_parameterlist,units_parameterlist) ]
		self.fit_paramlist_inputs = [ \
		FloatSpin(self,value=str(defval),increment=definc,size=(80,-1),digits=2) for defval,definc in zip(defaultvals_parameterlist,defaultvals_increments) ]

		# Don't need to bind SpinCtrl/ComboBox inputs to actions - the values are 
		# read whenever 'Compute' or 'Fit' buttons are pressed

		### Layout panel in sizers
		
		# main sizer for the panel
		panel_sizer = wx.BoxSizer(wx.VERTICAL)
		
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		label_sizer = wx.BoxSizer(wx.HORIZONTAL)
		label_sizer.Add((10,-1),0,wx.EXPAND)
		elem_label = wx.StaticText(self,wx.ID_ANY,"Select element and D-line")
		font = wx.Font(12,wx.DEFAULT, wx.NORMAL,wx.NORMAL)
		elem_label.SetFont(font)
		label_sizer.Add(elem_label)
		label_sizer.Add((10,-1),1,wx.EXPAND)
		panel_sizer.Add(label_sizer)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		# Add common elements first:
		
		for label, input in zip(fixed_paramlist_labels,self.fixed_paramlist_inputs):
			hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
			hor_sizer.Add((30,-1),0,wx.EXPAND)
			hor_sizer.Add(label,0,wx.EXPAND)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(input,0,wx.EXPAND)
			hor_sizer.Add((10,-1),0,wx.EXPAND)
			
			panel_sizer.Add(hor_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,3),0,wx.EXPAND)
		
		# vertical space
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		label_sizer = wx.BoxSizer(wx.HORIZONTAL)
		label_sizer.Add((10,-1),0,wx.EXPAND)
		if self.paneltype == 'Theory':
			parameter_label_text = 'Computation parameters'
		elif self.paneltype == 'Fit':
			parameter_label_text = 'Initial parameters for fit'
		paramlabeltextbox = wx.StaticText(self,wx.ID_ANY,parameter_label_text)
		font = wx.Font(12,wx.DEFAULT, wx.NORMAL,wx.NORMAL)
		paramlabeltextbox.SetFont(font)
		label_sizer.Add(paramlabeltextbox,0,wx.EXPAND)
		label_sizer.Add((10,-1),1,wx.EXPAND)
		if self.paneltype == 'Fit':
			floatlabel = wx.StaticText(self,wx.ID_ANY,'Float?')
			label_sizer.Add(floatlabel,0,wx.ALIGN_BOTTOM)
			label_sizer.Add((10,-1),0,wx.EXPAND)
		panel_sizer.Add(label_sizer,0,wx.EXPAND)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		
		for label,input,boolbox in zip(fit_paramlist_labels,self.fit_paramlist_inputs,self.fit_paramlist_bools):
			hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
			hor_sizer.Add((30,-1),0,wx.EXPAND)
			hor_sizer.Add(label,0,wx.EXPAND)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(input,0,wx.EXPAND)
			if self.paneltype == 'Fit':
				hor_sizer.Add((20,-1),0,wx.EXPAND)
				hor_sizer.Add(boolbox)
			hor_sizer.Add((10,-1),0,wx.EXPAND)
			panel_sizer.Add(hor_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,3),0,wx.EXPAND)
		
		if self.paneltype == 'Fit':
			# vertical space
			panel_sizer.Add((-1,5),0,wx.EXPAND)
			label_sizer = wx.BoxSizer(wx.HORIZONTAL)
			label_sizer.Add((20,-1),0,wx.EXPAND)
			parameter_label_text = 'Fit Algorithm'
			paramlabeltextbox = wx.StaticText(self,wx.ID_ANY,parameter_label_text)
			font = wx.Font(12,wx.DEFAULT, wx.NORMAL,wx.NORMAL)
			paramlabeltextbox.SetFont(font)
			label_sizer.Add(paramlabeltextbox,0,wx.EXPAND)
			label_sizer.Add((20,-1),1,wx.EXPAND)
			panel_sizer.Add(label_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,10),0,wx.EXPAND)
		
			for radiobtn in [fit_ML, fit_RR, fit_SA]:
				hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
				hor_sizer.Add((40,-1),0,wx.EXPAND)
				hor_sizer.Add(radiobtn,0,wx.EXPAND)
				hor_sizer.Add((20,-1),1,wx.EXPAND)
				panel_sizer.Add(hor_sizer,0,wx.EXPAND)
				panel_sizer.Add((-1,4),0,wx.EXPAND)		
		elif self.paneltype == 'Theory':
			# vertical space
			panel_sizer.Add((-1,5),0,wx.EXPAND)
			label_sizer = wx.BoxSizer(wx.HORIZONTAL)
			label_sizer.Add((20,-1),0,wx.EXPAND)
			parameter_label_text = 'Detuning Range'
			paramlabeltextbox = wx.StaticText(self,wx.ID_ANY,parameter_label_text)
			font = wx.Font(12,wx.DEFAULT, wx.NORMAL,wx.NORMAL)
			paramlabeltextbox.SetFont(font)
			label_sizer.Add(paramlabeltextbox,0,wx.EXPAND)
			label_sizer.Add((20,-1),1,wx.EXPAND)
			panel_sizer.Add(label_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,10),0,wx.EXPAND)
			
			for label, input in zip(Detuning_Labels, self.DetuningCtrl):
				hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
				hor_sizer.Add((40,-1),0,wx.EXPAND)
				hor_sizer.Add(label,0,wx.EXPAND)
				hor_sizer.Add((20,-1),1,wx.EXPAND)
				hor_sizer.Add(input,0,wx.EXPAND)
				hor_sizer.Add((10,-1),0,wx.EXPAND)
				panel_sizer.Add(hor_sizer,0,wx.EXPAND)
				panel_sizer.Add((-1,3),0,wx.EXPAND)				
			
			
			
		# vertical space
		panel_sizer.Add((-1,5),1,wx.EXPAND)
		
		hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
		hor_sizer.Add((30,-1),0,wx.EXPAND)
		hor_sizer.Add(ImportButton,0,wx.EXPAND)
		hor_sizer.Add((20,-1),1,wx.EXPAND)
		if self.paneltype == 'Theory':
			hor_sizer.Add(ComputeButton,0,wx.EXPAND)
		elif self.paneltype == 'Fit':
			hor_sizer.Add(RunFitButton,0,wx.EXPAND)
		hor_sizer.Add((30,-1),0,wx.EXPAND)
		panel_sizer.Add(hor_sizer,0,wx.EXPAND)
				
		panel_sizer.Add((-1,2),0,wx.EXPAND)
		
		self.SetSizer(panel_sizer)
		self.SetupScrolling()
		self.Layout()			
				
class PlotToolPanel(wx.Panel):
	def __init__(self, parent, mainwin, ID):
		""" mainwin is the main panel so we can bind buttons to actions in the main frame """
		wx.Panel.__init__(self, parent)
		
		self.fig = plt.figure(ID,facecolor=(240./255,240./255,240./255))
		
		#self.ax = self.fig.add_subplot(111)
		
		# create the wx objects to hold the figure
		self.canvas = FigureCanvasWxAgg(self, wx.ID_ANY, self.fig)
		self.toolbar = Toolbar(self.canvas) #matplotlib toolbar (pan, zoom, save etc)
		
		# Create vertical sizer to hold figure and toolbar - dynamically expand with window size
		plot_sizer = wx.BoxSizer(wx.VERTICAL)
		plot_sizer.Add(self.canvas, 1, wx.TOP|wx.LEFT|wx.GROW)
		plot_sizer.Add(self.toolbar, 0, wx.EXPAND)

		mainwin.figs.append(self.fig)
		mainwin.fig_IDs.append(ID) # use an ID number to keep track of figures
		mainwin.canvases.append(self.canvas)

		# display some text in the middle of the window to begin with
		self.fig.text(0.5,0.5,'ElecSus GUI\n\nVersion 0.9a\n\nTo get started, use the panel on the right\nto either Compute a spectrum or Import some data...', ha='center',va='center')
		self.fig.hold(False)
		
		self.SetSizer(plot_sizer)
		self.Fit()

class StatusPanel(scrolled.ScrolledPanel):
	def __init__(self, parent, mainwin, ID):
		""" mainwin is the main panel so we can bind buttons to actions in the main frame """
		scrolled.ScrolledPanel.__init__(self, parent)
	
		statusTitle = wx.StaticText(self, wx.ID_ANY, ID)
		font = wx.Font(12,wx.DEFAULT, wx.NORMAL,wx.NORMAL)
		statusTitle.SetFont(font)

		self.StatusTextBox = wx.TextCtrl(self,wx.ID_ANY,"",style=wx.TE_READONLY|wx.TE_MULTILINE)
		#self.StatusTextBox.Size.SetHeight(500)
		
		self.SaveBtn = wx.Button(self,wx.ID_ANY,"Save Text to file",size=(100,-1))
		self.Bind(wx.EVT_BUTTON,self.OnSave,self.SaveBtn)
		
		
		panel_sizer = wx.BoxSizer(wx.VERTICAL)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		panel_sizer.Add(statusTitle,0,wx.EXPAND|wx.LEFT,border=20)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		panel_sizer.Add(self.StatusTextBox,1,wx.EXPAND|wx.LEFT|wx.RIGHT,border=40)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		panel_sizer.Add(self.SaveBtn,0,wx.EXPAND|wx.LEFT|wx.RIGHT,border=40)
		panel_sizer.Add((-1,20),0,wx.EXPAND)
		

		self.SetSizer(panel_sizer)
		self.SetupScrolling()
		self.Layout()
	
	def OnSave(self,event):
		pass
	
	def write(self,textstring):
		self.StatusTextBox.AppendText(textstring)
		
class AdvancedFitOptions(wx.Dialog):
	def __init__(self,parent,title,id):
		""" Dialog box for selecting advanced fit options ... """
		
		wx.Dialog.__init__(self,parent,id,title,size=(300,300))
				
		# main sizer for this dialog box
		panel_sizer = wx.BoxSizer(wx.VERTICAL)
		
		panel_blurb = wx.StaticText(self,wx.ID_ANY,"These are the options given to scipy.curve_fit (leastsq) when fitting data. \n\nMost of the time these will not need to be \nchanged from their defaults. However, if fitting \nis returning errors then consult the scipy \ndocumentation and adjust these parameters accordingly.\n ",size=(260,-1),style=wx.ALIGN_CENTRE_HORIZONTAL)
		panel_blurb.Wrap(260)
		
		self.fitopt_labels = ['Max Iterations', 'Max Tolerance (x 1e-8)']
		self.fitopt_argnames = ['maxfev', 'ftol']
		fitopt_defaults = [5000, 15]
		fitopt_increments = [100, 1]
		self.fitopt_scaling = [1, 1e-9]
		fitopt_texts = [ wx.StaticText(self,wx.ID_ANY,label) for label in self.fitopt_labels ]
		self.fitopt_ctrl = [ wx.SpinCtrl(self,value=str(defval),size=(80,-1)) for defval,definc in zip(fitopt_defaults, fitopt_increments) ]
		
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		panel_sizer.Add(panel_blurb,0,wx.LEFT|wx.RIGHT,border=20)
		panel_sizer.Add((-1,15),0,wx.EXPAND)
		panel_sizer.Add(wx.StaticLine(self,-1,size=(-1,1),style=wx.LI_HORIZONTAL),0,wx.EXPAND|wx.LEFT|wx.RIGHT,border=20
		)
		panel_sizer.Add((-1,15),0,wx.EXPAND)

		for static,ctrl in zip(fitopt_texts, self.fitopt_ctrl):
			hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
			hor_sizer.Add(static,0,wx.EXPAND|wx.LEFT,border=20)
			hor_sizer.Add((10,-1),1,wx.EXPAND)
			hor_sizer.Add(ctrl,0,wx.EXPAND|wx.RIGHT,border=20)
			panel_sizer.Add(hor_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,5),0,wx.EXPAND)
			
		# ok and cancel buttons
		btnbar = self.CreateButtonSizer(wx.OK|wx.CANCEL)

		#panel_sizer.Add((-1,10),0,wx.EXPAND)
		#panel_sizer.Add(wx.StaticLine(self,-1,size=(-1,1),style=wx.LI_HORIZONTAL),0, wx.EXPAND|wx.LEFT|wx.RIGHT,border=20)

		panel_sizer.Add((-1,10),1,wx.EXPAND)
		panel_sizer.Add(btnbar,0,wx.ALIGN_CENTER)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		
		self.SetSizer(panel_sizer)
		self.Layout()
		
	def return_all_options(self):
		# get list of ctrl elements
		opt_dict = dict([(label,input.GetValue()*scaling) for label,input,scaling in zip(self.fitopt_argnames,self.fitopt_ctrl,self.fitopt_scaling)])
		print 'Advanced fit options dictionary: ',opt_dict
		return opt_dict
		
class ElecSus_GUI_Frame(wx.Frame):
	""" Main window """
	def __init__(self,parent,title):
		wx.Frame.__init__(self,None,title=title,size=(1200,800))
		
		# Set icons for top-left of frame, alt-tab window ...
		frame_icon = wx.IconBundle()
		frame_icon.AddIconFromFile('./images/elecsus_t_group.ico', wx.BITMAP_TYPE_ANY)
		self.SetIcons(frame_icon)

		#if the window is closed, exit
		self.Bind(wx.EVT_CLOSE,self.OnExit)

		self.panel = wx.Panel(self)
		
		self.panel.SetBackgroundColour(wx.Colour(240,240,240))
		
		self._init_default_values()	
		self._init_plot_defaults()
		self._init_panels()
		self._init_menus()
		
		# redirect stdout (command line text) to status box
		sys.stdout = self.StatusPanel
		sys.stderr = self.ErrorPanel
		
		# Create initially blank set of axes
		#self.OnCreateAxes(self.figs[0],self.canvases[0])
		
		## Bind the event EVT_FIT_COMPLETE to function
		## This executes in the main thread once the fitting thread (separate from the main thread) completes
		self.Bind(EVT_FIT_COMPLETE,self.OnFitCompleted)
		
		print 'here'
				
	def _init_default_values(self):
		""" Initialise default values for various things ... """
		
		self.figs = []
		self.canvases = []
		self.fig_IDs = []
		
		self.hold = False
		
		self.already_fitting = False

		# initialise advanced fit options dictionary
		dlg = AdvancedFitOptions(self,"Advanced Fit Options",wx.ID_ANY)
		self.advanced_fitoptions = dlg.return_all_options()
		dlg.Destroy()		
				
		
	def _init_plot_defaults(self):
		""" List of default values for all plots. Theory plots can be up to 7 panels, each of which can be turned on/off """
		
		self.plot_outputs = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix, Iy', 'Alpha Plus/Minus', 'N Plus/Minus']
		self.plot_ylabels = ['Transmission', 'S1', 'S2', 'S3', '$I_x$, $I_y$', r'$\alpha^{\pm}$', r'$n^\pm -1$']
		self.plot_output_indices = [OutputPlotTypes.index(po) for po in self.plot_outputs]
		
		self.xrange = [detuning_defaults[0],detuning_defaults[1]] # detuning range, in GHz
		self.npoints = detuning_defaults[2] # number of detuning points to calculate
		
		# default data for plots - blank lists
		self.x_array = np.linspace(self.xrange[0],self.xrange[1],self.npoints)
		self.y_arrays = [None, None, None, None, None, None, None, None]
		self.x_expt_arrays = [None, None, None, None, None, None, None, None]
		self.y_expt_arrays = [None, None, None, None, None, None, None, None]
		self.x_fit_array, self.y_fit_array = None,None
		
		# plot data visible when set to True
		self.display_theory_curves = [False,False,False,False,False,False,False]
		self.display_expt_curves = [False,False,False,False,False,False,False]
		
		# ...
	
	def _init_menus(self):
		""" Initialise menu bar items """
		
		# Create menuBar object
		menuBar = wx.MenuBar()
		
		# File
		fileMenu = wx.Menu()
		fM_open = fileMenu.Append(wx.ID_OPEN, "&Open\tCtrl+O", "Open file for plotting and/or fitting.")
		self.Bind(wx.EVT_MENU, self.OnFileOpen, fM_open)
		fileMenu.AppendSeparator()
		fm_saveplot = fileMenu.Append(wx.ID_ANY, "Save Plot as Image", "Save Data")
		self.Bind(wx.EVT_MENU, self.OnSaveFig, fm_saveplot)
		fm_savecsv = fileMenu.Append(wx.ID_SAVE, "E&xport CSV Data\tCtrl+S", "Save CSV Data")
		self.Bind(wx.EVT_MENU, self.OnSaveCSVData, fm_savecsv)
		fileMenu.AppendSeparator()
		fM_exit = fileMenu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", "Close window and exit program.")
		self.Bind(wx.EVT_MENU, self.OnExit, fM_exit)
		#
		
		# About
		aboutMenu = wx.Menu()
		aM_aboutthis = aboutMenu.Append(wx.ID_ABOUT, "&About this program", "About this program.")
		self.Bind(wx.EVT_MENU, self.OnAboutThis, aM_aboutthis)
		aM_aboutelecsus = aboutMenu.Append(wx.ID_ANY, "About &ElecSus", "About ElecSus")
		self.Bind(wx.EVT_MENU, self.OnAboutElecSus, aM_aboutelecsus)

		# View
		viewMenu = wx.Menu()
		vM_dummy = viewMenu.Append(wx.ID_ANY, "&Dummy", "Dummy")
		# Select plot types to display

		# Theory Plot
		theoryplotMenu = wx.Menu()
		tpM_plotholdon = theoryplotMenu.Append(wx.ID_ANY, "&Hold data on plot update", "Select whether to hold or clear current plot data on updating the figure", kind=wx.ITEM_CHECK)
		tpM_plotholdon.Check(self.hold)
		self.Bind(wx.EVT_MENU, self.OnPlotHold, tpM_plotholdon)
		tpM_clearplot = theoryplotMenu.Append(wx.ID_ANY, "&Clear current plot data", "Clear current plot data on all axes")
		self.Bind(wx.EVT_MENU, self.OnClearPlot, tpM_clearplot)
		
		tpM_grid = theoryplotMenu.Append(wx.ID_ANY, "&Grid on axes", "Grid", kind=wx.ITEM_CHECK)
		self.Bind(wx.EVT_MENU, self.OnGridToggle, tpM_grid)

		self.showTplotsSubMenu = wx.Menu()
		id_S0 = 2000
		id_S1 = 2001
		id_S2 = 2002
		id_S3 = 2003
		id_IxIy = 2004
		id_alpha = 2005
		id_n = 2006
		show_S0 = self.showTplotsSubMenu.AppendCheckItem(id_S0,"&Transmission (S0)")
		show_S1 = self.showTplotsSubMenu.AppendCheckItem(id_S1,"S&1")
		show_S2 = self.showTplotsSubMenu.AppendCheckItem(id_S2,"S&2")
		show_S3 = self.showTplotsSubMenu.AppendCheckItem(id_S3,"S&3")
		show_IxIy = self.showTplotsSubMenu.AppendCheckItem(id_IxIy,"&Ix/Iy")
		show_alpha = self.showTplotsSubMenu.AppendCheckItem(id_alpha,"&Alpha +/-")
		show_n = self.showTplotsSubMenu.AppendCheckItem(id_n,"&Refractive Index +/-")
		# bind event to check box selections
		for checkitem in [show_S0, show_S1, show_S2, show_S3, show_IxIy, show_alpha, show_n]:
			self.Bind(wx.EVT_MENU, self.OnShowTplots, checkitem)
		# add to parent menu
		tpM_showTplots = theoryplotMenu.AppendMenu(wx.ID_ANY,"&Theory Curves to Show", self.showTplotsSubMenu)
		
		self.showEplotsSubMenu = wx.Menu()
		id_S0 = 3000
		id_S1 = 3001
		id_S2 = 3002
		id_S3 = 3003
		id_IxIy = 3004
		id_alpha = 3005
		id_n = 3006
		show_S0 = self.showEplotsSubMenu.AppendCheckItem(id_S0,"&Transmission (S0)")
		show_S1 = self.showEplotsSubMenu.AppendCheckItem(id_S1,"S&1")
		show_S2 = self.showEplotsSubMenu.AppendCheckItem(id_S2,"S&2")
		show_S3 = self.showEplotsSubMenu.AppendCheckItem(id_S3,"S&3")
		show_IxIy = self.showEplotsSubMenu.AppendCheckItem(id_IxIy,"&Ix/Iy")
		show_alpha = self.showEplotsSubMenu.AppendCheckItem(id_alpha,"&Alpha +/-")
		show_n = self.showEplotsSubMenu.AppendCheckItem(id_n,"&Refractive Index +/-")
		# bind event to check box selections
		for checkitem in [show_S0, show_S1, show_S2, show_S3, show_IxIy, show_alpha, show_n]:
			self.Bind(wx.EVT_MENU, self.OnShowEplots, checkitem)
		# add to parent menu
		tpM_showEplots = theoryplotMenu.AppendMenu(wx.ID_ANY,"&Experimental Curves to Show", self.showEplotsSubMenu)

		
		# Experimental Plot
		exptplotMenu = wx.Menu()
		epM_plotholdon = exptplotMenu.Append(wx.ID_ANY, "&Hold data on plot update", "Select whether to hold or clear current plot data on updating the figure", kind=wx.ITEM_CHECK)
		epM_clearplot = exptplotMenu.Append(wx.ID_ANY, "&Clear current plot data", "Clear current plot data on all axes")
		self.Bind(wx.EVT_MENU, self.OnClearPlot,epM_clearplot)
		
		epM_grid = exptplotMenu.Append(wx.ID_ANY, "&Grid on axes", "Grid", kind=wx.ITEM_CHECK)
		self.Bind(wx.EVT_MENU, lambda evt, temp=self.figs[1]: self.OnGridToggle(evt, temp),epM_grid)
		
		

		
		# Fit
		fitMenu = wx.Menu()
		fitM_dummy = fitMenu.Append(wx.ID_ANY, "&Dummy", "Dummy")
		
		# Select fitting type
		self.fittypeSubMenu = wx.Menu()
		#define some id tags that we will use when items are selected
		id_ML = 1001
		id_RR = 1002
		id_SA = 1003
		fit_ML = self.fittypeSubMenu.AppendRadioItem(id_ML, "&Marquardt-Levenberg", "Use ML Fitting")
		fit_RR = self.fittypeSubMenu.AppendRadioItem(id_RR, "&Random Restart", "Use RR Fitting")
		fit_SA = self.fittypeSubMenu.AppendRadioItem(id_SA, "&Simulated Annealing", "Use SA Fitting")
		
		#Bind action to each item
		self.Bind(wx.EVT_MENU, self.OnFitTypeChange, fit_ML)
		self.Bind(wx.EVT_MENU, self.OnFitTypeChange, fit_RR)
		self.Bind(wx.EVT_MENU, self.OnFitTypeChange, fit_SA)
		
		# by default, select ML fitting
		self.fittypeSubMenu.Check(id_ML, True)
		self.OnFitTypeChange(1)

		# Add sub-menu to fit menu
		fitM_fittype = fitMenu.AppendMenu(wx.ID_ANY, "Fit &Type", self.fittypeSubMenu)
		
		fitM_advanced = fitMenu.Append(wx.ID_ANY, "&Advanced Fit Options...", "Advanced Fit Options")
		self.Bind(wx.EVT_MENU, self.OnAdvancedOptions, fitM_advanced)
		
		#initialise
		#
		
			
		#
		# ... other menu items to add...
		# 
		# Advanced fitting - setting tolerances, maximum processor count, ... etc
		# Data processing - binning and/or smoothing data before fitting
		# Residual analysis - fit gaussian to residuals etc
		#		-- Use normalised/raw residuals (requires importing data with errorbars)
		# Main plot - set axes limits
		#
		
		# Add menu items to menu bar
		menuBar.Append(fileMenu, "&File")
		menuBar.Append(viewMenu, "&View")
		menuBar.Append(theoryplotMenu, "&Main Plot")
		menuBar.Append(exptplotMenu, "&Residuals Plot")
		menuBar.Append(fitMenu, "F&it")
		menuBar.Append(aboutMenu, "&About")
	
		# Add Menu Bar to the main panel
		self.SetMenuBar(menuBar)
		
	def _init_panels(self):
		""" Initialise panel with matplotlib window, buttons, text boxes etc. Doesn't really need to be in a separate function, but makes it easier to find stuff... """
		
		# Sizer constructs are used for placement of all GUI elements. 
		# This makes the whole window scalable in a predictable way.


	## Create plot part of the window
		
		## create plot in a notebook-style for adding more tabs later on
		PlotTabs = wx.Notebook(self.panel)
		
		# The plot panels
		self.T_Panel = PlotToolPanel(PlotTabs,self,'Theory')
		self.E_Panel = PlotToolPanel(PlotTabs,self,'Experiment')
		
		# Text tabs - for stdout, stderr messages and fitting information generated after fitting data (optimum parameters etc)
		self.StatusPanel = StatusPanel(PlotTabs,self,'Status Information')
		self.FitInformation = StatusPanel(PlotTabs,self,'Fitting Information')
		self.ErrorPanel = StatusPanel(PlotTabs,self, 'Error Information')
		self.ErrorPanel.write("If this is the only thing displayed here, the program is working well...\n\n")
		
		# Add all tabs to the tab bar
		PlotTabs.AddPage(self.T_Panel, "Main Plot")
		PlotTabs.AddPage(self.E_Panel, "Fit / Residuals Plot")
		PlotTabs.AddPage(self.StatusPanel, "Status Panel")
		PlotTabs.AddPage(self.FitInformation, "Fitting Information")
		PlotTabs.AddPage(self.ErrorPanel, "Error Information")
		
		# Add the Tab bar to the main panel sizer
		plot_sizer = wx.BoxSizer(wx.VERTICAL)
		plot_sizer.Add(PlotTabs,1, wx.EXPAND)
		
	## Create button part of the window
		
		# elecsus and JQC logos at the top
		#elecsuslogo = wx.Image('images/elecsus.ico',wx.BITMAP_TYPE_ANY)
		#elecsuslogo.Rescale(108/2,138/2)
		#elecsus_bmp = wx.StaticBitmap(self.panel,wx.ID_ANY,wx.BitmapFromImage(elecsuslogo),size=(108/2,-1))
		jqclogo = wx.Image('images/jqc-logo.png',wx.BITMAP_TYPE_ANY)
		#jqclogo.Rescale(191*3/4,84*3/4)
		jqc_bmp = wx.StaticBitmap(self.panel,wx.ID_ANY,wx.BitmapFromImage(jqclogo),size=(191,-1))
		
		image_sizer = wx.BoxSizer(wx.HORIZONTAL)
		image_sizer.Add((10,-1),1,wx.EXPAND)
		image_sizer.Add((80,-1),0,wx.EXPAND)
		#image_sizer.Add(elecsus_bmp,0,wx.EXPAND)
		#image_sizer.Add((10,-1),0,wx.EXPAND)
		image_sizer.Add(jqc_bmp,0,wx.EXPAND)
		image_sizer.Add((80,-1),0,wx.EXPAND)
		image_sizer.Add((10,-1),1,wx.EXPAND)
		
	## Menus are tabbed into theory and fitting sections
	## use wx.Notebook for this
		TabPanel = wx.Notebook(self.panel)
		
		# add tabs
		self.ThyPanel = OptionsPanel(TabPanel,self,'Theory')
		self.FitPanel = OptionsPanel(TabPanel,self,'Fit')
		
		TabPanel.AddPage(self.ThyPanel, "Theory Settings")
		TabPanel.AddPage(self.FitPanel, "Fit Settings")
		
		tab_sizer = wx.BoxSizer(wx.HORIZONTAL)
		#tab_sizer.Add((40,-1),0,wx.EXPAND)
		tab_sizer.Add(TabPanel, 1, wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		#tab_sizer.Add((40,-1),0,wx.EXPAND)
		
		
	## Create sizer for right-hand side of panel
		button_sizer = wx.BoxSizer(wx.VERTICAL)
		button_sizer.Add((-1,10),0,wx.EXPAND)
		button_sizer.Add(image_sizer,0,wx.EXPAND)
		button_sizer.Add((-1,5),0,wx.EXPAND)
		button_sizer.Add(tab_sizer,1,wx.EXPAND)
		button_sizer.Add((-1,5),0,wx.EXPAND)
		
		
	## Put plot part and button part together
		main_sizer = wx.BoxSizer(wx.HORIZONTAL)
		main_sizer.Add(plot_sizer,1,wx.EXPAND)
		main_sizer.Add(wx.StaticLine(self.panel,-1,size=(1,-1),style=wx.LI_VERTICAL),0,wx.EXPAND)
		main_sizer.Add(button_sizer,0,wx.EXPAND)

		self.panel.SetSizer(main_sizer)
		self.panel.Layout()
		
		
		##
		self.T_Panel.Layout()
		self.T_Panel.Refresh()
		#self.E_Panel.Layout()
		self.update()
		#self.Refresh()
		#self.SendSizeEvent()
		#TabPanel.SendSizeEvent()
		
#		
## Actions for events
#
	
##
#### General Actions - Close window, Save Data etc...
##

	def OnExit(self,event):
		self.Destroy()
		# explicitly close all figures (bug with matplotlib and wx??)
		plt.close('all')
		app.ExitMainLoop()
		
	def OnAboutThis(self,event):
		## do stuff
		dlg = wx.MessageDialog(self, "Blurb about this program", "About", wx.OK)
		if dlg.ShowModal() == wx.ID_OK:
			dlg.Destroy()

	def OnAboutElecSus(self,event):
		## do stuff
		dlg = wx.MessageDialog(self, NOTICE.noticestring, "About", wx.OK)
		if dlg.ShowModal() == wx.ID_OK:
			dlg.Destroy()
					
	def OnSaveCSVData(self,event):
		SaveFileDialog = wx.FileDialog(self,"Save Output File", "./", "Outputs",
			"CSV files (*.csv)|*.csv", wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
		
		if SaveFileDialog.ShowModal() == wx.ID_OK:
			output_filename = SaveFileDialog.GetPath()
			SaveFileDialog.Destroy()
			
			### Select which data traces to save
			#
			#
			#
			
			dlg = wx.MessageDialog(self, "Save ability not yet implemented...", "No no no", wx.OK)
			if dlg.ShowModal() == wx.ID_OK:
				dlg.Destroy()
			
			profile_filename = output_filename + ".csv"
			#check for overwrite current files
			if os.path.isfile(profile_filename):
				OverwriteDialog = wx.MessageDialog(self,"Warning: file exists already! Overwrite?",\
					"Overwrite?",wx.YES_NO|wx.NO_DEFAULT)
				
				if OverwriteDialog.ShowModal() == wx.NO:
					OverwriteDialog.Destroy()
					return # exit without saving
				else:
					OverwriteDialog.Destroy()
			
			SaveMessage = wx.MessageDialog(self, \
				"Files created", "Files created", wx.OK|wx.ICON_INFORMATION)
			SaveMessage.ShowModal()
			SaveMessage.Destroy()

	def OnSaveFig(self,event):
		#widcards for file type selection
		wilds = "PDF (*.pdf)|*.pdf|" \
				"PNG (*.png)|*.png|" \
				"EPS (*.eps)|*.eps|" \
				"All files (*.*)|*.*"
		exts = ['.pdf','.png','.eps','.pdf'] # default to pdf
		SaveFileDialog = wx.FileDialog(self,message="Save Figure", defaultDir="./", defaultFile="figure",
			wildcard=wilds, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
		SaveFileDialog.SetFilterIndex(0)
		
		if SaveFileDialog.ShowModal() == wx.ID_OK:
			output_filename = SaveFileDialog.GetPath()
			if output_filename[-4:] == exts[SaveFileDialog.GetFilterIndex()]:
				output_filename = output_filename[:-4]
			
			#save all current figures
			for fig, id in zip(self.figs, self.fig_IDs):
				fig.savefig(output_filename+'_'+str(id)+exts[SaveFileDialog.GetFilterIndex()])
				
		SaveFileDialog.Destroy()				
##
####			
####### Plot-specific actions
####
##
	def OnFileOpen(self,event):
		self.dirname= ''
		dlg_choice = wx.SingleChoiceDialog(self,"Choose type of data to be imported","Data import",choices=OutputTypes)
		
		# wait for OK to be clicked
		if dlg_choice.ShowModal() == wx.ID_OK:
			choice = dlg_choice.GetSelection()
			# use the choice index to select which axes the data appears on	
			choice_index = self.plot_output_indices.index(choice)
			
			#set experimental display on, and update menus
			self.display_expt_curves[choice_index] = True
			self.showEplotsSubMenu.GetMenuItems()[choice_index].Check(True)
			
			dlg_choice.Destroy()
		
			dlg_open = wx.FileDialog(self,"Choose 2-column csv file (Detuning, Transmission)",
									self.dirname,"","*.csv",wx.OPEN)
			
			# if OK button clicked, open and read file
			if dlg_open.ShowModal() == wx.ID_OK:
				self.filename = dlg_open.GetFilename()
				self.dirname = dlg_open.GetDirectory()
				#call read
				self.x_expt_arrays[choice_index],self.y_expt_arrays[choice_index] = read_CSV(os.path.join(self.dirname,self.filename),spacing=0)
				
				#overwrite fit_array data - i.e. last data to be loaded
				self.x_fit_array = self.x_expt_arrays[choice_index]
				self.y_fit_array = self.y_expt_arrays[choice_index]
				
				# implicit that the fit type is the same as last data imported
				self.fit_datatype = OutputTypes[choice_index]
				
				## create main plot				
				self.OnCreateAxes(self.figs[0],self.canvases[0],clear_current=True)
				
				## if a similar theory curve exists for this type of experimental data, create a residual plot as well - interpolate to match the two arrays
				#
				#
				#
				
			dlg_open.Destroy()

	def OnCreateAxes(self,fig,canvas,clear_current=True):
		""" (Re)Create as many sets of axes in the main figure as are needed, and label them all """

		# calculate how many axes should be displayed - count if any elements of display_expt/theory_curves are true
		n_axes = [i|j for i,j in \
			zip(self.display_theory_curves,self.display_expt_curves)].count(True)
		
		# clear figure and start again from nothing if number of axes has changed, or if hold is off
		if clear_current:
			fig.clf()
		if n_axes != len(fig.axes):
			fig.clf()
		
		# create bare axes, all equal sizes
		if n_axes == 0:
			n_axes = 1
		
		fig.add_subplot(n_axes,1,1)
		i=0
		for i in range(1,n_axes):
			fig.add_subplot(n_axes,1,i+1,sharex=fig.axes[0])
		
		
		# Testing:
		#print self.display_expt_curves
		#print self.y_expt_arrays
		
		## PLOT BASED ON BOOL DISPLAY_x_CURVE
		i=0
		for displayT,displayE,yt,xe,ye,ylabel in zip(self.display_theory_curves,self.display_expt_curves,\
								self.y_arrays,self.x_expt_arrays,self.y_expt_arrays,self.plot_ylabels):
			try:
				#print i
				ax = fig.axes[i]
			except:
				pass
			if displayE and ye is not None:
				if isinstance(ye, (list,tuple)):
					for xi,yi in zip(xe,ye):
						ax.plot(xi,yi)
				else:
					ax.plot(xe,ye)
			if displayT and yt is not None:
				if isinstance(yt, (list,tuple)):
					for yi in yt:
						ax.plot(self.x_array,yi)
				else:
					ax.plot(self.x_array,yt)
			if displayT or displayE:
				ax.set_ylabel(ylabel)
				i += 1
			
		# set x axis label and rescale all axes to fit data
		fig.axes[-1].set_xlabel('Detuning (GHz)')
		for ax in fig.axes:
			ax.autoscale_view(tight=True)
			
		#fig.axes[-1].set_xlim(self.xrange)
		
		# remove the rest of the x tick labels from all but the bottom panel
		for ax in fig.axes[:-1]:
			plt.setp(ax.get_xticklabels(),visible=False)
		
		#print 'Created '+str(n_axes)+' axes'
		
		# update the plot window
		self.draw_fig(fig,canvas)
	
	def OnCreateResidualPlot(self,fig,canvas):
		""" Create a single-panel plot with residuals and histogram of residuals, using
			subplot2grid in matplotlib """
		
		yy = 8
		xx = 8
		ax_main = fig.subplot2grid((yy,xx),(0,0),colspan=xx-1,rowspan=yy-1)
		ax_residual = fig.subplot2grid((yy,xx),(yy-1,0),colspan=xx-1,sharex=ax_main)
		ax_hist = fig.subplot2grid((yy,xx), (yy-1,xx-1), sharey=ax_residual)
		
		ax_main.set_ylabel()
		ax_residual.set_xlabel('Detuning (GHz)')
		ax_residual.set_ylabel('Residuals (%)')
		self.draw_fig(fig,canvas)
	
	def OnCreateColorMapPlot(self,event):
		""" Create the plot for colormap/contour data ... """
		# to be implemented...
		pass
	
	def OnCreateResultsPlot(self,event):
		""" Create a new plot for fit results """
		# to be implemented...
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
	
	def OnGridToggle(self,event):
		for ax in self.fig.axes:
			ax.grid(bool(event.Checked()))
		self.draw_fig()
					
	def OnShowTplots(self,event):
		for ii, item in enumerate(self.showTplotsSubMenu.GetMenuItems()):
			if item.IsChecked():
				self.display_theory_curves[ii] = True
			else:
				self.display_theory_curves[ii] = False
		
		#redraw plot
		self.OnCreateAxes(self.figs[0],self.canvases[0])

	def OnShowEplots(self,event):
		for ii, item in enumerate(self.showEplotsSubMenu.GetMenuItems()):
			if item.IsChecked():
				self.display_expt_curves[ii] = True
			else:
				self.display_expt_curves[ii] = False
		
		#redraw plot
		self.OnCreateAxes(self.figs[0],self.canvases[0])
				
	def draw_fig(self,fig,canvas):
		""" shortcut method for redrawing figure and rearranging subplots to fill space """
		try:
			fig.tight_layout()
		except:
			pass
		self.update()
	
	def update(self):
		for can in self.canvases:
			can.draw()
		self.SendSizeEvent()

	## Fitting specific actions
	def OnFitTypeChange(self,event):
		for item in self.fittypeSubMenu.GetMenuItems():
			if item.IsChecked():
				self.fit_algorithm = item.GetLabel()
				print self.fit_algorithm
		### THIS NEEDS TO TALK TO THE RADIO BUTTONS ON THE FIT PANEL, AND VICE-VERSA

	def OnAdvancedOptions(self,event):
		dlg = AdvancedFitOptions(self,"Advanced Fit Options",wx.ID_ANY)
		
		# Show() rather than ShowModal() - doesn't halt program flow
		if dlg.ShowModal() == wx.ID_OK:
			self.advanced_fitoptions = dlg.return_all_options()
			print self.advanced_fitoptions
				
##
#### Main purpose of this - call elecsus with various settings!
##
	def OnFitButton(self,event):
		""" Call elecsus to fit data """
		if self.y_fit_array == None:
			#warn about no data present
			dlg = wx.MessageDialog(self, "No experimental data has been loaded, cannot proceed with fitting...", "No no no", wx.OK)
			
			if dlg.ShowModal() == wx.ID_OK:
				pass
		## if number of booleans > 3 and ML fitting selected, warn about fitting methods
		#
		#
		
		else:
			self.Call_ElecSus('Fit')
		
	def OnComputeButton(self,event):
		""" Call elecsus to compute spectrum """
		self.Call_ElecSus('Theory')

	def Call_ElecSus(self,calc_or_fit):
		""" call elecsus with passed arguments - single calculation or fit, depending on option """
		
		#get parameter list, xrange from either theory or fit panel
		if calc_or_fit == 'Theory':
			panel = self.ThyPanel
			xmin,xmax,npts = [input.GetValue() for input in panel.DetuningCtrl]
			self.x_array = np.linspace(xmin,xmax,npts)
			xrange = self.x_array
			
		elif calc_or_fit == 'Fit':
			panel = self.FitPanel
			xrange = self.x_fit_array	
			
		
		# get parameter list
		const_params = panel.fixed_paramlist_inputs
		elem = const_params[0].GetStringSelection()
		dline = const_params[1].GetStringSelection()
		constrain = const_params[2].IsChecked()
		
		#print elem, dline, constrain
		
		self.fit_params = [input.GetValue() for input in panel.fit_paramlist_inputs]
		#print xrange
		#print fit_params

		
		## convert parameter order into order required by elecsus...
		## Could change elecsus' order to make it match, but this might cause
		## backwards compatibility issues. This is a bit clunky and not at all
		## elegant, but it will do for now!
		self.params = [0]*14
		self.re_order_list = [0,1,2,8,7,5,6,3,4]
		
		self.params[0:2] = elem, dline # Elem, Dline
		self.params[11] = constrain # constrain doppler/number density
		self.params[12] = self.fit_params[9]   # % of K40
		self.params[13] = self.fit_params[10]  # % of K41
		for i, el in enumerate(self.re_order_list):
			self.params[i+2] = self.fit_params[el] # B, T, lcell
		
		## do the same rearranging for fit_bools
		if calc_or_fit=='Fit':
			self.fit_bools = [checkbox.IsChecked() for checkbox in panel.fit_paramlist_bools]
			print self.fit_bools
			self.fit_bools_ordered = [0]*len(self.fit_bools)
			self.re_order_bool_list = [0,1,2,8,7,5,6,3,4,9,10]
			for i, el in enumerate(self.re_order_bool_list):
				self.fit_bools_ordered[i] = self.fit_bools[el]
		
		
		## Finally, actually call ElecSus!
		if calc_or_fit == 'Theory':
			print '\n\n'
			print time.ctime()
			print 'Calling ElecSus for single calculation with parameters:'
			print params
			
			spectrum_data = elecsus_new.calculate(xrange,params)
			self.y_arrays[0:4] = spectrum_data[0:4] #S0,S1,S2,S3
			self.y_arrays[4] = [spectrum_data[4],spectrum_data[5]] # Ix,Iy
			self.y_arrays[5] = [spectrum_data[9],spectrum_data[10]] # alpha+,alpha-
			self.y_arrays[6] = [spectrum_data[6],spectrum_data[7]] # n-,n+
			self.y_arrays[7] = spectrum_data[8] # Faraday rotation angle phi
			
			#print self.y_arrays
			#for i, array in enumerate(self.y_arrays):
			#	print i, type(array)
			self.OnCreateAxes(self.figs[0],self.canvases[0])
			
		elif calc_or_fit == 'Fit':
			print '\n\n'
			print time.ctime()
			print 'Calling ElecSus for fitting data with initial parameters:'
			print self.params
			print self.fit_bools_ordered
			
			## log time, data set to be fitted, initial parameters on the 'Fit Info' notebook:
			font1 = wx.Font(10, wx.TELETYPE, wx.NORMAL, wx.NORMAL)#, False, u'Consolas')
			self.FitInformation.StatusTextBox.SetFont(font1)
			self.FitInformation.write('Fit started at:'.ljust(25)+time.ctime()+'\n')
			self.FitInformation.write('Data set to be fitted:'.ljust(25)+os.path.join(self.dirname,self.filename)+'\n')
			self.FitInformation.write('Experimental Data Type:'.ljust(25)+self.fit_datatype+'\n')
			self.FitInformation.write('Initial parameters (Floated):\n')
			[ self.FitInformation.write('\t'+param_name.ljust(30)+str(param_value).ljust(10)+'('+str(boo)+')\n') for param_name, param_value,boo in zip(fittable_parameterlist, self.fit_params, self.fit_bools) ]
			self.FitInformation.write('Using algorithm:'+self.fit_algorithm)

			self.FitInformation.write('\n\nRunning Fit...')

			##
			#### Run Fitting algorithm
			##
			
			### Use another thread for running elecsus fitting so the main panel doesn't stop responding
			if not self.already_fitting:
				#initialise thread for fitting
				fitThread = FittingThread(self)
				
				# start thread
				fitThread.start()
				
				# only allow one fit to run at any time
				self.already_fitting = True
				
				## show dialog box with status bar?
				self.fitting_dlg = wx.ProgressDialog("Fit in progress","Fitting in progress. Please wait for fitting to complete. This may take a while.",maximum=100,style=wx.PD_ELAPSED_TIME)
				self.fitting_dlg.Show(True)
				fitProgressThread = ProgressThread(self)
				fitProgressThread.start()
				
			else:
				pass
				## dialog box for already fitting...
			


	def OnFitCompleted(self,event):
		""" Task to complete after fitting completed """
		self.FitInformation.write('  ...Fit Completed\n\n')
			
		## reverse the re-ordering we did before...
		self.opt_params_new = self.fit_params
		for i, el in enumerate(self.re_order_list):
			self.opt_params_new[el] = self.opt_params[i+2]
		
		#self.opt_params = self.opt_params_new
		
		## Add information to the Fitting text box
		self.FitInformation.write('Optimised parameters are (Fitted):\n')
		[ self.FitInformation.write('\t'+param_name.ljust(30)+str(param_value).ljust(10)+'('+str(boo)+')\n') for param_name, param_value,boo in zip(fittable_parameterlist, self.opt_params_new, self.fit_bools) ]
		#
		# Add chi-squared (if data has errorbars), and/or RMS numbers
		#
		
		
		# use fitted parameters to calculate new theory arrays
		print '\n\n'
		print time.ctime()
		print 'Calling ElecSus for single calculation with optimised parameters:'
		print self.opt_params_new
		
		spectrum_data = elecsus_new.calculate(self.x_fit_array,self.opt_params)
		self.x_array = self.x_fit_array
		self.y_arrays[0:4] = spectrum_data[0:4] #S0,S1,S2,S3
		self.y_arrays[4] = [spectrum_data[4],spectrum_data[5]] # Ix,Iy
		self.y_arrays[5] = [spectrum_data[9],spectrum_data[10]] # alpha+,alpha-
		self.y_arrays[6] = [spectrum_data[6],spectrum_data[7]] # n-,n+
		self.y_arrays[7] = spectrum_data[8] # Faraday rotation angle phi
		
		# update main plot
		# turn on theory curve for the type of plot we just fitted..
		index_fit_datatype = OutputTypes_index[OutputTypes.index(self.fit_datatype)]
		self.display_theory_curves[index_fit_datatype] = True
		
		# update the figure
		self.OnCreateAxes(self.figs[0],self.canvases[0])			
		print 'Updating main plot...'
		
		# create/update residuals plot
		#
		#
		
		# Reset the 'already fitting' flag
		self.already_fitting = False
		
		# close the 'fitting in progress' dialog box
		self.fitting_dlg.Show(False)
		
		
##
## Other global functions
##

def read_CSV(filename,spacing=0,columns=2):
	""" Read an n-column csv file. Return data as n numpy arrays """
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

### Run the thing...
if __name__ == '__main__':
	app = wx.App(redirect=False)
	frame = ElecSus_GUI_Frame(None,"ElecSus GUI")
				
	frame.Show()
	app.MainLoop()