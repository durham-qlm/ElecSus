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
rc('font',**{'family':'serif', 'size':14})
rc('lines', linewidth=2)


import wx
import os
import sys
import csv

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg as Toolbar

###IMPORT ELECSUS MODULES
import elecsus_new
sys.path.append('../libs/')
import spectra
import NOTICE


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
OutputTypes = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix, Iy', 'Alpha Plus/Minus', 'N Plus/Minus']

fixed_parameterlist = ['Element','D-line','Constrain Dopp./Density Temps']
element_list = ['Na', 'K', 'Rb', 'Cs']
D_line_list = ['D1', 'D2']

fittable_parameterlist = ['Bfield','Temperature','Cell length','Shift','Additional-Broadening','Theta-0','LC Polarisation','Doppler-Temperature','Rb-85','K40','K41']
units_parameterlist = [' [G]',' [C]',' [mm]',' [MHz]',' [MHz]',' [deg]',' [%]',' [C]',' [%]',' [%]',' [%]']
defaultvals_parameterlist = [0,20,75,0,0,0,50,20,72.17,1,0]

class OptionsPanel(wx.Panel):
	def __init__(self, parent, mainwin, paneltype):
		""" 
		Most of the options are the same for theory and fitting, but there are some differences. 
		Hence, the same panel class with an argument 'paneltype' which will set what is displayed.
		""" 
		# mainwin is the main panel so we can bind buttons to actions in the main frame """
		
		wx.Panel.__init__(self, parent)
		
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
			fit_ML = wx.RadioButton(self,label="Marquardt-Levenburg", style=wx.RB_GROUP)
			fit_RR = wx.RadioButton(self,label="Random-Restart")
			fit_SA = wx.RadioButton(self,label="Simulated Annealing")
			
			self.Bind(wx.EVT_RADIOBUTTON, mainwin.OnFitTypeChange, fit_ML)
			self.Bind(wx.EVT_RADIOBUTTON, mainwin.OnFitTypeChange, fit_RR)
			self.Bind(wx.EVT_RADIOBUTTON, mainwin.OnFitTypeChange, fit_SA)
		
			# Run Fit Button
			RunFitButton = wx.Button(self,wx.ID_ANY, 'Run Fit', size=(-1,1.5*BtnSize))
			self.Bind(wx.EVT_BUTTON, mainwin.OnFitButton, RunFitButton)
			
		# Theory only:
		elif self.paneltype == 'Theory':

			# Calculate spectrum
			ComputeButton = wx.Button(self,wx.ID_ANY, 'Compute Spectrum', size=(-1,1.5*BtnSize))
			self.Bind(wx.EVT_BUTTON, mainwin.OnComputeButton, ComputeButton)
			#dummy variable
			self.fit_paramlist_bools = [0]*len(fittable_parameterlist)
			
		# Common options to both experiment and theory
		ImportButton = wx.Button(self,wx.ID_OPEN,label="Import Data", size=(-1,1.5*BtnSize))
		self.Bind(wx.EVT_BUTTON,mainwin.OnFileOpen,ImportButton)
			
		fixed_paramlist_labels = [ wx.StaticText(self,wx.ID_ANY,fixed_param) for fixed_param in fixed_parameterlist ]
		self.fixed_paramlist_inputs = [ \
		wx.ComboBox(self,wx.ID_ANY,choices=element_list,style=wx.CB_READONLY,size=(110,-1)), \
		wx.ComboBox(self,wx.ID_ANY,choices=D_line_list, style=wx.CB_READONLY,size=(110,-1)), \
		wx.CheckBox(self, label="")]
		
		# create list of parameters - labels and spin-control boxes
		fit_paramlist_labels = [ wx.StaticText(self,wx.ID_ANY,fit_param+unit) for fit_param,unit in zip(fittable_parameterlist,units_parameterlist) ]
		self.fit_paramlist_inputs = [ \
		wx.SpinCtrl(self,value=str(defval)) for defval in defaultvals_parameterlist ]

		# Don't need to bind SpinCtrl/ComboBox inputs to actions - the values are 
		# read whenever 'Compute' or 'Fit' buttons are pressed

		
		# Layout panel in sizers
		
		# main sizer for the panel
		panel_sizer = wx.BoxSizer(wx.VERTICAL)
		
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		label_sizer = wx.BoxSizer(wx.HORIZONTAL)
		label_sizer.Add((20,-1),0,wx.EXPAND)
		elem_label = wx.StaticText(self,wx.ID_ANY,"Select element and D-line")
		font = wx.Font(12,wx.DEFAULT, wx.NORMAL,wx.NORMAL)
		elem_label.SetFont(font)
		label_sizer.Add(elem_label)
		label_sizer.Add((20,-1),1,wx.EXPAND)
		panel_sizer.Add(label_sizer)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		# Add common elements first:
		
		for label, input in zip(fixed_paramlist_labels,self.fixed_paramlist_inputs):
			hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
			hor_sizer.Add((40,-1),0,wx.EXPAND)
			hor_sizer.Add(label,0,wx.EXPAND)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(input,0,wx.EXPAND)
			hor_sizer.Add((20,-1),0,wx.EXPAND)
			
			panel_sizer.Add(hor_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,3),0,wx.EXPAND)
		
		# vertical space
		panel_sizer.Add((-1,20),0,wx.EXPAND)
		label_sizer = wx.BoxSizer(wx.HORIZONTAL)
		label_sizer.Add((20,-1),0,wx.EXPAND)
		if self.paneltype == 'Theory':
			parameter_label_text = 'Computation parameters'
		elif self.paneltype == 'Fit':
			parameter_label_text = 'Initial parameters for fit'
		paramlabeltextbox = wx.StaticText(self,wx.ID_ANY,parameter_label_text)
		font = wx.Font(12,wx.DEFAULT, wx.NORMAL,wx.NORMAL)
		paramlabeltextbox.SetFont(font)
		label_sizer.Add(paramlabeltextbox,0,wx.EXPAND)
		label_sizer.Add((20,-1),1,wx.EXPAND)
		if self.paneltype == 'Fit':
			floatlabel = wx.StaticText(self,wx.ID_ANY,'Float?')
			label_sizer.Add(floatlabel,0,wx.ALIGN_BOTTOM)
			label_sizer.Add((15,-1),0,wx.EXPAND)
		panel_sizer.Add(label_sizer,0,wx.EXPAND)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		
		for label,input,boolbox in zip(fit_paramlist_labels,self.fit_paramlist_inputs,self.fit_paramlist_bools):
			hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
			hor_sizer.Add((40,-1),0,wx.EXPAND)
			hor_sizer.Add(label,0,wx.EXPAND)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(input,0,wx.EXPAND)
			if self.paneltype == 'Fit':
				hor_sizer.Add((25,-1),0,wx.EXPAND)
				hor_sizer.Add(boolbox)
			hor_sizer.Add((20,-1),0,wx.EXPAND)
			panel_sizer.Add(hor_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,3),0,wx.EXPAND)
		
		if self.paneltype == 'Fit':
			# vertical space
			panel_sizer.Add((-1,10),0,wx.EXPAND)
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
		
		# vertical space
		panel_sizer.Add((-1,20),1,wx.EXPAND)
		
		if self.paneltype == 'Theory':
			hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(ImportButton,0,wx.EXPAND)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(ComputeButton,0,wx.EXPAND)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			panel_sizer.Add(hor_sizer,0,wx.EXPAND)
		elif self.paneltype == 'Fit':
			hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(ImportButton)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(RunFitButton)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			panel_sizer.Add(hor_sizer,0,wx.EXPAND)
			
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		
		self.SetSizer(panel_sizer)
		self.Layout()
			
		
		
class PlotPanel(wx.Panel):
	def __init__(self, parent, mainwin, ID):
		""" mainwin is the main panel so we can bind buttons to actions in the main frame """
		wx.Panel.__init__(self, parent)
		
		self.fig = plt.figure(facecolor=(240./255,240./255,240./255))
		
		# display some text in the middle of the window to begin with
		#self.fig.text(0.5,0.5,"Welcome to the ElecSus GUI\n\nSelect some options to begin",ha='center',va='center')
		
		# create the wx objects to hold the figure
		self.canvas = FigureCanvasWxAgg(self, wx.ID_ANY, self.fig)
		self.toolbar = Toolbar(self.canvas) #matplotlib toolbar (pan, zoom, save etc)
		
		# Create vertical sizer to hold figure and toolbar - dynamically expand with window size
		plot_sizer = wx.BoxSizer(wx.VERTICAL)
		plot_sizer.Add(self.canvas, 1, wx.LEFT|wx.RIGHT|wx.GROW,border=0)
		plot_sizer.Add(self.toolbar, 0, wx.LEFT|wx.RIGHT|wx.EXPAND,border=0)

		mainwin.figs.append(self.fig)
		mainwin.fig_IDs.append(ID) # use an ID number to keep track of figures
		mainwin.canvases.append(self.canvas)

		self.fig.text(0.5,0.5,'Import some data to get started...', ha='center',va='center')
		self.fig.hold(False)
		
		self.SetSizer(plot_sizer)
		self.Layout()
		
class ElecSus_GUI_Frame(wx.Frame):
	""" Main window """
	def __init__(self,parent,title):
		wx.Frame.__init__(self,None,title=title,size=(1400,800))
		
		#if the window is closed, exit
		self.Bind(wx.EVT_CLOSE,self.OnExit)

		self.panel = wx.Panel(self)
		
		self.panel.SetBackgroundColour(wx.Colour(240,240,240))
		
		self._init_default_values()
		
		self._init_plot_defaults()
	
		self._init_panels()
		
		self._init_menus()
		
		self.OnCreateAxes(self.figs[0],self.canvases[0])
		
		
	def _init_default_values(self):
		""" Initialise default values for various things ... """
		
		self.figs = []
		self.canvases = []
		self.fig_IDs = []
		
		self.hold = False
		
	def _init_plot_defaults(self):
		""" List of default values for all plots. Theory plots can be up to 7 panels, each of which can be turned on/off """
		
		self.plot_outputs = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix, Iy', 'Alpha Plus/Minus', 'N Plus/Minus']
		self.plot_ylabels = ['Transmission', 'S1', 'S2', 'S3', '$I_x$, $I_y$', r'$\alpha^{\pm}$', r'$n^\pm -1$']
		self.plot_output_indices = [OutputTypes.index(po) for po in self.plot_outputs]
		
		self.xrange = [-10,10] # detuning range, in GHz
		self.npoints = 10000 # number of detuning points to calculate
		
		# default data for plots
		self.x_array = np.linspace(self.xrange[0],self.xrange[1],self.npoints)
		self.y_arrays = [None, None, None, None, None, None, None]
		self.x_expt_arrays = [None, None, None, None, None, None, None]
		self.y_expt_arrays = [None, None, None, None, None, None, None]
		
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
		fit_ML = self.fittypeSubMenu.AppendRadioItem(id_ML, "&Marquardt-Levenburg", "Use ML Fitting")
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
		#initialise
		#
		
			
		#
		# ... other menus
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
		
	## Create plot part of the window
		
		## create plot in a notebook-style for adding more tabs later on
		PlotTabs = wx.Notebook(self.panel)
		
		T_Panel = PlotPanel(PlotTabs,self,'Theory')
		E_Panel = PlotPanel(PlotTabs,self,'Experiment')
		
		PlotTabs.AddPage(T_Panel, "Main Plot")
		PlotTabs.AddPage(E_Panel, "Fit / Residuals Plot")

		plot_sizer = wx.BoxSizer(wx.VERTICAL)
		plot_sizer.Add(PlotTabs,1, wx.EXPAND)
		
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
		self.ThyPanel = OptionsPanel(TabPanel,self,'Theory')
		self.FitPanel = OptionsPanel(TabPanel,self,'Fit')
		
		TabPanel.AddPage(self.ThyPanel, "Theory Settings")
		TabPanel.AddPage(self.FitPanel, "Fit Settings")
		
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
				
#		
## Actions for events
#
	
##
#### General Actions - Close window, Save Data etc...
##

	def OnExit(self,event):
		self.Destroy()
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
				
				self.OnCreateAxes(self.figs[0],self.canvases[0],clear_current=True)
				
			dlg_open.Destroy()

	def OnCreateAxes(self,fig,canvas,clear_current=True):
		""" (Re)Create as many sets of axes in the current figure as are needed, and label them all """

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
		
		## PLOT BASED ON BOOL DISPLAY_x_CURVE
		
		print self.display_expt_curves
		print self.y_expt_arrays
		
		i=0
		for displayT,displayE,yt,xe,ye,ylabel in zip(self.display_theory_curves,self.display_expt_curves,\
								self.y_arrays,self.x_expt_arrays,self.y_expt_arrays,self.plot_ylabels):
			try:
				#print i
				ax = fig.axes[i]
			except:
				pass
			if displayT and yt is not None:
				if type(yt)=='list':
					for yi in yt:
						ax.plot(self.x_array,yi)
				else:
					ax.plot(self.x_array,yt)
			if displayE and ye is not None:
				if type(ye)=='list':
					for xi,yi in zip(xe,ye):
						ax.plot(xi,yi)
				else:
					ax.plot(xe,ye)
			if displayT or displayE:
				ax.set_ylabel(ylabel)
				i += 1
			
		fig.axes[-1].set_xlabel('Detuning (GHz)')
		for ax in fig.axes:
			ax.autoscale_view(tight=True)
			
		#fig.axes[-1].set_xlim(self.xrange)
		
		for ax in fig.axes[:-1]:
			plt.setp(ax.get_xticklabels(),visible=False)
		
		#print 'Created '+str(n_axes)+' axes'
		self.draw_fig(fig,canvas)
	
	def OnCreateResidualPlot(self,fig,canvas):
		""" Create a single-panel plot with residuals and histogram of residuals, using
			subplot2grid in matplotlib """
		
		yy = 8
		xx = 8
		ax_main = fig.subplot2grid((yy,xx),(0,0),colspan=xx-1,rowspan=yy-1)
		ax_residual = fig.subplot2grid((yy,xx),(yy-1,0),colspan=xx-1,sharex=ax_main)
		ax_hist = fig.subplot2grid((yy,xx), (yy-1,xx-1), sharey=ax_residual)
		
		self.draw_fig(fig,canvas)
	
	def OnCreateColorMapPlot(self,event):
		""" Create the plot for colormap/contour data ... """
		pass
	
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
		canvas.draw()

	## Fitting specific actions
	def OnFitTypeChange(self,event):
		for item in self.fittypeSubMenu.GetMenuItems():
			if item.IsChecked():
				self.fit_type = item.GetLabel()
				#print self.fit_type
		### THIS NEEDS TO TALK TO THE RADIO BUTTONS ON THE FIT PANEL, AND VICE-VERSA

##
#### Main purpose of this - call elecsus with various settings!
##
	def OnFitButton(self,event):
		""" Call elecsus to fit data """
		self.Call_ElecSus('Fit')
		
	def OnComputeButton(self,event):
		""" Call elecsus to compute spectrum """
		self.Call_ElecSus('Theory')

	def Call_ElecSus(self,calc_or_fit):
		""" call elecsus with passed arguments - single calculation or fit, depending on option """
		
		#get parameter list from either theory or fit panel
		if calc_or_fit == 'Theory':
			panel = self.ThyPanel
		elif calc_or_fit == 'Fit':
			panel = self.FitPanel
		
		# get parameter list
		const_params = panel.fixed_paramlist_inputs
		elem = const_params[0].GetCurrentSelection()
		dline = const_params[1].GetCurrentSelection()
		constrain = const_params[2].IsChecked()
		fit_params = [input.GetValue() for input in panel.fit_paramlist_inputs]
		print fit_params
		if calc_or_fit=='Fit':
			fit_bools = [checkbox.IsChecked() for checkbox in panel.fit_paramlist_bools]
			print fit_bools
		
		## call elecsus here...
		#
		#
		#
		#
		#
		
		
		
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
	#frame.Maximize()
	frame.Show()
app.MainLoop()