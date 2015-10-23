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

ElecSus GUI

v1.0.1 (2015-10-23)
	-- minor bug fix where the plot selection popups would not display the Phi plots
v1.0.0 (2015-09-03)


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
pylab.ioff() # set interactive mode off

import os
import sys
import csv
import time
import cPickle as pickle

import numpy as np
import matplotlib.pylab as plt

# use relative file paths
elecsus_dir = os.path.dirname(__file__)

try:
	import wx
except ImportError:
	print "wxPython 2.8 needs to be installed for this program to work! \n\
It is not currently possible to install this automatically through pip/easy_install.\n"
	if os.name == 'posix':
		print "For Ubuntu/Debian, wxPython is not supported in Enthought Canopy.\n\
Instead, use the system python distribution (/usr/bin/python) and install through apt:\n\n\
>    (sudo) apt-get install python-wxgtk2.8 python-wxtools wx2.8-i18n libwxgtk2.8-dev libgtk2.0-dev"
	else:
		print 'For Windows, recommended install is using Enthought Canopy'
	raise ImportError
	
wx.Log.SetLogLevel(0) # Ignore warnings
import wx.lib.scrolledpanel as scrolled
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN
from wx.lib.popupctl import PopupControl
# Matplotlib/wx integration
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg as Toolbar

## fits run in a separate thread so the GUI isn't blocked
import threading
# define a new event type which will be triggered when fitting is completed
import wx.lib.newevent
FitCompleteEvent, EVT_FIT_COMPLETE = wx.lib.newevent.NewEvent()

# import elecsus modules
import elecsus_methods as elecsus
from libs import spectra
from libs import NOTICE
from libs import data_proc
from libs.durhamcolours import *
from libs.durhamcolours import cols as durhamcols

#replace default matplotlib text and color sequence with durham colours
from matplotlib import rc
rc('text', usetex=False)
rc('font',**{'family':'serif', 'size':14})
rc('lines', linewidth=2)
rc('axes', color_cycle=durhamcols)

# preamble.py includes tooltip text, default values, labels...
from libs.preamble import *

class PlotSelectionPopUp(wx.PopupTransientWindow):
	""" Popup box to handle which plots are displayed. """
	def __init__(self,parent,style,mainwin,plottype):
		wx.PopupTransientWindow.__init__(self,parent,style)

		self.mainwin = mainwin
		self.plottype = plottype

		self.win = wx.Panel(self)#,wx.ID_ANY,pos=(0,0),size=(180,200),style=0)

		self.selection = wx.CheckListBox(self.win, wx.ID_ANY, choices = OutputPlotTypes, size=(150,-1))#,pos=(0,0))
		self.win.Bind(wx.EVT_CHECKLISTBOX, self.OnTicked, self.selection)
		
		if plottype == 'Theory':
			display_curves = self.mainwin.display_theory_curves
		else:
			display_curves = self.mainwin.display_expt_curves

		checked_items = []
		for i in range(len(display_curves)):
			if display_curves[i]:
				checked_items.append(i)

		self.selection.SetChecked(checked_items)
		#self.okbtn = wx.Button(self.win,wx.ID_OK,size=(120,BtnSize))
		#self.Bind(wx.EVT_BUTTON, self.OnOK, self.okbtn)

		self.SetSize(self.selection.GetSize()+(10,10))
		
		popup_sizer = wx.BoxSizer(wx.VERTICAL)
		popup_sizer.Add(self.selection,0,wx.EXPAND)
		#popup_sizer.Add((-1,5),1,wx.EXPAND)
		#popup_sizer.Add(self.okbtn,0,wx.EXPAND)

		#sz = popup_sizer.GetBestSize()
		#self.win.SetSize((sz.width+20, sz.height+20))
		win_sizer = wx.BoxSizer(wx.HORIZONTAL)
		win_sizer.Add(popup_sizer,0,wx.EXPAND|wx.ALL, border=5)

		self.win.SetSizer(win_sizer)
		self.win.Fit()

		wx.CallAfter(self.Refresh)

	def OnDismiss(self):
		""" 
		Action to perform when the popup loses focus and is closed. 
		
		Gets the tick box values, updates the main plot and changes 
		the menu items to match the popup box
		"""
		# re=plot the figure
		#self.OnTicked(1)

		#self.mainwin.Refresh()

		self.Dismiss()
	
	def OnTicked(self,event):
		""" 
		Action to perform when tick boxes are ticked 
		
		Gets the tick box values, updates the main plot and changes 
		the menu items to match the popup box
		"""
		if self.plottype == 'Theory':
			items = self.selection.GetChecked()
			self.mainwin.display_theory_curves = [False]*9
			for menuitem in self.mainwin.showTplotsSubMenu.GetMenuItems():
				menuitem.Check(False)
			for item in items:
				self.mainwin.display_theory_curves[item] = True
				self.mainwin.showTplotsSubMenu.GetMenuItems()[item].Check(True)
			print self.mainwin.display_theory_curves
		elif self.plottype == 'Fit':
			items = self.selection.GetChecked()
			self.mainwin.display_expt_curves = [False]*9
			for menuitem in self.mainwin.showEplotsSubMenu.GetMenuItems():
				menuitem.Check(False)
			for item in items:
				self.mainwin.display_expt_curves[item] = True
				self.mainwin.showEplotsSubMenu.GetMenuItems()[item].Check(True)
			print self.mainwin.display_expt_curves

		self.mainwin.OnCreateAxes(self.mainwin.figs[0],self.mainwin.canvases[0])
		
class ProgressBarFrame(wx.Frame):
	""" Not actually used at the moment - a custom alternative to the progressdialog in wx """
	def __init__(self, parent, title, pbrange = 100) :
		wx.Frame.__init__(self, parent = parent, title = title)
		self.range = pbrange
		self.createProgressbar()
		self.SetMinSize((400, 10))
		self.Centre()
		self.Show()
		self.t0 = time.time()
		self.elapsed_time_timer.Start(1000)

	def createProgressbar(self):
		""" Create a progress bar using wx.Gauge element """
		self.pb = wx.Gauge(self)
		self.pb.SetRange(range = self.range)

		self.elapsed_time_st  = wx.StaticText(self, label = 'Elapsed Time:')
		self.elapsed_time_val = wx.StaticText(self, label = '00:00:00')

		vbox_main = wx.BoxSizer(wx.VERTICAL)
		hbox_time = wx.BoxSizer(wx.HORIZONTAL)
		hbox_time.Add(self.elapsed_time_st,  0, wx.ALIGN_LEFT | wx.EXPAND | wx.ALL, 5)
		hbox_time.Add(self.elapsed_time_val, 0, wx.ALIGN_LEFT | wx.EXPAND | wx.ALL, 5)
		vbox_main.Add(self.pb,   0, wx.EXPAND | wx.ALL, 5)
		vbox_main.Add(hbox_time, 0, wx.EXPAND | wx.ALL, 5)

		self.SetSizerAndFit(vbox_main)

		self.elapsed_time_timer = wx.Timer(self)
		self.Bind(wx.EVT_TIMER, self.onTickTimer, self.elapsed_time_timer)

	def onTickTimer(self, event):
		""" update the time """
		fmt='%H:%M:%S'
		self.elapsed_time_val.SetLabel(time.strftime(fmt, time.gmtime(time.time()-self.t0)))

class AbortError(Exception):
	""" Custom exception for handling aborting fitting process -- not implemented yet."""
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class FittingThread(threading.Thread):
	""" 
	Fitting takes a long time, so we run the fit in another thread, leaving the GUI active in the meantime.
	When the fitting completes, this thread triggers an event in the GUI that then completes the rest of the
	fitting routine - i.e. updates the plot, writes text ... etc
	
	child class of the main threading.Thread class
	"""
	def __init__(self,parent):
		threading.Thread.__init__(self)
		
		self.mainwin = parent
		#self.start() ## this calls the run() method - called from the GUI
	
	def run(self):
		""" Run fitting thread - called by <>.start() method """
		
		# mainwin is the top-level window (the ElecSus_GUI_Frame instance)
		mainwin = self.mainwin
		
		
		## crop x and y arrays sent to fit_data() method, based on fit_bounds (if it's not [None, None])
		if None in mainwin.fit_bounds:
			# use full arrays if fit bounds not specified (default)
			x_array, y_array = mainwin.x_fit_array, mainwin.y_fit_array
		else:
			## crop data to specified range
			
			#create aliases for ease
			fb = mainwin.fit_bounds
			xfa = mainwin.x_fit_array
			yfa = mainwin.y_fit_array
			#
			y_array = yfa[(xfa>fb[0]) & (xfa<fb[1])]
			x_array = xfa[(xfa>fb[0]) & (xfa<fb[1])]
			
			fb, xfa, yfa = None, None, None
		
		print 'Fitting data in the detuning range (GHz):  ',x_array[0], ' to ',x_array[-1]
		mainwin.FitInformation.write('Fitting in the detuning range (GHz):  '+str(x_array[0])+' to '+str(x_array[-1]))
		mainwin.FitInformation.write('\n\n')

		## stick this in a try/except to have an abort control?
		#try:
		mainwin.opt_params, mainwin.rms = elecsus.fit_data((x_array,y_array),\
								mainwin.params,mainwin.fit_bools_ordered,\
								mainwin.fit_datatype,mainwin.fit_algorithm,\
								**mainwin.advanced_fitoptions)
		#print 'Fit completed'
			
			
		'''except RuntimeError: 
			## specify which kind of exception! - maybe add custom exception to cancel fitting?
			#self.fitting_dlg.Destroy()
			
			dlg = wx.MessageDialog(mainwin,"Fitting encountered an error! Check the errors panel for details.","Fitting Error",style=wx.OK|wx.CENTRE|wx.ICON_ERROR)
			
			if dlg.ShowModal() == wx.ID_OK:
				pass
			raise RuntimeError

		'''
		
		#post an event to the main panel to run the OnFitComplete() method
		evt = FitCompleteEvent()
		wx.PostEvent(mainwin, evt)
		## alternately, could use wx.CallAfter here ...		

class ProgressThread(threading.Thread):
	""" update the progress bar continually to show user something is happening ..."""
	def __init__(self,parent):
		""" Create the thread """
		threading.Thread.__init__(self)
		self.mainwin = parent
		self.pb = parent.fitting_dlg #.pb
		
	def run(self):
		""" Action to run when .start() method is called """
		while self.mainwin.already_fitting:
			for i in range(100):
				if not self.mainwin.already_fitting: break
				#print i,' ',
				#wx.CallAfter(self.pb.SetValue,i)
				wx.CallAfter(self.pb.Update,i)
				time.sleep(0.1)
		print 'Quitting progress bar update'
	
class OptionsPanel(scrolled.ScrolledPanel):
	def __init__(self, parent, mainwin, paneltype):
		""" 
		Panel which holds most of the control elements in this program. 
		Most of the options are the same for theory and fitting, but there are some differences.
		Hence, the same panel class with an argument 'paneltype' which will set what is displayed.
		
		This ScrolledPanel class will automatically add horizontal and vertical scroll bars when required, and dynamically turns them on and off as window is resized. Neat!
		
		##
		Note - this class has got a bit messy in it's old age, and could do with a rewrite to make it more clear!
		Another one for the To-Do list...
		##
		
		# 'mainwin' is the main application frame so we can bind buttons to actions in the main frame
		# [ it's the parent (frame) of the parent (tab notebook) of the options panel (self) ]
		
		"""
		
		scrolled.ScrolledPanel.__init__(self, parent)
		
		self.paneltype = paneltype
		
		# Fitting only:
		if self.paneltype == 'Fit':
			# Import data from csv
			
			# Booleans (check boxes) for fitting
			self.fit_paramlist_bools = \
				[ wx.CheckBox(self, label='') for parameter in fittable_parameterlist ]
						
			# Fitting algorithm selection
			# ties in with menu items
			self.fit_types = [wx.RadioButton(self,label="Marquardt-Levenberg", style=wx.RB_GROUP), \
						 wx.RadioButton(self,label="Random-Restart"),\
						 wx.RadioButton(self,label="Simulated Annealing") ]
			
			for fit_type, tt in zip(self.fit_types,fit_algorithm_tooltips):
				self.Bind(wx.EVT_RADIOBUTTON, mainwin.OnFitTypeChangePanel, fit_type)
				fit_type.SetToolTip(wx.ToolTip(tt))
			
			# Run Fit Button
			RunFitButton = wx.Button(self,wx.ID_ANY, 'Run Fit', size=(140,1.5*BtnSize))
			self.Bind(wx.EVT_BUTTON, mainwin.OnFitButton, RunFitButton)
			
		# Theory only:
		elif self.paneltype == 'Theory':
			# Detuning range selection
			Detuning_Labels = [ wx.StaticText(self,wx.ID_ANY,lab) for lab in ["Start [GHz]", "Stop [GHz]", "No. of points"]]
			self.DetuningCtrl = [ FloatSpin(self,value=str(defval),increment=definc,size=(80,-1),digits=2) for defval,definc in zip(detuning_defaults,detuning_increments) ]
			for ctrl in Detuning_Labels:
				ctrl.SetToolTip(wx.ToolTip(parameter_adjust_tooltip))
			
			# Calculate spectrum
			ComputeButton = wx.Button(self,wx.ID_ANY, 'Compute Spectrum', size=(140,1.5*BtnSize))
			self.Bind(wx.EVT_BUTTON, mainwin.OnComputeButton, ComputeButton)

			#dummy variable
			self.fit_paramlist_bools = [0]*len(fittable_parameterlist)
			
		# Common options to both experiment and theory
		ImportButton = wx.Button(self,wx.ID_OPEN,label="Import Data", size=(140,1.5*BtnSize))
		self.Bind(wx.EVT_BUTTON,mainwin.OnFileOpen,ImportButton)
		
		# Add tooltips describing each parameter
		fixed_paramlist_labels = [ wx.StaticText(self,wx.ID_ANY,fixed_param) for fixed_param in fixed_parameterlist ]
		for label, tt in zip(fixed_paramlist_labels,fixed_parameter_tooltips):
			label.SetToolTip(wx.ToolTip(tt))
		
		self.fixed_paramlist_inputs = [ \
		wx.ComboBox(self,wx.ID_ANY,choices=element_list,style=wx.CB_READONLY,size=(80,-1)), \
		wx.ComboBox(self,wx.ID_ANY,choices=D_line_list, style=wx.CB_READONLY,size=(80,-1)), \
		wx.CheckBox(self, label="")]
		self.fixed_paramlist_inputs[0].SetSelection(0)
		self.fixed_paramlist_inputs[1].SetSelection(0)
		self.fixed_paramlist_inputs[2].SetValue(True)
		
		# create list of parameters - labels and spin-control boxes
		fit_paramlist_labels = [ wx.StaticText(self,wx.ID_ANY,fit_param+unit) for fit_param,unit in zip(fittable_parameterlist,units_parameterlist) ]
		for label, tt in zip(fit_paramlist_labels, fit_parameter_tooltips):
			label.SetToolTip(wx.ToolTip(tt+parameter_adjust_tooltip))
		
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
			hor_sizer.Add(label,0,wx.ALIGN_CENTER_VERTICAL)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(input,0,wx.ALIGN_CENTER_VERTICAL)
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
			hor_sizer.Add(label,0,wx.ALIGN_CENTER_VERTICAL)
			hor_sizer.Add((20,-1),1,wx.EXPAND)
			hor_sizer.Add(input,0,wx.ALIGN_CENTER_VERTICAL)
			if self.paneltype == 'Fit':
				hor_sizer.Add((20,-1),0,wx.EXPAND)
				hor_sizer.Add(boolbox,0,wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
			hor_sizer.Add((10,-1),0,wx.EXPAND)
			panel_sizer.Add(hor_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,3),0,wx.EXPAND)
		
		if self.paneltype == 'Fit':
			# vertical space
			panel_sizer.Add((-1,5),0,wx.EXPAND)
			label_sizer = wx.BoxSizer(wx.HORIZONTAL)
			label_sizer.Add((10,-1),0,wx.EXPAND)
			parameter_label_text = 'Fit Algorithm'
			paramlabeltextbox = wx.StaticText(self,wx.ID_ANY,parameter_label_text)
			font = wx.Font(12,wx.DEFAULT, wx.NORMAL,wx.NORMAL)
			paramlabeltextbox.SetFont(font)
			label_sizer.Add(paramlabeltextbox,0,wx.EXPAND)
			label_sizer.Add((20,-1),1,wx.EXPAND)
			panel_sizer.Add(label_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,10),0,wx.EXPAND)
		
			for radiobtn in self.fit_types:
				hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
				hor_sizer.Add((30,-1),0,wx.EXPAND)
				hor_sizer.Add(radiobtn,0,wx.EXPAND)
				hor_sizer.Add((20,-1),1,wx.EXPAND)
				panel_sizer.Add(hor_sizer,0,wx.EXPAND)
				panel_sizer.Add((-1,4),0,wx.EXPAND)		
		elif self.paneltype == 'Theory':
			# vertical space
			panel_sizer.Add((-1,5),0,wx.EXPAND)
			label_sizer = wx.BoxSizer(wx.HORIZONTAL)
			label_sizer.Add((10,-1),0,wx.EXPAND)
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
				hor_sizer.Add((30,-1),0,wx.EXPAND)
				hor_sizer.Add(label,0,wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
				hor_sizer.Add((20,-1),1,wx.EXPAND)
				hor_sizer.Add(input,0,wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
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
	""" Panel to hold the matplotlib figure/canvas and toolbar """
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
		self.fig.text(0.5,0.5,'ElecSus GUI\n\nVersion 1.0.1\n\nTo get started, use the panel on the right\nto either Compute a spectrum or Import some data...', ha='center',va='center')
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
		""" 
		Save text file dump of all data in the panel - select filename, 
		check for overwrite and then pass to save_data() method
		"""
		SaveFileDialog = wx.FileDialog(self,"Save Output File", "./", "Outputs",
			"Text files (*.txt)|*.txt", wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
		
		if SaveFileDialog.ShowModal() == wx.ID_OK:
			
			output_filename = SaveFileDialog.GetPath()
			print output_filename
			#if output_filename[-4:] == exts[SaveFileDialog.GetFilterIndex()]:
			#	output_filename = output_filename[:-4]
			SaveFileDialog.Destroy()

			#check for overwrite current files
			if os.path.isfile(output_filename):
				OverwriteDialog = wx.MessageDialog(self,"Warning: file exists already! Overwrite?",\
					"Overwrite?",wx.YES_NO|wx.NO_DEFAULT)
				
				if OverwriteDialog.ShowModal() == wx.NO:
					OverwriteDialog.Destroy()
					return # exit without saving
				else:
					OverwriteDialog.Destroy()
			
			# do save
			self.save_data(output_filename)
		
	def save_data(self,filename):
		""" Save the data using wx.TextCtrl built-in method """
		success = self.StatusTextBox.SaveFile(filename) # returns true if no errors
		if not success:
			problem_dlg = wx.MessageDialog(self, "There was an error saving the data...", "Error saving", wx.OK|wx.ICON_ERROR)
			problem_dlg.ShowModal()

	def write(self,textstring):
		""" Append text to the text box """
		self.StatusTextBox.AppendText(textstring)
		
class DataProcessingDlg(wx.Dialog):
	def __init__(self,parent,title,id):
		""" Dialog box for smoothing and/or binning experimental data """
		
		wx.Dialog.__init__(self,parent,id,title,size=(360,450))
		
		self.parent = parent
		
		# make copy of non-binned / smoothed data in case reset button is pressed
		self.original_xdata = parent.x_fit_array
		self.original_ydata = parent.y_expt_arrays		
		
		panel_blurb = wx.StaticText(self,wx.ID_ANY,"Data smoothing and binning options.",size=(350,-1),style=wx.ALIGN_CENTRE_HORIZONTAL)
		panel_blurb.Wrap(350)
		
		bin_blurb = wx.StaticText(self,wx.ID_ANY,"Binning data is useful where the initial data size is very large. The output from this function is smaller by a specified factor n, which also removes some noise since every n data points are averaged into one.\n\n For fitting data, the number of experimental data points is a factor in how long the fit takes to run - making the array smaller makes the fit quicker.",size=(350,-1),style=wx.ALIGN_CENTRE_HORIZONTAL)
		bin_blurb.Wrap(350)
		
		self.cur_dsize_label = wx.StaticText(self,wx.ID_ANY,"Current data size:")
		self.cur_dsize = wx.TextCtrl(self,wx.ID_ANY,"",style=wx.TE_READONLY|wx.TE_CENTRE,size=(80,-1))
		self.cur_dsize.ChangeValue(str(len(parent.x_fit_array)))
		
		self.bin_size_label = wx.StaticText(self,wx.ID_ANY,"Bin size:")
		self.bin_size = wx.TextCtrl(self,wx.ID_ANY,"11",style=wx.TE_CENTRE,size=(80,-1))
		self.Bind(wx.EVT_TEXT, self.OnBinSizeChange, self.bin_size)
		
		self.new_dsize_label = wx.StaticText(self,wx.ID_ANY,"Data size after binning:")
		self.new_dsize = wx.TextCtrl(self,wx.ID_ANY,"",style=wx.TE_READONLY|wx.TE_CENTRE,size=(80,-1))
		self.new_dsize.ChangeValue(str(int(len(parent.x_fit_array)/float(self.bin_size.GetValue()))))
		
		bin_btn = wx.Button(self,label="Bin Data")
		self.Bind(wx.EVT_BUTTON,self.OnDoBin,bin_btn)
		
		smooth_blurb = wx.StaticText(self,wx.ID_ANY,"Smoothing is based on a moving average with a triangular weighting function of a specified width. For convenience, the output array is the same length as the input, but the smoothing is only applied over N - 2n data points (removing n from each end of the array), where N is the length of the data array and n is the width of the smoothing window.",size=(350,-1),style=wx.ALIGN_CENTRE_HORIZONTAL)
		smooth_blurb.Wrap(350)
		
		self.smooth_size_label = wx.StaticText(self,wx.ID_ANY,"Moving average window size:")
		self.smooth_size = wx.TextCtrl(self,wx.ID_ANY,"11",style=wx.TE_CENTRE,size=(80,-1))
				
		smooth_btn = wx.Button(self,label="Smooth Data")
		self.Bind(wx.EVT_BUTTON,self.OnDoSmooth,smooth_btn)
		
		reset_btn = wx.Button(self,label='Reset to original data')
		self.Bind(wx.EVT_BUTTON,self.OnReset,reset_btn)
		
		close_btn = wx.Button(self,wx.ID_OK,label='Close')
		
		## layout dialog box
		panel_sizer = wx.BoxSizer(wx.VERTICAL)
		
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		panel_sizer.Add(panel_blurb,0,wx.ALIGN_CENTRE_HORIZONTAL)
		
		bin_sizer = wx.BoxSizer(wx.HORIZONTAL)
		bin_sizer.Add((20,-1),0,wx.EXPAND)
		bin_sizer.Add(self.bin_size_label,0,wx.EXPAND)
		bin_sizer.Add((10,-1),1,wx.EXPAND)
		bin_sizer.Add(self.bin_size,0,wx.EXPAND)
		bin_sizer.Add((10,-1),0,wx.EXPAND)
		bin_sizer.Add(bin_btn,0,wx.EXPAND)
		bin_sizer.Add((10,-1),0,wx.EXPAND)

		cdata_sizer = wx.BoxSizer(wx.HORIZONTAL)
		cdata_sizer.Add((20,-1),0,wx.EXPAND)
		cdata_sizer.Add(self.cur_dsize_label,0,wx.EXPAND)
		cdata_sizer.Add((10,-1),1,wx.EXPAND)
		cdata_sizer.Add(self.cur_dsize,0,wx.EXPAND)
		cdata_sizer.Add((95,-1),0,wx.EXPAND)
		
		ndata_sizer = wx.BoxSizer(wx.HORIZONTAL)
		ndata_sizer.Add((20,-1),0,wx.EXPAND)
		ndata_sizer.Add(self.new_dsize_label,0,wx.EXPAND)
		ndata_sizer.Add((10,-1),1,wx.EXPAND)
		ndata_sizer.Add(self.new_dsize,0,wx.EXPAND)
		ndata_sizer.Add((95,-1),0,wx.EXPAND)
		
		smooth_sizer = wx.BoxSizer(wx.HORIZONTAL)
		smooth_sizer.Add((20,-1),0,wx.EXPAND)
		smooth_sizer.Add(self.smooth_size_label,0,wx.EXPAND)
		smooth_sizer.Add((10,-1),1,wx.EXPAND)
		smooth_sizer.Add(self.smooth_size,0,wx.EXPAND)
		smooth_sizer.Add((10,-1),0,wx.EXPAND)
		smooth_sizer.Add(smooth_btn,0,wx.EXPAND)
		smooth_sizer.Add((10,-1),0,wx.EXPAND)
		
		btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
		btn_sizer.Add((20,-1),1,wx.EXPAND)
		btn_sizer.Add(reset_btn,0,wx.EXPAND)
		btn_sizer.Add((20,-1),0,wx.EXPAND)
		btn_sizer.Add(close_btn,0,wx.EXPAND)
		btn_sizer.Add((20,-1),1,wx.EXPAND)
		
		
		panel_sizer.Add((-1,20),0,wx.EXPAND)		
		panel_sizer.Add(bin_blurb,0,wx.EXPAND)
		panel_sizer.Add((-1,20),0,wx.EXPAND)
		panel_sizer.Add(cdata_sizer,0,wx.EXPAND)
		panel_sizer.Add(bin_sizer,0,wx.EXPAND)
		panel_sizer.Add(ndata_sizer,0,wx.EXPAND)
		panel_sizer.Add((-1,20),0,wx.EXPAND)
		panel_sizer.Add(wx.StaticLine(self,-1,size=(-1,1),style=wx.LI_HORIZONTAL),
							0,wx.EXPAND|wx.LEFT|wx.RIGHT,border=20)
		panel_sizer.Add((-1,20),0,wx.EXPAND)	
		panel_sizer.Add(smooth_blurb,0,wx.EXPAND)
		panel_sizer.Add((-1,20),0,wx.EXPAND)		
		panel_sizer.Add(smooth_sizer,0,wx.EXPAND)
		panel_sizer.Add((-1,20),1,wx.EXPAND)		
		panel_sizer.Add(btn_sizer,0,wx.EXPAND)
		panel_sizer.Add((-1,10),0,wx.EXPAND)		
		
		self.SetSizer(panel_sizer)
		self.Layout()
		
		
	def OnBinSizeChange(self,event):
		""" update ui elements for value of new data length """
		#try:
		self.cur_dsize.ChangeValue(str(int(len(self.parent.x_fit_array))))
		self.new_dsize.ChangeValue(str(int(len(self.parent.x_fit_array)/float(self.bin_size.GetValue()))))
		#except:
		#	pass
		
	def OnDoBin(self,event):
		""" Action when bin button is pressed """
		bin_amnt = int(self.bin_size.GetValue())
		
		## only bin/smooth the most recently loaded data set
		self.parent.x_fit_array,self.parent.y_fit_array,ye = data_proc.bin_data(self.parent.x_fit_array,self.parent.y_fit_array,bin_amnt)
		
		#re-plot data
		self.update_plot()
		
		#update ui elements with new data size...
		self.OnBinSizeChange(1)
		
		
	def OnDoSmooth(self,event):
		""" Action when smooth button is pressed """
		smth_amnt = int(self.smooth_size.GetValue())
		
		self.parent.y_fit_array = data_proc.smooth_data(self.parent.y_fit_array,smth_amnt)
		
		#re-plot data
		self.update_plot()

	def OnReset(self,event):
		""" Reset to the original data, i.e. the data present when the dialog box was created """
		self.parent.x_fit_array = self.original_xdata
		self.parent.y_fit_array = self.original_ydata
		
		#re-plot data
		self.update_plot()
				
	def update_plot(self):
		""" Common method for updating the plot """
		#update plot arrays
		cindex = self.parent.choice_index
		self.parent.x_expt_arrays[cindex] = self.parent.x_fit_array
		self.parent.y_expt_arrays[cindex] = self.parent.y_fit_array

		self.parent.OnCreateAxes(self.parent.figs[0],self.parent.canvases[0])
						
class AdvancedFitOptions(wx.Dialog):
	def __init__(self,parent,title,id):
		""" Dialog box for selecting advanced fit options which are passed through to curve_fit / leastsq optimisation routines... """
		
		wx.Dialog.__init__(self,parent,id,title,size=(300,300))
				
		# main sizer for this dialog box
		panel_sizer = wx.BoxSizer(wx.VERTICAL)
		
		panel_blurb = wx.StaticText(self,wx.ID_ANY,"These are the options given to scipy.curve_fit (leastsq) when fitting data. \n\nMost of the time these will not need to be \nchanged from their defaults. However, if fitting \nis returning errors then consult the scipy \ndocumentation and adjust these parameters accordingly.\n ",size=(260,-1),style=wx.ALIGN_CENTRE_HORIZONTAL)
		panel_blurb.Wrap(260)
		
		self.fitopt_labels = ['Max Iterations (maxfev)', 'Max Tolerance [1e-8] (ftol)']
		self.fitopt_argnames = ['maxfev', 'ftol']
		fitopt_defaults = [1000, 15]
		fitopt_increments = [100, 1]
		self.fitopt_scaling = [1, 1e-9]
		fitopt_texts = [ wx.StaticText(self,wx.ID_ANY,label) for label in self.fitopt_labels ]
		self.fitopt_ctrl = [ wx.SpinCtrl(self,value=str(defval),size=(80,-1),min=0,max=10000,initial=defval) for defval,definc in zip(fitopt_defaults, fitopt_increments) ]
		
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
		""" get list of all ctrl elements """
		opt_dict = dict([(label,input.GetValue()*scaling) for label,input,scaling \
			in zip(self.fitopt_argnames,self.fitopt_ctrl,self.fitopt_scaling)])
		#print 'Advanced fit options dictionary: ',opt_dict
		return opt_dict

class FitBoundsDialog(wx.Dialog):
	def __init__(self,parent,title,id):
		""" Dialog box for selecting advanced fit options which are passed through to curve_fit / leastsq optimisation routines... """
		
		wx.Dialog.__init__(self,parent,id,title,size=(300,300))
		
		self.parent = parent

		# main sizer for this dialog box
		panel_sizer = wx.BoxSizer(wx.VERTICAL)
		
		panel_blurb = wx.StaticText(self,wx.ID_ANY,
			"Set bounds to fit data. By default, all data is fitted, but here the fitted data can be limited to a certain detuning range",
			size=(260,-1),style=wx.ALIGN_CENTRE_HORIZONTAL)
		panel_blurb.Wrap(260)
		
		self.fit_bounds = [wx.RadioButton(self,label="Use Full Data Range", style=wx.RB_GROUP), \
						 wx.RadioButton(self,label="Cropped Data Range") ]
		for btn in self.fit_bounds:
			self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelector, btn)
		
		self.fit_range_ctrl = [ wx.TextCtrl(self,wx.ID_ANY,"-5",style=wx.TE_CENTRE,size=(80,-1)),
							wx.TextCtrl(self,wx.ID_ANY,"5",style=wx.TE_CENTRE,size=(80,-1)) ]
		self.OnRadioSelector(1)
	
		fit_range_labels = [ wx.StaticText(self,wx.ID_ANY,"Min: "), wx.StaticText(self,wx.ID_ANY,"Max: ") ]
		
		panel_sizer.Add((-1,20),0,wx.EXPAND)
		panel_sizer.Add(panel_blurb,0,wx.ALIGN_CENTRE_HORIZONTAL)
		
		panel_sizer.Add((-1,20),0,wx.EXPAND)
		panel_sizer.Add(wx.StaticLine(self,-1,size=(-1,1),style=wx.LI_HORIZONTAL),
							0,wx.EXPAND|wx.LEFT|wx.RIGHT,border=20)
		panel_sizer.Add((-1,20),0,wx.EXPAND)
		
		panel_sizer.Add(self.fit_bounds[0],0,wx.EXPAND|wx.LEFT,border=30)
		panel_sizer.Add((-1,5),0,wx.EXPAND)
		panel_sizer.Add(self.fit_bounds[1],0,wx.EXPAND|wx.LEFT,border=30)
		
		panel_sizer.Add((-1,20),0,wx.EXPAND)

		for i in range(2):
			ctrl_sizer = wx.BoxSizer(wx.HORIZONTAL)
			ctrl_sizer.Add((45,-1),1,wx.EXPAND)
			ctrl_sizer.Add(fit_range_labels[i],0,wx.EXPAND)
			ctrl_sizer.Add((15,-1),0,wx.EXPAND)
			ctrl_sizer.Add(self.fit_range_ctrl[i],0,wx.EXPAND)
			ctrl_sizer.Add((45,-1),1,wx.EXPAND)
			
			panel_sizer.Add(ctrl_sizer,0,wx.EXPAND)
			panel_sizer.Add((-1,5),0,wx.EXPAND)
			
		panel_sizer.Add((-1,15),1,wx.EXPAND)
		
		btnbar = self.CreateButtonSizer(wx.OK|wx.CANCEL)
		btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
		btn_sizer.Add((20,-1),1,wx.EXPAND)
		btn_sizer.Add(btnbar,0,wx.EXPAND)
		btn_sizer.Add((20,-1),1,wx.EXPAND)
		
		panel_sizer.Add(btn_sizer,0,wx.EXPAND)
		panel_sizer.Add((-1,10),0,wx.EXPAND)
		
		self.SetSizer(panel_sizer)
		self.Layout()	
		
	def OnRadioSelector(self,event):
		## enable / disbale the text input box
		enabled = self.fit_bounds[1].GetValue()
		for ctrl in self.fit_range_ctrl:
			ctrl.Enable(enabled)
		
	def OnUpdateRange(self):
		if self.fit_bounds[0].GetValue():
			# if using full data range
			try:
				return self.parent.x_fit_array[0], self.parent.x_fit_array[-1]
			except: #no data loaded
				return [None, None]
		else:
			return [float(self.fit_range_ctrl[0].GetValue()), float(self.fit_range_ctrl[1].GetValue())]
		
		
class ElecSus_GUI_Frame(wx.Frame):
	""" Main class for this program - top-level window """
	def __init__(self,parent,title):
		""" Initialise main frame """
		wx.Frame.__init__(self,None,title=title,size=(1300,850))
		
		#ubuntu sizing:
		if os.name == 'posix':
			font = wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL)
			self.SetFont(font)
		
		# Set icons for top-left of frame, alt-tab window ...
		frame_icon = wx.IconBundle()
		frame_icon.AddIconFromFile(os.path.join(elecsus_dir,'images/elecsus_t_group.ico'), wx.BITMAP_TYPE_ANY)
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
						
	def _init_default_values(self):
		""" Initialise default values for various things ... """
		
		self.figs = []
		self.canvases = []
		self.fig_IDs = []
		
		self.hold = False
		
		self.fit_algorithm = 'Marquardt-Levenberg'
		self.already_fitting = False
		self.warnings = True

		# initialise advanced fit options dictionary
		dlg = AdvancedFitOptions(self,"Advanced Fit Options",wx.ID_ANY)
		self.advanced_fitoptions = dlg.return_all_options()
		dlg.Destroy()		
						
	def _init_plot_defaults(self):
		""" 
		List of default values for all plots. 
		Theory plots can be up to 7 panels, each of which can be turned on/off 
		"""
		
		self.plot_outputs = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix, Iy', 'Alpha Plus/Minus', 'N Plus/Minus', 'GI Plus/Minus', 'Phi']
		self.plot_ylabels = ['Transmission', 'S1', 'S2', 'S3', '$I_x$, $I_y$', r'$\alpha^{\pm}$', r'Refractive Index $n^\pm -1$', r'Group Index $\omega$ d$n^{\pm}$/d$\omega$', 'Rotation Angle $\phi$']
		self.plot_output_indices = [OutputPlotTypes.index(po) for po in self.plot_outputs]
		
		self.xrange = [detuning_defaults[0],detuning_defaults[1]] # detuning range, in GHz
		self.npoints = detuning_defaults[2] # number of detuning points to calculate
		
		# default data for plots - blank lists
		self.x_array = np.linspace(self.xrange[0],self.xrange[1],self.npoints)
		self.y_arrays = [None]*len(self.plot_outputs)
		self.x_expt_arrays = [None]*len(self.plot_outputs)
		self.y_expt_arrays = [None]*len(self.plot_outputs)
		self.x_fit_array, self.y_fit_array = None,None
		
		# plot data visible when set to True - default = view transmission (S0)
		self.display_theory_curves = [True,False,False,False,False,False,False,False,False]
		self.display_expt_curves = [False,False,False,False,False,False,False,False,False]
		
		# live plotting
		self.LiveEventsBound = False

		# residual plot settings
		self.normalised_residuals =  False
		self.residual_histogram = False
		
		#fit bounds
		self.fit_bounds = [None, None]
	
	def _init_menus(self):
		""" Initialise menu bar items """
		
		# Create menuBar object
		menuBar = wx.MenuBar()
		
		# File
		fileMenu = wx.Menu()
		fM_open = fileMenu.Append(wx.ID_OPEN, "&Open Experimental Data (csv)\tCtrl+O", "Open file for plotting and/or fitting.")
		self.Bind(wx.EVT_MENU, self.OnFileOpen, fM_open)
		fileMenu.AppendSeparator()
		fm_saveplot = fileMenu.Append(wx.ID_ANY, "Save Plot as Image", "Save Data")
		self.Bind(wx.EVT_MENU, self.OnSaveFig, fm_saveplot)
		fm_savecsv = fileMenu.Append(wx.ID_SAVE, "E&xport CSV Data\tCtrl+S", "Save CSV Data")
		self.Bind(wx.EVT_MENU, self.OnSaveCSVData, fm_savecsv)
		fm_saveconfig = fileMenu.Append(wx.ID_ANY, "Save Current Configuration", "Save Config")
		self.Bind(wx.EVT_MENU, self.OnSaveConfig, fm_saveconfig)
		fileMenu.AppendSeparator()
		fM_exit = fileMenu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", "Close window and exit program.")
		self.Bind(wx.EVT_MENU, self.OnExit, fM_exit)
		#
		
		# Edit
		editMenu = wx.Menu()
		eM_CopyTtoF = editMenu.Append(wx.ID_ANY, "Copy Parameters: &Theory to Fit", "Copy - T to F")
		eM_CopyFtoT = editMenu.Append(wx.ID_ANY, "Copy Parameters: &Fit to Theory", "Copy - F to T")
		
		self.Bind(wx.EVT_MENU, self.OnCopyParamsTtoF, eM_CopyTtoF)
		self.Bind(wx.EVT_MENU, self.OnCopyParamsFtoT, eM_CopyFtoT)
		
		# About
		aboutMenu = wx.Menu()
		aM_aboutthis = aboutMenu.Append(wx.ID_ABOUT, "&About this program", "About this program.")
		self.Bind(wx.EVT_MENU, self.OnAboutThis, aM_aboutthis)
		aM_aboutelecsus = aboutMenu.Append(wx.ID_ANY, "About &ElecSus", "About ElecSus")
		self.Bind(wx.EVT_MENU, self.OnAboutElecSus, aM_aboutelecsus)

		# View
		viewMenu = wx.Menu()
		vM_liveplot = viewMenu.Append(wx.ID_ANY, "&Live Plot", "Update Plot in real-time",kind=wx.ITEM_CHECK)
		self.Bind(wx.EVT_MENU, self.OnLivePlotting, vM_liveplot)
		# Select plot types to display

		# Theory Plot
		theoryplotMenu = wx.Menu()
		tpM_plotholdon = theoryplotMenu.Append(wx.ID_ANY, 
			"&Hold data on plot update", 
			"Select whether to hold or clear current plot data on updating the figure", 
			kind=wx.ITEM_CHECK)
		tpM_plotholdon.Check(self.hold)
		self.Bind(wx.EVT_MENU, self.OnPlotHold, tpM_plotholdon)
		#tpM_clearplot = theoryplotMenu.Append(wx.ID_ANY, "&Clear current plot data", "Clear current plot data on all axes")
		#self.Bind(wx.EVT_MENU, self.OnClearPlot, tpM_clearplot)
		
		tpM_grid = theoryplotMenu.Append(wx.ID_ANY, "&Grid on axes", "Grid", kind=wx.ITEM_CHECK)
		self.Bind(wx.EVT_MENU, self.OnGridToggleMain, tpM_grid)

		 
		## Plot selection on menus - superceded by panel buttons ##
		self.showTplotsSubMenu = wx.Menu()
		id_S0 = 2000
		id_S1 = 2001
		id_S2 = 2002
		id_S3 = 2003
		id_IxIy = 2004
		id_alpha = 2005
		id_n = 2006
		id_phi = 2007
		show_S0 = self.showTplotsSubMenu.AppendCheckItem(id_S0,"&Transmission (S0)")
		show_S1 = self.showTplotsSubMenu.AppendCheckItem(id_S1,"S&1")
		show_S2 = self.showTplotsSubMenu.AppendCheckItem(id_S2,"S&2")
		show_S3 = self.showTplotsSubMenu.AppendCheckItem(id_S3,"S&3")
		show_IxIy = self.showTplotsSubMenu.AppendCheckItem(id_IxIy,"&Ix/Iy")
		show_alpha = self.showTplotsSubMenu.AppendCheckItem(id_alpha,"&Alpha +/-")
		show_n = self.showTplotsSubMenu.AppendCheckItem(id_n,"&Refractive Index +/-")
		show_ng = self.showTplotsSubMenu.AppendCheckItem(id_n,"&Group Index +/-")
		show_phi = self.showTplotsSubMenu.AppendCheckItem(id_phi,"&Phi (Rotation Angle)")
		# bind event to check box selections
		for checkitem in [show_S0, show_S1, show_S2, show_S3, show_IxIy, show_alpha, show_n, show_ng, show_phi]:
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
		id_phi = 3007
		show_S0 = self.showEplotsSubMenu.AppendCheckItem(id_S0,"&Transmission (S0)")
		show_S1 = self.showEplotsSubMenu.AppendCheckItem(id_S1,"S&1")
		show_S2 = self.showEplotsSubMenu.AppendCheckItem(id_S2,"S&2")
		show_S3 = self.showEplotsSubMenu.AppendCheckItem(id_S3,"S&3")
		show_IxIy = self.showEplotsSubMenu.AppendCheckItem(id_IxIy,"&Ix/Iy")
		show_alpha = self.showEplotsSubMenu.AppendCheckItem(id_alpha,"&Alpha +/-")
		show_n = self.showEplotsSubMenu.AppendCheckItem(id_n,"&Refractive Index +/-")
		show_ng = self.showEplotsSubMenu.AppendCheckItem(id_n,"&Group Index +/-")
		show_phi = self.showEplotsSubMenu.AppendCheckItem(id_phi,"&Phi (Rotation Angle)")
		# bind event to check box selections
		for checkitem in [show_S0, show_S1, show_S2, show_S3, show_IxIy, show_alpha, show_n, show_ng, show_phi]:
			self.Bind(wx.EVT_MENU, self.OnShowEplots, checkitem)
		# add to parent menu
		tpM_showEplots = theoryplotMenu.AppendMenu(wx.ID_ANY,"&Experimental Curves to Show", self.showEplotsSubMenu)
		

		
		# Experimental Plot
		exptplotMenu = wx.Menu()
		epM_plotholdon = exptplotMenu.Append(wx.ID_ANY, "&Hold data on plot update", "Select whether to hold or clear current plot data on updating the figure", kind=wx.ITEM_CHECK)
		#epM_clearplot = exptplotMenu.Append(wx.ID_ANY, "&Clear current plot data", "Clear current plot data on all axes")
		#self.Bind(wx.EVT_MENU, self.OnClearPlot,epM_clearplot)
		
		epM_grid = exptplotMenu.Append(wx.ID_ANY, "&Grid on axes", "Grid", kind=wx.ITEM_CHECK)
		self.Bind(wx.EVT_MENU, self.OnGridToggleRes, epM_grid)

		epM_reshist = exptplotMenu.Append(wx.ID_ANY, "Show &Histogram of Residuals", "Res Hist", kind=wx.ITEM_CHECK)
		self.Bind(wx.EVT_MENU, self.OnResHist, epM_reshist)
		
	
		# Fit
		fitMenu = wx.Menu()
		#fitM_dummy = fitMenu.Append(wx.ID_ANY, "&Dummy", "Dummy")
		
		'''
		# Select fitting type
		self.fittypeSubMenu = wx.Menu()
		#define some id tags that we will use when items are selected
		id_ML = 1001
		id_RR = 1002
		id_SA = 1003
		self.fit_ids = [id_ML,id_RR,id_SA]
		fit_ML = self.fittypeSubMenu.AppendRadioItem(id_ML, "&Marquardt-Levenberg", "Use ML Fitting")
		fit_RR = self.fittypeSubMenu.AppendRadioItem(id_RR, "&Random Restart", "Use RR Fitting")
		fit_SA = self.fittypeSubMenu.AppendRadioItem(id_SA, "&Simulated Annealing", "Use SA Fitting")
		
		#Bind action to each item
		self.Bind(wx.EVT_MENU, self.OnFitTypeChangeMenu, fit_ML)
		self.Bind(wx.EVT_MENU, self.OnFitTypeChangeMenu, fit_RR)
		self.Bind(wx.EVT_MENU, self.OnFitTypeChangeMenu, fit_SA)
		
		# by default, select ML fitting
		self.fittypeSubMenu.Check(id_RR, True)
		self.OnFitTypeChangeMenu(1)

		# Add sub-menu to fit menu
		fitM_fittype = fitMenu.AppendMenu(wx.ID_ANY, "Fit &Type", self.fittypeSubMenu)
		'''

		fitM_dataproc = fitMenu.Append(wx.ID_ANY, "&Data Processing...", "Data Processing Options")
		self.Bind(wx.EVT_MENU, self.OnDataProcessing, fitM_dataproc)
		
		fitM_fitbounds = fitMenu.Append(wx.ID_ANY, "&Set Fit Bounds", "Set detuning bounds for fitting")
		self.Bind(wx.EVT_MENU, self.OnFitBounds, fitM_fitbounds)

		fitM_advanced = fitMenu.Append(wx.ID_ANY, "&Advanced Fit Options...", "Advanced Fit Options")
		self.Bind(wx.EVT_MENU, self.OnAdvancedOptions, fitM_advanced)
		
		fitM_warnings = fitMenu.Append(wx.ID_ANY, "&Warn about fit settings", "Warn about possible bad fit settings", kind=wx.ITEM_CHECK)
		
		self.Bind(wx.EVT_MENU, self.OnFitWarnings, fitM_warnings)
		fitM_warnings.Check(True)
		
		#initialise
		#
		
			
		#
		# ... other menu items to add...
		# 
		# Residual analysis - fit gaussian to residuals etc
		#		-- Use normalised/raw residuals (requires importing data with errorbars.. to do..)
		# Main plot - set axes limits
		#
		
		# Add menu items to menu bar
		menuBar.Append(fileMenu, "&File")
		menuBar.Append(editMenu, "&Edit")
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
		self.T_Panel = PlotToolPanel(PlotTabs,self,1)
		self.E_Panel = PlotToolPanel(PlotTabs,self,2)
		
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
		jqclogo = wx.Image(os.path.join(elecsus_dir,'images/jqc-logo.png'),wx.BITMAP_TYPE_ANY)
		#jqc_bmp = wx.StaticText(self.panel,wx.ID_ANY,"image") 
		jqc_bmp = wx.StaticBitmap(self.panel,wx.ID_ANY,wx.BitmapFromImage(jqclogo),size=(191,-1))
		
		image_sizer = wx.BoxSizer(wx.HORIZONTAL)
		image_sizer.Add((10,-1),1,wx.EXPAND)
		image_sizer.Add((90,-1),0,wx.EXPAND)
		#image_sizer.Add(elecsus_bmp,0,wx.EXPAND)
		#image_sizer.Add((10,-1),0,wx.EXPAND)
		image_sizer.Add(jqc_bmp,0,wx.EXPAND)
		image_sizer.Add((90,-1),0,wx.EXPAND)
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
		
		
		
		## Plot selection buttons
		plot_selection_text = wx.StaticText(self.panel,wx.ID_ANY,"Data types to show:")

		Tplot_selection_list = wx.Button(self.panel,wx.ID_ANY, 'Theory Plots', size=(100,BtnSize))
		self.Bind(wx.EVT_BUTTON, self.OnTheoryPlotSelection, Tplot_selection_list)

		Eplot_selection_list = wx.Button(self.panel,wx.ID_ANY, 'Fit Plots', size=(100,BtnSize))
		self.Bind(wx.EVT_BUTTON, self.OnFitPlotSelection, Eplot_selection_list)
		
		#self.Bind(EVT_)

		#for checkitem in [show_S0, show_S1, show_S2, show_S3, show_IxIy, show_alpha, show_n]:
		#	self.Bind(wx.EVT_MENU, self.OnShowTplots, checkitem)
		## add to parent menu
		#tpM_showTplots = theoryplotMenu.AppendMenu(wx.ID_ANY,"&Theory Curves to Show", self.showTplotsSubMenu)


		plot_selection_sizer = wx.BoxSizer(wx.HORIZONTAL)
		plot_selection_sizer.Add((10,-1),0,wx.EXPAND)
		plot_selection_sizer.Add(plot_selection_text,0, wx.ALIGN_CENTER_VERTICAL)
		plot_selection_sizer.Add((10,-1),1,wx.EXPAND)
		plot_selection_sizer.Add(Tplot_selection_list,0, wx.EXPAND)
		plot_selection_sizer.Add((20,-1),0,wx.EXPAND)
		plot_selection_sizer.Add(Eplot_selection_list,0, wx.EXPAND)
		plot_selection_sizer.Add((20,-1),1,wx.EXPAND)


	## Create sizer for right-hand side of panel
		button_sizer = wx.BoxSizer(wx.VERTICAL)
		button_sizer.Add((-1,10),0,wx.EXPAND)
		button_sizer.Add(image_sizer,0,wx.EXPAND)
		button_sizer.Add((-1,5),0,wx.EXPAND)
		button_sizer.Add(plot_selection_sizer,0,wx.EXPAND)#,0,wx.EXPAND)
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
		
				
#		
## Actions for events
#
	
##
#### General Actions - Close window, Save Data etc...
##

	def OnExit(self,event):
		""" What to do when the window is closed """
		## explicitly close all figures (bug with matplotlib and wx??)
		plt.close('all')
		self.Destroy()
		app.ExitMainLoop()
		
	def OnAboutThis(self,event):
		""" Show a message box about the program """
		dlg = wx.MessageDialog(self, "A graphical interface for ElecSus, a program to calculate the weak-probe electric susceptibility of alkali atoms.\n\n Written in python using wxPython and matplotlib.", "About", wx.OK)
		if dlg.ShowModal() == wx.ID_OK:
			dlg.Destroy()

	def OnAboutElecSus(self,event):
		""" Show a message box about ElecSus """
		dlg = wx.MessageDialog(self, NOTICE.noticestring, "About", wx.OK)
		if dlg.ShowModal() == wx.ID_OK:
			dlg.Destroy()
					
	def OnSaveCSVData(self,event):
		""" 
		Method to save the main plot data traces as a n-column csv file. 
		
		*All* data is saved, whether or not it is displayed 
		
		This method selects the file name, then passes that to SaveTheoryCurves()
		"""
		
		SaveFileDialog = wx.FileDialog(self,"Save Output File", "./", "Outputs",
			"CSV files (*.csv)|*.csv", wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
		
		if SaveFileDialog.ShowModal() == wx.ID_OK:
			output_filename = SaveFileDialog.GetPath()
			SaveFileDialog.Destroy()
			
			print output_filename
			
			# don't need this - FD_OVERWRITE_PROMPT does the same job
			#check for overwrite current files
			#if os.path.isfile(output_filename):
			#	OverwriteDialog = wx.MessageDialog(self,"Warning: file exists already! Overwrite?",\
			#		"Overwrite?",wx.YES_NO|wx.NO_DEFAULT)
			#	
			#	if OverwriteDialog.ShowModal() == wx.NO:
			#		OverwriteDialog.Destroy()
			#		return # exit without saving
			#	else:
			#		OverwriteDialog.Destroy()
			
			## do save
			self.SaveTheoryCurves(output_filename)

	def SaveTheoryCurves(self,filename):
		""" 
		Method to actually do the saving of xy data to csv file. 
		Separate to the OnSaveCSVData method in case a filename is automatically chosen
		- i.e. when fitting data, there is an option to autosave the results which calls this method
		"""
		xy_data = zip(self.x_array,self.y_arrays[0],self.y_arrays[1],self.y_arrays[2],self.y_arrays[3],\
			self.y_arrays[4][0], self.y_arrays[4][1], self.y_arrays[5][0], self.y_arrays[5][1], \
			self.y_arrays[6][0], self.y_arrays[6][1], self.y_arrays[7][0], self.y_arrays[7][1], self.y_arrays[8])
		success = write_CSV(xy_data, filename, titles=['Detuning']+OutputTypes)
		if not success:
			problem_dlg = wx.MessageDialog(self, "There was an error saving the data...", "Error saving", wx.OK|wx.ICON_ERROR)
			problem_dlg.ShowModal()
			
	def OnSaveConfig(self,event):
		""" 
		Save present configuration settings - plot types, data ranges, parameters ... for faster repeating of common tasks
		"""
		
		dlg = wx.MessageDialog(self, "Save current configuration of the program - theory/fit parameters, plot settings etc...\n\nNot implemented yet...", "No no no", wx.OK)
		dlg.ShowModal()

		data_dump = [0]
		
		# get plot settings
		# ...
		# get theory tab settings
		# ...
		# get fit tab settings
		# ...
		# get mist settings
		# ...
		
		## filedialog window for selecting filename/location
		filename = './elecsus_gui_config.dat'
		
		# save data in python-readable (binary) format using pickle module
		with open(filename,'wb') as file_obj:
			pickle.dump(data_dump, file_obj)
		
	def OnLoadConfig(self,event):
		"""
		Load previous configuration settings from pickle file, for picking up where you left off...
		"""
		dlg = wx.MessageDialog(self, "Load previous configuration of the program - theory/fit parameters, plot settings etc...\n\nNot implemented yet...", "No no no", wx.OK)
		dlg.ShowModal()
		
		## do the opposite of save config ...
		
	def OnSaveFig(self,event):
		""" Basically the same as saving the figure by the toolbar button. """
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

	def OnCopyParamsTtoF(self,event):
		self.CopyParams(1)

	def OnCopyParamsFtoT(self,event):
		self.CopyParams(2)

	def CopyParams(self,order):
		T_params = self.ThyPanel.fixed_paramlist_inputs + self.ThyPanel.fit_paramlist_inputs
		F_params = self.FitPanel.fixed_paramlist_inputs + self.FitPanel.fit_paramlist_inputs

		for fp, tp in zip(F_params,T_params):
			if order==1:
				# Theory to Fit
				fp.SetValue(tp.GetValue())
			elif order==2:
				tp.SetValue(fp.GetValue())
				
		print 'Parameters Copied'
		
##
####			
####### Plot-specific actions
####
##
	def OnTheoryPlotSelection(self,event):
		""" Select plots to show, via the popup button on the main panel """
		win = PlotSelectionPopUp(self.panel,wx.SIMPLE_BORDER,self,'Theory')
 		
		btn = event.GetEventObject()
		pos = btn.ClientToScreen( (0,0) )
		sz =  btn.GetSize()
		win.Position(pos, (0, sz[1]))
        
		win.Popup()

	def OnFitPlotSelection(self,event):
		""" Select plots to show, via the popup button on the main panel """
		win = PlotSelectionPopUp(self.panel,wx.SIMPLE_BORDER,self,'Fit')
 		
		btn = event.GetEventObject()
		pos = btn.ClientToScreen( (0,0) )
		sz =  btn.GetSize()
		win.Position(pos, (0, sz[1]))
        
		win.Popup()

	def OnLivePlotting(self,event):
		""" 
		Turn on/off live plotting when parameters are changed. 
		Main plot updates in pseudo-realtime, rather than having to click 'compute' every time
		
		This method binds or unbinds a call to the compute button each time a control is changed.
		"""
		LivePlotOn = bool(event.Checked())
		
		if LivePlotOn:
			for ctrl in self.ThyPanel.DetuningCtrl + self.ThyPanel.fit_paramlist_inputs:
				self.Bind(EVT_FLOATSPIN, self.OnComputeButton, ctrl)
			
			for ctrl in self.ThyPanel.fixed_paramlist_inputs[0:2]:
				self.Bind(wx.EVT_COMBOBOX, self.OnComputeButton, ctrl)
			
			self.LiveEventsBound = True
		else:
			# unbind the events, if currently bound
			if self.LiveEventsBound:
				# unbind
				for ctrl in self.ThyPanel.DetuningCtrl + self.ThyPanel.fit_paramlist_inputs:
					self.Unbind(EVT_FLOATSPIN, ctrl)

				for ctrl in self.ThyPanel.fixed_paramlist_inputs[0:2]:
					self.Unbind(EVT_FLOATSPIN, ctrl)
				
			self.LiveEventsBound = False
	
	def OnFileOpen(self,event):
		""" 
		Open a csv data file and plot the data. Detuning is assumed to be in GHz. 
		Vertical units are assumed to be the same as in the theory curves 
		"""
		self.dirname= ''
		dlg_choice = wx.SingleChoiceDialog(self,"Choose type of data to be imported","Data import",choices=OutputTypes)
		
		# wait for OK to be clicked
		if dlg_choice.ShowModal() == wx.ID_OK:
			choice = dlg_choice.GetSelection()
			#print 'Choice:', choice
			self.expt_type = OutputTypes[choice]
			# use the choice index to select which axes the data appears on	- may be different 
			# if axes order is rearranged later?
			self.choice_index = OutputTypes_index[choice]
			#print self.choice_index
			
			dlg_choice.Destroy()
		
			dlg_open = wx.FileDialog(self,"Choose 2-column csv file (Detuning, Transmission)",
									self.dirname,"","*.csv",wx.OPEN)
			
			# if OK button clicked, open and read file
			if dlg_open.ShowModal() == wx.ID_OK:
				#set experimental display on, and update menus
				self.display_expt_curves[self.choice_index] = True
				self.showEplotsSubMenu.GetMenuItems()[self.choice_index].Check(True)
			
				self.filename = dlg_open.GetFilename()
				self.dirname = dlg_open.GetDirectory()
				#call read
				self.x_expt_arrays[self.choice_index],self.y_expt_arrays[self.choice_index] = read_CSV(os.path.join(self.dirname,self.filename),spacing=0)
				
				#overwrite fit_array data - i.e. last data to be loaded
				self.x_fit_array = self.x_expt_arrays[self.choice_index]
				self.y_fit_array = self.y_expt_arrays[self.choice_index]
				
				# implicit that the fit type is the same as last data imported
				self.fit_datatype = self.expt_type
				
				## create main plot				
				self.OnCreateAxes(self.figs[0],self.canvases[0],clear_current=True)
				
				
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
				#if isinstance(ye, (list,tuple)):
				#	for xi,yi in zip(xe,ye):
				#		ax.plot(xi,yi,color=d_grey)
				#else:
				ax.plot(xe,ye,color=d_olive)
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
		"""
		Create a plot with main data, residuals and (optionally) histogram of residuals, using
		subplot2grid in matplotlib 
		"""
		
		fig.clf()
		
		print 'Debugging...'
		print self.fit_datatype
		print self.y_optimised
		
		#normalised_residuals = False
		if self.normalised_residuals:
			### not done yet! -- requires error bars in imported data
			residuals = 100*(self.y_fit_array-self.y_optimised)
		else:
			residuals = 100*(self.y_fit_array-self.y_optimised)
		
		fig = plt.figure(2)
		yy = 4
		xx = 6
		if self.residual_histogram:
			ax_main = plt.subplot2grid((yy,xx),(0,0),colspan=xx-1,rowspan=yy-1)
			ax_residual = plt.subplot2grid((yy,xx),(yy-1,0),colspan=xx-1,sharex=ax_main)
			ax_hist = plt.subplot2grid((yy,xx), (yy-1,xx-1), sharey=ax_residual)

			plt.setp(ax_hist.get_yticklabels(),visible=False)
			ax_hist.set_xticklabels([])

		else:
			ax_main = plt.subplot2grid((yy,xx),(0,0),colspan=xx,rowspan=yy-1)
			ax_residual = plt.subplot2grid((yy,xx),(yy-1,0),colspan=xx,sharex=ax_main)
			
		plt.setp(ax_main.get_xticklabels(),visible=False)
		
		ax_residual.set_xlabel('Detuning (GHz)')
		ax_residual.set_ylabel('Residuals (%)')

		ax_main.set_ylabel(self.expt_type)
		
		ax_main.plot(self.x_fit_array,self.y_fit_array,color=d_olive)
		print len(self.x_fit_array), len(self.y_optimised)
		ax_main.plot(self.x_fit_array,self.y_optimised)
		ax_residual.plot(self.x_fit_array,residuals,lw=1.25)
		ax_residual.axhline(0,color='k',linestyle='dashed')
		
		if self.residual_histogram:
			bins = 25
			ax_hist.hist(residuals, bins=bins, orientation='horizontal')

		ax_main.autoscale_view(tight=True)
		
		self.draw_fig(fig,canvas)

	def OnResHist(self,event):
		""" Turn histogram of residuals on/off """
		self.residual_histogram = bool(event.Checked())
		self.OnCreateResidualPlot(self.figs[1],self.canvases[1])
					
	def OnPlotHold(self,event):
		""" 
		Toggle plot hold (keep data on updating figure) on/off
		Allows multiple data sets to be shown on same plot
		"""
		#self.PlotHold = bool(event.Checked())
		#self.figs[0].hold(self.PlotHold)
		dlg = wx.MessageDialog(self, "Not implemented yet...", "No no no", wx.OK)
		dlg.ShowModal()
			
	#def OnClearPlot(self,event):
	#	""" Get list of all axes in figure and clear them all """
	#	for fig in self.figs[0]:
	#		for ax in fig.axes:
	#			ax.cla()		
	
	def OnPlotLegend(self,event):
		""" Toggle plot legend on/off """
		self.legendOn = bool(event.Checked())
	
	def OnGridToggleMain(self,event):
		""" Toggle axes grid on main plot """
		for ax in self.figs[0].axes:
			ax.grid(bool(event.Checked()))
		self.draw_fig(self.figs[0],self.canvases[0])
					
	def OnGridToggleRes(self,event):
		""" Toggle axes grid on residuals plot """
		for ax in self.figs[1].axes:
			ax.grid(bool(event.Checked()))
		self.draw_fig(self.figs[1],self.canvases[1])

	def OnShowTplots(self,event):
		""" Action when plot type is changed from the menu """
		for ii, item in enumerate(self.showTplotsSubMenu.GetMenuItems()):
			if item.IsChecked():
				self.display_theory_curves[ii] = True
			else:
				self.display_theory_curves[ii] = False
		
		#redraw plot
		self.OnCreateAxes(self.figs[0],self.canvases[0])

	def OnShowEplots(self,event):
		""" Action when plot type is changed from the menu """
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
	
		for can in self.canvases:
			can.draw()
		#self.SendSizeEvent()
##
#### Fitting specific actions
##
	def OnFitWarnings(self,event):
		""" 
		Warn about bad fitting practices - e.g. when there are many data points to fit, or many
		fit parameters where the fit algorithm could be improved.
		"""
		self.warnings = bool(event.Checked())
	
	def OnFitTypeChangePanel(self,event):
		""" Action to perform when fit type is changed """
		#print self.fittypeSubMenu.GetMenuItems()
		print ''
		for panelitem in self.FitPanel.fit_types:
			if panelitem.GetValue():
				self.fit_algorithm = panelitem.GetLabel()
				#self.fittypeSubMenu.Check(idfit,True)
				#self.OnFitTypeChangeMenu(1)
				print 'Fit Algorithm changed to:', self.fit_algorithm

			#self.fittypeSubMenu.Check(idfit,False)
			#print menuitem.IsChecked()

	def OnDataProcessing(self,event):
		""" Open the dialog box for binning / smoothing data """
		if self.x_fit_array is not None:
			dlg = DataProcessingDlg(self, "Data Processing", wx.ID_ANY)
			if dlg.Show() == wx.ID_OK:
				dlg.Destroy()	
		else:
			dlg = wx.MessageDialog(self,"Can't process data that hasn't been loaded...", "Nope.", wx.OK|wx.ICON_INFORMATION)
			dlg.ShowModal()
	
	def OnFitBounds(self,event):
		dlg = FitBoundsDialog(self,"Fit Bounds", wx.ID_ANY)
		
		if dlg.ShowModal() == wx.ID_OK:
			self.fit_bounds = dlg.OnUpdateRange()
			#print self.fit_bounds		
		
	def OnAdvancedOptions(self,event):
		""" Open the dialog box for adjusting advanced fit options """
		dlg = AdvancedFitOptions(self,"Advanced Fit Options",wx.ID_ANY)
		
		# Show() rather than ShowModal() - doesn't halt program flow
		if dlg.ShowModal() == wx.ID_OK:
			self.advanced_fitoptions = dlg.return_all_options()
			print self.advanced_fitoptions
				
##
#### Main purpose of this - call elecsus with various settings!
##
	def OnFitButton(self,event):
		""" Call elecsus to fit data. Some sanity checking takes place first. """

		## check for things that will prevent fitting from working, e.g. no data loaded - halt fitting if found
		if self.y_fit_array == None:
			#warn about no data present
			dlg = wx.MessageDialog(self, "No experimental data has been loaded, cannot proceed with fitting...", "No no no", wx.OK|wx.ICON_EXCLAMATION)
			dlg.ShowModal()
			return
		
		self.fit_bools = [checkbox.IsChecked() for checkbox in self.FitPanel.fit_paramlist_bools]
		if self.fit_bools.count(True) == 0:
			dlg = wx.MessageDialog(self, "No fit parameters are floating, cannot proceed with fitting...", "No no no", wx.OK|wx.ICON_EXCLAMATION)
			dlg.ShowModal()
			return
		
		## check for non-optimal conditions and warn user
		if self.warnings:
			## if number of booleans > 3 and ML fitting selected, warn about fitting methods
			if self.fit_bools.count(True) > 3 and self.fit_algorithm=='Marquardt-Levenberg':
				dlg = wx.MessageDialog(self, "The number of fitted parameters is large and the Marquardt-Levenberg algorithm is selected. There is a high probability that the fit will return a local minimum, rather than the global minimum. \n\nTo find the global minimum more reliably, consider changing to either Random-Restart or Simulated Annealing algorithms.\n\n Continue with fitting anyway?", "Warning", wx.YES|wx.NO|wx.ICON_WARNING)
				
				if dlg.ShowModal() == wx.ID_NO:
					return
				
			## if large number of data points
			if len(self.x_fit_array) > 5000:
				dlg = wx.MessageDialog(self, "The number of data points is quite high and fitting may be very slow, especially with Random-Restart or Simulated Annealing methods. \n\nConsider reducing the number of data points by binning data (Menu Fit --> Data Processing... --> Bin Data).\n\n Continue with fitting anyway?", "Warning", wx.YES|wx.NO|wx.ICON_WARNING)

				if dlg.ShowModal() == wx.ID_NO:
					return
			
		# If all conditions satisfied, start the fitting process
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
			
		elif calc_or_fit == 'Fit':
			panel = self.FitPanel
			
		
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
			print self.params
			
			spectrum_data = elecsus.calculate(self.x_array,self.params)
			self.y_arrays[0:4] = spectrum_data[0:4] #S0,S1,S2,S3
			self.y_arrays[4] = [spectrum_data[4],spectrum_data[5]] # Ix,Iy
			self.y_arrays[5] = [spectrum_data[9],spectrum_data[10]] # alpha+,alpha-
			self.y_arrays[6] = [spectrum_data[6],spectrum_data[7]] # n-,n+
			self.y_arrays[8] = spectrum_data[8] # Faraday rotation angle phi
			self.y_arrays[7] = [spectrum_data[11],spectrum_data[12]]
			
			#print self.y_arrays
			#for i, array in enumerate(self.y_arrays):
			#	print i, type(array)
			self.OnCreateAxes(self.figs[0],self.canvases[0])
			
		elif calc_or_fit == 'Fit':
			### Use another thread for running elecsus fitting so the main panel doesn't stop responding
			if not self.already_fitting:
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

				
				##
				#### Run Fitting algorithm
				##
			
				#initialise thread for fitting
				fitThread = FittingThread(self)
				
				# start thread
				fitThread.start()
				
				# only allow one fit to run at any time
				self.already_fitting = True
				
				## show dialog box with status bar?
				#self.fitting_dlg = ProgressBarFrame(self,"Fit in progress",range=100)
				#fitProgressThread = ProgressThread(self)
				#fitProgressThread.start()
				
				self.fitting_dlg = wx.ProgressDialog("Fit in progress","Fitting in progress. Please wait for fitting to complete. This may take a while.",maximum=100,style=wx.PD_ELAPSED_TIME)
				self.fitting_dlg.Show(True)
				fitProgressThread = ProgressThread(self)
				fitProgressThread.start()
				
				#time.sleep(0.5)
				#self.FitInformation.write('\n\nRunning Fit...')

				# After the fit completes, the thread triggers an EVT_FIT_COMPLETE which 
				# calls (in the main thread) the OnFitCompleted() method
			else:
				dlg = wx.MessageDialog(self,"Fit already in progress. Only one fit may be run at the same time.", "Patience required...",style=wx.OK|wx.ICON_ERROR)
				dlg.ShowModal()
				## dialog box for already fitting...
			


	def OnFitCompleted(self,event):
		""" Task to complete after fitting completed """
		self.FitInformation.write('\n============== Fit Completed ============ \n\n')
			
		# reverse the re-ordering we did before...
		self.opt_params_new = self.fit_params
		for i, el in enumerate(self.re_order_list):
			self.opt_params_new[el] = self.opt_params[i+2]
		
		#self.opt_params = self.opt_params_new
		
		## Add information to the Fitting text box
		self.FitInformation.write('Optimised parameters are (Fitted):\n')
		[ self.FitInformation.write('\t'+param_name.ljust(30)+str(param_value).ljust(10)+'('+str(boo)+')\n') for param_name, param_value,boo in zip(fittable_parameterlist, self.opt_params_new, self.fit_bools) ]
		
		# Add RMS info to the Fitting text box
		self.FitInformation.write('\nRMS error between theory and experiment (%): '+str(round(self.rms*100,2))+'\n')
		self.FitInformation.write('\n-------------------------------------------------------------------------------------\n')
		
		
		
		# use fitted parameters to calculate new theory arrays
		print '\n\n'
		print time.ctime()
		print 'Calling ElecSus for single calculation with optimised parameters:'
		print self.opt_params_new
		
		self.y_optimised = elecsus.calculate(self.x_fit_array,self.opt_params,OutputType=self.fit_datatype)
		spectrum_data = elecsus.calculate(self.x_array,self.opt_params)
		self.y_arrays[0:4] = spectrum_data[0:4] #S0,S1,S2,S3
		self.y_arrays[4] = [spectrum_data[4],spectrum_data[5]] # Ix,Iy
		self.y_arrays[5] = [spectrum_data[9],spectrum_data[10]] # alpha+,alpha-
		self.y_arrays[6] = [spectrum_data[6],spectrum_data[7]] # n-,n+
		self.y_arrays[7] = spectrum_data[8] # Faraday rotation angle phi
		
		# update main plot
		# turn on theory curve for the type of plot we just fitted..
		index_fit_datatype = OutputTypes_index[OutputTypes.index(self.fit_datatype)]
		self.display_theory_curves[index_fit_datatype] = True
		
		# Reset the 'already fitting' flag
		self.already_fitting = False
		
		# close the 'fitting in progress' dialog box
		self.fitting_dlg.Show(False)
		self.fitting_dlg.Destroy()

		# update the figure
		self.OnCreateAxes(self.figs[0],self.canvases[0])			
		print 'Updating main plot...'
		
		# ask whether to update the initial fitting parameters
		dlg = wx.MessageDialog(self, 
			"Fit completed successfully. Details can be found in the Fitting Information Tab.\
			\n\nReplace initial fit parameters with optimised parameters?", 
			"Fit Complete",
			style=wx.YES|wx.NO)
		if dlg.ShowModal() == wx.ID_YES:
			for i, value in enumerate(self.opt_params_new):
				self.FitPanel.fit_paramlist_inputs[i].SetValue(value)
		dlg.Destroy()

		# create/update residuals plot
		self.OnCreateResidualPlot(self.figs[1],self.canvases[1])				
		
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

def write_CSV(xy,filename,titles=None):
	""" Method for writing csv data with arbitrary
		number of columns to filename.
		takes in xy, which should be of the form [[x1,y1],[x2,y2] ...]
		this can be done by zipping arrays, e.g.
			xy = zip(x,y,z)
			where x, y and z are 1d arrays
		This creates a csv file that looks like
		x1, y1, z1, ...
		x2, y2, z2, ...
		...
		xN, yN, zN, ...
		
		The optional titles writes a single header row at the top of the csv file.
	"""	
	
	try:
		with open(filename, 'wb') as csvfile:
			csv_writer = csv.writer(csvfile,delimiter=',')
			if titles is not None:
				csv_writer.writerow(titles)
			for xy_line in xy:
				csv_writer.writerow(xy_line)
		return True
	except:
		return False
	
	
## Run the thing...
def start():
	print 'Starting ElecSus GUI...'
	app = wx.App(redirect=False)
	frame = ElecSus_GUI_Frame(None,"ElecSus GUI")
				
	frame.Centre()
	frame.Show()
	app.MainLoop()
	print '...Closed ElecSus GUI'
	
if __name__ == '__main__':
	start()
	