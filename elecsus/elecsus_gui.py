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

v3.1.0__dev (2021-10-25)
    -- elecsus_gui.py reformat to align with PEP directives (work in progress)

v3.0.7 (2019-10-22)
    -- Large speed improvement for electric field calculations.
    -- Fixed bug to allow data saving in python 3.
    -- Bug with relative paths fixed.
    -- Reduced number of initial points in RR fitting routine.
    -- Fixed compatibility issue with matplotlib versions > 3.
    -- Fixed the fitting test modules

v3.0.6 (2018-04-12)
    --	Bug in the data processing module (libs/data_proc.py) fixed
    -- Minor change to color cycler to support new version of matplotlib (v2.2)

v3.0.5 (2018-02-27)
    --	 Bug fix with some menu items not working properly

v3.0.4 (2018-02-19)
    -- Support for python 3.x added (maintains compatibility with python 2.7)

v3.0.3 (2017-12-06)
    -- Minor fixes to GUI for file input/output not working properly and an error that stopped fitting working

v3.0.2 (2017-11-14)
    -- Minor fix: changed a wx.OPEN to wx.FD_OPEN that affected newer versions of wx FileDialogs not opening properly

v3.0.1 (2017-09-25)
    -- Fix missing square-root in solve_dielectric.py
    -- Implement proper fine-structure constants properly for all alkalis

v3.0.0 (2017-08-18)
    -- MAJOR OVERHAUL to add many new program features and streamline a lot of others. This, unfortunately but
    necessarily, makes version 3 incompatible with previous versions. The changes are detailed below.
    The phyics are discussed in the new paper, a pre-print can be found at https://arxiv.org/abs/1708.05305

    -- First major change is to allow direct propagation of electric fields through the medium, which can enable
    simulating, for example: magnetic field gradients; polarising optical elements; birefringence or other
    optical imperfections. The consequence of this change is that polarisation is now defined in terms of the cartesian
    components of electric field and the phase difference between them, rather than the somewhat opaque notation
    in previous versions.

    -- Second major change is to allow arbitrary angles between the magnetic field axis and the light propagation
    axis (in 3D). This allows for a much broader set of physics to be investigated, since it removes a large
    constraint of previous versions. At zero degrees between the two fields we recover the standard Faraday geometry
    of previous ElecSus versions. At 90 degrees, we have the Voigt geometry (which also has analytic solutions,
    so is fast to run). We can also calculate for any angle, but BE WARNED, this is very computationally intensive,
    since the solutions for propagation are no longer analytic and must be calculated for each detuning point
    separately (with a corresponding slow-down in performance, typically of around a factor of 1000!).

    -- Finally, the fitting routines have been completely rewritten to use parameter dictionaries. We now utilise
    the lmfit module (https://lmfit.github.io/lmfit-py/). Performance is broadly similar (since lmfit also runs over
    the top of scipy.optimize modules), but there are many advantages of this model: Firstly, all parameters can now
    be selected to vary during a fit or be held constant, and bounds for these parameters can be given to prevent
    unphysical values from being returned. In addition, the differential_evolution solver is now availablle, which
    is very efficient at finding the global minimum for multi-parameter fits (we leave in random_restat and
    simulated_annealing for now, though these might disappear in future versions as it appears differential_evolution
    is better in all tested cases so far...).

    -- Other changes include better code documentation and minor bug fixes in many places.

v2.2.0 (2016-02-10)
    -- GUI version number and ElecSus version number are now the same
    -- Since Enthought Canopy now ships with wxPython version 3.0, GUI has been
        updated to work with this version of wxPython. All previous functionality should
        remain the same, but a few things have changed:
            -- Theory/Fit plot selections are no longer Transient Popups - they are now Dialog Boxes
            -- Default plot scaling may look a bit odd when viewing on small resolution monitors -
                not sure what the real fix for this is, but as a workaround, resizing the GUI window
                should fix this
    -- Added ability to use experimental detuning axis on theory curve,
        for direct comparison when outputting data (Menu > Edit > Use Experimental Detuning Axis)
    -- Added ability to turn off automatic scaling of axes (Menu > View > Autoscale)
    -- Fixed an issue where save files would not save correctly after performing a fit
    -- Minor fix for an issue where starting from the python interpreter would not exit properly on clicking the 'X'
    -- Corrected some incorrect tooltips
    -- Added show_versions() method, which shows the currently used version numbers of
        modules important to running this program

v1.0.1 (2015-10-23)
    -- minor bug fix where the plot selection popups would not display the Phi plots
v1.0.0 (2015-09-03)


A GUI based on wxpython for ElecSus, intended to augment/replace the 
runcard method of calling ElecSus.

See https://arxiv.org/abs/1708.05305 for the publication and more details on the physics of ElecSus.

Requirements:
    python, matplotlib, numpy, scipy

    tested on:
    Windows 8.1, Windows 10
        python 3.6 - 64-bit
        wxpython 4.1
        matplotlib 1.4.3, 1.5.1
        numpy 1.9.2
        scipy 0.15.1, 0.16.1
        lmfit 0.9.5

    Currently installed version numbers can be shown by running the show_versions() method in this module

    Requires
    --------

    wxpython version >2.8
        In newer versions, there are bugs with either matplotlib or wxpython that cause figures
        not to fill the entire panel - a temporary solution is to resize the panel.


LICENSE info: APACHE version 2

James Keaveney, Mark Zentile and co-authors
2017-21
"""
# py 2.7 compatibility
from __future__ import (division, print_function, absolute_import)

__version__ = '3.1.0'

# !/usr/bin/env python
import csv
import lmfit as lm
import os
import numpy as np
import pickle as pickle
import psutil
import scipy as sp
import sys
import threading
import time
from cycler import cycler
from matplotlib import rc

import matplotlib as mpl

# mpl.use('WxAgg')  # needs to be before pyplot import?
import matplotlib.pyplot as plt

# most important import - if wx doesn't work, the gui doesn't work!
# TODO: Update for newer version numbers
try:
    import wx
except ImportError:
    print("wxPython cannot be imported")
    print("wxPython >=2.8 needs to be installed for this program to work! \n\
    It is not currently possible to install this automatically through pip/easy_install.\n")
    if os.name == 'posix':
        print("For Ubuntu/Debian, wxPython is not supported in Enthought Canopy.\n\
        Instead, use the system python distribution (/usr/bin/python) and install through apt:\n\n\
        >   (sudo) apt-get install python-wxgtk2.8 python-wxtools wx2.8-i18n libwxgtk2.8-dev libgtk2.0-dev")
    else:
        print("For Windows, recommended install is using Enthought Canopy")
    raise ImportError

if 'phoenix' in wx.PlatformInfo:
    import wx.lib.agw.aui as aui
else:
    import wx.aui as aui

# Import lots of other wx stuff
import wx.lib.scrolledpanel as scrolled
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN
# Matplotlib/wx integration
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg as Toolbar
# Define a new event type which will be triggered when fitting is completed
import wx.lib.newevent

# Import local elecsus modules
from .elecsus_methods import calculate, fit_data
from .libs import NOTICE
from .libs import data_proc
from .libs.durhamcolours import *
from .libs.durhamcolours import cols as durhamcols
from .libs import polarisation_animation_mpl as pol_ani

# preamble.py includes tooltip text, default values, labels...
from .libs.preamble import *

# wx.Log.SetLogLevel(0)  # Ignore warnings  #TODO: check this is needed...

FitCompleteEvent, EVT_FIT_COMPLETE = wx.lib.newevent.NewEvent()

# use relative file paths
elecsus_dir = os.path.dirname(__file__)

# replace default matplotlib text and color sequence with durham colours
plt.rc('font', **{'family': 'Serif', 'serif': ['Times New Roman']})
params = {'axes.labelsize': 14, 'xtick.labelsize': 13, 'ytick.labelsize': 13, 'legend.fontsize': 12,
          'mathtext.fontset': 'cm', 'mathtext.rm': 'serif'}
plt.rcParams.update(params)
rc('text', usetex=False)
rc('font', **{'family': 'serif', 'size': 14})
rc('lines', linewidth=2)
rc('axes', prop_cycle=(cycler('color', durhamcols)))


def show_versions():
    """ Shows installed version numbers """
    print('Packages required for GUI: (this displays currently installed version numbers)')
    print('\tElecSus: ', __version__)
    print('\tWxPython: ', wx.__version__)
    print('\tNumpy: ', np.__version__)
    print('\tMatplotlib: ', mpl.__version__)
    print('Required for fitting (in addition to above):')
    print('\tScipy: ', sp.__version__)
    print('\tPSUtil: ', psutil.__version__)
    print('\tLMfit: ', lm.__version__)


class PlotSelectionPopUp(wx.PopupTransientWindow):
    """ Popup box to handle which plots are displayed. """

    def __init__(self, parent, style, mainwin, plottype):
        wx.PopupTransientWindow.__init__(self, parent, style)

        self.mainwin = mainwin
        self.plottype = plottype

        win = wx.Panel(self)  # ,wx.ID_ANY,pos=(0,0),size=(180,200),style=0)

        self.selection = wx.CheckListBox(win, wx.ID_ANY, choices=OutputPlotTypes, size=(150, -1))  # ,pos=(0,0))
        # self.win.Bind(wx.EVT_CHECKLISTBOX, self.OnTicked, self.selection)
        self.selection.Bind(wx.EVT_CHECKLISTBOX, self.OnTicked)

        if plottype == 'Theory':
            display_curves = self.mainwin.display_theory_curves
        else:
            display_curves = self.mainwin.display_expt_curves

        checked_items = []
        for i in range(len(display_curves)):
            if display_curves[i]:
                checked_items.append(i)
        self.selection.SetChecked(checked_items)

        self.SetSize(self.selection.GetSize() + (10, 10))
        self.SetMinSize(self.selection.GetSize() + (10, 10))

        popup_sizer = wx.BoxSizer(wx.VERTICAL)
        popup_sizer.Add(self.selection, 0, wx.ALL, 7)

        win.SetSizer(popup_sizer)
        popup_sizer.Fit(win)
        self.Layout()

    def OnDismiss(self):
        """
        Action to perform when the popup loses focus and is closed.

        Gets the tick box values, updates the main plot and changes
        the menu items to match the popup box
        """
        # re=plot the figure
        # self.OnTicked(1)

        # self.mainwin.Refresh()

        self.Dismiss()

    def OnTicked(self, event):
        """
        Action to perform when tick boxes are ticked

        Gets the tick box values, updates the main plot and changes
        the menu items to match the popup box
        """
        # print 'Ticked Event'

        if self.plottype == 'Theory':
            items = self.selection.GetChecked()
            self.mainwin.display_theory_curves = [False] * 9
            for menuitem in self.mainwin.showTplotsSubMenu.GetMenuItems():
                menuitem.Check(False)
            for item in items:
                self.mainwin.display_theory_curves[item] = True
                self.mainwin.showTplotsSubMenu.GetMenuItems()[item].Check(True)
            print(self.mainwin.display_theory_curves)
        elif self.plottype == 'Fit':
            items = self.selection.GetChecked()
            self.mainwin.display_expt_curves = [False] * 9
            for menuitem in self.mainwin.show_e_plots_submenu.GetMenuItems():
                menuitem.Check(False)
            for item in items:
                self.mainwin.display_expt_curves[item] = True
            self.mainwin.show_e_plots_submenu.GetMenuItems()[items[-1]].Check(True)
            print(self.mainwin.display_expt_curves)

        self.mainwin.OnCreateAxes(self.mainwin.figs[0])
        event.Skip()


class PlotSelectionDialog(wx.Dialog):
    """ Popup box to handle which plots are displayed. """

    def __init__(self, parent, mainwin, title, plottype, pos):
        wx.Dialog.__init__(self, parent, wx.ID_ANY, title, size=(400, 600), pos=pos)

        self.mainwin = mainwin
        self.plottype = plottype

        win = wx.Panel(self)  # ,wx.ID_ANY,pos=(0,0),size=(180,200),style=0)

        self.selection = wx.CheckListBox(win, wx.ID_ANY, choices=OutputPlotTypes, size=(120, -1))  # ,pos=(0,0))
        # self.win.Bind(wx.EVT_CHECKLISTBOX, self.OnTicked, self.selection)
        self.selection.Bind(wx.EVT_CHECKLISTBOX, self.OnTicked)

        if plottype == 'Theory':
            display_curves = self.mainwin.display_theory_curves
        else:
            display_curves = self.mainwin.display_expt_curves

        checked_items = []
        for i in range(len(display_curves)):
            if display_curves[i]:
                checked_items.append(i)

        self.selection.SetChecked(checked_items)
        # self.okbtn = wx.Button(self.win,wx.ID_OK,size=(120,BtnSize))
        # self.Bind(wx.EVT_BUTTON, self.OnOK, self.okbtn)

        self.SetSize(self.selection.GetSize() + (50, 50))
        self.SetMinSize(self.selection.GetSize() + (50, 50))

        popup_sizer = wx.BoxSizer(wx.VERTICAL)
        popup_sizer.Add(self.selection, 0, wx.EXPAND | wx.ALL, 7)
        # popup_sizer.Add((-1,5),1,wx.EXPAND)
        # popup_sizer.Add(self.okbtn,0,wx.EXPAND)

        # sz = popup_sizer.GetBestSize()
        # self.win.SetSize((sz.width+20, sz.height+20))

        self.SetSizer(popup_sizer)
        wx.CallAfter(self.Refresh)

    def OnTicked(self, event):
        """
        Action to perform when tick boxes are ticked

        Gets the tick box values, updates the main plot and changes
        the menu items to match the popup box
        """
        print('Ticked Event')

        if self.plottype == 'Theory':
            items = self.selection.GetChecked()
            self.mainwin.display_theory_curves = [False] * 9
            # for menuitem in self.mainwin.showTplotsSubMenu.GetMenuItems():
            #     menuitem.Check(False)
            for item in items:
                self.mainwin.display_theory_curves[item] = True
            # self.mainwin.showTplotsSubMenu.GetMenuItems()[item].Check(True)
            print(self.mainwin.display_theory_curves)
        elif self.plottype == 'Fit':
            items = self.selection.GetChecked()
            self.mainwin.display_expt_curves = [False] * 9
            # for menuitem in self.mainwin.showEplotsSubMenu.GetMenuItems():
            #    menuitem.Check(False)
            for item in items:
                self.mainwin.display_expt_curves[item] = True
            # self.mainwin.showEplotsSubMenu.GetMenuItems()[item].Check(True)
            print(self.mainwin.display_expt_curves)

        self.mainwin.OnCreateAxes(self.mainwin.figs[0])
        event.Skip()


class FittingThread(threading.Thread):
    """
    Fitting takes quite a bit of system resources, so we run the fit in another thread
    (and the processes spawned are low priority),
    leaving the GUI active in the meantime. When the fitting completes, this thread triggers an event in the
    GUI that then completes the rest of the fitting routine - i.e. updates the plot, writes text ... etc

    child class of the main threading.Thread class
    """

    def __init__(self, parent):
        threading.Thread.__init__(self)

        self.mainwin = parent

    # self.start() ## this calls the run() method - called from the GUI

    def run(self):
        """ Run fitting thread - called by <>.start() method """

        # mainwin is the top-level window (the ElecSus_GUI_Frame instance)
        mainwin = self.mainwin

        # crop x and y arrays sent to fit_data() method, based on fit_bounds (if it's not [None, None])
        if None in mainwin.detuning_fit_bounds:
            # use full arrays if fit bounds not specified (default)
            x_array, y_array = mainwin.x_fit_array, mainwin.y_fit_array
        else:
            # crop data to specified range

            # create aliases for ease
            fb = mainwin.detuning_fit_bounds
            xfa = mainwin.x_fit_array
            yfa = mainwin.y_fit_array
            #
            y_array = yfa[(xfa > fb[0]) & (xfa < fb[1])]
            x_array = xfa[(xfa > fb[0]) & (xfa < fb[1])]

            # fb, xfa, yfa = None, None, None

        print('Fitting data in the detuning range (GHz):  ', x_array[0], ' to ', x_array[-1])
        mainwin.FitInformation.write(
            'Fitting in the detuning range (GHz):  ' + str(x_array[0]) + ' to ' + str(x_array[-1]))
        mainwin.FitInformation.write('\n\n')

        # Check fit algorithms are ok
        if mainwin.fit_algorithm == 'Marquardt-Levenberg':
            fa = 'ML'
        elif mainwin.fit_algorithm == 'Random-Restart':
            fa = 'RR'
        elif mainwin.fit_algorithm == 'Simulated Annealing':
            fa = 'SA'
        elif mainwin.fit_algorithm == 'Differential Evolution':
            fa = 'DE'
        else:
            print("!! Fitting error - fit algorithm not defined")
            raise

        # rename data_type if S0:
        if mainwin.fit_datatype == 'Transmission (S0)':
            dt = 'S0'
        else:
            dt = mainwin.fit_datatype

        mainwin.opt_params, mainwin.rms, mainwin.fit_result = fit_data((x_array * 1e3, y_array),
                                                                       mainwin.params_dict,
                                                                       mainwin.params_dict_bools,
                                                                       p_dict_bounds=mainwin.params_dict_bounds,
                                                                       data_type=dt, fit_algorithm=fa)
        # **mainwin.advanced_fitoptions)

        # post an event to the main panel to run the OnFitComplete() method
        evt = FitCompleteEvent()
        wx.PostEvent(mainwin, evt)
    # alternately, could use wx.CallAfter here ...


class ProgressThread(threading.Thread):
    """ Update the progress bar continually to show user something is happening ..."""

    def __init__(self, parent):
        """ Create the thread """
        threading.Thread.__init__(self)
        self.mainwin = parent
        self.pb = parent.fitting_dlg  # .pb

    def run(self):
        """ Action to run when .start() method is called """
        while self.mainwin.already_fitting:
            for i in range(100):
                if not self.mainwin.already_fitting:
                    break
                # print i,' ',
                # wx.CallAfter(self.pb.SetValue,i)
                wx.CallAfter(self.pb.Update, i)
                time.sleep(0.1)
        print('Quitting progress bar update')


class OptionsPanel(scrolled.ScrolledPanel):
    def __init__(self, parent, mainwin, size, paneltype):
        """
        Panel which holds most of the control elements in this program.
        Most of the options are the same for theory and fitting, but there are some differences.
        Hence, the same panel class with an argument 'paneltype' which will set what is displayed.

        This ScrolledPanel class will automatically add horizontal and vertical scroll bars when required,
        and dynamically turns them on and off as window is resized. Neat!

        ##
        Note - this class has got a bit messy in it's old age, and could do with a rewrite to make it more clear!
        Another one for the To-Do list...
        ##

        # 'mainwin' is the main application frame so we can bind buttons to actions in the main frame
        # [ it's the parent (frame) of the parent (tab notebook) of the options panel (self) ]
        """

        self.mainwin = mainwin

        scrolled.ScrolledPanel.__init__(self, parent, size=size)

        self.paneltype = paneltype

        self.all_floatspin_inputs = []

        # Fitting only:
        if self.paneltype == 'Fit':
            # Import data from csv

            # Booleans (check boxes) for fitting
            self.main_paramlist_bools = \
                [wx.CheckBox(self, label='') for _ in main_parameterlist]
            self.magnet_paramlist_bools = \
                [wx.CheckBox(self, label='') for _ in magnet_parameterlist]

            for box in self.main_paramlist_bools + self.magnet_paramlist_bools:
                self.Bind(wx.EVT_CHECKBOX, self.OnFloatTicked, box)

            self.main_paramlist_usebounds = \
                [wx.CheckBox(self, label='') for _ in main_parameterlist]
            self.magnet_paramlist_usebounds = \
                [wx.CheckBox(self, label='') for _ in magnet_parameterlist]
            for boundsbox in self.main_paramlist_usebounds + self.magnet_paramlist_usebounds:
                self.Bind(wx.EVT_CHECKBOX, self.OnBoundsTicked, boundsbox)

            # Fitting algorithm selection
            # ties in with menu items
            self.fit_types = [wx.RadioButton(self, label="Marquardt-Levenberg", style=wx.RB_GROUP),
                              wx.RadioButton(self, label="Random-Restart"),
                              wx.RadioButton(self, label="Simulated Annealing"),
                              wx.RadioButton(self, label="Differential Evolution")]

            for fit_type, tt in zip(self.fit_types, fit_algorithm_tooltips):
                self.Bind(wx.EVT_RADIOBUTTON, mainwin.OnFitTypeChangePanel, fit_type)
                fit_type.SetToolTip(wx.ToolTip(tt))

            # Run Fit Button
            # RunFitButton = wx.Button(self,wx.ID_ANY, 'Run Fit', size=(140,1.5*BtnSize))
            # self.Bind(wx.EVT_BUTTON, mainwin.OnFitButton, RunFitButton)

            # Bounds - min/max boxes
            self.main_paramlist_mins = [FloatSpin(self, value=str(defval), increment=definc, size=(60, -1), digits=3)
                                        for defval, definc in
                                        zip(main_paramlist_mindefaults, defaultvals_main_increments)]
            self.main_paramlist_maxs = [FloatSpin(self, value=str(defval), increment=definc, size=(60, -1), digits=3)
                                        for defval, definc in
                                        zip(main_paramlist_maxdefaults, defaultvals_main_increments)]
            self.magnet_paramlist_mins = [FloatSpin(self, value=str(defval), increment=definc, size=(60, -1), digits=3)
                                          for defval, definc in
                                          zip(magnet_paramlist_mindefaults, defaultvals_magnet_increments)]
            self.magnet_paramlist_maxs = [FloatSpin(self, value=str(defval), increment=definc, size=(60, -1), digits=3)
                                          for defval, definc in
                                          zip(magnet_paramlist_maxdefaults, defaultvals_magnet_increments)]

            self.fit_polarisation_checkbox = wx.CheckBox(self, label='Fit Polarisation?')
            self.fit_polarisation_checkbox.SetToolTip(wx.ToolTip(fit_polarisation_tooltip))
            self.constrain_linear_checkbox = wx.CheckBox(self, label='Constrain Linear?')
            self.constrain_linear_checkbox.SetToolTip(wx.ToolTip(constrain_polarisation_tooltip))
            self.Bind(wx.EVT_CHECKBOX, self.OnConstrainLinear, self.constrain_linear_checkbox)

        # Theory only:
        elif self.paneltype == 'Theory':
            # Detuning range selection
            detuning_labels = [wx.StaticText(self, wx.ID_ANY, lab) for lab in
                               ["Start [GHz]", "Stop [GHz]", "No. of points"]]
            self.DetuningCtrl = [FloatSpin(self, value=str(defval), increment=definc, size=(70, -1), digits=3) for
                                 defval, definc in zip(detuning_defaults, detuning_increments)]
            for ctrl in detuning_labels:
                ctrl.SetToolTip(wx.ToolTip(parameter_adjust_tooltip))
            for ctrl in self.DetuningCtrl:
                self.all_floatspin_inputs.append(ctrl)

            self.main_paramlist_bools = [0] * len(main_parameterlist)
            self.main_paramlist_usebounds = [0] * len(main_parameterlist)
            self.main_paramlist_mins = [0] * len(main_parameterlist)
            self.main_paramlist_maxs = [0] * len(main_parameterlist)
            self.magnet_paramlist_bools = [0] * len(magnet_parameterlist)
            self.magnet_paramlist_usebounds = [0] * len(magnet_parameterlist)
            self.magnet_paramlist_mins = [0] * len(magnet_parameterlist)
            self.magnet_paramlist_maxs = [0] * len(magnet_parameterlist)

        # FIXED PARAMETERS
        fixed_paramlist_labels = [wx.StaticText(self, wx.ID_ANY, fixed_param) for
                                  fixed_param in fixed_parameterlist]
        # Add tooltips describing each parameter
        for label, tt in zip(fixed_paramlist_labels, fixed_parameter_tooltips):
            label.SetToolTip(wx.ToolTip(tt))
        # Add Combo Boxes + Check Box
        self.fixed_paramlist_inputs = [
            wx.ComboBox(self, wx.ID_ANY, choices=element_list, style=wx.CB_READONLY, size=(70, -1)),
            wx.ComboBox(self, wx.ID_ANY, choices=D_line_list, style=wx.CB_READONLY, size=(70, -1)),
            wx.CheckBox(self, label="")]
        self.fixed_paramlist_inputs[0].SetSelection(2)  # default selections
        self.fixed_paramlist_inputs[1].SetSelection(1)
        self.fixed_paramlist_inputs[2].SetValue(True)

        self.Bind(wx.EVT_COMBOBOX, self.OnElementSelect, self.fixed_paramlist_inputs[0])

        # MAIN PARAMETERS
        # create list of parameters - labels and spin-control boxes
        main_paramlist_labels = [wx.StaticText(self, wx.ID_ANY, fit_param + unit) for fit_param, unit in zip(
            main_parameterlist, mainunits_parameterlist)]
        # Tooltips
        for label, tt in zip(main_paramlist_labels, main_parameter_tooltips):
            label.SetToolTip(wx.ToolTip(tt + parameter_adjust_tooltip))
        # Input boxes
        self.main_paramlist_inputs = [
            FloatSpin(self, value=str(defval), increment=definc, size=(70, -1), digits=3) for defval, definc in
            zip(defaultvals_main_parameterlist, defaultvals_main_increments)]
        for item in self.main_paramlist_inputs:
            self.all_floatspin_inputs.append(item)
        # Don't need to bind SpinCtrl/ComboBox inputs to actions - the values are
        # read whenever 'Compute' or 'Fit' buttons are pressed

        # Layout panel in sizers

        # main sizer for the panel
        panel_sizer = wx.BoxSizer(wx.VERTICAL)

        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        label_sizer = wx.BoxSizer(wx.HORIZONTAL)
        label_sizer.Add((10, -1), 0, wx.EXPAND)
        elem_label = wx.StaticText(self, wx.ID_ANY, "Select element and D-line")
        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        elem_label.SetFont(font)
        label_sizer.Add(elem_label)
        label_sizer.Add((10, -1), 1, wx.EXPAND)
        panel_sizer.Add(label_sizer)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        # Add common elements first:

        for label, ipt in zip(fixed_paramlist_labels, self.fixed_paramlist_inputs):
            hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
            hor_sizer.Add((30, -1), 0, wx.EXPAND)
            hor_sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
            hor_sizer.Add((20, -1), 1, wx.EXPAND)
            hor_sizer.Add(ipt, 0, wx.ALIGN_CENTER_VERTICAL)
            hor_sizer.Add((10, -1), 0, wx.EXPAND)

            panel_sizer.Add(hor_sizer, 0, wx.EXPAND)
            panel_sizer.Add((-1, 3), 0, wx.EXPAND)

        # vertical space
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        label_sizer = wx.BoxSizer(wx.HORIZONTAL)
        label_sizer.Add((10, -1), 0, wx.EXPAND)
        parameter_label_text = 'Main parameters'
        mainparamlabeltextbox = wx.StaticText(self, wx.ID_ANY, parameter_label_text)
        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        mainparamlabeltextbox.SetFont(font)
        label_sizer.Add(mainparamlabeltextbox, 0, wx.EXPAND)
        label_sizer.Add((10, -1), 1, wx.EXPAND)
        if self.paneltype == 'Fit':
            floatlabel = wx.StaticText(self, wx.ID_ANY, 'Initial Value')
            label_sizer.Add(floatlabel, 0, wx.ALIGN_BOTTOM)
            label_sizer.Add((10, -1), 0, wx.EXPAND)
            floatlabel = wx.StaticText(self, wx.ID_ANY, 'Float?')
            label_sizer.Add(floatlabel, 0, wx.ALIGN_BOTTOM)
            label_sizer.Add((8, -1), 0, wx.EXPAND)
            boundslabel = wx.StaticText(self, wx.ID_ANY, 'Bounds?')
            label_sizer.Add(boundslabel, 0, wx.ALIGN_BOTTOM)
            label_sizer.Add((10, -1), 0, wx.EXPAND)
            minmaxlabel = wx.StaticText(self, wx.ID_ANY, 'Min Value')
            label_sizer.Add(minmaxlabel, 0, wx.ALIGN_BOTTOM)
            label_sizer.Add((13, -1), 0, wx.EXPAND)
            minmaxlabel = wx.StaticText(self, wx.ID_ANY, 'Max Value')
            label_sizer.Add(minmaxlabel, 0, wx.ALIGN_BOTTOM)
            label_sizer.Add((15, -1), 0, wx.EXPAND)
        panel_sizer.Add(label_sizer, 0, wx.EXPAND)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)

        for label, ipt, boolbox, boundsbox, minFS, maxFS in \
                zip(main_paramlist_labels, self.main_paramlist_inputs,
                    self.main_paramlist_bools, self.main_paramlist_usebounds,
                    self.main_paramlist_mins, self.main_paramlist_maxs):
            hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
            hor_sizer.Add((30, -1), 0, wx.EXPAND)
            hor_sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
            hor_sizer.Add((20, -1), 1, wx.EXPAND)
            hor_sizer.Add(ipt, 0, wx.ALIGN_CENTER_VERTICAL)
            if self.paneltype == 'Fit':
                hor_sizer.Add((20, -1), 0, wx.EXPAND)
                hor_sizer.Add(boolbox, 0, wx.ALIGN_CENTER_VERTICAL)
                hor_sizer.Add((20, -1), 0, wx.EXPAND)
                hor_sizer.Add(boundsbox, 0, wx.ALIGN_CENTER_VERTICAL)
                hor_sizer.Add((20, -1), 0, wx.EXPAND)
                hor_sizer.Add(minFS, 0, wx.ALIGN_CENTER_VERTICAL)
                hor_sizer.Add((10, -1), 0, wx.EXPAND)
                hor_sizer.Add(maxFS, 0, wx.ALIGN_CENTER_VERTICAL)
            hor_sizer.Add((10, -1), 0, wx.EXPAND)
            panel_sizer.Add(hor_sizer, 0, wx.EXPAND)
            panel_sizer.Add((-1, 3), 0, wx.EXPAND)

        # MAGNET PARAMETERS
        # create list of magnetic parameters - labels and spin-control boxes
        magnet_paramlist_labels = [wx.StaticText(self, wx.ID_ANY, param + unit) for param, unit in zip(
            magnet_parameterlist, magnetunits_parameterlist)]
        # Add Tooltips
        for label, tt in zip(magnet_paramlist_labels, magnet_parameter_tooltips):
            label.SetToolTip(wx.ToolTip(tt + parameter_adjust_tooltip))
        # Input boxes
        self.magnet_paramlist_inputs = [
            FloatSpin(self, value=str(defval), increment=definc, size=(70, -1), digits=3) for defval, definc in
            zip(defaultvals_magnet_parameterlist, defaultvals_magnet_increments)]
        for item in self.magnet_paramlist_inputs:
            self.all_floatspin_inputs.append(item)

        # Sizers
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        label_sizer = wx.BoxSizer(wx.HORIZONTAL)
        label_sizer.Add((10, -1), 0, wx.EXPAND)
        parameter_label_text = 'Magnet parameters'
        magparamlabeltextbox = wx.StaticText(self, wx.ID_ANY, parameter_label_text)
        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        magparamlabeltextbox.SetFont(font)
        label_sizer.Add(magparamlabeltextbox, 0, wx.EXPAND)
        label_sizer.Add((10, -1), 1, wx.EXPAND)
        panel_sizer.Add(label_sizer, 0, wx.EXPAND)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)

        for label, ipt, boolbox, boundsbox, minFS, maxFS in \
                zip(magnet_paramlist_labels, self.magnet_paramlist_inputs,
                    self.magnet_paramlist_bools, self.magnet_paramlist_usebounds,
                    self.magnet_paramlist_mins, self.magnet_paramlist_maxs):
            hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
            hor_sizer.Add((30, -1), 0, wx.EXPAND)
            hor_sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
            hor_sizer.Add((20, -1), 1, wx.EXPAND)
            hor_sizer.Add(ipt, 0, wx.ALIGN_CENTER_VERTICAL)
            if self.paneltype == 'Fit':
                hor_sizer.Add((20, -1), 0, wx.EXPAND)
                hor_sizer.Add(boolbox, 0, wx.ALIGN_CENTER_VERTICAL)
                hor_sizer.Add((20, -1), 0, wx.EXPAND)
                hor_sizer.Add(boundsbox, 0, wx.ALIGN_CENTER_VERTICAL)
                hor_sizer.Add((20, -1), 0, wx.EXPAND)
                hor_sizer.Add(minFS, 0, wx.ALIGN_CENTER_VERTICAL)
                hor_sizer.Add((10, -1), 0, wx.EXPAND)
                hor_sizer.Add(maxFS, 0, wx.ALIGN_CENTER_VERTICAL)

            hor_sizer.Add((10, -1), 0, wx.EXPAND)
            panel_sizer.Add(hor_sizer, 0, wx.EXPAND)
            panel_sizer.Add((-1, 3), 0, wx.EXPAND)

        # POLARISATION CONTROLS
        self.pol_radios = [wx.RadioButton(self, label="Linear", style=wx.RB_GROUP),
                           wx.RadioButton(self, label="Left CP"),
                           wx.RadioButton(self, label="Right CP"),
                           wx.RadioButton(self, label="Elliptical")
                           ]

        for pol_radio, tt in zip(self.pol_radios, polarisation_tooltips):
            self.Bind(wx.EVT_RADIOBUTTON, self.OnPolarisationTypeChange, pol_radio)
            pol_radio.SetToolTip(wx.ToolTip(tt))

        pol_control_labels = [wx.StaticText(self, wx.ID_ANY, pol_ctrl) for
                              pol_ctrl in polarisation_controls]
        for label, tt in zip(pol_control_labels, polarisation_parameter_tooltips):
            label.SetToolTip(wx.ToolTip(tt + parameter_adjust_tooltip))

        self.pol_control_inputs = [FloatSpin(self, value=str(defval), increment=definc, size=(70, -1), digits=3) for
                                   defval, definc in zip(defaultvals_pol_parameterlist, defaultvals_pol_increments)]
        for item in self.pol_control_inputs:
            self.all_floatspin_inputs.append(item)

        self.Bind(EVT_FLOATSPIN, self.OnTheta0, self.pol_control_inputs[0])
        self.Bind(EVT_FLOATSPIN, self.OnExCtrl, self.pol_control_inputs[1])
        self.Bind(EVT_FLOATSPIN, self.OnEyCtrl, self.pol_control_inputs[2])

        visualise_btn = wx.Button(self, label='Visualise Polarisation')
        self.Bind(wx.EVT_BUTTON, self.OnVisualisePol, visualise_btn)

        # initially, all but the first input are disabled
        for ctrl in self.pol_control_inputs[1:]:
            ctrl.Disable()

        if self.paneltype == 'Fit':
            pol_sizer_fit = wx.BoxSizer(wx.VERTICAL)
            pol_sizer_fit.Add((-1, 10), 1, wx.EXPAND)
            pol_sizer_fit.Add(self.fit_polarisation_checkbox, 0, wx.EXPAND)
            pol_sizer_fit.Add((-1, 10), 0, wx.EXPAND)
            pol_sizer_fit.Add(self.constrain_linear_checkbox, 0, wx.EXPAND)
            pol_sizer_fit.Add((-1, 10), 1, wx.EXPAND)

        # Sizers
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        label_sizer = wx.BoxSizer(wx.HORIZONTAL)
        label_sizer.Add((10, -1), 0, wx.EXPAND)
        parameter_label_text = 'Polarisation Parameters'
        polparamlabeltextbox = wx.StaticText(self, wx.ID_ANY, parameter_label_text)
        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        polparamlabeltextbox.SetFont(font)
        label_sizer.Add(polparamlabeltextbox, 0, wx.EXPAND)
        label_sizer.Add((10, -1), 1, wx.EXPAND)
        panel_sizer.Add(label_sizer, 0, wx.EXPAND)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)

        pol_sizer = wx.BoxSizer(wx.HORIZONTAL)
        pol_sizer_l = wx.BoxSizer(wx.VERTICAL)
        for radiobox in self.pol_radios:
            pol_sizer_l.Add(radiobox, 0, wx.EXPAND)
            pol_sizer_l.Add((-1, 5), 0, wx.EXPAND)
        pol_sizer_l.Add((-1, 5), 0, wx.EXPAND)
        pol_sizer_l.Add(visualise_btn, 0, wx.EXPAND)

        pol_sizer_r = wx.BoxSizer(wx.VERTICAL)
        for label, ipt in zip(pol_control_labels, self.pol_control_inputs):
            hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
            hor_sizer.Add((5, -1), 0, wx.EXPAND)
            hor_sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
            hor_sizer.Add((20, -1), 1, wx.EXPAND)
            hor_sizer.Add(ipt, 0, wx.ALIGN_CENTER_VERTICAL)

            hor_sizer.Add((10, -1), 0, wx.EXPAND)
            pol_sizer_r.Add(hor_sizer, 0, wx.EXPAND)
            pol_sizer_r.Add((-1, 3), 0, wx.EXPAND)

        pol_sizer.Add((30, -1), 0, wx.EXPAND)
        pol_sizer.Add(pol_sizer_l, 0, wx.EXPAND)
        pol_sizer.Add((5, -1), 1, wx.EXPAND)
        pol_sizer.Add(pol_sizer_r, 0, wx.EXPAND)
        if self.paneltype == 'Fit':
            pol_sizer.Add((5, -1), 0, wx.EXPAND)
            pol_sizer.Add(pol_sizer_fit, 0, wx.EXPAND)

        panel_sizer.Add(pol_sizer, 0, wx.EXPAND)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)

        if self.paneltype == 'Fit':
            # vertical space
            panel_sizer.Add((-1, 5), 0, wx.EXPAND)
            label_sizer = wx.BoxSizer(wx.HORIZONTAL)
            label_sizer.Add((10, -1), 0, wx.EXPAND)
            parameter_label_text = 'Fit Algorithm'
            paramlabeltextbox = wx.StaticText(self, wx.ID_ANY, parameter_label_text)
            font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
            paramlabeltextbox.SetFont(font)
            label_sizer.Add(paramlabeltextbox, 0, wx.EXPAND)
            label_sizer.Add((20, -1), 1, wx.EXPAND)
            panel_sizer.Add(label_sizer, 0, wx.EXPAND)
            panel_sizer.Add((-1, 10), 0, wx.EXPAND)

            for radiobtn in self.fit_types:
                hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
                hor_sizer.Add((30, -1), 0, wx.EXPAND)
                hor_sizer.Add(radiobtn, 0, wx.EXPAND)
                hor_sizer.Add((20, -1), 1, wx.EXPAND)
                panel_sizer.Add(hor_sizer, 0, wx.EXPAND)
                panel_sizer.Add((-1, 4), 0, wx.EXPAND)

        elif self.paneltype == 'Theory':
            panel_sizer.Add((-1, 5), 0, wx.EXPAND)
            label_sizer = wx.BoxSizer(wx.HORIZONTAL)
            label_sizer.Add((10, -1), 0, wx.EXPAND)
            parameter_label_text = 'Detuning Range'
            paramlabeltextbox = wx.StaticText(self, wx.ID_ANY, parameter_label_text)
            font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
            paramlabeltextbox.SetFont(font)
            label_sizer.Add(paramlabeltextbox, 0, wx.EXPAND)
            label_sizer.Add((20, -1), 1, wx.EXPAND)
            panel_sizer.Add(label_sizer, 0, wx.EXPAND)
            panel_sizer.Add((-1, 10), 0, wx.EXPAND)

            for label, ipt in zip(detuning_labels, self.DetuningCtrl):
                hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
                hor_sizer.Add((30, -1), 0, wx.EXPAND)
                hor_sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
                hor_sizer.Add((20, -1), 1, wx.EXPAND)
                hor_sizer.Add(ipt, 0, wx.ALIGN_CENTER_VERTICAL)
                hor_sizer.Add((10, -1), 0, wx.EXPAND)
                panel_sizer.Add(hor_sizer, 0, wx.EXPAND)
                panel_sizer.Add((-1, 3), 0, wx.EXPAND)

        panel_sizer.Add((-1, 5), 1, wx.EXPAND)
        panel_sizer.Add((-1, 2), 0, wx.EXPAND)

        # init the right combination of enabled boxes
        self.OnElementSelect()
        if self.paneltype == 'Fit':
            self.OnBoundsTicked()
            self.OnFloatTicked()

        self.SetSizer(panel_sizer)
        self.SetupScrolling()
        self.Layout()

    def OnConstrainLinear(self, event):
        if self.constrain_linear_checkbox.IsChecked():
            self.pol_control_inputs[3].SetValue(0)
        event.Skip()

    def OnElementSelect(self, event=None):
        """ Action when the element selection drop-down is clicked """
        selection = self.fixed_paramlist_inputs[0].GetValue()
        if selection in ('Cs', 'Na'):
            # disable Rb and K isotope selections
            for inputbox in self.main_paramlist_inputs[5:]:
                inputbox.Disable()
            if self.paneltype == 'Fit':
                for checkbox in self.main_paramlist_bools[5:]:
                    checkbox.SetValue(False)
                    checkbox.Disable()
        elif selection == 'Rb':
            # Enable Rb
            self.main_paramlist_inputs[5].Enable()
            if self.paneltype == 'Fit':
                self.main_paramlist_bools[5].Enable()
            # Disbale K
            for inputbox in self.main_paramlist_inputs[6:]:
                inputbox.Disable()
            if self.paneltype == 'Fit':
                for checkbox in self.main_paramlist_bools[6:]:
                    checkbox.Disable()
        elif selection == 'K':
            # Disabel Rb
            self.main_paramlist_inputs[5].Disable()
            if self.paneltype == 'Fit':
                self.main_paramlist_bools[5].Disable()
            # Disbale K
            for inputbox in self.main_paramlist_inputs[6:]:
                inputbox.Disable()
            if self.paneltype == 'Fit':
                for checkbox in self.main_paramlist_bools[6:]:
                    checkbox.Enable()
        if event is not None:
            event.Skip()

    def OnFloatTicked(self, event=None):
        """ Action when any of the 'float?' tick boxes are clicked """
        # Enable / Disable controls for varying parameters
        for i, box in enumerate(self.main_paramlist_bools):
            if not box.IsChecked():
                self.main_paramlist_usebounds[i].SetValue(False)
                self.OnBoundsTicked()
                self.main_paramlist_usebounds[i].Disable()
            else:
                self.main_paramlist_usebounds[i].Enable()

        for i, box in enumerate(self.magnet_paramlist_bools):
            if not box.IsChecked():
                self.magnet_paramlist_usebounds[i].SetValue(False)
                self.OnBoundsTicked()
                self.magnet_paramlist_usebounds[i].Disable()
            else:
                self.magnet_paramlist_usebounds[i].Enable()
        if event is not None:
            event.Skip()

    def OnBoundsTicked(self, event=None):
        """ Action when any of the fit bounds tick boxes are clicked """
        # enable / disable bounds floatspin controls if checkboxes are ticked
        for i, box in enumerate(self.main_paramlist_usebounds):
            if not box.IsChecked():
                self.main_paramlist_mins[i].Disable()
                self.main_paramlist_maxs[i].Disable()
            else:
                self.main_paramlist_mins[i].Enable()
                self.main_paramlist_maxs[i].Enable()

        for i, box in enumerate(self.magnet_paramlist_usebounds):
            if not box.IsChecked():
                self.magnet_paramlist_mins[i].Disable()
                self.magnet_paramlist_maxs[i].Disable()
            else:
                self.magnet_paramlist_mins[i].Enable()
                self.magnet_paramlist_maxs[i].Enable()
        if event is not None:
            event.Skip()

    def OnPolarisationTypeChange(self, event):
        """ Action when any of the polarisation radio buttons are clicked """

        # if Linear selected, make Ex,Ey,Phase inactive, set theta0 active and
        # set Phase to 0. On theta0 change, fill in Ex Ey as necessary
        if self.pol_radios[0].GetValue():
            for ctrl in self.pol_control_inputs[1:]:
                ctrl.Disable()
            self.pol_control_inputs[0].Enable()

            # set phase to 0 (linear)
            self.pol_control_inputs[3].SetValue(0)

        if self.pol_radios[1].GetValue():
            for ctrl in self.pol_control_inputs:
                ctrl.Disable()
            self.pol_control_inputs[1].SetValue(1. / np.sqrt(2))
            self.pol_control_inputs[2].SetValue(1. / np.sqrt(2))
            self.pol_control_inputs[3].SetValue(90)

        if self.pol_radios[2].GetValue():
            for ctrl in self.pol_control_inputs:
                ctrl.Disable()
            self.pol_control_inputs[1].SetValue(1. / np.sqrt(2))
            self.pol_control_inputs[2].SetValue(1. / np.sqrt(2))
            self.pol_control_inputs[3].SetValue(270)

        if self.pol_radios[3].GetValue():
            self.pol_control_inputs[0].Disable()
            for ctrl in self.pol_control_inputs[1:]:
                ctrl.Enable()
        event.Skip()

    def OnTheta0(self, event):
        """ Update values of Ex,Ey when theta0 is changed """
        self.pol_control_inputs[1].SetValue(np.cos(np.pi / 180 * self.pol_control_inputs[0].GetValue()))
        self.pol_control_inputs[2].SetValue(np.sin(np.pi / 180 * self.pol_control_inputs[0].GetValue()))
        event.Skip()

    def OnExCtrl(self, event):
        if self.pol_control_inputs[1].GetValue() <= 1:
            self.pol_control_inputs[2].SetValue(np.sqrt(1 - self.pol_control_inputs[1].GetValue() ** 2))
        else:
            scaling = np.sqrt(self.pol_control_inputs[1].GetValue() ** 2 + self.pol_control_inputs[2].GetValue() ** 2)
            self.pol_control_inputs[2].SetValue(self.pol_control_inputs[2].GetValue() / scaling)
            self.pol_control_inputs[1].SetValue(self.pol_control_inputs[1].GetValue() / scaling)

        self.pol_control_inputs[0].SetValue(
            180. / np.pi * np.arctan(self.pol_control_inputs[2].GetValue() / self.pol_control_inputs[1].GetValue()))
        event.Skip()

    def OnEyCtrl(self, event):
        if self.pol_control_inputs[2].GetValue() <= 1:
            self.pol_control_inputs[1].SetValue(np.sqrt(1 - self.pol_control_inputs[2].GetValue() ** 2))
        else:
            scaling = np.sqrt(self.pol_control_inputs[1].GetValue() ** 2 + self.pol_control_inputs[2].GetValue() ** 2)
            self.pol_control_inputs[2].SetValue(self.pol_control_inputs[2].GetValue() / scaling)
            self.pol_control_inputs[1].SetValue(self.pol_control_inputs[1].GetValue() / scaling)

        self.pol_control_inputs[0].SetValue(
            180. / np.pi * np.arctan(self.pol_control_inputs[2].GetValue() / self.pol_control_inputs[1].GetValue()))
        event.Skip()

    def OnVisualisePol(self, event):
        """ Run visualise polarisation animation in 3D plot"""
        e_x, e_y, phase = \
            self.pol_control_inputs[1].GetValue(), \
            self.pol_control_inputs[2].GetValue(), \
            self.pol_control_inputs[3].GetValue() * np.pi / 180

        plt.close(self.mainwin.figs[0])
        plt.close(self.mainwin.figs[1])

        pol_ani.animate_vectors(e_x, e_y, phase)
        event.Skip()


class PlotToolPanel(wx.Panel):
    """ Panel to hold the matplotlib figure/canvas and toolbar """

    def __init__(self, parent, mainwin, wx_id):
        """ mainwin is the main panel so we can bind buttons to actions in the main frame """
        wx.Panel.__init__(self, parent, id=-1)

        self.fig = plt.figure(wx_id, facecolor=(240. / 255, 240. / 255, 240. / 255), figsize=(12.9, 9.75), dpi=80)

        # self.ax = self.fig.add_subplot(111)

        # create the wx objects to hold the figure
        self.canvas = FigureCanvasWxAgg(self, -1, self.fig)
        self.toolbar = Toolbar(self.canvas)  # matplotlib toolbar (pan, zoom, save etc)
        # self.toolbar.Realize()

        # Create vertical sizer to hold figure and toolbar - dynamically expand with window size
        plot_sizer = wx.BoxSizer(wx.VERTICAL)
        plot_sizer.Add(self.canvas, 1, flag=wx.EXPAND | wx.ALL)  # wx.TOP|wx.LEFT|wx.GROW)
        plot_sizer.Add(self.toolbar, 0, wx.EXPAND)

        mainwin.figs.append(self.fig)
        mainwin.fig_IDs.append(wx_id)  # use an ID number to keep track of figures
        mainwin.canvases.append(self.canvas)

        # display some text in the middle of the window to begin with
        self.fig.text(0.5, 0.5,
                      'ElecSus GUI\n\nVersion ' + __version__ +
                      '\n\nTo get started, use the panel on the right\nto either Compute a spectrum or '
                      'Import some data...',
                      ha='center', va='center')
        # self.fig.hold(False)

        self.SetSizer(plot_sizer)


# self.Layout() #Fit()


class StatusPanel(scrolled.ScrolledPanel):
    """ Panel to hold status information (contents of the tab) """

    def __init__(self, parent, wx_id):
        """ mainwin is the main panel so we can bind buttons to actions in the main frame """
        scrolled.ScrolledPanel.__init__(self, parent)

        status_title = wx.StaticText(self, wx.ID_ANY, wx_id)
        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        status_title.SetFont(font)

        self.StatusTextBox = wx.TextCtrl(self, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_MULTILINE)
        # self.StatusTextBox.Size.SetHeight(500)

        self.SaveBtn = wx.Button(self, wx.ID_ANY, "Save Text to file", size=(100, -1))
        self.Bind(wx.EVT_BUTTON, self.OnSave, self.SaveBtn)

        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        panel_sizer.Add(status_title, 0, wx.EXPAND | wx.LEFT, border=20)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        panel_sizer.Add(self.StatusTextBox, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, border=40)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        panel_sizer.Add(self.SaveBtn, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, border=40)
        panel_sizer.Add((-1, 20), 0, wx.EXPAND)

        self.SetSizer(panel_sizer)
        self.SetupScrolling()
        self.Layout()

    def OnSave(self, event):
        """
        Save text file dump of all data in the panel - select filename,
        check for overwrite and then pass to save_data() method
        """
        save_file_dialog = wx.FileDialog(self, "Save Output File", "./", "Outputs",
                                         "Text files (*.txt)|*.txt", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if save_file_dialog.ShowModal() == wx.ID_OK:

            output_filename = save_file_dialog.GetPath()
            print(output_filename)
            # if output_filename[-4:] == exts[save_file_dialog.GetFilterIndex()]:
            #    output_filename = output_filename[:-4]
            save_file_dialog.Destroy()

            # check for overwrite current files
            if os.path.isfile(output_filename):
                overwrite_dialog = wx.MessageDialog(self, "Warning: file exists already! Overwrite?",
                                                    "Overwrite?", wx.YES_NO | wx.NO_DEFAULT)

                if overwrite_dialog.ShowModal() == wx.NO:
                    overwrite_dialog.Destroy()
                    return  # exit without saving
                else:
                    overwrite_dialog.Destroy()

            # do save
            self.save_data(output_filename)
        event.Skip()

    def save_data(self, filename):
        """ Save the data using wx.TextCtrl built-in method """
        success = self.StatusTextBox.SaveFile(filename)  # returns true if no errors
        if not success:
            problem_dlg = wx.MessageDialog(self, "There was an error saving the data...", "Error saving",
                                           wx.OK | wx.ICON_ERROR)
            problem_dlg.ShowModal()

    def write(self, textstring):
        """ Append text to the text box """
        self.StatusTextBox.AppendText(textstring)


class DataProcessingDlg(wx.Dialog):
    def __init__(self, parent, title, wx_id):
        """ Dialog box for smoothing and/or binning experimental data """

        wx.Dialog.__init__(self, parent, wx_id, title, size=(360, 450))

        self.SetMinSize((360, 450))

        self.parent = parent

        # make copy of non-binned / smoothed data in case reset button is pressed
        self.original_xdata = parent.x_fit_array
        self.original_ydata = parent.y_expt_arrays

        print('Length before binning: ', len(self.parent.x_fit_array))

        panel_blurb = wx.StaticText(self, wx.ID_ANY, "Data smoothing and binning options.", size=(350, -1),
                                    style=wx.ALIGN_CENTRE_HORIZONTAL)
        panel_blurb.Wrap(350)

        bin_blurb = wx.StaticText(self, wx.ID_ANY,
                                  "Binning data is useful where the initial data size is very large. "
                                  "The output from this function is smaller by a specified factor n, "
                                  "which also removes some noise since every n data points are averaged "
                                  "into one.\n\n For fitting data, the number of experimental data points "
                                  "is a factor in how long the fit takes to run - making the array smaller"
                                  " makes the fit quicker.",
                                  size=(350, -1), style=wx.ALIGN_CENTRE_HORIZONTAL)
        bin_blurb.Wrap(350)

        self.cur_dsize_label = wx.StaticText(self, wx.ID_ANY, "Current data size:")
        self.cur_dsize = wx.TextCtrl(self, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_CENTRE, size=(80, -1))
        self.cur_dsize.ChangeValue(str(len(parent.x_fit_array)))

        self.bin_size_label = wx.StaticText(self, wx.ID_ANY, "Bin size:")
        self.bin_size = wx.TextCtrl(self, wx.ID_ANY, "11", style=wx.TE_CENTRE, size=(80, -1))
        self.Bind(wx.EVT_TEXT, self.OnBinSizeChange, self.bin_size)

        self.new_dsize_label = wx.StaticText(self, wx.ID_ANY, "Data size after binning:")
        self.new_dsize = wx.TextCtrl(self, wx.ID_ANY, "", style=wx.TE_READONLY | wx.TE_CENTRE, size=(80, -1))
        self.new_dsize.ChangeValue(str(int(len(parent.x_fit_array) / float(self.bin_size.GetValue()))))

        bin_btn = wx.Button(self, label="Bin Data", size=(75, -1))
        self.Bind(wx.EVT_BUTTON, self.OnDoBin, bin_btn)

        smooth_blurb = wx.StaticText(self, wx.ID_ANY,
                                     "Smoothing is based on a moving average with a triangular "
                                     "weighting function of a specified width. For convenience, "
                                     "the output array is the same length as the input, but the "
                                     "smoothing is only applied over N - 2n data points (removing "
                                     "n from each end of the array), where N is the length of the data"
                                     " array and n is the width of the smoothing window.",
                                     size=(350, -1), style=wx.ALIGN_CENTRE_HORIZONTAL)
        smooth_blurb.Wrap(350)

        self.smooth_size_label = wx.StaticText(self, wx.ID_ANY, "Moving average window size:")
        self.smooth_size = wx.TextCtrl(self, wx.ID_ANY, "11", style=wx.TE_CENTRE, size=(80, -1))

        smooth_btn = wx.Button(self, label="Smooth Data")
        self.Bind(wx.EVT_BUTTON, self.OnDoSmooth, smooth_btn)

        reset_btn = wx.Button(self, label='Reset to original data')
        self.Bind(wx.EVT_BUTTON, self.OnReset, reset_btn)

        close_btn = wx.Button(self, wx.ID_OK, label='Close')

        # layout dialog box
        panel_sizer = wx.BoxSizer(wx.VERTICAL)

        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        panel_sizer.Add(panel_blurb, 0, wx.ALIGN_CENTRE_HORIZONTAL)

        bin_sizer = wx.BoxSizer(wx.HORIZONTAL)
        bin_sizer.Add((20, -1), 0, wx.EXPAND)
        bin_sizer.Add(self.bin_size_label, 0, wx.EXPAND)
        bin_sizer.Add((10, -1), 1, wx.EXPAND)
        bin_sizer.Add(self.bin_size, 0, wx.EXPAND)
        bin_sizer.Add((10, -1), 0, wx.EXPAND)
        bin_sizer.Add(bin_btn, 0, wx.EXPAND)
        bin_sizer.Add((10, -1), 0, wx.EXPAND)

        cdata_sizer = wx.BoxSizer(wx.HORIZONTAL)
        cdata_sizer.Add((20, -1), 0, wx.EXPAND)
        cdata_sizer.Add(self.cur_dsize_label, 0, wx.EXPAND)
        cdata_sizer.Add((10, -1), 1, wx.EXPAND)
        cdata_sizer.Add(self.cur_dsize, 0, wx.EXPAND)
        cdata_sizer.Add((95, -1), 0, wx.EXPAND)

        ndata_sizer = wx.BoxSizer(wx.HORIZONTAL)
        ndata_sizer.Add((20, -1), 0, wx.EXPAND)
        ndata_sizer.Add(self.new_dsize_label, 0, wx.EXPAND)
        ndata_sizer.Add((10, -1), 1, wx.EXPAND)
        ndata_sizer.Add(self.new_dsize, 0, wx.EXPAND)
        ndata_sizer.Add((95, -1), 0, wx.EXPAND)

        smooth_sizer = wx.BoxSizer(wx.HORIZONTAL)
        smooth_sizer.Add((20, -1), 0, wx.EXPAND)
        smooth_sizer.Add(self.smooth_size_label, 0, wx.EXPAND)
        smooth_sizer.Add((10, -1), 1, wx.EXPAND)
        smooth_sizer.Add(self.smooth_size, 0, wx.EXPAND)
        smooth_sizer.Add((10, -1), 0, wx.EXPAND)
        smooth_sizer.Add(smooth_btn, 0, wx.EXPAND)
        smooth_sizer.Add((10, -1), 0, wx.EXPAND)

        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_sizer.Add((20, -1), 1, wx.EXPAND)
        btn_sizer.Add(reset_btn, 0, wx.EXPAND)
        btn_sizer.Add((20, -1), 0, wx.EXPAND)
        btn_sizer.Add(close_btn, 0, wx.EXPAND)
        btn_sizer.Add((20, -1), 1, wx.EXPAND)

        panel_sizer.Add((-1, 20), 0, wx.EXPAND)
        panel_sizer.Add(bin_blurb, 0, wx.EXPAND)
        panel_sizer.Add((-1, 20), 0, wx.EXPAND)
        panel_sizer.Add(cdata_sizer, 0, wx.EXPAND)
        panel_sizer.Add(bin_sizer, 0, wx.EXPAND)
        panel_sizer.Add(ndata_sizer, 0, wx.EXPAND)
        panel_sizer.Add((-1, 20), 0, wx.EXPAND)
        panel_sizer.Add(wx.StaticLine(self, -1, size=(-1, 1), style=wx.LI_HORIZONTAL),
                        0, wx.EXPAND | wx.LEFT | wx.RIGHT, border=20)
        panel_sizer.Add((-1, 20), 0, wx.EXPAND)
        panel_sizer.Add(smooth_blurb, 0, wx.EXPAND)
        panel_sizer.Add((-1, 20), 0, wx.EXPAND)
        panel_sizer.Add(smooth_sizer, 0, wx.EXPAND)
        panel_sizer.Add((-1, 20), 1, wx.EXPAND)
        panel_sizer.Add(btn_sizer, 0, wx.EXPAND)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)

        self.SetSizerAndFit(panel_sizer)

        self.Layout()

    def OnBinSizeChange(self, event):
        """ update ui elements for value of new data length """
        # try:
        self.cur_dsize.ChangeValue(str(int(len(self.parent.x_fit_array))))
        self.new_dsize.ChangeValue(str(int(len(self.parent.x_fit_array) / float(self.bin_size.GetValue()))))
        event.Skip()

    def OnDoBin(self, event):
        """ Action when bin button is pressed """
        bin_amnt = int(self.bin_size.GetValue())

        # only bin/smooth the most recently loaded data set
        self.parent.x_fit_array, self.parent.y_fit_array, ye = data_proc.bin_data(self.parent.x_fit_array,
                                                                                  self.parent.y_fit_array, bin_amnt)

        print('Length ater binning: ', len(self.parent.x_fit_array))

        # re-plot data
        self.update_plot()

        # update ui elements with new data size...
        self.OnBinSizeChange(1)
        event.Skip()

    def OnDoSmooth(self, event):
        """ Action when smooth button is pressed """
        smth_amnt = int(self.smooth_size.GetValue())

        self.parent.y_fit_array = data_proc.smooth_data(self.parent.y_fit_array, smth_amnt)

        # re-plot data
        self.update_plot()
        event.Skip()

    def OnReset(self, event):
        """ Reset to the original data, i.e. the data present when the dialog box was created """
        self.parent.x_fit_array = self.original_xdata
        self.parent.y_fit_array = self.original_ydata

        # re-plot data
        self.update_plot()
        event.Skip()

    def update_plot(self):
        """ Common method for updating the plot """
        # update plot arrays
        cindex = self.parent.choice_index
        self.parent.x_expt_arrays[cindex] = self.parent.x_fit_array
        self.parent.y_expt_arrays[cindex] = self.parent.y_fit_array

        self.parent.OnCreateAxes(self.parent.figs[0], clear_current=True)


class AdvancedFitOptions(wx.Dialog):
    def __init__(self, parent, title, wx_id):
        """ Dialog box for selecting advanced fit options which are passed through to lmfit optimisation routines """

        wx.Dialog.__init__(self, parent, wx_id, title, size=(500, 900))

        # main sizer for this dialog box
        panel_sizer = wx.BoxSizer(wx.VERTICAL)

        panel_blurb = wx.StaticText(self, wx.ID_ANY,
                                    "These are the options given to lmfit when fitting data. "
                                    "\n\n BLURB TO UPDATE LATER",
                                    size=(460, -1), style=wx.ALIGN_CENTRE_HORIZONTAL)
        panel_blurb.Wrap(460)

        # UPDATE FOR DIFFERENTIAL EVOLUTION + BOUNDS ON FIT PARAMETERS ####

        self.fitopt_labels = ['Max Iterations (maxfev)', 'Max Tolerance [1e-8] (ftol)']
        self.fitopt_argnames = ['maxfev', 'ftol']
        fitopt_defaults = [1000, 15]
        fitopt_increments = [100, 1]
        self.fitopt_scaling = [1, 1e-9]
        fitopt_texts = [wx.StaticText(self, wx.ID_ANY, label) for label in self.fitopt_labels]
        self.fitopt_ctrl = [wx.SpinCtrl(self, value=str(defval), size=(80, -1), min=0, max=10000, initial=defval) for
                            defval, definc in zip(fitopt_defaults, fitopt_increments)]

        panel_sizer.Add((-1, 10), 0, wx.EXPAND)
        panel_sizer.Add(panel_blurb, 0, wx.LEFT | wx.RIGHT, border=20)
        panel_sizer.Add((-1, 15), 0, wx.EXPAND)
        panel_sizer.Add(wx.StaticLine(self, -1, size=(-1, 1), style=wx.LI_HORIZONTAL), 0,
                        wx.EXPAND | wx.LEFT | wx.RIGHT, border=20
                        )
        panel_sizer.Add((-1, 15), 0, wx.EXPAND)

        for static, ctrl in zip(fitopt_texts, self.fitopt_ctrl):
            hor_sizer = wx.BoxSizer(wx.HORIZONTAL)
            hor_sizer.Add(static, 0, wx.EXPAND | wx.LEFT, border=20)
            hor_sizer.Add((10, -1), 1, wx.EXPAND)
            hor_sizer.Add(ctrl, 0, wx.EXPAND | wx.RIGHT, border=20)
            panel_sizer.Add(hor_sizer, 0, wx.EXPAND)
            panel_sizer.Add((-1, 5), 0, wx.EXPAND)

        # ok and cancel buttons
        btnbar = self.CreateButtonSizer(wx.OK | wx.CANCEL)

        panel_sizer.Add((-1, 10), 1, wx.EXPAND)
        panel_sizer.Add(btnbar, 0, wx.ALIGN_CENTER)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)

        self.SetSizer(panel_sizer)
        self.Layout()

    def return_all_options(self):
        """ get list of all ctrl elements """
        opt_dict = dict([(label, ipt.GetValue() * scaling) for label, ipt, scaling
                         in zip(self.fitopt_argnames, self.fitopt_ctrl, self.fitopt_scaling)])
        # print 'Advanced fit options dictionary: ',opt_dict
        return opt_dict


class FitBoundsDialog(wx.Dialog):
    def __init__(self, parent, title, wx_id):
        """
        Dialog box for selecting advanced fit options which are passed through to
        curve_fit / leastsq optimisation routines...
        """

        wx.Dialog.__init__(self, parent, wx_id, title, size=(300, 300))

        self.parent = parent

        # main sizer for this dialog box
        panel_sizer = wx.BoxSizer(wx.VERTICAL)

        panel_blurb = wx.StaticText(self, wx.ID_ANY,
                                    "Set bounds to fit data. By default, all data is fitted, but "
                                    "here the fitted data can be limited to a certain detuning range",
                                    size=(260, -1), style=wx.ALIGN_CENTRE_HORIZONTAL)
        panel_blurb.Wrap(260)

        self.detuning_fit_bounds = [wx.RadioButton(self, label="Use Full Data Range", style=wx.RB_GROUP),
                                    wx.RadioButton(self, label="Cropped Data Range")]
        for btn in self.detuning_fit_bounds:
            self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelector, btn)

        self.fit_range_ctrl = [wx.TextCtrl(self, wx.ID_ANY, "-5", style=wx.TE_CENTRE, size=(80, -1)),
                               wx.TextCtrl(self, wx.ID_ANY, "5", style=wx.TE_CENTRE, size=(80, -1))]
        self.OnRadioSelector(1)

        fit_range_labels = [wx.StaticText(self, wx.ID_ANY, "Min: "), wx.StaticText(self, wx.ID_ANY, "Max: ")]

        panel_sizer.Add((-1, 20), 0, wx.EXPAND)
        panel_sizer.Add(panel_blurb, 0, wx.ALIGN_CENTRE_HORIZONTAL)

        panel_sizer.Add((-1, 20), 0, wx.EXPAND)
        panel_sizer.Add(wx.StaticLine(self, -1, size=(-1, 1), style=wx.LI_HORIZONTAL),
                        0, wx.EXPAND | wx.LEFT | wx.RIGHT, border=20)
        panel_sizer.Add((-1, 20), 0, wx.EXPAND)

        panel_sizer.Add(self.detuning_fit_bounds[0], 0, wx.EXPAND | wx.LEFT, border=30)
        panel_sizer.Add((-1, 5), 0, wx.EXPAND)
        panel_sizer.Add(self.detuning_fit_bounds[1], 0, wx.EXPAND | wx.LEFT, border=30)

        panel_sizer.Add((-1, 20), 0, wx.EXPAND)

        for i in range(2):
            ctrl_sizer = wx.BoxSizer(wx.HORIZONTAL)
            ctrl_sizer.Add((45, -1), 1, wx.EXPAND)
            ctrl_sizer.Add(fit_range_labels[i], 0, wx.EXPAND)
            ctrl_sizer.Add((15, -1), 0, wx.EXPAND)
            ctrl_sizer.Add(self.fit_range_ctrl[i], 0, wx.EXPAND)
            ctrl_sizer.Add((45, -1), 1, wx.EXPAND)

            panel_sizer.Add(ctrl_sizer, 0, wx.EXPAND)
            panel_sizer.Add((-1, 5), 0, wx.EXPAND)

        panel_sizer.Add((-1, 15), 1, wx.EXPAND)

        btnbar = self.CreateButtonSizer(wx.OK | wx.CANCEL)
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_sizer.Add((20, -1), 1, wx.EXPAND)
        btn_sizer.Add(btnbar, 0, wx.EXPAND)
        btn_sizer.Add((20, -1), 1, wx.EXPAND)

        panel_sizer.Add(btn_sizer, 0, wx.EXPAND)
        panel_sizer.Add((-1, 10), 0, wx.EXPAND)

        self.SetSizer(panel_sizer)
        self.Layout()

    def OnRadioSelector(self, event):
        # enable / disbale the text input box
        enabled = self.detuning_fit_bounds[1].GetValue()
        for ctrl in self.fit_range_ctrl:
            ctrl.Enable(enabled)
        event.Skip()

    def OnUpdateRange(self):
        if self.detuning_fit_bounds[0].GetValue():
            # if using full data range
            try:
                return self.parent.x_fit_array[0], self.parent.x_fit_array[-1]
            except:  # no data loaded TODO: define exception type!
                return [None, None]
        else:
            return [float(self.fit_range_ctrl[0].GetValue()), float(self.fit_range_ctrl[1].GetValue())]


class ElecSusGuiFrame(wx.Frame):
    """ Main class for this program - top-level window """
    # init class attributes
    rms = 0
    fit_result = None
    opt_params = None
    x_array = None
    y_optimised = None

    params_dict = {}
    params_dict_bools = {}
    params_dict_bounds = {}

    main_param_labels = ['T', 'lcell', 'shift', 'GammaBuf', 'DoppTemp',
                         'rb85frac', 'K40frac', 'K41frac']
    magnet_param_labels = ['Bfield', 'Btheta', 'Bphi']
    polarisation_param_labels = ['E_x', 'E_y', 'E_phase']

    main_params = None
    magnet_params = None
    polarisation_params = None

    already_fitting = False
    fit_bools = None
    detuning_fit_bounds = None
    fit_algorithm = "ML"
    fit_warnings = ""
    x_fit_array = None
    y_fit_array = None
    fit_datatype = None

    xlim = None
    ylim = None

    autoscale_axes = True
    live_events_bound = False
    legend_on = False

    expt_type = None
    choice_index = 0

    filename = ''
    dirname = ''

    residual_histogram = False
    use_exp_detuning = True

    def __init__(self, parent, title):
        """ Initialise main frame """
        wx.Frame.__init__(self, parent, title=title, size=(2000, 900))

        # ubuntu sizing:
        if os.name == 'posix':
            font = wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL)
            self.SetFont(font)

        # Set icons for top-left of frame, alt-tab window ...
        frame_icon = wx.IconBundle()
        try:
            frame_icon.AddIconFromFile(os.path.join(elecsus_dir, 'images/elecsus_t_group.ico'), wx.BITMAP_TYPE_ANY)
        except:  # TODO: define exception type!
            # new wx version
            frame_icon.AddIcon(os.path.join(elecsus_dir, 'images/elecsus_t_group.ico'), wx.BITMAP_TYPE_ANY)

        self.SetIcons(frame_icon)

        # if the window is closed, exit
        self.Bind(wx.EVT_CLOSE, self.OnExit)

        self.panel = wx.Panel(self)

        self.panel.SetBackgroundColour(wx.Colour(240, 240, 240))

        self._init_default_values()
        self._init_plot_defaults()
        self._init_menus()
        self._init_panels()

        # redirect stdout (command line text) to status box
        sys.stdout = self.StatusPanel
        sys.stderr = self.ErrorPanel

        # Create initially blank set of axes
        # self.OnCreateAxes(self.figs[0],self.canvases[0])

        # Bind the event EVT_FIT_COMPLETE to function
        # This executes in the main thread once the fitting thread
        # (separate from the main thread) completes
        self.Bind(EVT_FIT_COMPLETE, self.OnFitCompleted)

    def _init_default_values(self):
        """ Initialise default values for various things ... """

        self.figs = []
        self.canvases = []
        self.fig_IDs = []

        self.hold = False

        self.fit_algorithm = 'Marquardt-Levenberg'
        self.already_fitting = False
        self.fit_warnings = True

        # initialise advanced fit options dictionary
        # dlg = AdvancedFitOptions(self, "Advanced Fit Options", wx.ID_ANY)
        # self.advanced_fitoptions = dlg.return_all_options()
        # dlg.Destroy()

    def _init_plot_defaults(self):
        """
        List of default values for all plots.
        Theory plots can be up to 8 panels, each of which can be turned on/off
        """

        self.plot_outputs = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix, Iy', 'I+45 / I-45', 'Ircp / Ilcp',
                             'Alpha Plus/Minus/Z']
        self.plot_ylabels = ['Transmission, $S_0$', '$S_1$', '$S_2$', '$S_3$', '$I_x$, $I_y$', '$I_{+45}$, $I_{-45}$',
                             r'$I_{RCP}, I_{LCP}$', r'$\alpha^{+,-,Z}$']
        self.plot_output_indices = [OutputPlotTypes.index(po) for po in self.plot_outputs]

        self.xrange = [detuning_defaults[0], detuning_defaults[1]]  # detuning range, in GHz
        self.npoints = detuning_defaults[2]  # number of detuning points to calculate

        # default data for plots - blank lists
        self.x_array = np.linspace(self.xrange[0], self.xrange[1], self.npoints)
        self.y_arrays = [None] * len(self.plot_outputs)
        self.x_expt_arrays = [None] * len(self.plot_outputs)
        self.y_expt_arrays = [None] * len(self.plot_outputs)
        self.x_fit_array, self.y_fit_array = [], []

        # x and y-limits
        self.xlim = None
        self.ylim = [None] * len(OutputPlotTypes)

        # plot data visible when set to True - default = view transmission (S0)
        self.display_theory_curves = [False] * len(OutputPlotTypes)
        self.display_theory_curves[0] = True
        self.display_expt_curves = [False] * len(OutputPlotTypes)

        # live plotting
        self.live_events_bound = False

        # residual plot settings
        self.normalised_residuals = False
        self.residual_histogram = False

        # fit bounds
        self.detuning_fit_bounds = [None, None]

        # Theory curves are calculated with same detuning axis as experimental data
        self.use_exp_detuning = False

        # autoscale
        self.autoscale_axes = True

    def _init_menus(self):
        """ Initialise menu bar items """

        # Create menu_bar object
        menu_bar = wx.MenuBar()

        # File
        file_menu = wx.Menu()
        fm_open = file_menu.Append(wx.ID_OPEN, "&Open Experimental Data (csv)\tCtrl+O",
                                   "Open file for plotting and/or fitting.")
        self.Bind(wx.EVT_MENU, self.OnFileOpen, fm_open)
        file_menu.AppendSeparator()
        fm_saveplot = file_menu.Append(wx.ID_ANY, "Save Plot as Image", "Save Data")
        self.Bind(wx.EVT_MENU, self.OnSaveFig, fm_saveplot)
        fm_savecsv = file_menu.Append(wx.ID_SAVE, "E&xport CSV Data\tCtrl+S", "Save CSV Data")
        self.Bind(wx.EVT_MENU, self.OnSaveCSVData, fm_savecsv)
        # fm_saveconfig = file_menu.Append(wx.ID_ANY, "Save Current Configuration", "Save Config")
        # self.Bind(wx.EVT_MENU, self.OnSaveConfig, fm_saveconfig)
        file_menu.AppendSeparator()
        fm_exit = file_menu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", "Close window and exit program.")
        self.Bind(wx.EVT_MENU, self.OnExit, fm_exit)
        #

        # Edit
        edit_menu = wx.Menu()
        em_copy_t_to_f = edit_menu.Append(wx.ID_ANY, "Copy Parameters: &Theory to Fit", "Copy - T to F")
        em_copy_f_to_t = edit_menu.Append(wx.ID_ANY, "Copy Parameters: &Fit to Theory", "Copy - F to T")
        self.Bind(wx.EVT_MENU, self.OnCopyParamsTtoF, em_copy_t_to_f)
        self.Bind(wx.EVT_MENU, self.OnCopyParamsFtoT, em_copy_f_to_t)

        self.em_use_exp_detuning = edit_menu.Append(wx.ID_ANY, "&Use Experimental Detuning Axis", "Use Exp Axis",
                                                    kind=wx.ITEM_CHECK)
        self.Bind(wx.EVT_MENU, self.OnUseExpDetuning, self.em_use_exp_detuning)

        # About
        about_menu = wx.Menu()
        am_aboutthis = about_menu.Append(wx.ID_ABOUT, "&About this program", "About this program.")
        self.Bind(wx.EVT_MENU, self.OnAboutThis, am_aboutthis)
        am_aboutelecsus = about_menu.Append(wx.ID_ANY, "About &ElecSus", "About ElecSus")
        self.Bind(wx.EVT_MENU, self.OnAboutElecSus, am_aboutelecsus)

        # View
        view_menu = wx.Menu()
        vm_liveplot = view_menu.Append(wx.ID_ANY, "&Live Plot", "Update Plot in real-time", kind=wx.ITEM_CHECK)
        self.Bind(wx.EVT_MENU, self.OnLivePlotting, vm_liveplot)
        vm_autoscale = view_menu.Append(wx.ID_ANY, "&Autoscale when updating axes", "Autoscale", kind=wx.ITEM_CHECK)
        vm_autoscale.Check(True)
        self.Bind(wx.EVT_MENU, self.OnAutoscale, vm_autoscale)

        # Select plot types to display

        # Theory Plot
        theoryplot_menu = wx.Menu()
        tpm_plotholdon = theoryplot_menu.Append(wx.ID_ANY,
                                                "&Hold data on plot update",
                                                "Select whether to hold or clear current plot data on updating"
                                                " the figure",
                                                kind=wx.ITEM_CHECK)
        tpm_plotholdon.Check(self.hold)
        self.Bind(wx.EVT_MENU, self.OnPlotHold, tpm_plotholdon)
        # tpM_clearplot = theoryplot_menu.Append(wx.ID_ANY, "&Clear current plot data",
        #                                        "Clear current plot data on all axes")
        # self.Bind(wx.EVT_MENU, self.OnClearPlot, tpM_clearplot)

        tpm_grid = theoryplot_menu.Append(wx.ID_ANY, "&Grid on axes", "Grid", kind=wx.ITEM_CHECK)
        self.Bind(wx.EVT_MENU, self.OnGridToggleMain, tpm_grid)

        self.show_e_plots_submenu = wx.Menu()
        id_s0 = 3000
        id_s1 = 3001
        id_s2 = 3002
        id_s3 = 3003
        id_ix_iy = 3004
        id_alpha = 3005
        id_n = 3006
        id_phi = 3007
        show_s0 = self.show_e_plots_submenu.AppendCheckItem(id_s0, "&Transmission (S0)")
        show_s1 = self.show_e_plots_submenu.AppendCheckItem(id_s1, "S&1")
        show_s2 = self.show_e_plots_submenu.AppendCheckItem(id_s2, "S&2")
        show_s3 = self.show_e_plots_submenu.AppendCheckItem(id_s3, "S&3")
        show_ix_iy = self.show_e_plots_submenu.AppendCheckItem(id_ix_iy, "&Ix/Iy")
        show_alpha = self.show_e_plots_submenu.AppendCheckItem(id_alpha, "&Alpha +/-")
        show_n = self.show_e_plots_submenu.AppendCheckItem(id_n, "&Refractive Index +/-")
        show_ng = self.show_e_plots_submenu.AppendCheckItem(id_n, "&Group Index +/-")
        show_phi = self.show_e_plots_submenu.AppendCheckItem(id_phi, "&Phi (Rotation Angle)")
        # bind event to check box selections
        for checkitem in [show_s0, show_s1, show_s2, show_s3, show_ix_iy, show_alpha, show_n, show_ng, show_phi]:
            self.Bind(wx.EVT_MENU, self.OnShowEplots, checkitem)
        # add to parent menu
        tpm_show_e_plots = theoryplot_menu.AppendMenu(wx.ID_ANY, "&Experimental Curves to Show",
                                                      self.show_e_plots_submenu)

        # Experimental Plot
        exptplot_menu = wx.Menu()
        epm_plotholdon = exptplot_menu.Append(wx.ID_ANY, "&Hold data on plot update",
                                              "Select whether to hold or clear current plot data "
                                              "on updating the figure",
                                              kind=wx.ITEM_CHECK)
        # epM_clearplot = exptplot_menu.Append(wx.ID_ANY, "&Clear current plot data",
        #                                      "Clear current plot data on all axes")
        # self.Bind(wx.EVT_MENU, self.OnClearPlot,epM_clearplot)

        epm_grid = exptplot_menu.Append(wx.ID_ANY, "&Grid on axes", "Grid", kind=wx.ITEM_CHECK)
        self.Bind(wx.EVT_MENU, self.OnGridToggleRes, epm_grid)

        epm_reshist = exptplot_menu.Append(wx.ID_ANY, "Show &Histogram of Residuals", "Res Hist", kind=wx.ITEM_CHECK)
        self.Bind(wx.EVT_MENU, self.OnResHist, epm_reshist)

        # Fit
        fit_menu = wx.Menu()
        # fitM_dummy = fit_menu.Append(wx.ID_ANY, "&Dummy", "Dummy")

        fitm_dataproc = fit_menu.Append(wx.ID_ANY, "&Data Processing...", "Data Processing Options")
        self.Bind(wx.EVT_MENU, self.OnDataProcessing, fitm_dataproc)

        fitm_fitbounds = fit_menu.Append(wx.ID_ANY, "&Set Fit Bounds", "Set detuning bounds for fitting")
        self.Bind(wx.EVT_MENU, self.OnFitBounds, fitm_fitbounds)

        # fitM_advanced = fit_menu.Append(wx.ID_ANY, "&Advanced Fit Options...", "Advanced Fit Options")
        # self.Bind(wx.EVT_MENU, self.OnAdvancedOptions, fitM_advanced)

        fitm_warnings = fit_menu.Append(wx.ID_ANY, "&Warn about fit settings",
                                        "Warn about possible bad fit settings",
                                        kind=wx.ITEM_CHECK)

        self.Bind(wx.EVT_MENU, self.OnFitWarnings, fitm_warnings)
        fitm_warnings.Check(True)

        #
        # ... other menu items to add...
        #
        # Residual analysis - fit gaussian to residuals etc
        #       -- Use normalised/raw residuals (requires importing data with errorbars.. to do..)
        # Main plot - set axes limits
        #

        # Add menu items to menu bar
        menu_bar.Append(file_menu, "&File")
        menu_bar.Append(edit_menu, "&Edit")
        menu_bar.Append(view_menu, "&View")
        menu_bar.Append(theoryplot_menu, "&Main Plot")
        menu_bar.Append(exptplot_menu, "&Residuals Plot")
        menu_bar.Append(fit_menu, "F&it")
        menu_bar.Append(about_menu, "&About")

        # Add Menu Bar to the main panel
        self.SetMenuBar(menu_bar)

    def _init_panels(self):
        """
        Initialise panel with matplotlib window, buttons, text boxes etc.
        Doesn't really need to be in a separate function, but makes it easier to find stuff...
        """

        # Sizer constructs are used for placement of all GUI elements.
        # This makes the whole window scalable in a predictable way.

        #
        # Create plot part of the window
        #

        # create plot in a notebook-style for adding more tabs later on
        plot_tabs = wx.Notebook(self.panel)
        # plot_tabs = aui.AuiNotebook(self.panel)

        # The plot panels
        self.T_Panel = PlotToolPanel(plot_tabs, self, 1)
        self.E_Panel = PlotToolPanel(plot_tabs, self, 2)

        # Text tabs - for stdout, stderr messages and fitting information
        # generated after fitting data (optimum parameters etc)
        self.StatusPanel = StatusPanel(plot_tabs, 'Status Information')
        self.FitInformation = StatusPanel(plot_tabs, 'Fitting Information')
        self.ErrorPanel = StatusPanel(plot_tabs, 'Error Information')
        self.ErrorPanel.write("If this is the only thing displayed here, the program is working well...\n\n")

        # Add all tabs to the tab bar
        plot_tabs.AddPage(self.T_Panel, "Main Plot")
        plot_tabs.AddPage(self.E_Panel, "Fit / Residuals Plot")
        plot_tabs.AddPage(self.StatusPanel, "Status Panel")
        plot_tabs.AddPage(self.FitInformation, "Fitting Information")
        plot_tabs.AddPage(self.ErrorPanel, "Error Information")

        # Add the Tab bar to the main panel sizer
        plot_sizer = wx.BoxSizer(wx.VERTICAL)
        plot_sizer.Add(plot_tabs, 1, wx.EXPAND)

        #
        # Create button part of the window
        #

        # elecsus and JQC logos at the top
        # elecsuslogo = wx.Image('images/elecsus.ico',wx.BITMAP_TYPE_ANY)
        # elecsuslogo.Rescale(108/2,138/2)
        # elecsus_bmp = wx.StaticBitmap(self.panel,wx.ID_ANY,wx.BitmapFromImage(elecsuslogo),size=(108/2,-1))
        jqclogo = wx.Image(os.path.join(elecsus_dir, 'images/jqc-logo.png'), wx.BITMAP_TYPE_ANY)
        # jqc_bmp = wx.StaticText(self.panel,wx.ID_ANY,"image")
        jqc_bmp = wx.StaticBitmap(self.panel, wx.ID_ANY, wx.BitmapFromImage(jqclogo), size=(191, -1))

        spacer = 167
        image_sizer = wx.BoxSizer(wx.HORIZONTAL)
        image_sizer.Add((10, -1), 1, wx.EXPAND)
        image_sizer.Add((spacer, -1), 0, wx.EXPAND)
        # image_sizer.Add(elecsus_bmp,0,wx.EXPAND)
        # image_sizer.Add((10,-1),0,wx.EXPAND)
        image_sizer.Add(jqc_bmp, 0, wx.EXPAND)
        image_sizer.Add((spacer, -1), 0, wx.EXPAND)
        image_sizer.Add((10, -1), 1, wx.EXPAND)

        # Menus are tabbed into theory and fitting sections
        # use wx.Notebook for this
        tab_panel = wx.Notebook(self.panel)

        # add tabs
        self.thy_panel = wx.Panel(tab_panel)
        self.fit_panel = wx.Panel(tab_panel)
        self.thy_options = OptionsPanel(self.thy_panel, self, (600, -1), 'Theory')
        self.fit_options = OptionsPanel(self.fit_panel, self, (600, -1), 'Fit')

        # Run Fit Button
        run_fit_button = wx.Button(self.fit_panel, wx.ID_ANY, 'Run Fit', size=(140, 1.5 * BtnSize))
        self.Bind(wx.EVT_BUTTON, self.OnFitButton, run_fit_button)
        # Calculate spectrum
        compute_button = wx.Button(self.thy_panel, wx.ID_ANY, 'Compute Spectrum', size=(140, 1.5 * BtnSize))
        self.Bind(wx.EVT_BUTTON, self.OnComputeButton, compute_button)
        # Common options to both experiment and theory
        import_button_thy = wx.Button(self.thy_panel, wx.ID_OPEN, label="Import Data", size=(140, 1.5 * BtnSize))
        self.Bind(wx.EVT_BUTTON, self.OnFileOpen, import_button_thy)
        import_button_fit = wx.Button(self.fit_panel, wx.ID_OPEN, label="Import Data", size=(140, 1.5 * BtnSize))
        self.Bind(wx.EVT_BUTTON, self.OnFileOpen, import_button_fit)

        thy_btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        thy_btn_sizer.Add((30, -1), 1, wx.EXPAND)
        thy_btn_sizer.Add(import_button_thy, 0, wx.EXPAND)
        thy_btn_sizer.Add((30, -1), 0, wx.EXPAND)
        thy_btn_sizer.Add(compute_button, 0, wx.EXPAND)
        thy_btn_sizer.Add((30, -1), 1, wx.EXPAND)

        fit_btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        fit_btn_sizer.Add((30, -1), 1, wx.EXPAND)
        fit_btn_sizer.Add(import_button_fit, 0, wx.EXPAND)
        fit_btn_sizer.Add((30, -1), 0, wx.EXPAND)
        fit_btn_sizer.Add(run_fit_button, 0, wx.EXPAND)
        fit_btn_sizer.Add((30, -1), 1, wx.EXPAND)

        self.thy_sizer = wx.BoxSizer(wx.VERTICAL)
        self.thy_sizer.Add(self.thy_options, 1, wx.EXPAND)
        self.thy_sizer.Add((-1, 10), 0, wx.EXPAND)
        self.thy_sizer.Add(thy_btn_sizer, 0, wx.EXPAND)

        self.fit_sizer = wx.BoxSizer(wx.VERTICAL)
        self.fit_sizer.Add(self.fit_options, 1, wx.EXPAND)
        self.fit_sizer.Add((-1, 10), 0, wx.EXPAND)
        self.fit_sizer.Add(fit_btn_sizer, 0, wx.EXPAND)

        self.thy_panel.SetSizer(self.thy_sizer)
        self.fit_panel.SetSizer(self.fit_sizer)

        tab_panel.AddPage(self.thy_panel, "Theory Settings")
        tab_panel.AddPage(self.fit_panel, "Fit Settings")

        tab_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # tab_sizer.Add((40,-1),0,wx.EXPAND)
        tab_sizer.Add(tab_panel, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, border=5)
        # tab_sizer.Add((40,-1),0,wx.EXPAND)

        # Plot selection buttons
        plot_selection_text = wx.StaticText(self.panel, wx.ID_ANY, "Data types to show:")
        t_plot_selection_list = wx.Button(self.panel, wx.ID_ANY, 'Theory Plots', size=(100, BtnSize))
        self.Bind(wx.EVT_BUTTON, self.OnTheoryPlotSelection, t_plot_selection_list)
        e_plot_selection_list = wx.Button(self.panel, wx.ID_ANY, 'Fit Plots', size=(100, BtnSize))
        self.Bind(wx.EVT_BUTTON, self.OnFitPlotSelection, e_plot_selection_list)

        # for checkitem in [show_S0, show_S1, show_S2, show_S3, show_IxIy, show_alpha, show_n]:
        #   self.Bind(wx.EVT_MENU, self.OnShowTplots, checkitem)
        # Add to parent menu
        # tpM_showTplots = theoryplotMenu.AppendMenu(wx.ID_ANY,"&Theory Curves to Show", self.showTplotsSubMenu)

        plot_selection_sizer = wx.BoxSizer(wx.HORIZONTAL)
        plot_selection_sizer.Add((10, -1), 0, wx.EXPAND)
        plot_selection_sizer.Add(plot_selection_text, 0, wx.ALIGN_CENTER_VERTICAL)
        plot_selection_sizer.Add((10, -1), 1, wx.EXPAND)
        plot_selection_sizer.Add(t_plot_selection_list, 0, wx.EXPAND)
        plot_selection_sizer.Add((20, -1), 0, wx.EXPAND)
        plot_selection_sizer.Add(e_plot_selection_list, 0, wx.EXPAND)
        plot_selection_sizer.Add((20, -1), 1, wx.EXPAND)

        # Create sizer for right-hand side of panel
        button_sizer = wx.BoxSizer(wx.VERTICAL)
        button_sizer.Add((-1, 10), 0, wx.EXPAND)
        button_sizer.Add(image_sizer, 0, wx.EXPAND)
        button_sizer.Add((-1, 5), 0, wx.EXPAND)
        button_sizer.Add(plot_selection_sizer, 0, wx.EXPAND)  # ,0,wx.EXPAND)
        button_sizer.Add((-1, 5), 0, wx.EXPAND)
        button_sizer.Add(tab_sizer, 1, wx.EXPAND)
        button_sizer.Add((-1, 5), 0, wx.EXPAND)

        # Put plot part and button part together
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)
        main_sizer.Add(plot_sizer, 1, wx.EXPAND)
        main_sizer.Add(wx.StaticLine(self.panel, -1, size=(1, -1), style=wx.LI_VERTICAL), 0, wx.EXPAND)
        main_sizer.Add(button_sizer, 0, wx.EXPAND)

        self.panel.SetSizerAndFit(main_sizer)
        self.panel.Layout()
        self.SetSize((1300, 900))
        self.SetSize((1300, 850))

    # wx.CallAfter(self.canvases[0].SetSize,(1200,900))

    #
    # Actions for events
    #

    #
    # General Actions - Close window, Save Data etc...
    #

    def _draw_fig(self, fig):
        """ shortcut method for redrawing figure and rearranging subplots to fill space """
        try:
            fig.tight_layout()
        except:  # TODO: define exception type!
            pass

        for can in self.canvases:
            can.draw()

    # self.SendSizeEvent()

    def call_elecsus(self, calc_or_fit):
        """ Call elecsus with passed arguments - single calculation or fit, depending on option """

        # get parameter list, xrange from either theory or fit panel
        panel = None
        if calc_or_fit == 'Theory':
            panel = self.thy_options
            if not self.use_exp_detuning:
                # get detuning range from GUI controls
                xmin, xmax, npts = [ipt.GetValue() for ipt in panel.DetuningCtrl]
                self.x_array = np.linspace(xmin, xmax, int(npts))
            else:
                # get detuning from experimental data
                self.x_array = self.x_fit_array
        elif calc_or_fit == 'Fit':
            panel = self.fit_options

        # Dictionary of parameters, bools (floating or not) and bounds ((min,max) pairs)
        self.params_dict = {}
        self.params_dict_bools = {}
        self.params_dict_bounds = {}

        # get parameter list
        const_params = panel.fixed_paramlist_inputs
        self.params_dict['Elem'] = const_params[0].GetStringSelection()
        self.params_dict['Dline'] = const_params[1].GetStringSelection()
        self.params_dict['Constrain'] = const_params[2].IsChecked()

        # strings of dictionary entries as used by elecsus
        self.main_params = [ipt.GetValue() for ipt in panel.main_paramlist_inputs]
        for label, value in zip(self.main_param_labels, self.main_params):
            self.params_dict[label] = value

        # modfiy cell length - convert from mm to m
        self.params_dict['lcell'] *= 1e-3

        self.magnet_params = [ipt.GetValue() for ipt in panel.magnet_paramlist_inputs]
        for label, value in zip(self.magnet_param_labels, self.magnet_params):
            self.params_dict[label] = value
        # convert from degree to radians
        self.params_dict['Btheta'] *= np.pi / 180
        self.params_dict['Bphi'] *= np.pi / 180

        self.polarisation_params = [ipt.GetValue() for ipt in
                                    panel.pol_control_inputs[1:]]  # theta-0 is not needed in general
        for label, value in zip(self.polarisation_param_labels, self.polarisation_params):
            self.params_dict[label] = value

        # get Efield vector
        self.params_dict['E_phase'] *= np.pi / 180
        e_vec = np.array(
            [self.params_dict['E_x'], self.params_dict['E_y'] * np.exp(1.j * self.params_dict['E_phase']), 0])
        print(self.params_dict)

        if calc_or_fit == 'Fit':
            # Append floating parameters to params_dict_bools
            for i, checkbox in enumerate(panel.main_paramlist_bools):
                if checkbox.IsChecked():
                    self.params_dict_bools[self.main_param_labels[i]] = True

            for i, checkbox in enumerate(panel.magnet_paramlist_bools):
                if checkbox.IsChecked():
                    self.params_dict_bools[self.magnet_param_labels[i]] = True

            # Append parameter bounds to params_dict_bounds
            for i, checkbox in enumerate(panel.main_paramlist_usebounds):
                if checkbox.IsChecked():
                    minn = self.fit_options.main_paramlist_mins[i].GetValue()
                    maxx = self.fit_options.main_paramlist_maxs[i].GetValue()
                    self.params_dict_bounds[self.main_param_labels[i]] = (minn, maxx)

            for i, checkbox in enumerate(panel.magnet_paramlist_usebounds):
                if checkbox.IsChecked():
                    minn = self.fit_options.magnet_paramlist_mins[i].GetValue()
                    maxx = self.fit_options.magnet_paramlist_maxs[i].GetValue()
                    self.params_dict_bounds[self.magnet_param_labels[i]] = (minn, maxx)

            # POLARISATION FIT OPTIONS TO ADD

            if self.fit_options.fit_polarisation_checkbox.IsChecked():
                print('Fitting Polarisation...')

                self.params_dict_bools['E_x'] = True
                self.params_dict_bools['E_y'] = True

                self.params_dict_bounds['E_x'] = [-1., 1.]
                self.params_dict_bounds['E_y'] = [-1., 1.]

                if not self.fit_options.constrain_linear_checkbox.IsChecked():
                    print('Fitting Pol. Phase...')
                    self.params_dict_bools['E_phase'] = True
                #   self.params_dict_bounds['E_phase'] = (0.,np.pi)
                else:
                    self.params_dict_bools['E_phase'] = False

            print(self.params_dict_bools)

        if calc_or_fit == 'Theory':
            print('\n\n')
            print(time.ctime())
            print('Calling ElecSus for single calculation with parameters:')
            print(self.params_dict)

            spectrum_data = calculate(self.x_array * 1e3, e_vec, self.params_dict,
                                      outputs=['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45', 'Ir', 'Il',
                                               'alphaPlus', 'alphaMinus', 'alphaZ'])
            self.y_arrays[0:4] = spectrum_data[0:4]  # S0,S1,S2,S3
            self.y_arrays[4] = [spectrum_data[4], spectrum_data[5]]  # Ix,Iy
            self.y_arrays[5] = [spectrum_data[6], spectrum_data[7]]  # I+45,I-45
            self.y_arrays[6] = [spectrum_data[8], spectrum_data[9]]  # I_RCP, I_LCP
            self.y_arrays[7] = [spectrum_data[10], spectrum_data[11], spectrum_data[12]]  # alpha +/-/z

            # print self.y_arrays
            # for i, array in enumerate(self.y_arrays):
            #   print i, type(array)

            self.get_current_axes_limits(self.figs[0])
            self.OnCreateAxes(self.figs[0], clear_current=True)

        elif calc_or_fit == 'Fit':
            # Use another thread for running elecsus fitting so the main panel doesn't stop responding

            # Check if fit is already running, only proceed if not fitting already
            if not self.already_fitting:

                # Verbose output
                print('\n\n')
                print(time.ctime())
                print('Calling ElecSus for fitting data with initial parameters:')
                for key in self.params_dict:
                    print(key, self.params_dict[key])
                print('\nVarying the following parameters:')
                for key in self.params_dict_bools:
                    print(key)
                print('\nSubject to the following boundaries:')
                for key in self.params_dict_bounds:
                    val = self.params_dict_bounds[key]
                    print(key, val[0], 'to ', val[1])

                # log time, data set to be fitted, initial parameters on the 'Fit Info' notebook:
                font1 = wx.Font(10, wx.TELETYPE, wx.NORMAL, wx.NORMAL)  # , False, u'Consolas')
                self.FitInformation.StatusTextBox.SetFont(font1)
                self.FitInformation.write('Fit started at:'.ljust(25) + time.ctime() + '\n')
                self.FitInformation.write(
                    'Data set to be fitted:'.ljust(25) + os.path.join(self.dirname, self.filename) + '\n')
                self.FitInformation.write('Experimental Data Type:'.ljust(25) + self.fit_datatype + '\n')
                self.FitInformation.write('Initial parameters (Floated):\n')
                for key in self.params_dict:
                    val = self.params_dict[key]
                    self.FitInformation.write('\t' + key.ljust(30) + str(val).ljust(10))
                    if key in self.params_dict_bools:
                        self.FitInformation.write(' (Floating')
                        if key in self.params_dict_bounds:
                            self.FitInformation.write(
                                ', bounded between ' + str(self.params_dict_bounds[key][0]) + ' and' + str(
                                    self.params_dict_bounds[key][1]) + ' )')
                        else:
                            self.FitInformation.write(', unbounded)')

                self.FitInformation.write('Using algorithm: ' + self.fit_algorithm)

                ##
                # Run Fitting algorithm
                ##

                # initialise thread for fitting
                fit_thread = FittingThread(self)

                # start thread - thread grabs all information from main window, so start() is all that's needed
                fit_thread.start()
                # After the fit completes, the thread triggers an EVT_FIT_COMPLETE which
                # calls (in the main thread) the OnFitCompleted() method

                # set flag to only allow one fit to run at any time
                self.already_fitting = True

                # show dialog box with constantly cycling status bar
                self.fitting_dlg = wx.ProgressDialog("Fit in progress",
                                                "Fitting in progress. Please wait for fitting to complete. "
                                                "This may take a while.",
                                                maximum=100, style=wx.PD_ELAPSED_TIME)
                self.fitting_dlg.Show(True)
                fit_progress_thread = ProgressThread(self)
                fit_progress_thread.start()

            else:
                dlg = wx.MessageDialog(self, "Fit already in progress. Only one fit may be run at a time.",
                                       "Patience required...", style=wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()

    def copy_params(self, order):
        """ Copy current parameter values, either between theory and fit tabs, or other way around """
        f_params = \
            self.fit_options.fixed_paramlist_inputs + self.fit_options.main_paramlist_inputs \
            + self.fit_options.magnet_paramlist_inputs + self.fit_options.pol_control_inputs

        t_params = \
            self.thy_options.fixed_paramlist_inputs + self.thy_options.main_paramlist_inputs \
            + self.thy_options.magnet_paramlist_inputs + self.thy_options.pol_control_inputs

        for fp, tp in zip(f_params, t_params):
            if order == 1:
                # Theory to Fit
                fp.SetValue(tp.GetValue())
            elif order == 2:
                tp.SetValue(fp.GetValue())

        print('Parameters Copied')

    def get_current_axes_limits(self, fig):
        """ Get current axes limits (doesn't return, but updates self.<> variables"""
        displayed_axes = [i | j for i, j in
                          zip(self.display_theory_curves, self.display_expt_curves)]
        # n_axes = displayed_axes.count(True)

        # get current axes limits
        if len(fig.axes) > 0:
            self.xlim = fig.axes[0].get_xlim()
            j = 0
            for i, displayed_ax in enumerate(displayed_axes):
                if displayed_ax:
                    self.ylim[i] = fig.axes[j].get_ylim()
                    j += 1

    def OnAboutElecSus(self, event):
        """ Show a message box about ElecSus """
        dlg = wx.MessageDialog(self, NOTICE.noticestring, "About", wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            dlg.Destroy()
        event.Skip()

    def OnAboutThis(self, event):
        """ Show a message box about the program """
        dlg = wx.MessageDialog(self,
                               "ElecSus GUI\n\nA graphical interface for ElecSus, a program to calculate"
                               " the weak-probe electric susceptibility of alkali atoms.\n\n "
                               "Written in python using wxPython and matplotlib.\n\nCopyright "
                               "2016-21 James Keaveney and co-authors",
                               "About", wx.OK)
        if dlg.ShowModal() == wx.ID_OK:
            dlg.Destroy()
        event.Skip()

    def OnAdvancedOptions(self, event):
        """ Open the dialog box for adjusting advanced fit options """
        dlg = AdvancedFitOptions(self, "Advanced Fit Options", wx.ID_ANY)

        # Show() rather than ShowModal() - doesn't halt program flow
        if dlg.ShowModal() == wx.ID_OK:
            advanced_fitoptions = dlg.return_all_options()
            print(advanced_fitoptions)
        event.Skip()

    def OnAutoscale(self, event):
        """ Set axes auto-scaling """
        self.autoscale_axes = bool(event.IsChecked())

    def OnComputeButton(self, event):
        """ Call elecsus to compute spectrum """
        self.call_elecsus('Theory')
        event.Skip()

    def OnCopyParamsTtoF(self, event):
        """ Copy parameters from theory tab to fitting tab """
        self.copy_params(1)
        event.Skip()

    def OnCopyParamsFtoT(self, event):
        """ Copy parameters from fitting tab to theory tab """
        self.copy_params(2)
        event.Skip()

    def OnCreateAxes(self, fig, clear_current=True):
        """ (Re)Create as many sets of axes in the main figure as are needed, and label them all """

        # calculate how many axes should be displayed - count if any elements of display_expt/theory_curves are true
        displayed_axes = [i | j for i, j in zip(self.display_theory_curves, self.display_expt_curves)]
        n_axes = displayed_axes.count(True)

        # clear figure and start again from nothing if number of axes has changed, or if hold is off
        if clear_current:
            fig.clf()
        if n_axes != len(fig.axes):
            fig.clf()

        # create bare axes, all equal sizes
        if n_axes == 0:
            n_axes = 1

        fig.add_subplot(n_axes, 1, 1)
        for i in range(1, n_axes):
            fig.add_subplot(n_axes, 1, i + 1, sharex=fig.axes[0])

        # Testing:
        # print self.display_expt_curves
        # print self.y_expt_arrays

        # PLOT BASED ON BOOL DISPLAY_x_CURVE
        i = 0
        for displayT, displayE, yt, xe, ye, ylabel in zip(self.display_theory_curves, self.display_expt_curves,
                                                          self.y_arrays, self.x_expt_arrays, self.y_expt_arrays,
                                                          self.plot_ylabels):
            try:
                # print i
                ax = fig.axes[i]
            except:  # TODO: define exception type (index error?)
                ax = None
            if displayE and ye is not None:
                # if isinstance(ye, (list,tuple)):
                #   for xi,yi in zip(xe,ye):
                #       ax.plot(xi,yi,color=d_grey)
                # else:
                ax.plot(xe, ye, color=d_olive)
            if displayT and yt is not None:
                if isinstance(yt, (list, tuple)):
                    for yi in yt:
                        ax.plot(self.x_array, yi)
                else:
                    ax.plot(self.x_array, yt)
            if displayT or displayE:
                ax.set_ylabel(ylabel)
                i += 1

        # set x axis label and rescale all axes to fit data
        fig.axes[-1].set_xlabel('Detuning (GHz)')
        if self.autoscale_axes:
            for ax in fig.axes:
                ax.autoscale_view(tight=True)
        else:
            fig.axes[0].set_xlim(self.xlim)
            j = 0
            for i, displayed_ax in enumerate(displayed_axes):
                if displayed_ax:
                    fig.axes[j].set_ylim(self.ylim[i])
                    j += 1

        # fig.axes[-1].set_xlim(self.xrange)

        # remove the rest of the x tick labels from all but the bottom panel
        for ax in fig.axes[:-1]:
            plt.setp(ax.get_xticklabels(), visible=False)

        # print 'Created '+str(n_axes)+' axes'

        # update the plot window
        self._draw_fig(fig)

    def on_create_residual_plot(self, fig):
        """
        Create a plot with main data, residuals and (optionally) histogram of residuals, using
        subplot2grid in matplotlib
        """

        fig.clf()

        print('Debugging...')
        print(self.fit_datatype)
        print(self.y_optimised)

        # normalised_residuals = False
        if self.normalised_residuals:
            # not done yet! -- requires error bars in imported data
            residuals = 100 * (self.y_fit_array - self.y_optimised)
        else:
            residuals = 100 * (self.y_fit_array - self.y_optimised)

        fig = plt.figure(2)
        yy = 4
        xx = 6
        if self.residual_histogram:
            ax_main = plt.subplot2grid((yy, xx), (0, 0), colspan=xx - 1, rowspan=yy - 1)
            ax_residual = plt.subplot2grid((yy, xx), (yy - 1, 0), colspan=xx - 1, sharex=ax_main)
            ax_hist = plt.subplot2grid((yy, xx), (yy - 1, xx - 1), sharey=ax_residual)

            plt.setp(ax_hist.get_yticklabels(), visible=False)
            ax_hist.set_xticklabels([])

        else:
            ax_main = plt.subplot2grid((yy, xx), (0, 0), colspan=xx, rowspan=yy - 1)
            ax_residual = plt.subplot2grid((yy, xx), (yy - 1, 0), colspan=xx, sharex=ax_main)

        plt.setp(ax_main.get_xticklabels(), visible=False)

        ax_residual.set_xlabel('Detuning (GHz)')
        ax_residual.set_ylabel('Residuals (%)')

        ax_main.set_ylabel(self.expt_type)

        ax_main.plot(self.x_fit_array, self.y_fit_array, color=d_olive)
        print(len(self.x_fit_array), len(self.y_optimised))
        ax_main.plot(self.x_fit_array, self.y_optimised)
        ax_residual.plot(self.x_fit_array, residuals, lw=1.25)
        ax_residual.axhline(0, color='k', linestyle='dashed')

        if self.residual_histogram:
            bins = 25
            ax_hist.hist(residuals, bins=bins, orientation='horizontal')
            ax_hist.axhline(0, color='k', linestyle='dashed')

        ax_main.autoscale_view(tight=True)

        self._draw_fig(fig)

    def OnDataProcessing(self, event):
        """ Open the dialog box for binning / smoothing data """
        if self.x_fit_array is not None:
            dlg = DataProcessingDlg(self, "Data Processing", wx.ID_ANY)
            if dlg.Show() == wx.ID_OK:
                dlg.Destroy()
        else:
            dlg = wx.MessageDialog(self, "Can't process data that hasn't been loaded...", "Nope.",
                                   wx.OK | wx.ICON_INFORMATION)
            dlg.ShowModal()
        event.Skip()

    def OnExit(self, event):
        """ What to do when the window is closed """
        # explicitly close all figures (bug with matplotlib and wx??)
        plt.close('all')
        event.Skip()
        self.Destroy()

    def OnFileOpen(self, event):
        """
        Open a csv data file and plot the data. Detuning is assumed to be in GHz.
        Vertical units are assumed to be the same as in the theory curves
        """
        self.dirname = ''
        dlg_choice = wx.SingleChoiceDialog(self, "Choose type of data to be imported", "Data import",
                                           choices=OutputTypes)

        # wait for OK to be clicked
        if dlg_choice.ShowModal() == wx.ID_OK:
            choice = dlg_choice.GetSelection()
            # print 'Choice:', choice
            self.expt_type = OutputTypes[choice]
            # use the choice index to select which axes the data appears on	- may be different
            # if axes order is rearranged later?
            self.choice_index = OutputTypes_index[choice]
            # print self.choice_index

            dlg_choice.Destroy()

            dlg_open = wx.FileDialog(self, "Choose 2-column csv file (Detuning, Transmission)",
                                     self.dirname, "", "*.csv", wx.FD_OPEN)

            # if OK button clicked, open and read file
            if dlg_open.ShowModal() == wx.ID_OK:
                # set experimental display on, and update menus
                self.display_expt_curves[self.choice_index] = True
                # self.showEplotsSubMenu.GetMenuItems()[self.choice_index].Check(True)

                self.filename = dlg_open.GetFilename()
                self.dirname = dlg_open.GetDirectory()
                # call read
                self.x_expt_arrays[self.choice_index], self.y_expt_arrays[self.choice_index] = np.loadtxt(
                    os.path.join(self.dirname, self.filename), delimiter=',', usecols=[0, 1]).T

                # overwrite fit_array data - i.e. last data to be loaded
                self.x_fit_array = self.x_expt_arrays[self.choice_index]
                self.y_fit_array = self.y_expt_arrays[self.choice_index]

                # implicit that the fit type is the same as last data imported
                self.fit_datatype = self.expt_type

                # create main plot
                self.OnCreateAxes(self.figs[0], clear_current=True)

            dlg_open.Destroy()
        event.Skip()

    def OnFitBounds(self, event):
        """ Set fit bounds (for detuning-axis range) """
        dlg = FitBoundsDialog(self, "Fit Bounds", wx.ID_ANY)

        if dlg.ShowModal() == wx.ID_OK:
            self.detuning_fit_bounds = dlg.OnUpdateRange()
        event.Skip()

    def OnFitButton(self, event):
        """ Call elecsus to fit data. Some sanity checking takes place first. """

        # check for things that will prevent fitting from working, e.g. no data loaded - halt fitting if found
        if len(self.y_fit_array) == 0:
            # warn about no data present
            dlg = wx.MessageDialog(self, "No experimental data has been loaded, cannot proceed with fitting...",
                                   "No no no", wx.OK | wx.ICON_EXCLAMATION)
            dlg.ShowModal()
            return

        self.fit_bools = [checkbox.IsChecked() for checkbox in self.fit_options.main_paramlist_bools]
        self.fit_bools += [checkbox.IsChecked() for checkbox in self.fit_options.magnet_paramlist_bools]
        self.fit_bools += [self.fit_options.fit_polarisation_checkbox.IsChecked()]
        if self.fit_bools.count(True) == 0:
            dlg = wx.MessageDialog(self, "No fit parameters are floating, cannot proceed with fitting...", "No no no",
                                   wx.OK | wx.ICON_EXCLAMATION)
            dlg.ShowModal()
            return

        # check for non-optimal conditions and warn user
        if self.fit_warnings:
            # if number of booleans > 3 and ML fitting selected, warn about fitting methods
            if self.fit_bools.count(True) > 3 and self.fit_algorithm == 'Marquardt-Levenberg':
                dlg = wx.MessageDialog(self,
                                       "The number of fitted parameters is large and the Marquardt-Levenberg "
                                       "algorithm is selected. There is a high probability that the fit "
                                       "will return a local minimum, rather than the global minimum. \n\n"
                                       "To find the global minimum more reliably, consider changing to either "
                                       "Random-Restart or Simulated Annealing algorithms.\n\n Continue with "
                                       "fitting anyway?",
                                       "Warning", wx.YES | wx.NO | wx.ICON_WARNING)

                if dlg.ShowModal() == wx.ID_NO:
                    return

            # if large number of data points
            if len(self.x_fit_array) > 5000:
                dlg = wx.MessageDialog(self,
                                       "The number of data points is quite high and fitting may be very slow, "
                                       "especially with Random-Restart or Simulated Annealing methods. \n\n"
                                       "Consider reducing the number of data points by binning data "
                                       "(Menu Fit --> Data Processing... --> Bin Data).\n\n "
                                       "Continue with fitting anyway?",
                                       "Warning", wx.YES | wx.NO | wx.ICON_WARNING)

                if dlg.ShowModal() == wx.ID_NO:
                    return

        # Bounds on fit parameters must be enabled for differential evolution
        if self.fit_algorithm == 'Differential Evolution':
            # status of checkboxes for floated params must equal checkboxes ticked for bounds
            ticked_bools = [i for i, box in
                            enumerate(self.fit_options.main_paramlist_bools +
                                      self.fit_options.magnet_paramlist_bools) if
                            box.IsChecked()]
            print(ticked_bools)
            ticked_bounds = [i for i, box in enumerate(
                self.fit_options.main_paramlist_usebounds + self.fit_options.magnet_paramlist_usebounds) if
                             box.IsChecked()]
            print(ticked_bounds)

            if ticked_bools != ticked_bounds:
                dlg = wx.MessageDialog(self,
                                       "Bounds on fit parameters must be specified to use "
                                       "Differential Evolution algorithm",
                                       "Bounds not specified", style=wx.OK | wx.CENTRE | wx.ICON_ERROR)
                dlg.ShowModal()
                return

        # If all conditions satisfied, start the fitting process
        self.call_elecsus('Fit')
        event.Skip()

    def OnFitCompleted(self, event):
        """ Task to complete after fitting completed """
        panel = self.fit_options

        self.FitInformation.write('\n============== Fit Completed ============ \n\n')

        # self.opt_params = self.opt_params_new

        # Add information to the Fitting text box
        self.FitInformation.write(
            '\nRMS error between theory and experiment (%): ' + str(round(self.rms * 100, 2)) + '\n')
        self.FitInformation.write(
            '\n-------------------------------------------------------------------------------------\n')

        self.FitInformation.write('Fit report from LM Fit:\n')
        self.FitInformation.write(self.fit_result.fit_report())

        # Add RMS info to the Fitting text box

        # use fitted parameters to calculate new theory arrays
        print('\n\n')
        print(time.ctime())
        print('Calling ElecSus for single calculation with optimised parameters')
        print(self.opt_params)

        self.y_optimised = self.fit_result.best_fit
        # e_in_opt = [self.opt_params['Ex'],[self.opt_params['Ey'],self.opt_params['Phase']]]
        e_in_opt = np.array(
            [self.opt_params['E_x'], self.opt_params['E_y'] * np.exp(1.j * self.opt_params['E_phase']), 0])
        spectrum_data = calculate(self.x_array * 1e3, e_in_opt, self.opt_params,
                                  outputs=['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I_P45', 'I_M45', 'Ir', 'Il',
                                           'alphaPlus', 'alphaMinus', 'alphaZ'])

        self.y_arrays[0:4] = spectrum_data[0:4]  # S0,S1,S2,S3
        self.y_arrays[4] = [spectrum_data[4], spectrum_data[5]]  # Ix,Iy
        self.y_arrays[5] = [spectrum_data[6], spectrum_data[7]]  # I+45,I-45
        self.y_arrays[6] = [spectrum_data[8], spectrum_data[9]]  # I_RCP, I_LCP
        self.y_arrays[7] = [spectrum_data[10], spectrum_data[11], spectrum_data[12]]  # alpha +/-/z

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
        self.OnCreateAxes(self.figs[0])
        print('Updating main plot...')

        # ask whether to update the initial fitting parameters

        dlg = wx.MessageDialog(self,
                               "Fit completed successfully. Details can be found in the Fitting Information Tab."
                               "\n\nReplace initial fit parameters with optimised parameters?",
                               "Fit Complete",
                               style=wx.YES | wx.NO)
        if dlg.ShowModal() == wx.ID_YES:

            # set optimised parameters into GUI control boxes

            self.main_param_labels = ['T', 'lcell', 'shift', 'GammaBuf', 'DoppTemp',
                                      'rb85frac', 'K40frac', 'K41frac']
            # strings of dictionary entries as used by elecsus

            self.opt_params['lcell'] *= 1e3
            # main parameters
            for i in range(len(self.main_param_labels)):
                panel.main_paramlist_inputs[i].SetValue(self.opt_params[self.main_param_labels[i]])

            # magnet parameters
            self.opt_params['Btheta'] /= (np.pi / 180)
            self.opt_params['Bphi'] /= (np.pi / 180)
            self.opt_params['E_phase'] /= (np.pi / 180)
            for i in range(len(self.magnet_param_labels)):
                panel.magnet_paramlist_inputs[i].SetValue(self.opt_params[self.magnet_param_labels[i]])

            # polarisation parameters
            for i in range(len(self.polarisation_param_labels)):
                panel.pol_control_inputs[i + 1].SetValue(self.opt_params[self.polarisation_param_labels[i]])

        dlg.Destroy()

        # create/update residuals plot
        self.on_create_residual_plot(self.figs[1])
        event.Skip()

    def OnFitPlotSelection(self, event):
        """ Select plots to show, via the popup button on the main panel """
        btn = event.GetEventObject()
        pos = btn.ClientToScreen((0, 0))
        dlg = PlotSelectionDialog(self.panel, self, 'Fit Plots Shown:', 'Fit', pos)

        if dlg.Show() == wx.ID_OK:
            dlg.Destroy()

    def OnFitTypeChangePanel(self, event):
        """ Action to perform when fit type is changed """
        # print self.fittypeSubMenu.GetMenuItems()
        print('')
        for panelitem in self.fit_options.fit_types:
            if panelitem.GetValue():
                self.fit_algorithm = panelitem.GetLabel()
                # self.fittypeSubMenu.Check(idfit,True)
                # self.OnFitTypeChangeMenu(1)
                print('Fit Algorithm changed to:', self.fit_algorithm)
        event.Skip()

    def OnFitWarnings(self, event):
        """
        Warn about bad fitting practices - e.g. when there are many data points to fit, or many
        fit parameters where the fit algorithm could be improved.
        """
        self.fit_warnings = bool(event.IsChecked())

    def OnGridToggleMain(self, event):
        """ Toggle axes grid on main plot """
        for ax in self.figs[0].axes:
            ax.grid(bool(event.IsChecked()))
        self._draw_fig(self.figs[0])

    def OnGridToggleRes(self, event):
        """ Toggle axes grid on residuals plot """
        for ax in self.figs[1].axes:
            ax.grid(bool(event.IsChecked()))
        self._draw_fig(self.figs[1])

    def OnLivePlotting(self, event):
        """
        Turn on/off live plotting when parameters are changed.
        Main plot updates in pseudo-realtime, rather than having to click 'compute' every time

        This method binds or unbinds a call to the compute button each time a control is changed.
        """
        live_plot_on = bool(event.IsChecked())

        if live_plot_on:
            for ctrl in self.thy_options.all_floatspin_inputs:
                self.Bind(EVT_FLOATSPIN, self.OnComputeButton, ctrl)

            for ctrl in self.thy_options.fixed_paramlist_inputs[0:2]:
                self.Bind(wx.EVT_COMBOBOX, self.OnComputeButton, ctrl)

            self.live_events_bound = True
        else:
            # unbind the events, if currently bound
            if self.live_events_bound:
                # unbind
                for ctrl in self.thy_options.all_floatspin_inputs:
                    self.Unbind(EVT_FLOATSPIN, ctrl)

                for ctrl in self.thy_options.fixed_paramlist_inputs[0:2]:
                    self.Unbind(EVT_FLOATSPIN, ctrl)

            self.live_events_bound = False

    def OnLoadConfig(self, event):
        """
        Load previous configuration settings from pickle file, for picking up where you left off...
        """
        dlg = wx.MessageDialog(self,
                               "Load previous configuration of the program - theory/fit parameters, "
                               "plot settings etc...\n\nNot implemented yet...",
                               "No no no", wx.OK)
        dlg.ShowModal()
        event.Skip()

    def OnPlotHold(self, event):
        """
        Toggle plot hold (keep data on updating figure) on/off
        Allows multiple data sets to be shown on same plot
        """
        # self.PlotHold = bool(event.IsChecked())
        # self.figs[0].hold(self.PlotHold)
        dlg = wx.MessageDialog(self, "Not implemented yet...", "No no no", wx.OK)
        dlg.ShowModal()
        event.Skip()

    def OnPlotLegend(self, event):
        """ Toggle plot legend on/off """
        self.legend_on = bool(event.IsChecked())

    def OnResHist(self, event):
        """ Turn histogram of residuals on/off """
        self.residual_histogram = bool(event.IsChecked())
        self.on_create_residual_plot(self.figs[1])

    def OnSaveConfig(self, event):
        """
        Save present configuration settings - plot types, data ranges, parameters - for faster repeating of common tasks
        """

        dlg = wx.MessageDialog(
            self,
            "Save current configuration of the program - theory/fit parameters, "
            "plot settings etc...\n\nNot implemented yet...",
            "No no no", wx.OK
        )
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

        # File Dialog window for selecting filename/location
        filename = './elecsus_gui_config.dat'

        # Save data in python-readable (binary) format using pickle module
        with open(filename, 'wb') as file_obj:
            pickle.dump(data_dump, file_obj)
        event.Skip()

    def OnSaveCSVData(self, event):
        """
        Method to save the main plot data traces as a n-column csv file.

        *All* data is saved, whether or not it is displayed

        This method selects the file name, then passes that to SaveTheoryCurves()
        """

        save_file_dialog = wx.FileDialog(self, "Save Output File", "./", "Outputs",
                                         "CSV files (*.csv)|*.csv", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if save_file_dialog.ShowModal() == wx.ID_OK:
            output_filename = save_file_dialog.GetPath()
            save_file_dialog.Destroy()

            print(output_filename)
            self.SaveTheoryCurves(output_filename)
        event.Skip()

    def OnSaveFig(self, event):
        """ Basically the same as saving the figure by the toolbar button. """
        # widcards for file type selection
        wilds = "PDF (*.pdf)|*.pdf|" \
                "PNG (*.png)|*.png|" \
                "EPS (*.eps)|*.eps|" \
                "All files (*.*)|*.*"
        exts = ['.pdf', '.png', '.eps', '.pdf']  # default to pdf
        save_file_dialog = wx.FileDialog(self, message="Save Figure", defaultDir="./", defaultFile="figure",
                                         wildcard=wilds, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        save_file_dialog.SetFilterIndex(0)

        if save_file_dialog.ShowModal() == wx.ID_OK:
            output_filename = save_file_dialog.GetPath()
            if output_filename[-4:] == exts[save_file_dialog.GetFilterIndex()]:
                output_filename = output_filename[:-4]

            # save all current figures
            for fig, f_id in zip(self.figs, self.fig_IDs):
                fig.savefig(output_filename + '_' + str(f_id) + exts[save_file_dialog.GetFilterIndex()])

        save_file_dialog.Destroy()
        event.Skip()

    def OnShowEplots(self, event):
        """ Action when plot type is changed from the menu """
        for ii, item in enumerate(self.show_e_plots_submenu.GetMenuItems()):
            if item.IsChecked():
                self.display_expt_curves[ii] = True
            else:
                self.display_expt_curves[ii] = False

        # redraw plot
        self.OnCreateAxes(self.figs[0])
        event.Skip()

    def OnShowTplots(self, event):
        """ Action when plot type is changed from the menu """
        for ii, item in enumerate(self.showTplotsSubMenu.GetMenuItems()):
            if item.IsChecked():
                self.display_theory_curves[ii] = True
            else:
                self.display_theory_curves[ii] = False

        # redraw plot
        self.OnCreateAxes(self.figs[0])
        event.Skip()

    def OnTheoryPlotSelection(self, event):
        """ Select plots to show, via the popup button on the main panel """

        btn = event.GetEventObject()
        pos = btn.ClientToScreen((0, 0))
        dlg = PlotSelectionDialog(self.panel, self, 'Theory Plots Shown:', 'Theory', pos)

        if dlg.Show() == wx.ID_OK:
            dlg.Destroy()

    def OnUseExpDetuning(self, event):
        """ Action when Menu Item 'Use Experimental Detuning' is clicked """
        if self.x_fit_array is not None:
            # only do anything if data has been loaded already
            self.use_exp_detuning = bool(event.IsChecked())
            self.OnComputeButton(1)
        else:
            problem_dlg = wx.MessageDialog(self, "No experimental data has been loaded yet", "No data to use",
                                           wx.OK | wx.ICON_ERROR)
            problem_dlg.ShowModal()
            self.em_use_exp_detuning.Check(False)

    def SaveTheoryCurves(self, filename):
        """
        Method to actually do the saving of xy data to csv file.
        Separate to the OnSaveCSVData method in case a filename is automatically chosen
        - i.e. when fitting data, there is an option to autosave the results which calls this method
        """
        # print len(self.y_arrays)
        # print self.y_arrays

        try:
            xy_data = list(zip(self.x_array, self.y_arrays[0].real, self.y_arrays[1].real, self.y_arrays[2].real,
                               self.y_arrays[3].real,
                               self.y_arrays[4][0].real, self.y_arrays[4][1].real,
                               self.y_arrays[5][0].real, self.y_arrays[5][1].real,
                               self.y_arrays[6][0].real, self.y_arrays[6][1].real,
                               self.y_arrays[7][0].real, self.y_arrays[7][1].real, self.y_arrays[7][2].real))
            success = write_csv(xy_data, filename, titles=['Detuning'] + OutputTypes)
            if not success:
                raise AttributeError

        except AttributeError:
            problem_dlg = wx.MessageDialog(self, "There was an error saving the data...", "Error saving",
                                           wx.OK | wx.ICON_ERROR)
            problem_dlg.ShowModal()


def write_csv(xy, filename, titles=None):
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
        with open(filename, 'wt') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',')
            if titles is not None:
                csv_writer.writerow(titles)
            for xy_line in xy:
                csv_writer.writerow(xy_line)
        return True
    except csv.Error:
        return False


# Run the thing...
def start():
    """ Start the GUI """
    print('Starting ElecSus GUI...')
    # global app
    app = wx.App(redirect=False)
    frame = ElecSusGuiFrame(parent=None, title="ElecSus GUI - General Magneto Optics")

    frame.SetSize((1200, 900))
    frame.Centre()
    frame.Maximize()
    frame.Show()
    app.MainLoop()
    print('...Closed ElecSus GUI')


if __name__ == '__main__':
    """ Called when program started with 'python elecsus_gui.py' from the command line """
    show_versions()
    start()
