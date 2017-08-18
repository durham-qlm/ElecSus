# Copyright 2017 J. Keaveney

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" Variables used by the GUI. For formatting, default values etc."""

# Button (vertical) size
BtnSize = 30

## Master list of all output types that will be referenced for dynamic plotting
OutputTypes = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'I+45', 'I-45', 'Ircp', 'Ilcp', 'Alpha Plus', 'Alpha Minus', 'Alpha Z']
# some of the outputtypes are grouped when plotting
OutputTypes_index = [0,1,2,3,4,4,5,5,6,6,7,7,7]
OutputPlotTypes = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix, Iy', 'I+45 / I-45', 'Ircp / Ilcp', 'Alpha Plus/Minus/Z']

fixed_parameterlist = ['Element','D-line','Constrain Dopp./Density Temps']
element_list = ['Na', 'K', 'Rb', 'Cs']
D_line_list = ['D1', 'D2']
polarisation_options = ['Linear','LCP','RCP','Elliptical']
polarisation_controls = ['Theta-0','Ex','Ey','Phase [deg]']
magnet_parameterlist = ['Bfield', 'Btheta', 'Bphi']
magnetunits_parameterlist = [' [G]',' [deg]',' [deg]']
main_parameterlist = ['Temperature','Cell length','Shift','Additional-Broadening','Doppler-Temperature','Rb-85','K40','K41']
mainunits_parameterlist = [' [C]',' [mm]',' [MHz]',' [MHz]',' [C]',' [%]',' [%]',' [%]']

# default values and increments for the floatspin controls on the options panels
# order is the same as for the fittable_parameterlist list
defaultvals_magnet_parameterlist = [0, 0, 0]
defaultvals_magnet_increments = [1, 0.5, 0.5]
defaultvals_pol_parameterlist = [0, 1, 0, 0]
defaultvals_pol_increments = [1, 0.01, 0.01, 1]
defaultvals_main_parameterlist = [20, 75, 0, 0, 20, 72.17, 0.01, 6.73]
defaultvals_main_increments = [0.5, 0.5, 1, 1, 0.5, 1, 1, 1]
main_paramlist_mindefaults = [0, 0.001, -4000, 0, -273.15, 0, 0, 0]
main_paramlist_maxdefaults = [1000, 100, 4000, 1000, 1000, 100, 100, 100]
magnet_paramlist_mindefaults = [0, 0, 0]
magnet_paramlist_maxdefaults = [20000, 360, 360]
detuning_defaults = [-20, 20, 2500]
detuning_increments = [1,1,1000]


# Text strings for tooltips
magnet_parameter_tooltips = [ \
'Applied magnetic field in Gauss. Taken as absolute (only pos values allowed)', \
'Angle magnetic field axis makes with the z-axis (k-vector)', \
'Angle magnetic field axis makes with the x-axis' \
]

polarisation_tooltips = ['Linearly polarised light', \
'Left-hand circularly polarised light', \
'Right-hand circularly polarised light', \
'Elliptically polarised light - set ellipticity through Ex, Ey, and relative phase controls']

polarisation_parameter_tooltips = ['The angle the axis of\n\
the Electric-field oscillation makes (for linearly polarised light).\n\
0 degrees = polarised along x, 90 degrees = polarised along y.\n\
For circularly polarised light, this acts as a global \n\
phase factor and hence has no effect.', \
'Electric field amplitude in the x-axis', \
'Electric field amplitude in the y-axis', \
'Relative phase difference [in degrees] between the x and y components of the electric field']

fit_polarisation_tooltip = 'Fit the angle of polarisation. \nEx, Ey and Phase will all be allowed to vary. \nEx and Ey are constrained such that Ex^2 + Ey^2 = 1, \nand Phase is bounded between 0 and 360 degrees.'

constrain_polarisation_tooltip = 'Constrains the polarisation to be linearly polarised, \ni.e. only fit the angle Theta-0, with Phase = 0'

fixed_parameter_tooltips = [\
'Select the alkali element to use', \
'Select the D-line to use.\n\
D1 = nS_1/2 to nP_1/2\n\
D2 = nS_1/2 to nP_3/2', \
'Constrain the Doppler temperature (which sets Doppler width)\n\
and the temperature used to calculate number density. In most cases\n\
this should be checked, but in cells where the reservoir temperature\n\
is significantly different to the window temperature, allowing the two\n\
temperatures to vary independently may be a better idea.' ]
main_parameter_tooltips = [ \
'Vapour Temperature in degrees Celsius.\nThis is the temperature that sets the number density of the vapour',
'Cell length in mm.\nThis sets the amount of absorption/dispersion, via essentially the Beer-Lambert law',
'Shift the global line-centre by an amount in MHz. Negative numbers mean red-shift.',
'Additional homogeneous broadening, in MHz.\n\
This can be any source of homogeneous broadening, most commonly buffer gases.\n\
Note that the resonant alkali-alkali collisions (self-broadening or \n\
pressure-broadening) are already included in the model.',
'Doppler Temperature sets the Doppler width',
'Percentage of Rubidium-85 isotope. Default value is natural abundance. (72% Rb85, 28% Rb87)\n\
This has no effect when the element is set to something other than Rb',
'Percentage of Potassium-40 isotope. Default value is natural abundance. (93.28% K39, 0.01%K40, 6.73% K41)\n\
This has no effect when the element is set to something other than K',
'Percentage of Potassium-41 isotope. Default value is natural abundance. (93.28% K39, 0.01%K40, 6.73% K41)\n\
This has no effect when the element is set to something other than K'
]
parameter_adjust_tooltip = \
'\n\n\nShortcuts for faster adjustment of parameters:\n\
Shift + Arrow/MouseWheel = Adjust by 2x increment\n\
Ctrl + Arrow/MouseWheel = Adjust by 10x increment\n\
Alt + Arrow/MouseWheel = Adjust by 100x increment\n\
Escape = Set back to default value'
fit_algorithm_tooltips = [\
'See the documentation (GUI user manual and journal articles) for detailed\n\
information on fitting methods. For simple problems, we recommend \n\
Marquardt-Levenburg, but this is likely to find a local optimum in \n\
parameter space for complex problems. For complex problems we \n\
recommend using Differential Evolution, which is faster and more reliable \n\
than either the RR or SA methods.']*4