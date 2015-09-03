# Button (vertical) size
BtnSize = 30

## Master list of all output types that will be referenced for dynamic plotting
OutputTypes = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix', 'Iy', 'Alpha Plus', 'Alpha Minus', 'N Plus', 'N Minus','GI Plus', 'GI Minus','Phi']
# some of the outputtypes are grouped when plotting
OutputTypes_index = [0,1,2,3,4,4,5,5,6,6,7,7,8]
OutputPlotTypes = ['Transmission (S0)', 'S1', 'S2', 'S3', 'Ix, Iy', 'Alpha Plus/Minus', 'N Plus/Minus', 'GI Plus/Minus', 'Phi']

fixed_parameterlist = ['Element','D-line','Constrain Dopp./Density Temps']
element_list = ['Na', 'K', 'Rb', 'Cs']
D_line_list = ['D1', 'D2']
fittable_parameterlist = ['Bfield','Temperature','Cell length','Shift','Additional-Broadening','Theta-0','Sigma-Minus Polarisation','Doppler-Temperature','Rb-85','K40','K41']
units_parameterlist = [' [G]',' [C]',' [mm]',' [MHz]',' [MHz]',' [deg]',' [%]',' [C]',' [%]',' [%]',' [%]']

# default values and increments for the floatspin controls on the options panels
# order is the same as for the fittable_parameterlist list
defaultvals_parameterlist = [0,20,75,0,0,0,50,20,72.17,0.01,6.73]
defaultvals_increments = [1, 0.5, 0.5, 1, 1, 1, 1, 0.5, 1, 1, 1]
detuning_defaults = [-10, 10, 5000]
detuning_increments = [1,1,1000]


# Text strings for tooltips
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
fit_parameter_tooltips = [ \
'Applied magnetic field in Gauss. The field is currently only modeled along the axis of propagation.\n\
Positive numbers mean B and k are parallel, negative means anti-parallel',
'Vapour Temperature in degrees Celsius.\nThis is the temperature that sets the number density of the vapour',
'Cell length in mm.\nThis sets the amount of absorption/dispersion, via essentially the Beer-Lambert law',
'Shift the global line-centre by an amount in MHz. Negative numbers mean red-shift.',
'Additional homogeneous broadening, in MHz.\n\
This can be any source of homogeneous broadening, most commonly buffer gases.\n\
Note that the resonant alkali-alkali collisions (self-broadening or \n\
pressure-broadening) are already included in the model.',
'Initial axis of polarisation.\n\
If the light is completely linearly polarised, this is the axis of\n\
the Electric-field oscillation. 0 degrees = horizontal, 90 degrees = vertical.\n\
For circularly polarised light, this acts as a global \n\
phase factor and hence has no effect.',
'Percentage of the light that drives sigma-minus transitions. If the light is linearly polarised,\
in an axial magnetic field (which sets the quantisation axis) half the light drives sigma-plus, half sigma-minus',
'Doppler Temperature sets the Doppler width',
'Percentage of Rubidium-85 isotope. Default value is natural abundance. (72% Rb85, 28% Rb87)\n\
Obviously this has no effect when the element is set to something other than Rb',
'Percentage of Potassium-40 isotope. Default value is natural abundance. (99% Rb85, 28% Rb87)\n\
Obviously this has no effect when the element is set to something other than Rb',
]
parameter_adjust_tooltip = \
'\n\n\nShortcuts for faster adjustment of parameters:\n\
Shift + Arrow/MouseWheel = Adjust by 2x increment\n\
Ctrl + Arrow/MouseWheel = Adjust by 10x increment\n\
Alt + Arrow/MouseWheel = Adjust by 100x increment\n\
Escape = Set back to default value'
fit_algorithm_tooltips = [\
'Marquardt-Levenberg (ML) is the simplest fitting algorithm, and \n\
involves a simple downhill method to find a minimum in parameter space.\n\
Its simplicity makes it relatively fast, and a fit on a data set\n\
with around 5000 points takes around 15 seconds, depending on computer\n\
speed and choice of starting parameters. It runs on a single processing core.\n\
However, for large numbers of free parameters, it is known to find local \n\
minima and hence often returns non-optimal parameters.\n\
When the number of free parameters is more than 3, we recommend using either \n\
Random-Restart or Simulated Annealing instead.',
'Random-Restart takes a spread of starting parameters, centred\n\
around the initial parameters. The ML method is called separately\n\
on each of these starting parameters in order to find a global minimum.\n\
The speed of this method depends on the number of processing cores\n\
and number of free parameters in the fit.\n\
Typically this method takes 1-5 minutes.', \
'Simulated Annealing is the most complex algorithm included with ElecSus,\n\
and is therefore the most computationally intensive. However, it is \n\
most likely to return the global minimum, assuming a sensible starting condition.\n\
This fit uses all processor cores, but the method is slow, and should\n\
therefore only be used for the most complex fitting problems, where RR fails\n\
to find a good solution, or to check the results of an RR fit. \n\
Typical time to compute is around 10-15 minutes.' ]