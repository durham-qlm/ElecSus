- LUT: how to allow both change of power and diameter? or how to calculate corresponding intensity?
- beam_diameter is fixed, while power can be changed later... (transit time broadening)
- How to calculate isotope shifts? Currently taken from elecsus/libs/AtomConstants.py (different for D1/D2)
- If we remove 'Bfield':0.,'Btheta':0., 'Bphi':0.,'GammaBuf':0.,'shift':0 from p_dict (already non-bwf),
  we get "ValueError: array must not contain infs or Nans"

FITTING:
- Check so that lookup table does not has to be loaded multiple times
- Why does fit_data use strange E-field definition (
    params['E_x'].value = E_in[0]
	params['E_y'].value = E_in[1][0]
	params['E_phase'].value = E_in[1][1])