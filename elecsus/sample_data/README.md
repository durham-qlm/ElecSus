# ElecSus sample data

Small, real experimental spectra bundled so the tutorials and tests have
something genuine to fit. Each file is two columns with a header row.

## `rb_d2_doppler_40C.csv`

A Doppler-broadened **Rb D2** absorption spectrum recorded in **zero magnetic
field** through a **75 mm** natural-abundance vapour cell at a nominal cell
temperature of **40 C**.

| column | meaning |
| --- | --- |
| `Detuning (MHz)` | linear detuning about the Rb D2 line centre |
| `Transmission` | normalised transmission (0 = fully absorbed, 1 = clear) |

A 3000-point uniform subset of the original ~166k-point trace. ElecSus
reproduces it to about **0.9 % RMS** with `{'Elem':'Rb', 'Dline':'D2',
'Bfield':0, 'lcell':75e-3}` and the best-fit `T` recovering ~40 C. It is used
in `Notebooks/06_fitting.ipynb` and pinned by `elecsus/tests/test_sample_data.py`.

### Source and attribution

Data from the **scope2elecsus** tutorial dataset,
<https://github.com/durham-qlm/scope2elecsus>, accompanying:

> D. Pizzey, J. D. Briscoe, F. D. Logue, F. S. Ponciano-Ojeda,
> S. A. Wrathmall and I. G. Hughes, "Laser spectroscopy of hot atomic vapours:
> from 'scope to theoretical fit," *New J. Phys.* **24**, 125001 (2022).
> doi:10.1088/1367-2630/ac9cfe

Open access, CC BY 4.0. Reproduced here as a downsampled subset with
attribution; please cite the paper above if you use it.
