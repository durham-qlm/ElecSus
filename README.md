# ElecSus v4

> **ElecSus 4.0.0b1 is a beta.** Until 4.0.0, its interface and
> defaults may still change. Some of its results differ from ElecSus 3's on
> purpose: see [Backwards compatibility](#backwards-compatibility) and the
> [CHANGELOG](CHANGELOG.md). Please report problems on the
> [issue tracker](https://github.com/durham-qlm/ElecSus/issues).
>
> **ElecSus 3** (version 3.0.8) is on the [`elecsus3`](https://github.com/durham-qlm/ElecSus/tree/elecsus3)
> branch, tagged `v3.0.8`; its README there says how to install it. ElecSus 3
> and 4 both install a package called `elecsus`, so keep them in separate
> environments.

**ElecSus** calculates the electric susceptibility, refractive index, and transmission
spectra of thermal alkali-metal vapour cells on the D-lines. It includes least-squares
fitting routines for extracting experimental parameters from spectroscopic data.

This is **version 4**, a significant extension of
[ElecSus v3](https://doi.org/10.1016/j.cpc.2017.12.001) (Keaveney et al. 2018).
A companion paper describing the new physics is in preparation for
*Computer Physics Communications*.

> The v3 API (`calculate`, `fit_data`) keeps its arguments. Existing scripts run
> unchanged unless they relied on the absolute (un-normalised) values of S1, S2, S3,
> on the sense of `Bphi` or of a rotation (v3 computed the mirror image of its
> documented geometry), on input that v3 ignored and v4 now rejects, on
> simulated annealing (`fit_algorithm='SA'`), which has been removed, on
> v3's propagation in oblique fields (wrong by up to ~0.5 in S0; v4 solves
> for the non-orthogonal modes by default), or on spectra identical to v3's
> to better than ~3×10⁻³ (v3 truncated the line frequencies to whole MHz) — see [Backwards compatibility](#backwards-compatibility).

---

## What's new in v4

| Feature | Description |
|---|---|
| **Non-orthogonal propagation** | Correct biorthogonal mode decomposition for oblique magnetic fields (Briscoe et al., *Opt. Express* **33**, 32099, 2025). The default (`use_non_orthogonal=True`); `False` gives ES3's method, which is wrong in oblique fields by up to ~0.5 in S0. |
| **Fast oblique solver** | Webber's fully vectorised analytic eigenvalue/eigenvector solution replaces the ES3 SymPy + Python for-loop: ~130–200× faster, no unphysical artefacts at high fields. |
| **Dynamic buffer gas model** | `BufferGas`, `BufferPressure`, `BufferTemp_fill` parameters compute temperature- and pressure-dependent broadening and shift automatically. |
| **Cell-window Fresnel correction** | Full 4-interface model (air → glass → vapour → glass → air) applied per eigenmode. Opt in with `n_window` (1.51 for borosilicate); the default `None` is a windowless cell. |
| **Cascaded cells** | `calculate_cascade()` — passes the output electric field of cell 1 as the input to cell 2, correctly propagating polarisation state and intensity. |
| **Fractional Stokes parameters** | S1, S2, S3 now normalised by S0 (lying in [−1, +1]), with ES3's (Born & Wolf) sign for S3. |
| **The documented geometry** | `Btheta`/`Bphi` mean what the documents say: **B** along (sin θ cos φ, sin θ sin φ, cos θ). ES3 computed the mirror image of that geometry; see [Backwards compatibility](#backwards-compatibility). |
| **Validated parameter dictionary** | Unknown keys and unphysical values are rejected with a named error instead of being ignored. A mistyped `n_windows` or `use_nonorthogonal` used to switch the new physics off in silence. |
| **Real-valued intensity outputs** | `S0`–`S3`, `Ix`, `Iy`, `I_P45`, `I_M45`, `Il`, `Ir`, `psi` and `DoLP` are `float64`; they were `complex128` with a ~1e-17 imaginary residue. |

---

## Installation

ElecSus 4 needs Python 3.10 or newer. It is on PyPI as **`elecsus4`**, and is
imported as `elecsus`:

```bash
pip install "elecsus4[gui]"     # with the browser GUI; plain elecsus4 without it
```

Install it in a new environment (a conda environment or a `venv`). ElecSus 3
installs a package of the same name, `elecsus`, and the two must not share an
environment: if ElecSus 3 is installed, `pip uninstall ElecSus` first.

**From the source**, for the tests, the notebooks, or to change the code:

```bash
# 1. Get the source and enter it
git clone https://github.com/durham-qlm/ElecSus.git
cd ElecSus

# 2. Create and activate the conda environment
#    (installs Python 3.12, NumPy, SciPy, lmfit, Dash, and all other dependencies)
conda env create -f environment.yml
conda activate elecsus4

# 3. Install ElecSus in editable mode
pip install -e .
```

**Existing environment?**  Activate it and install ElecSus with the GUI
dependencies:

```bash
conda activate elecsus4         # or whatever your env is called
pip install -e ".[gui]"         # ElecSus plus the Dash GUI stack
```

---

## Quick start

### Basic transmission spectrum (ES3-compatible)

```python
import numpy as np
from elecsus import elecsus_methods as em

detuning = np.linspace(-10000, 10000, 1000)   # MHz

p = {
    'Elem':   'Rb',
    'Dline':  'D2',
    'T':      20.0,    # °C
    'lcell':  75e-3,   # m
    'Bfield': 100,     # Gauss
}

[S0] = em.calculate(detuning, [1, 0, 0], p, outputs=['S0'])
```

### Non-orthogonal propagation (oblique B field)

```python
p = {
    'Elem':   'Rb',   'Dline': 'D2',
    'T':      60.0,   'lcell': 75e-3,
    'Bfield': 300,    'Btheta': np.pi/4,   # radians from +z (the forward light's direction)
    'use_non_orthogonal': True,        # Briscoe et al. 2025 formalism (the default)
    'n_window': 1.51,                  # opt in to the window model (default None)
}

[S0, S1, S3] = em.calculate(detuning, [1, 0, 0], p, outputs=['S0', 'S1', 'S3'])
```

### Buffer gas model

```python
p = {
    'Elem':  'K',   'Dline': 'D1',
    'T':     120.0, 'lcell': 75e-3,
    'Bfield': 0,
    'BufferGas':      'Ne',
    'BufferPressure':  96.0,   # Torr fill pressure
    'BufferTemp_fill': 293.15, # K at filling
}

[S0] = em.calculate(detuning, [1, 0, 0], p, outputs=['S0'])
```

### Cascaded vapour cells

```python
p1 = {'Elem': 'Rb', 'Dline': 'D2', 'T': 20, 'lcell': 75e-3, 'Bfield': 100}
p2 = {'Elem': 'Rb', 'Dline': 'D2', 'T': 80, 'lcell': 50e-3, 'Bfield':   0}

(S0_cell1,), (S0_cell2,) = em.calculate_cascade(
    detuning, [1, 0, 0], p1, p2, outputs=['S0'])

S0_combined = S0_cell1 * S0_cell2   # total cascaded transmission
```

The output electric field from cell 1 feeds as the input to cell 2, preserving
the polarisation state and intensity attenuation. See
`elecsus/tests/test_cascade.py` for the combination rules each output type
obeys (identity, product and polarisation-coupling cases).

For cascaded magneto-optical filters (MOFs) an optional linear polariser can
be inserted between the cells — e.g. crossed with respect to an x-polarised
input, the standard MOF filter-section topology:

```python
import numpy as np

(S0_1,), (S0_total,) = em.calculate_cascade(
    detuning, [1, 0, 0], p1, p2, outputs=['S0'],
    polariser_angle=np.pi/2)          # radians; π/2 = crossed
```

With `polariser_angle` set, the cell-2 outputs are **instrument totals**
normalised to the original input (they already include the cell-1 attenuation
and the polariser loss), so no combination rule is needed. Per-stage
normalisation would be 0/0 in the crossed-polariser wings where almost no
light reaches cell 2.

For wing-selector instruments that add a quarter-wave plate after the
polariser — preparing a *fixed* circular state for cell 2 regardless of
wavelength, e.g. VAMOS (Oliviero, PhD thesis, Univ. of Naples, 1998) — use
`reset_polarisation` (`'left_cp'` / `'right_cp'`) together with
`polariser_angle` (without it, `reset_polarisation` raises: no passive optic
turns every polarisation into one circular state without the polariser):

```python
(S0_1,), (S0_total,) = em.calculate_cascade(
    detuning, [1, 0, 0], p1, p2, outputs=['S0'],
    polariser_angle=np.pi/2, reset_polarisation='left_cp')
```

Cell-2 outputs are instrument totals here too. In this mode the total equals
the simple product of two independent single-cell calculations.

### Fitting

```python
import numpy as np
from elecsus import elecsus_methods as em

raw = np.loadtxt('my_data.csv', delimiter=',')   # (N, 2) array: detuning, signal
data = [raw[:, 0], raw[:, 1]]                    # fit_data expects [x_array, y_array] (2×N), not raw

p0   = {'Elem': 'Rb', 'Dline': 'D2', 'T': 20, 'lcell': 75e-3, 'Bfield': 0}
vary = {'T': True, 'Bfield': False, 'lcell': False}

best_params, rms, result = em.fit_data(data, p0, vary, data_type='S0', fit_algorithm='ML')
# best_params  — dict of fitted parameter values
# rms          — RMS residual between data and best fit
# result       — lmfit MinimizerResult (result.best_fit, result.fit_report(), etc.)
```

`fit_algorithm` accepts `'ML'` (Marquardt-Levenberg, default), `'RR'` (random-restart),
or `'differential_evolution'` (`'DE'`). Simulated annealing (`'SA'`) has been removed;
use differential evolution for a global fit.

To fit the input polarisation, vary `E_x` or `E_y` and `E_phase`, with `E_in`
as the start. Not both amplitudes: every output is normalised by the input
intensity, so only their ratio is observable, and `fit_data` raises a
`ValueError` if both vary.

If your cell was filled with a buffer gas, pass `BufferGas`/`BufferPressure`/
`BufferTemp_fill` in `p0` as usual. The model's broadening then replaces
`GammaBuf`, and its pressure shift is added to `shift`, both evaluated at the
temperature of the gas in the beam (`T`, or `DoppTemp` with `Constrain=False`)
at each fit iteration. In `vary`:

- `BufferPressure` may vary: the fit finds the amount of gas in the cell, from
  a starting pressure above 0 Torr;
- `shift` may vary, as a frequency offset on top of the pressure shift;
- `GammaBuf` may not, as the model replaces it, and nor may `BufferTemp_fill`:
  the model depends on `BufferPressure / BufferTemp_fill` only, so the two are
  degenerate. Either raises a `ValueError`.

---

### Energy levels and transitions

`elecsus.energy_levels` gives what a spectrum is made of, for one isotope
(`Na`, `K39`, `K40`, `K41`, `Rb85`, `Rb87`, `Cs`) and D line in a field: the
levels (with their |mJ, mI⟩ components and the |F, mF⟩ state each joins at
zero field), the σ+, σ− and π lines between them with their detunings and
strengths (the ones the spectra use), the levels against field, and a figure
of all of it above the transmission spectrum, with a stick spectrum underneath
(every line's relative strength on a log scale, so the weak lines the field
mixes in can be read off).

```python
from elecsus import energy_levels as EL

levels = EL.energy_levels('Rb87', 'D1', 6000.)          # B in gauss
lines = EL.transitions('Rb87', 'D1', 6000., kinds=('sigma+', 'sigma-'))
fig = EL.level_diagram('Rb87', 'D1', 6000., T=75., lcell=2e-3)   # matplotlib
fig.savefig('rb87_d1_levels.pdf')
ax = EL.stick_spectrum('Na', 'D1', 2440.)                 # the sticks on their own
```

Or from the command line, which also prints the table of lines:

```bash
python -m elecsus.energy_levels Rb87 D1 6000 --T 75 --out rb87_d1_levels.pdf
python -m elecsus.energy_levels Cs D2 1000 --geometry voigt --labels FmF --out cs.png
```

`--geometry voigt` (field across the beam) adds the π lines; `--labels FmF`
names the levels by their zero-field states; `--min-strength` (default 0.01)
leaves out arrows for lines weaker than that fraction of the strongest (the
stick spectrum still shows them); `--no-sticks` leaves out the stick spectrum.
The GUI's *Energy Levels* tab draws the same figure.

## Notebooks

Tutorial Jupyter notebooks are in the [`Notebooks/`](Notebooks/) directory:

| Notebook | Topic |
|---|---|
| [`01_basic_transmission.ipynb`](Notebooks/01_basic_transmission.ipynb) | D1/D2 transmission across Rb, K, Cs; temperature and field effects |
| [`02_magnetic_field_and_polarimetry.ipynb`](Notebooks/02_magnetic_field_and_polarimetry.ipynb) | Stokes parameters, Faraday rotation, circular dichroism |
| [`03_non_orthogonal_propagation.ipynb`](Notebooks/03_non_orthogonal_propagation.ipynb) | Oblique B-field propagation, Briscoe 2025 correction vs standard |
| [`04_buffer_gas.ipynb`](Notebooks/04_buffer_gas.ipynb) | Dynamic buffer gas model; pressure broadening and shift (K/Ne) |
| [`05_cascaded_cells.ipynb`](Notebooks/05_cascaded_cells.ipynb) | Two-cell cascade with polarisation propagation; FADOF and the crossed-polariser MOF (`polariser_angle`) |
| [`06_fitting.ipynb`](Notebooks/06_fitting.ipynb) | Fitting measured spectra with `fit_data`: ML and differential-evolution, parameter uncertainties, and buffer-gas fits |

The notebooks are committed **with their executed output**, so the expected
plots render on GitHub without running anything. Launch to run them yourself
with:

```bash
jupyter notebook Notebooks/
```

---

## Key parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `Elem` | str | `'Rb'` | Element: `'Rb'`, `'Cs'`, `'K'` or `'Na'` |
| `Dline` | str | `'D2'` | Transition: `'D1'` or `'D2'` |
| `T` | float | `20` | Temperature (°C) — sets number density, and the buffer-gas broadening and shift |
| `lcell` | float | `75e-3` | Cell length (m) |
| `Bfield` | float | `0` | Magnetic field magnitude (Gauss) |
| `Btheta` | float | `0` | Angle of B from the +z axis, the forward pass's direction (**radians**, on a reversed pass too; the GUI inputs are in degrees and convert internally) |
| `Bphi` | float | `0` | Azimuthal angle of B (**radians**; the GUI inputs are in degrees and convert internally) |
| `GammaBuf` | float | `0` | Extra phenomenological Lorentzian broadening (MHz). Adds to the natural linewidth. **Replaced** by the buffer gas model's broadening when `BufferGas` is set (with a warning if it is not 0) — do not use both. |
| `shift` | float | `0` | Global frequency shift of all resonance lines (MHz). When `BufferGas` is set, the Pitz pressure shift is **added** to it, so it acts as an offset of the frequency axis (it used to be replaced). |
| `BufferGas` | str | — | Buffer gas species: `'He'`, `'Ne'`, `'Ar'`, `'Kr'`, `'Xe'`, `'N2'`, `'H2'`, `'CH4'`, `'C2H6'`, `'C3H8'`, `'C4H10'`; which of them have data depends on the element and line (58 entries, each with its source in `elecsus/libs/BufferGasData.py`). Activates the dynamic broadening and shift model of Pitz *et al.* and others, evaluated at the temperature of the gas in the beam (`T`, or `DoppTemp` with `Constrain=False`). **When set, `GammaBuf` is ignored, and the pressure shift is added to `shift`.** |
| `BufferPressure` | float | — | Fill pressure (Torr) — the pressure at which the sealed cell was filled at room temperature. The actual pressure at operating temperature is computed internally. A fit can vary it. |
| `BufferTemp_fill` | float | `293.15` | Temperature (K) at which the cell was sealed. Combined with `BufferPressure` to fix the number of buffer gas atoms. Default is 293.15 K (20 °C). A fit cannot vary it: only `BufferPressure / BufferTemp_fill` matters. |
| `use_non_orthogonal` | bool | `True` | Briscoe et al.'s biorthogonal propagation, exact in any field direction. `False` is ES3's orthogonal method, exact only for `Btheta` = 0, π/2 or π: in oblique fields it is off by up to ~0.5 in S0 (see [Backwards compatibility](#backwards-compatibility)) |
| `n_window` | float or None | `None` | Window glass index, e.g. 1.51 for borosilicate. The default `None` models a windowless cell, so both propagation paths agree on absolute transmission; setting it costs ~15.5% off resonance |
| `rb85frac` | float | `72.17` | ⁸⁵Rb fraction (%) |
| `DoppTemp` | float | `T` | Temperature of the vapour in the beam (°C): the Doppler width and the buffer-gas model. Used only with `Constrain=False`, for a cell whose body is hotter than its stem; `T` is then the stem's temperature, which sets the number density |
| `Constrain` | bool | `True` | Force DoppTemp = T |

All parameters have defaults; only those you need to change must be specified.

**Keys are validated.** A key ElecSus does not recognise raises `ValueError`
naming the closest valid parameter, and values outside the range where the model
is defined (a negative `lcell`, an isotope fraction outside 0–100, a temperature
below absolute zero) raise before anything is computed. Up to v3 these were
ignored in silence, so `n_windows` or `use_nonorthogonal` returned a plausible
spectrum with the v4 physics switched off. Pass `strict=False` to `calculate()`,
`calculate_cascade()` or `fit_data()` if a workflow genuinely needs to carry
extra keys through the dictionary; the range checks still apply. `fit_data()`
also checks that every key of `p_dict_bools` and `p_dict_bounds` is a fit
parameter. Settings that are valid but have no effect warn: `n_window` without
`use_non_orthogonal=True`, `BufferPressure` or `BufferTemp_fill` without
`BufferGas`, and `E_x`/`E_y`/`E_phase`, which only `fit_data()` reads.

```python
>>> calculate(x, E_in, {'Elem': 'Rb', 'lcel': 75e-3})
ValueError: unknown parameter key: 'lcel' (did you mean 'lcell'?). ...
```

---

## Outputs

Pass a list of strings to the `outputs` keyword of `calculate()`:

| Name | Description |
|---|---|
| `S0` | Total transmitted intensity (normalised to input) |
| `S1` | Fractional Stokes: (Ix − Iy) / S0 |
| `S2` | Fractional Stokes: (I₊₄₅ − I₋₄₅) / S0 |
| `S3` | Fractional Stokes: (Ir − Il) / S0, positive for right-circular light (Born & Wolf, as in ES3). Like every output, about the lab +z axis: on a reversed pass (`wavevector_dirn=-1`) S3 > 0 is left-circular for the returning beam |
| `Ix`, `Iy` | Intensities through x- and y-polarisers |
| `I_P45`, `I_M45` | Intensities through ±45° polarisers |
| `Il`, `Ir` | Left- and right-circular intensities (`[1,+i]/√2` and `[1,-i]/√2`) |
| `alphaPlus`, `alphaMinus`, `alphaZ` | Absorption coefficients (σ⁺, σ⁻, π) |
| `ChiPlus`, `ChiMinus`, `ChiZ` | Complex susceptibilities of the σ⁺, σ⁻, π manifolds — see the naming note below |
| `nPlus`, `nMinus` | Complex refractive indices √(1+χ) of the σ⁺ and σ⁻ modes (Faraday geometry only) |
| `TransLeft`, `TransRight` | Single-mode σ⁺ / σ⁻ transmissions (Faraday geometry only) |
| `phi`, `phiPlus`, `phiMinus` | Net Faraday rotation (from x towards y, for **B** along +z or −z) and the two single-mode phases (Faraday geometry only; **dispersion-only** — see `psi`) |
| `psi` | Polarisation-ellipse azimuth ½·atan2(S2, S1) from the full field propagation — the rotation angle a polarimeter actually reports. Equals `phi` (mod π) in Faraday geometry (`Btheta` = 0 or π); remains correct in *any* geometry, where `phi` does not. Principal value in (−π/2, π/2]; pair with `DoLP` for validity |
| `DoLP` | Degree of linear polarisation √(S1²+S2²) ∈ [0, 1]. Validity diagnostic for any rotation angle: DoLP → 0 where dichroism has circularised the light and rotation is unmeasurable |
| `E_out` | Output electric field vector, shape (3, N) |

If `outputs` is omitted, the default set (S0, S1, S2, S3, Ix, Iy, I_P45, I_M45,
alphaPlus, alphaMinus, alphaZ) is returned as a tuple.

### The `Plus`/`Minus` naming convention

`Plus`/`Minus` is the σ label of the atomic transitions:

| Output key | Is physically the… | Governs, with **B** along +z |
|---|---|---|
| `ChiPlus`, `nPlus`, `phiPlus`, `alphaPlus`, `TransLeft` | **σ⁺** (Δm = +1) manifold | `[1,+i]/√2` (left-circular) light |
| `ChiMinus`, `nMinus`, `phiMinus`, `alphaMinus`, `TransRight` | **σ⁻** (Δm = −1) manifold | `[1,-i]/√2` (right-circular) light |

**Changed from ES3.** ES3, and the 4.0.0 alpha pre-releases, returned the σ⁻
quantities under the `Plus`/`Left` names (all but `alphaPlus`, which crossed
over internally): the labelling of their mirrored geometry (see
[Backwards compatibility](#backwards-compatibility)). ES3's docstrings
described them as σ⁺, which they now are. ES3 code that read `nPlus`,
`ChiPlus`, `phiPlus` or `TransLeft` by name gets ES3's values from the
`Minus`/`Right` keys. Likewise `spectra.calc_chi` returns
`(ChiPlus, ChiMinus, ChiZ)` with σ⁺ first, where ES3 returned σ⁻ first, so
ES3 code that unpacked its first element as the left-circular (σ⁺)
susceptibility gets ES3's from the second. The mapping is pinned by
`elecsus/tests/test_faraday_outputs.py` and
`elecsus/tests/test_geometry_conventions.py`, and the GUI labels each output
by its σ.

---

## Web GUI

The interactive GUI runs in any browser. It needs the `gui` extra
(`pip install "elecsus4[gui]"`); `environment.yml` includes it.

```bash
# 1. Activate the environment ElecSus is installed in
conda activate elecsus4         # or source .venv/bin/activate, ...

# 2. Launch the GUI
python -m elecsus.elecsus_dash_gui

# 3. Open your browser at
#    http://127.0.0.1:8050
```

To stop the server press **Ctrl-C** in the terminal.

**Port conflict?**  If port 8050 is already in use (e.g. from a previous run),
pass `--port` to pick a free one:

```bash
python -m elecsus.elecsus_dash_gui --port 8051
# then open http://127.0.0.1:8051
```

To find and kill a stale process on port 8050:

```bash
lsof -ti tcp:8050 | xargs kill   # macOS / Linux
```

The GUI supports:

- Single-cell forward calculation for all output types (S0, S1, S2, S3, Ix, Iy, …),
  with outputs of different units (e.g. α in m⁻¹, ψ in rad) in separate panels
- Non-orthogonal propagation (Briscoe et al. 2025) and the optional cell-window
  Fresnel model, as switches in the parameter panel (the solver on and the
  windows off by default, as in the API); `S0_P` is available with the
  non-orthogonal solver
- Optional independent Doppler temperature, for dual-temperature cells
- Fit-to-data on the single-cell tab (CSV upload, ML or differential-evolution),
  with **Export fit CSV** (detuning, data, best fit, residual). **Load sample
  spectrum** loads the bundled Rb D2 spectrum and sets the panel to match, so a
  fit works with no data of your own. Fittable: T, B, GammaBuf, shift, the
  buffer-gas pressure, cell length, field angles, Doppler T, ⁸⁵Rb fraction and
  the linear polarisation angle. Files may be comma-, tab- or space-separated,
  with header lines; detuning is read as GHz unless the first column's header
  says MHz. A blank *Fit min* or *Fit max* leaves that end of the range open.
  The results panel says when a fit did not move from its start or a parameter
  is not determined by the data
- The core's warnings (a buffer gas at 0 Torr, no data for the element, …)
  appear in a box on the page
- Cascaded two-cell filters (same species in both cells, element/D-line locked
  to cell 1) with an optional intermediate linear polariser — default 90° from
  the x axis, which is crossed for x-polarised input (the default): the
  magneto-optical-filter (MOF) topology for laser input light.
  The plot shows the spectrum after cell 1 and the total after both filters.
- Buffer gas parameters (`BufferGas`, `BufferPressure`, `BufferTemp_fill`)
- An *Energy Levels* tab: for one isotope in a field, the levels against field,
  the levels at the field, an arrow for each σ+, σ− (and, with the field
  across the beam, π) line at its detuning, the transmission spectrum above
  and a log-scale stick spectrum below, with a table of the lines and
  **Export lines CSV**
- Interactive plots with zoom/pan and CSV export of the current spectrum
- Dark/light theme is fixed to the Durham light theme; the theme and font are
  bundled, so the GUI works with no network connection

---

## Backwards compatibility

- `calculate()` and `fit_data()` take the same positional arguments as in ES3.
  What they used to ignore now raises:
  - an unknown `p_dict` key (`strict=False` downgrades it to a warning);
  - a misspelt keyword argument to `fit_data()`;
  - an unknown `fit_algorithm`, which used to run random restart.
- **`fit_algorithm='SA'` (simulated annealing) has been removed** and raises
  `ValueError`. It had not run on Linux or macOS since ES3 3.0.0, and
  differential evolution (`'DE'`) does the same job faster.
- **The geometry is the documented one.** ES3 computed the mirror image
  (y → −y) of its documented geometry: the senses of `Bphi`, of a field along
  **k** and of a Faraday rotation were reversed. ES3's result for
  (`Btheta`, `Bphi`) is v4's for (π − `Btheta`, −`Bphi`), to 3×10⁻⁸ in
  `S0`–`S3` (10⁻⁴ in the two oblique test cases, where the solver also
  changed); the golden tests check this. Transmission of linear light in a zero or longitudinal
  field is the same at the same angles. The 4.0.0 alpha pre-releases computed
  the mirror image too.
- **Line frequencies are exact.** ES3 truncated every transition frequency to
  whole MHz, so each line sat up to 1 MHz too close to the line centre, and
  the spectrum changed in steps as `Bfield` changed. That stalled ML fits of
  `Bfield` near their start, and made the error bar on B from a DE fit ~10×
  too large. v4 keeps the exact frequencies. Against ES3, S0
  moves by up to ~3×10⁻³ and the fractional S1–S3 by up to ~2×10⁻²; with
  ES3's truncation restored (as the golden tests do), the two agree to ~10⁻⁹.
- **K D1 spectra** move slightly near the K41 lines: ES3 gave K41 4²P1/2 the
  ground state's hyperfine constant (127.007 MHz for a measured
  15.245 MHz). The golden tests restore ES3's value for their K D1 case.
- **S1, S2, S3** are now fractional (divided by S0), with ES3's signs (S3 > 0
  for right-circular light, as in Born & Wolf; the alpha pre-releases flipped
  it). ES3 returned absolute intensity differences. Code that compared
  S1/S2/S3 to experimental data will need updating.
- **`nPlus`/`ChiPlus`/`phiPlus`/`TransLeft` (and `…Minus`/`TransRight`)** are
  now the σ⁺ (and σ⁻) quantities, as ES3's docstrings said; ES3 returned the
  other one under each name. See the `Plus`/`Minus` naming note under
  [Outputs](#outputs).
- **The propagation in oblique fields.** The default,
  `use_non_orthogonal=True`, solves for the two non-orthogonal modes of the
  vapour. ES3 projected the field onto them as if they were orthogonal,
  which is exact only in Faraday and Voigt geometry: in Briscoe et al.'s
  Regime I cell (Rb D2, 2 mm, 120 °C, 239 G) its S0 is off by up to 0.03
  at Btheta = 30°, 0.07 at 45° and 0.48 at 80° for x-polarised light (0.07,
  0.13 and 0.53 for the worst linear polarisation). `use_non_orthogonal=False`
  still gives ES3's method. Faraday and Voigt results are the same either
  way.
- The `BufferGas` parameters are optional; `GammaBuf`/`shift` still work as before.

### `scope2elecsus` / Pizzey *et al.* (2022) tutorial

The widely-used [`scope2elecsus`](https://github.com/durham-qlm/scope2elecsus)
notebook accompanying *Laser spectroscopy of hot atomic vapours: from 'scope to
theoretical fit* (D. Pizzey *et al.*, *New J. Phys.* **24**, 125001, 2022) runs
on v4 **unchanged**. It uses ElecSus only to fit and generate transmission
(`S0`) spectra via `EM.fit_data(..., data_type='S0')` and
`sp.get_spectra(..., outputs=['S0'])` — both preserved with identical
signatures, and `S0` of linear light in a zero or longitudinal field is v3's
apart from the exact line frequencies (a difference of ≲ 3×10⁻³ in S0; see
above), so fitted parameters agree with v3's to within a correspondingly
small amount. None of the other breaking changes above apply (the tutorial
does not use Stokes parameters, oblique fields, the non-orthogonal path, or
the buffer-gas model). Just install v4 into the notebook's
environment; no code edits are required.

---

## Tests and validation

### Validate your installation (one command)

```bash
python -m elecsus.validate
```

The suites run under pytest, which a plain `pip install elecsus4` does not
install: use `pip install "elecsus4[test]"` (the command says so if pytest is
missing). This runs the packaged validation suites and writes overlay figures to
`./elecsus_validation/` showing spectra computed by **your installation**
on top of stored reference output from genuine ElecSus v3.0.8, for the
classic published-figure configurations (Siddons 2008, Zentile 2015,
Keaveney 2016, Rotondaro 2015). You see the agreement, not just an exit code.

It also runs the **convention-free physics suites**, whose expected values
come from first principles and the literature only, never from ElecSus's own
output (constants and sources in `elecsus/tests/literature.py`):

- `test_sum_rule.py`: the integrated absorption ∫ k Im(Tr χ/3) dν equals
  N π r_e c f, with the oscillator strength f from the measured lifetime, for
  every element and line, and does not change with the field, the Doppler
  temperature, broadening, shifts or the buffer gas; the optically thin limit
  ∫ (1 − S0) dν = N L π r_e c f for any polarisation and field direction;
- `test_kramers_kronig.py`: Re χ is the Hilbert transform of Im χ;
- `test_hyperfine_levels.py`: the ground levels follow Breit–Rabi, and the
  excited hyperfine levels match the measured A and B constants;
- `test_symmetries.py`: reversing B gives the mirror image, B = 0 is
  isotropic, S0 ≤ 1, and a pure state stays fully polarised.

A bug these find is marked as a strict expected failure (`xfail`) until it is
fixed, so the suite stays green and the fix has to remove the mark.

The core of this is a **golden-reference regression suite**
(`elecsus/tests/test_es3_golden.py`): 20 configurations spanning
Faraday/Voigt/oblique geometries, Rb/Cs/K on both D-lines, isotopically pure
cells, buffer-gas broadening and rotated polarisations, generated by running
the final release of v3 (3.0.8: commit `1933456`, tag `v3.0.8`; see
`elecsus/tests/golden/generate_es3_reference.py` for provenance). The tests
assert that, with v3's whole-MHz line frequencies, its K41 P₁/₂ constant and
its 273.16 in the Boltzmann populations restored for the comparison (v4 keeps
the exact values; `test_line_frequencies.py` bounds the difference), and with
v4 at the mirrored angles (π − `Btheta`, −`Bphi`) of the geometry paragraph
above:

- `S0`, `Ix`, `Iy` agree with v3 to 3×10⁻⁹ (8×10⁻⁶ in the two oblique cases,
  where v4's analytic solver replaces v3's SymPy solve);
- `S1`, `S2`, `S3` equal the v3 values divided by S0, with v3's sign of S3
  (Born & Wolf's), to 3×10⁻⁸ (10⁻⁴ in the oblique cases) at every detuning
  where the vapour transmits.

### Developer test commands

```bash
# Run the full auto-discovered pytest suite (a couple of minutes)
python -m pytest

# Non-orthogonal propagation against Briscoe et al. 2025 theory
python -m pytest elecsus/tests/briscoe_2025_theory_repro.py -v

# Cascade cell validation (identity, combination rules, polarisation coupling)
python -m pytest elecsus/tests/test_cascade.py -v

# Buffer gas model vs Alqarni et al. 2025 (K D1 / Ne)
python -m pytest elecsus/tests/test_buffer_gas.py -v

# Fitting against a BufferGas cell (T/Bfield recovery, GammaBuf-vary guard)
python -m pytest elecsus/tests/test_fitting_buffer_gas.py -v
```

`elecsus/tests/` also contains paper-figure-reproduction scripts (Rotondaro
2015, Siddons/Weller/Zentile/Keaveney spectroscopy examples, the field-gradient
fit from the 2017 ElecSus paper) converted to fast, non-blocking pytest checks.
Beyond physical-sanity assertions, these now pin quantitative anchors (peak
transmissions, line-centre depths) to the v3.0.8 golden reference, so a physics
regression fails loudly instead of needing a human to eyeball a plot. Set `ELECSUS_SHOW_PLOTS=1` before running
pytest to additionally pop up the original matplotlib figures interactively.
A couple of slow/data-dependent entry points (a ~15-minute differential-evolution
fit requiring an external CSV not in this repo, and a multi-algorithm fitting
comparison) are kept as manual, non-test-prefixed functions rather than
automated tests - see the module docstrings in `test_fieldgradient.py` and
`test_fitting.py`.

---

## Plot style

All scripts and notebooks in this repository use a shared house style defined in
[`elecsus/libs/plot_style.py`](elecsus/libs/plot_style.py).

```python
from elecsus.libs.plot_style import apply_style
apply_style()          # 14 pt — notebooks and standalone scripts
apply_style('paper')   # 11 pt — paper-width (≈ 6.7 inch) figures
```

This sets:

- **Typography** — Times New Roman throughout, including mathtext (`$...$`)
- **Colour cycle** — Durham University palette: Palatinate Purple → Durham Blue → Red → Yellow → Grey → Black
- **Resolution** — 300 dpi on save; 150 dpi for on-screen display
- **Font embedding** — `pdf.fonttype = 42` ensures fonts are embedded when saving to PDF

Individual colours are exported directly for use in explicit `color=` arguments:

```python
from elecsus.libs.plot_style import d_purple, d_blue, d_red, d_yellow, d_grey, d_black
```

The full colour definitions live in [`elecsus/libs/durhamcolours.py`](elecsus/libs/durhamcolours.py).

---

## Citation

If you use ElecSus v4, please cite the software, the original ElecSus papers,
and the v4 paper when it is available:

- **The software (ElecSus 4, every version):** S. A. Wrathmall, R. M. Potvliege,
  I. G. Hughes *et al.*, *ElecSus*, Zenodo. https://doi.org/10.5281/zenodo.21509735
  (each version also has its own DOI, on its Zenodo page).

- **ES3:** J. Keaveney et al., *Comput. Phys. Commun.* **224**, 311–324 (2018).
  https://doi.org/10.1016/j.cpc.2017.12.001

- **ES1/ES2:** M. A. Zentile et al., *Comput. Phys. Commun.* **189**, 162–174 (2015).
  https://doi.org/10.1016/j.cpc.2014.11.023

- **ES4:** Wrathmall et al., *Comput. Phys. Commun.* (in preparation, 2026).

The non-orthogonal propagation method is described in:

- J. D. Briscoe et al., *Opt. Express* **33**, 32099–32114 (2025).
  https://doi.org/10.1364/OE.567332

---

## License

Distributed under the Apache License, Version 2.0. See `LICENSE.txt`.
