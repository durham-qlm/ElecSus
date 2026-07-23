# ElecSus v4

**ElecSus** calculates the electric susceptibility, refractive index, and transmission
spectra of thermal alkali-metal vapour cells on the D-lines. It includes least-squares
fitting routines for extracting experimental parameters from spectroscopic data.

This is **version 4**, a significant extension of
[ElecSus v3](https://doi.org/10.1016/j.cpc.2017.12.001) (Keaveney et al. 2018).
A companion paper describing the new physics is in preparation for
*Computer Physics Communications*.

> ⚠️ **`v4.0.0-alpha` — pre-release.** APIs may still change before the stable
> v4.0.0 release; the current stable version is v3.

> The v3 API (`calculate`, `fit_data`) is fully preserved. Existing scripts run
> unchanged unless they relied on the absolute (un-normalised) values of S1, S2, S3
> — see [Backwards compatibility](#backwards-compatibility).

---

## What's new in v4

| Feature | Description |
|---|---|
| **Non-orthogonal propagation** | Correct biorthogonal mode decomposition for oblique magnetic fields (Briscoe et al., *Opt. Express* **33**, 32099, 2025). Enabled with `use_non_orthogonal=True`. |
| **Fast oblique solver** | Webber's fully vectorised analytic eigenvalue/eigenvector solution replaces the ES3 SymPy + Python for-loop: ~130–200× faster, no unphysical artefacts at high fields. |
| **Dynamic buffer gas model** | `BufferGas`, `BufferPressure`, `BufferTemp_fill` parameters compute temperature- and pressure-dependent broadening and shift automatically. |
| **Cell-window Fresnel correction** | Full 4-interface model (air → glass → vapour → glass → air) applied per eigenmode. Controlled by `n_window` (default `None`, off); set e.g. `1.51` for borosilicate glass to enable. |
| **Cascaded cells** | `calculate_cascade()` — passes the output electric field of cell 1 as the input to cell 2, correctly propagating polarisation state and intensity. |
| **Fractional Stokes parameters** | S1, S2, S3 now normalised by S0 (lying in [−1, +1]) and S3 sign-corrected. |

---

## Documentation

This README is the quick-start reference. The API and parameter conventions
(units, sigma labelling, non-orthogonal propagation, window Fresnel model,
buffer gas, Stokes conventions, cascades/MOFs) are documented in the module
docstrings and exercised by the test suite. A full user and theory manual is
in preparation for a later release.

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/durham-qlm/ElecSus.git
cd ElecSus

# 2. Create and activate the conda environment
#    (installs Python 3.12, NumPy, SciPy, lmfit, Dash, and all other dependencies)
conda env create -f environment.yml
conda activate elecsus4

# 3. Install ElecSus in editable mode
pip install -e .
```

**Existing environment?**  If you already have an ElecSus conda environment (it may
be named `elecsus` rather than `elecsus4`), activate it and install the v4 GUI
dependencies manually:

```bash
conda activate elecsus          # or whatever your env is called
pip install "dash>=2.17" "dash-bootstrap-components>=1.6" "plotly>=5.22"
pip install -e .
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
    'Bfield': 300,    'Btheta': np.pi/4,   # radians from k-vector
    'use_non_orthogonal': True,        # Briscoe et al. 2025 formalism
    'n_window': 1.51,                  # enable borosilicate windows (default None = off)
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
the polarisation state and intensity attenuation. The combination rules for all
output types are verified in the cascade tests under `elecsus/tests/`.

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
`reset_polarisation` (`'left_cp'` / `'right_cp'`), on its own or combined with
`polariser_angle`:

```python
(S0_1,), (S0_total,) = em.calculate_cascade(
    detuning, [1, 0, 0], p1, p2, outputs=['S0'],
    polariser_angle=np.pi/2, reset_polarisation='left_cp')
```

This also puts cell-2 outputs in instrument-totals mode; in this mode the total
equals the simple product of two independent single-cell calculations.

### Fitting

```python
import numpy as np
from elecsus import elecsus_methods as em

raw = np.loadtxt('my_data.csv', delimiter=',')   # (N, 2) array: detuning, signal
data = [raw[:, 0], raw[:, 1]]                    # fit_data expects [x_array, y_array]

p0   = {'Elem': 'Rb', 'Dline': 'D2', 'T': 20, 'lcell': 75e-3, 'Bfield': 0}
vary = {'T': True, 'Bfield': False, 'lcell': False}

best_params, rms, result = em.fit_data(data, p0, vary, data_type='S0', fit_algorithm='ML')
# best_params  — dict of fitted parameter values
# rms          — RMS residual between data and best fit
# result       — lmfit MinimizerResult (result.best_fit, result.fit_report(), etc.)
```

`fit_algorithm` accepts `'ML'` (Marquardt-Levenberg, default), `'RR'` (random-restart),
`'SA'` (simulated annealing), or `'differential_evolution'`.

If your cell was filled with a buffer gas, pass `BufferGas`/`BufferPressure`/
`BufferTemp_fill` in `p0` as usual — they are fixed cell-fabrication properties
known from how the cell was made, so leave them out of `vary`. GammaBuf/shift
are then derived from the buffer gas model at each fit iteration rather than
independently fitted; setting `vary['GammaBuf']` or `vary['shift']` to `True`
while `BufferGas` is set raises a `ValueError`, since the buffer gas model
would silently override them anyway.

---

## Notebooks

Tutorial Jupyter notebooks — basic transmission, magnetic field & polarimetry,
non-orthogonal propagation, buffer gas, cascaded cells, and fitting — are being
finalised and will ship in the next update.

---

## Key parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `Elem` | str | `'Rb'` | Element: `'Rb'`, `'Cs'`, `'K'`, `'Na'`, `'Li'` |
| `Dline` | str | `'D2'` | Transition: `'D1'` or `'D2'` |
| `T` | float | `20` | Temperature (°C) — sets number density |
| `lcell` | float | `75e-3` | Cell length (m) |
| `Bfield` | float | `0` | Magnetic field magnitude (Gauss) |
| `Btheta` | float | `0` | Angle of B from k-vector (**radians**; the GUI inputs are in degrees and convert internally) |
| `Bphi` | float | `0` | Azimuthal angle of B (**radians**; the GUI inputs are in degrees and convert internally) |
| `GammaBuf` | float | `0` | Extra phenomenological Lorentzian broadening (MHz). Adds to the natural linewidth. **Silently overridden** by the buffer gas model when `BufferGas` is set — do not use both. |
| `shift` | float | `0` | Global frequency shift of all resonance lines (MHz). **Silently overridden** by the Pitz pressure-shift when `BufferGas` is set — do not use both. |
| `BufferGas` | str | — | Buffer gas species: `'Ne'`, `'Ar'`, `'Kr'`, `'Xe'`, `'N2'`, `'He'`. Activates the Pitz 2014 dynamic broadening and shift model. **When set, `GammaBuf` and `shift` are ignored.** |
| `BufferPressure` | float | — | Fill pressure (Torr) — the pressure at which the sealed cell was filled at room temperature. The actual pressure at operating temperature is computed internally. |
| `BufferTemp_fill` | float | `293.15` | Temperature (K) at which the cell was sealed. Combined with `BufferPressure` to fix the number of buffer gas atoms. Default is 293.15 K (20 °C). |
| `use_non_orthogonal` | bool | `False` | Enable Briscoe et al. biorthogonal propagation |
| `n_window` | float or None | `None` | Window glass index for the 4-interface Fresnel model; `None` (default) disables it |
| `rb85frac` | float | `72.17` | ⁸⁵Rb fraction (%) |
| `DoppTemp` | float | — | Independent Doppler temperature (°C) |
| `Constrain` | bool | `True` | Force DoppTemp = T |

All parameters have defaults; only those you need to change must be specified.

---

## Outputs

Pass a list of strings to the `outputs` keyword of `calculate()`:

| Name | Description |
|---|---|
| `S0` | Total transmitted intensity (normalised to input) |
| `S1` | Fractional Stokes: (Ix − Iy) / S0 |
| `S2` | Fractional Stokes: (I₊₄₅ − I₋₄₅) / S0 |
| `S3` | Fractional Stokes: circular degree of polarisation |
| `Ix`, `Iy` | Intensities through x- and y-polarisers |
| `I_P45`, `I_M45` | Intensities through ±45° polarisers |
| `Il`, `Ir` | Circular component intensities (`[1,-i]/√2` and `[1,+i]/√2`) |
| `alphaPlus`, `alphaMinus`, `alphaZ` | Absorption coefficients (σ⁺, σ⁻, π) |
| `ChiPlus`, `ChiMinus`, `ChiZ` | Complex susceptibilities of the σ⁻, σ⁺, π manifolds — see naming note below |
| `nPlus`, `nMinus` | Complex refractive indices √(1+χ) of the two circular modes (Faraday geometry only) |
| `TransLeft`, `TransRight` | Single-mode σ⁻ / σ⁺ transmissions (Faraday geometry only) |
| `phi`, `phiPlus`, `phiMinus` | Net Faraday rotation and the two single-mode phases (Faraday geometry only; **dispersion-only** — see `psi`) |
| `psi` | Polarisation-ellipse azimuth ½·atan2(S2, S1) from the full field propagation — the rotation angle a polarimeter actually reports. Equals `phi` (mod π) in Faraday geometry; remains correct in *any* geometry, where `phi` does not. Principal value in (−π/2, π/2]; pair with `DoLP` for validity |
| `DoLP` | Degree of linear polarisation √(S1²+S2²) ∈ [0, 1]. Validity diagnostic for any rotation angle: DoLP → 0 where dichroism has circularised the light and rotation is unmeasurable |
| `E_out` | Output electric field vector, shape (3, N) |

If `outputs` is omitted, the default set (S0, S1, S2, S3, Ix, Iy, I_P45, I_M45,
alphaPlus, alphaMinus, alphaZ) is returned as a tuple.

### ⚠️ The `Plus`/`Minus` naming convention (an ES3 labelling error, corrected in v4)

The historical output names `ChiPlus` / `nPlus` / `phiPlus` (and `…Minus`)
carry a **circular-basis** suffix that is **opposite to the σ transition** they
describe:

| Output key | Is physically the… | Governs |
|---|---|---|
| `ChiPlus`, `nPlus`, `phiPlus` | **σ⁻** (Δm = −1) manifold | `[1,+i]/√2` light |
| `ChiMinus`, `nMinus`, `phiMinus` | **σ⁺** (Δm = +1) manifold | `[1,-i]/√2` light |

**What ES3 got wrong.** The computed spectra were correct — v4 reproduces them
exactly (golden tests) — but ES3's docstrings/comments described these
*secondary* outputs with the wrong σ (e.g. "`nPlus` = refractive index driving
σ⁺ transitions"). So an ES3 user who requested `nPlus` or `ChiPlus` **by name**,
expecting σ⁺, actually received the σ⁻ quantity. ES3 was also internally
inconsistent: `alphaPlus`/`alphaMinus` quietly cross over
(`alphaPlus = 2·Im(nMinus)·k`), so *they* alone carried the correct physical σ.

**What v4 does.** The true mapping was verified two independent ways — the *L*ᵤ
(`lz`) operator basis ordering **and** the high-field Zeeman-shift direction —
and is pinned by regression tests (`elecsus/tests/test_faraday_outputs.py`).
All code comments and docstrings are corrected. The **API key names are kept
unchanged** for backward compatibility, but the **GUI labels every output by
the physical σ** — the GUI shows `Trans σ⁺`, `n(σ⁺)`, `α σ⁺` … which map to the
API keys `TransRight`, `nMinus`, `alphaPlus`.

---

## Web GUI

The interactive GUI runs in any browser and requires no extra dependencies beyond
the standard install.

```bash
# 1. Activate the conda environment
#    (fresh install from environment.yml → elecsus4;
#     older installs may be named elecsus — use whichever name you have)
conda activate elecsus4

# 2. Launch the GUI
python -m elecsus.elecsus_dash_gui

# 3. Open your browser at
#    http://127.0.0.1:8050
```

The GUI dependencies (`dash`, `dash-bootstrap-components`, `plotly`) are already
included in `environment.yml`.  To stop the server press **Ctrl-C** in the terminal.

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

- Single-cell forward calculation for all output types (S0, S1, S2, S3, Ix, Iy, …)
- Fit-to-data on the single-cell tab (CSV upload, ML or differential-evolution),
  with **Export fit CSV** (detuning, data, best fit, residual)
- Cascaded two-cell filters (same species in both cells, element/D-line locked
  to cell 1) with an optional intermediate linear polariser — default crossed
  at 90°, the magneto-optical-filter (MOF) topology for laser input light.
  The plot shows the spectrum after cell 1 and the total after both filters.
- Buffer gas parameters (`BufferGas`, `BufferPressure`, `BufferTemp_fill`)
- Interactive plots with zoom/pan and CSV export of the current spectrum
- Dark/light theme is fixed to the Durham light theme

---

## Backwards compatibility

- `calculate()` and `fit_data()` signatures are unchanged.
- **S1, S2, S3** are now fractional (divided by S0) and S3 is sign-corrected.
  ES3 returned absolute intensity differences. Code that compared S1/S2/S3 to
  experimental data will need updating.
- **`nPlus`/`ChiPlus`/`phiPlus` (and `…Minus`)** are the σ∓ circular-basis
  quantities, *opposite* to what ES3's docstrings claimed — see the
  `Plus`/`Minus` naming note under [Outputs](#outputs). The values match ES3
  exactly; only the (now corrected) labels differ. Any ES3 code that read
  `nPlus`/`ChiPlus` by name expecting σ⁺ was reading the σ⁻ quantity.
- `use_non_orthogonal=False` (default) restores the ES3 propagation path exactly.
- The `BufferGas` parameters are optional; `GammaBuf`/`shift` still work as before.

### `scope2elecsus` / Pizzey *et al.* (2022) tutorial

The widely-used [`scope2elecsus`](https://github.com/durham-qlm/scope2elecsus)
notebook accompanying *Laser spectroscopy of hot atomic vapours: from 'scope to
theoretical fit* (D. Pizzey *et al.*, *New J. Phys.* **24**, 125001, 2022) runs
on v4 **unchanged**. It uses ElecSus only to fit and generate transmission
(`S0`) spectra via `EM.fit_data(..., data_type='S0')` and
`sp.get_spectra(..., outputs=['S0'])` — both preserved with identical
signatures, and `S0` is pinned bit-for-bit to ES3 v3.0.8 by the golden-reference
suite, so fitted parameters reproduce exactly. None of the breaking changes
above apply (the tutorial does not use Stokes parameters, the non-orthogonal
path, or the buffer-gas model). Just install v4 into the notebook's
environment; no code edits are required.

---

## Tests and validation

### Validate your installation (one command)

```bash
python -m elecsus.validate
```

This runs the packaged validation suites and writes overlay figures to
`./elecsus_validation/` showing spectra computed by **your installation**
on top of stored reference output from genuine ElecSus v3.0.8, for the
classic published-figure configurations (Siddons 2008, Zentile 2015,
Keaveney 2018, Rotondaro 2015). You see the agreement, not just an exit code.

The core of this is a **golden-reference regression suite**
(`elecsus/tests/test_es3_golden.py`): 20 configurations spanning
Faraday/Voigt/oblique geometries, Rb/Cs/K on both D-lines, isotopically pure
cells, buffer-gas broadening and rotated polarisations, generated by running
the final upstream v3.0.8 release (commit `065fa28`; see
`elecsus/tests/golden/generate_es3_reference.py` for provenance). The tests
assert that:

- `S0`, `Ix`, `Iy` are **unchanged** from v3 (agreement ≲ 10⁻⁹; oblique cases
  ≲ 10⁻⁵, where v4's analytic solver replaces v3's SymPy solve);
- `S1`, `S2`, `S3` match the v3 values under exactly the documented v4
  renormalisation (divided by S0, with the S3 sign convention corrected) at
  every detuning where the vapour transmits.

### Developer test commands

```bash
# Run the full auto-discovered pytest suite (~65 tests, a few seconds)
python -m pytest

# Non-orthogonal propagation against Briscoe et al. 2025 theory
python -m pytest elecsus/tests/briscoe_2025_theory_repro.py -v

# Cascade + buffer-gas validation with figures (one command)
python -m elecsus.validate

# Buffer gas model tests (K D1 / Ne vs Alqarni et al. 2025)
python -m pytest test_buffer_gas.py -v

# Fitting against a BufferGas cell (T/Bfield recovery, GammaBuf-vary guard)
python -m pytest test_fitting_buffer_gas.py -v
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

All plotting in this package uses a shared house style defined in
[`elecsus/libs/plot_style.py`](elecsus/libs/plot_style.py).

```python
from elecsus.libs.plot_style import apply_style
apply_style()          # 14 pt — default size
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

If you use ElecSus v4, please cite the original ElecSus papers and the v4 paper
(when available):

- **ES3:** J. Keaveney et al., *Comput. Phys. Commun.* **224**, 64–75 (2018).
  https://doi.org/10.1016/j.cpc.2017.12.001

- **ES1/ES2:** M. A. Zentile et al., *Comput. Phys. Commun.* **189**, 162–174 (2015).
  https://doi.org/10.1016/j.cpc.2014.11.023

- **ES4:** Wrathmall et al., *Comput. Phys. Commun.* (in preparation, 2026).

The non-orthogonal propagation method is described in:

- T. Briscoe et al., *Opt. Express* **33**, 32099–32114 (2025).
  https://doi.org/10.1364/OE.553296

---

## License

Distributed under the Apache License, Version 2.0. See `LICENSE.txt`.
