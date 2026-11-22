# Changelog

ElecSus 4 is installed from PyPI as `elecsus4` and imported as `elecsus`.
Until 4.0.0 its versions are pre-releases: betas (4.0.0b1, 4.0.0b2, …), then
release candidates. **Breaking** marks a change that can make existing code
fail or give different results.

## 4.0.0b1 — beta

The first beta of ElecSus 4, and its first release on PyPI. Since 4.0.0a1 the
physics, the fitting and the new GUI have been tested against first
principles, the literature and measured spectra (Briscoe et al. 2025, Rb D2 in
an oblique field; Alqarni et al. 2025, K D1 with neon), and the errors those
tests found have been fixed. Many of them are in ElecSus 3 too. The interface
is now frozen: until 4.0.0, a change that alters a result or removes a
parameter will be listed here, with its reason.

Coming from ElecSus 3? Read *Results that change* below, and *Backwards
compatibility* in the README.

### Results that change

- **The geometry is the documented one.** **Breaking** wherever a sense
  matters. ElecSus 3, and 4.0.0a1, computed the mirror image (y → −y) of the
  geometry they document. **B** is now along (sin θ cos φ, sin θ sin φ,
  cos θ), with `Btheta` = θ measured from +z, the direction of the light, and
  `Bphi` = φ from x towards y. ElecSus 3's result for (`Btheta`, `Bphi`) is
  ElecSus 4's for (π − `Btheta`, −`Bphi`).
  - A Faraday rotation has the opposite sense, and in a transverse field
    with `Bphi` ≠ 0, light polarised along the field's azimuth is now the π
    mode.
  - Transmission of linear light in a zero or longitudinal field is
    unchanged.
  - Briscoe et al.'s measured spectra are reproduced with the input
    polarisation of the experiment, [cos θ_E, sin θ_E]; 4.0.0a1 needed
    [cos θ_E, −sin θ_E]. (ElecSus 3, whose oblique propagation was also
    wrong, is closer with [cos θ_E, sin θ_E], but misses the data either way.)
- **`Plus` and `Minus` name σ⁺ and σ⁻.** **Breaking.** `ChiPlus`, `nPlus`,
  `phiPlus` and `TransLeft` are the σ⁺ quantities, as ElecSus 3's docstrings
  said, and `…Minus` and `TransRight` the σ⁻ ones; ElecSus 3 returned each
  under the other's name. `alphaPlus` and `alphaMinus` are unchanged. So
  does `elecsus.libs.spectra.calc_chi`: it returns (σ⁺, σ⁻, π), where
  ElecSus 3 returned the σ⁻ susceptibility first.
- **Fractional Stokes parameters.** **Breaking.** `S1`, `S2` and `S3` are
  divided by `S0`, so they lie in [−1, 1]; ElecSus 3 returned intensity
  differences. `S3` has ElecSus 3's sign, Born & Wolf's: `S3` > 0 for
  right-circular light, so `[1, 1j]/√2` light has `S3` = −1. 4.0.0a1 had
  reversed it.
- **Oblique fields: the non-orthogonal modes are the default.** **Breaking**
  for `Btheta` other than 0, π/2 and π. `use_non_orthogonal` now defaults to
  `True`. ElecSus 3 projected the light onto the vapour's two modes as if they
  were orthogonal, which is exact only in Faraday and Voigt geometry. In an
  oblique field its `S0` is wrong (Rb D2, 2 mm, 120 °C, 239 G) by up to
  0.03 at 30°, 0.07 at 45° and 0.48 at 80° for x-polarised light, and by up
  to 0.07, 0.13 and 0.53 for the worst linear polarisation. Faraday and Voigt results are
  unchanged, and `use_non_orthogonal=False` still gives ElecSus 3's method.
- **Exact line frequencies.** ElecSus 3 truncated every line's frequency to
  whole MHz. Every spectrum changes slightly (`S0` by up to ~3×10⁻³, the
  fractional `S1`–`S3` by up to ~2×10⁻²), and ML fits of `Bfield` no longer
  stall near their starting value.
- **⁴¹K D1.** The 4²P₁/₂ hyperfine constant of ⁴¹K was the ground state's
  (127.007 MHz, for a measured 15.245 MHz), so its D1 lines were up to
  ~140 MHz out. Natural potassium D1 spectra change slightly near them.
- **The reversed pass absorbs.** Since ElecSus 3.0.0, `wavevector_dirn=-1`
  made the vapour amplify light where it should absorb it (a transmission of
  ~10⁶² on resonance; the Faraday rotation came out right). The reversed pass
  is now a forward pass in the frame turned by π about x, and a round trip
  equals one pass through twice the length. `wavevector_dirn` takes only 1
  or −1.
- **B along −x** (`Btheta` = 3π/2 or −π/2) gave `S0` ~ 10¹⁰ in 4.0.0a1
  (ElecSus 3 was not affected). Angles a whole turn apart now give the same
  spectrum.
- **Buffer-gas data.** Each of the 58 entries was checked against its source,
  now named in `elecsus/libs/BufferGasData.py`. Since 4.0.0a1, eleven
  changed: broadenings by up to 15 % (Cs D1/He) and shifts by up to 24 %
  (Rb D2/Kr).
- **Buffer-gas model.** **Breaking.** `shift` is now added to the model's
  pressure shift, as an offset of the frequency axis; the model used to
  replace it. With `Constrain=False` the model is evaluated at `DoppTemp`, the
  temperature of the gas in the beam, not at `T` (2.8 % more broadening and
  shift for an 80 °C stem and a 100 °C body).
- **Cascades.** With `polariser_angle` set, the instrument-total `S3` of
  `calculate_cascade` had 4.0.0a1's sign; it now has `calculate`'s. Where an opaque first cell lets no
  light through, the second cell's outputs are 0, not NaN. **Breaking:**
  `reset_polarisation` needs `polariser_angle`. Without a polariser it made a
  circular state of the same intensity, which no passive optic can do.
- **Smaller changes:**
  - one Kelvin offset, 273.15 (the Boltzmann populations used 273.16), which
    moves `S0` by ~10⁻⁸;
  - a field below 0.01 G at an oblique angle is no longer taken to lie along
    ±z, which made its Faraday rotation too large, the more so the nearer
    the field is to transverse (6 times at 85°, 300 times at 89.9°);
  - `phi` has the sense of the measured rotation when **B** is along −z
    (4.0.0a1 had the opposite);
  - the light leaving the cell on the non-orthogonal path is transverse
    (`S0` counted the modes' longitudinal part, up to ~2×10⁻⁶ of it);
  - with `Constrain=False` and no `DoppTemp`, the Doppler width is taken at
    `T`; it was taken at 20 °C.

### Added

- **A browser GUI**: `python -m elecsus.elecsus_dash_gui`, after
  `pip install "elecsus4[gui]"`. Spectra, fits to your data or to a bundled
  sample spectrum, two-cell filters and energy-level diagrams, each with CSV
  export; plots against detuning or wavelength. It works offline.
- **Energy levels and transitions** (`elecsus.energy_levels`). For one isotope
  in a field: the levels, labelled by their |mJ, mI⟩ components and the
  |F, mF⟩ state each joins at zero field; the σ⁺, σ⁻ and π lines, with the
  detunings and strengths the spectra use; the levels against field; and
  `level_diagram`, which draws it all above the spectrum, with a stick
  spectrum below. `python -m elecsus.energy_levels ISOTOPE LINE B --out FILE`
  draws it and prints the lines.
- **Fitting the buffer-gas pressure**: `fit_data` can vary `BufferPressure`,
  with the buffer-gas model on.
- **Six tutorial notebooks** in `Notebooks/`, with their output.
- **Checks of the physics that take nothing from ElecSus**: the
  integrated-absorption sum rule against the measured oscillator strengths,
  Kramers–Kronig, Breit–Rabi and the measured hyperfine constants, and the
  symmetries (B → −B, B = 0, S0 ≤ 1, purity). `python -m elecsus.validate`
  runs them, with the comparison against ElecSus 3.
- **Checked parameters.** An unknown `p_dict` key raises `ValueError`, naming
  the closest valid one (`strict=False` makes it a warning): a mistyped
  `n_window` or `use_non_orthogonal` used to switch that feature off in
  silence. A value outside the model's range (a negative `lcell`, an isotope
  fraction outside 0–100, a temperature below absolute zero, `n_window` below
  1, an unsupported element or line) raises before anything is computed.
- `elecsus.__version__`.

### Changed

- **Packaging.** ElecSus 4 is on PyPI as `elecsus4`, and is still imported as
  `elecsus`. It is built from `pyproject.toml`, with the extras `gui`, `test`
  and `dev`, and needs Python 3.10, NumPy 1.24, SciPy 1.10, matplotlib 3.7
  and lmfit 1.3.4 or newer. `sympy` and `psutil` are no longer dependencies.
- **Breaking:** `elecsus/libs/` no longer puts itself on `sys.path`, where it
  registered `spectra`, `EigenSystem`, `AtomConstants` and eight other modules
  under top-level names. Import them as `from elecsus.libs import spectra`.
- Intensity-like outputs (`S0`–`S3`, `Ix`, `Iy`, `I_P45`, `I_M45`, `Il`, `Ir`,
  `psi`, `DoLP`) are `float64`. They were `complex128`, with an imaginary part
  of ~10⁻¹⁷.
- `n_window` on its own applies the window model; with
  `use_non_orthogonal=False` it warns that the windows are not modelled.
- `fit_data` refuses to vary what the data cannot determine, naming it: both
  `E_x` and `E_y` (only their ratio is observable), `BufferTemp_fill`
  (degenerate with the pressure), `BufferPressure` with the buffer-gas model
  off, and `GammaBuf` with it on.

### Fixed

Most of these are in ElecSus 3 too.

- **Fits use the model you ask for.** The ML fit, which is the default and the
  last stage of differential evolution, dropped `use_non_orthogonal`,
  `n_window` and `wavevector_dirn`, so it fitted the orthogonal, windowless
  model (a temperature 7 °C out in an oblique field). Random restart's last
  fit was a Rb D2 fit without the buffer gas, whatever was asked for.
  `BoltzmannFactor` now reaches the fit model.
- **Fits stay in the physical range.** Search windows, `p_dict_bounds` and
  random restart's starts are kept inside each parameter's range, so a
  negative width or pressure no longer stops a fit. An unphysical starting
  value is refused by name; a random-restart start whose fit leaves the range
  is dropped, with a warning.
- **Fits keep what they find for the input polarisation** (`E_x`, `E_y`,
  `E_phase`). After differential evolution, and in random restart, the last
  fit started again from the original `E_in`.
- **`fit_data` rejects what it used to ignore**: a misspelt keyword (so
  `datatype='S1'` fitted the S0 model), an unknown `fit_algorithm` (which ran
  random restart), and an unknown key in `p_dict`, `p_dict_bools` or
  `p_dict_bounds` (`strict=False` makes the last a warning, as in
  `calculate`).
- **`p_dict_bounds` reaches random restart** as the width of its search, as
  documented; it was ignored. Each entry is checked in the form its algorithm
  takes.
- **`fit_data`'s arguments and result.** `data` is `[x, y]`, 2×N, as the code
  always read it: the docstrings said N×2, and an N×2 array was fitted as two
  points. Anything else raises `ValueError`. `E_in` may be the Jones vector
  `[Ex, Ey, 0]` that `calculate` takes, and lmfit's bounds no longer turn
  right-circular light into linear light. `Elem`, `Dline`, `Bfield`, `T` and
  `lcell` may be left out, as in `calculate`. The best parameters returned can
  be passed straight back to `calculate`.
- **Parameters that crashed or were ignored** now raise or warn, by name:
  `Elem='Ideal'`, which has no susceptibility model; `K40frac` + `K41frac`
  above 100; switches given as strings (`'False'` switched the non-orthogonal
  solver on); `E_x`, `E_y` or `E_phase` given to `calculate`, which takes
  `E_in`; part of an ElecSus 3 input field (`E_x` without `E_y` and `E_phase`)
  in a fit's `p_dict`; and `BufferPressure` without `BufferGas`.
- A single detuning works: `calculate(5.0, ...)`. `calculate_cascade` takes
  `strict`.
- Importing ElecSus no longer silences every warning in your session,
  ElecSus's own included. Its warnings point at the line that called it.
- `python -m elecsus.validate` says when pytest is missing.
- The plot style falls back to another serif font where Times New Roman is
  not installed, instead of printing a warning for every label it draws.
- `RotationMatrices` and `JonesMatrices` no longer use the deprecated
  `np.matrix`.

### Removed

- **Simulated annealing** (`fit_algorithm='SA'`) raises `ValueError`. It had
  not run on Linux or macOS since ElecSus 3.0.0, and differential evolution
  (`'DE'`) does the same job faster.
- **ElecSus 3's wxPython GUI.** The browser GUI replaces it.

## 4.0.0a1 — alpha

The first ElecSus 4 pre-release, on GitHub only.

### Added

- **Non-orthogonal propagation** in oblique fields, validated against
  Briscoe et al. (2025): `use_non_orthogonal=True`, the default since 4.0.0b1.
- **An analytic solver for the modes in oblique fields**, ~130× faster than
  ElecSus 3's SymPy path, with corrected eigenmodes.
- **An optional model of the cell windows**, four Fresnel interfaces, set by
  `n_window`. The default, `None`, is a cell without windows, which reproduces
  Briscoe et al.'s Eq. 9; `n_window=1.51` models borosilicate glass.
- **A buffer-gas model**: `BufferGas`, `BufferPressure` and `BufferTemp_fill`,
  validated on K D1 with neon against Alqarni et al. (2025).
- **Two cells in series**: `calculate_cascade()`, with an optional polariser
  between them (magneto-optical filters).
- **New outputs**: the Stokes azimuth `psi` and the degree of linear
  polarisation `DoLP`; the refractive indices n(σ±), the Faraday rotation φ
  and α_π; and the Poynting-flux transmission `S0_P`.
- **Validation**: a regression test against stored output of ElecSus 3.0.8,
  `python -m elecsus.validate`, and continuous integration.

### Changed

- **Breaking:** `S1`, `S2` and `S3` are fractional (divided by `S0`), and the
  sign of `S3` was reversed; 4.0.0b1 restores ElecSus 3's.
