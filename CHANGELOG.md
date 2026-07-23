# Changelog

## 4.0.0a1 — pre-release

First public v4 (alpha) release. v4 is a substantial extension of ElecSus 3
(Keaveney et al. 2018) with new physics, a new GUI, and a full regression
suite. As a pre-release, APIs may still change before 4.0.0.

### Added
- **Analytic oblique (non-orthogonal) mode solver** — ~130× faster than the ES3
  SymPy oblique path, with corrected eigenmodes; validated against Briscoe
  et al. (2025) experimental data. Enable with `use_non_orthogonal=True`.
- **Optional cell-window Fresnel model** — a four-interface treatment controlled
  by `n_window`. The default is `None` (no window correction), which reproduces
  the published Briscoe et al. Eq. 9; set e.g. `n_window=1.51` for borosilicate
  glass to model absolute transmission.
- **Dynamic buffer-gas model** — `BufferGas`, `BufferPressure`, `BufferTemp_fill`
  parameters; K D1/Ne validated against al-Qarni et al. (2025).
- **Two-cell cascade** — `calculate_cascade()` with intermediate-polariser (MOF)
  support.
- **New outputs** — Stokes rotation azimuth `psi` and degree of linear
  polarisation `DoLP`; refractive indices n(σ±), Faraday rotation φ and α_π; and
  the Poynting-flux transmission `S0_P`.
- **Plotly Dash web GUI.**
- **Validation infrastructure** — golden-reference regression against ES3 v3.0.8,
  one-command validation via `python -m elecsus.validate`, and GitHub Actions CI.

### Changed
- Circular-basis outputs are labelled by physical σ in the GUI and documentation
  (API keys unchanged).
- **Breaking:** Stokes parameters `S1`/`S2`/`S3` are now fractional (normalised by
  `S0`, range [−1, 1]) and the `S3` sign convention is corrected. Code relying on
  the ES3 absolute Stokes values must be updated. `S0`, `Ix`, `Iy` are unchanged.

### Backwards compatibility
- `use_non_orthogonal=False` (default) reproduces ES3 propagation exactly:
  `S0`/`Ix`/`Iy` agree with ES3 v3.0.8 to ≤ 5e-9 (oblique cases ≤ 1e-5, where the
  analytic solver replaces the ES3 SymPy solve).
- The static `GammaBuf`/`shift` interface is unchanged; `BufferGas` is opt-in.
- The `scope2elecsus` / Pizzey (2022) fitting workflow runs unchanged.

### Cite
- M. A. Zentile et al., Comput. Phys. Commun. **189**, 162–174 (2015)
- J. Keaveney et al., Comput. Phys. Commun. **224**, 311–324 (2018)
