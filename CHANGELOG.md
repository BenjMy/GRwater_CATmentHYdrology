# Changelog

All notable changes to this repository are documented here. This project
is a **living repo** tied to ongoing GRwater work, so releases mark
snapshots worth citing (e.g. from a TFG report) rather than a stable
public API.

Format loosely follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Once the repo is connected to Zenodo, each tagged release below should
also get a DOI (§1 of the README).

---

## [Unreleased]


---

## [v0.1.0] — TFG Xela — 2026-09-06

First snapshot of the pipeline, corresponding to the state of the code
used for Xela Carracedo's 2026 TFG report.

### Added
- Full pyCATHY pipeline (`Agramon_withLAI_withETp.py`) driving CATHY with
  satellite ETp + rain forcing and monthly LAI-driven vegetation updates,
  for any of the 9 candidate Agramon DEM sub-plots.
- Monthly LAI → 5-class vegetation/root-depth reclassification
  (`lai_to_veg_map`), with salt-and-pepper patch smoothing
  (`smooth_veg_map`).
- Configurable boundary conditions: lateral outlet flux and optional
  free-bottom-drainage Neumann flux (`configure_boundary_conditions`).
- Optional pre-2016 spin-up by cycling one forcing year
  (`run_spinup`).
- Self-registering run log: every run appends its parameters to
  `simulation_log_withETp_withLAI.csv` (`sim_index`) and writes a
  per-run narrative `simulation_log.md` (`MarkdownLog`).
- Batch scripts `run_scenarii_withETp_withLAI.sh` and
  `plot_scenarii_withETp_withLAI.sh` for crossing dem-plot × PMIN ×
  outlet/bottom BC × spin-up and for post-processing scenario runs.

### Fixed
- ETa artefact floor (`apply_eta_artefact_floor`): patches the
  net-forcing ETa artefact so actual ET is floored at ETp when net rain
  ≥ ETp for a cell, without altering the raw `ACT. ETRA` column.


### Known limitations (see README §4.3)
- `Agramon_utils` and `geoutils.uCATHY` are external dependencies not
  included in this repository.

---

<!--
Template for the next entry:

## [vX.Y.Z] — <short name, e.g. TFG defence / plot review> — YYYY-MM-DD

### Added
### Changed
### Fixed
### Removed
-->
