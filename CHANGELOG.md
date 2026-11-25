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

### Added
- Run registry on top of the CSV run log (`agramon_registry.py`): status, seconds, n_months,
  relative `sim_dir`, `param_hash`, `git_commit`, `error_tail`; `backfill_registry.py` for old logs.
- `run_catalog.py` (summary table of runs) and `plot_compare_runs.py` (comparison figures from the log).
- Partial-run detection: runs that stopped early are never reused.
- `manifest.json` of a hypothesis study keeps the stderr tail (`error_tail`) of failed variants.

### Changed
- `--study` of `run_catalog.py`, `plot_compare_runs.py`, `clean_log.py` and the new `post_process.py --study`
  accepts the hypothesis name (`hyp4`, `hyp4_SW`) as well as a folder; resolved from the config by
  `hypotheses_common.resolve_study`.
- Hypothesis figures: variants of one group differ in shade, dash pattern and width (they were all one colour);
  the legend of `response` / `annual` / `eta_vs_EB` is a single wrapped block under the plots; `balance` wraps
  its y labels.
- `manifest.json` stores `sim_dir` relative to the study `prj/` (absolute paths still read).
- `scenario_runner.py` skip check ignores failed, running and partial runs.
- `clean_log.py` and `post_process.py` read relative manifest paths and skip failed/running runs.

### Added
- `lib/data_paths.toml`: every input path in one file (`AGRAMON_DATA_CONFIG` to use another site).
- Two-layer Ks scenario (`--permx-bottom`, `--ks-layer-depth`).
- Hypothesis framework: `hypotheses.toml`, `hyp1_layers_Ks.py` .. `hyp4_SW.py`, `run_all_hypotheses.py`
  (dry-run, validation vs EO energy-balance ETa, plots), `test/test_hypotheses.py`.
- `lib/hypothesis_tables.py`: annual / pre- / post-fire water-balance totals (rain, ETa, recharge, runoff,
  closure), per-variant means, and the recovery hook (`recovery_metrics`, active once a control run is
  configured); `study_report.md` per study; mean ETa / recharge / saturation columns in `summary_metrics.csv`.
- `lib/clean_log.py`: removes the run-log rows (and stale hypothesis-manifest entries) whose scenario
  folder no longer exists; `--dry-run`, backup of the CSV, `--study` / `--all-studies`, `--incomplete`.
- `lib/scenario_runner.py`: Python batch runner (REFACTORING 2.7), same flags as the old shell runner;
  `run_paper_scenarios.py --dry-run` now shows the real driver calls and the runs it would skip.

### Changed
- `run_scenarii_withETp_withLAI_ALL.sh` is a thin wrapper around `scenario_runner.py` (deprecated).
- Skip-if-done check: `--zroot-factors` are compared as the driver logs them (`2`, not `2.0`);
  before, runs with such factors were never skipped.
- `withLAI_results.py` -> `plot_run_results.py`, `results_plotter.py` -> `plot_results_shared.py`.
- Fig 8 (ETa/ETp vs psi) moved into `plot_eta.py`; vegetation colormap fixed (5 classes).

### Removed
- Batch `.log` writing, `to_run_paper.md`, `plot_fig8_ETa_ETp_vs_psi.py`.


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
