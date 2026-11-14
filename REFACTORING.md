# REFACTORING.md — code map, structure and refactoring tracker

Companion to `TO-DO.md` (science, experimental design, task lists, open questions). This file holds the **code map**, the CLI reference, what has been changed, what remains, and how each change was checked. (`README_paper.md` has been merged into these two files.)

Legend: ✅ done · 🟡 partly done · ⬜ not started · ❓ decision needed · 🐞 known bug

---

## 0. Branch and rules

**Working branch:** `refactor-agramon-modules` of `https://github.com/BenjMy/GRwater_CATmentHYdrology` (not `main`, which still has the 3,901-line original). It holds the WIP baseline commit and the Phase 1–2 refactor commit.

**Rule for Claude, every conversation:** at the start of every conversation about this repository, before reading or editing code, clone the repository and check out the working branch. Do not rely on `main`, on files from earlier chats, or on memory of the code.

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 --branch refactor-agramon-modules \
    https://github.com/BenjMy/GRwater_CATmentHYdrology
git log --oneline -3     # confirm the expected tip before working
```

If the tip is not the commit the user expects, ask before continuing. The session clone is read-only: changes are delivered as patches (`git am`), one commit per step, authored as BenjMy. `Agramon_utils.py` must be committed on the branch for Claude to import and smoke-test the main script.

1. **Behaviour-preserving first.** Phases 1 and 2 must not change any simulation result. Behaviour changes (Ks wiring, withIrr rain forcing, ERA5 removal) stay in Phase 3, one at a time.
2. **Move code verbatim**, then clean. A moved function must be AST-identical to the original unless listed below.
3. **Acceptance test before adopting** (§6): rerun one existing scenario and compare the output NetCDF files with the pre-refactoring run.

---

## 1. Code map

### 1.1 Repository contents
| File | Role |
|---|---|
| `Agramon_withLAI_withETp.py` | **Main simulation script.** Full pyCATHY pipeline: monthly LAI loop, ETp/rain forcing, BCs, ICs, spin-up. |
| `run_scenarii_withETp_withLAI_ALL.sh` | **Batch runner.** Crosses dem-plot × PMIN × outlet BC × bottom BC × spin-up × IC, forwards the other flags, skips completed scenarios, post-processes. |
| `to_run_paper.md` | Scratch file of run recipes, notes, known bugs. |
| `Agramon_utils.py` | Helpers: shapefile/treatment plots, DEM loading, EM sensors and stats, ERT profiles, ERA5 point series, notebook plot helpers. |
| `results_plotter.py` | Shared plotting module (hydro, ET, LAI, recharge, hydrograph, node maps, fire-window and weekly variants, Fig. 8-style stack). |
| `agramon_config.py`, `agramon_forcing.py`, `agramon_lai.py`, `agramon_checks.py`, `agramon_runlog.py`, `agramon_outputs.py` | **New modules** split from the main script (see §3.1). Must sit in the same folder as the main script (`lib/`). |
| `analysis_ET_compare.py` | Compares EO energy-balance ETa (TSEB) vs CATHY water-balance ETa, plus comparison against the treatments CSV for `etp_ts` runs. |
| `ks_fire_utils.py` | Post-fire Ks modes (`time`, `lai`, `burn_severity`) and presets; imported by the main script through `agramon_ks.py` (task 3.4, opt-in). Science in `TO-DO.md` §5.3. |

Referenced in code but **not in the upload** (check they are in the repo): `withLAI_results.py` (post-processing driver), `geoutils` (`uCATHY`), `plot_fig8_ETa_ETp_vs_psi.py`.

### 1.2 Simulation pipeline (`Agramon_withLAI_withETp.py`)
`run_pipeline` is being split into steps (§3.2). So far steps 1–3 are in `setup_inputs(args)` → `RunInputs`; step 4 in `setup_mesh_and_dem(args, inp)` → `MeshSetup`; step 5 in `build_atmbc(args, inp, ms)` → `AtmbcSetup`; combined outputs are written by `agramon_outputs.save_combined_outputs`. The rest is still inline.

1. Resolve the DEM source (GeoTIFF, `--dem-tif`), log simulation parameters (CSV + markdown).
2. Load shapefiles (13 micro-catchments, index `fid_1`, selected with `--dem-plot 1..13`).
3. Load forcing and LAI. Grids go to EPSG:25830; LAI/ETp/rain are natively EPSG:32630. Monthly LAI comes from EO biophysical studies. Forcing switch `--forcing {era5 | etp | etp_ts}`:
   - `era5`: ERA5 point series (deprecated).
   - `etp`: satellite ETp raster + satellite rain raster (default).
   - `etp_ts`: ETp from the Agramon ET CSV, one daily value uniform over the catchment (`--etp-ts-column ET_0_wet|ET_0_dry`, `--etp-ts-feature-id`) + satellite rain.
4. Clip, mask, optionally smooth the DEM (`--dem-smooth-sigma`), build the mesh, run the pyCATHY preprocessor.
5. Map ETp and rain to surface nodes, build ATMBC arrays (net flux = rain − ETp; alternative in `TO-DO.md` §5.2).
6. Configure soil, ICs, BCs; write the `etp_output.nc` artifact.
7. Optional spin-up: cycle one year of forcing N times.
8. **Main loop, one CATHY run per calendar month:** monthly LAI → 5 vegetation classes (bare/sparse/moderate/dense/very dense; `LAI_THRESHOLDS = [1.0, 1.5, 2.0, 2.5]`); class → root depth (`ZROOT = --zroot × factor`); update soil/veg table (`update_soil_for_active_zones`) with SCF + Feddes; run processor; apply ETa artefact floor and monthly sanity checks.
9. Aggregate and write NetCDF outputs (`psi_output.nc`, `sw_output.nc`, `et_output.nc`, `etp_output.nc`, …) into `scenario_<n>_<tag>/`.

Scenario folders are self-describing, e.g. `scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c_pmin-m50_ic-wt[_forcing-etpts-wet-fid1]`.

### 1.3 CLI reference
| Flag | Meaning |
|---|---|
| `--pmin` | ψmin threshold (m) for the atmospheric ↔ soil-limited BC switch |
| `--scf [0–1]` | 0 = PMIN, 1 = Feddes (default 1). Only 0 or 1 are meaningful. |
| `--zroot`, `--zroot-factors B S M D V` | Base root depth and 5 per-class multipliers |
| `--pcana --pcref --pcwlt --pz --omgc`, `--feddes-crop` | Feddes parameters (`FP_map`) |
| `--permx`, `--permy-ratio`, `--permz-ratio` | Saturated Ks override |
| `--outlet-flux`, `--bottom-flux` | Lateral outlet / bottom Neumann flux (0 closed, negative draining) |
| `--ic-mode wt\|pressure`, `--ic-wtposition`, `--ic-pressure-head` | Initial condition |
| `--no-spinup`, `--spinup-cycles`, `--spinup-year` | Spin-up |
| `--forcing`, `--etp-ts-column` | ETp source / wet-dry bound |
| `--eta-net-neg-mode {etp,gross,off}` | Handling of the "net = rain − ETp < 0" ETa artefact |
| `--debug` | Extra stdout diagnostics (LAI NaN checks, per-catchment LAI coverage). Off by default; not forwarded by the batch runner. |

### 1.4 Quick start
```bash
# single run
python Agramon_withLAI_withETp.py --dem-plot 3 --pmin -50 --scf 1 \
    --zroot-factors 0.3 0.6 2.0 3.0 5.0 --permx 1e-4 \
    --ic-mode pressure --ic-pressure-head -100 --outlet-flux 0 --bottom-flux 0

# ETp sensitivity (wet vs dry bound)
./run_scenarii_withETp_withLAI_ALL.sh --forcing etp_ts --etp-ts-column ET_0_wet --dem-plot 3 ...
./run_scenarii_withETp_withLAI_ALL.sh --forcing etp_ts --etp-ts-column ET_0_dry --dem-plot 3 ...

# ETa evaluation
python analysis_ET_compare.py --plots all --scenario <N>
```
`./run_scenarii_withETp_withLAI_ALL.sh --help` lists all flags; `--force` reruns completed scenarios.

---

## 2. Phase 1 — cleanup (no behaviour change)

| # | Change | Status |
|---|---|---|
| 1.1 | Removed duplicated definitions in `Agramon_withLAI_withETp.py` (`zroot_factors_from_args`, `add_datetime_col`, `check_monthly_sanity`, `N_VEG_CLASSES`, four sanity thresholds). Python was already using the *second* copy, so that one was kept (`add_datetime_col` differed; effective one kept). | ✅ |
| 1.2 | Rewrote the header docstring (wrong script name, outdated usage, forcing modes) | ✅ |
| 1.3 | Fixed `lai_to_veg_map` docstring (old 0–0.6 thresholds → `[1.0, 1.5, 2.0, 2.5]`) | ✅ |
| 1.4 | Removed commented-out `ZROOT_FACTORS` block and `sys.path` stub | ✅ |
| 1.5 | Removed unused imports and needless `f` prefixes | ✅ |
| 1.6 | Unused local `fig` in `setup_dem` (kept: needed for 3D axes) | ✅ left as is |

---

## 3. Phase 2 — structure

### 3.1 Module extraction ✅
Main script went from **3,901 to 2,774 lines** (the +164 since the module split is the new dataclasses and wrapper code; `run_pipeline` itself is down to 708 lines from 1,155). Moved verbatim into sibling modules:

| Module | Contents |
|---|---|
| `agramon_config.py` | Machine paths, CRS codes, `LAI_THRESHOLDS`, `N_VEG_CLASSES`, `ZROOT_FACTORS`, `SCF_DEFAULT`, `scf_type`, `scf_label_for`, sanity thresholds, `ETP_TS_COLUMN_DEFAULT` |
| `agramon_forcing.py` | `load_monthly_lai`, `load_era5` (legacy), `load_spatial_etp`, `load_rain`, `load_agramon_et`, `load_treatments_et`, `load_all`, `build_etp_ds_from_timeseries` + helpers |
| `agramon_lai.py` | `lai_to_veg_map`, `smooth_veg_map`, `update_soil_for_active_zones`, `iter_monthly_lai` |
| `agramon_checks.py` | `apply_eta_artefact_floor`, `add_datetime_col`, `check_monthly_sanity` |
| `agramon_runlog.py` | `MarkdownLog`, `log_simulation`, `build_scenario_dirname` |
| `agramon_outputs.py` | `save_combined_outputs` (was a closure; its ten captured names are now arguments) |
| `agramon_spinup.py` | `pick_default_spinup_year`, `run_spinup` (moved unchanged; 2.5) |
| `agramon_diagnostics.py` | `save_forcing_diagnostics` (moved unchanged; 2.6) |
| `agramon_mesh.py` | `load_dem_from_tif`, `initialise_project`, `smooth_dem`, `setup_dem`, `configure_boundary_conditions` (moved unchanged; 2.3) |
| `agramon_ks.py` | per-month post-fire Ks, run key (3.4) |
| `agramon_config.py` (added) / `agramon_forcing.py` (added) | `zroot_factors_from_args`, `feddes_params_from_args` (lazy pyCATHY import) / `subset_atmbc_for_month` |

The main script re-imports these names, so `from Agramon_withLAI_withETp import load_dem_from_tif` (used by `withLAI_results.py`) keeps working.

### 3.2 Still to do
| # | Task | Status |
|---|---|---|
| 2.2 | **Split `run_pipeline`** into `setup_inputs`, `setup_mesh_and_dem`, `build_atmbc`, `configure_model`, `run_monthly_loop`, `save_outputs`, passing state in a small dataclass instead of ~160 locals | 🟡 4 of 6 extracted (`save_outputs` closure, `setup_inputs`, `setup_mesh_and_dem`, `build_atmbc`) |
| 2.3 | Move DEM/mesh/BC code (`load_dem_from_tif`, `smooth_dem`, `setup_dem`, `configure_boundary_conditions`, `initialise_project`) to `agramon_mesh.py`; keep a re-export for `withLAI_results.py` | ✅ moved unchanged; re-exported (also used by the synthetic study). `_check_crs_alignment` / `load_shapefiles` stay in the main script. Checked: pyflakes, import under stubs, `smooth_dem` test; not run with CATHY. Base for the mesh-`.vtk` properties task (`TO-DO.md` 5.4) |
| 2.4 | Replace machine auto-detection in `agramon_config.py` (username check, hard-coded Nextcloud paths) with env vars or a config file | ⬜ |
| 2.5 | Move `run_spinup` / `pick_default_spinup_year` to `agramon_spinup.py` | ✅ moved unchanged with the helpers it needs; the main script re-exports the names. Checked: pyflakes, import under stubs, smoke test; not run with CATHY |
| 2.6 | Move `save_forcing_diagnostics` (350 lines of plotting) to a diagnostics module | ✅ `agramon_diagnostics.py`, re-exported; same checks as 2.5 |
| 2.7 | Scenario config files (TOML/YAML) instead of ~40 forwarded CLI flags in the 900-line shell runner; `to_run_paper.md` becomes a list of config files | 🟡 `lib/paper_scenarios.toml` + `lib/run_paper_scenarios.py` wrap the existing shell runner (flags unchanged); the shell runner itself and the Python CLI are not replaced yet |
| 2.8 | Single source for defaults: CLI default period (`--start-year 2022 --end-year 2022`) is empty while the batch runner uses 2016–2026; shell and Python defaults differ; default `ZROOT_FACTORS` (`0.05 0.25 1.0 1.5 2.0`) differ from paper runs (`0.3 0.6 2.0 3.0 5.0`) | ❓ choose the defaults |
| 2.9 | Typed containers (`Forcing`, `BCConfig`, `KsFireConfig`) instead of long tuples/kwargs; document shapes and units of `ETp_nodes`/`rain_nodes` | ⬜ |
| 2.10 | Replace `print` banners by `logging` | ⬜ |
| 2.11 | Split `results_plotter.py` (4,400 lines) by figure family (hydro, ET, LAI, recharge, GIFs) | ⬜ |
| 2.12 | Unit tests for pure functions (LAI → veg class, Ks recovery, ETa artefact floor, scenario dirname, sanity checks) + regression test for the `--dem-plot 2` empty-mesh crash | ⬜ |
| 2.13 | Global mutable state: `FORCING_MODE`, `USE_ETP_FORCING`, `USE_ETP_TS_FORCING` are module globals set by `apply_forcing_mode`; replace by a field of the run configuration | ✅ frozen `ForcingMode` (`agramon_config.py`: `name`, `label`, `use_etp`, `use_etp_ts`), built in `setup_inputs`, carried in `RunInputs.forcing`; `apply_forcing_mode` and the globals are gone; the `era5` DeprecationWarning is kept. Checked: pyflakes, import + CLI parse under stubs, unit test; not run with CATHY |
| 2.14 | Hard-coded machine paths (`Nextcloud/...`) → config/env (same as 2.4) | ⬜ |

### 3.3 `run_pipeline` split: progress map (static analysis)
Boundaries follow the `[n/9]` banners and numbered comments.

| Step | Covers | Status |
|---|---|---|
| `setup_inputs` | config/banner, `[1/9]`–`[3/9]`: shapefiles, forcing + LAI loading, run logging. Returns frozen `RunInputs` (11 fields). Stays in main script: reads `FORCING_MODE`/`USE_ETP_*` globals (2.13) | ✅ |
| `setup_mesh_and_dem` | `[4/9]`–`[7/9]`: DEM load/mask/smooth, CRS checks, project init, DEM → mesh, preprocessor, grid dims. `(args, inp)` → frozen `MeshSetup` (`simu`, `raster_DEM`, `raster_DEM_masked`, `xllcorner`, `yllcorner`, `res_x`, `res_y`, `grid3d`, `N`, `M`). Stays until `agramon_mesh.py` exists (2.3) | ✅ |
| `build_atmbc` | `[8/9]`: `ds_mesh`, ETp/rain on surface nodes, net flux, forcing diagnostics (~140 lines). `(args, inp, ms)` → frozen `AtmbcSetup` (`ds_mesh`, `n_surf`, `ETp_nodes`, `rain_nodes`, `net_flux`, `t_atmbc`, `valid_times`). Covers ETp branch and legacy ERA5 branch (there `ETp_nodes`, `rain_nodes`, `n_surf` are `None`). Stays: reads `USE_ETP_*` (2.13), calls `save_forcing_diagnostics` (2.6) | ✅ |
| `configure_model` | BCs and IC, LAI interpolation to the DEM grid, spin-up. Consumes `MeshSetup` + `AtmbcSetup`; `n_surf` and `ds_mesh` are *rebound inside the monthly loop* (grid refresh), so they must stay plain locals for `run_monthly_loop` | ⬜ |
| `run_monthly_loop` | one CATHY run per month (~455 lines; itself needs a `run_one_month` split) | ⬜ |
| `save_outputs` | `save_combined_outputs` now in `agramon_outputs.py`; `etp_output.nc` write and final prints still inline | 🟡 closure extracted |

Findings:
- `save_combined_outputs` was a closure over ten names, each bound exactly once in `run_pipeline` (lists/dicts mutated in place), so passing them as arguments is behaviour-preserving.
- The LAI-interpolation block had leftover debug prints (raw LAI patch shape, NaN counts, per-`fid_1` loop). ✅ Kept behind the new `--debug` flag (off by default); a commented-out variant was removed.
- The monthly loop assigns/reads ~60 locals; shared state is mainly `simu`, `ds_mesh`, `ETp_nodes`/`rain_nodes`, `t_atmbc`/`net_flux`, `valid_times`, the `*_records` lists and `out_dir` → candidates for a `RunState` dataclass.

**Per-step checks (all static / stubbed; a real simulation was not run):**
- `setup_inputs`: moved block (38 statements) AST-identical; rest of `run_pipeline` (79 statements) identical; only additions are `setup_inputs`/`RunInputs` and 12 alias lines. Run with loaders stubbed in the three forcing modes: `RunInputs` filled as expected. `df_et` stays local (reassigned before being read).
- `build_atmbc`: moved block (6 top-level statements) AST-identical; rest (52 statements) identical. The only added statements are `ETp_nodes = None`, `rain_nodes = None`, `n_surf = None` before the `if` (in ERA5 mode the original never bound them and every later read is behind `USE_ETP_FORCING`). Old and new `run_pipeline` run side by side with fakes (DEM loader, project, `setup_dem`, `mt.build_mesh_dataset`, `mt.map_grid_to_mesh`, `save_forcing_diagnostics`, ERA5 extractor) in all three forcing modes, with partly overlapping ETp/rain time axes (37 common steps) and node-dependent values: stdout, call log and hashed contents of `ETp_nodes`, `rain_nodes`, `net_flux`, `t_atmbc`, `valid_times`, `n_surf`, `ds_mesh` identical. The test fails when arrays are swapped on purpose.
- `setup_mesh_and_dem`: moved block (21 statements) AST-identical; rest (58 statements) identical. Side-by-side run with recording fakes, stopping at the start of `[8/9]`: stdout, call log and handed-over locals identical for the default case and `--dem-smooth-sigma 1.5`. Fails when `N`/`M` are swapped on purpose. `pyflakes` clean apart from `fig` in `setup_dem`.
- `save_outputs` (closure): original closure and `save_combined_outputs` run on the same synthetic data gave identical files (`xarray.testing.assert_identical` for NetCDF, byte-equal CSV, equal pickle), including the empty-state case.

**🐞 Finding (not changed):** the legacy ERA5 branch computes `t_atmbc = [float(t - _t0) / 1e9 for t in valid_times]`, valid only for `datetime64[ns]`. With pandas 3 / current xarray (`datetime64[us]`) it raises `TypeError: float() argument must be a string or a real number, not 'datetime.timedelta'` (reproduced with the original code). The ETp branch uses `/ np.timedelta64(1, "s")` and is unaffected. Not fixed (Phases 1–2 are behaviour-preserving); one more reason to remove ERA5 (3.3).

---

## 4. Phase 3 — variants (behaviour changes, one at a time)

### Should the script be split per CATHY version and per input?
**Recommendation: split by responsibility, not by copy.** Do not keep `Agramon_withLAI_withETp.py`, `..._withIrr.py`, `..._era5.py` as separate full scripts: ~80 % of the pipeline is identical (DEM, mesh, LAI loop, spin-up, outputs, logging), so copies would drift and bug fixes would be made several times (the Phase 1 duplicates are a small example). Isolate the *differences* behind two interfaces chosen by CLI/config flag:

1. **Forcing provider** (`--forcing`): returns ETp and rain on the mesh. `SatelliteForcing` (`etp`), `CsvEtpForcing` (`etp_ts`, wet/dry), legacy `Era5Forcing`. Each in its own module; the pipeline calls `provider.build(...)`.
2. **CATHY backend** (`--cathy-version`): wraps what differs between CATHY/pyCATHY builds (sketch; the `withIrr` side has not been read yet):

| Concern | `scf_lai` (new LAI version) | `with_irr` |
|---|---|---|
| Atmospheric file | ATMBC with net flux (rain − ETp), plus ETa artefact patch | ATMBC with **ETp only**; rain goes into the irrigation map file |
| Soil / veg table | LAI classes, `ZROOT`, Feddes, SCF | to be checked |
| Output reading / ETa | patched ETa | raw ETa (no patch) |

| # | Task | Status |
|---|---|---|
| 3.1 | Define `ForcingProvider` and move the three modes behind it | ⬜ |
| 3.2 | Define `CathyBackend`; implement `scf_lai` (current behaviour) first, then `with_irr` (rain via irrigation file, no ETa floor), then A/B test (`TO-DO.md` §5.2). Read the `CATHY_src_WithIrrigation` / wrapper code before fixing the interface | ⬜ |
| 3.3 | **ERA5:** deprecated (`DeprecationWarning` on `--forcing era5`; `etp` stays default). Removal list: `load_era5`, the `era5` entry of `FORCING_LABELS`/`--forcing`, `ds_era5` loading, the "ERA5 point time-series" branch in `build_atmbc` (and the `datetime64[ns]` assumption), ERA5 labels in plotting code. ERA5 could stay as a small test fixture | 🟡 deprecated, ❓ remove when you confirm |
| 3.4 | Wire `ks_fire_utils` into the monthly loop (`compute_permx_for_month` → `update_soil_for_active_zones(permx=...)`; CLI helpers `add_ks_fire_args`, `ks_fire_config_from_args` exist) and define the Ks scenario list (preset × mode). Decide surface-layer-only vs all zones (`TO-DO.md` §5.3) | 🟡 wired (opt-in) via `lib/agramon_ks.py`: `--ks-fire-mode off|time|lai|burn_severity` (+ `--ks-fire-*`) in the main script, per-month PERMX in the monthly soil update, unburned Ks before the fire date, Ks tag in the folder name, config validated at startup. All zones (like `--permx`); surface-only not implemented. The shell runner forwards `--ks-fire-*` (same for every run in the batch) and its skip-if-done check compares the new `ks_fire` column of the CSV log (rows without it = static Ks); checked on a mock log, shell runner not run end to end. Ks scenario list (preset × mode) still to define. Not run with CATHY |

---

## 5. Behaviour changes made so far (only these)
1. `apply_forcing_mode("era5")` emits a `DeprecationWarning`. The run is unchanged.
2. New `--debug` flag (default off): the LAI-interpolation diagnostics print only with it. No simulation input or output is affected (no later code reads those variables).

Everything else in Phases 1–2 is code moved or removed without changing what runs.

---

## 6. Verification

**Checked (static, no real run):**
- All files compile; `pyflakes` reports no undefined names.
- Every top-level function/class of the original is present exactly once after the split with identical AST (the *effective* definition where duplicates existed). Only three differ, as intended: `apply_forcing_mode` (deprecation warning), `lai_to_veg_map` (docstring), `configure_boundary_conditions` (two `f` prefixes removed).
- With `pyCATHY`, `geoutils`, `Agramon_utils`, `xarray`, `geopandas`, `rioxarray` stubbed, old and new scripts import, parse CLI arguments to identical results (default and a long example), and expose identical constants and paths.

**Not checked:** a real simulation (the pyCATHY/CATHY executable and Agramon data are not available in the session).

**Acceptance test to run on your machine** before replacing the old files:
1. Pick one short existing scenario (e.g. one year, `--dem-plot 3`, fixed `--pmin`/`--permx`).
2. Run it with the old and the refactored script into different `--path2prj` folders, same `--force`.
3. Compare `psi_output.nc`, `sw_output.nc`, `et_output.nc`, `etp_output.nc` and the hydrograph CSV (`xarray.testing.assert_identical` / `assert_allclose`). They should match.

---

## 7. Known bugs
- 🐞 **`--dem-plot 2` crashes** in `run_pipeline` at `ds_mesh.x.min()`: `ValueError: zero-size array to reduction operation minimum` (mesh/DEM empty for that plot: clip or mask returns no valid cells). Check footprint, CRS, NaN masking, `--dem-smooth-sigma`. Needs a regression test (2.12).
- 🐞 pyCATHY `map_prop2mesh` uses `~isinstance(...)` instead of `not isinstance(...)`, so per-node arrays are broadcast as scalars; a related stale `mesh_pv_attributes` cache can cause node-count mismatches.
- 🟡 `n_surf = int(simu.grid3d["nnod"])` used the stale attribute, giving "SURFACE NODE out of bounds" warnings. The script now reads `grid3d["nnod"]` from the re-read grid, so it looks fixed; confirm with a run (no warning) before closing.
- 🐞 ERA5 branch `datetime64[ns]` assumption (§3.3 finding).
- 🐞 Batch recipe: `--permx "${KS_VALUES}"` runs only the first Ks (see `TO-DO.md` §2.3).

---

## 8. Applying the patch (`agramon_refactor_phase1-2.patch`)

Made against the working copy of `lib/Agramon_withLAI_withETp.py` shared for this work. It touches `lib/` (main script changed, five `agramon_*.py` modules added) and adds this `REFACTORING.md` at the repo root. `refactored/` holds the new main script (same name) and the `agramon_*.py` modules; copy them into the folder with `Agramon_utils.py` (`lib/`).

```bash
cd /home/ben/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/

# 0. clean tree, new branch (easy to undo)
git status
git checkout -b refactor-agramon-modules

# 1. if REFACTORING.md was copied by hand, remove it (the patch creates it)
rm -f REFACTORING.md

# 2. dry run
git apply --check /path/to/agramon_refactor_phase1-2.patch

# 3. apply
git am /path/to/agramon_refactor_phase1-2.patch      # or: git apply ...

# 4. quick checks (sibling modules next to Agramon_utils.py)
cd lib
python -m py_compile Agramon_withLAI_withETp.py agramon_*.py
python Agramon_withLAI_withETp.py --help
```

**Undo:**
```bash
git reset --hard HEAD~1                      # if applied with 'git am'
git apply -R /path/to/agramon_refactor_phase1-2.patch   # if applied with 'git apply'
```

**If `--check` fails:** the main script in your tree is not the version the patch was built against (the GitHub copy differs by three lines: `MODULE_PATH` vs `EO_PATH` and the shapefile path). Do not force it; send the current `lib/Agramon_withLAI_withETp.py` and the patch will be rebuilt.

**Flat folder without `lib/`:** from the folder containing `Agramon_withLAI_withETp.py`:
```bash
patch -p2 --dry-run < /path/to/agramon_refactor_phase1-2.patch
patch -p2 < /path/to/agramon_refactor_phase1-2.patch
```
`patch` reports one skipped hunk (`REFACTORING.md`); copy it to the repo root by hand. (`git apply -p2` does not work for this layout.)

After applying, run the acceptance test (§6) before relying on the refactored code.
