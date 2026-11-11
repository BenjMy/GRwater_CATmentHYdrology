# REFACTORING.md — readability / structure tracker

Companion to `README.md`. It records what has been changed in the code, what still has to be done, and how each change was checked.

Legend: ✅ done · 🟡 partly done · ⬜ not started · ❓ decision needed

---

## 0. Branch and rules for this refactoring

**Working branch:** `refactor-agramon-modules` of `https://github.com/BenjMy/GRwater_CATmentHYdrology` (not `main`; `main` still has the 3,901-line original). It holds the WIP baseline commit and the Phase 1–2 refactor commit.

**Rule for Claude, every conversation:** at the start of every conversation about this repository, before reading or editing any code, clone the repository and check out the working branch. Do not rely on `main`, on files from earlier chats, or on memory of the code.

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 --branch refactor-agramon-modules \
    https://github.com/BenjMy/GRwater_CATmentHYdrology
git log --oneline -3     # confirm the expected tip before working
```

If the branch tip is not the commit the user expects, ask before continuing. The session clone is read-only: changes are delivered as patches (`git am`), one commit per step. `Agramon_utils.py` must be committed on the branch for Claude to import and smoke-test the main script.

1. **Behaviour-preserving first.** Phases 1 and 2 must not change any simulation result. Behaviour changes (Ks wiring, withIrr rain forcing, ERA5 removal) are kept in Phase 3 and are made one at a time.
2. **Move code verbatim**, then clean. A moved function must be AST-identical to the original unless the change is listed below.
3. **Acceptance test before adopting** (see §5): rerun one existing scenario with the new code and compare the output NetCDF files with the pre-refactoring run.

---

## 1. Phase 1 — cleanup (no behaviour change)

| # | Change | Status |
|---|---|---|
| 1.1 | Removed duplicated definitions in `Agramon_withLAI_withETp.py`: `zroot_factors_from_args`, `add_datetime_col`, `check_monthly_sanity`, `N_VEG_CLASSES`, and the four sanity thresholds (`SW_TOL`, `PSI_SANITY_MIN`, `ETA_TOL`, `WB_RESIDUAL_FRAC`). Python was already using the *second* copy of each function, so the second copy was kept. (`add_datetime_col` differed between the two copies; the effective one was kept.) | ✅ |
| 1.2 | Rewrote the file header docstring (wrong script name, outdated usage lines, forcing modes) | ✅ |
| 1.3 | Fixed the `lai_to_veg_map` docstring, which listed the old 0–0.6 LAI thresholds; it now matches `LAI_THRESHOLDS = [1.0, 1.5, 2.0, 2.5]` | ✅ |
| 1.4 | Removed a commented-out alternative `ZROOT_FACTORS` block and a commented `sys.path` stub | ✅ |
| 1.5 | Removed unused imports (`os`, `sys`, `mdates`, and the ones that became unused after the split), and `f` prefixes on strings without placeholders | ✅ |
| 1.6 | Unused local `fig` in `setup_dem` (kept: `plt.figure` is needed for the 3D axes) | ✅ left as is |

---

## 2. Phase 2 — structure

### 2.1 Module extraction ✅
`Agramon_withLAI_withETp.py` went from **3,901 to 2,774 lines** (2,683 after the first module split; 2,610 after the outputs extraction; the +164 since then is the new `RunInputs`/`MeshSetup`/`AtmbcSetup` dataclasses and the `setup_inputs`/`setup_mesh_and_dem`/`build_atmbc` wrapper code and docstrings; `run_pipeline` itself is down to 708 lines). Moved verbatim into sibling modules (same folder, imported by name like `Agramon_utils`):

| New module | Contents |
|---|---|
| `agramon_config.py` | Machine paths, CRS codes, `LAI_THRESHOLDS`, `N_VEG_CLASSES`, `ZROOT_FACTORS`, `SCF_DEFAULT`, `scf_type`, `scf_label_for`, sanity thresholds, `ETP_TS_COLUMN_DEFAULT` |
| `agramon_forcing.py` | `load_monthly_lai`, `load_era5` (legacy), `load_spatial_etp`, `load_rain`, `load_agramon_et`, `load_treatments_et`, `load_all`, `build_etp_ds_from_timeseries` and their helpers |
| `agramon_lai.py` | `lai_to_veg_map`, `smooth_veg_map`, `update_soil_for_active_zones`, `iter_monthly_lai` |
| `agramon_checks.py` | `apply_eta_artefact_floor`, `add_datetime_col`, `check_monthly_sanity` |
| `agramon_runlog.py` | `MarkdownLog`, `log_simulation`, `build_scenario_dirname` |
| `agramon_outputs.py` | `save_combined_outputs` (was a closure inside `run_pipeline`; its ten captured names are now arguments) |

The main script imports these names, so `from Agramon_withLAI_withETp import load_dem_from_tif` (used by `withLAI_results.py`) and the other existing entry points keep working.

### 2.2 Still to do
| # | Task | Status |
|---|---|---|
| 2.2 | **Split `run_pipeline`** (was 1,155 lines, now 708) into steps: `setup_inputs`, `setup_mesh_and_dem`, `build_atmbc`, `configure_model`, `run_monthly_loop`, `save_outputs`; pass state in a small dataclass instead of ~160 locals. See "2.2 progress" below | 🟡 map done, 4 of 6 steps extracted (`save_outputs` closure, `setup_inputs`, `setup_mesh_and_dem`, `build_atmbc`) |
| 2.3 | Move DEM / mesh / BC code (`load_dem_from_tif`, `smooth_dem`, `setup_dem`, `configure_boundary_conditions`, `initialise_project`) to `agramon_mesh.py`. Keep a re-export in the main script because `withLAI_results.py` imports `load_dem_from_tif` from it | ⬜ |
| 2.4 | Replace machine auto-detection in `agramon_config.py` (home-directory username check, hard-coded Nextcloud paths) with environment variables or a config file | ⬜ |
| 2.5 | Move `run_spinup` / `pick_default_spinup_year` to `agramon_spinup.py` | ⬜ |
| 2.6 | Move `save_forcing_diagnostics` (350 lines of plotting) to a diagnostics module | ⬜ |
| 2.7 | Scenario config files (TOML/YAML) instead of ~40 forwarded CLI flags in the 900-line shell runner. `to_run_paper.md` would become a list of config files | ⬜ |
| 2.8 | Single source for defaults (CLI default period is empty, shell and Python defaults differ, default `ZROOT_FACTORS` differ from those used in the paper runs) | ❓ choose the defaults |
| 2.9 | Typed containers (`Forcing`, `BCConfig`, `KsFireConfig`) instead of long tuples/keyword lists; document shapes and units of `ETp_nodes`/`rain_nodes` | ⬜ |
| 2.10 | Replace `print` banners by `logging` | ⬜ |
| 2.11 | Split `results_plotter.py` (4,400 lines) by figure family (hydro, ET, LAI, recharge, GIFs) | ⬜ |
| 2.12 | Unit tests for pure functions (LAI → veg class, Ks recovery, ETa artefact floor, scenario dirname, sanity checks) and a regression test for the `--dem-plot 2` empty-mesh crash | ⬜ |
| 2.13 | Global mutable state: `FORCING_MODE`, `USE_ETP_FORCING`, `USE_ETP_TS_FORCING` are module globals changed by `apply_forcing_mode`. Replace by a field of the run configuration | ⬜ |

#### 2.2 progress (map of `run_pipeline`, static analysis)
Boundaries follow the `[n/9]` banners and numbered comments in the code, so they survive edits.

| Step | Covers | Status |
|---|---|---|
| `setup_inputs` | config/banner, `[1/9]`–`[3/9]`: shapefiles, forcing + LAI loading, run logging. Returns a frozen `RunInputs` dataclass (11 fields). Stays in the main script for now: it reads `FORCING_MODE`/`USE_ETP_*`, module globals set by `apply_forcing_mode` (see 2.13) | ✅ |
| `setup_mesh_and_dem` | `[4/9]`–`[7/9]`: DEM load/mask/smooth, CRS checks, project init, DEM → mesh, preprocessor, grid dims. Takes `(args, inp)` and returns a frozen `MeshSetup` dataclass (10 fields: `simu`, `raster_DEM`, `raster_DEM_masked`, `xllcorner`, `yllcorner`, `res_x`, `res_y`, `grid3d`, `N`, `M`). Stays in the main script until `agramon_mesh.py` exists (2.3), because `load_dem_from_tif`/`setup_dem`/`initialise_project` still live there | ✅ |
| `build_atmbc` | `[8/9]`: `ds_mesh` build, ETp/rain mapped to surface nodes, net flux, forcing diagnostics (~140 lines). Takes `(args, inp, ms)` and returns a frozen `AtmbcSetup` dataclass (7 fields: `ds_mesh`, `n_surf`, `ETp_nodes`, `rain_nodes`, `net_flux`, `t_atmbc`, `valid_times`). Covers both the ETp branch and the legacy ERA5 branch. Stays in the main script because it reads `USE_ETP_FORCING`/`USE_ETP_TS_FORCING` (2.13) and calls `save_forcing_diagnostics` (2.6). In ERA5 mode `ETp_nodes`, `rain_nodes` and `n_surf` are `None` (see the check below) | ✅ |
| `configure_model` | BCs and IC, LAI interpolation to the DEM grid, spin-up. Next step: consumes `MeshSetup` + `AtmbcSetup`; `n_surf` and `ds_mesh` are *rebound inside the monthly loop* (grid refresh), so they must stay plain locals of `run_pipeline` for `run_monthly_loop` | ⬜ |
| `run_monthly_loop` | one CATHY run per month (~455 lines, itself needs a `run_one_month` split) | ⬜ |
| `save_outputs` | `save_combined_outputs` (now in `agramon_outputs.py`); `etp_output.nc` write and final prints still inline | 🟡 closure extracted |

Findings from the map:
- `save_combined_outputs` was a closure over ten names. Each is bound exactly once in `run_pipeline` (lists/dicts mutated in place), so passing them as arguments is behaviour-preserving.
- The LAI-interpolation block contains leftover debug prints (raw LAI patch shape, NaN counts, a per-`fid_1` loop). They only write to stdout. ✅ Decided: kept behind a new `--debug` flag (off by default); a commented-out variant of the same checks was removed. Note that `run_scenarii_withETp_withLAI_ALL.sh` does not forward `--debug`; call the Python script directly to use it.
- The monthly loop assigns/reads ~60 locals; most of the shared state is `simu`, `ds_mesh`, `ETp_nodes`/`rain_nodes`, `t_atmbc`/`net_flux`, `valid_times`, the `*_records` lists and `out_dir`. These are the candidates for the `RunState` dataclass.

Check for step `setup_inputs`: the moved block (38 statements) is AST-identical to the original, and so is the rest of `run_pipeline` (79 statements); the only additions are `setup_inputs`/`RunInputs` and 12 lines in `run_pipeline` (`inp = setup_inputs(args)` plus 11 aliases, so the rest of the function stays verbatim). `setup_inputs` was also run with the loaders stubbed in the three forcing modes (`etp`, `etp_ts`, `era5`): `RunInputs` is filled as expected (ERA5 mode has `ds_etp`/`ds_rain` = None, ETp modes have `ds_era5` = None). `df_et` stays local to `setup_inputs`: later code reassigns it before reading it. A real simulation was not run.

Check for step `build_atmbc`: the moved block (6 top-level statements, ~140 lines: from the `[8/9]` banner to the end of the `if USE_ETP_FORCING … else` ERA5 branch) is AST-identical to the original, and so is the rest of `run_pipeline` (52 statements); no other definition changed. Additions: `AtmbcSetup`, `build_atmbc` (10 aliases of `RunInputs`/`MeshSetup` at its top, same names as before, so the body is verbatim), and 8 lines in `run_pipeline` (`atm = build_atmbc(args, inp, ms)` plus 7 aliases). The `show_plots`, `ds_era5` and `etp_ts_fid` aliases in `run_pipeline` were dropped because only the moved block read them. Names that the loop overwrites before reading (`ds_ETp_mapped`, `ds_rain_mapped`, `ds_ETp_surf`, `ds_rain_surf`, and the ERA5 plot locals `fig`, `t`) stay local to `build_atmbc`.

**The only added statements in the moved code:** `ETp_nodes = None`, `rain_nodes = None`, `n_surf = None` before the `if`. In ERA5 mode the original never bound these names, and every later read is behind `USE_ETP_FORCING` (`apply_eta_artefact_floor`, `check_monthly_sanity` arguments, the `etp_output.nc` write, the n_surf refresh in the loop), so the new `None` values are never used; they exist so the dataclass can be returned. In the ETp modes they are overwritten by the branch.

Old (previous commit) and new `run_pipeline` were run side by side with the DEM loader, project, `setup_dem`, `mt.build_mesh_dataset`, `mt.map_grid_to_mesh`, `save_forcing_diagnostics` and the ERA5 extractor replaced by fakes, stopping at the call to `configure_boundary_conditions`, in all three forcing modes (`etp`, `etp_ts`, `era5` with `show_plots=True`). ETp and rain were given different, partly overlapping time axes (40 days each, 3 days apart → 37 common steps) and node-dependent values, so the time alignment, the `n_surf` slice and the mm/day → m/s conversion are exercised. Stdout, the full call log (including the arguments of `save_forcing_diagnostics` and the markdown log entries) and the hashed contents of `ETp_nodes`, `rain_nodes`, `net_flux`, `t_atmbc`, `valid_times`, `n_surf` and `ds_mesh` are identical (ERA5: the three names above are unbound in the old code and `None` in the new one, as intended). The test fails when `ETp_nodes`/`rain_nodes` or `t_atmbc`/`net_flux` are swapped on purpose. A real simulation was not run.

**Finding (not changed):** the legacy ERA5 branch computes `t_atmbc = [float(t - _t0) / 1e9 for t in valid_times]`, which only works for `datetime64[ns]` time axes. With pandas 3 / current xarray, which keep `datetime64[us]`, it raises `TypeError: float() argument must be a string or a real number, not 'datetime.timedelta'` (reproduced with the *original* code; the ETp branch uses `/ np.timedelta64(1, "s")` and is not affected). Not fixed because Phases 1–2 are behaviour-preserving; it is one more reason to remove ERA5 (3.3).

Check for step `setup_mesh_and_dem`: the moved block (21 statements, from the `[4/9]` banner to the `Grid:` log bullet) is AST-identical to the original, and so is the rest of `run_pipeline` (58 statements). The only additions are `MeshSetup`, `setup_mesh_and_dem` (6 aliases of `RunInputs` at its top, same names as before, so the body is verbatim) and 10 lines in `run_pipeline` (`ms = setup_mesh_and_dem(args, inp)` plus 8 aliases for the names read later: `simu`, `raster_DEM`, `raster_DEM_masked`, `res_x`, `res_y`, `grid3d`, `N`, `M`). The now-unused `dem_tif_path` alias in `run_pipeline` was dropped; `xllcorner`/`yllcorner` are in `MeshSetup` but not aliased because nothing reads them afterwards. `pyflakes` is clean apart from the known `fig` in `setup_dem` (1.6). Old and new `run_pipeline` were run side by side (original taken from git) with the DEM loader, project, `setup_dem`, preprocessor and `mt.build_mesh_dataset` replaced by recording fakes, stopping at the start of `[8/9]`: stdout and the full call log are identical for the default case and for `--dem-smooth-sigma 1.5` (smoothing + its two asserts), and the locals handed to the next step (`N`, `M`, `grid3d`, `raster_DEM`, `raster_DEM_masked`, `res_x`, `res_y`, `simu`) are identical. The test fails when `N`/`M` are swapped on purpose, so it does detect wiring mistakes. A real simulation was not run.

Check for step `save_outputs` (closure): the original closure (from the previous commit) and `agramon_outputs.save_combined_outputs` were run on the same synthetic data and produced identical files (`xarray.testing.assert_identical` for the NetCDF files, byte-equal CSV, equal pickle), including the empty-state case. A real simulation was not run.

---

## 3. Phase 3 — variants (behaviour changes, one at a time)

### Should the script be split per CATHY version and per input?
**Recommendation: split by responsibility, not by copy.** Do not keep `Agramon_withLAI_withETp.py`, `..._withIrr.py`, `..._era5.py` as separate full scripts: ~80 % of the pipeline is identical (DEM, mesh, LAI loop, spin-up, outputs, logging), so copies would drift and every bug fix would have to be made several times (the duplicated functions found in Phase 1 are a small example).

Instead isolate the *differences* behind two small interfaces, chosen by CLI/config flag:

1. **Forcing provider** (`--forcing`): returns ETp and rain on the mesh. Candidates: `SatelliteForcing` (`etp`), `CsvEtpForcing` (`etp_ts`, wet/dry bound), and the legacy `Era5Forcing`. Each lives in its own module/file; the pipeline only calls `provider.build(...)`.
2. **CATHY backend** (`--cathy-version`): wraps what differs between CATHY/pyCATHY builds. The sketch below is based on what is known from the code and your notes; the `withIrr` side has not been read yet.

| Concern | `scf_lai` (new LAI version) | `with_irr` |
|---|---|---|
| Atmospheric file | ATMBC with net flux (rain − ETp), plus the ETa artefact patch | ATMBC with **ETp only**; rain goes into the irrigation map file (see README §4.5) |
| Soil / veg table | LAI classes, `ZROOT`, Feddes, SCF | to be checked |
| Output reading / ETa | patched ETa | raw ETa (no patch needed) |

| # | Task | Status |
|---|---|---|
| 3.1 | Define the `ForcingProvider` interface and move the three modes behind it | ⬜ |
| 3.2 | Define the `CathyBackend` interface; implement `scf_lai` (current behaviour) first, then `with_irr` (rain via irrigation file, no ETa floor), then A/B test (README §4.5). Read the `CATHY_src_WithIrrigation` / wrapper code before fixing the interface | ⬜ |
| 3.3 | **ERA5:** deprecated now (a `DeprecationWarning` is raised when `--forcing era5` is used; `etp` stays the default). Removal list: `load_era5`, the `era5` entry of `FORCING_LABELS`/`--forcing`, the `ds_era5` loading and the "ERA5 point time-series" branch in `build_atmbc` (and the `datetime64[ns]` assumption noted in §2.2), and the ERA5 labels in the plotting code. ERA5 could be kept as a small test fixture instead | 🟡 deprecated, ❓ remove when you confirm |
| 3.4 | Wire `ks_fire_utils` into the monthly loop (`compute_permx_for_month` → `update_soil_for_active_zones(permx=...)`) and define the Ks scenario list (README §4.6) | ⬜ |

---

## 4. Behaviour changes made so far (only these)

1. `apply_forcing_mode("era5")` now emits a `DeprecationWarning`. It does not change the run.
2. New `--debug` flag (default off). The LAI-interpolation diagnostics (raw LAI patch NaN map, NaN counts inside/outside the DEM footprint, per-`fid_1` LAI coverage) are printed only with `--debug`. Without it they are skipped, so default runs print less. No simulation input or output is affected: no later code reads the variables they computed.

Everything else in Phases 1 and 2 is code moved or removed without changing what runs.

---

## 5. How the changes were checked, and what was not

**Checked (static, no real run):**
- All files compile; `pyflakes` reports no undefined names.
- Every top-level function and class of the original is present exactly once after the split, and its AST is identical to the original (the *effective* definition where duplicates existed). Only three differ, as intended: `apply_forcing_mode` (deprecation warning), `lai_to_veg_map` (docstring), `configure_boundary_conditions` (two `f` prefixes removed).
- With `pyCATHY`, `geoutils`, `Agramon_utils`, `xarray`, `geopandas` and `rioxarray` stubbed out, the old and new scripts import fine, parse the same CLI arguments to identical results (default and a long example), and expose identical constants and paths.

**Not checked:** a real simulation. The pyCATHY / CATHY executable and the Agramon data are not available here.

**Acceptance test to run on your machine** before replacing the old files:
1. Pick one short existing scenario (for example one year, `--dem-plot 3`, a fixed `--pmin`/`--permx`).
2. Run it with the old script and with the refactored one into different `--path2prj` folders, with the same `--force`.
3. Compare `psi_output.nc`, `sw_output.nc`, `et_output.nc`, `etp_output.nc` and the hydrograph CSV (for example with `xarray.testing.assert_identical`/`assert_allclose`). They should match.

---

## 6. Files

`refactored/` contains the new main script (same file name as before) and the five `agramon_*.py` modules. Copy all six into the same folder as `Agramon_utils.py` (`lib/` in the repository).

---

## 7. Applying the patch (`agramon_refactor_phase1-2.patch`)

The patch is a git patch made against the working copy of `lib/Agramon_withLAI_withETp.py` that was shared for this work. It touches `lib/` (main script changed, five `agramon_*.py` modules added) and adds this `REFACTORING.md` at the repo root.

```bash
cd /home/ben/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/

# 0. start from a clean tree, on a new branch, so it is easy to undo
git status
git checkout -b refactor-agramon-modules

# 1. if you already copied REFACTORING.md here by hand, remove it first
#    (the patch creates it, and git refuses to overwrite an existing file)
rm -f REFACTORING.md

# 2. dry run: changes nothing
git apply --check /path/to/agramon_refactor_phase1-2.patch

# 3. apply (keep the commit message with 'git am', or just apply the changes)
git am /path/to/agramon_refactor_phase1-2.patch
#   or: git apply /path/to/agramon_refactor_phase1-2.patch

# 4. quick checks (the sibling modules must sit next to Agramon_utils.py)
cd lib
python -m py_compile Agramon_withLAI_withETp.py agramon_*.py
python Agramon_withLAI_withETp.py --help
```

**Undo:**
```bash
git reset --hard HEAD~1                      # if applied with 'git am'
# or, if applied with 'git apply':
git apply -R /path/to/agramon_refactor_phase1-2.patch
```

**If `--check` fails:** the main script in your tree is not the version the patch was built against (the GitHub copy, for example, differs by three lines: `MODULE_PATH` vs `EO_PATH` and the shapefile path). Do not force it. Send the current `lib/Agramon_withLAI_withETp.py` and the patch will be rebuilt against it.

**Flat folder without `lib/` (for example a copy of the scripts elsewhere):** run it from the folder that contains `Agramon_withLAI_withETp.py`:
```bash
patch -p2 --dry-run < /path/to/agramon_refactor_phase1-2.patch
patch -p2 < /path/to/agramon_refactor_phase1-2.patch
```
`patch` will report one skipped hunk, which is `REFACTORING.md`; copy it to the repo root by hand. (`git apply -p2` does not work for this layout.)

After applying, run the acceptance test in section 5 before relying on the refactored code.
