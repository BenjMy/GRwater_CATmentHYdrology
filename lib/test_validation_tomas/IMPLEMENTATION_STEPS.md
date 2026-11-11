# Implementation tracker — CATHY LAI loop

Legend: ✅ written (not executed) · 🔲 to do / to verify by running · ❓ needs your input

## A. Analysis of the inputs
- ✅ Read `run_real_case_pycathy.py`, `run_scenarios_pycathy.py` (structure to clone: runner + scenario driver + comparison)
- ✅ Read `Agramon_withLAI_withETp.py` (monthly LAI loop, veg map, soil update, restart with `update_ic(INDP=1)`, ETa floor)
- ✅ Read `uCATHY.py` / `uET.py` (functions to reuse, see README §3)
- ✅ `synthetic_ksfire_LAI.py` and `synthetic_ksfire_LAI_1.py` are identical (`diff` empty); the Ks-fire part is not used here
- ✅ Inspect observations: `GWL_LANZHOT.xlsx` is a `;` CSV; EC ETa ends 2019-12-31; piezometers 2020-2024 (no overlap)
- ✅ Inspect NetCDF inventory (your `inspect-nc` output): dims `t,y,x`; LAI/ET on 502×613; meteo on 37×44; `TP-DD` units `?`
- ❓ Site mismatch (Česke Švýcarsko model data vs Lanžhot observations)

## B. Design
- ✅ One CATHY project per scenario; month-by-month loop; LAI → veg class → ZROOT; ETp product or ETr×Kc(LAI)
- ✅ Paired scenarios `obs_*` / `clim_*` (LAI drop vs no drop)
- ✅ Scenario table = CATHY version × LAI mode × extra options

## C. Code
- ✅ `lai_loop_utils.py`: `GridSource` (bbox-restricted lazy read), `CellAggregator`, `build_lai_cube`, `lai_climatology`,
  `lai_to_veg_map`, `kc_from_lai`, `water_table_depth`, `build_dem` (+ reprojection), obs readers, metrics
- ✅ `run_lai_loop_pycathy.py`: args → polygon/CRS → mesh (as runner) → LAI cube → ETp/rain → spin-up + monthly loop →
  harvest ETa / GWL / Sw / psi → resume backup → merged outputs + overview figure
- ✅ `scenarios_lai_pycathy.py`: scenario table, subprocess per scenario, then comparison
- ✅ `compare_obs_pycathy.py`: EC ETa, piezometers, TSEB; daily + monthly metrics; raw + anomaly GWL; LAI-effect figures
- ✅ `README.md` (objective, design, data findings, limitations, open questions)
- ✅ Python syntax check only (`py_compile`) — no execution

## D. First run — to do on your machine (suggested order)
1. 🔲 `python run_lai_loop_pycathy.py --no-run --months-limit 2 --dx 100` — checks polygon, CRS, DEM, mesh size, LAI cube, rain/ETp series, date overlap (prints annual rain totals → verify `TP-DD` units)
2. 🔲 Same without `--no-run`, `--spinup-cycles 0` — 2 CATHY months; check `eo_checks`-style prints, `outputs/timeseries_daily.csv`
3. 🔲 Verify the points I could not check without CATHY:
   - `simu.read_outputs("ET")` columns (`time_sec`, `ACT. ETRA`, `SURFACE NODE` or `X/Y`) and ETRA units (assumed m/s)
   - `simu.read_outputs("psi"/"sw")` columns = `nstr*nnod`, layer 0 = surface (the script stops if not)
   - `update_prepo_inputs(base=..., dr=...)`, `update_nansfdirbc(time=…, no_flow=True)` accepted by the chosen version
   - `update_ic(INDP=1, pressure_head_ini=...)` after a month on the new mesh (Agramon mesh-refresh note)
   - ETa ≤ ETp sanity and water-table depth vs mesh depth
4. 🔲 Full run 2016-2024, `obs_withIrr` + `clim_withIrr` (check `lai_monthly.csv`: the LAI drop must be visible)
5. 🔲 `compare_obs_pycathy.py`; check `eta_01`, `gwl_01`, `eta_03`, `gwl_03`
6. 🔲 Second CATHY version (`obs_v1`, `clim_v1`) once the version list is confirmed ❓
7. 🔲 Calibrate: LAI thresholds, `--zroot-factors`, `--pmin`, soil Ks, `--wt-height`, mesh depth

## E. Decisions to confirm ❓
- CATHY versions to include (your message: "see )" — link missing)
- Site: Česke Švýcarsko vs Lanžhot; polygon / DEM / piezometer coordinates
- LAI file (`s2-l2b` vs `lcd-l2b`), deforestation date, reference years
- Rain units; ETp choice (`ET_0_wet` vs `ET_0_dry`, `product` vs `etr_kc`)
