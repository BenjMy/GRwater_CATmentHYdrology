# CATHY model looping over LAI — deforestation experiment

pyCATHY wrapper (clone of `run_real_case_pycathy.py` / `run_scenarios_pycathy.py`, with the monthly LAI loop of
`Agramon_withLAI_withETp.py`) built with `geoutils.uCATHY` and `geoutils.uET`.

> **Status: written, NOT run.** Nothing here has been executed (as requested). Follow
> [IMPLEMENTATION_STEPS.md](IMPLEMENTATION_STEPS.md) for the first smoke test and the list of things to verify.

## 1. Objective

Deforestation shows up in the Sentinel-2 LAI as a sudden, lasting drop. Question:

> **Can CATHY, driven only by rain, ETp and a time-varying LAI, recover the changes in actual evapotranspiration (ETa)
> and groundwater level (GWL) that go with the LAI drop?**

Experimental design (two paired simulations per CATHY version):

| scenario | LAI given to CATHY | meaning |
|---|---|---|
| `obs_*`  | observed monthly LAI (the drop is in the data) | what CATHY makes of the deforestation |
| `clim_*` | monthly median LAI of the reference years (`--ref-years`, default 2016-2018) | counterfactual: **same** forcing, **no** LAI drop |

`obs − clim` isolates the effect of the LAI change on ETa and water table (figures `eta_03`, `gwl_03`).
Both are then compared with the ground truth (§5) and with different CATHY versions (§6).

## 2. How LAI enters CATHY

For every calendar month (one CATHY run per month, restarted from the previous psi with `update_ic(INDP=1)`):

1. **LAI cube** (`lai_loop_utils.build_lai_cube`): acquisitions of the month → pixel-wise median (clouds = NaN) → mean of
   the pixels inside each CATHY cell. Gaps are carried forward / back-filled / nearest-cell filled.
2. **Vegetation map**: LAI classes (`--lai-thresholds`, default `0.5 1.5 3.0 4.5` → 5 classes, fixed for the whole run)
   written with `simu.update_veg_map(veg, maxveg=K)`.
3. **Root depth per class**: `ZROOT = --zroot × --zroot-factors[class]` (default 1.5 m × `0.1 0.3 0.6 0.9 1.0`) in the
   Feddes table (`update_soil`). A clear-cut cell moves to class 1 → root zone shrinks → transpiration is cut.
4. **ETp (atmbc)**, `--etp-mode`:
   * `product` (default): gridded `ET_0_wet` of `*SW-ET-gf.nc`, per cell (already canopy-aware, it comes from the
     satellite energy balance). LAI acts through the root zone only.
   * `etr_kc`: `ETr(meteo) × Kc(LAI)`, `Kc = kc_min + (kc_max − kc_min)(1 − exp(−k·LAI))`: explicit LAI → ETp coupling.
5. **Rain**: `TP-DD` of `*daily-meteo.nc` (field mean, spread uniformly over the day).
6. atmbc = rain − ETp per node, daily steps, `HSPATM=0, IETO=1`; closed box (optional lateral outlet `--outlet-flux`).

## 3. Files

| file | role |
|---|---|
| `run_lai_loop_pycathy.py` | one CATHY project / scenario: mesh, LAI cube, monthly loop, outputs, resume |
| `scenarios_lai_pycathy.py` | scenario table (CATHY version × LAI mode), runs the runner, then the comparison |
| `compare_obs_pycathy.py` | CATHY vs Lanžhot EC ETa, piezometers (and TSEB satellite ETa of the polygon) |
| `lai_loop_utils.py` | helpers (bbox-restricted NetCDF reader, pixel→cell aggregation, LAI cube, water table from psi, DEM, obs readers, metrics) |
| `IMPLEMENTATION_STEPS.md` | step tracking / checklist |

geoutils functions re-used: `uCATHY.GridSpec / polygon_mask / cells_to_nodes / node_index_from_mesh / check_date_overlap /
log_simulation`, `uET.field_daily_series / norm_units / _grid_crs`.

## 4. Usage

```bash
cd lib/test_tomas/lai_loop          # (data dir default: .../lib/test_tomas)

# smoke test: 2 months, coarse mesh, one scenario
python scenarios_lai_pycathy.py --only obs_withIrr --scenarios obs_withIrr --no-compare --months-limit 2 --dx 100 --spinup-cycles 0

# real run: observed vs no-LAI-drop, 2016-2024, + comparison with observations
python scenarios_lai_pycathy.py --start 2016-01-01 --end 2024-12-31 --dx 50

# another CATHY version
python scenarios_lai_pycathy.py --scenarios obs_withIrr obs_v1      # edit SCENARIOS in the file first (python ... --list shows them)

python scenarios_lai_pycathy.py --compare-only --cmp-args "--piezo-xy X Y"
```

Runner outputs per scenario (`<project>/outputs/`): `timeseries_daily.csv` (ETa, ETp, rain, water table mean/p10/p90, Sw/psi
of 3 layers, LAI), `gwl_nodes.nc` (water-table depth per surface node and day), `eta_nodes_monthly.nc`,
`lai_cells_monthly.nc`, `lai_monthly.csv`, `run_config.json`, `fig_overview.png`, and `_monthly_backup/` (resume: re-running
with identical options skips the months already done; `--no-resume` wipes it).

Comparison outputs (`<base>/comparison_obs/`): `eta_01_timeseries`, `eta_02_seasonal_scatter`, `eta_03_lai_effect`,
`gwl_01_timeseries`, `gwl_02_scatter_anomaly`, `gwl_03_lai_effect`, `metrics_eta.csv`, `metrics_gwl.csv`, `summary.md`.

## 5. Data found and what it implies (checked on the files you attached)

| data | content | period |
|---|---|---|
| `*_s2-l2b.nc` / `*_lcd-l2b.nc` | LAI (+ leaf traits), 502×613 px | 2015-12 → 2026-01 (1384 steps) / 2016-01 → 2025-12 (898) |
| `*SW-ET-gf.nc` | ET_0 / T_0, dry & wet | 2016-01-01 → 2025-12-31, daily |
| `*TSEB-ET-gf.nc` | ETa, T | 2016-01-01 → 2025-12-31, daily |
| `*daily-meteo.nc` | ETr, Tmin/Tmax, SW, EA, WS, PA, TP-DD (coarse 37×44) | 2016-01-01 → 2025-12-31 |
| `Lanzhot_ETa_2015-2020.csv` | EC ETa `Date;LNZ`, mm/day, 1826 days | **2015-01-01 → 2019-12-31** (not 2020, despite the name) |
| `GWL_LANZHOT.xlsx` | 4 piezometers, 30 min, m rel. to ground | **2020-01-01 → 2024-12-31** |

Consequences:

* `GWL_LANZHOT.xlsx` is **not an Excel file**: it is a `;`-separated text file. `read_obs_gwl` reads both.
* **EC ETa (2015-2019) and piezometers (2020-2024) do not overlap.** Model period 2016-2024 is needed to cover both
  (ETa on 2016-2019, GWL on 2020-2024). 2015 EC data fall before the NetCDF forcing (starts 2016-01-01).
* **Site mismatch to confirm**: the NetCDF files are `CeskeSvycarsko` (Bohemian Switzerland) while the observations are
  from **Lanžhot** (South Moravia, ~300 km away). `compare_obs_pycathy.py` prints the distance and a warning, and also
  compares against the TSEB satellite ETa *of the modelled polygon* (consistent location). See open questions.
* Piezometer datum ≠ model datum, so GWL is compared both as absolute level and as **anomaly** (means removed).

## 6. CATHY versions

`--version` is passed to `pyCATHY.cathy_tools.CATHY(version=...)`. In the table of `scenarios_lai_pycathy.py`:
`withIrr` (the one in `run_real_case_pycathy.py`) and `1.0.0` (pyCATHY default, **placeholder**) — other versions
to be added once you tell me which. The ET output is read through `simu.read_outputs("ET")`, so it is version-independent as far as pyCATHY
abstracts the file (`fort.777` vs `fort.1178`); this must be verified for the non-withIrr versions.

## 7. Assumptions and known limitations

* **ETa artefact**: with atmbc = rain − ETp, CATHY cannot separate rain from ET when net flux ≥ 0. As in Agramon,
  ETa is floored at ETp for those node-days (`eta_mm_d`); the raw value is kept (`eta_raw_mm_d`; `--no-eta-floor`).
* **Rain units**: `TP-DD` has empty name and units `?` in the file. Assumed mm/day; the runner prints the annual totals
  (expect ~700-1000 mm). Use `--rain-factor` if not.
* **Mesh depth**: default template (~3 m). Piezometers are 0-2.5 m deep, so the water table may leave the mesh
  (reported as mesh base depth). `--base` / `--wt-height` control it.
* **Soil, Feddes and root parameters** are pyCATHY defaults (forest soil NOT calibrated); `ZROOT` factors and LAI
  thresholds are first guesses (Agramon values were for crops/shrubs, LAI 0-3).
* Monthly restart through `update_ic(INDP=1)`: mesh-refresh issue noted in Agramon (map_prop2mesh) is not handled here.
* One uniform rain series for the polygon (meteo grid is coarse); ETp is spatially variable.
* ETa/GWL aggregated over the polygon nodes; piezometer-specific comparison needs `--piezo-xy`.
* Runtime: 9 years = 108 monthly CATHY runs (+ spin-up) per scenario; start with `--dx 100` and `--months-limit`.

## 8. Open questions (see also the chat)

1. **CATHY versions**: your message ends with "(see )" — which versions / branches / links do you want in `SCENARIOS`?
2. **Site**: model files are Česke Švýcarsko, observations are Lanžhot — intended? If Lanžhot is the target I need its LAI/ET/meteo NetCDFs and polygon.
3. **Polygon / DEM**: `ceske_svycarsko.geojson` and the DEM (7 kB, very coarse?) were not attached; the DEM is used as-is (`--dem`).
4. Which **LAI product** (`s2-l2b` or `lcd-l2b`) and which **deforestation date** / reference years?
5. Rain variable `TP-DD` units?
