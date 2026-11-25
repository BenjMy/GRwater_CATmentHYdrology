# Agramon/Hellín Catchment Hydrology

> **⚠️ This is a living repository.** Code, data paths, and documentation are
> actively evolving alongside the 2026 TFG work — expect breaking changes,
> `[FILL IN]` placeholders, and hard-coded paths (§4.3) until the first
> tagged release. Check [`CHANGELOG.md`](./CHANGELOG.md) for what has
> changed between releases, and open an issue if something you relied on
> shifted under you.

Physically-based, long-term (multi-year) subsurface hydrology simulation of the
Agramon plots (Hellín, Albacete, Spain), driven by satellite-derived
potential evapotranspiration (ETp) and rainfall and by a monthly,
satellite-derived Leaf Area Index (LAI) that dynamically updates the
vegetation/root map during the run. Built on
[pyCATHY](https://github.com/BenjMy/pycathy_wrapper), the Python wrapper for
the CATHY (CATchment HYdrology) coupled surface/subsurface model.

This repository is part of the **GRwater** project
(<https://grwater.ica.csic.es/>) and its outputs feed into a 2026 TFG
(*Trabajo de Fin de Grado*) report `[FILL IN: title / student / degree]`.

---

## 1. Overview (Findable / citable)

| | |
|---|---|
| **Project** | GRwater — <https://grwater.ica.csic.es/> |
| **Study site** | Agramon plots, Hellín catchment (Albacete, Spain) — 9 candidate DEM sub-plots (`--dem-plot 1`–`9`) |
| **Model** | CATHY (Camporese, Paniconi, Putti & Orlandini, 2010, *Water Resour. Res.* 46, <https://doi.org/10.1029/2008WR007536>) via the [pyCATHY](https://github.com/BenjMy/pycathy_wrapper) Python wrapper (BSD-3-Clause) |
| **Institutions** | ICA-CSIC (Tech4Agro group); UCLM |
| **License** | GNU General Public License v3.0 — see [`LICENSE`](./LICENSE) |
| **Long-term archive** | [Zenodo](https://zenodo.org) — the repository is connected to Zenodo so that every GitHub release is automatically archived and issued a DOI. `[FILL IN once the first release is made: concept DOI badge/link, e.g. 10.5281/zenodo.XXXXXXX]` |
| **Funding** | Fundación BBVA (FBBVA) |

---

## 2. Contributors

| Name | Role | Affiliation | Contact | ORCID |
|---|---|---|---|---|
| **Hector N.** | PI — Tenured Scientist | Tech4Agro, ICA-CSIC | hector.nieto@ica.csic.es | [0000-0003-4250-6424](https://orcid.org/0000-0003-4250-6424) |
| **Benjamin M.** | Co-PI — Post-doc Researcher | Tech4Agro, ICA-CSIC | benjamin.mary@ica.csic.es | [0000-0001-7199-2885](https://orcid.org/0000-0001-7199-2885) — [website](https://benjmy.github.io/myCareerBook/intro.html) |
| **María B.** | Postdoc Researcher | Tech4Agro, ICA-CSIC | maria.burguet@ica.csic.es | [0000-0002-9748-8078](https://orcid.org/0000-0002-9748-8078) |
| **Manuel Esteban L.B.** | Professor | Universidad de Castilla-La Mancha (UCLM) | manuelesteban.lucas@uclm.es | [0000-0001-6270-8408](https://orcid.org/0000-0001-6270-8408) — [website](https://blog.uclm.es/manuelestebanlucas/) |


---

## 3. Repository structure

```
.
├── lib/                                     # simulation + post-processing code
│   ├── Agramon_withLAI_withETp.py           # main simulation driver (documented in §4)
│   ├── data_paths.toml                      # every INPUT path (ETp csv, LAI, EO ETa, DEM, shapefiles); see §7.1
│   ├── agramon_*.py                         # driver modules (config, forcing, LAI, Ks, mesh, spin-up, outputs, ...)
│   ├── plot_run_results.py                  # per-scenario post-processing (maps, hydro, monthly)
│   ├── plot_*.py                            # plotting modules (plot_results_shared = shared figure functions)
│   ├── analysis_ET_compare.py               # ETa vs EO energy-balance validation
│   ├── hypotheses.toml                      # the four hypotheses (variants, shared options)
│   ├── hyp1_layers_Ks.py .. hyp4_SW.py      # one runner per hypothesis (see §7.2)
│   ├── run_all_hypotheses.py                # runs / dry-runs / validates / plots all hypotheses
│   ├── scenario_runner.py                   # batch runner: crosses dem-plot × PMIN × outlet/bottom BC × spin-up
│   ├── run_scenarii_withETp_withLAI_ALL.sh  # deprecated wrapper around scenario_runner.py
│   └── plot_scenarii_withETp_withLAI.sh     # batch post-processor over sim_index / scenario folders
├── data/
│   ├── satellite/Agramon/
│   │   ├── 20161001_20241231_ET0.nc         # potential evapotranspiration (ETp), 2016-10-01–2024-12-31
│   │   └── 20161001_20241231_TP.nc          # rainfall (precipitation), same period
│   └── biophysical/agramon/
│       └── input/
│           └── lai_monthly.nc               # monthly Leaf Area Index (LAI)
│   └── microcuencas_13/agramon/
│   └── DTMplots/
└── README.md
```

---

## 4. Code description

### 4.1 Main driver — `Agramon_withLAI_withETp.py`

Runs the full pyCATHY pipeline for one scenario (one dem-plot × one set of
boundary conditions/spin-up/PMIN):

1. Load satellite ETp + rain NetCDFs (native EPSG:32630 → reprojected to
   EPSG:25830), or optionally ERA5 (`USE_ETP_FORCING = False` in-script).
2. Load the monthly LAI dataset.
3. Log the run's parameters to a shared CSV registry (`sim_index`).
4. Load and pre-process the DEM for the selected plot.
5. Initialise and mesh the CATHY model.
6. Map ETp/rain grids onto mesh nodes and build the full ATMBC (atmospheric
   boundary condition) time series.
7. Configure soil properties, initial conditions, boundary conditions
   (lateral outlet + optional free-bottom-drainage Neumann flux), and
   output settings.
8. For each calendar month in the simulation period: reclassify the LAI
   raster into a discrete vegetation/root-depth zone map, update CATHY's
   vegetation input, and advance the solver one month.
9. Aggregate and write outputs to NetCDF, plus a per-run Markdown log.

Every run is self-registering: parameters are appended to
`simulation_log_withETp_withLAI.csv` with an auto-incremented `sim_index`,
and outputs land in a self-describing folder name, e.g.:

```
scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c_pmin-m10
```

so scenarios can be told apart from folder names alone, without
cross-referencing the CSV.

**Key CLI arguments** (`python Agramon_withLAI_withETp.py --help` for the
full, current list):

| Group | Flag | Meaning |
|---|---|---|
| Simulation | `--start-year`, `--end-year` | forcing period (end exclusive) |
| Simulation | `--zroot` | base root depth (m), scaled per vegetation zone |
| Simulation | `--pmin` | minimum pressure-head threshold (m) |
| Boundary conditions | `--outlet-flux`, `--outlet-side` | Neumann outflow on the topographically lowest mesh side (m/s); `--outlet-side` overrides auto-detection |
| Boundary conditions | `--bottom-flux` | optional free-bottom-drainage-style Neumann outflow on the mesh base (m/s); `0.0` = closed/no-flow (default) |
| Spin-up | `--no-spinup`, `--spinup-cycles`, `--spinup-year` | pre-2016 equilibration by cycling one forcing year |
| DEM | `--dem-plot {1..9}`, `--dem-folder` | which Agramon sub-plot to simulate |
| Solver | `--dtmin`, `--dtmax`, `--deltat` | CATHY time-stepping controls |
| Output | `--path2prj`, `--no-plots` | where scenario folders are created; skip matplotlib figures |

**Root-depth / vegetation classes** (from monthly LAI, thresholds
calibrated to the observed Agramon LAI range 0–0.6 m² m⁻²):

| Class | LAI range | Root-depth factor |
|---|---|---|
| 1 – bare soil | ≤ 0.10 | 0.05 × `--zroot` |
| 2 – sparse | 0.10–0.25 | 0.25 × `--zroot` |
| 3 – moderate | 0.25–0.40 | 1.0 × `--zroot` |
| 4 – dense | 0.40–0.55 | 1.5 × `--zroot` |
| 5 – very dense | > 0.55 | 2.0 × `--zroot` |

### 4.2 Batch scripts

- **`scenario_runner.py`** (the old `run_scenarii_withETp_withLAI.sh` is a wrapper around it) — crosses dem-plot(s) × PMIN
  value(s) × outlet BC (free/closed) × bottom BC (free/closed) × spin-up
  (on/off), calling `Agramon_withLAI_withETp.py` once per combination.
  Does not write a batch log file (the per-run CSV log is kept). Run
  `--help` for the full flag reference and worked examples.
- **`plot_scenarii_withETp_withLAI.sh`** — iterates a hard-coded list of
  `sim_index` values and calls `plot_run_results.py --scenario N` for
  post-processing/plotting. **The `SCENARIOS=(...)` list must be kept in
  sync by hand** with `simulation_log_withETp_withLAI.csv` (or the output
  folder names) whenever more than one BC/spin-up/PMIN combination has
  been run per plot — see the comment block at the top of the script.

### 4.3 Known reuse caveats (please read before rerunning)

- **External code dependencies not in this repo**: `Agramon_utils`
  (`AgUtils`) and `geoutils.uCATHY` are imported from `MODULE_PATH`
  (`Tech4agro_org/GRwater_geophy`) and are **not included here** —
  `[FILL IN: are these published elsewhere / should they be vendored into
  ./lib or cited as a separate repository?]`. (Plot shapefiles no longer need `GRwater_geophy`: they are
  read from `[terrain]` in `lib/data_paths.toml`.)

### 4.4 Run registry and comparing runs

The CSV run log (`simulation_log*.csv`) is the single registry of runs. Besides the parameters it
now records, per run, `status` (running / ok / failed), `seconds`, `n_months`, `sim_dir` (relative to
`path2prj`, so the project can be moved), `param_hash`, `git_commit`, `error_tail` of a failed run and,
for hypothesis runs, the `hypothesis` / `variant` names. It is written by `agramon_registry.py` (called
from `agramon_runlog.log_simulation` and `run_pipeline`); a registry error never stops a run. Old logs
still load: `python backfill_registry.py --study <study>` fills `param_hash`, `sim_dir`, `status` and
`n_months` for them (a backup is written unless `--no-backup`).

- **`run_catalog.py`** — one table of the runs (status, duration, months, outputs on disk, duplicates,
  which parameters vary, optional ETa / recharge / closure): `python run_catalog.py --study <study> --metrics`,
  filters `--where "PMIN == -50"`, `--status ok`, `--unique`, `--csv runs_summary.csv`.
- **`plot_compare_runs.py`** — comparison figures (response, annual totals, water balance, optional ETa
  vs EO-EB with `--eo`) for runs picked by `--sim 2 4 6` or `--where`; labels show only the parameters that
  differ. Default selection: finished, not partial, latest of each parameter set.
- **Reuse**: `scenario_runner.py` and the hypothesis runner skip a run only when a matching one finished.
  Runs marked failed/running never count, and neither do **partial** runs (fewer than 90 % of the months of
  the longest run with the same period): the pipeline writes all NetCDFs even when months were skipped, so
  "outputs exist" alone is not proof that a run completed. `agramon_registry.find_run(log, params, path2prj)`
  does the same lookup by `param_hash` from Python.
- The hash covers the logged parameters only. It does not cover code version (see `git_commit`) or options
  that are not in the log (e.g. `--dem-smooth-sigma`), so identical hashes are a strong hint, not a proof.
- `python lib/test/test_registry.py` runs the registry smoke test (no CATHY needed).

### 4.5 Plot and compare results

Three levels, from one run to many. All commands run from `lib/`. A *study* is a hypothesis folder,
e.g. `<WB_PRJ_ROOT>/hypotheses/hyp1_layers_Ks` (contains `prj/`, `simulation_log.csv`, `manifest.json`).

`--study` (in `run_catalog.py`, `plot_compare_runs.py`, `post_process.py`, `clean_log.py`) takes either that
folder or just the hypothesis name: `--study hyp4` or `--study hyp4_SW`. The name is looked up under the
hypotheses `out_root` (`[general] out_root` in `hypotheses.toml`, else `<wb_prj_root>/hypotheses` from
`data_paths.toml`); a unique prefix is enough, and an ambiguous or unknown name prints the available studies.
Print the folder in use with `python -c "from hypotheses_common import default_out_root as d; print(d())"`.

**a) One run: figures per scenario**

```bash
python post_process.py --scenario 8          # maps, hydrographs, monthly, ETa vs EO-EB (Fig I)

# a run of a hypothesis study (log and prj/ come from the study)
python post_process.py --study hyp4 --list              # sim_index values of the study
python post_process.py --study hyp4 --scenario 0 1
```

Figures go to `../Analysis/ET_compare/<scenario_dir>/`. Needs the EO ETa file and the shapefile (§7.1).
Without `--study` the script reads `simulation_log_withETp_withLAI.csv` in the current folder. `--dem-res`
(default 5) must equal the run's `resample_resolution`; `--forcing etp|etp_ts` restricts the default
selection to one forcing (hyp4 mixes both).

**b) Several runs: side-by-side comparison from the run log**

```bash
# list what exists first (status, months, which parameters vary)
python run_catalog.py --study <study> --metrics

# compare everything that finished (latest of each parameter set)
python plot_compare_runs.py --study <study>

# pick runs, filter, add the EO energy-balance comparison (slower)
python plot_compare_runs.py --study <study> --sim 2 4 6 --out figs/ks --title "1 vs 2 layers"
python plot_compare_runs.py --study <study> --where "PMIN == -50" --fire-date 2022-08-01
python plot_compare_runs.py --study <study> --eo
```

Without a study: `--log-file simulation_log_withETp_withLAI.csv --path2prj ../prj`.
Output in `<study>/figures/compare_sim-<ids>/`:

| File | Content | Needs |
|---|---|---|
| `response.png` | ETa, soil moisture, recharge, psi per run | `monthly_stats.csv` |
| `annual.png` | annual ETa / recharge / runoff totals | monthly balance |
| `balance.png` | mean annual ETa, recharge, runoff, closure residual | monthly balance |
| `eta_eo.png`, `skill.png`, `metrics.csv` | ETa vs EO-EB, KGE / NSE / RMSE per run | `--eo` (shapefile + EO file) |
| `runs.csv` | runs plotted, with the labels used in the legends | always |

Labels show only the parameters that differ between the runs (`#4 permx=0.0006`). Failed, running and
partial runs are skipped unless `--include-partial` / `--sim` names them; `--all` keeps duplicates.

Figure layout: the legend is one block under the plots (long labels are wrapped, the figure grows to fit it);
in `balance.png` the labels are wrapped on the y axis. Variants of one `group` (hypotheses.toml) share a hue and
differ in shade, dash pattern and line width, so identical curves stay visible; with more groups than the 7 base
colours, the colours repeat and the dash pattern changes.

**c) A whole hypothesis**

```bash
python run_all_hypotheses.py --only hyp1_layers_Ks   # runs missing variants, validates, plots
python hyp1_layers_Ks.py --dry-run                   # variant matrix only, runs nothing
```

Writes `metrics.csv`, `water_balance_monthly.csv`, `annual_totals.csv`, `variant_summary.csv`,
`study_report.md` and `figures/` into the study folder (§7.3).

**Troubleshooting `--eo` / validation: `DataSourceError ... microcuencas_13.shp: No such file`.**
Shapefile paths come from `lib/data_paths.toml` (`[terrain] shapefile_dir`, relative to the repo root,
default `data/microcuencas_13`). Check with:

```bash
python -c "import agramon_config as C; print(C.SHAPEFILE_DIR, (C.SHAPEFILE_DIR/'microcuencas_13.shp').exists())"
```

If this prints `False`, fix `shapefile_dir` in your TOML (or point `AGRAMON_DATA_CONFIG` at your own
copy). Older versions looked in `GRwater_geophy/shapefiles` regardless of the TOML; that is fixed.

---

## 5. Data description

| Dataset | Path | Format | Native CRS | Coverage | Notes |
|---|---|---|---|---|---|
| Potential ET (ETp) | `data/satellite/Agramon/20161001_20241231_ET0.nc` | NetCDF | EPSG:32630 | 2016-10-01 – 2024-12-31 | satellite-derived potential evapotranspiration; reprojected to EPSG:25830 at load time |
| Rainfall | `data/satellite/Agramon/20161001_20241231_TP.nc` | NetCDF | EPSG:32630 | 2016-10-01 – 2024-12-31 | reprojected to EPSG:25830 at load time |
| Leaf Area Index (LAI) | `data/biophysical/agramon/input/lai_monthly.nc` | NetCDF | EPSG:32630 | `[FILL IN: date range]` | monthly composite; drives the vegetation/root-depth zone map (§4.1) |
| DEM (Agramon, 100 m) | `data/DTMplots/` (`[terrain] dem_tif` in `lib/data_paths.toml`) | GeoTIFF / `.adf` | EPSG:25830 | — | sub-plot selected via `--dem-plot 1`–`9` |
| Plot boundaries (`microcuencas_13`) | `data/microcuencas_13/microcuencas_13.shp` (`[terrain] shapefile_dir` / `shapefile_names` in `lib/data_paths.toml`) | Shapefile | EPSG:25830 | — | read through `agramon_config.SHAPEFILE_DIR`; no `GRwater_geophy` checkout needed |
| ERA5 (optional alternate forcing) | external, not in `./data` | NetCDF | — | `[FILL IN]` | Copernicus C3S; used only if `USE_ETP_FORCING = False` |

**Model outputs** written per scenario folder (not versioned in this
repository — regenerated by running the code):

| File | Content |
|---|---|
| `veg_map_history.pkl` / `veg_map_history.nc` | monthly vegetation-zone raster time series (1=bare … 5=very dense) |
| `et_output.nc` | actual evapotranspiration (ETa) per node/time, incl. an artefact-floor correction flag when `USE_ETP_FORCING=True` |
| `psi_output.nc` | pressure head (ψ) per node/time |
| `sw_output.nc` | soil water content per node/time |
| `recharge_output.nc` | recharge per node/time |
| `simulation_log.md` | human-readable per-run narrative log (config, per-month table, warnings, run summary) |
| `simulation_log_withETp_withLAI.csv` | shared registry of every run's parameters, keyed by `sim_index` |

---

## 6. Software requirements

No pinned environment file is provided; install into a fresh Python
environment (`[FILL IN: minimum Python version, e.g. ≥3.10]`, matching
pyCATHY's requirement):

- [pyCATHY](https://github.com/BenjMy/pycathy_wrapper) (BSD-3-Clause) —
  includes the CATHY Fortran core, no separate CATHY install needed
- `numpy`, `pandas`, `xarray`, `rioxarray`, `scipy`, `matplotlib`,
  `geopandas` (for shapefile loading)
- Internal/site-specific modules: `Agramon_utils`, `geoutils` (`uCATHY`)
  — see §4.3

---

## 7. How to reproduce a run

```bash
cd lib/

# Single scenario, default plot (3) and forcing period:
python Agramon_withLAI_withETp.py

# A specific plot, boundary conditions and PMIN:
python Agramon_withLAI_withETp.py --dem-plot 2 --outlet-flux=-1e-6 \
    --bottom-flux=0 --pmin=-10 --start-year 2019 --end-year 2022

# Batch: cross several plots/PMIN values and both outlet/spin-up settings
./run_scenarii_withETp_withLAI.sh --dem-plot 1 2 3 --pmin -5 -10 -20 \
    --outlet both --spinup both

# Post-process a known list of sim_index values
./plot_scenarii_withETp_withLAI.sh
```

Before running: repoint the input paths (§7.1) at the `./data` layout in §3,
or supply your own data at those locations.

### 7.1 Input paths: `lib/data_paths.toml`

All input data paths live in `lib/data_paths.toml` (read by `agramon_config.py`, imported by the
driver, the post-processing and the hypothesis runners): satellite ETp/rain rasters, the
**ETp time series CSV** (`forcing.etp_csv`, columns `ET_0_wet` / `ET_0_dry`, used by
`--forcing etp_ts`), the EO energy-balance ETa, LAI, DEM, shapefiles, CATHY project root.
To run another catchment copy the file, edit it and point to it:

```bash
export AGRAMON_DATA_CONFIG=/path/to/my_site.toml
```

Precedence: `AGRAMON_*` environment variables > data_paths.toml > built-in defaults.

### 7.2 Two-layer Ks

`--permx-bottom <m/s> --ks-layer-depth <m>` sets a two-layer soil: `--permx` above the depth
(CATHY layers whose centre is shallower), `--permx-bottom` below. Compatible with
`--ks-fire-mode` (the post-fire Ks drives the top layer).

### 7.3 Hypotheses

Each hypothesis is one Python file reading its table of `lib/hypotheses.toml`:

| file | question | variants |
|---|---|---|
| `hyp1_layers_Ks.py` | 1 vs 2 layers, static vs LAI-driven Ks | 9 |
| `hyp2_LAI_SCF.py` | SCF 0/1 (Feddes / PMIN threshold) x PMIN -15/-150/-1e35 | 6 |
| `hyp3_ic.py` | initial pressure head / water-table depth, with/without spin-up | 9 |
| `hyp4_SW.py` | ETp forcing: wet vs dry Shuttleworth-Wallace, satellite ETp | 3 |

```bash
python run_all_hypotheses.py --dry-run          # list every run + variant-matrix figure, runs nothing
python run_all_hypotheses.py --only hyp1_layers_Ks hyp3_ic    # run + validate + plot
python hyp2_LAI_SCF.py --dry-run                # one hypothesis, same options
```

Each hypothesis writes a single folder `<out_root>/<hyp>/` (default
`<WB_PRJ_ROOT>/hypotheses`) with `prj/` (all its CATHY runs), `simulation_log.csv`,
`manifest.json`, `metrics.csv` (bias, RMSE, r, NSE, KGE of each run vs the EO energy-balance
ETa) and `figures/`. Finished runs are skipped unless `--force`.
Independently of the EO reference (`--no-balance` to skip) it also writes, from each run's
`monthly_stats.csv` and `hgraph_output.csv`: `water_balance_monthly.csv`, `annual_totals.csv`,
`period_totals.csv` (all / pre-fire / post-fire totals and mm/yr of rain, ETa, recharge, runoff, closure
residual), `variant_summary.csv` (mean ETa, recharge, saturation per variant), balance figures and
`study_report.md` (failed or skipped items). `run_all_hypotheses.py` joins skill and means in
`summary_metrics.csv`. Runoff is in m3 unless `[validation] catchment_area_m2` is set; the closure residual
is flux-only (no storage change). `recovery_metrics.csv` and per-variant recovery figures appear once
`[validation] control_sim_dir` points to a control-catchment run; until then the report says it was skipped.
Limitation: only the LAI/SCF CATHY build exists in the pipeline, so hyp2 varies SCF and PMIN
but cannot run a classic CATHY without LAI.

---

## 8. How to cite

Once the first GitHub release/Zenodo DOI exists, cite this repository as:

`[FILL IN: full citation
once available]`

Please also cite the underlying model:

> Camporese, M., Paniconi, C., Putti, M., Orlandini, S. (2010).
> Surface-subsurface flow modeling with path-based runoff routing,
> boundary condition-based coupling, and assimilation of multisource
> observation data. *Water Resources Research*, 46.
> <https://doi.org/10.1029/2008WR007536>

and, if used, the associated 2026 TFG report: `[FILL IN: full citation
once available]`.

---

## 9. License

Code in this repository is released under the **GNU General Public
License v3.0** — see [`LICENSE`](./LICENSE). `[FILL IN: confirm the data
in ./data should carry the same GPL-3.0 license, or should it instead use
a data-appropriate license such as CC-BY-4.0?]`

---

## 10. Funding

This project has received funding from **Fundación BBVA (FBBVA)**.
`[FILL IN: grant/call name or reference number]`

---

