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
│   ├── withLAI_results.py                   # per-scenario post-processing ([FILL IN: confirm this lives in ./lib])
│   ├── run_scenarii_withETp_withLAI.sh      # batch runner: crosses dem-plot × PMIN × outlet/bottom BC × spin-up
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

- **`run_scenarii_withETp_withLAI.sh`** — crosses dem-plot(s) × PMIN
  value(s) × outlet BC (free/closed) × bottom BC (free/closed) × spin-up
  (on/off), calling `Agramon_withLAI_withETp.py` once per combination.
  Logs everything to a timestamped, argument-tagged `.log` file. Run
  `--help` for the full flag reference and worked examples.
- **`plot_scenarii_withETp_withLAI.sh`** — iterates a hard-coded list of
  `sim_index` values and calls `withLAI_results.py --scenario N` for
  post-processing/plotting. **The `SCENARIOS=(...)` list must be kept in
  sync by hand** with `simulation_log_withETp_withLAI.csv` (or the output
  folder names) whenever more than one BC/spin-up/PMIN combination has
  been run per plot — see the comment block at the top of the script.

### 4.3 Known reuse caveats (please read before rerunning)

- **External code dependencies not in this repo**: `Agramon_utils`
  (`AgUtils`) and `geoutils.uCATHY` are imported from `MODULE_PATH`
  (`Tech4agro_org/GRwater_geophy`) and are **not included here** —
  `[FILL IN: are these published elsewhere / should they be vendored into
  ./lib or cited as a separate repository?]`.

---

## 5. Data description

| Dataset | Path | Format | Native CRS | Coverage | Notes |
|---|---|---|---|---|---|
| Potential ET (ETp) | `data/satellite/Agramon/20161001_20241231_ET0.nc` | NetCDF | EPSG:32630 | 2016-10-01 – 2024-12-31 | satellite-derived potential evapotranspiration; reprojected to EPSG:25830 at load time |
| Rainfall | `data/satellite/Agramon/20161001_20241231_TP.nc` | NetCDF | EPSG:32630 | 2016-10-01 – 2024-12-31 | reprojected to EPSG:25830 at load time |
| Leaf Area Index (LAI) | `data/biophysical/agramon/input/lai_monthly.nc` | NetCDF | EPSG:32630 | `[FILL IN: date range]` | monthly composite; drives the vegetation/root-depth zone map (§4.1) |
| DEM (9 Agramon sub-plots) | `[FILL IN — see §3]` | `.adf` raster | EPSG:25830 | — | selected via `--dem-plot 1`–`9` |
| Plot boundaries (`microcuencas_13`) | `[FILL IN — see §3]` | Shapefile | EPSG:25830 | — | |
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

Before running: repoint the hard-coded paths described in §4.3 at the
`./data` layout in §3, or supply your own data at those locations.

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

