# Simulation Log — scenario_1_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2 (dem-plot 3)

_Started: 2026-11-25 10:21:39_


## Known fixes applied

- ETa net-forcing artefact floor: `apply_eta_artefact_floor()` adds an `ACT. ETRA_patched` column, flooring cells where CATHY reports `ACT. ETRA ≈ 0` because net forcing (rain − ETp) was ≥ 0 that day/node. Cells with net < 0 (ETp > rain) are set per `--eta-net-neg-mode` (= `etp`; `etp` → patched ETa = ETp). Original `ACT. ETRA` is left unmodified.

## Run configuration

- **Forcing source**: ETp time series (CSV) + satellite rain
- **Period**: 2016 – 2025
- **Target CRS**: EPSG:25830
- **ZROOT / PMIN / SCF**: 1.0 m / -50.0 m / 1 (Feddes)
- **ZROOT factors (bare/sparse/moderate/dense/very_dense)**: 0.3,0.6,2,3,5
- **PERMX / Ks override**: PERMX=1.000e-04 m/s, PERMY ratio=1.0, PERMZ ratio=1.0
- **Feddes (FP_map) override**: pyCATHY default
- **Initial condition**: uniform pressure head -100.0 m (INDP=0)
- **DEM source (tif)**: /home/z0272571a@campus.csic.es/Documentos/GitHub/GRwater_CATmentHYdrology/DTMplots/20250618_AGRAMON100m_micasense_dtm.tif
- **DEM plot (fid_1 clip)**: 3
- **Spin-up**: disabled
- **Outlet flux**: 0.0e+00 m/s
- **Outlet side**: auto-detect
- **Bottom flux**: disabled (closed)
- **Solver dtmin/dtmax/deltat**: 0.1 / 10000.0 / 1000.0 s
- **CSV run registry**: /home/z0272571a@campus.csic.es/Documentos/GitHub/GRwater_CATmentHYdrology/Agramon_SShydro/hypotheses/hyp1_layers_Ks/simulation_log.csv
- **sim_index**: 1

## Setup

- **DEM shape (rows x cols)**: (32, 9)
- **DEM CRS**: EPSG:25830
- **DEM pixel size (dx, dy)**: (5.000, 5.000) m
- **DEM origin (xllcorner, yllcorner)**: (615858.5, 4255338.2)
- **DEM valid area**: 288 px  ≈  0.72 ha (7200 m²)
- **DEM smoothing (gaussian sigma, px)**: none
- **LAI CRS**: EPSG:25830
- **Shapefile CRS**: EPSG:25830
