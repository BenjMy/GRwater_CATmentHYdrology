# Simulation Log — scenario_8_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2_ks2L-d0.5m-b1e-5 (dem-plot 3)

_Started: 2026-11-20 14:24:54_


## Known fixes applied

- ETa net-forcing artefact floor: `apply_eta_artefact_floor()` adds an `ACT. ETRA_patched` column, flooring cells where CATHY reports `ACT. ETRA ≈ 0` because net forcing (rain − ETp) was ≥ 0 that day/node. Cells with net < 0 (ETp > rain) are set per `--eta-net-neg-mode` (= `etp`; `etp` → patched ETa = ETp). Original `ACT. ETRA` is left unmodified.

## Run configuration

- **Forcing source**: ETp time series (CSV) + satellite rain
- **Period**: 2016 – 2025
- **Target CRS**: EPSG:25830
- **ZROOT / PMIN / SCF**: 1.0 m / -50.0 m / 1 (Feddes)
- **ZROOT factors (bare/sparse/moderate/dense/very_dense)**: 0.3,0.6,2,3,5
- **PERMX / Ks override**: PERMX=6.000e-04 m/s, PERMY ratio=1.0, PERMZ ratio=1.0
- **Feddes (FP_map) override**: pyCATHY default
- **Initial condition**: uniform pressure head -100.0 m (INDP=0)
- **DEM source (tif)**: /home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/DTMplots/20250618_AGRAMON100m_micasense_dtm.tif
- **DEM plot (fid_1 clip)**: 3
- **Spin-up**: disabled
- **Outlet flux**: 0.0e+00 m/s
- **Outlet side**: auto-detect
- **Bottom flux**: disabled (closed)
- **Solver dtmin/dtmax/deltat**: 0.1 / 10000.0 / 1000.0 s
- **CSV run registry**: /home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/Agramon_SShydro/hypotheses/hyp1_layers_Ks/simulation_log.csv
- **sim_index**: 8

## Setup

- **DEM shape (rows x cols)**: (32, 9)
- **DEM CRS**: EPSG:25830
- **DEM pixel size (dx, dy)**: (5.000, 5.000) m
- **DEM origin (xllcorner, yllcorner)**: (615858.5, 4255338.2)
- **DEM valid area**: 288 px  ≈  0.72 ha (7200 m²)
- **DEM smoothing (gaussian sigma, px)**: none
- **LAI CRS**: EPSG:25830
- **Shapefile CRS**: EPSG:25830
- Grid: 9 rows × 32 cols
- **Mesh x range**: [615858.5, 615903.5] m
- **Mesh y range**: [4255338.0, 4255498.0] m
- **Mesh CRS**: EPSG:25830

## Atmospheric forcing — provenance & units

- **ETp source**: /home/z0272571a@campus.csic.es/Nextcloud/GRwater/data/satellite/Agramon/agramon_treatments_ET.csv (feature_id=2, column=ET_0_wet; same value for every pixel)
- **Rain source file**: /home/z0272571a@campus.csic.es/Nextcloud/GRwater/data/satellite/Agramon/20161001_20241231_TP.nc
- **Native CRS (source)**: EPSG:32630
- **Target CRS (mesh)**: EPSG:25830
- **ETp native pixel size (dx, dy)**: (20.0, -20.0) m
- **Rain native pixel size (dx, dy)**: (20.0, -20.0) m
- **Surface mesh nodes forced**: 242
- **Forcing period**: 2016-10-01 → 2025-08-31 (3257 daily steps)
- **Units**: ETp/rain source: mm/day; ATMBC/CATHY forcing: m/s (spatial mean over mesh surface nodes, negative = ET loss)
- **ATMBC steps (ETp)**: 3257
- **Net flux range (m/s)**: [-1.25e-07, 6.84e-07]
- **ETp range (mm/day)**: [0.34, 10.78]
- **Rain range (mm/day)**: [0.00, 59.75]
- **Total rain over period (spatial mean)**: 3650.9 mm
- **Total ETp over period (spatial mean)**: 12795.4 mm
- **Resolved outlet side**: ymin
- **Outlet diagnostic node**: 239
- **Bottom flux**: disabled (closed)
- ds_LAI_raster shape: (120, 32, 9)

## Spin-up (warm-up)

- Disabled (`--no-spinup`) — starting cold from the uniform IC.

## Monthly simulation loop


| Month | LAI mean | LAI min–max | Net flux mean (mm/day) | Veg zones | Vegetated % | LAI updated | Solver |
|---|---|---|---|---|---|---|---|
| 2015-11 | — | — | — | — | — | — | skipped (no forcing) |
| 2015-12 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-01 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-02 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-03 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-04 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-05 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-06 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-07 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-08 | — | — | — | — | — | — | skipped (no forcing) |
| 2016-09 | — | — | — | — | — | — | skipped (no forcing) |
- 2016-10: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2016-10 stats — psi mean: -89.668 m, sw mean: 0.301, ETa mean: 1.117e-08 m/s, recharge mean: n/a
- 2016-11: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2016-11 stats — psi mean: -84.222 m, sw mean: 0.311, ETa mean: 4.437e-09 m/s, recharge mean: n/a
- 2016-12: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2016-12 stats — psi mean: -48.356 m, sw mean: 0.381, ETa mean: 6.428e-09 m/s, recharge mean: n/a
- 2017-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-01 stats — psi mean: -42.876 m, sw mean: 0.364, ETa mean: 1.430e-08 m/s, recharge mean: n/a
- 2017-02: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2017-02 stats — psi mean: -44.771 m, sw mean: 0.347, ETa mean: 1.541e-08 m/s, recharge mean: n/a
- 2017-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-03 stats — psi mean: -49.258 m, sw mean: 0.351, ETa mean: 2.317e-08 m/s, recharge mean: n/a
- 2017-04: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-04 stats — psi mean: -60.772 m, sw mean: 0.334, ETa mean: 1.265e-08 m/s, recharge mean: n/a
- 2017-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-05 stats — psi mean: -60.474 m, sw mean: 0.332, ETa mean: 2.241e-08 m/s, recharge mean: n/a
- 2017-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-06 stats — psi mean: -69.690 m, sw mean: 0.320, ETa mean: 1.454e-08 m/s, recharge mean: n/a
- 2017-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-07 stats — psi mean: -71.867 m, sw mean: 0.318, ETa mean: 1.131e-08 m/s, recharge mean: n/a
- 2017-08: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-08 stats — psi mean: -68.800 m, sw mean: 0.320, ETa mean: 9.138e-09 m/s, recharge mean: n/a
- 2017-09: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-09 stats — psi mean: -55.484 m, sw mean: 0.333, ETa mean: 2.101e-08 m/s, recharge mean: n/a
- 2017-10: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2017-10 stats — psi mean: -66.263 m, sw mean: 0.318, ETa mean: 7.771e-09 m/s, recharge mean: n/a
- 2017-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-11 stats — psi mean: -68.135 m, sw mean: 0.317, ETa mean: 3.274e-09 m/s, recharge mean: n/a
- 2017-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-12 stats — psi mean: -49.369 m, sw mean: 0.326, ETa mean: 9.777e-09 m/s, recharge mean: n/a
- 2018-01: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-01 stats — psi mean: -58.013 m, sw mean: 0.322, ETa mean: 5.670e-09 m/s, recharge mean: n/a
- 2018-02: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 4840 cell(s) set for net < 0 (mode=etp).
- 2018-02 stats — psi mean: -41.025 m, sw mean: 0.343, ETa mean: 1.267e-08 m/s, recharge mean: n/a
- 2018-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2018-03 stats — psi mean: -46.966 m, sw mean: 0.332, ETa mean: 1.472e-08 m/s, recharge mean: n/a
- 2018-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2018-04 stats — psi mean: -63.948 m, sw mean: 0.318, ETa mean: 9.856e-09 m/s, recharge mean: n/a
- 2018-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-05 stats — psi mean: -67.620 m, sw mean: 0.316, ETa mean: 9.808e-09 m/s, recharge mean: n/a
- 2018-06: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-06 stats — psi mean: -60.844 m, sw mean: 0.321, ETa mean: 1.697e-08 m/s, recharge mean: n/a
- 2018-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2018-07 stats — psi mean: -71.950 m, sw mean: 0.313, ETa mean: 9.668e-09 m/s, recharge mean: n/a
- 2018-08: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-08 stats — psi mean: -70.352 m, sw mean: 0.313, ETa mean: 7.626e-09 m/s, recharge mean: n/a
- 2018-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-09 stats — psi mean: -56.674 m, sw mean: 0.321, ETa mean: 1.459e-08 m/s, recharge mean: n/a
