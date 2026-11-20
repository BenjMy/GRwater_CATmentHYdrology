# Simulation Log — scenario_4_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2 (dem-plot 3)

_Started: 2026-11-20 12:53:51_


## Known fixes applied

- ETa net-forcing artefact floor: `apply_eta_artefact_floor()` adds an `ACT. ETRA_patched` column, flooring cells where CATHY reports `ACT. ETRA ≈ 0` because net forcing (rain − ETp) was ≥ 0 that day/node. Cells with net < 0 (ETp > rain) are set per `--eta-net-neg-mode` (= `etp`; `etp` → patched ETa = ETp). Original `ACT. ETRA` is left unmodified.

## Run configuration

- **Forcing source**: ETp time series (CSV) + satellite rain
- **Period**: 2016 – 2025
- **Target CRS**: EPSG:25830
- **ZROOT / PMIN / SCF**: 1.0 m / -50.0 m / 1 (Feddes)
- **ZROOT factors (bare/sparse/moderate/dense/very_dense)**: 0.3,0.6,2,3,5
- **PERMX / Ks override**: PERMX=1.000e-03 m/s, PERMY ratio=1.0, PERMZ ratio=1.0
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
- **sim_index**: 4

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
- 2016-10 stats — psi mean: -89.266 m, sw mean: 0.301, ETa mean: 1.112e-08 m/s, recharge mean: n/a
- 2016-11: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2016-11 stats — psi mean: -83.111 m, sw mean: 0.311, ETa mean: 4.388e-09 m/s, recharge mean: n/a
- 2016-12: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2016-12 stats — psi mean: -46.732 m, sw mean: 0.376, ETa mean: 6.432e-09 m/s, recharge mean: n/a
- 2017-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-01 stats — psi mean: -40.920 m, sw mean: 0.359, ETa mean: 1.419e-08 m/s, recharge mean: n/a
- 2017-02: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2017-02 stats — psi mean: -42.852 m, sw mean: 0.343, ETa mean: 1.523e-08 m/s, recharge mean: n/a
- 2017-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-03 stats — psi mean: -47.278 m, sw mean: 0.347, ETa mean: 2.280e-08 m/s, recharge mean: n/a
- 2017-04: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-04 stats — psi mean: -58.860 m, sw mean: 0.332, ETa mean: 1.210e-08 m/s, recharge mean: n/a
- 2017-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-05 stats — psi mean: -58.060 m, sw mean: 0.330, ETa mean: 2.202e-08 m/s, recharge mean: n/a
- 2017-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-06 stats — psi mean: -66.914 m, sw mean: 0.319, ETa mean: 1.448e-08 m/s, recharge mean: n/a
- 2017-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-07 stats — psi mean: -69.026 m, sw mean: 0.317, ETa mean: 1.124e-08 m/s, recharge mean: n/a
- 2017-08: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-08 stats — psi mean: -65.839 m, sw mean: 0.319, ETa mean: 9.120e-09 m/s, recharge mean: n/a
- 2017-09: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-09 stats — psi mean: -52.609 m, sw mean: 0.331, ETa mean: 2.066e-08 m/s, recharge mean: n/a
- 2017-10: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2017-10 stats — psi mean: -63.048 m, sw mean: 0.317, ETa mean: 7.839e-09 m/s, recharge mean: n/a
- 2017-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-11 stats — psi mean: -64.863 m, sw mean: 0.316, ETa mean: 3.306e-09 m/s, recharge mean: n/a
- 2017-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-12 stats — psi mean: -45.960 m, sw mean: 0.325, ETa mean: 9.790e-09 m/s, recharge mean: n/a
- 2018-01: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-01 stats — psi mean: -54.369 m, sw mean: 0.321, ETa mean: 5.779e-09 m/s, recharge mean: n/a
- 2018-02: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 4840 cell(s) set for net < 0 (mode=etp).
- 2018-02 stats — psi mean: -38.135 m, sw mean: 0.341, ETa mean: 1.248e-08 m/s, recharge mean: n/a
- 2018-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2018-03 stats — psi mean: -43.893 m, sw mean: 0.331, ETa mean: 1.454e-08 m/s, recharge mean: n/a
- 2018-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2018-04 stats — psi mean: -60.617 m, sw mean: 0.317, ETa mean: 9.880e-09 m/s, recharge mean: n/a
- 2018-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-05 stats — psi mean: -64.199 m, sw mean: 0.315, ETa mean: 9.858e-09 m/s, recharge mean: n/a
- 2018-06: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-06 stats — psi mean: -57.311 m, sw mean: 0.320, ETa mean: 1.701e-08 m/s, recharge mean: n/a
- 2018-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2018-07 stats — psi mean: -68.469 m, sw mean: 0.312, ETa mean: 9.708e-09 m/s, recharge mean: n/a
- 2018-08: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-08 stats — psi mean: -66.846 m, sw mean: 0.312, ETa mean: 7.691e-09 m/s, recharge mean: n/a
- 2018-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-09 stats — psi mean: -52.827 m, sw mean: 0.320, ETa mean: 1.462e-08 m/s, recharge mean: n/a
- 2018-10: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-10 stats — psi mean: -57.251 m, sw mean: 0.316, ETa mean: 8.305e-09 m/s, recharge mean: n/a
- 2018-11: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2018-11 stats — psi mean: -41.817 m, sw mean: 0.339, ETa mean: 9.231e-09 m/s, recharge mean: n/a
- 2018-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-12 stats — psi mean: -45.120 m, sw mean: 0.328, ETa mean: 1.198e-08 m/s, recharge mean: n/a
- 2019-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2019-01 stats — psi mean: -61.381 m, sw mean: 0.317, ETa mean: 6.749e-09 m/s, recharge mean: n/a
- 2019-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2019-02 stats — psi mean: -64.811 m, sw mean: 0.314, ETa mean: 5.447e-09 m/s, recharge mean: n/a
- 2019-03: ETa artefact floor patched 5 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-03 stats — psi mean: -66.856 m, sw mean: 0.313, ETa mean: 5.886e-09 m/s, recharge mean: n/a
- 2019-04: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2019-04 stats — psi mean: -52.234 m, sw mean: 0.337, ETa mean: 1.542e-08 m/s, recharge mean: n/a
- 2019-05: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-05 stats — psi mean: -54.933 m, sw mean: 0.330, ETa mean: 2.024e-08 m/s, recharge mean: n/a
- 2019-06: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-06 stats — psi mean: -64.371 m, sw mean: 0.320, ETa mean: 1.215e-08 m/s, recharge mean: n/a
- 2019-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2019-07 stats — psi mean: -66.083 m, sw mean: 0.318, ETa mean: 1.045e-08 m/s, recharge mean: n/a
- 2019-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-08 stats — psi mean: -65.521 m, sw mean: 0.317, ETa mean: 9.080e-09 m/s, recharge mean: n/a
- 2019-09: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2019-09 stats — psi mean: -47.913 m, sw mean: 0.347, ETa mean: 1.790e-08 m/s, recharge mean: n/a
- 2019-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2019-10 stats — psi mean: -48.549 m, sw mean: 0.340, ETa mean: 1.554e-08 m/s, recharge mean: n/a
- 2019-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2019-11 stats — psi mean: -47.433 m, sw mean: 0.332, ETa mean: 1.192e-08 m/s, recharge mean: n/a
- 2019-12: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2019-12 stats — psi mean: -34.266 m, sw mean: 0.358, ETa mean: 1.598e-08 m/s, recharge mean: n/a
- 2020-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2020-01 stats — psi mean: -40.825 m, sw mean: 0.352, ETa mean: 1.062e-08 m/s, recharge mean: n/a
- 2020-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2020-02 stats — psi mean: -44.298 m, sw mean: 0.342, ETa mean: 1.616e-08 m/s, recharge mean: n/a
- 2020-03: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2020-03 stats — psi mean: -45.740 m, sw mean: 0.349, ETa mean: 9.445e-09 m/s, recharge mean: n/a
- 2020-04: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2020-04 stats — psi mean: -36.595 m, sw mean: 0.356, ETa mean: 2.387e-08 m/s, recharge mean: n/a
- 2020-05: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2020-05 stats — psi mean: -50.537 m, sw mean: 0.337, ETa mean: 2.253e-08 m/s, recharge mean: n/a
- 2020-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2020-06 stats — psi mean: -61.881 m, sw mean: 0.326, ETa mean: 1.310e-08 m/s, recharge mean: n/a
- 2020-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2020-07 stats — psi mean: -64.264 m, sw mean: 0.323, ETa mean: 1.231e-08 m/s, recharge mean: n/a
- 2020-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2020-08 stats — psi mean: -64.789 m, sw mean: 0.321, ETa mean: 1.098e-08 m/s, recharge mean: n/a
- 2020-09: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2020-09 stats — psi mean: -64.676 m, sw mean: 0.320, ETa mean: 7.568e-09 m/s, recharge mean: n/a
- 2020-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2020-10 stats — psi mean: -62.473 m, sw mean: 0.320, ETa mean: 6.721e-09 m/s, recharge mean: n/a
- 2020-11: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2020-11 stats — psi mean: -42.590 m, sw mean: 0.337, ETa mean: 1.265e-08 m/s, recharge mean: n/a
- 2020-12: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2020-12 stats — psi mean: -37.446 m, sw mean: 0.336, ETa mean: 1.135e-08 m/s, recharge mean: n/a
- 2021-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6296 cell(s) set for net < 0 (mode=etp).
- 2021-01 stats — psi mean: -38.779 m, sw mean: 0.349, ETa mean: 1.506e-08 m/s, recharge mean: n/a
- 2021-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2021-02 stats — psi mean: -50.738 m, sw mean: 0.330, ETa mean: 8.778e-09 m/s, recharge mean: n/a
- 2021-03: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2021-03 stats — psi mean: -45.513 m, sw mean: 0.334, ETa mean: 1.449e-08 m/s, recharge mean: n/a
- 2021-04: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2021-04 stats — psi mean: -48.594 m, sw mean: 0.332, ETa mean: 9.671e-09 m/s, recharge mean: n/a
- 2021-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2021-05 stats — psi mean: -47.216 m, sw mean: 0.332, ETa mean: 2.549e-08 m/s, recharge mean: n/a
- 2021-06: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2021-06 stats — psi mean: -50.927 m, sw mean: 0.328, ETa mean: 2.364e-08 m/s, recharge mean: n/a
- 2021-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2021-07 stats — psi mean: -63.215 m, sw mean: 0.320, ETa mean: 1.178e-08 m/s, recharge mean: n/a
- 2021-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2021-08 stats — psi mean: -63.975 m, sw mean: 0.319, ETa mean: 9.711e-09 m/s, recharge mean: n/a
- 2021-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2021-09 stats — psi mean: -47.742 m, sw mean: 0.330, ETa mean: 2.000e-08 m/s, recharge mean: n/a
- 2021-10: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2021-10 stats — psi mean: -51.773 m, sw mean: 0.323, ETa mean: 8.818e-09 m/s, recharge mean: n/a
- 2021-11: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2021-11 stats — psi mean: -38.226 m, sw mean: 0.334, ETa mean: 1.191e-08 m/s, recharge mean: n/a
- 2021-12: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2021-12 stats — psi mean: -45.900 m, sw mean: 0.324, ETa mean: 7.942e-09 m/s, recharge mean: n/a
- 2022-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-01 stats — psi mean: -55.786 m, sw mean: 0.319, ETa mean: 4.286e-09 m/s, recharge mean: n/a
- 2022-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2022-02 stats — psi mean: -59.889 m, sw mean: 0.318, ETa mean: 4.178e-09 m/s, recharge mean: n/a
- 2022-03: ETa artefact floor patched 4598 cell(s) (net rain ≥ ETp); 3146 cell(s) set for net < 0 (mode=etp).
- 2022-03 stats — psi mean: -30.453 m, sw mean: 0.373, ETa mean: 8.637e-09 m/s, recharge mean: n/a
- 2022-04: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2022-04 stats — psi mean: -28.914 m, sw mean: 0.381, ETa mean: 3.111e-08 m/s, recharge mean: n/a
- 2022-05: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-05 stats — psi mean: -46.280 m, sw mean: 0.346, ETa mean: 2.382e-08 m/s, recharge mean: n/a
- 2022-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2022-06 stats — psi mean: -57.046 m, sw mean: 0.334, ETa mean: 1.823e-08 m/s, recharge mean: n/a
- 2022-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-07 stats — psi mean: -60.360 m, sw mean: 0.329, ETa mean: 1.257e-08 m/s, recharge mean: n/a
- 2022-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-08 stats — psi mean: -61.142 m, sw mean: 0.328, ETa mean: 1.126e-08 m/s, recharge mean: n/a
- 2022-09: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2022-09 stats — psi mean: -48.813 m, sw mean: 0.334, ETa mean: 1.612e-08 m/s, recharge mean: n/a
- 2022-10: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-10 stats — psi mean: -44.954 m, sw mean: 0.333, ETa mean: 1.395e-08 m/s, recharge mean: n/a
- 2022-11: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2022-11 stats — psi mean: -56.248 m, sw mean: 0.326, ETa mean: 5.876e-09 m/s, recharge mean: n/a
- 2022-12: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-12 stats — psi mean: -56.777 m, sw mean: 0.325, ETa mean: 3.759e-09 m/s, recharge mean: n/a
- 2023-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-01 stats — psi mean: -59.475 m, sw mean: 0.324, ETa mean: 3.423e-09 m/s, recharge mean: n/a
- 2023-02: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-02 stats — psi mean: -51.896 m, sw mean: 0.326, ETa mean: 4.705e-09 m/s, recharge mean: n/a
- 2023-03: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-03 stats — psi mean: -60.549 m, sw mean: 0.323, ETa mean: 6.367e-09 m/s, recharge mean: n/a
- 2023-04: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-04 stats — psi mean: -61.300 m, sw mean: 0.322, ETa mean: 8.387e-09 m/s, recharge mean: n/a
- 2023-05: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2023-05 stats — psi mean: -47.942 m, sw mean: 0.341, ETa mean: 9.370e-09 m/s, recharge mean: n/a
- 2023-06: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-06 stats — psi mean: -36.690 m, sw mean: 0.353, ETa mean: 2.341e-08 m/s, recharge mean: n/a
- 2023-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-07 stats — psi mean: -55.943 m, sw mean: 0.329, ETa mean: 1.094e-08 m/s, recharge mean: n/a
- 2023-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-08 stats — psi mean: -58.437 m, sw mean: 0.326, ETa mean: 8.703e-09 m/s, recharge mean: n/a
- 2023-09: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-09 stats — psi mean: -38.377 m, sw mean: 0.339, ETa mean: 1.868e-08 m/s, recharge mean: n/a
- 2023-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2023-10 stats — psi mean: -48.968 m, sw mean: 0.329, ETa mean: 9.160e-09 m/s, recharge mean: n/a
- 2023-11: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-11 stats — psi mean: -57.126 m, sw mean: 0.324, ETa mean: 4.780e-09 m/s, recharge mean: n/a
- 2023-12: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-12 stats — psi mean: -57.993 m, sw mean: 0.323, ETa mean: 3.787e-09 m/s, recharge mean: n/a
- 2024-01: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2024-01 stats — psi mean: -48.112 m, sw mean: 0.326, ETa mean: 6.348e-09 m/s, recharge mean: n/a
- 2024-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-02 stats — psi mean: -55.760 m, sw mean: 0.323, ETa mean: 6.053e-09 m/s, recharge mean: n/a
- 2024-03: ETa artefact floor patched 826 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-03 stats — psi mean: -53.769 m, sw mean: 0.324, ETa mean: 9.295e-09 m/s, recharge mean: n/a
- 2024-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2024-04 stats — psi mean: -57.762 m, sw mean: 0.322, ETa mean: 8.182e-09 m/s, recharge mean: n/a
- 2024-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2024-05 stats — psi mean: -57.544 m, sw mean: 0.322, ETa mean: 1.156e-08 m/s, recharge mean: n/a
- 2024-06: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2024-06 stats — psi mean: -50.427 m, sw mean: 0.327, ETa mean: 1.980e-08 m/s, recharge mean: n/a
- 2024-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2024-07 stats — psi mean: -59.595 m, sw mean: 0.320, ETa mean: 1.033e-08 m/s, recharge mean: n/a
- 2024-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2024-08 stats — psi mean: -60.293 m, sw mean: 0.320, ETa mean: 7.534e-09 m/s, recharge mean: n/a
- 2024-09: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-09 stats — psi mean: -55.838 m, sw mean: 0.321, ETa mean: 8.050e-09 m/s, recharge mean: n/a
- 2024-10: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2024-10 stats — psi mean: -52.894 m, sw mean: 0.327, ETa mean: 5.482e-09 m/s, recharge mean: n/a
- 2024-11: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2024-11 stats — psi mean: -25.770 m, sw mean: 0.362, ETa mean: 1.528e-08 m/s, recharge mean: n/a
- 2024-12: ETa artefact floor patched 1 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2024-12 stats — psi mean: -47.681 m, sw mean: 0.330, ETa mean: 7.680e-09 m/s, recharge mean: n/a
- 2025-01: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2025-01 stats — psi mean: -51.818 m, sw mean: 0.326, ETa mean: 3.918e-09 m/s, recharge mean: n/a
- 2025-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2025-02 stats — psi mean: -54.885 m, sw mean: 0.324, ETa mean: 4.278e-09 m/s, recharge mean: n/a
- 2025-03: ETa artefact floor patched 3872 cell(s) (net rain ≥ ETp); 3872 cell(s) set for net < 0 (mode=etp).
- 2025-03 stats — psi mean: -26.971 m, sw mean: 0.373, ETa mean: 1.705e-08 m/s, recharge mean: n/a
- 2025-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2025-04 stats — psi mean: -45.937 m, sw mean: 0.336, ETa mean: 1.397e-08 m/s, recharge mean: n/a
- 2025-05: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2025-05 stats — psi mean: -45.795 m, sw mean: 0.333, ETa mean: 1.879e-08 m/s, recharge mean: n/a
- 2025-06: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2025-06 stats — psi mean: -56.660 m, sw mean: 0.327, ETa mean: 9.769e-09 m/s, recharge mean: n/a
- 2025-07: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2025-07 stats — psi mean: -51.752 m, sw mean: 0.330, ETa mean: 1.427e-08 m/s, recharge mean: n/a
- 2025-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2025-08 stats — psi mean: -55.054 m, sw mean: 0.325, ETa mean: 9.606e-09 m/s, recharge mean: n/a
| 2025-09 | — | — | — | — | — | — | skipped (no forcing) |

## Run summary

- **Months run**: 107
- **Months skipped**: 12
- **Output directory**: /home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/Agramon_SShydro/hypotheses/hyp1_layers_Ks/prj/outputs/scenario_4_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2
- **ETa artefact cells patched (total)**: 90130
- **ETa net<0 cells set (total)**: 723100
- **Warnings logged**: 0
- **Elapsed**: 0:44:34

---

_Finished: 2026-11-20 13:38:26_
