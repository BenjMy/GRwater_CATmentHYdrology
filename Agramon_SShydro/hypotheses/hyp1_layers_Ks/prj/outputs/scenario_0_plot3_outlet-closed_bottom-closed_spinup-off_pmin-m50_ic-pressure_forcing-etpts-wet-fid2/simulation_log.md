# Simulation Log — scenario_0_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2 (dem-plot 3)

_Started: 2026-11-20 11:25:47_


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
- **DEM source (tif)**: /home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/DTMplots/20250618_AGRAMON100m_micasense_dtm.tif
- **DEM plot (fid_1 clip)**: 3
- **Spin-up**: disabled
- **Outlet flux**: 0.0e+00 m/s
- **Outlet side**: auto-detect
- **Bottom flux**: disabled (closed)
- **Solver dtmin/dtmax/deltat**: 0.1 / 10000.0 / 1000.0 s
- **CSV run registry**: /home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/Agramon_SShydro/hypotheses/hyp1_layers_Ks/simulation_log.csv
- **sim_index**: 0

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
- 2016-10 stats — psi mean: -91.017 m, sw mean: 0.304, ETa mean: 1.122e-08 m/s, recharge mean: n/a
- 2016-11: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2016-11 stats — psi mean: -87.729 m, sw mean: 0.314, ETa mean: 4.900e-09 m/s, recharge mean: n/a
- 2016-12: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2016-12 stats — psi mean: -52.527 m, sw mean: 0.397, ETa mean: 6.426e-09 m/s, recharge mean: n/a
- 2017-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-01 stats — psi mean: -45.646 m, sw mean: 0.378, ETa mean: 1.465e-08 m/s, recharge mean: n/a
- 2017-02: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2017-02 stats — psi mean: -46.389 m, sw mean: 0.358, ETa mean: 1.621e-08 m/s, recharge mean: n/a
- 2017-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-03 stats — psi mean: -51.275 m, sw mean: 0.361, ETa mean: 2.428e-08 m/s, recharge mean: n/a
- 2017-04: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-04 stats — psi mean: -63.036 m, sw mean: 0.339, ETa mean: 1.382e-08 m/s, recharge mean: n/a
- 2017-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-05 stats — psi mean: -63.485 m, sw mean: 0.335, ETa mean: 2.274e-08 m/s, recharge mean: n/a
- 2017-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-06 stats — psi mean: -72.759 m, sw mean: 0.322, ETa mean: 1.420e-08 m/s, recharge mean: n/a
- 2017-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-07 stats — psi mean: -74.932 m, sw mean: 0.319, ETa mean: 1.113e-08 m/s, recharge mean: n/a
- 2017-08: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-08 stats — psi mean: -72.031 m, sw mean: 0.322, ETa mean: 8.884e-09 m/s, recharge mean: n/a
- 2017-09: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-09 stats — psi mean: -58.531 m, sw mean: 0.336, ETa mean: 2.166e-08 m/s, recharge mean: n/a
- 2017-10: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2017-10 stats — psi mean: -70.056 m, sw mean: 0.319, ETa mean: 7.308e-09 m/s, recharge mean: n/a
- 2017-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-11 stats — psi mean: -71.792 m, sw mean: 0.318, ETa mean: 3.085e-09 m/s, recharge mean: n/a
- 2017-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-12 stats — psi mean: -53.328 m, sw mean: 0.329, ETa mean: 9.679e-09 m/s, recharge mean: n/a
- 2018-01: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-01 stats — psi mean: -62.346 m, sw mean: 0.324, ETa mean: 5.350e-09 m/s, recharge mean: n/a
- 2018-02: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 4840 cell(s) set for net < 0 (mode=etp).
- 2018-02 stats — psi mean: -43.207 m, sw mean: 0.349, ETa mean: 1.312e-08 m/s, recharge mean: n/a
- 2018-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2018-03 stats — psi mean: -50.313 m, sw mean: 0.334, ETa mean: 1.452e-08 m/s, recharge mean: n/a
- 2018-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2018-04 stats — psi mean: -67.734 m, sw mean: 0.319, ETa mean: 9.413e-09 m/s, recharge mean: n/a
- 2018-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-05 stats — psi mean: -70.960 m, sw mean: 0.317, ETa mean: 9.352e-09 m/s, recharge mean: n/a
- 2018-06: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-06 stats — psi mean: -64.229 m, sw mean: 0.323, ETa mean: 1.692e-08 m/s, recharge mean: n/a
- 2018-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2018-07 stats — psi mean: -74.913 m, sw mean: 0.314, ETa mean: 9.440e-09 m/s, recharge mean: n/a
- 2018-08: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-08 stats — psi mean: -73.382 m, sw mean: 0.314, ETa mean: 7.358e-09 m/s, recharge mean: n/a
- 2018-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-09 stats — psi mean: -60.930 m, sw mean: 0.324, ETa mean: 1.421e-08 m/s, recharge mean: n/a
- 2018-10: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-10 stats — psi mean: -65.070 m, sw mean: 0.319, ETa mean: 7.734e-09 m/s, recharge mean: n/a
- 2018-11: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2018-11 stats — psi mean: -49.622 m, sw mean: 0.346, ETa mean: 9.052e-09 m/s, recharge mean: n/a
- 2018-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-12 stats — psi mean: -49.793 m, sw mean: 0.334, ETa mean: 1.291e-08 m/s, recharge mean: n/a
- 2019-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2019-01 stats — psi mean: -67.174 m, sw mean: 0.318, ETa mean: 6.926e-09 m/s, recharge mean: n/a
- 2019-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2019-02 stats — psi mean: -71.341 m, sw mean: 0.315, ETa mean: 5.163e-09 m/s, recharge mean: n/a
- 2019-03: ETa artefact floor patched 12 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-03 stats — psi mean: -73.447 m, sw mean: 0.314, ETa mean: 5.665e-09 m/s, recharge mean: n/a
- 2019-04: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2019-04 stats — psi mean: -59.086 m, sw mean: 0.344, ETa mean: 1.519e-08 m/s, recharge mean: n/a
- 2019-05: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-05 stats — psi mean: -59.734 m, sw mean: 0.334, ETa mean: 2.269e-08 m/s, recharge mean: n/a
- 2019-06: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-06 stats — psi mean: -70.737 m, sw mean: 0.321, ETa mean: 1.211e-08 m/s, recharge mean: n/a
- 2019-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2019-07 stats — psi mean: -72.676 m, sw mean: 0.318, ETa mean: 1.026e-08 m/s, recharge mean: n/a
- 2019-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-08 stats — psi mean: -72.228 m, sw mean: 0.318, ETa mean: 8.749e-09 m/s, recharge mean: n/a
- 2019-09: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2019-09 stats — psi mean: -54.064 m, sw mean: 0.356, ETa mean: 1.897e-08 m/s, recharge mean: n/a
- 2019-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2019-10 stats — psi mean: -54.056 m, sw mean: 0.345, ETa mean: 1.683e-08 m/s, recharge mean: n/a
- 2019-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2019-11 stats — psi mean: -54.268 m, sw mean: 0.335, ETa mean: 1.225e-08 m/s, recharge mean: n/a
- 2019-12: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2019-12 stats — psi mean: -40.475 m, sw mean: 0.368, ETa mean: 1.622e-08 m/s, recharge mean: n/a
- 2020-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2020-01 stats — psi mean: -46.463 m, sw mean: 0.360, ETa mean: 1.131e-08 m/s, recharge mean: n/a
- 2020-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2020-02 stats — psi mean: -49.405 m, sw mean: 0.349, ETa mean: 1.769e-08 m/s, recharge mean: n/a
- 2020-03: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2020-03 stats — psi mean: -52.864 m, sw mean: 0.355, ETa mean: 9.677e-09 m/s, recharge mean: n/a
- 2020-04: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2020-04 stats — psi mean: -42.412 m, sw mean: 0.366, ETa mean: 2.520e-08 m/s, recharge mean: n/a
- 2020-05: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2020-05 stats — psi mean: -57.461 m, sw mean: 0.340, ETa mean: 2.342e-08 m/s, recharge mean: n/a
- 2020-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2020-06 stats — psi mean: -69.268 m, sw mean: 0.326, ETa mean: 1.288e-08 m/s, recharge mean: n/a
- 2020-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2020-07 stats — psi mean: -71.787 m, sw mean: 0.323, ETa mean: 1.213e-08 m/s, recharge mean: n/a
- 2020-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2020-08 stats — psi mean: -72.374 m, sw mean: 0.321, ETa mean: 1.082e-08 m/s, recharge mean: n/a
- 2020-09: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2020-09 stats — psi mean: -72.375 m, sw mean: 0.320, ETa mean: 7.424e-09 m/s, recharge mean: n/a
- 2020-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2020-10 stats — psi mean: -70.374 m, sw mean: 0.321, ETa mean: 6.364e-09 m/s, recharge mean: n/a
- 2020-11: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2020-11 stats — psi mean: -50.176 m, sw mean: 0.340, ETa mean: 1.325e-08 m/s, recharge mean: n/a
- 2020-12: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2020-12 stats — psi mean: -45.665 m, sw mean: 0.338, ETa mean: 1.151e-08 m/s, recharge mean: n/a
- 2021-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6298 cell(s) set for net < 0 (mode=etp).
- 2021-01 stats — psi mean: -46.513 m, sw mean: 0.354, ETa mean: 1.565e-08 m/s, recharge mean: n/a
- 2021-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2021-02 stats — psi mean: -58.015 m, sw mean: 0.330, ETa mean: 9.181e-09 m/s, recharge mean: n/a
- 2021-03: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2021-03 stats — psi mean: -53.959 m, sw mean: 0.334, ETa mean: 1.445e-08 m/s, recharge mean: n/a
- 2021-04: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2021-04 stats — psi mean: -57.780 m, sw mean: 0.333, ETa mean: 9.206e-09 m/s, recharge mean: n/a
- 2021-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2021-05 stats — psi mean: -55.573 m, sw mean: 0.333, ETa mean: 2.622e-08 m/s, recharge mean: n/a
- 2021-06: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2021-06 stats — psi mean: -60.914 m, sw mean: 0.327, ETa mean: 2.270e-08 m/s, recharge mean: n/a
- 2021-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2021-07 stats — psi mean: -71.872 m, sw mean: 0.318, ETa mean: 1.143e-08 m/s, recharge mean: n/a
- 2021-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2021-08 stats — psi mean: -72.577 m, sw mean: 0.318, ETa mean: 9.457e-09 m/s, recharge mean: n/a
- 2021-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2021-09 stats — psi mean: -57.388 m, sw mean: 0.331, ETa mean: 2.025e-08 m/s, recharge mean: n/a
- 2021-10: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2021-10 stats — psi mean: -61.867 m, sw mean: 0.322, ETa mean: 8.334e-09 m/s, recharge mean: n/a
- 2021-11: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2021-11 stats — psi mean: -48.842 m, sw mean: 0.335, ETa mean: 1.203e-08 m/s, recharge mean: n/a
- 2021-12: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2021-12 stats — psi mean: -56.786 m, sw mean: 0.324, ETa mean: 7.724e-09 m/s, recharge mean: n/a
- 2022-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-01 stats — psi mean: -66.142 m, sw mean: 0.318, ETa mean: 3.965e-09 m/s, recharge mean: n/a
- 2022-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2022-02 stats — psi mean: -69.357 m, sw mean: 0.317, ETa mean: 3.898e-09 m/s, recharge mean: n/a
- 2022-03: ETa artefact floor patched 4598 cell(s) (net rain ≥ ETp); 3146 cell(s) set for net < 0 (mode=etp).
- 2022-03 stats — psi mean: -39.987 m, sw mean: 0.385, ETa mean: 8.435e-09 m/s, recharge mean: n/a
- 2022-04: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2022-04 stats — psi mean: -36.485 m, sw mean: 0.401, ETa mean: 3.193e-08 m/s, recharge mean: n/a
- 2022-05: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-05 stats — psi mean: -53.607 m, sw mean: 0.353, ETa mean: 2.542e-08 m/s, recharge mean: n/a
- 2022-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2022-06 stats — psi mean: -65.188 m, sw mean: 0.337, ETa mean: 1.819e-08 m/s, recharge mean: n/a
- 2022-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-07 stats — psi mean: -68.515 m, sw mean: 0.332, ETa mean: 1.255e-08 m/s, recharge mean: n/a
- 2022-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-08 stats — psi mean: -69.378 m, sw mean: 0.330, ETa mean: 1.120e-08 m/s, recharge mean: n/a
- 2022-09: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2022-09 stats — psi mean: -58.253 m, sw mean: 0.338, ETa mean: 1.608e-08 m/s, recharge mean: n/a
- 2022-10: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-10 stats — psi mean: -55.180 m, sw mean: 0.336, ETa mean: 1.363e-08 m/s, recharge mean: n/a
- 2022-11: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2022-11 stats — psi mean: -65.011 m, sw mean: 0.328, ETa mean: 5.684e-09 m/s, recharge mean: n/a
- 2022-12: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-12 stats — psi mean: -65.399 m, sw mean: 0.327, ETa mean: 3.551e-09 m/s, recharge mean: n/a
- 2023-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-01 stats — psi mean: -67.940 m, sw mean: 0.326, ETa mean: 3.349e-09 m/s, recharge mean: n/a
- 2023-02: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-02 stats — psi mean: -61.354 m, sw mean: 0.329, ETa mean: 4.511e-09 m/s, recharge mean: n/a
- 2023-03: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-03 stats — psi mean: -68.997 m, sw mean: 0.324, ETa mean: 6.273e-09 m/s, recharge mean: n/a
- 2023-04: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-04 stats — psi mean: -69.723 m, sw mean: 0.324, ETa mean: 8.299e-09 m/s, recharge mean: n/a
- 2023-05: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2023-05 stats — psi mean: -57.190 m, sw mean: 0.347, ETa mean: 9.203e-09 m/s, recharge mean: n/a
- 2023-06: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-06 stats — psi mean: -44.027 m, sw mean: 0.363, ETa mean: 2.543e-08 m/s, recharge mean: n/a
- 2023-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-07 stats — psi mean: -64.296 m, sw mean: 0.330, ETa mean: 1.112e-08 m/s, recharge mean: n/a
- 2023-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-08 stats — psi mean: -67.498 m, sw mean: 0.326, ETa mean: 8.470e-09 m/s, recharge mean: n/a
- 2023-09: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-09 stats — psi mean: -47.825 m, sw mean: 0.342, ETa mean: 1.880e-08 m/s, recharge mean: n/a
- 2023-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2023-10 stats — psi mean: -59.312 m, sw mean: 0.330, ETa mean: 8.707e-09 m/s, recharge mean: n/a
- 2023-11: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-11 stats — psi mean: -66.534 m, sw mean: 0.324, ETa mean: 4.551e-09 m/s, recharge mean: n/a
- 2023-12: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-12 stats — psi mean: -67.273 m, sw mean: 0.324, ETa mean: 3.572e-09 m/s, recharge mean: n/a
- 2024-01: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2024-01 stats — psi mean: -58.094 m, sw mean: 0.328, ETa mean: 5.884e-09 m/s, recharge mean: n/a
- 2024-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-02 stats — psi mean: -64.992 m, sw mean: 0.324, ETa mean: 5.670e-09 m/s, recharge mean: n/a
- 2024-03: ETa artefact floor patched 825 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-03 stats — psi mean: -63.131 m, sw mean: 0.325, ETa mean: 8.881e-09 m/s, recharge mean: n/a
- 2024-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2024-04 stats — psi mean: -66.717 m, sw mean: 0.323, ETa mean: 7.930e-09 m/s, recharge mean: n/a
- 2024-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2024-05 stats — psi mean: -66.608 m, sw mean: 0.323, ETa mean: 1.127e-08 m/s, recharge mean: n/a
- 2024-06: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2024-06 stats — psi mean: -60.516 m, sw mean: 0.328, ETa mean: 1.947e-08 m/s, recharge mean: n/a
- 2024-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2024-07 stats — psi mean: -68.622 m, sw mean: 0.321, ETa mean: 1.025e-08 m/s, recharge mean: n/a
- 2024-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2024-08 stats — psi mean: -69.318 m, sw mean: 0.321, ETa mean: 7.312e-09 m/s, recharge mean: n/a
- 2024-09: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-09 stats — psi mean: -65.253 m, sw mean: 0.323, ETa mean: 7.673e-09 m/s, recharge mean: n/a
- 2024-10: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2024-10 stats — psi mean: -62.409 m, sw mean: 0.329, ETa mean: 5.279e-09 m/s, recharge mean: n/a
- 2024-11: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2024-11 stats — psi mean: -33.872 m, sw mean: 0.374, ETa mean: 1.601e-08 m/s, recharge mean: n/a
- 2024-12: ETa artefact floor patched 2 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2024-12 stats — psi mean: -55.093 m, sw mean: 0.332, ETa mean: 8.698e-09 m/s, recharge mean: n/a
- 2025-01: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2025-01 stats — psi mean: -60.660 m, sw mean: 0.326, ETa mean: 3.798e-09 m/s, recharge mean: n/a
- 2025-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2025-02 stats — psi mean: -64.014 m, sw mean: 0.324, ETa mean: 4.109e-09 m/s, recharge mean: n/a
- 2025-03: ETa artefact floor patched 3872 cell(s) (net rain ≥ ETp); 3872 cell(s) set for net < 0 (mode=etp).
- 2025-03 stats — psi mean: -35.991 m, sw mean: 0.385, ETa mean: 1.798e-08 m/s, recharge mean: n/a
- 2025-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2025-04 stats — psi mean: -53.399 m, sw mean: 0.339, ETa mean: 1.541e-08 m/s, recharge mean: n/a
- 2025-05: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2025-05 stats — psi mean: -55.139 m, sw mean: 0.334, ETa mean: 1.809e-08 m/s, recharge mean: n/a
- 2025-06: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2025-06 stats — psi mean: -65.637 m, sw mean: 0.326, ETa mean: 9.719e-09 m/s, recharge mean: n/a
- 2025-07: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2025-07 stats — psi mean: -61.089 m, sw mean: 0.330, ETa mean: 1.459e-08 m/s, recharge mean: n/a
- 2025-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2025-08 stats — psi mean: -64.855 m, sw mean: 0.324, ETa mean: 9.254e-09 m/s, recharge mean: n/a
| 2025-09 | — | — | — | — | — | — | skipped (no forcing) |

## Run summary

- **Months run**: 107
- **Months skipped**: 12
- **Output directory**: /home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/Agramon_SShydro/hypotheses/hyp1_layers_Ks/prj/outputs/scenario_0_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2
- **ETa artefact cells patched (total)**: 90137
- **ETa net<0 cells set (total)**: 723102
- **Warnings logged**: 0
- **Elapsed**: 0:43:29

---

_Finished: 2026-11-20 12:09:17_
