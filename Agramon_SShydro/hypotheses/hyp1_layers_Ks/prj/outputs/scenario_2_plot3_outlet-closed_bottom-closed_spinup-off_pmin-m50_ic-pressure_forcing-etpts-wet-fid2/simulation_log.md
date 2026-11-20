# Simulation Log — scenario_2_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2 (dem-plot 3)

_Started: 2026-11-20 12:09:26_


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
- **sim_index**: 2

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
- 2016-12 stats — psi mean: -48.038 m, sw mean: 0.381, ETa mean: 6.431e-09 m/s, recharge mean: n/a
- 2017-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-01 stats — psi mean: -41.980 m, sw mean: 0.363, ETa mean: 1.429e-08 m/s, recharge mean: n/a
- 2017-02: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2017-02 stats — psi mean: -43.855 m, sw mean: 0.345, ETa mean: 1.531e-08 m/s, recharge mean: n/a
- 2017-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-03 stats — psi mean: -48.212 m, sw mean: 0.350, ETa mean: 2.308e-08 m/s, recharge mean: n/a
- 2017-04: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-04 stats — psi mean: -59.791 m, sw mean: 0.333, ETa mean: 1.243e-08 m/s, recharge mean: n/a
- 2017-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-05 stats — psi mean: -59.249 m, sw mean: 0.331, ETa mean: 2.216e-08 m/s, recharge mean: n/a
- 2017-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-06 stats — psi mean: -68.255 m, sw mean: 0.320, ETa mean: 1.437e-08 m/s, recharge mean: n/a
- 2017-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-07 stats — psi mean: -70.312 m, sw mean: 0.317, ETa mean: 1.120e-08 m/s, recharge mean: n/a
- 2017-08: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-08 stats — psi mean: -67.152 m, sw mean: 0.320, ETa mean: 9.071e-09 m/s, recharge mean: n/a
- 2017-09: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-09 stats — psi mean: -53.852 m, sw mean: 0.332, ETa mean: 2.087e-08 m/s, recharge mean: n/a
- 2017-10: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2017-10 stats — psi mean: -64.585 m, sw mean: 0.318, ETa mean: 7.720e-09 m/s, recharge mean: n/a
- 2017-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-11 stats — psi mean: -66.425 m, sw mean: 0.317, ETa mean: 3.240e-09 m/s, recharge mean: n/a
- 2017-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-12 stats — psi mean: -47.531 m, sw mean: 0.326, ETa mean: 9.770e-09 m/s, recharge mean: n/a
- 2018-01: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-01 stats — psi mean: -56.214 m, sw mean: 0.321, ETa mean: 5.658e-09 m/s, recharge mean: n/a
- 2018-02: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 4840 cell(s) set for net < 0 (mode=etp).
- 2018-02 stats — psi mean: -39.197 m, sw mean: 0.343, ETa mean: 1.265e-08 m/s, recharge mean: n/a
- 2018-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2018-03 stats — psi mean: -45.283 m, sw mean: 0.332, ETa mean: 1.457e-08 m/s, recharge mean: n/a
- 2018-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2018-04 stats — psi mean: -62.095 m, sw mean: 0.318, ETa mean: 9.807e-09 m/s, recharge mean: n/a
- 2018-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-05 stats — psi mean: -65.678 m, sw mean: 0.315, ETa mean: 9.770e-09 m/s, recharge mean: n/a
- 2018-06: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-06 stats — psi mean: -58.806 m, sw mean: 0.321, ETa mean: 1.694e-08 m/s, recharge mean: n/a
- 2018-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2018-07 stats — psi mean: -69.915 m, sw mean: 0.313, ETa mean: 9.634e-09 m/s, recharge mean: n/a
- 2018-08: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-08 stats — psi mean: -68.306 m, sw mean: 0.313, ETa mean: 7.601e-09 m/s, recharge mean: n/a
- 2018-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-09 stats — psi mean: -54.596 m, sw mean: 0.321, ETa mean: 1.458e-08 m/s, recharge mean: n/a
- 2018-10: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-10 stats — psi mean: -59.078 m, sw mean: 0.317, ETa mean: 8.170e-09 m/s, recharge mean: n/a
- 2018-11: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2018-11 stats — psi mean: -43.635 m, sw mean: 0.340, ETa mean: 9.169e-09 m/s, recharge mean: n/a
- 2018-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-12 stats — psi mean: -46.241 m, sw mean: 0.330, ETa mean: 1.213e-08 m/s, recharge mean: n/a
- 2019-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2019-01 stats — psi mean: -62.671 m, sw mean: 0.317, ETa mean: 6.768e-09 m/s, recharge mean: n/a
- 2019-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2019-02 stats — psi mean: -66.250 m, sw mean: 0.315, ETa mean: 5.379e-09 m/s, recharge mean: n/a
- 2019-03: ETa artefact floor patched 7 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-03 stats — psi mean: -68.297 m, sw mean: 0.313, ETa mean: 5.836e-09 m/s, recharge mean: n/a
- 2019-04: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2019-04 stats — psi mean: -53.725 m, sw mean: 0.339, ETa mean: 1.537e-08 m/s, recharge mean: n/a
- 2019-05: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-05 stats — psi mean: -56.035 m, sw mean: 0.331, ETa mean: 2.071e-08 m/s, recharge mean: n/a
- 2019-06: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-06 stats — psi mean: -65.805 m, sw mean: 0.321, ETa mean: 1.213e-08 m/s, recharge mean: n/a
- 2019-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2019-07 stats — psi mean: -67.530 m, sw mean: 0.318, ETa mean: 1.041e-08 m/s, recharge mean: n/a
- 2019-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-08 stats — psi mean: -66.976 m, sw mean: 0.317, ETa mean: 9.013e-09 m/s, recharge mean: n/a
- 2019-09: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2019-09 stats — psi mean: -49.276 m, sw mean: 0.349, ETa mean: 1.813e-08 m/s, recharge mean: n/a
- 2019-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2019-10 stats — psi mean: -49.798 m, sw mean: 0.341, ETa mean: 1.579e-08 m/s, recharge mean: n/a
- 2019-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2019-11 stats — psi mean: -48.907 m, sw mean: 0.333, ETa mean: 1.203e-08 m/s, recharge mean: n/a
- 2019-12: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2019-12 stats — psi mean: -35.634 m, sw mean: 0.360, ETa mean: 1.605e-08 m/s, recharge mean: n/a
- 2020-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2020-01 stats — psi mean: -42.123 m, sw mean: 0.354, ETa mean: 1.075e-08 m/s, recharge mean: n/a
- 2020-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2020-02 stats — psi mean: -45.483 m, sw mean: 0.344, ETa mean: 1.646e-08 m/s, recharge mean: n/a
- 2020-03: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2020-03 stats — psi mean: -47.321 m, sw mean: 0.351, ETa mean: 9.490e-09 m/s, recharge mean: n/a
- 2020-04: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2020-04 stats — psi mean: -37.936 m, sw mean: 0.359, ETa mean: 2.414e-08 m/s, recharge mean: n/a
- 2020-05: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2020-05 stats — psi mean: -52.084 m, sw mean: 0.338, ETa mean: 2.271e-08 m/s, recharge mean: n/a
- 2020-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2020-06 stats — psi mean: -63.545 m, sw mean: 0.326, ETa mean: 1.303e-08 m/s, recharge mean: n/a
- 2020-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2020-07 stats — psi mean: -65.921 m, sw mean: 0.324, ETa mean: 1.226e-08 m/s, recharge mean: n/a
- 2020-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2020-08 stats — psi mean: -66.440 m, sw mean: 0.322, ETa mean: 1.094e-08 m/s, recharge mean: n/a
- 2020-09: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2020-09 stats — psi mean: -66.339 m, sw mean: 0.321, ETa mean: 7.531e-09 m/s, recharge mean: n/a
- 2020-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2020-10 stats — psi mean: -64.167 m, sw mean: 0.321, ETa mean: 6.638e-09 m/s, recharge mean: n/a
- 2020-11: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2020-11 stats — psi mean: -44.145 m, sw mean: 0.338, ETa mean: 1.281e-08 m/s, recharge mean: n/a
- 2020-12: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2020-12 stats — psi mean: -39.179 m, sw mean: 0.337, ETa mean: 1.137e-08 m/s, recharge mean: n/a
- 2021-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6297 cell(s) set for net < 0 (mode=etp).
- 2021-01 stats — psi mean: -40.416 m, sw mean: 0.350, ETa mean: 1.518e-08 m/s, recharge mean: n/a
- 2021-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2021-02 stats — psi mean: -52.292 m, sw mean: 0.330, ETa mean: 8.860e-09 m/s, recharge mean: n/a
- 2021-03: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2021-03 stats — psi mean: -47.300 m, sw mean: 0.334, ETa mean: 1.451e-08 m/s, recharge mean: n/a
- 2021-04: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2021-04 stats — psi mean: -50.527 m, sw mean: 0.333, ETa mean: 9.588e-09 m/s, recharge mean: n/a
- 2021-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2021-05 stats — psi mean: -48.966 m, sw mean: 0.333, ETa mean: 2.558e-08 m/s, recharge mean: n/a
- 2021-06: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2021-06 stats — psi mean: -52.972 m, sw mean: 0.328, ETa mean: 2.355e-08 m/s, recharge mean: n/a
- 2021-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2021-07 stats — psi mean: -65.090 m, sw mean: 0.320, ETa mean: 1.168e-08 m/s, recharge mean: n/a
- 2021-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2021-08 stats — psi mean: -65.819 m, sw mean: 0.319, ETa mean: 9.645e-09 m/s, recharge mean: n/a
- 2021-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2021-09 stats — psi mean: -49.790 m, sw mean: 0.331, ETa mean: 2.004e-08 m/s, recharge mean: n/a
- 2021-10: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2021-10 stats — psi mean: -53.976 m, sw mean: 0.323, ETa mean: 8.684e-09 m/s, recharge mean: n/a
- 2021-11: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2021-11 stats — psi mean: -40.401 m, sw mean: 0.335, ETa mean: 1.197e-08 m/s, recharge mean: n/a
- 2021-12: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2021-12 stats — psi mean: -48.211 m, sw mean: 0.325, ETa mean: 7.855e-09 m/s, recharge mean: n/a
- 2022-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-01 stats — psi mean: -58.101 m, sw mean: 0.319, ETa mean: 4.202e-09 m/s, recharge mean: n/a
- 2022-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2022-02 stats — psi mean: -62.021 m, sw mean: 0.318, ETa mean: 4.088e-09 m/s, recharge mean: n/a
- 2022-03: ETa artefact floor patched 4598 cell(s) (net rain ≥ ETp); 3146 cell(s) set for net < 0 (mode=etp).
- 2022-03 stats — psi mean: -32.559 m, sw mean: 0.375, ETa mean: 8.596e-09 m/s, recharge mean: n/a
- 2022-04: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2022-04 stats — psi mean: -30.642 m, sw mean: 0.385, ETa mean: 3.129e-08 m/s, recharge mean: n/a
- 2022-05: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-05 stats — psi mean: -48.043 m, sw mean: 0.348, ETa mean: 2.395e-08 m/s, recharge mean: n/a
- 2022-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2022-06 stats — psi mean: -58.957 m, sw mean: 0.335, ETa mean: 1.815e-08 m/s, recharge mean: n/a
- 2022-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-07 stats — psi mean: -62.259 m, sw mean: 0.330, ETa mean: 1.253e-08 m/s, recharge mean: n/a
- 2022-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-08 stats — psi mean: -63.051 m, sw mean: 0.328, ETa mean: 1.123e-08 m/s, recharge mean: n/a
- 2022-09: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2022-09 stats — psi mean: -50.893 m, sw mean: 0.335, ETa mean: 1.611e-08 m/s, recharge mean: n/a
- 2022-10: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-10 stats — psi mean: -47.350 m, sw mean: 0.334, ETa mean: 1.389e-08 m/s, recharge mean: n/a
- 2022-11: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2022-11 stats — psi mean: -58.351 m, sw mean: 0.327, ETa mean: 5.802e-09 m/s, recharge mean: n/a
- 2022-12: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-12 stats — psi mean: -58.818 m, sw mean: 0.326, ETa mean: 3.697e-09 m/s, recharge mean: n/a
- 2023-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-01 stats — psi mean: -61.462 m, sw mean: 0.324, ETa mean: 3.391e-09 m/s, recharge mean: n/a
- 2023-02: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-02 stats — psi mean: -54.120 m, sw mean: 0.327, ETa mean: 4.642e-09 m/s, recharge mean: n/a
- 2023-03: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-03 stats — psi mean: -62.479 m, sw mean: 0.323, ETa mean: 6.332e-09 m/s, recharge mean: n/a
- 2023-04: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-04 stats — psi mean: -63.193 m, sw mean: 0.323, ETa mean: 8.356e-09 m/s, recharge mean: n/a
- 2023-05: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2023-05 stats — psi mean: -50.046 m, sw mean: 0.343, ETa mean: 9.323e-09 m/s, recharge mean: n/a
- 2023-06: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-06 stats — psi mean: -38.294 m, sw mean: 0.356, ETa mean: 2.387e-08 m/s, recharge mean: n/a
- 2023-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-07 stats — psi mean: -57.774 m, sw mean: 0.329, ETa mean: 1.098e-08 m/s, recharge mean: n/a
- 2023-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-08 stats — psi mean: -60.389 m, sw mean: 0.326, ETa mean: 8.659e-09 m/s, recharge mean: n/a
- 2023-09: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-09 stats — psi mean: -40.291 m, sw mean: 0.340, ETa mean: 1.881e-08 m/s, recharge mean: n/a
- 2023-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2023-10 stats — psi mean: -51.110 m, sw mean: 0.330, ETa mean: 9.103e-09 m/s, recharge mean: n/a
- 2023-11: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-11 stats — psi mean: -59.162 m, sw mean: 0.325, ETa mean: 4.703e-09 m/s, recharge mean: n/a
- 2023-12: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-12 stats — psi mean: -59.967 m, sw mean: 0.324, ETa mean: 3.727e-09 m/s, recharge mean: n/a
- 2024-01: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2024-01 stats — psi mean: -50.254 m, sw mean: 0.327, ETa mean: 6.239e-09 m/s, recharge mean: n/a
- 2024-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-02 stats — psi mean: -57.697 m, sw mean: 0.324, ETa mean: 5.953e-09 m/s, recharge mean: n/a
- 2024-03: ETa artefact floor patched 825 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-03 stats — psi mean: -55.730 m, sw mean: 0.325, ETa mean: 9.190e-09 m/s, recharge mean: n/a
- 2024-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2024-04 stats — psi mean: -59.611 m, sw mean: 0.323, ETa mean: 8.108e-09 m/s, recharge mean: n/a
- 2024-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2024-05 stats — psi mean: -59.407 m, sw mean: 0.323, ETa mean: 1.150e-08 m/s, recharge mean: n/a
- 2024-06: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2024-06 stats — psi mean: -52.440 m, sw mean: 0.327, ETa mean: 1.977e-08 m/s, recharge mean: n/a
- 2024-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2024-07 stats — psi mean: -61.436 m, sw mean: 0.321, ETa mean: 1.030e-08 m/s, recharge mean: n/a
- 2024-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2024-08 stats — psi mean: -62.125 m, sw mean: 0.320, ETa mean: 7.475e-09 m/s, recharge mean: n/a
- 2024-09: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-09 stats — psi mean: -57.745 m, sw mean: 0.322, ETa mean: 7.978e-09 m/s, recharge mean: n/a
- 2024-10: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2024-10 stats — psi mean: -54.810 m, sw mean: 0.328, ETa mean: 5.446e-09 m/s, recharge mean: n/a
- 2024-11: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2024-11 stats — psi mean: -27.438 m, sw mean: 0.364, ETa mean: 1.543e-08 m/s, recharge mean: n/a
- 2024-12: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2024-12 stats — psi mean: -49.271 m, sw mean: 0.331, ETa mean: 7.846e-09 m/s, recharge mean: n/a
- 2025-01: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2025-01 stats — psi mean: -53.686 m, sw mean: 0.326, ETa mean: 3.879e-09 m/s, recharge mean: n/a
- 2025-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2025-02 stats — psi mean: -56.783 m, sw mean: 0.324, ETa mean: 4.233e-09 m/s, recharge mean: n/a
- 2025-03: ETa artefact floor patched 3872 cell(s) (net rain ≥ ETp); 3872 cell(s) set for net < 0 (mode=etp).
- 2025-03 stats — psi mean: -28.846 m, sw mean: 0.375, ETa mean: 1.724e-08 m/s, recharge mean: n/a
- 2025-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2025-04 stats — psi mean: -47.575 m, sw mean: 0.337, ETa mean: 1.419e-08 m/s, recharge mean: n/a
- 2025-05: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2025-05 stats — psi mean: -47.811 m, sw mean: 0.334, ETa mean: 1.865e-08 m/s, recharge mean: n/a
- 2025-06: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2025-06 stats — psi mean: -58.531 m, sw mean: 0.327, ETa mean: 9.742e-09 m/s, recharge mean: n/a
- 2025-07: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2025-07 stats — psi mean: -53.637 m, sw mean: 0.330, ETa mean: 1.435e-08 m/s, recharge mean: n/a
- 2025-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2025-08 stats — psi mean: -57.087 m, sw mean: 0.325, ETa mean: 9.472e-09 m/s, recharge mean: n/a
| 2025-09 | — | — | — | — | — | — | skipped (no forcing) |

## Run summary

- **Months run**: 107
- **Months skipped**: 12
- **Output directory**: /home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/Agramon_SShydro/hypotheses/hyp1_layers_Ks/prj/outputs/scenario_2_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2
- **ETa artefact cells patched (total)**: 90130
- **ETa net<0 cells set (total)**: 723101
- **Warnings logged**: 0
- **Elapsed**: 0:44:15

---

_Finished: 2026-11-20 12:53:42_
