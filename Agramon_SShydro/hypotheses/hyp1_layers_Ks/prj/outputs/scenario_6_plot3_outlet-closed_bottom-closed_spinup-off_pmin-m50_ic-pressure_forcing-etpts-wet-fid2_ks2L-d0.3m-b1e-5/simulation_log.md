# Simulation Log — scenario_6_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2_ks2L-d0.3m-b1e-5 (dem-plot 3)

_Started: 2026-11-20 13:38:35_


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
- **sim_index**: 6

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
- 2016-11 stats — psi mean: -84.292 m, sw mean: 0.311, ETa mean: 4.436e-09 m/s, recharge mean: n/a
- 2016-12: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2016-12 stats — psi mean: -49.698 m, sw mean: 0.384, ETa mean: 6.423e-09 m/s, recharge mean: n/a
- 2017-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-01 stats — psi mean: -44.415 m, sw mean: 0.371, ETa mean: 1.446e-08 m/s, recharge mean: n/a
- 2017-02: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2017-02 stats — psi mean: -45.568 m, sw mean: 0.355, ETa mean: 1.610e-08 m/s, recharge mean: n/a
- 2017-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-03 stats — psi mean: -49.720 m, sw mean: 0.355, ETa mean: 2.422e-08 m/s, recharge mean: n/a
- 2017-04: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-04 stats — psi mean: -61.648 m, sw mean: 0.336, ETa mean: 1.399e-08 m/s, recharge mean: n/a
- 2017-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-05 stats — psi mean: -62.619 m, sw mean: 0.332, ETa mean: 2.285e-08 m/s, recharge mean: n/a
- 2017-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-06 stats — psi mean: -71.841 m, sw mean: 0.320, ETa mean: 1.484e-08 m/s, recharge mean: n/a
- 2017-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2017-07 stats — psi mean: -74.256 m, sw mean: 0.317, ETa mean: 1.147e-08 m/s, recharge mean: n/a
- 2017-08: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2017-08 stats — psi mean: -71.297 m, sw mean: 0.319, ETa mean: 9.232e-09 m/s, recharge mean: n/a
- 2017-09: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-09 stats — psi mean: -58.206 m, sw mean: 0.332, ETa mean: 2.102e-08 m/s, recharge mean: n/a
- 2017-10: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2017-10 stats — psi mean: -68.970 m, sw mean: 0.317, ETa mean: 7.795e-09 m/s, recharge mean: n/a
- 2017-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2017-11 stats — psi mean: -70.837 m, sw mean: 0.316, ETa mean: 3.306e-09 m/s, recharge mean: n/a
- 2017-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2017-12 stats — psi mean: -52.134 m, sw mean: 0.325, ETa mean: 9.801e-09 m/s, recharge mean: n/a
- 2018-01: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-01 stats — psi mean: -60.699 m, sw mean: 0.321, ETa mean: 5.720e-09 m/s, recharge mean: n/a
- 2018-02: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 4840 cell(s) set for net < 0 (mode=etp).
- 2018-02 stats — psi mean: -43.936 m, sw mean: 0.342, ETa mean: 1.265e-08 m/s, recharge mean: n/a
- 2018-03: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2018-03 stats — psi mean: -49.979 m, sw mean: 0.331, ETa mean: 1.465e-08 m/s, recharge mean: n/a
- 2018-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2018-04 stats — psi mean: -66.784 m, sw mean: 0.317, ETa mean: 9.940e-09 m/s, recharge mean: n/a
- 2018-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-05 stats — psi mean: -70.541 m, sw mean: 0.315, ETa mean: 9.815e-09 m/s, recharge mean: n/a
- 2018-06: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2018-06 stats — psi mean: -63.675 m, sw mean: 0.320, ETa mean: 1.701e-08 m/s, recharge mean: n/a
- 2018-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2018-07 stats — psi mean: -74.829 m, sw mean: 0.312, ETa mean: 9.680e-09 m/s, recharge mean: n/a
- 2018-08: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-08 stats — psi mean: -73.252 m, sw mean: 0.312, ETa mean: 7.638e-09 m/s, recharge mean: n/a
- 2018-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-09 stats — psi mean: -59.687 m, sw mean: 0.320, ETa mean: 1.452e-08 m/s, recharge mean: n/a
- 2018-10: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2018-10 stats — psi mean: -64.090 m, sw mean: 0.316, ETa mean: 8.205e-09 m/s, recharge mean: n/a
- 2018-11: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2018-11 stats — psi mean: -48.710 m, sw mean: 0.340, ETa mean: 9.175e-09 m/s, recharge mean: n/a
- 2018-12: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2018-12 stats — psi mean: -51.092 m, sw mean: 0.329, ETa mean: 1.225e-08 m/s, recharge mean: n/a
- 2019-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2019-01 stats — psi mean: -67.556 m, sw mean: 0.316, ETa mean: 6.908e-09 m/s, recharge mean: n/a
- 2019-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2019-02 stats — psi mean: -71.348 m, sw mean: 0.314, ETa mean: 5.414e-09 m/s, recharge mean: n/a
- 2019-03: ETa artefact floor patched 13 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-03 stats — psi mean: -73.459 m, sw mean: 0.312, ETa mean: 5.858e-09 m/s, recharge mean: n/a
- 2019-04: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2019-04 stats — psi mean: -58.976 m, sw mean: 0.338, ETa mean: 1.538e-08 m/s, recharge mean: n/a
- 2019-05: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-05 stats — psi mean: -60.457 m, sw mean: 0.331, ETa mean: 2.160e-08 m/s, recharge mean: n/a
- 2019-06: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-06 stats — psi mean: -70.867 m, sw mean: 0.319, ETa mean: 1.236e-08 m/s, recharge mean: n/a
- 2019-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2019-07 stats — psi mean: -72.895 m, sw mean: 0.317, ETa mean: 1.051e-08 m/s, recharge mean: n/a
- 2019-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2019-08 stats — psi mean: -72.338 m, sw mean: 0.316, ETa mean: 9.090e-09 m/s, recharge mean: n/a
- 2019-09: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2019-09 stats — psi mean: -54.294 m, sw mean: 0.349, ETa mean: 1.852e-08 m/s, recharge mean: n/a
- 2019-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2019-10 stats — psi mean: -54.083 m, sw mean: 0.341, ETa mean: 1.660e-08 m/s, recharge mean: n/a
- 2019-11: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2019-11 stats — psi mean: -53.616 m, sw mean: 0.332, ETa mean: 1.243e-08 m/s, recharge mean: n/a
- 2019-12: ETa artefact floor patched 2178 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2019-12 stats — psi mean: -41.162 m, sw mean: 0.360, ETa mean: 1.605e-08 m/s, recharge mean: n/a
- 2020-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2020-01 stats — psi mean: -47.220 m, sw mean: 0.353, ETa mean: 1.102e-08 m/s, recharge mean: n/a
- 2020-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2020-02 stats — psi mean: -50.262 m, sw mean: 0.344, ETa mean: 1.712e-08 m/s, recharge mean: n/a
- 2020-03: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2020-03 stats — psi mean: -52.567 m, sw mean: 0.350, ETa mean: 9.868e-09 m/s, recharge mean: n/a
- 2020-04: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2020-04 stats — psi mean: -43.090 m, sw mean: 0.359, ETa mean: 2.467e-08 m/s, recharge mean: n/a
- 2020-05: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2020-05 stats — psi mean: -57.324 m, sw mean: 0.337, ETa mean: 2.349e-08 m/s, recharge mean: n/a
- 2020-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2020-06 stats — psi mean: -69.223 m, sw mean: 0.324, ETa mean: 1.350e-08 m/s, recharge mean: n/a
- 2020-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2020-07 stats — psi mean: -71.888 m, sw mean: 0.321, ETa mean: 1.257e-08 m/s, recharge mean: n/a
- 2020-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2020-08 stats — psi mean: -72.586 m, sw mean: 0.319, ETa mean: 1.116e-08 m/s, recharge mean: n/a
- 2020-09: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2020-09 stats — psi mean: -72.641 m, sw mean: 0.318, ETa mean: 7.695e-09 m/s, recharge mean: n/a
- 2020-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2020-10 stats — psi mean: -70.566 m, sw mean: 0.318, ETa mean: 6.762e-09 m/s, recharge mean: n/a
- 2020-11: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2020-11 stats — psi mean: -50.631 m, sw mean: 0.336, ETa mean: 1.281e-08 m/s, recharge mean: n/a
- 2020-12: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2020-12 stats — psi mean: -45.721 m, sw mean: 0.335, ETa mean: 1.139e-08 m/s, recharge mean: n/a
- 2021-01: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6297 cell(s) set for net < 0 (mode=etp).
- 2021-01 stats — psi mean: -46.859 m, sw mean: 0.348, ETa mean: 1.530e-08 m/s, recharge mean: n/a
- 2021-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2021-02 stats — psi mean: -58.457 m, sw mean: 0.328, ETa mean: 9.117e-09 m/s, recharge mean: n/a
- 2021-03: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2021-03 stats — psi mean: -53.787 m, sw mean: 0.332, ETa mean: 1.461e-08 m/s, recharge mean: n/a
- 2021-04: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2021-04 stats — psi mean: -57.075 m, sw mean: 0.331, ETa mean: 9.669e-09 m/s, recharge mean: n/a
- 2021-05: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2021-05 stats — psi mean: -55.584 m, sw mean: 0.330, ETa mean: 2.581e-08 m/s, recharge mean: n/a
- 2021-06: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2021-06 stats — psi mean: -59.760 m, sw mean: 0.326, ETa mean: 2.359e-08 m/s, recharge mean: n/a
- 2021-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2021-07 stats — psi mean: -72.041 m, sw mean: 0.317, ETa mean: 1.177e-08 m/s, recharge mean: n/a
- 2021-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2021-08 stats — psi mean: -72.855 m, sw mean: 0.316, ETa mean: 9.725e-09 m/s, recharge mean: n/a
- 2021-09: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2021-09 stats — psi mean: -56.848 m, sw mean: 0.328, ETa mean: 2.013e-08 m/s, recharge mean: n/a
- 2021-10: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2021-10 stats — psi mean: -60.933 m, sw mean: 0.321, ETa mean: 8.760e-09 m/s, recharge mean: n/a
- 2021-11: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2021-11 stats — psi mean: -47.586 m, sw mean: 0.332, ETa mean: 1.197e-08 m/s, recharge mean: n/a
- 2021-12: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2021-12 stats — psi mean: -55.286 m, sw mean: 0.322, ETa mean: 7.909e-09 m/s, recharge mean: n/a
- 2022-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-01 stats — psi mean: -65.083 m, sw mean: 0.317, ETa mean: 4.299e-09 m/s, recharge mean: n/a
- 2022-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2022-02 stats — psi mean: -69.164 m, sw mean: 0.316, ETa mean: 4.164e-09 m/s, recharge mean: n/a
- 2022-03: ETa artefact floor patched 4598 cell(s) (net rain ≥ ETp); 3146 cell(s) set for net < 0 (mode=etp).
- 2022-03 stats — psi mean: -39.855 m, sw mean: 0.373, ETa mean: 8.598e-09 m/s, recharge mean: n/a
- 2022-04: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5566 cell(s) set for net < 0 (mode=etp).
- 2022-04 stats — psi mean: -37.560 m, sw mean: 0.390, ETa mean: 3.153e-08 m/s, recharge mean: n/a
- 2022-05: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-05 stats — psi mean: -53.385 m, sw mean: 0.350, ETa mean: 2.652e-08 m/s, recharge mean: n/a
- 2022-06: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2022-06 stats — psi mean: -65.372 m, sw mean: 0.333, ETa mean: 1.914e-08 m/s, recharge mean: n/a
- 2022-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-07 stats — psi mean: -69.048 m, sw mean: 0.328, ETa mean: 1.311e-08 m/s, recharge mean: n/a
- 2022-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2022-08 stats — psi mean: -70.027 m, sw mean: 0.326, ETa mean: 1.164e-08 m/s, recharge mean: n/a
- 2022-09: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2022-09 stats — psi mean: -58.210 m, sw mean: 0.333, ETa mean: 1.630e-08 m/s, recharge mean: n/a
- 2022-10: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2022-10 stats — psi mean: -54.588 m, sw mean: 0.331, ETa mean: 1.410e-08 m/s, recharge mean: n/a
- 2022-11: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2022-11 stats — psi mean: -65.512 m, sw mean: 0.325, ETa mean: 6.037e-09 m/s, recharge mean: n/a
- ⚠️ **WARNING**: 2022-12: solver did not converge (PSI last index 259216.754 < expected 2592000.0) — keeping last known state.
- 2022-12: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 1210 cell(s) set for net < 0 (mode=etp).
- 2022-12 stats — psi mean: n/a, sw mean: 0.323, ETa mean: 1.225e-09 m/s, recharge mean: n/a
- 2023-01: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-01 stats — psi mean: -68.718 m, sw mean: 0.323, ETa mean: 3.624e-09 m/s, recharge mean: n/a
- 2023-02: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-02 stats — psi mean: -61.494 m, sw mean: 0.325, ETa mean: 4.823e-09 m/s, recharge mean: n/a
- 2023-03: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-03 stats — psi mean: -70.202 m, sw mean: 0.321, ETa mean: 6.579e-09 m/s, recharge mean: n/a
- 2023-04: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-04 stats — psi mean: -71.218 m, sw mean: 0.320, ETa mean: 8.561e-09 m/s, recharge mean: n/a
- 2023-05: ETa artefact floor patched 2420 cell(s) (net rain ≥ ETp); 5324 cell(s) set for net < 0 (mode=etp).
- 2023-05 stats — psi mean: -58.276 m, sw mean: 0.340, ETa mean: 9.426e-09 m/s, recharge mean: n/a
- 2023-06: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-06 stats — psi mean: -46.388 m, sw mean: 0.354, ETa mean: 2.428e-08 m/s, recharge mean: n/a
- 2023-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2023-07 stats — psi mean: -65.989 m, sw mean: 0.327, ETa mean: 1.119e-08 m/s, recharge mean: n/a
- 2023-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-08 stats — psi mean: -68.883 m, sw mean: 0.323, ETa mean: 8.752e-09 m/s, recharge mean: n/a
- 2023-09: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2023-09 stats — psi mean: -48.914 m, sw mean: 0.337, ETa mean: 1.880e-08 m/s, recharge mean: n/a
- 2023-10: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2023-10 stats — psi mean: -59.682 m, sw mean: 0.327, ETa mean: 9.171e-09 m/s, recharge mean: n/a
- 2023-11: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-11 stats — psi mean: -67.743 m, sw mean: 0.322, ETa mean: 4.812e-09 m/s, recharge mean: n/a
- 2023-12: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2023-12 stats — psi mean: -68.626 m, sw mean: 0.321, ETa mean: 3.815e-09 m/s, recharge mean: n/a
- 2024-01: ETa artefact floor patched 1694 cell(s) (net rain ≥ ETp); 6050 cell(s) set for net < 0 (mode=etp).
- 2024-01 stats — psi mean: -58.935 m, sw mean: 0.324, ETa mean: 6.310e-09 m/s, recharge mean: n/a
- 2024-02: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-02 stats — psi mean: -66.486 m, sw mean: 0.321, ETa mean: 6.059e-09 m/s, recharge mean: n/a
- 2024-03: ETa artefact floor patched 809 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-03 stats — psi mean: -64.656 m, sw mean: 0.322, ETa mean: 9.288e-09 m/s, recharge mean: n/a
- 2024-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2024-04 stats — psi mean: -68.669 m, sw mean: 0.320, ETa mean: 8.204e-09 m/s, recharge mean: n/a
- 2024-05: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2024-05 stats — psi mean: -68.579 m, sw mean: 0.320, ETa mean: 1.160e-08 m/s, recharge mean: n/a
- 2024-06: ETa artefact floor patched 968 cell(s) (net rain ≥ ETp); 6534 cell(s) set for net < 0 (mode=etp).
- 2024-06 stats — psi mean: -61.714 m, sw mean: 0.324, ETa mean: 1.982e-08 m/s, recharge mean: n/a
- 2024-07: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2024-07 stats — psi mean: -70.726 m, sw mean: 0.318, ETa mean: 1.037e-08 m/s, recharge mean: n/a
- 2024-08: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2024-08 stats — psi mean: -71.474 m, sw mean: 0.317, ETa mean: 7.541e-09 m/s, recharge mean: n/a
- 2024-09: ETa artefact floor patched 726 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2024-09 stats — psi mean: -67.112 m, sw mean: 0.319, ETa mean: 8.041e-09 m/s, recharge mean: n/a
- 2024-10: ETa artefact floor patched 1936 cell(s) (net rain ≥ ETp); 5808 cell(s) set for net < 0 (mode=etp).
- 2024-10 stats — psi mean: -64.227 m, sw mean: 0.325, ETa mean: 5.492e-09 m/s, recharge mean: n/a
- 2024-11: ETa artefact floor patched 1210 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2024-11 stats — psi mean: -36.700 m, sw mean: 0.362, ETa mean: 1.546e-08 m/s, recharge mean: n/a
- 2024-12: ETa artefact floor patched 4 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2024-12 stats — psi mean: -57.983 m, sw mean: 0.328, ETa mean: 8.106e-09 m/s, recharge mean: n/a
- 2025-01: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2025-01 stats — psi mean: -62.686 m, sw mean: 0.323, ETa mean: 3.961e-09 m/s, recharge mean: n/a
- 2025-02: ETa artefact floor patched 242 cell(s) (net rain ≥ ETp); 6776 cell(s) set for net < 0 (mode=etp).
- 2025-02 stats — psi mean: -65.927 m, sw mean: 0.321, ETa mean: 4.302e-09 m/s, recharge mean: n/a
- 2025-03: ETa artefact floor patched 3872 cell(s) (net rain ≥ ETp); 3872 cell(s) set for net < 0 (mode=etp).
- 2025-03 stats — psi mean: -38.144 m, sw mean: 0.373, ETa mean: 1.734e-08 m/s, recharge mean: n/a
- 2025-04: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7018 cell(s) set for net < 0 (mode=etp).
- 2025-04 stats — psi mean: -55.914 m, sw mean: 0.335, ETa mean: 1.492e-08 m/s, recharge mean: n/a
- 2025-05: ETa artefact floor patched 1452 cell(s) (net rain ≥ ETp); 6292 cell(s) set for net < 0 (mode=etp).
- 2025-05 stats — psi mean: -56.824 m, sw mean: 0.331, ETa mean: 1.886e-08 m/s, recharge mean: n/a
- 2025-06: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7502 cell(s) set for net < 0 (mode=etp).
- 2025-06 stats — psi mean: -67.904 m, sw mean: 0.323, ETa mean: 9.874e-09 m/s, recharge mean: n/a
- 2025-07: ETa artefact floor patched 484 cell(s) (net rain ≥ ETp); 7260 cell(s) set for net < 0 (mode=etp).
- 2025-07 stats — psi mean: -63.113 m, sw mean: 0.327, ETa mean: 1.442e-08 m/s, recharge mean: n/a
- 2025-08: ETa artefact floor patched 0 cell(s) (net rain ≥ ETp); 7744 cell(s) set for net < 0 (mode=etp).
- 2025-08 stats — psi mean: -66.494 m, sw mean: 0.322, ETa mean: 9.581e-09 m/s, recharge mean: n/a
| 2025-09 | — | — | — | — | — | — | skipped (no forcing) |

## Run summary

- **Months run**: 107
- **Months skipped**: 12
- **Output directory**: /home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/Agramon_SShydro/hypotheses/hyp1_layers_Ks/prj/outputs/scenario_6_plot3_outlet-closed_bottom-closed_spinup-off_pmin-m50_ic-pressure_forcing-etpts-wet-fid2_ks2L-d0.3m-b1e-5
- **ETa artefact cells patched (total)**: 89640
- **ETa net<0 cells set (total)**: 717293
- **Warnings logged**: 1
- **Elapsed**: 0:46:09

---

_Finished: 2026-11-20 14:24:44_
