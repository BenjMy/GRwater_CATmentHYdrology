"""
Agramon_withLAI_withETp.py — Long-term subsurface hydrology with LAI + ETp forcing (CATHY)
==========================================================================================
GRwater project — Agramon study site

Entry point of the pyCATHY pipeline: monthly Earth-Observation LAI updates
the vegetation / root-depth map during the run, and the ATMBC boundary is
forced from satellite (or CSV time-series) ETp and rain.

Pipeline (see ``run_pipeline``)
-------------------------------
  1.  Resolve DEM, shapefiles, forcing and LAI inputs      (agramon_forcing)
  2.  Log the simulation parameters (CSV + markdown)       (agramon_runlog)
  3.  Load and pre-process the DEM, build the mesh
  4.  Map ETp and rain onto mesh nodes; build full ATMBC arrays
  5.  Configure soil, ICs, BCs, and output settings
  6.  Optional spin-up (pre-2016 equilibration)
  7.  For each calendar month: update veg map from LAI -> run processor
                                                           (agramon_lai)
  8.  Monthly sanity checks / ETa artefact patch           (agramon_checks)
  9.  Aggregate and save outputs (NetCDF)

Related modules (same folder): agramon_config (paths + constants),
agramon_forcing (input loaders), agramon_lai (LAI -> veg -> soil),
agramon_checks (sanity checks, ETa patch), agramon_runlog (run logs,
scenario folder names), ks_fire_utils (post-fire Ks recovery, not yet
wired in). See REFACTORING.md for the refactoring plan.

Usage
-----
    python Agramon_withLAI_withETp.py
    python Agramon_withLAI_withETp.py --start-year 2019 --end-year 2022
    python Agramon_withLAI_withETp.py --zroot 0.5 --pmin -10
    python Agramon_withLAI_withETp.py --dem-plot 4
    python Agramon_withLAI_withETp.py --no-plots
    python Agramon_withLAI_withETp.py --zroot-factors 0.1 0.3 1.0 1.5 2.5
    python Agramon_withLAI_withETp.py --permx 1e-5 --permy-ratio 1.0 --permz-ratio 0.1

Forcing switch (``--forcing``, default FORCING_MODE_DEFAULT in agramon_config)
-----------------------------------------------------------------
      "etp"    : satellite ETp raster + satellite rain raster (default)
      "etp_ts" : ETp from the Agramon ET CSV time series (one value per day,
                 identical for every pixel; treatment picked from --dem-plot,
                 column from --etp-ts-column: ET_0_wet | ET_0_dry) + satellite rain
      "era5"   : DEPRECATED, quick-test only (AgUtils.extract_point_timeseries)
"""

import argparse
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import matplotlib.pyplot as plt
import rioxarray as rxr  # noqa: F401  (activates .rio accessor)
from scipy.ndimage import gaussian_filter
from pyCATHY import meshtools as mt


from utils import Agramon_utils as AgUtils
from utils import ks_fire_utils as KsFire
from agramon_ks import (FIRE_DATE_DEFAULT, ks_banner_label, ks_config_from_args,
                        ks_for_month, ks_key, ks_tag)
from pyCATHY import cathy_tools, cathy_utils
from pyCATHY.petro import Feddes_litterature  # literature-based Feddes lookup

from agramon_checks import (
    add_datetime_col,
    apply_eta_artefact_floor,
    check_monthly_sanity,
)
from agramon_config import (
    ETP_TS_COLUMN_DEFAULT,
    FORCING_LABELS,
    FORCING_MODE_DEFAULT,
    ForcingMode,
    EO_PATH,
    LAI_PATH,
    N_VEG_CLASSES,
    NATIVE_CRS,
    PathET,
    SCF_DEFAULT,
    TARGET_CRS,
    ZROOT_FACTORS,
    etp_path,
    rain_path,
    scf_label_for,
    feddes_params_from_args,
    scf_type,
    zroot_factors_from_args,
)
from agramon_mesh import (  # noqa: F401  (moved, REFACTORING 2.3)
    configure_boundary_conditions,
    initialise_project,
    load_dem_from_tif,
    setup_dem,
    smooth_dem,
)
from agramon_diagnostics import save_forcing_diagnostics  # noqa: F401  (moved, REFACTORING 2.6)
from agramon_spinup import pick_default_spinup_year, run_spinup  # noqa: F401  (moved, 2.5)
from agramon_forcing import (
    build_etp_ds_from_timeseries,
    load_all,
    load_era5,
    load_monthly_lai,
    load_rain,
    load_spatial_etp,
    subset_atmbc_for_month,
)
from agramon_lai import (
    iter_monthly_lai,
    lai_to_veg_map,
    smooth_veg_map,
    update_soil_for_active_zones,
)
from agramon_outputs import save_combined_outputs
from agramon_runlog import MarkdownLog, build_scenario_dirname, log_simulation


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Long-Term Subsurface Hydrology Simulation with LAI + ETp — pyCATHY",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    sim = p.add_argument_group("Simulation parameters")
    sim.add_argument("--start-year", type=int, default=2022, metavar="YEAR", #2019
                     help="First year of forcing to load")
    sim.add_argument("--end-year",   type=int, default=2022, metavar="YEAR", #2022
                     help="Last year (exclusive) of forcing")
    sim.add_argument("--zroot",  type=float, default=1.0,  metavar="M",
                     help="Base root depth (m) used as multiplier reference")
    sim.add_argument("--pmin",   type=float, default=-15, metavar="M",
                     help="Minimum pressure head threshold (m)")
    sim.add_argument("--scf",    type=scf_type, default=SCF_DEFAULT,
                     metavar="[0-1]",
                     help="Root water uptake stress compensation factor, a "
                          "float between 0 and 1: 1 = Feddes approach (ET ~ "
                          "all transpiration, the usual assumption), 0 = "
                          "PMIN-threshold approach (sparse/arid cover, bare "
                          "soil + plants). Values in between blend the two.")
    sim.add_argument("--log-file",
                     default="simulation_log_withETp_withLAI.csv",
                     metavar="PATH",
                     help="Shared CSV log file")

    veg = p.add_argument_group("Root depth mapping")
    veg.add_argument(
        "--zroot-factors", type=float, nargs=5, default=list(ZROOT_FACTORS.values()),
        metavar=("BARE", "SPARSE", "MODERATE", "DENSE", "VERY_DENSE"),
        help="Root-depth multiplier for each of the 5 fixed LAI-derived "
             "vegetation classes, applied per-class as "
             "ZROOT = --zroot * factor (replaces the hardcoded "
             "ZROOT_FACTORS mapping in update_soil_for_active_zones()). "
             "Must be given as 5 values in class order: bare soil, sparse, "
             "moderate, dense, very dense vegetation."
    )

    feddes = p.add_argument_group(
        "Feddes parameters (root water uptake — FP_map)")
    feddes.add_argument(
        "--feddes-crop", default=None, choices=["Corn", "Soybeans"],
        metavar="{Corn,Soybeans}",
        help="Infer PCANA/PCREF/PCWLT from literature values via "
             "pyCATHY.petro.Feddes_litterature(crop=...) instead of "
             "pyCATHY's hardcoded FP_map_default — pyCATHY's own tool for "
             "inferring Feddes parameters. Used as the base value for "
             "those three columns; --pcana/--pcref/--pcwlt below still "
             "override it if also given. Has no effect on PZ/OMGC (no "
             "literature equivalent) or on ZROOT (--zroot/--zroot-factors "
             "handle that separately). Default: unset, i.e. leave "
             "pyCATHY's own FP_map default untouched."
    )
    feddes.add_argument(
        "--pcana", type=float, default=None, metavar="M",
        help="Anaerobiosis point h2 (m, ≈ -0.5 m) — FP_map column PCANA. "
             "Applied uniformly across every row of df_fp (all veg "
             "zones), same broadcast pattern as --permx for SPP columns. "
             "Overrides --feddes-crop's literature value if both given. "
             "Default: unset (pyCATHY default, or --feddes-crop's value)."
    )
    feddes.add_argument(
        "--pcref", type=float, default=None, metavar="M",
        help="Field capacity / incipient water stress h3 (m, ≈ -4 m) — "
             "FP_map column PCREF. Same broadcast/override behaviour as "
             "--pcana. Default: unset (pyCATHY default, or "
             "--feddes-crop's value)."
    )
    feddes.add_argument(
        "--pcwlt", type=float, default=None, metavar="M",
        help="Wilting point h4 (m, ≈ -150 m) — FP_map column PCWLT. Same "
             "broadcast/override behaviour as --pcana. Default: unset "
             "(pyCATHY default, or --feddes-crop's value)."
    )
    feddes.add_argument(
        "--pz", type=float, default=None, metavar="EXP",
        help="Empirical root-water-uptake exponent PZ (FP_map column PZ; "
             "pyCATHY default 1.0). Applied uniformly across every row of "
             "df_fp. No literature equivalent in --feddes-crop. Default: "
             "unset, i.e. leave pyCATHY's own FP_map default untouched."
    )
    feddes.add_argument(
        "--omgc", type=float, default=None, metavar="0-1",
        help="Compensatory mechanism for root water uptake, 0<OMGC<1 "
             "(FP_map column OMGC; pyCATHY default 1.0). Applied "
             "uniformly across every row of df_fp. No literature "
             "equivalent in --feddes-crop. Default: unset, i.e. leave "
             "pyCATHY's own FP_map default untouched."
    )

    soil = p.add_argument_group("Soil hydraulic conductivity (PERMX / Ks)")
    soil.add_argument(
        "--permx", type=float, default=None, metavar="M/S",
        help="Override the saturated hydraulic conductivity Ks in the x "
             "direction (SOIL SPP column PERMX, m/s) for every soil zone. "
             "Applied uniformly across the SPP_map returned by "
             "simu.set_SOIL_defaults(SPP_map_default=True). Default: leave "
             "pyCATHY's own SPP default untouched."
    )
    soil.add_argument(
        "--permy-ratio", type=float, default=1.0, metavar="RATIO",
        help="Anisotropy ratio: PERMY = --permx * this value. Only takes "
             "effect when --permx is set."
    )
    soil.add_argument(
        "--permz-ratio", type=float, default=1.0, metavar="RATIO",
        help="Anisotropy ratio: PERMZ = --permx * this value. Only takes "
             "effect when --permx is set."
    )

    bc = p.add_argument_group("Boundary conditions")
    bc.add_argument("--outlet-flux", type=float, default=0.0, metavar="M/S",
                    help="Prescribed Neumann outflow (m/s, negative = leaving "
                         "the domain) applied on the topographically lowest "
                         "mesh side, so the outlet is open instead of a "
                         "closed box. Order of magnitude should reflect the "
                         "outlet-zone soil's saturated conductivity — check "
                         "this value against your Ks before trusting results. "
                         "Internally converted to a volumetric rate (m^3/s) "
                         "using the outlet side's average nodal face area "
                         "before being written to nansfneubc — see "
                         "configure_boundary_conditions() for details.")
    bc.add_argument("--outlet-side", default=None,
                    choices=["xmin_bound", "xmax_bound", "ymin_bound", "ymax_bound"],
                    help="Force which mesh side is treated as the outlet "
                         "instead of auto-detecting it from DEM elevation.")
    bc.add_argument("--bottom-flux", type=float, default=0.0, metavar="M/S",
                    help="Prescribed Neumann outflow (m/s, negative = leaving "
                         "the domain) applied to every node on the mesh "
                         "bottom — a free-bottom-drainage-style "
                         "approximation, same spirit as --outlet-flux but "
                         "for the base of the domain instead of the lateral "
                         "outlet. Default 0.0 keeps the bottom fully closed "
                         "(no-flow), i.e. the original behaviour. Converted "
                         "to a volumetric rate (m^3/s) per node using the "
                         "DEM pixel area (res_x * res_y) — see "
                         "configure_boundary_conditions() for details.")

    ic = p.add_argument_group("Initial conditions")
    ic.add_argument("--ic-mode", default="wt", choices=["wt", "pressure"],
                    help="How the initial pressure-head field is set: "
                         "'wt' places a hydrostatic water table at depth "
                         "--ic-wtposition below the surface (simu.update_ic("
                         "INDP=3, WTPOSITION=...)); 'pressure' sets a "
                         "uniform pressure head everywhere via "
                         "--ic-pressure-head (simu.update_ic(INDP=0, "
                         "pressure_head_ini=...)). Only affects the initial "
                         "condition of the very first run — subsequent "
                         "spin-up cycles / monthly steps still carry "
                         "forward the previous run's psi field regardless "
                         "of this setting.")
    ic.add_argument("--ic-wtposition", type=float, default=1.0, metavar="M",
                    help="Water-table depth (m) below the surface, used "
                         "when --ic-mode=wt")
    ic.add_argument("--ic-pressure-head", type=float, default=-10.0,
                    metavar="M",
                    help="Uniform initial pressure head (m), used when "
                         "--ic-mode=pressure")

    spinup = p.add_argument_group("Spin-up (pre-2016 equilibration)")
    spinup.add_argument("--no-spinup", action="store_true",
                        help="Disable pre-run spin-up (spin-up runs by default)")
    spinup.add_argument("--spinup-cycles", type=int, default=0, metavar="N",
                        help="Number of times to repeat the spin-up year")
    spinup.add_argument("--spinup-year", type=int, default=None, metavar="YEAR",
                        help="Calendar year of forcing to cycle for spin-up "
                             "(default: first complete year available — "
                             "satellite ETp/rain starts 2016-10-01, so 2016 "
                             "is skipped automatically)")

    sim.add_argument("--forcing", choices=list(FORCING_LABELS), default=FORCING_MODE_DEFAULT,
                     help="Atmospheric forcing: era5 = ERA5 point series | "
                          "etp = satellite ETp raster + rain raster | "
                          "etp_ts = ETp CSV time series (uniform over pixels, "
                          "see --etp-ts-*) + rain raster")

    sim.add_argument("--eta-net-neg-mode", choices=["etp", "gross", "off"],
                     default="etp",
                     help="'ACT. ETRA_patched' where net = rain - ETp < 0: "
                          "etp = set to ETp (default) | gross = min(ETp, "
                          "ETRA + rain), keeps soil-limited cells below ETp | "
                          "off = keep raw ACT. ETRA")

    etpts = p.add_argument_group(
        "ETp time-series forcing (only used with --forcing etp_ts)")
    etpts.add_argument("--etp-ts-column", default=ETP_TS_COLUMN_DEFAULT,
                       choices=["ET_0_wet", "ET_0_dry"],
                       help="ETp bound used as forcing (mm/day): ET_0_wet or "
                            "ET_0_dry column of agramon_treatments_ET.csv")
    etpts.add_argument("--etp-ts-feature-id", type=int, default=None, metavar="ID",
                       help="Treatment feature_id to use. Default: derived from "
                            "--dem-plot as feature_id = dem_plot - 1")

    dem = p.add_argument_group("DEM selection")
    dem.add_argument("--dem-tif", default=str(EO_PATH / "DTMplots/20250618_AGRAMON100m_micasense_dtm.tif"),
                     metavar="PATH",
                     help="Full path to the GeoTIFF DEM used to build the mesh "
                          "(default: 20250618_AGRAMON100m_micasense_dtm.tif)")
    dem.add_argument("--dem-plot",   type=int, default=3, choices=range(1, 14), metavar="N",
                     help="Plot index (1-13, matches gdf_Agramon's fid_1 column) used to "
                          "clip the DEM (loaded from --dem-tif) to that plot's footprint")
    dem.add_argument("--dem-folder", default=None, metavar="PATH",
                     help="Legacy override: full path to .adf raster folder "
                          "(unused now that --dem-tif is the DEM source)")
    dem.add_argument("--dem-smooth-sigma", type=float, default=0.0, metavar="SIGMA",
                     help="Gaussian smoothing sigma (in pixels) applied to the DEM "
                          "elevation after masking (default: 0.0 = no smoothing). "
                          "Array shape and NaN mask are preserved exactly — only "
                          "valid-pixel elevations are blended, so grid dims (N, M) "
                          "stay unaffected. Typical values: 0.5-2.0.")

    solver = p.add_argument_group("Solver parameters")
    solver.add_argument("--dtmin",  type=float, default=1e-1, metavar="S")
    solver.add_argument("--dtmax",  type=float, default=1e4,  metavar="S")
    solver.add_argument("--deltat", type=float, default=1e3,  metavar="S")

    out = p.add_argument_group("Output")
    out.add_argument(
        "--path2prj",
        default=str(Path("../SSHydro_withETp_withLAI/").resolve()),
        metavar="PATH",
        help="Parent directory for the CATHY project folder"
    )
    out.add_argument("--no-plots", action="store_true",
                     help="Skip all matplotlib visualisations")
    out.add_argument("--debug", action="store_true",
                     help="Print extra diagnostics (LAI interpolation NaN "
                          "checks and per-catchment LAI coverage). Off by default.")
    out.add_argument("--save-monthly", action="store_true",
                     help="Also write a standalone NetCDF file per month per "
                          "variable to outputs/<scenario>/monthly/ (psi, sw, "
                          "et, recharge). Off by default — the single running "
                          "psi_output.nc / sw_output.nc / et_output.nc / "
                          "recharge_output.nc files in outputs/<scenario>/ "
                          "are always refreshed after every month regardless.")

    # post-fire Ks recovery (opt-in; default 'off' = static --permx, unchanged behaviour)
    KsFire.add_ks_fire_args(p)
    p.set_defaults(ks_fire_date=FIRE_DATE_DEFAULT)
    return p.parse_args()


# ── Helper: CRS alignment check ───────────────────────────────────────────────

def _check_crs_alignment(label_a: str, crs_a, label_b: str, crs_b) -> None:
    if str(crs_a) != str(crs_b):
        print(f"  ⚠ CRS mismatch: {label_a}={crs_a} vs {label_b}={crs_b}")
    else:
        print(f"  ✓ CRS aligned ({label_a} & {label_b}): {crs_a}")


# ── Pipeline steps ────────────────────────────────────────────────────────────

def load_shapefiles():
    """Load Agramon plot shapefiles, ensuring TARGET_CRS."""
    #gdf = AgUtils.load_plot_shapefiles(MODULE_PATH / "shapefiles")
    gdf = AgUtils.load_plot_shapefiles(['microcuencas_13'],
                                               EO_PATH / "data/microcuencas_13")

    if gdf.crs is None:
        raise ValueError("Shapefile has no CRS — assign one before proceeding.")
    if str(gdf.crs) != TARGET_CRS:
        print(f"  Reprojecting shapefiles {gdf.crs} → {TARGET_CRS}")
        gdf = gdf.to_crs(TARGET_CRS)
    print(f"  Shapefile CRS : {gdf.crs}")
    return gdf


def plot_forcing(t_atmbc, net_flux, label="ETp forcing") -> None:
    """Bar chart of daily net flux forcing."""
    fig, ax = plt.subplots(figsize=(15, 4))
    ax.bar(range(len(net_flux)), net_flux, width=1.0, color="skyblue")
    ax.set_title(label)
    ax.set_ylabel("m/s")
    ax.set_xlabel("Time step")
    plt.tight_layout()


# ── Spin-up (pre-2016 equilibration) ────────────────────────────────────────


# ── Main pipeline ─────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class RunInputs:
    """Everything ``run_pipeline`` needs from ``setup_inputs``."""
    show_plots: bool
    dem_tif_path: Path
    feddes: dict
    gdf_Agramon: Any        # GeoDataFrame, micro-catchment polygons
    ds_LAI: Any             # xarray.Dataset, monthly LAI
    ds_rain: Any            # xarray.Dataset or None (ERA5 mode)
    ds_etp: Any             # xarray.Dataset or None (ERA5 mode)
    ds_era5: Any            # xarray.Dataset or None (ETp modes)
    etp_ts_fid: Any         # treatment id used in "etp_ts" mode
    forcing: ForcingMode    # immutable forcing choice (was module globals, 2.13)
    scenario_dirname: str
    md_log: MarkdownLog


def setup_inputs(args: argparse.Namespace) -> RunInputs:

    """
    Steps [1/9]-[3/9] of the pipeline: resolve labels and the run banner,
    load shapefiles, forcing and LAI, log the simulation (CSV + markdown).

    Moved verbatim out of ``run_pipeline``. The forcing choice is an immutable
    ``ForcingMode`` built here from ``--forcing`` and returned in ``RunInputs.forcing``
    (it replaces the former module globals, REFACTORING.md 2.13).
    """
    fm = ForcingMode.from_name(args.forcing)
    show_plots = not args.no_plots

    # ── 1. Resolve DEM source ─────────────────────────────────
    dem_tif_path = Path(args.dem_tif)

    # SCF is taken straight from --scf and reused everywhere below (soil
    # updates, CSV log, markdown log) so it can never drift between call
    # sites.
    scf = args.scf
    scf_label = scf_label_for(scf)

    # Feddes (FP_map) parameters are resolved once here too, same reasoning
    # as SCF above — --feddes-crop's literature lookup vs. explicit
    # --pcana/--pcref/--pcwlt/--pz/--omgc overrides — so the CSV log,
    # markdown log, and every update_soil_for_active_zones() call in the
    # main loop below all agree.
    feddes = feddes_params_from_args(args)
    feddes_label = ", ".join(
        f"{k.upper()}={v:g}" for k, v in feddes.items()
        if k in ("pcana", "pcref", "pcwlt", "pz", "omgc") and v is not None
    ) or "pyCATHY default"
    if feddes["feddes_crop"] is not None:
        feddes_label += f" (base: --feddes-crop {feddes['feddes_crop']})"

    # Banner labels are built here (not inline in the f-strings) so the code
    # stays valid on Python < 3.12, where an f-string replacement field can't
    # span lines or reuse the enclosing quote character.
    forcing_label = fm.label
    if args.permx is None:
        banner_permx_label = "pyCATHY default"
    else:
        banner_permx_label = (
            f"{args.permx:.3e} m/s"
            f"  |  PERMY ratio {args.permy_ratio}"
            f"  |  PERMZ ratio {args.permz_ratio}"
        )
    ks_cfg = ks_config_from_args(args)      # validates --ks-fire-* now, not at month 1
    spinup_banner_label = (
        "disabled" if args.no_spinup else f"{args.spinup_cycles} cycle(s)"
    )

    print("\n" + "=" * 60)
    print("  LT Subsurface Hydrology with LAI + ETp — pyCATHY")
    print(f"  Forcing source : {forcing_label}")
    print(f"  Period         : {args.start_year} – {args.end_year - 1}")
    print(f"  Target CRS     : {TARGET_CRS}")
    print(f"  ZROOT          : {args.zroot} m  |  PMIN : {args.pmin} m  |  SCF : {scf_label}")
    print(f"  ZROOT factors  : {args.zroot_factors} "
          f"(bare/sparse/moderate/dense/very_dense)")
    print(f"  Feddes (FP)    : {feddes_label}")
    print(f"  PERMX (Ks)     : {banner_permx_label}")
    print(f"  Post-fire Ks   : {ks_banner_label(ks_cfg)}")
    print(f"  DEM source     : {dem_tif_path}")
    print(f"  Log file       : {args.log_file}")
    print(f"  Spin-up        : {spinup_banner_label}")
    print("=" * 60)

    # ── 2. Shapefiles ─────────────────────────────────────────
    print("\n[1/9] Loading shapefiles …")
    gdf_Agramon = load_shapefiles()

    # ── 3. Forcing datasets + LAI ─────────────────────────────
    print("\n[2/9] Loading forcing datasets and LAI …")
    ds_LAI = load_monthly_lai(LAI_PATH)

    if fm.use_etp:
        ds_rain = load_rain()
        # Satellite forcing is defined from 2016-10-01 (see pick_default_spinup_year);
        # the rain file also contains earlier years the CSV ETp doesn't cover.
        ds_rain = ds_rain.sel(time=slice("2016-10-01", None))
        if fm.use_etp_ts:
            df_et = load_all(PathET)
            ds_etp = build_etp_ds_from_timeseries(
                df_et, ds_rain,
                plot_id=args.dem_plot,
                column=args.etp_ts_column,
                feature_id=args.etp_ts_feature_id,
            )
        else:
            ds_etp = load_spatial_etp()
        ds_era5 = None
    else:
        ds_era5 = load_era5(args.start_year, args.end_year)
        ds_etp  = None
        ds_rain = None

    # ── 4. Log simulation ─────────────────────────────────────
    print("\n[3/9] Logging simulation parameters …")
    log_path  = Path(".").resolve() / args.log_file
    # Treatment actually used in "etp_ts" mode (same rule as
    # build_etp_ds_from_timeseries: explicit --etp-ts-feature-id, else plot-1).
    etp_ts_fid = (args.etp_ts_feature_id if args.etp_ts_feature_id is not None
                  else args.dem_plot - 1)
    sim_index = log_simulation(
        log_path,
        start_year=args.start_year,
        end_year=args.end_year,
        zroot=args.zroot,
        pmin=args.pmin,
        scf=scf,
        watershed_nb=args.dem_plot,
        with_lai=1,
        outlet_flux=args.outlet_flux,
        bottom_flux=args.bottom_flux,
        spinup=("off" if args.no_spinup else f"on_{args.spinup_cycles}cycles"),
        ic_mode=args.ic_mode,
        ic_wtposition=args.ic_wtposition,
        ic_pressure_head=args.ic_pressure_head,
        zroot_factors=args.zroot_factors,
        permx=args.permx,
        permy_ratio=args.permy_ratio,
        permz_ratio=args.permz_ratio,
        pcana=feddes["pcana"], pcref=feddes["pcref"], pcwlt=feddes["pcwlt"],
        pz=feddes["pz"], omgc=feddes["omgc"],
        feddes_crop=feddes["feddes_crop"],
        forcing_mode=fm.name,
        etp_ts_column=args.etp_ts_column,
        etp_ts_feature_id=etp_ts_fid,
        ks_fire=ks_key(args),
    )

    # Single source of truth for this run's folder name (project dir,
    # outputs dir, markdown log, and — via withLAI_results.py — the
    # figures dir all reuse it), so the same run always lands in the
    # same, comparably-named place. See build_scenario_dirname() above.
    scenario_dirname = build_scenario_dirname(
        sim_index, args.dem_plot, args.outlet_flux, args.bottom_flux,
        args.no_spinup, args.spinup_cycles, args.pmin, args.ic_mode,
        forcing_tag="_".join(t for t in (
            (f"forcing-etpts-{args.etp_ts_column.removeprefix('ET_0_')}"
             f"-fid{etp_ts_fid}" if fm.use_etp_ts else ""),
            ks_tag(args)) if t),
    )

    # ── Markdown run log ───────────────────────────────────────
    # Same output directory the final NetCDFs land in (out_dir, defined
    # further below); computed again here from scenario_dirname/args.path2prj
    # so the log can start capturing steps immediately instead of waiting
    # for out_dir's later definition to be reordered.
    md_log_path = (Path(args.path2prj) / "outputs" / scenario_dirname
                    / "simulation_log.md")
    md_log = MarkdownLog(
        md_log_path,
        title=f"Simulation Log — {scenario_dirname} (dem-plot {args.dem_plot})",
    )
    if fm.use_etp:
        md_log.fixes_applied([
            "ETa net-forcing artefact floor: `apply_eta_artefact_floor()` adds "
            "an `ACT. ETRA_patched` column, flooring cells where CATHY reports "
            "`ACT. ETRA ≈ 0` because net forcing (rain − ETp) was ≥ 0 that "
            "day/node. Cells with net < 0 (ETp > rain) are set per "
            f"`--eta-net-neg-mode` (= `{args.eta_net_neg_mode}`; `etp` → "
            "patched ETa = ETp). Original `ACT. ETRA` is left unmodified.",
        ])
    md_log.section("Run configuration")
    md_log.meta({
        "Forcing source": fm.label,
        "Period": f"{args.start_year} – {args.end_year - 1}",
        "Target CRS": TARGET_CRS,
        "ZROOT / PMIN / SCF": f"{args.zroot} m / {args.pmin} m / {scf_label}",
        "ZROOT factors (bare/sparse/moderate/dense/very_dense)":
            ",".join(f"{v:g}" for v in args.zroot_factors),
        "PERMX / Ks override": (
            "pyCATHY default" if args.permx is None
            else f"PERMX={args.permx:.3e} m/s, PERMY ratio={args.permy_ratio}, "
                 f"PERMZ ratio={args.permz_ratio}"
        ),
        "Feddes (FP_map) override": feddes_label,
        "Initial condition": (f"WT depth {args.ic_wtposition} m (INDP=3)"
                               if args.ic_mode == "wt"
                               else f"uniform pressure head "
                                    f"{args.ic_pressure_head} m (INDP=0)"),
        "DEM source (tif)": str(dem_tif_path),
        "DEM plot (fid_1 clip)": args.dem_plot,
        "Spin-up": ("disabled" if args.no_spinup
                     else f"{args.spinup_cycles} cycle(s)"),
        "Outlet flux": f"{args.outlet_flux:.1e} m/s",
        "Outlet side": args.outlet_side or "auto-detect",
        "Bottom flux": (f"{args.bottom_flux:.1e} m/s"
                        if args.bottom_flux != 0.0 else "disabled (closed)"),
        "Solver dtmin/dtmax/deltat": f"{args.dtmin} / {args.dtmax} / {args.deltat} s",
        "CSV run registry": str(log_path),
        "sim_index": sim_index,
    })

    return RunInputs(
        show_plots=show_plots,
        dem_tif_path=dem_tif_path,
        feddes=feddes,
        gdf_Agramon=gdf_Agramon,
        ds_LAI=ds_LAI,
        ds_rain=ds_rain,
        ds_etp=ds_etp,
        ds_era5=ds_era5,
        etp_ts_fid=etp_ts_fid,
        forcing=fm,
        scenario_dirname=scenario_dirname,
        md_log=md_log,
    )


@dataclass(frozen=True)
class MeshSetup:
    """Everything ``run_pipeline`` needs from ``setup_mesh_and_dem``."""
    simu: Any               # pyCATHY.CATHY project (preprocessor already run)
    raster_DEM: Any         # xarray.DataArray, DEM clipped to the catchment bbox
    raster_DEM_masked: Any  # np.ndarray, DEM with NaN outside the catchment
    xllcorner: float
    yllcorner: float
    res_x: float
    res_y: float
    grid3d: dict            # freshly read grid3d (see the n_surf note in the ATMBC step)
    N: int                  # grid rows
    M: int                  # grid columns


def setup_mesh_and_dem(args: argparse.Namespace, inp: RunInputs) -> MeshSetup:
    """
    Steps [4/9]-[7/9] of the pipeline: load/mask/smooth the DEM, check CRS
    alignment, initialise the CATHY project, feed the DEM into it, and run
    the CATHY preprocessor.

    Moved verbatim out of ``run_pipeline``. It stays in this module because
    ``setup_dem`` / ``initialise_project`` / ``load_dem_from_tif`` live here
    until step 2.3 (``agramon_mesh.py``).
    """
    show_plots       = inp.show_plots
    dem_tif_path     = inp.dem_tif_path
    gdf_Agramon      = inp.gdf_Agramon
    ds_LAI           = inp.ds_LAI
    scenario_dirname = inp.scenario_dirname
    md_log           = inp.md_log

    # ── 5. DEM ────────────────────────────────────────────────
    print("\n[4/9] Loading and masking DEM …")
    raster_DEM, raster_DEM_masked, xllcorner, yllcorner, res_x, res_y = load_dem_from_tif(
        dem_tif_path,
        gdf_clip=gdf_Agramon,
        fid=args.dem_plot,
        resample_resolution=5
    )


    _check_crs_alignment("DEM",   raster_DEM.rio.crs,  "target", TARGET_CRS)
    _check_crs_alignment("LAI",   ds_LAI.rio.crs,       "target", TARGET_CRS)
    _check_crs_alignment("shape", gdf_Agramon.crs,       "target", TARGET_CRS)

    # ── 5b. Optional DEM smoothing (shape/NaN mask preserved) ──
    if args.dem_smooth_sigma > 0:
        shape_before = raster_DEM_masked.shape
        valid_before = int(np.isfinite(raster_DEM_masked).sum())
        raster_DEM_masked = smooth_dem(raster_DEM_masked, args.dem_smooth_sigma)
        assert raster_DEM_masked.shape == shape_before, (
            "smooth_dem changed DEM shape — this must never happen "
            f"({raster_DEM_masked.shape} != {shape_before})"
        )
        assert int(np.isfinite(raster_DEM_masked).sum()) == valid_before, (
            "smooth_dem changed the valid-pixel (NaN mask) footprint — "
            "this must never happen"
        )
        print(f"  DEM smoothed: gaussian sigma={args.dem_smooth_sigma} px "
              f"(shape unchanged: {shape_before})")

    md_log.section("Setup")
    dem_valid_px  = int(np.isfinite(raster_DEM_masked).sum())
    dem_area_m2   = dem_valid_px * abs(res_x) * abs(res_y)
    md_log.meta({
        "DEM shape (rows x cols)": f"{raster_DEM_masked.shape}",
        "DEM CRS": str(raster_DEM.rio.crs),
        "DEM pixel size (dx, dy)": f"({res_x:.3f}, {res_y:.3f}) m",
        "DEM origin (xllcorner, yllcorner)": f"({xllcorner:.1f}, {yllcorner:.1f})",
        "DEM valid area": f"{dem_valid_px} px  ≈  {dem_area_m2 / 1e4:.2f} ha "
                           f"({dem_area_m2:.0f} m²)",
        "DEM smoothing (gaussian sigma, px)": (
            args.dem_smooth_sigma if args.dem_smooth_sigma > 0 else "none"
        ),
        "LAI CRS": str(ds_LAI.rio.crs),
        "Shapefile CRS": str(gdf_Agramon.crs),
    })

    # ── 6. Initialise CATHY project ───────────────────────────
    print("\n[5/9] Initialising CATHY project …")
    simu = initialise_project(args.path2prj, args.start_year, args.end_year,
                              args.zroot, args.pmin, scenario_dirname)

    # ── 7. DEM → mesh ─────────────────────────────────────────
    print("\n[6/9] Feeding DEM into CATHY and building VTK mesh …")
    setup_dem(simu, raster_DEM_masked,
              xllcorner=xllcorner, yllcorner=yllcorner,
              delta_x=res_x, delta_y=res_y,
              show=show_plots)

    # ── 8. Preprocessor + grid dims ───────────────────────────
    print("\n[7/9] Running CATHY preprocessor …")
    simu.run_preprocessor(verbose=False)
    grid3d = simu.read_outputs("grid3d")
    N = int(simu.hapin["N"])
    M = int(simu.hapin["M"])
    print(f"  Grid : {N} rows × {M} cols")
    md_log.bullet(f"Grid: {N} rows × {M} cols")

    return MeshSetup(
        simu=simu,
        raster_DEM=raster_DEM,
        raster_DEM_masked=raster_DEM_masked,
        xllcorner=xllcorner,
        yllcorner=yllcorner,
        res_x=res_x,
        res_y=res_y,
        grid3d=grid3d,
        N=N,
        M=M,
    )


@dataclass(frozen=True)
class AtmbcSetup:
    """Everything ``run_pipeline`` needs from ``build_atmbc``."""
    ds_mesh: Any            # xarray.Dataset, mesh nodes (x, y, z) built from the DEM
    n_surf: Any             # int, surface node count; None in ERA5 mode
    ETp_nodes: Any          # ndarray (time, n_surf) in m/s; None in ERA5 mode
    rain_nodes: Any         # ndarray (time, n_surf) in m/s; None in ERA5 mode
    net_flux: Any           # ndarray (time,) scalar mean net flux (rain - ETp)
    t_atmbc: Any            # seconds since the first forcing step
    valid_times: Any        # datetime64 time axis of the forcing


def build_atmbc(args: argparse.Namespace, inp: RunInputs,
                ms: MeshSetup) -> AtmbcSetup:
    """
    Step [8/9] of the pipeline: build the mesh dataset, map ETp and rain to
    the surface nodes (or read the ERA5 point series), and build the ATMBC
    forcing arrays.

    Moved verbatim out of ``run_pipeline``. It stays in this module because
    it reads the run's ``ForcingMode`` (``inp.forcing``, REFACTORING.md 2.13) and calls
    ``save_forcing_diagnostics``.

    ``ETp_nodes``, ``rain_nodes`` and ``n_surf`` are only bound by the ETp
    branch. In ERA5 mode the original code left them unbound (and every later
    read is guarded by ``fm.use_etp``); here they are ``None`` so they
    can be returned.
    """
    fm                = inp.forcing
    show_plots        = inp.show_plots
    gdf_Agramon       = inp.gdf_Agramon
    ds_rain           = inp.ds_rain
    ds_etp            = inp.ds_etp
    ds_era5           = inp.ds_era5
    etp_ts_fid        = inp.etp_ts_fid
    md_log            = inp.md_log
    simu              = ms.simu
    grid3d            = ms.grid3d
    raster_DEM_masked = ms.raster_DEM_masked

    ETp_nodes  = None
    rain_nodes = None
    n_surf     = None

    # ── 9. Build ATMBC forcing arrays ─────────────────────────
    print("\n[8/9] Building ATMBC forcing arrays …")

    ds_mesh = mt.build_mesh_dataset(simu, raster_DEM_masked=raster_DEM_masked)
    print(f"  Mesh x : {float(ds_mesh.x.min()):.1f} – {float(ds_mesh.x.max()):.1f}")
    print(f"  Mesh y : {float(ds_mesh.y.min()):.1f} – {float(ds_mesh.y.max()):.1f}")
    md_log.meta({
        "Mesh x range": f"[{float(ds_mesh.x.min()):.1f}, {float(ds_mesh.x.max()):.1f}] m",
        "Mesh y range": f"[{float(ds_mesh.y.min()):.1f}, {float(ds_mesh.y.max()):.1f}] m",
        "Mesh CRS": TARGET_CRS,
    })

    if fm.use_etp:
        # ── Map ETp and rain to mesh surface nodes ────────────
        # FIX: `simu.grid3d` can be stale (only updated as a side effect of
        # certain internal calls, not by read_outputs()). Use the `grid3d`
        # dict we just freshly read a few lines above instead, so n_surf
        # matches the mesh actually built for THIS run.
        n_surf = int(grid3d["nnod"])

        ds_ETp_mapped  = mt.map_grid_to_mesh(ds_etp,  ds_mesh, variables=["ET_0-gf"])
        ds_rain_mapped = mt.map_grid_to_mesh(ds_rain, ds_mesh, variables=["TP-DD"])

        # Restrict to surface nodes only
        ds_ETp_surf  = ds_ETp_mapped.isel( node=slice(0, n_surf))
        ds_rain_surf = ds_rain_mapped.isel(node=slice(0, n_surf))

        # Align on common time axis
        common_time  = np.intersect1d(ds_rain_mapped.time, ds_ETp_mapped.time)
        ds_ETp_surf  = ds_ETp_surf.sel(time=common_time)
        ds_rain_surf = ds_rain_surf.sel(time=common_time)

        # Unit conversions: ETp mm/day → m/s (negative = loss);  rain mm/day → m/s
        ETp_nodes  = ds_ETp_surf["ET_0-gf"].values  * 1e-3 / 86400   # (time, nodes)
        rain_nodes = ds_rain_surf["TP-DD"].values   * 1e-3 / 86400   # (time, nodes)

        # Net flux per node (positive = input to soil)
        net_flux_2d = rain_nodes - ETp_nodes                          # (time, nodes)
        # net_flux_2d = - ETp_nodes                          # (time, nodes)

        # Scalar mean over nodes for ATMBC scalar mode
        net_flux_1d = net_flux_2d.mean(axis=1)                        # (time,)

        # Reference time: first common time step
        t0          = common_time[0]
        t_atmbc     = (common_time - t0) / np.timedelta64(1, "s")    # seconds
        valid_times = common_time

        print(f"  ATMBC steps (ETp) : {len(t_atmbc)}  |  flux range : "
              f"[{net_flux_1d.min():.2e}, {net_flux_1d.max():.2e}] m/s")

        # FIX: full provenance record for the forcing data — native/target
        # CRS, native pixel size, mesh node count, and total accumulated
        # depths over the whole forcing period (not just per-day ranges),
        # so unit/quantity mistakes are visible in the log without having
        # to reopen the source NetCDFs.
        etp_res_native  = ds_etp.rio.resolution()
        rain_res_native = ds_rain.rio.resolution()
        n_days_forcing  = len(valid_times)
        total_rain_mm   = float((rain_nodes.mean(axis=1) * 86400 * 1e3).sum())
        total_etp_mm    = float((ETp_nodes.mean(axis=1) * 86400 * 1e3).sum())
        md_log.section("Atmospheric forcing — provenance & units")
        md_log.meta({
            "ETp source": (
                f"{PathET / 'agramon_treatments_ET.csv'} "
                f"(feature_id={etp_ts_fid}, "
                f"column={args.etp_ts_column}; same value for every pixel)"
                if fm.use_etp_ts else str(etp_path)),
            "Rain source file": str(rain_path),
            "Native CRS (source)": NATIVE_CRS,
            "Target CRS (mesh)": TARGET_CRS,
            "ETp native pixel size (dx, dy)": (
                f"({etp_res_native[0]:.1f}, {etp_res_native[1]:.1f}) m"),
            "Rain native pixel size (dx, dy)": (
                f"({rain_res_native[0]:.1f}, {rain_res_native[1]:.1f}) m"),
            "Surface mesh nodes forced": n_surf,
            "Forcing period": (f"{pd.Timestamp(valid_times[0]):%Y-%m-%d} → "
                                f"{pd.Timestamp(valid_times[-1]):%Y-%m-%d} "
                                f"({n_days_forcing} daily steps)"),
            "Units": "ETp/rain source: mm/day; ATMBC/CATHY forcing: m/s "
                     "(spatial mean over mesh surface nodes, negative = ET loss)",
            "ATMBC steps (ETp)": len(t_atmbc),
            "Net flux range (m/s)": f"[{net_flux_1d.min():.2e}, {net_flux_1d.max():.2e}]",
            "ETp range (mm/day)": (f"[{ETp_nodes.min()*1e3*86400:.2f}, "
                                     f"{ETp_nodes.max()*1e3*86400:.2f}]"),
            "Rain range (mm/day)": (f"[{rain_nodes.min()*1e3*86400:.2f}, "
                                      f"{rain_nodes.max()*1e3*86400:.2f}]"),
            "Total rain over period (spatial mean)": f"{total_rain_mm:.1f} mm",
            "Total ETp over period (spatial mean)": f"{total_etp_mm:.1f} mm",
        })

        # net_flux passed to subset_atmbc_for_month (1-D scalar per time step)
        net_flux = net_flux_1d

        forcing_dir = Path("forcing_plots")
        forcing_dir.mkdir(exist_ok=True)
        # fig.savefig(forcing_dir / f"veg_map_{dt:%Y%m}.png", dpi=300, bbox_inches="tight")
        # plt.close(fig)

        # forcing_dir = (
        #         Path(args.path2prj)
        #         / "outputs"
        #         / f"scenario_{sim_index}"
        #         / "forcing_plots"
        #     )

        save_forcing_diagnostics(
            valid_times=valid_times,
            net_flux_1d=net_flux_1d,
            rain_nodes=rain_nodes,
            ETp_nodes=ETp_nodes,
            net_flux_2d=net_flux_2d,
            outdir=forcing_dir,
        )

    else:
        # ── ERA5 point time-series ────────────────────────────
        daily_ts = AgUtils.extract_point_timeseries(ds_era5, gdf_Agramon)
        if show_plots:
            pev_series = daily_ts["pev"].to_series()
            tp_series  = daily_ts["tp"].to_series()
            fig, ax = plt.subplots(figsize=(15, 4))
            ax.bar(pev_series.index, pev_series.values, width=1.0,
                   color="skyblue", label="PEV")
            ax.bar(tp_series.index,  tp_series.values,  width=1.0,
                   color="orange",  label="TP")
            ax.set_title("Daily PEV and TP forcing (ERA5)")
            ax.set_ylabel("mm/day")
            ax.legend()
            plt.tight_layout()

        _t0         = daily_ts["pev"].valid_time[0].values
        valid_times = daily_ts["pev"].valid_time.values
        t_atmbc     = [(float(t - _t0) / 1e9) for t in valid_times]  # ns → s
        net_flux    = np.array(daily_ts["tp"]) + np.array(daily_ts["pev"])
        print(f"  ATMBC steps (ERA5) : {len(t_atmbc)}  |  flux range : "
              f"[{net_flux.min():.2e}, {net_flux.max():.2e}] m")
        md_log.meta({
            "ATMBC steps (ERA5)": len(t_atmbc),
            "Net flux range (m)": f"[{net_flux.min():.2e}, {net_flux.max():.2e}]",
        })

    return AtmbcSetup(
        ds_mesh=ds_mesh,
        n_surf=n_surf,
        ETp_nodes=ETp_nodes,
        rain_nodes=rain_nodes,
        net_flux=net_flux,
        t_atmbc=t_atmbc,
        valid_times=valid_times,
    )


def run_pipeline(args: argparse.Namespace) -> None:
    """Execute the full pyCATHY + LAI + ETp simulation pipeline."""

    inp = setup_inputs(args)
    fm = inp.forcing

    ks_cfg = ks_config_from_args(args)   # per-month post-fire Ks (off by default)
    feddes           = inp.feddes
    gdf_Agramon      = inp.gdf_Agramon
    ds_LAI           = inp.ds_LAI
    ds_rain          = inp.ds_rain
    ds_etp           = inp.ds_etp
    scenario_dirname = inp.scenario_dirname
    md_log           = inp.md_log

    # ── [4/9]–[7/9] DEM, CATHY project, mesh, preprocessor ────
    ms = setup_mesh_and_dem(args, inp)
    simu              = ms.simu
    raster_DEM        = ms.raster_DEM
    raster_DEM_masked = ms.raster_DEM_masked
    res_x             = ms.res_x
    res_y             = ms.res_y
    grid3d            = ms.grid3d
    N                 = ms.N
    M                 = ms.M

    # ── [8/9] ATMBC forcing arrays ────────────────────────────
    atm = build_atmbc(args, inp, ms)
    ds_mesh     = atm.ds_mesh
    n_surf      = atm.n_surf
    ETp_nodes   = atm.ETp_nodes
    rain_nodes  = atm.rain_nodes
    net_flux    = atm.net_flux
    t_atmbc     = atm.t_atmbc
    valid_times = atm.valid_times

    # ── 10. ICs, BCs ──────────────────────────────────────────
    # FIX: set ICs once here; BCs are built once, up front.
    if args.ic_mode == "pressure":
        print(f"  IC             : uniform pressure head = "
              f"{args.ic_pressure_head} m (INDP=0)")
        simu.update_ic(INDP=0, pressure_head_ini=args.ic_pressure_head)
    else:
        print(f"  IC             : hydrostatic water table at "
              f"{args.ic_wtposition} m depth (INDP=3)")
        simu.update_ic(INDP=3, WTPOSITION=args.ic_wtposition)
    # FIX (sfbc EOF crash): configure_boundary_conditions used to be fed
    # list(simu.read_inputs("atmbc")["time"]) here - but at this point in
    # the pipeline simu.update_atmbc() has NEVER been called yet (it's
    # only called later, once per month, inside the main loop), so that
    # read back whatever short placeholder/template atmbc file pyCATHY
    # wrote during mesh setup (~1 day's worth of entries). update_sfbc()
    # writes exactly one BC block per entry in the time list it's given,
    # and unlike nansfdirbc/nansfneubc's reader (bcnxt.f, which has an
    # END= label and just freezes on EOF), the seepage-face reader
    # (sfvnxt.f) has no END= handling and hard-crashes with a Fortran
    # "End of file" runtime error once the solver's TIME runs past the
    # last declared block - which happened right around the 1-day mark,
    # matching the crash log exactly.
    #
    # Fix: build the BC files from a synthetic time axis instead of
    # reading (the wrong) atmbc file off disk. Content here is identical
    # at every block ("no-flow") regardless of spacing, so all that's
    # needed is an axis that outlasts any single month's local run - each
    # run_processor() call starts its own clock at 0 (see
    # subset_atmbc_for_month's re-zeroing), so this only needs to cover
    # one calendar month, not the whole multi-year forcing period.
    bc_time_axis = np.arange(0, 32 * 86400, 86400).tolist()  # 32 daily blocks
    outlet_side, outlet_node_id = configure_boundary_conditions(
        simu, bc_time_axis,
        outlet_flux=args.outlet_flux,
        outlet_side=args.outlet_side,
        bottom_flux=args.bottom_flux,
        res_x=res_x, res_y=res_y,
    )
    md_log.meta({
        "Resolved outlet side": outlet_side,
        "Outlet diagnostic node": outlet_node_id,
        "Bottom flux": (f"{args.bottom_flux:.1e} m/s"
                        if args.bottom_flux != 0.0 else "disabled (closed)"),
    })

    # Interpolate LAI onto DEM raster grid once, outside the loop
    print("  Interpolating LAI timeseries → DEM raster grid …")
    x_t, y_t = raster_DEM.x.values, raster_DEM.y.values
    lai_lin  = ds_LAI.interp(x=xr.DataArray(x_t, dims="x"),
                             y=xr.DataArray(y_t, dims="y"), method="linear")
    lai_near = ds_LAI.reindex(x=x_t, y=y_t, method="nearest")
    ds_LAI_raster = lai_lin.fillna(lai_near)

    # Diagnostics on the LAI interpolation (stdout only); off by default.
    if args.debug:
        x0, x1 = float(raster_DEM.x.min()), float(raster_DEM.x.max())
        y0, y1 = float(raster_DEM.y.min()), float(raster_DEM.y.max())
        raw = ds_LAI["LAI"].isel(time=1).where(
            (ds_LAI.x > x0-20) & (ds_LAI.x < x1+20) &
            (ds_LAI.y > y0-20) & (ds_LAI.y < y1+20), drop=True)
        print("raw LAI patch shape:", raw.shape, "| NaN fraction:", float(np.isnan(raw).mean()))
        print(np.isnan(raw.values).astype(int))        # 1 = NaN in the source LAI

        valid = raster_DEM_masked != -9999
        v = ds_LAI_raster["LAI"].isel(time=1).values
        print("NaN inside DEM footprint after interp:", int(np.isnan(v[valid]).sum()), "of", int(valid.sum()))
        print("NaN outside footprint:", int(np.isnan(v[~valid]).sum()), "of", int((~valid).sum()))


        print(gdf_Agramon[["fid_1"] + [c for c in ("PlotID",) if c in gdf_Agramon]].to_string())

        lai_any = ds_LAI["LAI"].notnull().any("time").astype(float)   # LAI exists in ≥1 month
        lai_all = ds_LAI["LAI"].notnull().all("time").astype(float)   # LAI exists in every month
        for fid, g in gdf_Agramon.groupby("fid_1"):
            a = lai_any.rio.clip(g.geometry, all_touched=True, drop=True)
            b = lai_all.rio.clip(g.geometry, all_touched=True, drop=True)
            print(f"fid_1={fid:>2}  bounds={[round(v) for v in g.total_bounds]}  "
                  f"LAI valid in ≥1 month: {np.nanmean(a):.0%} | in all months: {np.nanmean(b):.0%}")

    print(f"  ds_LAI_raster shape : {ds_LAI_raster['LAI'].shape}")
    md_log.bullet(f"ds_LAI_raster shape: {ds_LAI_raster['LAI'].shape}")

    # ── 10b. Spin-up (pre-2016 equilibration) ─────────────────
    first_iteration = True
    spinup_done     = False
    if not args.no_spinup:
        print("\n[Spin-up] Equilibrating initial pressure-head field "
              "(spin-up antes de 2016) …")
        first_iteration, spinup_done = run_spinup(
            simu=simu,
            args=args,
            t_atmbc=t_atmbc,
            net_flux=net_flux,
            valid_times=valid_times,
            ds_LAI_raster=ds_LAI_raster,
            first_iteration=first_iteration,
            outlet_node_id=outlet_node_id,
            md_log=md_log,
        )
    else:
        print("\n[Spin-up] Disabled (--no-spinup) — starting cold from "
              "the uniform IC.")
        md_log.section("Spin-up (warm-up)")
        md_log.bullet("Disabled (`--no-spinup`) — starting cold from the uniform IC.")

    # ── 11. Storage containers ────────────────────────────────
    df_psi = None
    psi_records       = []
    sw_records        = []
    et_records        = []
    all_recharge_recs = []
    all_hgraph_recs   = []
    veg_map_history   = {}
    lai_history       = {}
    stats_records     = []

    active_zones_prev = None

    out_dir = Path(args.path2prj) / "outputs" / scenario_dirname
    out_dir.mkdir(parents=True, exist_ok=True)

    # Persist the per-node ETp forcing actually used to drive ATMBC. This
    # is a static (time, node) array built once above (step 9) and never
    # touched again, so it's written straight away rather than accumulated
    # through the monthly loop like psi/sw/et. Downstream analysis (e.g.
    # ETa/ETp vs. shallow pressure head, Fig. 8 of Camporese et al. 2014)
    # needs ETp per node and otherwise has no way to recover it without
    # re-running the satellite-to-mesh mapping (mt.map_grid_to_mesh).
    # Only meaningful in the ETp-forcing branch — ERA5-forced runs never
    # build a per-node ETp array.
    if fm.use_etp and ETp_nodes is not None:
        etp_out_path = out_dir / "etp_output.nc"
        xr.DataArray(
            ETp_nodes,
            dims=["datetime", "node"],
            coords={"datetime": valid_times, "node": np.arange(ETp_nodes.shape[1])},
            name="ETp",
            attrs={
                "units": "m/s",
                "description": (
                    "Potential ET forcing per surface node, as mapped from "
                    + ("the uniform CSV ETp time series" if fm.use_etp_ts
                       else "the satellite ET0 raster")
                    + " onto the mesh (node index "
                    "matches the first n_surf rows of grid3d['mesh3d_nodes'])."
                ),
            },
        ).to_netcdf(etp_out_path)
        print(f"  ETp forcing saved -> {etp_out_path}")

    # Per-month outputs are written here as soon as each month finishes,
    # so a crash mid-run still leaves usable files on disk instead of
    # losing everything unwritten (the combined *_output.nc files at the
    # end of the loop are still produced too, for convenience).
    # Per-month outputs are written here as soon as each month finishes,
    # only when --save-monthly is passed (off by default). The single
    # running combined files below are always refreshed every month
    # regardless of this flag.
    monthly_dir = out_dir / "monthly"
    if args.save_monthly:
        monthly_dir.mkdir(parents=True, exist_ok=True)

    # ── 12. Main LAI loop — one CATHY run per calendar month ──
    print("\n── Starting LAI → vegetation-map loop ──────────────────────────────")
    md_log.section("Monthly simulation loop")
    md_log.table_header(["Month", "LAI mean", "LAI min–max", "Net flux mean (mm/day)",
                          "Veg zones", "Vegetated %", "LAI updated", "Solver"])
    n_months_run     = 0
    n_months_skipped = 0

    for dt, lai_2d in iter_monthly_lai(ds_LAI_raster.isel(time=slice(1, None))):
        # 12.1  LAI → consecutive vegetation map + remap dict
        veg_map_remapped, remap = lai_to_veg_map(simu, lai_2d)
        # 12.2  Smooth salt-and-pepper artefacts
        veg_map = smooth_veg_map(veg_map_remapped)

        fig, (ax_lai, ax_veg) = plt.subplots(1, 2, figsize=(8, 6))

        # LAI plot
        im_lai = ax_lai.imshow(lai_2d, cmap="YlGn")
        ax_lai.set_title(f"LAI - {dt:%Y-%m}")
        fig.colorbar(im_lai, ax=ax_lai, label="LAI")

        # Vegetation map plot
        im_veg = ax_veg.imshow(veg_map)
        ax_veg.set_title(f"Vegetation map - {dt:%Y-%m}")
        fig.colorbar(im_veg, ax=ax_veg)

        fig.suptitle(f"{dt:%Y-%m}", fontsize=13, fontweight="bold")
        fig.tight_layout()

        outdir_veg = Path("figures_lai_temp")
        outdir_veg.mkdir(exist_ok=True)
        fig.savefig(outdir_veg / f"veg_map_{dt:%Y%m}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        exclude_veg, veg_valid = simu._check_outside_DEM(veg_map)
        simu.update_veg_map(veg_map, maxveg=N_VEG_CLASSES)

        inv_remap             = {v: k for k, v in remap.items()}
        veg_map_semantic      = np.vectorize(inv_remap.get)(veg_map).astype(int)
        lai_history[dt]       = lai_2d.copy()
        veg_map_history[dt]   = veg_map_semantic

        # 12.3  Update soil only when active zones change
        active_zones = set(int(z) for z in veg_valid)
        lai_updated  = (active_zones != active_zones_prev)
        # if lai_updated:
        update_soil_for_active_zones(simu,
                                     active_zones,
                                     remap=remap,
                                     base_zroot=args.zroot,
                                     pmin=args.pmin,
                                     scf=args.scf,
                                     zroot_factors=zroot_factors_from_args(args),
                                     permx=ks_for_month(ks_cfg, dt, lai_2d, args.permx),
                                     permy_ratio=args.permy_ratio,
                                     permz_ratio=args.permz_ratio,
                                     pcana=feddes["pcana"], pcref=feddes["pcref"],
                                     pcwlt=feddes["pcwlt"], pz=feddes["pz"],
                                     omgc=feddes["omgc"],
                                     feddes_crop=feddes["feddes_crop"],
                                     )
        active_zones_prev = active_zones

        # 12.4  Subset ATMBC to current month
        t_month, nv_month = subset_atmbc_for_month(
            t_atmbc, net_flux, valid_times, dt.year, dt.month
        )
        if t_month is None:
            print(f"  {dt:%Y-%m-%d} | No forcing data — skipping")
            n_months_skipped += 1
            md_log.table_row([f"{dt:%Y-%m}", "—", "—", "—", "—", "—", "—",
                               "skipped (no forcing)"])
            continue

        simu.update_atmbc(HSPATM=1,
                          IETO=1,
                          time=t_month,
                          netValue=nv_month
                          )

        # Indices into ETp_nodes/rain_nodes' time axis for this month —
        # needed by both the ETa artefact floor (12.10) and the sanity
        # checks (12.9b). Computed once here so both reuse the same mask.
        idx_month = None
        if fm.use_etp:
            month_mask = np.array([
                pd.Timestamp(t).year == dt.year and pd.Timestamp(t).month == dt.month
                for t in valid_times
            ])
            idx_month = np.where(month_mask)[0]


        # Optional: plot the monthly forcing
        outdir_etp = Path("figures_ETp_temp")
        outdir_etp.mkdir(exist_ok=True)
        fig2, ax2 = plt.subplots(figsize=(8, 4))
        ax2.plot(nv_month)
        ax2.set_title(f"Net flux ({'ETp' if fm.use_etp else 'ERA5'}) - {dt:%Y-%m}")
        ax2.set_ylabel("m/s")
        fig2.savefig(outdir_etp / f"forcing_{dt:%Y%m}.png", dpi=150, bbox_inches="tight")
        plt.close(fig2)

        # 12.5  Update PARM print list
        df_atmbc_month = simu.read_inputs("atmbc")
        simu.update_parm(TIMPRTi=list(df_atmbc_month.time),
                         IPRT=4,
                         VTKF=2
                         )

        # 12.6  Preprocessor (first iteration only)
        # FIX: mirror the working withLAI logic exactly — run preprocessor
        # once on first iteration; BCs were already configured above from
        # the full time axis, not re-built here.
        if first_iteration:
            #simu.run_preprocessor(verbose=True)
            first_iteration = False

            # FIX: this preprocessor call can finalize a DIFFERENT (larger)
            # mesh than the one built earlier (setup_dem's create_mesh_vtk
            # + the very first run_preprocessor() before this loop). Two
            # things upstream were computed against that earlier, smaller
            # mesh and go stale the moment the grid actually changes size:
            #   1. `update_ic`'s internal map_prop2mesh() reuses the mesh
            #      cached by create_mesh_vtk() — never rebuilt afterwards,
            #      causing "'ic' has length N but the mesh has M points".
            #   2. `n_surf` / ETp_nodes / rain_nodes were sliced using the
            #      old grid3d["nnod"], causing the "SURFACE NODE out of
            #      bounds" warnings later when reading real ETa/recharge
            #      output (which is indexed on the true, larger grid).
            # Refresh both against the grid this run is actually using.
            simu.create_mesh_vtk(verbose=True)
            grid3d = simu.read_outputs("grid3d")
            n_surf_real = int(grid3d["nnod"])

            if fm.use_etp and n_surf_real != n_surf:
                print(f"  Refreshing ETp/rain node mapping: n_surf "
                      f"{n_surf} → {n_surf_real} (grid finalized by "
                      f"preprocessor)")
                n_surf = n_surf_real

                ds_mesh = mt.build_mesh_dataset(
                    simu, raster_DEM_masked=raster_DEM_masked
                )
                ds_ETp_mapped  = mt.map_grid_to_mesh(
                    ds_etp,  ds_mesh, variables=["ET_0-gf"]
                )
                ds_rain_mapped = mt.map_grid_to_mesh(
                    ds_rain, ds_mesh, variables=["TP-DD"]
                )

                ds_ETp_surf  = ds_ETp_mapped.isel( node=slice(0, n_surf))
                ds_rain_surf = ds_rain_mapped.isel(node=slice(0, n_surf))

                common_time_refresh = np.intersect1d(
                    ds_rain_mapped.time, ds_ETp_mapped.time
                )
                ds_ETp_surf  = ds_ETp_surf.sel(time=common_time_refresh)
                ds_rain_surf = ds_rain_surf.sel(time=common_time_refresh)

                ETp_nodes  = ds_ETp_surf["ET_0-gf"].values * 1e-3 / 86400
                rain_nodes = ds_rain_surf["TP-DD"].values  * 1e-3 / 86400

                # NOTE: t_atmbc / net_flux (the scalar, spatially-averaged
                # forcing already applied via update_atmbc for THIS month,
                # a few lines above) are not retroactively corrected here —
                # HSPATM=1 means CATHY only ever saw their spatial mean, so
                # this refresh mainly fixes the per-node ETp_nodes/
                # rain_nodes arrays used by the ETa/recharge sanity checks
                # from this month onward, plus (critically) the mesh used
                # by update_ic below.



        if len(psi_records) > 0 or spinup_done:
            df_psi = simu.read_outputs("psi").copy()
            psi_ini = df_psi.iloc[-1].values
            simu.update_ic(INDP=1,
                           pressure_head_ini=psi_ini
                           )


        # 12.7  Run processor
        print(f"\n  Running processor for {dt:%Y-%m} "
              f"(LAI updated: {lai_updated}) …")

        # MAXVEG already fixed to N_VEG_CLASSES once before the loop; no
        # need to re-set it (to a mismatched value of 10) every month.

        simu.run_processor(
            IPRT1=2,
            DTMIN=args.dtmin,
            DTMAX=args.dtmax,
            DELTAT=args.deltat,
            TRAFLAG=0,
            verbose=False,
        )

        # 12.8  PSI — check solver actually progressed; keep last known state if not
        month_stats = {"month": f"{dt:%Y-%m}"}
        psi_month_df = None
        sw_month_df  = None
        df_et        = None
        df_rec       = None
        solver_status = "✓ OK"
        try:
            df_psi = simu.read_outputs("psi").copy()
            last_psi_t   = df_psi.index[-1]
            last_month_t = t_month[-1]
            if last_psi_t < last_month_t:
                print(f"  WARNING: PSI last index ({last_psi_t}) is behind "
                      f"t_month[-1] ({last_month_t}) for {dt:%Y-%m} — "
                      f"solver did not converge; keeping last known state and continuing.")
                solver_status = f"⚠ not converged ({last_psi_t} < {last_month_t})"
                md_log.warn(f"{dt:%Y-%m}: solver did not converge "
                             f"(PSI last index {last_psi_t} < expected {last_month_t}) "
                             f"— keeping last known state.")
                # Do not append new psi_records; psi_ini will reuse previous last row
                # on the next iteration (len(psi_records) > 0 guard still holds).
            else:
                print(f"  time index ic: {last_psi_t} ✓ (t_month[-1]: {last_month_t})")
                df_psi["datetime"] = cathy_utils.change_x2date(df_psi.index, dt)
                psi_month_df = (
                    df_psi.reset_index(drop=True)
                          .melt(id_vars="datetime", var_name="node", value_name="psi")
                          .assign(node=lambda d: d["node"].astype(int))
                )
                psi_records.append(psi_month_df)
                month_stats["psi_mean"] = float(psi_month_df["psi"].mean())
                month_stats["psi_min"]  = float(psi_month_df["psi"].min())
                month_stats["psi_max"]  = float(psi_month_df["psi"].max())
                if args.save_monthly:
                    try:
                        (
                            psi_month_df
                            .groupby(["datetime", "node"], as_index=False)
                            .last()
                            .set_index(["datetime", "node"])
                            .to_xarray()
                            .to_netcdf(monthly_dir / f"psi_{dt:%Y%m}.nc")
                        )
                    except Exception as exc:
                        print(f"  Warning: could not save monthly PSI for {dt:%Y-%m} ({exc})")
        except Exception as exc:
            print(f"  WARNING: could not read PSI for {dt:%Y-%m} ({exc}) — "
                  f"keeping last known state and continuing.")
            solver_status = f"✗ PSI read failed: {exc}"
            md_log.warn(f"{dt:%Y-%m}: could not read PSI ({exc}) — "
                         f"keeping last known state.")


        # try:
        #     df_psi = simu.read_outputs("psi").copy()
        #     last_psi_t = df_psi.index[-1]
        #     last_month_t = t_month[-1]
        #
        #     if last_psi_t < last_month_t:
        #         print(f"  ERROR: PSI last index ({last_psi_t}) is behind "
        #               f"t_month[-1] ({last_month_t}) for {dt:%Y-%m} — "
        #               f"solver did not converge. STOPPING.")
        #         sys.exit(1)
        #
        #     print(f"  time index ic: {last_psi_t} ✓ (t_month[-1]: {last_month_t})")
        #     df_psi["datetime"] = cathy_utils.change_x2date(df_psi.index, dt)
        #
        #     psi_records.append(
        #         df_psi.reset_index(drop=True)
        #               .melt(id_vars="datetime", var_name="node", value_name="psi")
        #               .assign(node=lambda d: d["node"].astype(int))
        #     )
        #
        # except Exception as exc:
        #     print(f"  ERROR: could not read PSI for {dt:%Y-%m} ({exc}) — STOPPING.")
        #     sys.exit(1)


        # 12.9  SW
        try:
            df_sw, _ = simu.read_outputs("sw")
            df_sw = df_sw.copy()
            df_sw["datetime"] = cathy_utils.change_x2date(df_sw.index, dt)
            sw_month_df = (
                df_sw.reset_index(drop=True)
                     .melt(id_vars="datetime", var_name="node", value_name="sw")
                     .assign(node=lambda d: d["node"].astype(int))
            )
            sw_records.append(sw_month_df)
            month_stats["sw_mean"] = float(sw_month_df["sw"].mean())
            if args.save_monthly:
                try:
                    (
                        sw_month_df
                        .groupby(["datetime", "node"], as_index=False)
                        .last()
                        .set_index(["datetime", "node"])
                        .to_xarray()
                        .to_netcdf(monthly_dir / f"sw_{dt:%Y%m}.nc")
                    )
                except Exception as exc:
                    print(f"  Warning: could not save monthly SW for {dt:%Y-%m} ({exc})")
        except Exception as exc:
            print(f"  Warning: could not read SW for {dt:%Y-%m} ({exc})")

        # 12.10  ET
        try:
            df_et = simu.read_outputs("ET").copy()
            df_et = add_datetime_col(df_et, dt)

            # FIX: floor the net-forcing ETa artefact at ETp for this month
            # (adds 'ACT. ETRA_patched'; 'ACT. ETRA' stays unmodified).
            if fm.use_etp:
                df_et, n_eta_patched, n_eta_net_neg = apply_eta_artefact_floor(
                    df_et, t_month, ETp_nodes, rain_nodes, idx_month,
                    net_neg_mode=args.eta_net_neg_mode,
                )
                if n_eta_patched or n_eta_net_neg:
                    print(f"  ETa patch: {n_eta_patched} cells (net ≥ 0 artefact) + "
                          f"{n_eta_net_neg} cells (net < 0, mode={args.eta_net_neg_mode}) "
                          f"for {dt:%Y-%m}")
                md_log.bullet(f"{dt:%Y-%m}: ETa artefact floor patched "
                               f"{n_eta_patched} cell(s) (net rain ≥ ETp); "
                               f"{n_eta_net_neg} cell(s) set for net < 0 "
                               f"(mode={args.eta_net_neg_mode}).")
                md_log.add_metric("ETa artefact cells patched (total)", n_eta_patched)
                md_log.add_metric("ETa net<0 cells set (total)", n_eta_net_neg)

            et_xr = df_et.set_index(["X", "Y", "datetime"]).to_xarray()
            et_xr = et_xr.rio.set_spatial_dims(x_dim="X", y_dim="Y")
            et_xr = et_xr.rio.write_crs(TARGET_CRS)
            et_records.append(et_xr)
            if "ACT. ETRA" in df_et.columns:
                month_stats["eta_mean"] = float(df_et["ACT. ETRA"].mean())
            if "ACT. ETRA_patched" in df_et.columns:
                month_stats["eta_patched_mean"] = float(df_et["ACT. ETRA_patched"].mean())
            if args.save_monthly:
                try:
                    et_xr.to_netcdf(monthly_dir / f"et_{dt:%Y%m}.nc")
                except Exception as exc:
                    print(f"  Warning: could not save monthly ET for {dt:%Y-%m} ({exc})")
        except Exception as exc:
            print(f"  Warning: could not read ET for {dt:%Y-%m} ({exc})")

        # 12.11  Recharge
        try:
            df_rec = simu.read_outputs("recharge").copy()
            df_rec = add_datetime_col(df_rec, dt)
            xr_rec = df_rec.set_index(["X", "Y", "datetime"]).to_xarray()
            xr_rec = xr_rec.rio.set_spatial_dims(x_dim="X", y_dim="Y")
            xr_rec = xr_rec.rio.write_crs(TARGET_CRS)
            all_recharge_recs.append(xr_rec)
            rec_value_cols = [c for c in df_rec.select_dtypes(include=[np.number]).columns
                               if c not in ("X", "Y", "time_sec")]
            if rec_value_cols:
                month_stats["recharge_mean"] = float(df_rec[rec_value_cols].mean().mean())
            if args.save_monthly:
                try:
                    xr_rec.to_netcdf(monthly_dir / f"recharge_{dt:%Y%m}.nc")
                except Exception as exc:
                    print(f"  Warning: could not save monthly recharge for {dt:%Y-%m} ({exc})")
        except Exception as exc:
            print(f"  Warning: could not read recharge for {dt:%Y-%m} ({exc})")

        # 12.11a  Surface runoff hydrograph (hgraph) — discharge at the
        # outlet, one time series per monthly run (no X/Y spatial dims,
        # so it's accumulated/saved as CSV rather than NetCDF).
        try:
            df_hgraph = simu.read_outputs("hgraph").copy()
            df_hgraph = df_hgraph.rename(columns={"time": "time_sec"})
            df_hgraph = add_datetime_col(df_hgraph, dt)
            all_hgraph_recs.append(df_hgraph)
            if args.save_monthly:
                try:
                    df_hgraph.to_csv(monthly_dir / f"hgraph_{dt:%Y%m}.csv", index=False)
                except Exception as exc:
                    print(f"  Warning: could not save monthly hgraph for {dt:%Y-%m} ({exc})")
        except Exception as exc:
            print(f"  Warning: could not read hgraph for {dt:%Y-%m} ({exc})")

        # 12.11b  Monthly LAI raster backup (already kept in lai_history
        # for the combined lai_output.nc written by save_combined_outputs;
        # also dropped per-month here when --save-monthly is passed).
        if args.save_monthly:
            try:
                xr.DataArray(
                    lai_2d,
                    dims=["row", "col"],
                    coords={"row": np.arange(M), "col": np.arange(N)},
                    name="LAI",
                    attrs={"crs": TARGET_CRS},
                ).to_netcdf(monthly_dir / f"lai_{dt:%Y%m}.nc")
            except Exception as exc:
                print(f"  Warning: could not save monthly LAI for {dt:%Y-%m} ({exc})")

        # 12.11c  Sanity checks (physical bounds + approximate water balance)
        try:
            sanity_warnings = check_monthly_sanity(
                dt, month_stats,
                psi_month_df=psi_month_df,
                sw_month_df=sw_month_df,
                df_et=df_et,
                df_rec=df_rec,
                t_month=t_month,
                idx_month=idx_month,
                ETp_nodes=ETp_nodes if fm.use_etp else None,
                rain_nodes=rain_nodes if fm.use_etp else None,
            )
        except Exception as exc:
            sanity_warnings = []
            print(f"  Warning: sanity checks failed for {dt:%Y-%m} ({exc})")

        for w in sanity_warnings:
            print(f"  ⚠ sanity: {w}")
            md_log.warn(w)

        # 12.12  Progress summary
        semantic_zones = sorted(inv_remap.get(z, z) for z in active_zones)
        n_veg = int((veg_map > 1).sum())
        print(f"  {dt:%Y-%m-%d} | zones (semantic): {semantic_zones} | "
              f"vegetated px: {n_veg}/{N*M} ({100*n_veg/(N*M):.1f} %) | "
              f"LAI updated: {lai_updated}")

        # 12.12b  Monthly info stats (psi/sw/ETa/recharge means) — printed,
        # logged to the markdown run log, and appended to a single running
        # monthly_stats.csv in out_dir that's rewritten after every month.
        def _fmt(key: str, unit: str = "", sci: bool = False) -> str:
            v = month_stats.get(key)
            if v is None:
                return "n/a"
            return f"{v:.3e}{unit}" if sci else f"{v:.3f}{unit}"

        print(f"  stats     | psi mean: {_fmt('psi_mean', ' m')} "
              f"(min {_fmt('psi_min', ' m')}, max {_fmt('psi_max', ' m')}) | "
              f"sw mean: {_fmt('sw_mean')} | "
              f"ETa mean: {_fmt('eta_mean', ' m/s', sci=True)} | "
              f"recharge mean: {_fmt('recharge_mean', ' m/s', sci=True)}")
        print(f"  sanity    | SW out-of-bounds: {month_stats.get('sw_out_of_bounds_n', 'n/a')} | "
              f"PSI extreme: {month_stats.get('psi_extreme_n', 'n/a')} | "
              f"ETa>ETp: {month_stats.get('eta_exceeds_etp_n', 'n/a')} | "
              f"water balance residual: {_fmt('water_balance_residual_mm', ' mm')}")
        md_log.bullet(
            f"{dt:%Y-%m} stats — psi mean: {_fmt('psi_mean', ' m')}, "
            f"sw mean: {_fmt('sw_mean')}, "
            f"ETa mean: {_fmt('eta_mean', ' m/s', sci=True)}, "
            f"recharge mean: {_fmt('recharge_mean', ' m/s', sci=True)}"
        )
        month_stats["solver_status"] = solver_status
        stats_records.append(month_stats)
        try:
            pd.DataFrame(stats_records).to_csv(out_dir / "monthly_stats.csv", index=False)
        except Exception as exc:
            print(f"  Warning: could not refresh monthly_stats.csv after {dt:%Y-%m} ({exc})")

        # FIX: this counter was declared but never incremented before —
        # completed months now register in the run summary written by
        # md_log.close().
        n_months_run += 1

        # 12.13  Refresh the single running combined-output files with
        # everything completed so far (in addition to this month's
        # standalone files already saved above).
        try:
            save_combined_outputs(
                out_dir, veg_map_history, lai_history, et_records,
                psi_records, sw_records, all_recharge_recs,
                all_hgraph_recs, M, N,
            )
        except Exception as exc:
            print(f"  Warning: could not refresh combined outputs after {dt:%Y-%m} ({exc})")

    print("\n── Loop complete (partial runs saved if any months were skipped) ────────")

    # ── 13/14. Final refresh of the combined outputs (veg history, PSI,
    # SW, ET, recharge) — same helper called after every month in the
    # loop above, so this call just guarantees the very last month is
    # reflected too.
    save_combined_outputs(
        out_dir, veg_map_history, lai_history, et_records,
        psi_records, sw_records, all_recharge_recs,
        all_hgraph_recs, M, N,
    )
    print(f"veg_map_history saved → {out_dir / 'veg_map_history.pkl'}")
    if et_records:
        print(f"ET saved → {out_dir / 'et_output.nc'}")
    if psi_records:
        print(f"PSI saved → {out_dir / 'psi_output.nc'}")
    if sw_records:
        print(f"SW saved → {out_dir / 'sw_output.nc'}")
    if all_recharge_recs:
        print(f"Recharge saved → {out_dir / 'recharge_output.nc'}")
    if lai_history:
        print(f"LAI saved → {out_dir / 'lai_output.nc'}")
    if all_hgraph_recs:
        print(f"Hydrograph saved → {out_dir / 'hgraph_output.csv'}")

    print("\n" + "=" * 60)
    print("  SIMULATION COMPLETE (partial if any months were skipped)")
    print(f"  Project   : {simu.project_name}")
    print(f"  Workdir   : {simu.workdir}")
    print(f"  Outputs   : {out_dir}")
    print("=" * 60)

    # FIX: close() was previously never called — the markdown log had no
    # run-summary section, no elapsed time, and no final footer. Now every
    # run's simulation_log.md ends with months run/skipped, accumulated
    # metrics (e.g. total ETa artefact cells patched), warning count, and
    # elapsed wall time.
    md_log.close(summary={
        "Months run": n_months_run,
        "Months skipped": n_months_skipped,
        "Output directory": str(out_dir),
    })


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    args = _parse_args()
    run_pipeline(args)
