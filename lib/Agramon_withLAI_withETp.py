"""
LT_SShydro_run_withLAI_withETp.py — Long-Term Subsurface Hydrology Simulation
              with LAI + satellite ETp/rain forcing (CATHY)
================================================================================
GRwater project — Agramon study site

Runs the full pyCATHY pipeline with monthly LAI-driven vegetation updates,
forcing the ATMBC boundary from satellite ETp and rain (not ERA5).

  1.  Load satellite ETp + rain NetCDFs (EPSG:32630 → EPSG:25830)
  2.  Load monthly LAI dataset
  3.  Log simulation parameters
  4.  Load and pre-process the DEM
  5.  Initialise and mesh the CATHY model
  6.  Map ETp and rain grids onto mesh nodes; build full ATMBC arrays
  7.  Configure soil, ICs, BCs, and output settings
  8.  For each calendar month: update veg map from LAI → run processor
  9.  Aggregate and save outputs (NetCDF)

Usage
-----
    python Agramon_withLAI_withETp_fixed.py
    python Agramon_withLAI_withETp_fixed.py --start-year 2019 --end-year 2022
    python Agramon_withLAI_withETp_fixed.py --zroot 0.5 --pmin -10
    python Agramon_withLAI_withETp_fixed.py --dem-plot 4
    python Agramon_withLAI_withETp_fixed.py --no-plots
    python Agramon_withLAI_withETp_fixed.py --zroot-factors 0.1 0.3 1.0 1.5 2.5
    python Agramon_withLAI_withETp_fixed.py --permx 1e-5 --permy-ratio 1.0 --permz-ratio 0.1

Forcing switch
--------------
    ERA5  : set USE_ETP_FORCING = False  (uses AgUtils.extract_point_timeseries)
    ETp   : set USE_ETP_FORCING = True   (default — satellite ETp + rain)
"""

import argparse
import os
import pickle
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import rioxarray as rxr  # noqa: F401  (activates .rio accessor)
from scipy.ndimage import binary_dilation, label as sp_label, gaussian_filter
from scipy.stats import mode
from pyCATHY import meshtools as mt


# ── Forcing switch ─────────────────────────────────────────────────────────────
# True  → use satellite ETp + rain (PathET files)
# False → use ERA5 (original withLAI behaviour)
USE_ETP_FORCING = True


# ── Path setup ────────────────────────────────────────────────────────────────

MODULE_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
    #"/home/ben/Nextcloud/BenCSIC/Codes"
    "/Tech4agro_org/GRwater_geophy"
).resolve()

EO_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
    #"/home/ben/Nextcloud/BenCSIC/Codes"
    "/Tech4agro_org/GRwater_CATmentHYdrology"
).resolve()

LAI_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Training_Supervision"
    #"/home/ben/Nextcloud/BenCSIC/Training_Supervision"
    "/Supervision/Xela_Carracedo_Practicas_2026_data/301a-biophysical/agramon/input"
).resolve()

PathET    = Path("/home/z0272571a@campus.csic.es/Nextcloud/GRwater/data/satellite/Agramon")
#PathET    = Path("/home/ben/Nextcloud/GRwater/data/satellite/Agramon")
etp_path  = PathET / "20161001_20241231_ET0.nc"
rain_path = PathET / "20161001_20241231_TP.nc"

if str(MODULE_PATH) not in sys.path:
    sys.path.append(str(MODULE_PATH))

import Agramon_utils as AgUtils
from geoutils import uCATHY
from pyCATHY import cathy_tools, cathy_utils


# ── Constants ─────────────────────────────────────────────────────────────────

# All spatial datasets are reprojected to this CRS before any operation.
# DEM / shapefiles : EPSG:25830 (UTM zone 30N, ETRS89)
# LAI / ETp / rain : EPSG:32630 (UTM zone 30N, WGS84) — same zone, different datum
TARGET_CRS = "EPSG:25830"
NATIVE_CRS = "EPSG:32630"

# LAI class boundaries  [bare | sparse | moderate | dense | very dense]
# Calibrated to the observed Agramon LAI range (0 – 0.6 m² m⁻²).
LAI_THRESHOLDS = [0.1, 0.25, 0.4, 0.55]

# Fixed number of semantic vegetation classes (bare/sparse/moderate/dense/very
# dense). This must stay constant for the whole run and match MAXVEG passed
# to simu.update_cathyH(), since pyCATHY's SOIL vegetation table has a fixed
# row count for the whole project — it is NOT allowed to vary month to month
# with however many classes happen to be spatially present.
N_VEG_CLASSES = len(LAI_THRESHOLDS) + 1  # = 5

# Vegetation zone → root-depth multiplier
ZROOT_FACTORS = {
    1: 0.05,   # bare soil
    2: 0.25,   # sparse vegetation
    3: 1.0,    # moderate vegetation
    4: 1.5,    # dense vegetation
    5: 2.0,    # very dense vegetation
}


def zroot_factors_from_args(args: argparse.Namespace) -> dict[int, float]:
    """
    Build the semantic-zone → root-depth-multiplier mapping consumed by
    ``update_soil_for_active_zones`` from ``--zroot-factors`` (5 values,
    positional order: bare, sparse, moderate, dense, very dense —
    matching ``ZROOT_FACTORS``' key order 1-5). Centralised here so the
    spin-up and main-loop call sites build the exact same mapping from
    the exact same CLI values instead of risking drift between the two.
    """
    return {i + 1: v for i, v in enumerate(args.zroot_factors)}

# ── Sanity-check thresholds (monthly loop) ─────────────────────────────────
# Tolerances for the physical-bounds and water-balance checks run after
# every month (see `check_monthly_sanity`). Not physical hard limits —
# just flags for "this deserves a look", so kept generous on purpose.
SW_TOL           = 0.01     # SW allowed slightly outside [0, 1] before flagging
PSI_SANITY_MIN   = -50.0    # m — PSI below this is flagged as extreme drying
ETA_TOL          = 1e-8     # m/s — ETa allowed slightly above ETp before flagging
WB_RESIDUAL_FRAC = 0.5      # flag if |water balance residual| > this fraction of rain_in


# ── Markdown run log ─────────────────────────────────────────────────────────

class MarkdownLog:
    """
    Lightweight, append-as-you-go markdown logger for a single simulation
    run — one narrative file mirroring the console output (run config,
    per-step progress, spin-up diagnostics, and a per-month LAI/ETp/solver
    table), so a run can be reviewed later without re-parsing stdout.

    Writes to disk on every call (re-opened in append mode each time)
    rather than buffering in memory, so a crash mid-run still leaves a
    readable partial log instead of losing everything unwritten.
    """

    def __init__(self, path: Path, title: str) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._n_warnings = 0
        self._start_time = datetime.now()
        self._metrics: dict[str, float] = {}
        try:
            with open(self.path, "w") as fh:
                fh.write(f"# {title}\n\n")
                fh.write(f"_Started: {datetime.now():%Y-%m-%d %H:%M:%S}_\n\n")
        except OSError as exc:
            # Logging failures should never abort the simulation itself.
            print(f"  WARNING: could not open markdown log at {self.path} ({exc}) "
                  f"— continuing without it.")

    def _write(self, text: str) -> None:
        try:
            with open(self.path, "a") as fh:
                fh.write(text)
        except OSError as exc:
            print(f"  WARNING: markdown log write failed ({exc}).")

    def section(self, title: str) -> None:
        self._write(f"\n## {title}\n\n")

    def subsection(self, title: str) -> None:
        self._write(f"\n### {title}\n\n")

    def meta(self, items: dict) -> None:
        """Key/value bullet list — run parameters, step summaries, etc."""
        for k, v in items.items():
            self._write(f"- **{k}**: {v}\n")

    def line(self, text: str = "") -> None:
        self._write(f"{text}\n")

    def bullet(self, text: str) -> None:
        self._write(f"- {text}\n")

    def warn(self, text: str) -> None:
        self._n_warnings += 1
        self._write(f"- ⚠️ **WARNING**: {text}\n")

    def table_header(self, headers: list) -> None:
        self._write("\n| " + " | ".join(headers) + " |\n")
        self._write("|" + "|".join(["---"] * len(headers)) + "|\n")

    def table_row(self, cells: list) -> None:
        self._write("| " + " | ".join(str(c) for c in cells) + " |\n")

    def fixes_applied(self, fixes: list) -> None:
        """
        Document known pipeline patches active for this run — written once,
        near the top of the log, so a run's markdown self-describes which
        fixes were in effect without needing to diff the script against git.
        """
        self.section("Known fixes applied")
        for fx in fixes:
            self.bullet(fx)

    def add_metric(self, name: str, value: float) -> None:
        """
        Accumulate a running numeric total across the whole run (e.g. cells
        patched per month) without the caller needing to track it manually.
        Automatically rolled into the summary written by close().
        """
        self._metrics[name] = self._metrics.get(name, 0) + value

    def close(self, summary: dict | None = None) -> None:
        elapsed = datetime.now() - self._start_time
        self.section("Run summary")
        merged = dict(summary) if summary else {}
        merged.update(self._metrics)
        merged["Warnings logged"] = self._n_warnings
        merged["Elapsed"] = str(elapsed).split(".")[0]  # drop microseconds
        self.meta(merged)
        self._write(f"\n---\n\n_Finished: {datetime.now():%Y-%m-%d %H:%M:%S}_\n")


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
    sim.add_argument("--pmin",   type=float, default=-1e99, metavar="M",
                     help="Minimum pressure head threshold (m)")
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

    dem = p.add_argument_group("DEM selection")
    dem.add_argument("--dem-tif", default=str(MODULE_PATH / "DTMplots/20250618_AGRAMON100m_micasense_dtm.tif"),
                     metavar="PATH",
                     help="Full path to the GeoTIFF DEM used to build the mesh "
                          "(default: 20250618_AGRAMON100m_micasense_dtm.tif in MODULE_PATH)")
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
    out.add_argument("--save-monthly", action="store_true",
                     help="Also write a standalone NetCDF file per month per "
                          "variable to outputs/<scenario>/monthly/ (psi, sw, "
                          "et, recharge). Off by default — the single running "
                          "psi_output.nc / sw_output.nc / et_output.nc / "
                          "recharge_output.nc files in outputs/<scenario>/ "
                          "are always refreshed after every month regardless.")

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
                                               MODULE_PATH / "shapefiles")

    if gdf.crs is None:
        raise ValueError("Shapefile has no CRS — assign one before proceeding.")
    if str(gdf.crs) != TARGET_CRS:
        print(f"  Reprojecting shapefiles {gdf.crs} → {TARGET_CRS}")
        gdf = gdf.to_crs(TARGET_CRS)
    print(f"  Shapefile CRS : {gdf.crs}")
    return gdf
def load_dem_from_tif(
    tif_path: Path, target_crs: str = TARGET_CRS,
    gdf_clip: gpd.GeoDataFrame | None = None,
    fid: int | str | None = None,
    plot_id: str | None = None,
    gbuffer: float | None = None,
    show: bool = False,
    resample_resolution: float | None = None,
):
    """
    Load a DEM directly from a GeoTIFF, independent of ``AgUtils.load_dem``
    (which expects a folder of .adf grids). Mirrors that function's
    conventions exactly so it's a drop-in swap at the call site:

      - reprojects to ``target_crs`` (assumes EPSG:25830 if the GeoTIFF
        has no CRS embedded, same fallback as ``AgUtils.load_dem``);
      - by default, does NOT clip to any shapefile — it only crops to
        the DEM's own valid-data bounding box (trims all-NaN rows/cols
        at the edges), exactly like ``AgUtils.load_dem`` does (its
        shapefile-clip code path is present but disabled/dead there
        too). Pass ``gdf_clip`` to additionally clip to a shape's
        bounding box + buffer if the new DEM's extent turns out to be
        much larger than the catchment of interest;
      - fills nodata with ``-9999`` (CATHY's expected DEM no-data
        sentinel), NOT NaN — matching ``AgUtils.load_dem`` bit for bit,
        since the pyCATHY preprocessor/mesh writer expects this value;
      - sets ``xllcorner``/``yllcorner`` to the raw min(x)/min(y) pixel
        centre coordinates (no half-pixel edge offset), same as
        ``AgUtils.load_dem``.

    Parameters
    ----------
    tif_path : Path
        Full path to the GeoTIFF DEM.
    target_crs : str
        CRS every spatial dataset in the pipeline is aligned to.
    gdf_clip : gpd.GeoDataFrame, optional
        GeoDataFrame whose bounding box is used to clip the DEM.
        If provided, the DEM is clipped to this bounding box + buffer.
    fid : int or str, optional
        Feature ID to select a single polygon from ``gdf_clip`` by its
        ``fid_1`` column. Mutually exclusive with ``plot_id``.
    plot_id : str, optional
        Plot ID to select polygons from ``gdf_clip`` by its ``PlotID``
        column (e.g., ``args.dem_plot``). Mutually exclusive with ``fid``.
    gbuffer : float, optional
        Padding in metres around ``gdf_clip``'s bbox. Defaults to 1
        pixel of the current raster resolution.
    show : bool
        If True, plot the masked DEM for a quick visual sanity check.
    resample_resolution : float, optional
        Target resolution in metres for resampling (assumes square
        pixels). If provided, the DEM is resampled to this resolution
        after CRS transformation but before cropping to valid data.

    Returns
    -------
    raster_DEM : xr.DataArray
        DEM reprojected/cropped to its own valid-data extent, with
        ``.rio.crs`` set. Same x/y grid as ``raster_DEM_masked``.
    raster_DEM_masked : np.ndarray
        DEM elevations, with NaN → -9999 (CATHY no-data convention).
    xllcorner, yllcorner : float
        Lower-left corner of the grid (pixel-centre coordinate, not
        cell edge — matches ``AgUtils.load_dem``).
    res_x, res_y : float
        Pixel size (positive values).
    """
    tif_path = Path(tif_path)
    if not tif_path.is_file():
        raise FileNotFoundError(f"DEM GeoTIFF not found: {tif_path}")

    print(f"  DEM source (tif) : {tif_path}")
    raster_DEM = rxr.open_rasterio(tif_path, masked=True).isel(band=0)

    if raster_DEM.rio.crs is None:
        raster_DEM = raster_DEM.rio.write_crs("EPSG:25830")
        print("  DEM CRS was missing — assumed EPSG:25830")

    # Handle CRS transformation and resampling together
    needs_reproject = str(raster_DEM.rio.crs) != target_crs
    needs_resample = resample_resolution is not None

    if needs_reproject or needs_resample:
        dst_crs = target_crs if needs_reproject else raster_DEM.rio.crs
        kwargs = {"resolution": resample_resolution} if needs_resample else {}
        print(f"  Transforming DEM to {dst_crs} at {resample_resolution} m" if needs_resample
              else f"  Transforming DEM to {dst_crs}")
        raster_DEM = raster_DEM.rio.reproject(dst_crs, **kwargs)

    print(f"  DEM CRS : {raster_DEM.rio.crs}")

    res_x = abs(float(raster_DEM.rio.resolution()[0]))
    res_y = abs(float(raster_DEM.rio.resolution()[1]))
    print(f"  DEM resolution : {res_x:.2f} × {res_y:.2f} m")

    # Select and clip to gdf_clip bounding box + buffer if provided
    if gdf_clip is not None:
        # Select specific feature(s) from gdf_clip if fid or plot_id is provided
        if fid is not None and plot_id is not None:
            raise ValueError("Cannot specify both 'fid' and 'plot_id'. Use one or the other.")
        if fid is not None:
            gdf_clip = gdf_clip[gdf_clip["fid_1"] == fid].copy()
            print(f"  Selected fid={fid} from gdf_clip ({len(gdf_clip)} features)")
        elif plot_id is not None:
            gdf_clip = gdf_clip[gdf_clip["PlotID"] == plot_id].copy()
            print(f"  Selected PlotID={plot_id} from gdf_clip ({len(gdf_clip)} features)")

        if len(gdf_clip) == 0:
            raise ValueError(f"No features found in gdf_clip for {'fid' if fid is not None else 'PlotID'}={fid if fid is not None else plot_id}")

        if str(gdf_clip.crs) != target_crs:
            gdf_clip = gdf_clip.to_crs(target_crs)
            print(f"  Reprojected gdf_clip to {target_crs}")

        if gbuffer is None:
            gbuffer = res_x  # Default: 1 pixel of current resolution
        minx, miny, maxx, maxy = gdf_clip.total_bounds
        minx -= gbuffer
        miny -= gbuffer
        maxx += gbuffer
        maxy += gbuffer
        print(f"  Clipping DEM to gdf_clip bbox + {gbuffer:.2f} m buffer")
        #raster_DEM = raster_DEM.rio.clip_box(minx, miny, maxx, maxy)
        raster_DEM = raster_DEM.rio.clip(gdf_clip.geometry, all_touched=True)

    # Crop to valid-data bounding box (trim all-NaN rows/cols)
    valid_mask = ~np.isnan(raster_DEM)
    valid_x = raster_DEM["x"].where(valid_mask.any(dim="y"), drop=True)
    valid_y = raster_DEM["y"].where(valid_mask.any(dim="x"), drop=True)
    raster_DEM = (
        raster_DEM
        .where((raster_DEM["x"] >= float(valid_x.min())) &
               (raster_DEM["x"] <= float(valid_x.max())), drop=True)
        .where((raster_DEM["y"] >= float(valid_y.min())) &
               (raster_DEM["y"] <= float(valid_y.max())), drop=True)
    )

    xllcorner = float(raster_DEM["x"].min())
    yllcorner = float(raster_DEM["y"].min())

    raster_DEM_masked = np.where(np.isnan(raster_DEM), -9999, raster_DEM.values)

    print(f"  DEM shape      : {raster_DEM_masked.shape}")
    print(f"  xllcorner={xllcorner:.2f}, yllcorner={yllcorner:.2f}")

    if show:
        valid = raster_DEM_masked[raster_DEM_masked != -9999]
        fig, ax = plt.subplots(figsize=(8, 6))
        img = ax.imshow(raster_DEM_masked, cmap="terrain",
                        vmin=valid.min(), vmax=valid.max())
        plt.colorbar(img, ax=ax, label="Elevation (m)")
        ax.set_title(f"DEM — masked ({target_crs})")
        plt.tight_layout()

    return raster_DEM, raster_DEM_masked, xllcorner, yllcorner, res_x, res_y


def load_monthly_lai(lai_input_dir: Path) -> xr.Dataset:
    """
    Load monthly LAI dataset, assign its native CRS (EPSG:32630),
    and reproject to TARGET_CRS.
    """
    ds_lai = xr.open_dataset(lai_input_dir / "lai_monthly.nc")

    if ds_lai.rio.crs is None:
        ds_lai = ds_lai.rio.write_crs("EPSG:32630")
        print("  LAI CRS assigned : EPSG:32630")
    print(f"  LAI native CRS   : {ds_lai.rio.crs}")

    if str(ds_lai.rio.crs) != TARGET_CRS:
        print(f"  Reprojecting LAI {ds_lai.rio.crs} → {TARGET_CRS}")
        ds_lai = ds_lai.rio.reproject(TARGET_CRS)
    print(f"  LAI final CRS    : {ds_lai.rio.crs}")
    print(f"  LAI x range      : {float(ds_lai.x.min()):.1f} – {float(ds_lai.x.max()):.1f}")
    print(f"  LAI y range      : {float(ds_lai.y.min()):.1f} – {float(ds_lai.y.max()):.1f}")
    return ds_lai


def load_era5(start_year: int, end_year: int) -> xr.Dataset:
    """Load ERA5 climate forcing."""
    ds = AgUtils.load_era5_series(root_path=EO_PATH, start_year=start_year, end_year=end_year)
    if ds is not None:
        print(f"  Loaded ERA5 : {start_year} – {end_year - 1}")
    else:
        print("  No ERA5 datasets loaded — check EO_PATH.")
    return ds


def _finalise_ds(ds: xr.Dataset, label: str) -> xr.Dataset:
    """Write CRS + spatial dims after load or reproject."""
    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  {label} final CRS : {ds.rio.crs}")
    print(f"  {label} x range   : {float(ds.x.min()):.1f} – {float(ds.x.max()):.1f}")
    print(f"  {label} y range   : {float(ds.y.min()):.1f} – {float(ds.y.max()):.1f}")
    return ds


def _load_nc(path: Path, label: str) -> xr.Dataset:
    """Open NetCDF, assign native CRS (EPSG:32630), reproject to TARGET_CRS."""
    ds = (
        xr.open_dataset(path, decode_times=True, mask_and_scale=True)
          .drop_vars("spatial_ref", errors="ignore")
    )
    ds = ds.rio.write_crs(NATIVE_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  {label} native CRS : {ds.rio.crs}")

    if str(ds.rio.crs) != TARGET_CRS:
        print(f"  Reprojecting {label}  {ds.rio.crs} → {TARGET_CRS}")
        ds = ds.rio.reproject(TARGET_CRS)

    return _finalise_ds(ds, label)


def load_etp(path: Path = etp_path) -> xr.Dataset:
    """Load ETp dataset and reproject to TARGET_CRS."""
    return _load_nc(path, "ETp")


def load_rain(path: Path = rain_path) -> xr.Dataset:
    """Load rain dataset and reproject to TARGET_CRS."""
    return _load_nc(path, "rain")


def log_simulation(log_path: Path, start_year: int, end_year: int,
                   zroot: float, pmin: float, watershed_nb: int,
                   with_lai: int, outlet_flux: float = None,
                   bottom_flux: float = None, spinup: str = None,
                   ic_mode: str = None, ic_wtposition: float = None,
                   ic_pressure_head: float = None,
                   zroot_factors: list = None,
                   permx: float = None,
                   permy_ratio: float = None,
                   permz_ratio: float = None) -> int:
    """Register parameters in the shared CSV log and return the sim_index.

    outlet_flux / bottom_flux / spinup / ic_* are recorded too (not just
    used to tag the folder name) so `--list-scenarios` and the CSV itself
    are enough to compare runs without having to parse folder names — and
    so batch runners can tell whether a given parameter combination
    (including initial condition) has already been logged, e.g. to skip
    re-running scenarios that already completed successfully.

    zroot_factors / permx / permy_ratio / permz_ratio are logged the same
    way — a run's root-depth mapping and Ks override (incl. anisotropy)
    are as much a part of its identity as ZROOT or PMIN, so they belong
    in the CSV rather than only ever showing up in the console/markdown
    log.
    """
    params = {
        "start_year":   start_year,
        "end_year":     end_year,
        "ZROOT":        zroot,
        "PMIN":         pmin,
        "watershed_nb": watershed_nb,
        "withLAI":      with_lai,
        "outlet_flux":  outlet_flux,
        "bottom_flux":  bottom_flux,
        "spinup":       spinup,
        "ic_mode":          ic_mode,
        "ic_wtposition":    ic_wtposition,
        "ic_pressure_head": ic_pressure_head,
        "zroot_factors": (",".join(f"{v:g}" for v in zroot_factors)
                           if zroot_factors is not None else None),
        "permx":        permx,
        "permy_ratio":  permy_ratio if permx is not None else None,
        "permz_ratio":  permz_ratio if permx is not None else None,
    }
    sim_index = uCATHY.log_simulation(log_path, params)
    print(f"  sim_index = {sim_index}  (withLAI={with_lai})")
    return sim_index


def build_scenario_dirname(sim_index: int, dem_plot: int,
                           outlet_flux: float, bottom_flux: float,
                           no_spinup: bool, spinup_cycles: int,
                           pmin: float, ic_mode: str = "wt") -> str:
    """
    Human-readable, self-describing scenario folder name.

    ``scenario_<sim_index>`` stays the root/prefix of the name — every
    place that keys off the CSV-assigned sim_index (post-processing
    scripts, --scenario N, log lookups) keeps working unchanged — but a
    short tag is appended describing the axes that actually differ
    between batch runs (dem-plot, outlet/bottom BC, spin-up, PMIN,
    initial condition), so scenarios run with different args can be told
    apart, and compared, just by reading folder names side by side
    instead of cross-referencing the CSV log for every one of them.

    ``ic_mode`` is included so runs that only differ by --ic-mode (wt vs
    pressure) get distinct folders/tags instead of silently colliding —
    previously two such runs shared the same tag (differing only by the
    sim_index prefix), which also made them indistinguishable to any
    "has this combination already run?" check based on the tag alone.

    e.g. scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c_pmin-3.0_ic-wt
    """
    outlet_label = "closed" if outlet_flux == 0 else "free"
    bottom_label = "closed" if bottom_flux == 0 else "free"
    spinup_label = "off" if no_spinup else f"on{spinup_cycles}c"
    pmin_label = f"{pmin:g}".replace("-", "m").replace("+", "").replace(".", "p")
    tag = (f"plot{dem_plot}"
           f"_outlet-{outlet_label}"
           f"_bottom-{bottom_label}"
           f"_spinup-{spinup_label}"
           f"_pmin-{pmin_label}"
           f"_ic-{ic_mode}")
    return f"scenario_{sim_index}_{tag}"


def initialise_project(path2prj: str, start_year: int, end_year: int,
                       zroot: float, pmin: float, scenario_dirname: str):
    """Create and return the CATHY project object inside a sim-specific folder."""
    sim_dir = Path(path2prj) / scenario_dirname
    sim_dir.mkdir(parents=True, exist_ok=True)
    simu = cathy_tools.CATHY(dirName=str(sim_dir))
    print(f"  Working directory  : {sim_dir}")
    return simu


def smooth_dem(raster_DEM_masked: np.ndarray, sigma: float) -> np.ndarray:
    """
    Gaussian-smooth DEM elevations while leaving array shape and the
    NaN (masked/outside-catchment) footprint exactly unchanged.

    Plain ``scipy.ndimage.gaussian_filter`` on an array containing NaNs
    would leak NaN into every pixel within ~sigma of the mask edge (and
    vice-versa smear elevation into masked cells), silently growing or
    shrinking the effective valid-pixel footprint. That would change
    what "the DEM" is even though .shape stays (N, M) — the grid dims
    CATHY sees would still be right, but the terrain wouldn't be.

    Standard fix: zero-fill NaNs, smooth both the filled field and a
    0/1 validity mask with the same kernel, then divide the two
    (Gaussian-weighted normalisation) so only valid neighbours
    contribute to each valid pixel's new value. NaNs are then written
    back at the original locations, so the mask is bit-for-bit
    identical to the input — only masked-in elevations move.

    Parameters
    ----------
    raster_DEM_masked : np.ndarray
        DEM elevation grid, NaN outside the catchment/basin.
    sigma : float
        Gaussian sigma in pixels. <= 0 returns the input unchanged
        (no-op — smoothing was requested but effectively disabled).

    Returns
    -------
    np.ndarray
        Smoothed DEM, same shape and same NaN mask as the input.
    """
    if sigma is None or sigma <= 0:
        return raster_DEM_masked

    dem = np.asarray(raster_DEM_masked, dtype=float)
    nan_mask = ~np.isfinite(dem)

    filled = np.where(nan_mask, 0.0, dem)
    valid  = np.where(nan_mask, 0.0, 1.0)

    filled_smooth = gaussian_filter(filled, sigma=sigma, mode="nearest")
    valid_smooth  = gaussian_filter(valid,  sigma=sigma, mode="nearest")

    with np.errstate(invalid="ignore", divide="ignore"):
        dem_smooth = filled_smooth / valid_smooth

    # Re-impose the exact original mask — shape and valid-pixel count
    # (and therefore N, M and dem_valid_px used for area/logging) are
    # unchanged; only elevations at already-valid pixels move.
    dem_smooth[nan_mask] = np.nan
    return dem_smooth.astype(raster_DEM_masked.dtype, copy=False)


def setup_dem(simu, raster_DEM_masked: np.ndarray,
              xllcorner: float, yllcorner: float,
              delta_x: float, delta_y: float,
              show: bool = False) -> None:
    """Feed DEM into CATHY preprocessing inputs and optionally show 3-D terrain."""
    simu.update_prepo_inputs(
        xllcorner=xllcorner,
        yllcorner=yllcorner,
        DEM=raster_DEM_masked,
        delta_x=delta_x,
        delta_y=delta_y,
        ivert=1,
        #base=30,
    )
    if show:
        fig = plt.figure(figsize=(10, 6))
        ax  = plt.axes(projection="3d")
        simu.show_input(prop="dem", ax=ax)
        plt.tight_layout()
    simu.create_mesh_vtk(verbose=True)

def configure_boundary_conditions(
    simu, t_atmbc: list, outlet_flux: float = -1e-6, outlet_side: str = None,
    bottom_flux: float = 0.0, res_x: float = None, res_y: float = None,
) -> tuple[str, int]:
    """
    No-flow on the bottom and on 3 of the 4 lateral sides; a small
    prescribed outflow on the topographically lowest lateral side (the
    catchment outlet), so ponded/perched water can actually leave the
    domain instead of piling up against a fully closed box. Optionally
    also a free-bottom-drainage-style outflow on the mesh bottom.

    Thin orchestration layer only — the actual geometry (face areas),
    node identification, and file writing now live in pyCATHY itself
    (``CATHY.get_outlet_side``, ``CATHY.get_outlet_node``,
    ``CATHY.compute_nodal_face_areas``, ``CATHY.compute_bottom_face_area``,
    ``CATHY.update_nansfneubc``, ``CATHY.update_nansfdirbc``,
    ``CATHY.update_sfbc``, ``CATHY.check_neumann_flux_vs_ks``) rather than
    being reimplemented per-experiment here, so any other pyCATHY script
    needing an "open outlet" BC gets the same, already-reviewed logic
    instead of a fresh copy that can drift out of sync / re-introduce
    bugs already fixed once.

    ``outlet_flux`` / ``bottom_flux`` are in m/s, negative = water
    leaving the domain — tunable parameters, not physical constants. A
    reasonable starting point is a fraction of the saturated hydraulic
    conductivity (Ks) of the relevant soil zone; ``check_neumann_flux_vs_ks``
    below warns (does not block) if the magnitude exceeds Ks.

    UNITS NOTE: ``nansfneubc`` needs a volumetric flow rate Q [m^3/s]
    per node, not a flux density [m/s] like ``atmbc``. This is converted
    internally using each node's *own* face area
    (``compute_nodal_face_areas`` / ``compute_bottom_face_area``) —
    exact per node, not a side-wide average.

    CAVEAT (physical, not a code bug): this applies the same flux
    density to *every* node on the resolved side, from the DEM surface
    down to the mesh bottom — including nodes well above any realistic
    water table. That's an "always-open tap" approximation, not a true
    saturation-aware seepage face (pyCATHY's ``update_sfbc`` doesn't
    support non-trivial seepage faces yet — see its docstring). Treat
    behaviour near the shallow outlet-side nodes with more scepticism
    than behaviour near the water table / mesh bottom. If your
    catchment's real outlet is an interior pour point rather than a
    domain edge, treating a whole mesh side as the outlet is itself an
    approximation.

    Returns
    -------
    (outlet_side, outlet_node_id)
    """
    simu.create_mesh_bounds_df("nansfneubc", simu.grid3d["mesh3d_nodes"], t_atmbc)

    resolved_side = simu.get_outlet_side(outlet_side=outlet_side)
    outlet_node_id = simu.get_outlet_node()

    # Dirichlet & seepage-face channels stay fully closed for every time
    # block in t_atmbc (an explicit "no nodes imposed" record per block,
    # not an empty/0-byte file — see CATHY.update_nansfdirbc /
    # CATHY.update_sfbc docstrings for why time=[] used to be invalid).
    simu.update_nansfdirbc(time=t_atmbc, no_flow=True)
    simu.update_sfbc(time=t_atmbc, no_flow=True)

    # ── Neumann: per-node exact face-area conversion (m/s -> m^3/s) ──
    outlet_areas = simu.compute_nodal_face_areas(resolved_side)
    outlet_node_ids = np.array(sorted(outlet_areas.keys()), dtype=int)
    outlet_flux_m3s = np.array(
        [outlet_flux * outlet_areas[n] for n in outlet_node_ids]
    )
    simu.check_neumann_flux_vs_ks(outlet_flux, label="outlet_flux")

    bottom_active = bottom_flux != 0.0
    bottom_node_ids = np.array([], dtype=int)
    bottom_flux_m3s = np.array([])
    n_bottom_nodes = 0
    if bottom_active:
        if res_x is None or res_y is None:
            raise ValueError(
                "bottom_flux is non-zero but res_x/res_y (DEM pixel size) "
                "were not provided - needed to convert m/s -> m^3/s per "
                "bottom node."
            )
        df0 = simu.mesh_bound_cond_df[
            simu.mesh_bound_cond_df["time"] == simu.mesh_bound_cond_df["time"].iloc[0]
        ]
        bottom_node_ids = df0.loc[df0["bot_bound"], "id_node"].to_numpy(dtype=int)
        n_bottom_nodes = len(bottom_node_ids)
        bottom_nodal_area = simu.compute_bottom_face_area(res_x, res_y)
        bottom_flux_m3s = np.full(n_bottom_nodes, bottom_flux * bottom_nodal_area)
        simu.check_neumann_flux_vs_ks(bottom_flux, label="bottom_flux")

    # Combine both groups into ONE nansfneubc call. update_nansfneubc
    # rewrites the whole file on every call ("w+"), so calling it twice
    # (once per group) would silently lose the first group - merge
    # node/flux arrays here instead, exactly once.
    overlap_mask = np.isin(outlet_node_ids, bottom_node_ids)
    if overlap_mask.any():
        overlap_nodes = outlet_node_ids[overlap_mask]
        print(f"  WARNING: {len(overlap_nodes)} node(s) belong to both the "
              f"outlet side and the mesh bottom (e.g. a bottom corner also "
              f"on the lateral outlet edge): {sorted(overlap_nodes.tolist())}. "
              f"bottom_flux wins for these nodes (outlet contribution "
              f"dropped for them) since it's applied after the outlet "
              f"group below.")
    keep_outlet = ~overlap_mask
    all_nodes = np.concatenate([outlet_node_ids[keep_outlet], bottom_node_ids])
    all_flux = np.concatenate([outlet_flux_m3s[keep_outlet], bottom_flux_m3s])

    # Drop exactly-zero-flux nodes before writing. An explicit Neumann
    # node with flux=0 is physically identical to an unlisted node (the
    # solver's natural BC for any node absent from nansfneubc IS
    # zero-flux), but writing it explicitly makes the file noisier than
    # it needs to be, and needlessly inflates NQMAX via update_cathyH
    # below. When outlet_flux=0 and bottom_flux=0 (fully closed run) this
    # empties the list entirely, so we fall back to the same no_flow=True
    # convention update_nansfdirbc/update_sfbc already use above, instead
    # of writing a spurious "explicit list of N nodes, all flux=0" block.
    nonzero_mask = all_flux != 0.0
    all_nodes = all_nodes[nonzero_mask]
    all_flux = all_flux[nonzero_mask]

    if len(all_nodes) == 0:
        simu.update_nansfneubc(time=t_atmbc, no_flow=True)
    else:
        simu.update_nansfneubc(time=t_atmbc, nodes=all_nodes, flux=all_flux)

    outlet_total_area = sum(outlet_areas.values())
    if len(all_nodes) == 0:
        print(f"  BCs: Dirichlet + seepage-face + Neumann all closed "
              f"(no_flow) - outlet_flux and bottom_flux both zero.")
    else:
        print(f"  BCs: Dirichlet + seepage-face closed; Neumann outflow on "
              f"outlet side '{resolved_side}': {outlet_flux:.1e} m/s x exact "
              f"per-node face area (side total {outlet_total_area:.3g} m^2 "
              f"over {len(outlet_node_ids)} nodes). Outlet diagnostic node: "
              f"{outlet_node_id}.")
    if bottom_active:
        print(f"  Free bottom drainage: {bottom_flux:.1e} m/s x nodal "
              f"area {bottom_nodal_area:.3g} m^2, over {n_bottom_nodes} "
              f"bottom nodes.")
    else:
        print(f"  Bottom: closed (no-flow).")
    return resolved_side, outlet_node_id




def plot_forcing(t_atmbc, net_flux, label="ETp forcing") -> None:
    """Bar chart of daily net flux forcing."""
    fig, ax = plt.subplots(figsize=(15, 4))
    ax.bar(range(len(net_flux)), net_flux, width=1.0, color="skyblue")
    ax.set_title(label)
    ax.set_ylabel("m/s")
    ax.set_xlabel("Time step")
    plt.tight_layout()


# ── LAI helpers ───────────────────────────────────────────────────────────────

def lai_to_veg_map(simu, lai_2d: np.ndarray) -> tuple[np.ndarray, dict[int, int]]:
    """
    Convert a 2-D LAI field to an integer vegetation-class map with
    *consecutive* classes starting at 1, as required by pyCATHY.

    Raw semantic classes (calibrated to Agramon observed range 0 – 0.6 m² m⁻²)
    ---------------------------------------------------------------------------
    1 : LAI ≤ 0.10           (bare soil)
    2 : 0.10 < LAI ≤ 0.25   (sparse)
    3 : 0.25 < LAI ≤ 0.40   (moderate)
    4 : 0.40 < LAI ≤ 0.55   (dense)
    5 : LAI > 0.55           (very dense)

    Returns
    -------
    veg_remapped : (M, N) int array — consecutive classes 1 … K
    remap        : {original_semantic_class → remapped_class}
    """
    N = int(simu.hapin["N"])
    M = int(simu.hapin["M"])
    veg = np.ones((M, N), dtype=int)
    for i, thresh in enumerate(LAI_THRESHOLDS, start=2):
        veg[lai_2d > thresh] = i

    # NOTE: remap is intentionally the IDENTITY mapping over the fixed set of
    # semantic classes 1..N_VEG_CLASSES — NOT derived from np.unique(veg).
    # Re-deriving it per month (as before) renumbered whichever classes were
    # spatially present that month to consecutive 1..K, so the same veg-index
    # meant a different semantic class from one month to the next, and the
    # SOIL table size (and MAXVEG) mismatched whenever fewer than 5 classes
    # were present. Keeping it fixed means the soil table always has
    # N_VEG_CLASSES rows and row i always corresponds to semantic class i,
    # whether or not that class actually appears in the DEM this month.
    remap: dict[int, int] = {c: c for c in range(1, N_VEG_CLASSES + 1)}
    veg_remapped = veg  # already using the fixed semantic classes 1..N_VEG_CLASSES
    return veg_remapped, remap


def smooth_veg_map(
    veg: np.ndarray,
    min_patch_pixels: int = 4,
    n_classes: int | None = None,
) -> np.ndarray:
    """
    Remove isolated / salt-and-pepper pixels from a vegetation-class map
    using connected-component minimum-patch removal.
    """
    if n_classes is None:
        n_classes = int(veg.max())

    out = veg.copy()

    for cls in range(1, n_classes + 1):
        binary = out == cls
        labeled, n_features = sp_label(binary)

        for feat_id in range(1, n_features + 1):
            patch = labeled == feat_id
            if patch.sum() >= min_patch_pixels:
                continue

            border = binary_dilation(patch, iterations=1) & ~patch
            if not border.any():
                continue

            neighbour_vals = out[border]
            replacement = int(mode(neighbour_vals, keepdims=True).mode[0])
            out[patch] = replacement

    return out


def update_soil_for_active_zones(
    simu,
    active_zones: set[int],
    remap: dict[int, int],
    base_zroot: float,
    pmin: float,
    zroot_factors: dict[int, float] | None = None,
    permx: float | None = None,
    permy_ratio: float = 1.0,
    permz_ratio: float = 1.0,
    vgncell: float | None = None,
    vgrmccell: float | None = None,
    vgpsatcell: float | None = None,
) -> None:
    """
    ``zroot_factors`` overrides the module-level ZROOT_FACTORS mapping
    (semantic veg class 1-5 → root-depth multiplier) when given —
    normally built from ``--zroot-factors`` — falling back to
    ZROOT_FACTORS itself when None, so existing callers that don't pass
    it keep the original hardcoded behaviour.

    ``permx`` (m/s), when given — normally from ``--permx`` — overrides
    PERMX uniformly across every row of the SPP_map returned by
    ``simu.set_SOIL_defaults(SPP_map_default=True)``; PERMY/PERMZ follow
    via ``permy_ratio``/``permz_ratio`` (anisotropy). Left untouched
    (pyCATHY's own SPP default) when permx is None.

    ``vgncell`` / ``vgrmccell`` / ``vgpsatcell``, when given, override the
    van Genuchten retention-curve columns (VGNCELL = n, VGRMCCELL =
    residual water content θr, VGPSATCELL = 1/alpha) uniformly across
    every row of df_spp, same broadcast pattern as permx above. Each is
    left at pyCATHY's own SPP default when None. Note VGRMCCELL must
    stay below POROS — simu.update_soil() itself warns (doesn't raise)
    if VGRMCCELL >= POROS for any row.
    """
    inv_remap: dict[int, int] = {v: k for k, v in remap.items()}
    if zroot_factors is None:
        zroot_factors = ZROOT_FACTORS

    # Always build the SOIL vegetation table for the FULL fixed set of
    # N_VEG_CLASSES zones — matching MAXVEG — regardless of how many of
    # them are spatially present in the DEM this particular month. Zones
    # not currently active still get a physically-meaningful ZROOT (their
    # fixed zroot_factors value), they just aren't occupying any pixels
    # this month; they must not be dropped from the table, or the row
    # count stops matching MAXVEG and/or rows shift relative to the
    # semantic class each veg-index in the raster refers to.
    n_zones = N_VEG_CLASSES

    df_fp  = simu.set_SOIL_defaults(FP_map_default=True, nveg=n_zones)
    df_spp = simu.set_SOIL_defaults(SPP_map_default=True)

    template = df_fp.iloc[[0]]
    df_fp = pd.concat([template] * n_zones, ignore_index=True)
    # pd.concat(..., ignore_index=True) resets the index to a plain 0-based
    # RangeIndex(0, n_zones). But set_SOIL_defaults()/init_soil_FP_map_df()
    # index this table 1..n_zones ("Veg nb"), and cathy_tools.py's
    # _prepare_SOIL_vegetation_tb() reads it back with
    # FP_map[sfp].loc[iveg] for iveg in range(1, MAXVEG + 1) — so losing
    # the 1-based labels here means the last class (iveg == n_zones) is
    # never found, raising KeyError: n_zones. Re-attach the expected
    # labels; row order is already class 1..n_zones, so this is just
    # restoring the index, not reordering anything.
    df_fp.index = pd.RangeIndex(1, n_zones + 1)
    df_fp.index.name = "Veg nb"
    for key in df_fp.columns:
        if key == "ZROOT":
            df_fp[key] = [
                base_zroot * zroot_factors.get(inv_remap.get(z, z), 1.0)
                for z in range(1, n_zones + 1)
            ]

    permx_label = "pyCATHY default"
    if permx is not None:
        overrides = {
            "PERMX": permx,
            "PERMY": permx * permy_ratio,
            "PERMZ": permx * permz_ratio,
        }
        for col, val in overrides.items():
            if col in df_spp.columns:
                df_spp[col] = val
            else:
                print(f"    WARNING: SPP_map has no '{col}' column — "
                      f"--permx/--permy-ratio/--permz-ratio ignored for it.")
        permx_label = (f"PERMX={permx:.3e}, PERMY={overrides['PERMY']:.3e}, "
                       f"PERMZ={overrides['PERMZ']:.3e} m/s")

    # Retention curve (van Genuchten) overrides — VGNCELL/VGRMCCELL/
    # VGPSATCELL — same broadcast-to-every-row pattern as permx above.
    vg_overrides = {
        "VGNCELL": vgncell,
        "VGRMCCELL": vgrmccell,
        "VGPSATCELL": vgpsatcell,
    }
    vg_label_parts = []
    for col, val in vg_overrides.items():
        if val is None:
            continue
        if col in df_spp.columns:
            df_spp[col] = val
            vg_label_parts.append(f"{col}={val:.4g}")
        else:
            print(f"    WARNING: SPP_map has no '{col}' column — override ignored for it.")
    if (df_spp["VGRMCCELL"] >= df_spp["POROS"]).any():
        print("    WARNING: VGRMCCELL >= POROS for at least one zone — "
              "residual water content exceeds porosity, fix before running.")
    vg_label = ", ".join(vg_label_parts) if vg_label_parts else "pyCATHY default"

    simu.update_soil(FP_map=df_fp, SPP_map=df_spp, PMIN=pmin)
    semantic_active = sorted(inv_remap.get(z, z) for z in active_zones)
    print(f"    → Soil updated | table rows: {n_zones} (all zones, fixed) "
          f"| spatially active this month: {sorted(active_zones)} "
          f"(semantic: {semantic_active}) | Ks: {permx_label} | "
          f"retention curve: {vg_label}")


def iter_monthly_lai(lai_ds: xr.Dataset):
    """Yield the first (datetime, lai_2d) of each calendar month."""
    seen = set()
    for t in lai_ds["time"].values:
        dt  = pd.Timestamp(t).to_pydatetime()
        key = (dt.year, dt.month)
        if key in seen:
            continue
        seen.add(key)
        yield dt, lai_ds["LAI"].sel(time=t).values


def subset_atmbc_for_month(t_atmbc, net_flux: np.ndarray,
                           valid_times: np.ndarray,
                           year: int, month: int):
    """
    Subset the full ATMBC vectors to a single calendar month
    and re-zero the time axis so it starts at 0.

    net_flux may be 1-D (n_steps,) or 2-D (n_steps, n_nodes).
    Returns (None, None) when no data exist for that month.
    """
    mask = np.array([
        pd.Timestamp(t).year == year and pd.Timestamp(t).month == month
        for t in valid_times
    ])
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return None, None
    t_sub  = np.array(t_atmbc)[idx]
    nv_sub = net_flux[idx]        # works for both 1-D and 2-D
    t_sub  = (t_sub - t_sub[0]).tolist()
    return t_sub, nv_sub

def save_forcing_diagnostics(
    valid_times,
    net_flux_1d,
    rain_nodes=None,
    ETp_nodes=None,
    net_flux_2d=None,
    outdir="forcing_plots",
):
    """
    Save forcing diagnostic figures.

    Parameters
    ----------
    valid_times : array-like
        Time axis.
    net_flux_1d : ndarray (nt,)
        Mean net flux [m/s].
    rain_nodes : ndarray (nt, nnodes), optional
        Rainfall forcing [m/s].
    ETp_nodes : ndarray (nt, nnodes), optional
        ETp forcing [m/s].
    net_flux_2d : ndarray (nt, nnodes), optional
        Spatially distributed net flux [m/s].
    outdir : str or Path
        Output directory.
    """

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # Net flux timeseries
    # --------------------------------------------------
    net_mmday = net_flux_1d * 86400 * 1000

    fig, ax = plt.subplots(figsize=(15, 4))

    ax.plot(
        valid_times,
        net_mmday,
        color="k",
        lw=1.5,
        label="Rain - ETp",
    )

    ax.axhline(0, color="0.5", ls="--", lw=0.8)

    ax.fill_between(
        valid_times,
        net_mmday,
        0,
        where=(net_mmday > 0),
        color="royalblue",
        alpha=0.4,
        label="Net input",
    )

    ax.fill_between(
        valid_times,
        net_mmday,
        0,
        where=(net_mmday < 0),
        color="firebrick",
        alpha=0.4,
        label="Net ET loss",
    )

    ax.set_title("Net Atmospheric Flux")
    ax.set_ylabel("mm/day")
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / "net_flux_timeseries.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --------------------------------------------------
    # Components
    # --------------------------------------------------
    if rain_nodes is not None and ETp_nodes is not None:

        mean_rain = rain_nodes.mean(axis=1) * 86400 * 1000
        mean_etp = ETp_nodes.mean(axis=1) * 86400 * 1000

        fig, ax = plt.subplots(figsize=(15, 5))

        ax.plot(valid_times, mean_rain,
                color="blue", label="Rain")

        ax.plot(valid_times, mean_etp,
                color="orange", label="ETp")

        ax.plot(valid_times, net_mmday,
                color="black", lw=2,
                label="Rain - ETp")

        ax.axhline(0, color="0.5", ls="--")

        ax.set_title("Atmospheric Forcing Components")
        ax.set_ylabel("mm/day")
        ax.grid(alpha=0.3)
        ax.legend()

        fig.tight_layout()
        fig.savefig(outdir / "forcing_components.png",
                    dpi=300, bbox_inches="tight")
        plt.close(fig)

    # --------------------------------------------------
    # Spatial variability
    # --------------------------------------------------
    if net_flux_2d is not None:

        p10 = np.percentile(net_flux_2d, 10, axis=1) * 86400 * 1000
        p50 = np.percentile(net_flux_2d, 50, axis=1) * 86400 * 1000
        p90 = np.percentile(net_flux_2d, 90, axis=1) * 86400 * 1000

        fig, ax = plt.subplots(figsize=(15, 4))

        ax.plot(valid_times, p50,
                color="k", label="Median")

        ax.fill_between(
            valid_times,
            p10,
            p90,
            alpha=0.3,
            label="10–90 %",
        )

        ax.set_title("Spatial Distribution of Net Flux")
        ax.set_ylabel("mm/day")
        ax.grid(alpha=0.3)
        ax.legend()

        fig.tight_layout()
        fig.savefig(
            outdir / "net_flux_spatial_variability.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)

        # --------------------------------------------------
        # Annual water balance
        # --------------------------------------------------
        if rain_nodes is not None and ETp_nodes is not None:

            import pandas as pd

            # Mean forcing over surface nodes
            rain_mmday = rain_nodes.mean(axis=1) * 86400 * 1000
            etp_mmday  = ETp_nodes.mean(axis=1) * 86400 * 1000
            net_mmday  = rain_mmday - etp_mmday

            df = pd.DataFrame(
                {
                    "rain": rain_mmday,
                    "etp": etp_mmday,
                    "net": net_mmday,
                },
                index=pd.to_datetime(valid_times),
            )

            # Annual totals [mm/year]
            annual = df.resample("Y").sum()
            annual.index = annual.index.year

            fig, ax = plt.subplots(figsize=(10, 5))

            x = np.arange(len(annual))
            w = 0.28

            ax.bar(
                x - w,
                annual["rain"],
                width=w,
                color="royalblue",
                label="Rain",
            )

            ax.bar(
                x,
                annual["etp"],
                width=w,
                color="darkorange",
                label="ETp",
            )

            ax.bar(
                x + w,
                annual["net"],
                width=w,
                color="forestgreen",
                label="Net",
            )

            ax.axhline(0, color="k", lw=0.8)

            ax.set_xticks(x)
            ax.set_xticklabels(annual.index.astype(str))

            ax.set_ylabel("Annual total [mm]")
            ax.set_title("Annual Water Balance")
            ax.legend()

            fig.tight_layout()

            fig.savefig(
                outdir / "annual_water_balance.png",
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(fig)

            years = np.unique(pd.to_datetime(valid_times).year)

            ncols = 3
            nrows = int(np.ceil(len(years) / ncols))

            fig, axes = plt.subplots(
                nrows,
                ncols,
                figsize=(4*ncols, 3*nrows),
                sharex=True,
                sharey=True,
            )

            axes = np.atleast_1d(axes).ravel()

            for ax, year in zip(axes, years):

                mask = pd.to_datetime(valid_times).year == year

                ax.hist(
                    net_mmday[mask],
                    bins=40,
                    color="steelblue",
                    alpha=0.8,
                )

                ax.axvline(0, color="k", ls="--", lw=1)
                ax.set_title(str(year))

            for ax in axes[len(years):]:
                ax.set_visible(False)

            fig.suptitle("Distribution of Daily Net Flux by Year")
            fig.tight_layout()

            fig.savefig(
                outdir / "net_flux_histogram_by_year.png",
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(fig)

            # ==================================================
            # Monthly climatology
            # ==================================================

            monthly = (
                df.groupby(df.index.month)
                  .sum()
            )

            month_labels = [
                "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
            ]

            fig, ax = plt.subplots(figsize=(12, 5))

            x = np.arange(12)
            w = 0.28

            ax.bar(
                x - w,
                monthly["rain"],
                width=w,
                color="royalblue",
                label="Rain",
            )

            ax.bar(
                x,
                monthly["etp"],
                width=w,
                color="darkorange",
                label="ETp",
            )

            ax.bar(
                x + w,
                monthly["net"],
                width=w,
                color="forestgreen",
                label="Net",
            )

            ax.axhline(0, color="k", lw=0.8)

            ax.set_xticks(x)
            ax.set_xticklabels(month_labels)

            ax.set_ylabel("Accumulated forcing [mm]")
            ax.set_title("Monthly Water Balance (all years combined)")
            ax.legend()

            fig.tight_layout()

            fig.savefig(
                outdir / "monthly_water_balance.png",
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(fig)

            fig, ax = plt.subplots(figsize=(12, 5))

            data = [
                df.loc[df.index.month == m, "net"].values
                for m in range(1, 13)
            ]

            ax.boxplot(
                data,
                labels=month_labels,
                showfliers=False,
            )

            ax.axhline(0, color="k", ls="--")

            ax.set_ylabel("Net flux [mm/day]")
            ax.set_title("Monthly Distribution of Net Flux")

            fig.tight_layout()

            fig.savefig(
                outdir / "monthly_net_flux_boxplot.png",
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(fig)


# ── ETa "net-forcing" artefact patch ────────────────────────────────────────
# Known issue: netValue = rain - ETp is passed to CATHY as a single combined
# flux (HSPATM=1, IETO=1). CATHY has no way to separate "net infiltration"
# from "ET still happening at the same time" — so on any day/node where
# rain >= ETp, the reported ACT. ETRA collapses to 0, even though real ET was
# likely still occurring near the potential rate. This floors ONLY those
# artefact cells (ACT. ETRA ~ 0 AND rain >= ETp) at ETp; genuinely
# soil-moisture-limited zero-ETa days (dry soil, no rain) are left untouched.
def apply_eta_artefact_floor(
    df_et: pd.DataFrame,
    t_month: list,
    ETp_nodes: np.ndarray,
    rain_nodes: np.ndarray,
    idx_month: np.ndarray,
    eta_eps: float = 1e-9,
) -> tuple[pd.DataFrame, int]:
    """
    Add an 'ACT. ETRA_patched' column to df_et, flooring the net-forcing
    artefact at ETp. The original 'ACT. ETRA' column is left untouched.

    Parameters
    ----------
    df_et : DataFrame from simu.read_outputs("ET") for one month, with
        columns 'time_sec', 'SURFACE NODE', 'ACT. ETRA'.
    t_month : time axis (seconds from month start) passed to
        simu.update_atmbc() for this month.
    ETp_nodes, rain_nodes : ndarray (n_time_full, n_surf_nodes), m/s —
        the full-run forcing arrays built in step 9.
    idx_month : ndarray
        Indices into ETp_nodes/rain_nodes' time axis for this month
        (same mask used to build t_month/nv_month via
        subset_atmbc_for_month).
    eta_eps : float
        ACT. ETRA magnitude (m/s) at/below which a value is treated as
        CATHY's "collapsed to 0" artefact rather than a real reading.

    Returns
    -------
    (df_et, n_patched) : patched copy of df_et, and count of cells patched.
    """
    if len(df_et) == 0 or len(t_month) == 0:
        df_et = df_et.copy()
        df_et["ACT. ETRA_patched"] = df_et.get("ACT. ETRA", pd.Series(dtype=float))
        return df_et, 0

    t_month_arr = np.asarray(t_month, dtype=float)
    tsec        = df_et["time_sec"].to_numpy(dtype=float)

    # Nearest-time match (defensive against float round-trip through CATHY I/O)
    i_right  = np.clip(np.searchsorted(t_month_arr, tsec), 0, len(t_month_arr) - 1)
    i_left   = np.clip(i_right - 1, 0, len(t_month_arr) - 1)
    use_left = np.abs(t_month_arr[i_left] - tsec) < np.abs(t_month_arr[i_right] - tsec)
    i_time   = np.where(use_left, i_left, i_right)

    node_idx = df_et["SURFACE NODE"].to_numpy(dtype=int)
    n_surf   = ETp_nodes.shape[1]
    # CATHY node numbering is 1-based (Fortran) — auto-detect and correct
    if node_idx.max() == n_surf:
        node_idx = node_idx - 1
    assert node_idx.min() >= 0 and node_idx.max() < n_surf, (
        f"SURFACE NODE out of bounds for ETp_nodes (n_surf={n_surf}); "
        f"got range [{node_idx.min()}, {node_idx.max()}]"
    )

    rows      = idx_month[i_time]
    etp_vals  = ETp_nodes[rows, node_idx]    # m/s
    rain_vals = rain_nodes[rows, node_idx]   # m/s
    net_vals  = rain_vals - etp_vals

    eta_raw     = df_et["ACT. ETRA"].to_numpy(dtype=float)
    is_artefact = (np.abs(eta_raw) <= eta_eps) & (net_vals >= 0)
    eta_patched = np.where(is_artefact, etp_vals, eta_raw)

    df_et = df_et.copy()
    df_et["ACT. ETRA_patched"] = eta_patched
    return df_et, int(is_artefact.sum())


def add_datetime_col(df: pd.DataFrame, base_dt) -> pd.DataFrame:
    """
    Add a 'datetime' column as base_dt + time_sec (returns df).

    Guards against a pyCATHY output-reading quirk where 'time_sec' is
    sometimes already a Timedelta dtype (instead of a plain float number
    of seconds) — calling pd.to_timedelta(..., unit="s") on data that's
    already Timedelta makes pandas try Timedelta+float arithmetic
    internally and raises "unsupported operand type(s) for +: 'Timedelta'
    and 'float'". When time_sec is already a Timedelta we use it as-is;
    otherwise we coerce to float first and apply unit="s" as before.
    """
    ts = df["time_sec"]
    if pd.api.types.is_timedelta64_dtype(ts):
        df["datetime"] = base_dt + ts
    else:
        df["datetime"] = base_dt + pd.to_timedelta(ts.astype(float), unit="s")
    return df


def check_monthly_sanity(
    dt,
    month_stats: dict,
    psi_month_df,
    sw_month_df,
    df_et,
    df_rec,
    t_month,
    idx_month,
    ETp_nodes,
    rain_nodes,
    sw_tol: float = SW_TOL,
    psi_sanity_min: float = PSI_SANITY_MIN,
    eta_tol: float = ETA_TOL,
    wb_residual_frac: float = WB_RESIDUAL_FRAC,
) -> list[str]:
    """
    Run cheap physical-bounds and (approximate) water-balance sanity
    checks on one month's outputs. Fills diagnostic fields into
    `month_stats` in place and returns a list of human-readable warning
    strings (empty if nothing looked off). Never raises — a check that
    can't be computed for this month (missing data, unexpected columns)
    is silently skipped rather than aborting the run.

    Water-balance note
    -------------------
    The residual is rain_in − ETa_out − recharge_out only: it does NOT
    account for the change in subsurface storage (ΔS), because that
    would require per-node porosity and control-volume information not
    readily available here. So this is not a strict mass-balance
    closure — treat a large residual as "worth checking the psi/sw
    trend for this month", not as proof of a bug on its own. It also
    assumes the 'recharge' output column(s) are in the same units
    (m/s) as ETa/rain; if your CATHY build reports recharge in
    different units, `recharge_out_mm` below will be wrong and should
    be ignored.
    """
    warnings: list[str] = []

    # ── Physical bounds: SW should stay within [0, 1] ──
    if sw_month_df is not None and len(sw_month_df):
        sw_vals = sw_month_df["sw"].to_numpy(dtype=float)
        n_bad = int(((sw_vals < -sw_tol) | (sw_vals > 1 + sw_tol)).sum())
        month_stats["sw_out_of_bounds_n"] = n_bad
        if n_bad:
            warnings.append(
                f"{dt:%Y-%m}: {n_bad} SW value(s) outside [0,1] "
                f"(range [{sw_vals.min():.3f}, {sw_vals.max():.3f}])."
            )

    # ── Physical bounds: PSI shouldn't be sitting at extreme drying values ──
    if psi_month_df is not None and len(psi_month_df):
        psi_vals = psi_month_df["psi"].to_numpy(dtype=float)
        n_extreme = int((psi_vals < psi_sanity_min).sum())
        month_stats["psi_extreme_n"] = n_extreme
        if n_extreme:
            warnings.append(
                f"{dt:%Y-%m}: {n_extreme} PSI value(s) below sanity threshold "
                f"{psi_sanity_min:.0f} m (min {psi_vals.min():.2f} m)."
            )

    # ── Physical bounds: ETa should not exceed ETp ──
    if (df_et is not None and len(df_et) and t_month is not None and len(t_month)
            and ETp_nodes is not None and idx_month is not None and len(idx_month)):
        try:
            t_month_arr = np.asarray(t_month, dtype=float)
            tsec = df_et["time_sec"].to_numpy(dtype=float)
            i_right  = np.clip(np.searchsorted(t_month_arr, tsec), 0, len(t_month_arr) - 1)
            i_left   = np.clip(i_right - 1, 0, len(t_month_arr) - 1)
            use_left = np.abs(t_month_arr[i_left] - tsec) < np.abs(t_month_arr[i_right] - tsec)
            i_time   = np.where(use_left, i_left, i_right)

            node_idx = df_et["SURFACE NODE"].to_numpy(dtype=int)
            n_surf   = ETp_nodes.shape[1]
            if node_idx.max() == n_surf:
                node_idx = node_idx - 1

            rows     = idx_month[i_time]
            etp_vals = ETp_nodes[rows, node_idx]
            eta_col  = "ACT. ETRA_patched" if "ACT. ETRA_patched" in df_et.columns else "ACT. ETRA"
            eta_vals = df_et[eta_col].to_numpy(dtype=float)

            n_exceed = int((eta_vals > etp_vals + eta_tol).sum())
            month_stats["eta_exceeds_etp_n"] = n_exceed
            if n_exceed:
                max_excess = float((eta_vals - etp_vals).max())
                warnings.append(
                    f"{dt:%Y-%m}: {n_exceed} ETa value(s) exceed ETp "
                    f"(max excess {max_excess:.3e} m/s)."
                )
        except Exception:
            pass  # sanity check itself must never break the run

    # ── Water balance: rain_in vs ETa_out vs recharge_out (flux-only) ──
    rain_in_mm = eta_out_mm = recharge_out_mm = None

    if rain_nodes is not None and idx_month is not None and len(idx_month):
        rain_in_mm = float((rain_nodes[idx_month].mean(axis=1) * 86400 * 1000).sum())
        month_stats["rain_in_mm"] = rain_in_mm

    if df_et is not None and len(df_et):
        eta_col = "ACT. ETRA_patched" if "ACT. ETRA_patched" in df_et.columns else "ACT. ETRA"
        if eta_col in df_et.columns:
            eta_by_time = df_et.groupby("time_sec")[eta_col].mean()
            eta_out_mm = float((eta_by_time * 86400 * 1000).sum())
            month_stats["eta_out_mm"] = eta_out_mm

    if df_rec is not None and len(df_rec):
        rec_value_cols = [c for c in df_rec.select_dtypes(include=[np.number]).columns
                           if c not in ("X", "Y", "time_sec")]
        if rec_value_cols:
            rec_by_time = df_rec.groupby("time_sec")[rec_value_cols].mean().mean(axis=1)
            recharge_out_mm = float((rec_by_time * 86400 * 1000).sum())
            month_stats["recharge_out_mm"] = recharge_out_mm

    if rain_in_mm is not None and eta_out_mm is not None and recharge_out_mm is not None:
        residual_mm = rain_in_mm - eta_out_mm - recharge_out_mm
        month_stats["water_balance_residual_mm"] = residual_mm
        if rain_in_mm > 1e-6 and abs(residual_mm) > wb_residual_frac * rain_in_mm:
            warnings.append(
                f"{dt:%Y-%m}: water balance residual {residual_mm:.1f} mm "
                f"(rain {rain_in_mm:.1f}, ETa {eta_out_mm:.1f}, "
                f"recharge {recharge_out_mm:.1f}) — large relative to rain "
                f"input. Excludes storage change (ΔS); also check the "
                f"psi/sw trend for this month before assuming a bug."
            )

    return warnings


# ── Spin-up (pre-2016 equilibration) ────────────────────────────────────────

def pick_default_spinup_year(valid_times: np.ndarray) -> int:
    """
    Pick a full calendar year of forcing to use as the repeating spin-up
    cycle. Skips the first calendar year present in the forcing record if
    it is incomplete (satellite ETp/rain start 2016-10-01, so 2016 only
    has Oct–Dec and cannot be cycled as a full year).
    """
    times      = pd.to_datetime(valid_times)
    first_year = int(times.year.min())
    months_in_first_year = times[times.year == first_year].month
    if len(np.unique(months_in_first_year)) < 12:
        return first_year + 1
    return first_year


def run_spinup(
    simu,
    args: argparse.Namespace,
    t_atmbc: list,
    net_flux: np.ndarray,
    valid_times: np.ndarray,
    ds_LAI_raster: xr.Dataset,
    first_iteration: bool,
    outlet_node_id: int = None,
    md_log: "MarkdownLog | None" = None,
) -> tuple[bool, bool]:
    """
    Equilibrate the pressure-head field before the real simulation period
    ("spin-up antes de 2016"): since satellite ETp/rain forcing only
    starts 2016-10-01, there is no real forcing to run before 2016, so
    the spin-up instead cycles one full calendar year of forcing
    (default: the first complete year available) ``args.spinup_cycles``
    times, warm-starting psi between runs, to bring the model to a
    quasi-equilibrium state before the real simulation loop begins.

    Vegetation is held static during spin-up (the vegetation state of the
    first month of the real record) — only ATMBC forcing and psi evolve.
    Spin-up runs are not written to the output NetCDFs; they only leave
    behind the final psi field on disk, which the main loop's existing
    "warm start from simu.read_outputs('psi')" logic then picks up.

    Parameters
    ----------
    first_iteration : bool
        Whether the CATHY preprocessor has not yet been run. If True,
        the first spin-up run triggers it, and the returned value is
        False so the main loop does not re-run it.
    outlet_node_id : int, optional
        Node id to track psi at across spin-up cycles (see
        ``get_outlet_node``). If given, a diagnostic CSV
        ("spinup_outlet_psi_diagnostic.csv") and plot
        ("spinup_outlet_psi_diagnostic.png") are written at the end of
        spin-up, one line per cycle per month. A monotonically rising
        trend toward 0 across cycles (rather than flattening out) means
        spin-up is drifting toward saturation because of closed
        boundaries, not equilibrating — the diagnostic this issue asked
        for.

    Returns
    -------
    (first_iteration, spinup_done) : tuple[bool, bool]
        Updated ``first_iteration`` flag and whether at least one
        spin-up run completed (used by the main loop to decide whether
        to warm-start ICs from disk even before its own first run).
    """
    spinup_year = args.spinup_year or pick_default_spinup_year(valid_times)
    print(f"  Spin-up cycle year : {spinup_year}  |  cycles : {args.spinup_cycles}")
    if md_log is not None:
        md_log.section("Spin-up (warm-up)")
        md_log.meta({
            "Spin-up year": spinup_year,
            "Cycles": args.spinup_cycles,
            "Outlet diagnostic node": outlet_node_id if outlet_node_id is not None else "n/a",
        })

    # Static vegetation for spin-up: first available month of the real record.
    try:
        dt0, lai0 = next(iter(
            iter_monthly_lai(ds_LAI_raster.isel(time=slice(1, None)))
        ))
    except StopIteration:
        print("  WARNING: no LAI data available — skipping spin-up.")
        if md_log is not None:
            md_log.warn("No LAI data available — spin-up skipped entirely.")
        return first_iteration, False

    veg_map_remapped, remap = lai_to_veg_map(simu, lai0)
    veg_map = smooth_veg_map(veg_map_remapped)
    _, veg_valid = simu._check_outside_DEM(veg_map)
    simu.update_veg_map(veg_map, maxveg=N_VEG_CLASSES)
    active_zones = set(int(z) for z in veg_valid)
    update_soil_for_active_zones(
        simu, active_zones, remap=remap,
        base_zroot=args.zroot, pmin=args.pmin,
        zroot_factors=zroot_factors_from_args(args),
        permx=args.permx,
        permy_ratio=args.permy_ratio,
        permz_ratio=args.permz_ratio,
    )
    if md_log is not None:
        md_log.bullet(f"Static spin-up vegetation from {dt0:%Y-%m}: "
                       f"zones {sorted(active_zones)}")
        md_log.table_header(["Cycle", "Month", "LAI mean", "LAI min–max",
                              "ψ outlet (m)", "Status"])

    spinup_ran = False
    mesh_refreshed = False     # see note at the first update_ic below
    outlet_psi_by_cycle = []   # diagnostic: psi at outlet node, per cycle/month
    for cycle in range(1, args.spinup_cycles + 1):
        for month in range(1, 13):
            t_month, nv_month = subset_atmbc_for_month(
                t_atmbc, net_flux, valid_times, spinup_year, month
            )
            if t_month is None:
                continue

            simu.update_atmbc(HSPATM=1, IETO=1, time=t_month, netValue=nv_month)

            df_atmbc_month = simu.read_inputs("atmbc")
            simu.update_parm(TIMPRTi=list(df_atmbc_month.time), IPRT=4, VTKF=2)

            #if first_iteration:
                #simu.run_preprocessor(verbose=True)
                #first_iteration = False

            if spinup_ran:
                # FIX: same stale-mesh issue the main loop already handles
                # (see "this preprocessor call can finalize a DIFFERENT
                # (larger) mesh" note there). update_ic -> map_prop2mesh
                # validates against the mesh cached by create_mesh_vtk()
                # in setup_dem, which predates the processor finalizing the
                # 3D grid. psi has one value per real 3D node, so without
                # this refresh we get "'ic' has length N but the mesh has M
                # points". Rebuild once, after the first spin-up run.
                if not mesh_refreshed:
                    try:
                        simu.create_mesh_vtk(verbose=False)
                        mesh_refreshed = True
                    except Exception as exc:
                        print(f"  WARNING: could not refresh mesh before "
                              f"spin-up IC warm start ({exc}).")
                        if md_log is not None:
                            md_log.warn(f"Mesh refresh failed in spin-up: {exc}")

                df_psi   = simu.read_outputs("psi").copy()
                psi_ini  = df_psi.iloc[-1].values
                simu.update_ic(INDP=1, pressure_head_ini=psi_ini)

            print(f"  [spin-up {cycle}/{args.spinup_cycles}] "
                  f"{spinup_year}-{month:02d} …")
            # MAXVEG already fixed to N_VEG_CLASSES before the loop; no need
            # to re-set it here every cycle (it was previously re-set to a
            # mismatched value of 10 on every iteration).

            try:
                simu.run_processor(
                    IPRT1=2,
                    DTMIN=args.dtmin,
                    DTMAX=args.dtmax,
                    DELTAT=args.deltat,
                    TRAFLAG=0,
                    verbose=False,
                )
            except Exception as exc:
                print(f"  WARNING: spin-up run failed for "
                      f"{spinup_year}-{month:02d} ({exc}) — stopping spin-up early.")
                if md_log is not None:
                    md_log.table_row([cycle, f"{spinup_year}-{month:02d}",
                                       f"{np.nanmean(lai0):.3f}",
                                       f"{np.nanmin(lai0):.3f}–{np.nanmax(lai0):.3f}",
                                       "n/a", f"✗ FAILED: {exc}"])
                break
            spinup_ran = True

            # Diagnostic: track psi at the outlet node this cycle/month.
            psi_outlet_str = "n/a"
            if outlet_node_id is not None:
                try:
                    df_psi_diag = simu.read_outputs("psi")
                    psi_outlet = float(df_psi_diag.iloc[-1][outlet_node_id])
                    outlet_psi_by_cycle.append(
                        {"cycle": cycle, "month": month, "psi_outlet": psi_outlet}
                    )
                    psi_outlet_str = f"{psi_outlet:.4f}"
                except Exception as exc:
                    print(f"  WARNING: could not read outlet psi for "
                          f"{spinup_year}-{month:02d} ({exc})")
                    psi_outlet_str = "read failed"

            if md_log is not None:
                md_log.table_row([cycle, f"{spinup_year}-{month:02d}",
                                   f"{np.nanmean(lai0):.3f}",
                                   f"{np.nanmin(lai0):.3f}–{np.nanmax(lai0):.3f}",
                                   psi_outlet_str, "✓ OK"])

    if outlet_psi_by_cycle:
        df_diag = pd.DataFrame(outlet_psi_by_cycle)
        df_diag.to_csv("spinup_outlet_psi_diagnostic.csv", index=False)

        fig, ax = plt.subplots(figsize=(8, 4))
        for c, grp in df_diag.groupby("cycle"):
            ax.plot(grp["month"], grp["psi_outlet"], marker="o", label=f"cycle {c}")
        ax.axhline(0, color="k", ls="--", lw=1)
        ax.set_xlabel("Month (within spin-up year)")
        ax.set_ylabel("\u03c8 outlet (m)")
        ax.legend()
        ax.set_title(f"Spin-up convergence check — outlet node {outlet_node_id}")
        fig.savefig("spinup_outlet_psi_diagnostic.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # Cheap monotonic-drift check: compare each cycle's December value
        # (or last available month) to the previous cycle's.
        last_by_cycle = df_diag.sort_values("month").groupby("cycle")["psi_outlet"].last()
        if len(last_by_cycle) > 1:
            diffs = last_by_cycle.diff().dropna()
            if (diffs > 0).all():
                verdict = ("end-of-cycle outlet psi rises monotonically toward 0 "
                            "across every spin-up cycle "
                            f"({list(last_by_cycle.round(4))}) — this looks like "
                            "accumulation from a closed boundary, not equilibration. "
                            "Re-check the outlet BC.")
                print(f"  DIAGNOSTIC WARNING: {verdict}")
                if md_log is not None:
                    md_log.warn(verdict)
            else:
                verdict = ("end-of-cycle outlet psi across cycles: "
                            f"{list(last_by_cycle.round(4))} — not monotonically "
                            "rising, consistent with equilibration rather than "
                            "closed-boundary drift.")
                print(f"  DIAGNOSTIC: {verdict}")
                if md_log is not None:
                    md_log.bullet(f"✓ {verdict}")
        print("  Spin-up outlet-psi diagnostic saved → "
              "spinup_outlet_psi_diagnostic.{csv,png}")
        if md_log is not None:
            md_log.bullet("Diagnostic detail saved → "
                           "`spinup_outlet_psi_diagnostic.csv` / `.png`")

    if not spinup_ran:
        print(f"  WARNING: no forcing found for spin-up year {spinup_year} "
              f"— spin-up skipped.")
        if md_log is not None:
            md_log.warn(f"No forcing found for spin-up year {spinup_year} "
                         f"— spin-up skipped.")
        return first_iteration, False

    print(f"  Spin-up complete ({args.spinup_cycles} cycle(s) of {spinup_year}).")
    if md_log is not None:
        md_log.bullet(f"Spin-up complete: {args.spinup_cycles} cycle(s) of "
                       f"{spinup_year}.")
    return first_iteration, True


# ── Main pipeline ─────────────────────────────────────────────────────────────

def run_pipeline(args: argparse.Namespace) -> None:
    """Execute the full pyCATHY + LAI + ETp simulation pipeline."""

    show_plots = not args.no_plots

    # ── 1. Resolve DEM source ─────────────────────────────────
    dem_tif_path = Path(args.dem_tif)

    print("\n" + "=" * 60)
    print("  LT Subsurface Hydrology with LAI + ETp — pyCATHY")
    print(f"  Forcing source : {'satellite ETp + rain' if USE_ETP_FORCING else 'ERA5'}")
    print(f"  Period         : {args.start_year} – {args.end_year - 1}")
    print(f"  Target CRS     : {TARGET_CRS}")
    print(f"  ZROOT          : {args.zroot} m  |  PMIN : {args.pmin} m")
    print(f"  ZROOT factors  : {args.zroot_factors} "
          f"(bare/sparse/moderate/dense/very_dense)")
    print(f"  PERMX (Ks)     : "
          f"{'pyCATHY default' if args.permx is None else f'{args.permx:.3e} m/s'}"
          f"{'' if args.permx is None else f'  |  PERMY ratio {args.permy_ratio}'
                                            f'  |  PERMZ ratio {args.permz_ratio}'}")
    print(f"  DEM source     : {dem_tif_path}")
    print(f"  Log file       : {args.log_file}")
    print(f"  Spin-up        : "
          f"{'disabled' if args.no_spinup else f'{args.spinup_cycles} cycle(s)'}")
    print("=" * 60)

    # ── 2. Shapefiles ─────────────────────────────────────────
    print("\n[1/9] Loading shapefiles …")
    gdf_Agramon = load_shapefiles()

    # ── 3. Forcing datasets + LAI ─────────────────────────────
    print("\n[2/9] Loading forcing datasets and LAI …")
    ds_LAI = load_monthly_lai(LAI_PATH)

    if USE_ETP_FORCING:
        ds_etp  = load_etp()
        ds_rain = load_rain()
        ds_era5 = None
    else:
        ds_era5 = load_era5(args.start_year, args.end_year)
        ds_etp  = None
        ds_rain = None

    # ── 4. Log simulation ─────────────────────────────────────
    print("\n[3/9] Logging simulation parameters …")
    log_path  = Path(".").resolve() / args.log_file
    sim_index = log_simulation(
        log_path,
        start_year=args.start_year,
        end_year=args.end_year,
        zroot=args.zroot,
        pmin=args.pmin,
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
    )

    # Single source of truth for this run's folder name (project dir,
    # outputs dir, markdown log, and — via withLAI_results.py — the
    # figures dir all reuse it), so the same run always lands in the
    # same, comparably-named place. See build_scenario_dirname() above.
    scenario_dirname = build_scenario_dirname(
        sim_index, args.dem_plot, args.outlet_flux, args.bottom_flux,
        args.no_spinup, args.spinup_cycles, args.pmin, args.ic_mode,
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
    if USE_ETP_FORCING:
        md_log.fixes_applied([
            "ETa net-forcing artefact floor: `apply_eta_artefact_floor()` adds "
            "an `ACT. ETRA_patched` column, flooring cells where CATHY reports "
            "`ACT. ETRA ≈ 0` because net forcing (rain − ETp) was ≥ 0 that "
            "day/node. Original `ACT. ETRA` is left unmodified.",
        ])
    md_log.section("Run configuration")
    md_log.meta({
        "Forcing source": "satellite ETp + rain" if USE_ETP_FORCING else "ERA5",
        "Period": f"{args.start_year} – {args.end_year - 1}",
        "Target CRS": TARGET_CRS,
        "ZROOT / PMIN": f"{args.zroot} m / {args.pmin} m",
        "ZROOT factors (bare/sparse/moderate/dense/very_dense)":
            ",".join(f"{v:g}" for v in args.zroot_factors),
        "PERMX / Ks override": (
            "pyCATHY default" if args.permx is None
            else f"PERMX={args.permx:.3e} m/s, PERMY ratio={args.permy_ratio}, "
                 f"PERMZ ratio={args.permz_ratio}"
        ),
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

    # ── 5. DEM ────────────────────────────────────────────────
    print("\n[4/9] Loading and masking DEM …")
    raster_DEM, raster_DEM_masked, xllcorner, yllcorner, res_x, res_y = load_dem_from_tif(
        dem_tif_path,
        gdf_clip=gdf_Agramon,
        fid=args.dem_plot,
        resample_resolution=2
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

    if USE_ETP_FORCING:
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
            "ETp source file": str(etp_path),
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
    ds_LAI_raster = ds_LAI.interp(
        x=xr.DataArray(raster_DEM.x.values, dims="x"),
        y=xr.DataArray(raster_DEM.y.values, dims="y"),
        method="linear",
    )
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
    if USE_ETP_FORCING and ETp_nodes is not None:
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
                    "the satellite ET0 raster onto the mesh (node index "
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

    def save_combined_outputs() -> None:
        """
        (Re)write the single running output files in out_dir from whatever
        has been accumulated in the *_records lists / veg_map_history so
        far. Safe to call after every month: each call overwrites the
        previous version of these files with the latest cumulative data,
        so — unlike the per-month files in monthly_dir — there's always
        exactly one psi_output.nc / sw_output.nc / et_output.nc /
        recharge_output.nc / veg_map_history.nc reflecting all months
        completed up to now, without waiting for the loop to finish.
        """
        if veg_map_history:
            veg_history_path = out_dir / "veg_map_history.pkl"
            with open(veg_history_path, "wb") as fh:
                pickle.dump(veg_map_history, fh)

            sorted_dates = sorted(veg_map_history)
            veg_arr = np.stack([veg_map_history[d] for d in sorted_dates], axis=0)
            veg_xr = xr.DataArray(
                veg_arr,
                dims=["time", "row", "col"],
                coords={"time": sorted_dates, "row": np.arange(M), "col": np.arange(N)},
                name="veg_map",
                attrs={
                    "description": "Vegetation zone map (1=bare, 2=sparse, 3=moderate, 4=dense, 5=very dense)",
                    "crs": TARGET_CRS,
                },
            )
            veg_xr.to_netcdf(out_dir / "veg_map_history.nc")

        if et_records:
            ET_xr_all = xr.concat(et_records, dim="datetime")
            ET_xr_all.time.attrs.pop("dtype", None)
            ET_xr_all.to_netcdf(out_dir / "et_output.nc",
                                 encoding={"time": {"dtype": "float64"}})

        if psi_records:
            df_psi_all = pd.concat(psi_records, ignore_index=True)
            df_psi_all = (
                df_psi_all
                .groupby(["datetime", "node"], as_index=False)
                .last()
            )
            df_psi_all.set_index(["datetime", "node"]).to_xarray().to_netcdf(
                out_dir / "psi_output.nc"
            )

        if sw_records:
            df_sw_all = pd.concat(sw_records, ignore_index=True)
            df_sw_all = (
                df_sw_all
                .groupby(["datetime", "node"], as_index=False)
                .last()
            )
            df_sw_all.set_index(["datetime", "node"]).to_xarray().to_netcdf(
                out_dir / "sw_output.nc"
            )

        if all_recharge_recs:
            xr.concat(all_recharge_recs, dim="datetime").to_netcdf(
                out_dir / "recharge_output.nc")

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
                                     zroot_factors=zroot_factors_from_args(args),
                                     permx=args.permx,
                                     permy_ratio=args.permy_ratio,
                                     permz_ratio=args.permz_ratio,
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
        if USE_ETP_FORCING:
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
        ax2.set_title(f"Net flux ({'ETp' if USE_ETP_FORCING else 'ERA5'}) - {dt:%Y-%m}")
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

            if USE_ETP_FORCING and n_surf_real != n_surf:
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
            if USE_ETP_FORCING:
                df_et, n_eta_patched = apply_eta_artefact_floor(
                    df_et, t_month, ETp_nodes, rain_nodes, idx_month,
                )
                if n_eta_patched:
                    print(f"  ETa artefact floor: {n_eta_patched} cells patched for {dt:%Y-%m}")
                md_log.bullet(f"{dt:%Y-%m}: ETa artefact floor patched "
                               f"{n_eta_patched} cell(s) (net rain ≥ ETp).")
                md_log.add_metric("ETa artefact cells patched (total)", n_eta_patched)

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

        # 12.11b  Sanity checks (physical bounds + approximate water balance)
        try:
            sanity_warnings = check_monthly_sanity(
                dt, month_stats,
                psi_month_df=psi_month_df,
                sw_month_df=sw_month_df,
                df_et=df_et,
                df_rec=df_rec,
                t_month=t_month,
                idx_month=idx_month,
                ETp_nodes=ETp_nodes if USE_ETP_FORCING else None,
                rain_nodes=rain_nodes if USE_ETP_FORCING else None,
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
            save_combined_outputs()
        except Exception as exc:
            print(f"  Warning: could not refresh combined outputs after {dt:%Y-%m} ({exc})")

    print("\n── Loop complete (partial runs saved if any months were skipped) ────────")

    # ── 13/14. Final refresh of the combined outputs (veg history, PSI,
    # SW, ET, recharge) — same helper called after every month in the
    # loop above, so this call just guarantees the very last month is
    # reflected too.
    save_combined_outputs()
    print(f"veg_map_history saved → {out_dir / 'veg_map_history.pkl'}")
    if et_records:
        print(f"ET saved → {out_dir / 'et_output.nc'}")
    if psi_records:
        print(f"PSI saved → {out_dir / 'psi_output.nc'}")
    if sw_records:
        print(f"SW saved → {out_dir / 'sw_output.nc'}")
    if all_recharge_recs:
        print(f"Recharge saved → {out_dir / 'recharge_output.nc'}")

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
