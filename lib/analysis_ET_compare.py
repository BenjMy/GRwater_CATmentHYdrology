"""
analyse_ET_compare.py
=====================
GRwater project — Agramon study site

Compares two ETa sources side-by-side:
  A) Earth-Observation / Energy-Balance (EO-EB)
       — monthly NetCDF from the TSEB/biophysical pipeline
         (same file used by analyse_ET.py: 20161001_20241231_ET.nc)
  B) Water-Balance model (WB)
       — et_output.nc produced by LT_SShydro_run_withLAI.py
         (variable 'ACT. ETRA', units m/s, dims X/Y/datetime)

Figures produced
----------------
Saved under ``<fig_path>/<scenario_dirname>/``, where scenario_dirname is
the resolved (possibly tagged) scenario_<sim_index>[_tag] folder name for
--scenario (e.g. scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c) —
same naming convention as withLAI_results.py's figures/ tree, so each
scenario's ET-compare plots land in their own folder instead of
overwriting one another in a flat directory.

  compare_ET_monthly_mean_by_year_<label>.png
      Spatial-mean monthly ETa, one line per year — EO-EB vs WB on the same axes.

  compare_ET_annual_total_<label>.png
      Annual domain-mean ETa bar chart — EO-EB (blue) vs WB (orange).

  compare_ET_monthly_total_<label>.png
      Monthly domain-mean ETa bar chart — EO-EB vs WB grouped bars.

  compare_ET_scatter_<label>.png
      Scatter plot of coincident monthly spatial-mean values:
      EO-EB (x) vs WB (y), coloured by year.

  compare_ET_spatial_monthly_<year>_<label>.png  [one per source]
      Monthly spatial maps (mosaic) for a selected year — EO-EB and WB
      rendered with a shared colour scale for direct comparison.

  compare_ET_monthly_mean_by_year_per_zone_<label>.png
      Same as the domain-mean monthly-by-year figure, but faceted into
      LAI zones (default: Bare / Moderate veg. / Dense veg., tercile-split
      on the PRE-FIRE mean LAI so zone boundaries reflect pre-disturbance
      structure). One panel per zone, EO-EB (solid) vs WB (dashed).

  compare_ET_annual_bias_per_zone_<label>.png
      Annual mean bias (EO-EB minus WB), one line per LAI zone, fire date
      marked. Isolates whether disagreement concentrates in high-LAI
      zones pre-fire and collapses toward the bare-soil-zone behaviour
      post-fire (mechanism check, not just domain-mean skill).

<<<<<<< HEAD
  compare_ETa_vs_csv_<label>.png   [only for forcing=etp_ts runs]
      Simulated WB ETa vs the Agramon treatments CSV (agramon_treatments_ET.csv)
      for the run's feature_id: real ET ("ET"), ETp dry bound ("ET_0_dry") and
      ETp wet bound ("ET_0_wet"). Panels: daily series, monthly means (with
      EO-EB), daily scatter vs real ET. A merged CSV with the same label is
      written next to it and RMSE / bias / r are printed. The forcing column
      (dry|wet) and feature_id are parsed from the scenario folder name
      (..._forcing-etpts-<dry|wet>-fid<K>); override with --csv-fid, disable
      with --no-csv, point to another file with --et-csv.

  Wet / dry matching rule (all EO-EB-vs-WB and CSV figures)
      A run forced with the WET (DRY) ETp bound is compared only with the
      WET (DRY) reference: Fig H draws just the matching ETp bound
      (ET_0_wet / ET_0_dry) plus the real-ET column for that condition
      (ET_wet / ET_dry if the CSV has them, else the single "ET" column);
      Figs A-G use the EO-EB variable for that condition (ET_wet / ET_dry)
      when the EO file provides one, else its first variable. Condition is
      parsed from the folder name (..._forcing-etpts-<wet|dry>-fid<K>).

=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
Usage
-----
    python analyse_ET_compare.py                          # default plot 0, scenario 0
    python analyse_ET_compare.py --plots 0 2              # by index
    python analyse_ET_compare.py --plots Plot4Mulching    # by name
    python analyse_ET_compare.py --plots all              # all plots
    python analyse_ET_compare.py --scenario 3             # WB scenario index
    python analyse_ET_compare.py --year 2021              # spatial-map year
    python analyse_ET_compare.py --n-zones 3               # LAI-zone count
    python analyse_ET_compare.py --no-zones                # skip zone figures
<<<<<<< HEAD
    python analyse_ET_compare.py --scenario 7 --path2prj ../SSHydro_withETp_withLAI \
                                  --log-file simulation_log_withETp_withLAI.csv
        # point --scenario's resolution at a specific project tree/log —
        # same values plot_scenarii_withETp_withLAI.sh passes to
        # withLAI_results.py for that experiment. Without these, --scenario
        # is resolved against the WB_PRJ_ROOT/LOG_FILE module defaults.
=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import geopandas as gpd
from shapely.geometry import box
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401
import xarray as xr

# =============================================================================
# Paths & constants
# =============================================================================

EO_ET_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/GRwater/data/satellite/Agramon"
    "/20161001_20241231_ET.nc"
)

# Root of CATHY project — simulation artefacts live under <WB_PRJ_ROOT>/sim<N>/outputs/
# WB_PRJ_ROOT = Path(
#     "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
#     "/Tech4agro_org/GRwater_CATmentHYdrology"
# )
WB_PRJ_ROOT = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Training_Supervision/"
    "Supervision/Xela_Carracedo_Practicas_2026_data/SSHydro_withETp_withLAI"
)

# Shared CSV log produced by LT_SShydro_run_withLAI.py — same file used by
# withLAI_results.py to resolve a --scenario index into the actual (possibly
# tagged) scenario_<N>[_tag] output folder. Resolved under
# <WB_PRJ_ROOT>/../Codigos_HydroModel/<LOG_FILE> via _resolve_log().
LOG_FILE = "simulation_log_withETp_withLAI.csv"


MODULE_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
    "/Tech4agro_org/GRwater_geophy"
).resolve()
if str(MODULE_PATH) not in sys.path:
    sys.path.append(str(MODULE_PATH))

import Agramon_utils as AgUtils

TARGET_CRS  = "EPSG:25830"
NATIVE_CRS  = "EPSG:32630"

#WB_ET_VAR   = "ACT. ETRA"
WB_ET_VAR   = "ACT. ETRA_patched"
WB_ET_SCALE = 1e3 * 86400        # m/s  → mm/day

# Agramon ET CSV (per-treatment, long format) — same folder as the EO NetCDF.
# Columns used: ET (real ET), ET_0_dry, ET_0_wet (mm/day), feature_id, t.
ET_CSV_PATH = EO_ET_PATH.parent / "agramon_treatments_ET.csv"
COLOR_CSV_ET  = "black"
COLOR_CSV_DRY = "#e08214"
COLOR_CSV_WET = "#2c7fb8"

CMAP_EO = "Blues"
CMAP_WB = "Oranges"
COLOR_EO = "#2166ac"
COLOR_WB = "#d6604d"

fig_path = Path().cwd() / "../Analysis/ET_compare"

# -----------------------------------------------------------------------------
# LAI-zone comparison (per-zone ETa, not just domain-mean)
# -----------------------------------------------------------------------------
# Same monthly LAI product used by withLAI_results.py / Agramon_withLAI_
# withETp.py (LAI_PATH / "lai_monthly.nc") — NOT alongside EO_ET_PATH.
# (Previously guessed as EO_ET_PATH.parent / "20161001_20241231_LAI.nc",
# which doesn't exist, silently disabling Fig F/G every run.)
LAI_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Training_Supervision"
    #"/home/ben/Nextcloud/BenCSIC/Training_Supervision"
    "/Supervision/Xela_Carracedo_Practicas_2026_data/301a-biophysical/agramon/input"
    "/lai_monthly.nc"
)

# Zones are tercile-split (by default) on the PRE-FIRE time-mean LAI, so
# zone boundaries are fixed at the pre-disturbance vegetation structure —
# any post-fire convergence between EO-EB and WB inside a formerly
# high-LAI zone then shows up as a change in that zone's bias over time,
# rather than as a re-shuffling of which pixels belong to which zone.
N_LAI_ZONES  = 3
ZONE_LABELS  = ("Bare", "Moderate veg.", "Dense veg.")
CMAP_ZONES   = mcolors.ListedColormap(["#d4b483", "#78c679", "#006837"])

# Duplicated locally (rather than imported from results_plotter, which
# pulls in pyCATHY) so this script keeps its own lighter dependency set.
FIRE_DATE = pd.Timestamp("2020-07-01")


def add_fire_vline(
    ax,
    d0: pd.Timestamp | None = None,
    d1: pd.Timestamp | None = None,
    *,
    label: str = "Fire (Jul 2020)",
    color: str = "#B22222",
) -> None:
    """Draw a dashed vertical line marking the Agramon fire date, unless
    the supplied (d0, d1) window falls entirely outside FIRE_DATE."""
    if d0 is not None and FIRE_DATE < pd.Timestamp(d0):
        return
    if d1 is not None and FIRE_DATE > pd.Timestamp(d1):
        return
    ax.axvline(FIRE_DATE, color=color, linestyle="--", linewidth=1.3,
               alpha=0.85, zorder=4, label=label)




# =============================================================================
# Log helpers
# =============================================================================

def _resolve_log(path2prj: Path, log_file: str) -> Path:
    lp = Path(log_file)
    if not lp.is_absolute():
        # .resolve() elimina los '..' y convierte la ruta en absoluta
        lp = (path2prj / 'lib' / lp).resolve()
    else:
        lp = lp.resolve()

    if not lp.exists():
        raise FileNotFoundError(f"Log file not found: {lp}")
    return lp


def _find_scenario_dir(base_dir: Path, sim_index: int) -> Path:
    """
    Locate the folder for *sim_index* under *base_dir*.

    Simulation runs (Agramon_withLAI_withETp.py) name scenario folders
    ``scenario_<sim_index>`` (older runs, before this tagging existed) or
    ``scenario_<sim_index>_<tag>``, where <tag> describes the dem-plot /
    outlet / bottom / spin-up combination used for that run — see
    build_scenario_dirname() in Agramon_withLAI_withETp.py. Matching is
    anchored on sim_index as a whole underscore-delimited segment so e.g.
    sim_index=2 never matches a folder belonging to sim_index=20.

    Raises FileNotFoundError if nothing (or, after a warning, more than
    one candidate) matches.
    """
    if not base_dir.exists():
        raise FileNotFoundError(f"{base_dir} does not exist.")

    pattern = re.compile(rf"^scenario_{sim_index}(_.+)?$")
    matches = sorted(
        p for p in base_dir.iterdir() if p.is_dir() and pattern.match(p.name)
    )

    if not matches:
        available = sorted(p.name for p in base_dir.iterdir() if p.is_dir())
        raise FileNotFoundError(
            f"No folder matching 'scenario_{sim_index}' (optionally "
            f"'scenario_{sim_index}_<tag>') found under {base_dir}.\n"
            f"  Available: {available}"
        )
    if len(matches) > 1:
        print(f"  WARNING: multiple folders match scenario_{sim_index} under "
              f"{base_dir} — using '{matches[0].name}'. "
              f"All matches: {[m.name for m in matches]}")
    return matches[0]


def _load_log_row(path2prj: Path, log_file: str, scenario: int | None) -> tuple[int, pd.Series]:
    """
    Return (sim_index, log_row) for the requested withLAI scenario.

    If *scenario* is None, the first row where ``withLAI == 1`` is used.
    Raises ValueError when no matching row is found.
    """
    log_path = _resolve_log(path2prj, log_file)
    # log_path = _resolve_log(path2prj/args.log_file)

    df_log   = pd.read_csv(log_path)

    # Filter to withLAI=1 runs when the column exists
    lai_col = "withLAI" if "withLAI" in df_log.columns else None
    df_lai  = df_log[df_log[lai_col] == 1] if lai_col else df_log

    if df_lai.empty:
        raise ValueError(
            f"No withLAI=1 rows found in {log_path}. "
            "Run LT_SShydro_run_withLAI.py first."
        )

    sim_index = scenario if scenario is not None else int(df_lai["sim_index"].iloc[0])

    row = df_lai[df_lai["sim_index"] == sim_index]
    if row.empty:
        raise ValueError(
            f"Scenario {sim_index} not found (or not withLAI=1) in {log_path}. "
            f"Available withLAI scenarios: {list(df_lai['sim_index'])}"
        )
    return sim_index, row.iloc[0]



# =============================================================================
# Shared utilities (mirrors analyse_ET.py helpers)
# =============================================================================

def load_shapefiles() -> gpd.GeoDataFrame:
    gdf = AgUtils.load_plot_shapefiles(['microcuencas_13'],
                                       MODULE_PATH / "shapefiles")
    if gdf.crs is None:
        raise ValueError("Shapefile has no CRS.")
    if str(gdf.crs) != TARGET_CRS:
        gdf = gdf.to_crs(TARGET_CRS)
    print(f"  Shapefile CRS : {gdf.crs}")
    for i, pid in enumerate(gdf["PlotID"]):
        print(f"    [{i}] {pid}")
    return gdf

def resolve_plot_selection(
    gdf: gpd.GeoDataFrame,
    selection: list[str | int] | None,
) -> gpd.GeoDataFrame:
    if selection is None or selection == ["all"]:
        return gdf.copy()

    # Ensure selection is a list
    if not isinstance(selection, list):
        selection = [selection]

    resolved_indices = []
    plot_ids = list(gdf["PlotID"])

    for item in selection:
        if isinstance(item, int) or (isinstance(item, str) and item.isdigit()):
            idx = int(item)
            if not (0 <= idx < len(gdf)):
                raise IndexError(f"Index {idx} out of range.")
            resolved_indices.append(idx)  # Store DataFrame index, not PlotID value
        else:
            if item not in plot_ids:
                raise ValueError(f"PlotID '{item}' not found. Available: {plot_ids}")
            matches = gdf[gdf["PlotID"] == item].index.tolist()
            resolved_indices.extend(matches)

    # Deduplicate indices
    seen: set[int] = set()
    unique_indices = [i for i in resolved_indices if not (i in seen or seen.add(i))]

    return gdf.iloc[unique_indices].copy().reset_index(drop=True)

def selection_label(gdf: gpd.GeoDataFrame, selection: list[str | int] | None) -> str:
    if selection is None or selection == ["all"]:
        return "all"

    if not isinstance(selection, list):
        selection = [selection]

    names: list[str] = []
    for item in selection:
        if isinstance(item, int) or (isinstance(item, str) and item.isdigit()):
            names.append(str(int(item)))  # Use the index directly in the label
        else:
            names.append(item)

    seen: set[str] = set()
    unique = [n for n in names if not (n in seen or seen.add(n))]
    return "_".join(unique)


def resolve_plot_for_dataset(gdf: gpd.GeoDataFrame, ds: xr.Dataset) -> gpd.GeoDataFrame:
    """Select the shapefile row(s) that actually overlap *ds*'s own extent.

    Used as a fallback when the user hasn't explicitly picked a plot via
    --plots: each scenario's WB grid (ds_wb) is already cropped to a single
    dem-plot by the solver, and that plot doesn't necessarily correspond to
    shapefile row 0 (the --plots default). Clipping ds_wb against an
    unrelated microcuenca polygon finds zero overlap and rioxarray raises
    NoDataInBounds — so instead of trusting row 0, pick whichever row(s)
    intersect ds's bounding box.
    """
    bounds = ds.rio.bounds()
    hits = gdf[gdf.intersects(box(*bounds))]
    if hits.empty:
        raise ValueError(
            f"No shapefile polygon overlaps dataset bounds {bounds} — "
            "check CRS / plot registration."
        )
    return hits.reset_index(drop=True)


def _clip_ds(ds: xr.Dataset, gdf: gpd.GeoDataFrame) -> xr.Dataset:
    clipped: dict[str, xr.DataArray] = {}
    for var in ds.data_vars:
        da = (
            ds[var]
            .rio.write_crs(TARGET_CRS, inplace=False)
            .rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
        )
        clipped[var] = da.rio.clip(gdf.geometry, gdf.crs, drop=True, all_touched=True)
    return xr.Dataset(clipped, attrs=ds.attrs)


# =============================================================================
# EO-EB dataset loader  (mirrors analyse_ET._load_nc)
# =============================================================================

def load_eo_et(path: Path = EO_ET_PATH) -> xr.Dataset:
    """Load EO ETa NetCDF, reproject to TARGET_CRS.  Returns dataset in mm/day."""
    ds = (
        xr.open_dataset(path, decode_times=True, mask_and_scale=True)
          .drop_vars("spatial_ref", errors="ignore")
    )
    ds = ds.rio.write_crs(NATIVE_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    if str(ds.rio.crs) != TARGET_CRS:
        ds = ds.rio.reproject(TARGET_CRS)
    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  EO-ET CRS : {ds.rio.crs}   vars: {list(ds.data_vars)}")
    return ds


def load_lai_ref(path: Path = LAI_PATH) -> xr.Dataset:
    """Load the LAI NetCDF used to build vegetation zones, reproject to
    TARGET_CRS. Mirrors load_eo_et(); returns None (with a warning) if
    *path* doesn't exist so callers can skip the zone figures gracefully
    instead of crashing the whole comparison run."""
    if not path.exists():
        print(f"  WARNING: LAI file not found at {path} — "
              "per-zone comparison will be skipped. Set LAI_PATH.")
        return None
    ds = (
        xr.open_dataset(path, decode_times=True, mask_and_scale=True)
          .drop_vars("spatial_ref", errors="ignore")
    )
    ds = ds.rio.write_crs(NATIVE_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    if str(ds.rio.crs) != TARGET_CRS:
        ds = ds.rio.reproject(TARGET_CRS)
    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  LAI CRS : {ds.rio.crs}   vars: {list(ds.data_vars)}")
    return ds


# =============================================================================
# WB dataset loader  (reads et_output.nc from simulation artefacts)
# =============================================================================

def _find_sim_dir(scenario: int, log_file: str | None = None) -> Path:
    """
    Locate the WB simulation output folder for *scenario*.

    Mirrors withLAI_results.py's scenario-selection pattern instead of
    assuming a flat <WB_PRJ_ROOT>/sim<N>/outputs/ layout (which is not how
    LT_SShydro_run_withLAI.py actually writes its outputs):
      1. Resolve *scenario* to a sim_index via the shared log
         (<WB_PRJ_ROOT>/../Codigos_HydroModel/<log_file>).
      2. Find the matching (possibly tagged) scenario_<sim_index>[_tag]
         folder under <WB_PRJ_ROOT>/outputs/.

    *log_file* defaults to the module-level LOG_FILE, looked up here
    (rather than as a ``= LOG_FILE`` parameter default) so that a
    ``--log-file`` CLI override — which reassigns the LOG_FILE global in
    ``__main__`` — actually takes effect. A default bound at function
    -definition time would freeze on the module's original value.
    """
    if log_file is None:
        log_file = LOG_FILE
    sim_index, log_row = _load_log_row(WB_PRJ_ROOT / "..", log_file, scenario)
    sim_dir = _find_scenario_dir(WB_PRJ_ROOT / "outputs", sim_index)
    print(f"  WB scenario  : sim_index={sim_index}  → {sim_dir.name}")
    return sim_dir


<<<<<<< HEAD
def _safe_open_wb_dataset(nc_path: Path, max_attempts: int = 5,
                           retry_wait_s: float = 2.0) -> xr.Dataset:
    """
    Open et_output.nc defensively while it may still be written by
    LT_SShydro_run_withLAI.py (a separate, unmodified script).

    We never touch the live file with a long-lived handle: instead we
    copy it to a scratch location, load it fully into memory, and
    delete the scratch copy. If the copy lands mid-write (truncated /
    inconsistent HDF5 metadata), opening will raise — in that case we
    wait and retry, since the writer is presumably still going and a
    later copy should be complete.
    """
    import os
    import shutil
    import tempfile
    import time

    last_err: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        tmp_dir = Path(tempfile.mkdtemp(prefix="et_output_snapshot_"))
        tmp_path = tmp_dir / "et_output.nc"
        try:
            # shutil.copy2 can itself race with an in-progress write, but
            # any resulting inconsistency is confined to tmp_path, never
            # the original file.
            shutil.copy2(nc_path, tmp_path)
            with xr.open_dataset(tmp_path, decode_times=True,
                                  mask_and_scale=True,
                                  engine="netcdf4") as ds_tmp:
                # .load() pulls all data into memory now, so we can close
                # and remove tmp_path immediately afterwards regardless
                # of what happens to the original file next.
                return ds_tmp.load()
        except Exception as exc:  # noqa: BLE001 - deliberately broad: any
            # failure here means "this snapshot wasn't a clean read",
            # not a fatal problem with the pipeline.
            last_err = exc
            print(f"  [snapshot read attempt {attempt}/{max_attempts} failed: "
                  f"{exc!r} — retrying in {retry_wait_s}s]")
            time.sleep(retry_wait_s)
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    raise RuntimeError(
        f"Could not get a clean read of {nc_path} after {max_attempts} "
        f"attempts (last error: {last_err!r}). The writer process may "
        f"still be running — try again once it finishes."
    ) from last_err


=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
def load_wb_et(scenario: int = 0) -> tuple[xr.Dataset, Path]:
    """
    Load et_output.nc from a CATHY simulation, normalise to mm/day,
    and align spatial dims so rioxarray can clip it.

    The raw file has dims (X, Y, datetime) and values in m/s.
    Returns (ds, sim_dir):
      • ds      — Dataset with variable name preserved as WB_ET_VAR
                  ('ACT. ETRA'), values scaled to mm/day, dims renamed
                  to (x, y, time), CRS written as TARGET_CRS (25830)
      • sim_dir — the resolved (possibly tagged) scenario_<N>[_tag]
                  folder the data came from, e.g.
                  scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c
                  — handed back so callers can mirror it as a figures
                  subfolder name (see __main__).
    """
    sim_dir = _find_sim_dir(scenario)
    nc_path = sim_dir / "et_output.nc"
    if not nc_path.exists():
        raise FileNotFoundError(f"et_output.nc not found in {sim_dir}")
    print(f"Attempting to open file at: {nc_path}")
    ds = _safe_open_wb_dataset(nc_path)
    print(f"  WB-ET raw dims : {dict(ds.dims)}   vars: {list(ds.data_vars)}")

    if WB_ET_VAR not in ds:
        raise KeyError(
            f"Variable '{WB_ET_VAR}' not in WB dataset. "
            f"Available: {list(ds.data_vars)}"
        )

    # Scale m/s → mm/day
    ds = ds.assign({WB_ET_VAR: ds[WB_ET_VAR] * WB_ET_SCALE})

    # Normalise dimension names → (x, y, time)
    rename_map: dict[str, str] = {}
    if "X" in ds.dims:
        rename_map["X"] = "x"
    if "Y" in ds.dims:
        rename_map["Y"] = "y"
    if "datetime" in ds.dims and "time" not in ds.dims:
        # Drop any existing 'time' coord OR data variable that would
        # conflict with renaming the 'datetime' dim to 'time'. The raw
        # CATHY et_output.nc carries its own 'time' data var (seconds/
        # timestamp helper column) alongside the 'datetime' dim — using
        # ds.variables (not just ds.coords) catches that case too.
        if "time" in ds.variables:
            ds = ds.drop_vars("time")
        rename_map["datetime"] = "time"
    if rename_map:
        ds = ds.rename(rename_map)

    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  WB-ET final dims : {dict(ds.dims)}   CRS: {ds.rio.crs}")
    return ds, sim_dir


# =============================================================================
# Build common monthly spatial-mean time series
# =============================================================================

def _spatial_mean_monthly(ds: xr.Dataset, var: str, time_dim: str = "time") -> pd.Series:
    """Spatial mean (x,y) → resample to monthly mean → pd.Series."""
    sp_dims = [d for d in ds[var].dims if d != time_dim]
    da_mean = ds[var].mean(dim=sp_dims, skipna=True)
    df = pd.DataFrame(
        {var: da_mean.values},
        index=pd.DatetimeIndex(ds[time_dim].values),
    )
    return df[var].resample("ME").mean()


def build_common_series(
    ds_eo: xr.Dataset, var_eo: str,
    ds_wb: xr.Dataset,
    s_eo_override: pd.Series | None = None,
) -> tuple[pd.Series, pd.Series]:
    """
    Return (s_eo, s_wb): coincident monthly spatial-mean ETa series (mm/day).
    Only months present in BOTH datasets are kept. *s_eo_override* replaces
    the NetCDF-derived reference (e.g. the Agramon CSV real-ET series).
    """
    s_eo = (s_eo_override if s_eo_override is not None
            else _spatial_mean_monthly(ds_eo, var_eo, time_dim="time"))
    s_wb = _spatial_mean_monthly(ds_wb, WB_ET_VAR, time_dim="time")

    # Align to the period common to both
    idx_common = s_eo.index.intersection(s_wb.index)
    if len(idx_common) == 0:
        raise ValueError(
            "EO-EB and WB datasets share no overlapping months.\n"
            f"EO range: {s_eo.index[0]} – {s_eo.index[-1]}\n"
            f"WB range: {s_wb.index[0]} – {s_wb.index[-1]}"
        )
    print(f"  Common period : {idx_common[0].date()} – {idx_common[-1].date()}  "
          f"({len(idx_common)} months)")
    return s_eo.loc[idx_common], s_wb.loc[idx_common]


# =============================================================================
# LAI-zone classification & masking
# =============================================================================
# Mechanism check: does WB-vs-EO-EB disagreement concentrate in the
# high-LAI zones pre-fire, and do both series funnel toward the
# bare-soil-zone behaviour post-fire almost everywhere? Domain-mean series
# (build_common_series above) can't show this — it averages the effect
# away. The functions below split ETa into fixed vegetation zones (based
# on PRE-FIRE mean LAI) and rebuild the EO-EB/WB comparison inside each.

def compute_lai_zone_reference(
    da_lai:        xr.DataArray,
    reference_end: pd.Timestamp = FIRE_DATE,
    time_dim:      str = "time",
) -> xr.DataArray:
    """
    Time-mean LAI map using only records before *reference_end* (default:
    FIRE_DATE). This static, pre-disturbance map is what zone boundaries
    are drawn from — so a pixel keeps its zone label even after the fire
    changes its actual LAI.
    """
    mask = da_lai[time_dim] < pd.Timestamp(reference_end)
    if not bool(mask.any()):
        raise ValueError(
            f"No LAI records before {pd.Timestamp(reference_end).date()} "
            "— cannot build a pre-fire zone reference."
        )
    return da_lai.sel({time_dim: mask}).mean(time_dim, skipna=True)


def classify_lai_zones(
    da_lai_ref: xr.DataArray,
    n_zones:    int = N_LAI_ZONES,
) -> tuple[xr.DataArray, list[float]]:
    """
    Classify a static LAI reference map into *n_zones* classes by
    within-domain quantiles (terciles by default): zone 1 = lowest LAI
    ('Bare'), zone n_zones = highest LAI ('Dense veg.').

    Returns (da_zone, bounds):
      da_zone : same shape as da_lai_ref, values in {1, ..., n_zones}
                (NaN where da_lai_ref is NaN)
      bounds  : the (n_zones - 1) percentile thresholds used
    """
    vals = da_lai_ref.values
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        raise ValueError("LAI reference map has no finite values to classify.")

    qs     = np.linspace(0, 100, n_zones + 1)[1:-1]
    bounds = [float(b) for b in np.percentile(finite, qs)]
    edges  = [-np.inf] + bounds + [np.inf]

    da_zone = xr.full_like(da_lai_ref, np.nan, dtype=float)
    for z in range(1, n_zones + 1):
        m = (da_lai_ref >= edges[z - 1]) & (da_lai_ref < edges[z])
        da_zone = xr.where(m, z, da_zone)
    da_zone = da_zone.rename("lai_zone")
    print(f"  LAI zone bounds (m²/m²): {['%.2f' % b for b in bounds]}")
    return da_zone, bounds


def _zone_mask_on_grid(
    da_zone_ref: xr.DataArray,
    ds_target:   xr.Dataset,
    var_target:  str,
) -> xr.DataArray:
    """
    Resample the categorical zone map onto *ds_target*'s own (x, y) grid
    via nearest-neighbour matching, so the same pre-fire zone definition
    can mask both the EO-EB grid and the (differently-resolved) WB grid.
    """
    from rasterio.enums import Resampling

    target_da = ds_target[var_target]
    if "time" in target_da.dims:
        target_da = target_da.isel(time=0, drop=True)
    target_da = (
        target_da.rio.write_crs(TARGET_CRS, inplace=False)
                  .rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    )
    zone_src = (
        da_zone_ref.rio.write_crs(TARGET_CRS, inplace=False)
                   .rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    )
    return zone_src.rio.reproject_match(target_da, resampling=Resampling.nearest)


def _spatial_mean_monthly_masked(
    ds:       xr.Dataset,
    var:      str,
    zone_da:  xr.DataArray,
    zone_id:  int,
    time_dim: str = "time",
) -> pd.Series:
    """Like _spatial_mean_monthly, but restricted to pixels where
    zone_da == zone_id."""
    da_masked = ds[var].where(zone_da == zone_id)
    sp_dims   = [d for d in da_masked.dims if d != time_dim]
    da_mean   = da_masked.mean(dim=sp_dims, skipna=True)
    df = pd.DataFrame(
        {var: da_mean.values},
        index=pd.DatetimeIndex(ds[time_dim].values),
    )
    return df[var].resample("ME").mean()


def build_zone_series(
    ds_eo:      xr.Dataset,
    var_eo:     str,
    ds_wb:      xr.Dataset,
    da_lai_ref: xr.DataArray,
    n_zones:    int = N_LAI_ZONES,
) -> tuple[dict[int, tuple[pd.Series, pd.Series]], list[float]]:
    """
    Build coincident monthly spatial-mean EO-EB/WB series per LAI zone.

    Returns (zone_series, bounds) where zone_series maps
    zone_id -> (s_eo_zone, s_wb_zone), each pair aligned to their common
    months (same convention as build_common_series).
    """
    da_zone_ref, bounds = classify_lai_zones(da_lai_ref, n_zones=n_zones)
    zone_on_eo = _zone_mask_on_grid(da_zone_ref, ds_eo, var_eo)
    zone_on_wb = _zone_mask_on_grid(da_zone_ref, ds_wb, WB_ET_VAR)

    zone_series: dict[int, tuple[pd.Series, pd.Series]] = {}
    for z in range(1, n_zones + 1):
        s_eo_z = _spatial_mean_monthly_masked(ds_eo, var_eo, zone_on_eo, z)
        s_wb_z = _spatial_mean_monthly_masked(ds_wb, WB_ET_VAR, zone_on_wb, z)
        idx_common = s_eo_z.index.intersection(s_wb_z.index)
        if len(idx_common) == 0:
            print(f"  Zone {z}: no overlapping months (zone may be empty "
                  "in one of the two grids) — skipped.")
            continue
        zone_series[z] = (s_eo_z.loc[idx_common], s_wb_z.loc[idx_common])
    return zone_series, bounds


def print_zone_bias_summary(
    zone_series: dict[int, tuple[pd.Series, pd.Series]],
    zone_labels: tuple[str, ...] = ZONE_LABELS,
) -> None:
    """
    Print mean (EO-EB − WB) bias per zone, pre- vs post-fire — the direct
    numeric check for "disagreement concentrates in high-LAI zones
    pre-fire, converges post-fire": a shrinking |Δ| for the pre-fire
    high-LAI zone, converging toward the (typically small) bare-soil-zone
    bias, is the signature to look for.
    """
    print("\n  Per-zone EO-EB − WB bias (mm/day):")
    print(f"  {'Zone':<16}{'Pre-fire':>12}{'Post-fire':>12}{'Δ post-pre':>14}")
    for z in sorted(zone_series):
        s_eo, s_wb = zone_series[z]
        bias = s_eo - s_wb
        pre  = bias[bias.index < FIRE_DATE].mean()
        post = bias[bias.index >= FIRE_DATE].mean()
        delta = post - pre
        zl = zone_labels[z - 1] if z - 1 < len(zone_labels) else f"Zone {z}"
        print(f"  {zl:<16}{pre:>12.3f}{post:>12.3f}{delta:>14.3f}")


# =============================================================================
# Figure A — monthly mean by year (EO vs WB on same axes)
# =============================================================================

def plot_compare_monthly_by_year(
    s_eo:     pd.Series,
    s_wb:     pd.Series,
    fig_path: Path,
    label:    str  = "",
    dpi:      int  = 150,
    show:     bool = False,
) -> None:
    """
    Spatial-mean monthly ETa (mm/day), one line per year.
    Solid lines = EO-EB, dashed lines = WB.  Matching colours per year.
    """
    years = sorted(set(s_eo.index.year) | set(s_wb.index.year))
    cmap_lines = plt.get_cmap("tab10", len(years))

    fig, ax = plt.subplots(figsize=(10, 5))
    months = np.arange(1, 13)
    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    for i, year in enumerate(years):
        color = cmap_lines(i)

        def _monthly_profile(s: pd.Series) -> pd.Series:
            sub = s[s.index.year == year]
            return sub.groupby(sub.index.month).mean()

        eo_m = _monthly_profile(s_eo)
        wb_m = _monthly_profile(s_wb)

        if not eo_m.empty:
            ax.plot(eo_m.index, eo_m.values,
                    marker="o", lw=2, color=color,
                    label=f"{year} {EO_LABEL}")
        if not wb_m.empty:
            ax.plot(wb_m.index, wb_m.values,
                    marker="s", lw=2, ls="--", color=color,
                    label=f"{year} WB")

    ax.set_xticks(months)
    ax.set_xticklabels(month_labels)
    ax.set_xlabel("Month")
    ax.set_ylabel("ETa (mm/day)")
    ax.set_title(f"Spatial-mean monthly ETa — {EO_LABEL} (solid) vs WB (dashed)",
                 fontweight="bold")
    ax.legend(ncol=2, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_monthly_mean_by_year{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure B — annual total bar chart (grouped EO vs WB)
# =============================================================================

def plot_compare_annual_total(
    s_eo:     pd.Series,
    s_wb:     pd.Series,
    fig_path: Path,
    label:    str  = "",
    dpi:      int  = 150,
    show:     bool = False,
) -> None:
    """Grouped annual bar chart: EO-EB (blue) and WB (orange)."""
    # Monthly series are mean rates (mm/day): weight by days in month before
    # summing, otherwise "mm/year" is ~12 x the mean daily rate (~1/30 of the
    # true annual total).
    annual_eo = (s_eo * s_eo.index.days_in_month).resample("YE").sum()
    annual_wb = (s_wb * s_wb.index.days_in_month).resample("YE").sum()

    # Align to common years
    idx = annual_eo.index.intersection(annual_wb.index)
    annual_eo = annual_eo.loc[idx]
    annual_wb = annual_wb.loc[idx]
    years = annual_eo.index.year

    x   = np.arange(len(years))
    w   = 0.35
    fig, ax = plt.subplots(figsize=(max(8, len(years) * 1.1), 5))

    bars_eo = ax.bar(x - w / 2, annual_eo.values, w, label=EO_LABEL,
                     color=COLOR_EO, edgecolor="white")
    bars_wb = ax.bar(x + w / 2, annual_wb.values, w, label="WB",
                     color=COLOR_WB, edgecolor="white")

    for bars in (bars_eo, bars_wb):
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + max(annual_eo.max(), annual_wb.max()) * 0.01,
                    f"{h:.0f}", ha="center", va="bottom", fontsize=8)

    if FIRE_DATE.year in list(years):
        xf = (list(years).index(FIRE_DATE.year)
              + (FIRE_DATE.dayofyear - 1) / 365.0 - 0.5)
        ax.axvline(xf, color="#B22222", linestyle="--", linewidth=1.3,
                   alpha=0.85, zorder=4, label="Fire (Jul 2020)")
    ax.set_xticks(x)
    ax.set_xticklabels(years)
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual ETa sum (mm/year)")
    ax.set_title(f"Annual domain-mean ETa — {EO_LABEL} vs Water-Balance model",
                 fontweight="bold")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_annual_total{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure C — monthly total bar chart (grouped EO vs WB)
# =============================================================================

def plot_compare_monthly_total(
    s_eo:     pd.Series,
    s_wb:     pd.Series,
    fig_path: Path,
    label:    str  = "",
    dpi:      int  = 150,
    show:     bool = False,
) -> None:
    """Grouped monthly bar chart for the full common record."""
    time_labels = s_eo.index.strftime("%Y-%m")
    x = np.arange(len(s_eo))
    w = 0.35

    fig, ax = plt.subplots(figsize=(max(12, len(s_eo) * 0.22), 5))
    ax.bar(x - w / 2, s_eo.values, w, label=EO_LABEL, color=COLOR_EO, edgecolor="none")
    ax.bar(x + w / 2, s_wb.values, w, label="WB",    color=COLOR_WB, edgecolor="none")

    # Year-boundary ticks only
    jan_pos    = [i for i, lbl in enumerate(time_labels) if lbl.endswith("-01")]
    jan_labels = [lbl[:4] for lbl in time_labels if lbl.endswith("-01")]
    ax.set_xticks(jan_pos)
    ax.set_xticklabels(jan_labels, rotation=0, ha="center", fontsize=9)
    for pos in jan_pos:
        ax.axvline(pos - 0.5 - w / 2, color="gray", lw=0.5, alpha=0.4, zorder=0)

    _fire_lbl = f"{FIRE_DATE:%Y-%m}"
    if _fire_lbl in list(time_labels):
        ax.axvline(list(time_labels).index(_fire_lbl) - 0.5, color="#B22222",
                   linestyle="--", linewidth=1.3, alpha=0.85, zorder=4,
                   label="Fire (Jul 2020)")
    ax.set_ylabel("Monthly mean ETa (mm/day)")
    ax.set_title(f"Monthly domain-mean ETa — {EO_LABEL} vs Water-Balance model",
                 fontweight="bold")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_monthly_total{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure D — scatter EO-EB vs WB (coincident monthly means)
# =============================================================================

def plot_compare_scatter(
    s_eo:     pd.Series,
    s_wb:     pd.Series,
    fig_path: Path,
    label:    str  = "",
    dpi:      int  = 150,
    show:     bool = False,
) -> None:
    """Scatter of coincident monthly spatial-mean ETa, coloured by year."""
    years      = s_eo.index.year
    uniq_years = sorted(set(years))
    cmap_sc    = plt.get_cmap("tab10", len(uniq_years))
    year_color = {y: cmap_sc(i) for i, y in enumerate(uniq_years)}

    fig, ax = plt.subplots(figsize=(6, 6))
    for year in uniq_years:
        mask = years == year
        ax.scatter(s_eo[mask], s_wb[mask], color=year_color[year],
                   alpha=0.75, s=40, label=str(year))

    # 1:1 line
    all_vals = np.concatenate([s_eo.values, s_wb.values])
    lim = (0, np.nanmax(all_vals) * 1.05)
    ax.plot(lim, lim, "k--", lw=1, alpha=0.5, label="1:1")

    # Linear regression
    finite = np.isfinite(s_eo.values) & np.isfinite(s_wb.values)
    if finite.sum() > 2:
        m, b = np.polyfit(s_eo.values[finite], s_wb.values[finite], 1)
        corr = np.corrcoef(s_eo.values[finite], s_wb.values[finite])[0, 1]
        rmse = np.sqrt(np.mean((s_eo.values[finite] - s_wb.values[finite]) ** 2))
        xs   = np.array(lim)
        ax.plot(xs, m * xs + b, color="#555", lw=1.5,
                label=f"Reg.  R={corr:.2f}  RMSE={rmse:.2f}")
        ax.text(0.05, 0.93,
                f"R = {corr:.3f}\nRMSE = {rmse:.3f} mm/day\ny = {m:.2f}x + {b:.2f}",
                transform=ax.transAxes, fontsize=9,
                verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel(f"{EO_LABEL} ETa (mm/day)")
    ax.set_ylabel("WB ETa (mm/day)")
    ax.set_title(f"Monthly spatial-mean ETa: {EO_LABEL} vs Water-Balance",
                 fontweight="bold")
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_aspect("equal")
    ax.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_scatter{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure E — spatial monthly maps for a selected year (EO and WB, shared cscale)
# =============================================================================

def _spatial_maps_mosaic(
    ds:           xr.Dataset,
    var:          str,
    time_dim:     str,
    selected_year: int,
    vmin:         float,
    vmax:         float,
    cmap:         str,
    title_prefix: str,
    fig_path:     Path,
    fname:        str,
    dpi:          int  = 150,
    show:         bool = False,
) -> None:
    """Generic monthly spatial mosaic for one dataset."""
    da_year = ds[var].sel(
        {time_dim: ds[time_dim].dt.year == selected_year}
    )
    if da_year.sizes[time_dim] == 0:
        print(f"  {title_prefix}: no data for {selected_year} — skipped.")
        return

    months   = np.unique(da_year[time_dim].dt.month.values)
    n_months = len(months)
    ncols    = min(4, n_months)
    nrows    = int(np.ceil(n_months / ncols))

    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    sp_dims = [d for d in da_year.dims if d != time_dim]
    x_dim, y_dim = sp_dims[0], sp_dims[1]

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows),
                             squeeze=False)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap_obj = plt.get_cmap(cmap)

    for idx, month in enumerate(months):
        ax = axes[idx // ncols][idx % ncols]
        da_m = da_year.sel(
            {time_dim: da_year[time_dim].dt.month == month}
        ).mean(time_dim)

        im = ax.imshow(
            da_m.values,
            origin="upper",
            cmap=cmap_obj,
            norm=norm,
            extent=[float(da_m.x.min()), float(da_m.x.max()),
                    float(da_m.y.min()), float(da_m.y.max())],
            aspect="auto",
        )
        ax.set_title(month_names[month - 1], fontsize=9)
        ax.set_xlabel(x_dim, fontsize=7)
        ax.set_ylabel(y_dim, fontsize=7)
        ax.tick_params(labelsize=7)
        fig.colorbar(im, ax=ax, shrink=0.8).set_label("mm/day", fontsize=7)

    for idx in range(n_months, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(f"{title_prefix} — monthly ETa maps {selected_year}",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    out = fig_path / fname
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


def plot_compare_spatial_maps(
    ds_eo:         xr.Dataset,
    var_eo:        str,
    ds_wb:         xr.Dataset,
    selected_year: int,
    fig_path:      Path,
    label:         str  = "",
    dpi:           int  = 150,
    show:          bool = False,
) -> None:
    """
    Produce two spatial-map mosaics for *selected_year* (EO-EB and WB)
    using a shared colour scale so maps can be compared directly.
    """
    suffix = f"_{label}" if label else ""

    # Determine shared colour limits from both datasets for that year
    def _year_vals(ds: xr.Dataset, var: str, td: str) -> np.ndarray:
        mask = ds[td].dt.year == selected_year
        da   = ds[var].sel({td: mask})
        v    = da.values.ravel()
        return v[np.isfinite(v)]

    td_eo = "time"
    td_wb = "time"

    v_eo = _year_vals(ds_eo, var_eo, td_eo)
    v_wb = _year_vals(ds_wb, WB_ET_VAR, td_wb)

    if v_eo.size == 0 and v_wb.size == 0:
        print(f"  No data for {selected_year} in either dataset — spatial maps skipped.")
        return

    all_vals = np.concatenate([v for v in (v_eo, v_wb) if v.size > 0])
    vmin     = float(np.nanpercentile(all_vals, 2))
    vmax     = float(np.nanpercentile(all_vals, 98))

    _spatial_maps_mosaic(
        ds_eo, var_eo, td_eo, selected_year, vmin, vmax, CMAP_EO,
        "EO-EB ETa", fig_path,
        f"compare_ET_spatial_EO_{selected_year}{suffix}.png",
        dpi=dpi, show=show,
    )
    _spatial_maps_mosaic(
        ds_wb, WB_ET_VAR, td_wb, selected_year, vmin, vmax, CMAP_WB,
        "WB ETa", fig_path,
        f"compare_ET_spatial_WB_{selected_year}{suffix}.png",
        dpi=dpi, show=show,
    )


# =============================================================================
# Figure F — monthly mean by year, faceted per LAI zone
# =============================================================================

def plot_compare_monthly_by_year_per_zone(
    zone_series: dict[int, tuple[pd.Series, pd.Series]],
    fig_path:    Path,
    label:       str = "",
    zone_labels: tuple[str, ...] = ZONE_LABELS,
    dpi:         int  = 150,
    show:        bool = False,
) -> None:
    """
    One panel per LAI zone, each reproducing Fig A (spatial-mean monthly
    ETa, one line per year, EO-EB solid / WB dashed) restricted to that
    zone's pixels. Lets disagreement that's invisible in the domain mean
    (Fig A) show up if it's concentrated in a particular zone.
    """
    zones = sorted(zone_series)
    if not zones:
        print("  plot_compare_monthly_by_year_per_zone: no zones with data — skipped.")
        return

    fig, axes = plt.subplots(1, len(zones), figsize=(5.5 * len(zones), 5),
                             sharey=True, squeeze=False)
    axes = axes[0]
    months       = np.arange(1, 13)
    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    for zi, z in enumerate(zones):
        ax = axes[zi]
        s_eo, s_wb = zone_series[z]
        years = sorted(set(s_eo.index.year) | set(s_wb.index.year))
        cmap_lines = plt.get_cmap("tab10", max(len(years), 1))

        for i, year in enumerate(years):
            color = cmap_lines(i)
            eo_y = s_eo[s_eo.index.year == year]
            wb_y = s_wb[s_wb.index.year == year]
            eo_m = eo_y.groupby(eo_y.index.month).mean()
            wb_m = wb_y.groupby(wb_y.index.month).mean()
            if not eo_m.empty:
                ax.plot(eo_m.index, eo_m.values, marker="o", lw=1.6, color=color,
                        label=f"{year} EO-EB" if zi == 0 else None)
            if not wb_m.empty:
                ax.plot(wb_m.index, wb_m.values, marker="s", lw=1.6, ls="--", color=color,
                        label=f"{year} WB" if zi == 0 else None)

        ax.set_xticks(months)
        ax.set_xticklabels(month_labels, fontsize=7)
        zl = zone_labels[z - 1] if z - 1 < len(zone_labels) else f"Zone {z}"
        ax.set_title(zl, fontsize=10, fontweight="bold")
        ax.set_xlabel("Month")
        ax.grid(True, linestyle="--", alpha=0.3)

    axes[0].set_ylabel("ETa (mm/day)")
    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_, ncol=min(len(labels_), 6), fontsize=7,
               loc="upper center", bbox_to_anchor=(0.5, 1.1))
    fig.suptitle("Spatial-mean monthly ETa by LAI zone — EO-EB (solid) vs WB (dashed)",
                 fontsize=12, fontweight="bold", y=1.18)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_monthly_mean_by_year_per_zone{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure G — annual bias per LAI zone (mechanism diagnostic)
# =============================================================================

def plot_compare_annual_bias_per_zone(
    zone_series: dict[int, tuple[pd.Series, pd.Series]],
    fig_path:    Path,
    label:       str = "",
    zone_labels: tuple[str, ...] = ZONE_LABELS,
    dpi:         int  = 150,
    show:        bool = False,
) -> None:
    """
    Annual mean bias (EO-EB minus WB), one line per LAI zone, fire date
    marked. The mechanism signature to look for: the pre-fire high-LAI
    zone's |bias| shrinks after the fire and converges toward the
    (typically smaller) bare-soil-zone bias, rather than all zones moving
    in lockstep.
    """
    zones = sorted(zone_series)
    if not zones:
        print("  plot_compare_annual_bias_per_zone: no zones with data — skipped.")
        return

    fig, ax = plt.subplots(figsize=(9, 5))
    zone_cmap = plt.get_cmap("viridis", len(zones))

    idx_all = []
    for zi, z in enumerate(zones):
        s_eo, s_wb = zone_series[z]
        bias = (s_eo - s_wb).resample("YE").mean()
        idx_all.append(bias.index)
        zl = zone_labels[z - 1] if z - 1 < len(zone_labels) else f"Zone {z}"
        ax.plot(bias.index, bias.values, marker="o", lw=2, color=zone_cmap(zi), label=zl)

    ax.axhline(0, color="black", lw=0.8, alpha=0.6)
    if idx_all:
        d0 = min(idx.min() for idx in idx_all)
        d1 = max(idx.max() for idx in idx_all)
        add_fire_vline(ax, d0, d1)

    ax.set_ylabel("Annual mean bias, EO-EB − WB (mm/day)")
    ax.set_xlabel("Year")
    ax.set_title("Per-zone ETa bias — EO-EB vs Water-Balance model", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.3)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_annual_bias_per_zone{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
<<<<<<< HEAD
# Figure H — simulated WB ETa vs Agramon CSV (real ET, ETp dry / wet bounds)
# =============================================================================

_FORCING_RE = re.compile(r"forcing-etpts-(dry|wet)-fid(\d+)")


def parse_etpts_forcing(sim_dir: Path) -> tuple[str, int] | None:
    """
    Parse (column, feature_id) from a scenario folder name such as
    ``..._forcing-etpts-dry-fid2`` -> ("ET_0_dry", 2). Returns None for runs
    that were not forced with the ETp time series (era5 / satellite etp).
    """
    m = _FORCING_RE.search(sim_dir.name)
    if m is None:
        return None
    return f"ET_0_{m.group(1)}", int(m.group(2))


# Label of the "EO" reference series in Figs A-D. run_comparison() switches it
# from the NetCDF to the Agramon CSV for etp_ts (wet / dry) runs.
EO_LABEL = "EO-EB"


def csv_monthly_series(df_csv: pd.DataFrame, cond: str | None) -> tuple[pd.Series, str]:
    """Monthly-mean real-ET series (mm/day) from the daily Agramon CSV frame
    for the WB forcing condition; returns (series, column_used)."""
    col = real_et_col(df_csv, cond)
    return df_csv[col].dropna().resample("ME").mean().dropna(), col


def forcing_condition(sim_dir: Path) -> str | None:
    """'wet' / 'dry' ETp bound the WB run was forced with (None if not etp_ts)."""
    parsed = parse_etpts_forcing(sim_dir)
    return None if parsed is None else parsed[0].removeprefix("ET_0_")


def real_et_col(df: pd.DataFrame, cond: str | None) -> str:
    """
    Real-ET column to use as reference for a WB run forced with *cond*
    ('wet' | 'dry'): ``ET_<cond>`` when the CSV provides condition-specific
    real-ET columns, otherwise the single real-ET column ``ET`` (the
    agramon_treatments_ET.csv currently has only ``ET``).
    """
    if cond and f"ET_{cond}" in df.columns:
        return f"ET_{cond}"
    return "ET"


def pick_eo_var(ds_eo: xr.Dataset, cond: str | None) -> str:
    """
    EO-EB ET variable matching the WB forcing condition. If the EO file
    carries condition-specific ET variables (ET_wet / ET_dry, ET-wet, ...),
    use the one matching *cond*; otherwise fall back to the first variable.
    """
    names = list(ds_eo.data_vars)
    if cond:
        wanted = {f"et_{cond}", f"et-{cond}", f"et{cond}", f"et_{cond}-gf",
                  f"et-{cond}-gf"}
        for n in names:
            if n.lower() in wanted:
                print(f"  EO-EB variable matched to WB forcing '{cond}': {n}")
                return n
        print(f"  No '{cond}'-specific variable in EO-EB file {names} — "
              f"using '{names[0]}'.")
    return names[0]


def load_csv_et(csv_path: Path, feature_id: int) -> pd.DataFrame | None:
    """Daily ET / ET_0_dry / ET_0_wet (mm/day) of one treatment, indexed by date."""
    if not csv_path.exists():
        print(f"  WARNING: ET CSV not found: {csv_path} — Fig H skipped.")
        return None
    df = pd.read_csv(csv_path, parse_dates=["t"])
    sel = df[df["feature_id"] == feature_id]
    if sel.empty:
        print(f"  WARNING: feature_id={feature_id} not in {csv_path.name} "
              f"(available: {sorted(int(v) for v in df['feature_id'].unique())}) "
              "— Fig H skipped.")
        return None
    keep = [c for c in ("ET", "ET_dry", "ET_wet", "ET_0_dry", "ET_0_wet")
            if c in sel.columns]
    sel = sel.set_index("t").sort_index()[keep]
    sel.index = sel.index.normalize()
    if sel["ET"].notna().sum() == 0:
        print(f"  WARNING: 'ET' is all-NaN for feature_id={feature_id}; "
              "only the ETp bounds will be drawn.")
    return sel


def _spatial_mean_daily(ds: xr.Dataset, var: str, time_dim: str = "time") -> pd.Series:
    """Spatial mean (x, y) -> one value per calendar day (mm/day)."""
    sp_dims = [d for d in ds[var].dims if d != time_dim]
    vals = ds[var].mean(dim=sp_dims, skipna=True).values
    s = pd.Series(vals, index=pd.DatetimeIndex(ds[time_dim].values).normalize())
    return s.groupby(level=0).mean().sort_index()


def _skill(sim: pd.Series, ref: pd.Series) -> dict:
    d = pd.concat([sim, ref], axis=1, join="inner").dropna()
    if len(d) < 3:
        return {"n": len(d)}
    e = d.iloc[:, 0] - d.iloc[:, 1]
    return {"n": len(d), "RMSE": round(float(np.sqrt((e ** 2).mean())), 3),
            "bias": round(float(e.mean()), 3),
            "r": round(float(d.corr().iloc[0, 1]), 3)}


def plot_compare_wb_vs_csv(
    ds_wb:       xr.Dataset,
    s_eo:        pd.Series | None,
    df_csv:      pd.DataFrame,
    forcing_col: str,
    feature_id:  int,
    fig_path:    Path,
    label:       str  = "",
    dpi:         int  = 150,
    show:        bool = False,
) -> None:
    """
    Simulated WB ETa (domain mean, mm/day) against the CSV used as forcing:
    real ET, ETp dry and ETp wet. The bound that forced this run
    (*forcing_col*) is drawn thicker. *s_eo* (monthly EO-EB mean) is added to
    the monthly panel when available.
    """
    sim = _spatial_mean_daily(ds_wb, WB_ET_VAR)
    ref = df_csv.loc[sim.index.min():sim.index.max()]
    # Condition rule: a run forced with the wet (dry) ETp bound is compared
    # ONLY against the wet (dry) reference — matching ETp bound + the real-ET
    # column for that condition; the other bound is not drawn.
    cond    = forcing_col.removeprefix("ET_0_")          # 'wet' | 'dry'
    et_col  = real_et_col(ref, cond)
    print(f"  Period : {sim.index.min():%Y-%m-%d} -> {sim.index.max():%Y-%m-%d} "
          f"({len(sim)} sim days) | feature_id={feature_id} | forcing={forcing_col}")
    print("  Simulated ETa vs CSV columns (mm/day):")
    for col in dict.fromkeys((et_col, forcing_col)):
        print(f"    vs {col:9s}: {_skill(sim, ref[col])}")

    colors = {"ET": COLOR_CSV_ET, "ET_dry": COLOR_CSV_ET, "ET_wet": COLOR_CSV_ET,
              "ET_0_dry": COLOR_CSV_DRY, "ET_0_wet": COLOR_CSV_WET}
    names  = {"ET": "Real ET (CSV)",
              "ET_dry": "Real ET dry (CSV)", "ET_wet": "Real ET wet (CSV)",
              "ET_0_dry": "ETp dry (CSV)", "ET_0_wet": "ETp wet (CSV)"}
    shown  = list(dict.fromkeys((forcing_col, et_col)))   # matching series only

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))

    ax = axes[0]
    for col in shown:
        lw = 1.8 if col == forcing_col else 1.0
        tag = names[col] + (" — forcing" if col == forcing_col else "")
        ax.plot(ref.index, ref[col], color=colors[col], lw=lw, label=tag)
    ax.plot(sim.index, sim.values, color=COLOR_WB, lw=1.2, label="ETa WB (sim)")
    add_fire_vline(ax, sim.index.min(), sim.index.max())
    ax.set_ylabel("mm/day")
    ax.set_title(f"Daily — WB ETa ({cond} ETp forcing) vs Agramon CSV "
                 f"[{cond} reference] (feature_id={feature_id})",
                 fontweight="bold")
    ax.legend(ncol=4, fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.35)

    ax = axes[1]
    mon = pd.concat([sim.rename("ETa WB (sim)"), ref], axis=1).dropna(
        subset=["ETa WB (sim)"]).resample("ME").mean()
    ax.plot(mon.index, mon["ETa WB (sim)"], "-o", ms=3, color=COLOR_WB, label="ETa WB (sim)")
    for col in shown:
        ax.plot(mon.index, mon[col], "--" if col.startswith("ET_0") else "-",
                color=colors[col], lw=1.6 if col == forcing_col else 1.0,
                label=names[col])
    if s_eo is not None and len(s_eo):
        ax.plot(s_eo.index, s_eo.values, "-s", ms=3, color=COLOR_EO,
                label=f"ETa EO-EB ({cond})")
    add_fire_vline(ax, mon.index.min(), mon.index.max())
    ax.set_ylabel("mm/day")
    ax.set_title("Monthly mean", fontweight="bold")
    ax.legend(ncol=3, fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.35)

    ax = axes[2]
    d = pd.concat([ref[et_col].rename("ET"), sim.rename("sim")],
                  axis=1, join="inner").dropna()
    if len(d):
        ax.scatter(d["ET"], d["sim"], s=8, alpha=0.5, color=COLOR_WB)
        lim = [0, max(float(d.max().max()), 1.0)]
        ax.plot(lim, lim, "k--", lw=0.8)
        k = _skill(d["sim"], d["ET"])
        ax.text(0.03, 0.95, f"n={k.get('n')}  RMSE={k.get('RMSE')}  "
                            f"bias={k.get('bias')}  r={k.get('r')}",
                transform=ax.transAxes, va="top", fontsize=9)
    ax.set_xlabel(f"{names[et_col]}, mm/day")
    ax.set_ylabel("ETa WB (sim), mm/day")
    ax.set_title(f"Daily scatter vs {names[et_col].lower()}", fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.35)

    plt.tight_layout()
    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ETa_vs_csv{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    csv_out = fig_path / f"compare_ETa_vs_csv{suffix}.csv"
    pd.concat([sim.rename("ETa_WB_sim"), ref], axis=1).to_csv(csv_out)
    print(f"  Saved → {csv_out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
# Orchestrator
# =============================================================================

def run_comparison(
    ds_eo:         xr.Dataset,
    var_eo:        str,
    ds_wb:         xr.Dataset,
    label:         str,
    selected_year: int = 2020,
    dpi:           int = 150,
    show:          bool = False,
    fig_path:      Path = fig_path,
    da_lai_ref:    xr.DataArray | None = None,
    n_zones:       int = N_LAI_ZONES,
<<<<<<< HEAD
    df_csv:        pd.DataFrame | None = None,
    csv_forcing:   tuple[str, int] | None = None,
=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
) -> None:
    """
    Run all comparison figures. Writes under *fig_path*.

<<<<<<< HEAD
    *df_csv* / *csv_forcing* (daily Agramon CSV frame and the
    (column, feature_id) of the etp_ts run) enable Fig H; None skips it.

=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
    *da_lai_ref* is an optional clipped LAI DataArray (with a 'time' dim)
    used to build the per-zone figures (F, G); pass None to skip them
    (e.g. when LAI_PATH doesn't resolve to a real file).
    """
    fig_path.mkdir(parents=True, exist_ok=True)

    global EO_LABEL
    print("\n── Building common monthly series ───────────────────────────────")
    s_eo_nc, s_wb = build_common_series(ds_eo, var_eo, ds_wb)   # NetCDF EO-EB
    s_eo = s_eo_nc
    EO_LABEL = "EO-EB"
    if df_csv is not None and csv_forcing is not None:
        # etp_ts (wet / dry) run: the EO reference of Figs A-D is the Agramon
        # CSV real-ET column for that condition, not the NetCDF.
        cond = csv_forcing[0].removeprefix("ET_0_")
        s_csv, csv_col = csv_monthly_series(df_csv, cond)
        s_eo, s_wb = build_common_series(ds_eo, var_eo, ds_wb, s_eo_override=s_csv)
        EO_LABEL = f"Agramon CSV {csv_col} ({cond})"
        print(f"  EO reference for Figs A-D: CSV column '{csv_col}' "
              f"(feature_id={csv_forcing[1]}, condition={cond}) — NetCDF kept "
              f"for spatial maps (E) and LAI-zone figures (F/G).")

    print("\n── Fig A : monthly mean by year ─────────────────────────────────")
    plot_compare_monthly_by_year(s_eo, s_wb, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig B : annual total bar chart ───────────────────────────────")
    plot_compare_annual_total(s_eo, s_wb, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig C : monthly total bar chart ──────────────────────────────")
    plot_compare_monthly_total(s_eo, s_wb, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig D : scatter EO-EB vs WB ──────────────────────────────────")
    plot_compare_scatter(s_eo, s_wb, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig E : spatial maps ─────────────────────────────────────────")
    plot_compare_spatial_maps(
        ds_eo, var_eo, ds_wb, selected_year,
        fig_path, label=label, dpi=dpi, show=show,
    )

<<<<<<< HEAD
    if df_csv is not None and csv_forcing is not None:
        print("\n── Fig H : WB ETa vs Agramon CSV (real ET, dry/wet) ─────────────")
        plot_compare_wb_vs_csv(
            ds_wb, s_eo_nc, df_csv, csv_forcing[0], csv_forcing[1],
            fig_path, label=label, dpi=dpi, show=show,
        )
    else:
        print("\n── Fig H : CSV comparison skipped (not an etp_ts run / no CSV) ──")

=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
    if da_lai_ref is None:
        print("\n── Fig F/G : per-zone comparison skipped (no LAI reference) ─────")
        return

    print("\n── Building per-zone series (pre-fire LAI terciles) ──────────────")
    zone_series, _bounds = build_zone_series(
        ds_eo, var_eo, ds_wb, da_lai_ref, n_zones=n_zones
    )
    if not zone_series:
        print("  No zones produced overlapping data — Fig F/G skipped.")
        return

    print_zone_bias_summary(zone_series)

    print("\n── Fig F : monthly mean by year, per LAI zone ────────────────────")
    plot_compare_monthly_by_year_per_zone(zone_series, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig G : annual bias per LAI zone ──────────────────────────────")
    plot_compare_annual_bias_per_zone(zone_series, fig_path, label=label, dpi=dpi, show=show)


# =============================================================================
# CLI argument parsing
# =============================================================================

<<<<<<< HEAD
def parse_args() -> tuple[list[str], int, int, str, int, bool, str, str, Path, int | None, bool]:
=======
def parse_args() -> tuple[list[str], int, int, str]:
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
    parser = argparse.ArgumentParser(
        description="Compare EO-EB vs Water-Balance ETa for Agramon"
    )
    parser.add_argument(
        "--plots", nargs="+", metavar="PLOT", default=["0"],
        help="Plots: row index, PlotID name, 'all', or mixed. Default: 0.",
    )
    parser.add_argument(
        "--year", type=int, default=2020,
        help="Year for spatial-map mosaics (default: 2020).",
    )
    parser.add_argument(
        "--scenario", type=int, default=6,
        help="sim_index to load from the simulation log (default: 0).",
    )
    parser.add_argument(
        "--wb-et-var", default=WB_ET_VAR,
        choices=["ACT. ETRA", "ACT. ETRA_patched"],
        help=(
            "Which WB ET variable to compare: 'ACT. ETRA' (raw solver "
            "output) or 'ACT. ETRA_patched' (artefact-floored). "
            f"Default: '{WB_ET_VAR}' (the WB_ET_VAR module constant)."
        ),
    )
    parser.add_argument(
        "--n-zones", type=int, default=N_LAI_ZONES,
        help=f"Number of LAI zones for Fig F/G (default: {N_LAI_ZONES}).",
    )
    parser.add_argument(
        "--no-zones", action="store_true",
        help="Skip the per-zone comparison (Fig F/G) even if LAI_PATH resolves.",
    )
<<<<<<< HEAD
    parser.add_argument(
        "--path2prj", default=str(WB_PRJ_ROOT),
        help=(
            "Root of the CATHY project outputs — mirrors withLAI_results.py's "
            "--path2prj. scenario_<N>[_tag] folders are resolved under "
            f"<path2prj>/outputs/, and the log CSV under "
            f"<path2prj>/../Codigos_HydroModel/. Default: {WB_PRJ_ROOT}"
        ),
    )
    parser.add_argument(
        "--log-file", default=LOG_FILE,
        help=(
            "Simulation log CSV — mirrors withLAI_results.py's --log-file. "
            f"Default: {LOG_FILE}"
        ),
    )

    parser.add_argument(
        "--et-csv", default=str(ET_CSV_PATH),
        help=f"Agramon treatments ET CSV for Fig H. Default: {ET_CSV_PATH}",
    )
    parser.add_argument(
        "--csv-fid", type=int, default=None,
        help=("feature_id in the CSV for Fig H. Default: parsed from the "
              "scenario folder name (..._forcing-etpts-<dry|wet>-fid<K>)."),
    )
    parser.add_argument(
        "--no-csv", action="store_true",
        help="Skip Fig H (WB ETa vs Agramon CSV).",
    )

    known_flags = {"--plots", "--year", "--scenario", "--wb-et-var", "--n-zones",
                    "--path2prj", "--log-file", "--et-csv", "--csv-fid"}
    bool_flags  = {"--no-zones", "--no-csv"}
=======

    known_flags = {"--plots", "--year", "--scenario", "--wb-et-var", "--n-zones"}
    bool_flags  = {"--no-zones"}
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
    clean: list[str] = []
    i = 1
    while i < len(sys.argv):
        tok = sys.argv[i]
        if tok in bool_flags:
            clean.append(tok)
            i += 1
        elif tok in known_flags:
            clean.append(tok)
            i += 1
            while i < len(sys.argv) and not sys.argv[i].startswith("--"):
                clean.append(sys.argv[i])
                i += 1
        else:
            i += 1

    args = parser.parse_args(clean)
<<<<<<< HEAD
    return (args.plots, args.year, args.scenario, args.wb_et_var, args.n_zones,
            args.no_zones, args.path2prj, args.log_file,
            Path(args.et_csv), args.csv_fid, args.no_csv)
=======
    return args.plots, args.year, args.scenario, args.wb_et_var, args.n_zones, args.no_zones
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
<<<<<<< HEAD
    (selection, selected_year, scenario, WB_ET_VAR, n_zones, no_zones,
     path2prj, log_file, et_csv, csv_fid, no_csv) = parse_args()

    # Mirror the --path2prj / --log-file that plot_scenarii_withETp_withLAI.sh
    # already passes to withLAI_results.py for this scenario, instead of the
    # module's hardcoded defaults — otherwise --scenario resolves against the
    # wrong project tree and every load below fails to find the ET data.
    WB_PRJ_ROOT = Path(path2prj)
    LOG_FILE    = log_file
=======
    selection, selected_year, scenario, WB_ET_VAR, n_zones, no_zones = parse_args()
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4

    print("\n── Loading shapefiles ──────────────────────────────────────────")
    gdf_proj  = load_shapefiles()

    print(f"\n── Loading WB ETa (scenario {scenario}) ────────────────────────")
    ds_wb, sim_dir = load_wb_et(scenario)

    # --plots defaults to ["0"] (shapefile row 0), which only matches the
    # right dem-plot by coincidence. Each scenario's WB grid is already
    # cropped to its own plot by the solver, so when the user hasn't
    # explicitly overridden --plots, derive the clip polygon from the WB
    # grid's own extent instead of blindly trusting row 0 — otherwise
    # scenarios whose plot != row 0 fail with rioxarray's NoDataInBounds.
    if selection == ["0"]:
        gdf_clip = resolve_plot_for_dataset(gdf_proj, ds_wb)
        label    = f"auto_{sim_dir.name}"
    else:
        gdf_clip = resolve_plot_selection(gdf_proj, selection)
        label    = selection_label(gdf_proj, selection)

    print("\n── Loading EO-EB ETa ───────────────────────────────────────────")
    ds_eo     = load_eo_et()
    ds_eo_c   = _clip_ds(ds_eo, gdf_clip)
    # Condition rule: WB forced with wet (dry) ETp -> compare with the wet
    # (dry) EO-EB variable when the file provides one (all EB-vs-WB figures).
    wb_cond   = forcing_condition(sim_dir)
    var_eo    = pick_eo_var(ds_eo_c, wb_cond)

    ds_wb_c   = _clip_ds(ds_wb, gdf_clip)

    da_lai_ref = None
    if not no_zones:
        print("\n── Loading LAI (for per-zone comparison) ────────────────────────")
        ds_lai = load_lai_ref()
        if ds_lai is not None:
            ds_lai_c  = _clip_ds(ds_lai, gdf_clip)
            var_lai   = next(iter(ds_lai_c.data_vars))
            da_lai_ref = compute_lai_zone_reference(ds_lai_c[var_lai])

    # Mirror withLAI_results.py: reuse the resolved (possibly tagged)
    # scenario_<N>[_tag] folder name as the figures subfolder, so ET
    # compare plots from different scenarios never overwrite each other
    # and stay easy to match back to their source outputs/ folder.
    scenario_dirname = sim_dir.name
    scenario_fig_path = fig_path / scenario_dirname

<<<<<<< HEAD
    # Fig H inputs: only for etp_ts runs (forcing column + feature_id are in
    # the scenario folder name); --csv-fid overrides the feature_id.
    df_csv, csv_forcing = None, None
    if not no_csv:
        parsed = parse_etpts_forcing(sim_dir)
        if parsed is None:
            print("\n  Not an etp_ts run (no 'forcing-etpts-<dry|wet>-fid<K>' in "
                  f"'{sim_dir.name}') — Fig H skipped.")
        else:
            col, fid = parsed
            fid = csv_fid if csv_fid is not None else fid
            print(f"\n── Loading Agramon ET CSV (feature_id={fid}, forcing={col}) ──")
            df_csv = load_csv_et(et_csv, fid)
            csv_forcing = (col, fid) if df_csv is not None else None

=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
    print(f"\n── Running comparison ──────────────────────────────────────────")
    print(f"  EO-EB var    : {var_eo}")
    print(f"  WB condition : {wb_cond or 'n/a (not an etp_ts run)'}")
    print(f"  WB var       : {WB_ET_VAR}")
    print(f"  Plot label   : {label}")
    print(f"  Spatial year : {selected_year}")
    print(f"  Scenario dir : {scenario_dirname}")
    print(f"  Output dir   : {scenario_fig_path}\n")

    run_comparison(
        ds_eo_c, var_eo,
        ds_wb_c,
        label=label,
        selected_year=selected_year,
        show=False,
        fig_path=scenario_fig_path,
        da_lai_ref=da_lai_ref,
        n_zones=n_zones,
<<<<<<< HEAD
        df_csv=df_csv,
        csv_forcing=csv_forcing,
=======
>>>>>>> 1d78c26b5b57db0ec851d536dacab6ca80d61ed4
    )

    print("\nDone.")
