"""
analyse_ET.py
=============
GRwater project — Agramon study site

Loads ET and ET0 NetCDF files and produces four ETa diagnostic figures:

    plot_eta_monthly_by_year   — spatial-mean ETa per month, one line per year
    plot_eta_spatial_maps_year — monthly spatial maps for a chosen year (mosaic)
    plot_eta_annual_hist       — histogram of per-pixel annual ETa sums
    plot_eta_annual_total      — bar chart of domain-total (spatial sum) annual ETa

Usage
-----
    python analyse_ET.py                         # default: plot index 0
    python analyse_ET.py --plots 0 2 4           # by index
    python analyse_ET.py --plots Plot4Mulching   # by name
    python analyse_ET.py --plots 1 Plot6Control  # mixed
    python analyse_ET.py --plots all             # all plots
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401 — registers .rio accessor
import xarray as xr

# =============================================================================
# Paths & constants
# =============================================================================

PathET     = Path("/home/z0272571a@campus.csic.es/Nextcloud/GRwater/data/satellite/Agramon")
TARGET_CRS = "EPSG:25830"
NATIVE_CRS = "EPSG:32630"

MODULE_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
    "/Tech4agro_org/GRwater_geophy"
).resolve()
if str(MODULE_PATH) not in sys.path:
    sys.path.append(str(MODULE_PATH))

import Agramon_utils as AgUtils
import results_plotter
from results_plotter import save_fig, maybe_show

fig_path = Path().cwd() / "../Analysis"
et_path  = PathET / "20161001_20241231_ET.nc"
etp_path = PathET / "20161001_20241231_ET0.nc"

# =============================================================================
# Shapefile loading
# =============================================================================

def load_shapefiles() -> gpd.GeoDataFrame:
    """Load Agramon plot shapefiles, reproject to TARGET_CRS."""
    gdf = AgUtils.load_plot_shapefiles(['microcuencas_13'], 
                                               MODULE_PATH / "shapefiles")
    if gdf.crs is None:
        raise ValueError("Shapefile has no CRS — assign one before proceeding.")
    if str(gdf.crs) != TARGET_CRS:
        print(f"  Reprojecting shapefiles {gdf.crs} → {TARGET_CRS}")
        gdf = gdf.to_crs(TARGET_CRS)
    print(f"  Shapefile CRS : {gdf.crs}")
    print("  Available plots (index : PlotID):")
    for i, pid in enumerate(gdf["PlotID"]):
        print(f"    [{i}] {pid}")
    return gdf


# =============================================================================
# Plot selection helpers
# =============================================================================

def resolve_plot_selection(
    gdf: gpd.GeoDataFrame,
    selection: list[str | int] | None,
) -> gpd.GeoDataFrame:
    """
    Filter gdf by PlotID, accepting indices, names, or a mix.

    Parameters
    ----------
    gdf       : GeoDataFrame with a 'PlotID' column
    selection : None / ["all"]     → all plots
                [0, 2]             → by row index
                ["Plot4Mulching"]  → by name
                [1, "Plot6Control"] → mixed
    """
    plot_ids = list(gdf["PlotID"])

    if selection is None or selection == ["all"]:
        print(f"  Selected plots : all ({len(plot_ids)} plots)")
        return gdf.copy()

    resolved: list[str] = []
    for item in selection:
        if isinstance(item, int) or (isinstance(item, str) and item.isdigit()):
            idx = int(item)
            if idx < 0 or idx >= len(plot_ids):
                raise IndexError(
                    f"Index {idx} out of range — valid range 0–{len(plot_ids) - 1}."
                )
            resolved.append(plot_ids[idx])
        else:
            if item not in plot_ids:
                raise ValueError(
                    f"PlotID '{item}' not found.\nAvailable: {plot_ids}"
                )
            resolved.append(item)

    # deduplicate, preserving order
    seen: set[str] = set()
    unique = [p for p in resolved if not (p in seen or seen.add(p))]

    filtered = gdf[gdf["PlotID"].isin(unique)].copy().reset_index(drop=True)
    print(f"  Selected plots : {unique}")
    return filtered

def selection_label(gdf: gpd.GeoDataFrame, selection: list[str | int] | None) -> str:
    """
    Build a filename-safe label from the plot selection using PlotID names.

    None / ["all"]    → "all"
    [0]               → "Plot2PostFirefascine"
    [0, 2, 4]         → "Plot2PostFirefascine_Plot4Mulching_Plot6Control"
    ["Plot4Mulching"] → "Plot4Mulching"
    """
    plot_ids = list(gdf["PlotID"])

    if selection is None or selection == ["all"]:
        return "all"

    names: list[str] = []
    for item in selection:
        if isinstance(item, int) or (isinstance(item, str) and item.isdigit()):
            names.append(plot_ids[int(item)])
        else:
            names.append(item)

    # deduplicate, preserving order
    seen: set[str] = set()
    unique = [n for n in names if not (n in seen or seen.add(n))]

    return f'watershedNb_{selection}_{"_".join(unique)}'

# =============================================================================
# Dataset loaders
# =============================================================================

def _finalise_ds(ds: xr.Dataset, label: str) -> xr.Dataset:
    """Write CRS + spatial dims after load or reproject."""
    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  {label} final CRS : {ds.rio.crs}")
    print(f"  {label} x range   : {float(ds.x.min()):.1f} – {float(ds.x.max()):.1f}")
    print(f"  {label} y range   : {float(ds.y.min()):.1f} – {float(ds.y.max()):.1f}")
    return ds


def _load_nc(path: Path, label: str) -> xr.Dataset:
    """Open NetCDF, assign native CRS, reproject to TARGET_CRS."""
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


def load_et(path: Path = et_path) -> xr.Dataset:
    """Load ET dataset and reproject to TARGET_CRS."""
    return _load_nc(path, "ET")


def load_etp(path: Path = etp_path) -> xr.Dataset:
    """Load ETp dataset and reproject to TARGET_CRS."""
    return _load_nc(path, "ETp")


# =============================================================================
# Clipping
# =============================================================================

def clip_dataset(
    ds: xr.Dataset,
    gdf: gpd.GeoDataFrame,
    selection: list[str | int] | None = None,
) -> xr.Dataset:
    """
    Clip dataset to the (optionally filtered) GeoDataFrame extent.

    Variables are clipped individually to avoid rioxarray CRS-strip bugs.
    """
    gdf_clip = resolve_plot_selection(gdf, selection)
    clipped: dict[str, xr.DataArray] = {}
    for var in ds.data_vars:
        da = (
            ds[var]
            .rio.write_crs(TARGET_CRS, inplace=False)
            .rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
        )
        clipped[var] = da.rio.clip(
            gdf_clip.geometry, gdf_clip.crs, drop=True, all_touched=True
        )
    return xr.Dataset(clipped, attrs=ds.attrs)


# =============================================================================
# Plotting — thin wrappers that inject the label into filenames
#
# Root cause of the blank-figure bug:
#   results_plotter functions call maybe_show(fig, show) internally, which
#   calls plt.close(fig) even when show=False.  plt.gcf() after that returns
#   a NEW blank figure.
#
# Fix: pass show=False so the figure is saved inside results_plotter with its
#   default hardcoded name, then immediately rename / re-save it with our
#   label-aware name using the still-open figure object retrieved via
#   plt.get_fignums() before the close happens.
#
# Cleanest approach: create the figure, let results_plotter draw into it (it
#   always calls plt.tight_layout() before save_fig), capture it before close.
#   Since results_plotter always creates its own fig internally we instead
#   patch around the close by saving plt.gcf() *before* calling the plotter,
#   then retrieving the last figure number after the call.
# =============================================================================

def _last_open_fig() -> plt.Figure | None:
    """Return the most-recently created still-open Figure, or None."""
    nums = plt.get_fignums()
    return plt.figure(nums[-1]) if nums else None

def run_plots(
    ds_crop: xr.Dataset,
    et_var: str,
    label: str,
    selected_year: int = 2020,
    product_name: str = "ETa",
    show: bool = False,
) -> None:
    """
    Run all diagnostic figures for either ETa or ETp.
    """

    out_path = fig_path / product_name
    out_path.mkdir(parents=True, exist_ok=True)

    # Monthly mean by year
    results_plotter.plot_eta_monthly_by_year(
        ds_crop,
        et_var=et_var,
        fig_path=out_path,
        fname=f"{product_name}_monthly_mean_by_year_{label}.png",
        show=show,
        product_name=product_name,
    )

    # Spatial maps
    _orig_save = results_plotter.save_fig

    def _labelled_save(fig, fp, name, dpi=150):
        return _orig_save(
            fig,
            fp,
            f"{product_name}_spatial_maps_{selected_year}_{label}.png",
            dpi,
        )

    results_plotter.save_fig = _labelled_save

    try:
        results_plotter.plot_eta_spatial_maps_year(
            ds_crop,
            et_var=et_var,
            selected_year=selected_year,
            fig_path=out_path,
            show=show,
            product_name=product_name,
        )
    finally:
        results_plotter.save_fig = _orig_save

    # Annual histogram
    results_plotter.plot_eta_annual_hist(
        ds_crop,
        et_var=et_var,
        fig_path=out_path,
        fname=f"{product_name}_annual_sum_histogram_{label}.png",
        show=show,
        product_name=product_name,
    )

    # Monthly + annual totals
    results_plotter.plot_eta_total(
        ds_crop,
        fig_path=out_path,
        label=label,
        product_name=product_name,
        dpi=150,
        show=show,
    )

# =============================================================================
# CLI argument parsing
# =============================================================================

def parse_args() -> list[str] | None:
    """
    Parse --plots from the command line.  Default is index 0.

    Examples
    --------
    python analyse_ET.py                         # → [0]  (default)
    python analyse_ET.py --plots all             # → all plots
    python analyse_ET.py --plots 0 2             # → indices
    python analyse_ET.py --plots Plot4Mulching   # → name
    python analyse_ET.py --plots 1 Plot6Control  # → mixed
    """
    parser = argparse.ArgumentParser(description="ETa diagnostics for Agramon")
    parser.add_argument(
        "--plots",
        nargs="+",
        metavar="PLOT",
        default=["0"],          # ← default: plot index 0
        help=(
            "Plots to process: row index (0-based), PlotID name, 'all', or a mix. "
            "Default: 0 (first plot)."
        ),
    )
    parser.add_argument(
        "--year",
        type=int,
        default=2020,
        help="Year for spatial-map mosaic (default: 2020).",
    )

    # sys.argv can contain injected garbage from conda launchers, Spyder, or
    # other wrappers (e.g. path fragments like "Install" from
    # .../envs/Install/bin/python3).  Keep only tokens belonging to our two
    # known flags and their values; discard everything else.
    known_flags = {"--plots", "--year"}
    clean: list[str] = []
    i = 1  # skip argv[0] (script name)
    while i < len(sys.argv):
        tok = sys.argv[i]
        if tok in known_flags:
            clean.append(tok)
            i += 1
            while i < len(sys.argv) and not sys.argv[i].startswith("--"):
                clean.append(sys.argv[i])
                i += 1
        else:
            i += 1  # silently discard unrecognised tokens

    args = parser.parse_args(clean)

    # "all" is a special sentinel — pass through as-is for resolve_plot_selection
    return args.plots, args.year


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    selection, selected_year = parse_args()

    print("\n── Loading shapefiles ──────────────────────────────────────────")
    gdf_proj = load_shapefiles()

    print("\n── Loading ET datasets ─────────────────────────────────────────")
    ds_et  = load_et()
    ds_etp = load_etp()

    print("\n── Clipping ────────────────────────────────────────────────────")
    ds_et_crop  = clip_dataset(ds_et,  gdf_proj, selection)
    ds_etp_crop = clip_dataset(ds_etp, gdf_proj, selection)

    label  = selection_label(gdf_proj, selection)
    et_var = next(iter(ds_et_crop.data_vars))

    print(f"\n── Plotting ────────────────────────────────────────────────────")
    print(f"  ET variable  : {et_var}")
    print(f"  Plot label   : {label}")
    print(f"  Spatial year : {selected_year}")
    print(f"  Output dir   : {fig_path}\n")

    run_plots(ds_et_crop, et_var=et_var, label=label,
              selected_year=selected_year, show=False)

    et_var = next(iter(ds_etp_crop.data_vars))

    run_plots(ds_etp_crop, et_var=et_var, label=label,
              selected_year=selected_year, product_name="ETp", show=False)

    print("\nDone.")
