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
from matplotlib.colors import TwoSlopeNorm
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

# -----------------------------------------------------------------------------
# CLI defaults — edit here so the script runs with no arguments.
# (diagnose_ET.py shares parse_args(), so these apply to it as well.)
# -----------------------------------------------------------------------------
DEFAULT_PLOTS = ["all"]       # row indices, PlotID names, or "all"
DEFAULT_YEAR = 2019           # spatial-map year; also the default month-by-month year
DEFAULT_RATIO_YEARS = None    # None = every year present in the clipped data
DEFAULT_RATIO_MONTH_YEARS = None  # None = [DEFAULT_YEAR] (or --year if given)

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
    n = len(plot_ids)

    if selection is None or selection == ["all"]:
        print(f"  Selected plots : all ({n} plots)")
        return gdf.copy()

    # Resolve every selection item to a ROW POSITION rather than a PlotID
    # name. PlotID is not guaranteed unique (e.g. a shapefile where every
    # feature shares the same layer name) — filtering by name in that case
    # would silently match every row instead of just the one requested.
    positions: list[int] = []
    for item in selection:
        if isinstance(item, int) or (isinstance(item, str) and item.isdigit()):
            idx = int(item)
            if idx < 0 or idx >= n:
                raise IndexError(
                    f"Index {idx} out of range — valid range 0–{n - 1}."
                )
            positions.append(idx)
        else:
            matches = [i for i, pid in enumerate(plot_ids) if pid == item]
            if not matches:
                raise ValueError(
                    f"PlotID '{item}' not found.\nAvailable: {plot_ids}"
                )
            if len(matches) > 1:
                print(
                    f"  WARNING: PlotID '{item}' matches {len(matches)} rows "
                    f"(indices {matches}) — PlotID is not unique in this "
                    f"shapefile, selecting all of them. Use a row index "
                    f"instead (e.g. --plots {matches[0]}) to pick just one."
                )
            positions.extend(matches)

    # deduplicate, preserving order
    seen: set[int] = set()
    unique_positions = [p for p in positions if not (p in seen or seen.add(p))]

    filtered = gdf.iloc[unique_positions].copy().reset_index(drop=True)
    print(
        f"  Selected plots : row index {unique_positions} "
        f"→ PlotID {[plot_ids[i] for i in unique_positions]}"
    )
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

def plot_eta_etp_timeseries(
    ds_et_crop: xr.Dataset,
    ds_etp_crop: xr.Dataset,
    et_var: str,
    etp_var: str,
    fig_path: Path,
    label: str,
    show: bool = False,
) -> None:
    """
    Time-series diagnostic comparing ETa and ETp for the current plot selection.

    Top panel    : spatial-mean ETa and ETp per time step [mm]
    Bottom panel : ETa / ETp ratio per time step (dashed line at 1.0 for reference)

    NOTE: this is a standalone plot (not routed through results_plotter.py,
    which wasn't available when this was added) — it reuses the same
    save_fig / maybe_show helpers so behaviour and output layout stay
    consistent with the other diagnostics. If results_plotter.py already has
    a house style (colors, fonts, dpi) this doesn't pick it up automatically —
    happy to align it once that file's shared.
    """
    out_path = fig_path / "ETa_ETp_ratio"
    out_path.mkdir(parents=True, exist_ok=True)

    eta_ts = ds_et_crop[et_var].mean(dim=("x", "y"), skipna=True)
    etp_ts = ds_etp_crop[etp_var].mean(dim=("x", "y"), skipna=True)

    # Align on common time steps in case ET / ETp have slightly different coverage
    eta_ts, etp_ts = xr.align(eta_ts, etp_ts, join="inner")

    if eta_ts.sizes.get("time", 0) == 0:
        print("  [plot_eta_etp_timeseries] No overlapping time steps between ET and ETp — skipping.")
        return

    ratio_ts = eta_ts / etp_ts.where(etp_ts != 0)
    time = eta_ts["time"].values

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 6), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
    )

    ax1.plot(time, eta_ts.values, label="ETa", color="tab:blue", lw=1.2)
    ax1.plot(time, etp_ts.values, label="ETp", color="tab:orange", lw=1.2)
    ax1.set_ylabel("ET [mm]")
    ax1.legend(loc="upper right")
    ax1.set_title(f"ETa & ETp time series — {label}")
    ax1.grid(alpha=0.3)

    ax2.plot(time, ratio_ts.values, color="tab:green", lw=1.0)
    ax2.axhline(1.0, color="grey", lw=0.8, ls="--")
    ax2.set_ylabel("ETa / ETp")
    ax2.set_xlabel("Date")
    ax2.grid(alpha=0.3)

    fig.tight_layout()

    fname = f"ETa_ETp_ratio_timeseries_{label}.png"
    save_fig(fig, out_path, fname, dpi=150)
    maybe_show(fig, show)


def plot_ratio_spatial_by_year(
    ds_et_crop: xr.Dataset,
    ds_etp_crop: xr.Dataset,
    et_var: str,
    etp_var: str,
    years: list[int],
    fig_path: Path,
    label: str,
    gdf: "gpd.GeoDataFrame | None" = None,
    selection: list[str | int] | None = None,
    show: bool = False,
) -> None:
    """
    Spatial map of the annual ETa/ETp ratio (per-pixel annual ETa sum divided
    by per-pixel annual ETp sum), one panel per requested year, mosaicked
    side by side with a shared diverging colorbar centered at 1.0.

    A value > 1 means that pixel's ETa exceeded its ETp for that year;
    < 1 means the opposite.

    ds_et_crop / ds_etp_crop are expected to already be clipped (see
    clip_dataset). Because clip_dataset crops to the selected polygon's
    BOUNDING BOX (drop=True), an irregularly-shaped subcatchment can still
    look "unclipped" in the rendered rectangle even though pixels outside
    the polygon are correctly NaN. To make the clip visually unambiguous —
    and to guard against the two datasets not sharing an identical grid —
    this function:
      1. Re-clips both DataArrays to `gdf`/`selection` again if provided
         (redundant with clip_dataset, but catches any upstream leak).
      2. Draws the selected polygon's boundary on top of each panel.
      3. Prints the realised x/y extent so it can be checked against the
         expected subcatchment bounds.
    """
    out_path = fig_path / "ETa_ETp_ratio"
    out_path.mkdir(parents=True, exist_ok=True)

    eta_da = ds_et_crop[et_var]
    etp_da = ds_etp_crop[etp_var]

    gdf_clip = None
    if gdf is not None:
        gdf_clip = resolve_plot_selection(gdf, selection)
        eta_da = eta_da.rio.clip(gdf_clip.geometry, gdf_clip.crs, drop=True, all_touched=True)
        etp_da = etp_da.rio.clip(gdf_clip.geometry, gdf_clip.crs, drop=True, all_touched=True)

    eta_da, etp_da = xr.align(eta_da, etp_da, join="inner")

    print(
        f"  [plot_ratio_spatial_by_year] realised extent — "
        f"x: {float(eta_da.x.min()):.1f}–{float(eta_da.x.max()):.1f}, "
        f"y: {float(eta_da.y.min()):.1f}–{float(eta_da.y.max()):.1f} "
        f"({eta_da.sizes.get('x', 0)}x{eta_da.sizes.get('y', 0)} px)"
    )

    ratio_maps: dict[int, xr.DataArray] = {}
    for yr in years:
        eta_yr = eta_da.sel(time=eta_da["time"].dt.year == yr).sum(dim="time", skipna=True)
        etp_yr = etp_da.sel(time=etp_da["time"].dt.year == yr).sum(dim="time", skipna=True)
        if eta_yr.size == 0 or etp_yr.size == 0:
            print(f"  [plot_ratio_spatial_by_year] No data for year {yr} — skipping.")
            continue
        ratio_maps[yr] = eta_yr / etp_yr.where(etp_yr != 0)

    if not ratio_maps:
        print("  [plot_ratio_spatial_by_year] No years had data — nothing to plot.")
        return

    all_vals = np.concatenate(
        [r.values.ravel()[~np.isnan(r.values.ravel())] for r in ratio_maps.values()]
    )
    if all_vals.size == 0:
        print("  [plot_ratio_spatial_by_year] All-NaN ratio maps — nothing to plot.")
        return

    # Diverging colormap centred on 1.0; bounds from data (2nd/98th pct),
    # but always straddling 1.0 so the reference colour stays meaningful.
    vmin = min(float(np.nanpercentile(all_vals, 2)), 0.5)
    vmax = max(float(np.nanpercentile(all_vals, 98)), 1.5)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    yrs_ok = list(ratio_maps.keys())
    ncols = min(len(yrs_ok), 4)
    nrows = int(np.ceil(len(yrs_ok) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False
    )

    im = None
    for ax, yr in zip(axes.ravel(), yrs_ok):
        im = ratio_maps[yr].plot(
            ax=ax, cmap="RdBu_r", norm=norm, add_colorbar=False
        )
        if gdf_clip is not None:
            gdf_clip.boundary.plot(ax=ax, color="black", linewidth=1.2)
        ax.set_title(str(yr))
        ax.set_xlabel("")
        ax.set_ylabel("")

    for ax in axes.ravel()[len(yrs_ok):]:
        ax.axis("off")

    fig.suptitle(f"Annual ETa / ETp ratio — {label}")
    fig.subplots_adjust(right=0.88, top=0.85)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.65])
    fig.colorbar(im, cax=cbar_ax, label="ETa / ETp")

    fname = f"ETa_ETp_ratio_spatial_{'-'.join(str(y) for y in yrs_ok)}_{label}.png"
    save_fig(fig, out_path, fname, dpi=150)
    maybe_show(fig, show)


_MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def plot_ratio_spatial_by_month(
    ds_et_crop: xr.Dataset,
    ds_etp_crop: xr.Dataset,
    et_var: str,
    etp_var: str,
    year: int,
    fig_path: Path,
    label: str,
    gdf: "gpd.GeoDataFrame | None" = None,
    selection: list[str | int] | None = None,
    show: bool = False,
) -> None:
    """
    Spatial map of the ETa/ETp ratio for each month of a single year
    (per-pixel monthly ETa sum divided by per-pixel monthly ETp sum), one
    panel per month present, mosaicked with a shared diverging colorbar
    centered at 1.0.

    Same re-clip / boundary-overlay / extent-print safeguards as
    plot_ratio_spatial_by_year — see that function's docstring for why.
    """
    out_path = fig_path / "ETa_ETp_ratio"
    out_path.mkdir(parents=True, exist_ok=True)

    eta_da = ds_et_crop[et_var]
    etp_da = ds_etp_crop[etp_var]

    gdf_clip = None
    if gdf is not None:
        gdf_clip = resolve_plot_selection(gdf, selection)
        eta_da = eta_da.rio.clip(gdf_clip.geometry, gdf_clip.crs, drop=True, all_touched=True)
        etp_da = etp_da.rio.clip(gdf_clip.geometry, gdf_clip.crs, drop=True, all_touched=True)

    eta_da, etp_da = xr.align(eta_da, etp_da, join="inner")

    print(
        f"  [plot_ratio_spatial_by_month] realised extent — "
        f"x: {float(eta_da.x.min()):.1f}–{float(eta_da.x.max()):.1f}, "
        f"y: {float(eta_da.y.min()):.1f}–{float(eta_da.y.max()):.1f} "
        f"({eta_da.sizes.get('x', 0)}x{eta_da.sizes.get('y', 0)} px)"
    )

    eta_yr = eta_da.sel(time=eta_da["time"].dt.year == year)
    etp_yr = etp_da.sel(time=etp_da["time"].dt.year == year)

    if eta_yr.sizes.get("time", 0) == 0:
        print(f"  [plot_ratio_spatial_by_month] No data for year {year} — skipping.")
        return

    months_present = sorted(set(pd.to_datetime(eta_yr["time"].values).month))

    ratio_maps: dict[int, xr.DataArray] = {}
    for m in months_present:
        eta_m = eta_yr.sel(time=eta_yr["time"].dt.month == m).sum(dim="time", skipna=True)
        etp_m = etp_yr.sel(time=etp_yr["time"].dt.month == m).sum(dim="time", skipna=True)
        ratio_maps[m] = eta_m / etp_m.where(etp_m != 0)

    all_vals = np.concatenate(
        [r.values.ravel()[~np.isnan(r.values.ravel())] for r in ratio_maps.values()]
    )
    if all_vals.size == 0:
        print(f"  [plot_ratio_spatial_by_month] All-NaN ratio maps for {year} — nothing to plot.")
        return

    vmin = min(float(np.nanpercentile(all_vals, 2)), 0.5)
    vmax = max(float(np.nanpercentile(all_vals, 98)), 1.5)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    ncols = min(len(months_present), 4)
    nrows = int(np.ceil(len(months_present) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False
    )

    im = None
    for ax, m in zip(axes.ravel(), months_present):
        im = ratio_maps[m].plot(
            ax=ax, cmap="RdBu_r", norm=norm, add_colorbar=False
        )
        if gdf_clip is not None:
            gdf_clip.boundary.plot(ax=ax, color="black", linewidth=1.2)
        ax.set_title(_MONTH_NAMES[m - 1])
        ax.set_xlabel("")
        ax.set_ylabel("")

    for ax in axes.ravel()[len(months_present):]:
        ax.axis("off")

    fig.suptitle(f"Monthly ETa / ETp ratio — {year} — {label}")
    fig.subplots_adjust(right=0.88, top=0.9)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.65])
    fig.colorbar(im, cax=cbar_ax, label="ETa / ETp")

    fname = f"ETa_ETp_ratio_spatial_monthly_{year}_{label}.png"
    save_fig(fig, out_path, fname, dpi=150)
    maybe_show(fig, show)


# =============================================================================
# CLI argument parsing
# =============================================================================

def parse_args() -> tuple[list[str], int, list[int] | None, list[int] | None]:
    """
    Parse --plots from the command line.  Default is index 0.

    Examples
    --------
    python analyse_ET.py                         # → DEFAULT_PLOTS / DEFAULT_YEAR
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
        default=DEFAULT_PLOTS,  # ← edit DEFAULT_PLOTS at the top of the file
        help=(
            "Plots to process: row index (0-based), PlotID name, 'all', or a mix. "
            f"Default: {DEFAULT_PLOTS}."
        ),
    )
    parser.add_argument(
        "--year",
        type=int,
        default=DEFAULT_YEAR,
        help=f"Year for spatial-map mosaic (default: {DEFAULT_YEAR}).",
    )
    parser.add_argument(
        "--ratio-years",
        nargs="+",
        type=int,
        default=DEFAULT_RATIO_YEARS,
        metavar="YEAR",
        help=(
            "One or more years to map for the spatial ETa/ETp ratio mosaic "
            "(e.g. --ratio-years 2019 2020 2023). Default: every year present "
            "in the clipped data."
        ),
    )
    parser.add_argument(
        "--ratio-month-years",
        nargs="+",
        type=int,
        default=DEFAULT_RATIO_MONTH_YEARS,
        metavar="YEAR",
        help=(
            "One or more years to break out month-by-month in a spatial "
            "ETa/ETp ratio mosaic (e.g. --ratio-month-years 2022 2023). "
            "Each year gets its own 12-panel figure. Default: the year given "
            "by --year."
        ),
    )

    # sys.argv can contain injected garbage from conda launchers, Spyder, or
    # other wrappers (e.g. path fragments like "Install" from
    # .../envs/Install/bin/python3).  Keep only tokens belonging to our
    # known flags and their values; discard everything else.
    known_flags = {"--plots", "--year", "--ratio-years", "--ratio-month-years"}
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
    return args.plots, args.year, args.ratio_years, args.ratio_month_years


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    selection, selected_year, ratio_years, ratio_month_years = parse_args()

    print("\n── Loading shapefiles ──────────────────────────────────────────")
    gdf_proj = load_shapefiles()

    print("\n── Loading ET datasets ─────────────────────────────────────────")
    ds_et  = load_et()
    ds_etp = load_etp()

    print("\n── Clipping ────────────────────────────────────────────────────")
    ds_et_crop  = clip_dataset(ds_et,  gdf_proj, selection)
    ds_etp_crop = clip_dataset(ds_etp, gdf_proj, selection)

    label   = selection_label(gdf_proj, selection)
    et_var  = next(iter(ds_et_crop.data_vars))
    etp_var = next(iter(ds_etp_crop.data_vars))

    if ratio_years is None:
        ratio_years = sorted(set(pd.to_datetime(ds_et_crop["time"].values).year))
    if ratio_month_years is None:
        ratio_month_years = [selected_year]

    print(f"\n── Plotting ────────────────────────────────────────────────────")
    print(f"  ET variable      : {et_var}")
    print(f"  ETp variable     : {etp_var}")
    print(f"  Plot label       : {label}")
    print(f"  Spatial year     : {selected_year}")
    print(f"  Ratio years      : {ratio_years}")
    print(f"  Ratio month-years: {ratio_month_years}")
    print(f"  Output dir       : {fig_path}\n")

    run_plots(ds_et_crop, et_var=et_var, label=label,
              selected_year=selected_year, show=False)

    run_plots(ds_etp_crop, et_var=etp_var, label=label,
              selected_year=selected_year, product_name="ETp", show=False)

    plot_eta_etp_timeseries(
        ds_et_crop, ds_etp_crop,
        et_var=et_var, etp_var=etp_var,
        fig_path=fig_path, label=label,
        show=False,
    )

    plot_ratio_spatial_by_year(
        ds_et_crop, ds_etp_crop,
        et_var=et_var, etp_var=etp_var,
        years=ratio_years,
        fig_path=fig_path, label=label,
        gdf=gdf_proj, selection=selection,
        show=False,
    )

    for yr in ratio_month_years:
        plot_ratio_spatial_by_month(
            ds_et_crop, ds_etp_crop,
            et_var=et_var, etp_var=etp_var,
            year=yr,
            fig_path=fig_path, label=label,
            gdf=gdf_proj, selection=selection,
            show=False,
        )

    print("\nDone.")
