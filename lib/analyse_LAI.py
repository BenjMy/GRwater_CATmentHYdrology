"""
analyse_LAI.py
==============
GRwater project — Agramon study site

Loads the LAI NetCDF produced by the TSEB/biophysical pipeline and produces
four LAI diagnostic figures:

    plot_lai_monthly_by_year   — spatial-mean LAI per month, one line per year
    plot_lai_spatial_maps_year — monthly spatial maps for a chosen year (mosaic)
    plot_lai_annual_hist       — histogram of per-pixel annual mean LAI
    plot_lai_annual_mean       — bar chart of domain-mean annual LAI

Usage
-----
    python analyse_LAI.py                         # default: plot index 0
    python analyse_LAI.py --plots 0 2 4           # by index
    python analyse_LAI.py --plots Plot4Mulching   # by name
    python analyse_LAI.py --plots 1 Plot6Control  # mixed
    python analyse_LAI.py --plots all             # all plots
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

# LAI NetCDF — location taken from LAI_TSEB.py (lai_monthly.nc input)
LAI_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Training_Supervision"
    "/Supervision/Xela_Carracedo_Practicas_2026_data/301a-biophysical/agramon/input"
).resolve()

TARGET_CRS = "EPSG:25830"
NATIVE_CRS = "EPSG:32630"

LAI_VAR   = "LAI"          # variable name inside lai_monthly.nc
LAI_UNITS = "m² m⁻²"
LAI_VMAX  = 2.0            # colour-scale ceiling (adjust if needed)
CMAP_LAI  = "YlGn"

MODULE_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
    "/Tech4agro_org/GRwater_geophy"
).resolve()
if str(MODULE_PATH) not in sys.path:
    sys.path.append(str(MODULE_PATH))

import Agramon_utils as AgUtils
import results_plotter
from results_plotter import save_fig, maybe_show

fig_path = Path().cwd() / "../Analysis/LAI"
lai_file = LAI_PATH / "lai_monthly.nc"

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
    selection : None / ["all"]      → all plots
                [0, 2]              → by row index
                ["Plot4Mulching"]   → by name
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
# Dataset loader
# =============================================================================

def _finalise_ds(ds: xr.Dataset, label: str) -> xr.Dataset:
    """Write CRS + spatial dims after load or reproject."""
    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  {label} final CRS : {ds.rio.crs}")
    print(f"  {label} x range   : {float(ds.x.min()):.1f} – {float(ds.x.max()):.1f}")
    print(f"  {label} y range   : {float(ds.y.min()):.1f} – {float(ds.y.max()):.1f}")
    return ds


def load_lai(path: Path = lai_file) -> xr.Dataset:
    """
    Open lai_monthly.nc, assign native CRS and reproject to TARGET_CRS.

    The file uses a 'time' dimension (monthly) and spatial dims (x, y).
    If the file already carries EPSG:25830 metadata it is used directly;
    otherwise we assume EPSG:32630 (same convention as the ET files).
    """
    ds = (
        xr.open_dataset(path, decode_times=True, mask_and_scale=True)
          .drop_vars("spatial_ref", errors="ignore")
    )

    # Detect or assign native CRS
    detected = ds.rio.crs
    if detected is None:
        print(f"  LAI: no CRS metadata found — assuming {NATIVE_CRS}")
        ds = ds.rio.write_crs(NATIVE_CRS, inplace=False)
    else:
        print(f"  LAI native CRS : {detected}")

    # Ensure spatial dims are registered
    # lai_monthly.nc may use 'x'/'y' or 'lon'/'lat' — handle both
    x_dim = "x" if "x" in ds.dims else ("lon" if "lon" in ds.dims else None)
    y_dim = "y" if "y" in ds.dims else ("lat" if "lat" in ds.dims else None)
    if x_dim is None or y_dim is None:
        raise ValueError(
            f"Cannot identify spatial dimensions in {path}. "
            f"Available dims: {list(ds.dims)}"
        )
    if x_dim != "x" or y_dim != "y":
        ds = ds.rename({x_dim: "x", y_dim: "y"})
        print(f"  Renamed dims: ({x_dim},{y_dim}) → (x,y)")

    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)

    if str(ds.rio.crs) != TARGET_CRS:
        print(f"  Reprojecting LAI  {ds.rio.crs} → {TARGET_CRS}")
        ds = ds.rio.reproject(TARGET_CRS)

    return _finalise_ds(ds, "LAI")


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
# Plotting functions
# =============================================================================

def plot_lai_monthly_by_year(
    ds: xr.Dataset,
    lai_var: str = LAI_VAR,
    fig_path: Path = fig_path,
    fname: str = "LAI_monthly_mean_by_year.png",
    show: bool = False,
) -> None:
    """
    Spatial-mean LAI per calendar month, one line per year.
    Mirrors results_plotter.plot_eta_monthly_by_year in style.
    """
    da = ds[lai_var]

    # Spatial mean → time series
    mean_ts = da.mean(dim=[d for d in da.dims if d not in ("time",)])
    ts = mean_ts.to_series()

    years = sorted(ts.index.year.unique())
    cmap  = plt.get_cmap("tab10", len(years))

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, yr in enumerate(years):
        yr_ts = ts[ts.index.year == yr]
        monthly = yr_ts.groupby(yr_ts.index.month).mean()
        ax.plot(monthly.index, monthly.values,
                marker="o", lw=2, color=cmap(i), label=str(yr))

    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    ax.set(title="Spatial-mean LAI — monthly average by year",
           xlabel="Month", ylabel=f"LAI ({LAI_UNITS})")
    ax.legend(title="Year", bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.grid(True, alpha=0.4)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    out = fig_path / fname
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


def plot_lai_spatial_maps_year(
    ds: xr.Dataset,
    selected_year: int = 2020,
    lai_var: str = LAI_VAR,
    fig_path: Path = fig_path,
    fname: str | None = None,
    show: bool = False,
) -> None:
    """
    Monthly spatial maps (mosaic) for *selected_year*.
    Mirrors results_plotter.plot_eta_spatial_maps_year in style.
    """
    da = ds[lai_var]
    da_year = da.sel(time=da["time"].dt.year == selected_year)

    if da_year.sizes["time"] == 0:
        print(f"  Skipped — no LAI data for {selected_year}.")
        return

    times = [pd.Timestamp(t) for t in da_year["time"].values]
    ncols  = 4
    nrows  = (len(times) + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5 * ncols, 4 * nrows),
        squeeze=False,
    )

    for idx, t in enumerate(times):
        r, c = divmod(idx, ncols)
        ax   = axes[r][c]
        frame = da_year.sel(time=t).values
        im    = ax.imshow(frame, origin="upper", cmap=CMAP_LAI,
                          vmin=0.0, vmax=LAI_VMAX, aspect="auto")
        plt.colorbar(im, ax=ax, label=LAI_UNITS, fraction=0.046, pad=0.04)
        ax.set_title(t.strftime("%B %Y"), fontsize=10)
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")

    # hide unused axes
    for idx in range(len(times), nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    fig.suptitle(f"LAI monthly maps — {selected_year}", fontsize=13, fontweight="bold")
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    if fname is None:
        fname = f"LAI_spatial_maps_{selected_year}.png"
    out = fig_path / fname
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


def plot_lai_annual_hist(
    ds: xr.Dataset,
    lai_var: str = LAI_VAR,
    fig_path: Path = fig_path,
    fname: str = "LAI_annual_mean_histogram.png",
    show: bool = False,
) -> None:
    """
    Histogram of per-pixel annual mean LAI (all years combined).
    Mirrors results_plotter.plot_eta_annual_hist in style.
    """
    da = ds[lai_var]

    # Annual mean per pixel
    annual_mean = da.groupby("time.year").mean("time")  # (year, y, x)
    pixel_vals  = annual_mean.values.ravel()
    pixel_vals  = pixel_vals[np.isfinite(pixel_vals)]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(pixel_vals, bins=40, color="#78c679", edgecolor="white", alpha=0.85)
    ax.axvline(np.nanmedian(pixel_vals), color="#006837", lw=2,
               linestyle="--", label=f"Median = {np.nanmedian(pixel_vals):.2f}")
    ax.set(title="Distribution of per-pixel annual mean LAI",
           xlabel=f"LAI ({LAI_UNITS})", ylabel="Pixel count")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.4)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    out = fig_path / fname
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


def plot_lai_annual_mean(
    ds: xr.Dataset,
    lai_var: str = LAI_VAR,
    fig_path: Path = fig_path,
    label: str = "",
    fname: str | None = None,
    dpi: int = 150,
    show: bool = False,
) -> None:
    """
    Bar chart of domain-mean annual LAI with an overlaid trend line.
    Mirrors results_plotter.plot_eta_total in style.
    """
    da = ds[lai_var]

    # Spatial mean → annual mean
    spatial_mean = da.mean(dim=[d for d in da.dims if d not in ("time",)])
    annual       = spatial_mean.groupby("time.year").mean("time").to_series()

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.bar(annual.index, annual.values,
                  color="#78c679", edgecolor="white", width=0.6)
    ax.plot(annual.index, annual.values,
            marker="o", color="#006837", lw=2, zorder=5)

    # Annotate bar tops
    for bar, val in zip(bars, annual.values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.02 * LAI_VMAX,
                f"{val:.2f}", ha="center", va="bottom", fontsize=9)

    ax.set(title=f"Domain-mean annual LAI{' — ' + label if label else ''}",
           xlabel="Year", ylabel=f"LAI ({LAI_UNITS})")
    ax.set_xticks(annual.index)
    ax.grid(True, axis="y", alpha=0.4)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    if fname is None:
        suffix = f"_{label}" if label else ""
        fname  = f"LAI_annual_mean{suffix}.png"
    out = fig_path / fname
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Orchestrator
# =============================================================================

def run_plots(
    ds_lai_crop: xr.Dataset,
    lai_var: str,
    label: str,
    selected_year: int = 2020,
    show: bool = False,
) -> None:
    """Run all four LAI figures, saving each with the plot-selection label."""
    fig_path.mkdir(parents=True, exist_ok=True)

    # ── Fig A : monthly mean by year ─────────────────────────────────────────
    plot_lai_monthly_by_year(
        ds_lai_crop, lai_var=lai_var, fig_path=fig_path,
        fname=f"LAI_monthly_mean_by_year_{label}.png", show=show,
    )

    # ── Fig B : spatial maps for selected_year ───────────────────────────────
    plot_lai_spatial_maps_year(
        ds_lai_crop, selected_year=selected_year, lai_var=lai_var,
        fig_path=fig_path,
        fname=f"LAI_spatial_maps_{selected_year}_{label}.png", show=show,
    )

    # ── Fig C : annual histogram (per-pixel) ─────────────────────────────────
    plot_lai_annual_hist(
        ds_lai_crop, lai_var=lai_var, fig_path=fig_path,
        fname=f"LAI_annual_mean_histogram_{label}.png", show=show,
    )

    # ── Fig D : domain-mean annual bar chart ─────────────────────────────────
    plot_lai_annual_mean(
        ds_lai_crop, lai_var=lai_var, fig_path=fig_path,
        label=label, dpi=150, show=show,
    )


# =============================================================================
# CLI argument parsing
# =============================================================================

def parse_args() -> tuple[list[str], int]:
    """
    Parse --plots and --year from the command line.  Default is index 0.

    Examples
    --------
    python analyse_LAI.py                         # → [0]  (default)
    python analyse_LAI.py --plots all             # → all plots
    python analyse_LAI.py --plots 0 2             # → indices
    python analyse_LAI.py --plots Plot4Mulching   # → name
    python analyse_LAI.py --plots 1 Plot6Control  # → mixed
    """
    parser = argparse.ArgumentParser(description="LAI diagnostics for Agramon")
    parser.add_argument(
        "--plots",
        nargs="+",
        metavar="PLOT",
        default=["0"],
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

    # Strip injected garbage tokens (conda launchers, Spyder, etc.)
    known_flags = {"--plots", "--year"}
    clean: list[str] = []
    i = 1
    while i < len(sys.argv):
        tok = sys.argv[i]
        if tok in known_flags:
            clean.append(tok)
            i += 1
            while i < len(sys.argv) and not sys.argv[i].startswith("--"):
                clean.append(sys.argv[i])
                i += 1
        else:
            i += 1

    args = parser.parse_args(clean)
    return args.plots, args.year


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    selection, selected_year = parse_args()

    print("\n── Loading shapefiles ──────────────────────────────────────────")
    gdf_proj = load_shapefiles()

    print("\n── Loading LAI dataset ─────────────────────────────────────────")
    ds_lai = load_lai()

    if LAI_VAR not in ds_lai:
        raise KeyError(
            f"Variable '{LAI_VAR}' not found in {lai_file}.\n"
            f"Available variables: {list(ds_lai.data_vars)}"
        )

    print("\n── Clipping ────────────────────────────────────────────────────")
    ds_lai_crop = clip_dataset(ds_lai, gdf_proj, selection)

    label   = selection_label(gdf_proj, selection)
    lai_var = LAI_VAR

    print(f"\n── Plotting ────────────────────────────────────────────────────")
    print(f"  LAI variable : {lai_var}")
    print(f"  Plot label   : {label}")
    print(f"  Spatial year : {selected_year}")
    print(f"  Output dir   : {fig_path}\n")

    run_plots(ds_lai_crop, lai_var=lai_var, label=label,
              selected_year=selected_year, show=False)

    print("\nDone.")
