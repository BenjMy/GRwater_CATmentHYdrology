"""
results_plotter.py — Shared Plotting Module
============================================
GRwater project — Agramon study site

Centralises every figure-generation function used by both:
  - test_withLAI_results.py   (LAI loop post-processing)
  - LT_SShydro_results.py     (long-term subsurface hydrology post-processing)

All functions are pure: they receive data, produce a matplotlib Figure,
save it to *fig_path*, and optionally display it.  No global state is
mutated.  Callers are responsible for loading artefacts and passing them in.

Public API
----------
Shared helpers
    save_fig(fig, fig_path, name, dpi)
    maybe_show(fig, show)
    time_label(t_val)
    et_global_limits(et_dataset, var)
    veg_colorbar(ax, cmap, labels)
    date_window(date_start, date_end, df_dates)
    add_fire_vline(ax, d0, d1, label, color)   # marks FIRE_DATE = 2020-07-01
    plot_catchment_et_freq(ET_xr, ET_xr_mean, time_dates, freq, fig_path, dpi, show)

LAI-loop figures (test_withLAI_results)
    plot_recharge_bar(recharge_monthly, fig_path, dpi, show)
    plot_et_monthly_mean(ET_xr_all, et_var, et_scale, fig_path, dpi, show)
    plot_et_by_veg_type(ET_xr_all, veg_map_history, et_var, et_scale, fig_path, dpi, show)
    plot_veg_map_samples(veg_map_history, cmap_veg, fig_path, dpi, show)
    plot_veg_map_year(veg_map_history, selected_year, cmap_veg, fig_path, dpi, show)
    plot_lai_map_year(ds_lai, var_name, selected_year, cmap, vmin, vmax, fig_path, dpi, show)
    make_veg_gif(veg_history, cmap_veg, out_path, duration)
    make_et_gif(et_dataset, var, scale, out_path, cmap, duration)

LAI diagnostic figures (analyse_LAI / withLAI_results)
    plot_lai_monthly_by_year(ds, lai_var, fig_path, dpi, show, fname)
    plot_lai_spatial_maps_year(ds, lai_var, selected_year, fig_path, cmap, vmin, vmax, dpi, show, fname)
    plot_lai_annual_hist(ds, lai_var, fig_path, bins, dpi, show, fname)
    plot_lai_annual_mean(ds, lai_var, fig_path, label, dpi, show, fname)
    plot_lai_annual_median(ds, lai_var, fig_path, label, dpi, show, fname)

SSHydro figures (LT_SShydro_results)
    plot_hydro(data, nodes_dict, args, fig_path)
    plot_et_sw(data, nodes_dict, args, fig_path, colors)
    plot_era5(data, args, fig_path)
    plot_catchment_et(data, fig_path, dpi, show)
    plot_hgraph(data, fig_path, dpi, show, hgraph_fname)
    plot_node_locations(grid3d, nodes_dict, fig_path, dem_da, colors_dict, scenario_id, dpi, show)
    plot_recharge(data, fig_path, dpi, show, scenario_id)
    plot_recharge_xr(xr_recharge, fig_path, dpi, show, scenario_id)

LAI-loop xarray-based hydro figures (test_withLAI_results — Fig 6 / 7)
    plot_hydro_xr(psi_xr, sw_xr, pev_series, tp_series, nodes_dict,
                  fig_path, date_start, date_end, ...)
    plot_et_sw_xr(sw_xr, ET_xr_all, et_var, et_scale, nodes_dict,
                  fig_path, date_start, date_end, ...)
"""

from __future__ import annotations

import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import imageio
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from pyCATHY import meshtools as mt  # local import — optional dependency
from mpl_toolkits.axes_grid1 import make_axes_locatable


# ── Colour constants (callers may override via kwargs) ────────────────────────

CMAP_ET  = "YlGnBu"
CMAP_VEG = mcolors.ListedColormap(["#d4b483", "#78c679", "#006837"])

COLORS_HYDRO: dict[str, str] = {
    "uphill_surface":   "#0072B2",
    "uphill_1m":        "#56B4E9",
    "downhill_surface": "#D55E00",
    "downhill_1m":      "#E69F00",
    "ET_mean":          "#009E73",
}

# ── Fire event marker ──────────────────────────────────────────────────────

FIRE_DATE = pd.Timestamp("2020-07-01")


def add_fire_vline(
    ax: plt.Axes,
    d0: pd.Timestamp | None = None,
    d1: pd.Timestamp | None = None,
    *,
    label: str = "Fire (Jul 2020)",
    color: str = "#B22222",
) -> None:
    """Draw a dashed vertical line marking the Agramon fire date (2020-07).

    Intended for any axis whose x-values are real datetimes. If *d0*/*d1*
    (the plotted date window) are supplied and ``FIRE_DATE`` falls outside
    that window, nothing is drawn — this keeps scenarios that don't cover
    July 2020 from having their x-axis distorted by an out-of-range line.
    """
    if d0 is not None and FIRE_DATE < pd.Timestamp(d0):
        return
    if d1 is not None and FIRE_DATE > pd.Timestamp(d1):
        return
    ax.axvline(FIRE_DATE, color=color, linestyle="--", linewidth=1.3,
               alpha=0.85, zorder=4, label=label)


# =============================================================================
# Shared helpers
# =============================================================================

def save_fig(fig: plt.Figure, fig_path: Path, name: str, dpi: int = 150) -> Path:
    """Save *fig* to *fig_path / name* and print the destination."""
    out = fig_path / name
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    return out


def maybe_show(fig: plt.Figure, show: bool) -> None:
    """Optionally display *fig*, then close it."""
    if show:
        plt.show()
    plt.close(fig)


def time_label(t_val: Any) -> str:
    """Human-readable string for a time coordinate (timedelta64 or datetime64)."""
    if np.issubdtype(np.array(t_val).dtype, np.timedelta64):
        t_days = float(t_val / np.timedelta64(1, "s")) / 86400
        return f"t = {t_days:.1f} d"
    return pd.Timestamp(t_val).strftime("%Y-%m-%d")


def et_global_limits(et_dataset: xr.Dataset, var: str) -> tuple[float, float]:
    """Return (vmin, vmax) across the whole ET dataset for consistent colour axes."""
    vals = et_dataset[var].values
    return float(np.nanmin(vals)), float(np.nanmax(vals))


def veg_colorbar(
    ax: plt.Axes,
    cmap: mcolors.Colormap,
    labels: tuple[str, ...] = ("Bare", "Moderate veg.", "Dense veg."),
) -> plt.colorbar:
    """Attach a discrete 3-class colourbar to *ax*."""
    bounds = [0.5, 1.5, 2.5, 3.5]
    norm   = mcolors.BoundaryNorm(bounds, cmap.N)
    sm     = cm.ScalarMappable(cmap=cmap, norm=norm)
    cb     = plt.colorbar(sm, ax=ax, ticks=[1, 2, 3])
    cb.ax.set_yticklabels(labels)
    return cb


def date_window(
    date_start: str | None,
    date_end:   str | None,
    df_dates:   pd.Series | pd.DatetimeIndex,
) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Return (d0, d1) as Timestamps, falling back to series bounds."""
    d0 = pd.Timestamp(date_start) if date_start else pd.Timestamp(df_dates.min())
    d1 = pd.Timestamp(date_end)   if date_end   else pd.Timestamp(df_dates.max())
    return d0, d1


# =============================================================================
# LAI-loop figures   (used by test_withLAI_results.py)
# =============================================================================

def plot_recharge_bar(
    recharge_monthly: pd.Series,
    fig_path: Path,
    *,
    dpi:  int  = 150,
    show: bool = False,
    fname: str = "recharge_LAI_loop.png",
) -> None:
    """Fig 1 — Monthly recharge bar + cumulative line."""
    if recharge_monthly is None or len(recharge_monthly) == 0:
        print("  plot_recharge_bar: skipped (no recharge data).")
        return

    fig, (ax_bar, ax_cum) = plt.subplots(2, 1, figsize=(12, 7), sharex=False)

    recharge_monthly.plot.bar(ax=ax_bar, color="#4C72B0", edgecolor="white")
    ax_bar.set_title("Monthly recharge (LAI loop)", fontsize=14, fontweight="bold")
    ax_bar.set_ylabel("Recharge (mm / month)")
    ax_bar.set_xlabel("")
    ax_bar.tick_params(axis="x", rotation=45)

    recharge_monthly.cumsum().plot(ax=ax_cum, color="#C44E52", lw=2, marker="o", ms=4)
    ax_cum.set_title("Cumulative recharge", fontsize=14)
    ax_cum.set_ylabel("Cumulative recharge (mm)")
    ax_cum.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    idx = recharge_monthly.index
    add_fire_vline(ax_cum, idx.min(), idx.max())
    ax_cum.legend(fontsize=8)

    plt.tight_layout()
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_et_monthly_mean(
    ET_xr_all: xr.Dataset,
    et_var:    str,
    et_scale:  float,
    fig_path:  Path,
    *,
    dpi:  int  = 150,
    show: bool = False,
    fname: str = "ET_monthly_mean.png",
) -> None:
    """Fig 2 — Spatial mean ET per month, one line per year."""
    if et_var not in ET_xr_all:
        print(f"  plot_et_monthly_mean: skipped — '{et_var}' not in ET dataset.")
        return

    et_mean = ET_xr_all[et_var].mean(dim=["X", "Y"]) * et_scale
    years   = np.unique([pd.Timestamp(t.item()).year for t in et_mean.datetime])
    months  = np.arange(1, 13)

    fig, ax = plt.subplots(figsize=(10, 5))
    cmap_lines = plt.get_cmap("tab10", len(years))

    for i, year in enumerate(years):
        mask = np.array([
            pd.Timestamp(t.item()).year == year for t in et_mean.datetime
        ])
        et_year_monthly = et_mean.isel(datetime=mask).groupby("datetime.month").mean()
        ax.plot(
            et_year_monthly.month,
            et_year_monthly.values,
            marker="o", linewidth=2,
            color=cmap_lines(i), label=str(year),
        )

    ax.set_title("Mean Monthly ET — by year")
    ax.set_xlabel("Month")
    ax.set_ylabel("ET (mm/day)")
    ax.set_xticks(months)
    ax.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    ax.legend(title="Year", bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.grid(True)
    plt.tight_layout()
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_et_by_veg_type(
    ET_xr_all:       xr.Dataset,
    veg_map_history: dict,
    et_var:          str,
    et_scale:        float,
    fig_path:        Path,
    *,
    dpi:  int  = 150,
    show: bool = False,
    fname: str = "ET_by_vegetation_trend.png",
) -> None:
    """Fig 3 — ET time series grouped by vegetation type."""

    ET = ET_xr_all[et_var] * et_scale
    ET_monthly_mean = ET.resample(datetime="1MS").mean(dim="datetime")

    time_vals = sorted(veg_map_history.keys())
    veg_types = np.unique(list(veg_map_history.values())[0])
    veg_types = veg_types[~np.isnan(veg_types)]
    veg_et_series: dict[Any, list[float]] = {v: [] for v in veg_types}
    times: list = []

    for t in time_vals:
        try:
            et_t = ET_monthly_mean.sel(datetime=pd.to_datetime(t)).values
        except KeyError:
            print(f"  plot_et_by_veg_type: missing datetime {t}, skipping.")
            continue

        veg_map = veg_map_history[t].T
        veg_map_nodes = mt.map_cells_to_nodes(
            veg_map,
            grid3d_shape=(veg_map.shape[0] + 1, veg_map.shape[1] + 1),
        )
        times.append(t)
        for v in veg_types:
            mask = veg_map_nodes == v
            veg_et_series[v].append(np.nanmean(et_t[mask]) if np.any(mask) else np.nan)

    fig, ax = plt.subplots(figsize=(10, 5))
    for v in veg_types:
        ax.plot(times, veg_et_series[v], marker="o", label=f"Veg {v}")
    ax.set_title("ET trend by vegetation type")
    ax.set_xlabel("Time")
    ax.set_ylabel("ET (mm/day)")
    ax.grid(True)
    if times:
        add_fire_vline(ax, min(times), max(times))
    ax.legend()
    plt.tight_layout()
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)

    veg_et_series[2]
    veg_et_series[1]
     
     
def plot_veg_map_samples(
    veg_map_history: dict,
    cmap_veg: mcolors.Colormap,
    fig_path: Path,
    *,
    dpi:  int  = 150,
    show: bool = False,
    fname: str = "veg_map_samples.png",
) -> None:
    """Fig 5 — One vegetation map sample per year (mosaic layout)."""
    by_year: dict[int, list] = defaultdict(list)
    for d in sorted(veg_map_history):
        by_year[d.year].append(d)

    sample_dates = []
    for yr in sorted(by_year):
        day1 = [d for d in by_year[yr] if d.day == 1]
        sample_dates.append(day1[0] if day1 else by_year[yr][0])

    if not sample_dates:
        print("  plot_veg_map_samples: skipped (no veg_map_history entries).")
        return

    ncols = 2
    nrows = (len(sample_dates) + ncols - 1) // ncols
    mosaic = [
        [str(sample_dates[r * ncols + c].year) if r * ncols + c < len(sample_dates) else "."
         for c in range(ncols)]
        for r in range(nrows)
    ]

    fig, axes = plt.subplot_mosaic(mosaic, figsize=(5 * ncols, 4 * nrows), empty_sentinel=".")
    for d in sample_dates:
        ax = axes[str(d.year)]
        ax.imshow(veg_map_history[d], origin="upper", cmap=cmap_veg,
                  vmin=0.5, vmax=3.5, aspect="auto")
        veg_colorbar(ax, cmap_veg)
        ax.set_title(f"{d:%Y-%m}", fontsize=10)
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")

    plt.suptitle("Vegetation zone maps — one sample per year", fontsize=13, fontweight="bold")
    plt.tight_layout()
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_veg_map_year(
    veg_map_history: dict,
    selected_year:   int,
    cmap_veg:        mcolors.Colormap,
    fig_path:        Path,
    *,
    dpi:  int  = 150,
    show: bool = False,
) -> None:
    """Fig 5b — All monthly vegetation maps for a single year."""
    monthly_dates = sorted(d for d in veg_map_history if d.year == selected_year)
    if not monthly_dates:
        print(f"  plot_veg_map_year: skipped — no entries for {selected_year}.")
        return

    ncols = 4
    nrows = (len(monthly_dates) + ncols - 1) // ncols
    mosaic = [
        [f"{monthly_dates[r * ncols + c]:%Y-%m}"
         if r * ncols + c < len(monthly_dates) else "."
         for c in range(ncols)]
        for r in range(nrows)
    ]

    fig, axes = plt.subplot_mosaic(mosaic, figsize=(5 * ncols, 4 * nrows), empty_sentinel=".")
    for d in monthly_dates:
        ax  = axes[f"{d:%Y-%m}"]
        ax.imshow(veg_map_history[d], origin="upper", cmap=cmap_veg,
                  vmin=0.5, vmax=3.5, aspect="auto")
        veg_colorbar(ax, cmap_veg)
        ax.set_title(f"{d:%B %Y}", fontsize=10)
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")

    plt.suptitle(f"Vegetation zone maps — {selected_year} (monthly)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    save_fig(fig, fig_path, f"veg_map_{selected_year}_monthly.png", dpi)
    maybe_show(fig, show)


def plot_lai_map_year(
    ds_lai:        xr.Dataset,
    var_name:      str,
    selected_year: int,
    cmap:          str,
    vmin:          float,
    vmax:          float,
    fig_path:      Path,
    *,
    time_coord: str  = "time",
    dpi:        int  = 150,
    show:       bool = False,
) -> None:
    """Fig 5c — All monthly LAI maps for a single year."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable  # optional dep

    da      = ds_lai[var_name]
    da_year = da.sel({time_coord: da[time_coord].dt.year == selected_year})

    if da_year.sizes[time_coord] == 0:
        print(f"  plot_lai_map_year: skipped — no entries for {selected_year}.")
        return

    times = da_year[time_coord].values
    ncols = 4
    nrows = (len(times) + ncols - 1) // ncols
    mosaic = [
        [pd.Timestamp(times[r * ncols + c]).strftime("%Y-%m")
         if r * ncols + c < len(times) else "."
         for c in range(ncols)]
        for r in range(nrows)
    ]

    fig, axes = plt.subplot_mosaic(mosaic, figsize=(5 * ncols, 4 * nrows), empty_sentinel=".")

    for t in times:
        ts  = pd.Timestamp(t)
        key = ts.strftime("%Y-%m")
        ax  = axes[key]
        frame = da_year.sel({time_coord: t})

        im = ax.imshow(frame.values, origin="upper", cmap=cmap,
                       vmin=vmin, vmax=vmax, aspect="auto")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax.set_title(ts.strftime("%B %Y"), fontsize=10)
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")

    plt.suptitle(f"LAI monthly maps — {selected_year}", fontsize=13, fontweight="bold")
    plt.tight_layout()
    save_fig(fig, fig_path, f"lai_fromEO_{selected_year}_monthly.png", dpi)
    maybe_show(fig, show)


# =============================================================================
# LAI diagnostic figures  (used by analyse_LAI.py / withLAI_results.py)
# =============================================================================

LAI_UNITS = "m² m⁻²"
LAI_VMAX  = 2.0
CMAP_LAI  = "YlGn"


def lai_history_to_dataset(
    lai_history: dict,
    lai_var: str = "LAI",
) -> xr.Dataset:
    """
    Convert a ``{datetime: 2-D np.ndarray}`` LAI history dict to an
    ``xr.Dataset`` with dimensions ``(time, row, col)``.

    This is the canonical adapter used by all four LAI diagnostic functions
    so they can consume the in-memory artefact loaded from
    ``veg_map_history.pkl`` without touching any NetCDF file.

    Parameters
    ----------
    lai_history : dict mapping ``datetime`` → 2-D float array (row, col)
    lai_var     : name to give the data variable (default: 'LAI')

    Returns
    -------
    xr.Dataset with one variable *lai_var* and dims (time, row, col).
    Returns an empty Dataset (no time steps) when *lai_history* is empty.
    """
    if not lai_history:
        return xr.Dataset()

    sorted_dates = sorted(lai_history.keys())
    arrays       = np.stack([lai_history[d] for d in sorted_dates], axis=0)  # (T, row, col)

    da = xr.DataArray(
        arrays,
        dims=["time", "row", "col"],
        coords={"time": pd.DatetimeIndex(sorted_dates)},
        name=lai_var,
    )
    return da.to_dataset()


def plot_lai_monthly_by_year(
    ds: xr.Dataset,
    lai_var: str,
    fig_path: Path,
    *,
    dpi:  int  = 150,
    show: bool = False,
    fname: str = "LAI_monthly_mean_by_year.png",
) -> None:
    """
    Spatial-mean LAI per calendar month, one coloured line per year.

    Mirrors ``plot_eta_monthly_by_year`` in style.

    Parameters
    ----------
    ds       : xr.Dataset containing *lai_var* with a 'time' dimension
    lai_var  : variable name inside *ds* (e.g. 'LAI')
    fig_path : output directory (created if absent)
    dpi      : figure resolution
    show     : display interactively before closing
    fname    : output filename
    """
    da       = _da(ds, lai_var)
    time_dim = _infer_time_dim(ds)
    sp_dims  = _spatial_dims(da, time_dim)

    mean_ts = da.mean(dim=sp_dims, skipna=True)
    ts      = mean_ts.to_series()
    years   = sorted(ts.index.year.unique())
    cmap    = plt.get_cmap("tab10", len(years))

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, yr in enumerate(years):
        yr_ts   = ts[ts.index.year == yr]
        monthly = yr_ts.groupby(yr_ts.index.month).mean()
        ax.plot(monthly.index, monthly.values,
                marker="o", lw=2, color=cmap(i), label=str(yr))

    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    ax.set(title=f"Spatial-mean LAI — monthly average by year  [{lai_var}]",
           xlabel="Month", ylabel=f"LAI ({LAI_UNITS})")
    ax.legend(title="Year", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_lai_spatial_maps_year(
    ds: xr.Dataset,
    lai_var: str,
    selected_year: int,
    fig_path: Path,
    *,
    cmap:  str        = CMAP_LAI,
    vmin:  float      = 0.0,
    vmax:  float      = LAI_VMAX,
    dpi:   int        = 150,
    show:  bool       = False,
    fname: str | None = None,
) -> None:
    """
    Monthly spatial maps (4-column mosaic) for *selected_year*.

    Mirrors ``plot_eta_spatial_maps_year`` in style.

    Parameters
    ----------
    ds            : xr.Dataset containing *lai_var*
    lai_var       : variable name (e.g. 'LAI')
    selected_year : calendar year to display
    fig_path      : output directory (created if absent)
    cmap          : matplotlib colourmap (default: 'YlGn')
    vmin / vmax   : colour-scale limits
    dpi           : figure resolution
    show          : display interactively before closing
    fname         : output filename; auto-generated from year if None
    """
    da       = _da(ds, lai_var)
    time_dim = _infer_time_dim(ds)

    da_year = da.sel({time_dim: da[time_dim].dt.year == selected_year})
    if da_year.sizes[time_dim] == 0:
        print(f"  plot_lai_spatial_maps_year: no data for {selected_year} — skipped.")
        return

    times = [pd.Timestamp(t) for t in da_year[time_dim].values]
    ncols = 4
    nrows = (len(times) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(5 * ncols, 4 * nrows), squeeze=False)

    for idx, t in enumerate(times):
        r, c  = divmod(idx, ncols)
        ax    = axes[r][c]
        frame = da_year.sel({time_dim: t}).values
        im    = ax.imshow(frame, origin="upper", cmap=cmap,
                          vmin=vmin, vmax=vmax, aspect="auto")
        plt.colorbar(im, ax=ax, label=LAI_UNITS, fraction=0.046, pad=0.04)
        ax.set_title(t.strftime("%B %Y"), fontsize=10)
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")

    for idx in range(len(times), nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    fig.suptitle(f"LAI monthly spatial maps — {selected_year}  [{lai_var}]",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    if fname is None:
        fname = f"LAI_spatial_maps_{selected_year}.png"
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_lai_annual_hist(
    ds: xr.Dataset,
    lai_var: str,
    fig_path: Path,
    *,
    bins:  int  = 40,
    dpi:   int  = 150,
    show:  bool = False,
    fname: str  = "LAI_annual_mean_histogram.png",
) -> None:
    """
    Histogram of per-pixel annual mean LAI (all years combined).

    One coloured series per year, matching the style of
    ``plot_eta_annual_hist``.

    Parameters
    ----------
    ds       : xr.Dataset containing *lai_var*
    lai_var  : variable name (e.g. 'LAI')
    fig_path : output directory (created if absent)
    bins     : number of histogram bins
    dpi      : figure resolution
    show     : display interactively before closing
    fname    : output filename
    """
    da       = _da(ds, lai_var)
    time_dim = _infer_time_dim(ds)
    times    = pd.DatetimeIndex(da[time_dim].values)
    years    = sorted(set(times.year))

    fig, ax = plt.subplots(figsize=(9, 5))
    cmap_lines = plt.get_cmap("tab10", len(years))

    for i, year in enumerate(years):
        mask       = times.year == year
        da_year    = da.isel({time_dim: mask})
        annual_mean = da_year.mean(dim=time_dim, skipna=True).values.ravel()
        annual_mean = annual_mean[np.isfinite(annual_mean)]
        if annual_mean.size == 0:
            continue
        ax.hist(annual_mean, bins=bins, alpha=0.55,
                color=cmap_lines(i), edgecolor="none",
                label=str(year), density=False)

    ax.set_xlabel(f"Annual mean LAI ({LAI_UNITS})", fontsize=11)
    ax.set_ylabel("Pixel count", fontsize=11)
    ax.set_title(f"Per-pixel annual mean LAI distribution  [{lai_var}]",
                 fontsize=13, fontweight="bold")
    ax.legend(title="Year", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_lai_annual_mean(
    ds: xr.Dataset,
    lai_var: str,
    fig_path: Path,
    *,
    label: str       = "",
    dpi:   int       = 150,
    show:  bool      = False,
    fname: str | None = None,
) -> None:
    """
    Bar chart of domain-mean annual LAI with an overlaid trend line.

    Mirrors ``plot_eta_total`` in style.  Produces one figure:
        LAI_annual_mean[_<label>].png

    Parameters
    ----------
    ds       : xr.Dataset containing *lai_var*
    lai_var  : variable name (e.g. 'LAI')
    fig_path : output directory (created if absent)
    label    : plot-selection label appended to the filename
    dpi      : figure resolution
    show     : display interactively before closing
    fname    : explicit output filename; auto-generated if None
    """
    da       = _da(ds, lai_var)
    time_dim = _infer_time_dim(ds)
    sp_dims  = _spatial_dims(da, time_dim)

    spatial_mean = da.mean(dim=sp_dims, skipna=True)
    annual       = spatial_mean.groupby(f"{time_dim}.year").mean(time_dim).to_series()

    norm   = mcolors.Normalize(vmin=float(annual.min()), vmax=float(annual.max()))
    cmap   = cm.get_cmap("YlGn")
    colors = [cmap(norm(v)) for v in annual.values]

    fig, ax = plt.subplots(figsize=(max(8, len(annual) * 0.8), 5))
    bars = ax.bar(annual.index.astype(str), annual.values,
                  color=colors, edgecolor="white", width=0.6, zorder=2)

    ax.plot(range(len(annual)), annual.values,
            marker="o", color="#006837", lw=2, zorder=5)

    mean_val = float(annual.mean())
    ax.axhline(mean_val, color="#C44E52", lw=1.5, linestyle="--",
               label=f"mean  {mean_val:.2f}", zorder=3)

    # colourbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.01, shrink=0.85)
    cb.set_label(f"LAI ({LAI_UNITS})", fontsize=9)

    for bar, val in zip(bars, annual.values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + LAI_VMAX * 0.015,
                f"{val:.2f}", ha="center", va="bottom", fontsize=9)

    ax.set(title=f"Domain-mean annual LAI{' — ' + label if label else ''}  [{lai_var}]",
           xlabel="Year", ylabel=f"LAI ({LAI_UNITS})")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(True, axis="y", linestyle="--", alpha=0.4, zorder=1)
    ax.legend(fontsize=9)
    ax.set_ylim(0, float(annual.max()) * 1.18)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    if fname is None:
        suffix = f"_{label}" if label else ""
        fname  = f"LAI_annual_mean{suffix}.png"
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_lai_annual_median(
    ds: xr.Dataset,
    lai_var: str,
    fig_path: Path,
    *,
    label: str        = "",
    dpi:   int        = 150,
    show:  bool       = False,
    fname: str | None = None,
) -> None:
    """
    Yearly time series of domain-median LAI (one point per year).

    Companion to ``plot_lai_annual_mean``: uses the median instead of the
    mean, which is less sensitive to outlier pixels (e.g. bare/burned
    patches skewing a spatial average). Drawn on a real datetime x-axis
    (each point anchored at 1 July of its year) so the fire date can be
    marked directly on the series via ``add_fire_vline``.

    Produces one figure:
        LAI_annual_median[_<label>].png

    Parameters
    ----------
    ds       : xr.Dataset containing *lai_var*
    lai_var  : variable name (e.g. 'LAI')
    fig_path : output directory (created if absent)
    label    : plot-selection label appended to the filename
    dpi      : figure resolution
    show     : display interactively before closing
    fname    : explicit output filename; auto-generated if None
    """
    da       = _da(ds, lai_var)
    time_dim = _infer_time_dim(ds)
    sp_dims  = _spatial_dims(da, time_dim)

    spatial_median = da.median(dim=sp_dims, skipna=True)
    annual = spatial_median.groupby(f"{time_dim}.year").median(time_dim).to_series()

    if annual.empty:
        print("  plot_lai_annual_median: skipped (no data).")
        return

    years   = annual.index.astype(int)
    x_dates = pd.to_datetime([f"{y}-07-01" for y in years])  # mid-year anchor

    fig, ax = plt.subplots(figsize=(max(8, len(annual) * 0.8), 5))
    ax.plot(x_dates, annual.values, marker="o", color="#004529", lw=2,
            zorder=5, label=f"{lai_var} median")

    median_val = float(annual.mean())
    ax.axhline(median_val, color="#C44E52", lw=1.5, linestyle="--",
               label=f"mean of yearly medians  {median_val:.2f}", zorder=3)

    add_fire_vline(ax, x_dates.min(), x_dates.max())

    for x, val in zip(x_dates, annual.values):
        ax.text(x, val, f"{val:.2f}", ha="center", va="bottom", fontsize=9)

    ax.set(title=f"Domain-median annual LAI{' — ' + label if label else ''}  [{lai_var}]",
           xlabel="Year", ylabel=f"LAI ({LAI_UNITS})")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(fontsize=9)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    if fname is None:
        suffix = f"_{label}" if label else ""
        fname  = f"LAI_annual_median{suffix}.png"
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def make_veg_gif(
    veg_history: dict,
    cmap_veg:    mcolors.Colormap,
    out_path:    Path,
    *,
    duration: float = 0.6,
) -> None:
    """GIF A — Animate monthly vegetation zone maps."""
    dates   = sorted(veg_history.keys())
    tmp_dir = Path(tempfile.mkdtemp())
    frames  = []

    for i, d in enumerate(dates):
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.imshow(veg_history[d], origin="upper", cmap=cmap_veg,
                  vmin=0.5, vmax=3.5, aspect="auto")
        veg_colorbar(ax, cmap_veg)
        ax.set_title(f"Vegetation map — {d:%Y-%m}", fontsize=11)
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")
        plt.tight_layout()

        fp = tmp_dir / f"veg_{i:04d}.png"
        fig.savefig(fp, dpi=100, bbox_inches="tight")
        plt.close(fig)
        frames.append(fp)

    with imageio.get_writer(str(out_path), mode="I", duration=duration) as writer:
        for fp in frames:
            writer.append_data(imageio.imread(fp))

    for fp in frames:
        fp.unlink()
    tmp_dir.rmdir()
    print(f"  Saved → {out_path}")


def make_et_gif(
    et_dataset: xr.Dataset,
    var:        str,
    scale:      float,
    out_path:   Path,
    *,
    cmap:     str   = CMAP_ET,
    duration: float = 0.4,
) -> None:
    """GIF B — Animate ET spatial maps from an (X, Y, datetime) DataArray."""
    data_var = et_dataset[var]
    vmin = float(data_var.min()) * scale
    vmax = float(data_var.max()) * scale

    tmp_dir = Path(tempfile.mkdtemp())
    frames  = []

    for t_idx in range(len(data_var["datetime"])):
        et_da = data_var.isel(datetime=t_idx) * scale
        t_val = data_var["datetime"].isel(datetime=t_idx).values

        fig, ax = plt.subplots(figsize=(7, 5))
        et_da.plot.imshow(
            ax=ax, x="X", y="Y",
            cmap=cmap, vmin=vmin, vmax=vmax,
            cbar_kwargs={"label": "ET (mm/day)"},
            add_colorbar=True,
        )
        ax.set_title(time_label(t_val), fontsize=11)
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        plt.tight_layout()

        fp = tmp_dir / f"frame_{t_idx:04d}.png"
        fig.savefig(fp, dpi=100, bbox_inches="tight")
        plt.close(fig)
        frames.append(fp)

    with imageio.get_writer(str(out_path), mode="I", duration=duration) as writer:
        for fp in frames:
            writer.append_data(imageio.imread(fp))

    for fp in frames:
        fp.unlink()
    tmp_dir.rmdir()
    print(f"  Saved → {out_path}")


# =============================================================================
# SSHydro figures   (used by LT_SShydro_results.py)
# =============================================================================

def plot_hydro(
    data:       dict,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    date_start:  str | None = None,
    date_end:    str | None = None,
    show_rain:   bool = True,
    show_psi:    bool = True,
    show_sw:     bool = True,
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    dpi:         int  = 300,
    show:        bool = False,
    fname:       str  = "hydro_timeseries.png",
) -> None:
    """Plot Rain/ETp + ψ + sw time series via AgUtils."""
    import Agramon_utils as AgUtils  # project-local import

    colors_dict = colors_dict or COLORS_HYDRO
    d0, d1 = date_window(date_start, date_end, data["df_dates"])
    print(f"  Plotting hydro time-series ({d0.date()} → {d1.date()}) …")

    fig = AgUtils.plot_hydro_notebook(
        date_start=d0,
        date_end=d1,
        node_data=nodes_dict,
        pev_series=data["pev_series"],
        tp_series=data["tp_series"],
        df_psi=data["df_psi"],
        df_sw=data["df_sw"],
        df_dates=data["df_dates"],
        df_dates_sw=data["df_dates_sw"],
        show_rain=show_rain,
        show_psi=show_psi,
        show_sw=show_sw,
        show_uh_s=show_uh_s,
        show_uh_1=show_uh_1,
        show_dh_s=show_dh_s,
        show_dh_1=show_dh_1,
        colors_dict=colors_dict,
        scenario_id=data["scenario"],
    )
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_et_sw(
    data:       dict,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    date_start:  str | None = None,
    date_end:    str | None = None,
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    dpi:         int  = 300,
    show:        bool = False,
    fname:       str  = "et_timeseries.png",
) -> None:
    """Plot soil water content + actual ET time series."""
    colors_dict = colors_dict or COLORS_HYDRO
    d0, d1 = date_window(date_start, date_end, data["df_dates_sw"])
    print(f"  Plotting ET+sw time-series ({d0.date()} → {d1.date()}) …")

    df_sw       = data["df_sw"]
    df_dates_sw = data["df_dates_sw"]
    ET_data     = data["ET_xr_mean"]["ACT. ETRA"].values
    time_dates  = data["time_dates"]
    scenario_id = data["scenario"]
    depth_label = nodes_dict.get("depth_label", 1)

    def get_node_id(key: str) -> int:
        return nodes_dict[key][0][0]

    uh_s = get_node_id("mid_z-0m")
    uh_1 = get_node_id("mid_z-1m")
    dh_s = get_node_id("outlet_z-0m")
    dh_1 = get_node_id("outlet_z-1m")

    sw_m = (df_dates_sw >= d0) & (df_dates_sw <= d1)
    et_m = (time_dates  >= d0) & (time_dates  <= d1)

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

    if show_uh_s:
        ax0.plot(df_dates_sw[sw_m], df_sw.loc[sw_m, uh_s],
                 label=f"Mid surf n={uh_s}",
                 color=colors_dict["uphill_surface"], marker="+", linestyle="-")
    if show_uh_1:
        ax0.plot(df_dates_sw[sw_m], df_sw.loc[sw_m, uh_1],
                 label=f"Mid −{depth_label}m n={uh_1}",
                 color=colors_dict["uphill_1m"], linestyle="--")
    if show_dh_s:
        ax0.plot(df_dates_sw[sw_m], df_sw.loc[sw_m, dh_s],
                 label=f"Outlet surf n={dh_s}",
                 color=colors_dict["downhill_surface"], marker="o", linestyle="-", markersize=3)
    if show_dh_1:
        ax0.plot(df_dates_sw[sw_m], df_sw.loc[sw_m, dh_1],
                 label=f"Outlet −{depth_label}m n={dh_1}",
                 color=colors_dict["downhill_1m"], linestyle="--")

    ax0.set_ylabel("sw (−)")
    ax0.set_title("Soil Water Content", fontsize=10, fontweight="bold")
    add_fire_vline(ax0, d0, d1)
    ax0.legend(loc="upper right", fontsize=8)
    ax0.grid(True, linestyle="--", alpha=0.3)

    ax1.scatter(time_dates[et_m], ET_data[et_m],
                c=colors_dict["ET_mean"], s=14, label="ETa mean")
    ax1.set_ylabel("ETa (mm/day)")
    ax1.set_xlabel("Date")
    ax1.set_title("Actual Evapotranspiration (spatial mean)", fontsize=10, fontweight="bold")
    add_fire_vline(ax1, d0, d1)
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(True, linestyle="--", alpha=0.3)
    ax1.xaxis.set_major_formatter(mdates.ConciseDateFormatter(mdates.AutoDateLocator()))

    fig.suptitle(f"Sc.{scenario_id} | sw & ETa | {d0.date()} → {d1.date()}",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_era5(
    data:      dict,
    fig_path:  Path,
    *,
    date_start: str | None = None,
    date_end:   str | None = None,
    dpi:        int  = 300,
    show:       bool = False,
    fname:      str  = "era5_forcing.png",
) -> None:
    """ERA5 forcing bar chart (ETp + Rain)."""
    pev_series = data["pev_series"]
    tp_series  = data["tp_series"]
    d0, d1 = date_window(date_start, date_end, pev_series.index.to_series())
    print(f"  Plotting ERA5 forcing ({d0.date()} → {d1.date()}) …")

    mask = (pev_series.index >= d0) & (pev_series.index <= d1)
    fig, ax = plt.subplots(figsize=(14, 3))
    ax.bar(pev_series[mask].index, pev_series[mask].values,
           width=1.0, color="grey", label="ETp", alpha=0.5)
    ax.bar(tp_series[mask].index, tp_series[mask].values,
           width=1.0, color="skyblue", label="Rain")
    ax.invert_yaxis()
    ax.set_ylabel("mm/day")
    ax.set_title(f"ERA5 Forcing  |  {d0.date()} → {d1.date()}", fontweight="bold")
    add_fire_vline(ax, d0, d1)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.4)
    locator   = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    fig.autofmt_xdate()
    plt.tight_layout()
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_catchment_et_freq(
    ET_xr:      xr.Dataset,
    ET_xr_mean: xr.Dataset,
    time_dates: np.ndarray,
    freq:       str,
    fig_path:   Path,
    *,
    dpi:  int  = 300,
    show: bool = False,
) -> None:
    """Catchment ETa aggregated at *freq* frequency ('ME' = monthly, 'YE' = annual)."""
    # 1. Catchment area
    res_x, res_y    = ET_xr.rio.resolution()
    pixel_area       = abs(float(res_x) * float(res_y))
    variable_name    = list(ET_xr.data_vars)[0]
    _t               = _infer_time_dim(ET_xr)
    valid_pixels     = int(ET_xr[variable_name].notnull().isel({_t: 0}).sum().item())
    area_m2          = valid_pixels * pixel_area

    # 2. Resample
    df_et = pd.DataFrame(
        {"ETa_mm": ET_xr_mean["ACT. ETRA"].values},
        index=pd.to_datetime(time_dates),
    )
    df_resampled = df_et.resample(freq).sum()
    df_resampled["ETa_m3"] = (df_resampled["ETa_mm"] / 1000.0) * area_m2

    # 3. Labels
    p_str       = "monthly" if freq == "ME" else "yearly"
    time_labels = (
        df_resampled.index.strftime("%Y-%m") if freq == "ME"
        else df_resampled.index.strftime("%Y")
    )
    values = df_resampled["ETa_mm"]

    # 4. Figure
    n_bars = len(values)
    fig_w  = max(12, n_bars * 0.18)   # ~0.18 in per bar, minimum 12
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    bars = ax.bar(time_labels, values, color="skyblue", edgecolor="navy", alpha=0.7)
    ax.set_ylabel("ETa accumulated (mm)")
    ax.set_title(
        f"Catchment ETa — {p_str.capitalize()} (Area: {area_m2:,.1f} m²)",
        fontweight="bold",
    )
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    # ── Fire marker ────────────────────────────────────────────────────────
    # x is categorical (bar labels), so locate the matching label's position
    # rather than using a real datetime value.
    fire_label = "2020-07" if freq == "ME" else "2020"
    all_labels_list = list(time_labels)
    if fire_label in all_labels_list:
        fire_pos = all_labels_list.index(fire_label)
        ax.axvline(fire_pos, color="#B22222", linestyle="--", linewidth=1.3,
                   alpha=0.85, zorder=4, label="Fire (Jul 2020)")
        ax.legend(fontsize=8, loc="upper right")

    # ── X-axis: show only January ticks labelled as the year ─────────────────
    if freq == "ME":
        all_labels     = list(time_labels)
        tick_positions = [i for i, lbl in enumerate(all_labels) if lbl.endswith("-01")]
        tick_labels    = [lbl[:4] for lbl in all_labels if lbl.endswith("-01")]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=0, ha="center", fontsize=10)
        # light vertical rule at each year boundary
        for pos in tick_positions:
            ax.axvline(pos - 0.5, color="gray", lw=0.5, alpha=0.4, zorder=0)
    else:
        ax.tick_params(axis="x", rotation=45)

    # ── Bar value labels (skip zeros) ─────────────────────────────────────────
    for bar in bars:
        yval = bar.get_height()
        if yval > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, yval,
                    f"{yval:,.1f}", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()


def plot_catchment_et(
    data:     dict,
    fig_path: Path,
    *,
    dpi:  int  = 300,
    show: bool = False,
) -> None:
    """Monthly + annual catchment ET summaries (wrapper)."""
    print("  Computing monthly catchment ET …")
    plot_catchment_et_freq(
        data["ET_xr"], data["ET_xr_mean"], data["time_dates"],
        freq="ME", fig_path=fig_path, dpi=dpi, show=show,
    )
    print("  Computing annual catchment ET …")
    plot_catchment_et_freq(
        data["ET_xr"], data["ET_xr_mean"], data["time_dates"],
        freq="YE", fig_path=fig_path, dpi=dpi, show=show,
    )


def plot_hgraph(
    data:     dict,
    fig_path: Path,
    *,
    dpi:   int  = 300,
    show:  bool = False,
    fname: str  = "hgraph.png",
) -> None:
    """Hydraulic graph from pyCATHY simu.show()."""
    print("  Plotting hydraulic graph …")
    fig, ax = plt.subplots()
    data["simu"].show(prop="hgraph", ax=ax)
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_node_locations(
    grid3d:      dict,
    nodes_dict:  dict,
    fig_path:    Path,
    *,
    dem_da:      xr.DataArray | None = None,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig_node_locations.png",
) -> None:
    """
    Fig 6b — 2-D map of the reference-node locations where the psi and sw
    time-series (Fig 6 / Fig 7) are extracted.

    Nodes sharing the same site (e.g. every depth resolved at "outlet" or
    "mid") sit at the same (X, Y) and are plotted together, annotated with
    their depth tag and node index.

    Parameters
    ----------
    grid3d      : dict returned by ``simu.read_outputs("grid3d")``; must
                  contain "mesh3d_nodes" (array of X, Y, Z per node index).
    nodes_dict  : dict returned by ``resolve_nodes()``, i.e.
                  {"outlet_z-0m": <node id>, "mid_z-1m": <node id>, …}.
    fig_path    : output directory.
    dem_da      : optional xr.DataArray (dims "x"/"y" or "X"/"Y") used as a
                  DEM basemap. Skipped silently if it can't be plotted.
    colors_dict : override COLORS_HYDRO; "uphill_surface" colours "mid"
                  sites, "downhill_surface" colours "outlet" sites.
    scenario_id : label for the figure title.
    """
    colors_dict = colors_dict or COLORS_HYDRO
    nodes = grid3d["mesh3d_nodes"]

    # ── Group node labels by site (everything before "_z-") ───────────────────
    sites: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for label in nodes_dict:
        if "_z-" not in label:
            continue
        site = label.split("_z-")[0]
        sites[site].append((label, _node_id(nodes_dict, label)))

    if not sites:
        print("  plot_node_locations: skipped — no depth-tagged nodes in nodes_dict.")
        return

    site_color_map = {
        "mid":    colors_dict.get("uphill_surface",   "#0072B2"),
        "outlet": colors_dict.get("downhill_surface", "#D55E00"),
    }
    fallback_cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(7, 7))

    # ── Optional DEM basemap ────────────────────────────────────────────────
    if dem_da is not None:
        try:
            xdim = "x" if "x" in dem_da.dims else "X"
            ydim = "y" if "y" in dem_da.dims else "Y"
            dem_da.plot.imshow(
                ax=ax, x=xdim, y=ydim, cmap="terrain", alpha=0.85,
                add_colorbar=True,
                cbar_kwargs={"label": "Elevation (m)", "shrink": 0.75},
            )
            ax.set_title("")
        except Exception as exc:
            print(f"  plot_node_locations: DEM basemap skipped ({exc}).")

    # ── Reference-node markers ──────────────────────────────────────────────
    for i, (site, entries) in enumerate(sites.items()):
        color = site_color_map.get(site, fallback_cmap(i % 10))
        xs, ys = [], []
        for label, node_id in entries:
            x, y, _z = nodes[node_id]
            xs.append(x)
            ys.append(y)
            depth_tag = label.split("_z-")[1]
            ax.annotate(
                f"{depth_tag}\n(n={node_id})",
                (x, y), textcoords="offset points", xytext=(6, 6),
                fontsize=7, color=color,
            )
        ax.scatter(
            xs, ys, s=70, color=color, edgecolor="black",
            linewidth=0.6, zorder=3, label=site.capitalize(),
        )

    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_aspect("equal")
    ax.set_title(
        f"Sc.{scenario_id} | Reference node locations (psi / sw extraction)",
        fontsize=11, fontweight="bold",
    )
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.3)

    plt.tight_layout()
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_recharge(
    data:        dict,
    fig_path:    Path,
    *,
    dpi:         int  = 300,
    show:        bool = False,
    scenario_id: int | str | None = None,
) -> None:
    """Monthly + annual recharge bar charts from CATHY recharge output."""
    print("  Plotting recharge …")

    simu        = data["simu"]
    start_date  = str(data["start_date"].date())
    scenario_id = scenario_id if scenario_id is not None else data.get("scenario", "?")

    # ── Load & process ──────────────────────────────────────────────────────
    df_recharge  = simu.read_outputs("recharge")
    xr_recharge  = df_recharge.set_index(["time", "X", "Y"]).to_xarray()
    recharge_ms  = xr_recharge["recharge"].mean(dim=["X", "Y"])

    time_sec      = xr_recharge.time.values / np.timedelta64(1, "s")
    dt            = np.diff(time_sec, prepend=0)
    t0            = pd.Timestamp(start_date)
    time_datetime = t0 + pd.to_timedelta(xr_recharge.time.values)

    recharge_m          = pd.Series(recharge_ms.values * dt, index=time_datetime)
    recharge_monthly_mm = recharge_m.resample("ME").sum() * 1000
    recharge_yearly_mm  = recharge_m.resample("YE").sum() * 1000

    # ── Monthly figure ──────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(12, 6))
    recharge_monthly_mm.plot.bar(ax=axes[0])
    axes[0].set_title("Monthly recharge")
    axes[0].set_ylabel("Recharge (mm/month)")
    axes[0].set_xlabel("")
    recharge_monthly_mm.cumsum().plot(ax=axes[1])
    axes[1].set_title("Cumulative recharge")
    axes[1].set_ylabel("Cumulative recharge (mm)")
    idx_m = recharge_monthly_mm.index
    add_fire_vline(axes[1], idx_m.min(), idx_m.max())
    axes[1].legend(fontsize=8)
    fig.suptitle(f"Scenario: {scenario_id}", fontweight="bold")
    plt.tight_layout()
    save_fig(fig, fig_path, "recharge.png", dpi)
    maybe_show(fig, show)

    # ── Annual figure ───────────────────────────────────────────────────────
    fig2, axes2 = plt.subplots(2, 1, figsize=(10, 6))
    year_labels = [str(t.year) for t in recharge_yearly_mm.index]

    axes2[0].bar(range(len(recharge_yearly_mm)), recharge_yearly_mm.values,
                 color="steelblue", edgecolor="white")
    axes2[0].set_xticks(range(len(recharge_yearly_mm)))
    axes2[0].set_xticklabels(year_labels, rotation=45, ha="right")
    axes2[0].set_title("Annual recharge")
    axes2[0].set_ylabel("Recharge (mm/year)")
    axes2[0].grid(True, axis="y", linestyle="--", alpha=0.4)

    axes2[1].plot(recharge_yearly_mm.index, recharge_yearly_mm.cumsum().values,
                  marker="o", color="steelblue")
    axes2[1].set_title("Cumulative annual recharge")
    axes2[1].set_ylabel("Cumulative recharge (mm)")
    axes2[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    add_fire_vline(axes2[1], recharge_yearly_mm.index.min(), recharge_yearly_mm.index.max())
    axes2[1].legend(fontsize=8)
    axes2[1].grid(True, linestyle="--", alpha=0.4)

    fig2.suptitle(f"Scenario: {scenario_id} — Annual recharge", fontweight="bold")
    plt.tight_layout()
    save_fig(fig2, fig_path, "recharge_yearly.png", dpi)
    maybe_show(fig2, show)

def plot_recharge_xr(
    xr_recharge: xr.Dataset,
    fig_path:    Path,
    *,
    dpi:         int  = 300,
    show:        bool = False,
    scenario_id: int | str | None = None,
    start_date:  str | pd.Timestamp | None = None,
) -> None:
    """Monthly + annual recharge bar charts from a pre-saved recharge_output.nc.

    This function is the NetCDF-based counterpart of ``plot_recharge`` and
    replaces ``simu.read_outputs("recharge")`` with a direct xarray load.

    The dataset is expected to have:
      - a ``datetime`` (or ``time``) dimension with absolute datetimes,
      - ``X`` and ``Y`` spatial dimensions,
      - a recharge variable (the first data variable is used if the name is
        not one of the recognised candidates).

    Parameters
    ----------
    xr_recharge : xr.Dataset
        Dataset loaded from ``recharge_output.nc`` (saved by
        ``test_Agramon_withLAI.py`` via ``xr.concat(all_recharge_recs, ...)``.
    fig_path    : Path to the output figure directory.
    dpi         : Figure resolution.
    show        : If True, call ``plt.show()`` after saving.
    scenario_id : Label for the figure suptitle.
    """
    print("  Plotting recharge from recharge_output.nc …")
    scenario_id = scenario_id if scenario_id is not None else "?"

    # ── Identify variable and time dimension ──────────────────────────────────
    # Prefer a variable whose name contains "recharge"; fall back to the first.
    rec_var = next(
        (v for v in xr_recharge.data_vars if "recharge" in v.lower()),
        list(xr_recharge.data_vars)[0],
    )
    time_dim = next(
        (d for d in xr_recharge.dims if d.lower() in ("datetime", "time", "t")),
        list(xr_recharge.dims)[0],
    )

    da = xr_recharge[rec_var]   # (time, X, Y)

    # ── Spatial mean → 1-D time series (m/s) ─────────────────────────────────
    spatial_dims = [d for d in da.dims if d != time_dim]
    recharge_ms  = da.mean(dim=spatial_dims)

    # ── Build a DatetimeIndex from the time coordinate ────────────────────────
    time_vals = xr_recharge[time_dim].values
    if np.issubdtype(time_vals.dtype, np.timedelta64):
        # Relative timedeltas — convert to seconds and integrate Δt
        time_sec = time_vals.astype("timedelta64[s]").astype(float)
        dt_s     = np.diff(time_sec, prepend=0.0)
        # Resolve absolute origin: explicit arg > dataset attribute > warn & use 1970
        if start_date is not None:
            origin = pd.Timestamp(start_date)
        elif "start_date" in xr_recharge.attrs:
            origin = pd.Timestamp(xr_recharge.attrs["start_date"])
        elif "time_origin" in xr_recharge.attrs:
            origin = pd.Timestamp(xr_recharge.attrs["time_origin"])
        else:
            import warnings
            warnings.warn(
                "plot_recharge_xr: time coordinate is relative (timedelta64) but no "
                "start_date was provided and none found in dataset attributes. "
                "Dates will be anchored to 1970-01-01. Pass start_date=<your_date> "
                "to fix this.",
                UserWarning, stacklevel=2,
            )
            origin = pd.Timestamp("1970-01-01")
        time_datetime = pd.DatetimeIndex(
            [origin + pd.to_timedelta(td) for td in time_vals]
        )
    else:
        time_datetime = pd.DatetimeIndex(time_vals)
        time_sec = (time_datetime - time_datetime[0]).total_seconds().values
        dt_s     = np.diff(time_sec, prepend=0.0)

    # ── Integrate: recharge volume [m] per time step ──────────────────────────
    recharge_m          = pd.Series(recharge_ms.values * dt_s, index=time_datetime)
    recharge_monthly_mm = recharge_m.resample("ME").sum() * 1000   # → mm/month
    recharge_yearly_mm  = recharge_m.resample("YE").sum() * 1000   # → mm/year

    # ── Monthly figure — all years plotted by calendar month ────────────────
    month_names  = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    all_years    = sorted(recharge_monthly_mm.index.year.unique())
    cmap_years   = plt.get_cmap("tab10", len(all_years))

    fig, axes = plt.subplots(2, 1, figsize=(10, 6))

    for i, yr in enumerate(all_years):
        col     = cmap_years(i)
        mask    = recharge_monthly_mm.index.year == yr
        yr_data = recharge_monthly_mm[mask]
        if yr_data.empty:
            continue
        months = yr_data.index.month
        axes[0].plot(months, yr_data.values, marker="o", lw=2, color=col, label=str(yr))

    axes[0].set_title("Monthly recharge — all years")
    axes[0].set_ylabel("Recharge (mm/month)")
    axes[0].set_xlabel("")
    axes[0].set_xticks(range(1, 13))
    axes[0].set_xticklabels(month_names, fontsize=9)
    axes[0].grid(axis="y", linestyle="--", alpha=0.3)
    axes[0].legend(title="Year", fontsize=9, bbox_to_anchor=(1.01, 1), loc="upper left")

    for i, yr in enumerate(all_years):
        col     = cmap_years(i)
        mask    = recharge_monthly_mm.index.year == yr
        yr_data = recharge_monthly_mm[mask]
        if yr_data.empty:
            continue
        axes[1].plot(yr_data.index.month, yr_data.cumsum().values,
                     marker="o", lw=2, color=col, label=str(yr))

    axes[1].set_title("Cumulative recharge (within year)")
    axes[1].set_ylabel("Cumulative recharge (mm)")
    axes[1].set_xticks(range(1, 13))
    axes[1].set_xticklabels(month_names, fontsize=9)
    axes[1].grid(linestyle="--", alpha=0.3)
    axes[1].legend(title="Year", fontsize=9, bbox_to_anchor=(1.01, 1), loc="upper left")

    fig.suptitle(f"Scenario: {scenario_id}", fontweight="bold")
    plt.tight_layout()
    save_fig(fig, fig_path, "recharge_monthly_all_years.png", dpi)
    maybe_show(fig, show)
    fig2, axes2 = plt.subplots(2, 1, figsize=(10, 6))
    year_labels = [str(t.year) for t in recharge_yearly_mm.index]

    axes2[0].bar(range(len(recharge_yearly_mm)), recharge_yearly_mm.values,
                 color="steelblue", edgecolor="white")
    axes2[0].set_xticks(range(len(recharge_yearly_mm)))
    axes2[0].set_xticklabels(year_labels, rotation=45, ha="right")
    axes2[0].set_title("Annual recharge")
    axes2[0].set_ylabel("Recharge (mm/year)")
    axes2[0].grid(True, axis="y", linestyle="--", alpha=0.4)

    axes2[1].plot(recharge_yearly_mm.index, recharge_yearly_mm.cumsum().values,
                  marker="o", color="steelblue")
    axes2[1].set_title("Cumulative annual recharge")
    axes2[1].set_ylabel("Cumulative recharge (mm)")
    axes2[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    add_fire_vline(axes2[1], recharge_yearly_mm.index.min(), recharge_yearly_mm.index.max())
    axes2[1].legend(fontsize=8)
    axes2[1].grid(True, linestyle="--", alpha=0.4)

    fig2.suptitle(f"Scenario: {scenario_id} — Annual recharge", fontweight="bold")
    plt.tight_layout()
    save_fig(fig2, fig_path, "recharge_yearly.png", dpi)
    maybe_show(fig2, show)



# =============================================================================
# LAI-loop xarray-based hydro figures  (Fig 6 / 7 in test_withLAI_results.py)
# =============================================================================

def _xr_time_to_dates(da: xr.DataArray, time_dim: str) -> pd.DatetimeIndex:
    """Convert any xarray time coordinate to a pandas DatetimeIndex."""
    raw = da[time_dim].values
    # timedelta64 relative to first value → absolute if offset stored as attr
    if np.issubdtype(raw.dtype, np.timedelta64):
        raw = raw.astype("timedelta64[s]")
        # Fall back: treat t=0 as origin (caller should supply pev_series.index[0])
        origin = pd.Timestamp("1970-01-01")
        return origin + pd.to_timedelta(raw)
    return pd.DatetimeIndex(raw)


def _infer_time_dim(ds: xr.Dataset | xr.DataArray) -> str:
    """Return the name of the time dimension from a dataset or data-array."""
    dims = ds.dims if isinstance(ds, xr.DataArray) else list(ds.dims)
    for candidate in ("datetime", "time", "t"):
        if candidate in dims:
            return candidate
    return list(dims)[0]


def _node_id(nodes_dict: dict, key: str) -> int:
    """Extract a scalar node index from the nodes_dict regardless of nesting."""
    val = nodes_dict[key]
    # Support int, (int,), ((int,),), np.ndarray …
    while hasattr(val, "__iter__") and not isinstance(val, (str, np.integer)):
        val = next(iter(val))
    return int(val)


def plot_hydro_xr(
    psi_xr:     xr.Dataset,
    sw_xr:      xr.Dataset,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    pev_series:  pd.Series | None = None,
    tp_series:   pd.Series | None = None,
    date_start:  str | None = None,
    date_end:    str | None = None,
    show_rain:   bool = True,
    show_psi:    bool = True,
    show_sw:     bool = True,
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig6_hydro_timeseries.png",
) -> None:
    """Fig 6 — Rain/ETp + ψ + sw time-series built entirely from xarray Datasets.

    Parameters
    ----------
    psi_xr      : xr.Dataset with dimension ``(datetime|time, node)`` and at
                  least one variable containing pressure-head values (m).
    sw_xr       : xr.Dataset with the same structure for soil-water content (−).
    nodes_dict  : dict returned by ``resolve_nodes()``.
    fig_path    : output directory.
    pev_series  : optional ERA5 potential ET series (mm/day, DatetimeIndex).
    tp_series   : optional ERA5 precipitation series (mm/day, DatetimeIndex).
    date_start/end : ISO-date strings for the x-axis window.
    show_rain/psi/sw : toggle the three main panels.
    show_uh_s/uh_1/dh_s/dh_1 : toggle individual node series within psi/sw.
    colors_dict : override COLORS_HYDRO.
    scenario_id : label for the figure title.
    """
    colors_dict = colors_dict or COLORS_HYDRO

    # ── Identify psi variable and time dimension ──────────────────────────────
    psi_time_dim = _infer_time_dim(psi_xr)
    sw_time_dim  = _infer_time_dim(sw_xr)
    psi_var      = [v for v in psi_xr.data_vars][0]
    sw_var       = [v for v in sw_xr.data_vars][0]

    psi_dates = _xr_time_to_dates(psi_xr[psi_var], psi_time_dim)
    sw_dates  = _xr_time_to_dates(sw_xr[sw_var],   sw_time_dim)

    d0, d1 = date_window(date_start, date_end, psi_dates.to_series())
    print(f"  Plotting hydro_xr time-series ({d0.date()} → {d1.date()}) …")

    # ── Node indices ──────────────────────────────────────────────────────────
    uh_s_id = _node_id(nodes_dict, "mid_z-0m")
    uh_1_id = _node_id(nodes_dict, "mid_z-1m")
    dh_s_id = _node_id(nodes_dict, "outlet_z-0m")
    dh_1_id = _node_id(nodes_dict, "outlet_z-1m")

    psi_mask = (psi_dates >= d0) & (psi_dates <= d1)
    sw_mask  = (sw_dates  >= d0) & (sw_dates  <= d1)

    # ── psi arrays at reference nodes  (shape: T) ────────────────────────────
    psi_da = psi_xr[psi_var]   # (datetime, node)
    sw_da  = sw_xr[sw_var]

    def _psi_at(node_id: int) -> np.ndarray:
        return psi_da.isel({psi_time_dim: psi_mask}).sel(node=node_id).values

    def _sw_at(node_id: int) -> np.ndarray:
        return sw_da.isel({sw_time_dim: sw_mask}).sel(node=node_id).values

    # ── Count active panels ───────────────────────────────────────────────────
    n_panels = sum([show_rain, show_psi, show_sw])
    if n_panels == 0:
        print("  plot_hydro_xr: all panels disabled — nothing to plot.")
        return

    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(14, 3.2 * n_panels),
        sharex=True,
        constrained_layout=True,
    )
    if n_panels == 1:
        axes = [axes]
    ax_iter = iter(axes)

    # ── Panel 1: Rain / ETp ───────────────────────────────────────────────────
    if show_rain:
        ax = next(ax_iter)
        if pev_series is not None and tp_series is not None and len(pev_series):
            rain_m = (pev_series.index >= d0) & (pev_series.index <= d1)
            ax.bar(pev_series[rain_m].index, pev_series[rain_m].values,
                   width=1.0, color="grey",   label="ETp", alpha=0.5)
            ax.bar(tp_series[rain_m].index,  tp_series[rain_m].values,
                   width=1.0, color="skyblue", label="Rain")
            ax.invert_yaxis()
        else:
            ax.text(0.5, 0.5, "ERA5 forcing not available",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="#888888")
        ax.set_ylabel("mm/day")
        ax.set_title("Rain / ETp (ERA5)", fontsize=10, fontweight="bold")
        add_fire_vline(ax, d0, d1)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, linestyle="--", alpha=0.3)

    # ── Panel 2: ψ at reference nodes ─────────────────────────────────────────
    if show_psi:
        ax = next(ax_iter)
        psi_dates_w = psi_dates[psi_mask]
        if show_uh_s:
            ax.plot(psi_dates_w, _psi_at(uh_s_id),
                    label=f"Mid surf n={uh_s_id}",
                    color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
        if show_uh_1:
            ax.plot(psi_dates_w, _psi_at(uh_1_id),
                    label=f"Mid −1m n={uh_1_id}",
                    color=colors_dict["uphill_1m"], linestyle="--")
        if show_dh_s:
            ax.plot(psi_dates_w, _psi_at(dh_s_id),
                    label=f"Outlet surf n={dh_s_id}",
                    color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
        if show_dh_1:
            ax.plot(psi_dates_w, _psi_at(dh_1_id),
                    label=f"Outlet −1m n={dh_1_id}",
                    color=colors_dict["downhill_1m"], linestyle="--")
        ax.set_ylabel("ψ (m)")
        ax.set_title("Pressure Head ψ at reference nodes", fontsize=10, fontweight="bold")
        add_fire_vline(ax, d0, d1)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.axhline(0, color="black", linewidth=0.8, linestyle=":")

    # ── Panel 3: sw at reference nodes ────────────────────────────────────────
    if show_sw:
        ax = next(ax_iter)
        sw_dates_w = sw_dates[sw_mask]
        if show_uh_s:
            ax.plot(sw_dates_w, _sw_at(uh_s_id),
                    label=f"Mid surf n={uh_s_id}",
                    color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
        if show_uh_1:
            ax.plot(sw_dates_w, _sw_at(uh_1_id),
                    label=f"Mid −1m n={uh_1_id}",
                    color=colors_dict["uphill_1m"], linestyle="--")
        if show_dh_s:
            ax.plot(sw_dates_w, _sw_at(dh_s_id),
                    label=f"Outlet surf n={dh_s_id}",
                    color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
        if show_dh_1:
            ax.plot(sw_dates_w, _sw_at(dh_1_id),
                    label=f"Outlet −1m n={dh_1_id}",
                    color=colors_dict["downhill_1m"], linestyle="--")
        ax.set_ylabel("sw (−)")
        ax.set_title("Soil Water Content at reference nodes", fontsize=10, fontweight="bold")
        add_fire_vline(ax, d0, d1)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, linestyle="--", alpha=0.3)

    # ── Shared x-axis formatting ──────────────────────────────────────────────
    axes[-1].xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(mdates.AutoDateLocator())
    )
    axes[-1].set_xlabel("Date")
    fig.suptitle(
        f"Sc.{scenario_id} | Hydro time-series | {d0.date()} → {d1.date()}",
        fontsize=11, fontweight="bold",
    )
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_et_sw_xr(
    sw_xr:      xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    date_start:  str | None = None,
    date_end:    str | None = None,
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig7_et_sw_timeseries.png",
) -> None:
    """Fig 7 — Soil water content + spatial-mean actual ETa, from xarray Datasets.

    Parameters
    ----------
    sw_xr      : xr.Dataset with dimension ``(datetime|time, node)`` for sw.
    ET_xr_all  : xr.Dataset with spatial ET (X, Y, datetime|time), variable
                 name given by *et_var*.
    et_var     : name of the ET variable inside *ET_xr_all*.
    et_scale   : multiplicative factor (e.g. 1e3*86400 for m/s → mm/day).
    nodes_dict : dict returned by ``resolve_nodes()``.
    fig_path   : output directory.
    date_start/end : ISO-date strings for the x-axis window.
    show_uh_s/uh_1/dh_s/dh_1 : toggle individual sw node series.
    colors_dict : override COLORS_HYDRO.
    scenario_id : label for the figure title.
    """
    colors_dict = colors_dict or COLORS_HYDRO

    sw_time_dim = _infer_time_dim(sw_xr)
    et_time_dim = _infer_time_dim(ET_xr_all)
    sw_var      = list(sw_xr.data_vars)[0]

    sw_dates = _xr_time_to_dates(sw_xr[sw_var], sw_time_dim)
    et_dates = _xr_time_to_dates(ET_xr_all[et_var], et_time_dim)

    d0, d1 = date_window(date_start, date_end, sw_dates.to_series())
    print(f"  Plotting et_sw_xr time-series ({d0.date()} → {d1.date()}) …")

    sw_mask = (sw_dates >= d0) & (sw_dates <= d1)
    et_mask = (et_dates >= d0) & (et_dates <= d1)

    sw_da = sw_xr[sw_var]

    # Spatial-mean ET (all non-time dims averaged)
    spatial_dims = [d for d in ET_xr_all[et_var].dims
                    if d.lower() not in ("datetime", "time", "t")]
    et_mean = ET_xr_all[et_var].mean(dim=spatial_dims) * et_scale  # (time,)
    et_vals = et_mean.values[et_mask]
    et_dates_w = et_dates[et_mask]

    # Node IDs
    uh_s_id = _node_id(nodes_dict, "mid_z-0m")
    uh_1_id = _node_id(nodes_dict, "mid_z-1m")
    dh_s_id = _node_id(nodes_dict, "outlet_z-0m")
    dh_1_id = _node_id(nodes_dict, "outlet_z-1m")

    sw_dates_w = sw_dates[sw_mask]

    def _sw_at(node_id: int) -> np.ndarray:
        return sw_da.isel({sw_time_dim: sw_mask}).sel(node=node_id).values

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(14, 6), sharex=True, constrained_layout=True
    )

    # ── Top: sw at reference nodes ────────────────────────────────────────────
    if show_uh_s:
        ax0.plot(sw_dates_w, _sw_at(uh_s_id),
                 label=f"Mid surf n={uh_s_id}",
                 color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
    if show_uh_1:
        ax0.plot(sw_dates_w, _sw_at(uh_1_id),
                 label=f"Mid −1m n={uh_1_id}",
                 color=colors_dict["uphill_1m"], linestyle="--")
    if show_dh_s:
        ax0.plot(sw_dates_w, _sw_at(dh_s_id),
                 label=f"Outlet surf n={dh_s_id}",
                 color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
    if show_dh_1:
        ax0.plot(sw_dates_w, _sw_at(dh_1_id),
                 label=f"Outlet −1m n={dh_1_id}",
                 color=colors_dict["downhill_1m"], linestyle="--")

    ax0.set_ylabel("sw (−)")
    ax0.set_title("Soil Water Content at reference nodes", fontsize=10, fontweight="bold")
    add_fire_vline(ax0, d0, d1)
    ax0.legend(loc="upper right", fontsize=8)
    ax0.grid(True, linestyle="--", alpha=0.3)

    # ── Bottom: spatial-mean ETa ──────────────────────────────────────────────
    if et_var in ET_xr_all and len(et_vals):
        ax1.scatter(et_dates_w, et_vals,
                    c=colors_dict["ET_mean"], s=14, label="ETa spatial mean")
    else:
        ax1.text(0.5, 0.5, f"'{et_var}' not found in ET dataset",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=9, color="#888888")

    ax1.set_ylabel("ETa (mm/day)")
    ax1.set_xlabel("Date")
    ax1.set_title(
        f"Actual Evapotranspiration — spatial mean ({et_var})",
        fontsize=10, fontweight="bold",
    )
    add_fire_vline(ax1, d0, d1)
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(True, linestyle="--", alpha=0.3)
    ax1.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(mdates.AutoDateLocator())
    )

    fig.suptitle(
        f"Sc.{scenario_id} | sw & ETa | {d0.date()} → {d1.date()}",
        fontsize=11, fontweight="bold",
    )
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)
    
    

# =============================================================================
# Colour / label constants
# =============================================================================

CMAP_ET    = "YlGnBu"
MONTH_ABBR = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


# =============================================================================
# Internal helpers
# =============================================================================

def _infer_time_dim(ds: xr.Dataset | xr.DataArray) -> str:
    """Return the time dimension name, regardless of capitalisation."""
    dims = list(ds.dims)
    for candidate in ("time", "datetime", "t", "date", "TIME", "Date"):
        if candidate in dims:
            return candidate
    for d in dims:
        if d in ds.coords and np.issubdtype(ds[d].dtype, np.datetime64):
            return d
    raise ValueError(f"Cannot detect time dimension in dims={dims}")


def _spatial_dims(da: xr.DataArray, time_dim: str) -> list[str]:
    """Return all non-time dimension names."""
    return [d for d in da.dims if d != time_dim]


def _da(ds: xr.Dataset, et_var: str) -> xr.DataArray:
    if et_var not in ds:
        raise KeyError(f"Variable '{et_var}' not found. Available: {list(ds.data_vars)}")
    return ds[et_var]


# =============================================================================
# Fig A — Monthly spatial-mean ETa, one line per year
# =============================================================================

def plot_eta_monthly_by_year(
    ds_et:        xr.Dataset,
    et_var:       str,
    fig_path:     Path,
    *,
    et_scale:     float = 1.0,
    dpi:          int   = 150,
    show:         bool  = False,
    fname:        str   = "ETa_monthly_mean_by_year.png",
    product_name: str   = "ETa",
) -> None:
    """
    Spatial-mean ETa per month, one coloured line per year.

    Parameters
    ----------
    ds_et    : xr.Dataset containing *et_var*
    et_var   : name of the ET variable (e.g. "ET-gf")
    fig_path : output directory
    et_scale : multiplicative scale factor (default 1.0 — assumes mm/day already)
    dpi      : figure resolution
    show     : display interactively before closing
    fname    : output filename
    """
    da       = _da(ds_et, et_var)
    time_dim = _infer_time_dim(ds_et)
    sp_dims  = _spatial_dims(da, time_dim)

    et_mean = da.mean(dim=sp_dims, skipna=True) * et_scale
    times   = pd.DatetimeIndex(et_mean[time_dim].values)
    years   = sorted(set(times.year))
    months  = np.arange(1, 13)

    cmap_lines = plt.get_cmap("tab10", len(years))

    fig, ax = plt.subplots(figsize=(10, 5))

    for i, year in enumerate(years):
        mask    = times.year == year
        ts_year = et_mean.isel({time_dim: mask})
        mo_mean = ts_year.groupby(f"{time_dim}.month").mean()
        ax.plot(
            mo_mean.month.values,
            mo_mean.values,
            marker="o", linewidth=2,
            color=cmap_lines(i), label=str(year),
        )

    ax.set_title(f"Monthly Spatial-Mean {product_name} — by year  [{et_var}]",
                 fontsize=13, fontweight="bold")
    ax.set_xlabel("Month")
    ax.set_ylabel(f"{product_name} (mm/day)" + (f"  ×{et_scale}" if et_scale != 1 else ""))
    ax.set_xticks(months)
    ax.set_xticklabels(MONTH_ABBR)
    ax.legend(title="Year", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()

    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


# =============================================================================
# Fig B — Monthly spatial maps for a chosen year  (mosaic, 4-column layout)
# =============================================================================

def plot_eta_spatial_maps_year(
    ds_et:         xr.Dataset,
    et_var:        str,
    selected_year: int,
    fig_path:      Path,
    *,
    et_scale:     float       = 1.0,
    cmap:         str         = CMAP_ET,
    vmin:         float | None = None,
    vmax:         float | None = None,
    dpi:          int         = 150,
    show:         bool        = False,
    product_name: str         = "ETa",
) -> None:
    """
    12-panel mosaic of monthly-mean ETa spatial maps for *selected_year*.

    Parameters
    ----------
    ds_et         : xr.Dataset containing *et_var*
    et_var        : name of the ET variable
    selected_year : calendar year to plot (e.g. 2020)
    fig_path      : output directory
    et_scale      : multiplicative scale factor
    cmap          : matplotlib colourmap
    vmin / vmax   : colourmap limits; auto = 2nd/98th percentile of the year
    dpi           : figure resolution
    show          : display interactively before closing
    """
    da       = _da(ds_et, et_var) * et_scale
    time_dim = _infer_time_dim(ds_et)

    da_year = da.sel({time_dim: da[time_dim].dt.year == selected_year})
    if da_year.sizes[time_dim] == 0:
        print(f"  plot_eta_spatial_maps_year: no data for {selected_year} — skipped.")
        return

    da_monthly = da_year.resample({time_dim: "1MS"}).mean()
    times      = da_monthly[time_dim].values

    if len(times) == 0:
        print(f"  plot_eta_spatial_maps_year: no monthly slices for {selected_year}.")
        return

    if vmin is None:
        vmin = float(np.nanpercentile(da_year.values, 2))
    if vmax is None:
        vmax = float(np.nanpercentile(da_year.values, 98))

    ncols  = 4
    nrows  = (len(times) + ncols - 1) // ncols
    keys   = [pd.Timestamp(t).strftime("%Y-%m") for t in times]
    mosaic = [
        [keys[r * ncols + c] if r * ncols + c < len(keys) else "."
         for c in range(ncols)]
        for r in range(nrows)
    ]

    fig, axes = plt.subplot_mosaic(mosaic, figsize=(5 * ncols, 4 * nrows),
                                   empty_sentinel=".")

    sp_orig      = _spatial_dims(da, time_dim)
    y_dim, x_dim = sp_orig[0], sp_orig[1]
    extent = [
        float(da[x_dim].min()), float(da[x_dim].max()),
        float(da[y_dim].min()), float(da[y_dim].max()),
    ]

    for t, key in zip(times, keys):
        ax    = axes[key]
        frame = da_monthly.sel({time_dim: t})

        im = ax.imshow(frame.values, origin="upper", extent=extent,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       aspect="equal", interpolation="nearest")

        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="4%", pad=0.05)
        cb      = plt.colorbar(im, cax=cax)
        cb.set_label("mm/day", fontsize=7)
        cb.ax.tick_params(labelsize=7)

        ax.set_title(pd.Timestamp(t).strftime("%B %Y"), fontsize=10)
        ax.set_xlabel(x_dim, fontsize=7)
        ax.set_ylabel(y_dim, fontsize=7)
        ax.tick_params(labelsize=7)

    plt.suptitle(f"Monthly {product_name} spatial maps — {selected_year}  [{et_var}]",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()

    save_fig(fig, fig_path, f"{product_name}_spatial_maps_{selected_year}.png", dpi)
    maybe_show(fig, show)


# =============================================================================
# Fig C — Histogram of per-pixel annual ETa sums
# =============================================================================

def plot_eta_annual_hist(
    ds_et:        xr.Dataset,
    et_var:       str,
    fig_path:     Path,
    *,
    et_scale:     float = 1.0,
    bins:         int   = 40,
    dpi:          int   = 150,
    show:         bool  = False,
    fname:        str   = "ETa_annual_sum_histogram.png",
    product_name: str   = "ETa",
) -> None:
    """
    Histogram of per-pixel annual ETa sums.

    For every calendar year in the dataset, each spatial pixel is summed
    over time to obtain its annual ETa total.  All per-pixel annual values
    are collected into a single histogram so the spatial distribution of
    annual ETa is visible.  One coloured histogram series is drawn per year.

    Parameters
    ----------
    ds_et    : xr.Dataset containing *et_var*
    et_var   : name of the ET variable (e.g. "ET-gf")
    fig_path : output directory
    et_scale : multiplicative scale factor (default 1.0)
    bins     : number of histogram bins
    dpi      : figure resolution
    show     : display interactively before closing
    fname    : output filename
    """
    da       = _da(ds_et, et_var) * et_scale
    time_dim = _infer_time_dim(ds_et)
    times    = pd.DatetimeIndex(da[time_dim].values)
    years    = sorted(set(times.year))

    fig, ax = plt.subplots(figsize=(9, 5))
    cmap_lines = plt.get_cmap("tab10", len(years))

    for i, year in enumerate(years):
        mask    = times.year == year
        da_year = da.isel({time_dim: mask})
        # sum over time → per-pixel annual total; flatten to 1-D
        annual  = da_year.sum(dim=time_dim, skipna=True).values.ravel()
        annual  = annual[~np.isnan(annual)]
        if annual.size == 0:
            continue
        ax.hist(
            annual, bins=bins, alpha=0.55,
            color=cmap_lines(i), edgecolor="none",
            label=str(year), density=False,
        )

    ax.set_xlabel(f"Annual {product_name} sum (mm/year)", fontsize=11)
    ax.set_ylabel("Pixel count", fontsize=11)
    ax.set_title(
        f"Per-pixel annual {product_name} distribution  [{et_var}]",
        fontsize=13, fontweight="bold",
    )
    ax.legend(title="Year", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()

    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


# =============================================================================
# Private core  (Fig D helper)
# =============================================================================
 
def _plot_eta_total_freq(
    ds_et:        xr.Dataset,
    freq:         str,          # "ME" or "YE"
    fig_path:     Path,
    *,
    label:        str  = "",
    dpi:          int  = 150,
    show:         bool = False,
    product_name: str  = "ETa",
) -> None:
    """
    Bar chart of domain-total ETa aggregated at *freq* frequency.
 
    Bar colour encodes ETa intensity (value-normalised colourmap YlGnBu),
    so drier periods are yellow and wetter periods are blue-green.
 
    Parameters
    ----------
    ds_et    : clipped xr.Dataset with one ET variable
    freq     : "ME" (monthly) or "YE" (yearly)
    fig_path : output directory
    label    : plot-selection label appended to filename (e.g. "p00")
    dpi      : figure resolution
    show     : display interactively before closing
    """
    # ── 1. Variable & spatial mean → DataFrame ────────────────────────────────
    et_var  = list(ds_et.data_vars)[0]
    da      = ds_et[et_var]
    sp_dims = [d for d in da.dims if d != "time"]
    et_mean = da.mean(dim=sp_dims, skipna=True)
 
    df = pd.DataFrame(
        {f"{product_name}_mm": et_mean.values},
        index=pd.DatetimeIndex(da["time"].values),
    ).resample(freq).sum()
 
    # ── 2. Labels ─────────────────────────────────────────────────────────────
    p_str = "monthly" if freq == "ME" else "annual"
    time_labels = (
        df.index.strftime("%Y-%m") if freq == "ME"
        else df.index.strftime("%Y")
    )
 
    # ── 3. Colour by ETa intensity (value-normalised) ─────────────────────────
    values   = df[f"{product_name}_mm"].values
    norm     = mcolors.Normalize(vmin=np.nanmin(values), vmax=np.nanmax(values))
    cmap     = cm.get_cmap("YlGnBu")
    colors   = [cmap(norm(v)) for v in values]
    mean_val = float(np.nanmean(values))
 
    # ── 4. Figure ─────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(10, len(values) * 0.65), 5))
 
    bars = ax.bar(time_labels, values, color=colors,
                  edgecolor="white", linewidth=0.5, zorder=2)
 
    # colourbar (ETa intensity scale)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.01, shrink=0.85)
    cb.set_label(f"{product_name} intensity (mm)", fontsize=9)
 
    # mean line
    ax.axhline(mean_val, color="#C44E52", linewidth=1.5, linestyle="--",
               label=f"mean  {mean_val:,.1f} mm", zorder=3)
 
    # value labels on bars
    for bar, v in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(values) * 0.01,
            f"{v:,.1f}", ha="center", va="bottom", fontsize=7, color="#222222",
        )
 
    ax.set_ylabel(f"{product_name} accumulated (mm)", fontsize=10)
    ax.set_title(
        f"Domain {product_name} — {p_str.capitalize()}  [{et_var}]",
        fontsize=12, fontweight="bold",
    )
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", linestyle="--", alpha=0.3, zorder=1)
    ax.legend(fontsize=9)
    ax.set_ylim(0, max(values) * 1.14)
 
    plt.tight_layout()
 
    suffix = f"_{label}" if label else ""
    out = fig_path / f"{product_name}_total_{p_str}{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
 
    if show:
        plt.show()
    plt.close(fig)
 
 
# =============================================================================
# Public wrapper
# =============================================================================
 
def plot_eta_total(
    ds_et:        xr.Dataset,
    fig_path:     Path,
    *,
    label:        str  = "",
    dpi:          int  = 150,
    show:         bool = False,
    product_name: str  = "ETa",
) -> None:
    """
    Monthly and annual domain ETa bar charts.
 
    Produces two figures:
        ETa_total_monthly_<label>.png
        ETa_total_annual_<label>.png
 
    Parameters
    ----------
    ds_et    : clipped xr.Dataset (one ET variable, 'time' dimension)
    fig_path : output directory
    label    : plot-selection label (e.g. "p00", "all")
    dpi      : figure resolution
    show     : display interactively before closing
    """
    print(f"  Computing monthly domain {product_name} …")
    _plot_eta_total_freq(ds_et, "ME", fig_path, label=label, dpi=dpi, show=show, product_name=product_name)
 
    print(f"  Computing annual domain {product_name} …")
    _plot_eta_total_freq(ds_et, "YE", fig_path, label=label, dpi=dpi, show=show, product_name=product_name)