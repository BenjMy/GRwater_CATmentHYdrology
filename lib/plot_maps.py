"""
plot_maps.py - Vegetation / LAI / ET maps, histograms, annual statistics and GIFs.
================================================================================
Split unchanged out of results_plotter.py (REFACTORING.md 2.11); results_plotter re-exports every name.
"""
from __future__ import annotations

from typing import Any
from pathlib import Path
from collections import defaultdict
import tempfile

import matplotlib.cm as cm
import imageio
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from pyCATHY import meshtools as mt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr

from plot_common import CMAP_ET, CMAP_LAI, LAI_UNITS, LAI_VMAX, _da, _infer_time_dim, _spatial_dims, add_fire_vline, maybe_show, save_fig, time_label, veg_colorbar


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
                  vmin=0.5, vmax=3.5, aspect="equal")
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

    n     = len(monthly_dates)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols

    _h, _w = np.asarray(veg_map_history[monthly_dates[0]]).shape[:2]
    _map_aspect = float(np.clip(_h / max(_w, 1), 0.4, 2.5))   # height / width

    # constrained_layout + one shared discrete colourbar (a colourbar per
    # panel used to shrink every map and leave uneven gaps).
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.4 * ncols + 1.4, (3.4 * _map_aspect + 0.30) * nrows + 0.6),
        squeeze=False,
        constrained_layout=True,
    )
    fig.get_layout_engine().set(w_pad=0.04, h_pad=0.05, wspace=0.03, hspace=0.06)

    used = []
    for k, d in enumerate(monthly_dates):
        ax = axes[k // ncols][k % ncols]
        used.append(ax)
        ax.imshow(veg_map_history[d], origin="upper", cmap=cmap_veg,
                  vmin=0.5, vmax=3.5, aspect="equal")
        ax.set_title(f"{d:%B}", fontsize=10, pad=3)
        ax.tick_params(labelsize=7, length=2.5, pad=1.5,
                       labelleft=(k % ncols == 0),
                       labelbottom=(k + ncols >= n))
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].set_visible(False)

    bounds = [0.5, 1.5, 2.5, 3.5]
    sm = cm.ScalarMappable(cmap=cmap_veg, norm=mcolors.BoundaryNorm(bounds, cmap_veg.N))
    cb = fig.colorbar(sm, ax=used, ticks=[1, 2, 3], shrink=0.7, aspect=25, pad=0.015)
    cb.ax.set_yticklabels(("Bare", "Moderate veg.", "Dense veg."))
    fig.supxlabel("Column", fontsize=9)
    fig.supylabel("Row", fontsize=9)
    fig.suptitle(f"Vegetation zone maps — {selected_year} (monthly)",
                 fontsize=13, fontweight="bold")
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
                       vmin=vmin, vmax=vmax, aspect="equal")
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
                          vmin=vmin, vmax=vmax, aspect="equal")
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
                  vmin=0.5, vmax=3.5, aspect="equal")
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
