"""
plot_eta.py - ETa figures: monthly by year, spatial maps, annual histograms, totals, ETa/ETp ratios.
==================================================================================================
Split unchanged out of plot_results_shared.py (REFACTORING.md 2.11); plot_results_shared re-exports every name.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr

from scipy.odr import ODR, Model, RealData   # Fig. 8 total-least-squares trend
from plot_common import unwrap_node_id
from plot_common import CMAP_ET, COLORS_HYDRO, MONTH_ABBR, _da, _fmt_node_id, _infer_time_dim, _node_ids, _spatial_dims, _values_at_nodes, _xr_time_to_dates, add_fire_vline, date_window, maybe_show, save_fig


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


def plot_lai_et_sw_psi_ratio_xr(
    lai_ds:     xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    ETp_xr:     xr.Dataset,
    etp_var:    str,
    etp_scale:  float,
    sw_xr:      xr.Dataset,
    psi_xr:     xr.Dataset,
    nodes_dict: dict,
    grid3d:     dict,
    fig_path:   Path,
    *,
    lai_var:     str  = "LAI",
    date_start:  str | None = None,
    date_end:    str | None = None,
    show_uh:     bool = True,
    show_dh:     bool = True,
    show_uh_1:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig8_lai_et_sw_psi_ratio_timeseries.png",
    et_unit:     str  = "mm/day",
) -> None:
    """Fig 8 -- LAI, ETa, sw, psi and the ETa/ETp ratio, stacked on one
    shared time axis, all read off the same reference nodes used by
    Fig 6/7 (``resolve_nodes()`` -> ``nodes_dict``).

    Five panels, top to bottom: LAI, ETa, sw, psi, ETa/ETp.

    Parameters
    ----------
    lai_ds   : xr.Dataset with a time dimension and variable *lai_var*,
               carrying real ``x``/``y`` coordinates (e.g. ``ds_LAI_raster``
               before it is stacked into ``lai_history``, or the
               ``lai_history_to_dataset()`` output with DEM x/y coordinates
               assigned onto its row/col dims). LAI lives on the satellite
               raster grid, not on the CATHY mesh, so reference-node LAI is
               a nearest-pixel lookup at each node's (x, y) rather than a
               ``.sel(node=...)`` like the other panels.
    ET_xr_all, et_var, et_scale : actual-ET dataset, same (time, node)
               convention as *sw_xr* / *psi_xr* (i.e. already CATHY
               surface-node output, not the gridded satellite product).
    ETp_xr, etp_var, etp_scale : potential-ET dataset mapped onto the same
               mesh nodes as *ET_xr_all* (e.g. via
               ``pyCATHY.meshtools.map_grid_to_mesh(ds_etp, ds_mesh, ...)``
               in ``Agramon_withLAI_withETp.py``) -- needed so the ETa/ETp
               ratio panel is computed node-by-node on a common time axis.
    sw_xr, psi_xr : same (time, node) Datasets as in ``plot_hydro_xr``.
    nodes_dict : dict from ``resolve_nodes()`` (mid_z-*m / outlet_z-*m).
    grid3d     : dict from ``simu.read_outputs("grid3d")`` -- used to look
               up each reference node's (x, y) for the LAI pixel lookup.
    fig_path   : output directory.
    show_uh/show_dh : toggle the mid-uphill / outlet reference nodes.
    show_uh_1/show_dh_1 : toggle the -1 m depth variants (sw & psi panels
               only -- LAI, ETa and the ratio have no vertical structure,
               so only the surface position of each site is plotted there).
    colors_dict : override COLORS_HYDRO.
    scenario_id : label for the figure title.
    et_unit    : y-axis unit label for the ETa panel.

    Notes
    -----
    - LAI and the ETa/ETp ratio are only meaningful at the *surface*
      reference nodes (mid_z-0m / outlet_z-0m): depth variants share the
      same (x, y), so plotting them again would just duplicate the line.
    - If *ET_xr_all* and *ETp_xr* turn out not to share the same time
      axis (different native time steps), the ratio panel shows a
      placeholder instead of silently misaligning the two series --
      resample one onto the other before calling this function.
    """
    colors_dict = colors_dict or COLORS_HYDRO

    # ── Time dims / dates for every input ─────────────────────────────────────
    psi_time_dim = _infer_time_dim(psi_xr)
    sw_time_dim  = _infer_time_dim(sw_xr)
    et_time_dim  = _infer_time_dim(ET_xr_all)
    etp_time_dim = _infer_time_dim(ETp_xr)
    lai_time_dim = _infer_time_dim(lai_ds)

    psi_var = list(psi_xr.data_vars)[0]
    sw_var  = list(sw_xr.data_vars)[0]

    psi_dates = _xr_time_to_dates(psi_xr[psi_var],   psi_time_dim)
    sw_dates  = _xr_time_to_dates(sw_xr[sw_var],     sw_time_dim)
    et_dates  = _xr_time_to_dates(ET_xr_all[et_var], et_time_dim)
    etp_dates = _xr_time_to_dates(ETp_xr[etp_var],   etp_time_dim)
    lai_dates = _xr_time_to_dates(lai_ds[lai_var],   lai_time_dim)

    d0, d1 = date_window(date_start, date_end, sw_dates.to_series())
    print(f"  Plotting lai_et_sw_psi_ratio_xr time-series ({d0.date()} -> {d1.date()}) ...")

    psi_mask = (psi_dates >= d0) & (psi_dates <= d1)
    sw_mask  = (sw_dates  >= d0) & (sw_dates  <= d1)
    et_mask  = (et_dates  >= d0) & (et_dates  <= d1)
    etp_mask = (etp_dates >= d0) & (etp_dates <= d1)
    lai_mask = (lai_dates >= d0) & (lai_dates <= d1)

    psi_da = psi_xr[psi_var]
    sw_da  = sw_xr[sw_var]
    et_da  = ET_xr_all[et_var] * et_scale
    etp_da = ETp_xr[etp_var]   * etp_scale

    # ── Reference node ids (band-preserving, see _node_ids) ───────────────────
    uh_s_id = _node_ids(nodes_dict, "mid_z-0m")
    uh_1_id = _node_ids(nodes_dict, "mid_z-1m")
    dh_s_id = _node_ids(nodes_dict, "outlet_z-0m")
    dh_1_id = _node_ids(nodes_dict, "outlet_z-1m")

    def _psi_at(ids):
        return _values_at_nodes(psi_da, psi_time_dim, psi_mask, ids)

    def _sw_at(ids):
        return _values_at_nodes(sw_da, sw_time_dim, sw_mask, ids)

    def _et_at(ids):
        """Nearest-pixel ETa series at the (x, y) of a reference node/band.
        ET_xr_all lives on the (X, Y) raster grid (et_output.nc's own
        convention — see the module docstring and plot_et_sw_xr, which
        takes a spatial mean rather than a node .sel() for the same
        reason), not on CATHY mesh nodes. This mirrors _lai_at below
        rather than routing through _values_at_nodes()."""
        x, y = _xy_of(ids)
        da_pt = et_da.sel(X=x, Y=y, method="nearest")
        return da_pt.isel({et_time_dim: et_mask}).values

    def _etp_at(ids):
        return _values_at_nodes(etp_da, etp_time_dim, etp_mask, ids)

    # ── (x, y) of a reference node/band, for the LAI pixel lookup ─────────────
    mesh_nodes = grid3d["mesh3d_nodes"]

    def _xy_of(ids) -> tuple[float, float]:
        """Mean (x, y) of a scalar node id or a band of ids. Depth doesn't
        change (x, y), so this is exact for a single id and a sensible
        summary for a band."""
        if isinstance(ids, (list, tuple, np.ndarray)):
            idx = np.asarray(ids, dtype=int)
            return float(mesh_nodes[idx, 0].mean()), float(mesh_nodes[idx, 1].mean())
        x, y, _z = mesh_nodes[int(ids)]
        return float(x), float(y)

    def _lai_at(ids) -> np.ndarray:
        """Nearest-pixel LAI series at the (x, y) of a reference node/band."""
        x, y = _xy_of(ids)
        da_pt = lai_ds[lai_var].sel(x=x, y=y, method="nearest")
        return da_pt.isel({lai_time_dim: lai_mask}).values

    lai_dates_w = lai_dates[lai_mask]
    psi_dates_w = psi_dates[psi_mask]
    sw_dates_w  = sw_dates[sw_mask]
    et_dates_w  = et_dates[et_mask]
    etp_dates_w = etp_dates[etp_mask]

    # ── Figure: LAI / ETa / sw / psi / ETa÷ETp ────────────────────────────────
    fig, (ax_lai, ax_et, ax_sw, ax_psi, ax_ratio) = plt.subplots(
        5, 1, figsize=(14, 13.5), sharex=True, constrained_layout=True
    )

    def _panel(ax, get_fn, dates_w, ylabel, title, *, with_depth=True):
        if show_uh:
            ax.plot(dates_w, get_fn(uh_s_id),
                    label=f"Mid surf n={_fmt_node_id(uh_s_id)}",
                    color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
        if with_depth and show_uh_1:
            ax.plot(dates_w, get_fn(uh_1_id),
                    label=f"Mid -1m n={_fmt_node_id(uh_1_id)}",
                    color=colors_dict["uphill_1m"], linestyle="--")
        if show_dh:
            ax.plot(dates_w, get_fn(dh_s_id),
                    label=f"Outlet surf n={_fmt_node_id(dh_s_id)}",
                    color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
        if with_depth and show_dh_1:
            ax.plot(dates_w, get_fn(dh_1_id),
                    label=f"Outlet -1m n={_fmt_node_id(dh_1_id)}",
                    color=colors_dict["downhill_1m"], linestyle="--")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10, fontweight="bold")
        add_fire_vline(ax, d0, d1)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, linestyle="--", alpha=0.3)

    # Panel 1 — LAI (surface positions only — no vertical structure)
    _panel(ax_lai, _lai_at, lai_dates_w, "LAI (m² m⁻²)",
           "Leaf Area Index at reference nodes", with_depth=False)

    # Panel 2 — ETa
    _panel(ax_et, _et_at, et_dates_w, f"ETa ({et_unit})",
           f"Actual Evapotranspiration at reference nodes ({et_var})", with_depth=False)

    # Panel 3 — sw
    _panel(ax_sw, _sw_at, sw_dates_w, "sw (−)",
           "Soil Water Content at reference nodes")

    # Panel 4 — psi
    _panel(ax_psi, _psi_at, psi_dates_w, "ψ (m)",
           "Pressure Head ψ at reference nodes")
    ax_psi.axhline(0, color="black", linewidth=0.8, linestyle=":")

    # Panel 5 — ETa / ETp ratio (surface only)
    # ETa (CATHY model output, et_dates_w) and ETp (satellite ET0 mapped
    # onto mesh nodes, etp_dates_w) don't share native time steps or
    # coverage, so resample ETa onto ETp's time axis (nearest within 1
    # day) before dividing, rather than requiring the two masks to
    # already line up.
    def _ratio_at(ids):
        eta_vals = np.asarray(_et_at(ids), dtype=float)
        etp_vals = np.asarray(_etp_at(ids), dtype=float)
        eta_series = (
            pd.Series(eta_vals, index=pd.DatetimeIndex(et_dates_w))
              .groupby(level=0).mean()
              .sort_index()
        )
        eta_on_etp = eta_series.reindex(
            pd.DatetimeIndex(etp_dates_w),
            method="nearest",
            tolerance=pd.Timedelta("1D"),
        ).to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(etp_vals > 0, eta_on_etp / etp_vals, np.nan)

    _panel(ax_ratio, _ratio_at, etp_dates_w, "ETa / ETp (−)",
           "ETa \u00f7 ETp ratio at reference nodes", with_depth=False)
    ax_ratio.axhline(1.0, color="black", linewidth=0.8, linestyle=":")
    ylim_top = ax_ratio.get_ylim()[1]
    ax_ratio.set_ylim(0, max(1.2, ylim_top))

    ax_ratio.set_xlabel("Date")
    ax_ratio.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(mdates.AutoDateLocator())
    )

    fig.suptitle(
        f"Sc.{scenario_id} | LAI, ETa, sw, \u03c8 & ETa\u00f7ETp | {d0.date()} \u2192 {d1.date()}",
        fontsize=11, fontweight="bold",
    )
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


# =============================================================================
# Catchment-scale ETa/ETp ratio (domain mean & median)
# =============================================================================
def plot_catchment_eta_etp_ratio(
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    etp_mean_series:   pd.Series,
    etp_median_series: pd.Series,
    fig_path:   Path,
    *,
    date_start:  str | None = None,
    date_end:    str | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig_catchment_eta_etp_ratio.png",
    freq:        str | None = None,
) -> None:
    """Catchment-scale ETa/ETp ratio time series — domain **mean** and
    domain **median**, both computed over the whole spatial extent,
    unlike ``plot_lai_et_sw_psi_ratio_xr``'s ratio panel which is per
    reference node.

    Parameters
    ----------
    ET_xr_all, et_var, et_scale : the same raw spatial ET dataset used
        everywhere else (e.g. Fig 7/7h) — dims (X, Y, datetime|time).
        Domain mean/median are taken directly across the (X, Y) extent,
        the same convention ``plot_et_sw_xr`` uses for its spatial-mean
        ETa line.
    etp_mean_series, etp_median_series : pd.Series (DatetimeIndex) —
        domain-mean / domain-median satellite ETp already reduced to the
        catchment polygon (``_domain_mean_daily_series(..., stat=...)``
        in plot_run_results.py).
    fig_path : output directory.
    date_start/end : ISO-date strings for the x-axis window; defaults to
        the mean-ETp series' own coverage.
    scenario_id : label for the figure title.

    Notes
    -----
    ETa (CATHY model output) and ETp (satellite-derived) don't share
    native time steps, so each ETa series (mean, median) is resampled
    onto its matching ETp series' own dates (nearest within 1 day)
    before the ratio is taken — the same approach used for the
    node-level ratio panel in ``plot_lai_et_sw_psi_ratio_xr``.
    """
    et_time_dim = _infer_time_dim(ET_xr_all[et_var])
    et_dates    = _xr_time_to_dates(ET_xr_all[et_var], et_time_dim)

    spatial_dims = [d for d in ET_xr_all[et_var].dims
                    if d.lower() not in ("datetime", "time", "t")]
    eta_mean_da   = ET_xr_all[et_var].mean(dim=spatial_dims)   * et_scale
    eta_median_da = ET_xr_all[et_var].median(dim=spatial_dims) * et_scale

    d0, d1 = date_window(date_start, date_end, etp_mean_series.index)
    print(f"  Plotting catchment_eta_etp_ratio time-series "
          f"({d0.date()} → {d1.date()}) …")

    et_mask         = (et_dates >= d0) & (et_dates <= d1)
    etp_mean_mask   = (etp_mean_series.index   >= d0) & (etp_mean_series.index   <= d1)
    etp_median_mask = (etp_median_series.index >= d0) & (etp_median_series.index <= d1)

    def _to_clean_series(vals: np.ndarray, dates: pd.DatetimeIndex) -> pd.Series:
        """DataArray values → sorted, de-duplicated pd.Series, ready for
        a nearest-tolerance reindex (duplicate timestamps averaged)."""
        return (
            pd.Series(np.asarray(vals, dtype=float), index=pd.DatetimeIndex(dates))
              .groupby(level=0).mean()
              .sort_index()
        )

    eta_mean_series   = _to_clean_series(eta_mean_da.values[et_mask],   et_dates[et_mask])
    eta_median_series = _to_clean_series(eta_median_da.values[et_mask], et_dates[et_mask])

    etp_mean_w   = etp_mean_series.loc[etp_mean_mask]
    etp_median_w = etp_median_series.loc[etp_median_mask]

    eta_mean_on_etp = eta_mean_series.reindex(
        etp_mean_w.index, method="nearest", tolerance=pd.Timedelta("1D")
    )
    eta_median_on_etp = eta_median_series.reindex(
        etp_median_w.index, method="nearest", tolerance=pd.Timedelta("1D")
    )

    # Optional temporal aggregation (e.g. freq="MS" -> monthly). ETa and ETp
    # are first paired day by day (so both sums cover exactly the same days),
    # then each is summed per period and the ratio of the sums is taken —
    # the correct way to aggregate a flux ratio (not a mean of daily ratios).
    def _agg(eta_on_etp: pd.Series, etp_w: pd.Series):
        if not freq:
            return etp_w.index, eta_on_etp.values, etp_w.values
        d = pd.concat([eta_on_etp.rename("eta"), etp_w.rename("etp")],
                      axis=1).dropna()
        m = d.resample(freq).sum()
        return m.index, m["eta"].to_numpy(), m["etp"].to_numpy()

    idx_mean,   eta_mean_v,   etp_mean_v   = _agg(eta_mean_on_etp,   etp_mean_w)
    idx_median, eta_median_v, etp_median_v = _agg(eta_median_on_etp, etp_median_w)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_mean = np.where(etp_mean_v > 0, eta_mean_v / etp_mean_v, np.nan)
        ratio_median = np.where(etp_median_v > 0, eta_median_v / etp_median_v, np.nan)

    fig, ax = plt.subplots(figsize=(14, 4.5), constrained_layout=True)
    _mk = "o" if freq else None
    ax.plot(idx_mean, ratio_mean,
            label="Mean (domain-mean ETa ÷ domain-mean ETp)",
            color="#1b7837", linewidth=1.2, marker=_mk, ms=3)
    ax.plot(idx_median, ratio_median,
            label="Median (domain-median ETa ÷ domain-median ETp)",
            color="#762a83", linewidth=1.2, linestyle="--", marker=_mk, ms=3)
    ax.axhline(1.0, color="black", linewidth=0.8, linestyle=":")
    ax.set_ylabel("ETa / ETp (−)")
    ax.set_xlabel("Date")
    ax.set_title(
        f"Sc.{scenario_id} | Catchment-scale ETa ÷ ETp"
        f"{' (monthly)' if freq in ('MS', 'ME', 'M') else ''} | {d0.date()} → {d1.date()}",
        fontsize=11, fontweight="bold",
    )
    add_fire_vline(ax, d0, d1)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(mdates.AutoDateLocator()))
    ylim_top = ax.get_ylim()[1]
    ax.set_ylim(0, max(1.2, ylim_top))

    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


# =============================================================================
# Fig. 8 - ETa/ETp vs shallow pressure head (moved here from plot_fig8_ETa_ETp_vs_psi.py)
# =============================================================================
# Reproduces Fig. 8 of Camporese et al. (2014, Adv. Water Resour. 69:95-105): daily ETa/ETp
# against pressure head averaged over the shallowest ~10 cm, one row per scenario (psi_min),
# two columns (outlet vs mid-uphill), with a total-least-squares (orthogonal) trend per panel.
# Needs etp_output.nc (per-surface-node ETp forcing written by the driver next to et_output.nc);
# older runs must be rerun with --forcing etp / etp_ts. "Shallowest 10 cm" = mean of psi over the
# --shallow-depths node set (default 0 and 0.1 m): a mesh-discretisation-dependent proxy.
# ETa/ETp is clipped to [0, 1.05]; days with ETp < etp_min_mmday are dropped (ratio blow-up).
# The pipeline calls fig8_from_artefacts() for every scenario (plot_run_results.py); the
# multi-scenario figure (hypothesis 2, PMIN) feeds several entries to plot_fig8().

LOCATIONS = ("outlet", "mid")
LOCATION_LABELS = {
    "outlet": "Near outlet (creek)",
    "mid":    "Upslope (mid-elevation)",
}
LOCATION_COLORS = {
    "outlet": "#2166ac",
    "mid":    "#d6604d",
}


# =============================================================================
# Artefact loading
# =============================================================================

def load_etp_forcing(out_dir: Path) -> xr.Dataset:
    """Load etp_output.nc (per-surface-node ETp forcing, dims datetime x node).

    Written by the ETp-persistence patch in Agramon_withLAI_withETp.py,
    right after out_dir is created — a static (time, node) artefact, not
    accumulated through the monthly loop like psi/sw/et.
    """
    etp_nc = out_dir / "etp_output.nc"
    if not etp_nc.exists():
        raise FileNotFoundError(
            f"'etp_output.nc' not found in {out_dir}.\n"
            "  -> This artefact is written by the ETp-persistence patch in "
            "Agramon_withLAI_withETp.py (right after out_dir is created). "
            "Rerun this scenario with --forcing etp or etp_ts to generate it."
        )
    return xr.open_dataset(etp_nc)


# =============================================================================
# Reference locations — mirrors resolve_nodes() in plot_run_results.py, but
# also keeps the raw (X, Y) coordinates (needed to sample the gridded ETa
# raster) alongside the shallow node ids (needed to sample psi/ETp).
# =============================================================================

def resolve_shallow_nodes_and_xy(
    simu,
    grid3d: dict,
    *,
    outlet_xy: tuple[float, float] | None = None,
    uphill_xy: tuple[float, float] | None = None,
    shallow_depths: tuple[float, ...] = (0.0, 0.1),
) -> dict:
    """
    Resolve outlet / mid-uphill locations at each of *shallow_depths*.

    Returns a dict with:
      "outlet_xy", "mid_xy"       — (x, y) tuples, for sampling the ETa raster
      "outlet_nodes", "mid_nodes" — list (one entry per shallow_depths, in
                                     order) of node id (outlet) or list of
                                     node ids (mid — a band, mirroring
                                     resolve_nodes()'s mid-elevation handling)
    """
    nnod          = int(grid3d["nnod"])
    nodes         = grid3d["mesh3d_nodes"]
    surface_nodes = nodes[:nnod]

    if outlet_xy:
        ox, oy = outlet_xy
        dists  = np.hypot(surface_nodes[:, 0] - ox, surface_nodes[:, 1] - oy)
        outlet_idx = int(dists.argmin())
    else:
        outlet_idx = int(surface_nodes[:, 2].argmin())
    outlet_x, outlet_y, outlet_z = surface_nodes[outlet_idx]

    mid_band_idx = None
    if uphill_xy:
        ux, uy  = uphill_xy
        dists   = np.hypot(surface_nodes[:, 0] - ux, surface_nodes[:, 1] - uy)
        mid_idx = int(dists.argmin())
        mid_x, mid_y, mid_z = surface_nodes[mid_idx]
    else:
        from plot_hydro_xr import find_representative_mid_node   # lazy: avoids an import cycle
        mid_band_idx = find_representative_mid_node(surface_nodes, return_band=True)
        mid_idx      = find_representative_mid_node(surface_nodes, return_band=False)
        mid_x, mid_y, mid_z = surface_nodes[mid_idx]

    print(f"  Outlet: X={outlet_x:.2f}, Y={outlet_y:.2f}, Z={outlet_z:.2f}")
    print(f"  Mid   : X={mid_x:.2f}, Y={mid_y:.2f}, Z={mid_z:.2f}"
          + (f"  (band of {len(mid_band_idx)} nodes)" if mid_band_idx is not None else "  (single point)"))

    outlet_nodes_by_depth: list = []
    mid_nodes_by_depth:    list = []
    for depth in shallow_depths:
        outlet_nodes_by_depth.append(
            unwrap_node_id(simu.find_nearest_node([outlet_x, outlet_y, outlet_z - depth]))
        )
        if mid_band_idx is not None:
            if depth == 0:
                # Surface band indices are already node indices at depth 0.
                mid_nodes_by_depth.append([int(i) for i in mid_band_idx])
            else:
                mid_nodes_by_depth.append([
                    unwrap_node_id(
                        simu.find_nearest_node(
                            [surface_nodes[i, 0], surface_nodes[i, 1], surface_nodes[i, 2] - depth]
                        )
                    )
                    for i in mid_band_idx
                ])
        else:
            mid_nodes_by_depth.append(
                unwrap_node_id(simu.find_nearest_node([mid_x, mid_y, mid_z - depth]))
            )

    return dict(
        outlet_xy=(float(outlet_x), float(outlet_y)),
        mid_xy=(float(mid_x), float(mid_y)),
        outlet_nodes=outlet_nodes_by_depth,
        mid_nodes=mid_nodes_by_depth,
    )


# =============================================================================
# Series extraction
# =============================================================================

def _mean_over_nodes(da: xr.DataArray, node_ids) -> np.ndarray:
    """Per-timestep mean across one or more node ids (scalar id or a band)."""
    if isinstance(node_ids, (list, tuple, np.ndarray)):
        ids = np.atleast_1d(np.asarray(node_ids, dtype=int))
        return da.sel(node=ids).mean(dim="node").values
    return da.sel(node=int(node_ids)).values


def shallow_psi_series(psi_xr: xr.Dataset, nodes_by_depth: list) -> pd.Series:
    """Mean psi across the requested shallow depths — the paper's 'pressure
    head averaged over the shallowest 10 cm of soil' (see module docstring
    for the discretisation caveat)."""
    da       = psi_xr["psi"]
    time_dim = _infer_time_dim(psi_xr)
    per_depth = [_mean_over_nodes(da, ids) for ids in nodes_by_depth]
    psi_mean  = np.mean(np.stack(per_depth, axis=0), axis=0)
    return pd.Series(psi_mean, index=pd.DatetimeIndex(da[time_dim].values))


def etp_series_at_surface(etp_xr: xr.Dataset, surface_node_ids) -> pd.Series:
    """ETp [m/s] at the surface node(s) for one location (the depth=0 entry
    of resolve_shallow_nodes_and_xy's node lists)."""
    da       = etp_xr["ETp"]
    time_dim = _infer_time_dim(etp_xr)
    vals     = _mean_over_nodes(da, surface_node_ids)
    return pd.Series(vals, index=pd.DatetimeIndex(da[time_dim].values))


def eta_series_at_xy(ET_xr_all: xr.Dataset, et_var: str, xy: tuple[float, float]) -> pd.Series:
    """ETa [m/s] at the nearest et_output.nc raster pixel to *xy*."""
    x, y     = xy
    time_dim = _infer_time_dim(ET_xr_all)
    da       = ET_xr_all[et_var].sel(X=x, Y=y, method="nearest")
    return pd.Series(da.values, index=pd.DatetimeIndex(da[time_dim].values))


def build_location_series(
    *,
    psi_xr: xr.Dataset,
    etp_xr: xr.Dataset,
    ET_xr_all: xr.Dataset,
    et_var: str,
    resolved: dict,
    location: str,
    etp_min_mmday: float,
) -> pd.DataFrame:
    """Coincident daily (psi, ETa/ETp) pairs for one location."""
    nodes_by_depth   = resolved[f"{location}_nodes"]
    surface_node_ids = nodes_by_depth[0]           # the depth=0 entry
    xy               = resolved[f"{location}_xy"]

    psi_s = shallow_psi_series(psi_xr, nodes_by_depth)
    etp_s = etp_series_at_surface(etp_xr, surface_node_ids)
    eta_s = eta_series_at_xy(ET_xr_all, et_var, xy)

    df = pd.concat({"psi": psi_s, "etp": etp_s, "eta": eta_s}, axis=1).dropna()

    etp_min_ms = etp_min_mmday * 1e-3 / 86400   # mm/day -> m/s
    df = df[df["etp"] > etp_min_ms].copy()
    df["ratio"] = (df["eta"] / df["etp"]).clip(lower=0.0, upper=1.05)
    return df[["psi", "ratio"]]


# =============================================================================
# Fitting — "standard total least squares", as in the paper's Fig. 8
# =============================================================================

def fit_tls(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Orthogonal-distance (total-least-squares) linear fit y = a*x + b."""
    def _linear(beta, x):
        return beta[0] * x + beta[1]

    p0  = np.polyfit(x, y, 1)   # OLS slope/intercept as the ODR initial guess
    odr = ODR(RealData(x, y), Model(_linear), beta0=p0)
    out = odr.run()
    return float(out.beta[0]), float(out.beta[1])


# =============================================================================
# Scenario loading
# =============================================================================


def fig8_entry_from_artefacts(
    simu, psi_xr: xr.Dataset, ET_xr_all: xr.Dataset, out_dir: Path, *,
    pmin: float, et_var: str,
    outlet_xy=None, uphill_xy=None,
    shallow_depths: tuple[float, ...] = (0.0, 0.1),
    etp_min_mmday: float = 0.5,
) -> dict:
    """One ``scenario_data`` entry for plot_fig8 from already-loaded artefacts (no log lookup)."""
    if simu is None:
        raise RuntimeError("CATHY project could not be loaded - node resolution needs it.")
    grid3d = simu.read_outputs("grid3d")
    etp_xr = load_etp_forcing(Path(out_dir))
    resolved = resolve_shallow_nodes_and_xy(simu, grid3d, outlet_xy=outlet_xy,
                                            uphill_xy=uphill_xy, shallow_depths=shallow_depths)
    entry: dict = {"pmin": float(pmin), "out_dir": Path(out_dir)}
    for location in LOCATIONS:
        entry[location] = build_location_series(
            psi_xr=psi_xr, etp_xr=etp_xr, ET_xr_all=ET_xr_all, et_var=et_var,
            resolved=resolved, location=location, etp_min_mmday=etp_min_mmday)
    return entry


def plot_fig8(
    scenario_data: dict,
    fig_dir: Path,
    *,
    dpi: int = 150,
    show: bool = True,
) -> Path:
    sim_indices = sorted(scenario_data)
    n_rows = len(sim_indices)
    fig, axes = plt.subplots(
        n_rows, 2, figsize=(9, 3.2 * n_rows), squeeze=False, sharey=True
    )

    for row, sim_index in enumerate(sim_indices):
        entry = scenario_data[sim_index]
        pmin  = entry["pmin"]
        for col, location in enumerate(LOCATIONS):
            ax = axes[row][col]
            df = entry[location]

            if df.empty:
                ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
                ax.set_axis_off()
                continue

            ax.scatter(
                df["psi"], df["ratio"], s=10, alpha=0.35,
                color=LOCATION_COLORS[location], edgecolors="none",
            )

            if len(df) >= 5:
                slope, intercept = fit_tls(df["psi"].to_numpy(), df["ratio"].to_numpy())
                x_line = np.linspace(df["psi"].min(), df["psi"].max(), 50)
                ax.plot(x_line, slope * x_line + intercept, "k--", lw=1.3)
                ax.text(
                    0.03, 0.93, f"TLS slope = {slope:.3g}",
                    transform=ax.transAxes, fontsize=7.5, va="top",
                )

            ax.axhline(0, color="0.75", lw=0.6)
            ax.axhline(1, color="0.75", lw=0.6)
            ax.set_ylim(-0.05, 1.1)

            if row == 0:
                ax.set_title(LOCATION_LABELS[location], fontsize=10, fontweight="bold")
            if row == n_rows - 1:
                ax.set_xlabel(r"$\psi$ (m)")
            if col == 0:
                ax.set_ylabel(f"$\\psi_{{min}}$ = {pmin:g} m\n" + r"$ET_a / ET_p$")

    fig.suptitle(
        r"$ET_a/ET_p$ vs. shallow ($\leq$10 cm) pressure head"
        " — cf. Fig. 8, Camporese et al. (2014)",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    fig_dir.mkdir(parents=True, exist_ok=True)
    out = save_fig(fig, fig_dir, "fig8_ETa_ETp_vs_psi.png", dpi=dpi)
    maybe_show(fig, show)
    return out


