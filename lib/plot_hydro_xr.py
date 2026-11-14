"""
plot_hydro_xr.py - xarray-based hydrology figures at reference nodes (psi, sw, ETa, ETp, LAI) and their fire-window / weekly variants.
====================================================================================================================================
Split unchanged out of results_plotter.py (REFACTORING.md 2.11); results_plotter re-exports every name.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr

from plot_common import COLORS_HYDRO, _fire_month_window, _fire_year_window, _fmt_node_id, _infer_time_dim, _node_ids, _values_at_nodes, _xr_time_to_dates, add_fire_vline, date_window, maybe_show, save_fig


def find_representative_mid_node(
    nodes: np.ndarray,
    margin: float = 0.05,
    dz_pct: float = 0.02,
    return_band: bool = False,
) -> int | np.ndarray:
    """
    Pick a representative node (or node band) from the middle elevation range.

    Instead of picking the single node closest to the median elevation —
    which can easily land right on the mesh edge — this function:
    1. Filters out boundary nodes (using ``margin``)
    2. Finds the median elevation of non-boundary nodes
    3. Identifies all non-boundary nodes within ±``dz_pct`` of that median
       elevation (a *tiny* band around mid-elevation, not a single point)
    4. Either:
       - ``return_band=False`` (default): computes the median (X, Y, Z)
         position of that band and returns the single node closest to it
         — i.e. a *recentred* mid-elevation node, guaranteed off the
         boundary and not just whichever node happens to be nearest in Z.
       - ``return_band=True``: returns the indices of every node in the
         band, so the caller can take the per-timestep *median* across the
         whole band (see ``_values_at_nodes`` / ``_node_ids``) rather than
         relying on any single node.

    Parameters
    ----------
    nodes : np.ndarray
        Array of shape (N, 3) with X, Y, Z coordinates for each node
    margin : float
        Fraction of X/Y extent to exclude as boundary (default: 0.05)
    dz_pct : float
        Elevation-band half-width, as a fraction of the median elevation
        (default: 0.02, i.e. ±2% — a tiny range around mid-elevation)
    return_band : bool
        If True, return all node indices within the band instead of
        collapsing to the single closest-to-median-position node.

    Returns
    -------
    int or np.ndarray
        Index of the representative mid-elevation node, or (if
        ``return_band``) an array of node indices making up the band.
    """
    # Extract coordinates
    x_coords = nodes[:, 0]
    y_coords = nodes[:, 1]
    z_coords = nodes[:, 2]

    # Filter out boundary nodes
    x_min, x_max = np.nanmin(x_coords), np.nanmax(x_coords)
    y_min, y_max = np.nanmin(y_coords), np.nanmax(y_coords)
    x_range = x_max - x_min
    y_range = y_max - y_min

    non_boundary = (
        (x_coords >= x_min + margin * x_range) &
        (x_coords <= x_max - margin * x_range) &
        (y_coords >= y_min + margin * y_range) &
        (y_coords <= y_max - margin * y_range)
    )

    # Work with non-boundary nodes only
    elevations_nb = z_coords[non_boundary]
    median_elev = np.nanmedian(elevations_nb)

    # Define elevation band: ±dz_pct of median elevation
    elev_tolerance = dz_pct * median_elev
    elev_min = median_elev - elev_tolerance
    elev_max = median_elev + elev_tolerance

    # Find nodes within this elevation band AND not on boundary
    elev_mask = (
        (z_coords >= elev_min) &
        (z_coords <= elev_max) &
        non_boundary
    )

    cluster_nodes = nodes[elev_mask]

    if len(cluster_nodes) == 0:
        # Fallback: use original method (closest to median elevation)
        print(f"  WARNING: No nodes found in ±{dz_pct:.1%} elevation band, "
              "falling back to median elevation node")
        all_elevations = z_coords[non_boundary]
        median_elev_all = np.nanmedian(all_elevations)
        distances = np.abs(all_elevations - median_elev_all)
        closest_idx_local = np.argmin(distances)
        non_boundary_indices = np.where(non_boundary)[0]
        fallback_idx = non_boundary_indices[closest_idx_local]
        return np.array([fallback_idx]) if return_band else int(fallback_idx)

    if return_band:
        return np.where(elev_mask)[0]

    # Compute median position of the cluster
    median_x = np.nanmedian(cluster_nodes[:, 0])
    median_y = np.nanmedian(cluster_nodes[:, 1])
    median_z = np.nanmedian(cluster_nodes[:, 2])
    median_pos = np.array([median_x, median_y, median_z])

    # Find node closest to this median position (from non-boundary nodes)
    non_boundary_nodes = nodes[non_boundary]
    distances = np.linalg.norm(non_boundary_nodes - median_pos, axis=1)
    closest_idx_local = np.argmin(distances)

    # Map back to original index
    non_boundary_indices = np.where(non_boundary)[0]
    return int(non_boundary_indices[closest_idx_local])


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
    bar_width:   float = 1.0,
    rain_unit:   str   = "mm/day",
) -> None:
    """Fig 6 — Rain/ETp + ψ + sw time-series built entirely from xarray Datasets.

    Parameters
    ----------
    psi_xr      : xr.Dataset with dimension ``(datetime|time, node)`` and at
                  least one variable containing pressure-head values (m).
    sw_xr       : xr.Dataset with the same structure for soil-water content (−).
    nodes_dict  : dict returned by ``resolve_nodes()``.
    fig_path    : output directory.
    pev_series  : optional potential ET series (mm/day, DatetimeIndex) —
                  satellite ETp (ET0) by default; ERA5 also accepted.
    tp_series   : optional precipitation series (mm/day, DatetimeIndex) —
                  satellite rain (TP) by default; ERA5 also accepted.
    date_start/end : ISO-date strings for the x-axis window.
    show_rain/psi/sw : toggle the three main panels.
    show_uh_s/uh_1/dh_s/dh_1 : toggle individual node series within psi/sw.
    colors_dict : override COLORS_HYDRO.
    scenario_id : label for the figure title.
    bar_width   : width (in day units) of the rain/ETp bars — widen this for
                  coarser (e.g. weekly) resampled data so bars stay visible.
    rain_unit   : y-axis unit label for the rain/ETp panel (e.g. "mm/week").
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
    # "mid_z-*" may be a single node or a whole mid-elevation band (a list of
    # node ids within ±dz_pct of the target elevation) — see
    # find_representative_mid_node(..., return_band=True) in resolve_nodes().
    uh_s_id = _node_ids(nodes_dict, "mid_z-0m")
    uh_1_id = _node_ids(nodes_dict, "mid_z-1m")
    dh_s_id = _node_ids(nodes_dict, "outlet_z-0m")
    dh_1_id = _node_ids(nodes_dict, "outlet_z-1m")

    psi_mask = (psi_dates >= d0) & (psi_dates <= d1)
    sw_mask  = (sw_dates  >= d0) & (sw_dates  <= d1)

    # ── psi arrays at reference nodes  (shape: T) ────────────────────────────
    psi_da = psi_xr[psi_var]   # (datetime, node)
    sw_da  = sw_xr[sw_var]

    def _psi_at(node_id: int | list[int]) -> np.ndarray:
        return _values_at_nodes(psi_da, psi_time_dim, psi_mask, node_id)

    def _sw_at(node_id: int | list[int]) -> np.ndarray:
        return _values_at_nodes(sw_da, sw_time_dim, sw_mask, node_id)

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
            pev_m = (pev_series.index >= d0) & (pev_series.index <= d1)
            tp_m  = (tp_series.index  >= d0) & (tp_series.index  <= d1)
            ax.bar(pev_series[pev_m].index, pev_series[pev_m].values,
                   width=bar_width, color="grey",   label="ETp", alpha=0.5)
            ax.bar(tp_series[tp_m].index,  tp_series[tp_m].values,
                   width=bar_width, color="skyblue", label="Rain")
            ax.invert_yaxis()
        else:
            ax.text(0.5, 0.5, "Rain/ETp forcing not available",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="#888888")
        ax.set_ylabel(rain_unit)
        ax.set_title("Rain / ETp (satellite)", fontsize=10, fontweight="bold")
        add_fire_vline(ax, d0, d1)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, linestyle="--", alpha=0.3)

    # ── Panel 2: ψ at reference nodes ─────────────────────────────────────────
    if show_psi:
        ax = next(ax_iter)
        psi_dates_w = psi_dates[psi_mask]
        if show_uh_s:
            ax.plot(psi_dates_w, _psi_at(uh_s_id),
                    label=f"Mid surf n={_fmt_node_id(uh_s_id)}",
                    color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
        if show_uh_1:
            ax.plot(psi_dates_w, _psi_at(uh_1_id),
                    label=f"Mid −1m n={_fmt_node_id(uh_1_id)}",
                    color=colors_dict["uphill_1m"], linestyle="--")
        if show_dh_s:
            ax.plot(psi_dates_w, _psi_at(dh_s_id),
                    label=f"Outlet surf n={_fmt_node_id(dh_s_id)}",
                    color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
        if show_dh_1:
            ax.plot(psi_dates_w, _psi_at(dh_1_id),
                    label=f"Outlet −1m n={_fmt_node_id(dh_1_id)}",
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
                    label=f"Mid surf n={_fmt_node_id(uh_s_id)}",
                    color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
        if show_uh_1:
            ax.plot(sw_dates_w, _sw_at(uh_1_id),
                    label=f"Mid −1m n={_fmt_node_id(uh_1_id)}",
                    color=colors_dict["uphill_1m"], linestyle="--")
        if show_dh_s:
            ax.plot(sw_dates_w, _sw_at(dh_s_id),
                    label=f"Outlet surf n={_fmt_node_id(dh_s_id)}",
                    color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
        if show_dh_1:
            ax.plot(sw_dates_w, _sw_at(dh_1_id),
                    label=f"Outlet −1m n={_fmt_node_id(dh_1_id)}",
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
    et_unit:     str  = "mm/day",
    et_as_line:  bool = False,
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
    et_unit    : y-axis unit label for the ETa panel (e.g. "mm/week" for
                 resampled data).
    et_as_line : draw the ETa panel as a connected line+marker instead of a
                 bare scatter — clearer for coarser (e.g. weekly) series.
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

    # Node IDs — band-preserving: "mid_z-*" may be a single node or a whole
    # mid-elevation band (see find_representative_mid_node(..., return_band=True)
    # in resolve_nodes()); _node_ids keeps that band intact instead of
    # collapsing it to a single id.
    uh_s_id = _node_ids(nodes_dict, "mid_z-0m")
    uh_1_id = _node_ids(nodes_dict, "mid_z-1m")
    dh_s_id = _node_ids(nodes_dict, "outlet_z-0m")
    dh_1_id = _node_ids(nodes_dict, "outlet_z-1m")

    sw_dates_w = sw_dates[sw_mask]

    def _sw_at(node_id: int | list[int]) -> np.ndarray:
        return _values_at_nodes(sw_da, sw_time_dim, sw_mask, node_id)

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(14, 6), sharex=True, constrained_layout=True
    )

    # ── Top: sw at reference nodes ────────────────────────────────────────────
    if show_uh_s:
        ax0.plot(sw_dates_w, _sw_at(uh_s_id),
                 label=f"Mid surf n={_fmt_node_id(uh_s_id)}",
                 color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
    if show_uh_1:
        ax0.plot(sw_dates_w, _sw_at(uh_1_id),
                 label=f"Mid −1m n={_fmt_node_id(uh_1_id)}",
                 color=colors_dict["uphill_1m"], linestyle="--")
    if show_dh_s:
        ax0.plot(sw_dates_w, _sw_at(dh_s_id),
                 label=f"Outlet surf n={_fmt_node_id(dh_s_id)}",
                 color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
    if show_dh_1:
        ax0.plot(sw_dates_w, _sw_at(dh_1_id),
                 label=f"Outlet −1m n={_fmt_node_id(dh_1_id)}",
                 color=colors_dict["downhill_1m"], linestyle="--")

    ax0.set_ylabel("sw (−)")
    ax0.set_title("Soil Water Content at reference nodes", fontsize=10, fontweight="bold")
    add_fire_vline(ax0, d0, d1)
    ax0.legend(loc="upper right", fontsize=8)
    ax0.grid(True, linestyle="--", alpha=0.3)

    # ── Bottom: spatial-mean ETa ──────────────────────────────────────────────
    if et_var in ET_xr_all and len(et_vals):
        if et_as_line:
            ax1.plot(et_dates_w, et_vals, color=colors_dict["ET_mean"],
                     marker="o", ms=4, linewidth=1.3, label="ETa spatial mean")
        else:
            ax1.scatter(et_dates_w, et_vals,
                        c=colors_dict["ET_mean"], s=14, label="ETa spatial mean")
    else:
        ax1.text(0.5, 0.5, f"'{et_var}' not found in ET dataset",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=9, color="#888888")

    ax1.set_ylabel(f"ETa ({et_unit})")
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


def plot_hydro_xr_fire_year(
    psi_xr:     xr.Dataset,
    sw_xr:      xr.Dataset,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    pev_series:  pd.Series | None = None,
    tp_series:   pd.Series | None = None,
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
    fname:       str  = "fig6c_hydro_timeseries_fire_years.png",
) -> None:
    """Fig 6c -- Fig 6 zoomed to the fire year + the following year
    (e.g. 2020-01-01 -> 2021-12-31 for ``FIRE_DATE = 2020-07-01``)."""
    d0, d1 = _fire_year_window()
    plot_hydro_xr(
        psi_xr, sw_xr, nodes_dict, fig_path,
        pev_series=pev_series, tp_series=tp_series,
        date_start=str(d0.date()), date_end=str(d1.date()),
        show_rain=show_rain, show_psi=show_psi, show_sw=show_sw,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
    )


def plot_hydro_xr_fire_month(
    psi_xr:     xr.Dataset,
    sw_xr:      xr.Dataset,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    pev_series:  pd.Series | None = None,
    tp_series:   pd.Series | None = None,
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
    fname:       str  = "fig6d_hydro_timeseries_fire_month.png",
) -> None:
    """Fig 6d -- Fig 6 zoomed to just the calendar month of the fire
    (e.g. July 2020 for ``FIRE_DATE = 2020-07-01``)."""
    d0, d1 = _fire_month_window()
    plot_hydro_xr(
        psi_xr, sw_xr, nodes_dict, fig_path,
        pev_series=pev_series, tp_series=tp_series,
        date_start=str(d0.date()), date_end=str(d1.date()),
        show_rain=show_rain, show_psi=show_psi, show_sw=show_sw,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
    )


def plot_hydro_xr_weekly(
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
    fname:       str  = "fig6e_hydro_timeseries_weekly.png",
) -> None:
    """Fig 6e -- Weekly-resampled version of Fig 6.

    psi and sw (state variables) are resampled with a weekly *mean*; rain and
    ETp (fluxes) are resampled with a weekly *sum*, so the bars read as
    "mm accumulated that week". Falls back to the full record if
    *date_start*/*date_end* are not given.
    """
    psi_time_dim = _infer_time_dim(psi_xr)
    sw_time_dim  = _infer_time_dim(sw_xr)
    psi_weekly = psi_xr.resample({psi_time_dim: "W"}).mean()
    sw_weekly  = sw_xr.resample({sw_time_dim: "W"}).mean()

    pev_weekly = pev_series.resample("W").sum() if pev_series is not None and len(pev_series) else pev_series
    tp_weekly  = tp_series.resample("W").sum()  if tp_series  is not None and len(tp_series)  else tp_series

    plot_hydro_xr(
        psi_weekly, sw_weekly, nodes_dict, fig_path,
        pev_series=pev_weekly, tp_series=tp_weekly,
        date_start=date_start, date_end=date_end,
        show_rain=show_rain, show_psi=show_psi, show_sw=show_sw,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        bar_width=6.0, rain_unit="mm/week",
    )


def plot_et_sw_xr_fire_year(
    sw_xr:      xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig7c_et_sw_timeseries_fire_years.png",
) -> None:
    """Fig 7c -- Fig 7 zoomed to the fire year + the following year."""
    d0, d1 = _fire_year_window()
    plot_et_sw_xr(
        sw_xr, ET_xr_all, et_var, et_scale, nodes_dict, fig_path,
        date_start=str(d0.date()), date_end=str(d1.date()),
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
    )


def plot_et_sw_xr_fire_month(
    sw_xr:      xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig7d_et_sw_timeseries_fire_month.png",
) -> None:
    """Fig 7d -- Fig 7 zoomed to just the calendar month of the fire."""
    d0, d1 = _fire_month_window()
    plot_et_sw_xr(
        sw_xr, ET_xr_all, et_var, et_scale, nodes_dict, fig_path,
        date_start=str(d0.date()), date_end=str(d1.date()),
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
    )


def plot_et_sw_xr_weekly(
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
    fname:       str  = "fig7e_et_sw_timeseries_weekly.png",
) -> None:
    """Fig 7e -- Weekly-resampled version of Fig 7.

    sw (state variable) is resampled with a weekly *mean*; ETa (a flux) is
    resampled with a weekly *sum*, so the panel reads as "mm accumulated
    that week". Falls back to the full record if *date_start*/*date_end*
    are not given.
    """
    sw_time_dim = _infer_time_dim(sw_xr)
    et_time_dim = _infer_time_dim(ET_xr_all)
    sw_weekly = sw_xr.resample({sw_time_dim: "W"}).mean()
    et_weekly = ET_xr_all.resample({et_time_dim: "W"}).sum()

    plot_et_sw_xr(
        sw_weekly, et_weekly, et_var, et_scale, nodes_dict, fig_path,
        date_start=date_start, date_end=date_end,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        et_unit="mm/week", et_as_line=True,
    )


def plot_et_psi_xr(
    psi_xr:     xr.Dataset,
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
    fname:       str  = "fig7g_et_psi_timeseries.png",
    et_unit:     str  = "mm/day",
    et_as_line:  bool = False,
) -> None:
    """Fig 7g -- same as ``plot_et_sw_xr`` (Fig 7), but the top panel plots
    pressure head (psi) at the reference nodes instead of soil water
    content (sw).

    Parameters
    ----------
    psi_xr     : xr.Dataset with dimension ``(datetime|time, node)`` for psi.
    ET_xr_all  : xr.Dataset with spatial ET (X, Y, datetime|time), variable
                 name given by *et_var*.
    et_var     : name of the ET variable inside *ET_xr_all*.
    et_scale   : multiplicative factor (e.g. 1e3*86400 for m/s → mm/day).
    nodes_dict : dict returned by ``resolve_nodes()``.
    fig_path   : output directory.
    date_start/end : ISO-date strings for the x-axis window.
    show_uh_s/uh_1/dh_s/dh_1 : toggle individual psi node series.
    colors_dict : override COLORS_HYDRO.
    scenario_id : label for the figure title.
    et_unit    : y-axis unit label for the ETa panel (e.g. "mm/week" for
                 resampled data).
    et_as_line : draw the ETa panel as a connected line+marker instead of a
                 bare scatter — clearer for coarser (e.g. weekly) series.
    """
    colors_dict = colors_dict or COLORS_HYDRO

    psi_time_dim = _infer_time_dim(psi_xr)
    et_time_dim  = _infer_time_dim(ET_xr_all)
    psi_var      = list(psi_xr.data_vars)[0]

    psi_dates = _xr_time_to_dates(psi_xr[psi_var], psi_time_dim)
    et_dates  = _xr_time_to_dates(ET_xr_all[et_var], et_time_dim)

    d0, d1 = date_window(date_start, date_end, psi_dates.to_series())
    print(f"  Plotting et_psi_xr time-series ({d0.date()} → {d1.date()}) …")

    psi_mask = (psi_dates >= d0) & (psi_dates <= d1)
    et_mask  = (et_dates  >= d0) & (et_dates  <= d1)

    psi_da = psi_xr[psi_var]

    # Spatial-mean ET (all non-time dims averaged)
    spatial_dims = [d for d in ET_xr_all[et_var].dims
                    if d.lower() not in ("datetime", "time", "t")]
    et_mean = ET_xr_all[et_var].mean(dim=spatial_dims) * et_scale  # (time,)
    et_vals = et_mean.values[et_mask]
    et_dates_w = et_dates[et_mask]

    # Node IDs — band-preserving: "mid_z-*" may be a single node or a whole
    # mid-elevation band (see find_representative_mid_node(..., return_band=True)
    # in resolve_nodes()); _node_ids keeps that band intact instead of
    # collapsing it to a single id.
    uh_s_id = _node_ids(nodes_dict, "mid_z-0m")
    uh_1_id = _node_ids(nodes_dict, "mid_z-1m")
    dh_s_id = _node_ids(nodes_dict, "outlet_z-0m")
    dh_1_id = _node_ids(nodes_dict, "outlet_z-1m")

    psi_dates_w = psi_dates[psi_mask]

    def _psi_at(node_id: int | list[int]) -> np.ndarray:
        return _values_at_nodes(psi_da, psi_time_dim, psi_mask, node_id)

    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(14, 6), sharex=True, constrained_layout=True
    )

    # ── Top: psi at reference nodes ───────────────────────────────────────────
    if show_uh_s:
        ax0.plot(psi_dates_w, _psi_at(uh_s_id),
                 label=f"Mid surf n={_fmt_node_id(uh_s_id)}",
                 color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
    if show_uh_1:
        ax0.plot(psi_dates_w, _psi_at(uh_1_id),
                 label=f"Mid −1m n={_fmt_node_id(uh_1_id)}",
                 color=colors_dict["uphill_1m"], linestyle="--")
    if show_dh_s:
        ax0.plot(psi_dates_w, _psi_at(dh_s_id),
                 label=f"Outlet surf n={_fmt_node_id(dh_s_id)}",
                 color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
    if show_dh_1:
        ax0.plot(psi_dates_w, _psi_at(dh_1_id),
                 label=f"Outlet −1m n={_fmt_node_id(dh_1_id)}",
                 color=colors_dict["downhill_1m"], linestyle="--")

    ax0.set_ylabel("ψ (m)")
    ax0.set_title("Pressure Head ψ at reference nodes", fontsize=10, fontweight="bold")
    add_fire_vline(ax0, d0, d1)
    ax0.legend(loc="upper right", fontsize=8)
    ax0.grid(True, linestyle="--", alpha=0.3)
    ax0.axhline(0, color="black", linewidth=0.8, linestyle=":")

    # ── Bottom: spatial-mean ETa ──────────────────────────────────────────────
    if et_var in ET_xr_all and len(et_vals):
        if et_as_line:
            ax1.plot(et_dates_w, et_vals, color=colors_dict["ET_mean"],
                     marker="o", ms=4, linewidth=1.3, label="ETa spatial mean")
        else:
            ax1.scatter(et_dates_w, et_vals,
                        c=colors_dict["ET_mean"], s=14, label="ETa spatial mean")
    else:
        ax1.text(0.5, 0.5, f"'{et_var}' not found in ET dataset",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=9, color="#888888")

    ax1.set_ylabel(f"ETa ({et_unit})")
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
        f"Sc.{scenario_id} | psi & ETa | {d0.date()} → {d1.date()}",
        fontsize=11, fontweight="bold",
    )
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_et_psi_xr_fire_year(
    psi_xr:     xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig7gc_et_psi_timeseries_fire_years.png",
) -> None:
    """Fig 7gc -- Fig 7g (psi + ETa) zoomed to the fire year + the following
    year — the psi counterpart of ``plot_et_sw_xr_fire_year`` (Fig 7c)."""
    d0, d1 = _fire_year_window()
    plot_et_psi_xr(
        psi_xr, ET_xr_all, et_var, et_scale, nodes_dict, fig_path,
        date_start=str(d0.date()), date_end=str(d1.date()),
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
    )


def plot_et_sw_lai_xr(
    sw_xr:      xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    lai_ds:     xr.Dataset | None,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    lai_var:     str = "LAI",
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
    fname:       str  = "fig7b_et_sw_lai_timeseries.png",
    et_unit:     str  = "mm/day",
    et_as_line:  bool = False,
) -> None:
    """Fig 7b -- Fig 7 (sw + spatial-mean ETa) with a third panel for the
    domain-mean LAI trend, so vegetation state and the water-balance
    response can be read off the same time axis.

    Parameters
    ----------
    sw_xr, ET_xr_all, et_var, et_scale, nodes_dict, fig_path : same as
        ``plot_et_sw_xr``.
    lai_ds   : xr.Dataset with a time dimension and variable *lai_var*
               (e.g. from ``lai_history_to_dataset(art["lai_history"])``).
               May be ``None`` or empty — the LAI panel then shows a
               "not available" placeholder instead of failing the whole
               figure.
    lai_var  : name of the LAI variable inside *lai_ds* (default 'LAI').
    et_unit / et_as_line : same as ``plot_et_sw_xr`` (used by the weekly
               variant to switch the ETa panel to mm/week + line style).
    """
    colors_dict = colors_dict or COLORS_HYDRO

    sw_time_dim = _infer_time_dim(sw_xr)
    et_time_dim = _infer_time_dim(ET_xr_all)
    sw_var      = list(sw_xr.data_vars)[0]

    sw_dates = _xr_time_to_dates(sw_xr[sw_var], sw_time_dim)
    et_dates = _xr_time_to_dates(ET_xr_all[et_var], et_time_dim)

    d0, d1 = date_window(date_start, date_end, sw_dates.to_series())
    print(f"  Plotting et_sw_lai_xr time-series ({d0.date()} → {d1.date()}) …")

    sw_mask = (sw_dates >= d0) & (sw_dates <= d1)
    et_mask = (et_dates >= d0) & (et_dates <= d1)

    sw_da = sw_xr[sw_var]

    # Spatial-mean ET (all non-time dims averaged)
    spatial_dims = [d for d in ET_xr_all[et_var].dims
                    if d.lower() not in ("datetime", "time", "t")]
    et_mean = ET_xr_all[et_var].mean(dim=spatial_dims) * et_scale  # (time,)
    et_vals = et_mean.values[et_mask]
    et_dates_w = et_dates[et_mask]

    # Node IDs — band-preserving: "mid_z-*" may be a single node or a whole
    # mid-elevation band (see find_representative_mid_node(..., return_band=True)
    # in resolve_nodes()); _node_ids keeps that band intact instead of
    # collapsing it to a single id.
    uh_s_id = _node_ids(nodes_dict, "mid_z-0m")
    uh_1_id = _node_ids(nodes_dict, "mid_z-1m")
    dh_s_id = _node_ids(nodes_dict, "outlet_z-0m")
    dh_1_id = _node_ids(nodes_dict, "outlet_z-1m")

    sw_dates_w = sw_dates[sw_mask]

    def _sw_at(node_id: int | list[int]) -> np.ndarray:
        return _values_at_nodes(sw_da, sw_time_dim, sw_mask, node_id)

    # Domain-mean LAI trend (guarded — lai_ds is optional / may be empty)
    lai_available = (
        lai_ds is not None and len(lai_ds.dims) > 0 and lai_var in lai_ds
    )
    lai_vals: np.ndarray = np.array([])
    lai_dates_w = None
    if lai_available:
        lai_time_dim = _infer_time_dim(lai_ds)
        lai_dates = _xr_time_to_dates(lai_ds[lai_var], lai_time_dim)
        lai_mask  = (lai_dates >= d0) & (lai_dates <= d1)
        lai_spatial_dims = [d for d in lai_ds[lai_var].dims if d != lai_time_dim]
        lai_mean = lai_ds[lai_var].mean(dim=lai_spatial_dims)
        lai_vals = lai_mean.values[lai_mask]
        lai_dates_w = lai_dates[lai_mask]

    fig, (ax0, ax1, ax2) = plt.subplots(
        3, 1, figsize=(14, 8.5), sharex=True, constrained_layout=True
    )

    # ── Top: sw at reference nodes ────────────────────────────────────────────
    if show_uh_s:
        ax0.plot(sw_dates_w, _sw_at(uh_s_id),
                 label=f"Mid surf n={_fmt_node_id(uh_s_id)}",
                 color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
    if show_uh_1:
        ax0.plot(sw_dates_w, _sw_at(uh_1_id),
                 label=f"Mid −1m n={_fmt_node_id(uh_1_id)}",
                 color=colors_dict["uphill_1m"], linestyle="--")
    if show_dh_s:
        ax0.plot(sw_dates_w, _sw_at(dh_s_id),
                 label=f"Outlet surf n={_fmt_node_id(dh_s_id)}",
                 color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
    if show_dh_1:
        ax0.plot(sw_dates_w, _sw_at(dh_1_id),
                 label=f"Outlet −1m n={_fmt_node_id(dh_1_id)}",
                 color=colors_dict["downhill_1m"], linestyle="--")

    ax0.set_ylabel("sw (−)")
    ax0.set_title("Soil Water Content at reference nodes", fontsize=10, fontweight="bold")
    add_fire_vline(ax0, d0, d1)
    ax0.legend(loc="upper right", fontsize=8)
    ax0.grid(True, linestyle="--", alpha=0.3)

    # ── Middle: spatial-mean ETa ──────────────────────────────────────────────
    if et_var in ET_xr_all and len(et_vals):
        if et_as_line:
            ax1.plot(et_dates_w, et_vals, color=colors_dict["ET_mean"],
                     marker="o", ms=4, linewidth=1.3, label="ETa spatial mean")
        else:
            ax1.scatter(et_dates_w, et_vals,
                        c=colors_dict["ET_mean"], s=14, label="ETa spatial mean")
    else:
        ax1.text(0.5, 0.5, f"'{et_var}' not found in ET dataset",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=9, color="#888888")

    ax1.set_ylabel(f"ETa ({et_unit})")
    ax1.set_title(
        f"Actual Evapotranspiration — spatial mean ({et_var})",
        fontsize=10, fontweight="bold",
    )
    add_fire_vline(ax1, d0, d1)
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(True, linestyle="--", alpha=0.3)

    # ── Bottom: domain-mean LAI trend ─────────────────────────────────────────
    if lai_available and len(lai_vals):
        ax2.plot(lai_dates_w, lai_vals,
                 color=colors_dict.get("LAI", "#4daf4a"),
                 marker="s", ms=4, linewidth=1.3, label="LAI domain mean")
    else:
        ax2.text(0.5, 0.5, f"'{lai_var}' not available", ha="center", va="center",
                 transform=ax2.transAxes, fontsize=9, color="#888888")

    ax2.set_ylabel("LAI (m² m⁻²)")
    ax2.set_xlabel("Date")
    ax2.set_title("Leaf Area Index — domain mean", fontsize=10, fontweight="bold")
    add_fire_vline(ax2, d0, d1)
    ax2.legend(loc="upper right", fontsize=8)
    ax2.grid(True, linestyle="--", alpha=0.3)
    ax2.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(mdates.AutoDateLocator())
    )

    fig.suptitle(
        f"Sc.{scenario_id} | sw, ETa & LAI | {d0.date()} → {d1.date()}",
        fontsize=11, fontweight="bold",
    )
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_et_sw_lai_xr_fire_year(
    sw_xr:      xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    lai_ds:     xr.Dataset | None,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    lai_var:     str = "LAI",
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig7bc_et_sw_lai_timeseries_fire_years.png",
) -> None:
    """Fig 7bc -- Fig 7b zoomed to the fire year + the following year."""
    d0, d1 = _fire_year_window()
    plot_et_sw_lai_xr(
        sw_xr, ET_xr_all, et_var, et_scale, lai_ds, nodes_dict, fig_path,
        lai_var=lai_var,
        date_start=str(d0.date()), date_end=str(d1.date()),
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
    )


def plot_et_sw_lai_xr_fire_month(
    sw_xr:      xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    lai_ds:     xr.Dataset | None,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    lai_var:     str = "LAI",
    show_uh_s:   bool = True,
    show_uh_1:   bool = True,
    show_dh_s:   bool = True,
    show_dh_1:   bool = True,
    colors_dict: dict | None = None,
    scenario_id: int | str = "?",
    dpi:         int  = 150,
    show:        bool = False,
    fname:       str  = "fig7bd_et_sw_lai_timeseries_fire_month.png",
) -> None:
    """Fig 7bd -- Fig 7b zoomed to just the calendar month of the fire."""
    d0, d1 = _fire_month_window()
    plot_et_sw_lai_xr(
        sw_xr, ET_xr_all, et_var, et_scale, lai_ds, nodes_dict, fig_path,
        lai_var=lai_var,
        date_start=str(d0.date()), date_end=str(d1.date()),
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
    )


def plot_et_sw_lai_xr_weekly(
    sw_xr:      xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    lai_ds:     xr.Dataset | None,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    lai_var:     str = "LAI",
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
    fname:       str  = "fig7be_et_sw_lai_timeseries_weekly.png",
) -> None:
    """Fig 7be -- Weekly-resampled version of Fig 7b.

    sw and LAI (state variables) are resampled with a weekly *mean*; ETa
    (a flux) is resampled with a weekly *sum* — same convention as
    ``plot_et_sw_xr_weekly``.
    """
    sw_time_dim = _infer_time_dim(sw_xr)
    et_time_dim = _infer_time_dim(ET_xr_all)
    sw_weekly = sw_xr.resample({sw_time_dim: "W"}).mean()
    et_weekly = ET_xr_all.resample({et_time_dim: "W"}).sum()

    lai_weekly = lai_ds
    if lai_ds is not None and len(lai_ds.dims) > 0:
        lai_time_dim = _infer_time_dim(lai_ds)
        lai_weekly = lai_ds.resample({lai_time_dim: "W"}).mean()

    plot_et_sw_lai_xr(
        sw_weekly, et_weekly, et_var, et_scale, lai_weekly, nodes_dict, fig_path,
        lai_var=lai_var,
        date_start=date_start, date_end=date_end,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        et_unit="mm/week", et_as_line=True,
    )


def plot_et_psi_lai_xr(
    psi_xr:     xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    lai_ds:     xr.Dataset | None,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    lai_var:     str = "LAI",
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
    fname:       str  = "fig7f_et_psi_lai_timeseries.png",
    et_unit:     str  = "mm/day",
    et_as_line:  bool = False,
) -> None:
    """Fig 7f -- same as ``plot_et_sw_lai_xr`` (Fig 7b), but the top panel
    plots pressure head (psi) at the reference nodes instead of soil water
    content (sw).

    Parameters
    ----------
    psi_xr, ET_xr_all, et_var, et_scale, nodes_dict, fig_path : same as
        ``plot_hydro_xr``/``plot_et_sw_xr``, but for the psi dataset.
    lai_ds   : xr.Dataset with a time dimension and variable *lai_var*
               (e.g. from ``lai_history_to_dataset(art["lai_history"])``).
               May be ``None`` or empty — the LAI panel then shows a
               "not available" placeholder instead of failing the whole
               figure.
    lai_var  : name of the LAI variable inside *lai_ds* (default 'LAI').
    et_unit / et_as_line : same as ``plot_et_sw_xr`` (used by the weekly
               variant to switch the ETa panel to mm/week + line style).
    """
    colors_dict = colors_dict or COLORS_HYDRO

    psi_time_dim = _infer_time_dim(psi_xr)
    et_time_dim  = _infer_time_dim(ET_xr_all)
    psi_var      = list(psi_xr.data_vars)[0]

    psi_dates = _xr_time_to_dates(psi_xr[psi_var], psi_time_dim)
    et_dates  = _xr_time_to_dates(ET_xr_all[et_var], et_time_dim)

    d0, d1 = date_window(date_start, date_end, psi_dates.to_series())
    print(f"  Plotting et_psi_lai_xr time-series ({d0.date()} → {d1.date()}) …")

    psi_mask = (psi_dates >= d0) & (psi_dates <= d1)
    et_mask  = (et_dates >= d0) & (et_dates <= d1)

    psi_da = psi_xr[psi_var]

    # Spatial-mean ET (all non-time dims averaged)
    spatial_dims = [d for d in ET_xr_all[et_var].dims
                    if d.lower() not in ("datetime", "time", "t")]
    et_mean = ET_xr_all[et_var].mean(dim=spatial_dims) * et_scale  # (time,)
    et_vals = et_mean.values[et_mask]
    et_dates_w = et_dates[et_mask]

    # Node IDs — band-preserving: "mid_z-*" may be a single node or a whole
    # mid-elevation band (see find_representative_mid_node(..., return_band=True)
    # in resolve_nodes()); _node_ids keeps that band intact instead of
    # collapsing it to a single id.
    uh_s_id = _node_ids(nodes_dict, "mid_z-0m")
    uh_1_id = _node_ids(nodes_dict, "mid_z-1m")
    dh_s_id = _node_ids(nodes_dict, "outlet_z-0m")
    dh_1_id = _node_ids(nodes_dict, "outlet_z-1m")

    psi_dates_w = psi_dates[psi_mask]

    def _psi_at(node_id: int | list[int]) -> np.ndarray:
        return _values_at_nodes(psi_da, psi_time_dim, psi_mask, node_id)

    # Domain-mean LAI trend (guarded — lai_ds is optional / may be empty)
    lai_available = (
        lai_ds is not None and len(lai_ds.dims) > 0 and lai_var in lai_ds
    )
    lai_vals: np.ndarray = np.array([])
    lai_dates_w = None
    if lai_available:
        lai_time_dim = _infer_time_dim(lai_ds)
        lai_dates = _xr_time_to_dates(lai_ds[lai_var], lai_time_dim)
        lai_mask  = (lai_dates >= d0) & (lai_dates <= d1)
        lai_spatial_dims = [d for d in lai_ds[lai_var].dims if d != lai_time_dim]
        lai_mean = lai_ds[lai_var].mean(dim=lai_spatial_dims)
        lai_vals = lai_mean.values[lai_mask]
        lai_dates_w = lai_dates[lai_mask]

    fig, (ax0, ax1, ax2) = plt.subplots(
        3, 1, figsize=(14, 8.5), sharex=True, constrained_layout=True
    )

    # ── Top: psi at reference nodes ───────────────────────────────────────────
    if show_uh_s:
        ax0.plot(psi_dates_w, _psi_at(uh_s_id),
                 label=f"Mid surf n={_fmt_node_id(uh_s_id)}",
                 color=colors_dict["uphill_surface"], marker="+", linestyle="-", ms=3)
    if show_uh_1:
        ax0.plot(psi_dates_w, _psi_at(uh_1_id),
                 label=f"Mid −1m n={_fmt_node_id(uh_1_id)}",
                 color=colors_dict["uphill_1m"], linestyle="--")
    if show_dh_s:
        ax0.plot(psi_dates_w, _psi_at(dh_s_id),
                 label=f"Outlet surf n={_fmt_node_id(dh_s_id)}",
                 color=colors_dict["downhill_surface"], marker="o", linestyle="-", ms=3)
    if show_dh_1:
        ax0.plot(psi_dates_w, _psi_at(dh_1_id),
                 label=f"Outlet −1m n={_fmt_node_id(dh_1_id)}",
                 color=colors_dict["downhill_1m"], linestyle="--")

    ax0.set_ylabel("ψ (m)")
    ax0.set_title("Pressure Head at reference nodes", fontsize=10, fontweight="bold")
    add_fire_vline(ax0, d0, d1)
    ax0.legend(loc="upper right", fontsize=8)
    ax0.grid(True, linestyle="--", alpha=0.3)

    # ── Middle: spatial-mean ETa ──────────────────────────────────────────────
    if et_var in ET_xr_all and len(et_vals):
        if et_as_line:
            ax1.plot(et_dates_w, et_vals, color=colors_dict["ET_mean"],
                     marker="o", ms=4, linewidth=1.3, label="ETa spatial mean")
        else:
            ax1.scatter(et_dates_w, et_vals,
                        c=colors_dict["ET_mean"], s=14, label="ETa spatial mean")
    else:
        ax1.text(0.5, 0.5, f"'{et_var}' not found in ET dataset",
                 ha="center", va="center", transform=ax1.transAxes,
                 fontsize=9, color="#888888")

    ax1.set_ylabel(f"ETa ({et_unit})")
    ax1.set_title(
        f"Actual Evapotranspiration — spatial mean ({et_var})",
        fontsize=10, fontweight="bold",
    )
    add_fire_vline(ax1, d0, d1)
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(True, linestyle="--", alpha=0.3)

    # ── Bottom: domain-mean LAI trend ─────────────────────────────────────────
    if lai_available and len(lai_vals):
        ax2.plot(lai_dates_w, lai_vals,
                 color=colors_dict.get("LAI", "#4daf4a"),
                 marker="s", ms=4, linewidth=1.3, label="LAI domain mean")
    else:
        ax2.text(0.5, 0.5, f"'{lai_var}' not available", ha="center", va="center",
                 transform=ax2.transAxes, fontsize=9, color="#888888")

    ax2.set_ylabel("LAI (m² m⁻²)")
    ax2.set_xlabel("Date")
    ax2.set_title("Leaf Area Index — domain mean", fontsize=10, fontweight="bold")
    add_fire_vline(ax2, d0, d1)
    ax2.legend(loc="upper right", fontsize=8)
    ax2.grid(True, linestyle="--", alpha=0.3)
    ax2.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(mdates.AutoDateLocator())
    )

    fig.suptitle(
        f"Sc.{scenario_id} | psi, ETa & LAI | {d0.date()} → {d1.date()}",
        fontsize=11, fontweight="bold",
    )
    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_et_psi_lai_xr_weekly(
    psi_xr:     xr.Dataset,
    ET_xr_all:  xr.Dataset,
    et_var:     str,
    et_scale:   float,
    lai_ds:     xr.Dataset | None,
    nodes_dict: dict,
    fig_path:   Path,
    *,
    lai_var:     str = "LAI",
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
    fname:       str  = "fig7fe_et_psi_lai_timeseries_weekly.png",
) -> None:
    """Fig 7fe -- Weekly-resampled version of Fig 7f.

    psi and LAI (state variables) are resampled with a weekly *mean*; ETa
    (a flux) is resampled with a weekly *sum* — same convention as
    ``plot_et_sw_lai_xr_weekly``/``plot_et_sw_xr_weekly``.
    """
    psi_time_dim = _infer_time_dim(psi_xr)
    et_time_dim  = _infer_time_dim(ET_xr_all)
    psi_weekly = psi_xr.resample({psi_time_dim: "W"}).mean()
    et_weekly  = ET_xr_all.resample({et_time_dim: "W"}).sum()

    lai_weekly = lai_ds
    if lai_ds is not None and len(lai_ds.dims) > 0:
        lai_time_dim = _infer_time_dim(lai_ds)
        lai_weekly = lai_ds.resample({lai_time_dim: "W"}).mean()

    plot_et_psi_lai_xr(
        psi_weekly, et_weekly, et_var, et_scale, lai_weekly, nodes_dict, fig_path,
        lai_var=lai_var,
        date_start=date_start, date_end=date_end,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        et_unit="mm/week", et_as_line=True,
    )
