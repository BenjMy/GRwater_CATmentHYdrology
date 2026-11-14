"""
plot_hydro_legacy.py - Pandas-based hydrology figures: hydrographs, recharge, node locations, ET/sw.
==================================================================================================
Split unchanged out of results_plotter.py (REFACTORING.md 2.11); results_plotter re-exports every name.
"""
from __future__ import annotations

from pathlib import Path
from collections import defaultdict

import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr

from plot_common import COLORS_HYDRO, DEFAULT_ET_VAR, _infer_time_dim, _node_id, add_fire_vline, date_window, maybe_show, save_fig


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
    et_var:      str  = DEFAULT_ET_VAR,
) -> None:
    """Plot soil water content + actual ET time series.

    et_var : name of the ET variable inside data["ET_xr_mean"] — e.g.
        "ACT. ETRA" (raw) or "ACT. ETRA_patched" (artefact-floored).
        Defaults to DEFAULT_ET_VAR; pass explicitly to match whichever
        ET_VAR the caller used upstream.
    """
    colors_dict = colors_dict or COLORS_HYDRO
    d0, d1 = date_window(date_start, date_end, data["df_dates_sw"])
    print(f"  Plotting ET+sw time-series ({d0.date()} → {d1.date()}) …")

    ET_xr_mean = data["ET_xr_mean"]
    if et_var not in ET_xr_mean:
        raise KeyError(
            f"plot_et_sw: '{et_var}' not found in ET_xr_mean. "
            f"Available: {list(ET_xr_mean.data_vars)}"
        )

    df_sw       = data["df_sw"]
    df_dates_sw = data["df_dates_sw"]
    ET_data     = ET_xr_mean[et_var].values
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

    mask_pev = (pev_series.index >= d0) & (pev_series.index <= d1)
    mask_tp  = (tp_series.index  >= d0) & (tp_series.index  <= d1)
    fig, ax = plt.subplots(figsize=(14, 3))
    ax.bar(pev_series[mask_pev].index, pev_series[mask_pev].values,
           width=1.0, color="grey", label="ETp", alpha=0.5)
    ax.bar(tp_series[mask_tp].index, tp_series[mask_tp].values,
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
    dpi:    int  = 300,
    show:   bool = False,
    et_var: str  = DEFAULT_ET_VAR,
) -> None:
    """Catchment ETa aggregated at *freq* frequency ('ME' = monthly, 'YE' = annual).

    et_var : name of the ET variable inside ET_xr_mean — e.g. "ACT. ETRA"
        (raw) or "ACT. ETRA_patched" (artefact-floored). Defaults to
        DEFAULT_ET_VAR; pass explicitly to match whichever ET_VAR the
        caller used upstream.
    """
    if et_var not in ET_xr_mean:
        raise KeyError(
            f"plot_catchment_et_freq: '{et_var}' not found in ET_xr_mean. "
            f"Available: {list(ET_xr_mean.data_vars)}"
        )

    # 1. Catchment area
    res_x, res_y    = ET_xr.rio.resolution()
    pixel_area       = abs(float(res_x) * float(res_y))
    variable_name    = list(ET_xr.data_vars)[0]
    _t               = _infer_time_dim(ET_xr)
    valid_pixels     = int(ET_xr[variable_name].notnull().isel({_t: 0}).sum().item())
    area_m2          = valid_pixels * pixel_area

    # 2. Resample
    df_et = pd.DataFrame(
        {"ETa_mm": ET_xr_mean[et_var].values},
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
    dpi:    int  = 300,
    show:   bool = False,
    et_var: str  = DEFAULT_ET_VAR,
) -> None:
    """Monthly + annual catchment ET summaries (wrapper).

    et_var : forwarded to plot_catchment_et_freq — name of the ET
        variable inside data["ET_xr_mean"] (e.g. "ACT. ETRA" or
        "ACT. ETRA_patched").
    """
    print("  Computing monthly catchment ET …")
    plot_catchment_et_freq(
        data["ET_xr"], data["ET_xr_mean"], data["time_dates"],
        freq="ME", fig_path=fig_path, dpi=dpi, show=show, et_var=et_var,
    )
    print("  Computing annual catchment ET …")
    plot_catchment_et_freq(
        data["ET_xr"], data["ET_xr_mean"], data["time_dates"],
        freq="YE", fig_path=fig_path, dpi=dpi, show=show, et_var=et_var,
    )


def plot_hgraph(
    data: dict,
    fig_path: Path,
    *,
    dpi: int = 300,
    show: bool = False,
    fname: str = "hgraph.png",
) -> None:
    """Hydraulic graph from pyCATHY simu.show() — single-run only.

    Limited to whatever one CATHY run ``data["simu"]`` currently holds in
    memory (one calendar month, in the LAI-loop pipeline). Prefer
    ``plot_hgraph_csv`` when a multi-month ``hgraph_output.csv`` backup is
    available — this function remains for callers with only a live
    ``simu`` object and no such backup (e.g. LT_SShydro_results.py).
    """
    print("  Plotting hydraulic graph …")
    fig, ax = plt.subplots()

    data["simu"].show(prop="hgraph", ax=ax)

    # Extract y-data from all lines rendered by simu.show()
    y_values = np.concatenate([line.get_ydata() for line in ax.get_lines()])

    # Calculate 5th and 95th percentiles of the data
    y_min, y_max = np.percentile(y_values, [5, 95])

    # Apply percentile bounds to y-axis limits
    ax.set_ylim(y_min, y_max)

    save_fig(fig, fig_path, fname, dpi)
    maybe_show(fig, show)


def plot_hgraph_csv(
    df_hgraph:   pd.DataFrame,
    fig_path:    Path,
    *,
    dpi:         int  = 300,
    show:        bool = False,
    scenario_id: int | str | None = None,
    date_start:  str | pd.Timestamp | None = None,
    date_end:    str | pd.Timestamp | None = None,
    fname:       str  = "hgraph.png",
) -> None:
    """Full-record streamflow hydrograph from a pre-saved hgraph_output.csv.

    This is the multi-month counterpart of ``plot_hgraph``: instead of
    ``simu.show(prop="hgraph")`` — which only ever has access to whichever
    single month's run is currently loaded into ``simu`` — this reads the
    ``all_hgraph_recs`` backup written every month by
    ``Agramon_withLAI_withETp.py`` (``hgraph_output.csv``), so the plotted
    hydrograph spans the *entire* simulation instead of one month.

    Parameters
    ----------
    df_hgraph : pd.DataFrame
        Loaded from ``hgraph_output.csv`` — one row per output step, with
        at least ``datetime`` and ``streamflow`` columns (as written by
        ``read_hgraph`` + ``add_datetime_col`` in the LAI-loop pipeline;
        the ``Unnamed0/1/2`` columns from the raw CATHY hgraph file are
        ignored here).
    fig_path    : Path to the output figure directory.
    dpi         : Figure resolution.
    show        : If True, call ``plt.show()`` after saving.
    scenario_id : Label for the figure suptitle.
    date_start / date_end : optional window (see ``date_window``); by
        default the full backed-up record is plotted.
    """
    print("  Plotting hydrograph from hgraph_output.csv …")
    scenario_id = scenario_id if scenario_id is not None else "?"

    if df_hgraph.empty:
        print("  WARNING: hgraph_output.csv is empty — nothing to plot.")
        return

    # ── Clean up: consecutive monthly runs can each re-emit a step at
    # their shared boundary (this month's last write-out timestep vs.
    # next month's first) — keep the later month's value and sort by
    # absolute time so the concatenated record reads as one continuous
    # series instead of a sawtooth at every month edge.
    df = (
        df_hgraph
        .dropna(subset=["datetime", "streamflow"])
        .sort_values("datetime")
        .drop_duplicates(subset="datetime", keep="last")
        .reset_index(drop=True)
    )
    if df.empty:
        print("  WARNING: hgraph_output.csv has no valid datetime/streamflow rows — skipping.")
        return

    d0, d1 = date_window(date_start, date_end, df["datetime"])
    df = df[(df["datetime"] >= d0) & (df["datetime"] <= d1)]
    if df.empty:
        print(f"  WARNING: no hgraph data in window {d0.date()} – {d1.date()} — skipping.")
        return

    q = pd.Series(df["streamflow"].values, index=pd.DatetimeIndex(df["datetime"]))

    # Robust y-limits: a handful of solver-transient spikes shouldn't
    # flatten the rest of the record, but don't clip so hard that real
    # peak-flow events (the whole point of a hydrograph) disappear.
    y_min, y_max = np.nanpercentile(q.values, [0.5, 99.5])
    pad = 0.05 * (y_max - y_min) if y_max > y_min else 1.0

    # ── Full-record hydrograph + cumulative discharge ─────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=False)

    axes[0].plot(q.index, q.values, color="#0072B2", linewidth=0.8)
    axes[0].set_ylim(y_min - pad, y_max + pad)
    axes[0].set_title("Streamflow at the outlet — full record")
    axes[0].set_ylabel("Streamflow")
    axes[0].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    add_fire_vline(axes[0], d0, d1)
    axes[0].legend(fontsize=8)
    axes[0].grid(True, linestyle="--", alpha=0.4)

    axes[1].plot(q.index, q.cumsum().values, color="#0072B2")
    axes[1].set_title("Cumulative discharge")
    axes[1].set_ylabel("Cumulative streamflow")
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    add_fire_vline(axes[1], d0, d1)
    axes[1].legend(fontsize=8)
    axes[1].grid(True, linestyle="--", alpha=0.4)

    fig.suptitle(f"Scenario: {scenario_id} — Streamflow hydrograph "
                 f"({d0.date()} – {d1.date()})", fontweight="bold")
    fig.tight_layout()
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
