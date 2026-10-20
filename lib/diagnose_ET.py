"""
diagnose_ET.py
==============
GRwater project — Agramon study site

Diagnostics for why the ETa/ETp ratio is frequently > 1. Reuses the
loading/clipping helpers from analyse_ET.py so the data going into these
plots is identical to what analyse_ET.py itself uses.

Each hypothesis gets its own figure(s), saved under <fig_path>/Diagnostics:

    H1  Units / temporal-resolution mismatch
            plot_units_check()        histograms of raw ETa vs ETp + printed
                                       attrs (units) and time-step spacing
    H2  Spatial / masking misalignment
            plot_coverage_maps()      per-pixel valid-time-step count, ETa vs ETp
            plot_scatter_eta_etp()    pooled per-pixel-per-time scatter vs 1:1 line
    H3  Seasonality / advection ("oasis effect")
            plot_ratio_by_month()     boxplot of ratio grouped by calendar month
            plot_ratio_rolling()      raw + rolling-mean ratio time series
    H4  Per-plot / land-cover comparison
            plot_ratio_by_plot()      ratio time series + boxplot, one per PlotID
    H5  Data-quality / gap-filling artifacts
            plot_ratio_flagged()      ratio coloured by QC/fill flag, if one exists
            plot_ratio_vs_doy()       ratio vs day-of-year, all years overlaid
    H6  Annual-total sanity check
            plot_annual_totals()      paired bar chart of annual ETa vs ETp sums
    H7  Satellite vs Agramon CSV (dry / wet)
            plot_csv_comparison()     satellite ETa/ETp vs the CSV's ET,
                                       ET_0_dry and ET_0_wet (daily, monthly,
                                       annual totals, ratios), one set of
                                       figures per selected plot

A failure in any one diagnostic (e.g. no QC flag variable present) is caught
and reported, so one missing ingredient doesn't stop the rest of the run.

Usage
-----
    python diagnose_ET.py                     # default: plot index 0
    python diagnose_ET.py --plots all         # all plots
    python diagnose_ET.py --plots 0 2 4       # by index
    python diagnose_ET.py --plots Plot4Mulching
    python diagnose_ET.py --plots 0 --et-csv /path/agramon_treatments_ET.csv
    python diagnose_ET.py --plots 0 --csv-fid 3   # force CSV feature_id (single plot)

H7 extra options (ignored by analyse_ET.parse_args, parsed here):
    --et-csv PATH   agramon_treatments_ET.csv (default: <PathET>/agramon_treatments_ET.csv)
    --fire-date D   fire date (YYYY-MM-DD): vertical line on all time-series figures
    --csv-fid ID    CSV feature_id to use. Default: the shapefile row index of
                    each selected plot (same rule as the simulation script:
                    feature_id = dem_plot - 1).
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from analyse_ET import (
    fig_path,
    PathET,
    load_shapefiles,
    load_et,
    load_etp,
    clip_dataset,
    resolve_plot_selection,
    selection_label,
    parse_args,
)
from results_plotter import save_fig, maybe_show

RNG = np.random.default_rng(0)

# Date of the fire, drawn as a vertical line on every time-axis figure.
# Set it here (e.g. "2019-07-27") or pass --fire-date YYYY-MM-DD.
# None = no line.
FIRE_DATE: str | None = None


def _mark_fire(ax) -> None:
    """Vertical line on the fire date (only if it falls inside the x-range)."""
    if FIRE_DATE is None:
        return
    x = mdates.date2num(pd.Timestamp(FIRE_DATE))
    lo, hi = ax.get_xlim()
    if lo <= x <= hi:
        ax.axvline(x, color="crimson", ls="-.", lw=1.4, label=f"Fire ({FIRE_DATE})",
                   zorder=5)

# =============================================================================
# H1 — Units / temporal-resolution mismatch
# =============================================================================

def plot_units_check(
    ds_et: xr.Dataset, ds_etp: xr.Dataset,
    et_var: str, etp_var: str,
    out_path: Path, label: str, show: bool = False,
) -> None:
    """Histograms of raw ETa vs ETp + printed units/time-step spacing."""
    print("\n  [H1] Units / temporal-resolution check")
    for name, ds, var in (("ETa", ds_et, et_var), ("ETp", ds_etp, etp_var)):
        units = ds[var].attrs.get("units", "unknown")
        t = pd.to_datetime(ds["time"].values)
        steps = np.unique(np.diff(t)[:20]) if len(t) > 1 else []
        print(f"    {name}: units='{units}', n_time={len(t)}, "
              f"sample time-step(s)={list(steps)}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, name, ds, var in (
        (axes[0], "ETa", ds_et, et_var),
        (axes[1], "ETp", ds_etp, etp_var),
    ):
        vals = ds[var].values.astype(float).ravel()
        vals = vals[~np.isnan(vals)]
        ax.hist(vals, bins=60, color="tab:blue" if name == "ETa" else "tab:orange")
        ax.set_title(f"{name} raw distribution\n(units: {ds[var].attrs.get('units', '?')})")
        ax.set_xlabel(name)
    fig.suptitle(f"H1 — Units / distribution check — {label}")
    fig.tight_layout()
    save_fig(fig, out_path, f"H1_units_check_{label}.png", dpi=150)
    maybe_show(fig, show)


# =============================================================================
# H2 — Spatial / masking misalignment
# =============================================================================

def plot_coverage_maps(
    ds_et_crop: xr.Dataset, ds_etp_crop: xr.Dataset,
    et_var: str, etp_var: str,
    out_path: Path, label: str, show: bool = False,
) -> None:
    """Per-pixel valid-time-step count for ETa vs ETp, side by side."""
    valid_et = ds_et_crop[et_var].notnull().sum(dim="time")
    valid_etp = ds_etp_crop[etp_var].notnull().sum(dim="time")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    valid_et.plot(ax=axes[0], cmap="viridis")
    axes[0].set_title("ETa — valid time-step count")
    valid_etp.plot(ax=axes[1], cmap="viridis")
    axes[1].set_title("ETp — valid time-step count")
    fig.suptitle(f"H2 — Coverage / masking check — {label}")
    fig.tight_layout()
    save_fig(fig, out_path, f"H2_coverage_maps_{label}.png", dpi=150)
    maybe_show(fig, show)


def plot_scatter_eta_etp(
    ds_et_crop: xr.Dataset, ds_etp_crop: xr.Dataset,
    et_var: str, etp_var: str,
    out_path: Path, label: str, show: bool = False,
    max_points: int = 20000,
) -> None:
    """Pooled per-pixel-per-time scatter of ETa vs ETp against the 1:1 line."""
    eta_da, etp_da = xr.align(ds_et_crop[et_var], ds_etp_crop[etp_var], join="inner")
    eta = eta_da.values.ravel()
    etp = etp_da.values.ravel()
    mask = ~np.isnan(eta) & ~np.isnan(etp)
    eta, etp = eta[mask], etp[mask]

    if eta.size == 0:
        print("  [H2] No overlapping valid ETa/ETp pixels after alignment — skipping scatter.")
        return

    if eta.size > max_points:
        idx = RNG.choice(eta.size, size=max_points, replace=False)
        eta, etp = eta[idx], etp[idx]

    lims = [0, max(np.nanpercentile(eta, 99.5), np.nanpercentile(etp, 99.5))]

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(etp, eta, s=4, alpha=0.15, color="tab:blue")
    ax.plot(lims, lims, color="grey", ls="--", lw=1, label="1:1")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("ETp")
    ax.set_ylabel("ETa")
    ax.set_title(f"H2 — ETa vs ETp (pooled pixels) — {label}\n"
                 f"{100 * np.mean(eta > etp):.1f}% of points above 1:1 line")
    ax.legend()
    fig.tight_layout()
    save_fig(fig, out_path, f"H2_scatter_eta_vs_etp_{label}.png", dpi=150)
    maybe_show(fig, show)


# =============================================================================
# H3 — Seasonality / advection ("oasis effect")
# =============================================================================

def _spatial_mean_ratio(
    ds_et_crop: xr.Dataset, ds_etp_crop: xr.Dataset, et_var: str, etp_var: str,
) -> xr.DataArray:
    eta_ts = ds_et_crop[et_var].mean(dim=("x", "y"), skipna=True)
    etp_ts = ds_etp_crop[etp_var].mean(dim=("x", "y"), skipna=True)
    eta_ts, etp_ts = xr.align(eta_ts, etp_ts, join="inner")
    return eta_ts / etp_ts.where(etp_ts != 0)


def plot_ratio_by_month(
    ds_et_crop: xr.Dataset, ds_etp_crop: xr.Dataset,
    et_var: str, etp_var: str,
    out_path: Path, label: str, show: bool = False,
) -> None:
    """Boxplot of the ETa/ETp ratio grouped by calendar month, all years pooled."""
    ratio = _spatial_mean_ratio(ds_et_crop, ds_etp_crop, et_var, etp_var)
    if ratio.sizes.get("time", 0) == 0:
        print("  [H3] No overlapping time steps — skipping monthly boxplot.")
        return

    months = pd.to_datetime(ratio["time"].values).month
    data = [ratio.values[months == m] for m in range(1, 13)]
    data = [d[~np.isnan(d)] for d in data]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.boxplot(data, positions=range(1, 13), showfliers=False)
    ax.axhline(1.0, color="grey", ls="--", lw=1)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
    ax.set_ylabel("ETa / ETp")
    ax.set_title(f"H3 — Ratio by month (all years pooled) — {label}")
    fig.tight_layout()
    save_fig(fig, out_path, f"H3_ratio_by_month_{label}.png", dpi=150)
    maybe_show(fig, show)


def plot_ratio_rolling(
    ds_et_crop: xr.Dataset, ds_etp_crop: xr.Dataset,
    et_var: str, etp_var: str,
    out_path: Path, label: str, show: bool = False,
    window: int = 30,
) -> None:
    """Raw + rolling-mean ETa/ETp ratio time series."""
    ratio = _spatial_mean_ratio(ds_et_crop, ds_etp_crop, et_var, etp_var)
    if ratio.sizes.get("time", 0) == 0:
        print("  [H3] No overlapping time steps — skipping rolling time series.")
        return

    s = pd.Series(ratio.values, index=pd.to_datetime(ratio["time"].values))
    rolling = s.rolling(window=window, min_periods=max(1, window // 3), center=True).mean()

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(s.index, s.values, color="tab:green", alpha=0.3, lw=0.8, label="raw")
    ax.plot(rolling.index, rolling.values, color="tab:green", lw=1.6,
            label=f"{window}-step rolling mean")
    ax.axhline(1.0, color="grey", ls="--", lw=1)
    ax.set_ylabel("ETa / ETp")
    ax.set_xlabel("Date")
    ax.set_title(f"H3 — Ratio time series (raw + rolling mean) — {label}")
    _mark_fire(ax)
    ax.legend()
    fig.tight_layout()
    save_fig(fig, out_path, f"H3_ratio_rolling_{label}.png", dpi=150)
    maybe_show(fig, show)


# =============================================================================
# H4 — Per-plot / land-cover comparison
# =============================================================================

def plot_ratio_by_plot(
    ds_et: xr.Dataset, ds_etp: xr.Dataset, gdf,
    et_var: str, etp_var: str,
    out_path: Path, show: bool = False,
) -> None:
    """Ratio time series and boxplot, one line/box per subcatchment row."""
    plot_ids = list(gdf["PlotID"])
    name_counts = {pid: plot_ids.count(pid) for pid in plot_ids}
    ratios: dict[str, xr.DataArray] = {}

    for i, pid in enumerate(plot_ids):
        # PlotID may not be unique (e.g. every row sharing a layer name), so
        # select by row position and only fall back to the bare name in the
        # label when it's actually unique.
        key = f"[{i}] {pid}" if name_counts[pid] > 1 else pid
        try:
            et_clip = clip_dataset(ds_et, gdf, [i])
            etp_clip = clip_dataset(ds_etp, gdf, [i])
            ratios[key] = _spatial_mean_ratio(et_clip, etp_clip, et_var, etp_var)
        except Exception as exc:  # noqa: BLE001 — keep going for other plots
            print(f"  [H4] Skipping plot index {i} ('{pid}'): {exc}")

    if not ratios:
        print("  [H4] No plots could be processed — skipping per-plot comparison.")
        return

    fig, ax = plt.subplots(figsize=(11, 4.5))
    for pid, ratio in ratios.items():
        s = pd.Series(ratio.values, index=pd.to_datetime(ratio["time"].values))
        ax.plot(s.index, s.rolling(30, min_periods=5, center=True).mean(),
                lw=1.2, label=pid)
    ax.axhline(1.0, color="grey", ls="--", lw=1)
    ax.set_ylabel("ETa / ETp (30-step rolling mean)")
    ax.set_title("H4 — Ratio time series by plot")
    _mark_fire(ax)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    save_fig(fig, out_path, "H4_ratio_timeseries_by_plot.png", dpi=150)
    maybe_show(fig, show)

    fig2, ax2 = plt.subplots(figsize=(9, 4.5))
    labels = list(ratios.keys())
    data = [r.values[~np.isnan(r.values)] for r in ratios.values()]
    ax2.boxplot(data, labels=labels, showfliers=False)
    ax2.axhline(1.0, color="grey", ls="--", lw=1)
    ax2.set_ylabel("ETa / ETp")
    ax2.set_title("H4 — Ratio distribution by plot")
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right")
    fig2.tight_layout()
    save_fig(fig2, out_path, "H4_ratio_boxplot_by_plot.png", dpi=150)
    maybe_show(fig2, show)


# =============================================================================
# H5 — Data-quality / gap-filling artifacts
# =============================================================================

def plot_ratio_flagged(
    ds_et_crop: xr.Dataset, ds_etp_crop: xr.Dataset,
    et_var: str, etp_var: str,
    out_path: Path, label: str, show: bool = False,
) -> None:
    """
    Colour the ratio time series by a QC/gap-fill flag variable, if one exists
    in either dataset. Looks for any data_var whose name contains 'flag',
    'qc', or 'fill' (case-insensitive).
    """
    def _find_flag_var(ds: xr.Dataset) -> str | None:
        for v in ds.data_vars:
            if any(k in v.lower() for k in ("flag", "qc", "fill")):
                return v
        return None

    flag_var = _find_flag_var(ds_et_crop) or _find_flag_var(ds_etp_crop)
    if flag_var is None:
        print("  [H5] No QC/flag/fill variable found in ET or ETp datasets — "
              "skipping flagged-ratio plot. If a flag exists under a different "
              "name, point plot_ratio_flagged() at it directly.")
        return

    src_ds = ds_et_crop if flag_var in ds_et_crop.data_vars else ds_etp_crop
    flag_ts = src_ds[flag_var].mean(dim=("x", "y"), skipna=True)
    ratio = _spatial_mean_ratio(ds_et_crop, ds_etp_crop, et_var, etp_var)
    ratio, flag_ts = xr.align(ratio, flag_ts, join="inner")

    fig, ax = plt.subplots(figsize=(11, 4.5))
    sc = ax.scatter(pd.to_datetime(ratio["time"].values), ratio.values,
                     c=flag_ts.values, cmap="autumn_r", s=10)
    ax.axhline(1.0, color="grey", ls="--", lw=1)
    fig.colorbar(sc, ax=ax, label=flag_var)
    ax.set_ylabel("ETa / ETp")
    ax.set_title(f"H5 — Ratio coloured by '{flag_var}' — {label}")
    _mark_fire(ax)
    fig.tight_layout()
    save_fig(fig, out_path, f"H5_ratio_flagged_{label}.png", dpi=150)
    maybe_show(fig, show)


def plot_ratio_vs_doy(
    ds_et_crop: xr.Dataset, ds_etp_crop: xr.Dataset,
    et_var: str, etp_var: str,
    out_path: Path, label: str, show: bool = False,
) -> None:
    """Ratio vs day-of-year, all years overlaid and coloured by year."""
    ratio = _spatial_mean_ratio(ds_et_crop, ds_etp_crop, et_var, etp_var)
    if ratio.sizes.get("time", 0) == 0:
        print("  [H5] No overlapping time steps — skipping ratio-vs-DOY plot.")
        return

    t = pd.to_datetime(ratio["time"].values)
    doy = t.dayofyear
    year = t.year

    fig, ax = plt.subplots(figsize=(9, 4.5))
    sc = ax.scatter(doy, ratio.values, c=year, cmap="viridis", s=10)
    ax.axhline(1.0, color="grey", ls="--", lw=1)
    fig.colorbar(sc, ax=ax, label="Year")
    ax.set_xlabel("Day of year")
    ax.set_ylabel("ETa / ETp")
    ax.set_title(f"H5 — Ratio vs day-of-year, all years overlaid — {label}")
    fig.tight_layout()
    save_fig(fig, out_path, f"H5_ratio_vs_doy_{label}.png", dpi=150)
    maybe_show(fig, show)


# =============================================================================
# H6 — Annual-total sanity check
# =============================================================================

def plot_annual_totals(
    ds_et_crop: xr.Dataset, ds_etp_crop: xr.Dataset,
    et_var: str, etp_var: str,
    out_path: Path, label: str, show: bool = False,
) -> None:
    """Paired bar chart of annual (spatial-mean) ETa vs ETp totals."""
    eta_ts = ds_et_crop[et_var].mean(dim=("x", "y"), skipna=True)
    etp_ts = ds_etp_crop[etp_var].mean(dim=("x", "y"), skipna=True)
    eta_ts, etp_ts = xr.align(eta_ts, etp_ts, join="inner")

    eta_annual = eta_ts.resample(time="1YE").sum()
    etp_annual = etp_ts.resample(time="1YE").sum()
    years = pd.to_datetime(eta_annual["time"].values).year

    x = np.arange(len(years))
    width = 0.35

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(x - width / 2, eta_annual.values, width, label="ETa", color="tab:blue")
    ax.bar(x + width / 2, etp_annual.values, width, label="ETp", color="tab:orange")
    ax.set_xticks(x)
    ax.set_xticklabels(years)
    ax.set_ylabel("Annual total [mm]")
    ax.set_title(f"H6 — Annual ETa vs ETp totals — {label}")
    ax.legend()
    fig.tight_layout()
    save_fig(fig, out_path, f"H6_annual_totals_{label}.png", dpi=150)
    maybe_show(fig, show)

    for yr, a, p in zip(years, eta_annual.values, etp_annual.values):
        flag = " <-- ETa > ETp" if a > p else ""
        print(f"    {yr}: ETa={a:.1f}  ETp={p:.1f}  ratio={a / p:.2f}{flag}")


# =============================================================================
# H7 — Satellite ETa/ETp vs Agramon CSV (dry / wet)
# =============================================================================

CSV_COLS = {
    "ET":       ("CSV ET (real)",  "black"),
    "ET_0_dry": ("CSV ETp dry",    "tab:red"),
    "ET_0_wet": ("CSV ETp wet",    "tab:green"),
}


def _selection_positions(gdf, selection) -> list[int]:
    """Row positions of the selected plots (mirrors resolve_plot_selection)."""
    plot_ids = list(gdf["PlotID"])
    if selection is None or selection == ["all"]:
        return list(range(len(plot_ids)))
    pos: list[int] = []
    for item in selection:
        if isinstance(item, int) or (isinstance(item, str) and item.isdigit()):
            pos.append(int(item))
        else:
            pos.extend(i for i, pid in enumerate(plot_ids) if pid == item)
    return list(dict.fromkeys(pos))


def _parse_csv_args() -> tuple[Path, int | None]:
    """--et-csv / --csv-fid / --fire-date. Unknown flags (--plots, ...) are left alone."""
    global FIRE_DATE
    ap = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    ap.add_argument("--et-csv", default=str(PathET / "agramon_treatments_ET.csv"))
    ap.add_argument("--csv-fid", type=int, default=None)
    ap.add_argument("--fire-date", default=None, metavar="YYYY-MM-DD")
    ns, _ = ap.parse_known_args(sys.argv[1:])
    if ns.fire_date:
        FIRE_DATE = str(pd.Timestamp(ns.fire_date).date())
    return Path(ns.et_csv), ns.csv_fid


def plot_csv_comparison(
    ds_et: xr.Dataset, ds_etp: xr.Dataset, gdf,
    selection, et_var: str, etp_var: str,
    csv_path: Path, csv_fid: int | None,
    out_path: Path, show: bool = False,
) -> None:
    """Satellite ETa/ETp vs CSV ET / ET_0_dry / ET_0_wet, one set per plot."""
    if not Path(csv_path).is_file():
        print(f"  [H7] CSV not found: {csv_path} — skipping (use --et-csv).")
        return
    df = pd.read_csv(csv_path, parse_dates=["t"])
    df["t"] = df["t"].dt.normalize()
    avail = sorted(int(v) for v in df["feature_id"].unique())
    plot_ids = list(gdf["PlotID"])
    positions = _selection_positions(gdf, selection)
    if csv_fid is not None and len(positions) != 1:
        print("  [H7] --csv-fid only applies to a single selected plot; "
              "using row index for each plot instead.")
        csv_fid = None

    for pos in positions:
        pid = plot_ids[pos]
        fid = csv_fid if csv_fid is not None else pos
        if fid not in avail:
            print(f"  [H7] Plot [{pos}] {pid}: feature_id {fid} not in CSV "
                  f"{avail} — skipping (try --csv-fid).")
            continue
        print(f"\n  [H7] Plot [{pos}] {pid}  ↔  CSV feature_id {fid}")

        et_clip = clip_dataset(ds_et, gdf, [pos])
        etp_clip = clip_dataset(ds_etp, gdf, [pos])
        sat = pd.DataFrame({
            "ETa_sat": et_clip[et_var].mean(dim=("x", "y"), skipna=True).to_series(),
            "ETp_sat": etp_clip[etp_var].mean(dim=("x", "y"), skipna=True).to_series(),
        })
        sat.index = pd.to_datetime(sat.index).normalize()
        csv = (df[df["feature_id"] == fid].set_index("t")[list(CSV_COLS)]
               .sort_index())
        both = sat.join(csv, how="inner")
        if both.empty:
            print("    No overlapping dates between satellite data and CSV — skipping.")
            continue
        print(f"    Overlap: {both.index[0]:%Y-%m-%d} … {both.index[-1]:%Y-%m-%d} "
              f"({len(both)} days)")
        tag = f"{pid}_fid{fid}"

        # --- (a) daily series (7-d rolling) + (b) monthly means --------------
        fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=False)
        ax = axes[0]
        ax.plot(both.index, both["ETa_sat"].rolling(7, min_periods=2, center=True).mean(),
                color="tab:blue", lw=1.4, label="Satellite ETa")
        ax.plot(both.index, both["ETp_sat"].rolling(7, min_periods=2, center=True).mean(),
                color="tab:orange", lw=1.4, label="Satellite ETp")
        for col, (lab, colr) in CSV_COLS.items():
            ax.plot(both.index, both[col].rolling(7, min_periods=2, center=True).mean(),
                    color=colr, lw=1.0, ls="--" if col != "ET" else "-", label=lab)
        ax.set_ylabel("mm/day (7-day rolling mean)")
        ax.set_title(f"H7 — Satellite vs CSV (dry/wet) — {tag}")
        _mark_fire(ax)
        ax.legend(fontsize=8, ncol=3)

        ax = axes[1]
        monthly = both.resample("MS").mean()
        for col, lab, colr in (
            ("ETa_sat", "Sat ETa", "tab:blue"),
            ("ETp_sat", "Sat ETp", "tab:orange"),
            *[(c, v[0], v[1]) for c, v in CSV_COLS.items()],
        ):
            ax.plot(monthly.index, monthly[col], marker="o", ms=3, lw=1.0,
                    color=colr, label=lab)
        ax.set_ylabel("mm/day (monthly mean)")
        _mark_fire(ax)
        ax.legend(fontsize=8, ncol=3)
        fig.tight_layout()
        save_fig(fig, out_path, f"H7_sat_vs_csv_timeseries_{tag}.png", dpi=150)
        maybe_show(fig, show)

        # --- (c) annual totals + (d) scatter ETa_sat vs CSV ------------------
        annual = both.resample("YE").sum(min_count=1)
        annual.index = annual.index.year
        fig2, axes2 = plt.subplots(1, 2, figsize=(12, 4.8))
        annual.plot.bar(ax=axes2[0],
                        color=["tab:blue", "tab:orange", "black", "tab:red", "tab:green"])
        axes2[0].set_ylabel("Annual total [mm] (overlap days only)")
        axes2[0].set_xlabel("")
        axes2[0].set_title("Annual totals")
        axes2[0].legend(["Sat ETa", "Sat ETp", *[v[0] for v in CSV_COLS.values()]],
                        fontsize=7)
        ax = axes2[1]
        for col, (lab, colr) in CSV_COLS.items():
            ax.scatter(both[col], both["ETa_sat"], s=8, alpha=0.5, color=colr, label=lab)
        lim = np.nanmax(both[["ETa_sat", *CSV_COLS]].values)
        ax.plot([0, lim], [0, lim], color="grey", ls="--", lw=1)
        ax.set_xlabel("CSV value [mm/day]")
        ax.set_ylabel("Satellite ETa [mm/day]")
        ax.set_title("Satellite ETa vs CSV")
        ax.legend(fontsize=8)
        fig2.suptitle(f"H7 — {tag}")
        fig2.tight_layout()
        save_fig(fig2, out_path, f"H7_sat_vs_csv_annual_scatter_{tag}.png", dpi=150)
        maybe_show(fig2, show)

        # --- (e) ratios --------------------------------------------------------
        fig3, ax = plt.subplots(figsize=(11, 4.5))
        ratios = {
            "Sat ETa / Sat ETp":   both["ETa_sat"] / both["ETp_sat"].where(both["ETp_sat"] != 0),
            "Sat ETa / CSV dry":   both["ETa_sat"] / both["ET_0_dry"].where(both["ET_0_dry"] != 0),
            "Sat ETa / CSV wet":   both["ETa_sat"] / both["ET_0_wet"].where(both["ET_0_wet"] != 0),
            "CSV ET / CSV dry":    both["ET"] / both["ET_0_dry"].where(both["ET_0_dry"] != 0),
            "CSV ET / CSV wet":    both["ET"] / both["ET_0_wet"].where(both["ET_0_wet"] != 0),
        }
        for lab, r in ratios.items():
            ax.plot(r.index, r.rolling(30, min_periods=5, center=True).mean(), lw=1.2, label=lab)
        ax.axhline(1.0, color="grey", ls="--", lw=1)
        ax.set_ylabel("Ratio (30-day rolling mean)")
        ax.set_title(f"H7 — Ratios vs dry/wet ETp — {tag}")
        _mark_fire(ax)
        ax.legend(fontsize=8, ncol=2)
        fig3.tight_layout()
        save_fig(fig3, out_path, f"H7_sat_vs_csv_ratios_{tag}.png", dpi=150)
        maybe_show(fig3, show)

        print("    Annual totals [mm]:")
        print(annual.round(1).to_string().replace("\n", "\n      "))
        print("    Mean daily [mm/day]: "
              + "  ".join(f"{c}={both[c].mean():.2f}" for c in both.columns))


# =============================================================================
# Entry point
# =============================================================================

def main() -> None:
    selection, *_ = parse_args()  # remaining CLI values unused here
    csv_path, csv_fid = _parse_csv_args()

    print("\n── Loading shapefiles ──────────────────────────────────────────")
    gdf_proj = load_shapefiles()

    print("\n── Loading ET datasets ─────────────────────────────────────────")
    ds_et = load_et()
    ds_etp = load_etp()

    print("\n── Clipping ────────────────────────────────────────────────────")
    ds_et_crop = clip_dataset(ds_et, gdf_proj, selection)
    ds_etp_crop = clip_dataset(ds_etp, gdf_proj, selection)

    label = selection_label(gdf_proj, selection)
    et_var = next(iter(ds_et_crop.data_vars))
    etp_var = next(iter(ds_etp_crop.data_vars))

    out_path = fig_path / "Diagnostics"
    out_path.mkdir(parents=True, exist_ok=True)

    print(f"\n── Running diagnostics ─────────────────────────────────────────")
    print(f"  ET variable  : {et_var}")
    print(f"  ETp variable : {etp_var}")
    print(f"  Plot label   : {label}")
    print(f"  ET CSV       : {csv_path}")
    print(f"  Fire date    : {FIRE_DATE or 'not set (use --fire-date or edit FIRE_DATE)'}")
    print(f"  Output dir   : {out_path}\n")

    diagnostics = [
        ("H1 units check", plot_units_check,
         (ds_et_crop, ds_etp_crop, et_var, etp_var, out_path, label)),
        ("H2 coverage maps", plot_coverage_maps,
         (ds_et_crop, ds_etp_crop, et_var, etp_var, out_path, label)),
        ("H2 scatter", plot_scatter_eta_etp,
         (ds_et_crop, ds_etp_crop, et_var, etp_var, out_path, label)),
        ("H3 ratio by month", plot_ratio_by_month,
         (ds_et_crop, ds_etp_crop, et_var, etp_var, out_path, label)),
        ("H3 ratio rolling", plot_ratio_rolling,
         (ds_et_crop, ds_etp_crop, et_var, etp_var, out_path, label)),
        ("H4 ratio by plot", plot_ratio_by_plot,
         (ds_et, ds_etp, gdf_proj, et_var, etp_var, out_path)),
        ("H5 ratio flagged", plot_ratio_flagged,
         (ds_et_crop, ds_etp_crop, et_var, etp_var, out_path, label)),
        ("H5 ratio vs doy", plot_ratio_vs_doy,
         (ds_et_crop, ds_etp_crop, et_var, etp_var, out_path, label)),
        ("H6 annual totals", plot_annual_totals,
         (ds_et_crop, ds_etp_crop, et_var, etp_var, out_path, label)),
        ("H7 satellite vs CSV dry/wet", plot_csv_comparison,
         (ds_et, ds_etp, gdf_proj, selection, et_var, etp_var,
          csv_path, csv_fid, out_path)),
    ]

    for name, func, args in diagnostics:
        try:
            print(f"  → {name}")
            func(*args, show=False)
        except Exception:  # noqa: BLE001 — one failure shouldn't stop the rest
            print(f"    [FAILED] {name}")
            traceback.print_exc()

    print("\nDone.")


if __name__ == "__main__":
    main()
