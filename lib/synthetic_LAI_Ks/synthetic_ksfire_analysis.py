"""
synthetic_ksfire_analysis.py — area-mean diagnostics and scenario comparison
=============================================================================
Companion of synthetic_ksfire_LAI.py.

(1) PER RUN (called by synthetic_ksfire_LAI.py after the simulation)
    * every quantity is the MEAN OVER THE FULL AREA (all 20x20 surface
      columns), not the value at the outlet:
        - ETa and recharge: mean over X, Y of the surface fields;
        - Sw and psi at 0, 1 and 2 m depth: mean over the nodes found at that
          depth below every surface node.
    * writes  outputs/areamean_timeseries.csv  (daily),
              outputs/areamean_monthly.csv, areamean_periods.csv, run_meta.json
    * figures: fig_areamean_hydro.png, fig_areamean_year_compare.png

(2) ACROSS SCENARIOS (standalone):
        python synthetic_ksfire_analysis.py --root ../synthetic_LAI_Ks
    reads  <root>/<name>/<scenario>/outputs/areamean_timeseries.csv  of every
    run and writes <root>/_comparison/ :
        fig_cmp_monthly_eta_recharge.png   monthly ETa / recharge, all scenarios
        fig_cmp_cumulative.png             cumulative ETa / recharge
        fig_cmp_delta_vs_control.png       scenario - matching control (monthly, cumulative)
        fig_cmp_periods.png                half-year totals (before/after fire) + delta vs control
        fig_cmp_ks_response.png            delta ETa / recharge vs the Ks applied that month
        fig_cmp_soil_water.png             Sw and psi at 0/1/2 m, all scenarios
        summary_periods.csv, summary_annual.csv
    The "matching control" of a run is the mode-off run with the same number of
    layers (control for 1 layer, control_2L for 2 layers). Because the weather
    is identical, scenario - control isolates the effect of Ks.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

DEPTHS = (0.0, 1.0, 2.0)
FLUX_COLS = ("eta_mm", "recharge_mm", "rain_mm", "etp_mm")


# ══ 1. Area means ═════════════════════════════════════════════════════════════

def area_node_sets(grid3d, depths=DEPTHS) -> dict:
    """{depth: {"idx": node indices (0-based) at that depth under EVERY surface
    node, "actual": mean depth really found (m)}}. If the mesh is shallower than
    a requested depth, the deepest node is used and `actual` shows it."""
    from scipy.spatial import cKDTree
    nodes = np.asarray(grid3d["mesh3d_nodes"], float)
    nnod = int(grid3d["nnod"])
    surf = nodes[:nnod]
    tree = cKDTree(nodes)
    out = {}
    for d in depths:
        pts = np.column_stack([surf[:, 0], surf[:, 1], surf[:, 2] - d])
        _, idx = tree.query(pts)
        actual = float(np.mean(surf[:, 2] - nodes[idx, 2]))
        out[float(d)] = {"idx": np.unique(idx), "actual": actual}
        if abs(actual - d) > max(0.25 * d, 0.1):
            print(f"  WARNING: requested {d:g} m but the mesh gives {actual:.2f} m on average")
    return out


def _mean_over_nodes(da, idx, n_total) -> pd.Series:
    """Mean of a (datetime, node) DataArray over the node indices `idx`.
    Handles 0-based and 1-based node labels."""
    nc = np.asarray(da["node"].values)
    lab = idx + 1 if (nc.min() == 1 and nc.max() == n_total) else idx
    lab = lab[np.isin(lab, nc)]
    if lab.size == 0:
        raise ValueError("no node of the requested depth found in the output")
    return da.sel(node=lab).mean("node").to_series()


def _daily(s: pd.Series) -> pd.Series:
    s = s.copy()
    s.index = pd.DatetimeIndex(s.index)
    return s.sort_index().resample("D").mean()


def pick_recharge_var(ds, name=None) -> str:
    """Recharge variable of the dataset: --rech-var, else a name containing
    'rech', else the first data variable that is not a coordinate/time column."""
    if name:
        return name
    cands = [v for v in ds.data_vars if "rech" in v.lower()]
    if cands:
        return cands[0]
    skip = ("time", "datetime", "x", "y")
    cands = [v for v in ds.data_vars if v.lower() not in skip]
    if not cands:
        raise ValueError("no usable variable in the recharge dataset")
    return cands[0]


def _area_mean_surface(ds, var, scale) -> pd.Series:
    da = ds[var] * scale
    dims = [d for d in da.dims if d != "datetime"]
    return _daily(da.mean(dim=dims, skipna=True).to_series())


def build_timeseries(psi_xr, sw_xr, et_xr, rech_xr, grid3d, *, et_var, et_scale,
                     rech_var=None, rech_scale=None, rain_mmd=None, etp_mmd=None,
                     lai_hist=None, ks_df=None, depths=DEPTHS):
    """Daily full-area-mean table + info dict (depths really used, recharge var)."""
    cols, info = {}, {}
    cols["eta_mm_d"] = _area_mean_surface(et_xr, et_var, et_scale)
    if rech_xr is not None:
        rv = pick_recharge_var(rech_xr, rech_var)
        info["rech_var"] = rv
        cols["recharge_mm_d"] = _area_mean_surface(
            rech_xr, rv, et_scale if rech_scale is None else rech_scale)
    n_total = len(np.asarray(grid3d["mesh3d_nodes"]))
    sets = area_node_sets(grid3d, depths)
    info["depth_actual_m"] = {f"{d:g}": v["actual"] for d, v in sets.items()}
    for d, v in sets.items():
        cols[f"sw_z{d:g}m"] = _daily(_mean_over_nodes(sw_xr["sw"], v["idx"], n_total))
        cols[f"psi_z{d:g}m"] = _daily(_mean_over_nodes(psi_xr["psi"], v["idx"], n_total))
    ts = pd.concat(cols, axis=1).sort_index()
    ix = ts.index
    if rain_mmd is not None:
        ts["rain_mm_d"] = rain_mmd.reindex(ix)
    if etp_mmd is not None:
        ts["etp_mm_d"] = etp_mmd.reindex(ix)
    if lai_hist:
        lai = pd.Series({pd.Timestamp(k): float(np.nanmean(v)) for k, v in lai_hist.items()}).sort_index()
        ts["lai_mean"] = lai.reindex(ix, method="ffill")
    if ks_df is not None and len(ks_df):
        k = ks_df.set_index(pd.to_datetime(ks_df["month"])).sort_index()
        for c in ("permx_mm_hr", "permx_deep_mm_hr", "lai_driver"):
            if c in k:
                ts[c] = k[c].astype(float).reindex(ix, method="ffill")
    return ts, info


def monthly_table(ts: pd.DataFrame) -> pd.DataFrame:
    """Monthly means (mm/d) and totals (mm) of the daily full-area means."""
    m = ts.resample("MS").mean()
    dim = m.index.days_in_month.to_numpy()
    for c in ("eta_mm_d", "recharge_mm_d", "rain_mm_d", "etp_mm_d"):
        if c in m:
            m[c.replace("_mm_d", "_mm")] = m[c] * dim
    return m


def period_table(mon: pd.DataFrame, fire_date=None, freq: str = "half") -> pd.DataFrame:
    """Totals (mm) per half-year (Jan-Jun / Jul-Dec) or calendar year; only
    complete periods. 'phase' = pre / post / fire (period straddles the fire)."""
    cols = [c for c in FLUX_COLS if c in mon]
    n_need = 6 if freq == "half" else 12
    if freq == "half":
        key = pd.Series([f"{i.year}-H{1 if i.month <= 6 else 2}" for i in mon.index],
                        index=mon.index)
    else:
        key = pd.Series([str(i.year) for i in mon.index], index=mon.index)
    g = mon[cols].groupby(key)
    out = g.sum(min_count=1)
    out["n_months"] = g.size()
    out = out[out["n_months"] == n_need].drop(columns="n_months")
    starts, phases = [], []
    for k in out.index:
        y = int(k[:4])
        s = pd.Timestamp(y, 1 if (freq == "year" or k.endswith("H1")) else 7, 1)
        e = s + pd.DateOffset(months=n_need)
        starts.append(s)
        if fire_date is None:
            phases.append("")
        else:
            f = pd.Timestamp(fire_date)
            phases.append("pre" if e <= f else "post" if s >= f else "fire")
    out["start"], out["phase"] = starts, phases
    if "eta_mm" in out and "etp_mm" in out:
        out["eta_over_etp"] = out["eta_mm"] / out["etp_mm"]
    if "recharge_mm" in out and "rain_mm" in out:
        out["recharge_over_rain"] = out["recharge_mm"] / out["rain_mm"]
    return out


# ══ 2. Per-run outputs ════════════════════════════════════════════════════════

def _plt():
    import matplotlib
    import matplotlib.pyplot as plt
    return plt


def _fire_line(ax, fire_date, label=None):
    if fire_date is not None:
        ax.axvline(pd.Timestamp(fire_date), color="red", ls="--", lw=1.2, alpha=0.8,
                   label=label)


def plot_run_hydro(ts, fig_dir, fire_date, title, info, dpi=150):
    plt = _plt()
    z = info.get("depth_actual_m", {})
    zl = lambda d: f"{d} m" if abs(z.get(d, float(d)) - float(d)) < 0.05 else \
        f"{d} m (mesh: {z[d]:.2f} m)"
    dz = [c[4:-1] for c in ts if c.startswith("sw_z")]
    fig, ax = plt.subplots(6, 1, figsize=(11, 14), sharex=True, constrained_layout=True)
    a = ax[0]
    if "rain_mm_d" in ts:
        a.bar(ts.index, ts["rain_mm_d"], width=1.0, color="#4575b4", alpha=0.6, label="rain")
    if "etp_mm_d" in ts:
        a.plot(ts.index, ts["etp_mm_d"], color="#d94801", lw=0.8, label="ETp")
    a.set_ylabel("mm d⁻¹")
    a.set_title(f"{title} — full-area means", fontsize=10, fontweight="bold")
    a = ax[1]
    a.plot(ts.index, ts["eta_mm_d"], color="#fdae61", lw=0.6, label="ETa daily")
    a.plot(ts.index, ts["eta_mm_d"].rolling(15, center=True, min_periods=1).mean(),
           color="#b2182b", lw=1.6, label="ETa 15-d mean")
    if "etp_mm_d" in ts:
        a.plot(ts.index, ts["etp_mm_d"], color="grey", lw=0.7, ls=":", label="ETp")
    a.set_ylabel("ETa (mm d⁻¹)")
    a = ax[2]
    if "recharge_mm_d" in ts:
        a.plot(ts.index, ts["recharge_mm_d"], color="#9ecae1", lw=0.6, label="recharge daily")
        a.plot(ts.index, ts["recharge_mm_d"].rolling(15, center=True, min_periods=1).mean(),
               color="#08519c", lw=1.6, label="recharge 15-d mean")
    a.set_ylabel("recharge (mm d⁻¹)")
    cols = ["#1b7837", "#762a83", "#b35806", "#555555"]
    for i, d in enumerate(dz):
        if f"sw_z{d}m" in ts:
            ax[3].plot(ts.index, ts[f"sw_z{d}m"], color=cols[i % 4], lw=1.2, label=zl(d))
        if f"psi_z{d}m" in ts:
            ax[4].plot(ts.index, ts[f"psi_z{d}m"], color=cols[i % 4], lw=1.2, label=zl(d))
    ax[3].set_ylabel("Sw (–)  area mean")
    ax[4].set_ylabel("ψ (m)  area mean")
    a = ax[5]
    if "lai_mean" in ts:
        a.plot(ts.index, ts["lai_mean"], color="#006837", lw=1.5, label="LAI mean")
    a.set_ylabel("LAI (m² m⁻²)")
    if "permx_mm_hr" in ts:
        b = a.twinx()
        b.plot(ts.index, ts["permx_mm_hr"], color="#762a83", lw=1.3, label="Ks top")
        if "permx_deep_mm_hr" in ts and ts["permx_deep_mm_hr"].notna().any():
            b.plot(ts.index, ts["permx_deep_mm_hr"], color="#1b7837", ls="--", lw=1.0,
                   label="Ks deep")
        b.set_ylabel("PERMX (mm hr⁻¹)")
        h1, l1 = a.get_legend_handles_labels()
        h2, l2 = b.get_legend_handles_labels()
        a.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)
    for i, a in enumerate(ax):
        _fire_line(a, fire_date, "fire" if i == 0 else None)
        a.grid(True, ls="--", alpha=0.3)
        if i < 5:
            a.legend(loc="upper right", fontsize=8, ncol=3)
    fig.savefig(Path(fig_dir) / "fig_areamean_hydro.png", dpi=dpi)
    plt.close(fig)


def plot_run_year_compare(mon, per, fig_dir, fire_date, title, dpi=150):
    """Top: ETa / recharge per calendar month, one line per year (same season,
    different years; filled marker = after the fire). Bottom: half-year totals."""
    plt = _plt()
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    f = pd.Timestamp(fire_date) if fire_date is not None else None
    years = sorted(set(mon.index.year))
    cmap = plt.get_cmap("viridis")
    for j, (col, lab) in enumerate((("eta_mm", "ETa (mm month⁻¹)"),
                                    ("recharge_mm", "recharge (mm month⁻¹)"))):
        a = ax[0, j]
        if col not in mon:
            continue
        for k, y in enumerate(years):
            m = mon[mon.index.year == y]
            c = cmap(k / max(len(years) - 1, 1))
            a.plot(m.index.month, m[col], "-", color=c, lw=1.2, label=str(y))
            post = np.array([f is not None and d >= f for d in m.index])
            a.scatter(m.index.month[post], m[col][post], color=c, s=28, zorder=3)
            a.scatter(m.index.month[~post], m[col][~post], facecolors="none",
                      edgecolors=c, s=28, zorder=3)
        a.set_xticks(range(1, 13))
        a.set_xlabel("calendar month")
        a.set_ylabel(lab)
        a.grid(True, ls="--", alpha=0.3)
        if j == 0:
            a.legend(title="year (● after fire, ○ before)", fontsize=8)
    pc = {"pre": "#4575b4", "post": "#d94801", "fire": "#999999", "": "#4575b4"}
    for j, (col, lab) in enumerate((("eta_mm", "ETa (mm per half-year)"),
                                    ("recharge_mm", "recharge (mm per half-year)"))):
        a = ax[1, j]
        if col not in per or per.empty:
            continue
        a.bar(range(len(per)), per[col], color=[pc[p] for p in per["phase"]])
        a.set_xticks(range(len(per)))
        a.set_xticklabels([f"{k}\n{p}" for k, p in zip(per.index, per["phase"])], fontsize=8)
        a.set_ylabel(lab)
        a.grid(True, axis="y", ls="--", alpha=0.3)
    fig.suptitle(f"{title} — same season, different years", fontsize=10, fontweight="bold")
    fig.savefig(Path(fig_dir) / "fig_areamean_year_compare.png", dpi=dpi)
    plt.close(fig)


def save_run_analysis(out_dir, fig_dir, ts, info, meta, fire_date, plots=True, dpi=150):
    """Write the area-mean CSVs / metadata and (optionally) the per-run figures."""
    out_dir, fig_dir = Path(out_dir), Path(fig_dir)
    ts.to_csv(out_dir / "areamean_timeseries.csv", index_label="datetime")
    mon = monthly_table(ts)
    mon.to_csv(out_dir / "areamean_monthly.csv", index_label="month")
    per = period_table(mon, fire_date, "half")
    per.to_csv(out_dir / "areamean_periods.csv", index_label="period")
    meta = dict(meta, fire_date=str(fire_date), **info)
    (out_dir / "run_meta.json").write_text(json.dumps(meta, indent=2, default=str))
    if plots:
        title = meta.get("scenario", "run")
        plot_run_hydro(ts, fig_dir, fire_date, title, info, dpi)
        plot_run_year_compare(mon, per, fig_dir, fire_date, title, dpi)
    return mon, per


# ══ 3. Comparison across scenarios ════════════════════════════════════════════

def load_runs(root, exclude=()) -> dict:
    runs = {}
    for csv in sorted(Path(root).glob("*/*/outputs/areamean_timeseries.csv")):
        label = csv.parents[2].name
        if label.startswith("_") or label in exclude:
            continue
        if label in runs:
            label = f"{label}:{csv.parents[1].name}"
        ts = pd.read_csv(csv, index_col=0, parse_dates=True)
        mp = csv.parent / "run_meta.json"
        meta = json.loads(mp.read_text()) if mp.exists() else {}
        kp = csv.parent / "ks_fire_monthly.csv"
        ks = pd.read_csv(kp, parse_dates=["month"]).set_index("month") if kp.exists() else None
        if "two_layers" not in meta:
            meta["two_layers"] = label.endswith("2L")
        if "mode" not in meta:
            meta["mode"] = "off" if label.startswith("control") else "lai"
        runs[label] = {"ts": ts, "meta": meta, "ks": ks, "mon": monthly_table(ts)}
    return runs


def matching_control(runs: dict) -> dict:
    """label -> label of the mode-off run with the same number of layers (or None)."""
    ctrl = {}
    for lab, r in runs.items():
        if r["meta"]["mode"] == "off":
            ctrl.setdefault(bool(r["meta"]["two_layers"]), lab)
    return {lab: ctrl.get(bool(r["meta"]["two_layers"])) for lab, r in runs.items()}


def _style(runs):
    plt = _plt()
    cmap = plt.get_cmap("tab10")
    return {lab: {"color": cmap(i % 10), "ls": "--" if r["meta"]["two_layers"] else "-"}
            for i, (lab, r) in enumerate(runs.items())}


def _delta(runs, ctrl, lab, col):
    c = ctrl.get(lab)
    if c is None or c == lab or col not in runs[lab]["mon"] or col not in runs[c]["mon"]:
        return None
    return (runs[lab]["mon"][col] - runs[c]["mon"][col]).dropna()


def _longest(runs):
    return max(runs.values(), key=lambda r: len(r["ts"]))["ts"]


def cmp_monthly(runs, st, fire, out, dpi):
    plt = _plt()
    ref = _longest(runs)
    refm = monthly_table(ref)
    fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True, constrained_layout=True)
    for lab, r in runs.items():
        for a, col in zip(ax, ("eta_mm", "recharge_mm")):
            if col in r["mon"]:
                a.plot(r["mon"].index, r["mon"][col], marker="o", ms=3, lw=1.3, label=lab, **st[lab])
    if "etp_mm" in refm:
        ax[0].plot(refm.index, refm["etp_mm"], ":", color="grey", label="ETp")
    if "rain_mm" in refm:
        ax[1].bar(refm.index, refm["rain_mm"], width=20, color="#4575b4", alpha=0.2, label="rain")
    ax[0].set_ylabel("ETa (mm month⁻¹)")
    ax[1].set_ylabel("recharge (mm month⁻¹)")
    ax[0].set_title("Full-area monthly ETa and recharge — all scenarios", fontsize=10,
                    fontweight="bold")
    for a in ax:
        _fire_line(a, fire)
        a.grid(True, ls="--", alpha=0.3)
        a.legend(fontsize=7, ncol=3, loc="upper right")
    fig.savefig(out / "fig_cmp_monthly_eta_recharge.png", dpi=dpi)
    plt.close(fig)


def cmp_cumulative(runs, st, fire, out, dpi):
    plt = _plt()
    fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True, constrained_layout=True)
    for lab, r in runs.items():
        for a, col in zip(ax, ("eta_mm", "recharge_mm")):
            if col in r["mon"]:
                a.plot(r["mon"].index, r["mon"][col].cumsum(), lw=1.6, label=lab, **st[lab])
    ax[0].set_ylabel("cumulative ETa (mm)")
    ax[1].set_ylabel("cumulative recharge (mm)")
    ax[0].set_title("Cumulative full-area ETa and recharge since the start of the run",
                    fontsize=10, fontweight="bold")
    for a in ax:
        _fire_line(a, fire)
        a.grid(True, ls="--", alpha=0.3)
        a.legend(fontsize=7, ncol=3, loc="upper left")
    fig.savefig(out / "fig_cmp_cumulative.png", dpi=dpi)
    plt.close(fig)


def cmp_delta(runs, ctrl, st, fire, out, dpi):
    plt = _plt()
    fig, ax = plt.subplots(2, 2, figsize=(13, 8), sharex=True, constrained_layout=True)
    drawn = False
    for lab in runs:
        for i, col in enumerate(("eta_mm", "recharge_mm")):
            d = _delta(runs, ctrl, lab, col)
            if d is None:
                continue
            drawn = True
            ax[i, 0].plot(d.index, d.values, marker="o", ms=3, lw=1.2, label=lab, **st[lab])
            ax[i, 1].plot(d.index, d.cumsum().values, lw=1.6, label=lab, **st[lab])
    if not drawn:
        plt.close(fig)
        print("  no control run found -> delta figure skipped")
        return
    names = ("ETa", "recharge")
    for i in range(2):
        ax[i, 0].set_ylabel(f"Δ {names[i]} (mm month⁻¹)\nscenario − control")
        ax[i, 1].set_ylabel(f"cumulative Δ {names[i]} (mm)")
        for j in range(2):
            ax[i, j].axhline(0, color="k", lw=0.8)
            _fire_line(ax[i, j], fire)
            ax[i, j].grid(True, ls="--", alpha=0.3)
    ax[0, 0].set_title("monthly difference", fontsize=10)
    ax[0, 1].set_title("cumulative difference", fontsize=10)
    ax[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle("Effect of Ks on ETa and recharge (same weather: scenario − matching control)",
                 fontsize=10, fontweight="bold")
    fig.savefig(out / "fig_cmp_delta_vs_control.png", dpi=dpi)
    plt.close(fig)


def summary_tables(runs, ctrl, fire, out):
    res = {}
    for freq, name in (("half", "summary_periods.csv"), ("year", "summary_annual.csv")):
        rows, pt = [], {lab: period_table(r["mon"], fire, freq) for lab, r in runs.items()}
        for lab, p in pt.items():
            c = ctrl.get(lab)
            for k, row in p.iterrows():
                rec = {"scenario": lab, "period": k, "phase": row["phase"],
                       **{c_: row[c_] for c_ in p.columns
                          if c_ not in ("start", "phase")}}
                if c and c != lab and k in pt[c].index:
                    for col, tag in (("eta_mm", "eta"), ("recharge_mm", "rech")):
                        if col in p:
                            base = pt[c].loc[k, col]
                            rec[f"d_{tag}_mm"] = row[col] - base
                            rec[f"d_{tag}_pct"] = 100 * (row[col] - base) / base if base else np.nan
                rows.append(rec)
        df = pd.DataFrame(rows)
        df.to_csv(out / name, index=False)
        res[freq] = (df, pt)
    return res


def cmp_periods(res, runs, ctrl, st, out, dpi):
    plt = _plt()
    df, pt = res["half"]
    if df.empty:
        print("  no complete half-year -> period figure skipped")
        return
    phase = df.drop_duplicates("period").set_index("period")["phase"]
    order = sorted(phase.index)
    fig, ax = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    for j, (col, lab) in enumerate((("eta_mm", "ETa"), ("recharge_mm", "recharge"))):
        piv = df.pivot(index="period", columns="scenario", values=col).reindex(order)
        piv = piv[[c for c in runs if c in piv.columns]]
        piv.plot.bar(ax=ax[0, j], color=[st[c]["color"] for c in piv.columns], width=0.85)
        ax[0, j].set_ylabel(f"{lab} (mm per half-year)")
        ax[0, j].set_title(f"{lab}: before / after the fire", fontsize=10)
        dcol = "d_eta_mm" if col == "eta_mm" else "d_rech_mm"
        if dcol in df:
            pd_ = df.dropna(subset=[dcol]).pivot(index="period", columns="scenario",
                                                 values=dcol).reindex(order)
            pd_ = pd_[[c for c in runs if c in pd_.columns]]
            pd_.plot.bar(ax=ax[1, j], color=[st[c]["color"] for c in pd_.columns], width=0.85)
        ax[1, j].axhline(0, color="k", lw=0.8)
        ax[1, j].set_ylabel(f"Δ {lab} (mm) scenario − control")
        for a in (ax[0, j], ax[1, j]):
            a.set_xticklabels([f"{p}\n{phase[p]}" for p in order], rotation=0, fontsize=8)
            a.set_xlabel("")
            a.grid(True, axis="y", ls="--", alpha=0.3)
            a.legend(fontsize=7, ncol=2)
    fig.savefig(out / "fig_cmp_periods.png", dpi=dpi)
    plt.close(fig)


def cmp_ks_response(runs, ctrl, st, fire, out, dpi):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    f = pd.Timestamp(fire)
    sc = None
    for lab, r in runs.items():
        if r["ks"] is None:
            continue
        for a, col, nm in ((ax[0], "eta_mm", "ETa"), (ax[1], "recharge_mm", "recharge")):
            d = _delta(runs, ctrl, lab, col)
            if d is None:
                continue
            ks = r["ks"]["permx_mm_hr"].reindex(d.index)
            ok = ks.notna()
            mo = ((d.index.year - f.year) * 12 + d.index.month - f.month).to_numpy()
            sc = a.scatter(ks[ok], d[ok], c=mo[ok.to_numpy()], cmap="coolwarm", s=30,
                           marker="os^vDp"[list(runs).index(lab) % 6], edgecolors="k",
                           linewidths=0.4, vmin=-6, vmax=30, label=lab)
            a.set_xlabel("Ks top applied that month (mm hr⁻¹)")
            a.set_ylabel(f"Δ {nm} (mm month⁻¹)  scenario − control")
    for a in ax:
        a.axhline(0, color="k", lw=0.8)
        a.grid(True, ls="--", alpha=0.3)
        a.legend(fontsize=7)
    if sc is None:
        plt.close(fig)
        print("  no Ks/control pair -> Ks-response figure skipped")
        return
    fig.colorbar(sc, ax=ax, label="months since fire")
    fig.suptitle("Monthly ETa / recharge response vs Ks", fontsize=10, fontweight="bold")
    fig.savefig(out / "fig_cmp_ks_response.png", dpi=dpi)
    plt.close(fig)


def cmp_soil(runs, st, fire, out, dpi):
    plt = _plt()
    dz = sorted({c[4:-1] for r in runs.values() for c in r["ts"] if c.startswith("sw_z")},
                key=float)
    if not dz:
        return
    fig, ax = plt.subplots(len(dz), 2, figsize=(14, 3.2 * len(dz)), sharex=True,
                           constrained_layout=True, squeeze=False)
    for i, d in enumerate(dz):
        for j, (pre, lab) in enumerate((("sw", "Sw (–)"), ("psi", "ψ (m)"))):
            a = ax[i, j]
            for k, r in runs.items():
                c = f"{pre}_z{d}m"
                if c in r["ts"]:
                    s = r["ts"][c].rolling(7, center=True, min_periods=1).mean()
                    a.plot(s.index, s.values, lw=1.2, label=k, **st[k])
            a.set_ylabel(f"{lab}  @ {d} m")
            _fire_line(a, fire)
            a.grid(True, ls="--", alpha=0.3)
    ax[0, 0].legend(fontsize=7, ncol=2)
    ax[0, 0].set_title("full-area mean saturation (7-d mean)", fontsize=10)
    ax[0, 1].set_title("full-area mean pressure head (7-d mean)", fontsize=10)
    fig.savefig(out / "fig_cmp_soil_water.png", dpi=dpi)
    plt.close(fig)


def explain_missing_runs(root, runs) -> None:
    """Say exactly what is required and what each scenario folder is missing."""
    root = Path(root)
    print(f"\nCannot compare: need >= 2 runs, found {len(runs)} usable under\n  {root}")
    print("A run is usable only if this file exists (written by synthetic_ksfire_LAI.py "
          "at the end of a successful simulation):")
    print("  <root>/<scenario_name>/<synth_ks-..._YYYY-YYYY>/outputs/areamean_timeseries.csv")
    print("(optional, used if present: run_meta.json, ks_fire_monthly.csv)\n")
    if not root.is_dir():
        print(f"  !! root folder does not exist: {root}\n"
              f"     -> pass the right one with --root (default ../synthetic_LAI_Ks, "
              f"relative to the CURRENT directory)")
        return
    sims = sorted(d for d in root.glob("*/*") if d.is_dir() and not d.parent.name.startswith("_"))
    if not sims:
        print("  !! no <scenario_name>/<scenario>/ folder inside: nothing was simulated yet.\n"
              "     -> run ./run_synthetic_scenarios.sh first (OUT_ROOT must match --root)")
        return
    for sim in sims:
        out = sim / "outputs"
        label = f"{sim.parent.name}/{sim.name}"
        if (out / "areamean_timeseries.csv").exists():
            status = "OK (usable)"
        elif not out.is_dir():
            status = "no outputs/ folder -> the run died before saving anything"
        elif (out / "psi_output.nc").exists():
            status = ("simulated (psi/sw/et .nc exist) but areamean_timeseries.csv is missing "
                      "-> the area-mean step failed: see failures.csv ('area-mean analysis')")
        else:
            n = "?"
            fp = out / "failures.csv"
            if fp.exists():
                try:
                    n = max(len(pd.read_csv(fp)), 0)
                except Exception:
                    pass
            status = (f"no simulation output (psi_output.nc absent); {n} failure(s) logged in "
                      f"outputs/failures.csv" if fp.exists() else
                      "no simulation output and no failures.csv (run crashed early; "
                      "check its console/log)")
        print(f"  - {label}: {status}")
    print("\nFix the failed scenario(s) (or leave them out with --exclude <name>) and re-run.")


def compare_scenarios(root, out=None, exclude=(), fire_date=None, dpi=150) -> Path | None:
    import matplotlib
    matplotlib.use("Agg")
    runs = load_runs(root, exclude)
    if len(runs) < 2:
        explain_missing_runs(root, runs)
        return None
    out = Path(out) if out else Path(root) / "_comparison"
    out.mkdir(parents=True, exist_ok=True)
    fire = fire_date or next((r["meta"].get("fire_date") for r in runs.values()
                              if r["meta"].get("fire_date")), "2020-07-01")
    ctrl = matching_control(runs)
    st = _style(runs)
    print(f"Comparing {len(runs)} runs: {', '.join(runs)}")
    print("  matching control:", {k: v for k, v in ctrl.items()})
    res = summary_tables(runs, ctrl, fire, out)
    for fn in (lambda: cmp_monthly(runs, st, fire, out, dpi),
               lambda: cmp_cumulative(runs, st, fire, out, dpi),
               lambda: cmp_delta(runs, ctrl, st, fire, out, dpi),
               lambda: cmp_periods(res, runs, ctrl, st, out, dpi),
               lambda: cmp_ks_response(runs, ctrl, st, fire, out, dpi),
               lambda: cmp_soil(runs, st, fire, out, dpi)):
        try:
            fn()
        except Exception as exc:
            print(f"  WARNING: {type(exc).__name__}: {exc}")
    print(f"Comparison written to {out}")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description="Compare the synthetic Ks/fire scenarios")
    p.add_argument("--root", default="../synthetic_LAI_Ks",
                   help="folder holding <name>/<scenario>/outputs/ of every run")
    p.add_argument("--out", default=None, help="default: <root>/_comparison")
    p.add_argument("--exclude", nargs="*", default=[],
                   help="scenario names to leave out (e.g. deep_wt: other initial condition)")
    p.add_argument("--fire-date", default=None, help="default: read from run_meta.json")
    p.add_argument("--dpi", type=int, default=150)
    a = p.parse_args(argv)
    compare_scenarios(Path(a.root).expanduser().resolve(), a.out, a.exclude, a.fire_date, a.dpi)


if __name__ == "__main__":
    main()
