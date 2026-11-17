"""
plot_hypotheses.py — validation and hypothesis-interpretation figures
=====================================================================
Pure matplotlib/pandas functions (headless, no CATHY): each takes tidy
pandas inputs and returns the ``Figure``; pass ``out`` to also save a PNG.

* ``plot_run_summary``     one-page health + response check of a run from its
                           ``monthly_stats.csv`` (ETa, SW, psi, recharge, solver flags)
* ``plot_validation_panel`` model vs observation: series, 1:1 scatter with
                           skill (bias, RMSE, r, NSE, KGE), monthly-mean bias
* ``plot_water_balance``   monthly P, ETa, recharge, runoff and the closure residual
* ``plot_recovery_panels`` burned vs control, ratio with tolerance band and the
                           recovery time (see ``recovery_metrics``)

Hypothesis studies (hyp1..hyp4, run_all_hypotheses.py; one study folder per hypothesis):

* ``plot_variant_matrix``      the runs of a hypothesis as a table (dry-run overview)
* ``plot_hypothesis_eta``      monthly WB ETa of every variant vs the EO-EB reference + seasonal cycle
* ``plot_hypothesis_skill``    bias / RMSE / KGE of every variant vs EO-EB (one bar per variant)
* ``plot_hypothesis_response`` ETa, soil moisture, recharge and psi of every variant (monthly_stats.csv)
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import recovery_metrics as RM

C_SIM, C_OBS, C_CTRL, C_BURN = "#0072B2", "#222222", "#009E73", "#D55E00"
C_FIRE, C_WARN = "#CC0000", "#E69F00"


def _finish(fig, out: Optional[str], dpi: int = 150):
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
    return fig


def _fire_line(ax, fire_date, label=True):
    if fire_date is not None:
        ax.axvline(pd.Timestamp(fire_date), color=C_FIRE, ls="--", lw=1,
                   label="fire" if label else None)


# ── skill ────────────────────────────────────────────────────
def skill_metrics(sim: pd.Series, obs: pd.Series) -> dict:
    """bias (sim-obs), RMSE, Pearson r, NSE, KGE(2009) on the common index (NaN dropped)."""
    c = pd.concat([sim.rename("s"), obs.rename("o")], axis=1).dropna()
    nan = dict(n=len(c), bias=np.nan, rmse=np.nan, r=np.nan, nse=np.nan, kge=np.nan)
    if len(c) < 2:
        return nan
    s, o = c["s"].to_numpy(float), c["o"].to_numpy(float)
    out = dict(n=len(c), bias=float(np.mean(s - o)), rmse=float(np.sqrt(np.mean((s - o) ** 2))))
    if np.std(s) > 0 and np.std(o) > 0:
        r = float(np.corrcoef(s, o)[0, 1])
        out["r"] = r
        out["kge"] = float(1 - np.sqrt((r - 1) ** 2 + (s.std() / o.std() - 1) ** 2
                                       + (s.mean() / o.mean() - 1) ** 2)) if o.mean() else np.nan
    else:
        out["r"] = out["kge"] = np.nan
    den = np.sum((o - o.mean()) ** 2)
    out["nse"] = float(1 - np.sum((s - o) ** 2) / den) if den > 0 else np.nan
    return out


# ── run summary ──────────────────────────────────────────────
def plot_run_summary(stats: pd.DataFrame, fire_date=None, title: str = "",
                     out: Optional[str] = None):
    """
    Four stacked panels from ``monthly_stats.csv`` (columns: month, eta_mean,
    eta_patched_mean, permx_mm_hr, sw_mean, psi_mean/min/max, recharge_mean, solver_status;
    missing columns are skipped). Months whose solver status is not OK are marked.
    """
    if "month" not in stats.columns:
        raise ValueError("stats needs a 'month' column")
    d = stats.copy()
    d["t"] = pd.to_datetime(d["month"].astype(str))
    d = d.sort_values("t")
    bad = (~d["solver_status"].astype(str).str.contains("OK")) if "solver_status" in d else \
        pd.Series(False, index=d.index)
    panels = [
        ("Ks (mm/hr)", [("permx_mm_hr", "PERMX applied", C_BURN)]),
        ("ETa (mm/day)", [("eta_mean", "ETa", C_SIM), ("eta_patched_mean", "ETa patched", C_BURN)]),
        ("Saturation (-)", [("sw_mean", "SW mean", C_CTRL)]),
        ("Pressure head (m)", [("psi_mean", "mean", C_OBS), ("psi_min", "min", C_BURN),
                               ("psi_max", "max", C_SIM)]),
        ("Recharge", [("recharge_mean", "recharge", C_SIM)]),
    ]
    panels = [(yl, [s for s in ser if s[0] in d.columns]) for yl, ser in panels]
    panels = [p for p in panels if p[1]]
    if not panels:
        raise ValueError("no plottable column in stats")
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 2.2 * len(panels) + 0.8),
                             sharex=True, squeeze=False)
    for ax, (ylab, series) in zip(axes[:, 0], panels):
        for col, lab, colr in series:
            ax.plot(d["t"], d[col], "-o", ms=3, color=colr, label=lab)
        if bad.any():
            first = series[0][0]
            ax.plot(d.loc[bad, "t"], d.loc[bad, first], "x", color=C_WARN, ms=9, mew=2,
                    label="solver not OK")
        _fire_line(ax, fire_date)
        if ylab.startswith("Ks"):
            ax.set_yscale("log")
        ax.set_ylabel(ylab); ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=4, loc="best")
    axes[0, 0].set_title(title or "Run summary (monthly)")
    fig.autofmt_xdate()
    return _finish(fig, out)


# ── validation ───────────────────────────────────────────────
def plot_validation_panel(sim: pd.Series, obs: pd.Series, label: str = "ETa",
                          units: str = "mm/day", out: Optional[str] = None):
    """Series, 1:1 scatter with skill text, and mean bias per calendar month."""
    c = pd.concat([sim.rename("sim"), obs.rename("obs")], axis=1).dropna()
    if c.empty:
        raise ValueError("sim and obs share no valid time steps")
    sk = skill_metrics(c["sim"], c["obs"])
    fig = plt.figure(figsize=(11, 6.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1])
    a1 = fig.add_subplot(gs[0, :]); a2 = fig.add_subplot(gs[1, 0]); a3 = fig.add_subplot(gs[1, 1])
    a1.plot(c.index, c["obs"], color=C_OBS, lw=1.2, label="observed")
    a1.plot(c.index, c["sim"], color=C_SIM, lw=1.2, label="simulated")
    a1.set_ylabel(f"{label} ({units})"); a1.legend(fontsize=8); a1.grid(alpha=0.3)
    lim = [min(c.min()), max(c.max())]
    a2.scatter(c["obs"], c["sim"], s=14, alpha=0.6, color=C_SIM)
    a2.plot(lim, lim, "k--", lw=1)
    a2.set_xlabel(f"observed ({units})"); a2.set_ylabel(f"simulated ({units})")
    a2.text(0.04, 0.96, "n={n}\nbias={bias:.3g}\nRMSE={rmse:.3g}\nr={r:.2f}\nNSE={nse:.2f}\n"
            "KGE={kge:.2f}".format(**sk), transform=a2.transAxes, va="top", fontsize=8,
            bbox=dict(fc="white", ec="0.7"))
    a2.grid(alpha=0.3)
    mb = (c["sim"] - c["obs"]).groupby(c.index.month).mean()
    a3.bar(mb.index, mb.values, color=np.where(mb.values >= 0, C_BURN, C_SIM))
    a3.axhline(0, color="k", lw=0.8)
    a3.set_xlabel("month"); a3.set_ylabel(f"mean bias sim-obs ({units})"); a3.set_xticks(range(1, 13))
    a3.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    return _finish(fig, out)


# ── water balance ────────────────────────────────────────────
def plot_water_balance(df: pd.DataFrame, out: Optional[str] = None, title: str = ""):
    """
    ``df`` indexed by month, columns in mm/month: ``precip`` and any of
    ``eta``, ``recharge``, ``runoff``. Bars: outputs (negative) vs precip
    (positive); line: closure residual ``precip - eta - recharge - runoff``
    (= storage change if the balance closes).
    """
    if "precip" not in df.columns:
        raise ValueError("df needs a 'precip' column")
    sinks = [c for c in ("eta", "recharge", "runoff") if c in df.columns]
    if not sinks:
        raise ValueError("df needs at least one of eta / recharge / runoff")
    colors = {"eta": C_CTRL, "recharge": C_SIM, "runoff": C_BURN}
    x = np.arange(len(df)); fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.bar(x, df["precip"], color="#9ecae1", label="precip")
    bottom = np.zeros(len(df))
    for c in sinks:
        v = df[c].fillna(0).to_numpy()
        ax.bar(x, -v, bottom=-bottom, color=colors[c], label=c); bottom += v
    resid = df["precip"] - df[sinks].fillna(0).sum(axis=1)
    ax.plot(x, resid, "k-o", ms=3, label="residual (P - outputs)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x[::max(1, len(x) // 12)])
    ax.set_xticklabels([pd.Timestamp(i).strftime("%Y-%m") for i in df.index[::max(1, len(x) // 12)]],
                       rotation=45, ha="right")
    ax.set_ylabel("mm / month"); ax.set_title(title or "Monthly water balance")
    ax.legend(fontsize=8, ncol=5); ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    return _finish(fig, out)


# ── recovery ─────────────────────────────────────────────────
def plot_recovery_panels(burned: pd.DataFrame, control: pd.DataFrame, fire_date,
                         tol: float = 0.1, sustain: int = 3, out: Optional[str] = None):
    """
    Per variable, left: burned vs control; right: ratio burned/control with
    the ``1 ± tol`` band, the fire date and the recovery date (start of the
    first ``sustain``-month run inside the band), titled with the recovery time.
    """
    cols = [c for c in burned.columns if c in control.columns]
    if not cols:
        raise ValueError("burned and control share no variable")
    fig, axes = plt.subplots(len(cols), 2, figsize=(11, 2.6 * len(cols)), squeeze=False)
    for i, v in enumerate(cols):
        ratio = RM.relative_to_control(burned[v], control[v])
        rec = RM.recovery_time(ratio, fire_date, tol, sustain)
        a, b = axes[i]
        a.plot(control.index, control[v], color=C_CTRL, label="control")
        a.plot(burned.index, burned[v], color=C_BURN, label="burned")
        _fire_line(a, fire_date); a.set_ylabel(v); a.legend(fontsize=7); a.grid(alpha=0.3)
        b.axhspan(1 - tol, 1 + tol, color=C_CTRL, alpha=0.15, label=f"1 ± {tol:g}")
        b.plot(ratio.index, ratio, "-o", ms=3, color=C_OBS)
        _fire_line(b, fire_date, label=False)
        if np.isfinite(rec):
            t0 = pd.Timestamp(fire_date) + pd.DateOffset(months=int(rec))
            hit = ratio.index[(ratio.index.year == t0.year) & (ratio.index.month == t0.month)]
            t_rec = hit[0] if len(hit) else t0       # the plotted sample of that month
            b.axvline(t_rec, color=C_CTRL, ls=":", lw=1.5)
            b.set_title(f"{v}: recovered after {int(rec)} months", fontsize=9)
        else:
            b.set_title(f"{v}: no sustained recovery in the series", fontsize=9)
        b.set_ylabel("burned / control"); b.grid(alpha=0.3); b.legend(fontsize=7)
    fig.autofmt_xdate(); fig.tight_layout()
    return _finish(fig, out)


# ── hypothesis studies ───────────────────────────────────────
_GROUP_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7A7A7A"]


def _colors(names, groups=None):
    """One colour per group (variants of a group share it, lightness varies); per variant if no groups."""
    groups = groups or {}
    order = list(dict.fromkeys(groups.get(n, "") or n for n in names))
    col = {g: _GROUP_COLORS[i % len(_GROUP_COLORS)] for i, g in enumerate(order)}
    return {n: col[groups.get(n, "") or n] for n in names}


def plot_variant_matrix(rows, title: str = "", out: Optional[str] = None):
    """Table figure of the runs of a hypothesis (rows = list of dicts: variant, group, varying options)."""
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(max(6, 1.6 * len(df.columns)), 0.38 * (len(df) + 2) + 0.6))
    ax.axis("off")
    cell = df.astype(str).values.tolist()
    tab = ax.table(cellText=cell, colLabels=list(df.columns), loc="center", cellLoc="left")
    tab.auto_set_font_size(False)
    tab.set_fontsize(8)
    tab.auto_set_column_width(list(range(len(df.columns))))
    cols = _colors(list(df["variant"]), dict(zip(df["variant"], df.get("group", [""] * len(df)))))
    for i, name in enumerate(df["variant"], start=1):
        tab[(i, 0)].set_facecolor(cols[name] + "33")
    for j in range(len(df.columns)):
        tab[(0, j)].set_text_props(fontweight="bold")
        tab[(0, j)].set_facecolor("#dddddd")
    ax.set_title(f"{title}  ({len(df)} runs)", fontsize=10, fontweight="bold")
    return _finish(fig, out)


def plot_hypothesis_eta(wb: dict, eo: dict, groups: Optional[dict] = None, title: str = "",
                        out: Optional[str] = None):
    """
    Monthly spatial-mean ETa (mm/day): every variant (``wb``: name -> Series) against the EO-EB
    reference(s) (``eo``: label -> Series, black), left as a time series, right as the mean
    seasonal cycle (calendar months).
    """
    col = _colors(list(wb), groups)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.2), gridspec_kw={"width_ratios": [2.2, 1]})
    for name, s in wb.items():
        a1.plot(s.index, s.values, color=col[name], lw=1.3, alpha=0.9, label=name)
        c = s.groupby(s.index.month).mean()
        a2.plot(c.index, c.values, color=col[name], lw=1.3, marker="o", ms=3)
    for label, s in eo.items():
        a1.plot(s.index, s.values, color=C_OBS, lw=2.0, ls="--", label=f"EO-EB ({label})")
        c = s.groupby(s.index.month).mean()
        a2.plot(c.index, c.values, color=C_OBS, lw=2.0, ls="--", marker="s", ms=3)
    a1.set_ylabel("ETa (mm/day)")
    a1.set_title("Monthly mean")
    a2.set_xticks(range(1, 13))
    a2.set_xlabel("month")
    a2.set_title("Seasonal cycle")
    a1.legend(fontsize=7, ncol=2, frameon=False)
    for a in (a1, a2):
        a.grid(alpha=0.25)
    fig.suptitle(title, fontsize=10, fontweight="bold")
    fig.tight_layout()
    return _finish(fig, out)


def plot_hypothesis_skill(metrics: pd.DataFrame, title: str = "", out: Optional[str] = None):
    """
    Skill of every variant vs EO-EB ETa from the metrics table (columns variant, bias, rmse, kge,
    optional group): three bar panels sorted by KGE (best on top); bias is sim-obs (mm/day).
    """
    d = metrics.sort_values("kge", ascending=True, na_position="first").reset_index(drop=True)
    col = _colors(list(d["variant"]), dict(zip(d["variant"], d["group"])) if "group" in d else None)
    fig, axes = plt.subplots(1, 3, figsize=(12, 0.34 * len(d) + 1.8), sharey=True)
    for ax, (m, lab) in zip(axes, (("kge", "KGE (-)"), ("rmse", "RMSE (mm/day)"), ("bias", "bias sim-obs (mm/day)"))):
        ax.barh(d["variant"], d[m], color=[col[v] for v in d["variant"]])
        ax.set_xlabel(lab)
        ax.grid(axis="x", alpha=0.25)
        if m in ("kge", "bias"):
            ax.axvline(0, color="0.4", lw=0.8)
    axes[0].tick_params(axis="y", labelsize=7)
    fig.suptitle(title, fontsize=10, fontweight="bold")
    fig.tight_layout()
    return _finish(fig, out)


def plot_hypothesis_response(stats: dict, groups: Optional[dict] = None, title: str = "",
                             out: Optional[str] = None):
    """
    Hydrological response of every variant from its monthly_stats.csv (``stats``: name -> DataFrame
    with month, eta_mean / eta_patched_mean, sw_mean, recharge_mean, psi_mean; missing columns skipped).
    """
    panels = [("eta_patched_mean", "eta_mean", "ETa (mm/day)"), ("sw_mean", None, "soil moisture (-)"),
              ("recharge_mean", None, "recharge"), ("psi_mean", None, "mean psi (m)")]
    col = _colors(list(stats), groups)
    avail = [p for p in panels if any(p[0] in d or (p[1] and p[1] in d) for d in stats.values())]
    if not avail:
        return None
    fig, axes = plt.subplots(len(avail), 1, figsize=(11, 2.3 * len(avail) + 0.8), sharex=True, squeeze=False)
    for ax, (c1, c2, lab) in zip(axes[:, 0], avail):
        for name, d in stats.items():
            c = c1 if c1 in d else c2
            if not c or "month" not in d:
                continue
            ax.plot(pd.to_datetime(d["month"].astype(str)), d[c], color=col[name], lw=1.2, label=name)
        ax.set_ylabel(lab, fontsize=8)
        ax.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=7, ncol=3, frameon=False)
    fig.suptitle(title, fontsize=10, fontweight="bold")
    fig.tight_layout()
    return _finish(fig, out)
