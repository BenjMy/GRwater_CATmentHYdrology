"""
hypothesis_tables.py - water-balance and recovery tables of a hypothesis study
===============================================================================
Pure pandas/numpy (no CATHY, no geopandas): everything is read from the files a finished run
already writes in its scenario folder, so a study can be summarised long after the runs.

Inputs of one run (``sim_dir``)
    monthly_stats.csv    rain_in_mm, eta_out_mm, recharge_out_mm (mm/month), sw_mean
                         (agramon_checks.check_monthly_sanity; flux-only, storage change NOT included)
    hgraph_output.csv    datetime, streamflow at the outlet (CATHY hgraph, assumed m3/s)

Monthly frame (``load_monthly_balance``), indexed by month start, columns
    precip_mm, eta_mm, recharge_mm      monthly sums, mm
    runoff_m3                           outlet volume = mean streamflow x seconds in the month
    runoff_mm                           runoff_m3 / catchment area (NaN unless the area is given)
    residual_mm                         precip - eta - recharge (what agramon_checks reports)
    residual_runoff_mm                  residual_mm - runoff_mm (NaN without runoff_mm)
    eta_mm_day, recharge_mm_day         monthly sums / days in the month
    sw_mean                             mean saturation (-)

Tables: ``annual_totals`` (one row per calendar year), ``period_totals`` (all / pre_fire /
post_fire, totals and mm/year), ``variant_summary`` (one row per variant: means and post-fire
mm/year, merged into summary_metrics.csv) and ``recovery_by_variant`` (burned/control ratio,
recovery time, deficit; needs a control run, see recovery_metrics.py).

Nothing here raises on missing or odd files: problems are returned as a list of ``notes`` so the
study report can show what was skipped and why.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import recovery_metrics as RM

FIRE_DATE_FALLBACK = "2020-07-01"
SUM_COLS = ("precip_mm", "eta_mm", "recharge_mm", "runoff_mm", "residual_mm", "residual_runoff_mm")
RECOVERY_VARS = ("eta_mm_day", "recharge_mm_day", "sw_mean")


def default_fire_date() -> str:
    """Fire date of the pipeline (agramon_ks.FIRE_DATE_DEFAULT), with a local fallback."""
    try:
        from agramon_ks import FIRE_DATE_DEFAULT
        return str(FIRE_DATE_DEFAULT)
    except Exception:
        return FIRE_DATE_FALLBACK


# ── one run -> monthly frame ─────────────────────────────────────────────────
def _month_start(values) -> pd.DatetimeIndex:
    return pd.DatetimeIndex(pd.to_datetime(pd.Series(values).astype(str), errors="coerce")
                            .dt.to_period("M").dt.to_timestamp())


def _runoff_volume(path: Path) -> pd.Series:
    """Monthly outlet volume (m3) from hgraph_output.csv: mean streamflow x seconds in the month."""
    hg = pd.read_csv(path, usecols=lambda c: c in ("datetime", "streamflow"))
    hg["datetime"] = pd.to_datetime(hg["datetime"], errors="coerce")
    hg["streamflow"] = pd.to_numeric(hg["streamflow"], errors="coerce")
    # consecutive monthly runs re-emit the step at their shared boundary: keep the later one
    hg = hg.dropna().sort_values("datetime").drop_duplicates("datetime", keep="last")
    q = hg.set_index("datetime")["streamflow"]
    mean_q = q.groupby(q.index.to_period("M")).mean()
    mean_q.index = mean_q.index.to_timestamp()
    return mean_q * mean_q.index.days_in_month * 86400.0


def load_monthly_balance(sim_dir, area_m2: Optional[float] = None) -> tuple:
    """Monthly water-balance frame of one run -> ``(DataFrame, notes)``; empty frame if unreadable."""
    sim_dir, notes = Path(sim_dir), []
    p = sim_dir / "monthly_stats.csv"
    if not p.exists():
        return pd.DataFrame(), [f"monthly_stats.csv not found in {sim_dir.name}"]
    st = pd.read_csv(p)
    if "month" not in st.columns:
        return pd.DataFrame(), ["monthly_stats.csv has no 'month' column"]
    st.index = _month_start(st["month"])
    st = st[st.index.notna()]
    st = st[~st.index.duplicated(keep="last")].sort_index()
    out = pd.DataFrame(index=st.index)
    out.index.name = "month"
    for dst, src in (("precip_mm", "rain_in_mm"), ("eta_mm", "eta_out_mm"),
                     ("recharge_mm", "recharge_out_mm"), ("sw_mean", "sw_mean")):
        if src in st.columns:
            out[dst] = pd.to_numeric(st[src], errors="coerce")
        else:
            out[dst] = np.nan
            notes.append(f"column {src} missing in monthly_stats.csv")
    out["runoff_m3"] = np.nan
    hp = sim_dir / "hgraph_output.csv"
    if hp.exists():
        try:
            out["runoff_m3"] = _runoff_volume(hp).reindex(out.index)
        except Exception as exc:
            notes.append(f"hgraph_output.csv unreadable ({exc}) - no runoff")
    else:
        notes.append("hgraph_output.csv not found - no runoff")
    out["runoff_mm"] = out["runoff_m3"] / area_m2 * 1000.0 if area_m2 else np.nan
    if not area_m2 and out["runoff_m3"].notna().any():
        notes.append("catchment_area_m2 not set: runoff is in m3 only (runoff_mm empty)")
    out["residual_mm"] = out["precip_mm"] - out["eta_mm"] - out["recharge_mm"]
    out["residual_runoff_mm"] = out["residual_mm"] - out["runoff_mm"]
    days = out.index.days_in_month
    out["eta_mm_day"] = out["eta_mm"] / days
    out["recharge_mm_day"] = out["recharge_mm"] / days
    return out, notes


# ── tables ───────────────────────────────────────────────────────────────────
def _closure_pct(residual, precip) -> float:
    return float(100.0 * residual / precip) if pd.notna(residual) and pd.notna(precip) and precip > 0 else np.nan


def annual_totals(monthly: pd.DataFrame) -> pd.DataFrame:
    """One row per calendar year: flux sums (mm), mean saturation, months present, closure %.
    Years with fewer than 12 months (run start / end) are flagged by ``n_months``."""
    if monthly is None or monthly.empty:
        return pd.DataFrame()
    g = monthly.groupby(monthly.index.year)
    out = g[[c for c in SUM_COLS if c in monthly]].sum(min_count=1)
    out["sw_mean"] = g["sw_mean"].mean()
    out["n_months"] = g.size()
    out["closure_pct"] = [_closure_pct(r, p) for r, p in zip(out["residual_mm"], out["precip_mm"])]
    out.index.name = "year"
    return out.reset_index()


def period_totals(monthly: pd.DataFrame, fire_date) -> pd.DataFrame:
    """Rows ``all`` / ``pre_fire`` / ``post_fire`` (a period with no month is left out).
    ``total_<flux>`` is the sum, ``annual_<flux>`` the mm/year rate over the months with data."""
    if monthly is None or monthly.empty:
        return pd.DataFrame()
    fire = pd.Timestamp(fire_date).to_period("M").to_timestamp()
    parts = {"all": monthly, "pre_fire": monthly[monthly.index < fire],
             "post_fire": monthly[monthly.index >= fire]}
    rows = []
    for name, d in parts.items():
        if d.empty:
            continue
        row = {"period": name, "n_months": len(d)}
        for c in SUM_COLS:
            n_ok = int(d[c].notna().sum())
            tot = d[c].sum(min_count=1)
            row[f"total_{c}"] = tot
            row[f"annual_{c}"] = float(tot / (n_ok / 12.0)) if n_ok else np.nan
        row["mean_eta_mm_day"] = d["eta_mm_day"].mean()
        row["mean_recharge_mm_day"] = d["recharge_mm_day"].mean()
        row["mean_sw"] = d["sw_mean"].mean()
        row["closure_pct"] = _closure_pct(row["total_residual_mm"], row["total_precip_mm"])
        rows.append(row)
    return pd.DataFrame(rows)


def tables_by_variant(monthly_by_variant: dict, fire_date, groups: Optional[dict] = None) -> dict:
    """Stack the per-variant tables: ``{'annual', 'period', 'summary'}`` DataFrames, ``variant`` first."""
    groups = groups or {}
    annual, period, periods = [], [], {}
    for name, m in monthly_by_variant.items():
        a, p = annual_totals(m), period_totals(m, fire_date)
        if len(a):
            annual.append(a.assign(variant=name))
        if len(p):
            period.append(p.assign(variant=name))
            periods[name] = p

    def stack(frames):
        if not frames:
            return pd.DataFrame()
        df = pd.concat(frames, ignore_index=True)
        return df[["variant"] + [c for c in df.columns if c != "variant"]]
    return {"annual": stack(annual), "period": stack(period),
            "summary": variant_summary(periods, groups)}


def variant_summary(period_by_variant: dict, groups: Optional[dict] = None) -> pd.DataFrame:
    """One row per variant: mean ETa / recharge / saturation over the record, and the same plus mm/year
    totals and closure % after the fire (``post_`` columns, NaN if the record has no post-fire month)."""
    groups = groups or {}
    rows = []
    for name, p in period_by_variant.items():
        row = {"variant": name, "group": groups.get(name, "")}
        for per, tag, full in (("all", "", False), ("post_fire", "post_", True)):
            r = p[p["period"] == per]
            r = r.iloc[0] if len(r) else None
            get = (lambda k: np.nan) if r is None else (lambda k: r[k])
            row[f"{tag}n_months"] = get("n_months")
            row[f"{tag}mean_eta_mm_day"] = get("mean_eta_mm_day")
            row[f"{tag}mean_recharge_mm_day"] = get("mean_recharge_mm_day")
            row[f"{tag}mean_sw"] = get("mean_sw")
            if full:
                row["post_eta_mm_yr"] = get("annual_eta_mm")
                row["post_recharge_mm_yr"] = get("annual_recharge_mm")
                row["post_runoff_mm_yr"] = get("annual_runoff_mm")
                row["post_closure_pct"] = get("closure_pct")
        rows.append(row)
    return pd.DataFrame(rows)


# ── recovery against a control catchment ─────────────────────────────────────
def recovery_by_variant(monthly_by_variant: dict, control: Optional[pd.DataFrame], fire_date,
                        tol: float = 0.1, sustain: int = 3, variables=RECOVERY_VARS) -> tuple:
    """
    burned/control recovery table of every variant -> ``(DataFrame, notes)``.
    Columns: variant, variable, recovery_months, min_ratio, final_ratio, deficit (recovery_metrics).
    Without a control frame (the control catchment is not delineated yet) the table is empty and the
    note says so; nothing else is skipped silently.
    """
    if control is None or control.empty:
        return pd.DataFrame(), ["recovery: no control catchment run configured "
                                "([validation] control_sim_dir) - skipped"]
    rows, notes = [], []
    for name, burned in monthly_by_variant.items():
        cols = [v for v in variables if v in burned and v in control
                and burned[v].notna().any() and control[v].notna().any()]
        if not cols:
            notes.append(f"recovery: {name} shares no usable variable with the control run")
            continue
        try:
            t = RM.recovery_table(burned, control, fire_date, tol, sustain, variables=cols)
        except Exception as exc:                      # one bad variant must not stop the study
            notes.append(f"recovery: {name} failed ({exc})")
            continue
        t = t.rename_axis("variable").reset_index()
        t.insert(0, "variant", name)
        rows.append(t)
    return (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()), notes
