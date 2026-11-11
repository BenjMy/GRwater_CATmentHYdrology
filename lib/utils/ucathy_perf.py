"""
ucathy_perf.py -- generic numerical-performance / QC evaluation for CATHY runs.

Works on anything that has a `read_outputs(name) -> DataFrame` method (a
pyCATHY `CATHY` object, or `FrameSource` wrapping DataFrames you loaded
yourself). Nothing is hard-wired to one script: column names, thresholds and
the set of checks are all configurable, and every check is guarded so a
missing column gives a note, never an exception.
"""
from __future__ import annotations

import contextlib
import time
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
# Logical name -> candidate column names, tried in order in every loaded frame.
# Extend this (never edit the checks) when a build uses different headers.
DEFAULT_ALIASES = {
    "dt":        ["Deltat", "DELTAT", "dt"],
    "nlin":      ["NLIN", "NL-a", "NL"],
    "lin":       ["LIN", "LINEAR", "LIN-IT"],          # unverified: check df.columns
    "back":      ["Back", "BACK"],
    "mbe":       ["MBE"],
    "rel_mbe":   ["REL. MBE", "REL_MBE", "RELMBE"],
    "atmact":    ["Atmact-vf", "ATMACT-VF"],
    "et_time":   ["time_sec", "TIME", "time"],
    "et_tran":   ["ACT. TRAN", "ACT_TRAN"],
    "et_evap":   ["ACT. EVAP", "ACT_EVAP"],
    "et_etra":   ["ACT. ETRA", "ACT_ETRA"],
}

# metric name -> value above which a warning is raised.
DEFAULT_THRESHOLDS = {
    "nlin_max": 10,
    "nlin_mean": 8,
    "back_per_1000": 20,
    "frac_at_dtmin": 0.10,
    "rel_mbe_max_pct": 1.0,
    "closure_pct": 2.0,
    "et_t0_evap": 1e-12,          # ACT. EVAP at the first output time must be ~0
    "et_etra_sum_err": 1e-9,      # |ETRA - (TRAN + EVAP)|
}


@dataclass
class PerfConfig:
    outputs: tuple = ("mbeconv", "dtcoupling", "ET")   # read_outputs() names to try
    aliases: dict = field(default_factory=lambda: {k: list(v) for k, v in DEFAULT_ALIASES.items()})
    thresholds: dict = field(default_factory=lambda: dict(DEFAULT_THRESHOLDS))
    checks: Optional[tuple] = None    # None -> run every registered check
    dtmin: Optional[float] = None     # enables frac_at_dtmin
    dtmax: Optional[float] = None     # enables frac_at_dtmax


@dataclass
class RunPerf:
    label: str
    metrics: dict
    warnings: list
    notes: list


class FrameSource:
    """Adapter: use DataFrames you loaded yourself instead of a pyCATHY object."""
    def __init__(self, **frames):
        self.frames = frames

    def read_outputs(self, name):
        return self.frames[name]


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def pick(frames: dict, key: str, aliases: dict) -> Optional[pd.Series]:
    """First column matching any alias of `key`, in any frame, as numbers."""
    for cand in aliases.get(key, [key]):
        for df in frames.values():
            if cand in df.columns:
                return pd.to_numeric(df[cand], errors="coerce")
    return None


def load_outputs(src, names) -> tuple[dict, list]:
    frames, notes = {}, []
    for n in names:
        try:
            frames[n] = src.read_outputs(n)
        except Exception as exc:                       # missing output is not fatal
            notes.append(f"{n} unreadable ({type(exc).__name__}: {exc})")
    return frames, notes


@contextlib.contextmanager
def timed(store: dict, key: str = "cpu_s"):
    """`with timed(extra): simu.run_processor(...)` -> extra['cpu_s']."""
    t0 = time.perf_counter()
    try:
        yield
    finally:
        store[key] = time.perf_counter() - t0


# --------------------------------------------------------------------------
# Check registry: each check gets (frames, cfg, extra) and returns {metric: value}
# --------------------------------------------------------------------------
CHECKS: dict[str, Callable] = {}


def register(name: str):
    def deco(fn):
        CHECKS[name] = fn
        return fn
    return deco


@register("stepping")
def _stepping(frames, cfg, extra):
    out = {}
    dt = pick(frames, "dt", cfg.aliases)
    if dt is not None:
        dt = dt.dropna()
        out.update(n_steps=len(dt), dt_min=dt.min(), dt_max=dt.max(), dt_mean=dt.mean())
        if cfg.dtmin:
            out["frac_at_dtmin"] = float((dt <= cfg.dtmin * 1.001).mean())
        if cfg.dtmax:
            out["frac_at_dtmax"] = float((dt >= cfg.dtmax * 0.999).mean())
    back = pick(frames, "back", cfg.aliases)
    if back is not None:
        out["back_steps"] = int(back.sum())
        if dt is not None and len(dt):
            out["back_per_1000"] = 1000.0 * out["back_steps"] / len(dt)
    return out


@register("iterations")
def _iterations(frames, cfg, extra):
    out = {}
    nl = pick(frames, "nlin", cfg.aliases)
    if nl is not None:
        out.update(nlin_mean=float(nl.mean()), nlin_max=int(nl.max()))
    li = pick(frames, "lin", cfg.aliases)
    if li is not None:
        out["lin_mean"] = float(li.mean())
    return out


@register("mass_balance")
def _mass_balance(frames, cfg, extra):
    out = {}
    mbe = pick(frames, "mbe", cfg.aliases)
    if mbe is not None:
        out["mbe_cum_abs"] = float(mbe.abs().sum())
    rel = pick(frames, "rel_mbe", cfg.aliases)
    if rel is not None:
        out.update(rel_mbe_max_pct=float(rel.abs().max()),
                   rel_mbe_mean_pct=float(rel.abs().mean()))
    return out


@register("flux_integral")
def _flux_integral(frames, cfg, extra):
    """Cumulative actual ET (m3) from the domain-integrated flux; negative = leaving."""
    dt, flux = pick(frames, "dt", cfg.aliases), pick(frames, "atmact", cfg.aliases)
    if dt is None or flux is None:
        return {}
    dt, flux = dt.to_numpy(float), flux.to_numpy(float)
    n = min(len(dt), len(flux))
    dt, flux = dt[:n], flux[:n]
    fwd = -np.sum(flux * dt)                      # flux of a step applied to that step
    bwd = -np.sum(flux[1:] * dt[:-1])             # flux of a step applied to the previous dt
    out = {"cum_flux_m3_forward": fwd, "cum_flux_m3_backward": bwd}
    if max(abs(fwd), abs(bwd)) > 0:
        out["cum_flux_rule_gap_pct"] = 100 * abs(fwd - bwd) / max(abs(fwd), abs(bwd))
    return out


@register("closure")
def _closure(frames, cfg, extra):
    """
    Independent closure: delta_storage + cumulative outflow ~ 0.
    Needs extra['storage_m3'] (series or [S0, S_end]); optional
    extra['other_fluxes_m3'] = net inflow through other boundaries (+ in).
    Only meaningful when everything except the atmospheric flux is accounted for.
    """
    s = extra.get("storage_m3")
    if s is None:
        return {}
    dt, flux = pick(frames, "dt", cfg.aliases), pick(frames, "atmact", cfg.aliases)
    if dt is None or flux is None:
        return {}
    n = min(len(dt), len(flux))
    outflow = -np.sum(flux.to_numpy(float)[:n] * dt.to_numpy(float)[:n])   # +ve = lost
    d_storage = float(np.asarray(s, float)[-1] - np.asarray(s, float)[0])
    residual = d_storage + outflow - float(extra.get("other_fluxes_m3", 0.0))
    out = {"closure_residual_m3": residual}
    if outflow:
        out["closure_pct"] = 100 * abs(residual) / abs(outflow)
    return out


@register("et_partition")
def _et_partition(frames, cfg, extra):
    """Consistency of an E/T output (fort.777-like). Skipped if columns are absent."""
    et = frames.get("ET")
    if et is None:
        return {}
    a = cfg.aliases
    sub = {"ET": et}
    t, tr, ev, tot = (pick(sub, k, a) for k in ("et_time", "et_tran", "et_evap", "et_etra"))
    if t is None or ev is None:
        return {}
    d = pd.DataFrame({"t": t, "tran": tr if tr is not None else 0.0, "evap": ev,
                      "tot": tot if tot is not None else np.nan})
    g = d.groupby("t").mean()                         # domain mean per output time
    out = {"et_t0_evap": float(abs(g["evap"].iloc[0]))}
    peak_t = (g["tran"] + g["evap"]).abs().idxmax()
    out["et_peak_at_t0"] = bool(peak_t == g.index[0] and len(g) > 1)
    if tot is not None:
        out["et_etra_sum_err"] = float((d["tot"] - (d["tran"] + d["evap"])).abs().max())
    return out


# --------------------------------------------------------------------------
# Evaluation, tables, comparison, sweeps
# --------------------------------------------------------------------------
def evaluate(src, label: str, cfg: Optional[PerfConfig] = None,
             extra: Optional[dict] = None) -> RunPerf:
    """Never raises: problems become entries in `notes`."""
    cfg, extra = cfg or PerfConfig(), dict(extra or {})
    frames, notes = load_outputs(src, cfg.outputs)
    metrics = {}
    if "cpu_s" in extra:
        metrics["cpu_s"] = extra["cpu_s"]
    for name in (cfg.checks or tuple(CHECKS)):
        try:
            metrics.update(CHECKS[name](frames, cfg, extra))
        except Exception as exc:
            notes.append(f"check '{name}' failed ({type(exc).__name__}: {exc})")
    warnings = []
    for key, limit in cfg.thresholds.items():
        v = metrics.get(key)
        if v is not None and not isinstance(v, bool) and v > limit:
            warnings.append(f"{key} = {v:.4g} exceeds {limit:g}")
    if metrics.get("et_peak_at_t0"):
        warnings.append("ET peak at first output time (possible t=0 artifact)")
    return RunPerf(label, metrics, warnings, notes)


def perf_table(runs: list[RunPerf]) -> pd.DataFrame:
    return pd.DataFrame({r.label: r.metrics for r in runs}).T


def table_to_markdown(df: pd.DataFrame, floatfmt: str = ".4g") -> str:
    cols = ["run"] + list(df.columns)
    rows = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for idx, r in df.iterrows():
        cells = [f"{v:{floatfmt}}" if isinstance(v, (float, np.floating)) else str(v) for v in r]
        rows.append("| " + " | ".join([str(idx)] + cells) + " |")
    return "\n".join(rows)


def compare(runs: list[RunPerf], a: str, b: str,
            metric: str = "cum_flux_m3_backward") -> Optional[float]:
    """Relative difference (%) of `metric` between two runs, relative to `a`."""
    m = {r.label: r.metrics for r in runs}
    va, vb = m.get(a, {}).get(metric), m.get(b, {}).get(metric)
    if va in (None, 0) or vb is None:
        return None
    return 100 * (vb - va) / abs(va)


def sweep(build_and_run: Callable, param: str, values, cfg: Optional[PerfConfig] = None) -> list[RunPerf]:
    """
    build_and_run(**{param: value}) -> (src, extra_dict); extra may hold 'cpu_s'
    and 'storage_m3'. Returns one RunPerf per value, labelled 'param=value'.
    """
    runs = []
    for v in values:
        src, extra = build_and_run(**{param: v})
        runs.append(evaluate(src, f"{param}={v}", cfg, extra))
    return runs


def time_axis(frames: dict, cfg: PerfConfig) -> Optional[np.ndarray]:
    """Cumulative time from the step sizes (same convention as show_dtcoupling)."""
    dt = pick(frames, "dt", cfg.aliases)
    return None if dt is None else np.concatenate([[dt.iloc[0]], dt.iloc[:-1].cumsum() + dt.iloc[0]])


# --------------------------------------------------------------------------
# Aggregation of a run made of several consecutive segments (e.g. monthly restarts)
# --------------------------------------------------------------------------
def aggregate_months(df: pd.DataFrame) -> dict:
    """One row of run-level metrics from per-segment rows (columns as in RunPerf.metrics).

    Sums: cpu_s, n_steps, back_steps, mbe_cum_abs, cum_flux_m3_*. Extremes: dt_min, dt_max,
    nlin_max, rel_mbe_max_pct. Step-weighted means: dt_mean, nlin_mean, lin_mean, frac_at_*.
    Missing columns are skipped. Adds n_segments, n_warn_segments, n_note_segments.
    """
    out = {"n_segments": len(df)}
    w = pd.to_numeric(df["n_steps"], errors="coerce") if "n_steps" in df else None
    for c in ("cpu_s", "n_steps", "back_steps", "mbe_cum_abs",
              "cum_flux_m3_forward", "cum_flux_m3_backward"):
        if c in df and df[c].notna().any():
            out[c] = float(pd.to_numeric(df[c], errors="coerce").sum())
    for c, fn in (("dt_min", "min"), ("dt_max", "max"), ("nlin_max", "max"),
                  ("rel_mbe_max_pct", "max")):
        if c in df and df[c].notna().any():
            out[c] = float(getattr(pd.to_numeric(df[c], errors="coerce"), fn)())
    for c in ("dt_mean", "nlin_mean", "lin_mean", "frac_at_dtmin", "frac_at_dtmax"):
        if c in df and w is not None:
            v = pd.to_numeric(df[c], errors="coerce")
            ok = v.notna() & w.notna() & (w > 0)
            if ok.any():
                out[c] = float((v[ok] * w[ok]).sum() / w[ok].sum())
    if out.get("n_steps"):
        out["back_per_1000"] = 1000.0 * out.get("back_steps", 0.0) / out["n_steps"]
    for src, dst in (("warnings", "n_warn_segments"), ("notes", "n_note_segments")):
        if src in df:
            out[dst] = int((df[src].fillna("").astype(str) != "").sum())
    return out
