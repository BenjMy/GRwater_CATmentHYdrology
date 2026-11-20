"""
agramon_registry.py — one run registry on top of the CSV run log
================================================================
GRwater project — Agramon study site

The CSV run log (``simulation_log*.csv``, one row per run, ``sim_index`` + the run parameters) is
the single registry. This module adds the columns that make runs reusable and comparable:

    created_at, finished_at   ISO timestamps (seconds)
    status                    running | ok | failed   (empty in legacy rows = unknown)
    seconds                   wall time of the run
    sim_dir                   scenario folder, relative to ``path2prj`` (portable)
    param_hash                identity of the parameters (see ``param_hash``)
    hypothesis, variant       set by the hypothesis runner through env vars
                              AGRAMON_HYPOTHESIS / AGRAMON_VARIANT
    git_commit                short hash of the code that produced the run
    error_tail                last lines of the traceback of a failed run

Old CSVs without these columns load fine (``load_registry`` adds them empty). Pure pandas/stdlib:
no CATHY, no geopandas. Registry failures never stop a run (callers wrap them, see ``safe``).

    start_run(log, idx, params)         after the row is logged (agramon_runlog.log_simulation)
    finish_run(log, idx, status, ...)   at the end of the pipeline (ok / failed)
    find_run(log, params, path2prj)     sim_index of a finished run with the same parameters
    resolve_sim_dir(row, path2prj)      scenario folder of a row (relative, absolute or by index)
"""
from __future__ import annotations

import hashlib
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd

REQUIRED = ("et_output.nc", "psi_output.nc", "sw_output.nc")   # written only by a finished run

# parameters that define a run (same names as the CSV columns written by log_simulation)
IDENTITY_COLS = (
    "start_year", "end_year", "ZROOT", "PMIN", "SCF", "watershed_nb", "withLAI",
    "outlet_flux", "bottom_flux", "spinup", "ic_mode", "ic_wtposition", "ic_pressure_head",
    "zroot_factors", "permx", "permy_ratio", "permz_ratio",
    "PCANA", "PCREF", "PCWLT", "PZ", "OMGC", "feddes_crop",
    "forcing_mode", "etp_ts_column", "etp_ts_feature_id", "ks_fire",
)
REGISTRY_COLS = ("created_at", "finished_at", "status", "seconds", "sim_dir", "param_hash",
                 "n_months", "hypothesis", "variant", "git_commit", "error_tail")
PARTIAL_BELOW = 0.9   # a run with fewer months than this x the longest run of its period is partial


# ── identity ──────────────────────────────────────────────────────────────────
def _norm(v) -> str:
    """Canonical text of one value: empty/NaN -> '', numbers -> '%g' (so 2, 2.0 and '2' agree)."""
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return ""
    s = str(v).strip()
    if s == "":
        return ""
    try:
        return format(float(s), "g")
    except ValueError:
        return s


def param_hash(params) -> str:
    """12-char hash of the identity parameters of a run (dict or Series with the CSV column names).

    Legacy rows are normalised the way the pipeline reads them: empty forcing_mode means "etp",
    empty ks_fire means "off", and etp_ts_* only count for "etp_ts" runs.
    """
    p = {c: _norm(params.get(c)) for c in IDENTITY_COLS}
    p["forcing_mode"] = p["forcing_mode"] or "etp"
    p["ks_fire"] = p["ks_fire"] or "off"
    if p["forcing_mode"] != "etp_ts":
        p["etp_ts_column"] = p["etp_ts_feature_id"] = ""
    text = "|".join(f"{c}={p[c]}" for c in IDENTITY_COLS)
    return hashlib.sha1(text.encode()).hexdigest()[:12]


# ── read / write ──────────────────────────────────────────────────────────────
def load_registry(log_path) -> pd.DataFrame:
    """The run log with every registry column present; missing status reads as 'unknown'."""
    df = pd.read_csv(log_path)
    for c in REGISTRY_COLS:
        if c not in df.columns:
            df[c] = pd.Series([None] * len(df), dtype=object)
    df["seconds"] = pd.to_numeric(df["seconds"], errors="coerce")
    df["n_months"] = pd.to_numeric(df["n_months"], errors="coerce")
    df["status"] = df["status"].where(df["status"].notna(), "unknown")
    return df


def count_months(sim_dir) -> Optional[int]:
    """Months a run simulated = rows of its monthly_stats.csv (None if unreadable)."""
    try:
        return int(len(pd.read_csv(Path(sim_dir) / "monthly_stats.csv", usecols=["month"])))
    except Exception:
        return None


def partial_rows(df: pd.DataFrame) -> pd.Series:
    """True for runs that stopped early: ``n_months`` below PARTIAL_BELOW x the longest run with
    the same start/end year. A finished-looking folder (all NetCDFs) can still be partial, because
    the pipeline writes the outputs even when months were skipped. Unknown n_months -> False."""
    if "n_months" not in df.columns:
        return pd.Series(False, index=df.index)
    n = pd.to_numeric(df["n_months"], errors="coerce")
    longest = n.groupby([df["start_year"], df["end_year"]]).transform("max")
    return (n < PARTIAL_BELOW * longest).fillna(False)


def _write_atomic(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def _update_row(log_path, sim_index: int, **fields) -> bool:
    """Set ``fields`` on the row of ``sim_index`` (last one if repeated); False if not found."""
    path = Path(log_path)
    df = pd.read_csv(path)
    hit = df.index[df["sim_index"] == int(sim_index)]
    if len(hit) == 0:
        return False
    for c, v in fields.items():
        if c not in df.columns:
            df[c] = pd.Series([None] * len(df), dtype=object)
        df[c] = df[c].astype(object)
        df.loc[hit[-1], c] = v
    _write_atomic(df, path)
    return True


def _git_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", str(Path(__file__).resolve().parent),
                              "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def start_run(log_path, sim_index: int, params) -> bool:
    """Mark a freshly logged run as running and stamp its identity, time and code version."""
    return _update_row(
        log_path, sim_index,
        created_at=_now(), status="running", param_hash=param_hash(params),
        hypothesis=os.environ.get("AGRAMON_HYPOTHESIS", ""),
        variant=os.environ.get("AGRAMON_VARIANT", ""),
        git_commit=_git_commit(),
    )


def finish_run(log_path, sim_index: int, status: str, seconds: Optional[float] = None,
               sim_dir=None, path2prj=None, error_tail: str = "") -> bool:
    """Close a run: status 'ok' or 'failed', wall time, scenario folder (relative to path2prj)."""
    if status not in ("ok", "failed"):
        raise ValueError(f"status must be 'ok' or 'failed', got {status!r}")
    rel = ""
    if sim_dir is not None:
        sim_dir = Path(sim_dir)
        try:
            rel = sim_dir.resolve().relative_to(Path(path2prj).resolve()).as_posix() \
                if path2prj is not None else sim_dir.as_posix()
        except ValueError:
            rel = sim_dir.as_posix()
    tail = re.sub(r"\s*\n\s*", " | ", (error_tail or "").strip())[-500:]
    return _update_row(
        log_path, sim_index,
        finished_at=_now(), status=status,
        seconds=None if seconds is None else round(float(seconds), 1),
        sim_dir=rel, error_tail=tail,
        n_months=None if sim_dir is None else count_months(sim_dir),
    )


# ── lookup / reuse ────────────────────────────────────────────────────────────
def outputs_ok(sim_dir) -> bool:
    """True when the scenario folder holds the NetCDFs a finished run writes."""
    return sim_dir is not None and all((Path(sim_dir) / f).exists() for f in REQUIRED)


def resolve_sim_dir(row, path2prj) -> Optional[Path]:
    """Scenario folder of a registry row.

    Order: ``sim_dir`` relative to path2prj, ``sim_dir`` absolute (legacy manifests), then the
    folder ``<path2prj>/outputs/scenario_<sim_index>[_<tag>]``. None if nothing exists.
    """
    base = Path(path2prj)
    sd = row.get("sim_dir") if hasattr(row, "get") else None
    if isinstance(sd, str) and sd.strip():
        for cand in (base / sd, Path(sd)):
            if cand.is_dir():
                return cand
    out = base / "outputs"
    if out.is_dir():
        pat = re.compile(rf"^scenario_{int(row['sim_index'])}(_.+)?$")
        hits = sorted(d for d in out.iterdir() if d.is_dir() and pat.match(d.name))
        if hits:
            return hits[0]
    return None


def find_run(log_path, params, path2prj=None, require_outputs: bool = True) -> Optional[int]:
    """sim_index of the latest finished run with the same parameters, or None.

    Rows marked failed/running, and partial runs (see ``partial_rows``), never match. Legacy rows (status unknown) match when their
    outputs are on disk, so ``path2prj`` is needed unless ``require_outputs`` is False.
    """
    path = Path(log_path)
    if not path.exists() or path.stat().st_size == 0:
        return None
    df = load_registry(path)
    partial = partial_rows(df)
    h = param_hash(params)
    for i, row in df.iloc[::-1].iterrows():
        if row["status"] in ("failed", "running") or partial[i]:
            continue
        rh = row["param_hash"] if isinstance(row["param_hash"], str) and row["param_hash"] \
            else param_hash(row)
        if rh != h:
            continue
        if not require_outputs:
            return int(row["sim_index"])
        if path2prj is not None and outputs_ok(resolve_sim_dir(row, path2prj)):
            return int(row["sim_index"])
    return None


def safe(fn, *args, **kwargs):
    """Call a registry function; a failure is logged and swallowed (never stops a run)."""
    try:
        return fn(*args, **kwargs)
    except Exception as exc:   # noqa: BLE001
        try:
            from agramon_log import get_logger
            get_logger(__name__).warning(f"  registry: {fn.__name__} failed ({exc})")
        except Exception:
            print(f"  registry: {fn.__name__} failed ({exc})")
        return None
