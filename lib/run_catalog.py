#!/usr/bin/env python3
"""
run_catalog.py - summary of the runs in a CSV run log
=====================================================
One table with a row per run: registry columns (status, seconds, folder, hash), whether the
outputs are on disk, duplicates (same parameters as an earlier run), which parameters actually
differ between the selected runs, and optionally water-balance metrics read from each run's
``monthly_stats.csv`` (hypothesis_tables). Use it to pick runs for plot_compare_runs.py.

    python run_catalog.py --study <out_root>/hyp1_layers_Ks
    python run_catalog.py --study <study> --metrics                    # + ETa, recharge, closure
    python run_catalog.py --study <study> --where "permx == 0.0006"    # pandas query on log columns
    python run_catalog.py --study <study> --unique --status ok --csv runs_summary.csv
    python run_catalog.py --log-file my.csv --path2prj ../prj

In Python::

    from run_catalog import load_runs, varying_params, run_labels
    runs = load_runs(log, prj, where="PMIN == -50", unique=True, metrics=True)

Works on old logs (no registry columns): status shows 'unknown', hashes are computed on the fly.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

from agramon_registry import (IDENTITY_COLS, _norm, load_registry, outputs_ok, param_hash,
                              partial_rows, resolve_sim_dir)

METRIC_COLS = ("annual_precip_mm", "annual_eta_mm", "annual_recharge_mm", "mean_sw", "closure_pct")


def load_runs(log, prj, *, where: Optional[str] = None, status: Optional[str] = None,
              unique: bool = False, metrics: bool = False, fire_date=None) -> pd.DataFrame:
    """Registry rows + ``sim_path``, ``outputs_ok``, ``dup_of`` (+ metrics), sorted by sim_index.

    where    pandas query on the log columns, e.g. ``"PMIN == -50 and permx > 1e-4"``
    status   keep one status (ok / failed / running / unknown)
    unique   keep only the latest run of each parameter set
    metrics  add METRIC_COLS from monthly_stats.csv (whole-period values, mm/year)
    """
    df = load_registry(log)
    df["param_hash"] = [h if isinstance(h, str) and h else param_hash(r)
                        for h, (_, r) in zip(df["param_hash"], df.iterrows())]
    first = df.drop_duplicates("param_hash").set_index("param_hash")["sim_index"]
    df["dup_of"] = [int(first[h]) if int(first[h]) != int(i) else pd.NA
                    for h, i in zip(df["param_hash"], df["sim_index"])]
    paths = [resolve_sim_dir(r, prj) for _, r in df.iterrows()]
    df["sim_path"] = paths
    df["outputs_ok"] = [outputs_ok(p) for p in paths]
    df["partial"] = partial_rows(df)
    if where:
        df = df.query(where)
    if status:
        df = df[df["status"] == status]
    if unique:
        df = df.drop_duplicates("param_hash", keep="last")
    df = df.sort_values("sim_index").reset_index(drop=True)
    if metrics:
        df = _add_metrics(df, fire_date)
    return df


def _add_metrics(df: pd.DataFrame, fire_date=None) -> pd.DataFrame:
    import hypothesis_tables as HT
    fire = fire_date or HT.default_fire_date()
    rows = []
    for p in df["sim_path"]:
        row = dict.fromkeys(METRIC_COLS, float("nan"))
        if p is not None:
            monthly, _ = HT.load_monthly_balance(p)
            tot = HT.period_totals(monthly, fire)
            tot = tot[tot["period"] == "all"]
            if not tot.empty:
                t = tot.iloc[0]
                row.update(annual_precip_mm=t["annual_precip_mm"], annual_eta_mm=t["annual_eta_mm"],
                           annual_recharge_mm=t["annual_recharge_mm"], mean_sw=t["mean_sw"],
                           closure_pct=t["closure_pct"])
        rows.append(row)
    return pd.concat([df, pd.DataFrame(rows, index=df.index)], axis=1)


def varying_params(df: pd.DataFrame) -> list:
    """Identity columns whose (normalised) value differs between the runs of ``df``."""
    return [c for c in IDENTITY_COLS
            if c in df.columns and df[c].map(_norm).nunique() > 1]


def run_labels(df: pd.DataFrame, varying: Optional[list] = None) -> dict:
    """``{sim_index: "#3 permx=0.0006, ks_fire=off|L2=0.3:1e-05"}`` from the varying parameters."""
    varying = varying_params(df) if varying is None else varying
    out = {}
    for _, r in df.iterrows():
        parts = [f"{c}={_norm(r[c]) or '-'}" for c in varying]
        out[int(r["sim_index"])] = f"#{int(r['sim_index'])}" + (" " + ", ".join(parts) if parts else "")
    return out


def format_table(df: pd.DataFrame, varying: list) -> str:
    cols = ["sim_index", "status", "seconds", "n_months", "outputs_ok", "partial", "dup_of", *varying,
            *[c for c in METRIC_COLS if c in df.columns]]
    t = df[cols].copy()
    t["outputs_ok"] = t["outputs_ok"].map({True: "yes", False: "NO"})
    t["partial"] = t["partial"].map({True: "PARTIAL", False: ""})
    t = t.rename(columns={"outputs_ok": "outputs"})
    return t.to_string(index=False, na_rep="-", float_format=lambda v: f"{v:.4g}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--log-file", default="simulation_log_withETp_withLAI.csv", help="CSV run log")
    ap.add_argument("--path2prj", default="../SSHydro_withETp_withLAI", help="project folder holding outputs/")
    ap.add_argument("--study", help="hypothesis study folder (uses its simulation_log.csv and prj/)")
    ap.add_argument("--where", help='pandas query on the log columns, e.g. "PMIN == -50"')
    ap.add_argument("--status", help="keep one status: ok, failed, running, unknown")
    ap.add_argument("--unique", action="store_true", help="only the latest run of each parameter set")
    ap.add_argument("--metrics", action="store_true", help="add ETa / recharge / closure from monthly_stats.csv")
    ap.add_argument("--csv", help="also write the table to this CSV")
    a = ap.parse_args(argv)
    log, prj = (Path(a.study) / "simulation_log.csv", Path(a.study) / "prj") if a.study \
        else (Path(a.log_file), Path(a.path2prj))
    if not log.exists():
        print(f"log not found: {log}")
        return 1
    runs = load_runs(log, prj, where=a.where, status=a.status, unique=a.unique, metrics=a.metrics)
    if runs.empty:
        print("no run matches")
        return 0
    varying = varying_params(runs)
    print(f"{log}: {len(runs)} run(s), {runs['param_hash'].nunique()} distinct parameter set(s)")
    print(f"parameters that vary: {', '.join(varying) or '(none)'}\n")
    print(format_table(runs, varying))
    bad = runs[(runs["status"] == "failed") | runs["partial"]
               | (~runs["outputs_ok"] & (runs["status"] != "running"))]
    if len(bad):
        print(f"\ncheck: sim_index {', '.join(map(str, bad['sim_index']))} failed, partial or "
              f"without outputs (not reused)")
    if a.csv:
        runs.drop(columns=["sim_path"]).to_csv(a.csv, index=False)
        print(f"\nwritten: {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
