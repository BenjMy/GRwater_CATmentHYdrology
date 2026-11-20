#!/usr/bin/env python3
"""
backfill_registry.py - add the registry columns to runs logged before they existed
==================================================================================
For every row of a CSV run log that lacks them, fills (never overwrites)

    param_hash   from the parameter columns
    sim_dir      scenario folder relative to the project (found by ``scenario_<sim_index>[_tag]``)
    status       'ok' when the folder holds et/psi/sw_output.nc, else left 'unknown'
    n_months     months simulated (rows of monthly_stats.csv); runs far shorter than the longest
                 run of their period are then treated as partial and never reused

Not filled on purpose: seconds / created_at / finished_at (nothing reliable to read them from:
file times change with a git checkout) and anything taken from manifest.json, whose sim_dir
indexes can disagree with the CSV after clean_log.py reused sim_index values.

    python backfill_registry.py --dry-run --study <out_root>/hyp1_layers_Ks
    python backfill_registry.py --study <out_root>/hyp1_layers_Ks
    python backfill_registry.py --log-file my.csv --path2prj ../prj

A timestamped backup of the CSV is written first (``--no-backup`` to skip). Rows that already have
a value keep it, so the script can be run again safely.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

from agramon_registry import (REGISTRY_COLS, _write_atomic, count_months, load_registry,
                              outputs_ok, param_hash, resolve_sim_dir)


def _empty(v) -> bool:
    return v is None or (not isinstance(v, str) and pd.isna(v)) or v in ("", "unknown")


def backfill(log: Path, prj: Path, *, dry_run: bool = False, backup: bool = True) -> dict:
    """Fill missing registry values; returns counts and the rows that have no scenario folder."""
    df = load_registry(log)
    df = df.astype({c: object for c in REGISTRY_COLS if c != "seconds"})
    n = {"hash": 0, "sim_dir": 0, "ok": 0, "months": 0}
    no_dir, no_outputs = [], []
    for i, row in df.iterrows():
        if _empty(row["param_hash"]):
            df.at[i, "param_hash"] = param_hash(row)
            n["hash"] += 1
        d = resolve_sim_dir(row, prj)
        if d is None:
            no_dir.append(int(row["sim_index"]))
            continue
        if _empty(row["sim_dir"]):
            df.at[i, "sim_dir"] = d.resolve().relative_to(prj.resolve()).as_posix()
            n["sim_dir"] += 1
        if _empty(row["n_months"]):
            m = count_months(d)
            if m is not None:
                df.at[i, "n_months"] = m
                n["months"] += 1
        if _empty(row["status"]):
            if outputs_ok(d):
                df.at[i, "status"] = "ok"
                n["ok"] += 1
            else:
                no_outputs.append(int(row["sim_index"]))
    res = {"rows": len(df), **n, "no_folder": no_dir, "no_outputs": no_outputs}
    if not dry_run and any(n.values()):
        if backup:
            bak = log.with_name(f"{log.name}.bak_{datetime.now():%Y%m%d-%H%M%S}")
            shutil.copy2(log, bak)
            res["backup"] = str(bak)
        _write_atomic(df, log)
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--log-file", default="simulation_log_withETp_withLAI.csv", help="CSV run log")
    ap.add_argument("--path2prj", default="../SSHydro_withETp_withLAI", help="project folder holding outputs/")
    ap.add_argument("--study", help="hypothesis study folder (uses its simulation_log.csv and prj/)")
    ap.add_argument("--dry-run", action="store_true", help="report only, change nothing")
    ap.add_argument("--no-backup", action="store_true", help="do not write the .bak copy of the CSV")
    a = ap.parse_args(argv)
    log, prj = (Path(a.study) / "simulation_log.csv", Path(a.study) / "prj") if a.study \
        else (Path(a.log_file), Path(a.path2prj))
    if not log.exists():
        print(f"log not found: {log}")
        return 1
    res = backfill(log, prj, dry_run=a.dry_run, backup=not a.no_backup)
    print(f"{log}  ({res['rows']} rows)")
    print(f"  param_hash filled: {res['hash']} | sim_dir filled: {res['sim_dir']} | status ok: {res['ok']} | n_months: {res['months']}")
    if res["no_folder"]:
        print(f"  no scenario folder for sim_index: {res['no_folder']}")
    if res["no_outputs"]:
        print(f"  folder without et/psi/sw_output.nc (status left unknown): {res['no_outputs']}")
    if res.get("backup"):
        print(f"  backup: {res['backup']}")
    print("Dry run: nothing changed." if a.dry_run else "Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
