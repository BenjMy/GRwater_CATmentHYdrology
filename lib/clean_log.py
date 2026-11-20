#!/usr/bin/env python3
"""
clean_log.py - drop the rows of the CSV run log whose scenario folder no longer exists
======================================================================================
Every run of Agramon_withLAI_withETp.py appends one row (``sim_index`` + its parameters) to the
CSV log and writes ``<path2prj>/outputs/scenario_<sim_index>_<tag>/``. When a scenario folder is
deleted by hand the row stays and later confuses ``--list-scenarios``, the skip-if-done check of
scenario_runner.py and the summaries. This script finds those rows and removes them (a timestamped
backup of the CSV is written first).

    python clean_log.py --dry-run                         # default log + project of the driver
    python clean_log.py                                   # clean it
    python clean_log.py --log-file my.csv --path2prj ../prj
    python clean_log.py --study  <out_root>/hyp1_layers_Ks   # a hypothesis study folder
    python clean_log.py --all-studies [<out_root>]        # every study folder of the out_root
    python clean_log.py --incomplete                      # also drop runs that never finished
                                                          # (no et/psi/sw_output.nc)

A *study folder* (run_all_hypotheses.py) holds ``simulation_log.csv``, ``prj/`` and
``manifest.json``; the manifest entries whose scenario folder is gone are removed too, so the
variant is rerun next time.

Notes: the next sim_index is max(sim_index) + 1, so after removing the LAST rows their numbers are
used again (no folder carries them, so nothing collides). Scenario folders without a log row are
only reported, never deleted. Exit code 0 if nothing was wrong or everything was cleaned, 2 on
``--dry-run`` when rows would be removed.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_LOG = "simulation_log_withETp_withLAI.csv"      # defaults of the driver (relative to lib/)
DEFAULT_PRJ = "../SSHydro_withETp_withLAI"
REQUIRED = ("et_output.nc", "psi_output.nc", "sw_output.nc")   # written only by a finished run


@dataclass
class CleanResult:
    log: Path
    kept: int = 0
    removed: list = field(default_factory=list)       # sim_index of the removed rows
    incomplete: list = field(default_factory=list)    # logged, folder exists, outputs missing
    unreadable: int = 0                               # rows without a usable sim_index (kept)
    orphans: list = field(default_factory=list)       # scenario folders with no log row
    manifest_removed: list = field(default_factory=list)
    backup: Path | None = None
    written: bool = False


def _abs(p, base: Path = HERE) -> Path:
    p = Path(p).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def scenario_dirs(prj: Path, sim_index: int) -> list[Path]:
    """Folders of one run: scenario_<idx>_<tag> (and the untagged legacy scenario_<idx>)."""
    out = prj / "outputs"
    if not out.exists():
        return []
    found = [d for d in out.glob(f"scenario_{sim_index}_*") if d.is_dir()]
    legacy = out / f"scenario_{sim_index}"
    return found + ([legacy] if legacy.is_dir() else [])


def is_finished(folder: Path) -> bool:
    return all((folder / f).exists() for f in REQUIRED)


def clean_log(log: Path, prj: Path, *, incomplete: bool = False, dry_run: bool = False,
              backup: bool = True, manifest: Path | None = None) -> CleanResult:
    """Remove the rows of ``log`` whose run folder is gone (and, with ``incomplete``, the ones
    whose folder lacks the final NetCDFs). ``manifest``: a study's manifest.json, cleaned the
    same way. Nothing is written when ``dry_run`` or when no row has to go."""
    import pandas as pd
    res = CleanResult(log=Path(log))
    if not res.log.exists() or res.log.stat().st_size == 0:
        raise FileNotFoundError(f"log file not found or empty: {res.log}")
    df = pd.read_csv(res.log, dtype=str, keep_default_na=False)    # text in, same text out
    if "sim_index" not in df.columns:
        raise ValueError(f"{res.log}: no 'sim_index' column")
    idx = pd.to_numeric(df["sim_index"], errors="coerce")
    drop = pd.Series(False, index=df.index)
    for i, v in idx.items():
        if pd.isna(v):
            res.unreadable += 1                       # cannot tell: keep the row
            continue
        folders = scenario_dirs(Path(prj), int(v))
        if not folders:
            drop[i] = True
        elif incomplete and not any(is_finished(f) for f in folders):
            drop[i] = True
            res.incomplete.append(int(v))
    res.removed = sorted(int(v) for v in idx[drop])
    res.kept = int((~drop).sum())

    logged = {int(v) for v in idx.dropna()}
    out = Path(prj) / "outputs"
    if out.exists():
        for d in sorted(out.glob("scenario_*")):
            head = d.name.split("_")[1] if d.name.count("_") >= 1 else ""
            if d.is_dir() and head.isdigit() and int(head) not in logged:
                res.orphans.append(d.name)

    if manifest is not None and Path(manifest).exists():
        man = json.loads(Path(manifest).read_text())
        stale = [k for k, e in man.items()
                 if e.get("sim_dir") and not Path(e["sim_dir"]).exists()]
        res.manifest_removed = stale
    if dry_run:
        return res
    if res.removed:
        if backup:
            res.backup = res.log.with_name(f"{res.log.name}.bak-{time.strftime('%Y%m%d-%H%M%S')}")
            shutil.copy2(res.log, res.backup)
        tmp = res.log.with_name(res.log.name + ".tmp")
        df[~drop].to_csv(tmp, index=False)
        tmp.replace(res.log)                          # atomic: never a half-written log
        res.written = True
    if manifest is not None and res.manifest_removed:
        for k in res.manifest_removed:
            man.pop(k, None)
        Path(manifest).write_text(json.dumps(man, indent=2, default=str))
    return res


def study_paths(study: Path) -> tuple[Path, Path, Path]:
    """(log, prj, manifest) of a hypothesis study folder."""
    return study / "simulation_log.csv", study / "prj", study / "manifest.json"


def report(res: CleanResult, prj: Path, dry_run: bool) -> None:
    verb = "would remove" if dry_run else "removed"
    print(f"\n{res.log}\n  project : {prj}")
    print(f"  rows kept {res.kept}, {verb} {len(res.removed)}"
          + (f" (sim_index {', '.join(map(str, res.removed))})" if res.removed else ""))
    if res.incomplete:
        print(f"  of which never finished (no et/psi/sw_output.nc): {', '.join(map(str, res.incomplete))}")
    if res.unreadable:
        print(f"  {res.unreadable} row(s) without a numeric sim_index were kept")
    if res.manifest_removed:
        print(f"  manifest.json: {verb} {len(res.manifest_removed)} stale variant(s): "
              f"{', '.join(res.manifest_removed)}")
    if res.orphans:
        print(f"  scenario folder(s) with no log row (not touched): {', '.join(res.orphans)}")
    if res.backup:
        print(f"  backup  : {res.backup}")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--log-file", default=DEFAULT_LOG, help="CSV run log")
    p.add_argument("--path2prj", default=DEFAULT_PRJ, help="project folder holding outputs/")
    p.add_argument("--study", type=Path, default=None,
                   help="hypothesis study folder (uses its simulation_log.csv, prj/, manifest.json)")
    p.add_argument("--all-studies", nargs="?", const="", default=None, metavar="OUT_ROOT",
                   help="every study folder under OUT_ROOT (default: the hypotheses out_root)")
    p.add_argument("--incomplete", action="store_true",
                   help="also remove runs whose folder has no finished et/psi/sw_output.nc")
    p.add_argument("--dry-run", action="store_true", help="report only, change nothing")
    p.add_argument("--no-backup", action="store_true", help="do not write the .bak copy of the CSV")
    a = p.parse_args(argv)

    targets = []                                      # (log, prj, manifest)
    if a.all_studies is not None:
        if a.all_studies:
            root = Path(a.all_studies)
        else:
            from hypotheses_common import default_out_root
            root = default_out_root()
        targets = [study_paths(d) for d in sorted(root.iterdir())
                   if d.is_dir() and (d / "simulation_log.csv").exists()]
        if not targets:
            print(f"No study folder with a simulation_log.csv under {root}")
            return 1
    elif a.study:
        targets = [study_paths(a.study)]
    else:
        targets = [(_abs(a.log_file), _abs(a.path2prj), None)]

    pending, rc = 0, 0
    for log, prj, man in targets:
        try:
            res = clean_log(log, prj, incomplete=a.incomplete, dry_run=a.dry_run,
                            backup=not a.no_backup, manifest=man)
        except (FileNotFoundError, ValueError) as exc:
            print(f"\n!! {log}: {exc}")
            rc = 1
            continue
        report(res, prj, a.dry_run)
        pending += len(res.removed) + len(res.manifest_removed)
    if a.dry_run and pending:
        print("\nDry run: nothing changed. Run again without --dry-run to clean.")
        return 2
    return rc


if __name__ == "__main__":
    sys.exit(main())
