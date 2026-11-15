#!/usr/bin/env python3
"""
post_process.py - post-process Agramon scenarios (replaces plot_scenarii_withETp_withLAI.sh)
===========================================================================================
For each sim_index: withLAI_results.py (maps, hydro, monthly) then analysis_ET_compare.py
(ETa vs EO-EB / CSV, validation panel). Also called at the end of a run with
``Agramon_withLAI_withETp.py --post-plots``.

    python post_process.py                                  # every sim_index in the log CSV
    python post_process.py --scenario 7 8 --forcing etp_ts --dem-res 5
    python post_process.py --scenario 7 --et-csv <ETA_CSV> --monthly-ts
    python post_process.py --list                           # only print the sim_index values

Exit code 0 if every step succeeded, 1 otherwise.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENT = "withETp_withLAI"
DEFAULT_LOG = f"simulation_log_{EXPERIMENT}.csv"
DEFAULT_PRJ = f"../SSHydro_{EXPERIMENT}"
FALLBACK_SCENARIOS = list(range(1, 14))
FORCINGS = ("etp", "etp_ts", "era5")


def list_logged_scenarios(log_csv, forcing: str = "") -> list[int]:
    """sim_index values in the run log; rows without forcing_mode count as 'etp'."""
    import pandas as pd
    df = pd.read_csv(log_csv)
    if forcing:
        fm = df["forcing_mode"].fillna("etp").astype(str) if "forcing_mode" in df.columns \
            else pd.Series("etp", index=df.index)
        df = df[fm == forcing]
    return sorted(int(i) for i in df["sim_index"].unique())


def build_commands(idx: int, log_csv, path2prj, dem_res=5, monthly_ts=False,
                   et_csv=None, no_csv=False) -> dict[str, list[str]]:
    """The two subprocess command lines for one scenario (pure, testable)."""
    res = [sys.executable, str(HERE / "withLAI_results.py"), "--scenario", str(idx),
           "--log-file", str(log_csv), "--path2prj", str(path2prj), "--dem-res", str(dem_res)]
    if monthly_ts:
        res.append("--monthly-ts")
    cmp_ = [sys.executable, str(HERE / "analysis_ET_compare.py"), "--scenario", str(idx),
            "--log-file", str(log_csv), "--path2prj", str(path2prj)]
    if et_csv:
        cmp_ += ["--et-csv", str(et_csv)]
    if no_csv:
        cmp_.append("--no-csv")
    return {"results": res, "et_compare": cmp_}


def post_process_scenario(idx: int, log_csv=DEFAULT_LOG, path2prj=DEFAULT_PRJ, dem_res=5,
                          monthly_ts=False, et_csv=None, no_csv=False, cwd=None) -> dict[str, bool]:
    """Run both steps for one scenario; the ET compare is skipped if results failed."""
    cmds = build_commands(idx, log_csv, path2prj, dem_res, monthly_ts, et_csv, no_csv)
    ok = {}
    for step, cmd in cmds.items():
        if step == "et_compare" and not ok.get("results", False):
            print(f"  ⊘ ET compare skipped for scenario {idx} (results step failed)")
            ok[step] = False
            continue
        t0 = time.time()
        rc = subprocess.run(cmd, cwd=str(cwd or HERE)).returncode
        ok[step] = rc == 0
        print(f"  {'✓' if rc == 0 else '✗'} {step} scenario {idx} ({time.time() - t0:.0f}s)")
    return ok


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1], formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--scenario", nargs="+", type=int, default=None, metavar="N")
    p.add_argument("--forcing", choices=FORCINGS, default="")
    p.add_argument("--dem-res", type=float, default=5, help="m; must equal the run's resample_resolution")
    p.add_argument("--et-csv", default=None, metavar="PATH")
    p.add_argument("--no-csv", action="store_true")
    p.add_argument("--monthly-ts", action="store_true")
    p.add_argument("--log-file", default=DEFAULT_LOG)
    p.add_argument("--path2prj", default=DEFAULT_PRJ)
    p.add_argument("--list", action="store_true")
    a = p.parse_args(argv)
    if a.et_csv and not Path(a.et_csv).is_file():
        p.error(f"--et-csv file not found: {a.et_csv}")
    scen = a.scenario
    if scen is None:
        if Path(a.log_file).is_file():
            scen = list_logged_scenarios(a.log_file, a.forcing)
            if not scen:
                p.error(f"no scenarios in {a.log_file}" + (f" for forcing '{a.forcing}'" if a.forcing else ""))
        else:
            print(f"WARNING: {a.log_file} not found, falling back to {FALLBACK_SCENARIOS}", file=sys.stderr)
            scen = FALLBACK_SCENARIOS
    if a.list:
        print(*scen)
        return 0
    failed = []
    for idx in scen:
        print("-" * 60 + f"\n  Scenario {idx} / {scen[-1]}  —  {time.ctime()}\n" + "-" * 60)
        ok = post_process_scenario(idx, a.log_file, a.path2prj, a.dem_res, a.monthly_ts, a.et_csv, a.no_csv)
        if not all(ok.values()):
            failed.append(idx)
    print("=" * 60)
    print("  Status   : all scenarios post-processed" if not failed else f"  FAILED   : {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
