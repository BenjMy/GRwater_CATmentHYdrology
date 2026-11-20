"""
Scenarios of the LAI-loop experiment: the SAME model run with different CATHY versions and / or LAI treatments,
then compared with the Lanzhot eddy-covariance ETa and piezometers (compare_obs_pycathy.py).

    python scenarios_lai_pycathy.py --list
    python scenarios_lai_pycathy.py                                   # default scenarios: obs_withIrr + clim_withIrr
    python scenarios_lai_pycathy.py --scenarios obs_withIrr obs_v1    # compare two CATHY versions on the observed LAI
    python scenarios_lai_pycathy.py --scenarios all --start 2016-01-01 --end 2024-12-31 --dx 50
    python scenarios_lai_pycathy.py --compare-only                    # runs already done
    python scenarios_lai_pycathy.py --only obs_withIrr --months-limit 2 --dx 100      # smoke test of one scenario

Each scenario = one CATHY project, built and run by run_lai_loop_pycathy.py in  <base>/<scenario>/ .
    obs_*   : observed LAI (the LAI drop is in the model)              -> what CATHY makes of the deforestation
    clim_*  : counterfactual LAI = monthly median of --ref-years      -> same model WITHOUT the LAI drop
    The difference obs - clim isolates the effect of the LAI change on ETa and groundwater level.

Any option not listed here (--dx, --lai-file, --etp-mode, --zroot, --pmin, --months-limit, ...) is forwarded unchanged to
run_lai_loop_pycathy.py. Edit SCENARIOS below to add a scenario (CATHY version, LAI mode, extra runner options).
"""
import argparse
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------------------------------------------
# SCENARIO TABLE.  version = value given to pyCATHY.cathy_tools.CATHY(version=...).
#   "withIrr" : CATHY_src_WithIrrigation (the one of run_real_case_pycathy.py)           [known to work in your setup]
#   "1.0.0"   : pyCATHY default (bitbucket cathy master)                                 [PLACEHOLDER: confirm / replace]
# Put here the other CATHY versions you want to try (TO CONFIRM with you, see README "Open questions").
# ---------------------------------------------------------------------------------------------------------------
SCENARIOS = {
    "obs_withIrr":   dict(version="withIrr", lai_mode="observed",    extra=[], desc="CATHY withIrr, observed LAI (with the LAI drop)"),
    "clim_withIrr":  dict(version="withIrr", lai_mode="climatology", extra=[], desc="CATHY withIrr, LAI climatology (no LAI drop)"),
    "obs_v1":        dict(version="1.0.0",   lai_mode="observed",    extra=[], desc="CATHY 1.0.0, observed LAI"),
    "clim_v1":       dict(version="1.0.0",   lai_mode="climatology", extra=[], desc="CATHY 1.0.0, LAI climatology"),
    "obs_withIrr_etrkc": dict(version="withIrr", lai_mode="observed", extra=["--etp-mode", "etr_kc"],
                              desc="CATHY withIrr, observed LAI, ETp = ETr x Kc(LAI) (explicit LAI -> ETp coupling)"),
}
DEFAULT = ["obs_withIrr", "clim_withIrr"]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", default=str(HERE / "CS_LAI_scenarios"), help="folder holding one CATHY project per scenario")
ap.add_argument("--scenarios", nargs="+", default=DEFAULT, help="scenario names (see --list) or 'all'")
ap.add_argument("--list", action="store_true", help="print the scenario table and exit")
ap.add_argument("--only", default=None, help="run only this scenario (still compared with the others already run)")
ap.add_argument("--compare-only", action="store_true", help="skip the runs, only build the comparison")
ap.add_argument("--no-compare", action="store_true", help="skip the comparison")
ap.add_argument("--cmp-args", default="", help='extra options of compare_obs_pycathy.py, one string, e.g. "--piezo-xy 123 456"')
ap.add_argument("--keep-going", action="store_true", help="continue with the next scenario if one fails")
a, forward = ap.parse_known_args()

if a.list:
    for k, v in SCENARIOS.items():
        print(f"{k:20s} version={v['version']:8s} lai_mode={v['lai_mode']:12s} {' '.join(v['extra']):24s} {v['desc']}")
    sys.exit(0)

names = list(SCENARIOS) if a.scenarios == ["all"] else a.scenarios
bad = [n for n in names if n not in SCENARIOS]
if bad:
    sys.exit(f"unknown scenario(s) {bad}; --list shows: {list(SCENARIOS)}")
base = Path(a.base)

if not a.compare_only:
    base.mkdir(parents=True, exist_ok=True)
    failed = []
    for n in names:
        if a.only and n != a.only:
            continue
        s = SCENARIOS[n]
        cmd = [sys.executable, str(HERE / "run_lai_loop_pycathy.py"), "--project", str(base / n) + "/",
               "--version", s["version"], "--lai-mode", s["lai_mode"], *s["extra"], *forward]
        print(f"\n===== scenario '{n}': {s['desc']} =====\n  " + " ".join(shlex.quote(c) for c in cmd), flush=True)
        if subprocess.run(cmd, cwd=HERE).returncode != 0:
            failed.append(n)
            print(f"!!! scenario '{n}' FAILED (see the messages above)", flush=True)
            if not a.keep_going:
                sys.exit(f"STOPPED at '{n}' (--keep-going to continue with the next scenario)")
    if failed:
        print(f"failed scenarios: {failed}")

if not a.no_compare:
    ran = [n for n in names if any((base / n).rglob("timeseries_daily.csv"))]
    if not ran:
        sys.exit("no scenario output found to compare")
    cmd = [sys.executable, str(HERE / "compare_obs_pycathy.py"), "--base", str(base), "--scenarios", *ran, *shlex.split(a.cmp_args)]
    print("\n===== comparison with observations =====\n  " + " ".join(shlex.quote(c) for c in cmd), flush=True)
    sys.exit(subprocess.run(cmd, cwd=HERE).returncode)
