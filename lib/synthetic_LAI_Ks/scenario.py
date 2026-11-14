#!/usr/bin/env python3
"""
scenario.py - run a batch of synthetic Ks/fire scenarios and compare them
=========================================================================
GRwater project - Agramon (synthetic variant)

Driver for synthetic_ksfire_LAI.py. Each scenario is a name plus a list of CLI
arguments (1 vs 2 layers, Ks mode, LAI reference, topsoil depth, initial
condition, ...). Every scenario is run in its own subprocess with

    --path2prj <root>/<name>

so the folder layout is exactly what synthetic_ksfire_analysis.py expects:

    <root>/<name>/<synth_ks-..._YYYY-YYYY>/outputs/areamean_timeseries.csv

and the comparison (<root>/_comparison/) can be built at the end.

Usage
-----
    python scenario.py --list                        # show scenarios and tags
    python scenario.py --dry-run                     # print the commands only
    python scenario.py                               # run the default set (tag "default")
    python scenario.py --only control control_2L lai_2L
    python scenario.py --only ks                     # a tag selects all scenarios with it
    python scenario.py --only all --skip deep_wt
    python scenario.py --start-year 2020 --end-year 2024 --root ../synthetic_LAI_Ks
    python scenario.py --jobs 3                      # run 3 scenarios in parallel
    python scenario.py --no-plots --extra --dtmin 1e-4 --deltat 10
    python scenario.py --config my_scenarios.toml    # add/override scenarios (see below)
    python scenario.py --compare-only                # skip the runs, only compare
    python scenario.py --wet --only control control_2L lai_2L   # wetter forcing + IC
    python scenario.py --quick --wet                 # 5x5 grid, fast sanity run

Custom scenarios (--config, .toml or .json)
-------------------------------------------
    [lai_2L_deep5e-5]
    args = ["--ks-fire-mode", "lai", "--permx-deep", "5e-5"]
    tags = ["custom"]
    description = "lai mode, stiffer subsoil"

A name that already exists is replaced. Matching controls for the comparison
are picked automatically by synthetic_ksfire_analysis.py (mode-off run with the
same number of layers), so keep at least one control per layer setup in the set.

Exit code: 0 all OK, 1 at least one scenario failed, 3 some were partial.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN_SCRIPT = HERE / "synthetic_ksfire_LAI.py"
ANALYSIS_SCRIPT = HERE / "synthetic_ksfire_analysis.py"

# "wetter world" preset (--wet): wetter forcing and initial condition, so that rain
# actually reaches the soil and ETa can be limited by soil water. Applied BEFORE the
# scenario's own args and --extra, so both can override it.
# IC: hydrostatic with the water table 2.5 m below the surface (mesh is 3 m deep). The
# former uniform head of -1 m saturated the whole column (Sw = 1 at 1 m and 2 m).
WET_ARGS = ["--rain-scale", "3", "--ic-wt-depth", "2.5", "--pmin", "-15",
            "--eta-net-neg-mode", "gross"]
# quick-test preset (--quick): 5x5 synthetic DEM (no spatial heterogeneity, much faster)
QUICK_ARGS = ["--n", "5"]
# Results live OUTSIDE lib/: <repo>/results/synthetic_LAI_Ks[_wet][_quick]
# (override the base with $SYNTH_KS_ROOT or --root).
RESULTS_BASE = Path(os.environ.get("SYNTH_KS_ROOT") or
                    HERE.parents[1] / "results" / "synthetic_LAI_Ks")
DEFAULT_ROOT = RESULTS_BASE
WET_ROOT = RESULTS_BASE.with_name(RESULTS_BASE.name + "_wet")

EXIT_PARTIAL = 3     # same codes as synthetic_ksfire_LAI.py
EXIT_NO_OUTPUT = 2


# -- Scenario definitions ------------------------------------------------------

@dataclass
class Scenario:
    name: str
    args: list[str]
    tags: set[str] = field(default_factory=set)
    description: str = ""


def _s(name, args, tags, description=""):
    return Scenario(name, list(args), set(tags), description)


L1 = ["--no-two-layers"]           # one homogeneous layer
L2 = ["--two-layers"]              # burned topsoil over stable subsoil (script default)

BUILTIN = [
    # controls: Ks constant before/after the fire (vegetation still burns)
    _s("control",    ["--ks-fire-mode", "off", *L1], {"default", "control", "layers1", "ks"},
       "constant Ks, 1 layer"),
    _s("control_2L", ["--ks-fire-mode", "off", *L2], {"default", "control", "layers2", "ks"},
       "constant Ks, 2 layers"),

    # LAI-driven Ks (Thomas et al. 2021, Eq. 11)
    _s("lai_1L", ["--ks-fire-mode", "lai", *L1], {"default", "lai", "layers1", "ks"},
       "Ks from LAI recovery ratio, 1 layer"),
    _s("lai_2L", ["--ks-fire-mode", "lai", *L2], {"default", "lai", "layers2", "ks"},
       "Ks from LAI recovery ratio, 2 layers"),
    _s("lai_raw_2L", ["--ks-fire-mode", "lai", "--lai-ref", "raw", *L2],
       {"lai", "layers2", "ks"}, "Ks from raw (seasonal) LAI, 2 layers"),

    # time-since-fire logistic Ks (Ebel & Martin 2017 / Thomas Eq. 9)
    _s("time_1L", ["--ks-fire-mode", "time", "--ks-fire-preset", "overall_average", *L1],
       {"default", "time", "layers1", "ks"}, "logistic recovery (overall_average), 1 layer"),
    _s("time_2L", ["--ks-fire-mode", "time", "--ks-fire-preset", "overall_average", *L2],
       {"default", "time", "layers2", "ks"}, "logistic recovery (overall_average), 2 layers"),
    _s("time_median_2L", ["--ks-fire-mode", "time", "--ks-fire-preset", "median", *L2],
       {"time", "layers2", "ks"}, "logistic recovery, fast median preset, 2 layers"),
    _s("time_lower_2L", ["--ks-fire-mode", "time", "--ks-fire-preset", "lower_bound", *L2],
       {"time", "layers2", "ks"}, "logistic recovery, lower-bound preset, 2 layers"),
    _s("time_upper_2L", ["--ks-fire-mode", "time", "--ks-fire-preset", "upper_bound", *L2],
       {"time", "layers2", "ks"}, "logistic recovery, upper-bound preset, 2 layers"),

    # static burn-severity step (Atchley et al. 2018-style)
    _s("burnsev_1L", ["--ks-fire-mode", "burn_severity", "--ks-fire-severity", "0.6", *L1],
       {"burnsev", "layers1", "ks"}, "static Ks reduction, severity 0.6, 1 layer"),
    _s("burnsev_2L", ["--ks-fire-mode", "burn_severity", "--ks-fire-severity", "0.6", *L2],
       {"burnsev", "layers2", "ks"}, "static Ks reduction, severity 0.6, 2 layers"),

    # two-layer geometry / subsoil sensitivity (lai mode)
    _s("lai_2L_top0.5", ["--ks-fire-mode", "lai", "--top-depth", "0.5"],
       {"lai", "layers2", "geometry"}, "burned topsoil 0.5 m thick"),
    _s("lai_2L_deepK5e-5", ["--ks-fire-mode", "lai", "--permx-deep", "5e-5"],
       {"lai", "layers2", "geometry"}, "stiffer stable subsoil (5e-5 m/s)"),

    # other initial condition: NOT comparable with the others -> own tag, and
    # excluded from the automatic comparison (see EXCLUDE_FROM_COMPARISON)
    _s("deep_wt", ["--ks-fire-mode", "lai", "--ic-wtposition", "8"],
       {"ic"}, "hydrostatic IC, water table 8 m above the grid base"),
]

EXCLUDE_FROM_COMPARISON = {"deep_wt"}


def load_config(path: Path) -> list[Scenario]:
    """Scenarios from a .toml or .json file: {name: {args, tags?, description?}}."""
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text())
    else:
        try:
            import tomllib
        except ImportError:                                   # Python < 3.11
            import tomli as tomllib                           # type: ignore
        data = tomllib.loads(path.read_text())
    out = []
    for name, spec in data.items():
        if "args" not in spec:
            raise ValueError(f"{path}: scenario '{name}' has no 'args' list")
        out.append(Scenario(name, [str(a) for a in spec["args"]],
                            set(spec.get("tags", ["custom"])),
                            str(spec.get("description", ""))))
    return out


def select(scenarios: dict[str, Scenario], only, skip) -> list[Scenario]:
    """Pick scenarios by name or tag (default: tag 'default'), then drop `skip`."""
    tokens = list(only) if only else ["default"]
    chosen: dict[str, Scenario] = {}
    for tok in tokens:
        if tok == "all":
            chosen.update(scenarios)
            continue
        if tok in scenarios:
            chosen[tok] = scenarios[tok]
            continue
        tagged = {n: s for n, s in scenarios.items() if tok in s.tags}
        if not tagged:
            names = ", ".join(sorted(scenarios))
            tags = ", ".join(sorted({t for s in scenarios.values() for t in s.tags}))
            raise SystemExit(f"Unknown scenario or tag '{tok}'.\n  names: {names}\n  tags : {tags}")
        chosen.update(tagged)
    for tok in skip or []:
        for n in [n for n, s in chosen.items() if n == tok or tok in s.tags]:
            del chosen[n]
    return list(chosen.values())


# -- Running -------------------------------------------------------------------

def build_command(sc: Scenario, root: Path, a: argparse.Namespace) -> list[str]:
    cmd = [sys.executable, str(RUN_SCRIPT),
           "--path2prj", str(root / sc.name),
           "--start-year", str(a.start_year), "--end-year", str(a.end_year),
           "--dpi", str(a.dpi)]
    if a.no_plots:
        cmd.append("--no-plots")
    if a.no_resume:
        cmd.append("--no-resume")
    if a.debug:
        cmd.append("--debug")
    if a.wet:
        cmd += WET_ARGS     # wetter forcing + IC (scenario args and --extra can override)
    if a.quick:
        cmd += QUICK_ARGS   # 5x5 grid
    if a.slope is not None:
        cmd += ["--slope", str(a.slope)]
    cmd += sc.args          # scenario-specific
    cmd += a.extra          # global overrides last, so they win (argparse: last value wins)
    return cmd


def run_one(sc: Scenario, root: Path, a: argparse.Namespace) -> dict:
    cmd = build_command(sc, root, a)
    sim_root = root / sc.name
    sim_root.mkdir(parents=True, exist_ok=True)
    log = sim_root / "run.log"
    t0 = time.time()
    print(f"[start] {sc.name:<18} log: {log}", flush=True)
    with open(log, "w") as fh:
        fh.write("# " + " ".join(cmd) + "\n\n")
        fh.flush()
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=HERE).returncode
    dt = time.time() - t0
    status = {0: "OK", EXIT_PARTIAL: "PARTIAL", EXIT_NO_OUTPUT: "NO OUTPUT"}.get(rc, f"FAILED({rc})")
    print(f"[done ] {sc.name:<18} {status:<10} {dt / 60:5.1f} min", flush=True)
    return {"name": sc.name, "rc": rc, "status": status, "minutes": dt / 60, "log": log}


def run_all(scs: list[Scenario], root: Path, a: argparse.Namespace) -> list[dict]:
    results = []
    if a.jobs <= 1:
        for sc in scs:
            r = run_one(sc, root, a)
            results.append(r)
            if r["rc"] not in (0, EXIT_PARTIAL) and a.stop_on_error:
                print("Stopping (--stop-on-error).")
                break
        return results
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_one, sc, root, a) for sc in scs]
        for f in as_completed(futs):
            results.append(f.result())
    order = {sc.name: i for i, sc in enumerate(scs)}
    return sorted(results, key=lambda r: order[r["name"]])


def compare(root: Path, a: argparse.Namespace, ran: list[str] | None = None) -> int:
    exclude = sorted(EXCLUDE_FROM_COMPARISON | set(a.compare_exclude))
    cmd = [sys.executable, str(ANALYSIS_SCRIPT), "--root", str(root), "--dpi", str(a.dpi)]
    if exclude:
        cmd += ["--exclude", *exclude]
    print("\n[compare] " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=HERE).returncode


def print_summary(results: list[dict]) -> None:
    print("\n" + "=" * 60)
    print(f"{'scenario':<20}{'status':<12}{'minutes':>8}  log")
    for r in results:
        print(f"{r['name']:<20}{r['status']:<12}{r['minutes']:>8.1f}  {r['log']}")
    print("=" * 60)


# -- CLI -----------------------------------------------------------------------

def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run and compare the synthetic Ks/fire scenarios",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--only", nargs="*", default=None, metavar="NAME_OR_TAG",
                   help="scenario names and/or tags ('all' = everything); default: tag 'default'")
    p.add_argument("--skip", nargs="*", default=[], metavar="NAME_OR_TAG")
    p.add_argument("--list", action="store_true", help="list scenarios and exit")
    p.add_argument("--dry-run", action="store_true", help="print the commands, run nothing")
    p.add_argument("--config", type=Path, default=None,
                   help=".toml/.json file with extra or overriding scenarios")
    p.add_argument("--root", type=Path, default=None,
                   help=f"root folder: one sub-folder per scenario, plus _comparison/ "
                        f"(default {DEFAULT_ROOT}, or {WET_ROOT} with --wet; _quick appended with --quick)")
    p.add_argument("--wet", action="store_true",
                   help="wetter forcing and initial condition: "
                        + " ".join(WET_ARGS) + ". Uses its own root so the controls "
                        "match (never mix wet and normal runs in one root)")
    p.add_argument("--quick", action="store_true",
                   help="quick test: 5x5 synthetic DEM (" + " ".join(QUICK_ARGS) + "). Uses its "
                        "own root (suffix _quick) so it never mixes with the 20x20 runs")
    p.add_argument("--slope", type=float, default=None, metavar="M/M",
                   help="DEM slope (m/m) for every run, e.g. 0.05 for 5 %%. Uses its own root "
                        "(suffix _slope<M/M>) so slopes never mix; compare slopes by running "
                        "once per slope, e.g. for s in 0 0.05 0.15; do python scenario.py "
                        "--slope $s ...; done")
    p.add_argument("--start-year", type=int, default=2020)
    p.add_argument("--end-year", type=int, default=2024,
                   help="exclusive (2024 -> simulates 2020-2023, 3.5 yr after the fire)")
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--no-resume", action="store_true",
                   help="ignore the monthly backups and re-simulate everything")
    p.add_argument("--debug", action="store_true",
                   help="each run stops at the first month the solver cannot complete")
    p.add_argument("--jobs", type=int, default=1, help="scenarios run in parallel")
    p.add_argument("--stop-on-error", action="store_true",
                   help="sequential mode: stop at the first scenario that fails")
    p.add_argument("--no-compare", action="store_true", help="do not build _comparison/")
    p.add_argument("--compare-only", action="store_true",
                   help="skip the simulations and only build _comparison/")
    p.add_argument("--compare-exclude", nargs="*", default=[], metavar="NAME",
                   help="extra scenario names to leave out of the comparison")
    p.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                   help="everything after --extra is appended to EVERY run "
                        "(e.g. --extra --dtmin 1e-4 --deltat 10)")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = _parse_args(argv)
    scenarios = {s.name: s for s in BUILTIN}
    if a.config:
        for s in load_config(a.config):
            scenarios[s.name] = s

    if a.list:
        print(f"{'name':<20}{'tags':<34}description")
        for s in scenarios.values():
            print(f"{s.name:<20}{','.join(sorted(s.tags)):<34}{s.description}")
            print(f"{'':<20}args: {' '.join(s.args)}")
        return 0

    base = a.root or (WET_ROOT if a.wet else DEFAULT_ROOT)
    if a.quick and not a.root:
        base = base.with_name(base.name + "_quick")
    if a.slope is not None and not a.root:
        base = base.with_name(base.name + f"_slope{a.slope:g}")
    root = base.expanduser().resolve()

    if a.compare_only:
        return compare(root, a)

    scs = select(scenarios, a.only, a.skip)
    if not scs:
        print("Nothing selected.")
        return 1
    print(f"{len(scs)} scenario(s): {', '.join(s.name for s in scs)}")
    print(f"root: {root} | years {a.start_year}-{a.end_year - 1} | jobs {a.jobs}")

    if a.dry_run:
        for sc in scs:
            print("\n# " + sc.name + (f" - {sc.description}" if sc.description else ""))
            print(" ".join(build_command(sc, root, a)))
        return 0

    if not RUN_SCRIPT.exists():
        raise SystemExit(f"Cannot find {RUN_SCRIPT}; scenario.py must sit next to it.")

    root.mkdir(parents=True, exist_ok=True)
    results = run_all(scs, root, a)
    print_summary(results)

    if not a.no_compare:
        usable = [r for r in results if r["rc"] in (0, EXIT_PARTIAL)]
        if len(usable) >= 2 or (root / "_comparison").exists() or len(list(root.glob("*/*/outputs"))) >= 2:
            compare(root, a)
        else:
            print("Comparison skipped: fewer than 2 usable runs.")

    if any(r["rc"] not in (0, EXIT_PARTIAL) for r in results):
        return 1
    if any(r["rc"] == EXIT_PARTIAL for r in results):
        return EXIT_PARTIAL
    return 0


if __name__ == "__main__":
    sys.exit(main())
