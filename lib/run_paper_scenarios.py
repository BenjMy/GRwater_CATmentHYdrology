#!/usr/bin/env python3
"""
run_paper_scenarios.py - run the Agramon paper scenarios listed in paper_scenarios.toml
=======================================================================================
Reads the TOML (one table per scenario, see its header), expands [name.loop] into one
call of run_scenarii_withETp_withLAI_ALL.sh per combination and runs them.

    python run_paper_scenarios.py --list
    python run_paper_scenarios.py --dry-run                  # default set: tag "paper"
    python run_paper_scenarios.py --only paper_good ks       # names and/or tags ('all')
    python run_paper_scenarios.py --only all --skip experimental --stop-on-error
    python run_paper_scenarios.py --extra --start-year 2018  # appended to every call

Exit code 0 if every call returned 0, 1 otherwise. Stdlib only (tomllib: Python >= 3.11).
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "paper_scenarios.toml"
DEFAULT_RUNNER = HERE / "run_scenarii_withETp_withLAI_ALL.sh"
DEFAULT_TAG = "paper"


@dataclass
class Scenario:
    name: str
    args: dict
    loop: dict = field(default_factory=dict)
    tags: set = field(default_factory=set)
    description: str = ""


def load_config(path: Path) -> dict[str, Scenario]:
    import tomllib
    data = tomllib.loads(Path(path).read_text())
    out = {}
    for name, spec in data.items():
        if not isinstance(spec, dict) or "args" not in spec:
            raise ValueError(f"{path}: scenario '{name}' needs an [{name}.args] table")
        loop = spec.get("loop", {})
        for k, v in loop.items():
            if not isinstance(v, list) or not v:
                raise ValueError(f"{path}: {name}.loop.{k} must be a non-empty list")
        both = set(loop) & set(spec["args"])
        if both:
            raise ValueError(f"{path}: {name}: {sorted(both)} in both args and loop")
        out[name] = Scenario(name, spec["args"], loop, set(spec.get("tags", [])),
                             str(spec.get("description", "")))
    return out


def _flag(key: str, value) -> list[str]:
    """'permx', '1e-4' -> ['--permx', '1e-4'] ; list -> several values ; True -> flag only."""
    if value is True:
        return [f"--{key}"]
    if value is False or value is None:
        return []
    vals = value if isinstance(value, list) else [value]
    return [f"--{key}", *[str(v) for v in vals]]


def expand(sc: Scenario) -> list[list[str]]:
    """Runner argument lists: one per combination of the loop options (one if no loop)."""
    base = [t for k, v in sc.args.items() for t in _flag(k, v)]
    if not sc.loop:
        return [base]
    keys = list(sc.loop)
    return [base + [t for k, v in zip(keys, combo) for t in _flag(k, v)]
            for combo in itertools.product(*(sc.loop[k] for k in keys))]


def select(scenarios: dict, only, skip) -> list[Scenario]:
    chosen: dict[str, Scenario] = {}
    for tok in (only or [DEFAULT_TAG]):
        if tok == "all":
            chosen.update(scenarios)
        elif tok in scenarios:
            chosen[tok] = scenarios[tok]
        else:
            tagged = {n: s for n, s in scenarios.items() if tok in s.tags}
            if not tagged:
                tags = ", ".join(sorted({t for s in scenarios.values() for t in s.tags}))
                raise SystemExit(f"Unknown scenario or tag '{tok}'.\n  names: "
                                 f"{', '.join(scenarios)}\n  tags : {tags}")
            chosen.update(tagged)
    for tok in skip or []:
        for n in [n for n, s in chosen.items() if n == tok or tok in s.tags]:
            del chosen[n]
    return list(chosen.values())


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    p.add_argument("--only", nargs="*", default=None, metavar="NAME_OR_TAG",
                   help=f"names and/or tags ('all' = everything); default: tag '{DEFAULT_TAG}'")
    p.add_argument("--skip", nargs="*", default=[], metavar="NAME_OR_TAG")
    p.add_argument("--list", action="store_true")
    p.add_argument("--dry-run", action="store_true", help="print the commands, run nothing")
    p.add_argument("--stop-on-error", action="store_true")
    p.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                   help="appended to every runner call (put last)")
    a = p.parse_args(argv)

    scenarios = load_config(a.config)
    if a.list:
        for s in scenarios.values():
            print(f"{s.name:<28}{','.join(sorted(s.tags)):<26}{len(expand(s))} call(s)  {s.description}")
        return 0
    scs = select(scenarios, a.only, a.skip)
    if not scs:
        print("Nothing selected.")
        return 1
    calls = [(s.name, ["bash", str(a.runner), *c, *a.extra]) for s in scs for c in expand(s)]
    print(f"{len(scs)} scenario(s), {len(calls)} runner call(s)")
    failed = []
    for name, cmd in calls:
        print(f"\n# {name}\n" + " ".join(cmd), flush=True)
        if a.dry_run:
            continue
        rc = subprocess.run(cmd, cwd=HERE).returncode
        if rc:
            failed.append((name, rc))
            print(f"!! {name} exited with {rc}", flush=True)
            if a.stop_on_error:
                break
    if failed:
        print("\nFailed: " + ", ".join(f"{n} ({rc})" for n, rc in failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
