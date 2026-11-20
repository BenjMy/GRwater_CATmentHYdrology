#!/usr/bin/env python3
"""
scenario_runner.py - Python batch runner for Agramon_withLAI_withETp.py
=======================================================================
Replaces ``run_scenarii_withETp_withLAI_ALL.sh`` (kept as a 3-line wrapper around this file).
Same flags, same run order, same skip-if-done rule; the options can now also come from a TOML
table (``run_paper_scenarios.py`` / ``paper_scenarios.toml``) instead of a long command line.

Crosses dem-plot x PMIN x outlet BC x bottom BC x spin-up x initial condition (one driver call
per combination), forwards the remaining options to every call, and skips a combination whose
parameters are already in the CSV run log AND whose scenario folder holds et/psi/sw_output.nc.

    python scenario_runner.py --dem-plot 2 --pmin -5 -10 -20 --outlet free --bottom closed \\
        --spinup on --spinup-cycles 3
    python scenario_runner.py --dem-plot 3 --pmin -50 -1e35 --ic both --permx 6e-4 --dry-run
    python scenario_runner.py --help

Values are given as separate tokens (``--pmin -1e35``); ``--pmin=-1e35`` also works. They reach
the driver in the ``--flag=value`` form, so negative values in scientific notation are not
mistaken for options by argparse.
"""
from __future__ import annotations

import itertools
import math
import re
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
DRIVER = "Agramon_withLAI_withETp.py"

NUM = re.compile(r"^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$")
KS_FIRE_FLAGS = ("ks-fire-mode", "ks-fire-preset", "ks-fire-date", "ks-fire-ll", "ks-fire-ul",
                 "ks-fire-rate", "ks-fire-lai-ll", "ks-fire-lai-ul", "ks-fire-severity",
                 "ks-fire-severity-reduction")

# flag -> destination (dest names are the keys of BatchConfig.opts)
SCALARS = {
    "forcing": "forcing", "etp-ts-column": "etp_ts_column",
    "etp-ts-feature-id": "etp_ts_feature_id", "dem-smooth-sigma": "dem_smooth_sigma",
    "start-year": "start_year", "end-year": "end_year",
    "outlet": "outlet", "outlet-flux": "outlet_flux",
    "bottom": "bottom", "bottom-flux": "bottom_flux",
    "spinup": "spinup", "spinup-cycles": "spinup_cycles",
    "ic": "ic", "ic-mode": "ic", "ic-wtposition": "ic_wtposition",
    "ic-pressure-head": "ic_pressure_head",
    "zroot": "zroot", "permx": "permx", "permy-ratio": "permy_ratio",
    "permz-ratio": "permz_ratio", "permx-bottom": "permx_bottom",
    "ks-layer-depth": "ks_layer_depth", "scf": "scf",
    "feddes-crop": "feddes_crop", "pcana": "pcana", "pcref": "pcref", "pcwlt": "pcwlt",
    "pz": "pz", "omgc": "omgc",
    "path2prj": "path2prj", "log-file": "log_file",
    # forwarded to the driver, not part of a run's identity
    "post-et-csv": "post_et_csv", "log-level": "log_level", "cathy-backend": "cathy_backend",
    **{f: f.replace("-", "_") for f in KS_FIRE_FLAGS},
}
LISTS = {"dem-plot": "dem_plot", "pmin": "pmin", "zroot-factors": "zroot_factors"}
BOOLS = {"force": "force", "dry-run": "dry_run", "stop-on-error": "stop_on_error",
         "post-plots": "post_plots", "no-plots": "no_plots", "debug": "debug"}

DEFAULTS = {
    "forcing": "etp_ts", "etp_ts_column": "ET_0_wet", "etp_ts_feature_id": "",
    "dem_smooth_sigma": "0.0", "start_year": "2016", "end_year": "2026",
    "outlet": "closed", "outlet_flux": "-1e-6", "bottom": "closed", "bottom_flux": "-1e-6",
    "spinup": "off", "spinup_cycles": "3", "ic": "both", "ic_wtposition": "1",
    "ic_pressure_head": "-10", "zroot": "1.0", "permx": "", "permy_ratio": "1.0",
    "permz_ratio": "1.0", "permx_bottom": "", "ks_layer_depth": "", "scf": "1.0",
    "feddes_crop": "", "pcana": "", "pcref": "", "pcwlt": "", "pz": "", "omgc": "",
    # the driver's own defaults (relative to lib/): the skip check must read the same files
    "path2prj": "../SSHydro_withETp_withLAI", "log_file": "simulation_log_withETp_withLAI.csv",
    "post_et_csv": "", "log_level": "", "cathy_backend": "",
    **{f.replace("-", "_"): "" for f in KS_FIRE_FLAGS},
}
LIST_DEFAULTS = {"dem_plot": ["2"], "pmin": ["-1e35"], "zroot_factors": []}


@dataclass
class BatchConfig:
    opts: dict = field(default_factory=lambda: dict(DEFAULTS))
    lists: dict = field(default_factory=lambda: {k: list(v) for k, v in LIST_DEFAULTS.items()})
    flags: dict = field(default_factory=lambda: {v: False for v in BOOLS.values()})

    def __getitem__(self, k):
        return self.opts[k] if k in self.opts else self.lists[k] if k in self.lists else self.flags[k]


@dataclass(frozen=True)
class Run:
    plot: str
    pmin: str
    outlet: str
    bottom: str
    spinup: str
    ic: str
    outlet_flux: str
    bottom_flux: str
    spinup_str: str          # string the driver logs: "off" | "on_<N>cycles"

    @property
    def label(self) -> str:
        return (f"dem-plot {self.plot} | pmin={self.pmin} | outlet={self.outlet} | "
                f"bottom={self.bottom} | spinup={self.spinup} | ic={self.ic}")


# ── parsing ──────────────────────────────────────────────────────────────────
def parse_args(argv) -> BatchConfig:
    """Shell-style parsing: list options take every following token up to the next ``--flag``."""
    cfg = BatchConfig()
    toks = [str(t) for t in argv]
    i = 0
    while i < len(toks):
        tok = toks[i]
        if tok in ("-h", "--help"):
            print(__doc__)
            raise SystemExit(0)
        if not tok.startswith("--"):
            raise SystemExit(f"ERROR: unexpected value '{tok}' (values must follow their --flag)")
        name, eq, inline = tok[2:].partition("=")
        i += 1
        if name in BOOLS:
            cfg.flags[BOOLS[name]] = True
        elif name in SCALARS:
            if eq:
                val = inline
            elif i < len(toks):
                val, i = toks[i], i + 1
            else:
                raise SystemExit(f"ERROR: --{name} needs a value")
            cfg.opts[SCALARS[name]] = val
        elif name in LISTS:
            vals = [inline] if eq else []
            while i < len(toks) and not toks[i].startswith("--"):
                vals.append(toks[i])
                i += 1
            if not vals:
                raise SystemExit(f"ERROR: --{name} needs at least one value")
            cfg.lists[LISTS[name]] = vals
        else:
            raise SystemExit(f"ERROR: unknown option '--{name}'")
    validate(cfg)
    return cfg


def _num(name, val):
    if not NUM.match(str(val)):
        raise SystemExit(f"ERROR: --{name} must be numeric (got '{val}')")


def validate(cfg: BatchConfig) -> None:
    o, l = cfg.opts, cfg.lists
    if o["forcing"] not in ("era5", "etp", "etp_ts"):
        raise SystemExit(f"ERROR: --forcing must be era5|etp|etp_ts (got '{o['forcing']}')")
    if o["etp_ts_column"] not in ("ET_0_wet", "ET_0_dry"):
        raise SystemExit(f"ERROR: --etp-ts-column must be ET_0_wet or ET_0_dry "
                         f"(got '{o['etp_ts_column']}')")
    if not re.match(r"^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$", o["dem_smooth_sigma"]):
        raise SystemExit(f"ERROR: --dem-smooth-sigma must be a non-negative number "
                         f"(got '{o['dem_smooth_sigma']}')")
    if o["etp_ts_feature_id"] and not o["etp_ts_feature_id"].isdigit():
        raise SystemExit(f"ERROR: --etp-ts-feature-id must be a non-negative integer "
                         f"(got '{o['etp_ts_feature_id']}')")
    for key, allowed in (("outlet", ("free", "closed", "both")),
                         ("bottom", ("free", "closed", "both")),
                         ("spinup", ("on", "off", "both")), ("ic", ("wt", "pressure", "both"))):
        if o[key] not in allowed:
            raise SystemExit(f"ERROR: --{key} must be {'|'.join(allowed)} (got '{o[key]}')")
    for p in l["dem_plot"]:
        if not (p.isdigit() and 1 <= int(p) <= 13):
            raise SystemExit(f"ERROR: --dem-plot values must be 1-13 (got '{p}')")
    for p in l["pmin"]:
        _num("pmin", p)
    if o["scf"]:
        _num("scf", o["scf"])
        if not 0 <= float(o["scf"]) <= 1:
            raise SystemExit(f"ERROR: --scf must be between 0 and 1 (got '{o['scf']}')")
    if o["feddes_crop"] and o["feddes_crop"] not in ("Corn", "Soybeans"):
        raise SystemExit(f"ERROR: --feddes-crop must be Corn or Soybeans (got '{o['feddes_crop']}')")
    for k in ("pcana", "pcref", "pcwlt", "pz", "omgc", "permx", "permx_bottom", "ks_layer_depth"):
        if o[k]:
            _num(k.replace("_", "-"), o[k])
    if l["zroot_factors"] and len(l["zroot_factors"]) != 5:
        raise SystemExit("ERROR: --zroot-factors needs exactly 5 values "
                         f"(bare sparse moderate dense very_dense), got {len(l['zroot_factors'])}")
    if bool(o["permx_bottom"]) != bool(o["ks_layer_depth"]):
        raise SystemExit("ERROR: --permx-bottom and --ks-layer-depth must be given together")


# ── expansion and commands ───────────────────────────────────────────────────
def _labels(mode: str, both: tuple) -> tuple:
    return both if mode == "both" else (mode,)


def expand_runs(cfg: BatchConfig) -> list[Run]:
    """One Run per combination, in the order of the old shell loops (plot, pmin, outlet,
    bottom, spin-up, ic)."""
    o = cfg.opts
    out = []
    for plot, pmin, outlet, bottom, spin, ic in itertools.product(
            cfg.lists["dem_plot"], cfg.lists["pmin"],
            _labels(o["outlet"], ("free", "closed")), _labels(o["bottom"], ("free", "closed")),
            _labels(o["spinup"], ("on", "off")), _labels(o["ic"], ("wt", "pressure"))):
        out.append(Run(plot, pmin, outlet, bottom, spin, ic,
                       o["outlet_flux"] if outlet == "free" else "0",
                       o["bottom_flux"] if bottom == "free" else "0",
                       f"on_{o['spinup_cycles']}cycles" if spin == "on" else "off"))
    return out


def driver_command(cfg: BatchConfig, run: Run) -> list[str]:
    """Driver call of one run; same options and order as the shell runner."""
    o = cfg.opts
    c = [sys.executable, str(HERE / DRIVER),
         "--start-year", o["start_year"], "--end-year", o["end_year"]]
    c += ["--spinup-cycles", o["spinup_cycles"]] if run.spinup == "on" else ["--no-spinup"]
    c += [f"--outlet-flux={run.outlet_flux}", f"--bottom-flux={run.bottom_flux}"]
    if run.ic == "wt":
        c += ["--ic-mode", "wt", f"--ic-wtposition={o['ic_wtposition']}"]
    else:
        c += ["--ic-mode", "pressure", f"--ic-pressure-head={o['ic_pressure_head']}"]
    if cfg.lists["zroot_factors"]:
        c += ["--zroot-factors", *cfg.lists["zroot_factors"]]
    if o["zroot"] != DEFAULTS["zroot"]:
        c += [f"--zroot={o['zroot']}"]
    if o["permx"]:
        c += [f"--permx={o['permx']}", f"--permy-ratio={o['permy_ratio']}",
              f"--permz-ratio={o['permz_ratio']}"]
    if o["permx_bottom"]:
        c += [f"--permx-bottom={o['permx_bottom']}", f"--ks-layer-depth={o['ks_layer_depth']}"]
    c += ["--scf", o["scf"]]
    c += [f"--{f}={o[f.replace('-', '_')]}" for f in KS_FIRE_FLAGS if o[f.replace("-", "_")]]
    if o["feddes_crop"]:
        c += ["--feddes-crop", o["feddes_crop"]]
    c += [f"--{k}={o[k]}" for k in ("pcana", "pcref", "pcwlt", "pz", "omgc") if o[k]]
    c += ["--forcing", o["forcing"], "--etp-ts-column", o["etp_ts_column"],
          "--dem-smooth-sigma", o["dem_smooth_sigma"]]
    if o["etp_ts_feature_id"]:
        c += [f"--etp-ts-feature-id={o['etp_ts_feature_id']}"]
    c += ["--dem-plot", run.plot, f"--pmin={run.pmin}"]
    if o["path2prj"] != DEFAULTS["path2prj"]:
        c += ["--path2prj", o["path2prj"]]
    if o["log_file"] != DEFAULTS["log_file"]:
        c += ["--log-file", o["log_file"]]
    c += [f"--{k.replace('_', '-')}={o[k]}" for k in ("post_et_csv", "log_level", "cathy_backend")
          if o[k]]
    c += [f"--{k.replace('_', '-')}" for k in ("post_plots", "no_plots", "debug") if cfg.flags[k]]
    return c


# ── skip-if-done ─────────────────────────────────────────────────────────────
def _ks_key(cfg: BatchConfig) -> str:
    """agramon_ks.ks_key of the batch's Ks options (post-fire modes + two-layer Ks)."""
    import argparse
    import agramon_ks as AK
    p = argparse.ArgumentParser()
    AK.KsFire.add_ks_fire_args(p)
    p.set_defaults(ks_fire_date=AK.FIRE_DATE_DEFAULT)
    ns = p.parse_args([f"--{f}={cfg.opts[f.replace('-', '_')]}" for f in KS_FIRE_FLAGS
                       if cfg.opts[f.replace("-", "_")]])
    pb, ld = cfg.opts["permx_bottom"], cfg.opts["ks_layer_depth"]
    ns.permx_bottom = float(pb) if pb else None
    ns.ks_layer_depth = float(ld) if ld else None
    return AK.ks_key(ns)


def check_already_done(cfg: BatchConfig, run: Run) -> str:
    """'DONE <folder>' or 'MISSING (<why>)'. Port of the shell runner's check_already_done:
    the run must be in the CSV log with every parameter that defines its identity AND its
    scenario folder must hold the three NetCDFs written only at the end of a finished run."""
    import pandas as pd
    o = cfg.opts
    log_path = (HERE / o["log_file"]).resolve() if not Path(o["log_file"]).is_absolute() \
        else Path(o["log_file"])
    if not log_path.exists():
        return "MISSING (no log file yet)"
    df = pd.read_csv(log_path)

    def close(series, value):
        return series.astype(float).apply(lambda v: math.isclose(v, value, rel_tol=1e-9, abs_tol=1e-9))

    mask = ((df["start_year"] == int(o["start_year"])) & (df["end_year"] == int(o["end_year"]))
            & close(df["ZROOT"], float(o["zroot"])) & close(df["PMIN"], float(run.pmin))
            & (df["watershed_nb"] == int(run.plot))
            & close(df["outlet_flux"], float(run.outlet_flux))
            & close(df["bottom_flux"], float(run.bottom_flux))
            & (df["spinup"].astype(str) == run.spinup_str))
    if "withLAI" in df.columns:
        mask &= df["withLAI"] == 1
    forcing = o["forcing"]
    if "forcing_mode" in df.columns:        # rows without it are satellite-ETp runs
        mask &= df["forcing_mode"].fillna("etp").astype(str) == forcing
        if forcing == "etp_ts":
            fid = int(o["etp_ts_feature_id"]) if o["etp_ts_feature_id"] else int(run.plot) - 1
            if "etp_ts_column" in df.columns:
                mask &= df["etp_ts_column"].fillna("").astype(str) == o["etp_ts_column"]
            if "etp_ts_feature_id" in df.columns:
                mask &= pd.to_numeric(df["etp_ts_feature_id"], errors="coerce") == fid
    elif forcing != "etp":
        mask &= False
    if "ic_mode" in df.columns:
        mask &= df["ic_mode"].astype(str) == run.ic
    if cfg.lists["zroot_factors"] and "zroot_factors" in df.columns:
        # the driver logs them with "{v:g}" ("2", not "2.0"): normalise the same way, otherwise
        # factors typed as 2.0 never match (the old shell check had this bug)
        mask &= df["zroot_factors"].astype(str) == ",".join(f"{float(v):g}" for v in cfg.lists["zroot_factors"])
    if "permx" in df.columns:
        if o["permx"]:
            mask &= close(df["permx"], float(o["permx"]))
            for col, key in (("permy_ratio", "permy_ratio"), ("permz_ratio", "permz_ratio")):
                if col in df.columns:
                    mask &= close(df[col], float(o[key]))
        else:
            mask &= df["permx"].isna()
    if "SCF" in df.columns:
        mask &= close(df["SCF"], float(o["scf"]) if o["scf"] else 1.0)
    expected_ks = _ks_key(cfg)
    if "ks_fire" in df.columns:
        mask &= df["ks_fire"].fillna("off").astype(str) == expected_ks
    elif expected_ks != "off":
        mask &= False

    lit = {}
    if o["feddes_crop"]:
        try:
            from pyCATHY.petro import Feddes_litterature
            f = Feddes_litterature(crop=o["feddes_crop"])
            lit = {c: float(f[c]) for c in ("PCANA", "PCREF", "PCWLT")}
        except Exception:
            pass                             # pyCATHY not importable here: fall back to isna
    for col in ("PCANA", "PCREF", "PCWLT", "PZ", "OMGC"):
        if col not in df.columns:
            continue
        explicit = o[col.lower()]
        if explicit:
            mask &= close(df[col], float(explicit))
        elif col in lit:
            mask &= close(df[col], lit[col])
        else:
            mask &= df[col].isna()
    if "feddes_crop" in df.columns:
        mask &= (df["feddes_crop"].astype(str) == o["feddes_crop"]) if o["feddes_crop"] \
            else df["feddes_crop"].isna()

    hits = df[mask]
    if hits.empty:
        return "MISSING (no matching logged run)"
    sim_index = int(hits.iloc[-1]["sim_index"])
    prj = Path(o["path2prj"])
    out_root = (prj if prj.is_absolute() else HERE / prj) / "outputs"
    required = ("et_output.nc", "psi_output.nc", "sw_output.nc")
    for d in sorted(out_root.glob(f"scenario_{sim_index}_*")) + [out_root / f"scenario_{sim_index}"]:
        if d.is_dir() and all((d / f).exists() for f in required):
            return f"DONE {d}"
    return f"MISSING (sim_index={sim_index} logged but outputs incomplete)"


def skip_status(cfg: BatchConfig, run: Run) -> str | None:
    """None when the run must be (re)run, else the folder of the finished scenario.
    The check is off with --force and with --dem-smooth-sigma != 0 (the log has no sigma column,
    so a smoothed run could not be told apart from an unsmoothed one)."""
    if cfg.flags["force"] or float(cfg.opts["dem_smooth_sigma"]) != 0:
        return None
    status = check_already_done(cfg, run)
    return status[5:] if status.startswith("DONE") else None


# ── main ─────────────────────────────────────────────────────────────────────
def describe(cfg: BatchConfig, runs: list) -> None:
    o, l = cfg.opts, cfg.lists
    print(f"Started       : {time.strftime('%c')}")
    print(f"dem-plot(s)   : {' '.join(l['dem_plot'])}")
    print(f"PMIN value(s) : {' '.join(l['pmin'])} m")
    print(f"outlet / bottom / spinup / ic : {o['outlet']} / {o['bottom']} / {o['spinup']} / {o['ic']}")
    print(f"forcing       : {o['forcing']}" + (f" (column={o['etp_ts_column']}, feature_id="
                                              f"{o['etp_ts_feature_id'] or 'auto from plot'})"
                                              if o["forcing"] == "etp_ts" else ""))
    print(f"zroot-factors : {' '.join(l['zroot_factors']) or '(driver default)'}")
    print(f"PERMX (Ks)    : {o['permx'] + ' m/s' if o['permx'] else '(pyCATHY default)'}")
    print(f"SCF           : {o['scf']}")
    print(f"skip-existing : {'disabled (--force)' if cfg.flags['force'] else 'enabled'}")
    print(f"runs          : {len(runs)}\n")


def main(argv=None) -> int:
    cfg = parse_args(sys.argv[1:] if argv is None else argv)
    runs = expand_runs(cfg)
    describe(cfg, runs)
    failed, skipped = [], 0
    for run in runs:
        print("=" * 60 + f"\n  {run.label} ...  {time.strftime('%c')}\n" + "=" * 60, flush=True)
        try:
            done = skip_status(cfg, run)
        except ImportError as exc:           # pandas / agramon_ks missing: run everything
            print(f"  (skip check unavailable: {exc})")
            done = None
        if done:
            skipped += 1
            print(f"  skip: already completed ({done}). Use --force to re-run.\n")
            continue
        cmd = driver_command(cfg, run)
        print("  " + " ".join(shlex.quote(c) for c in cmd), flush=True)
        if cfg.flags["dry_run"]:
            continue
        rc = subprocess.run(cmd, cwd=HERE).returncode
        if rc:
            failed.append((run.label, rc))
            print(f"!! {run.label} failed (rc={rc})", flush=True)
            if cfg.flags["stop_on_error"]:
                break
    print(f"\nFinished: {time.strftime('%c')}  ({len(runs)} run(s), {skipped} skipped, "
          f"{len(failed)} failed)")
    for label, rc in failed:
        print(f"  failed: {label} (rc={rc})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
