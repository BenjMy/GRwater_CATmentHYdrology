"""
agramon_runlog.py — Run logging and scenario naming
===================================================
GRwater project — Agramon study site

  * ``MarkdownLog`` — per-scenario markdown run log.
  * ``log_simulation`` — shared CSV log (assigns the sim_index).
  * ``build_scenario_dirname`` — self-describing scenario folder name.

Moved verbatim out of Agramon_withLAI_withETp.py (REFACTORING.md, step 2.1).
"""

from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from geoutils import uCATHY


class MarkdownLog:
    """
    Lightweight, append-as-you-go markdown logger for a single simulation
    run — one narrative file mirroring the console output (run config,
    per-step progress, spin-up diagnostics, and a per-month LAI/ETp/solver
    table), so a run can be reviewed later without re-parsing stdout.

    Writes to disk on every call (re-opened in append mode each time)
    rather than buffering in memory, so a crash mid-run still leaves a
    readable partial log instead of losing everything unwritten.
    """

    def __init__(self, path: Path, title: str) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._n_warnings = 0
        self._start_time = datetime.now()
        self._metrics: dict[str, float] = {}
        try:
            with open(self.path, "w") as fh:
                fh.write(f"# {title}\n\n")
                fh.write(f"_Started: {datetime.now():%Y-%m-%d %H:%M:%S}_\n\n")
        except OSError as exc:
            # Logging failures should never abort the simulation itself.
            print(f"  WARNING: could not open markdown log at {self.path} ({exc}) "
                  f"— continuing without it.")

    def _write(self, text: str) -> None:
        try:
            with open(self.path, "a") as fh:
                fh.write(text)
        except OSError as exc:
            print(f"  WARNING: markdown log write failed ({exc}).")

    def section(self, title: str) -> None:
        self._write(f"\n## {title}\n\n")

    def subsection(self, title: str) -> None:
        self._write(f"\n### {title}\n\n")

    def meta(self, items: dict) -> None:
        """Key/value bullet list — run parameters, step summaries, etc."""
        for k, v in items.items():
            self._write(f"- **{k}**: {v}\n")

    def line(self, text: str = "") -> None:
        self._write(f"{text}\n")

    def bullet(self, text: str) -> None:
        self._write(f"- {text}\n")

    def warn(self, text: str) -> None:
        self._n_warnings += 1
        self._write(f"- ⚠️ **WARNING**: {text}\n")

    def table_header(self, headers: list) -> None:
        self._write("\n| " + " | ".join(headers) + " |\n")
        self._write("|" + "|".join(["---"] * len(headers)) + "|\n")

    def table_row(self, cells: list) -> None:
        self._write("| " + " | ".join(str(c) for c in cells) + " |\n")

    def fixes_applied(self, fixes: list) -> None:
        """
        Document known pipeline patches active for this run — written once,
        near the top of the log, so a run's markdown self-describes which
        fixes were in effect without needing to diff the script against git.
        """
        self.section("Known fixes applied")
        for fx in fixes:
            self.bullet(fx)

    def add_metric(self, name: str, value: float) -> None:
        """
        Accumulate a running numeric total across the whole run (e.g. cells
        patched per month) without the caller needing to track it manually.
        Automatically rolled into the summary written by close().
        """
        self._metrics[name] = self._metrics.get(name, 0) + value

    def close(self, summary: dict | None = None) -> None:
        elapsed = datetime.now() - self._start_time
        self.section("Run summary")
        merged = dict(summary) if summary else {}
        merged.update(self._metrics)
        merged["Warnings logged"] = self._n_warnings
        merged["Elapsed"] = str(elapsed).split(".")[0]  # drop microseconds
        self.meta(merged)
        self._write(f"\n---\n\n_Finished: {datetime.now():%Y-%m-%d %H:%M:%S}_\n")


def log_simulation(log_path: Path, start_year: int, end_year: int,
                   zroot: float, pmin: float, watershed_nb: int,
                   with_lai: int, outlet_flux: float = None,
                   bottom_flux: float = None, spinup: str = None,
                   ic_mode: str = None, ic_wtposition: float = None,
                   ic_pressure_head: float = None,
                   zroot_factors: list = None,
                   permx: float = None,
                   permy_ratio: float = None,
                   permz_ratio: float = None,
                   scf: float = None,
                   pcana: float = None,
                   pcref: float = None,
                   pcwlt: float = None,
                   pz: float = None,
                   omgc: float = None,
                   feddes_crop: str = None,
                   forcing_mode: str = None,
                   etp_ts_column: str = None,
                   etp_ts_feature_id: int = None,
                   ks_fire: str = None) -> int:
    """Register parameters in the shared CSV log and return the sim_index.

    outlet_flux / bottom_flux / spinup / ic_* are recorded too (not just
    used to tag the folder name) so `--list-scenarios` and the CSV itself
    are enough to compare runs without having to parse folder names — and
    so batch runners can tell whether a given parameter combination
    (including initial condition) has already been logged, e.g. to skip
    re-running scenarios that already completed successfully.

    zroot_factors / permx / permy_ratio / permz_ratio are logged the same
    way — a run's root-depth mapping and Ks override (incl. anisotropy)
    are as much a part of its identity as ZROOT or PMIN, so they belong
    in the CSV rather than only ever showing up in the console/markdown
    log.

    pcana/pcref/pcwlt/pz/omgc/feddes_crop are the remaining Feddes (FP_map)
    parameters, logged the same way for the same reason — resolved once via
    ``feddes_params_from_args(args)`` at the call site, not recomputed here.

    ``scf`` (normally just ``args.scf``) is logged too, since it changes
    which ET-partitioning approach was actually used.

    ``forcing_mode`` ("era5" | "etp" | "etp_ts") is part of a run's identity
    too. ``etp_ts_column`` / ``etp_ts_feature_id`` are only logged for
    "etp_ts" runs (None otherwise). Rows logged before these columns existed
    were all satellite-ETp runs, i.e. an empty/NaN forcing_mode means "etp".
    """
    params = {
        "start_year":   start_year,
        "end_year":     end_year,
        "ZROOT":        zroot,
        "PMIN":         pmin,
        "SCF":          scf,
        "watershed_nb": watershed_nb,
        "withLAI":      with_lai,
        "outlet_flux":  outlet_flux,
        "bottom_flux":  bottom_flux,
        "spinup":       spinup,
        "ic_mode":          ic_mode,
        "ic_wtposition":    ic_wtposition,
        "ic_pressure_head": ic_pressure_head,
        "zroot_factors": (",".join(f"{v:g}" for v in zroot_factors)
                           if zroot_factors is not None else None),
        "permx":        permx,
        "permy_ratio":  permy_ratio if permx is not None else None,
        "permz_ratio":  permz_ratio if permx is not None else None,
        "PCANA":        pcana,
        "PCREF":        pcref,
        "PCWLT":        pcwlt,
        "PZ":           pz,
        "OMGC":         omgc,
        "feddes_crop":  feddes_crop,
        "forcing_mode":      forcing_mode,
        "etp_ts_column":     etp_ts_column if forcing_mode == "etp_ts" else None,
        "etp_ts_feature_id": etp_ts_feature_id if forcing_mode == "etp_ts" else None,
        # post-fire Ks setup (agramon_ks.ks_key); rows logged before this column
        # existed were all static-Ks runs, i.e. an empty value means "off"
        "ks_fire":           ks_fire,
    }
    sim_index = uCATHY.log_simulation(log_path, params)
    print(f"  sim_index = {sim_index}  (withLAI={with_lai})")
    return sim_index


def build_scenario_dirname(sim_index: int, dem_plot: int,
                           outlet_flux: float, bottom_flux: float,
                           no_spinup: bool, spinup_cycles: int,
                           pmin: float, ic_mode: str = "wt",
                           forcing_tag: str = "") -> str:
    """
    Human-readable, self-describing scenario folder name.

    ``scenario_<sim_index>`` stays the root/prefix of the name — every
    place that keys off the CSV-assigned sim_index (post-processing
    scripts, --scenario N, log lookups) keeps working unchanged — but a
    short tag is appended describing the axes that actually differ
    between batch runs (dem-plot, outlet/bottom BC, spin-up, PMIN,
    initial condition), so scenarios run with different args can be told
    apart, and compared, just by reading folder names side by side
    instead of cross-referencing the CSV log for every one of them.

    ``ic_mode`` is included so runs that only differ by --ic-mode (wt vs
    pressure) get distinct folders/tags instead of silently colliding —
    previously two such runs shared the same tag (differing only by the
    sim_index prefix), which also made them indistinguishable to any
    "has this combination already run?" check based on the tag alone.

    ``forcing_tag`` is appended only for non-default forcings (currently the
    "etp_ts" runs, e.g. ``forcing-etpts-wet-fid1``), so legacy satellite-ETp
    folder names stay exactly as they were.

    e.g. scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c_pmin-3.0_ic-wt
    """
    outlet_label = "closed" if outlet_flux == 0 else "free"
    bottom_label = "closed" if bottom_flux == 0 else "free"
    spinup_label = "off" if no_spinup else f"on{spinup_cycles}c"
    pmin_label = f"{pmin:g}".replace("-", "m").replace("+", "").replace(".", "p")
    tag = (f"plot{dem_plot}"
           f"_outlet-{outlet_label}"
           f"_bottom-{bottom_label}"
           f"_spinup-{spinup_label}"
           f"_pmin-{pmin_label}"
           f"_ic-{ic_mode}")
    if forcing_tag:
        tag += f"_{forcing_tag}"
    return f"scenario_{sim_index}_{tag}"
