"""
hypotheses_common.py - engine shared by hyp1_layers_Ks.py ... hyp4_SW.py
========================================================================
A *hypothesis* is a set of runs (variants) of Agramon_withLAI_withETp.py that differ in one
aspect (layers/Ks, SCF/PMIN, initial condition, ETp forcing). It is described in
``hypotheses.toml`` and handled in four steps, all inside ONE study folder
``<out_root>/<hyp>/`` (prj/, simulation_log.csv, manifest.json, metrics.csv, figures/):

    plan      list the runs (``--dry-run`` also draws the variant matrix)
    run       one driver call per variant, skipped when its et_output.nc already exists
    validate  WB ETa vs the EO energy-balance NetCDF (analysis_ET_compare.validate_against_eb)
    plot      hypothesis figures (plot_hypotheses.py)

and, independent of the EO reference, ``balance_study``: annual / post-fire water-balance totals
(ETa, recharge, runoff, closure), per-variant means and (when a control run is configured) the
recovery metrics, all from hypothesis_tables.py, plus a ``study_report.md`` listing what failed or
was skipped.

Heavy imports (analysis_ET_compare: geopandas/rioxarray, matplotlib) are done lazily, so planning
and ``--dry-run`` need only the standard library (+ matplotlib for the figure).
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "hypotheses.toml"
MANIFEST = "manifest.json"
RESERVED = ("name", "group")


# ── configuration ────────────────────────────────────────────────────────────
@dataclass
class Variant:
    name: str
    args: dict                       # merged driver options (common < base < variant)
    group: str = ""
    varying: tuple = ()              # option keys that differ between the variants of the hypothesis


@dataclass
class Hypothesis:
    key: str
    description: str
    variants: list
    validation: dict = field(default_factory=dict)
    driver: str = "Agramon_withLAI_withETp.py"
    out_root: Path | None = None


def _read_toml(path: Path) -> dict:
    from agramon_config import _read_toml as read
    return read(path)


def _short(value) -> str:
    return ",".join(map(str, value)) if isinstance(value, list) else str(value)


def _expand(key: str, spec: dict, common: dict) -> list[Variant]:
    base = {**common, **spec.get("base", {})}
    out: list[Variant] = []
    grid = spec.get("grid", {})
    for k, v in grid.items():
        if not isinstance(v, list) or not v:
            raise ValueError(f"{key}.grid.{k} must be a non-empty list")
    if grid:
        keys = list(grid)
        for combo in itertools.product(*(grid[k] for k in keys)):
            name = "_".join(f"{k}{_short(v)}" for k, v in zip(keys, combo))
            out.append(Variant(name, {**base, **dict(zip(keys, combo))}, group="", varying=tuple(keys)))
    for var in spec.get("variants", []):
        if "name" not in var:
            raise ValueError(f"{key}: every variant needs a name (got {var})")
        opts = {k: v for k, v in var.items() if k not in RESERVED}
        out.append(Variant(str(var["name"]), {**base, **opts}, str(var.get("group", ""))))
    if not out:
        raise ValueError(f"{key}: needs a [grid] table or a variants list")
    names = [v.name for v in out]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise ValueError(f"{key}: duplicate variant names {sorted(dup)}")
    if not grid:      # explicit variants: the varying options are those that differ
        keys = sorted({k for v in out for k in v.args if len({_short(w.args.get(k)) for w in out}) > 1})
        for v in out:
            v.varying = tuple(keys)
    return out


def load_hypotheses(path: Path = DEFAULT_CONFIG) -> dict[str, Hypothesis]:
    data = _read_toml(Path(path))
    common = data.get("common", {})
    general = data.get("general", {})
    val = data.get("validation", {})
    skip = {"general", "validation", "common"}
    out = {}
    for key, spec in data.items():
        if key in skip:
            continue
        out[key] = Hypothesis(key, str(spec.get("description", "")), _expand(key, spec, common),
                              dict(val), general.get("driver", "Agramon_withLAI_withETp.py"),
                              Path(general["out_root"]).expanduser() if general.get("out_root") else None)
    return out


def default_out_root() -> Path:
    from agramon_config import WB_PRJ_ROOT
    return Path(WB_PRJ_ROOT) / "hypotheses"


def resolve_study(arg, config: Path = DEFAULT_CONFIG) -> Path:
    """
    Study folder from ``--study``: either a path, or just the hypothesis name / prefix (``hyp4``,
    ``hyp4_SW``), looked up under the hypotheses out_root of the config (hypotheses.toml
    ``[general] out_root``, else ``<wb_prj_root>/hypotheses`` from data_paths.toml). An existing folder,
    or anything that looks like a path, is used as given. Exits with the list of choices if the name
    matches nothing or several studies.
    """
    text = str(arg)
    path = Path(text).expanduser()
    if path.is_dir() or any(c in text for c in "/\\") or text.startswith((".", "~")):
        return path
    try:
        hyps = load_hypotheses(config)
    except Exception:                                  # unreadable config: fall back to the folders on disk
        hyps = {}
    root = Path(next((h.out_root for h in hyps.values() if h.out_root), None) or default_out_root())
    on_disk = [d.name for d in sorted(root.iterdir()) if d.is_dir() and (d / "simulation_log.csv").exists()] \
        if root.is_dir() else []
    names = list(dict.fromkeys(on_disk + list(hyps)))
    hits = [n for n in names if n == text] or [n for n in names if n.startswith(text)]
    if len(hits) == 1:
        return root / hits[0]
    choices = ", ".join(names) or "none found"
    raise SystemExit(f"--study {text!r}: " + ("ambiguous" if hits else "no study of that name")
                     + f" under {root}\n  available: {choices}")


def study_dir(h: Hypothesis, out_root: Path | None = None) -> Path:
    return Path(out_root or h.out_root or default_out_root()) / h.key


# ── commands ─────────────────────────────────────────────────────────────────
def _flag(key: str, value) -> list[str]:
    if value is True:
        return [f"--{key}"]
    if value is False or value is None:
        return []
    vals = [str(v) for v in (value if isinstance(value, list) else [value])]
    if len(vals) == 1 and vals[0].startswith("-"):      # "--pmin=-1e35": argparse may read "-1e35" as a flag
        return [f"--{key}={vals[0]}"]
    return [f"--{key}", *vals]


def variant_command(h: Hypothesis, v: Variant, study: Path, extra=(), post_plots=False) -> list[str]:
    """Driver call of one variant; every run of the hypothesis shares prj/ and the log CSV."""
    cmd = [sys.executable, str(HERE / h.driver)]
    for k, val in v.args.items():
        cmd += _flag(k, val)
    cmd += ["--path2prj", str(study / "prj"), "--log-file", str(study / "simulation_log.csv")]
    if post_plots:
        cmd += ["--post-plots"]
    else:
        cmd += ["--no-plots"]
    return cmd + list(extra)


def plan_table(h: Hypothesis):
    """Tidy table (list of dicts): one row per variant with its varying options."""
    rows = []
    for v in h.variants:
        row = {"variant": v.name, "group": v.group}
        row.update({k: _short(v.args.get(k, "-")) for k in v.varying})
        rows.append(row)
    return rows


def print_plan(h: Hypothesis, study: Path, extra=(), post_plots=False, show_cmd=True) -> None:
    rows = plan_table(h)
    cols = list(rows[0])
    w = {c: max(len(c), *(len(str(r[c])) for r in rows)) for c in cols}
    print(f"\n=== {h.key}: {h.description}\n    study folder: {study}\n    {len(rows)} run(s)")
    print("    " + "  ".join(c.ljust(w[c]) for c in cols))
    for r in rows:
        print("    " + "  ".join(str(r[c]).ljust(w[c]) for c in cols))
    if show_cmd:
        for v in h.variants:
            print(f"\n# {v.name}\n" + " ".join(variant_command(h, v, study, extra, post_plots)))


# ── manifest / run ───────────────────────────────────────────────────────────
def load_manifest(study: Path) -> dict:
    p = study / MANIFEST
    return json.loads(p.read_text()) if p.exists() else {}


def save_manifest(study: Path, manifest: dict) -> None:
    study.mkdir(parents=True, exist_ok=True)
    (study / MANIFEST).write_text(json.dumps(manifest, indent=2, default=str))


def entry_sim_dir(study: Path, entry: dict):
    """Scenario folder of a manifest entry, or None. A relative sim_dir is read against
    <study>/prj; an absolute one (manifests written before the registry) is used as stored."""
    sd = (entry or {}).get("sim_dir")
    if not sd:
        return None
    p = Path(sd)
    return p if p.is_absolute() else study / "prj" / p


def _scenario_dirs(prj: Path) -> set:
    outs = prj / "outputs"
    return set(outs.glob("scenario_*")) if outs.exists() else set()


def _run_tee(cmd, cwd, env, keep: int = 15) -> tuple:
    """Run ``cmd``; stdout is inherited, stderr is echoed live and its last ``keep`` lines are
    returned as one line: ``(returncode, stderr_tail)`` (what the manifest keeps for a failed run)."""
    import collections
    tail = collections.deque(maxlen=keep)
    p = subprocess.Popen(cmd, cwd=cwd, env=env, stderr=subprocess.PIPE, text=True, errors="replace")
    for line in p.stderr:
        sys.stderr.write(line)
        sys.stderr.flush()
        tail.append(line.rstrip())
    return p.wait(), " | ".join(t for t in tail if t)[-500:]


def run_variants(h: Hypothesis, study: Path, *, only=None, force=False, extra=(),
                 post_plots=False, stop_on_error=False) -> dict:
    """Run the variants (skipping finished ones); returns the updated manifest."""
    manifest = load_manifest(study)
    study.mkdir(parents=True, exist_ok=True)
    for v in h.variants:
        if only and v.name not in only:
            continue
        done = manifest.get(v.name, {})
        sim_dir = entry_sim_dir(study, done)
        if not force and done.get("status") == "ok" and sim_dir and (sim_dir / "et_output.nc").exists():
            print(f"[{h.key}] {v.name}: already done -> {sim_dir.name} (use --force to rerun)")
            continue
        cmd = variant_command(h, v, study, extra, post_plots)
        print(f"\n[{h.key}] {v.name}\n  " + " ".join(cmd), flush=True)
        before = _scenario_dirs(study / "prj")
        t0 = time.time()
        # the driver stamps these in its registry row (agramon_registry.start_run)
        env = {**os.environ, "AGRAMON_HYPOTHESIS": h.key, "AGRAMON_VARIANT": v.name}
        rc, err_tail = _run_tee(cmd, HERE, env)
        new = sorted(_scenario_dirs(study / "prj") - before, key=lambda p: p.stat().st_mtime)
        ok = rc == 0 and bool(new)
        manifest[v.name] = {"group": v.group, "args": v.args, "returncode": rc,
                            "status": "ok" if ok else "failed",
                            "sim_dir": new[-1].relative_to(study / "prj").as_posix() if new else None,
                            "seconds": round(time.time() - t0)}
        if not ok:
            manifest[v.name]["error_tail"] = err_tail or "no stderr (rc=0 but no new scenario folder)"
        save_manifest(study, manifest)
        if not ok:
            print(f"!! [{h.key}] {v.name} failed (rc={rc})", flush=True)
            if stop_on_error:
                break
    return manifest


# ── validation vs EO-EB ──────────────────────────────────────────────────────
def validate_variants(h: Hypothesis, study: Path, *, validator=None) -> tuple:
    """
    Score every finished variant against the EO energy-balance ETa.
    Returns (metrics DataFrame, {variant: monthly WB series}, {eo_var: monthly EO series}).
    ``validator(sim_dir) -> dict`` defaults to analysis_ET_compare.validate_against_eb with
    the shapefiles and EO file loaded once (they are slow to open).
    """
    import pandas as pd
    manifest = load_manifest(study)
    if validator is None:
        import analysis_ET_compare as AE
        gdf, ds_eo = AE.load_shapefiles(), AE.load_eo_et()
        full = bool(h.validation.get("full_figures", False))
        year = int(h.validation.get("selected_year", 2020))

        def validator(sim_dir):
            return AE.validate_against_eb(sim_dir, gdf_proj=gdf, ds_eo=ds_eo, full_figures=full,
                                          fig_dir=study / "figures" / "per_run" / Path(sim_dir).name,
                                          selected_year=year)
    rows, wb, eo = [], {}, {}
    for v in h.variants:
        e = manifest.get(v.name)
        if not e or e.get("status") != "ok" or not e.get("sim_dir"):
            print(f"[{h.key}] {v.name}: no finished run, skipped in validation")
            continue
        try:
            r = validator(entry_sim_dir(study, e))
        except Exception as exc:                          # one bad run must not stop the study
            print(f"[{h.key}] {v.name}: validation failed - {exc}")
            continue
        wb[v.name] = r["s_wb"]
        eo[r["eo_var"]] = r["s_eo"]
        rows.append({"variant": v.name, "group": v.group, "eo_var": r["eo_var"],
                     "condition": r.get("condition"), **r["metrics"]})
    df = pd.DataFrame(rows)
    study.mkdir(parents=True, exist_ok=True)
    df.to_csv(study / "metrics.csv", index=False)
    return df, wb, eo


# ── water balance, totals and recovery ───────────────────────────────────────
def _opt_float(value):
    """float(value), or None when unset ('' / None) or not a number."""
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def balance_study(h: Hypothesis, study: Path, *, loader=None) -> dict:
    """
    Water-balance tables of every finished variant (hypothesis_tables.py), written in the study folder:
    water_balance_monthly.csv, annual_totals.csv, period_totals.csv (all / pre_fire / post_fire),
    variant_summary.csv and, when [validation] control_sim_dir points to a control run,
    recovery_metrics.csv. A table with no data removes its stale file instead of leaving it behind.

    Returns {'monthly': {variant: frame}, 'tables': {annual, period, summary}, 'recovery': frame,
    'control': frame | None, 'notes': [what was skipped and why], 'written': [paths], 'fire_date',
    'tol', 'sustain'}. ``loader(sim_dir, area_m2) -> (frame, notes)`` defaults to
    hypothesis_tables.load_monthly_balance (tests pass a fake). Never raises on a bad run.
    """
    import pandas as pd
    import hypothesis_tables as HT
    v = h.validation
    fire = str(v.get("fire_date") or HT.default_fire_date())
    area = _opt_float(v.get("catchment_area_m2"))
    tol = _opt_float(v.get("recovery_tol"))
    tol = 0.1 if tol is None else tol
    sustain = int(v.get("recovery_sustain") or 3)
    load = loader or HT.load_monthly_balance
    manifest = load_manifest(study)
    groups = {x.name: x.group for x in h.variants}
    notes, monthly = [], {}
    for x in h.variants:
        e = manifest.get(x.name)
        if not e or e.get("status") != "ok" or not e.get("sim_dir"):
            notes.append(f"{x.name}: no finished run, left out of the water-balance tables")
            continue
        try:
            frame, frame_notes = load(entry_sim_dir(study, e), area)
        except Exception as exc:                          # one bad run must not stop the study
            notes.append(f"{x.name}: outputs unreadable ({exc})")
            continue
        notes += [f"{x.name}: {m}" for m in frame_notes]
        if frame is None or frame.empty:
            notes.append(f"{x.name}: empty monthly balance, left out of the tables")
            continue
        monthly[x.name] = frame
    out = {"monthly": monthly, "tables": {}, "recovery": pd.DataFrame(), "control": None,
           "notes": notes, "written": [], "fire_date": fire, "tol": tol, "sustain": sustain}
    if not monthly:
        return out
    control_dir = str(v.get("control_sim_dir") or "").strip()
    if control_dir:
        try:
            control, control_notes = load(Path(control_dir).expanduser(), area)
            notes += [f"control: {m}" for m in control_notes]
            if control is None or control.empty:
                notes.append(f"control: no monthly balance in {control_dir}")
            else:
                out["control"] = control
        except Exception as exc:
            notes.append(f"control: outputs unreadable ({exc})")
    out["recovery"], recovery_notes = HT.recovery_by_variant(monthly, out["control"], fire, tol, sustain)
    notes += recovery_notes
    out["tables"] = HT.tables_by_variant(monthly, fire, groups)
    study.mkdir(parents=True, exist_ok=True)
    long = pd.concat(monthly, names=["variant"]).reset_index()
    for name, frame in (("water_balance_monthly.csv", long), ("annual_totals.csv", out["tables"]["annual"]),
                        ("period_totals.csv", out["tables"]["period"]),
                        ("variant_summary.csv", out["tables"]["summary"]),
                        ("recovery_metrics.csv", out["recovery"])):
        path = study / name
        if len(frame):
            frame.to_csv(path, index=False)
            out["written"].append(path)
        else:
            path.unlink(missing_ok=True)
    return out


def write_report(h: Hypothesis, study: Path, manifest: dict, metrics=None, bal=None) -> Path:
    """study_report.md: run status, how many variants were scored, files written, and every issue
    (failed or missing run, finished but unscored, skipped table) so nothing is lost in the console."""
    names = [x.name for x in h.variants]
    ok = [n for n in names if manifest.get(n, {}).get("status") == "ok"]
    issues = []
    for n in names:
        e = manifest.get(n)
        if not e:
            issues.append(f"{n}: not run")
        elif e.get("status") != "ok":
            issues.append(f"{n}: run failed (returncode {e.get('returncode')})")
    lines = [f"# {h.key}", "", h.description, "", "## Runs", "", "| variant | status | seconds |", "|---|---|---|"]
    for n in names:
        e = manifest.get(n) or {}
        lines.append(f"| {n} | {e.get('status', 'not run')} | {e.get('seconds', '')} |")
    lines += ["", "## Validation vs EO-EB ETa", ""]
    if metrics is None:
        lines.append("not run (--no-validate)")
    else:
        scored = set(metrics["variant"]) if len(metrics) else set()
        lines.append(f"{len(scored)} of {len(ok)} finished variants scored")
        issues += [f"{n}: finished but not scored against EO-EB (validation failed or skipped)"
                   for n in ok if n not in scored]
    lines += ["", "## Water balance", ""]
    if bal is None:
        lines.append("not run (--no-balance)")
    else:
        lines.append(f"fire date {bal['fire_date']}; tables: " +
                     (", ".join(p.name for p in bal["written"]) or "none written"))
        issues += bal["notes"]
    lines += ["", "## Issues", ""] + ([f"- {m}" for m in issues] or ["none"])
    path = study / "study_report.md"
    study.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


# ── plots ────────────────────────────────────────────────────────────────────
def monthly_stats_by_variant(h: Hypothesis, study: Path) -> dict:
    import pandas as pd
    out = {}
    for name, e in load_manifest(study).items():
        p = entry_sim_dir(study, e) / "monthly_stats.csv" if e.get("sim_dir") else None
        if p is not None and p.exists():
            out[name] = pd.read_csv(p)
    return out


def _save(fn, *args, out, **kw):
    """Call a plot_hypotheses function, close the figure, return the written path (None if nothing drawn)."""
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
    fig = fn(*args, out=out, **kw)
    if fig is None:
        return None
    plt.close(fig)
    return out


def plot_study(h: Hypothesis, study: Path, metrics, wb: dict, eo: dict) -> list:
    """Hypothesis figures -> study/figures; returns the written paths."""
    import plot_hypotheses as PH
    fig_dir = study / "figures"
    groups = {v.name: v.group for v in h.variants}
    written = []
    if wb:
        written.append(_save(PH.plot_hypothesis_eta, wb, eo, groups=groups,
                             title=f"{h.key}: {h.description}", out=fig_dir / f"{h.key}_eta_vs_EB.png"))
    if len(metrics):
        written.append(_save(PH.plot_hypothesis_skill, metrics, title=f"{h.key}: skill vs EO-EB ETa",
                             out=fig_dir / f"{h.key}_skill.png"))
    stats = monthly_stats_by_variant(h, study)
    if stats:
        written.append(_save(PH.plot_hypothesis_response, stats, groups=groups,
                             title=f"{h.key}: hydrological response", out=fig_dir / f"{h.key}_response.png"))
    return [w for w in written if w]


def plot_balance(h: Hypothesis, study: Path, bal: dict) -> list:
    """Water-balance figures -> study/figures: annual totals, post-fire budget per variant and, when a
    control run was loaded, one burned-vs-control recovery figure per variant; returns the paths."""
    import plot_hypotheses as PH
    fig_dir = study / "figures"
    groups = {x.name: x.group for x in h.variants}
    tables = bal.get("tables") or {}
    written = []
    if len(tables.get("annual", ())):
        written.append(_save(PH.plot_hypothesis_annual, tables["annual"], groups=groups,
                             fire_date=bal["fire_date"], title=f"{h.key}: annual water-balance totals",
                             out=fig_dir / f"{h.key}_annual_totals.png"))
    if len(tables.get("period", ())):
        written.append(_save(PH.plot_hypothesis_balance, tables["period"], groups=groups,
                             title=f"{h.key}: water balance per variant (mm/yr)",
                             out=fig_dir / f"{h.key}_water_balance.png"))
    control, recovery = bal.get("control"), bal.get("recovery")
    if control is not None and recovery is not None and len(recovery):
        for name, rows in recovery.groupby("variant", sort=False):
            cols = list(rows["variable"])
            written.append(_save(PH.plot_recovery_panels, bal["monthly"][name][cols], control[cols],
                                 bal["fire_date"], bal["tol"], bal["sustain"],
                                 out=fig_dir / "recovery" / f"{name}.png"))
    return [w for w in written if w]


def plot_plan(h: Hypothesis, study: Path):
    import plot_hypotheses as PH
    return _save(PH.plot_variant_matrix, plan_table(h), title=f"{h.key}: {h.description}",
                 out=study / "figures" / f"{h.key}_variants.png")


# ── CLI shared by the hyp*.py files and run_all_hypotheses.py ────────────────
def add_common_args(p: argparse.ArgumentParser,
                    only_help: str = "restrict to these variants") -> argparse.ArgumentParser:
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--out-root", type=Path, default=None,
                   help="parent of the study folders (default: general.out_root or <WB_PRJ_ROOT>/hypotheses)")
    p.add_argument("--dry-run", action="store_true",
                   help="print the runs and draw the variant matrix; run nothing")
    p.add_argument("--only", nargs="*", default=None, metavar="NAME", help=only_help)
    p.add_argument("--force", action="store_true", help="rerun variants that already finished")
    p.add_argument("--no-run", action="store_true", help="skip the runs (validate / plot existing ones)")
    p.add_argument("--no-validate", action="store_true", help="skip the EO-EB validation (and its plots)")
    p.add_argument("--no-balance", action="store_true",
                   help="skip the water-balance / post-fire totals and recovery tables (and their figures)")
    p.add_argument("--no-figs", action="store_true")
    p.add_argument("--post-plots", action="store_true",
                   help="per-run figures (plot_run_results / ET compare) at the end of each run")
    p.add_argument("--stop-on-error", action="store_true")
    p.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                   help="appended to every driver call, e.g. --extra --end-year 2018 (put last)")
    return p


def run_hypothesis(h: Hypothesis, a) -> int:
    """Plan -> (run) -> (validate) -> (plot) one hypothesis; returns 0 / 1."""
    study = study_dir(h, a.out_root)
    only = set(a.only) if a.only else None
    if only and (unknown := only - {v.name for v in h.variants}):
        raise SystemExit(f"{h.key}: unknown variant(s) {sorted(unknown)}; "
                         f"available: {[v.name for v in h.variants]}")
    if a.dry_run:
        print_plan(h, study, a.extra, a.post_plots)
        if not a.no_figs:
            try:
                print(f"\n  variant matrix figure: {plot_plan(h, study)}")
            except Exception as exc:
                print(f"\n  (variant matrix figure skipped: {exc})")
        return 0
    manifest = load_manifest(study) if a.no_run else run_variants(
        h, study, only=only, force=a.force, extra=a.extra, post_plots=a.post_plots,
        stop_on_error=a.stop_on_error)
    failed = [n for n, e in manifest.items() if e.get("status") != "ok"]
    metrics = bal = None
    if not a.no_validate:
        metrics, wb, eo = validate_variants(h, study)
        if not a.no_figs:
            for p in plot_study(h, study, metrics, wb, eo):
                print(f"  figure: {p}")
        if len(metrics):
            print("\n" + metrics.drop(columns=["eo_var"], errors="ignore").round(3).to_string(index=False))
    if not getattr(a, "no_balance", False):       # getattr: namespaces built before the flag existed
        bal = balance_study(h, study)
        if not a.no_figs and bal["monthly"]:
            try:
                for p in plot_balance(h, study, bal):
                    print(f"  figure: {p}")
            except Exception as exc:                          # figures never fail the study
                bal["notes"].append(f"water-balance figures failed ({exc})")
        summ = bal["tables"].get("summary") if bal["tables"] else None
        if summ is not None and len(summ):
            cols = ["variant", "post_eta_mm_yr", "post_recharge_mm_yr", "post_runoff_mm_yr", "post_closure_pct"]
            print("\npost-fire water balance (mm/yr)\n" +
                  summ[[c for c in cols if c in summ]].round(1).to_string(index=False))
        for note in bal["notes"]:
            print(f"  note: {note}")
    print(f"  report: {write_report(h, study, manifest, metrics, bal)}")
    if failed:
        print(f"\n[{h.key}] failed variants: {', '.join(failed)}")
    return 1 if failed else 0


def hypothesis_run(key: str, a) -> int:
    """Run one hypothesis by key with already-parsed arguments (used by run_all_hypotheses.py)."""
    hyps = load_hypotheses(a.config)
    if key not in hyps:
        raise SystemExit(f"'{key}' is not defined in {a.config}")
    return run_hypothesis(hyps[key], a)


def hypothesis_main(key: str, argv=None) -> int:
    """Entry point of a hyp*.py file (all of them are this call with their own key)."""
    mod = sys.modules.get("__main__")
    doc = (getattr(mod, "__doc__", None) or key).strip().split("\n\n")[0]
    p = add_common_args(argparse.ArgumentParser(description=doc,
                                                formatter_class=argparse.ArgumentDefaultsHelpFormatter))
    return hypothesis_run(key, p.parse_args(argv))
