#!/usr/bin/env python3
"""
plot_compare_runs.py - comparison figures for runs picked from a CSV run log
============================================================================
Selects runs with run_catalog.load_runs, labels each one with the parameters that actually differ
between them (``#4 permx=0.0006``), and draws the figures of plot_hypotheses.py. No study folder or
manifest is needed, only the run log and the scenario folders.

By default the selection is the runs that finished (status ok), are not partial, and are the
latest of their parameter set; ``--all`` keeps duplicates, ``--include-partial`` keeps runs that
stopped early.

    python plot_compare_runs.py --study <out_root>/hyp1_layers_Ks
    python plot_compare_runs.py --study <study> --sim 2 4 6 --out figs/ks
    python plot_compare_runs.py --study <study> --where "PMIN == -50" --fire-date 2022-08-01
    python plot_compare_runs.py --log-file my.csv --path2prj ../prj --eo     # + ETa vs EO-EB

Figures written to ``<study>/figures/compare_<tag>/`` (or ``--out``):

    response.png   ETa, soil moisture, recharge and psi (monthly_stats.csv)        [always]
    annual.png     annual ETa / recharge / runoff totals                              [if data]
    balance.png    mean annual ETa, recharge, runoff and closure residual            [if data]
    eta_eo.png, skill.png   ETa vs the EO energy-balance product                      [--eo]
    runs.csv       the runs of the figures with their labels and metrics

Pure matplotlib/pandas apart from ``--eo`` (analysis_ET_compare: geopandas, rioxarray).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import pandas as pd

import hypothesis_tables as HT
import plot_hypotheses as PH
from run_catalog import load_runs, run_labels, varying_params


def select_runs(log, prj, *, sims=None, where=None, keep_all=False, include_partial=False,
                metrics=True, fire_date=None) -> pd.DataFrame:
    """Runs to compare (see the module docstring for the default selection)."""
    runs = load_runs(log, prj, where=where, metrics=metrics, fire_date=fire_date)
    if sims:
        runs = runs[runs["sim_index"].isin([int(s) for s in sims])]
        missing = sorted(set(int(s) for s in sims) - set(runs["sim_index"]))
        if missing:
            print(f"  not in the log (or removed by --where): {missing}")
    else:
        runs = runs[runs["status"].isin(["ok", "unknown"])]
    runs = runs[runs["outputs_ok"]]
    if not include_partial:
        dropped = runs[runs["partial"]]["sim_index"].tolist()
        if dropped:
            print(f"  skipped partial runs (stopped early): {dropped}  (--include-partial to keep)")
        runs = runs[~runs["partial"]]
    if not keep_all and not sims:
        runs = runs.drop_duplicates("param_hash", keep="last")
    return runs.sort_values("sim_index").reset_index(drop=True)


def compare(runs: pd.DataFrame, out: Path, *, fire_date=None, title: str = "", eo: bool = False,
            dpi: int = 150) -> dict:
    """Draw the figures for ``runs``; returns {figure name: path or reason it was skipped}."""
    out.mkdir(parents=True, exist_ok=True)
    fire_date = fire_date or HT.default_fire_date()
    labels = run_labels(runs)
    done = {}

    stats, monthly = {}, {}
    for _, r in runs.iterrows():
        lab = labels[int(r["sim_index"])]
        p = Path(r["sim_path"])
        if (p / "monthly_stats.csv").exists():
            stats[lab] = pd.read_csv(p / "monthly_stats.csv")
        frame, _ = HT.load_monthly_balance(p)
        if not frame.empty:
            monthly[lab] = frame

    def save(name, fig_fn, *args, **kw):
        path = out / f"{name}.png"
        try:
            fig = fig_fn(*args, out=str(path), **kw)
        except Exception as exc:                       # one broken figure must not stop the rest
            done[name] = f"failed: {exc}"
            return
        done[name] = str(path) if fig is not None else "no data"
        if fig is not None:
            matplotlib.pyplot.close(fig)

    if stats:
        save("response", PH.plot_hypothesis_response, stats, title=title)
    if monthly:
        tabs = HT.tables_by_variant(monthly, fire_date)
        save("annual", PH.plot_hypothesis_annual, tabs["annual"], fire_date=fire_date, title=title)
        pname = "post_fire" if (tabs["period"]["period"] == "post_fire").any() else "all"
        save("balance", PH.plot_hypothesis_balance, tabs["period"], title=title, period_name=pname)
    if eo:
        done.update(_compare_eo(runs, labels, out, title))

    table = runs.drop(columns=["sim_path"]).assign(label=[labels[int(i)] for i in runs["sim_index"]])
    table.to_csv(out / "runs.csv", index=False)
    done["runs.csv"] = str(out / "runs.csv")
    return done


def _compare_eo(runs, labels, out: Path, title: str) -> dict:
    """ETa of every run vs the EO energy-balance product (slow: opens the shapefiles and EO file)."""
    try:
        import analysis_ET_compare as AE
        gdf, ds_eo = AE.load_shapefiles(), AE.load_eo_et()
    except Exception as exc:
        return {"eta_eo": f"skipped: {exc}", "skill": f"skipped: {exc}"}
    rows, wb, eo = [], {}, {}
    for _, r in runs.iterrows():
        lab = labels[int(r["sim_index"])]
        try:
            res = AE.validate_against_eb(Path(r["sim_path"]), gdf_proj=gdf, ds_eo=ds_eo,
                                         full_figures=False, fig_dir=out / "per_run" / Path(r["sim_path"]).name)
        except Exception as exc:
            print(f"  {lab}: validation failed - {exc}")
            continue
        wb[lab] = res["s_wb"]
        eo[res["eo_var"]] = res["s_eo"]
        rows.append({"variant": lab, **res["metrics"]})
    done = {}
    if wb:
        fig = PH.plot_hypothesis_eta(wb, eo, title=title, out=str(out / "eta_eo.png"))
        done["eta_eo"] = str(out / "eta_eo.png") if fig is not None else "no data"
        fig = PH.plot_hypothesis_skill(pd.DataFrame(rows), title=title, out=str(out / "skill.png"))
        done["skill"] = str(out / "skill.png") if fig is not None else "no data"
        pd.DataFrame(rows).to_csv(out / "metrics.csv", index=False)
    else:
        done["eta_eo"] = done["skill"] = "no run could be validated"
    return done


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--log-file", default="simulation_log_withETp_withLAI.csv", help="CSV run log")
    ap.add_argument("--path2prj", default="../SSHydro_withETp_withLAI", help="project folder holding outputs/")
    ap.add_argument("--study", help="hypothesis study: folder, or just its name from the config (hyp4, hyp4_SW); uses its simulation_log.csv and prj/")
    ap.add_argument("--sim", nargs="+", type=int, help="sim_index values to compare (default: automatic selection)")
    ap.add_argument("--where", help='pandas query on the log columns, e.g. "PMIN == -50"')
    ap.add_argument("--all", action="store_true", help="keep runs that repeat an earlier parameter set")
    ap.add_argument("--include-partial", action="store_true", help="keep runs that stopped early")
    ap.add_argument("--fire-date", help="fire date for the pre/post split (default: hypothesis_tables.default_fire_date)")
    ap.add_argument("--eo", action="store_true", help="also compare ETa with the EO energy-balance product (slow)")
    ap.add_argument("--out", help="output folder (default: <study or log folder>/figures/compare_<tag>)")
    ap.add_argument("--title", default="", help="figure title")
    ap.add_argument("--dpi", type=int, default=150)
    a = ap.parse_args(argv)
    if a.study:
        from hypotheses_common import resolve_study
        a.study = resolve_study(a.study)
    log, prj = (Path(a.study) / "simulation_log.csv", Path(a.study) / "prj") if a.study \
        else (Path(a.log_file), Path(a.path2prj))
    if not log.exists():
        print(f"log not found: {log}")
        return 1
    runs = select_runs(log, prj, sims=a.sim, where=a.where, keep_all=a.all,
                       include_partial=a.include_partial, fire_date=a.fire_date)
    if len(runs) < 1:
        print("no run to compare")
        return 1
    tag = "sim-" + "-".join(map(str, runs["sim_index"]))
    out = Path(a.out) if a.out else log.parent / "figures" / f"compare_{tag}"
    print(f"comparing {len(runs)} run(s); varying: {', '.join(varying_params(runs)) or '(none)'}")
    for lab in run_labels(runs).values():
        print(f"  {lab}")
    done = compare(runs, out, fire_date=a.fire_date, title=a.title, eo=a.eo, dpi=a.dpi)
    print(f"\nfigures in {out}")
    for k, v in done.items():
        print(f"  {k:10s} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
