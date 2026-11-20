#!/usr/bin/env python3
"""
test_hypotheses.py - assertion-based test of the hypothesis framework (no simulation, no data)
==============================================================================================
Checks variant counts, command building, dry-run, run/skip/manifest with a fake driver,
validation with a fake validator and the plotters with synthetic series.

    python lib/test/test_hypotheses.py          # exits 0 and prints "hypotheses test OK"
"""
import argparse
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import numpy as np
import pandas as pd

LIB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LIB))
sys.path.insert(0, str(LIB / "utils"))

import matplotlib
matplotlib.use("Agg")

import hypotheses_common as H

EXPECTED = {"hyp1_layers_Ks": 9, "hyp2_LAI_SCF": 6, "hyp3_ic": 9, "hyp4_SW": 3}
hyps = H.load_hypotheses()
assert {k: len(h.variants) for k, h in hyps.items()} == EXPECTED, hyps.keys()

# unique names, shared options come from [common]
for h in hyps.values():
    names = [v.name for v in h.variants]
    assert len(set(names)) == len(names), h.key
    assert all(v.args.get("forcing") for v in h.variants), h.key
    assert all(v.varying for v in h.variants), h.key

# ── command building ──
tmp = Path(tempfile.mkdtemp())
h2 = hyps["hyp2_LAI_SCF"]
cmd = H.variant_command(h2, next(v for v in h2.variants if v.args["pmin"] == "-1e35"), tmp)
assert "--pmin=-1e35" in cmd and "--no-plots" in cmd
assert cmd[cmd.index("--path2prj") + 1] == str(tmp / "prj")
assert "--no-spinup" in H.variant_command(hyps["hyp3_ic"], hyps["hyp3_ic"].variants[0], tmp)
sp = next(v for v in hyps["hyp3_ic"].variants if v.name == "wt1_spinup3")
assert "--no-spinup" not in H.variant_command(hyps["hyp3_ic"], sp, tmp)
two = next(v for v in hyps["hyp1_layers_Ks"].variants if v.name == "2L_d0.5_b1e-5")
c2 = H.variant_command(hyps["hyp1_layers_Ks"], two, tmp)
assert c2[c2.index("--permx-bottom") + 1] == "1e-5" and c2[c2.index("--ks-layer-depth") + 1] == "0.5"
assert len(H.plan_table(h2)) == 6

# ── dry-run (no driver call) ──
calls = []
real_run = subprocess.run
H.subprocess.run = lambda *a, **k: calls.append(a) or types.SimpleNamespace(returncode=0)
ns = dict(config=H.DEFAULT_CONFIG, out_root=tmp / "dry", only=None, force=False, extra=[], post_plots=False,
          stop_on_error=False, dry_run=True, no_figs=False, no_run=False, no_validate=False)
assert H.run_hypothesis(hyps["hyp4_SW"], argparse.Namespace(**ns)) == 0 and not calls
assert (tmp / "dry" / "hyp4_SW" / "figures" / "hyp4_SW_variants.png").exists()

# ── run with a fake driver: manifest, skip-if-done, failure ──
def fake_run(cmd, cwd=None, **kw):
    prj = Path(cmd[cmd.index("--path2prj") + 1])
    name = "scenario_" + "_".join(c for c in cmd[2:] if c.startswith("--etp-ts-column") or not c.startswith("--"))[:60]
    d = prj / "outputs" / name
    d.mkdir(parents=True, exist_ok=True)
    if "ET_0_dry" in cmd:
        return types.SimpleNamespace(returncode=1)
    (d / "et_output.nc").write_text("x")
    return types.SimpleNamespace(returncode=0)

H.subprocess.run = fake_run
study = tmp / "run" / "hyp4_SW"
man = H.run_variants(hyps["hyp4_SW"], study)
assert man["SW_wet"]["status"] == "ok" and man["SW_dry"]["status"] == "failed"
assert man["satellite_etp"]["status"] == "ok"
assert (study / "manifest.json").exists()
n_before = len(H.load_manifest(study))
H.run_variants(hyps["hyp4_SW"], study, only={"SW_wet"})          # already done -> skipped
assert len(H.load_manifest(study)) == n_before
H.subprocess.run = real_run

# ── validation with a fake validator + plots with synthetic series ──
idx = pd.date_range("2017-01-01", periods=24, freq="MS")
rng = np.random.default_rng(0)
eo_s = pd.Series(30 + 20 * np.sin(np.arange(24) / 12 * 2 * np.pi), idx)

def validator(sim_dir):
    s = eo_s * (1.0 if "wet" in Path(sim_dir).name else 0.8) + rng.normal(0, 2, 24)
    return {"s_wb": s, "s_eo": eo_s, "eo_var": "ET", "condition": "wet",
            "metrics": {"bias": float((s - eo_s).mean()), "rmse": float(((s - eo_s) ** 2).mean() ** .5),
                        "r": float(np.corrcoef(s, eo_s)[0, 1]), "nse": 0.5, "kge": 0.6}}

metrics, wb, eo = H.validate_variants(hyps["hyp4_SW"], study, validator=validator)
assert len(metrics) == 2 and set(wb) == {"SW_wet", "satellite_etp"} and (study / "metrics.csv").exists()
paths = H.plot_study(hyps["hyp4_SW"], study, metrics, wb, eo)
assert len(paths) >= 2 and all(Path(p).exists() for p in paths), paths

# ── water-balance tables, post-fire totals, recovery hook, cross-hypothesis summary ──
import dataclasses
import json

import hypothesis_tables as HT
import run_all_hypotheses as RA


def fake_sim(root, name, eta_pre, eta_post, rain=50.0, rech=10.0, q=0.01):
    """monthly_stats.csv + hgraph_output.csv of a 2019-2022 run (fire 2020-07: 18 pre / 30 post months)."""
    d = root / name
    d.mkdir(parents=True)
    months = pd.date_range("2019-01-01", "2022-12-01", freq="MS")
    eta = np.where(months < pd.Timestamp("2020-07-01"), eta_pre, eta_post)
    pd.DataFrame({"month": months.strftime("%Y-%m"), "rain_in_mm": rain, "eta_out_mm": eta,
                  "recharge_out_mm": rech, "sw_mean": 0.4}).to_csv(d / "monthly_stats.csv", index=False)
    pd.DataFrame({"datetime": pd.date_range("2019-01-01", "2022-12-31", freq="D"),
                  "streamflow": q}).to_csv(d / "hgraph_output.csv", index=False)
    return d


bal_root = tmp / "bal"
h4 = hyps["hyp4_SW"]
sims = {"SW_wet": fake_sim(bal_root / "sims", "wet", 30.0, 20.0),
        "satellite_etp": fake_sim(bal_root / "sims", "sat", 30.0, 30.0)}          # SW_dry left unrun
bstudy = bal_root / "hyp4_SW"
H.save_manifest(bstudy, {n: {"status": "ok", "sim_dir": str(p), "group": "", "returncode": 0, "seconds": 1}
                         for n, p in sims.items()})

# monthly frame: units and runoff volume (0.01 m3/s over January = 26784 m3)
m, notes = HT.load_monthly_balance(sims["SW_wet"], 1e6)
assert notes == [] and np.isclose(m.loc["2019-01-01", "runoff_m3"], 0.01 * 31 * 86400)
assert np.isclose(m.loc["2019-01-01", "runoff_mm"], 26.784)
assert np.isclose(m.loc["2019-01-01", "eta_mm_day"], 30.0 / 31)
m0, n0 = HT.load_monthly_balance(sims["SW_wet"])                            # no area: runoff in m3 only
assert m0["runoff_mm"].isna().all() and any("catchment_area_m2" in n for n in n0)
assert HT.load_monthly_balance(bal_root / "nope")[0].empty                  # missing run -> empty, no raise

# annual totals / closure, and the pre / post-fire split
an = HT.annual_totals(m)
y19 = an[an["year"] == 2019].iloc[0]
assert y19["n_months"] == 12 and np.isclose(y19["eta_mm"], 360) and np.isclose(y19["closure_pct"], 20.0)
pt = HT.period_totals(m, "2020-07-01").set_index("period")
assert list(pt.index) == ["all", "pre_fire", "post_fire"]
assert pt.loc["pre_fire", "n_months"] == 18 and pt.loc["post_fire", "n_months"] == 30
assert np.isclose(pt.loc["post_fire", "total_eta_mm"], 30 * 20) and np.isclose(pt.loc["post_fire", "annual_eta_mm"], 240)
assert np.isclose(pt.loc["pre_fire", "annual_eta_mm"], 360) and np.isclose(pt.loc["post_fire", "closure_pct"], 40.0)

# balance_study: files, missing variant reported, no control -> recovery skipped with a note
b = H.balance_study(dataclasses.replace(h4, validation={**h4.validation, "catchment_area_m2": 1e6}), bstudy)
for f in ("water_balance_monthly.csv", "annual_totals.csv", "period_totals.csv", "variant_summary.csv"):
    assert (bstudy / f).exists(), f
assert not (bstudy / "recovery_metrics.csv").exists() and b["recovery"].empty
assert any("SW_dry" in n and "no finished run" in n for n in b["notes"]), b["notes"]
assert any("control" in n for n in b["notes"]), b["notes"]
sm = pd.read_csv(bstudy / "variant_summary.csv").set_index("variant")
assert set(sm.index) == {"SW_wet", "satellite_etp"}
assert np.isclose(sm.loc["SW_wet", "post_eta_mm_yr"], 240) and np.isclose(sm.loc["satellite_etp", "post_eta_mm_yr"], 360)
assert sm.loc["SW_wet", "post_eta_mm_yr"] < sm.loc["SW_wet", "mean_eta_mm_day"] * 365.25      # post-fire drop visible
assert np.isclose(sm.loc["SW_wet", "mean_sw"], 0.4) and np.isclose(sm.loc["SW_wet", "post_runoff_mm_yr"], 315.36, rtol=0.01)
figs = H.plot_balance(h4, bstudy, b)
assert len(figs) == 2 and all(Path(p).exists() for p in figs), figs

# with a control run: recovery table + one figure per variant (control burns less, so ratio < 1 after the fire)
ctrl = fake_sim(bal_root / "sims", "control", 30.0, 30.0)
hc = dataclasses.replace(h4, validation={**h4.validation, "control_sim_dir": str(ctrl),
                                         "recovery_tol": 0.05, "recovery_sustain": 2})
bc = H.balance_study(hc, bstudy)
rec = pd.read_csv(bstudy / "recovery_metrics.csv")
assert (bstudy / "recovery_metrics.csv").exists() and set(rec["variant"]) == {"SW_wet", "satellite_etp"}
wet = rec[(rec["variant"] == "SW_wet") & (rec["variable"] == "eta_mm_day")].iloc[0]
sat = rec[(rec["variant"] == "satellite_etp") & (rec["variable"] == "eta_mm_day")].iloc[0]
assert np.isnan(wet["recovery_months"]) and wet["deficit"] > 0                # never recovers, loses water
assert sat["recovery_months"] == 0 and np.isclose(sat["deficit"], 0.0)       # identical to the control
assert len(H.plot_balance(hc, bstudy, bc)) == 4                               # 2 + recovery per variant
b_off = H.balance_study(h4, bstudy)                                           # control switched off again
assert not (bstudy / "recovery_metrics.csv").exists()                         # stale file removed

# a broken run is reported, not raised
def bad_loader(sim_dir, area):
    raise OSError("disk gone")
bb = H.balance_study(h4, bstudy, loader=bad_loader)
assert bb["monthly"] == {} and any("unreadable" in n for n in bb["notes"])

# study report lists failures and the unscored variant
metrics_fake = pd.DataFrame({"variant": ["SW_wet"], "kge": [0.5]})
rp = H.write_report(h4, bstudy, H.load_manifest(bstudy), metrics_fake, b).read_text()
assert "1 of 2 finished variants scored" in rp and "satellite_etp: finished but not scored" in rp
assert "SW_dry: not run" in rp and "## Issues" in rp

# summary_metrics.csv: skill and means joined, outer join keeps a variant without skill
metrics_fake = pd.DataFrame({"variant": ["SW_wet"], "group": ["SW"], "eo_var": ["ET"], "kge": [0.5],
                             "rmse": [1.0], "bias": [0.1], "r": [0.9], "nse": [0.4]})
metrics_fake.to_csv(bstudy / "metrics.csv", index=False)
H.balance_study(h4, bstudy)
RA.summary([h4], bal_root)
sumry = pd.read_csv(bal_root / "summary_metrics.csv").set_index("variant")
assert {"kge", "mean_eta_mm_day", "mean_recharge_mm_day", "mean_sw", "post_eta_mm_yr"} <= set(sumry.columns)
assert np.isclose(sumry.loc["SW_wet", "kge"], 0.5) and np.isnan(sumry.loc["satellite_etp", "kge"])
assert np.isclose(sumry.loc["satellite_etp", "post_recharge_mm_yr"], 120)
assert (bal_root / "summary_period_totals.csv").exists() and (bal_root / "summary_skill.png").exists()

# --no-balance flag exists, and a namespace without it still works (getattr)
assert H.add_common_args(argparse.ArgumentParser()).parse_args(["--no-balance"]).no_balance

# end to end through run_hypothesis with the fake driver replaced by the finished study (no-run)
ns2 = dict(ns, out_root=bal_root, dry_run=False, no_run=True, no_validate=True, no_figs=False)
assert H.run_hypothesis(h4, argparse.Namespace(**ns2)) == 0
assert (bstudy / "study_report.md").exists() and (bstudy / "figures" / "hyp4_SW_water_balance.png").exists()
json.loads((bstudy / "manifest.json").read_text())

# ── every hyp file exposes HYPOTHESIS + run ──
for mod in ("hyp1_layers_Ks", "hyp2_LAI_SCF", "hyp3_ic", "hyp4_SW", "run_all_hypotheses"):
    m = __import__(mod)
    if mod != "run_all_hypotheses":
        assert m.HYPOTHESIS in hyps and callable(m.run), mod

print("hypotheses test OK")
