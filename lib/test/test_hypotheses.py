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

# ── every hyp file exposes HYPOTHESIS + run ──
for mod in ("hyp1_layers_Ks", "hyp2_LAI_SCF", "hyp3_ic", "hyp4_SW", "run_all_hypotheses"):
    m = __import__(mod)
    if mod != "run_all_hypotheses":
        assert m.HYPOTHESIS in hyps and callable(m.run), mod

print("hypotheses test OK")
