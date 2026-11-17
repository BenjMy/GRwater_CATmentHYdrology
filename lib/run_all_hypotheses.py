#!/usr/bin/env python3
"""
run_all_hypotheses.py - run every hypothesis (hyp1..hyp4), validate vs EO-EB ETa and plot
=======================================================================================
Calls hyp1_layers_Ks.py ... hyp4_SW.py (all defined in hypotheses.toml). Each hypothesis keeps its
runs in its own study folder <out_root>/<hyp>/ ; a cross-hypothesis summary is written to
<out_root>/summary_metrics.csv and summary_skill.png.

    python run_all_hypotheses.py --dry-run                 # every run + one variant-matrix figure each
    python run_all_hypotheses.py                           # run + validate + plot everything
    python run_all_hypotheses.py --only hyp1 hyp3          # prefix of the hypothesis key is enough
    python run_all_hypotheses.py --no-run                  # only validate / plot existing runs
    python run_all_hypotheses.py --extra --end-year 2018   # shorten every run (put last)

Exit code 0 if every run of every selected hypothesis succeeded, 1 otherwise.
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

from hypotheses_common import add_common_args, default_out_root, load_hypotheses, study_dir

HERE = Path(__file__).resolve().parent
MODULES = ("hyp1_layers_Ks", "hyp2_LAI_SCF", "hyp3_ic", "hyp4_SW")   # one file per hypothesis


def select(hyps: dict, only) -> list:
    if not only:
        return list(hyps.values())
    chosen = []
    for tok in only:
        hit = [h for k, h in hyps.items() if k == tok or k.startswith(tok)]
        if not hit:
            raise SystemExit(f"Unknown hypothesis '{tok}'. Available: {', '.join(hyps)}")
        chosen += [h for h in hit if h not in chosen]
    return chosen


def summary(hyps: list, out_root: Path) -> None:
    """Concatenate the per-hypothesis metrics.csv and draw one skill figure."""
    import pandas as pd
    frames = []
    for h in hyps:
        p = study_dir(h, out_root) / "metrics.csv"
        if p.exists() and p.stat().st_size > 2:
            frames.append(pd.read_csv(p).assign(hypothesis=h.key))
    if not frames:
        print("\nNo metrics to summarise.")
        return
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(out_root / "summary_metrics.csv", index=False)
    import plot_hypotheses as PH
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
    plt.close(PH.plot_hypothesis_skill(df.assign(variant=df["hypothesis"] + " / " + df["variant"]),
                                       title="All hypotheses: skill vs EO-EB ETa",
                                       out=out_root / "summary_skill.png"))
    print(f"\nSummary: {out_root / 'summary_metrics.csv'}\n         {out_root / 'summary_skill.png'}")


def main(argv=None) -> int:
    p = add_common_args(argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                                formatter_class=argparse.ArgumentDefaultsHelpFormatter),
                        only_help="hypothesis keys or prefixes (hyp1, hyp2_LAI_SCF, ...); default: all")
    a = p.parse_args(argv)
    hyps = select(load_hypotheses(a.config), a.only)
    a.only = None                                 # here --only picked hypotheses, not variants
    out_root = Path(a.out_root or hyps[0].out_root or default_out_root())
    a.out_root = out_root
    rc = 0
    for h in hyps:
        mod = importlib.import_module(h.key)      # hyp1_layers_Ks.py ... : one file per hypothesis
        assert mod.HYPOTHESIS == h.key, (mod.__name__, h.key)
        rc |= mod.run(a)
    if not a.dry_run and not a.no_figs and not a.no_validate:
        summary(hyps, Path(out_root))
    return rc


if __name__ == "__main__":
    sys.exit(main())
