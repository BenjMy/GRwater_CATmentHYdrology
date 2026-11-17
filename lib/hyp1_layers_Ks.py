#!/usr/bin/env python3
"""
Hypothesis 1 - vertical layering and Ks handling: one vs two soil layers, static vs LAI-driven (post-fire) Ks.

Question: does a tighter layer below the top soil (and a Ks that follows the observed LAI) improve ETa
and the water balance against the EO energy-balance ETa? Variants: hypotheses.toml [hyp1_layers_Ks]
(--permx, --permx-bottom, --ks-layer-depth, --ks-fire-*).

Every run of the hypothesis lives in ONE folder, <out_root>/hyp1_layers_Ks/ (prj/, simulation_log.csv,
manifest.json, metrics.csv, figures/). Validation vs the EO energy-balance ETa (NetCDF) is part of
the hypothesis: see analysis_ET_compare.validate_against_eb and hypotheses_common.validate_variants.

    python hyp1_layers_Ks.py --dry-run                 # list the runs + variant matrix figure, run nothing
    python hyp1_layers_Ks.py                           # run, validate vs EO-EB, plot
    python hyp1_layers_Ks.py --only <variant> --force  # rerun one variant
    python hyp1_layers_Ks.py --no-run                  # validate / plot runs that already exist
    python hyp1_layers_Ks.py --extra --end-year 2018   # shorten every run
"""
import sys

from hypotheses_common import hypothesis_main, hypothesis_run

HYPOTHESIS = "hyp1_layers_Ks"


def run(args) -> int:
    """Called by run_all_hypotheses.py with the parsed common options."""
    return hypothesis_run(HYPOTHESIS, args)


if __name__ == "__main__":
    sys.exit(hypothesis_main(HYPOTHESIS))
