#!/usr/bin/env python3
"""
Hypothesis 3 - initial condition (pressure head or water-table depth), with and without spin-up.

Question: how long does the model remember its initial state, and does spin-up remove the dependence?
Variants: hypotheses.toml [hyp3_ic] (--ic-mode, --ic-pressure-head, --ic-wtposition, --spinup-cycles).

Every run of the hypothesis lives in ONE folder, <out_root>/hyp3_ic/ (prj/, simulation_log.csv,
manifest.json, metrics.csv, figures/). Validation vs the EO energy-balance ETa (NetCDF) is part of
the hypothesis: see analysis_ET_compare.validate_against_eb and hypotheses_common.validate_variants.

    python hyp3_ic.py --dry-run                 # list the runs + variant matrix figure, run nothing
    python hyp3_ic.py                           # run, validate vs EO-EB, plot
    python hyp3_ic.py --only <variant> --force  # rerun one variant
    python hyp3_ic.py --no-run                  # validate / plot runs that already exist
    python hyp3_ic.py --extra --end-year 2018   # shorten every run
"""
import sys

from hypotheses_common import hypothesis_main, hypothesis_run

HYPOTHESIS = "hyp3_ic"


def run(args) -> int:
    """Called by run_all_hypotheses.py with the parsed common options."""
    return hypothesis_run(HYPOTHESIS, args)


if __name__ == "__main__":
    sys.exit(hypothesis_main(HYPOTHESIS))
