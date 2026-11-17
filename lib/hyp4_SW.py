#!/usr/bin/env python3
"""
Hypothesis 4 - ETp forcing: wet vs dry Shuttleworth-Wallace bounds, and the gridded satellite ETp.

Question: which ETp bound reproduces the EO energy-balance ETa best? Each run is compared with the EO-EB
variable of the same condition when the EO file provides one. Variants: hypotheses.toml [hyp4_SW]
(--forcing, --etp-ts-column).

Every run of the hypothesis lives in ONE folder, <out_root>/hyp4_SW/ (prj/, simulation_log.csv,
manifest.json, metrics.csv, figures/). Validation vs the EO energy-balance ETa (NetCDF) is part of
the hypothesis: see analysis_ET_compare.validate_against_eb and hypotheses_common.validate_variants.

    python hyp4_SW.py --dry-run                 # list the runs + variant matrix figure, run nothing
    python hyp4_SW.py                           # run, validate vs EO-EB, plot
    python hyp4_SW.py --only <variant> --force  # rerun one variant
    python hyp4_SW.py --no-run                  # validate / plot runs that already exist
    python hyp4_SW.py --extra --end-year 2018   # shorten every run
"""
import sys

from hypotheses_common import hypothesis_main, hypothesis_run

HYPOTHESIS = "hyp4_SW"


def run(args) -> int:
    """Called by run_all_hypotheses.py with the parsed common options."""
    return hypothesis_run(HYPOTHESIS, args)


if __name__ == "__main__":
    sys.exit(hypothesis_main(HYPOTHESIS))
