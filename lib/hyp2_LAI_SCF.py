#!/usr/bin/env python3
"""
Hypothesis 2 - LAI-driven root-water-uptake stress (SCF 0/1) crossed with PMIN (-15/-150/-1e35 m).

Question: how much do the stress formulation (Feddes vs PMIN threshold) and the surface pressure-head
floor control ETa? Grid: hypotheses.toml [hyp2_LAI_SCF.grid]. Only the scf_lai CATHY build is available
in the pipeline (agramon_backend.BACKENDS); a classic CATHY without LAI/SCF is not runnable yet.

Every run of the hypothesis lives in ONE folder, <out_root>/hyp2_LAI_SCF/ (prj/, simulation_log.csv,
manifest.json, metrics.csv, figures/). Validation vs the EO energy-balance ETa (NetCDF) is part of
the hypothesis: see analysis_ET_compare.validate_against_eb and hypotheses_common.validate_variants.

    python hyp2_LAI_SCF.py --dry-run                 # list the runs + variant matrix figure, run nothing
    python hyp2_LAI_SCF.py                           # run, validate vs EO-EB, plot
    python hyp2_LAI_SCF.py --only <variant> --force  # rerun one variant
    python hyp2_LAI_SCF.py --no-run                  # validate / plot runs that already exist
    python hyp2_LAI_SCF.py --extra --end-year 2018   # shorten every run
"""
import sys

from hypotheses_common import hypothesis_main, hypothesis_run

HYPOTHESIS = "hyp2_LAI_SCF"


def run(args) -> int:
    """Called by run_all_hypotheses.py with the parsed common options."""
    return hypothesis_run(HYPOTHESIS, args)


if __name__ == "__main__":
    sys.exit(hypothesis_main(HYPOTHESIS))
