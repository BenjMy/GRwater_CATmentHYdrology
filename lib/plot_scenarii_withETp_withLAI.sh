#!/usr/bin/env bash
# =============================================================================
# GRwater — Agramon LAI hydrology batch post-processor
# Post-processes withETp_withLAI experiment scenarios by sim_index.
#
# NOTE: SCENARIOS=(1..9) below assumes one sim_index per dem-plot, which
# only holds if run_scenarii_withETp_withLAI.sh was run with a single
# outlet/bottom/spinup combo. Once you run multiple BC/spin-up combos per
# plot, sim_index no longer maps 1:1 to dem-plot — check
# simulation_log_${EXPERIMENT}.csv (or list ../SSHydro_${EXPERIMENT}/
# outputs/, folders are now self-describing, e.g.
# scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c) for the actual
# sim_index list and update SCENARIOS accordingly.
# withLAI_results.py resolves the tagged folder from sim_index for you —
# only this list needs to be kept in sync.
# analysis_ET_compare.py now mirrors that: its figures are written to
# ../Analysis/ET_compare/<scenario_dirname>/ (the tagged scenario folder
# name), so each scenario's ET-compare plots get their own subfolder
# instead of overwriting one another in a flat directory.
# =============================================================================
set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SCRIPT_RES="withLAI_results.py"
SCRIPT_CMP="analysis_ET_compare.py"
EXPERIMENT="withETp_withLAI"
PRJ_PATH="../SSHydro_${EXPERIMENT}"
LOG_CSV="simulation_log_${EXPERIMENT}.csv"
SCENARIOS=(1 2 3 4 5 6 7 8 9 10 11 12 13)   # default list; override with --scenario

usage() {
    cat <<EOF
Usage: $(basename "$0") [--scenario N [N ...]]

Post-processes withETp_withLAI experiment scenarios by sim_index: runs
${SCRIPT_RES} then ${SCRIPT_CMP} for each scenario in turn.

  --scenario N [N ...]  sim_index value(s) to post-process (default:
                         ${SCENARIOS[*]}). Overrides the built-in default
                         list entirely, e.g. --scenario 7 runs only
                         scenario 7. sim_index no longer maps 1:1 to
                         dem-plot once multiple outlet/bottom/spinup combos
                         are run per plot — check ${LOG_CSV} or list
                         ${PRJ_PATH}/outputs/ (folders are self-describing,
                         e.g. scenario_7_plot2_outlet-free_bottom-closed_
                         spinup-on3c) for the actual sim_index list.
  --help, -h             Show this help and exit

Examples:
  $(basename "$0")
      # post-processes the default scenario list: ${SCENARIOS[*]}

  $(basename "$0") --scenario 7
      # post-processes only scenario 7

  $(basename "$0") --scenario 12 13 14 15
      # post-processes scenarios 12, 13, 14, 15
EOF
}

# -----------------------------------------------------------------------------
# Parse args
# -----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            usage; exit 0 ;;
        --scenario)
            shift
            SCENARIOS=()
            while [[ $# -gt 0 && "$1" != --* ]]; do
                SCENARIOS+=("$1")
                shift
            done
            if [ ${#SCENARIOS[@]} -eq 0 ]; then
                echo "ERROR: --scenario needs at least one value" >&2
                exit 1
            fi
            ;;
        *)
            echo "ERROR: unknown option '$1'" >&2
            usage
            exit 1
            ;;
    esac
done

for scenario in "${SCENARIOS[@]}"; do
    if ! [[ "$scenario" =~ ^[0-9]+$ ]]; then
        echo "ERROR: --scenario values must be non-negative integers (got '${scenario}')" >&2
        exit 1
    fi
done

echo "============================================================"
echo "  GRwater — Agramon batch post-processor"
echo "  Experiment : ${EXPERIMENT}"
echo "  Scenarios  : ${SCENARIOS[*]}"
echo "  Started    : $(date)"
echo "============================================================"
echo

# -----------------------------------------------------------------------------
# Post-processing loop
# -----------------------------------------------------------------------------
failed=()
failed_cmp=()
skipped_cmp=()

for scenario in "${SCENARIOS[@]}"; do
    echo "------------------------------------------------------------"
    echo "  Scenario ${scenario} / ${SCENARIOS[-1]}  —  $(date)"
    echo "------------------------------------------------------------"
    t_start=$(date +%s)

    if python "${SCRIPT_RES}" \
            --scenario "${scenario}" \
            --log-file  "${LOG_CSV}" \
            --path2prj  "${PRJ_PATH}"; then
        t_end=$(date +%s)
        echo "  ✓ Scenario ${scenario} completed in $(( t_end - t_start ))s"
    else
        t_end=$(date +%s)
        echo "  ✗ Scenario ${scenario} FAILED after $(( t_end - t_start ))s" >&2
        failed+=("${scenario}")
        skipped_cmp+=("${scenario}")
        echo "  ⊘ ET compare skipped for scenario ${scenario} (main step failed)"
        echo
        continue
    fi
    echo

    echo "  -- ET compare (EO-EB vs WB) — scenario ${scenario} --------------"
    t_start=$(date +%s)

    if python "${SCRIPT_CMP}" \
            --scenario "${scenario}"; then
        t_end=$(date +%s)
        echo "  ✓ ET compare ${scenario} completed in $(( t_end - t_start ))s"
    else
        t_end=$(date +%s)
        echo "  ✗ ET compare ${scenario} FAILED after $(( t_end - t_start ))s" >&2
        failed_cmp+=("${scenario}")
    fi
    echo
done

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "============================================================"
echo "  Finished : $(date)"

if [[ ${#failed[@]} -eq 0 ]]; then
    echo "  Status   : all ${#SCENARIOS[@]} scenarios completed successfully"
else
    echo "  Status   : ${#failed[@]} scenario(s) FAILED — ${failed[*]}"
fi

if [[ ${#failed_cmp[@]} -eq 0 ]]; then
    echo "  ET compare : all $(( ${#SCENARIOS[@]} - ${#skipped_cmp[@]} )) attempted scenario(s) completed successfully"
else
    echo "  ET compare : ${#failed_cmp[@]} scenario(s) FAILED — ${failed_cmp[*]}"
fi
if [[ ${#skipped_cmp[@]} -gt 0 ]]; then
    echo "  ET compare : ${#skipped_cmp[@]} scenario(s) SKIPPED (main step failed) — ${skipped_cmp[*]}"
fi

if [[ ${#failed[@]} -gt 0 || ${#failed_cmp[@]} -gt 0 ]]; then
    exit 1
fi
echo "============================================================"
