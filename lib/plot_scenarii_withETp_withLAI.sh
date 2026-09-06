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
# =============================================================================
set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SCRIPT_RES="withLAI_results.py"
EXPERIMENT="withETp_withLAI"
PRJ_PATH="../SSHydro_${EXPERIMENT}"
LOG_CSV="simulation_log_${EXPERIMENT}.csv"
SCENARIOS=(1 2 3 4 5 6 7 8 9 10 11 12 13)

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
    exit 1
fi
echo "============================================================"
