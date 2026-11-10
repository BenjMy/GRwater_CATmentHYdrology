#!/usr/bin/env bash
# =============================================================================
# GRwater — Agramon LAI hydrology batch post-processor
# Post-processes withETp_withLAI experiment scenarios by sim_index.
#
# By default the scenario list is now read from simulation_log_${EXPERIMENT}.csv
# (every logged sim_index, optionally narrowed with --forcing), so it no
# longer has to be kept in sync by hand: sim_index does not map 1:1 to
# dem-plot once several BC/spin-up/forcing combos run per plot. If the CSV is
# missing, it falls back to the old fixed list (1..13). Pass --scenario N [N ...]
# to post-process specific sim_index values instead.
# The --forcing filter uses the forcing_mode column the simulation script now
# logs (etp | etp_ts | era5); rows logged before that column existed were all
# satellite-ETp runs, so they count as "etp". etp_ts scenario folders carry an
# extra tag, e.g. scenario_12_plot2_..._forcing-etpts-dry-fid1, but still start
# with scenario_<sim_index>, so withLAI_results.py (which resolves the tagged
# folder from sim_index) needs no change.
# analysis_ET_compare.py now mirrors that: its figures are written to
# ../Analysis/ET_compare/<scenario_dirname>/ (the tagged scenario folder
# name), so each scenario's ET-compare plots get their own subfolder
# instead of overwriting one another in a flat directory. It also now
# accepts --path2prj/--log-file (this script passes PRJ_PATH/LOG_CSV, same
# as for withLAI_results.py) instead of resolving --scenario against its
# own hardcoded WB_PRJ_ROOT/LOG_FILE defaults.
#
# Added with the DEM-resolution / ETa-vs-CSV update:
#   * --dem-res (default 5, = resample_resolution of Agramon_withLAI_withETp.py)
#     is forwarded to withLAI_results.py. Results must rebuild the DEM at the
#     SAME resolution as the run, otherwise the DEM mask no longer matches
#     veg_map_history (e.g. (79,22) vs (32,9) broadcast error).
#   * analysis_ET_compare.py now also draws "WB ETa vs Agramon CSV" (real ET,
#     ETp dry/wet) automatically for etp_ts scenarios (the dry|wet column and
#     feature_id come from the scenario folder name; other forcings skip it).
#     --et-csv PATH / --no-csv are forwarded to it.
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
FALLBACK_SCENARIOS=(1 2 3 4 5 6 7 8 9 10 11 12 13)   # used only if the log CSV is missing
SCENARIOS=()          # filled from LOG_CSV below unless --scenario is given
FORCING_FILTER=""     # empty = all forcings; else etp | etp_ts | era5
DEM_RES="5"           # m; must match resample_resolution of the simulation run
ET_CSV=""             # empty = analysis_ET_compare.py default (agramon_treatments_ET.csv)
NO_CSV=0              # 1 = skip the WB-vs-CSV figure in analysis_ET_compare.py
MONTHLY_TS=0          # 1 = also draw the optional monthly time-series figures

# Print the sim_index values logged in LOG_CSV (one per line), optionally
# restricted to one forcing_mode. Legacy rows without forcing_mode = "etp".
list_logged_scenarios() {
    LST_LOG="${LOG_CSV}" LST_FORCING="${FORCING_FILTER}" python3 - <<'PYEOF'
import os
import pandas as pd

df = pd.read_csv(os.environ["LST_LOG"])
want = os.environ.get("LST_FORCING", "")
if want:
    if "forcing_mode" in df.columns:
        fm = df["forcing_mode"].fillna("etp").astype(str)
    else:
        fm = pd.Series("etp", index=df.index)
    df = df[fm == want]
for i in sorted(df["sim_index"].astype(int).unique()):
    print(i)
PYEOF
}

usage() {
    cat <<EOF
Usage: $(basename "$0") [--scenario N [N ...]] [--forcing {etp|etp_ts|era5}]
                       [--dem-res M] [--et-csv PATH] [--no-csv] [--monthly-ts]

Post-processes withETp_withLAI experiment scenarios by sim_index: runs
${SCRIPT_RES} then ${SCRIPT_CMP} for each scenario in turn.

  --scenario N [N ...]  sim_index value(s) to post-process. Default: every
                         sim_index logged in ${LOG_CSV} (fallback if the
                         file is missing: ${FALLBACK_SCENARIOS[*]}).
                         Overrides that list entirely, e.g. --scenario 7
                         runs only scenario 7. sim_index no longer maps 1:1
                         to dem-plot once multiple outlet/bottom/spinup/
                         forcing combos are run per plot — check
                         ${LOG_CSV} or list ${PRJ_PATH}/outputs/ (folders
                         are self-describing, e.g. scenario_7_plot2_outlet-
                         free_bottom-closed_spinup-on3c) for the actual
                         sim_index list.
  --forcing MODE         Only post-process scenarios logged with this
                         forcing: etp (satellite ETp raster), etp_ts (ETp
                         CSV time series) or era5. Applies to the default
                         list only; ignored if --scenario is given.
  --dem-res M            DEM resample resolution in metres passed to
                         ${SCRIPT_RES} (default: ${DEM_RES}). Must equal the
                         resolution used by the simulation script.
  --et-csv PATH          Agramon treatments ET CSV for the WB-vs-CSV figure
                         of ${SCRIPT_CMP} (default: its built-in path).
  --no-csv               Skip the WB-vs-CSV figure (it is only produced for
                         etp_ts scenarios anyway).
  --monthly-ts           Also draw the optional monthly-resampled time-series
                         figures in ${SCRIPT_RES} (off by default).
  --help, -h             Show this help and exit

Examples:
  $(basename "$0")
      # post-processes every scenario logged in ${LOG_CSV}

  $(basename "$0") --forcing etp_ts
      # only the scenarios run with the ETp time-series forcing

  $(basename "$0") --forcing etp_ts --et-csv /path/agramon_treatments_ET.csv
      # etp_ts scenarios, ETa compared with a specific CSV

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
        --forcing)
            FORCING_FILTER="${2:-}"
            case "${FORCING_FILTER}" in
                etp|etp_ts|era5) ;;
                *) echo "ERROR: --forcing must be etp|etp_ts|era5 (got '${FORCING_FILTER}')" >&2; exit 1 ;;
            esac
            shift 2 ;;
        --dem-res)
            DEM_RES="${2:-}"
            if ! [[ "${DEM_RES}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                echo "ERROR: --dem-res must be a positive number (got '${DEM_RES}')" >&2; exit 1
            fi
            shift 2 ;;
        --et-csv)
            ET_CSV="${2:-}"
            if [ -z "${ET_CSV}" ] || [ ! -f "${ET_CSV}" ]; then
                echo "ERROR: --et-csv file not found: '${ET_CSV}'" >&2; exit 1
            fi
            shift 2 ;;
        --no-csv)
            NO_CSV=1
            shift ;;
        --monthly-ts)
            MONTHLY_TS=1
            shift ;;
        --scenario)
            shift
            SCENARIOS=()
            SCENARIOS_GIVEN=1
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

# Default scenario list: everything logged in LOG_CSV (optionally filtered by
# --forcing); only used when --scenario was not given.
if [ "${SCENARIOS_GIVEN:-0}" != "1" ]; then
    if [ -f "${LOG_CSV}" ]; then
        mapfile -t SCENARIOS < <(list_logged_scenarios)
        if [ ${#SCENARIOS[@]} -eq 0 ]; then
            echo "ERROR: no scenarios found in ${LOG_CSV}${FORCING_FILTER:+ for forcing '${FORCING_FILTER}'}" >&2
            exit 1
        fi
    else
        echo "WARNING: ${LOG_CSV} not found — falling back to scenarios ${FALLBACK_SCENARIOS[*]}" >&2
        SCENARIOS=("${FALLBACK_SCENARIOS[@]}")
    fi
fi

for scenario in "${SCENARIOS[@]}"; do
    if ! [[ "$scenario" =~ ^[0-9]+$ ]]; then
        echo "ERROR: --scenario values must be non-negative integers (got '${scenario}')" >&2
        exit 1
    fi
done

echo "============================================================"
echo "  GRwater — Agramon batch post-processor"
echo "  Experiment : ${EXPERIMENT}"
echo "  DEM res    : ${DEM_RES} m"
if [ "${NO_CSV}" = "1" ]; then
    echo "  ETa vs CSV : disabled (--no-csv)"
else
    echo "  ETa vs CSV : ${ET_CSV:-<default agramon_treatments_ET.csv>} (etp_ts scenarios only)"
fi
if [ "${SCENARIOS_GIVEN:-0}" != "1" ]; then
    echo "  Scenarios  : ${SCENARIOS[*]}${FORCING_FILTER:+  (forcing: ${FORCING_FILTER})}"
else
    echo "  Scenarios  : ${SCENARIOS[*]}"
fi
echo "  Started    : $(date)"
echo "============================================================"
echo

# -----------------------------------------------------------------------------
# Post-processing loop
# -----------------------------------------------------------------------------
failed=()
failed_cmp=()
skipped_cmp=()

# Optional extra args for analysis_ET_compare.py (kept in an array so empty
# values never become stray "" arguments).
RES_EXTRA=()
[ "${MONTHLY_TS}" = "1" ] && RES_EXTRA+=(--monthly-ts)
CMP_EXTRA=()
[ -n "${ET_CSV}" ] && CMP_EXTRA+=(--et-csv "${ET_CSV}")
[ "${NO_CSV}" = "1" ] && CMP_EXTRA+=(--no-csv)

for scenario in "${SCENARIOS[@]}"; do
    echo "------------------------------------------------------------"
    echo "  Scenario ${scenario} / ${SCENARIOS[-1]}  —  $(date)"
    echo "------------------------------------------------------------"
    t_start=$(date +%s)

    if python "${SCRIPT_RES}" \
            --scenario "${scenario}" \
            --log-file  "${LOG_CSV}" \
            --path2prj  "${PRJ_PATH}" \
            --dem-res   "${DEM_RES}" \
            "${RES_EXTRA[@]+"${RES_EXTRA[@]}"}"; then
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

    echo "  -- ET compare (EO-EB vs WB, + WB vs CSV for etp_ts) — scenario ${scenario} --"
    t_start=$(date +%s)

    if python "${SCRIPT_CMP}" \
            --scenario "${scenario}" \
            --log-file  "${LOG_CSV}" \
            --path2prj  "${PRJ_PATH}" \
            "${CMP_EXTRA[@]+"${CMP_EXTRA[@]}"}"; then
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
