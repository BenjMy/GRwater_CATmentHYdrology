#!/usr/bin/env bash
# =============================================================================
# run_all_plots.sh
# GRwater project — Agramon study site
#
# Runs both analyse_ET.py and analyse_LAI.py for every plot combination.
#
# Usage
# -----
#   chmod +x run_all_plots.sh
#   ./run_all_plots.sh              # spatial maps default year (2020)
#   ./run_all_plots.sh 2021         # override spatial-map year
#   ./run_all_plots.sh 2021 ET      # only ETa figures
#   ./run_all_plots.sh 2021 LAI     # only LAI figures
#
# Environment
# -----------
#   PYTHON=/path/to/python ./run_all_plots.sh
# =============================================================================
set -euo pipefail

PYTHON="${PYTHON:-python3}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

YEAR="${1:-2020}"
YEAR_ARGS=(--year "$YEAR")

MODE="${2:-ALL}"   # ET | LAI | ALL

# ── Known PlotIDs ─────────────────────────────────────────────────────────────
ALL_PLOTS=(
    Plot2PostFirefascine
    Plot3PostFirefascine
    Plot4Mulching
    Plot5Mulching
    Plot6Control
    Plot7Mulching
    Plot8Mulching
    Plot9Control
)

MULCHING=(Plot4Mulching Plot5Mulching Plot7Mulching Plot8Mulching)
CONTROL=(Plot6Control Plot9Control)
FASCINE=(Plot2PostFirefascine Plot3PostFirefascine)

# ── Helper ────────────────────────────────────────────────────────────────────
run() {
    local script="$1"; local desc="$2"; shift 2
    echo ""
    echo "──────────────────────────────────────────────────────────────"
    echo "  [$script]  $desc"
    echo "──────────────────────────────────────────────────────────────"
    "$PYTHON" "$DIR/$script" "${YEAR_ARGS[@]}" "$@"
}

run_for_script() {
    local script="$1"

    # 1. Individual plots
    for i in "${!ALL_PLOTS[@]}"; do
        run "$script" "Plot index $i — ${ALL_PLOTS[$i]}" --plots "$i"
    done

    # 2. Keyword groups
    run "$script" "Group: Mulching" --plots "${MULCHING[@]}"
    run "$script" "Group: Control"  --plots "${CONTROL[@]}"
    run "$script" "Group: Fascine"  --plots "${FASCINE[@]}"

    # 3. All plots combined
    run "$script" "All plots" --plots all
}

# =============================================================================
# Dispatch
# =============================================================================
case "${MODE^^}" in
    ET)
        echo "Running ETa figures only (year=$YEAR) …"
        run_for_script "analyse_ET.py"
        ;;
    LAI)
        echo "Running LAI figures only (year=$YEAR) …"
        run_for_script "analyse_LAI.py"
        ;;
    ALL)
        echo "Running ETa + LAI figures (year=$YEAR) …"
        run_for_script "analyse_ET.py"
        run_for_script "analyse_LAI.py"
        ;;
    *)
        echo "ERROR: unknown mode '$MODE'. Use ET, LAI, or ALL." >&2
        exit 1
        ;;
esac

echo ""
echo "============================================================"
echo "  All figures generated successfully."
echo "============================================================"
