#!/usr/bin/env bash
# =============================================================================
# GRwater — Agramon LAI hydrology batch runner
# Simulates dem-plots x PMIN values, crossed with up to three independent
# BC/spin-up axes:
#   - outlet boundary condition : free drainage (open) vs closed (no outflow)
#   - bottom boundary condition : free drainage (open) vs closed (no-flow)
#   - spin-up                   : on vs off
# then post-processes.
#
# Every flag below maps 1:1 onto an Agramon_withLAI_withETp.py CLI argument;
# nothing here is edit-the-script-only. Run with --help for details.
#
# Example call:
#   ./run_scenarii_withETp_withLAI.sh --dem-plot 2 --pmin -5 -10 -20 \
#       --outlet free --bottom closed --spinup on --spinup-cycles 3
#   # -> 3 runs (one per PMIN value), dem-plot 2, outlet free only,
#   #    bottom closed only, spin-up on with 3 cycles
# =============================================================================
set -uo pipefail

# -- Defaults ------------------------------------------------------------
START_YEAR=2016
END_YEAR=2026
PLOTS=(2)                     # dem-plot(s) to run; one 1-9 value per sim run
PMINS=(-10)                   # PMIN value(s) in m; one or more, crossed like --dem-plot
OUTLET_MODE="both"            # free | closed | both
SPINUP_MODE="both"            # on | off | both
BOTTOM_MODE="closed"          # free | closed | both  (default "closed" =
                               # bottom stays no-flow, i.e. unchanged from
                               # before --bottom-flux existed)
OUTLET_FREE_VAL="-1e-6"       # m/s, used when outlet_label=free
BOTTOM_FREE_VAL="-1e-6"       # m/s, used when bottom_label=free
SPINUP_CYCLES=3

SCRIPT_SIM="Agramon_withLAI_withETp.py"
SCRIPT_RES="withLAI_results.py"

# -- Help ----------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Runs Agramon_withLAI_withETp.py once per combination of dem-plot x pmin x
outlet mode x bottom mode x spin-up mode, then (optionally) post-processes.
Every option below is a real flag understood by this script — it forwards
straight to the matching Agramon_withLAI_withETp.py argument.

Simulation period:
  --start-year YYYY     First year of forcing (default: ${START_YEAR})
  --end-year YYYY       Last year, exclusive   (default: ${END_YEAR})

DEM plot(s):
  --dem-plot N [N ...]  One or more dem-plot indices, 1-9 (default: ${PLOTS[*]}).
                         Each value becomes its own simulation run.

PMIN (--pmin under the hood):
  --pmin N [N ...]      One or more PMIN values in m (default: ${PMINS[*]}).
                         Crossed with dem-plot/outlet/bottom/spinup just
                         like --dem-plot, i.e. each value becomes its own
                         simulation run.

Outlet boundary condition (lateral, --outlet-flux under the hood):
  --outlet {free|closed|both}
                         free   : open outlet, ${OUTLET_FREE_VAL} m/s
                         closed : outlet-flux = 0 (fully closed box)
                         both   : run both (default)
  --outlet-flux VAL      Override the "free" outlet flux value in m/s
                         (default: ${OUTLET_FREE_VAL})

Bottom boundary condition (--bottom-flux under the hood):
  --bottom {free|closed|both}
                         free   : free bottom drainage, ${BOTTOM_FREE_VAL} m/s
                         closed : bottom-flux = 0 (no-flow, original default)
                         both   : run both
                         Default: closed only, so existing batches are
                         unaffected unless you opt in.
  --bottom-flux VAL       Override the "free" bottom flux value in m/s
                         (default: ${BOTTOM_FREE_VAL})

Spin-up (--spinup-cycles / --no-spinup under the hood):
  --spinup {on|off|both} on   : spin-up enabled, N cycles
                         off  : --no-spinup
                         both : run both (default)
  --spinup-cycles N      Number of spin-up cycles when on (default: ${SPINUP_CYCLES})

Other:
  --help, -h             Show this help and exit

Runs = (#dem-plot) x (#pmin) x (1 or 2 outlet labels) x (1 or 2 bottom
       labels) x (1 or 2 spinup labels). All output is logged to
       batch_run_YYYYMMDD_HHMMSS_plot<N>_pmin<M>_outlet-<mode>_bottom-<mode>_spinup-<mode>.log
       (tee'd to the terminal too), so different argument combos land in
       distinguishable files instead of only differing by timestamp.

Examples:
  $(basename "$0")
      # default: dem-plot 2, outlet free+closed, bottom closed only,
      # spinup on+off -> 4 runs

  $(basename "$0") --dem-plot 1 2 3 --outlet free --spinup on
      # 3 runs: one per dem-plot, outlet free only, spin-up on only

  $(basename "$0") --pmin -5 -10 -20
      # 3 runs (times whatever --outlet/--bottom/--spinup give): one per
      # PMIN value, dem-plot unchanged

  $(basename "$0") --bottom both --bottom-flux -5e-7
      # adds the free-bottom-drainage axis at a custom magnitude,
      # doubling the run count on top of whatever --outlet/--spinup give
EOF
}

# -- Parse args ------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            usage; exit 0 ;;
        --start-year)
            START_YEAR="$2"; shift 2 ;;
        --end-year)
            END_YEAR="$2"; shift 2 ;;
        --dem-plot)
            shift
            PLOTS=()
            while [[ $# -gt 0 && "$1" != --* ]]; do
                PLOTS+=("$1")
                shift
            done
            if [ ${#PLOTS[@]} -eq 0 ]; then
                echo "ERROR: --dem-plot needs at least one value" >&2
                exit 1
            fi
            ;;
        --pmin)
            shift
            PMINS=()
            while [[ $# -gt 0 && "$1" != --* ]]; do
                PMINS+=("$1")
                shift
            done
            if [ ${#PMINS[@]} -eq 0 ]; then
                echo "ERROR: --pmin needs at least one value" >&2
                exit 1
            fi
            ;;
        --outlet)
            OUTLET_MODE="$2"; shift 2 ;;
        --outlet-flux)
            OUTLET_FREE_VAL="$2"; shift 2 ;;
        --bottom)
            BOTTOM_MODE="$2"; shift 2 ;;
        --bottom-flux)
            BOTTOM_FREE_VAL="$2"; shift 2 ;;
        --spinup)
            SPINUP_MODE="$2"; shift 2 ;;
        --spinup-cycles)
            SPINUP_CYCLES="$2"; shift 2 ;;
        *)
            echo "ERROR: unknown option '$1'" >&2
            usage
            exit 1
            ;;
    esac
done

# -- Validate + expand modes into label lists ------------------------------
case "$OUTLET_MODE" in
    free)   OUTLET_LABELS=(free) ;;
    closed) OUTLET_LABELS=(closed) ;;
    both)   OUTLET_LABELS=(free closed) ;;
    *) echo "ERROR: --outlet must be free|closed|both (got '${OUTLET_MODE}')" >&2; exit 1 ;;
esac

case "$BOTTOM_MODE" in
    free)   BOTTOM_LABELS=(free) ;;
    closed) BOTTOM_LABELS=(closed) ;;
    both)   BOTTOM_LABELS=(free closed) ;;
    *) echo "ERROR: --bottom must be free|closed|both (got '${BOTTOM_MODE}')" >&2; exit 1 ;;
esac

case "$SPINUP_MODE" in
    on)   SPINUP_LABELS=(on) ;;
    off)  SPINUP_LABELS=(off) ;;
    both) SPINUP_LABELS=(on off) ;;
    *) echo "ERROR: --spinup must be on|off|both (got '${SPINUP_MODE}')" >&2; exit 1 ;;
esac

for plot in "${PLOTS[@]}"; do
    if ! [[ "$plot" =~ ^[1-9]$ ]]; then
        echo "ERROR: --dem-plot values must be 1-9 (got '${plot}')" >&2
        exit 1
    fi
done

for pmin in "${PMINS[@]}"; do
    if ! [[ "$pmin" =~ ^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]; then
        echo "ERROR: --pmin values must be numeric (got '${pmin}')" >&2
        exit 1
    fi
done

YEARS="--start-year ${START_YEAR} --end-year ${END_YEAR}"

# Build a filename tag from the resolved plot list + mode selections, so
# different argument combos land in distinguishable log files instead of
# only differing by timestamp.
PLOTS_TAG=$(IFS=-; echo "${PLOTS[*]}")
PMINS_TAG=$(IFS=-; echo "${PMINS[*]}")
LOG_TAG="plot${PLOTS_TAG}_pmin${PMINS_TAG}_outlet-${OUTLET_MODE}_bottom-${BOTTOM_MODE}_spinup-${SPINUP_MODE}"
LOG_FILE="batch_run_$(date +%Y%m%d_%H%M%S)_${LOG_TAG}.log"

# Redirect all stdout and stderr to log file AND terminal
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "Log file      : ${LOG_FILE}"
echo "Started       : $(date)"
echo "dem-plot(s)   : ${PLOTS[*]}"
echo "PMIN value(s) : ${PMINS[*]} m"
echo "outlet modes  : ${OUTLET_LABELS[*]} (free = ${OUTLET_FREE_VAL} m/s)"
echo "bottom modes  : ${BOTTOM_LABELS[*]} (free = ${BOTTOM_FREE_VAL} m/s)"
echo "spinup modes  : ${SPINUP_LABELS[*]} (cycles = ${SPINUP_CYCLES})"
echo

# -- Simulations --------------------------------------------------------------
for plot in "${PLOTS[@]}"; do
  for pmin in "${PMINS[@]}"; do
    for outlet_label in "${OUTLET_LABELS[@]}"; do
        if [ "${outlet_label}" = "free" ]; then
            # FIX: "=" form, not a separate token - argparse's negative-number
            # heuristic (^-\d+$|^-\d*\.\d+$) doesn't match scientific notation
            # (e.g. -1e-6, -1e99), so "--outlet-flux -1e-6" gets misread as an
            # unrecognized "-1e-6" option instead of this flag's value. "=" sidesteps
            # the ambiguity entirely, regardless of the value's notation.
            OUTLET_ARGS="--outlet-flux=${OUTLET_FREE_VAL}"
        else
            OUTLET_ARGS="--outlet-flux 0"
        fi

        for bottom_label in "${BOTTOM_LABELS[@]}"; do
            if [ "${bottom_label}" = "free" ]; then
                # FIX: same "=" reasoning as OUTLET_ARGS above.
                BOTTOM_ARGS="--bottom-flux=${BOTTOM_FREE_VAL}"
            else
                BOTTOM_ARGS="--bottom-flux 0"
            fi

            for spinup_label in "${SPINUP_LABELS[@]}"; do
                if [ "${spinup_label}" = "on" ]; then
                    SPINUP_ARGS="--spinup-cycles ${SPINUP_CYCLES}"
                else
                    SPINUP_ARGS="--no-spinup"
                fi

                echo "============================================================"
                echo "  dem-plot ${plot} | pmin=${pmin} | outlet=${outlet_label} | bottom=${bottom_label} | spinup=${spinup_label} ...  $(date)"
                echo "============================================================"
                # FIX: --pmin="${pmin}", not a separate token - same argparse
                # negative-number ambiguity as OUTLET_ARGS/BOTTOM_ARGS above
                # (e.g. --pmin -1e99 was misread as an unrecognized option).
                python "${SCRIPT_SIM}" ${YEARS} ${SPINUP_ARGS} ${OUTLET_ARGS} ${BOTTOM_ARGS} --dem-plot "${plot}" --pmin="${pmin}"
            done
        done
    done
  done
done

# -- Post-processing ----------------------------------------------------------
# NOTE: sim_index is auto-incremented per run (via log_simulation into the
# shared CSV), so it no longer maps 1:1 to dem-plot number now that each
# plot can produce up to 8 runs per PMIN value (outlet x bottom x spinup,
# times however many --pmin values were given). Each run's folder is now
# named "scenario_<sim_index>_plot<N>_outlet-<label>_bottom-<label>_
# spinup-<label>_pmin-<label>" (see build_scenario_dirname() in
# Agramon_withLAI_withETp.py), so you CAN tell combos apart just by
# listing ../SSHydro_withETp_withLAI/outputs/ — but --scenario still wants
# the bare sim_index. Check simulation_log_withETp_withLAI.csv (or the
# folder names themselves) for the sim_index of each combo before
# re-enabling this loop, and drive it off that list instead of `seq 1 9`.
#
# for scenario in $(seq 1 9); do
#     echo "============================================================"
#     echo "  Post-processing scenario ${scenario} ...  $(date)"
#     echo "============================================================"
#     python "${SCRIPT_RES}" --scenario "${scenario}"
# done

echo
echo "Finished: $(date)"
echo "All done. Log saved to ${LOG_FILE}"
