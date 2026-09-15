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
PMINS=(-1e35)                 # PMIN value(s) in m; one or more, crossed like --dem-plot
OUTLET_MODE="closed"          # free | closed | both  (default "closed" so a
                               # bare run is a fully-closed box, matching
                               # BOTTOM_MODE's default below)
SPINUP_MODE="off"             # on | off | both  (default "off", matching
                               # OUTLET_MODE/BOTTOM_MODE now defaulting to
                               # closed — a bare run does the plainest case)
BOTTOM_MODE="closed"          # free | closed | both  (default "closed" =
                               # bottom stays no-flow, i.e. unchanged from
                               # before --bottom-flux existed)
OUTLET_FREE_VAL="-1e-6"       # m/s, used when outlet_label=free
BOTTOM_FREE_VAL="-1e-6"       # m/s, used when bottom_label=free
SPINUP_CYCLES=3
IC_MODE="both"                # wt | pressure | both
IC_WTPOSITION="1"             # m, used when ic_label=wt
IC_PRESSURE_HEAD="-10"        # m, used when ic_label=pressure

# Kept in sync with Agramon_withLAI_withETp.py's own defaults for the
# flags this script never overrides — needed so the "already succeeded?"
# check below looks at the same CSV log / project path the simulation
# script itself will actually use.
ZROOT=1.0                             # matches --zroot default
LOG_FILE_CHECK="simulation_log_withETp_withLAI.csv"   # matches --log-file default
PRJ_PATH_CHECK="../SSHydro_withETp_withLAI"            # matches --path2prj default

# Root-depth mapping / Ks (PERMX) overrides — forwarded as-is to
# Agramon_withLAI_withETp.py. Empty PERMX means "leave pyCATHY's own
# SPP default untouched" (same as the Python script's own default).
ZROOT_FACTORS=()                      # 5 values, empty = use Python script's default
PERMX=""                              # m/s, empty = pyCATHY default (not overridden)
PERMY_RATIO="1.0"
PERMZ_RATIO="1.0"

FORCE_RERUN=0                 # --force: run even if a matching completed
                               # scenario is already found on disk

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
                         both   : run both
                         Default: closed only, so existing batches are
                         unaffected unless you opt in.
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
                         both : run both
                         Default: off only, so existing batches are
                         unaffected unless you opt in.
  --spinup-cycles N      Number of spin-up cycles when on (default: ${SPINUP_CYCLES})

Initial condition (--ic-mode under the hood):
  --ic {wt|pressure|both}
                         wt       : hydrostatic water table, depth
                                    ${IC_WTPOSITION} m (INDP=3, --ic-wtposition)
                         pressure : uniform pressure head, ${IC_PRESSURE_HEAD} m
                                    (INDP=0, --ic-pressure-head)
                         both     : run both (default)
  --ic-wtposition VAL    Override the water-table depth in m, used when
                         ic_label=wt (default: ${IC_WTPOSITION})
  --ic-pressure-head VAL Override the uniform initial pressure head in m,
                         used when ic_label=pressure (default: ${IC_PRESSURE_HEAD})

Root depth mapping (--zroot-factors under the hood):
  --zroot-factors B S M D V
                         5 root-depth multipliers (bare, sparse, moderate,
                         dense, very dense), applied as ZROOT = --zroot *
                         factor. Default: whatever
                         Agramon_withLAI_withETp.py itself defaults to
                         (not overridden by this script unless given).
                         Same value used for every run in the batch —
                         it is not crossed like --dem-plot/--pmin.

Soil hydraulic conductivity (--permx / --permy-ratio / --permz-ratio):
  --permx VAL            Override saturated Ks in the x direction (m/s,
                          SOIL SPP column PERMX) for every soil zone.
                          Default: unset, i.e. leave pyCATHY's own SPP
                          default untouched. Same value used for every
                          run in the batch — it is not crossed.
  --permy-ratio VAL      PERMY = --permx * this value (default: ${PERMY_RATIO}).
                          Only takes effect when --permx is set.
  --permz-ratio VAL      PERMZ = --permx * this value (default: ${PERMZ_RATIO}).
                          Only takes effect when --permx is set.

Other:
  --force                Re-run a scenario even if a completed one with the
                          exact same parameters (incl. initial condition,
                          root-depth factors, and PERMX/Ks) is already
                          found on disk. Default: skip it.
  --help, -h             Show this help and exit

Runs = (#dem-plot) x (#pmin) x (1 or 2 outlet labels) x (1 or 2 bottom
       labels) x (1 or 2 spinup labels) x (1 or 2 ic labels). Scenarios
       that already completed successfully (same params, incl. initial
       condition) are skipped automatically — pass --force to re-run
       them anyway. All output is logged to
       batch_run_YYYYMMDD_HHMMSS_plot<N>_pmin<M>_outlet-<mode>_bottom-<mode>_spinup-<mode>_ic-<mode>.log
       (tee'd to the terminal too), so different argument combos land in
       distinguishable files instead of only differing by timestamp.

Examples:
  $(basename "$0")
      # default: dem-plot 2, outlet closed only, bottom closed only,
      # spinup off only -> 1 run

  $(basename "$0") --dem-plot 1 2 3 --outlet free --spinup on
      # 3 runs: one per dem-plot, outlet free only, spin-up on only

  $(basename "$0") --pmin -5 -10 -20
      # 3 runs (times whatever --outlet/--bottom/--spinup give): one per
      # PMIN value, dem-plot unchanged

  $(basename "$0") --bottom both --bottom-flux -5e-7
      # adds the free-bottom-drainage axis at a custom magnitude,
      # doubling the run count on top of whatever --outlet/--spinup give

  $(basename "$0") --pmin -15 --outlet free --spinup on
      # PMIN=-15 m, outlet free (draining), spin-up on -> 1 run

  $(basename "$0") --pmin -15 --outlet closed --spinup off
      # PMIN=-15 m, outlet closed, no spin-up -> 1 run

  $(basename "$0") --pmin -1e35 --outlet free --spinup on
      # PMIN=-1e35 m (effectively unbounded suction), outlet free,
      # spin-up on -> 1 run

  $(basename "$0") --pmin -1e35 --outlet closed --spinup off
      # same PMIN, outlet closed, no spin-up -> 1 run

  $(basename "$0") --pmin -15 -25 -1e35 --outlet both --spinup both
      # full 3x2x2 sweep: PMIN in {-15,-25,-1e35} x outlet {free,closed}
      # x spinup {on,off} -> 12 runs

  $(basename "$0") --ic wt --ic-wtposition 0.1
      # water-table IC at 0.1 m depth only (times whatever other axes give)

  $(basename "$0") --ic pressure --ic-pressure-head -10
      # uniform pressure-head IC of -10 m only

  $(basename "$0") --ic both
      # runs both IC options (wt @ ${IC_WTPOSITION} m and pressure @
      # ${IC_PRESSURE_HEAD} m) -> doubles the run count on top of whatever
      # --outlet/--bottom/--spinup give

  for ks in 1e-4 1e-5 1e-6; do
      $(basename "$0") --pmin -1e34 --outlet closed --bottom closed \\
          --ic both --zroot-factors 0.1 0.3 1.0 1.5 2.5 --permx "\${ks}"
  done
      # varying root depth (--zroot-factors) x varying Ks (--permx, one
      # value per loop iteration since it is NOT crossed like --pmin),
      # PMIN=-1e34, ic both, outlet+bottom both closed -> 4 runs per Ks
      # value (12 runs total). NOTE: --pmin/--permx here are THIS script's
      # own flags and take a space-separated value ("--pmin -1e34", NOT
      # "--pmin=-1e34") — the "=" form is only used internally, when this
      # script forwards values on to Agramon_withLAI_withETp.py's argparse
      # CLI (see FIX comments near OUTLET_ARGS below).
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
        --ic)
            IC_MODE="$2"; shift 2 ;;
        --ic-wtposition)
            IC_WTPOSITION="$2"; shift 2 ;;
        --ic-pressure-head)
            IC_PRESSURE_HEAD="$2"; shift 2 ;;
        --zroot-factors)
            shift
            ZROOT_FACTORS=()
            while [[ $# -gt 0 && "$1" != --* ]]; do
                ZROOT_FACTORS+=("$1")
                shift
            done
            if [ ${#ZROOT_FACTORS[@]} -ne 5 ]; then
                echo "ERROR: --zroot-factors needs exactly 5 values (bare sparse moderate dense very_dense), got ${#ZROOT_FACTORS[@]}" >&2
                exit 1
            fi
            ;;
        --permx)
            PERMX="$2"; shift 2 ;;
        --permy-ratio)
            PERMY_RATIO="$2"; shift 2 ;;
        --permz-ratio)
            PERMZ_RATIO="$2"; shift 2 ;;
        --force)
            FORCE_RERUN=1; shift ;;
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

case "$IC_MODE" in
    wt)       IC_LABELS=(wt) ;;
    pressure) IC_LABELS=(pressure) ;;
    both)     IC_LABELS=(wt pressure) ;;
    *) echo "ERROR: --ic must be wt|pressure|both (got '${IC_MODE}')" >&2; exit 1 ;;
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

# Root-depth / Ks overrides, forwarded verbatim if given (same value for
# every run in the batch — not crossed like --dem-plot/--pmin).
SOIL_ARGS=""
if [ ${#ZROOT_FACTORS[@]} -eq 5 ]; then
    SOIL_ARGS="${SOIL_ARGS} --zroot-factors ${ZROOT_FACTORS[*]}"
fi
if [ -n "${PERMX}" ]; then
    # "=" form, same negative-number-notation reasoning as
    # OUTLET_ARGS/BOTTOM_ARGS/--pmin above (e.g. --permx -1e-5).
    SOIL_ARGS="${SOIL_ARGS} --permx=${PERMX} --permy-ratio=${PERMY_RATIO} --permz-ratio=${PERMZ_RATIO}"
fi

# Build a filename tag from the resolved plot list + mode selections, so
# different argument combos land in distinguishable log files instead of
# only differing by timestamp.
PLOTS_TAG=$(IFS=-; echo "${PLOTS[*]}")
PMINS_TAG=$(IFS=-; echo "${PMINS[*]}")
LOG_TAG="plot${PLOTS_TAG}_pmin${PMINS_TAG}_outlet-${OUTLET_MODE}_bottom-${BOTTOM_MODE}_spinup-${SPINUP_MODE}_ic-${IC_MODE}"
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
echo "ic modes      : ${IC_LABELS[*]} (wt depth = ${IC_WTPOSITION} m, pressure = ${IC_PRESSURE_HEAD} m)"
if [ ${#ZROOT_FACTORS[@]} -eq 5 ]; then
    echo "zroot-factors : ${ZROOT_FACTORS[*]}"
else
    echo "zroot-factors : (script default)"
fi
if [ -n "${PERMX}" ]; then
    echo "PERMX (Ks)    : ${PERMX} m/s  (PERMY ratio ${PERMY_RATIO}, PERMZ ratio ${PERMZ_RATIO})"
else
    echo "PERMX (Ks)    : (pyCATHY default)"
fi
if [ "${FORCE_RERUN}" = "1" ]; then
    echo "skip-existing : disabled (--force given, all combos will run)"
else
    echo "skip-existing : enabled (already-completed combos will be skipped)"
fi
echo

# -----------------------------------------------------------------------------
# check_already_done — has this exact parameter combination already run to
# completion?
#
# Looks up ${LOG_FILE_CHECK} for a logged row matching every parameter that
# feeds into a scenario's identity (years, ZROOT, PMIN, plot, outlet/bottom
# flux, spinup, and initial condition), then confirms the corresponding
# scenario_<sim_index>_* folder under ${PRJ_PATH_CHECK}/outputs actually
# contains the three NetCDFs Agramon_withLAI_withETp.py only ever writes at
# the very end of a successful run (et/psi/sw_output.nc). Being logged is
# not enough on its own — a run can be interrupted after logging but before
# finishing, so both checks have to pass.
#
# Prints "DONE" (and the matched folder) or "MISSING" on stdout; callers
# should test the *last line*.
# -----------------------------------------------------------------------------
check_already_done() {
    CHK_LOG_FILE="${LOG_FILE_CHECK}" \
    CHK_PRJ_PATH="${PRJ_PATH_CHECK}" \
    CHK_START_YEAR="${START_YEAR}" \
    CHK_END_YEAR="${END_YEAR}" \
    CHK_ZROOT="${ZROOT}" \
    CHK_PMIN="$1" \
    CHK_PLOT="$2" \
    CHK_OUTLET_FLUX="$3" \
    CHK_BOTTOM_FLUX="$4" \
    CHK_SPINUP_STR="$5" \
    CHK_IC_MODE="$6" \
    CHK_ZROOT_FACTORS="$(IFS=,; echo "${ZROOT_FACTORS[*]}")" \
    CHK_PERMX="${PERMX}" \
    CHK_PERMY_RATIO="${PERMY_RATIO}" \
    CHK_PERMZ_RATIO="${PERMZ_RATIO}" \
    python3 - <<'PYEOF'
import math
import os
from pathlib import Path

import pandas as pd

log_path = Path(os.environ["CHK_LOG_FILE"])
if not log_path.exists():
    print("MISSING (no log file yet)")
    raise SystemExit

start_year  = int(os.environ["CHK_START_YEAR"])
end_year    = int(os.environ["CHK_END_YEAR"])
zroot       = float(os.environ["CHK_ZROOT"])
pmin        = float(os.environ["CHK_PMIN"])
plot        = int(os.environ["CHK_PLOT"])
outlet_flux = float(os.environ["CHK_OUTLET_FLUX"])
bottom_flux = float(os.environ["CHK_BOTTOM_FLUX"])
spinup_str  = os.environ["CHK_SPINUP_STR"]
ic_mode     = os.environ["CHK_IC_MODE"]
zroot_factors_str = os.environ.get("CHK_ZROOT_FACTORS", "")
permx       = os.environ.get("CHK_PERMX", "")
permy_ratio = os.environ.get("CHK_PERMY_RATIO", "")
permz_ratio = os.environ.get("CHK_PERMZ_RATIO", "")

df = pd.read_csv(log_path)

def close(series: pd.Series, value: float) -> pd.Series:
    return series.astype(float).apply(lambda v: math.isclose(v, value, rel_tol=1e-9, abs_tol=1e-9))

mask = (
    (df["start_year"] == start_year)
    & (df["end_year"] == end_year)
    & close(df["ZROOT"], zroot)
    & close(df["PMIN"], pmin)
    & (df["watershed_nb"] == plot)
    & close(df["outlet_flux"], outlet_flux)
    & close(df["bottom_flux"], bottom_flux)
    & (df["spinup"].astype(str) == spinup_str)
)
if "withLAI" in df.columns:
    mask &= (df["withLAI"] == 1)
if "ic_mode" in df.columns:
    # Rows logged before ic_mode existed won't match a specific ic_mode —
    # treated as "not confirmed for this axis", so they're re-run rather
    # than risk skipping the wrong initial condition.
    mask &= (df["ic_mode"].astype(str) == ic_mode)

# --zroot-factors / --permx (+ ratios): only enforced when this batch run
# actually overrides them (non-empty env var). Empty means "use whatever
# Agramon_withLAI_withETp.py itself defaults to" here too, so we don't
# require a specific logged value in that case — but a batch run that DID
# override Ks/root-depth must not be satisfied by a logged row that used
# a different (or no) override, or vice versa.
if zroot_factors_str and "zroot_factors" in df.columns:
    mask &= (df["zroot_factors"].astype(str) == zroot_factors_str)
if permx:
    if "permx" in df.columns:
        mask &= close(df["permx"], float(permx))
        if "permy_ratio" in df.columns:
            mask &= close(df["permy_ratio"], float(permy_ratio))
        if "permz_ratio" in df.columns:
            mask &= close(df["permz_ratio"], float(permz_ratio))
else:
    if "permx" in df.columns:
        mask &= df["permx"].isna()

matches = df[mask]
if matches.empty:
    print("MISSING (no matching logged run)")
    raise SystemExit

sim_index = int(matches.iloc[-1]["sim_index"])
out_root  = Path(os.environ["CHK_PRJ_PATH"]) / "outputs"
required  = ("et_output.nc", "psi_output.nc", "sw_output.nc")

for d in sorted(out_root.glob(f"scenario_{sim_index}_*")):
    if d.is_dir() and all((d / f).exists() for f in required):
        print(f"DONE {d}")
        raise SystemExit
# also allow the untagged legacy folder name
legacy = out_root / f"scenario_{sim_index}"
if legacy.is_dir() and all((legacy / f).exists() for f in required):
    print(f"DONE {legacy}")
    raise SystemExit

print(f"MISSING (sim_index={sim_index} logged but outputs incomplete)")
PYEOF
}

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
            OUTLET_FLUX_VAL="${OUTLET_FREE_VAL}"
            OUTLET_ARGS="--outlet-flux=${OUTLET_FLUX_VAL}"
        else
            OUTLET_FLUX_VAL="0"
            OUTLET_ARGS="--outlet-flux 0"
        fi

        for bottom_label in "${BOTTOM_LABELS[@]}"; do
            if [ "${bottom_label}" = "free" ]; then
                # FIX: same "=" reasoning as OUTLET_ARGS above.
                BOTTOM_FLUX_VAL="${BOTTOM_FREE_VAL}"
                BOTTOM_ARGS="--bottom-flux=${BOTTOM_FLUX_VAL}"
            else
                BOTTOM_FLUX_VAL="0"
                BOTTOM_ARGS="--bottom-flux 0"
            fi

            for spinup_label in "${SPINUP_LABELS[@]}"; do
                if [ "${spinup_label}" = "on" ]; then
                    SPINUP_ARGS="--spinup-cycles ${SPINUP_CYCLES}"
                    SPINUP_STR="on_${SPINUP_CYCLES}cycles"   # must match the
                                                              # string logged by
                                                              # Agramon_withLAI_withETp.py
                else
                    SPINUP_ARGS="--no-spinup"
                    SPINUP_STR="off"
                fi

                for ic_label in "${IC_LABELS[@]}"; do
                    if [ "${ic_label}" = "wt" ]; then
                        IC_ARGS="--ic-mode wt --ic-wtposition=${IC_WTPOSITION}"
                    else
                        IC_ARGS="--ic-mode pressure --ic-pressure-head=${IC_PRESSURE_HEAD}"
                    fi

                    echo "============================================================"
                    echo "  dem-plot ${plot} | pmin=${pmin} | outlet=${outlet_label} | bottom=${bottom_label} | spinup=${spinup_label} | ic=${ic_label} ...  $(date)"
                    echo "============================================================"

                    if [ "${FORCE_RERUN}" != "1" ]; then
                        status=$(check_already_done \
                            "${pmin}" "${plot}" "${OUTLET_FLUX_VAL}" \
                            "${BOTTOM_FLUX_VAL}" "${SPINUP_STR}" "${ic_label}")
                        if [[ "${status}" == DONE* ]]; then
                            echo "  ⏭  Already completed — skipping (${status#DONE }). Use --force to re-run."
                            echo
                            continue
                        else
                            echo "  → ${status} — running."
                        fi
                    fi

                    # FIX: --pmin="${pmin}", not a separate token - same argparse
                    # negative-number ambiguity as OUTLET_ARGS/BOTTOM_ARGS above
                    # (e.g. --pmin -1e99 was misread as an unrecognized option).
                    python "${SCRIPT_SIM}" ${YEARS} ${SPINUP_ARGS} ${OUTLET_ARGS} ${BOTTOM_ARGS} ${IC_ARGS} ${SOIL_ARGS} --dem-plot "${plot}" --pmin="${pmin}"
                done
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

# =============================================================================
# Extra scenario batch: PMIN in {-15, -25, -1e35} crossed with outlet
# {free,closed} and spinup {on,off} (bottom left at the script's default:
# closed).
# This block is OPT-IN: it only runs if you launch this script with
# RUN_EXTRA_PMIN_SCENARIOS=1, e.g.:
#   RUN_EXTRA_PMIN_SCENARIOS=1 ./run_scenarii_withETp_withLAI.sh
# It is independent of the flag-parsing above and does NOT respect any
# --pmin/--outlet/--spinup values you passed on the command line — it always
# runs exactly the 12 combinations below. Comment out any calls you don't
# want, or copy a single line out to run it standalone.
# =============================================================================
if [ "${RUN_EXTRA_PMIN_SCENARIOS:-0}" = "1" ]; then
    echo
    echo "############################################################"
    echo "  Extra scenario batch: PMIN {-15,-25,-1e35} x outlet x spinup"
    echo "############################################################"

    ./"$(basename "$0")" --pmin -15    --outlet free   --spinup on
    ./"$(basename "$0")" --pmin -15    --outlet free   --spinup off
    ./"$(basename "$0")" --pmin -15    --outlet closed --spinup on
    ./"$(basename "$0")" --pmin -15    --outlet closed --spinup off
    ./"$(basename "$0")" --pmin -25    --outlet free   --spinup on
    ./"$(basename "$0")" --pmin -25    --outlet free   --spinup off
    ./"$(basename "$0")" --pmin -25    --outlet closed --spinup on
    ./"$(basename "$0")" --pmin -25    --outlet closed --spinup off
    ./"$(basename "$0")" --pmin -1e35  --outlet free   --spinup on
    ./"$(basename "$0")" --pmin -1e35  --outlet free   --spinup off
    ./"$(basename "$0")" --pmin -1e35  --outlet closed --spinup on
    ./"$(basename "$0")" --pmin -1e35  --outlet closed --spinup off

    echo
    echo "Extra scenario batch finished: $(date)"
fi
