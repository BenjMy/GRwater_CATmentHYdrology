#!/usr/bin/env bash
# Typical test scenarios for synthetic_ksfire_LAI.py
# Usage: ./run_synthetic_scenarios.sh [all|smoke|control|lai|fast|slow|dry]   (time*/severity* disabled)
#
# TWO layers (burned top 0.1 m + stable subsoil) are the DEFAULT of the Python script and of
# every scenario below. Scenario folders: $OUT_ROOT/<name>/ (2 layers) or <name>_1L (1 layer).
#   control : Ks constant before/after the fire (LAI still burns)   [alias: baseline]
#   lai     : LAI-driven Ks recovery (default scenario)
#   fast    : fast LAI recovery   (--lai-recovery-rate 2.0)
#   slow    : slow LAI recovery   (--lai-recovery-rate 0.3, one more year)
#   dry     : deeper water table, no plots (quick run)
#   all     : smoke + control + lai + fast + slow + dry
# Layers (env LAYERS):
#   LAYERS=2     (default) two-layer soil: top (Ks from ks_fire_utils) / deep (stable Ks)
#   LAYERS=1     one homogeneous layer  (--no-two-layers), folders get the suffix _1L
#   LAYERS=both  2-layer set first, then the 1-layer set
# Legacy targets kept: control2l|lai2l|twolayer (= control / lai / control+lai, 2 layers)
#                      control1l|lai1l|onelayer (same, 1 layer)
# Env: PY (python), ETP_CSV (Agramon ETp csv), ETP_COL, RAIN_FILE (override forcing)
#      TOP_DEPTH (m, default: the script's 0.1), PERMX_DEEP (m/s, default = unburned Ks)
#      OUT_ROOT (default ../synthetic_LAI_Ks): each scenario goes to $OUT_ROOT/<name>/
#      SAVE_EVERY (default: script's 1) merged outputs are rewritten after every N simulated
#                 months, so nothing is lost if a run stops; 0 = only at the end of the run
#      NO_RESUME=1 re-simulate everything (ignore the monthly backup)
# NOTE: old 1-layer results made before this refactor sit in <name>/ (control, lai, ...); the
#       new 2-layer runs use the same folder names. Move the old ones to <name>_1L (or use a new
#       OUT_ROOT) before running synthetic_ksfire_analysis.py, so the comparison labels stay clear.
set -u
PY="${PY:-python}"
OUT_ROOT="${OUT_ROOT:-../synthetic_LAI_Ks}"
SCRIPT="synthetic_ksfire_LAI.py"
LAYERS="${LAYERS:-2}"
EXTRA=()
[ -n "${ETP_CSV:-}" ]    && EXTRA+=(--etp-csv "$ETP_CSV")
[ -n "${ETP_COL:-}" ]    && EXTRA+=(--etp-ts-column "$ETP_COL")   # ET_0_wet (default) | ET_0_dry
[ -n "${RAIN_FILE:-}" ]  && EXTRA+=(--rain-file "$RAIN_FILE")
[ -n "${NO_RESUME:-}" ]  && EXTRA+=(--no-resume)                  # NO_RESUME=1: re-simulate everything
[ -n "${SAVE_EVERY:-}" ] && EXTRA+=(--save-every "$SAVE_EVERY")
FAILED=(); PARTIAL=(); OK=()
# A scenario that fails or only partly finishes NEVER stops the batch: the next one
# still runs. Exit codes of the python script: 0 = complete, 3 = partial (some months
# not simulated, outputs of the others saved), anything else = failed.

# set_layers 2|1 -> LAYER_ARGS (flags for the python script) and SFX (scenario-name suffix)
set_layers() {
  case "$1" in
    2) SFX=""
       LAYER_ARGS=(--two-layers)                                   # python default, explicit here
       [ -n "${TOP_DEPTH:-}" ]  && LAYER_ARGS+=(--top-depth "$TOP_DEPTH")
       [ -n "${PERMX_DEEP:-}" ] && LAYER_ARGS+=(--permx-deep "$PERMX_DEEP");;
    1) SFX="_1L"
       LAYER_ARGS=(--no-two-layers);;
  esac
}

run() {  # run <name> <args...>   (name gets the layer suffix; LAYER_ARGS are added)
  local name="$1$SFX"; shift
  echo "=== [$name] ($([ -n "$SFX" ] && echo '1 layer' || echo '2 layers')) ==="
  "$PY" "$SCRIPT" --path2prj "${OUT_ROOT}/${name}" ${LAYER_ARGS[@]+"${LAYER_ARGS[@]}"} "$@" \
       ${EXTRA[@]+"${EXTRA[@]}"}
  local rc=$?
  case $rc in
    0) OK+=("$name");;
    3) PARTIAL+=("$name"); echo "!!! [$name] PARTIAL (some months not simulated) - continuing";;
    *) FAILED+=("$name");  echo "!!! [$name] FAILED (exit $rc) - continuing with the next scenario";;
  esac
  return 0
}

"$PY" -m py_compile "$SCRIPT" || { echo "py_compile failed"; exit 1; }

# sc_smoke() { "$PY" "$SCRIPT" --smoke-test && OK+=(smoke) || FAILED+=(smoke); return 0; }

# 1. CONTROL: Ks constant (= unburned Ks) before AND after the fire; LAI still burns
sc_control() { run control --ks-fire-mode off --start-year 2020 --end-year 2022; }
# 2. LAI-driven recovery (default scenario)
sc_lai()     { run lai --ks-fire-mode lai --ks-fire-date 2020-07-01 --start-year 2020 --end-year 2022; }
# 3. Time-driven logistic recovery (Thomas et al. 2021 preset)
# (disabled) sc_time()     { run time --ks-fire-mode time --ks-fire-preset overall_average --ks-fire-date 2020-07-01 --start-year 2020 --end-year 2022; }
# 4. Burn-severity reduction (moderate-high)
# (disabled) sc_severity() { run severity --ks-fire-mode burn_severity --ks-fire-severity 0.6 --ks-fire-date 2020-07-01 --start-year 2020 --end-year 2022; }
# 5. Fast vs slow LAI recovery
sc_fast()    { run fast_recovery --ks-fire-mode lai --lai-recovery-rate 2.0 --start-year 2020 --end-year 2022; }
sc_slow()    { run slow_recovery --ks-fire-mode lai --lai-recovery-rate 0.3 --start-year 2020 --end-year 2023; }
# 6. Deeper water table / no plots (quick run)
sc_dry()     { run deep_wt --ks-fire-mode lai --ic-wtposition 8.0 --no-plots --start-year 2020 --end-year 2021; }

# time / severity scenarios are disabled (commented) for now

# each scenario function runs once per requested layer set (2 first, then 1)
for_layers() {  # for_layers <layer-set: 2|1|both> <scenario function...>
  local sets="$1"; shift
  [ "$sets" = both ] && sets="2 1"
  local L f
  for L in $sets; do
    set_layers "$L"
    for f in "$@"; do "$f"; done
  done
}

case "$LAYERS" in
  2|1|both) ;;
  *) echo "LAYERS must be 2, 1 or both (got: $LAYERS)"; exit 2;;
esac

TARGET="${1:-all}"
case "$TARGET" in
  smoke)                  sc_smoke;;
  control|baseline)       for_layers "$LAYERS" sc_control;;
  lai)                    for_layers "$LAYERS" sc_lai;;
  fast)                   for_layers "$LAYERS" sc_fast;;
  slow)                   for_layers "$LAYERS" sc_slow;;
  dry)                    for_layers "$LAYERS" sc_dry;;
  #time)                  for_layers "$LAYERS" sc_time;;
  #severity)              for_layers "$LAYERS" sc_severity;;
  control2l)              for_layers 2 sc_control;;
  lai2l)                  for_layers 2 sc_lai;;
  twolayer)               for_layers 2 sc_control sc_lai;;
  control1l)              for_layers 1 sc_control;;
  lai1l)                  for_layers 1 sc_lai;;
  onelayer)               for_layers 1 sc_control sc_lai;;
  all)                    sc_smoke
                          for_layers "$LAYERS" sc_control sc_lai sc_fast sc_slow sc_dry;;
                          # + sc_time sc_severity when re-enabled
  *) echo "unknown scenario: $TARGET"; exit 2;;
esac

echo
echo "================ SUMMARY ================"
[ ${#OK[@]}      -gt 0 ] && echo "  complete : ${OK[*]}"
[ ${#PARTIAL[@]} -gt 0 ] && echo "  PARTIAL  : ${PARTIAL[*]}   (see <scenario>/<run>/outputs/failures.csv)"
[ ${#FAILED[@]}  -gt 0 ] && echo "  FAILED   : ${FAILED[*]}   (see <scenario>/<run>/outputs/failures.csv)"
[ ${#FAILED[@]} -gt 0 ] && exit 1
[ ${#PARTIAL[@]} -gt 0 ] && exit 0
echo "All requested scenarios completed."
