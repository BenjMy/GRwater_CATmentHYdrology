# Running the scenarios and plotting them (GRwater / Agramon)

All commands run from `lib/` unless stated. Placeholders: `<LIKE_THIS>` = fill in; items marked **[MISSING]** do not exist yet in the repo.
Do **not** use `--forcing era5` for the paper runs (agreed). **Default is now `--forcing etp_ts`** (CSV ETp series, uniform over pixels, + rain raster; options `--etp-ts-column ET_0_wet|ET_0_dry`, `--etp-ts-feature-id <ID>`); `--forcing etp` (satellite ET0 + TP rasters) stays available.

## 0. Setup (once per machine)

```bash
export AGRAMON_ROOT=<NEXTCLOUD_ROOT>              # contains GRwater/, BenCSIC/ ; default = old auto-detection
export AGRAMON_PATH_ET=<DIR_WITH_ET0_TP_NC>       # optional override (…/GRwater/data/satellite/Agramon)
export AGRAMON_LAI_PATH=<DIR_WITH_LAI_INPUTS>     # optional override
export AGRAMON_MODULE_PATH=<PATH>/GRwater_geophy   # optional, analysis scripts
export AGRAMON_WB_PRJ_ROOT=<PATH>/SSHydro_withETp_withLAI   # optional, analysis_ET_compare.py
pip install -e <PATH>/pycathy_wrapper             # install pyCATHY from source (not pip)
pip install tomli                                 # only if Python < 3.11 (run_paper_scenarios.py, scenario.py --config)
python test/test_smoke_agramon.py                 # must end with "smoke test OK"
```

## A. Real catchment (Agramon plots), run + plot

**A1. List / dry-run the paper scenarios (`paper_scenarios.toml`)**
```bash
python run_paper_scenarios.py --list
python run_paper_scenarios.py --dry-run                    # default tag "paper"
```

**A2. Run them** (each call = `run_scenarii_withETp_withLAI_ALL.sh`, which loops dem-plot × PMIN × outlet × bottom × spin-up)
```bash
python run_paper_scenarios.py --only paper_good            # one scenario or tag
python run_paper_scenarios.py --only all --skip experimental --stop-on-error
python run_paper_scenarios.py --only paper_good --extra --forcing etp_ts --etp-ts-column ET_0_dry --etp-ts-feature-id <FEATURE_ID> --start-year 2018
```
Single run, directly:
```bash
python Agramon_withLAI_withETp.py --forcing etp_ts --etp-ts-feature-id <FEATURE_ID> --dem-plot 3 --pmin -50 --outlet-flux 0 --permx 1e-4 \
    --eta-net-neg-mode etp|gross|off --save-monthly --mesh-props --log-level INFO
```
`--eta-net-neg-mode gross` is the recommended mode to expose water limitation (see TO-DO §5.2; decision pending).

Each run writes in `<path2prj>/outputs/<scenario_dir>/`: `*_output.nc`, `monthly_stats.csv`, `run_summary.png`, `water_balance.png`, `simulation_log.md`. One row per run goes to `simulation_log_withETp_withLAI.csv` (`sim_index`).

**A3. Post-process / plot** (also available inside the run: add `--post-plots [--post-et-csv <ETA_CSV>]` to `Agramon_withLAI_withETp.py`, or `--extra --post-plots` with `run_paper_scenarios.py`; `run_summary.png` / `water_balance.png` (incl. the applied-Ks panel) are refreshed every month during the run anyway)
```bash
python post_process.py                                  # every sim_index in the log CSV
python post_process.py --scenario 7 8 --forcing etp_ts --dem-res 5
python post_process.py --scenario 7 --et-csv <ETA_CSV_FROM_HECTOR> --monthly-ts
```
(`./plot_scenarii_withETp_withLAI.sh` is now a thin wrapper around `post_process.py`.)
`post_process.py` calls `plot_run_results.py` (maps, hydro, monthly) and `analysis_ET_compare.py` (ETa/ETp, validation panel Fig I) per scenario; figures go to `../Analysis/ET_compare/<scenario_dir>/`.
`<ETA_CSV_FROM_HECTOR>` **[MISSING]**: Hector's ETa data, `--et-csv` default is `agramon_treatments_ET.csv`.

Fig. 8 (ETa/ETp vs shallow psi) is now produced by `plot_run_results.py` for every scenario (Fig 7k, code in `plot_eta.py`). `plot_eta.py`, `plot_maps.py`, `plot_monthly.py` are libraries used by `plot_run_results.py`, not CLIs.

## B. Post-fire Ks curves (theory, no run needed)

`plot_ks_scenarios.py` reads **no simulation output**: it only draws the Ks(t) equations (`ks_fire_utils`) for the dates you give, so it is not part of A3 (which post-processes finished runs). Use it to choose/justify the scenarios before running them.

```bash
python plot_ks_scenarios.py --fire-date <YYYY-MM-DD> --start <YYYY-MM-DD> --end <YYYY-MM-DD> --out ks_scenarios.png
```
Draws Ks(t) for all modes (time / lai / burn_severity) and presets (upper_bound, median, lower_bound, overall_average).
**[MISSING]**: final scenario list and LAI anchors (TO-DO §4), control catchment `<CONTROL_FEATURE_ID>`, recovery thresholds (`recovery_metrics.py`: tol 0.1 / sustain 3 are placeholders).

## C. Synthetic Ks/fire study (no data needed), run + plot

```bash
cd synthetic_LAI_Ks
python synthetic_ksfire_LAI.py --smoke-test                 # one tiny run
python scenario.py --list
python scenario.py --dry-run
python scenario.py --root <OUT_ROOT> --jobs 3                       # default set + comparison
python scenario.py --root <OUT_ROOT> --only control_2L control_ks_median_2L time_median_2L control_ks_lower_2L time_lower_2L lai_2L --wet   # each time_* preset needs its control_ks_* (same unburned Ks)
python scenario.py --root <OUT_ROOT> --quick --wet                  # 5x5 grid sanity run
python scenario.py --root <OUT_ROOT> --config <MY_SCENARIOS.toml>   # custom, see scenario.py docstring
```
Plot / compare (also run automatically at the end of `scenario.py` unless `--no-compare`):
```bash
python synthetic_ksfire_analysis.py --root <OUT_ROOT> --out <OUT_ROOT>/_comparison --fire-date <YYYY-MM-DD> --exclude deep_wt
python synthetic_ksfire_analysis.py --slopes <ROOT_SLOPE1> <ROOT_SLOPE2> --out <OUT_SLOPES>
python scenario.py --root <OUT_ROOT> --compare-only
```
Output `<OUT_ROOT>/_comparison/`: `fig_cmp_*.png` (incl. `fig_cmp_eta_variants.png`: net / raw / gross ETa), `summary_*.csv`.
`scenarios.toml` is mentioned in `synthetic_LAI_Ks/README.md` but is **[MISSING]**: create `<MY_SCENARIOS.toml>` if you need extras.

## D. Rain-via-irrigation A/B test (withIrr build)

```bash
cd lib/ab_withirr
PYTHONPATH=<PATH_TO_WITHIRR_CHECKOUT> python run_ab.py A <OUT>/A
PYTHONPATH=<PATH_TO_WITHIRR_CHECKOUT> python run_ab.py B <OUT>/B
python compare_ab.py <OUT>/A <OUT>/B
```
Result (TO-DO log, 2026-11-15): B reproduces A, so it does not remove the ETa artefact.
