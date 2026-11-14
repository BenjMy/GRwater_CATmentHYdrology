# Synthetic Ks / fire scenarios (GRwater - Agramon)

Batch driver for the synthetic 20x20 pyCATHY test of fire-driven saturated hydraulic
conductivity (Ks) recovery.

| File | Role |
|---|---|
| `synthetic_ksfire_LAI.py` | One simulation (one scenario, one set of CLI arguments) |
| `synthetic_ksfire_analysis.py` | Area-mean diagnostics per run + comparison across runs |
| `ks_fire_utils.py` | Ks recovery equations (`time`, `lai`, `burn_severity`); imported by the run script |
| `scenario.py` | Runs many scenarios (each in its own subprocess) and builds the comparison |
| `scenarios.toml` | Extra / custom scenarios loaded with `--config` |

`scenario.py`, `scenarios.toml`, `synthetic_ksfire_LAI.py` and `synthetic_ksfire_analysis.py`
must sit in the same folder. `ks_fire_utils.py` is imported by `synthetic_ksfire_LAI.py` as
`from utils import ks_fire_utils`, so it must be reachable under that path (see the
`sys.path.insert` at the top of that script).

## Quick start

```bash
python scenario.py --list                  # what exists, with tags and arguments
python scenario.py --dry-run               # print the commands, run nothing
python scenario.py                         # run the default set, then compare
```

Before a long batch, check the single-run script works:

```bash
python synthetic_ksfire_LAI.py --smoke-test
```

## What a run does

Each scenario runs

```
python synthetic_ksfire_LAI.py --path2prj <root>/<name> --start-year ... --end-year ... <scenario args>
```

Outputs land in:

```
<root>/                                  (default <repo>/results/synthetic_LAI_Ks, outside lib/)
  <name>/
    run.log                              full console output of that scenario
    _monthly_backup/<scenario>/          one file per simulated month, enables restart
    synth_ks-<mode>[_2L]_<y0>-<y1>/
      outputs/                           psi/sw/et/recharge .nc, ks_fire_monthly.csv,
                                         areamean_timeseries.csv, areamean_monthly.csv,
                                         areamean_periods.csv, run_meta.json, failures.csv
      figures/                           per-run figures
    ...
  _comparison/                           cross-scenario figures and summary tables
```

## Built-in scenarios

Tag `default` is what runs when you give no `--only`.

| Name | Tags | Ks | Soil |
|---|---|---|---|
| `control` | default, control | constant (`off`) | 1 layer |
| `control_2L` | default, control | constant (`off`) | 2 layers |
| `lai_1L` | default, lai | from LAI recovery ratio | 1 layer |
| `lai_2L` | default, lai | from LAI recovery ratio | 2 layers |
| `lai_raw_2L` | lai | from raw (seasonal) LAI | 2 layers |
| `time_1L` | default, time | logistic, `overall_average` preset | 1 layer |
| `time_2L` | default, time | logistic, `overall_average` preset | 2 layers |
| `time_median_2L` | time | logistic, `median` preset | 2 layers |
| `time_lower_2L` | time | logistic, `lower_bound` preset | 2 layers |
| `time_upper_2L` | time | logistic, `upper_bound` preset | 2 layers |
| `burnsev_1L` | burnsev | static, severity 0.6 | 1 layer |
| `burnsev_2L` | burnsev | static, severity 0.6 | 2 layers |
| `lai_2L_top0.5` | lai, geometry | LAI | 2 layers, 0.5 m burned topsoil |
| `lai_2L_deepK5e-5` | lai, geometry | LAI | 2 layers, deep Ks 5e-5 m/s |
| `deep_wt` | ic | LAI | hydrostatic IC (not compared) |

Other tags: `ks` (all Ks-mode scenarios except geometry/ic), `layers1`, `layers2`.

Ks modes (see `ks_fire_utils.py`):

- `off`: Ks constant before and after the fire. Vegetation (LAI) still burns. This is the control.
- `lai`: Ks linearly rescaled from the monthly LAI (Thomas et al. 2021, Eq. 11). By default
  the LAI is expressed as a recovery ratio relative to the unburned LAI of the same calendar
  month, so seasonality cancels (`--lai-ref raw` keeps the seasonal LAI).
- `time`: Ks follows a logistic recovery in time since the fire (Ebel & Martin 2017;
  Thomas et al. 2021, Eq. 9).
- `burn_severity`: single static Ks reduction scaled by a severity fraction (placeholder
  linear relation in the spirit of Atchley et al. 2018).

Soil layers:

- 2 layers (default): a thin burned topsoil (`--top-depth`, default 0.1 m) whose Ks follows
  the fire model every month, over a stable subsoil (`--permx-deep`, default = unburned Ks).
- 1 layer (`--no-two-layers`): one homogeneous soil.

## Hypotheses: what each scenario should do (theory)

The synthetic fire removes LAI on 2020-07-01 (root-water uptake and ETa drop at once, in all
scenarios including the controls); LAI then regrows (`--lai-recovery-rate`). The scenarios differ
only in how Ks responds. Expected sign of scenario − matching control:

| Scenario | Ks after the fire | Expected ΔETa | Expected Δrecharge |
|---|---|---|---|
| `control`, `control_2L` | constant | 0 (reference; shows the pure vegetation effect) | 0 |
| `lai_1L` / `lai_2L` | drops to the post-fire minimum (16.2 mm/hr, 37 % of unburned), recovers with the LAI recovery ratio | small decrease or none: less infiltration, more ponding and runoff, drier topsoil in storms | decrease right after the fire, converging to 0 as LAI regrows. Stronger in 1L (whole 3 m column slowed) than 2L (only 0.1 m topsoil: it throttles infiltration but not deep drainage) |
| `lai_raw_2L` | follows the seasonal absolute LAI | as `lai`, plus a seasonal Ks oscillation (artefact of not normalising by the unburned LAI) | as `lai` with seasonal noise |
| `time_*` | logistic recovery in time since fire, independent of vegetation (`overall_average`, `median` fast, `lower_bound`/`upper_bound` extremes) | as `lai`, but the recovery speed is set by the preset, not by regrowth | as `lai`; `lower_bound` recovers slowest, so the longest deficit |
| `burnsev_*` | one static reduction (severity 0.6), no recovery | persistent small decrease | persistent decrease, no convergence to 0 |
| `lai_2L_top0.5` | as `lai_2L`, 0.5 m burned layer | slightly larger effect than 0.1 m | larger decrease (thicker low-Ks layer) |
| `lai_2L_deepK5e-5` | as `lai_2L`, stiffer subsoil | negligible | decrease in recharge (deep drainage limited), independent of the fire |

**Caveat.** The forcing is daily: even 3x rain stays well below the lowest post-fire Ks (16.2 mm/hr
≈ 390 mm/day), so infiltration capacity is rarely exceeded and Ks acts mainly through
redistribution and recharge, not through infiltration-excess runoff. Expect small Δ unless the
soil is wet or Ks is lowered further (`--ks-fire-ll`).

## Selecting scenarios

`--only` takes names and/or tags; `--skip` removes names and/or tags afterwards.

```bash
python scenario.py --only control control_2L lai_2L          # three named scenarios
python scenario.py --only time                               # every scenario tagged "time"
python scenario.py --only layers2 --skip geometry            # all 2-layer runs except geometry tests
python scenario.py --only all --skip ic                      # everything except the IC test
python scenario.py --only ks --skip burnsev                  # all Ks modes, no burn severity
```

An unknown name or tag stops the script and prints the valid names and tags.

## Common options

```bash
# period (end year is exclusive: 2020-2023 simulates 2020, 2021, 2022, 2023)
python scenario.py --start-year 2020 --end-year 2024

# different output root
python scenario.py --root /data/synthetic_LAI_Ks     # or: export SYNTH_KS_ROOT=/data/synthetic_LAI_Ks

# run 3 scenarios in parallel (each is a separate CATHY run; mind CPU/RAM/disk)
python scenario.py --only all --skip ic --jobs 3

# quick test: 5x5 grid (no spatial heterogeneity), own root results/synthetic_LAI_Ks_quick
python scenario.py --quick --only control control_2L lai_2L

# quick test: no figures
python scenario.py --only control control_2L lai_2L --no-plots

# stop at the first scenario that fails (sequential mode only)
python scenario.py --stop-on-error

# do not build _comparison/ at the end
python scenario.py --no-compare

# only rebuild _comparison/ from runs already on disk
python scenario.py --compare-only
python scenario.py --compare-only --compare-exclude lai_2L_wetIC deep_wt
```

### Wetter forcing and initial condition (`--wet`)

With the default forcing (rain ~1-3 mm/day against ETp 1-4.5 mm/day) and a uniform initial
head of -15 m, the soil stays dry and Ks has almost nothing to act on. `--wet` appends

```
--rain-scale 3 --ic-wt-depth 2.5 --pmin -15 --eta-net-neg-mode gross
```

to every run (before the scenario's own arguments and `--extra`, so both can override it) and
writes to its own root, `<repo>/results/synthetic_LAI_Ks_wet`, so the controls match. The IC is
hydrostatic with the water table 2.5 m below the surface (the mesh is 3 m deep); the earlier
uniform head of -1 m saturated the whole column (Sw = 1 at 1 m and 2 m). `--ic-wt-depth M`
sets that depth in a single run (`--ic-wtposition` is the same thing measured from the grid base).
Note the mesh is only 3 m deep, so "2 m" in the soil figures is really ~1.7 m (the figure
label says so). Never mix wet and
normal runs in the same root.

```bash
python scenario.py --wet --only control control_2L lai_2L time_2L
python scenario.py --wet --extra --rain-scale 5 --ic-head -0.5     # even wetter
```

`--rain-scale X` / `--etp-scale X` (single-run script) multiply the daily forcing; the run
prints the resulting annual rain and ETp. The default period is now 2020-2023
(`--end-year 2024`, 3.5 years after the 2020-07-01 fire).

### DEM slope (`--slope`)

The default DEM is flat (outlet = one cell 1 mm lower). `--slope S` (m/m, 0.1 = 10 %) tilts it
along the row index: row 0 is the lowest side and the single outlet; the elevation step per
cell is `S * --res`. Each slope has its own root (`..._slope<S>`), so slopes never mix:

```bash
for s in 0 0.05 0.15; do python scenario.py --wet --slope $s --only control_2L lai_2L; done
python scenario.py --quick --slope 0.1 --only control --dry-run     # check the command
```

Caveats: the hydrostatic IC (`--ic-wt-depth`) puts a horizontal water table at (column
thickness - depth) above the grid base, so on a slope it is that deep only on average; the
column thickness assumes a terrain-following mesh (checked only for the flat DEM). Compare
slopes by pointing `synthetic_ksfire_analysis.py --root` at each slope root.

### Arguments for every run (`--extra`)

Everything after `--extra` is appended to the command of **every** scenario. It goes last,
so it overrides the scenario's own value for the same option. Put `--extra` at the end of
the command line.

```bash
# tighter solver settings for all runs
python scenario.py --only lai_2L control_2L --extra --dtmin 1e-4 --deltat 10

# different fire date and thicker topsoil for all runs
python scenario.py --only layers2 --extra --ks-fire-date 2020-09-01 --top-depth 0.2

# user-defined static Ks for the controls (m/s)
python scenario.py --only control control_2L --extra --permx 1e-5
```

Changing a simulation setting (including through `--extra`) changes the run fingerprint,
which invalidates that scenario's monthly backup (see Restart below).

## Custom scenarios (`--config`)

Add scenarios in a `.toml` (or `.json`) file. Each table is one scenario:

```toml
[lai_2L_lowKs]
args = ["--ks-fire-mode", "lai", "--ks-fire-ll", "1.4e-6", "--ks-fire-ul", "1.2e-5",
        "--two-layers"]
tags = ["custom", "anchors"]
description = "LAI mode, strong post-fire Ks loss"
```

- `args`: list of strings passed to `synthetic_ksfire_LAI.py`; write numbers as strings.
- `tags`: optional, default `["custom"]`.
- `description`: optional, shown by `--list`.
- A name that already exists in `scenario.py` is replaced by the one in the file.
- Ks values are in m/s (1 mm/hr = 2.778e-7 m/s).

The provided `scenarios.toml` contains severity, Ks-anchor, recovery-rate, vegetation-speed
and topsoil-geometry variants:

```bash
python scenario.py --config scenarios.toml --list
python scenario.py --config scenarios.toml --only custom --dry-run
python scenario.py --config scenarios.toml --only custom --skip ic --jobs 3
python scenario.py --config scenarios.toml --only control control_2L sev_low_2L sev_high_2L
```

## Useful options of `synthetic_ksfire_LAI.py` (pass through `args` or `--extra`)

| Option | Meaning |
|---|---|
| `--ks-fire-mode {off,time,lai,burn_severity}` | how Ks is computed each month |
| `--ks-fire-preset {overall_average,median,lower_bound,upper_bound}` | Thomas et al. 2021 field fits (time mode) |
| `--ks-fire-date YYYY-MM-DD` | fire date (default 2020-07-01) |
| `--ks-fire-ll`, `--ks-fire-ul` | burned / unburned Ks anchors (m/s) |
| `--ks-fire-rate` | logistic recovery rate (yr^-1) |
| `--ks-fire-severity`, `--ks-fire-severity-reduction` | burn-severity mode |
| `--lai-ref {ratio,raw}` | LAI fed to the lai mode |
| `--lai-recovery-rate` | speed of the synthetic LAI regrowth (yr^-1) |
| `--two-layers` / `--no-two-layers` | 2-layer or homogeneous soil |
| `--top-depth M` / `--top-layers N` | burned topsoil thickness / number of strata |
| `--permx M/S` | static Ks (control only) |
| `--permx-deep M/S` | Ks of the stable subsoil |
| `--slope M/M` | DEM slope, single outlet (0 = flat) |
| `--ic-head M` / `--ic-wtposition M` | uniform or hydrostatic initial condition |
| `--dtmin`, `--dtmax`, `--deltat` | CATHY solver time-step settings |
| `--max-consecutive-fail N` | stop after N months in a row that fail |
| `--save-every N` | rewrite merged outputs every N simulated months |

Full list: `python synthetic_ksfire_LAI.py --help`.

## Comparison (`_comparison/`)

Built by `synthetic_ksfire_analysis.py` from every run that has an
`outputs/areamean_timeseries.csv`. All quantities are means over the full 20x20 area.

| File | Content |
|---|---|
| `fig_cmp_atmospheric.png` | atmospheric BCs shared by all runs: rain, ETp, cumulative rain vs ETp, LAI |
| `fig_cmp_monthly_eta_recharge.png` | monthly ETa and recharge, all scenarios |
| `fig_cmp_cumulative.png` | cumulative ETa and recharge, plus the cumulative ΔETa to the control (zoomed) |
| `fig_cmp_delta_vs_control.png` | scenario minus matching control (monthly with rain/ETp behind, cumulative) |
| `fig_cmp_periods.png` | half-year totals before/after the fire, and difference to control |
| `fig_cmp_soil_water.png` | Sw and psi at 0, 1, 2 m depth |

Every figure has a LAI panel at the bottom. Half-year labels read e.g. `2020 / Jul–Dec /
after fire` (the keys `2020-H2` are only used in the CSVs: H1 = Jan–Jun, H2 = Jul–Dec).
`fig_cmp_ks_response.png` was dropped: it showed nothing that `fig_cmp_delta_vs_control.png` does not.

**Why can the cumulative curves look identical?** ETa stays close to ETp, so absolute
cumulative ETa of different scenarios differs by far less than its own magnitude; read the
cumulative Δ row instead. If ETa == ETp in every month (dry default forcing, soil never
limiting) the comparison prints a NOTE and the Δ is exactly 0: Ks cannot act. Use `--wet`.
| `summary_periods.csv`, `summary_annual.csv` | totals and differences per period / year |

**Matching control:** for each run the baseline is the `--ks-fire-mode off` run with the
same number of layers (`control` for 1 layer, `control_2L` for 2 layers). The weather is
identical, so scenario minus control isolates the effect of Ks. Always include both
controls in a batch. Do not add further `off` scenarios without excluding them
(`--compare-exclude`), since the first one found per layer setup becomes the baseline.

Only compare runs that share the same period and initial condition. `deep_wt` is excluded
automatically; use `--compare-exclude` for any other run with a different setup.

## Restart and resuming

- Every month that runs cleanly is backed up under `<root>/<name>/_monthly_backup/<scenario>/`
  (`<name>` is the scenario name used by `scenario.py`).
- Re-running a scenario with **identical settings** (same arguments, same period, same
  forcing files) reloads those months and simulates only the missing ones, e.g. after a
  crash or a Ctrl-C.
- Any change to a simulation setting invalidates the backup of that scenario, including
  `--start-year` / `--end-year`: extending the period re-simulates everything.
- `--no-resume` forces a full re-simulation.
- Outputs are written after each simulated month, so an interrupted run keeps all completed months.

## Exit codes and troubleshooting

`scenario.py` returns 0 if all scenarios succeeded, 1 if at least one failed, 3 if some
were only partly simulated. Each scenario's own code is shown in the final table:
`OK`, `PARTIAL` (some months failed), `NO OUTPUT` (not a single month could be simulated),
`FAILED(n)`.

- Always read `<root>/<name>/run.log` first. Its first line is the exact command used.
- `outputs/failures.csv` lists every failure the run tracked.
- **Every month fails from the first one:** this is a model setup problem (initial
  condition, boundary conditions, atmospheric forcing), not a Ks problem. Try
  `--ic-wtposition <other depth>`, `--dtmin 1e-4 --deltat 10`, or `--debug` to stop at the
  first failing month.
- **"Cannot compare: need >= 2 runs":** the comparison needs at least two runs with an
  `areamean_timeseries.csv`. The message lists what each scenario folder is missing.
- **`--root` is relative to the current directory:** default `../synthetic_LAI_Ks`. Use an
  absolute path to avoid surprises.
- **Parallel runs:** with `--jobs N` the console only shows start/finish lines; progress is in
  each `run.log`.
- Check the plumbing first: `python synthetic_ksfire_LAI.py --smoke-test` and
  `python scenario.py --only <name> --dry-run`.

## Typical workflow

```bash
# 1. sanity checks
python synthetic_ksfire_LAI.py --smoke-test
python scenario.py --dry-run

# 2. minimal comparison: both controls and the LAI scenarios
python scenario.py --only control control_2L lai_1L lai_2L

# 3. all Ks modes, 1 and 2 layers
python scenario.py --only default --jobs 3

# 4. sensitivity tests from the TOML file
python scenario.py --config scenarios.toml --only control control_2L severity anchors vegetation

# 5. rebuild the comparison after adding or excluding runs
python scenario.py --config scenarios.toml --compare-only --compare-exclude lai_2L_wetIC
```
