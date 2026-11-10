# Agramon — LAI-driven CATHY for post-fire surface/subsurface hydrology in semi-arid conditions

## Repositories

| Repository | Role |
|---|---|
| [BenjMy/pycathy_wrapper](https://github.com/BenjMy/pycathy_wrapper) | `pyCATHY` Python wrapper used to build, run and post-process CATHY models |
| [CATHY-Org/CATHY_SCF_LAI](https://github.com/CATHY-Org/CATHY_SCF_LAI) | CATHY source version with SCF / LAI (the new version tested in this paper) |
| [CATHY-Org/CATHY_src_WithIrrigation](https://github.com/CATHY-Org/CATHY_src_WithIrrigation) | CATHY source version with irrigation (`withIrr`), used for the rain-via-irrigation-map trick (§4.5) |
| [BenjMy/GRwater_CATmentHYdrology](https://github.com/BenjMy/GRwater_CATmentHYdrology) | Project repository with the catchment hydrology codes (the scripts in this README) |

> **Tip:** `pycathy_wrapper` and `GRwater_CATmentHYdrology` can be **cloned** into the working environment to ease the process and testing (running the pipeline, reproducing bugs, testing the `withIrr` forcing).

> Working README for the peer-reviewed article. It doubles as a **code map** and a **progress tracker**.
> Legend: ✅ done / implemented · 🟡 partly done or needs checking · ⬜ not started · 🐞 known bug

---

## 1. Purpose of the paper

The article has two intertwined goals.

**Methodological.** Test a new version of CATHY (via `pyCATHY`) that is better suited to semi-arid conditions by adding a **time-varying LAI** to the root-water-uptake / ET partitioning.

- Standard CATHY has no LAI. Vegetation cover enters only through the soil cover fraction **SCF**, which is binary:
  - `SCF = 0`: the **ψmin threshold** approach. ETa is limited by `PMIN` and the surface BC switches between atmosphere-controlled and soil-limited. This is appropriate for sparse or arid cover, where bare-soil evaporation matters.
  - `SCF = 1`: the **Feddes** root-water-uptake approach. ET is assumed to be essentially all transpiration.
- Camporese et al. (2014) suggest that `SCF = 0` can be realistic for semi-arid systems too.
- Here, monthly LAI maps drive a vegetation-class map. Each class gets its own root depth (`ZROOT`) and Feddes parameters, so the canopy can change through time (burn, recovery).

**Scientific.** Apply the model to pre- and post-fire surface/subsurface hydrology at Agramon.

- Compare a **burned catchment** with a nearby **intact (control) catchment**.
- Quantify the **rate of hydrological recovery** (LAI → ET → soil moisture).
- Quantify **water redistribution** after the fire: runoff/streamflow vs infiltration vs recharge.
- The fire date is **2020-07-01** (confirmed). Geophysical data (**ERT lines + EM**) were collected after the fire. They are used to infer **bedrock/regolith depth**, which constrains the model domain. Erosion may have changed it between the fire and the survey, and the inversion did not work very well (see §8).

**Post-fire Ks scenarios.** Fire is reported in the literature to change soil hydraulic conductivity (e.g. water repellency / hydrophobicity, crusting, sealing), and this directly controls infiltration vs runoff. Rather than assume one value, the post-fire Ks is varied in a set of scenarios built from the reported literature changes (see §4.2).

**Validation (limited by design).** At Agramon there are no measured ET or discharge data to validate against. The only available check is to compare the **water-balance ETa from CATHY** with **ETa from energy-balance (EB) models (pyTSEB)**. EB ETa is *not* the truth, but it is a separate estimate and useful as a plausibility check (see §4.4). The procedure could potentially also be validated against a **deforested forest site in the Czech Republic**, but this is uncertain: that site is not water-limited like Agramon, so it would test the model in a different regime.

**Additional sensitivity.** Run the model with different ETp inputs (`ET_0_wet` vs `ET_0_dry`, from the SW wet/dry bounds) to show that **using the right ETp input matters** for CATHY.

---

## 2. Repository contents

| File | Role |
|---|---|
| `Agramon_withLAI_withETp.py` | **Main simulation script.** Full pyCATHY pipeline with monthly LAI loop, ETp/rain forcing, BCs, ICs, spin-up. |
| `run_scenarii_withETp_withLAI_ALL.sh` | **Batch runner.** Crosses dem-plot × PMIN × outlet BC × bottom BC × spin-up × IC, forwards the other flags, skips scenarios that already completed, and post-processes. |
| `to_run_paper.md` | Scratch file of run recipes (`for` loops of the batch runner), notes and known bugs. |
| `Agramon_utils.py` | Helpers: plot shapefiles/treatments, DEM loading, EM conductivity sensors and stats, ERT profile extraction, ERA5 point series, notebook plotting helpers. |
| `results_plotter.py` | Shared plotting module (hydro, ET, LAI, recharge, hydrograph, node maps, fire-window and weekly variants, Fig. 8-style stack). |
| `agramon_config.py`, `agramon_forcing.py`, `agramon_lai.py`, `agramon_checks.py`, `agramon_runlog.py`, `agramon_outputs.py` | **New modules** split out of the main script (paths/constants, input loaders, LAI → veg → soil, sanity checks + ETa patch, run logs, running combined NetCDF/CSV outputs). They must sit in the same folder as the main script (`lib/`). See `REFACTORING.md`. |
| `REFACTORING.md` | Tracker of the code restructuring (what changed, how it was checked, what is left). Work happens on the branch `refactor-agramon-modules`. |
| `analysis_ET_compare.py` | Compares **EO energy-balance ETa** (TSEB/biophysical pipeline) vs **water-balance (WB) ETa** from CATHY, plus a comparison against the treatments CSV for `etp_ts` runs. |

Referenced in the code but **not in this upload** (check they are in the repo): `withLAI_results.py` (post-processing driver), `geoutils` (`uCATHY`), `plot_fig8_ETa_ETp_vs_psi.py`.

---

## 3. Simulation pipeline (`Agramon_withLAI_withETp.py`)

> **Code structure (work in progress).** `run_pipeline` is being split into steps (`REFACTORING.md` §2.2). So far: steps 1–3 below are in `setup_inputs(args)` (returns a `RunInputs` dataclass) and the running output files are written by `agramon_outputs.save_combined_outputs`. The behaviour is unchanged; the other steps are still inline in `run_pipeline`.

1. Resolve the DEM source (GeoTIFF, `--dem-tif`) and log the simulation parameters (CSV + markdown log).
2. Load the shapefiles (13 micro-catchments, index = `fid_1`, selected with `--dem-plot 1..13`).
3. Load the forcing and LAI. All grids go to EPSG:25830; LAI/ETp/rain are natively EPSG:32630. The **monthly LAI comes from Earth-Observation biophysical studies** (not from a model), so it is an observed canopy forcing that drives the vegetation classes.
   - Forcing switch `--forcing {era5 | etp | etp_ts}`:
     - `era5`: ERA5 point series.
     - `etp`: satellite ETp raster + satellite rain raster (default).
     - `etp_ts`: ETp from the Agramon ET CSV, one daily value uniform over the catchment (`--etp-ts-column ET_0_wet|ET_0_dry`, `--etp-ts-feature-id`) + satellite rain.
4. Clip, mask and optionally smooth the DEM (`--dem-smooth-sigma`), build the mesh, and run the pyCATHY preprocessor.
5. Map ETp and rain to surface nodes and build the full ATMBC arrays (currently net flux = rain − ETp; a planned alternative that passes rain separately is described in §4.5).
6. Configure soil, ICs and BCs, and write the `etp_output.nc` artifact.
7. Optional spin-up: cycle one year of forcing N times before the run.
8. **Main loop, one CATHY run per calendar month:**
   - monthly LAI → 5 vegetation classes (bare / sparse / moderate / dense / very dense; thresholds `LAI_THRESHOLDS = [1.0, 1.5, 2.0, 2.5]`);
   - class → root depth (`ZROOT = --zroot × factor`);
   - update the soil/veg table (`update_soil_for_active_zones`) with SCF + Feddes parameters;
   - run the processor, then apply the ETa artefact floor and the monthly sanity checks.
9. Aggregate and write NetCDF outputs (`psi_output.nc`, `sw_output.nc`, `et_output.nc`, `etp_output.nc`, …) into `scenario_<n>_<tag>/`.

Scenario folders are self-describing, e.g. `scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c_pmin-m50_ic-wt[_forcing-etpts-wet-fid1]`.

### Key parameters (CLI → meaning)

| Flag | Meaning |
|---|---|
| `--pmin` | ψmin threshold (m) for the atmospheric ↔ soil-limited BC switch |
| `--scf [0–1]` | 0 = PMIN approach, 1 = Feddes (default 1). Only 0 or 1 are meaningful (see §8). |
| `--zroot`, `--zroot-factors B S M D V` | Base root depth and the 5 per-class multipliers |
| `--pcana --pcref --pcwlt --pz --omgc`, `--feddes-crop` | Feddes parameters (`FP_map`) |
| `--permx`, `--permy-ratio`, `--permz-ratio` | Saturated Ks override |
| `--outlet-flux`, `--bottom-flux` | Lateral outlet and bottom Neumann flux (0 = closed, negative = draining) |
| `--ic-mode wt|pressure`, `--ic-wtposition`, `--ic-pressure-head` | Initial condition |
| `--no-spinup`, `--spinup-cycles`, `--spinup-year` | Spin-up |
| `--forcing`, `--etp-ts-column` | ETp source / wet-dry bound |
| `--eta-net-neg-mode {etp,gross,off}` | Handling of the "net = rain − ETp < 0" ETa artefact |
| `--debug` | Extra stdout diagnostics (LAI interpolation NaN checks, per-catchment LAI coverage). Off by default. Not forwarded by the batch runner. |

---

## 4. Experimental design

### 4.1 Catchments
| Role | `--dem-plot` (fid_1) | Status |
|---|---|---|
| Burned catchment | all current `--dem-plot` values (the runs mainly use 3, with 2 and 1 in a few) | ✅ all current catchments are burned |
| Intact / control catchment (nearby) | **to be extracted** (hydrologically delineated, coming soon) | ⬜ not available yet |

The control catchment is not one of the existing `fid_1` plots: it has to be delineated from the DEM (hydrological extraction), and then added as a new DEM/plot input with its own LAI, ETp and rain forcing.

### 4.2 Scenario axes
- **ET partitioning:** `SCF ∈ {0, 1}` × `PMIN ∈ {−15, −50, −150, −1e34/−1e35}` × Feddes `PCANA`
- **Root depth:** baseline `0.1 0.3 1.0 1.5 2.5`; deeper `0.3 0.6 2.0 3.0 5.0` (the one used in most recent runs); uniform `2.5` for all classes
- **Ks (`PERMX`):** `1e-6, 1e-4, 6e-4, 1e-3` m/s
- **Post-fire Ks change (repellency etc.):** scenarios where Ks of the burned catchment is changed relative to the pre-fire / control value, with magnitudes taken from the post-fire literature. The literature values and equations are already coded in `lib/ks_fire_utils.py` of the `GRwater_CATmentHYdrology` repo (see §4.6), but **not yet wired into `Agramon_withLAI_withETp.py`**, where `--permx` is still one constant value for the whole run (see §6).
- **BCs:** outlet and bottom ∈ {closed, free}
- **IC:** water table vs uniform pressure head (`−100 m` in the latest runs)
- **Spin-up:** off vs 3 cycles
- **ETp input:** satellite ETp vs `ET_0_wet` vs `ET_0_dry` (`--forcing etp_ts`)

### 4.3 Recipe currently marked "Good one" in `to_run_paper.md`
```bash
KS_VALUES=(1e-4 6e-4 1e-3)
ZROOT_SETS=("0.3 0.6 2.0 3.0 5.0")
./run_scenarii_withETp_withLAI_ALL.sh \
    --dem-plot 3 --pmin -50 -1e35 \
    --outlet both --bottom both --ic both --ic-pressure-head -100 \
    --permx "${KS_VALUES}" --zroot-factors ${ZROOT_SETS}
```
> Check: `--permx "${KS_VALUES}"` expands an array to its **first element only**, so as written it runs only `1e-4`. Use a `for ks in "${KS_VALUES[@]}"` loop (as in the first blocks of `to_run_paper.md`) to cover all Ks values.

### 4.4 Validation strategy
- **Agramon (Spain):** no direct ET or discharge measurements. The check is **CATHY WB ETa vs EB ETa (pyTSEB)**, implemented in `analysis_ET_compare.py` (monthly, annual, spatial, per-LAI-zone, bias). EB ETa is itself a model product with its own uncertainty, so agreement is read as *consistency*, not proof, and disagreement does not by itself say which is wrong. The wet/dry matching rule is applied (wet run vs wet reference, dry run vs dry reference).
- **LAI input:** EO biophysical LAI drives CATHY and is also an input to the EB estimate, so the two ETa products are not fully independent. This should be acknowledged in the paper.
- **Czech Republic deforested forest (possible):** a candidate site for validating the procedure after canopy removal. ⬜ Not confirmed: the site is not water-limited, so it would mainly test the energy-limited regime, not the semi-arid behaviour (SCF/ψmin) that motivates the study. Data availability and what would be compared are to be defined.

### 4.5 Planned improvement: separate rain from ETp (withIrr wrapper trick)
**Problem.** The current setup passes `netValue = rain − ETp` to CATHY as a single combined atmospheric flux (`HSPATM=1`, `IETO=1`). CATHY cannot tell "net infiltration" from "ET still happening at the same time", so on days or nodes where rain ≥ ETp the reported `ACT. ETRA` collapses to 0. This is currently handled by an after-the-fact **ETa artefact floor** (`apply_eta_artefact_floor`, column `ACT. ETRA_patched`), which only floors ETa at ETp where ETa ≈ 0 and rain ≥ ETp. It is a post-processing patch, not a fix to the model forcing.

**Idea.** Use the **`withIrr` version of the pyCATHY wrapper** and supply the **rain through the irrigation map file** (the irrigation input is added as a water source at the surface), while the atmospheric forcing carries **ETp only**. Rain is then no longer netted against ETp, so CATHY computes ETa against the full ETp every day, and the artefact (and the need for the patch) should disappear.

**To verify when implementing:**
- the units and time-stepping of the irrigation map file vs the ATMBC time axis (rain is currently daily, in m/s per surface node);
- that rain entering as an irrigation-type source behaves like rainfall at the surface (infiltration, ponding, runoff, and the ψmin/Feddes switching) and is not treated differently in the water balance outputs, so that recharge and hydrograph diagnostics stay comparable;
- that the monthly loop (`subset_atmbc_for_month`, one CATHY run per month) works with the irrigation file split per month;
- an A/B test: same scenario with the net-forcing + artefact floor vs the withIrr rain forcing. Compare ETa (against EB ETa, §4.4), recharge and outlet flow, to show how much the artefact patch was biasing the results.

**Impact on the paper:** if it works, the artefact floor can be dropped (or kept only as a comparison) and the methods text becomes simpler. It matters most for the semi-arid ETa/ETp-vs-ψ analysis, where ETa on rainy days is otherwise underestimated.

### 4.6 Post-fire Ks: literature and existing code (`lib/ks_fire_utils.py`)
**Reference (corrected).** Thomas et al. (2021), *JGR: Earth Surface*, 126, e2021JF006091. The repo gives only the journal, volume and article number, so verify the full author list and title before citing. It builds on Ebel & Martin (2017), *Hydrological Processes* 31(21):3682-3696 (logistic Ks recovery with time since fire) and uses a burned:unburned Ks ratio of 0.37 from Ebel & Moody (2020). Atchley et al. (2018), *Vadose Zone Journal* 17(1) is the basis of the burn-severity mode.

**What the module provides** (all values in m/s internally; literature values in mm/hr are converted):

| Mode | Idea | Source |
|---|---|---|
| `time` | Ks recovers as a logistic function of time since fire (`ks_recovery_time_based`) | Thomas et al. 2021 Eq. 9; Ebel & Martin 2017 |
| `lai` | Ks linearly rescaled from the observed monthly LAI between a post-fire minimum and a recovered maximum (`ks_from_lai_linear`) | Thomas et al. 2021 Eq. 11 |
| `burn_severity` | One static Ks reduction scaled by burn severity (dNBR). The dNBR→Ks relation in the code is a **placeholder** and must be replaced with a calibrated one | Atchley et al. 2018 |

Field presets from Thomas et al. 2021 (Table 1, 2016 Fish Fire, loamy sand, mm/hr): `upper_bound`, `median`, `lower_bound`, `overall_average`. These give a ready-made set of Ks recovery scenarios. The satellite-LAI anchors are 16.2 mm/hr (post-fire minimum) and 43.8 mm/hr (unburned, Rosetta loamy sand).

**Caveats for Agramon.** The presets come from a different site and soil (San Gabriel Mountains, loamy sand), so they are literature-based scenarios, not calibrated Agramon values. The `lai` mode is the natural fit here because the pipeline already loads an observed monthly LAI, but the LAI thresholds/anchors (`lai_ll`, `lai_ul`) must come from the Agramon LAI record. Check that the Ks applied is meant for the surface layer only; Atchley et al. (2018) reduced Ks only in the top 2 cm of soil, whereas `--permx` currently overrides Ks for every soil zone.

**Integration (to do).** The module is written to be called once per month inside the loop (`compute_permx_for_month` → `update_soil_for_active_zones(permx=...)`) and has CLI helpers (`add_ks_fire_args`, `ks_fire_config_from_args`). In the uploaded main script and in the repo copy, nothing imports it yet.

---

## 5. Outputs and figures available

- **Hydrology at reference nodes:** ψ, sw, ETa, ETp, rain, LAI stacks, with fire-year, fire-month and weekly zooms (`plot_hydro_xr*`, `plot_et_sw_xr*`, `plot_et_sw_lai_xr*`, `plot_lai_et_sw_psi_ratio_xr`).
- **Water balance:** hydrograph, recharge (bar / xarray), catchment ET.
- **LAI diagnostics:** monthly by year, spatial maps, annual histogram/mean/median.
- **Vegetation maps:** history, yearly, GIFs.
- **ETa evaluation** (`analysis_ET_compare.py`):
  - EO-EB vs WB monthly-by-year, annual, monthly totals, scatter, spatial maps;
  - the same by LAI zone (tercile split on pre-fire LAI), with annual bias per zone;
  - for `etp_ts` runs: WB ETa vs treatments CSV (real ET, `ET_0_dry`, `ET_0_wet`), with RMSE/bias/r. A wet run is compared only with the wet reference, and a dry run only with the dry one.
- **Camporese et al. (2014) Fig. 8 reproduction:** ETa/ETp vs shallow ψ by scenario and location (`plot_fig8_ETa_ETp_vs_psi.py`; requires `etp_output.nc`, so older runs need a rerun).
- The fire date is marked on the figures (`FIRE_DATE = 2020-07-01` in `results_plotter.py`).

---

## 6. Progress tracker

### Model / code
- [x] ✅ Monthly LAI → vegetation-class map → root depth / soil table update loop
- [x] ✅ Fixed set of 5 vegetation classes (no month-to-month renumbering)
- [x] ✅ SCF as a plain CLI parameter, ψmin switch, Feddes parameters (`--pcana/--pcref/--pcwlt/--pz/--omgc`, `--feddes-crop`)
- [x] ✅ Root-depth factors and Ks / anisotropy overrides
- [x] ✅ Outlet and bottom flux BCs (closed / free drainage)
- [x] ✅ IC modes (water table / uniform pressure head) and spin-up cycles
- [x] ✅ Three forcing modes (`era5`, `etp`, `etp_ts` with wet/dry ETp)
- [x] ✅ ETa "net-forcing" artefact floor (`ACT. ETRA_patched`), as a temporary patch
- [ ] ⬜ Replace the net-forcing approach: rain via the irrigation map file with the `withIrr` pyCATHY wrapper, ETp only in ATMBC (§4.5), then A/B test against the patched version
- [x] ✅ Self-describing scenario folders, skip-if-already-done logic, markdown + CSV run logs
- [x] ✅ Batch runner with all axes crossed
- [x] ✅ Comparison against EO-EB ETa and treatments CSV
- [x] ✅ `etp_output.nc` artifact and Camporese Fig. 8 reproduction
- [ ] 🟡 `pyCATHY` wrapper fixes for BCs (`update_nansfdirbc`, `update_nansfneubc`, `update_sfbc`): in progress
- [ ] 🟡 Validate `pyCATHY` fixes for the `map_prop2mesh` bitwise-not bug and the stale `mesh_pv_attributes` cache
- [ ] ⬜ Hydrological recovery metrics (see §7)
- [x] ✅ Post-fire Ks equations and Thomas et al. 2021 presets coded in `lib/ks_fire_utils.py` (`time`, `lai`, `burn_severity` modes)
- [ ] ⬜ Wire `ks_fire_utils` into `Agramon_withLAI_withETp.py` (monthly `PERMX` override via CLI) and define the Ks scenario list (preset × mode)
- [ ] ⬜ Replace the placeholder dNBR→Ks relation if the burn-severity mode is used

### Simulations (to confirm which are completed and checked)
- [ ] 🟡 Burned catchment, baseline: `SCF=1`, `PMIN=−1e34`, closed, `Ks=1e-4`, both root-depth sets
- [ ] 🟡 Burned catchment, `PMIN` sweep (−15, −50, −150, −1e34/−1e35)
- [ ] 🟡 `SCF=0` vs `SCF=1`
- [ ] 🟡 Ks sweep (1e-6 … 1e-3)
- [ ] 🟡 Open vs closed BC combinations, `ic both`, `ic-pressure-head −100`
- [ ] 🟡 Spin-up (3 cycles) vs no spin-up
- [ ] ⬜ **Extract the control catchment hydrologically** (DEM delineation), then run the full selected configuration on it
- [ ] ⬜ **Control catchment simulations**: full repeat of the selected configuration
- [ ] ⬜ **Post-fire Ks scenarios** (burned catchment, several Ks reductions) vs the control
- [ ] ⬜ **ETp sensitivity**: `etp_ts` with `ET_0_wet` vs `ET_0_dry` vs satellite ETp, compared with EO-EB and CSV ET
- [ ] ⬜ Final "paper" configuration frozen and rerun on both catchments

### Validation
- [x] ✅ Tooling for CATHY WB ETa vs EB (pyTSEB) ETa comparison (`analysis_ET_compare.py`)
- [ ] 🟡 Run the comparison for the final burned and control scenarios and summarise skill (RMSE, bias, r)
- [ ] ⬜ Decide whether to include the Czech Republic deforested site (uncertain; not water-limited)

### Geophysics
- [ ] 🟡 ERT / EM inversion for bedrock depth (ran, but it did not work very well)
- [ ] ⬜ Decide how bedrock depth enters the CATHY mesh (flat bottom vs variable depth)
- [ ] ⬜ Discuss post-fire erosion (time between the 2020-07-01 fire and the geophysical survey; survey date to add)

### Paper
- [ ] ⬜ Methods text (LAI → ET partitioning, Feddes/ψmin/SCF)
- [ ] ⬜ Result figures chosen (burned vs control: ETa, sw, ψ, recharge, runoff)
- [ ] ⬜ Discussion: SCF 0 vs 1 in semi-arid conditions, ETp sensitivity
- [ ] ⬜ Write-up

---

## 7. Analyses still to define (for the scientific question)

- Recovery rate: e.g. time for LAI, ETa/ETp and soil moisture in the burned catchment to return to the control catchment's level.
- Water redistribution, burned vs control, by year: ETa, recharge (bottom flux / ΔS), outlet hydrograph and infiltration/runoff partition.
- Effect of root depth (the root-depth sets are the proxy for rooting recovery) and of Ks (post-fire repellency/crusting is explored through the Ks scenarios in §4.2, cf. Thomas et al. 2021, "Ks and LAI").
- Effect of `hmin` (surface-routing threshold; ponding vs saturation, roughness/microtopography) as a possible improvement (see the notes in `to_run_paper.md`).

---

## 8. Known issues and open questions

### Bugs
- 🐞 **`--dem-plot 2` crashes** in `run_pipeline` at `ds_mesh.x.min()`. The error is `ValueError: zero-size array to reduction operation minimum`, i.e. the mesh/DEM is empty for that plot (clip or mask returns no valid cells). Check the footprint, CRS, NaN masking and `--dem-smooth-sigma`.
- 🐞 pyCATHY `map_prop2mesh` uses `~isinstance(...)` instead of `not isinstance(...)`, so per-node arrays are broadcast as scalars. A related stale-`mesh_pv_attributes` cache can cause node-count mismatches.
- 🟡 `n_surf = int(simu.grid3d["nnod"])` used the stale attribute instead of the re-read `grid3d`, giving "SURFACE NODE out of bounds" warnings. The current script reads `grid3d["nnod"]` from the re-read grid, so this looks fixed; confirm with a run (no warning) before closing it.

### Code hygiene (before releasing the code with the paper)
- Code-hygiene items are now tracked in `REFACTORING.md` (duplicate definitions and stale docstrings are already cleaned in the refactored copy).
- ✅ The `lai_to_veg_map` docstring listed the old 0–0.6 LAI thresholds; it now matches the code (`[1.0, 1.5, 2.0, 2.5]`).
- The CLI defaults (`--start-year 2022 --end-year 2022`) give an empty period. The batch runner overrides them with 2016–2026, but make this consistent.
- The default `ZROOT_FACTORS` (`0.05 0.25 1.0 1.5 2.0`) differ from the sets used in the recipes (`0.3 0.6 2.0 3.0 5.0`). Pick one as the documented default.
- Hard-coded machine paths (`Nextcloud/...`) should become a config file or env vars.

### Scientific caveats (from discussion with the CATHY developer, noted in `to_run_paper.md`)
- **ETa never reaching 0 in published plots:** in the earlier paper, simulated ETa was interpolated in time to compare with measured ETa. Don't expect ETa → 0 in raw CATHY output, and compare like with like.
- **SCF:** designed as 0 (PMIN approach) or 1 (Feddes). Intermediate values were never properly tested and may cause numerical errors. Treat `0 < SCF < 1` as experimental and do not use it for the main results.
- ψmin alone cannot split ETa into soil evaporation vs plant transpiration (it is a lumped surface flux).
- The net-forcing formulation (rain − ETp as one flux) makes ETa collapse to 0 on rain ≥ ETp days. This is handled by the artefact floor for now, but it should be documented in the methods. The planned `withIrr` approach (§4.5) would remove the cause rather than patch the output.

### Open questions
- **Timing:** the fire date is **2020-07-01** (confirmed; matches `FIRE_DATE` in the plotting code). The earlier "12 years after the fire" figure for the geophysical survey does not fit with a 2020 fire, so the survey date needs to be stated.
- Is a stationary ψ at the outlet during spin-up physically expected? **Unknown.** What is known: there is no permanent river at the outlet. This favours a design where the outlet is not a perennial water source or sink. Test it by comparing outlet ψ and flow for closed vs free outlet in the spin-up (no stationary ψ should be forced by the BC itself).
- Which BC set is right for open vs closed domain? **Not settled.** What is known: to evaluate **recharge**, start from an initial condition with a **water table** (`--ic wt`), and the water table is thought to be deep. Still to check: the water-table depth value, and whether a closed bottom is consistent with a bedrock-limited domain (the bedrock depth inversion is uncertain).
- Why is pre-fire sw at mid-hillslope/outlet higher than post-fire when only LAI/`ZROOT` change? Less vegetation should leave more water, not less, so this is counterintuitive and needs a diagnosis.
- How to define and delineate the control catchment (§4.1), and whether its soil, geology and aspect are comparable to the burned ones.
- Which Ks scenarios (preset, mode) are representative for Agramon soils? Thomas et al. 2021 values are for a loamy sand at another site.
- Is the Czech Republic deforested site worth including given that it is not water-limited?
- How independent are CATHY ETa and pyTSEB ETa, given the shared EO LAI input?

---

## 9. Quick start

```bash
# single run
python Agramon_withLAI_withETp.py --dem-plot 3 --pmin -50 --scf 1 \
    --zroot-factors 0.3 0.6 2.0 3.0 5.0 --permx 1e-4 \
    --ic-mode pressure --ic-pressure-head -100 --outlet-flux 0 --bottom-flux 0

# ETp sensitivity (wet vs dry bound)
./run_scenarii_withETp_withLAI_ALL.sh --forcing etp_ts --etp-ts-column ET_0_wet --dem-plot 3 ...
./run_scenarii_withETp_withLAI_ALL.sh --forcing etp_ts --etp-ts-column ET_0_dry --dem-plot 3 ...

# ETa evaluation
python analysis_ET_compare.py --plots all --scenario <N>
```

Use `./run_scenarii_withETp_withLAI_ALL.sh --help` for the full flag list. Pass `--force` to rerun scenarios already completed.

---

## 10. References to cite
- Camporese et al. (2014), *Adv. Water Resour.*: ψmin boundary-condition switching, ETa/ETp vs ψ (Fig. 8).
- Thomas et al. (2021), *JGR: Earth Surface* 126, e2021JF006091: post-fire Ks recovery (time- and LAI-based). Verify author list and title before citing.
- Ebel & Martin (2017), *Hydrological Processes* 31(21):3682-3696; Ebel & Moody (2020); Atchley et al. (2018), *Vadose Zone Journal* 17(1).
- CATHY / pyCATHY documentation, plus the Feddes root-water-uptake model.
