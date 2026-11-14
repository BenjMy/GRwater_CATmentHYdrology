# TO-DO.md — Post-fire hydrology with CATHY (Agramon / GRwater)

Science, experimental design and task tracker. Code structure, CLI and refactoring status live in `REFACTORING.md` (single companion file; `README_paper.md` has been merged into these two files).

Legend: ✅ done · 🟡 partly done / needs checking · ⬜ not started · 🐞 known bug · ❓ decision needed

---

## 0. Context and rules for AI sessions

New methodological CATHY improvement and scientific questions from real post-fire affected field sites.

```
root_lib_folder = /home/ben/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/lib/
root_synthetic_Ks_folder = <root_lib_folder>/synthetic_LAI_Ks/
```

**Rules**
- Use geoutils (`uCATHY`, `uET`) and CATHY modules as much as possible. Ask for geoUtils files whenever needed.
- Clone the repo on branch `refactor-agramon-modules` (procedure in `REFACTORING.md` §0).
- Commits / pushes / patches / bundles are allowed, with **BenjMy as author**.
- The task list below can be reordered, renamed or improved.
- Keep track of what has been done (see the log in §9) for future AI requests.
- The repository may be refactored freely (file renames for readability, new folder structure).

**Repositories**

| Repository | Role |
|---|---|
| [BenjMy/pycathy_wrapper](https://github.com/BenjMy/pycathy_wrapper) | `pyCATHY` Python wrapper to build, run and post-process CATHY models |
| [CATHY-Org/CATHY_SCF_LAI](https://github.com/CATHY-Org/CATHY_SCF_LAI) | CATHY source with SCF / LAI (new version tested in the paper) |
| [CATHY-Org/CATHY_src_WithIrrigation](https://github.com/CATHY-Org/CATHY_src_WithIrrigation) | CATHY source with irrigation (`withIrr`), for the rain-via-irrigation-map trick (§5.2) |
| [BenjMy/GRwater_CATmentHYdrology](https://github.com/BenjMy/GRwater_CATmentHYdrology) | Project repo (catchment hydrology codes) |
| [CATHY-Org](https://github.com/CATHY-Org/) | Other CATHY repos of interest |

`pycathy_wrapper` and `GRwater_CATmentHYdrology` can be cloned into the working environment to ease testing (running the pipeline, reproducing bugs, testing the `withIrr` forcing).

---

## 1. Purpose of the paper

Two intertwined goals.

**Methodological.** Test a new CATHY version (via `pyCATHY`) better suited to semi-arid conditions, by adding a **time-varying LAI** to root-water-uptake / ET partitioning.
- Standard CATHY has no LAI. Vegetation enters only through the soil cover fraction **SCF**, which is binary:
  - `SCF = 0`: **ψmin threshold** approach. ETa limited by `PMIN`; the surface BC switches between atmosphere-controlled and soil-limited. Appropriate for sparse/arid cover where bare-soil evaporation matters.
  - `SCF = 1`: **Feddes** root-water-uptake. ET is essentially all transpiration.
- Camporese et al. (2014) suggest `SCF = 0` can be realistic for semi-arid systems too.
- Here, monthly LAI maps (from Earth-Observation biophysical studies, not a model) drive a vegetation-class map. Each class has its own root depth (`ZROOT`) and Feddes parameters, so the canopy can change through time (burn, recovery).

**Scientific.** Apply the model to pre- and post-fire surface/subsurface hydrology at Agramon.
- Compare a **burned catchment** with a nearby **intact (control) catchment**.
- Quantify the **rate of hydrological recovery** (LAI → ET → soil moisture).
- Quantify **water redistribution** after the fire: runoff/streamflow vs infiltration vs recharge.
- Fire date: **2020-07-01** (confirmed; `FIRE_DATE` in `results_plotter.py`). Geophysical data (ERT lines + EM) were collected after the fire and used to infer **bedrock/regolith depth**, which constrains the model domain. Erosion may have changed it between fire and survey, and the inversion did not work very well.

**Post-fire Ks scenarios.** Fire changes soil hydraulic conductivity (water repellency, crusting, sealing), which controls infiltration vs runoff. Rather than assuming one value, post-fire Ks is varied in scenarios built from literature changes (§5.3).

**Additional sensitivity.** Run with different ETp inputs (`ET_0_wet` vs `ET_0_dry`, from the SW wet/dry bounds) to show that the right ETp input matters for CATHY.

---

## 2. Experimental design

### 2.1 Catchments
| Role | `--dem-plot` (fid_1) | Status |
|---|---|---|
| Burned | all current values (runs mainly use 3, a few 2 and 1) | ✅ all current catchments are burned |
| Intact / control (nearby) | **to be extracted** (hydrological delineation) | ⬜ |

The control catchment is not one of the existing `fid_1` plots: delineate it from the DEM, then add it as a new DEM/plot input with its own LAI, ETp and rain forcing.

### 2.2 Scenario axes
- **ET partitioning:** `SCF ∈ {0, 1}` × `PMIN ∈ {−15, −50, −150, −1e34/−1e35}` × Feddes `PCANA`
- **Root depth:** baseline `0.1 0.3 1.0 1.5 2.5`; deeper `0.3 0.6 2.0 3.0 5.0` (most recent runs); uniform `2.5`
- **Ks (`PERMX`):** `1e-6, 1e-4, 6e-4, 1e-3` m/s
- **Post-fire Ks change:** Ks of the burned catchment changed relative to pre-fire/control value, magnitudes from the literature. Coded in `lib/ks_fire_utils.py` (§5.3) but **not yet wired** into the main script, where `--permx` is one constant for the whole run.
- **BCs:** outlet and bottom ∈ {closed, free}
- **IC:** water table vs uniform pressure head (`−100 m` in latest runs)
- **Spin-up:** off vs 3 cycles
- **ETp input:** satellite ETp vs `ET_0_wet` vs `ET_0_dry` (`--forcing etp_ts`)

### 2.3 Recipe marked "Good one" in `to_run_paper.md`
```bash
KS_VALUES=(1e-4 6e-4 1e-3)
ZROOT_SETS=("0.3 0.6 2.0 3.0 5.0")
./run_scenarii_withETp_withLAI_ALL.sh \
    --dem-plot 3 --pmin -50 -1e35 \
    --outlet both --bottom both --ic both --ic-pressure-head -100 \
    --permx "${KS_VALUES}" --zroot-factors ${ZROOT_SETS}
```
> `--permx "${KS_VALUES}"` expands an array to its **first element only** (runs only `1e-4`). Use `for ks in "${KS_VALUES[@]}"` as in the first blocks of `to_run_paper.md`.

### 2.4 Validation strategy
- **Agramon:** no measured ET or discharge. Check is **CATHY water-balance ETa vs energy-balance ETa (pyTSEB)**, implemented in `analysis_ET_compare.py` (monthly, annual, spatial, per-LAI-zone, bias). EB ETa is a model product, so agreement means *consistency*, not proof. Wet/dry matching rule: wet run vs wet reference, dry run vs dry reference.
- **Independence caveat:** EO biophysical LAI drives CATHY and is also an input to the EB estimate, so the two ETa products are not fully independent. Acknowledge in the paper.
- **Czech Republic deforested forest (possible):** candidate for validating the procedure after canopy removal. ⬜ Not confirmed: not water-limited, so it tests the energy-limited regime, not the semi-arid behaviour (SCF/ψmin) that motivates the study.
- **Ultimately** all scenarios should be validated with ETa from an EC tower where one exists (test_validation_tomas).

---

## 3. Task list by study

### 3.1 Synthetic Ks fire study (`lib/synthetic_LAI_Ks/synthetic_ksfire_LAI.py`)
- [x] ✅ Results relocated outside `lib/`: default root `<repo>/results/synthetic_LAI_Ks[_wet][_quick]` (override `$SYNTH_KS_ROOT` / `--root`); committed runs moved there
- [x] ✅ `_comparison`: new `fig_cmp_atmospheric.png` (rain, ETp, cumulative balance, LAI)
- [x] ✅ README: hypotheses table per scenario + daily-forcing caveat (to review scientifically)
- [x] ✅ `fig_cmp_ks_response` dropped
- [x] ✅ `fig_cmp_delta_vs_control`: rain (bars) and ETp (dotted) on twin axes of the monthly-difference subplots
- [x] ✅ H1/H2 relabelled `Jan–Jun` / `Jul–Dec` + before/after fire (CSV keys unchanged)
- [x] ✅ LAI panel in all `_comparison` figures (half-year mean LAI in `fig_cmp_periods`); per-run hydro figure already had it
- [x] 🟡 `--wet` IC now hydrostatic, water table 2.5 m below surface (`--ic-wt-depth`), was uniform −1 m (saturated column). **Rerun the wet batch to check**; mesh is 3 m deep so "2 m" is ~1.7 m (labelled)
- [x] 🟡 Diagnosed: cumulative ETa ≈ cumulative ETp (hundreds of mm) hides Δ of a few mm; added a cumulative-Δ row. In the committed dry runs ETa == ETp every month (Δ = 0), now flagged by a NOTE. Not verified on the wet runs (not in the repo)
- [x] 🟡 `scenario.py --quick` (`--n 5`, root suffix `_quick`); not run end-to-end here (no pyCATHY in the session)
- [ ] New study: different DEM slopes (single outlet)

### 3.2 `test_validation_tomas` (validation site, `lib/test_validation_tomas/`)
State: **waiting for actual ETa data from Hector** (existing data do not cover the right area). Inputs from the EB model are `.nc` (or `.csv` if a homogeneous spatial time series).
- [ ] Quick smoke-test files
- [ ] Report `.md` with hypotheses tested
- [ ] Plots showing validation against EC-tower ETa (if not already existing)

### 3.3 Agramon (real field site, `lib/`)
Test: Ks-vs-LAI scenarios (as in the synthetic case); CATHY code versions (withLAI vs classic, different Fortran sources); ETp inputs (computed from SW, wet/dry). Validate with EB-model ETa (pyTSEB).
- [ ] Improve `to_run_paper.md` (→ becomes a list of scenario config files, see `REFACTORING.md` 2.7)
- [ ] Plots showing validation (if not already existing)
- [ ] New nearby **control** subcatchment unaffected by fire; compare ETa/recharge with the affected ones
- [ ] Quick smoke-test files
- [ ] Report `.md` with hypotheses tested

### 3.4 Input data format (EB model; inspected with `inspect-nc` in `geoUtils/data_tests/agramon`)
| File | Content | Grid |
|---|---|---|
| `agramon_ERA5.nc` | ERA5 hourly (wind, T, dewpoint, precip, pressure, SW/LW down, …), 2015-12 → 2026-01 | 8×8 |
| `agramon_daily-meteo.nc` | daily ETr, Tmin/Tmax, SW-IN, EA, WS, PA, TP, 2016 → 2025 (3653 steps) | 13×13 |
| `agramon_dms.nc`, `agramon_dms-corr.nc` | LST (4000 steps) | 147×126 |
| `agramon_lcd-l2a.nc` | S2 bands B2–B7 + LST (313 steps) | 147×126 |
| `agramon_s2-l2a.nc` | S2 (1391 steps); `inspect-nc` fails on it (`isfinite` TypeError, non-numeric dtype) | 147×126 |
| `agramon_s2-lr.nc`, `agramon_s2-squared.nc` | S2, same TypeError | 33×33 |
| `agramon_s3-meteo.nc` | satellite-overpass TA, EA, WS, PA, AOT, TCWV, PAR/NIR direct+diffuse, LW-IN | 13×13 |
| `agramon_s3.nc` | S3 LST, VZA, SZA | 33×33 |

Note: some variables lack units/names (`TP-DD`, `LW-IN`, `VZA`, `SZA`); the dms LST files report `valid=-1` (check the NaN/fill handling).

---

## 4. Analyses still to define (scientific questions)
- **Recovery rate:** time for LAI, ETa/ETp and soil moisture in the burned catchment to return to the control level.
- **Water redistribution** (burned vs control, by year): ETa, recharge (bottom flux / ΔS), outlet hydrograph, infiltration/runoff partition.
- **Root depth** (proxy for rooting recovery) and **Ks** (post-fire repellency/crusting via the Ks scenarios, cf. Thomas et al. 2021: Ks and LAI).
- **`hmin`** (surface-routing threshold; ponding vs saturation, roughness/microtopography) as a possible improvement (notes in `to_run_paper.md`).

---

## 5. Method design notes

### 5.1 Scientific caveats (CATHY developer discussion, noted in `to_run_paper.md`)
- **ETa never reaching 0** in published plots: in the earlier paper, simulated ETa was interpolated in time to compare with measured ETa. Don't expect ETa → 0 in raw CATHY output; compare like with like.
- **SCF:** designed as 0 (PMIN) or 1 (Feddes). Intermediate values never properly tested and may cause numerical errors; treat `0 < SCF < 1` as experimental, not for main results.
- ψmin alone cannot split ETa into soil evaporation vs plant transpiration (lumped surface flux).
- **Net-forcing formulation** (rain − ETp as one flux) makes ETa collapse to 0 on rain ≥ ETp days. Handled by the artefact floor for now; document in the methods. The `withIrr` approach (§5.2) would remove the cause.

### 5.2 Planned improvement: separate rain from ETp (withIrr trick)
**Problem.** Currently `netValue = rain − ETp` is passed as one atmospheric flux (`HSPATM=1`, `IETO=1`). CATHY cannot tell net infiltration from ET happening at the same time, so where rain ≥ ETp the reported `ACT. ETRA` collapses to 0. Currently patched after the fact by the **ETa artefact floor** (`apply_eta_artefact_floor`, column `ACT. ETRA_patched`), which floors ETa at ETp where ETa ≈ 0 and rain ≥ ETp. It is a post-processing patch, not a fix.

**Idea.** Use the **`withIrr` pyCATHY wrapper**: rain through the **irrigation map file** (surface water source), ATMBC carrying **ETp only**. Rain is no longer netted against ETp, so CATHY computes ETa against the full ETp daily and the artefact (and the patch) should disappear.

**To verify:**
- units and time-stepping of the irrigation map file vs the ATMBC time axis (rain is daily, m/s per surface node);
- that irrigation-type rain behaves like rainfall (infiltration, ponding, runoff, ψmin/Feddes switching) and is treated the same in water-balance outputs (recharge, hydrograph stay comparable);
- that the monthly loop (`subset_atmbc_for_month`) works with the irrigation file split per month;
- **A/B test:** same scenario with net forcing + artefact floor vs withIrr rain forcing; compare ETa (vs EB ETa), recharge and outlet flow to quantify the bias of the patch.

**Impact:** if it works, drop the floor (or keep it only as comparison); simpler methods text. Matters most for the semi-arid ETa/ETp-vs-ψ analysis. Code-side tasks: `REFACTORING.md` 3.2.

### 5.3 Post-fire Ks: literature and existing code (`lib/ks_fire_utils.py`)
**Reference.** Thomas et al. (2021), *JGR: Earth Surface* 126, e2021JF006091 (verify full author list and title before citing). Builds on Ebel & Martin (2017), *Hydrol. Processes* 31(21):3682-3696 (logistic Ks recovery with time since fire); uses burned:unburned Ks ratio 0.37 from Ebel & Moody (2020). Atchley et al. (2018), *Vadose Zone J.* 17(1) is the basis of the burn-severity mode.

| Mode | Idea | Source |
|---|---|---|
| `time` | Ks recovers logistically with time since fire (`ks_recovery_time_based`) | Thomas 2021 Eq. 9; Ebel & Martin 2017 |
| `lai` | Ks linearly rescaled from observed monthly LAI between post-fire minimum and recovered maximum (`ks_from_lai_linear`) | Thomas 2021 Eq. 11 |
| `burn_severity` | One static Ks reduction scaled by dNBR (**placeholder** dNBR→Ks relation, to be replaced by a calibrated one) | Atchley 2018 |

Presets (Thomas 2021 Table 1, 2016 Fish Fire, loamy sand, mm/hr): `upper_bound`, `median`, `lower_bound`, `overall_average`. Satellite-LAI anchors: 16.2 mm/hr (post-fire min), 43.8 mm/hr (unburned, Rosetta loamy sand). Values are m/s internally.

**Caveats for Agramon.** Presets come from another site/soil (San Gabriel Mountains), so they are literature scenarios, not calibrated. The `lai` mode fits best (observed monthly LAI already loaded) but `lai_ll`/`lai_ul` must come from the Agramon LAI record. Atchley et al. reduced Ks only in the top 2 cm, whereas `--permx` overrides Ks for every soil zone: decide whether Ks should apply to the surface layer only.

---

## 6. Paper progress tracker

### Model / science
- [x] ✅ Monthly LAI → vegetation class → root depth / soil table loop; fixed 5 classes
- [x] ✅ SCF, ψmin switch, Feddes parameters, root-depth factors, Ks/anisotropy overrides
- [x] ✅ Outlet/bottom flux BCs, IC modes, spin-up cycles
- [x] ✅ Three forcing modes (`era5` deprecated, `etp`, `etp_ts` wet/dry)
- [x] ✅ ETa net-forcing artefact floor (temporary patch)
- [ ] ⬜ Rain via irrigation file with `withIrr` + A/B test (§5.2)
- [x] ✅ Self-describing scenario folders, skip-if-done, markdown + CSV run logs, batch runner
- [x] ✅ Comparison against EO-EB ETa and treatments CSV; `etp_output.nc` + Camporese Fig. 8
- [x] ✅ Post-fire Ks equations and presets in `ks_fire_utils.py`
- [ ] ⬜ Wire `ks_fire_utils` into main script + define Ks scenario list (preset × mode) (`REFACTORING.md` 3.4)
- [ ] ⬜ Replace placeholder dNBR→Ks relation if burn-severity mode is used
- [ ] 🟡 `pyCATHY` wrapper fixes for BCs (`update_nansfdirbc`, `update_nansfneubc`, `update_sfbc`): in progress
- [ ] 🟡 Validate `pyCATHY` fixes for the `map_prop2mesh` bug and the stale `mesh_pv_attributes` cache
- [ ] ⬜ Hydrological recovery metrics (§4)

### Simulations (confirm which are completed and checked)
- [ ] 🟡 Burned, baseline: `SCF=1`, `PMIN=−1e34`, closed, `Ks=1e-4`, both root-depth sets
- [ ] 🟡 `PMIN` sweep (−15, −50, −150, −1e34/−1e35)
- [ ] 🟡 `SCF=0` vs `SCF=1`
- [ ] 🟡 Ks sweep (1e-6 … 1e-3)
- [ ] 🟡 Open vs closed BCs, `ic both`, `ic-pressure-head −100`
- [ ] 🟡 Spin-up (3 cycles) vs none
- [ ] ⬜ **Extract the control catchment** (DEM delineation), then run the selected configuration on it
- [ ] ⬜ Control catchment simulations
- [ ] ⬜ Post-fire Ks scenarios (burned, several reductions) vs control
- [ ] ⬜ ETp sensitivity: `etp_ts` wet vs dry vs satellite ETp, vs EO-EB and CSV ET
- [ ] ⬜ Freeze the final "paper" configuration and rerun on both catchments

### Validation
- [x] ✅ Tooling for CATHY WB ETa vs EB (pyTSEB) ETa (`analysis_ET_compare.py`)
- [ ] 🟡 Run for final burned and control scenarios; summarise skill (RMSE, bias, r)
- [ ] ⬜ Decide on the Czech Republic deforested site

### Geophysics
- [ ] 🟡 ERT/EM inversion for bedrock depth (ran, did not work very well)
- [ ] ⬜ Decide how bedrock depth enters the mesh (flat vs variable bottom)
- [ ] ⬜ Discuss post-fire erosion (add the survey date; time between 2020-07-01 and survey)

### Paper
- [ ] ⬜ Methods text (LAI → ET partitioning, Feddes/ψmin/SCF)
- [ ] ⬜ Result figures (burned vs control: ETa, sw, ψ, recharge, runoff)
- [ ] ⬜ Discussion: SCF 0 vs 1 in semi-arid conditions, ETp sensitivity
- [ ] ⬜ Write-up

---

## 7. Outputs and figures available
- **Hydrology at reference nodes:** ψ, sw, ETa, ETp, rain, LAI stacks, with fire-year, fire-month and weekly zooms (`plot_hydro_xr*`, `plot_et_sw_xr*`, `plot_et_sw_lai_xr*`, `plot_lai_et_sw_psi_ratio_xr`).
- **Water balance:** hydrograph, recharge (bar / xarray), catchment ET.
- **LAI diagnostics:** monthly by year, spatial maps, annual histogram/mean/median. **Vegetation maps:** history, yearly, GIFs.
- **ETa evaluation** (`analysis_ET_compare.py`): EO-EB vs WB (monthly-by-year, annual, scatter, spatial maps); same by LAI zone (terciles of pre-fire LAI) with annual bias; for `etp_ts` runs, WB ETa vs treatments CSV (real ET, `ET_0_dry`, `ET_0_wet`) with RMSE/bias/r.
- **Camporese et al. (2014) Fig. 8:** ETa/ETp vs shallow ψ by scenario and location (`plot_fig8_ETa_ETp_vs_psi.py`; needs `etp_output.nc`, so older runs need a rerun).
- Fire date marked on figures (`FIRE_DATE = 2020-07-01`).

---

## 8. Open questions
- **Timing:** fire is 2020-07-01. The earlier "12 years after the fire" for the geophysical survey does not fit; state the real survey date.
- **Stationary ψ at the outlet during spin-up:** physically expected? Unknown. No permanent river at the outlet, which favours an outlet that is not a perennial source/sink. Test closed vs free outlet (outlet ψ and flow).
- **BC set (open vs closed domain):** not settled. For **recharge**, start from a water-table IC (`--ic wt`); the water table is thought to be deep. Check the water-table depth value and whether a closed bottom is consistent with a bedrock-limited domain (bedrock inversion uncertain).
- **Counter-intuitive result:** why is pre-fire sw at mid-hillslope/outlet higher than post-fire when only LAI/`ZROOT` change? Less vegetation should leave more water. Needs diagnosis.
- How to delineate the control catchment, and are its soil, geology and aspect comparable?
- Which Ks scenarios (preset, mode) are representative for Agramon soils?
- Is the Czech deforested site worth including?
- How independent are CATHY ETa and pyTSEB ETa given the shared EO LAI?
- ❓ Defaults to document (CLI, shell, `ZROOT_FACTORS`): see `REFACTORING.md` 2.8.

---

## 9. Log of what has been done (for future AI requests)
- Refactoring Phase 1–2 partly done on branch `refactor-agramon-modules`: details and checks in `REFACTORING.md`.
- `README_paper.md` merged into this file and `REFACTORING.md` (2026-11-13); keep only these two files.
- 2026-11-13: §3.1 implemented (see ticks); only "new study: different DEM slopes" remains. Checked: py_compile, `--smoke-test`, comparison rebuilt on committed dry runs. Not checked: full CATHY runs (`--quick`, `--wet`).

## 10. References to cite
- Camporese et al. (2014), *Adv. Water Resour.*: ψmin BC switching, ETa/ETp vs ψ (Fig. 8).
- Thomas et al. (2021), *JGR: Earth Surface* 126, e2021JF006091 (verify authors/title).
- Ebel & Martin (2017), *Hydrol. Processes* 31(21):3682-3696; Ebel & Moody (2020); Atchley et al. (2018), *Vadose Zone J.* 17(1).
- CATHY / pyCATHY documentation; Feddes root-water-uptake model.
