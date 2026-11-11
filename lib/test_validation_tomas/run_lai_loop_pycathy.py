"""
CATHY model looping over LAI (deforestation experiment) - pyCATHY wrapper.

One CATHY project = one scenario. The simulation is split in calendar months (like Agramon_withLAI_withETp.py):
for every month the LAI field is turned into a vegetation-class map (+ root depth ZROOT per class), the atmbc
(rain - ETp, per node, daily) of that month is written, CATHY is run restarting from the last pressure head of the
previous month (update_ic INDP=1), and the compact outputs (ETa, water-table depth, Sw, psi) are harvested.

    python run_lai_loop_pycathy.py --project CS_LAI/obs/ --start 2016-01-01 --end 2024-12-31
    python run_lai_loop_pycathy.py ... --lai-mode climatology --ref-years 2016 2017 2018   # counterfactual: NO LAI drop
    python run_lai_loop_pycathy.py ... --lai-mode observed --lai-scale 0.5                  # sensitivity: halve the LAI
    python run_lai_loop_pycathy.py ... --version withIrr                                    # CATHY version (see scenarios_lai_pycathy.py)
    python run_lai_loop_pycathy.py ... --no-run                                             # INPUT CHECKS ONLY: mesh + LAI cube + forcing, figures, no CATHY run
    python run_lai_loop_pycathy.py ... --months-limit 3                                     # first 3 months only

Data (all in --data-dir, see README.md): CeskeSvycarsko_{s2-l2b,lcd-l2b,SW-ET-gf,TSEB-ET-gf,daily-meteo}.nc,
CeskeSvycarsko_dem.tif, ceske_svycarsko.geojson.

geoutils re-used: uCATHY (GridSpec, polygon_mask, cells_to_nodes, node_index_from_mesh, check_date_overlap,
log_simulation), uET (field_daily_series, norm_units). Helpers specific to this experiment: lai_loop_utils.py.

Checks and plots (lai_loop_plots.py, written to <project>/.../outputs/checks and outputs/figures):
    before the run : check_01..07 (domain, LAI series / acquisitions / maps, rain-ETp-ETr, ETp maps), month_inputs_YYYYMM.png
                     (exactly what is written to CATHY for the first month, --check-every N for more), native_* (pyCATHY
                     show_input / show_bc plots), checks_report.md (OK / INFO / WARN list);
    after the run  : fig_01..06 (forcing vs ETa, water balance, annual budget, layers, LAI response, seasonal), map_01..05
                     (yearly ETa, water-table depth, LAI maps and their change), annual_budget.csv.
    python lai_loop_plots.py --outputs <outputs folder>   redraws the post-run figures without re-running.

NOT TESTED when written: run with --months-limit 2 --dx 100 first (see IMPLEMENTATION_STEPS.md).
"""
import argparse
import hashlib
import json
import pickle
import shutil
import sys
import time
from pathlib import Path

import matplotlib

HERE = Path(__file__).resolve().parent
DATA_DIR = Path("/home/ben/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/lib/test_tomas")

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
g = ap.add_argument_group("project / data")
g.add_argument("--project", default=str(HERE / "CS_LAI" / "obs") + "/", help="CATHY project folder (dirName)")
g.add_argument("--data-dir", default=str(DATA_DIR))
g.add_argument("--geojson", default=None, help="field / catchment polygon (default <data-dir>/ceske_svycarsko.geojson)")
g.add_argument("--lai-file", default="CeskeSvycarsko_s2-l2b.nc", help="LAI NetCDF (or CeskeSvycarsko_lcd-l2b.nc)")
g.add_argument("--lai-var", default="LAI")
g.add_argument("--etp-file", default="CeskeSvycarsko_SW-ET-gf.nc", help="gridded ETp NetCDF (etp-mode product)")
g.add_argument("--etp-var", default="ET_0_wet", help="ET_0_wet | ET_0_dry (| T_0_wet / T_0_dry)")
g.add_argument("--meteo-file", default="CeskeSvycarsko_daily-meteo.nc", help="daily meteo NetCDF: rain (TP-DD) and ETr")
g.add_argument("--rain-var", default="TP-DD")
g.add_argument("--rain-factor", type=float, default=1.0, help="multiplies the rain variable (its units are '?' in the file: "
               "CHECK; 1 = mm/day, 1000 if it is m/day)")
g.add_argument("--etr-var", default="ETr", help="FAO56 reference ET of the meteo file (etp-mode etr_kc)")
g.add_argument("--dem", default="auto", help="GeoTIFF DEM (auto = <data-dir>/CeskeSvycarsko_dem.tif if it exists, none = plane)")
g.add_argument("--crs", default="auto", help="CRS of the CATHY project (auto = CRS of the LAI file, else EPSG:32633)")
g.add_argument("--version", default="1.0.0", help='pyCATHY version: "withIrr" | "1.0.0" | ... (see scenarios_lai_pycathy.py)')
g = ap.add_argument_group("period / mesh")
g.add_argument("--start", default="2016-01-01")
g.add_argument("--end", default="2024-12-31")
g.add_argument("--months-limit", type=int, default=None, help="run only the first N months (smoke test)")
g.add_argument("--dx", type=float, default=50.0, help="cell size (m); the LAI / ET pixels are averaged inside each cell")
g.add_argument("--margin", type=float, default=0.0, help="margin around the polygon bounding box (m)")
g.add_argument("--max-nodes", type=int, default=20000, help="abort if the surface mesh has more nodes (CATHY runtime)")
g.add_argument("--base", type=float, default=None, help="mesh base elevation parameter passed to update_prepo_inputs(base=...) (default template, ~3 m deep)")
g.add_argument("--slope", type=float, default=0.01)
g.add_argument("--outlet", default="sw", choices=["sw", "se", "nw", "ne", "none"])
g.add_argument("--outlet-drop", type=float, default=1e-3)
g.add_argument("--skip-prepro", action="store_true", help="mesh already built in the project")
g = ap.add_argument_group("LAI scenario")
g.add_argument("--lai-mode", default="observed", choices=["observed", "climatology"],
               help="observed = LAI of the EO file; climatology = counterfactual WITHOUT LAI drop (median of --ref-years per calendar month)")
g.add_argument("--ref-years", type=int, nargs="+", default=[2016, 2017, 2018], help="years of the climatology (before the LAI drop)")
g.add_argument("--lai-scale", type=float, default=1.0, help="multiplies the LAI (any mode): sensitivity runs")
g.add_argument("--lai-min-obs", type=int, default=1, help="min acquisitions per month to use it (else carried forward)")
g.add_argument("--lai-thresholds", type=float, nargs="+", default=[0.5, 1.5, 3.0, 4.5],
               help="LAI class limits -> K = len+1 vegetation classes")
g.add_argument("--zroot", type=float, default=1.5, help="base root depth (m)")
g.add_argument("--zroot-factors", type=float, nargs="+", default=[0.1, 0.3, 0.6, 0.9, 1.0],
               help="root-depth factor of each class (len = number of classes); ZROOT = --zroot x factor")
g = ap.add_argument_group("forcing")
g.add_argument("--etp-mode", default="product", choices=["product", "etr_kc"],
               help="product = gridded ETp of --etp-file (already canopy-aware); etr_kc = ETr(meteo) x Kc(LAI) per cell")
g.add_argument("--kc", type=float, default=1.0, help="global factor on ETp")
g.add_argument("--kc-min", type=float, default=0.2, help="etr_kc: Kc at LAI = 0 (bare soil)")
g.add_argument("--kc-max", type=float, default=1.1, help="etr_kc: Kc at full cover")
g.add_argument("--k-ext", type=float, default=0.5, help="etr_kc: extinction coefficient, Kc = kcmin + (kcmax-kcmin)(1-exp(-k LAI))")
g.add_argument("--no-eta-floor", action="store_true", help="do not patch ETa = ETp where rain >= ETp (net flux >= 0 artefact)")
g = ap.add_argument_group("soil / initial & boundary conditions")
g.add_argument("--pmin", type=float, default=-10.0, help="PMIN (m); runner of the Maach case used -1e35 (no lower bound)")
g.add_argument("--scf", type=float, default=1.0)
g.add_argument("--permx", type=float, default=None, help="PERMX (m/s); default pyCATHY SPP default")
g.add_argument("--permy-ratio", type=float, default=1.0)
g.add_argument("--permz-ratio", type=float, default=1.0)
for k in ("pcana", "pcref", "pcwlt", "pz", "omgc"):
    g.add_argument(f"--{k}", type=float, default=None, help=f"Feddes {k.upper()} (default pyCATHY)")
g.add_argument("--wt-height", type=float, default=2.0, help="initial water table height above the mesh base (ic INDP=3)")
g.add_argument("--spinup-year", type=int, default=None, help="year replayed as spin-up (default first year of --start)")
g.add_argument("--spinup-cycles", type=int, default=1, help="number of spin-up cycles of that year (0 = none)")
g.add_argument("--outlet-flux", type=float, default=0.0, help="m/s, negative = lateral outflow on the lowest side (0 = closed box)")
g = ap.add_argument_group("checks / plots")
g.add_argument("--no-checks", action="store_true", help="skip the input checks (figures + checks_report.md) before the run")
g.add_argument("--check-every", type=int, default=0, help="besides the first month, plot the inputs written to CATHY every N months (0 = first only)")
g.add_argument("--no-native", action="store_true", help="skip the pyCATHY native plots (show_input / show_bc / show)")
g.add_argument("--no-plots", action="store_true", help="skip the post-run figures")
g = ap.add_argument_group("run control")
g.add_argument("--dtmin", type=float, default=None)
g.add_argument("--dtmax", type=float, default=None)
g.add_argument("--deltat", type=float, default=None)
g.add_argument("--on-fail", default="stop", choices=["stop", "skip"], help="solver did not reach the end of the month")
g.add_argument("--no-resume", action="store_true", help="ignore / wipe the monthly backup and restart")
g.add_argument("--save-every", type=int, default=12, help="rewrite the merged outputs every N recorded months")
g.add_argument("--no-run", action="store_true")
g.add_argument("--show", action="store_true")
a = ap.parse_args()
if not a.show:
    matplotlib.use("Agg")

import geopandas as gpd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import shapely  # noqa: E402
import xarray as xr  # noqa: E402

for up in HERE.parents:                                    # fallback if geoutils is not pip-installed
    for root in (up / "geoUtils", up):
        if (root / "src" / "geoutils").exists():
            sys.path[:0] = [str(root / "src")]
            break
    else:
        continue
    break
from pyCATHY import cathy_tools  # noqa: E402
from geoutils import uCATHY as eo  # noqa: E402
from geoutils import uET  # noqa: E402
import lai_loop_utils as U  # noqa: E402
import lai_loop_plots as PL  # noqa: E402

K = len(a.lai_thresholds) + 1
if len(a.zroot_factors) != K:
    sys.exit(f"--zroot-factors needs {K} values (one per LAI class), got {len(a.zroot_factors)}")
data = Path(a.data_dir)
prj = Path(a.project)
t_start, t_end = pd.Timestamp(a.start), pd.Timestamp(a.end)
spin_year = a.spinup_year or t_start.year
log = print


def P(name):
    p = Path(name)
    return p if p.is_absolute() or p.exists() else data / name


# =============================================================================================
# 1. field polygon, CRS, domain
# =============================================================================================
gj = Path(a.geojson) if a.geojson else data / "ceske_svycarsko.geojson"
if not gj.exists():
    sys.exit(f"polygon not found: {gj} (--geojson)")
gdf = gpd.read_file(gj)
if gdf.crs is None:
    gdf = gdf.set_crs("EPSG:4326")
lai_path = P(a.lai_file)
with xr.open_dataset(lai_path) as _d:
    _crs = uET._grid_crs(_d)
crs = (_crs.to_string() if _crs is not None else "EPSG:32633") if a.crs == "auto" else a.crs
gdf = gdf.to_crs(crs)
poly = gdf.union_all() if hasattr(gdf, "union_all") else gdf.unary_union
log(f"field   : {gj.name}  area {gdf.area.sum() / 1e6:.2f} km2  crs {crs}")

dem_path = (data / "CeskeSvycarsko_dem.tif" if (data / "CeskeSvycarsko_dem.tif").exists() else None) if a.dem == "auto" \
    else (None if a.dem.lower() == "none" else Path(a.dem))

# =============================================================================================
# 2. CATHY object + mesh (same calls as run_real_case_pycathy.py)
# =============================================================================================
simu = cathy_tools.CATHY(dirName=a.project, prj_name="", version=a.version)
if not a.skip_prepro:
    dem, x_ll, y_ll = U.build_dem(gdf.total_bounds, crs, a.dx, a.margin, a.slope, dem_path)
    if dem.shape[0] * dem.shape[1] > a.max_nodes:
        sys.exit(f"mesh {dem.shape[1]}x{dem.shape[0]} cells = {dem.shape[0] * dem.shape[1]} surface nodes > --max-nodes "
                 f"{a.max_nodes}: increase --dx (now {a.dx} m) or --max-nodes")
    dem = U.add_fake_outlet(dem, a.outlet, a.outlet_drop)
    log(f"DEM     : {dem.shape[1]}x{dem.shape[0]} cells, d={a.dx} m, ll=({x_ll:.1f}, {y_ll:.1f}), z=[{dem.min():.1f}, {dem.max():.1f}] "
        f"({'raster ' + dem_path.name if dem_path else 'synthetic plane'})")
    kw = dict(DEM=dem, delta_x=a.dx, delta_y=a.dx, xllcorner=x_ll, yllcorner=y_ll, dr=a.dx)
    if a.base is not None:
        kw["base"] = a.base
    simu.update_prepo_inputs(**kw)
    hap = Path(simu._prepro_path("hap.in"))          # withIrr template: drop the Strahler line mpar.f90 does not read
    mpar = hap.parent / "src" / "mpar.f90"
    if hap.exists() and mpar.exists() and "strahler" not in mpar.read_text(errors="ignore").lower():
        lines = hap.read_text().splitlines(keepends=True)
        keep = [ln for ln in lines if "strahler" not in ln.lower()]
        if len(keep) != len(lines):
            hap.write_text("".join(keep))
            log("hap.in  : removed the 'Strahler Stream Order' line (not read by the prepro mpar.f90)")
    simu.update_ic(INDP=3, IPOND=0, WTPOSITION=a.wt_height)
    log(f"ic      : INDP=3, water table {a.wt_height:g} m above the mesh base")
    simu.run_preprocessor(verbose=False)
simu.update_veg_map()
simu.update_soil(PMIN=a.pmin)
simu.run_processor(IPRT1=3, verbose=True)            # mesh + grid3d only

grid = eo.GridSpec.from_hapin(simu.hapin, crs=crs)
g3 = getattr(simu, "grid3d", None) or {}
if "mesh3d_nodes" not in g3:
    g3 = simu.read_outputs("grid3d")
mesh = np.asarray(g3["mesh3d_nodes"])
nnod = int(g3["nnod"])
nstr = len(mesh) // nnod
z3 = mesh[:, 2].reshape(nstr, nnod)
if z3[0].mean() < z3[-1].mean():
    sys.exit("unexpected node order: layer 0 is not the surface (z[0] < z[-1]); adapt the reshape of psi / sw")
node_idx = eo.node_index_from_mesh(grid, mesh[:nnod, :2])         # CATHY node -> flat (M+1)x(N+1) node grid
nx_, ny_ = grid.node_coords()
node_x, node_y = nx_.ravel()[node_idx], ny_.ravel()[node_idx]
node_mask = np.asarray(shapely.contains_xy(poly, node_x, node_y))
cell_mask = eo.polygon_mask(gdf, grid)
if not node_mask.any() or not cell_mask.any():
    sys.exit("the polygon contains no mesh node / cell (check --dx, --crs)")
bbox = (grid.xll, grid.yll, grid.xll + grid.N * grid.dx, grid.yll + grid.M * grid.dy)
xc, yc = grid.cell_centers()
chk_dir = Path(simu.workdir) / simu.project_name / "outputs" / "checks"
chk_dir.mkdir(parents=True, exist_ok=True)
rep_ = PL.Report(log)
log(f"mesh    : {grid.M}x{grid.N} cells, {nnod} surface nodes x {nstr} layers; polygon = {int(cell_mask.sum())} cells, "
    f"{int(node_mask.sum())} nodes; mesh depth {z3[0].mean() - z3[-1].mean():.2f} m")

if not a.no_checks:
    log("\nchecks  : domain")
    PL.plot_domain(node_x, node_y, z3[0], float(z3[-1].mean()), node_mask, gdf, chk_dir, rep_, a.dx)

# =============================================================================================
# 3. EO inputs on the mesh: LAI cube, ETp, rain
# =============================================================================================
all_months = U.month_starts(min(t_start, pd.Timestamp(spin_year, 1, 1)), t_end)
lai_src = U.GridSource(lai_path, a.lai_var, bbox, crs)
agg_lai = U.CellAggregator(lai_src.x, lai_src.y, lai_src.crs, grid)
lai_stats = {}
lai_obs, n_obs = U.build_lai_cube(lai_src, agg_lai, all_months, min_obs=a.lai_min_obs, log=log, stats=lai_stats)
lai_src.close()
lai_cube = U.lai_climatology(lai_obs, all_months, a.ref_years) if a.lai_mode == "climatology" else lai_obs
lai_cube = np.clip(lai_cube * a.lai_scale, 0.0, None)
log(f"LAI     : mode {a.lai_mode}{' (ref years ' + str(a.ref_years) + ')' if a.lai_mode == 'climatology' else ''}, "
    f"scale {a.lai_scale:g}; polygon-mean LAI {np.nanmean(lai_cube[:, cell_mask]):.2f} "
    f"[{np.nanmin(lai_cube[:, cell_mask].mean(axis=1)):.2f} - {np.nanmax(lai_cube[:, cell_mask].mean(axis=1)):.2f}]")

meteo = P(a.meteo_file)
rain_s, rain_info = uET.field_daily_series(meteo, gdf, var=a.rain_var, role="rain")
rain_s = rain_s * a.rain_factor
yr_rain = rain_s.groupby(rain_s.index.year).sum()
log(f"rain    : {a.rain_var} [{rain_info['units']}] {rain_info['selection']}; annual totals mm: "
    + ", ".join(f"{y}: {v:.0f}" for y, v in yr_rain.items()) + "   <- CHECK units (--rain-factor)")
etr_s = etp_fb = etp_src = agg_etp = None
if a.etp_mode == "etr_kc":
    etr_s, _ = uET.field_daily_series(meteo, gdf, var=a.etr_var, role="etp")
    etp_series = etr_s
else:
    etp_file = P(a.etp_file)
    etp_src = U.GridSource(etp_file, a.etp_var, bbox, crs)
    agg_etp = U.CellAggregator(etp_src.x, etp_src.y, etp_src.crs, grid)
    etp_fb, _ = uET.field_daily_series(etp_file, gdf, var=a.etp_var, role="etp")        # fallback / report series
    etp_series = etp_fb
log(f"ETp     : mode {a.etp_mode}, {a.etp_var if a.etp_mode == 'product' else a.etr_var}, global factor {a.kc:g}; "
    f"mean {etp_series.mean():.2f} mm/day over the field")
try:
    eo.check_date_overlap((t_start, t_end), {"rain": rain_s, "ETp": etp_series}, strict=True, log=log)
except eo.DateOverlapError as e_:
    sys.exit(f"STOPPED: {e_}")


def etp_cells_for(month, lai_c):
    """(days, ETp cells (nd, M, N) in mm/day) of the calendar month."""
    me = month + pd.offsets.MonthBegin(1)
    days = pd.date_range(month, me - pd.Timedelta(days=1), freq="D")
    nd, M, N = len(days), grid.M, grid.N
    if a.etp_mode == "etr_kc":
        et = etr_s.reindex(days).interpolate(limit_area="inside").bfill().ffill().fillna(0.0).to_numpy()
        kc = U.kc_from_lai(lai_c, a.kc_min, a.kc_max, a.k_ext)
        return days, et[:, None, None] * kc[None] * a.kc
    times, arr = etp_src.read(month, me)
    cells = np.full((nd, M, N), np.nan)
    if len(times):
        cc = agg_etp(arr * etp_src.to_mm_day)
        pos = ((pd.DatetimeIndex(times).normalize() - month).days).to_numpy()
        ok = (pos >= 0) & (pos < nd)
        cells[pos[ok]] = cc[ok]
    fb = etp_fb.reindex(days).interpolate(limit_area="inside").bfill().ffill().fillna(0.0).to_numpy()
    for d in range(nd):                                    # gaps: cell NaN <- field mean of the day <- field series
        c = cells[d]
        if not np.isfinite(c).any():
            cells[d] = fb[d]
        elif not np.isfinite(c).all():
            cells[d] = U.fill_nearest(c, fb[d])
    return days, cells * a.kc


# ---- input checks (figures + report) ---------------------------------------------------------
months_rec = list(U.month_starts(t_start, t_end))
if a.months_limit:
    months_rec = months_rec[: a.months_limit]
if not a.no_checks:
    log("\nchecks  : LAI")
    PL.plot_lai_inputs(all_months, lai_obs, lai_cube, n_obs, cell_mask, xc, yc, gdf, chk_dir, rep_, lai_stats,
                       a.lai_thresholds, a.lai_mode, a.ref_years, a.lai_min_obs, a.lai_scale)
    log("checks  : ETp / rain (reads the ETp of every month once; --no-checks to skip)")
    etp_maps, etp_daily = [], []
    for m_ in months_rec:
        d_, c_ = etp_cells_for(m_, lai_cube[int(np.where(all_months == m_)[0][0])])
        etp_maps.append(np.nanmean(c_, axis=0))
        etp_daily.append(pd.Series(c_[:, cell_mask].mean(axis=1), index=d_))
    etp_used = pd.concat(etp_daily)
    etr_chk = etr_s
    if etr_chk is None:
        try:
            etr_chk, _ = uET.field_daily_series(meteo, gdf, var=a.etr_var, role="etp")
        except Exception as e_:                            # noqa: BLE001
            rep_.info(f"ETr not read for the comparison ({type(e_).__name__}: {e_})")
    lai_poly_m = pd.Series(lai_cube[:, cell_mask].mean(axis=1), index=all_months)
    PL.plot_forcing_inputs(rain_s, etp_used, etr_chk, lai_poly_m, a.etp_mode,
                           dict(kc_min=a.kc_min, kc_max=a.kc_max, k_ext=a.k_ext), a.lai_thresholds, chk_dir, rep_,
                           window=(months_rec[0], months_rec[-1] + pd.offsets.MonthEnd(1)), rain_factor=a.rain_factor)
    PL.plot_etp_maps(pd.DatetimeIndex(months_rec), np.stack(etp_maps), xc, yc, gdf, chk_dir, rep_)
    rep_.write(chk_dir / "checks_report.md", "Input checks (before the run)")


# =============================================================================================
# 4. soil table (ZROOT per vegetation class) and boundary conditions
# =============================================================================================
FEDDES = dict(PCANA=a.pcana, PCREF=a.pcref, PCWLT=a.pcwlt, PZ=a.pz, OMGC=a.omgc)


def update_soil_zones():
    """SOIL table for the K LAI classes (fixed row count = MAXVEG); ZROOT = base x class factor (Agramon logic)."""
    df_fp = simu.set_SOIL_defaults(FP_map_default=True, nveg=K)
    df_spp = simu.set_SOIL_defaults(SPP_map_default=True)
    df_fp = pd.concat([df_fp.iloc[[0]]] * K, ignore_index=True)
    df_fp.index = pd.RangeIndex(1, K + 1)
    df_fp.index.name = "Veg nb"
    df_fp["ZROOT"] = [a.zroot * f for f in a.zroot_factors]
    for col, val in FEDDES.items():
        if val is not None and col in df_fp.columns:
            df_fp[col] = val
    if a.permx is not None:
        for col, val in (("PERMX", a.permx), ("PERMY", a.permx * a.permy_ratio), ("PERMZ", a.permx * a.permz_ratio)):
            if col in df_spp.columns:
                df_spp[col] = val
    simu.update_soil(FP_map=df_fp, SPP_map=df_spp, PMIN=a.pmin, SCF=a.scf)


def set_boundaries(t):
    """closed box (no Dirichlet / seepage face), optional Neumann outflow on the lowest side (Agramon logic)."""
    t = [float(x) for x in t]
    if a.outlet_flux == 0.0:
        simu.update_nansfdirbc(time=t, no_flow=True)
        simu.update_nansfneubc(time=t, no_flow=True)
        simu.update_sfbc(time=t, no_flow=True)
        return
    simu.create_mesh_bounds_df("nansfneubc", mesh, t)
    side = simu.get_outlet_side(outlet_side=None)
    simu.get_outlet_node()
    simu.update_nansfdirbc(time=t, no_flow=True)
    simu.update_sfbc(time=t, no_flow=True)
    areas = simu.compute_nodal_face_areas(side)
    ids = np.array(sorted(areas), dtype=int)
    simu.update_nansfneubc(time=t, nodes=ids, flux=np.array([a.outlet_flux * areas[n] for n in ids]))


# =============================================================================================
# 5. restart / resume bookkeeping
# =============================================================================================
out_dir = Path(simu.workdir) / simu.project_name / "outputs"
state_dir = out_dir / "_monthly_backup"
out_dir.mkdir(parents=True, exist_ok=True)
SKIP_KEYS = {"show", "no_run", "no_checks", "check_every", "no_native", "no_plots", "no_resume", "months_limit", "save_every", "end", "skip_prepro", "on_fail"}
cfg = {k: v for k, v in vars(a).items() if k not in SKIP_KEYS}
cfg_hash = hashlib.md5(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()
if a.no_resume and state_dir.exists():
    shutil.rmtree(state_dir)
if state_dir.exists() and (state_dir / "config.md5").exists() and (state_dir / "config.md5").read_text() != cfg_hash:
    log("resume  : settings changed since the backup -> backup wiped (use the same options to resume)")
    shutil.rmtree(state_dir)
state_dir.mkdir(parents=True, exist_ok=True)
(state_dir / "config.md5").write_text(cfg_hash)
done = {p.stem.split("_")[1] for p in state_dir.glob("rec_*.pkl")}
psi_last = np.load(state_dir / "psi_last.npy") if (state_dir / "psi_last.npy").exists() and done else None
if done:
    log(f"resume  : {len(done)} month(s) already simulated and reloaded from {state_dir}")

params = dict(version=a.version, lai_mode=a.lai_mode, lai_scale=a.lai_scale, etp_mode=a.etp_mode, dx=a.dx, zroot=a.zroot,
              pmin=a.pmin, scf=a.scf, start=a.start, end=a.end, wt_height=a.wt_height, nnod=nnod, nstr=nstr)
try:
    sim_index = eo.log_simulation(prj / "simulation_log.csv", params)
except Exception as e_:                                    # noqa: BLE001
    sim_index = -1
    log("simulation log:", e_)

# =============================================================================================
# 6. one month
# =============================================================================================
mask_cells_all = np.ones((grid.M, grid.N), bool)
lay_ids = sorted({1, nstr // 2, nstr - 1})
tree = None


def node_index_of_et(et):
    """0-based surface node of every row of the ET table."""
    global tree
    if "SURFACE NODE" in et.columns:
        return et["SURFACE NODE"].to_numpy(int) - 1
    from scipy.spatial import cKDTree
    tree = tree or cKDTree(mesh[:nnod, :2])
    dist, idx = tree.query(et[["X", "Y"]].to_numpy(float))
    if dist.max() > 1e-3 * a.dx:
        raise ValueError("ET table X/Y do not match the mesh surface nodes")
    return idx


def prepare_month(month, lai_c):
    """Write every CATHY input of one calendar month (veg map, soil / ZROOT, atmbc, boundaries, TIMPRT); returns what was sent."""
    veg = U.lai_to_veg_map(lai_c, a.lai_thresholds)
    simu.update_veg_map(veg, maxveg=K)
    update_soil_zones()

    days, etp_cells = etp_cells_for(month, lai_c)
    nd = len(days)
    rain_mm = rain_s.reindex(days).fillna(0.0).clip(lower=0.0).to_numpy()
    etp_n = U.cells_to_node_flat(etp_cells)[:, node_idx] / 1e3 / 86400.0          # (nd, nnod) m/s, CATHY node order
    net_n = rain_mm[:, None] / 1e3 / 86400.0 - etp_n                              # rain - ETp
    t = np.arange(nd + 1) * 86400.0
    simu.update_atmbc(HSPATM=0, IETO=1, time=[float(x) for x in t],
                      netValue=list(np.vstack([net_n, net_n[-1:]])), show=False)
    set_boundaries(t)
    simu.update_parm(TIMPRTi=[int(x) for x in t], IPRT=4, VTKF=0)                # after update_atmbc (resets TIMPRT)
    return dict(veg=veg, days=days, etp_cells=etp_cells, rain_mm=rain_mm, etp_n=etp_n, net_n=net_n, t=t)


def check_month(month, tag, lai_c, pm):
    """Figures of the inputs just written to CATHY (own plots + pyCATHY native plots) and the cells -> nodes mass check."""
    PL.plot_month_inputs(month, tag, lai_c, pm["veg"], pm["etp_cells"], pm["rain_mm"], pm["net_n"], node_x, node_y, node_mask,
                         cell_mask, xc, yc, gdf, a.lai_thresholds, [a.zroot * f for f in a.zroot_factors], chk_dir, rep_)
    if not a.no_native:
        PL.native_checks(simu, chk_dir, f"{month:%Y%m}", rep_)
    rep_.write(chk_dir / "checks_report.md", "Input checks")


def run_month(month, lai_c, record, tag):
    """Run CATHY for one calendar month; returns the record dict (None for spin-up / failed months)."""
    global psi_last
    t0 = time.time()
    pm = prepare_month(month, lai_c)
    veg, days, etp_cells, rain_mm, etp_n, net_n, t = (pm[k] for k in ("veg", "days", "etp_cells", "rain_mm", "etp_n", "net_n", "t"))
    nd = len(days)
    if record and not a.no_checks and month in chk_months:
        check_month(month, tag, lai_c, pm)
    if psi_last is not None:
        simu.update_ic(INDP=1, pressure_head_ini=psi_last)
    kw = {k: v for k, v in dict(DTMIN=a.dtmin, DTMAX=a.dtmax, DELTAT=a.deltat).items() if v is not None}
    simu.run_processor(IPRT1=2, TRAFLAG=0, verbose=False, **kw)

    psi_df = simu.read_outputs("psi")
    psi_df = psi_df[0] if isinstance(psi_df, tuple) else psi_df
    if float(psi_df.index[-1]) < t[-1] - 1.0:
        msg = f"{month:%Y-%m}: solver stopped at {float(psi_df.index[-1]):.0f} s < {t[-1]:.0f} s"
        if a.on_fail == "stop":
            raise RuntimeError(msg + " (--on-fail skip keeps the last state and continues)")
        log("  WARNING", msg, "-> month skipped, state kept")
        return None
    Pm = psi_df.to_numpy(float)
    if Pm.shape[1] != nstr * nnod:
        raise ValueError(f"psi has {Pm.shape[1]} columns, expected nstr*nnod = {nstr * nnod}")
    psi_last = Pm[-1].copy()
    np.save(state_dir / "psi_last.npy", psi_last)
    if not record:
        log(f"  [{tag}] {month:%Y-%m} spin-up done in {time.time() - t0:.0f} s")
        return None

    # ---- harvest ---------------------------------------------------------------------------
    ts_psi = np.asarray(psi_df.index, float)
    keep = ts_psi > 0
    dates_state = U.to_datetime(ts_psi[keep], month).normalize()
    P3 = Pm[keep].reshape(-1, nstr, nnod)
    gwl_depth = U.water_table_depth(P3, z3)                                          # (nt, nnod) m below surface
    sw_df, _ = simu.read_outputs("sw")
    S3 = sw_df.to_numpy(float)[np.asarray(sw_df.index, float) > 0].reshape(-1, nstr, nnod)
    state = pd.DataFrame({"date": dates_state, "gwl_mean_m": -gwl_depth[:, node_mask].mean(axis=1),
                          "gwl_p10_m": -np.percentile(gwl_depth[:, node_mask], 90, axis=1),
                          "gwl_p90_m": -np.percentile(gwl_depth[:, node_mask], 10, axis=1)})
    for L in lay_ids:
        state[f"sw_L{L}"] = S3[: len(state), L][:, node_mask].mean(axis=1)
        state[f"psi_L{L}"] = P3[:, L][:, node_mask].mean(axis=1)

    et = simu.read_outputs("ET").copy()
    ts = et["time_sec"]
    ts = ts.dt.total_seconds() if pd.api.types.is_timedelta64_dtype(ts) else ts.astype(float)
    et = et[ts.to_numpy() > 0].copy()
    ts = ts[ts.to_numpy() > 0].to_numpy()
    nid = node_index_of_et(et)
    eta_raw = np.abs(et["ACT. ETRA"].to_numpy(float))
    dd = np.clip(np.floor((ts - 1.0) / 86400.0).astype(int), 0, nd - 1)
    eta_fix = eta_raw if a.no_eta_floor else np.where(net_n[dd, nid] >= 0.0, etp_n[dd, nid], eta_raw)
    inm = node_mask[nid]
    et_df = pd.DataFrame({"date": month + pd.to_timedelta(dd, unit="D"), "eta": eta_fix * 1e3 * 86400.0,
                          "eta_raw": eta_raw * 1e3 * 86400.0, "nid": nid, "in": inm})
    eta_field = et_df[et_df["in"]].groupby("date")[["eta", "eta_raw"]].mean().rename(
        columns={"eta": "eta_mm_d", "eta_raw": "eta_raw_mm_d"})
    eta_map = np.full(nnod, np.nan)
    gm = et_df.groupby("nid")["eta"].mean()
    eta_map[gm.index.to_numpy()] = gm.to_numpy()

    forc = pd.DataFrame({"etp_mm_d": etp_cells[:, cell_mask].mean(axis=1), "rain_mm_d": rain_mm}, index=days)
    forc.index.name = "date"
    veg_frac = np.bincount(veg[cell_mask], minlength=K + 1)[1:] / cell_mask.sum()
    rec = dict(month=month, state=state, eta=eta_field, forcing=forc, gwl_nodes=gwl_depth.astype(np.float32),
               gwl_dates=dates_state, eta_map=eta_map.astype(np.float32), lai=lai_c.astype(np.float32),
               lai_mean=float(lai_c[cell_mask].mean()), veg_frac=veg_frac)
    dm = eta_field["eta_mm_d"].mean()
    log(f"  [{tag}] {month:%Y-%m}  LAI {rec['lai_mean']:.2f}  veg {np.round(veg_frac, 2).tolist()}  "
        f"rain {rain_mm.sum():.0f} mm  ETp {forc['etp_mm_d'].sum():.0f} mm  ETa {eta_field['eta_mm_d'].sum():.0f} mm "
        f"({dm:.2f} mm/d)  GWL {state['gwl_mean_m'].iloc[-1]:+.2f} m  [{time.time() - t0:.0f} s]")
    return rec


# =============================================================================================
# 7. outputs
# =============================================================================================
def write_outputs(recs):
    recs = sorted(recs, key=lambda r: r["month"])
    st = pd.concat([r["state"] for r in recs]).groupby("date").last()
    eta = pd.concat([r["eta"] for r in recs]).groupby("date").mean()
    fo = pd.concat([r["forcing"] for r in recs]).groupby("date").mean()
    lai_m = pd.Series({r["month"]: r["lai_mean"] for r in recs})
    ts_ = st.join(eta, how="outer").join(fo, how="outer")
    ts_["lai_mean"] = lai_m.reindex(ts_.index.to_period("M").to_timestamp()).to_numpy()
    ts_.index.name = "date"
    ts_.to_csv(out_dir / "timeseries_daily.csv")
    pd.DataFrame({"month": lai_m.index, "lai_mean": lai_m.values,
                  **{f"veg_class{i + 1}_frac": [r["veg_frac"][i] for r in recs] for i in range(K)}}).to_csv(
        out_dir / "lai_monthly.csv", index=False)
    dts = np.concatenate([np.asarray(r["gwl_dates"]) for r in recs])
    G = np.concatenate([r["gwl_nodes"] for r in recs])
    _, last = np.unique(dts[::-1], return_index=True)
    sel = np.sort(len(dts) - 1 - last)
    ds_g = xr.Dataset({"gwl_depth_m": (("date", "node"), G[sel])},
                      coords=dict(date=dts[sel], node=np.arange(nnod), x=("node", node_x), y=("node", node_y),
                                  in_polygon=("node", node_mask)))
    ds_g["gwl_depth_m"].attrs.update(long_name="water table depth below the surface (>=0); piezometers = -depth")
    ds_g.to_netcdf(out_dir / "gwl_nodes.nc")
    xr.Dataset({"eta_mm_d": (("month", "node"), np.stack([r["eta_map"] for r in recs]))},
               coords=dict(month=[r["month"] for r in recs], node=np.arange(nnod), x=("node", node_x), y=("node", node_y))
               ).to_netcdf(out_dir / "eta_nodes_monthly.nc")
    xc, yc = grid.cell_centers()
    xr.Dataset({"LAI": (("month", "row", "col"), np.stack([r["lai"] for r in recs]))},
               coords=dict(month=[r["month"] for r in recs], x=(("row", "col"), xc), y=(("row", "col"), yc),
                           in_polygon=(("row", "col"), cell_mask))).to_netcdf(out_dir / "lai_cells_monthly.nc")
    return ts_


# =============================================================================================
# 8. loop
# =============================================================================================
chk_months = {months_rec[0]} | (set(months_rec[:: a.check_every]) if a.check_every > 0 else set())
tasks = []
if a.spinup_cycles > 0 and not done:
    ms = [m for m in all_months if m.year == spin_year]
    tasks += [(m, False, f"spin-up {c + 1}") for c in range(a.spinup_cycles) for m in ms]
tasks += [(m, True, "run") for m in months_rec if f"{m:%Y%m}" not in done]
log(f"\nloop    : {len(tasks)} CATHY month(s) to run ({sum(1 for t_ in tasks if t_[1])} recorded), version {a.version}, "
    f"LAI mode {a.lai_mode}, project {prj}")

cfg_out = dict(cfg, crs=crs, geojson=str(gj), grid=dict(M=grid.M, N=grid.N, dx=grid.dx, xll=grid.xll, yll=grid.yll),
               nnod=nnod, nstr=nstr, layer_depth_m=(z3[0].mean() - z3[:, :].mean(axis=1)).round(3).tolist(), lay_ids=lay_ids,
               sim_index=int(sim_index), outputs=str(out_dir))
(out_dir / "run_config.json").write_text(json.dumps(cfg_out, indent=2, default=str))

if a.no_run:
    if not a.no_checks:                                   # write + plot the inputs of the first month, no CATHY run
        m0 = months_rec[0]
        check_month(m0, "no-run", lai_cube[int(np.where(all_months == m0)[0][0])],
                    prepare_month(m0, lai_cube[int(np.where(all_months == m0)[0][0])]))
        log(f"--no-run: inputs checked, figures in {chk_dir}; {rep_.n_warn} warning(s) (checks_report.md); nothing simulated")
    else:
        log("--no-run: inputs prepared, nothing simulated")
    sys.exit(0)

n_rec = 0
for month, record, tag in tasks:
    i_cube = int(np.where(all_months == month)[0][0])
    rec = run_month(month, lai_cube[i_cube], record, tag)
    if rec is not None:
        pickle.dump(rec, open(state_dir / f"rec_{month:%Y%m}.pkl", "wb"))
        done.add(f"{month:%Y%m}")
        n_rec += 1
        if a.save_every and n_rec % a.save_every == 0:
            recs = [pickle.load(open(p, "rb")) for p in state_dir.glob("rec_*.pkl")]
            write_outputs(recs)

recs = [pickle.load(open(p, "rb")) for p in sorted(state_dir.glob("rec_*.pkl"))]
if not recs:
    sys.exit("no month recorded")
ts_all = write_outputs(recs)
log(f"\noutputs : {out_dir}  (timeseries_daily.csv, gwl_nodes.nc, eta_nodes_monthly.nc, lai_cells_monthly.nc, lai_monthly.csv)")

fig, ax = plt.subplots(4, 1, figsize=(11, 10), sharex=True, layout="constrained")
ax[0].plot(ts_all.index, ts_all["lai_mean"], c="tab:green"); ax[0].set_ylabel("LAI (polygon mean)")
ax[1].plot(ts_all.index, ts_all["rain_mm_d"], c="tab:blue", lw=0.5, label="rain")
ax[1].plot(ts_all.index, ts_all["etp_mm_d"], c="gray", lw=0.8, label="ETp"); ax[1].set_ylabel("mm/day"); ax[1].legend()
ax[2].plot(ts_all.index, ts_all["eta_mm_d"], c="tab:red", lw=0.8); ax[2].set_ylabel("CATHY ETa (mm/day)")
ax[3].fill_between(ts_all.index, ts_all["gwl_p10_m"], ts_all["gwl_p90_m"], color="tab:cyan", alpha=0.3)
ax[3].plot(ts_all.index, ts_all["gwl_mean_m"], c="tab:blue"); ax[3].set_ylabel("CATHY water table (m)")
fig.suptitle(f"{prj.name or prj.parent.name}: LAI mode {a.lai_mode}, version {a.version}")
fig.savefig(out_dir / "fig_overview.png", dpi=150)
if not a.no_plots:
    try:
        PL.make_postrun_plots(out_dir, gdf=gdf, log=log)
        if not a.no_native:
            PL.native_output_checks(simu, out_dir / "figures", rep_)
    except Exception as e_:                               # noqa: BLE001  (plots never invalidate a finished run)
        log(f"post-run plots failed: {type(e_).__name__}: {e_}")
if not a.no_checks:
    rep_.write(chk_dir / "checks_report.md", "Input checks and run notes")
if a.show:
    plt.show()
