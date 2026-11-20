"""
Compare the CATHY LAI-loop scenarios with ground truth:
    * eddy-covariance ETa of Lanzhot      Lanzhot_ETa_2015-2020.csv   (daily, mm/day, 2015-01-01 .. 2019-12-31)
    * piezometers GWL1..GWL4 of Lanzhot   GWL_LANZHOT.xlsx            (30 min, m rel. to ground, 2020-01-01 .. 2024-12-31;
                                                                         the "xlsx" is a semicolon CSV, both formats are read)
    * (optional, same location as the model) satellite TSEB ETa of the modelled polygon (CeskeSvycarsko_TSEB-ET-gf.nc)

    python compare_obs_pycathy.py --base CS_LAI_scenarios --scenarios obs_withIrr clim_withIrr
    python compare_obs_pycathy.py ... --piezo-xy 462350 5632100        # also the model node nearest to the piezometers (model CRS)
    python compare_obs_pycathy.py ... --no-tseb --show

Reads <base>/<scenario>/**/outputs/{timeseries_daily.csv, gwl_nodes.nc, run_config.json} written by run_lai_loop_pycathy.py.
Writes <base>/comparison_obs/: eta_01..04, gwl_01..04 figures, metrics_eta.csv, metrics_gwl.csv, summary.md.

IMPORTANT (see README): the observations are at Lanzhot (South Moravia) while the model files are named CeskeSvycarsko
(Bohemian Switzerland). The script prints the distance between the model polygon and --obs-lonlat and the figures say
which data is compared; absolute piezometer levels are also compared after removing the mean ("anomaly" metrics), because
the model datum / mesh depth is not that of the piezometer wells.

Metrics: bias (model - obs), RMSE, Pearson r, NSE, KGE on the common days (daily) and on monthly means.
"""
import argparse
import sys
from pathlib import Path

import matplotlib

HERE = Path(__file__).resolve().parent
DATA_DIR = Path("/home/ben/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/lib/test_tomas")

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", default=str(HERE / "CS_LAI_scenarios"))
ap.add_argument("--scenarios", nargs="+", required=True)
ap.add_argument("--data-dir", default=str(DATA_DIR))
ap.add_argument("--eta-obs", default="Lanzhot_ETa_2015-2020.csv")
ap.add_argument("--gwl-obs", default="GWL_LANZHOT.xlsx")
ap.add_argument("--piezos", nargs="+", default=["GWL1", "GWL2", "GWL3", "GWL4"])
ap.add_argument("--piezo-xy", type=float, nargs=2, default=None, metavar=("X", "Y"),
                help="coordinates (model CRS) of the point of the model compared with the piezometers; default = polygon mean")
ap.add_argument("--tseb-nc", default="CeskeSvycarsko_TSEB-ET-gf.nc", help="satellite ETa NetCDF of the modelled polygon")
ap.add_argument("--tseb-var", default="ET")
ap.add_argument("--no-tseb", action="store_true")
ap.add_argument("--obs-lonlat", type=float, nargs=2, default=[16.95, 48.69], metavar=("LON", "LAT"),
                help="approximate location of the Lanzhot observations (only used to report the distance to the model polygon)")
ap.add_argument("--out", default=None)
ap.add_argument("--show", action="store_true")
a = ap.parse_args()
if not a.show:
    matplotlib.use("Agg")

import json  # noqa: E402

import geopandas as gpd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402

for up in HERE.parents:
    for root in (up / "geoUtils", up):
        if (root / "src" / "geoutils").exists():
            sys.path[:0] = [str(root / "src")]
            break
    else:
        continue
    break
from geoutils import uET  # noqa: E402
import lai_loop_utils as U  # noqa: E402

base, data = Path(a.base), Path(a.data_dir)
out = Path(a.out) if a.out else base / "comparison_obs"
out.mkdir(parents=True, exist_ok=True)
P = lambda n: Path(n) if Path(n).exists() else data / n          # noqa: E731
COL = {n: c for n, c in zip(a.scenarios, plt.rcParams["axes.prop_cycle"].by_key()["color"] * 3)}


def fmt_dates(fig):
    import matplotlib.dates as mdates
    for ax_ in fig.axes:
        try:
            if ax_.xaxis.units is not None:
                loc = mdates.AutoDateLocator(minticks=3, maxticks=8)
                ax_.xaxis.set_major_locator(loc)
                ax_.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
        except Exception:       # noqa: BLE001
            pass


# =============================================================================================
# 1. load
# =============================================================================================
S, CFG, GW = {}, {}, {}
for n in a.scenarios:
    f = sorted((base / n).rglob("timeseries_daily.csv"), key=lambda p: len(p.parts))
    if not f:
        print(f"scenario {n}: no timeseries_daily.csv under {base / n} -> skipped")
        continue
    S[n] = pd.read_csv(f[0], index_col=0, parse_dates=True)
    CFG[n] = json.loads((f[0].parent / "run_config.json").read_text()) if (f[0].parent / "run_config.json").exists() else {}
    if (f[0].parent / "gwl_nodes.nc").exists():
        GW[n] = xr.open_dataset(f[0].parent / "gwl_nodes.nc")
if not S:
    sys.exit("no scenario output to compare")
first = next(iter(S))
print(f"scenarios: {list(S)}  period {S[first].index.min():%Y-%m-%d} -> {S[first].index.max():%Y-%m-%d}")

eta_ec = U.read_obs_eta(P(a.eta_obs))
gwl_raw = U.read_obs_gwl(P(a.gwl_obs))
piezos = [c for c in a.piezos if c in gwl_raw.columns]
gwl_obs = gwl_raw[piezos].resample("D").mean()
print(f"EC ETa   : {eta_ec.index.min():%Y-%m-%d} -> {eta_ec.index.max():%Y-%m-%d}, {len(eta_ec)} days, mean {eta_ec.mean():.2f} mm/d")
print(f"piezos   : {piezos}, {gwl_obs.index.min():%Y-%m-%d} -> {gwl_obs.index.max():%Y-%m-%d}, means "
      + ", ".join(f"{c} {gwl_obs[c].mean():+.2f} m" for c in piezos))

# distance model polygon <-> observation site
tseb, tseb_info = None, None
gdf = None
gj = CFG[first].get("geojson")
if gj and Path(gj).exists():
    gdf = gpd.read_file(gj)
    cen = gdf.to_crs("EPSG:4326").geometry.union_all().centroid if hasattr(gdf.geometry, "union_all") else gdf.to_crs("EPSG:4326").unary_union.centroid
    lon0, lat0 = np.radians(a.obs_lonlat)
    lon1, lat1 = np.radians([cen.x, cen.y])
    dist = 6371.0 * 2 * np.arcsin(np.sqrt(np.sin((lat1 - lat0) / 2) ** 2 + np.cos(lat0) * np.cos(lat1) * np.sin((lon1 - lon0) / 2) ** 2))
    print(f"model polygon centroid ({cen.x:.3f}E, {cen.y:.3f}N) is {dist:.0f} km from the observation site {a.obs_lonlat}")
    if dist > 20:
        print("WARNING: the observations are NOT at the modelled location: the comparison is a consistency check of "
              "magnitude / seasonality / dynamics, not a validation at the site (see README, open questions)")
    if not a.no_tseb and P(a.tseb_nc).exists():
        try:
            tseb, tseb_info = uET.field_daily_series(P(a.tseb_nc), gdf.to_crs(CFG[first].get("crs", gdf.crs)), var=a.tseb_var, role="eta")
            print(f"TSEB ETa : {tseb_info['var']} {tseb.index.min():%Y-%m-%d} -> {tseb.index.max():%Y-%m-%d}, {len(tseb)} days, "
                  f"{tseb_info['n_pixels']} px, mean {tseb.mean():.2f} mm/d")
        except Exception as e_:     # noqa: BLE001
            print("TSEB ETa not usable:", e_)

lai = pd.read_csv(next((base / first).rglob("lai_monthly.csv")), parse_dates=["month"]).set_index("month")["lai_mean"] \
    if any((base / first).rglob("lai_monthly.csv")) else None

# model GWL: polygon mean (timeseries_daily.csv) or node nearest to --piezo-xy
model_gwl = {n: S[n]["gwl_mean_m"] for n in S}
if a.piezo_xy:
    for n, ds in GW.items():
        j = int(np.argmin((ds.x.values - a.piezo_xy[0]) ** 2 + (ds.y.values - a.piezo_xy[1]) ** 2))
        model_gwl[n] = pd.Series(-ds["gwl_depth_m"].isel(node=j).to_series().values, index=pd.DatetimeIndex(ds.date.values))
        print(f"{n}: model node {j} at ({float(ds.x[j]):.0f}, {float(ds.y[j]):.0f}) used for the piezometers")

# =============================================================================================
# 2. ETa
# =============================================================================================
rows = []
for n in S:
    m = S[n]["eta_mm_d"].dropna()
    for ref_name, ref in (("EC Lanzhot", eta_ec), ("TSEB (model polygon)", tseb)):
        if ref is None:
            continue
        for scale, mm_, rr_ in (("daily", m, ref), ("monthly", m.resample("MS").mean(), ref.resample("MS").mean())):
            r = U.metrics(mm_, rr_)
            rows.append(dict(scenario=n, variable="ETa mm/d", reference=ref_name, scale=scale, **r))
met_eta = pd.DataFrame(rows)
met_eta.to_csv(out / "metrics_eta.csv", index=False)

fig, ax = plt.subplots(3 if lai is not None else 2, 1, figsize=(12, 9), sharex=True, layout="constrained")
k = 0
if lai is not None:
    ax[0].plot(lai.index, lai.values, c="tab:green"); ax[0].set_ylabel("LAI polygon mean (model input)"); k = 1
ax[k].plot(eta_ec.index, eta_ec.values, ".", c="k", ms=2, label="EC Lanzhot (daily)")
if tseb is not None:
    ax[k].plot(tseb.index, tseb.values, ".", c="tab:gray", ms=2, label="TSEB satellite (polygon)")
for n in S:
    ax[k].plot(S[n].index, S[n]["eta_mm_d"].rolling(15, center=True, min_periods=5).mean(), c=COL[n], lw=1.2, label=f"CATHY {n} (15-d mean)")
ax[k].set_ylabel("ETa (mm/day)"); ax[k].legend(ncol=3, fontsize=7)
for n in S:
    ax[k + 1].plot(S[n].index, S[n]["eta_mm_d"].cumsum() / 1.0, c=COL[n], label=f"CATHY {n}")
ax[k + 1].set_ylabel("cumulative ETa since start (mm)"); ax[k + 1].legend(fontsize=7)
fmt_dates(fig); fig.savefig(out / "eta_01_timeseries.png", dpi=150)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
mo = lambda s: s.groupby(s.index.month).mean()          # noqa: E731
ax[0].plot(range(1, 13), mo(eta_ec), "k-o", label="EC Lanzhot")
if tseb is not None:
    ax[0].plot(range(1, 13), mo(tseb), "-o", c="tab:gray", label="TSEB")
for n in S:
    ax[0].plot(range(1, 13), mo(S[n]["eta_mm_d"].dropna()), "-o", c=COL[n], label=n)
ax[0].set_xlabel("month"); ax[0].set_ylabel("mean ETa (mm/day)"); ax[0].set_title("mean seasonal cycle (each on its own years)"); ax[0].legend(fontsize=7)
for n in S:
    c = eta_ec.index.intersection(S[n]["eta_mm_d"].dropna().index)
    if len(c):
        ax[1].scatter(eta_ec.loc[c], S[n].loc[c, "eta_mm_d"], s=6, alpha=0.4, c=COL[n], label=f"{n} (n={len(c)})")
lim = [0, max(6.0, float(eta_ec.max()))]
ax[1].plot(lim, lim, c="gray", lw=0.8); ax[1].set_xlabel("EC ETa (mm/day)"); ax[1].set_ylabel("CATHY ETa (mm/day)"); ax[1].legend(fontsize=7)
fig.savefig(out / "eta_02_seasonal_scatter.png", dpi=150)

# effect of the LAI change: scenarios sharing a version, observed vs climatology
pairs = [(n, n.replace("obs_", "clim_", 1)) for n in S if n.startswith("obs_") and n.replace("obs_", "clim_", 1) in S]
if pairs:
    fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True, layout="constrained")
    for no, nc in pairs:
        d = (S[no]["eta_mm_d"] - S[nc]["eta_mm_d"]).dropna()
        ax[0].plot(d.resample("MS").mean().index, d.resample("MS").mean().values, c=COL[no], label=f"{no} - {nc}")
        ax[1].plot(d.index, d.cumsum(), c=COL[no], label=f"{no} - {nc}")
    ax[0].axhline(0, c="gray", lw=0.8); ax[0].set_ylabel("monthly mean dETa (mm/day)"); ax[0].legend(fontsize=8)
    ax[0].set_title("effect of the LAI change on CATHY ETa (observed LAI - LAI climatology)")
    ax[1].set_ylabel("cumulative dETa (mm)")
    fmt_dates(fig); fig.savefig(out / "eta_03_lai_effect.png", dpi=150)

# =============================================================================================
# 3. groundwater level
# =============================================================================================
rows = []
for n in S:
    mg = model_gwl[n].dropna()
    for pz in piezos:
        o = gwl_obs[pz].dropna()
        for scale, mm_, oo_ in (("daily", mg, o), ("monthly", mg.resample("MS").mean(), o.resample("MS").mean())):
            for dem in (False, True):
                r = U.metrics(mm_, oo_, demean=dem)
                rows.append(dict(scenario=n, piezometer=pz, scale=scale, anomaly=dem, **r))
    ob = gwl_obs[piezos].mean(axis=1).dropna()
    for dem in (False, True):
        rows.append(dict(scenario=n, piezometer="mean of piezos", scale="daily", anomaly=dem, **U.metrics(mg, ob, demean=dem)))
met_gwl = pd.DataFrame(rows)
met_gwl.to_csv(out / "metrics_gwl.csv", index=False)

fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True, layout="constrained")
for pz in piezos:
    ax[0].plot(gwl_obs.index, gwl_obs[pz], lw=0.7, label=f"obs {pz}")
for n in S:
    ax[0].plot(model_gwl[n].index, model_gwl[n].values, c=COL[n], lw=1.4, label=f"CATHY {n}")
    if not a.piezo_xy and "gwl_p10_m" in S[n]:
        ax[0].fill_between(S[n].index, S[n]["gwl_p10_m"], S[n]["gwl_p90_m"], color=COL[n], alpha=0.15, lw=0)
ax[0].set_ylabel("water table (m rel. to the ground)"); ax[0].legend(ncol=4, fontsize=7)
ax[0].set_title("absolute levels (band = p10-p90 of the model nodes of the polygon)")
ob = gwl_obs[piezos].mean(axis=1)
ax[1].plot(ob.index, (ob - ob.mean()).values, c="k", lw=1, label="obs, mean of piezometers (anomaly)")
for n in S:
    g_ = model_gwl[n]
    ax[1].plot(g_.index, (g_ - g_[g_.index.isin(ob.dropna().index)].mean()).values, c=COL[n], lw=1.2, label=f"CATHY {n} (anomaly)")
ax[1].set_ylabel("anomaly (m) = level - mean over the common period"); ax[1].legend(ncol=3, fontsize=7)
fmt_dates(fig); fig.savefig(out / "gwl_01_timeseries.png", dpi=150)

fig, ax = plt.subplots(1, len(S), figsize=(5 * len(S), 4.5), squeeze=False, layout="constrained")
for axi, n in zip(ax[0], S):
    mg = model_gwl[n]
    c = ob.dropna().index.intersection(mg.dropna().index)
    axi.scatter(ob.loc[c] - ob.loc[c].mean(), mg.loc[c] - mg.loc[c].mean(), s=5, alpha=0.4, c=COL[n])
    lim = [min(ob.loc[c].min() - ob.loc[c].mean(), mg.loc[c].min() - mg.loc[c].mean()) if len(c) else -1,
           max(ob.loc[c].max() - ob.loc[c].mean(), mg.loc[c].max() - mg.loc[c].mean()) if len(c) else 1]
    axi.plot(lim, lim, c="gray", lw=0.8)
    m_ = U.metrics(mg, ob, demean=True)
    axi.set_title(f"{n}: anomaly r={m_['r']:.2f} RMSE={m_['rmse']:.2f} m (n={m_['n']})", fontsize=9)
    axi.set_xlabel("obs anomaly (m)"); axi.set_ylabel("model anomaly (m)")
fig.savefig(out / "gwl_02_scatter_anomaly.png", dpi=150)

if pairs:
    fig, ax = plt.subplots(figsize=(12, 4), layout="constrained")
    for no, nc in pairs:
        d = (model_gwl[no] - model_gwl[nc]).dropna()
        ax.plot(d.index, d.values, c=COL[no], label=f"{no} - {nc}")
    ax.axhline(0, c="gray", lw=0.8); ax.set_ylabel("dGWL (m): + = higher water table")
    ax.set_title("effect of the LAI change on the CATHY water table"); ax.legend(fontsize=8)
    fmt_dates(fig); fig.savefig(out / "gwl_03_lai_effect.png", dpi=150)

# =============================================================================================
# 4. summary
# =============================================================================================
def tab(df, cols):
    try:
        return df[cols].round(3).to_markdown(index=False)          # needs `tabulate`
    except ImportError:
        return "```\n" + df[cols].round(3).to_string(index=False) + "\n```"


try:
    md = ["# LAI-loop CATHY vs observations", "",
          f"Scenarios: {', '.join(S)}", "",
          "## ETa (model - reference)", tab(met_eta, ["scenario", "reference", "scale", "n", "bias", "rmse", "r", "nse", "kge"]), "",
          "## Groundwater level (daily; anomaly = means removed)",
          tab(met_gwl[met_gwl["scale"] == "daily"], ["scenario", "piezometer", "anomaly", "n", "bias", "rmse", "r", "nse", "kge"]), ""]
    (out / "summary.md").write_text("\n".join(md))
except Exception as e_:         # noqa: BLE001  (tabulate missing -> to_string fallback already used)
    print("summary.md:", e_)
print("\nETa metrics\n", met_eta[met_eta["scale"] == "daily"].round(3).to_string(index=False))
print("\nGWL metrics (daily, anomaly)\n",
      met_gwl[(met_gwl["scale"] == "daily") & met_gwl["anomaly"]].round(3).to_string(index=False))
print(f"\nfigures / tables in {out}")
if a.show:
    plt.show()
