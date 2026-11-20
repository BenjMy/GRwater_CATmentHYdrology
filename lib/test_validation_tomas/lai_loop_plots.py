"""
lai_loop_plots.py - input checks and output plots of the CATHY "loop over LAI" experiment
=========================================================================================
Used by run_lai_loop_pycathy.py (input checks before the run, output plots after it) and compare_obs_pycathy.py
(observation QC, extra comparison plots, scenario-to-scenario plots).  Can also be run on its own on finished outputs:

    python lai_loop_plots.py --outputs CS_LAI_scenarios/obs_withIrr/<...>/outputs       # re-draw all post-run figures

Re-used (nothing re-implemented when it exists):
    pyCATHY.plotters.cathy_plots via simu.show_input("dem" | "root_map" | "soil" | "atmbc"), simu.show_bc(),
        simu.show("WTD" | "spatialET" | "hgraph" | "cumflowvol" | "wtdepth")        -> native_*.png
    geoutils.uET.check_temporal_quality                                              -> gaps / cadence of every input series

Every figure goes to a folder given by the caller; every numeric check is appended to a Report (OK / INFO / WARN)
that is printed and written to checks_report.md.  Nothing here needs CATHY to be run (except the native_* helpers).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

MONTHS = "JFMAMJJASOND"


# =============================================================================================
# 0. small helpers
# =============================================================================================
class Report:
    """Collects the numeric checks (level OK / INFO / WARN) -> console + markdown."""

    def __init__(self, log=print):
        self.items, self.log = [], log

    def add(self, level, msg):
        self.items.append((level, msg))
        self.log(f"  [{level:4s}] {msg}")

    def ok(self, msg): self.add("OK", msg)
    def info(self, msg): self.add("INFO", msg)
    def warn(self, msg): self.add("WARN", msg)

    @property
    def n_warn(self):
        return sum(1 for lv, _ in self.items if lv == "WARN")

    def write(self, path, title="Input checks"):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        lines = [f"# {title}", "", f"{self.n_warn} warning(s)", ""] + [f"- **{lv}** {m}" for lv, m in self.items]
        path.write_text("\n".join(lines) + "\n")


def save(fig, outdir, name, dpi=130):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / f"{name}.png"
    fig.savefig(p, dpi=dpi)
    plt.close(fig)
    return p


def capture_figs(outdir, prefix, fn, *args, **kw):
    """Run ``fn`` (e.g. a pyCATHY plotter that opens its own figure), save every new figure as <prefix>[_i].png."""
    before = set(plt.get_fignums())
    try:
        fn(*args, **kw)
    except Exception as e:                       # noqa: BLE001  (a native plot must never stop the run)
        print(f"  [skip] {prefix}: {type(e).__name__}: {e}")
    new = sorted(set(plt.get_fignums()) - before)
    for i, n in enumerate(new):
        fig = plt.figure(n)
        fig.savefig(Path(outdir) / f"{prefix}{'' if len(new) == 1 else f'_{i + 1}'}.png", dpi=120)
        plt.close(fig)
    return len(new)


def _outline(ax, gdf):
    if gdf is not None:
        try:
            gdf.boundary.plot(ax=ax, color="k", lw=1.0)
        except Exception:                        # noqa: BLE001
            pass


def _axes(n, ncol=3, size=(4.3, 3.8)):
    nrow = int(np.ceil(n / ncol))
    fig, ax = plt.subplots(nrow, min(ncol, n), figsize=(size[0] * min(ncol, n), size[1] * nrow), squeeze=False,
                           layout="constrained")
    ax = ax.ravel()
    for a_ in ax[n:]:
        a_.set_visible(False)
    return fig, ax


def _map(ax, X, Y, Z, gdf=None, title="", cmap="viridis", vmin=None, vmax=None, label="", cbar=True):
    m = ax.pcolormesh(X, Y, np.ma.masked_invalid(Z), cmap=cmap, vmin=vmin, vmax=vmax, shading="nearest")
    _outline(ax, gdf)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=6)
    ax.ticklabel_format(style="plain", useOffset=False)
    if cbar:
        ax.figure.colorbar(m, ax=ax, shrink=0.8, label=label)
    return m


def node_grid(x, y, v):
    """Scattered structured-mesh nodes (x, y, v) -> 2-D arrays (X, Y, V) for pcolormesh."""
    x, y = np.round(np.asarray(x, float), 3), np.round(np.asarray(y, float), 3)
    xs, ys = np.unique(x), np.unique(y)
    V = np.full((len(ys), len(xs)), np.nan)
    V[np.searchsorted(ys, y), np.searchsorted(xs, x)] = v
    X, Y = np.meshgrid(xs, ys)
    return X, Y, V


def _fmt(fig):
    for a_ in fig.axes:
        try:
            if a_.xaxis.units is not None:
                loc = mdates.AutoDateLocator(minticks=3, maxticks=8)
                a_.xaxis.set_major_locator(loc)
                a_.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
        except Exception:                        # noqa: BLE001
            pass


def temporal_check(series: pd.Series, name: str, rep: Report, freq="D"):
    """geoutils.uET.check_temporal_quality on a pandas series: reports cadence, gaps (> 1 step) and NaN."""
    try:
        from geoutils import uET
        s = series.dropna()
        da = xr.DataArray(s.values, coords={"time": pd.DatetimeIndex(s.index)}, dims="time", name=name)
        r = uET.check_temporal_quality(da, time_dim="time", expected_freq=freq, report=False, plot=False)
        gaps = r["gap_df"]
        n_gap = 0 if gaps is None else len(gaps)
        miss = int(gaps["missing_steps"].sum()) if n_gap else 0
        msg = f"{name}: {s.index.min():%Y-%m-%d} -> {s.index.max():%Y-%m-%d}, {len(s)} values, {n_gap} gap(s), {miss} missing step(s)"
        if n_gap:
            big = gaps.sort_values("missing_steps", ascending=False).iloc[0]
            msg += f"; longest {int(big['missing_steps'])} step(s) from {pd.Timestamp(big['gap_start']):%Y-%m-%d}"
        (rep.warn if miss > 0.05 * max(len(s), 1) else rep.ok)(msg)
        return r
    except Exception as e:                       # noqa: BLE001
        rep.info(f"{name}: temporal check not available ({type(e).__name__}: {e})")


def lai_class(lai, thresholds):
    """vegetation class 1..K, same rule as lai_loop_utils.lai_to_veg_map."""
    return 1 + np.digitize(lai, thresholds, right=True)


def _monthly(ts: pd.DataFrame, sum_cols=(), mean_cols=()):
    """monthly table: totals (mean daily rate x days in month) for sum_cols, means for mean_cols."""
    out = {}
    for c in list(sum_cols) + list(mean_cols):
        if c in ts:
            m = ts[c].resample("MS").mean()
            out[c] = m * m.index.days_in_month if c in sum_cols else m
    return pd.DataFrame(out)


def annual_budget(ts: pd.DataFrame, min_days=300) -> pd.DataFrame:
    """Calendar-year budget from the daily series: rain / ETp / ETa (mm), ETa/ETp, LAI, water-table level."""
    rows = []
    for y, g in ts.groupby(ts.index.year):
        n = int(g["eta_mm_d"].notna().sum()) if "eta_mm_d" in g else len(g)
        ndy = 366 if pd.Timestamp(y, 12, 31).dayofyear == 366 else 365
        r = dict(year=y, days=n, partial=n < min_days)
        for c, k in (("rain_mm_d", "rain_mm"), ("etp_mm_d", "etp_mm"), ("eta_mm_d", "eta_mm"), ("eta_raw_mm_d", "eta_raw_mm")):
            if c in g:
                r[k] = g[c].mean() * ndy
        r["eta_over_etp"] = r.get("eta_mm", np.nan) / r["etp_mm"] if r.get("etp_mm") else np.nan
        for c, k in (("lai_mean", "lai_mean"), ("gwl_mean_m", "gwl_mean_m")):
            if c in g:
                r[k] = g[c].mean()
        if "gwl_mean_m" in g:
            r["gwl_min_m"], r["gwl_max_m"] = g["gwl_mean_m"].min(), g["gwl_mean_m"].max()
        rows.append(r)
    return pd.DataFrame(rows).set_index("year")


# =============================================================================================
# 1. INPUT CHECKS (before the run)
# =============================================================================================
def plot_domain(node_x, node_y, z_surf, z_bottom, node_mask, gdf, outdir, rep: Report, dx=None):
    """Mesh surface elevation, nodes inside the polygon, elevation histogram."""
    X, Y, Z = node_grid(node_x, node_y, z_surf)
    _, _, Mk = node_grid(node_x, node_y, node_mask.astype(float))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), layout="constrained")
    _map(ax[0], X, Y, Z, gdf, "mesh surface elevation", "terrain", label="z (m)")
    _map(ax[1], X, Y, Mk, gdf, f"surface nodes inside polygon ({int(node_mask.sum())} / {len(node_mask)})", "Greys", 0, 1, cbar=False)
    ax[2].hist(z_surf[node_mask], bins=30, color="tab:brown", alpha=0.8, label="polygon nodes")
    ax[2].hist(z_surf[~node_mask], bins=30, color="gray", alpha=0.4, label="outside")
    ax[2].set_xlabel("surface elevation (m)"); ax[2].legend(fontsize=8)
    save(fig, outdir, "check_01_domain")
    rng = float(z_surf[node_mask].max() - z_surf[node_mask].min())
    rep.info(f"surface elevation in polygon {z_surf[node_mask].min():.1f} - {z_surf[node_mask].max():.1f} m (range {rng:.1f} m); "
             f"mesh depth {np.mean(z_surf) - z_bottom:.2f} m")
    if rng < 1e-6:
        rep.warn("flat surface (plane DEM or DEM not read): no topographic gradient in the model")
    if dx and rng / max(dx, 1) > 1.0:
        rep.warn(f"very steep surface (range {rng:.0f} m over cells of {dx:g} m): check the DEM CRS / units")


def plot_lai_inputs(months, lai_obs, lai_cube, n_obs, cell_mask, xc, yc, gdf, outdir, rep: Report, stats=None,
                    thresholds=(), lai_mode="observed", ref_years=(), min_obs=1, lai_scale=1.0):
    """LAI checks: raw acquisitions vs monthly series, acquisitions per month, class fractions, year x month heatmap,
    JJA pixel distribution per year, maps."""
    months = pd.DatetimeIndex(months)
    po, pc = lai_obs[:, cell_mask], lai_cube[:, cell_mask]
    m_o, m_c = po.mean(axis=1), pc.mean(axis=1)
    stats = stats or {}
    ref = [y for y in ref_years if (months.year == y).any()]

    # ---- fig 1: time series ------------------------------------------------------------------
    fig, ax = plt.subplots(3, 1, figsize=(13, 9), sharex=True, layout="constrained", height_ratios=[3, 1.3, 1.7])
    if "acq_times" in stats and len(stats["acq_times"]):
        ax[0].plot(stats["acq_times"], stats["acq_mean"], ".", c="gray", ms=2, alpha=0.5, label="single acquisitions (bbox mean, raw)")
    ax[0].fill_between(months, np.percentile(po, 10, axis=1), np.percentile(po, 90, axis=1), color="tab:green", alpha=0.2, lw=0,
                       label="p10-p90 of the cells (observed monthly)")
    ax[0].plot(months, m_o, c="tab:green", label="observed, monthly (carried forward if no data)")
    if lai_mode != "observed" or lai_scale != 1.0:
        ax[0].plot(months, m_c, c="tab:red", lw=1.5, label=f"LAI used by CATHY ({lai_mode}, x{lai_scale:g})")
    for th in thresholds:
        ax[0].axhline(th, c="k", ls=":", lw=0.7)
    for y in ref:
        ax[0].axvspan(pd.Timestamp(y, 1, 1), pd.Timestamp(y, 12, 31), color="gold", alpha=0.12, lw=0)
    ax[0].set_ylabel("LAI polygon mean"); ax[0].legend(fontsize=7, ncol=2)
    ax[0].set_title("LAI input (dotted = class thresholds; yellow = reference years of the climatology)")
    if "n_raw" in stats:
        ax[1].bar(months, stats["n_raw"], width=20, color="lightgray", label="acquisitions in the month")
    ax[1].bar(months, n_obs, width=20, color="tab:green", label="with valid pixels")
    bad = n_obs < min_obs
    ax[1].plot(months[bad], np.zeros(bad.sum()), "rv", label=f"< {min_obs} valid -> carried forward")
    ax[1].set_ylabel("acquisitions / month"); ax[1].legend(fontsize=7, ncol=3)
    cls = lai_class(pc, thresholds) if len(thresholds) else np.ones_like(pc, int)
    K = int(cls.max()) if len(thresholds) == 0 else len(thresholds) + 1
    fr = np.stack([(cls == k).mean(axis=1) for k in range(1, K + 1)])
    ax[2].stackplot(months, fr, labels=[f"class {k}" for k in range(1, K + 1)], colors=plt.cm.YlGn(np.linspace(0.1, 0.95, K)))
    ax[2].set_ylabel("cell fraction"); ax[2].set_ylim(0, 1); ax[2].legend(fontsize=7, ncol=K, loc="upper center")
    _fmt(fig)
    save(fig, outdir, "check_02_lai_timeseries")

    # ---- fig 2: heatmap + JJA distributions ----------------------------------------------------
    yrs = sorted(set(months.year))
    H = np.full((len(yrs), 12), np.nan)
    for i, mo in enumerate(months):
        H[yrs.index(mo.year), mo.month - 1] = m_o[i]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6), layout="constrained")
    im = ax[0].imshow(H, aspect="auto", cmap="YlGn", origin="upper")
    ax[0].set_xticks(range(12), list(MONTHS)); ax[0].set_yticks(range(len(yrs)), yrs)
    fig.colorbar(im, ax=ax[0], label="observed LAI (polygon mean)")
    ax[0].set_title("year x month (a pre-/post-drop contrast shows up here)")
    jja = np.isin(months.month, (6, 7, 8))
    data, labs = [], []
    for y in yrs:
        sel = jja & (months.year == y)
        if sel.any():
            data.append(po[sel].ravel()); labs.append(str(y))
    ax[1].boxplot(data, tick_labels=labs, showfliers=False)
    ax[1].set_ylabel("observed LAI, cells x JJA months"); ax[1].set_title("summer LAI distribution per year")
    ax[1].tick_params(axis="x", rotation=45)
    if ref:
        a_ = po[jja & np.isin(months.year, ref)].ravel()
        b_ = po[jja & (months.year > max(ref))].ravel()
        if len(b_):
            bins = np.linspace(0, max(np.percentile(a_, 99.5), np.percentile(b_, 99.5), 1), 30)
            ax[2].hist(a_, bins, alpha=0.6, color="tab:green", density=True, label=f"ref years {ref} (JJA)")
            ax[2].hist(b_, bins, alpha=0.6, color="tab:brown", density=True, label=f"after {max(ref)} (JJA)")
            ax[2].set_xlabel("LAI"); ax[2].legend(fontsize=8); ax[2].set_title("LAI before / after")
    else:
        ax[2].axis("off")
    save(fig, outdir, "check_03_lai_heatmap_distribution")

    # ---- fig 3: maps ------------------------------------------------------------------------------
    sel_years = yrs if len(yrs) <= 6 else [yrs[i] for i in np.unique(np.linspace(0, len(yrs) - 1, 6).round().astype(int))]
    idx = []
    for y in sel_years:
        c = np.where((months.year == y) & (months.month == 7))[0]
        idx.append(int(c[0]) if c.size else int(np.where(months.year == y)[0][0]))
    fig, ax = _axes(len(idx))
    vmax = float(np.nanpercentile(lai_cube[idx][:, cell_mask], 99)) if len(idx) else 1.0
    for a_, i in zip(ax, idx):
        Z = np.where(cell_mask, lai_cube[i], lai_cube[i])
        _map(a_, xc, yc, Z, gdf, f"LAI used, {months[i]:%Y-%m}  (polygon mean {pc[i].mean():.2f})", "YlGn", 0, vmax, "LAI")
    save(fig, outdir, "check_04_lai_maps")

    # ---- numeric checks ------------------------------------------------------------------------------
    rep.info(f"LAI (observed) polygon mean {m_o.min():.2f} - {m_o.max():.2f}; cells {po.min():.2f} - {po.max():.2f}")
    if po.max() > 8:
        rep.warn(f"LAI max {po.max():.1f} > 8: unrealistic for Sentinel-2 LAI products (units / scaling?)")
    if m_o.max() < 1:
        rep.warn("LAI never exceeds 1: check the file / variable (--lai-var) or a scale factor")
    if "n_raw" in stats:
        nz = int((np.asarray(stats["n_raw"]) == 0).sum())
        (rep.warn if nz else rep.ok)(f"{nz} month(s) without any acquisition in the file")
    nb = int(bad.sum())
    (rep.warn if nb > 0.2 * len(months) else rep.info)(f"{nb} / {len(months)} month(s) with < {min_obs} valid acquisition(s) (carried forward)")
    if "cloud_frac" in stats:
        rep.info(f"mean fraction of invalid (NaN) pixels in the acquisitions: {np.nanmean(stats['cloud_frac']):.0%}")
    if ref:
        a_m = m_o[jja & np.isin(months.year, ref)].mean()
        late = jja & (months.year > max(ref))
        if late.any():
            b_m = m_o[late].mean()
            rep.info(f"JJA LAI: ref years {a_m:.2f} -> later years {b_m:.2f} ({(b_m - a_m) / a_m:+.0%})")
    if lai_mode == "climatology":
        spread = np.abs(lai_cube[months.year == yrs[0]] - lai_cube[months.year == yrs[-1]]).max() if len(yrs) > 1 else 0
        (rep.ok if spread < 1e-6 else rep.warn)(f"climatology cube identical in every year (max difference {spread:.2g})")
    z_frac = float((pc < 0.05).mean())
    if z_frac > 0.3:
        rep.warn(f"{z_frac:.0%} of the (cell, month) LAI values are ~0 (bare / filled with 0?)")


def plot_forcing_inputs(rain_s, etp_s, etr_s, lai_monthly, etp_mode, kc_params, thresholds, outdir, rep: Report,
                        window=None, rain_factor=1.0):
    """Rain / ETp (product) / ETr checks: daily rain, monthly water balance, annual totals, seasonal cycle,
    ETp-to-ETr ratio against LAI (is the product canopy-aware?), Kc(LAI) curve."""
    if window is not None:
        t0, t1 = pd.Timestamp(window[0]), pd.Timestamp(window[1])
        rain_s, etp_s = rain_s[t0:t1], etp_s[t0:t1]
        etr_s = etr_s[t0:t1] if etr_s is not None else None
    rm = rain_s.resample("MS").sum(min_count=1)
    em = etp_s.resample("MS").sum(min_count=1)
    xm = etr_s.resample("MS").sum(min_count=1) if etr_s is not None else None

    fig, ax = plt.subplots(3, 1, figsize=(13, 9.5), sharex=True, layout="constrained")
    ax[0].bar(rain_s.index, rain_s.values, width=1.0, color="tab:blue")
    ax[0].set_ylabel("daily rain (mm/d)")
    ax[0].set_title(f"rain (x{rain_factor:g}), ETp ({etp_mode}) and ETr: daily rain, monthly totals, cumulative P - ETp")
    ax[1].bar(rm.index, rm.values, width=20, color="tab:blue", alpha=0.6, label="rain")
    ax[1].plot(em.index, em.values, c="k", label="ETp used")
    if xm is not None:
        ax[1].plot(xm.index, xm.values, c="tab:orange", lw=1, label="ETr (FAO56 reference)")
    ax[1].set_ylabel("monthly total (mm)"); ax[1].legend(fontsize=8, ncol=3)
    cum = (rain_s - etp_s.reindex(rain_s.index)).fillna(0).cumsum()
    ax[2].plot(cum.index, cum.values, c="tab:purple")
    ax[2].axhline(0, c="gray", lw=0.8); ax[2].set_ylabel("cumulative P - ETp (mm)")
    _fmt(fig)
    save(fig, outdir, "check_05_forcing_timeseries")

    # annual + seasonal + diagnostics
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), layout="constrained")
    yr = pd.DataFrame({"rain": rain_s.groupby(rain_s.index.year).sum(), "ETp": etp_s.groupby(etp_s.index.year).sum()})
    if etr_s is not None:
        yr["ETr"] = etr_s.groupby(etr_s.index.year).sum()
    cnt = rain_s.groupby(rain_s.index.year).count()
    yr.plot.bar(ax=ax[0], color=["tab:blue", "k", "tab:orange"][: yr.shape[1]], width=0.8)
    ax[0].set_ylabel("annual total (mm)"); ax[0].set_title("annual totals (years with < 300 days are partial)")
    for i, y in enumerate(yr.index):
        if cnt.get(y, 0) < 300:
            ax[0].text(i, 0, "partial", rotation=90, ha="center", va="bottom", fontsize=7, color="red")
    sc = pd.DataFrame({"rain": rm.groupby(rm.index.month).mean(), "ETp": em.groupby(em.index.month).mean()})
    if xm is not None:
        sc["ETr"] = xm.groupby(xm.index.month).mean()
    sc.plot(ax=ax[1], style=["-o", "-o", "-o"][: sc.shape[1]], color=["tab:blue", "k", "tab:orange"][: sc.shape[1]])
    ax[1].set_xlabel("month"); ax[1].set_ylabel("mean monthly total (mm)"); ax[1].set_title("seasonal cycle")
    if etp_mode == "etr_kc" and kc_params:
        L = np.linspace(0, 7, 100)
        kc = kc_params["kc_min"] + (kc_params["kc_max"] - kc_params["kc_min"]) * (1 - np.exp(-kc_params["k_ext"] * L))
        ax[2].plot(L, kc, c="tab:green")
        for th in thresholds:
            ax[2].axvline(th, c="k", ls=":", lw=0.7)
        ax[2].set_xlabel("LAI"); ax[2].set_ylabel("Kc = ETp / ETr"); ax[2].set_title("Kc(LAI) used (etr_kc); dotted = class limits")
    elif xm is not None and lai_monthly is not None:
        ratio = (em / xm).replace([np.inf, -np.inf], np.nan)
        j = pd.concat([ratio.rename("r"), lai_monthly.rename("lai")], axis=1).dropna()
        s_ = ax[2].scatter(j["lai"], j["r"], c=j.index.month, cmap="twilight", s=14)
        fig.colorbar(s_, ax=ax[2], label="month")
        ax[2].set_xlabel("LAI (polygon mean)"); ax[2].set_ylabel("ETp product / ETr (monthly)")
        ax[2].set_title("does the ETp product follow the canopy?")
    else:
        ax[2].axis("off")
    save(fig, outdir, "check_06_forcing_annual_seasonal")

    # numeric checks
    ann = rain_s.groupby(rain_s.index.year).sum()
    full = ann[cnt.reindex(ann.index) >= 300]
    rep.info("annual rain (mm): " + ", ".join(f"{y}: {v:.0f}" for y, v in ann.items()))
    if len(full) and (full.min() < 250 or full.max() > 1500):
        rep.warn(f"annual rain {full.min():.0f} - {full.max():.0f} mm outside 250-1500 mm: check units (--rain-factor, now {rain_factor:g})")
    elif len(full):
        rep.ok(f"annual rain within a plausible range ({full.min():.0f} - {full.max():.0f} mm)")
    if rain_s.max() > 150:
        rep.warn(f"daily rain max {rain_s.max():.0f} mm/d ({rain_s.idxmax():%Y-%m-%d}): unit problem or a real extreme?")
    if (rain_s < 0).any():
        rep.warn(f"{int((rain_s < 0).sum())} negative rain value(s) (clipped to 0 in the run)")
    ea = etp_s.groupby(etp_s.index.year).sum()
    ea = ea[etp_s.groupby(etp_s.index.year).count() >= 300]
    rep.info(f"ETp mean {etp_s.mean():.2f} mm/d, max {etp_s.max():.1f}; annual (full years): " + ", ".join(f"{y}: {v:.0f}" for y, v in ea.items()))
    if etp_s.max() > 12 or (len(ea) and (ea.min() < 250 or ea.max() > 1300)):
        rep.warn("ETp outside the usual range for Central Europe (max 12 mm/d, 250-1300 mm/yr): check units")
    if (etp_s < 0).any():
        rep.warn(f"{int((etp_s < 0).sum())} negative ETp value(s)")
    if etr_s is not None and etp_mode == "product":
        r_ = float(etp_s.mean() / etr_s.mean())
        rep.info(f"ETp product / ETr = {r_:.2f} on average (< 1 => the product is already canopy-reduced; Kc = 1 would be pure FAO56 grass)")
    temporal_check(rain_s, "rain", rep)
    temporal_check(etp_s, "ETp", rep)
    if etr_s is not None:
        temporal_check(etr_s, "ETr", rep)


def plot_etp_maps(months, etp_month_mean, xc, yc, gdf, outdir, rep: Report):
    """Spatial ETp (mm/day) on the CATHY cells: mean per year (small multiples) + spatial dispersion."""
    months = pd.DatetimeIndex(months)
    yrs = sorted(set(months.year))
    fig, ax = _axes(len(yrs) + 1)
    vmax = float(np.nanpercentile(etp_month_mean, 99))
    for a_, y in zip(ax, yrs):
        _map(a_, xc, yc, np.nanmean(etp_month_mean[months.year == y], axis=0), gdf, f"mean ETp {y} (mm/d)", "YlOrBr", 0, vmax, "mm/d")
    pm = np.nanmean(etp_month_mean, axis=0)
    _map(ax[len(yrs)], xc, yc, pm - np.nanmean(pm), gdf, "period-mean ETp minus field mean (mm/d)", "RdBu_r", label="mm/d")
    save(fig, outdir, "check_07_etp_maps")
    spread = np.nanmax(np.nanmean(etp_month_mean, axis=0)) - np.nanmin(np.nanmean(etp_month_mean, axis=0))
    rep.info(f"spatial spread of the period-mean ETp over the cells: {spread:.2f} mm/d "
             f"({'nearly uniform: the EO ETp pixel is coarser than the mesh' if spread < 0.05 else 'spatially variable'})")


def plot_month_inputs(month, tag, lai_c, veg, etp_cells, rain_mm, net_node, node_x, node_y, node_mask, cell_mask, xc, yc,
                      gdf, thresholds, zroot, outdir, rep: Report):
    """One month, exactly what is written to CATHY: LAI map, vegetation class map, ETp map, net flux map (rain - ETp
    on the nodes), daily series, ZROOT table; plus a mass check cells -> nodes."""
    days = pd.date_range(month, periods=len(rain_mm), freq="D")
    etp_poly = etp_cells[:, cell_mask].mean(axis=1)
    net_mm_d = net_node.mean(axis=0) * 1e3 * 86400.0
    K = len(thresholds) + 1
    fig, ax = plt.subplots(2, 3, figsize=(16, 9), layout="constrained")
    ax = ax.ravel()
    _map(ax[0], xc, yc, lai_c, gdf, f"LAI {month:%Y-%m}", "YlGn", 0, None, "LAI")
    cm = plt.matplotlib.colors.ListedColormap(plt.cm.YlGn(np.linspace(0.15, 0.95, K)))
    m = ax[1].pcolormesh(xc, yc, veg, cmap=cm, vmin=0.5, vmax=K + 0.5, shading="nearest")
    _outline(ax[1], gdf); ax[1].set_aspect("equal"); ax[1].tick_params(labelsize=6)
    cb = fig.colorbar(m, ax=ax[1], ticks=range(1, K + 1), shrink=0.8)
    cb.ax.set_yticklabels([f"{k}  (ZROOT {zroot[k - 1]:.2f} m)" for k in range(1, K + 1)])
    ax[1].set_title("vegetation class sent to CATHY", fontsize=9)
    _map(ax[2], xc, yc, etp_cells.mean(axis=0), gdf, "ETp month mean (cells, mm/d)", "YlOrBr", 0, None, "mm/d")
    X, Y, V = node_grid(node_x, node_y, net_mm_d)
    _map(ax[3], X, Y, V, gdf, "net flux rain - ETp on the nodes (mm/d, month mean)", "RdBu", -np.nanmax(np.abs(V)) or -1, np.nanmax(np.abs(V)) or 1, "mm/d")
    ax[4].bar(days, rain_mm, color="tab:blue", label="rain (mm/d)")
    ax[4].plot(days, etp_poly, c="k", label="ETp polygon mean")
    ax[4].plot(days, rain_mm - etp_poly, c="tab:purple", lw=0.8, label="rain - ETp")
    ax[4].axhline(0, c="gray", lw=0.6); ax[4].legend(fontsize=7); ax[4].tick_params(axis="x", rotation=45, labelsize=7)
    ax[4].set_title("daily forcing of the month")
    fr = [(veg[cell_mask] == k).mean() for k in range(1, K + 1)]
    ax[5].bar(range(1, K + 1), zroot, color=plt.cm.YlGn(np.linspace(0.15, 0.95, K)))
    for k, f in enumerate(fr, start=1):
        ax[5].text(k, zroot[k - 1], f"{f:.0%}\nof cells", ha="center", va="bottom", fontsize=8)
    ax[5].set_xlabel("class"); ax[5].set_ylabel("ZROOT (m)"); ax[5].set_title("root depth per class and polygon area share")
    fig.suptitle(f"inputs written to CATHY - {tag} {month:%Y-%m}")
    save(fig, outdir, f"month_inputs_{month:%Y%m}")

    exp_mm = float(rain_mm.sum() - etp_poly.sum())
    got_mm = float(net_node[:, node_mask].mean(axis=1).sum() * 1e3 * 86400.0)
    d = got_mm - exp_mm
    msg = f"{month:%Y-%m} net water input over the polygon: series {exp_mm:.1f} mm, CATHY nodes {got_mm:.1f} mm (diff {d:+.1f} mm)"
    (rep.warn if abs(d) > max(2.0, 0.05 * abs(exp_mm)) else rep.ok)(msg)
    if not np.isfinite(net_node).all():
        rep.warn(f"{month:%Y-%m}: NaN / Inf in the atmbc fluxes")


def native_checks(simu, outdir, tag, rep: Report):
    """Native pyCATHY plots of the inputs currently written in the project (DEM, vegetation map, ZROOT, PERMX, atmbc, BC)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    jobs = [("dem", lambda: simu.show_input("dem")),
            ("veg_map", lambda: simu.show_input("root_map")),
            ("zroot", lambda: simu.show_input("soil", yprop="ZROOT")),
            ("permx", lambda: simu.show_input("soil", yprop="PERMX")),
            ("atmbc", lambda: simu.show_input("atmbc")),
            ("bc", lambda: simu.show_bc())]
    n = sum(capture_figs(outdir, f"native_{tag}_{name}", fn) for name, fn in jobs)
    rep.info(f"pyCATHY native input plots ({tag}): {n} figure(s) in {outdir}")


def native_output_checks(simu, outdir, rep: Report):
    """Native pyCATHY plots of the outputs of the LAST simulated month (water-table depth map, spatial ET, hgraph...)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    n = sum(capture_figs(outdir, f"native_out_{p}", simu.show, p) for p in ("WTD", "spatialET", "hgraph", "cumflowvol"))
    rep.info(f"pyCATHY native output plots (last month): {n} figure(s)")


# =============================================================================================
# 2. OUTPUT PLOTS (after the run; only need the files written by write_outputs)
# =============================================================================================
def _periods(ts, ref_years):
    ref = [y for y in ref_years if y in set(ts.index.year)]
    if not ref:
        return None, None
    return ts[ts.index.year.isin(ref)], ts[ts.index.year > max(ref)]


def plot_timeseries_set(ts: pd.DataFrame, outdir, title="", ref_years=()):
    """fig_01 forcing vs ETa, fig_02 water balance, fig_03 annual budget, fig_04 state, fig_05 LAI response, fig_06 seasonal."""
    outdir = Path(outdir)
    ab = annual_budget(ts)
    ab.round(2).to_csv(outdir / "annual_budget.csv")
    mo = _monthly(ts, sum_cols=("rain_mm_d", "etp_mm_d", "eta_mm_d", "eta_raw_mm_d"), mean_cols=("lai_mean", "gwl_mean_m"))
    mo.columns = [c.replace("_mm_d", "_mm") for c in mo.columns]

    # 01 forcing vs ETa
    fig, ax = plt.subplots(4, 1, figsize=(13, 11), sharex=True, layout="constrained")
    ax[0].plot(mo.index, mo["lai_mean"], c="tab:green"); ax[0].set_ylabel("LAI (polygon mean)")
    ax[1].bar(mo.index, mo["rain_mm"], width=20, color="tab:blue", alpha=0.6, label="rain")
    ax[1].plot(mo.index, mo["etp_mm"], c="gray", label="ETp"); ax[1].plot(mo.index, mo["eta_mm"], c="tab:red", label="ETa CATHY")
    ax[1].set_ylabel("monthly total (mm)"); ax[1].legend(ncol=3, fontsize=8)
    ratio = mo["eta_mm"] / mo["etp_mm"]
    ax[2].plot(ratio.index, ratio.values, c="tab:red", marker=".", ms=3)
    ax[2].axhline(1, c="gray", lw=0.7); ax[2].set_ylabel("ETa / ETp (monthly)"); ax[2].set_ylim(0, 1.3)
    ax[3].fill_between(ts.index, ts["gwl_p10_m"], ts["gwl_p90_m"], color="tab:cyan", alpha=0.3, lw=0, label="p10-p90 nodes")
    ax[3].plot(ts.index, ts["gwl_mean_m"], c="tab:blue", label="mean"); ax[3].set_ylabel("water table (m)"); ax[3].legend(fontsize=8)
    for a_ in ax:
        for y in ref_years:
            a_.axvspan(pd.Timestamp(y, 1, 1), pd.Timestamp(y, 12, 31), color="gold", alpha=0.1, lw=0)
    fig.suptitle(title); _fmt(fig)
    save(fig, outdir, "fig_01_forcing_vs_eta")

    # 02 water balance
    fig, ax = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True, layout="constrained")
    for lab, s, c in (("P - ETp", ts["rain_mm_d"] - ts["etp_mm_d"], "gray"), ("P - ETa (CATHY)", ts["rain_mm_d"] - ts["eta_mm_d"], "tab:red")):
        ax[0].plot(ts.index, s.fillna(0).cumsum(), c=c, label=lab)
    ax[0].set_ylabel("cumulative (mm)"); ax[0].legend(fontsize=8)
    ax[0].set_title("cumulative water input minus demand / actual ET (the drift is what recharges the groundwater or leaves the box)")
    ax[1].plot(ts.index, ts["gwl_mean_m"] - ts["gwl_mean_m"].iloc[0], c="tab:blue")
    ax[1].set_ylabel("water-table change since start (m)")
    _fmt(fig); save(fig, outdir, "fig_02_water_balance")

    # 03 annual budget
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), layout="constrained")
    ab[[c for c in ("rain_mm", "etp_mm", "eta_mm") if c in ab]].plot.bar(ax=ax[0], color=["tab:blue", "gray", "tab:red"], width=0.8)
    ax[0].set_ylabel("mm / year"); ax[0].set_title("annual budget (partial years flagged in annual_budget.csv)")
    ab["eta_over_etp"].plot.bar(ax=ax[1], color="tab:red"); ax[1].set_ylabel("ETa / ETp"); ax[1].set_title("transpiration stress index")
    ab["lai_mean"].plot(ax=ax[2], marker="o", c="tab:green"); ax[2].set_ylabel("LAI", color="tab:green")
    a2 = ax[2].twinx(); ab["gwl_mean_m"].plot(ax=a2, marker="s", c="tab:blue"); a2.set_ylabel("mean water table (m)", color="tab:blue")
    ax[2].set_title("LAI and water table, annual mean")
    save(fig, outdir, "fig_03_annual_budget")

    # 04 state layers
    cols_sw = [c for c in ts if c.startswith("sw_L")]
    cols_psi = [c for c in ts if c.startswith("psi_L")]
    if cols_sw:
        fig, ax = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True, layout="constrained")
        for c in cols_sw:
            ax[0].plot(ts.index, ts[c], lw=0.9, label=c)
        ax[0].set_ylabel("saturation Sw"); ax[0].legend(fontsize=8, ncol=3)
        for c in cols_psi:
            ax[1].plot(ts.index, ts[c], lw=0.9, label=c)
        ax[1].set_ylabel("pressure head (m)"); ax[1].legend(fontsize=8, ncol=3)
        ax[1].set_title("layer means over the polygon (L1 = shallow ... last = deep)")
        _fmt(fig); save(fig, outdir, "fig_04_state_layers")

    # 05 LAI response + patch
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), layout="constrained")
    j = mo.dropna(subset=["lai_mean", "eta_mm", "etp_mm"])
    s_ = ax[0].scatter(j["lai_mean"], j["eta_mm"], c=j.index.year, cmap="viridis", s=18)
    fig.colorbar(s_, ax=ax[0], label="year"); ax[0].set_xlabel("LAI"); ax[0].set_ylabel("monthly ETa (mm)"); ax[0].set_title("ETa vs LAI")
    s_ = ax[1].scatter(j["lai_mean"], j["eta_mm"] / j["etp_mm"], c=j.index.month, cmap="twilight", s=18)
    fig.colorbar(s_, ax=ax[1], label="month"); ax[1].set_xlabel("LAI"); ax[1].set_ylabel("ETa / ETp"); ax[1].set_title("stress vs LAI")
    if "eta_raw_mm_d" in ts:
        d = (ts["eta_mm_d"] - ts["eta_raw_mm_d"]).dropna()
        ax[2].plot(d.resample("MS").mean().index, d.resample("MS").mean().values, c="tab:orange")
        ax[2].set_ylabel("ETa patched - raw (mm/d, monthly mean)")
        ax[2].set_title(f"effect of the 'ETa = ETp where net flux >= 0' patch (days changed: {(d.abs() > 1e-9).mean():.0%})")
    save(fig, outdir, "fig_05_lai_response")

    # 06 seasonal, ref years vs later
    ref, late = _periods(ts, ref_years)
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.3), layout="constrained")
    for a_, (c, lab) in zip(ax, (("lai_mean", "LAI"), ("eta_mm_d", "ETa (mm/d)"), ("etp_mm_d", "ETp (mm/d)"), ("gwl_mean_m", "water table (m)"))):
        if ref is not None and late is not None and len(late):
            a_.plot(range(1, 13), ref[c].groupby(ref.index.month).mean(), "-o", c="tab:green", label=f"ref {sorted(set(ref.index.year))}")
            a_.plot(range(1, 13), late[c].groupby(late.index.month).mean(), "-o", c="tab:brown", label=f"after {max(ref_years)}")
        else:
            a_.plot(range(1, 13), ts[c].groupby(ts.index.month).mean(), "-o", c="k", label="all years")
        a_.set_xlabel("month"); a_.set_ylabel(lab); a_.legend(fontsize=7)
    fig.suptitle("seasonal cycle before / after the reference period")
    save(fig, outdir, "fig_06_seasonal_before_after")
    return ab


def _year_maps(outdir, name, gridder, V, dates, ref_years, title, cmap, label, gdf, diff_label, vsym=None, months_sel=None):
    """small multiples of the yearly mean of V (time, node) + (later - ref) difference."""
    dates = pd.DatetimeIndex(dates)
    sel = np.ones(len(dates), bool) if months_sel is None else np.isin(dates.month, months_sel)
    yrs = sorted(set(dates[sel].year))
    fields = {y: np.nanmean(V[sel & (dates.year == y)], axis=0) for y in yrs if (sel & (dates.year == y)).any()}
    fig, ax = _axes(len(fields) + 1)
    allv = np.concatenate([f[np.isfinite(f)] for f in fields.values()])
    vmin, vmax = np.percentile(allv, 1), np.percentile(allv, 99)
    for a_, (y, f) in zip(ax, fields.items()):
        X, Y, Z = gridder(f)
        _map(a_, X, Y, Z, gdf, f"{title} {y}", cmap, vmin, vmax, label)
    ref = [y for y in ref_years if y in fields]
    late = [y for y in fields if ref and y > max(ref)]
    if ref and late:
        d = np.nanmean([fields[y] for y in late], axis=0) - np.nanmean([fields[y] for y in ref], axis=0)
        v = np.nanmax(np.abs(d)) or 1.0
        X, Y, Z = gridder(d)
        _map(ax[len(fields)], X, Y, Z, gdf, f"{diff_label}: after {max(ref)} - ref {ref}", "RdBu_r", -v, v, label)
    else:
        ax[len(fields)].axis("off")
    save(fig, outdir, name)


def plot_maps_set(outputs: Path, outdir, gdf=None, ref_years=()):
    """Spatial maps from the NetCDF outputs: yearly mean ETa, water-table depth, LAI (JJA) and their changes."""
    outputs, outdir = Path(outputs), Path(outdir)
    f = outputs / "eta_nodes_monthly.nc"
    if f.exists():
        with xr.open_dataset(f) as ds:
            x, y, V, m = ds["x"].values, ds["y"].values, ds["eta_mm_d"].values, ds["month"].values
        X, Y, _ = node_grid(x, y, V[0])
        gr = lambda v: node_grid(x, y, v)                                   # noqa: E731
        _year_maps(outdir, "map_01_eta_yearly", gr, V, m, ref_years, "ETa", "YlGnBu", "mm/d", gdf, "ETa change")
    f = outputs / "gwl_nodes.nc"
    if f.exists():
        with xr.open_dataset(f) as ds:
            x, y, V, d = ds["x"].values, ds["y"].values, ds["gwl_depth_m"].values, ds["date"].values
        gr = lambda v: node_grid(x, y, v)                                   # noqa: E731
        _year_maps(outdir, "map_02_gwl_depth_yearly", gr, V, d, ref_years, "WT depth", "viridis_r", "depth (m)", gdf, "WT depth change")
        _year_maps(outdir, "map_03_gwl_depth_summer", gr, V, d, ref_years, "WT depth JAS", "viridis_r", "depth (m)", gdf,
                   "WT depth change", months_sel=(7, 8, 9))
        dd = pd.DatetimeIndex(d)
        ok = np.isfinite(V).all(axis=0)
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.3), layout="constrained")
        for p, c in ((10, "tab:cyan"), (50, "tab:blue"), (90, "navy")):
            ax[0].plot(dd, np.percentile(V[:, ok], p, axis=1), c=c, label=f"p{p}")
        ax[0].invert_yaxis(); ax[0].set_ylabel("water-table depth (m)"); ax[0].legend(fontsize=8); ax[0].set_title("distribution over the nodes")
        ax[1].hist(V[:, ok].mean(axis=0), bins=30, color="tab:blue")
        ax[1].set_xlabel("time-mean water-table depth (m)"); ax[1].set_title("spatial spread (nodes)")
        _fmt(fig); save(fig, outdir, "map_04_gwl_distribution")
    f = outputs / "lai_cells_monthly.nc"
    if f.exists():
        with xr.open_dataset(f) as ds:
            xc, yc, L, m = ds["x"].values, ds["y"].values, ds["LAI"].values, ds["month"].values
        _year_maps_cells(outdir, "map_05_lai_jja_yearly", xc, yc, L, m, ref_years, gdf)


def _year_maps_cells(outdir, name, xc, yc, L, months, ref_years, gdf):
    months = pd.DatetimeIndex(months)
    sel = np.isin(months.month, (6, 7, 8))
    yrs = sorted(set(months[sel].year))
    fields = {y: L[sel & (months.year == y)].mean(axis=0) for y in yrs}
    fig, ax = _axes(len(fields) + 1)
    vmax = max(np.nanmax(v) for v in fields.values())
    for a_, (y, v) in zip(ax, fields.items()):
        _map(a_, xc, yc, v, gdf, f"LAI JJA {y}", "YlGn", 0, vmax, "LAI")
    ref = [y for y in ref_years if y in fields]
    late = [y for y in fields if ref and y > max(ref)]
    if ref and late:
        d = np.mean([fields[y] for y in late], axis=0) - np.mean([fields[y] for y in ref], axis=0)
        v = np.nanmax(np.abs(d)) or 1.0
        _map(ax[len(fields)], xc, yc, d, gdf, f"LAI change: after {max(ref)} - ref", "RdBu", -v, v, "LAI")
    else:
        ax[len(fields)].axis("off")
    save(fig, outdir, name)


def make_postrun_plots(outputs, outdir=None, gdf=None, log=print):
    """All post-run figures from an ``outputs`` folder of run_lai_loop_pycathy.py. Returns the annual budget table."""
    outputs = Path(outputs)
    outdir = Path(outdir) if outdir else outputs / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((outputs / "run_config.json").read_text()) if (outputs / "run_config.json").exists() else {}
    if gdf is None and cfg.get("geojson") and Path(cfg["geojson"]).exists():
        import geopandas as gpd
        gdf = gpd.read_file(cfg["geojson"])
        gdf = gdf.set_crs("EPSG:4326") if gdf.crs is None else gdf
        gdf = gdf.to_crs(cfg.get("crs", gdf.crs))
    ts = pd.read_csv(outputs / "timeseries_daily.csv", index_col=0, parse_dates=True)
    ref = cfg.get("ref_years", [])
    title = f"{cfg.get('version', '')} | LAI {cfg.get('lai_mode', '')} | ETp {cfg.get('etp_mode', '')}"
    ab = plot_timeseries_set(ts, outdir, title, ref)
    try:
        plot_maps_set(outputs, outdir, gdf, ref)
    except Exception as e:                       # noqa: BLE001
        log(f"  maps skipped: {type(e).__name__}: {e}")
    log(f"post-run figures: {outdir}")
    return ab


# =============================================================================================
# 3. OBSERVATION QC + EXTRA COMPARISON PLOTS (used by compare_obs_pycathy.py)
# =============================================================================================
def plot_obs_qc(eta_ec: pd.Series, gwl_raw: pd.DataFrame, piezos, outdir, rep: Report):
    """Quality of the observations before comparing: EC ETa coverage / distribution, piezometer coverage, spikes, inter-piezo correlation."""
    yrs = sorted(set(eta_ec.index.year))
    H = np.full((len(yrs), 366), np.nan)
    for d, v in eta_ec.items():
        H[yrs.index(d.year), d.dayofyear - 1] = v
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6), layout="constrained")
    im = ax[0].imshow(H, aspect="auto", cmap="YlGnBu", origin="upper", vmin=0, vmax=max(6, np.nanpercentile(H, 99)))
    ax[0].set_yticks(range(len(yrs)), yrs); ax[0].set_xlabel("day of year"); fig.colorbar(im, ax=ax[0], label="EC ETa (mm/d)")
    ax[0].set_title("EC ETa coverage (white = missing)")
    cnt = eta_ec.groupby([eta_ec.index.year, eta_ec.index.month]).count().unstack()
    im = ax[1].imshow(cnt.values, aspect="auto", cmap="Greens", vmin=0, vmax=31)
    ax[1].set_xticks(range(cnt.shape[1]), cnt.columns); ax[1].set_yticks(range(len(cnt.index)), cnt.index)
    fig.colorbar(im, ax=ax[1], label="valid days"); ax[1].set_title("valid days per month")
    ax[2].hist(eta_ec.values, bins=40, color="tab:blue"); ax[2].set_xlabel("EC ETa (mm/d)"); ax[2].set_title("distribution")
    save(fig, outdir, "obs_01_eta_ec_qc")
    rep.info(f"EC ETa: {len(eta_ec)} days, mean {eta_ec.mean():.2f}, max {eta_ec.max():.1f} mm/d")
    if (eta_ec < 0).any():
        rep.warn(f"EC ETa: {int((eta_ec < 0).sum())} negative day(s)")
    if eta_ec.max() > 12:
        rep.warn(f"EC ETa max {eta_ec.max():.1f} mm/d is unusually high")
    temporal_check(eta_ec, "EC ETa", rep)

    daily = gwl_raw[piezos].resample("D").mean()
    fig, ax = plt.subplots(3, 1, figsize=(13, 10), layout="constrained", height_ratios=[3, 1.2, 2])
    for pz in piezos:
        ax[0].plot(gwl_raw.index, gwl_raw[pz], lw=0.4, label=pz)
    ax[0].set_ylabel("level (m rel. ground), raw 30 min"); ax[0].legend(ncol=4, fontsize=8)
    for i, pz in enumerate(piezos):
        v = daily[pz].notna().values
        ax[1].fill_between(daily.index, i, i + 0.8, where=v, step="mid", color=f"C{i}", lw=0)
    ax[1].set_yticks(np.arange(len(piezos)) + 0.4, piezos); ax[1].set_title("data coverage")
    cr = daily.corr()
    im = ax[2].imshow(cr.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax[2].set_xticks(range(len(piezos)), piezos); ax[2].set_yticks(range(len(piezos)), piezos)
    for i in range(len(piezos)):
        for j in range(len(piezos)):
            ax[2].text(j, i, f"{cr.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax[2], label="correlation (daily means)", shrink=0.7)
    save(fig, outdir, "obs_02_piezometers_qc")
    for pz in piezos:
        d1 = daily[pz].diff().dropna()
        n_sp = int((d1.abs() > 5 * d1.std()).sum()) if d1.std() > 0 else 0
        rep.info(f"{pz}: mean {daily[pz].mean():+.2f} m, range {daily[pz].min():+.2f} .. {daily[pz].max():+.2f}, "
                 f"{daily[pz].notna().mean():.0%} of days, {n_sp} day-to-day jump(s) > 5 sigma")
        if daily[pz].max() > 0.5:
            rep.warn(f"{pz}: levels above the ground surface (> +0.5 m): sign convention?")
        temporal_check(daily[pz], pz, rep)


def plot_eta_extra(S, eta_ec, tseb, COL, out, rep: Report | None = None):
    """eta_04 annual means + monthly bias; eta_05 CDF + residuals per month; eta_06 ETa/ETp."""
    out = Path(out)
    ann = {"EC Lanzhot": eta_ec.groupby(eta_ec.index.year).mean()}
    if tseb is not None:
        ann["TSEB"] = tseb.groupby(tseb.index.year).mean()
    for n in S:
        ann[n] = S[n]["eta_mm_d"].dropna().groupby(S[n]["eta_mm_d"].dropna().index.year).mean()
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.6), layout="constrained")
    pd.DataFrame(ann).plot.bar(ax=ax[0], width=0.85, color=["k", "tab:gray"][: 2 if tseb is not None else 1] + [COL[n] for n in S])
    ax[0].set_ylabel("mean ETa (mm/d)"); ax[0].set_title("annual mean ETa (each source on its own coverage)"); ax[0].legend(fontsize=7)
    w = 0.8 / max(len(S), 1)
    for k, n in enumerate(S):
        c = eta_ec.index.intersection(S[n]["eta_mm_d"].dropna().index)
        if len(c):
            b = (S[n].loc[c, "eta_mm_d"] - eta_ec.loc[c])
            bm = b.groupby(b.index.month).mean()
            ax[1].bar(bm.index + (k - (len(S) - 1) / 2) * w, bm.values, width=w, color=COL[n], label=n)
    ax[1].axhline(0, c="gray", lw=0.8); ax[1].set_xlabel("month"); ax[1].set_ylabel("bias model - EC (mm/d)")
    ax[1].set_title("monthly bias on common days"); ax[1].legend(fontsize=7)
    save(fig, out, "eta_04_annual_and_monthly_bias")

    fig, ax = plt.subplots(1, 2, figsize=(14, 4.6), layout="constrained")
    for name, s, c in [("EC", eta_ec, "k")] + ([("TSEB", tseb, "tab:gray")] if tseb is not None else []) + [(n, S[n]["eta_mm_d"].dropna(), COL[n]) for n in S]:
        v = np.sort(s.values)
        ax[0].plot(v, np.linspace(0, 1, len(v)), c=c, label=name)
    ax[0].set_xlabel("daily ETa (mm/d)"); ax[0].set_ylabel("CDF"); ax[0].legend(fontsize=7); ax[0].set_title("distribution of daily ETa")
    pos = 0
    for k, n in enumerate(S):
        c = eta_ec.index.intersection(S[n]["eta_mm_d"].dropna().index)
        if len(c):
            r = S[n].loc[c, "eta_mm_d"] - eta_ec.loc[c]
            bp = ax[1].boxplot([r[r.index.month == m_].values for m_ in range(1, 13)], positions=np.arange(1, 13) + (k - (len(S) - 1) / 2) * w,
                               widths=w * 0.9, showfliers=False, patch_artist=True)
            for b_ in bp["boxes"]:
                b_.set_facecolor(COL[n]); b_.set_alpha(0.6)
    ax[1].axhline(0, c="gray", lw=0.8); ax[1].set_xticks(range(1, 13), list(MONTHS)); ax[1].set_ylabel("residual model - EC (mm/d)")
    ax[1].set_title("daily residuals per month")
    save(fig, out, "eta_05_cdf_residuals")

    fig, ax = plt.subplots(2, 1, figsize=(13, 7), sharex=True, layout="constrained")
    for n in S:
        e, p = S[n]["eta_mm_d"].resample("MS").mean(), S[n]["etp_mm_d"].resample("MS").mean()
        ax[0].plot(e.index, e / p, c=COL[n], label=n)
        ax[1].plot(p.index, p.values, c=COL[n], label=f"ETp {n}")
    ec_m = eta_ec.resample("MS").mean()
    first = next(iter(S))
    p0 = S[first]["etp_mm_d"].resample("MS").mean()
    j = (ec_m / p0).dropna()
    ax[0].plot(j.index, j.values, "k.", label=f"EC / ETp({first})")
    ax[0].axhline(1, c="gray", lw=0.7); ax[0].set_ylabel("ETa / ETp (monthly)"); ax[0].legend(fontsize=7, ncol=3)
    ax[1].plot(ec_m.index, ec_m.values, "k.", label="EC ETa"); ax[1].set_ylabel("mm/d"); ax[1].legend(fontsize=7, ncol=3)
    ax[1].set_title("ETp forcing of each scenario and EC ETa")
    _fmt(fig); save(fig, out, "eta_06_eta_over_etp")


def plot_gwl_extra(S, model_gwl, gwl_obs, piezos, COL, out):
    """gwl_04 seasonal anomaly cycle + level distributions; gwl_05 monthly water-table change vs net input."""
    out = Path(out)
    ob = gwl_obs[piezos].mean(axis=1).dropna()
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6), layout="constrained")
    an = ob - ob.mean()
    ax[0].plot(range(1, 13), an.groupby(an.index.month).mean(), "k-o", label="obs (mean of piezometers)")
    for n in S:
        g = model_gwl[n].dropna()
        g = g - g[g.index.isin(ob.index)].mean() if g.index.isin(ob.index).any() else g - g.mean()
        ax[0].plot(range(1, 13), g.groupby(g.index.month).mean(), "-o", c=COL[n], label=n)
    ax[0].set_xlabel("month"); ax[0].set_ylabel("anomaly (m)"); ax[0].legend(fontsize=7); ax[0].set_title("mean seasonal cycle of the anomaly")
    for pz in piezos:
        ax[1].hist(gwl_obs[pz].dropna(), bins=30, alpha=0.4, label=pz)
    for n in S:
        ax[1].axvline(model_gwl[n].mean(), c=COL[n], lw=2, label=f"model {n}")
    ax[1].set_xlabel("level (m rel. ground)"); ax[1].legend(fontsize=6); ax[1].set_title("level distributions (vertical = model mean)")
    for n in S:
        net = (S[n]["rain_mm_d"] - S[n]["etp_mm_d"]).resample("MS").sum(min_count=1)
        d = model_gwl[n].resample("MS").mean().diff()
        j = pd.concat([net.rename("net"), d.rename("d")], axis=1).dropna()
        ax[2].scatter(j["net"], j["d"], s=10, alpha=0.5, c=COL[n], label=f"model {n}")
    first = next(iter(S))
    net = (S[first]["rain_mm_d"] - S[first]["etp_mm_d"]).resample("MS").sum(min_count=1)
    j = pd.concat([net.rename("net"), ob.resample("MS").mean().diff().rename("d")], axis=1).dropna()
    ax[2].scatter(j["net"], j["d"], s=14, c="k", marker="x", label="obs vs forcing of the first scenario")
    ax[2].axhline(0, c="gray", lw=0.7); ax[2].axvline(0, c="gray", lw=0.7)
    ax[2].set_xlabel("monthly rain - ETp (mm)"); ax[2].set_ylabel("monthly change of the level (m)"); ax[2].legend(fontsize=6)
    ax[2].set_title("level response to the net input")
    save(fig, out, "gwl_04_seasonal_distribution_response")


def plot_metric_bars(met_eta, met_gwl, COL, out):
    """Bars of r / RMSE / KGE per scenario for ETa (daily) and the water-table anomaly (mean of piezos)."""
    out = Path(out)
    fig, ax = plt.subplots(2, 3, figsize=(15, 7), layout="constrained")
    e = met_eta[met_eta["scale"] == "daily"]
    g = met_gwl[(met_gwl["piezometer"] == "mean of piezos") & met_gwl["anomaly"]]
    for row, (df, lab) in enumerate(((e, "ETa daily"), (g, "GWL anomaly daily"))):
        for col, m in enumerate(("r", "rmse", "kge")):
            if df.empty:
                ax[row, col].axis("off")
                continue
            piv = df.pivot_table(index="scenario", columns="reference" if row == 0 else "piezometer", values=m)
            piv.plot.bar(ax=ax[row, col], width=0.8, legend=(col == 0))
            ax[row, col].set_title(f"{lab}: {m}"); ax[row, col].tick_params(axis="x", rotation=30)
    save(fig, out, "metrics_bars")


def plot_scenarios_inputs(S, DIRS, out):
    """What differs between scenarios: LAI actually used, vegetation classes, ETp forcing, annual budget table of all scenarios."""
    out = Path(out)
    COL = {n: c for n, c in zip(S, plt.rcParams["axes.prop_cycle"].by_key()["color"] * 3)}
    fig, ax = plt.subplots(3, 1, figsize=(13, 9.5), sharex=True, layout="constrained")
    rows = []
    for n in S:
        lai_f = DIRS[n] / "lai_monthly.csv"
        if lai_f.exists():
            l = pd.read_csv(lai_f, parse_dates=["month"]).set_index("month")["lai_mean"]
            ax[0].plot(l.index, l.values, c=COL[n], label=n)
        ax[1].plot(S[n].index, S[n]["etp_mm_d"].rolling(30, center=True, min_periods=10).mean(), c=COL[n], label=n)
        ax[2].plot(S[n].index, S[n]["eta_mm_d"].rolling(30, center=True, min_periods=10).mean(), c=COL[n], label=n)
        ab = annual_budget(S[n]); ab.insert(0, "scenario", n); rows.append(ab.reset_index())
    ax[0].set_ylabel("LAI used (polygon mean)"); ax[0].legend(fontsize=7, ncol=3)
    ax[1].set_ylabel("ETp forcing (mm/d, 30-d mean)"); ax[2].set_ylabel("CATHY ETa (mm/d, 30-d mean)")
    _fmt(fig); save(fig, out, "scn_01_inputs_and_eta")
    tab = pd.concat(rows, ignore_index=True)
    tab.round(2).to_csv(out / "scenarios_annual_budget.csv", index=False)
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6), layout="constrained")
    for k, (c, lab) in enumerate((("eta_mm", "annual ETa (mm)"), ("eta_over_etp", "ETa / ETp"), ("gwl_mean_m", "mean water table (m)"))):
        piv = tab.pivot(index="year", columns="scenario", values=c)
        piv.plot.bar(ax=ax[k], width=0.85, color=[COL[n] for n in piv.columns], legend=(k == 0))
        ax[k].set_ylabel(lab)
    save(fig, out, "scn_02_annual_comparison")


def plot_diff_maps(DIRS, pairs, out, ref_years=(), gdf=None):
    """Maps of the LAI effect (obs - clim) on ETa and water-table depth: mean after the reference years."""
    out = Path(out)
    ref = list(ref_years)
    for no, nc in pairs:
        for fname, var, tdim, lab, cmap, name in (("eta_nodes_monthly.nc", "eta_mm_d", "month", "dETa (mm/d)", "RdBu_r", "eta"),
                                                   ("gwl_nodes.nc", "gwl_depth_m", "date", "d(WT depth) (m)", "RdBu_r", "gwl")):
            fo, fc = DIRS[no] / fname, DIRS[nc] / fname
            if not (fo.exists() and fc.exists()):
                continue
            with xr.open_dataset(fo) as a_, xr.open_dataset(fc) as b_:
                t = pd.DatetimeIndex(a_[tdim].values)
                common = t.intersection(pd.DatetimeIndex(b_[tdim].values))
                if not len(common):
                    continue
                late = common[common.year > max(ref)] if ref else common
                if not len(late):
                    late = common
                d = a_[var].sel({tdim: late}).mean(tdim).values - b_[var].sel({tdim: late}).mean(tdim).values
                x, y = a_["x"].values, a_["y"].values
                ins = a_["in_polygon"].values if "in_polygon" in a_ else np.ones(len(x), bool)
            X, Y, Z = node_grid(x, y, d)
            v = np.nanmax(np.abs(Z)) or 1.0
            fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), layout="constrained")
            _map(ax[0], X, Y, Z, gdf, f"{lab}: {no} - {nc}\n({late.min():%Y} - {late.max():%Y} mean)", cmap, -v, v, lab)
            ax[1].hist(d[ins], bins=30, color="tab:purple"); ax[1].axvline(0, c="gray")
            ax[1].set_xlabel(lab); ax[1].set_title(f"nodes in the polygon: mean {np.nanmean(d[ins]):+.3f}")
            save(fig, out, f"scn_map_{name}_{no}_minus_{nc}")


# =============================================================================================
# stand-alone
# =============================================================================================
if __name__ == "__main__":
    import argparse
    import matplotlib
    matplotlib.use("Agg")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", required=True, help="outputs folder of one scenario (contains timeseries_daily.csv)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    make_postrun_plots(a.outputs, a.out)
