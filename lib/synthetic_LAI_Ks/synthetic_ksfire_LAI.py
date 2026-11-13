"""
synthetic_ksfire_LAI.py — Synthetic 20x20 test of fire-driven Ks recovery
==========================================================================
GRwater project — Agramon (synthetic variant)

Minimal pyCATHY test case that reuses the Agramon pipeline building blocks
(Agramon_withLAI_withETp.py) but with:

  * a synthetic 20x20 planar DEM (no GeoTIFF / shapefile needed);
  * the SAME ATMBC forcing as the Agramon run (satellite rain - ETp,
    HSPATM=1) but spatially averaged -> one scalar series, no per-pixel ETp;
  * a synthetic monthly LAI field (seasonal cycle + fire on --ks-fire-date
    where LAI falls to 0, then logistic post-fire recovery);
  * PERMX (Ks) updated EVERY MONTH through ks_fire_utils from that month's
    domain-mean LAI (--ks-fire-mode lai, default), time since fire
    (time) or a static burn severity (burn_severity);
  * TWO-LAYER soil, the DEFAULT (--no-two-layers gives one homogeneous layer):
    a thin burned topsoil (--top-depth, default 0.1 m) whose PERMX follows
    ks_fire_utils every month, over a stable subsoil with a fixed PERMX
    (--permx-deep, default = unburned Ks). The layers are selected from the
    CATHY mesh strata: --top-depth (m) or --top-layers (count).
    Vertical heterogeneity goes through the (zone, layer) SPP_map of
    pyCATHY's update_soil, so no extra soil zones are needed;
  * CONTROL scenarios: --ks-fire-mode off keeps Ks constant before and after
    the fire (the vegetation/LAI still burns), in 1 or 2 layers.

Results are plotted with results_plotter / withLAI_results helpers, plus a
LAI vs Ks figure.

Usage
-----
    python synthetic_ksfire_LAI.py --smoke-test        # pure-python checks
    python synthetic_ksfire_LAI.py                     # lai mode, 2 layers (default)
    python synthetic_ksfire_LAI.py --no-two-layers     # lai mode, 1 homogeneous layer
    python synthetic_ksfire_LAI.py --top-depth 0.5 --permx-deep 5e-5
    python synthetic_ksfire_LAI.py --ks-fire-mode time --ks-fire-preset overall_average
    python synthetic_ksfire_LAI.py --ks-fire-mode off                # control, constant Ks, 2 layers
    python synthetic_ksfire_LAI.py --ks-fire-mode off --no-two-layers  # control, 1 layer
    python synthetic_ksfire_LAI.py --ks-fire-mode off --permx 1e-5   # control, user Ks
    python synthetic_ksfire_LAI.py --start-year 2020 --end-year 2023

Restart: every month that ran cleanly is backed up to
<path2prj>/_monthly_backup/<scenario>/. Re-running with the same settings
reloads those months and simulates only the missing ones (e.g. after a crash
or after extending --end-year). Any change of a simulation setting invalidates
the backup; --no-resume forces a full re-simulation.

Outputs are saved DURING the run, not only at the end: after every simulated
month (--save-every N, 0 = only at the end) the merged psi/sw/ET/recharge
NetCDF files, hgraph_output.csv, ks_fire_monthly.csv, failures.csv and the
area-mean CSVs are rewritten (atomically) in <path2prj>/<scenario>/outputs/,
so a crash or a Ctrl-C never loses the months already simulated. The figures
are made at the end of the run.

Console: for every month the script prints how ks_fire_utils turns that month's
LAI into a Ks and how the soil file is rewritten from it (driver, formula,
Ks before -> after, top/deep strata, read-back); a table of all Ks updates is
printed at the end of the monthly loop.
"""

import argparse
import hashlib
import json
import os
import pickle
import sys
from datetime import date, datetime
from pathlib import Path

import sys
sys.path.insert(0, "/home/ben/Nextcloud/BenCSIC/Codes/Tech4agro_org/GRwater_CATmentHYdrology/lib/")

import numpy as np
import pandas as pd

from utils import ks_fire_utils as KsFire

ET_VAR = "ACT. ETRA_patched"
ET_SCALE = 1e3 * 86400            # m/s -> mm/day
SYN_CRS = "EPSG:25830"
XLL, YLL = 450000.0, 4120000.0    # arbitrary UTM-like origin of the synthetic DEM
PATH_ET = Path.home() / "Nextcloud" / "GRwater" / "data" / "satellite" / "Agramon"

# ── Pure helpers (no pyCATHY needed -> usable in the smoke test) ─────────────

def synthetic_dem(n: int = 20, res: float = 10.0, seed: int = 0) -> np.ndarray:
    """Flat DEM at elevation 1 m; the (0, 0) corner is 1 - 1e-3 to fake an outlet cell.
    (res and seed kept for call compatibility.)"""
    z = np.ones((n, n), dtype=float)
    z[0, 0] = 1.0 - 1e-3
    return z


_KS_EQUATIONS = {
    "off": "static --permx (ks_fire_utils returns None)",
    "time": "Ks(t) = ks_ll + (ks_ul-ks_ll)*logistic(t_since_fire, r_ks)   [ks_recovery_time_based]",
    "lai": "Ks = ks_ll + (ks_ul-ks_ll)*(LAI-lai_ll)/(lai_ul-lai_ll), clipped [ks_from_lai_linear]",
    "burn_severity": "Ks = ks_unburned*(1 - max_reduction*severity)   [ks_burn_severity_step]",
}


def describe_ks_cfg(cfg, args) -> None:
    """Console summary of how ks_fire_utils is wired in and parametrised."""
    mm = lambda v: "n/a" if v is None else f"{KsFire.m_s_to_mm_hr(v):.2f} mm/hr ({v:.3e} m/s)"
    print("  ── ks_fire_utils usage ──")
    drv = ("LAI / unburned LAI of the same month (recovery ratio; LAI_UL=1)"
           if args.lai_ref == "ratio" else "raw domain-mean LAI (seasonal)")
    print(f"  per month: KsFire.compute_permx_for_month(cfg, date, lai_value=<{drv}>)")
    print(f"             -> AG.update_soil_for_active_zones(permx=...) "
          f"[PERMX; PERMY=x{args.permy_ratio:g}, PERMZ=x{args.permz_ratio:g}]")
    print(f"  mode       : {cfg.mode}  ->  {_KS_EQUATIONS.get(cfg.mode, '?')}")
    print(f"  fire_date  : {cfg.fire_date}")
    print(f"  ks_ll      : {mm(cfg.ks_ll)}")
    print(f"  ks_ul      : {mm(cfg.ks_ul)}")
    print(f"  r_ks       : {cfg.r_ks} yr^-1 (preset: {args.ks_fire_preset})")
    print(f"  lai_ll/ul  : {cfg.lai_ll} / {cfg.lai_ul}")
    print(f"  ks_unburned: {mm(cfg.ks_unburned)} | severity: {cfg.severity_fraction} "
          f"| max_reduction: {cfg.max_reduction}")
    if cfg.mode == "off":
        print(f"  static permx: {args.permx} (CONTROL: Ks constant before/after fire)")
    if args.two_layers:
        deep = resolve_deep_permx(args, cfg)
        print(f"  soil       : TWO layers -> top (burned, ks_fire_utils) / deep (stable) "
              f"{mm(deep)}")
    else:
        print("  soil       : ONE homogeneous layer (--no-two-layers)")
    print("  soil update, every month:")
    print("    LAI field -> vegetation zones (root depth / Feddes tables follow the zones)")
    if cfg.mode == "off":
        print("    no Ks override: PERMX stays at the static value in every update")
    else:
        print(f"    LAI-driver -> ks_fire_utils ({cfg.mode}) -> Ks -> update_soil: "
              f"PERMX = Ks, PERMY = x{args.permy_ratio:g}, PERMZ = x{args.permz_ratio:g}"
              + (" (top strata only; deep strata keep their constant Ks)"
                 if args.two_layers else " (all strata)"))
    print("    porosity / retention parameters are never modified")
    curve = ks_transfer_curve(cfg, args)
    if curve:
        print("  Ks transfer (ks_fire_utils output, mm/hr): "
              + " | ".join(f"{lab} -> {KsFire.m_s_to_mm_hr(k):.1f}" for lab, k in curve))


def soil_change_report(prev, new, rtol: float = 1e-6) -> dict:
    """Compare previous vs new PERMX (m/s). prev=None -> first iteration."""
    if new is None:
        return {"changed": False, "delta_pct": np.nan, "label": "no Ks override"}
    if prev is None:
        return {"changed": True, "delta_pct": np.nan, "label": "initial"}
    changed = not np.isclose(new, prev, rtol=rtol, atol=0.0)
    d = 100.0 * (new - prev) / prev if prev else np.nan
    return {"changed": bool(changed), "delta_pct": d,
            "label": "CHANGED" if changed else "unchanged"}


def ks_transfer_curve(cfg, args, n_points: int = 5) -> list:
    """[(driver label, Ks m/s)] returned by ks_fire_utils for a few driver values
    (lai mode: LAI_LL ... LAI_UL; time mode: 0-36 months after the fire).
    Display only; [] when not applicable or if ks_fire_utils refuses."""
    out = []
    try:
        fire = pd.Timestamp(args.ks_fire_date)
        if cfg.mode == "lai" and cfg.lai_ll is not None and cfg.lai_ul is not None:
            lab = "LAI/LAI_unburned" if args.lai_ref == "ratio" else "LAI"
            for x in np.linspace(cfg.lai_ll, cfg.lai_ul, n_points):
                out.append((f"{lab} {x:.2f}",
                            KsFire.compute_permx_for_month(cfg, fire, lai_value=float(x))))
        elif cfg.mode == "time":
            for m in (0, 3, 6, 12, 24, 36):
                out.append((f"+{m} mo",
                            KsFire.compute_permx_for_month(cfg, fire + pd.DateOffset(months=m))))
    except Exception:
        return []
    return [(lab, k) for lab, k in out if k is not None]


def ks_update_explanation(cfg, args, dt, lai_mean, lai_unb, driver, permx, prev_permx,
                          years_since_fire, permx_deep=None, layers=None, zones=None) -> list:
    """Console lines explaining ONE monthly soil update: what drives Ks, which
    ks_fire_utils equation gives it (with the numbers), how Ks changed since the
    previous month and which strata of the soil file receive it."""
    mm = lambda v: "n/a" if v is None else f"{KsFire.m_s_to_mm_hr(v):.2f}"
    rep_ = soil_change_report(prev_permx, permx)
    dtxt = "" if np.isnan(rep_["delta_pct"]) else f" ({rep_['delta_pct']:+.1f}%)"
    lines = []
    if cfg.mode == "off":
        lines.append(f"ks_fire_utils [off]: no override -> static PERMX {mm(permx)} mm/hr, "
                     f"identical every month (LAI {lai_mean:.2f} only changes the vegetation map)")
    else:
        rng = cfg.ks_ul - cfg.ks_ll
        frac = (permx - cfg.ks_ll) / rng if rng else np.nan
        if cfg.mode == "lai":
            if args.lai_ref == "ratio":
                drv = (f"LAI {lai_mean:.2f} / unburned LAI {lai_unb:.2f} "
                       f"= recovery ratio {driver:.3f}")
            else:
                drv = f"raw domain-mean LAI {driver:.2f}"
            lines.append(f"ks_fire_utils [lai]: driver = {drv}")
            span = cfg.lai_ul - cfg.lai_ll
            f_lai = 1.0 if not span else float(np.clip((driver - cfg.lai_ll) / span, 0.0, 1.0))
            expected = cfg.ks_ll + rng * f_lai
            lines.append(f"  Ks = Ks_LL + (Ks_UL-Ks_LL)*clip((x-LAI_LL)/(LAI_UL-LAI_LL)) = "
                         f"{mm(cfg.ks_ll)} + ({mm(cfg.ks_ul)}-{mm(cfg.ks_ll)}) * "
                         f"clip(({driver:.3f}-{cfg.lai_ll:.3f})/({cfg.lai_ul:.3f}-{cfg.lai_ll:.3f})) "
                         f"= {mm(expected)} mm/hr")
            if permx is not None and not np.isclose(expected, permx, rtol=1e-3):
                lines.append(f"  (value returned by ks_fire_utils: {mm(permx)} mm/hr)")
        elif cfg.mode == "time":
            lines.append(f"ks_fire_utils [time]: driver = {years_since_fire:+.2f} yr since the "
                         f"fire ({args.ks_fire_date})")
            lines.append(f"  Ks = Ks_LL + (Ks_UL-Ks_LL)*logistic(t, r_ks={cfg.r_ks}) -> "
                         f"fraction {frac:.3f} -> {mm(permx)} mm/hr")
        else:
            lines.append(f"ks_fire_utils [{cfg.mode}]: driver = burn severity "
                         f"{cfg.severity_fraction}")
            lines.append(f"  Ks = Ks_unburned*(1 - max_reduction*severity) = {mm(cfg.ks_unburned)}"
                         f"*(1 - {cfg.max_reduction}*{cfg.severity_fraction}) -> {mm(permx)} mm/hr")
        if np.isfinite(frac):
            lines.append(f"  recovery: {100 * frac:.0f}% of the way from burned Ks_LL "
                         f"({mm(cfg.ks_ll)}) to unburned Ks_UL ({mm(cfg.ks_ul)})")
    prev_txt = "initial value" if prev_permx is None else f"{mm(prev_permx)} -> {mm(permx)} mm/hr"
    lines.append(f"soil Ks vs previous month: {prev_txt}{dtxt}  [{rep_['label']}]")
    zt = f" | veg zones {zones}" if zones is not None else ""
    ratios = f"PERMY = x{args.permy_ratio:g}, PERMZ = x{args.permz_ratio:g}"
    if args.two_layers and layers is not None:
        lines.append(f"soil file (update_soil / SPP_map): strata 1-{layers['n_top']} "
                     f"(0-{layers['z_top']:.2f} m, burned) PERMX = {mm(permx)} mm/hr "
                     f"[{rep_['label']}]{zt}")
        lines.append(f"                                   strata {layers['n_top'] + 1}-{layers['nstr']} "
                     f"({layers['z_deep0']:.2f}-{layers['z_bot']:.2f} m, stable) "
                     f"PERMX = {mm(permx_deep)} mm/hr [constant]; {ratios}")
    else:
        lines.append(f"soil file (update_soil / SPP_map): PERMX = {mm(permx)} mm/hr in every "
                     f"stratum [{rep_['label']}]; {ratios}{zt}")
    return lines


def print_ks_update_summary(df, two_layers: bool = False, fire_date=None) -> None:
    """Table of every monthly soil update (what ks_fire_utils did over time)."""
    if df is None or len(df) == 0:
        return
    k = df["permx_mm_hr"].to_numpy(float)
    prev = np.r_[np.nan, k[:-1]]
    with np.errstate(divide="ignore", invalid="ignore"):
        dpc = 100.0 * (k - prev) / prev
    changed = np.isfinite(prev) & ~np.isclose(k, prev, rtol=1e-6)
    fire = None if fire_date is None else pd.Timestamp(fire_date)
    print("  ── ks_fire_utils: soil Ks updates over time "
          f"({'top layer' if two_layers else 'whole soil'}) ──")
    print("  month    LAI   driver  Ks mm/hr   change   soil file")
    for i, r in enumerate(df.itertuples(index=False)):
        m = pd.Timestamp(r.month)
        state = "initial" if i == 0 else ("CHANGED" if changed[i] else "unchanged")
        d = "" if not np.isfinite(dpc[i]) else f"{dpc[i]:+7.1f}%"
        mark = "  <- fire" if fire is not None and m.to_period("M") == fire.to_period("M") else ""
        print(f"  {m:%Y-%m}  {r.lai_mean:5.2f}  {r.lai_driver:6.3f}  {k[i]:8.2f}  {d:>8}   "
              f"{state}{mark}")
    n = int(changed.sum())
    print(f"  Ks changed in {n} of {max(len(df) - 1, 0)} monthly updates "
          f"(range {np.nanmin(k):.2f} - {np.nanmax(k):.2f} mm/hr)")
    if two_layers and "permx_deep_mm_hr" in df and df["permx_deep_mm_hr"].notna().any():
        print(f"  deep layer: constant {df['permx_deep_mm_hr'].dropna().iloc[0]:.2f} mm/hr "
              f"in every update")


def assemble_outputs(psi_rec, sw_rec, et_rec, rec_rec, hg_rec) -> dict:
    """Merge the per-month records into the run datasets (psi / sw need a
    (datetime, node) index; duplicated times at month boundaries keep the last)."""
    import xarray as xr
    key = ["datetime", "node"]
    return {
        "psi": (pd.concat(psi_rec, ignore_index=True).groupby(key, as_index=False)
                .last().set_index(key).to_xarray()),
        "sw": (pd.concat(sw_rec, ignore_index=True).groupby(key, as_index=False)
               .last().set_index(key).to_xarray()),
        "et": xr.concat(et_rec, dim="datetime"),
        "rech": xr.concat(rec_rec, dim="datetime") if rec_rec else None,
        "hg": pd.concat(hg_rec, ignore_index=True) if hg_rec else None,
    }


def _atomic_write(path: Path, writer) -> None:
    """writer(tmp_path) then rename: an interrupted save never corrupts the
    previous (complete) file."""
    path = Path(path)
    tmp = path.with_name(f"{path.stem}.tmp{path.suffix}")
    writer(tmp)
    os.replace(tmp, path)


def write_outputs(out_dir, asm: dict) -> None:
    """psi/sw/et/recharge NetCDF + hgraph CSV of the months simulated so far."""
    out_dir = Path(out_dir)
    for name, key in (("psi_output.nc", "psi"), ("sw_output.nc", "sw"),
                      ("et_output.nc", "et"), ("recharge_output.nc", "rech")):
        if asm.get(key) is not None:
            _atomic_write(out_dir / name, asm[key].to_netcdf)
    if asm.get("hg") is not None:
        _atomic_write(out_dir / "hgraph_output.csv",
                      lambda p: asm["hg"].to_csv(p, index=False))


def _get_spp_map(simu):
    """SPP_map DataFrame (index: zone, layer) held by the simulation, or None."""
    for getter in (lambda: simu.soil_SPP["SPP_map"],
                   lambda: simu.soil_SPP[0] if isinstance(simu.soil_SPP, (list, tuple)) else None,
                   lambda: simu.read_inputs("soil")["SPP_map"]):
        try:
            df = getter()
            if df is not None and "PERMX" in df:
                return df
        except Exception:
            continue
    return None


def read_back_permx(simu):
    """PERMX actually held by the simulation object after update_soil (m/s),
    or None if it cannot be read (pyCATHY attribute layout varies)."""
    df = _get_spp_map(simu)
    if df is None:
        return None
    v = np.asarray(df["PERMX"], dtype=float)
    return float(np.nanmean(v)), float(np.nanmin(v)), float(np.nanmax(v))


def read_back_permx_layers(simu, n_top):
    """Mean PERMX (m/s) of the top n_top strata and of the deeper ones, as
    held by the simulation after update_soil. None if it cannot be read."""
    df = _get_spp_map(simu)
    if df is None or "layer" not in df.index.names:
        return None
    layer = np.asarray(df.index.get_level_values("layer"))
    v = np.asarray(df["PERMX"], dtype=float)
    top, deep = v[layer <= n_top], v[layer > n_top]
    if top.size == 0 or deep.size == 0:
        return None
    return float(np.nanmean(top)), float(np.nanmean(deep))


# ── Two-layer soil helpers (pure) ─────────────────────────────────────────────

def n_top_layers(layer_tops, layer_bots, top_depth=0.5, top_layers=None) -> int:
    """Number of mesh strata that belong to the burned topsoil.

    top_layers (count) wins when given; otherwise a stratum is 'top' when its
    mid-depth is <= top_depth (m). At least 1 stratum is always top, and at
    least 1 must stay below (ValueError otherwise)."""
    n = len(layer_tops)
    if top_layers is not None:
        k = int(top_layers)
    else:
        mids = 0.5 * (np.asarray(layer_tops, float) + np.asarray(layer_bots, float))
        k = max(int(np.sum(mids <= float(top_depth))), 1)
    if not 1 <= k < n:
        raise ValueError(f"Two-layer soil needs 1 <= n_top < nstr (got n_top={k}, nstr={n}); "
                         f"adjust --top-depth / --top-layers or the mesh strata.")
    return k


def two_layer_spp(spp_map, n_top, permx_top, permx_deep,
                  permy_ratio: float = 1.0, permz_ratio: float = 1.0):
    """Copy of an SPP_map (index levels 'zone','layer') with PERMX/PERMY/PERMZ
    set to the burned-topsoil Ks for strata <= n_top and to the stable
    subsoil Ks below. All other soil properties are left untouched."""
    df = spp_map.copy()
    layer = np.asarray(df.index.get_level_values("layer"))
    is_top = layer <= n_top
    for col, ratio in (("PERMX", 1.0), ("PERMY", permy_ratio), ("PERMZ", permz_ratio)):
        df[col] = np.where(is_top, permx_top * ratio, permx_deep * ratio).astype(float)
    return df


def ic_settings(args) -> dict:
    """kwargs of simu.update_ic: hydrostatic (INDP=3) when --ic-wtposition is given,
    else the uniform pressure head --ic-head (INDP=0, default -15 m)."""
    if args.ic_wtposition is not None:
        return {"INDP": 3, "WTPOSITION": args.ic_wtposition}
    return {"INDP": 0, "pressure_head_ini": float(args.ic_head)}


def resolve_deep_permx(args, cfg):
    """Ks of the stable subsoil: --permx-deep, else the unburned Ks (ks_ul),
    else (control, mode off) the same constant Ks as the topsoil."""
    if args.permx_deep is not None:
        return args.permx_deep
    if cfg.mode == "off":
        return args.permx
    return cfg.ks_ul



def synthetic_lai(dt, shape=(20, 20), fire_date=None, rate: float = 0.8) -> np.ndarray:
    """Monthly LAI field in [0, ~3] m2/m2.

    Pre-fire: 1.3 + 1.4*season (peak in May) times a mild W->E factor.
    Fire: LAI * recov(t), i.e. LAI drops to 0 everywhere at the fire date
    (total canopy loss) and recovers toward the unburned field with
    recov(t) = 0 at the fire, -> 1 logistically (rate in yr^-1).
    fire_date=None gives the unburned reference field.
    """
    ny, nx = shape
    month = pd.Timestamp(dt).month
    season = 0.5 * (1.0 + np.sin(2 * np.pi * (month - 2) / 12.0))
    xfac = np.tile(np.linspace(1.1, 0.9, nx), (ny, 1))
    lai = (1.3 + 1.4 * season) * xfac
    if fire_date is not None:
        t = KsFire.years_since_fire(pd.Timestamp(dt).to_pydatetime(), fire_date)
        if t >= 0:
            recov = 2.0 * (KsFire.logistic_recovery(t, 0.0, 1.0, rate) - 0.5)
            lai = lai * recov          # severity 1: LAI = 0 at the fire
    return np.clip(lai, 0.0, None)


def lai_recovery_ratio(lai_post, lai_unburned) -> float:
    """Domain-mean LAI relative to the unburned reference of the SAME calendar
    month. Seasonality cancels: 1 = fully recovered/unburned, 0 = no canopy."""
    ref = float(np.nanmean(lai_unburned))
    return float(np.nanmean(lai_post)) / ref if ref > 0 else 1.0


def prefire_monthly_climatology(lai_by_date: dict, fire_date) -> dict:
    """Calendar month (1-12) -> mean PRE-fire LAI field, from an observed LAI
    series {date: 2-D array}. Use it as the 'unburned reference' for real data
    (Thomas et al. 2021 compare post-fire LAI with the long-term pre-fire
    monthly average). Months with no pre-fire observation are absent."""
    fire = pd.Timestamp(fire_date)
    acc: dict = {}
    for d, f in lai_by_date.items():
        if pd.Timestamp(d) < fire:
            acc.setdefault(pd.Timestamp(d).month, []).append(np.asarray(f, float))
    return {m: np.nanmean(np.stack(v), axis=0) for m, v in acc.items()}


def lai_driver_values(lai_ref: str, months, lai_all: dict, lai_unburned: dict) -> dict:
    """The LAI value that feeds ks_fire_utils each month: the recovery ratio
    (default, 'ratio') or the raw domain-mean LAI ('raw', seasonal)."""
    if lai_ref == "ratio":
        return {m: lai_recovery_ratio(lai_all[m], lai_unburned[m]) for m in months}
    return {m: float(np.nanmean(lai_all[m])) for m in months}


def finalize_ks_cfg(cfg: KsFire.KsFireConfig, lai_means: np.ndarray,
                    lai_ref: str = "raw") -> KsFire.KsFireConfig:
    """Fill unset anchors with Thomas et al. (2021) satellite anchors and
    the LAI limits, so every mode runs with no extra flags.

    lai_ref='ratio': lai_means are recovery ratios; LAI_UL = 1 (unburned) and
    LAI_LL = lowest post-fire ratio (Thomas et al. 2021: LAI_LL = LAI right
    after the fire, LAI_UL = unburned LAI).
    lai_ref='raw'  : limits = min / max of the (seasonal) LAI means."""
    if cfg.mode == "off":
        return cfg
    anc = KsFire.THOMAS2021_SATELLITE_KS_ANCHORS_MM_HR
    if cfg.ks_ll is None:
        cfg.ks_ll = KsFire.mm_hr_to_m_s(anc["ks_ll"])
    if cfg.ks_ul is None:
        cfg.ks_ul = KsFire.mm_hr_to_m_s(anc["ks_ul"])
    if cfg.mode == "time":
        pre = KsFire.get_field_preset_m_s("overall_average")
        cfg.r_ks = pre["r"] if cfg.r_ks is None else cfg.r_ks
    if cfg.mode == "lai":
        cfg.lai_ll = float(np.min(lai_means)) if cfg.lai_ll is None else cfg.lai_ll
        if lai_ref == "ratio":
            cfg.lai_ul = 1.0 if cfg.lai_ul is None else cfg.lai_ul
        else:
            cfg.lai_ul = float(np.max(lai_means)) if cfg.lai_ul is None else cfg.lai_ul
    if cfg.mode == "burn_severity":
        cfg.ks_unburned = cfg.ks_ul
        cfg.severity_fraction = 0.6 if cfg.severity_fraction is None else cfg.severity_fraction
    return cfg


def resolve_static_permx(args, cfg) -> float:
    """Constant Ks of the CONTROL (--ks-fire-mode off): --permx if given, else
    the unburned Ks anchor, so the control is directly comparable to the
    fire runs (whose Ks recovers toward that same upper limit)."""
    if args.permx is not None:
        return args.permx
    return KsFire.mm_hr_to_m_s(KsFire.THOMAS2021_SATELLITE_KS_ANCHORS_MM_HR["ks_ul"])


def planned_ks_series(cfg, args, months, lai_means, lai_driver, permx_deep=None) -> pd.DataFrame:
    """Ks (top, and deep if two layers) the monthly loop WILL apply, computed
    from the config and the synthetic LAI only -- no simulation needed.
    lai_means: domain-mean LAI; lai_driver: the value fed to ks_fire_utils
    (recovery ratio or raw LAI, see --lai-ref)."""
    rows = []
    for dt, lai, drv in zip(months, lai_means, lai_driver):
        k = KsFire.compute_permx_for_month(cfg, dt, lai_value=float(drv))
        k = k if k is not None else args.permx
        rows.append({"month": dt, "lai_mean": float(lai), "lai_driver": float(drv),
                     "permx_top_mm_hr": KsFire.m_s_to_mm_hr(k) if k else np.nan,
                     "permx_deep_mm_hr": (KsFire.m_s_to_mm_hr(permx_deep)
                                          if permx_deep else np.nan)})
    return pd.DataFrame(rows)


def plot_pre_simulation(plt, plotter, fig_dir, df_plan, lai_all, etp_mmd, rain_mmd,
                        fire_date, cfg, two_layers, dpi=150, lai_ref="raw"):
    """Input trends shown BEFORE the run: LAI, ETp (+ rain) and Ks vs time."""
    months = list(df_plan["month"])
    lai_min = [lai_all[m].min() for m in months]
    lai_max = [lai_all[m].max() for m in months]
    t0, t1 = min(months), max(months) + pd.offsets.MonthEnd(0)

    fig, (a0, a1, a2) = plt.subplots(3, 1, figsize=(11, 9), sharex=True,
                                     constrained_layout=True)
    a0.fill_between(months, lai_min, lai_max, color="#4daf4a", alpha=0.25,
                    label="LAI min–max")
    a0.plot(months, df_plan["lai_mean"], "s-", color="#006837", label="LAI domain mean")
    a0.set_ylabel("LAI (m² m⁻²)")
    if lai_ref == "ratio":
        a0b = a0.twinx()
        a0b.plot(months, df_plan["lai_driver"], "d--", color="#b35806", ms=4,
                 label="LAI / unburned LAI (drives Ks)")
        a0b.set_ylabel("recovery ratio (–)")
        a0b.set_ylim(-0.05, 1.1)
        hh, ll = a0.get_legend_handles_labels()
        h2, l2 = a0b.get_legend_handles_labels()
        a0.legend(hh + h2, ll + l2, loc="lower right", fontsize=8)
    a0.set_title("Inputs before simulation: LAI, ETp / rain, Ks", fontsize=10,
                 fontweight="bold")

    a1.plot(etp_mmd.index, etp_mmd.values, color="#fdae61", lw=0.7, alpha=0.8,
            label="ETp daily")
    etp_m = etp_mmd.resample("MS").mean()
    a1.plot(etp_m.index, etp_m.values, "o-", color="#d94801", ms=3, label="ETp monthly mean")
    a1.set_ylabel("ETp (mm d⁻¹)")
    a1b = a1.twinx()
    rain_m = rain_mmd.resample("MS").sum()
    a1b.bar(rain_m.index, rain_m.values, width=20, color="#4575b4", alpha=0.35,
            label="Rain monthly total")
    a1b.set_ylabel("Rain (mm month⁻¹)")
    h1, l1 = a1.get_legend_handles_labels()
    h2, l2 = a1b.get_legend_handles_labels()
    a1.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    ks_lbl = "PERMX top (burned)" if two_layers else "PERMX"
    a2.plot(months, df_plan["permx_top_mm_hr"], "o-", color="#762a83", label=ks_lbl)
    if two_layers:
        a2.plot(months, df_plan["permx_deep_mm_hr"], "s--", color="#1b7837",
                label="PERMX deep (stable)")
    a2.set_ylabel("PERMX (mm/hr)")
    a2.set_title(f"Ks (mode={cfg.mode}{', control: constant Ks' if cfg.mode == 'off' else ''})",
                 fontsize=9)

    for a in (a0, a1, a2):
        plotter.add_fire_vline(a, t0, t1)
        a.grid(True, linestyle="--", alpha=0.3)
    if lai_ref != "ratio":
        a0.legend(loc="upper right", fontsize=8)
    a2.legend(loc="upper right", fontsize=8)
    plotter.save_fig(fig, fig_dir, "fig_inputs_lai_etp_ks.png", dpi)
    plt.close(fig)


# ── Monthly backup / restart (pure file helpers) ──────────────────────────────
# After every month that ran cleanly, its results and the end-of-month pressure
# head (initial condition of the next month) are pickled to
# <path2prj>/_monthly_backup/<scenario>/. A later run with the SAME settings
# reloads those months instead of re-running CATHY and continues from the first
# missing one. Backups are local files written by this script: only reload
# folders you produced yourself (pickle).

BACKUP_VERSION = 1
_FP_IGNORE = {"path2prj", "dpi", "no_plots", "smoke_test", "no_resume", "save_every", "debug"}


def run_fingerprint(args, months) -> str:
    """Hash of everything that changes the simulated numbers (not plotting /
    paths). A different fingerprint invalidates the backup."""
    d = {k: v for k, v in vars(args).items() if k not in _FP_IGNORE}
    files = {}
    for f in (args.etp_file, args.rain_file):
        try:
            st = os.stat(f)
            files[str(f)] = [st.st_size, int(st.st_mtime)]
        except OSError:
            files[str(f)] = None
    payload = {"v": BACKUP_VERSION, "args": d, "files": files,
               "months": [str(pd.Timestamp(m).date()) for m in months]}
    return hashlib.sha1(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def _month_file(bdir: Path, dt) -> Path:
    return Path(bdir) / f"state_{pd.Timestamp(dt):%Y-%m}.pkl"


def prepare_backup(bdir: Path, fingerprint: str, resume: bool = True) -> bool:
    """Create the backup folder and decide whether existing months are usable.
    True -> reload them; False -> start from scratch (stale month files of a
    run with other settings are removed, the meta is rewritten)."""
    bdir = Path(bdir)
    bdir.mkdir(parents=True, exist_ok=True)
    meta = bdir / "meta.json"
    old = None
    if meta.exists():
        try:
            old = json.loads(meta.read_text()).get("fingerprint")
        except Exception:
            old = None
    existing = sorted(bdir.glob("state_*.pkl"))
    if resume and old == fingerprint:
        print(f"  Backup: {len(existing)} month(s) found in {bdir} -> will be reused")
        return True
    if existing:
        why = "--no-resume" if not resume else "settings changed since the backup"
        print(f"  Backup: ignoring {len(existing)} existing month(s) ({why}); "
              f"they will be overwritten")
        for f in existing:
            f.unlink()
    meta.write_text(json.dumps({"fingerprint": fingerprint, "version": BACKUP_VERSION}))
    return False


def save_month_backup(bdir: Path, dt, bundle: dict) -> None:
    """Atomic write (tmp file + rename) so an interrupted run never leaves a
    half-written month behind."""
    f = _month_file(bdir, dt)
    tmp = f.with_suffix(".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(bundle, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, f)


def load_month_backup(bdir: Path, dt):
    """Bundle of that month, or None if missing / unreadable."""
    f = _month_file(bdir, dt)
    if not f.exists():
        return None
    try:
        with open(f, "rb") as fh:
            return pickle.load(fh)
    except Exception as exc:
        print(f"  WARNING: unreadable backup {f.name}: {exc}")
        return None


def load_forcing_series(start_year, end_year, etp_path, rain_path):
    """Domain-mean daily ETp and rain (mm/day) over the whole satellite
    raster extent (no spatial detail kept), restricted to the run window."""
    import xarray as xr

    def _mean(path, var):
        ds = xr.open_dataset(path, decode_times=True, mask_and_scale=True)
        da = ds[var]
        sp = [d for d in da.dims if d != "time"]
        s = da.mean(dim=sp, skipna=True).to_series()
        s.index = pd.DatetimeIndex(s.index)
        return s.sort_index()

    etp = _mean(etp_path, "ET_0-gf")
    rain = _mean(rain_path, "TP-DD")
    idx = etp.index.intersection(rain.index)
    idx = idx[(idx >= f"{start_year}-01-01") & (idx <= f"{end_year - 1}-12-31")]
    if len(idx) == 0:
        raise ValueError(f"No common ETp/rain days in {start_year}-{end_year - 1}.")
    etp, rain = etp.loc[idx], rain.loc[idx]
    if etp.isna().any() or rain.isna().any():
        raise ValueError("NaN in domain-mean ETp/rain series.")
    return etp, rain


def resolve_nodes_simple(grid3d, plotter, depths=(0, 1)):
    """Outlet (lowest surface node) + mid-elevation band, at each depth."""
    from scipy.spatial import cKDTree
    nodes = grid3d["mesh3d_nodes"]
    surf = nodes[:int(grid3d["nnod"])]
    tree = cKDTree(nodes)
    o = int(surf[:, 2].argmin())
    ox, oy, oz = surf[o]
    band = plotter.find_representative_mid_node(surf, return_band=True)
    out = {}
    for d in depths:
        _, oi = tree.query([ox, oy, oz - d])
        out[f"outlet_z-{d}m"] = int(oi)
        pts = np.column_stack([surf[band, 0], surf[band, 1], surf[band, 2] - d])
        _, bi = tree.query(pts)
        out[f"mid_z-{d}m"] = [int(i) for i in bi]
    return out


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Synthetic 20x20 pyCATHY test: monthly LAI-driven fire Ks recovery",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--start-year", type=int, default=2020)
    p.add_argument("--end-year", type=int, default=2022, help="exclusive")
    p.add_argument("--n", type=int, default=20, help="synthetic DEM size (n x n)")
    p.add_argument("--res", type=float, default=1.0, help="DEM pixel size (m); 20x20 -> 20 m x 20 m")
    p.add_argument("--lai-recovery-rate", type=float, default=0.8, metavar="YR^-1",
                   help="rate of the synthetic post-fire LAI recovery")
    p.add_argument("--lai-ref", choices=["ratio", "raw"], default="ratio",
                   help="LAI fed to the 'lai' Ks mode. ratio: LAI relative to the unburned/"
                        "pre-fire LAI of the same calendar month (seasonality cancels; "
                        "LAI_UL=1 unburned, LAI_LL=post-fire minimum; Thomas et al. 2021). "
                        "raw: absolute seasonal LAI (Ks then follows the season)")
    p.add_argument("--zroot", type=float, default=1.0)
    p.add_argument("--pmin", type=float, default=-1e35)
    p.add_argument("--scf", type=float, default=1.0)
    p.add_argument("--permx", type=float, default=None, metavar="M/S",
                   help="static PERMX, used only when --ks-fire-mode off "
                        "(default there: the unburned Ks anchor ks_ul)")
    p.add_argument("--two-layers", action=argparse.BooleanOptionalAction, default=True,
                   help="two-layer soil (DEFAULT): burned topsoil (Ks from ks_fire_utils, "
                        "monthly) over a stable subsoil (--permx-deep). --no-two-layers: "
                        "one homogeneous layer")
    p.add_argument("--top-depth", type=float, default=0.1, metavar="M",
                   help="thickness of the burned topsoil (m); strata with mid-depth "
                        "<= this are 'top'. Ignored if --top-layers is given")
    p.add_argument("--top-layers", type=int, default=None, metavar="N",
                   help="number of mesh strata in the burned topsoil (overrides --top-depth)")
    p.add_argument("--permx-deep", type=float, default=None, metavar="M/S",
                   help="constant PERMX of the stable subsoil "
                        "(default: unburned Ks ks_ul; in mode off: same as --permx)")
    p.add_argument("--permy-ratio", type=float, default=1.0)
    p.add_argument("--permz-ratio", type=float, default=1.0)
    p.add_argument("--outlet-flux", type=float, default=0.0, metavar="M/S")
    p.add_argument("--bottom-flux", type=float, default=0.0, metavar="M/S")
    p.add_argument("--ic-head", type=float, default=-15.0, metavar="M",
                   help="DEFAULT initial condition: uniform (fixed) pressure head (m) at "
                        "all nodes (INDP=0)")
    p.add_argument("--ic-wtposition", type=float, default=None, metavar="M",
                   help="if given, overrides --ic-head: hydrostatic IC with the water table "
                        "at this height above the base of the grid (INDP=3)")
    p.add_argument("--eta-net-neg-mode", choices=["etp", "gross", "off"], default="etp")
    p.add_argument("--etp-file", default=str(PATH_ET / "20161001_20241231_ET0.nc"))
    p.add_argument("--rain-file", default=str(PATH_ET / "20161001_20241231_TP.nc"))
    p.add_argument("--rain-scale", type=float, default=1.0, metavar="X",
                   help="multiply the daily rain forcing by X (wetter/drier test; 1 = as read)")
    p.add_argument("--etp-scale", type=float, default=1.0, metavar="X",
                   help="multiply the daily ETp forcing by X (1 = as read)")
    p.add_argument("--dtmin", type=float, default=1e-2)
    p.add_argument("--dtmax", type=float, default=1e4)
    p.add_argument("--deltat", type=float, default=1e3)
    p.add_argument("--max-consecutive-fail", type=int, default=3, metavar="N",
                   help="stop the monthly loop after N months in a row that could not be "
                        "simulated (the months done so far are still saved); 0 = never stop")
    p.add_argument("--path2prj", default="../synthetic_LAI_Ks",
                   help="root folder of the run (project, outputs/ and figures/ "
                        "are created under <path2prj>/<scenario>/)")
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--debug", action="store_true",
                   help="stop the monthly loop at the FIRST month the solver cannot "
                        "complete (outputs of the earlier months are still saved)")
    p.add_argument("--no-resume", action="store_true",
                   help="ignore the monthly backup and re-simulate every month "
                        "(the backup is then overwritten). By default months already "
                        "simulated with the SAME settings are reloaded, not re-run")
    p.add_argument("--save-every", type=int, default=1, metavar="N",
                   help="rewrite the merged outputs (psi/sw/ET/recharge NetCDF, hgraph, "
                        "area-mean CSVs) after every N simulated months; ks_fire_monthly.csv "
                        "and failures.csv are saved after every month. 0 = only at the end")
    p.add_argument("--smoke-test", action="store_true",
                   help="run assertion-based checks of the pure helpers and exit")
    KsFire.add_ks_fire_args(p)
    p.set_defaults(ks_fire_mode="lai", ks_fire_date="2020-07-01")
    return p.parse_args(argv)


# ── Smoke test ────────────────────────────────────────────────────────────────

def smoke_test() -> None:
    args = _parse_args([])
    fd = datetime.strptime(args.ks_fire_date, "%Y-%m-%d").date()

    dem = synthetic_dem(20, 10.0)
    assert dem.shape == (20, 20) and np.isfinite(dem).all()
    assert dem[0, 0] == 1.0 - 1e-3 and np.isclose(dem.ravel()[1:], 1.0).all(), "flat DEM + outlet corner"
    assert dem.argmin() == 0, "outlet must be the unique lowest cell"

    months = pd.date_range("2020-01-01", "2022-12-01", freq="MS")
    burn = np.array([synthetic_lai(m, fire_date=fd).mean() for m in months])
    ref = np.array([synthetic_lai(m).mean() for m in months])
    assert synthetic_lai(months[0]).shape == (20, 20)
    assert np.all(burn <= ref + 1e-12), "burn must never raise LAI"
    i_f, i_late = list(months).index(pd.Timestamp("2020-08-01")), list(months).index(pd.Timestamp("2022-08-01"))
    assert burn[i_f] / ref[i_f] < burn[i_late] / ref[i_late] < 1.0 + 1e-9, "LAI must recover"
    assert np.allclose(burn[:6], ref[:6]), "pre-fire months unaffected"
    lai_fire = synthetic_lai(pd.Timestamp("2020-07-01"), fire_date=fd)
    assert np.allclose(lai_fire, 0.0), "LAI must fall to 0 at the fire date"
    assert synthetic_lai(pd.Timestamp("2020-07-01")).min() > 0.5, "unburned LAI stays > 0"

    cfg = finalize_ks_cfg(KsFire.ks_fire_config_from_args(args), burn)
    ks = np.array([KsFire.compute_permx_for_month(cfg, m, lai_value=l)
                   for m, l in zip(months, burn)])
    assert np.all(ks >= cfg.ks_ll - 1e-18) and np.all(ks <= cfg.ks_ul + 1e-18)
    assert ks[i_f] < ks[i_late], "Ks must recover with LAI"
    assert abs(KsFire.m_s_to_mm_hr(KsFire.mm_hr_to_m_s(16.2)) - 16.2) < 1e-9

    args_t = _parse_args(["--ks-fire-mode", "time"])
    cfg_t = finalize_ks_cfg(KsFire.ks_fire_config_from_args(args_t), burn)
    k_t = [KsFire.compute_permx_for_month(cfg_t, m) for m in months]
    assert k_t[0] == cfg_t.ks_ll and k_t[-1] > k_t[i_f]
    r0 = soil_change_report(None, ks[0]); r1 = soil_change_report(ks[0], ks[0])
    r2 = soil_change_report(ks[0], ks[-1])
    assert r0["label"] == "initial" and not r1["changed"] and r2["changed"]
    assert soil_change_report(1.0, None)["label"] == "no Ks override"
    assert np.isclose(r2["delta_pct"], 100 * (ks[-1] - ks[0]) / ks[0])
    assert np.sum(np.abs(np.diff(ks)) > 0) >= 3, "Ks must vary month to month in lai mode"

    # ── two-layer soil + control ──
    tops, bots = [0.0, 0.1, 0.3, 0.7, 1.5], [0.1, 0.3, 0.7, 1.5, 3.0]
    assert n_top_layers(tops, bots, top_depth=0.5) == 3          # mids .05 .2 .5 -> 3 strata
    assert n_top_layers(tops, bots, top_depth=0.01) == 1         # always >= 1 top stratum
    assert n_top_layers(tops, bots, top_layers=2) == 2           # count overrides depth
    for bad in (dict(top_layers=0), dict(top_layers=5), dict(top_depth=100.0)):
        try:
            n_top_layers(tops, bots, **bad)
            raise AssertionError(f"n_top_layers should reject {bad}")
        except ValueError:
            pass
    idx = pd.MultiIndex.from_product([[1], [1, 2, 3, 4, 5]], names=["zone", "layer"])
    spp0 = pd.DataFrame({"PERMX": 1.0, "PERMY": 1.0, "PERMZ": 1.0, "POROS": 0.5}, index=idx)
    spp2 = two_layer_spp(spp0, 2, 4e-6, 1e-5, permy_ratio=1.0, permz_ratio=0.5)
    assert np.allclose(spp2["PERMX"].values, [4e-6, 4e-6, 1e-5, 1e-5, 1e-5])
    assert np.allclose(spp2["PERMZ"].values, [2e-6, 2e-6, 5e-6, 5e-6, 5e-6])
    assert np.allclose(spp0["PERMX"].values, 1.0) and (spp2["POROS"] == 0.5).all(), \
        "input must not be modified, other properties untouched"
    # deep Ks must stay constant while the top Ks follows the fire in lai mode
    top_series = [two_layer_spp(spp0, 2, k, cfg.ks_ul).loc[(1, 5), "PERMX"] for k in ks]
    assert np.allclose(top_series, cfg.ks_ul), "deep Ks must not follow the fire"
    assert resolve_deep_permx(args, cfg) == cfg.ks_ul
    args_d = _parse_args(["--permx-deep", "2e-5"])
    assert resolve_deep_permx(args_d, cfg) == 2e-5
    args_c = _parse_args(["--ks-fire-mode", "off", "--two-layers"])
    cfg_c = finalize_ks_cfg(KsFire.ks_fire_config_from_args(args_c), burn)
    args_c.permx = resolve_static_permx(args_c, cfg_c)
    assert np.isclose(args_c.permx, KsFire.mm_hr_to_m_s(
        KsFire.THOMAS2021_SATELLITE_KS_ANCHORS_MM_HR["ks_ul"])), "control Ks defaults to ks_ul"
    assert resolve_deep_permx(args_c, cfg_c) == args_c.permx, "control: deep == top"
    ctrl = [KsFire.compute_permx_for_month(cfg_c, m, lai_value=l) for m, l in zip(months, burn)]
    assert all(c is None for c in ctrl), "mode off must not override Ks (constant control)"
    args_u = _parse_args(["--ks-fire-mode", "off", "--permx", "1e-5"])
    assert resolve_static_permx(args_u, cfg_c) == 1e-5

    # ── LAI relative to the unburned state (default --lai-ref ratio) ──
    ml = list(months)
    unb = {m: synthetic_lai(m) for m in ml}
    bur = {m: synthetic_lai(m, fire_date=fd) for m in ml}
    drv = lai_driver_values("ratio", ml, bur, unb)
    r = np.array([drv[m] for m in ml])
    i_fire = ml.index(pd.Timestamp("2020-07-01"))
    assert np.allclose(r[:i_fire], 1.0), "ratio = 1 before the fire (seasonality cancels)"
    assert np.isclose(r[i_fire], 0.0), "ratio = 0 at the fire"
    assert r.max() <= 1.0 + 1e-9 and r[i_late] > r[i_f], "ratio recovers toward 1"
    assert args.lai_ref == "ratio", "ratio must be the default"
    cfg_r = finalize_ks_cfg(KsFire.ks_fire_config_from_args(args), r, "ratio")
    assert cfg_r.lai_ul == 1.0 and np.isclose(cfg_r.lai_ll, 0.0)
    ks_r = np.array([KsFire.compute_permx_for_month(cfg_r, m, lai_value=v)
                     for m, v in zip(ml, r)])
    assert np.allclose(ks_r[:i_fire], cfg_r.ks_ul), "unburned Ks (ks_ul) before the fire"
    assert np.isclose(ks_r[i_fire], cfg_r.ks_ll), "Ks = ks_ll right after the fire"
    assert ks_r[i_late] > ks_r[i_f], "Ks recovers"
    assert np.ptp(ks[:i_fire]) > 0 and np.ptp(ks_r[:i_fire]) == 0, \
        "raw LAI makes Ks seasonal before the fire; ratio does not"
    assert np.isclose(lai_recovery_ratio(np.full((3, 3), 1.0), np.full((3, 3), 2.0)), 0.5)
    assert lai_recovery_ratio(np.ones((2, 2)), np.zeros((2, 2))) == 1.0
    obs = {pd.Timestamp("2018-03-01"): np.full((2, 2), 1.0),
           pd.Timestamp("2019-03-01"): np.full((2, 2), 3.0),
           pd.Timestamp("2019-08-01"): np.full((2, 2), 0.5),
           pd.Timestamp("2021-03-01"): np.full((2, 2), 9.0)}       # post-fire: ignored
    clim = prefire_monthly_climatology(obs, "2020-07-01")
    assert set(clim) == {3, 8} and np.allclose(clim[3], 2.0) and np.allclose(clim[8], 0.5)
    # ── 2 layers by default; console explanation; incremental saving ──
    assert args.two_layers is True and args.top_depth == 0.1, "2 layers (0.1 m) must be the default"
    assert _parse_args(["--no-two-layers"]).two_layers is False
    assert _parse_args(["--two-layers"]).two_layers is True
    layers = {"n_top": 2, "nstr": 5, "z_top": 0.1, "z_deep0": 0.1, "z_bot": 3.0}
    txt = "\n".join(ks_update_explanation(
        cfg_r, args, ml[i_fire + 2], 0.5, 2.0, r[i_fire + 2], ks_r[i_fire + 2], ks_r[i_fire + 1],
        1 / 6, permx_deep=cfg_r.ks_ul, layers=layers, zones=[1, 2]))
    assert "recovery ratio" in txt and "strata 1-2" in txt and "strata 3-5" in txt \
        and "CHANGED" in txt and "constant" in txt, txt
    txt1 = "\n".join(ks_update_explanation(
        cfg_r, _parse_args(["--no-two-layers"]), ml[0], 2.0, 2.0, 1.0, ks_r[0], None, -0.5))
    assert "every stratum" in txt1 and "initial value" in txt1, txt1
    txt_off = "\n".join(ks_update_explanation(
        cfg_c, args_c, ml[0], 2.0, 2.0, 1.0, args_c.permx, args_c.permx, -0.5))
    assert "no override" in txt_off and "unchanged" in txt_off, txt_off
    assert len(ks_transfer_curve(cfg_r, args)) == 5
    ks_tab = pd.DataFrame({"month": ml, "lai_mean": [0.0] * len(ml), "lai_driver": r,
                           "permx_mm_hr": [KsFire.m_s_to_mm_hr(v) for v in ks_r],
                           "permx_deep_mm_hr": KsFire.m_s_to_mm_hr(cfg_r.ks_ul)})
    print_ks_update_summary(ks_tab, True, fd)
    try:
        import tempfile
        import xarray as xr
    except ImportError:
        print("  (xarray missing: incremental-save check skipped)")
    else:
        t = pd.date_range("2020-01-01", periods=3)
        mk = lambda off: pd.DataFrame({"datetime": np.repeat(t + pd.Timedelta(days=off), 2),
                                       "node": np.tile([1, 2], 3), "psi": 1.0})
        recs = [mk(0), mk(2)]                       # 2nd month repeats one boundary date
        et = [xr.Dataset({"v": (("datetime", "X"), np.ones((3, 2)))},
                         coords={"datetime": t + pd.Timedelta(days=3 * k), "X": [0, 1]})
              for k in range(2)]
        with tempfile.TemporaryDirectory() as td:
            for n_m in (1, 2):                      # save after month 1, then after month 2
                asm = assemble_outputs(recs[:n_m], [x.rename(columns={"psi": "sw"}) for x in recs[:n_m]],
                                       et[:n_m], et[:n_m], [pd.DataFrame({"a": [1]})] * n_m)
                write_outputs(td, asm)
                got = xr.open_dataset(Path(td) / "psi_output.nc")
                assert got.sizes["datetime"] == (3 if n_m == 1 else 5), got.sizes   # shared boundary day kept once
                got.close()
            assert not list(Path(td).glob("*.tmp*")), "no temp file left behind"
            assert {p.name for p in Path(td).iterdir()} == {
                "psi_output.nc", "sw_output.nc", "et_output.nc", "recharge_output.nc",
                "hgraph_output.csv"}
    assert ic_settings(args) == {"INDP": 0, "pressure_head_ini": -15.0}, "default IC: -15 m uniform"
    assert ic_settings(_parse_args(["--ic-head", "-5"]))["pressure_head_ini"] == -5.0
    assert ic_settings(_parse_args(["--ic-wtposition", "8"])) == {"INDP": 3, "WTPOSITION": 8.0}
    describe_ks_cfg(cfg, args)
    print("smoke test OK  (Ks range mm/hr: "
          f"{KsFire.m_s_to_mm_hr(ks.min()):.1f} - {KsFire.m_s_to_mm_hr(ks.max()):.1f})")


# ── Pipeline ──────────────────────────────────────────────────────────────────

EXIT_PARTIAL = 3     # some months failed, outputs of the others were saved
EXIT_NO_OUTPUT = 2   # not a single month could be simulated


def _parm_solver_report(simu, args) -> str:
    """Solver settings as held by the simulation and as written in the parm file
    CATHY reads. '' if neither can be read. Flags any value that differs from args."""
    want = {"DTMIN": args.dtmin, "DTMAX": args.dtmax, "DELTAT": args.deltat}
    parts = []
    parm = getattr(simu, "parm", None)
    if isinstance(parm, dict) and all(k in parm for k in want):
        bad = [k for k, v in want.items() if not np.isclose(float(parm[k]), v, rtol=1e-6)]
        parts.append("simu.parm " + ", ".join(f"{k}={float(parm[k]):g}" for k in want)
                     + (f"  MISMATCH vs requested {want}" if bad else "  OK"))
    try:
        f = Path(simu.workdir) / simu.project_name / simu.input_dirname / "parm"
        rows = [ln.strip() for ln in f.read_text().splitlines()
                if any(ln.lstrip().startswith(k) or f" {k}" in ln for k in want)]
        if rows:
            parts.append("parm file: " + " | ".join(r.split("!")[0].strip() for r in rows[:6]))
    except Exception:
        pass
    return " || ".join(parts)


def run_month_solver(simu, args, t_end) -> tuple:
    """Run CATHY once for the current month (fixed DTMIN/DTMAX/DELTAT, no retry).
    The solver settings are written to the parm file explicitly before the run and
    read back, because run_processor kwargs were seen NOT to reach the parm file.
    Never raises. Returns (ok, message)."""
    try:
        simu.update_parm(DTMIN=args.dtmin, DTMAX=args.dtmax, DELTAT=args.deltat)
        rep = _parm_solver_report(simu, args)
        if rep:
            print(f"      solver settings: {rep}")
        simu.run_processor(IPRT1=2, DTMIN=args.dtmin, DTMAX=args.dtmax, DELTAT=args.deltat,
                           TRAFLAG=0, verbose=False)
        last = float(simu.read_outputs("psi").index[-1])
        if last >= t_end:
            return True, ""
        msg = f"solver stopped at t={last:g} s < {t_end:g} s"
    except Exception as exc:
        msg = f"{type(exc).__name__}: {exc}"
    print(f"    ! {msg}")
    return False, msg


def run_pipeline(args: argparse.Namespace) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import xarray as xr
    import rioxarray  # noqa: F401
    from pyCATHY import cathy_utils
    import Agramon_withLAI_withETp as AG
    import results_plotter as plotter

    failures: list[str] = []

    def _safe(label, fn, *a, **k):
        try:
            return fn(*a, **k)
        except Exception as exc:  # tracked, surfaced at the end
            msg = f"{label}: {type(exc).__name__}: {exc}"
            print(f"  WARNING: {msg}")
            failures.append(msg)
            return None

    fire_date = datetime.strptime(args.ks_fire_date, "%Y-%m-%d").date()
    layer_tag = "_2L" if args.two_layers else ""
    scenario = f"synth_ks-{args.ks_fire_mode}{layer_tag}_{args.start_year}-{args.end_year - 1}"
    # Absolute path: pyCATHY changes the working directory while it compiles /
    # runs CATHY, so a relative path like ../synthetic_LAI_Ks would then point
    # somewhere else (outputs/ "not found" when saving the CSV).
    args.path2prj = str(Path(args.path2prj).expanduser().resolve())
    sim_dir = Path(args.path2prj) / scenario
    out_dir = sim_dir / "outputs"
    fig_dir = sim_dir / "figures"
    for d in (out_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    # 1. Forcing (same ATMBC as Agramon: rain - ETp, scalar) ------------------
    print("[1/6] Forcing …")
    etp_mmd, rain_mmd = load_forcing_series(
        args.start_year, args.end_year, args.etp_file, args.rain_file)
    etp_mmd, rain_mmd = etp_mmd * args.etp_scale, rain_mmd * args.rain_scale
    _yrs = max((etp_mmd.index[-1] - etp_mmd.index[0]).days + 1, 1) / 365.25
    print(f"  Forcing: {len(etp_mmd)} days | rain x{args.rain_scale:g} = "
          f"{rain_mmd.sum() / _yrs:.0f} mm/yr | ETp x{args.etp_scale:g} = "
          f"{etp_mmd.sum() / _yrs:.0f} mm/yr | max daily rain {rain_mmd.max():.1f} mm/d")
    valid_times = etp_mmd.index.to_numpy()
    t_atmbc = ((valid_times - valid_times[0]) / np.timedelta64(1, "s")).astype(float)
    etp_ms = etp_mmd.to_numpy() * 1e-3 / 86400
    rain_ms = rain_mmd.to_numpy() * 1e-3 / 86400
    net_flux = rain_ms - etp_ms

    # 2. Synthetic LAI + Ks config -------------------------------------------
    shape = (args.n, args.n)
    forcing_months = {(t.year, t.month) for t in pd.DatetimeIndex(valid_times)}
    months = [m for m in pd.date_range(f"{args.start_year}-01-01",
                                       f"{args.end_year - 1}-12-01", freq="MS")
              if (m.year, m.month) in forcing_months]
    lai_all = {m: synthetic_lai(m, shape, fire_date, args.lai_recovery_rate) for m in months}
    # unburned / pre-fire reference of the same calendar month (for real data:
    # prefire_monthly_climatology(observed_LAI, fire_date)[month])
    lai_unburned = {m: synthetic_lai(m, shape) for m in months}
    lai_driver = lai_driver_values(args.lai_ref, months, lai_all, lai_unburned)
    cfg = finalize_ks_cfg(KsFire.ks_fire_config_from_args(args),
                          np.array([lai_driver[m] for m in months]), args.lai_ref)
    if cfg.mode == "off":   # control: constant Ks, same value all run long
        args.permx = resolve_static_permx(args, cfg)
    permx_deep = resolve_deep_permx(args, cfg) if args.two_layers else None
    print(f"  Ks mode: {cfg.mode} | {len(months)} months | "
          f"{'2 layers' if args.two_layers else '1 layer'}")
    describe_ks_cfg(cfg, args)

    # 3. Synthetic DEM -> mesh -------------------------------------------------
    print("[2/6] Synthetic DEM + mesh …")
    dem = synthetic_dem(args.n, args.res)
    simu = AG.initialise_project(args.path2prj, args.start_year, args.end_year,
                                 args.zroot, args.pmin, scenario)
    # (re)create outputs/figures: pyCATHY may have reset the project folder
    for d in (out_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    # 2b. Input trends BEFORE the simulation (LAI, ETp/rain, Ks) -------------
    plotter.FIRE_DATE = pd.Timestamp(fire_date)
    df_plan = planned_ks_series(cfg, args, months,
                                [lai_all[m].mean() for m in months],
                                [lai_driver[m] for m in months], permx_deep)
    df_plan.to_csv(out_dir / "inputs_planned_monthly.csv", index=False)
    if not args.no_plots:
        print("  Plotting input trends (LAI, ETp, Ks) before the run …")
        _safe("fig inputs LAI/ETp/Ks", plot_pre_simulation, plt, plotter, fig_dir, df_plan,
              lai_all, etp_mmd, rain_mmd, fire_date, cfg, args.two_layers, args.dpi,
              args.lai_ref)

    AG.setup_dem(simu, dem, xllcorner=XLL, yllcorner=YLL,
                 delta_x=args.res, delta_y=args.res, show=False)
    simu.run_preprocessor(verbose=False)
    grid3d = simu.read_outputs("grid3d")
    n_surf = int(grid3d["nnod"])
    print(f"  Grid {int(simu.hapin['N'])}x{int(simu.hapin['M'])} | surface nodes {n_surf}")

    n_top, layers_info = None, None
    if args.two_layers:
        from pyCATHY import meshtools as mt
        l_tops, l_bots = mt.get_layer_depths(simu.dem_parameters)
        n_top = n_top_layers(l_tops, l_bots, args.top_depth, args.top_layers)
        layers_info = {"n_top": n_top, "nstr": len(l_tops),
                       "z_top": float(l_bots[n_top - 1]), "z_deep0": float(l_tops[n_top]),
                       "z_bot": float(l_bots[-1])}
        print(f"  Two-layer soil: nstr={len(l_tops)} | top = strata 1-{n_top} "
              f"(0 - {l_bots[n_top - 1]:.2f} m, burned) | deep = strata {n_top + 1}-{len(l_tops)} "
              f"({l_tops[n_top]:.2f} - {l_bots[-1]:.2f} m, stable)")

    def _nodes(n):
        return (np.repeat(etp_ms[:, None], n, axis=1),
                np.repeat(rain_ms[:, None], n, axis=1))
    ETp_nodes, rain_nodes = _nodes(n_surf)

    # 4. ICs / BCs -------------------------------------------------------------
    ic_kw = ic_settings(args)
    print("  Initial condition: " + (f"hydrostatic, water table {ic_kw['WTPOSITION']:g} m "
          f"above the grid base" if ic_kw["INDP"] == 3
          else f"uniform pressure head {ic_kw['pressure_head_ini']:g} m"))
    simu.update_ic(**ic_kw)
    bc_time_axis = np.arange(0, 32 * 86400, 86400).tolist()
    AG.configure_boundary_conditions(
        simu, bc_time_axis, outlet_flux=args.outlet_flux,
        bottom_flux=args.bottom_flux, res_x=args.res, res_y=args.res)

    # 5. Monthly loop ------------------------------------------------------------
    print("[3/6] Monthly loop …")
    psi_rec, sw_rec, et_rec, rec_rec, hg_rec = [], [], [], [], []
    veg_hist, lai_hist, ks_rows = {}, {}, []
    first_iteration = True
    prev_permx, prev_zones = None, None
    state = {"psi_last": None}     # end-of-month pressure head -> next month's IC
    failed_months: list = []       # months the solver could not finish
    consec_fail = 0

    # ── saving DURING the run ──
    # _area_analysis: area-mean CSVs (+ figures only at the end of the run).
    # _checkpoint: after every month rewrite everything that exists so far, so a
    # crash / Ctrl-C keeps all the months already simulated (atomic writes).
    def _area_analysis(asm, plots, complete):
        # REQUIRED by synthetic_ksfire_analysis.py for the cross-scenario comparison
        import synthetic_ksfire_analysis as SA
        ts, info = SA.build_timeseries(
            asm["psi"], asm["sw"], asm["et"], asm["rech"], grid3d, et_var=ET_VAR,
            et_scale=ET_SCALE, rain_mmd=rain_mmd, etp_mmd=etp_mmd, lai_hist=lai_hist,
            ks_df=pd.DataFrame(ks_rows))
        meta = {"scenario": scenario, "mode": cfg.mode, "two_layers": bool(args.two_layers),
                "failed_months": failed_months, "complete": complete,
                "months_simulated": len(psi_rec), "months_planned": len(months)}
        SA.save_run_analysis(out_dir, fig_dir, ts, info, meta, pd.Timestamp(fire_date),
                             plots=plots, dpi=args.dpi)

    def _checkpoint(full=True):
        """Light (every month): ks_fire_monthly.csv + failures.csv.
        full: also the merged NetCDF/CSV outputs and the area-mean CSVs."""
        try:
            out_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(ks_rows).to_csv(out_dir / "ks_fire_monthly.csv", index=False)
            pd.Series(failures, name="failure", dtype=str).to_csv(
                out_dir / "failures.csv", index=False)
            what = "ks_fire_monthly.csv, failures.csv"
            if full and psi_rec and sw_rec and et_rec:
                asm = assemble_outputs(psi_rec, sw_rec, et_rec, rec_rec, hg_rec)
                write_outputs(out_dir, asm)
                what = "psi/sw/ET/recharge/hgraph, " + what
                try:
                    _area_analysis(asm, False, False)
                    what += ", area-mean CSVs"
                except Exception as exc:
                    what += f" (area-mean skipped: {type(exc).__name__}: {exc})"
            print(f"      saved [{len(psi_rec)}/{len(months)} months]: {what}")
        except Exception as exc:     # never let a save problem kill the simulation
            print(f"  WARNING: progress save failed: {type(exc).__name__}: {exc}")

    sim_count = 0

    # monthly backup: reuse months already simulated with the same settings
    bdir = Path(args.path2prj) / "_monthly_backup" / scenario
    resume_ok = prepare_backup(bdir, run_fingerprint(args, months), resume=not args.no_resume)
    n_restored = 0

    try:
        for dt in months:
            if resume_ok:
                bk = load_month_backup(bdir, dt)
                if bk is not None:
                    for lst, key in ((psi_rec, "psi"), (sw_rec, "sw"), (et_rec, "et"),
                                     (rec_rec, "rec"), (hg_rec, "hg")):
                        if bk.get(key) is not None:
                            lst.append(bk[key])
                    ks_rows.append(bk["ks_row"])
                    veg_hist[dt], lai_hist[dt] = bk["veg_map"], lai_all[dt].copy()
                    prev_permx, prev_zones = bk["permx"], bk["zones"]
                    state["psi_last"] = bk["psi_last"]
                    n_restored += 1
                    print(f"  {dt:%Y-%m} | restored from backup (not re-simulated) | "
                          f"LAI {bk['ks_row']['lai_mean']:.2f} | "
                          f"Ks {bk['ks_row']['permx_mm_hr']:.1f} mm/hr")
                    continue
                resume_ok = False   # a gap breaks the initial-condition chain
                print(f"  {dt:%Y-%m} | no backup -> simulating from here "
                      f"({n_restored} month(s) restored)")
            n_fail0 = len(failures)
            lai_2d = lai_all[dt]
            veg_remapped, remap = AG.lai_to_veg_map(simu, lai_2d)
            veg_map = AG.smooth_veg_map(veg_remapped)
            _, veg_valid = simu._check_outside_DEM(veg_map)
            simu.update_veg_map(veg_map, maxveg=AG.N_VEG_CLASSES)
            veg_hist[dt], lai_hist[dt] = veg_map.copy(), lai_2d.copy()

            # ── Ks from ks_fire_utils (the point of this test) ──
            permx_fire = KsFire.compute_permx_for_month(
                cfg, dt, lai_value=lai_driver[dt])
            permx = permx_fire if permx_fire is not None else args.permx
            ks_rows.append({"month": dt, "lai_mean": float(lai_2d.mean()),
                            "lai_driver": lai_driver[dt],
                            "lai_min": float(lai_2d.min()), "lai_max": float(lai_2d.max()),
                            "permx_m_s": permx,
                            "permx_mm_hr": KsFire.m_s_to_mm_hr(permx) if permx else np.nan,
                            "permx_deep_m_s": permx_deep if args.two_layers else np.nan,
                            "permx_deep_mm_hr": (KsFire.m_s_to_mm_hr(permx_deep)
                                                 if args.two_layers else np.nan)})
            AG.update_soil_for_active_zones(
                simu, set(int(z) for z in veg_valid), remap=remap,
                base_zroot=args.zroot, pmin=args.pmin, scf=args.scf,
                zroot_factors=AG.ZROOT_FACTORS, permx=permx,
                permy_ratio=args.permy_ratio, permz_ratio=args.permz_ratio)
            if args.two_layers:
                # AG wrote the same PERMX in every stratum (and set the Feddes/veg
                # tables). Re-apply the SPP_map with the burned top / stable deep
                # split; FP_map is re-used as is, no extra kwargs (they would be
                # stored in the soil/parm files).
                spp_2l = two_layer_spp(simu.soil_SPP["SPP_map"], n_top, permx, permx_deep,
                                       args.permy_ratio, args.permz_ratio)
                simu.update_soil(SPP_map=spp_2l, FP_map=simu.soil_FP["FP_map"])

            # ── soil-change check for this iteration ──
            rep = soil_change_report(prev_permx, permx)
            zones = sorted(int(z) for z in veg_valid)
            zones_changed = prev_zones is not None and zones != prev_zones
            rb = read_back_permx(simu)
            rb_top = rb_deep = np.nan
            if args.two_layers:
                rbl = read_back_permx_layers(simu, n_top)
                if rbl is None:
                    rb_txt = "layer read-back n/a"
                else:
                    rb_top, rb_deep = rbl
                    ok_top = permx is None or np.isclose(rb_top, permx, rtol=1e-6)
                    ok_deep = np.isclose(rb_deep, permx_deep, rtol=1e-6)
                    rb_txt = (f"read-back top {KsFire.m_s_to_mm_hr(rb_top):.2f} mm/hr "
                              f"{'OK' if ok_top else 'MISMATCH'} | deep "
                              f"{KsFire.m_s_to_mm_hr(rb_deep):.2f} mm/hr "
                              f"{'OK (stable)' if ok_deep else 'MISMATCH'}")
                    if not ok_top:
                        failures.append(f"{dt:%Y-%m}: top-layer PERMX read-back {rb_top:.3e} "
                                        f"!= expected {permx:.3e}")
                    if not ok_deep:
                        failures.append(f"{dt:%Y-%m}: deep-layer PERMX read-back {rb_deep:.3e} "
                                        f"!= expected (stable) {permx_deep:.3e}")
            elif rb is None:
                rb_txt = "read-back n/a"
            else:
                ok = permx is None or np.isclose(rb[0], permx, rtol=1e-6)
                rb_txt = (f"read-back PERMX {KsFire.m_s_to_mm_hr(rb[0]):.2f} mm/hr "
                          + ("OK" if ok else "MISMATCH"))
                if not ok:
                    failures.append(f"{dt:%Y-%m}: soil PERMX read-back {rb[0]:.3e} "
                                    f"!= expected {permx:.3e}")
            dtxt = "" if np.isnan(rep["delta_pct"]) else f" ({rep['delta_pct']:+.1f}%)"
            lbl = f"top Ks {rep['label']}" if args.two_layers else f"Ks {rep['label']}"
            deep_txt = (f" (top) / {ks_rows[-1]['permx_deep_mm_hr']:.1f} mm/hr (deep)"
                        if args.two_layers else "")
            print(f"  {dt:%Y-%m} | LAI {lai_2d.mean():.2f} "
                  f"(min {lai_2d.min():.2f}, max {lai_2d.max():.2f}) | "
                  f"Ks {ks_rows[-1]['permx_mm_hr']:.1f} mm/hr{deep_txt}")
            try:   # how ks_fire_utils turned this month's LAI into the soil update
                for ln in ks_update_explanation(
                        cfg, args, dt, float(lai_2d.mean()), float(np.nanmean(lai_unburned[dt])),
                        lai_driver[dt], permx, prev_permx,
                        KsFire.years_since_fire(dt.to_pydatetime(), fire_date),
                        permx_deep=permx_deep, layers=layers_info, zones=zones):
                    print(f"      {ln}")
            except Exception as exc:
                print(f"      (Ks explanation unavailable: {type(exc).__name__}: {exc})")
            print(f"      soil check: {lbl}{dtxt} | veg zones {zones}"
                  f"{' (changed)' if zones_changed else ''} | {rb_txt}")
            ks_rows[-1].update({"ks_changed": rep["changed"], "ks_delta_pct": rep["delta_pct"],
                                "veg_zones": str(zones), "veg_zones_changed": zones_changed,
                                "readback_permx_m_s": rb[0] if rb else np.nan,
                                "readback_top_m_s": rb_top, "readback_deep_m_s": rb_deep})
            prev_permx, prev_zones = permx, zones

            t_month, nv_month = AG.subset_atmbc_for_month(
                t_atmbc, net_flux, valid_times, dt.year, dt.month)
            if t_month is None:
                failures.append(f"{dt:%Y-%m}: no forcing — skipped")
                continue
            simu.update_atmbc(HSPATM=1, IETO=1, time=t_month, netValue=nv_month)
            idx_month = np.where([pd.Timestamp(t).year == dt.year and
                                  pd.Timestamp(t).month == dt.month
                                  for t in valid_times])[0]
            simu.update_parm(TIMPRTi=list(simu.read_inputs("atmbc").time), IPRT=4, VTKF=2)

            if first_iteration:
                first_iteration = False
                simu.create_mesh_vtk(verbose=False)
                grid3d = simu.read_outputs("grid3d")
                if int(grid3d["nnod"]) != n_surf:
                    n_surf = int(grid3d["nnod"])
                    ETp_nodes, rain_nodes = _nodes(n_surf)
            if state["psi_last"] is not None:
                simu.update_ic(INDP=1, pressure_head_ini=state["psi_last"])
            n_before = [len(l) for l in (psi_rec, sw_rec, et_rec, rec_rec, hg_rec)]

            print("      running CATHY …")
            ok_run, why = run_month_solver(simu, args, float(t_month[-1]))
            if not ok_run:
                # do NOT let one bad month kill the run: log it, skip it, go on.
                # (the next month restarts from the last COMPLETED month's psi)
                failed_months.append(f"{dt:%Y-%m}")
                failures.append(f"psi {dt:%Y-%m}: {why} -> month skipped")
                consec_fail += 1
                print(f"  {dt:%Y-%m} | FAILED, month skipped ({consec_fail} in a row)")
                _checkpoint(full=False)
                if args.debug:
                    print(f"  >>> --debug: {dt:%Y-%m} did not complete ({why}); "
                          f"stopping the monthly loop. Saving what exists.")
                    break
                if args.max_consecutive_fail and consec_fail >= args.max_consecutive_fail:
                    print(f"  >>> {consec_fail} months in a row failed: stopping the monthly "
                          f"loop (see --max-consecutive-fail / "
                          f"--ic-wtposition / --dtmin). Saving what exists.")
                    break
                continue
            consec_fail = 0

            def _psi():
                df = simu.read_outputs("psi").copy()
                if df.index[-1] < t_month[-1]:
                    raise RuntimeError(f"solver stopped at {df.index[-1]} < {t_month[-1]}")
                state["psi_last"] = df.iloc[-1].to_numpy(copy=True)   # before 'datetime' col
                df["datetime"] = cathy_utils.change_x2date(df.index, dt)
                psi_rec.append(df.reset_index(drop=True)
                                 .melt(id_vars="datetime", var_name="node", value_name="psi")
                                 .assign(node=lambda d: d["node"].astype(int)))
            _safe(f"psi {dt:%Y-%m}", _psi)

            def _sw():
                df, _ = simu.read_outputs("sw")
                df = df.copy()
                df["datetime"] = cathy_utils.change_x2date(df.index, dt)
                sw_rec.append(df.reset_index(drop=True)
                                .melt(id_vars="datetime", var_name="node", value_name="sw")
                                .assign(node=lambda d: d["node"].astype(int)))
            _safe(f"sw {dt:%Y-%m}", _sw)

            def _et():
                df = AG.add_datetime_col(simu.read_outputs("ET").copy(), dt)
                df, _, _ = AG.apply_eta_artefact_floor(
                    df, t_month, ETp_nodes, rain_nodes, idx_month,
                    net_neg_mode=args.eta_net_neg_mode)
                x = df.set_index(["X", "Y", "datetime"]).to_xarray()
                x = x.rio.set_spatial_dims(x_dim="X", y_dim="Y").rio.write_crs(SYN_CRS)
                et_rec.append(x)
            _safe(f"ET {dt:%Y-%m}", _et)

            def _rec():
                df = AG.add_datetime_col(simu.read_outputs("recharge").copy(), dt)
                x = df.set_index(["X", "Y", "datetime"]).to_xarray()
                x = x.rio.set_spatial_dims(x_dim="X", y_dim="Y").rio.write_crs(SYN_CRS)
                rec_rec.append(x)
            _safe(f"recharge {dt:%Y-%m}", _rec)

            def _hg():
                df = simu.read_outputs("hgraph").copy().rename(columns={"time": "time_sec"})
                hg_rec.append(AG.add_datetime_col(df, dt))
            _safe(f"hgraph {dt:%Y-%m}", _hg)

            # ── backup this month if everything worked (psi/sw/ET are mandatory) ──
            lists = (psi_rec, sw_rec, et_rec, rec_rec, hg_rec)
            added = [len(l) - n for l, n in zip(lists, n_before)]
            if len(failures) == n_fail0 and all(a == 1 for a in added[:3]):
                _safe(f"backup {dt:%Y-%m}", save_month_backup, bdir, dt, {
                    "psi": psi_rec[-1], "sw": sw_rec[-1], "et": et_rec[-1],
                    "rec": rec_rec[-1] if added[3] == 1 else None,
                    "hg": hg_rec[-1] if added[4] == 1 else None,
                    "psi_last": state["psi_last"], "ks_row": ks_rows[-1],
                    "veg_map": veg_map.copy(), "permx": permx, "zones": zones})

            # ── save everything simulated so far (do not wait for the end) ──
            sim_count += 1
            _checkpoint(full=args.save_every > 0 and sim_count % args.save_every == 0)
    finally:
        # never let a missing folder mask the real exception of the loop
        out_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(ks_rows).to_csv(out_dir / "ks_fire_monthly.csv", index=False)

    # 6. Assemble + save ---------------------------------------------------------
    _safe("Ks update table", print_ks_update_summary, pd.DataFrame(ks_rows),
          args.two_layers, fire_date)
    print("[4/6] Assembling outputs …")
    if not (psi_rec and sw_rec and et_rec):
        pd.Series(failures, name="failure", dtype=str).to_csv(out_dir / "failures.csv", index=False)
        print(f"\n  !! NO month could be simulated ({len(failed_months)} failed: "
              f"{failed_months[0] if failed_months else '-'} ... "
              f"{failed_months[-1] if failed_months else '-'}).")
        print(f"  First failure : {failures[0] if failures else 'n/a'}")
        print(f"  All failures  : {out_dir / 'failures.csv'}")
        print("  Every month failing from the first one is a MODEL SETUP problem (initial "
              "condition / boundary conditions / Ks / atmbc), not a Ks-recovery one. Try e.g.\n"
              "    --ic-wtposition <other depth>   --dtmin 1e-4 --deltat 10   --ks-fire-mode off")
        sys.exit(EXIT_NO_OUTPUT)
    asm = assemble_outputs(psi_rec, sw_rec, et_rec, rec_rec, hg_rec)
    psi_xr, sw_xr, ET_xr_all = asm["psi"], asm["sw"], asm["et"]
    xr_rech, df_hg = asm["rech"], asm["hg"]
    _safe("save netcdf", write_outputs, out_dir, asm)

    # area-mean analysis (+ per-run figures) -> outputs/areamean_timeseries.csv
    _safe("area-mean analysis", _area_analysis, asm, not args.no_plots, True)

    # 7. Plots -------------------------------------------------------------------
    if not args.no_plots:
        print("[5/6] Plotting …")
        dpi = args.dpi
        plotter.FIRE_DATE = pd.Timestamp(fire_date)
        df_ks = pd.DataFrame(ks_rows)

        def _ks_fig():
            fig, (a0, a1) = plt.subplots(2, 1, figsize=(11, 6), sharex=True,
                                         constrained_layout=True)
            a0.fill_between(df_ks["month"], df_ks["lai_min"], df_ks["lai_max"],
                            color="#4daf4a", alpha=0.25, label="LAI min–max")
            a0.plot(df_ks["month"], df_ks["lai_mean"], "s-", color="#006837",
                    label="LAI domain mean")
            a0.set_ylabel("LAI (m² m⁻²)")
            a1.plot(df_ks["month"], df_ks["permx_mm_hr"], "o-", color="#762a83",
                    label="PERMX top (burned)" if args.two_layers else "PERMX")
            if args.two_layers:
                a1.plot(df_ks["month"], df_ks["permx_deep_mm_hr"], "s--", color="#1b7837",
                        label="PERMX deep (stable)")
            a1.set_ylabel("PERMX (mm/hr)")
            a1.set_title(f"Ks driven by ks_fire_utils (mode={cfg.mode}, "
                         f"{'2 layers' if args.two_layers else '1 layer'}"
                         f"{', control: constant Ks' if cfg.mode == 'off' else ''})",
                         fontsize=10, fontweight="bold")
            for a in (a0, a1):
                plotter.add_fire_vline(a, df_ks["month"].min(), df_ks["month"].max())
                a.grid(True, linestyle="--", alpha=0.3)
                a.legend(loc="upper right", fontsize=8)
            plotter.save_fig(fig, fig_dir, "fig_lai_vs_ks.png", dpi)
            plt.close(fig)
        _safe("fig LAI vs Ks", _ks_fig)

        nodes_dict = _safe("resolve nodes", resolve_nodes_simple, grid3d, plotter)
        pev = etp_mmd.copy()
        tp = rain_mmd.copy()
        lai_ds = plotter.lai_history_to_dataset(lai_hist, lai_var="LAI")
        kw = dict(dpi=dpi, show=False, scenario_id=scenario)
        d0, d1 = f"{args.start_year}-01-01", f"{args.end_year - 1}-12-31"

        _safe("ET monthly mean", plotter.plot_et_monthly_mean,
              ET_xr_all, ET_VAR, ET_SCALE, fig_dir, dpi=dpi)
        if nodes_dict is not None:
            _safe("hydro", plotter.plot_hydro_xr, psi_xr, sw_xr, nodes_dict, fig_dir,
                  pev_series=pev, tp_series=tp, date_start=d0, date_end=d1, **kw)
            _safe("et+sw+lai", plotter.plot_et_sw_lai_xr, sw_xr, ET_xr_all, ET_VAR,
                  ET_SCALE, lai_ds, nodes_dict, fig_dir, date_start=d0, date_end=d1, **kw)
            _safe("et+psi+lai", plotter.plot_et_psi_lai_xr, psi_xr, ET_xr_all, ET_VAR,
                  ET_SCALE, lai_ds, nodes_dict, fig_dir, date_start=d0, date_end=d1, **kw)
            _safe("node map", plotter.plot_node_locations, grid3d, nodes_dict, fig_dir,
                  scenario_id=scenario, dpi=dpi)
        if df_hg is not None:
            _safe("hgraph", plotter.plot_hgraph_csv, df_hg, fig_dir, dpi=dpi,
                  scenario_id=scenario)
        if xr_rech is not None:
            _safe("recharge", plotter.plot_recharge_xr, xr_rech, fig_dir, dpi=dpi,
                  scenario_id=scenario)

        # withLAI_results helpers (spatial maps / LAI diagnostics)
        try:
            import withLAI_results as WR
        except Exception as exc:
            failures.append(f"import withLAI_results: {exc}")
            WR = None
        if WR is not None:
            years = sorted({m.year for m in months})
            for yr in years:
                _safe(f"ET maps {yr}", WR.plot_et_spatial_monthly, ET_xr_all, ET_VAR,
                      ET_SCALE, yr, fig_dir, cmap=WR.CMAP_ET, dpi=dpi)
            ds_sw = _safe("sw surface grid", WR.build_sw_surface_grid, sw_xr, grid3d, ET_xr_all)
            if ds_sw is not None:
                for yr in years:
                    _safe(f"sw maps {yr}", WR.plot_sw_spatial_monthly, ds_sw, "sw", yr,
                          fig_dir, cmap=WR.CMAP_SW, dpi=dpi)
            _safe("LAI monthly by year", WR.plot_index_monthly_by_year, lai_hist,
                  index_var="LAI", out_dir=fig_dir, scenario_id=scenario, dpi=dpi)
            for yr in years:
                _safe(f"LAI maps {yr}", WR.plot_index_spatial_maps_year, lai_hist,
                      index_var="LAI", selected_year=yr, out_dir=fig_dir,
                      scenario_id=scenario, dpi=dpi)
            _safe("LAI annual mean", WR.plot_index_annual_mean, lai_hist,
                  index_var="LAI", out_dir=fig_dir, scenario_id=scenario, dpi=dpi)

    # 8. Report ------------------------------------------------------------------
    print("[6/6] Done.")
    pd.Series(failures, name="failure", dtype=str).to_csv(out_dir / "failures.csv", index=False)
    print(f"  Outputs : {out_dir}\n  Figures : {fig_dir}")
    print(f"  Failures tracked: {len(failures)}" + ("" if not failures else
          "\n    - " + "\n    - ".join(failures)))
    if failed_months:
        print(f"  PARTIAL run: {len(failed_months)} month(s) not simulated: "
              f"{', '.join(failed_months)}")
        sys.exit(EXIT_PARTIAL)


if __name__ == "__main__":
    _args = _parse_args()
    if _args.smoke_test:
        smoke_test()
        sys.exit(0)
    run_pipeline(_args)
