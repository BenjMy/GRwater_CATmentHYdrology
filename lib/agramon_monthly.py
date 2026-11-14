"""
agramon_monthly.py — read one month of CATHY outputs (REFACTORING.md 2.2)
=========================================================================
Bodies moved verbatim from the monthly loop of ``Agramon_withLAI_withETp``
(steps 12.9-12.11a). Each function appends to the caller's record list,
fills ``month_stats`` and returns the monthly frame used by the sanity
checks (``None`` if the read failed). Errors are reported, never raised.
``monthly_dir`` is the directory for per-month files, or ``None`` when
``--save-monthly`` is off.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyCATHY import cathy_utils

from agramon_checks import add_datetime_col, apply_eta_artefact_floor
from agramon_config import N_VEG_CLASSES, TARGET_CRS
from agramon_lai import lai_to_veg_map, smooth_veg_map


def read_sw_month(simu, dt, sw_records, month_stats, monthly_dir):
    """Step 12.9: saturation (SW)."""
    sw_month_df = None
    try:
        df_sw, _ = simu.read_outputs("sw")
        df_sw = df_sw.copy()
        df_sw["datetime"] = cathy_utils.change_x2date(df_sw.index, dt)
        sw_month_df = (
            df_sw.reset_index(drop=True)
                 .melt(id_vars="datetime", var_name="node", value_name="sw")
                 .assign(node=lambda d: d["node"].astype(int))
        )
        sw_records.append(sw_month_df)
        month_stats["sw_mean"] = float(sw_month_df["sw"].mean())
        if monthly_dir is not None:
            try:
                (
                    sw_month_df
                    .groupby(["datetime", "node"], as_index=False)
                    .last()
                    .set_index(["datetime", "node"])
                    .to_xarray()
                    .to_netcdf(monthly_dir / f"sw_{dt:%Y%m}.nc")
                )
            except Exception as exc:
                print(f"  Warning: could not save monthly SW for {dt:%Y-%m} ({exc})")
    except Exception as exc:
        print(f"  Warning: could not read SW for {dt:%Y-%m} ({exc})")
    return sw_month_df


def read_et_month(simu, dt, et_records, month_stats, monthly_dir, md_log, *,
                  use_etp, t_month, ETp_nodes, rain_nodes, idx_month,
                  eta_net_neg_mode):
    """Step 12.10: actual ET (with the net-forcing artefact floor in ETp mode)."""
    df_et = None
    try:
        df_et = simu.read_outputs("ET").copy()
        df_et = add_datetime_col(df_et, dt)

        # FIX: floor the net-forcing ETa artefact at ETp for this month
        # (adds 'ACT. ETRA_patched'; 'ACT. ETRA' stays unmodified).
        if use_etp:
            df_et, n_eta_patched, n_eta_net_neg = apply_eta_artefact_floor(
                df_et, t_month, ETp_nodes, rain_nodes, idx_month,
                net_neg_mode=eta_net_neg_mode,
            )
            if n_eta_patched or n_eta_net_neg:
                print(f"  ETa patch: {n_eta_patched} cells (net ≥ 0 artefact) + "
                      f"{n_eta_net_neg} cells (net < 0, mode={eta_net_neg_mode}) "
                      f"for {dt:%Y-%m}")
            md_log.bullet(f"{dt:%Y-%m}: ETa artefact floor patched "
                           f"{n_eta_patched} cell(s) (net rain ≥ ETp); "
                           f"{n_eta_net_neg} cell(s) set for net < 0 "
                           f"(mode={eta_net_neg_mode}).")
            md_log.add_metric("ETa artefact cells patched (total)", n_eta_patched)
            md_log.add_metric("ETa net<0 cells set (total)", n_eta_net_neg)

        et_xr = df_et.set_index(["X", "Y", "datetime"]).to_xarray()
        et_xr = et_xr.rio.set_spatial_dims(x_dim="X", y_dim="Y")
        et_xr = et_xr.rio.write_crs(TARGET_CRS)
        et_records.append(et_xr)
        if "ACT. ETRA" in df_et.columns:
            month_stats["eta_mean"] = float(df_et["ACT. ETRA"].mean())
        if "ACT. ETRA_patched" in df_et.columns:
            month_stats["eta_patched_mean"] = float(df_et["ACT. ETRA_patched"].mean())
        if monthly_dir is not None:
            try:
                et_xr.to_netcdf(monthly_dir / f"et_{dt:%Y%m}.nc")
            except Exception as exc:
                print(f"  Warning: could not save monthly ET for {dt:%Y-%m} ({exc})")
    except Exception as exc:
        print(f"  Warning: could not read ET for {dt:%Y-%m} ({exc})")
    return df_et


def read_recharge_month(simu, dt, all_recharge_recs, month_stats, monthly_dir):
    """Step 12.11: recharge."""
    df_rec = None
    try:
        df_rec = simu.read_outputs("recharge").copy()
        df_rec = add_datetime_col(df_rec, dt)
        xr_rec = df_rec.set_index(["X", "Y", "datetime"]).to_xarray()
        xr_rec = xr_rec.rio.set_spatial_dims(x_dim="X", y_dim="Y")
        xr_rec = xr_rec.rio.write_crs(TARGET_CRS)
        all_recharge_recs.append(xr_rec)
        rec_value_cols = [c for c in df_rec.select_dtypes(include=[np.number]).columns
                           if c not in ("X", "Y", "time_sec")]
        if rec_value_cols:
            month_stats["recharge_mean"] = float(df_rec[rec_value_cols].mean().mean())
        if monthly_dir is not None:
            try:
                xr_rec.to_netcdf(monthly_dir / f"recharge_{dt:%Y%m}.nc")
            except Exception as exc:
                print(f"  Warning: could not save monthly recharge for {dt:%Y-%m} ({exc})")
    except Exception as exc:
        print(f"  Warning: could not read recharge for {dt:%Y-%m} ({exc})")
    return df_rec


def read_hgraph_month(simu, dt, all_hgraph_recs, monthly_dir):
    """Step 12.11a: surface-runoff hydrograph at the outlet (CSV, no X/Y dims)."""
    try:
        df_hgraph = simu.read_outputs("hgraph").copy()
        df_hgraph = df_hgraph.rename(columns={"time": "time_sec"})
        df_hgraph = add_datetime_col(df_hgraph, dt)
        all_hgraph_recs.append(df_hgraph)
        if monthly_dir is not None:
            try:
                df_hgraph.to_csv(monthly_dir / f"hgraph_{dt:%Y%m}.csv", index=False)
            except Exception as exc:
                print(f"  Warning: could not save monthly hgraph for {dt:%Y-%m} ({exc})")
    except Exception as exc:
        print(f"  Warning: could not read hgraph for {dt:%Y-%m} ({exc})")


def read_psi_month(simu, dt, psi_records, month_stats, monthly_dir, md_log, t_month):
    """
    Step 12.8: pressure head (PSI). Checks the solver actually reached the end
    of the month; if not, keeps the last known state. Returns
    ``(psi_month_df, solver_status)`` (frame is ``None`` if nothing was read).
    """
    psi_month_df = None
    solver_status = "✓ OK"
    try:
        df_psi = simu.read_outputs("psi").copy()
        last_psi_t   = df_psi.index[-1]
        last_month_t = t_month[-1]
        if last_psi_t < last_month_t:
            print(f"  WARNING: PSI last index ({last_psi_t}) is behind "
                  f"t_month[-1] ({last_month_t}) for {dt:%Y-%m} — "
                  f"solver did not converge; keeping last known state and continuing.")
            solver_status = f"⚠ not converged ({last_psi_t} < {last_month_t})"
            md_log.warn(f"{dt:%Y-%m}: solver did not converge "
                         f"(PSI last index {last_psi_t} < expected {last_month_t}) "
                         f"— keeping last known state.")
            # Do not append new psi_records; psi_ini will reuse previous last row
            # on the next iteration (len(psi_records) > 0 guard still holds).
        else:
            print(f"  time index ic: {last_psi_t} ✓ (t_month[-1]: {last_month_t})")
            df_psi["datetime"] = cathy_utils.change_x2date(df_psi.index, dt)
            psi_month_df = (
                df_psi.reset_index(drop=True)
                      .melt(id_vars="datetime", var_name="node", value_name="psi")
                      .assign(node=lambda d: d["node"].astype(int))
            )
            psi_records.append(psi_month_df)
            month_stats["psi_mean"] = float(psi_month_df["psi"].mean())
            month_stats["psi_min"]  = float(psi_month_df["psi"].min())
            month_stats["psi_max"]  = float(psi_month_df["psi"].max())
            if monthly_dir is not None:
                try:
                    (
                        psi_month_df
                        .groupby(["datetime", "node"], as_index=False)
                        .last()
                        .set_index(["datetime", "node"])
                        .to_xarray()
                        .to_netcdf(monthly_dir / f"psi_{dt:%Y%m}.nc")
                    )
                except Exception as exc:
                    print(f"  Warning: could not save monthly PSI for {dt:%Y-%m} ({exc})")
    except Exception as exc:
        print(f"  WARNING: could not read PSI for {dt:%Y-%m} ({exc}) — "
              f"keeping last known state and continuing.")
        solver_status = f"✗ PSI read failed: {exc}"
        md_log.warn(f"{dt:%Y-%m}: could not read PSI ({exc}) — "
                     f"keeping last known state.")
    return psi_month_df, solver_status


def report_month_stats(dt, month_stats, solver_status, md_log, stats_records, out_dir):
    """
    Step 12.12b: print and log the monthly stats, then rewrite the running
    ``monthly_stats.csv`` in ``out_dir``.
    """
    # 12.12b  Monthly info stats (psi/sw/ETa/recharge means) — printed,
    # logged to the markdown run log, and appended to a single running
    # monthly_stats.csv in out_dir that's rewritten after every month.
    def _fmt(key: str, unit: str = "", sci: bool = False) -> str:
        v = month_stats.get(key)
        if v is None:
            return "n/a"
        return f"{v:.3e}{unit}" if sci else f"{v:.3f}{unit}"

    print(f"  stats     | psi mean: {_fmt('psi_mean', ' m')} "
          f"(min {_fmt('psi_min', ' m')}, max {_fmt('psi_max', ' m')}) | "
          f"sw mean: {_fmt('sw_mean')} | "
          f"ETa mean: {_fmt('eta_mean', ' m/s', sci=True)} | "
          f"recharge mean: {_fmt('recharge_mean', ' m/s', sci=True)}")
    print(f"  sanity    | SW out-of-bounds: {month_stats.get('sw_out_of_bounds_n', 'n/a')} | "
          f"PSI extreme: {month_stats.get('psi_extreme_n', 'n/a')} | "
          f"ETa>ETp: {month_stats.get('eta_exceeds_etp_n', 'n/a')} | "
          f"water balance residual: {_fmt('water_balance_residual_mm', ' mm')}")
    md_log.bullet(
        f"{dt:%Y-%m} stats — psi mean: {_fmt('psi_mean', ' m')}, "
        f"sw mean: {_fmt('sw_mean')}, "
        f"ETa mean: {_fmt('eta_mean', ' m/s', sci=True)}, "
        f"recharge mean: {_fmt('recharge_mean', ' m/s', sci=True)}"
    )
    month_stats["solver_status"] = solver_status
    stats_records.append(month_stats)
    try:
        pd.DataFrame(stats_records).to_csv(out_dir / "monthly_stats.csv", index=False)
    except Exception as exc:
        print(f"  Warning: could not refresh monthly_stats.csv after {dt:%Y-%m} ({exc})")



def prepare_month_vegetation(simu, lai_2d, dt):
    """
    Steps 12.1-12.2: LAI -> consecutive vegetation map (+ remap dict), smooth
    it, save the LAI/veg preview figure and push the map to ``simu``.
    Returns ``(veg_map, remap, veg_valid, veg_map_semantic, inv_remap)``.
    """
    veg_map_remapped, remap = lai_to_veg_map(simu, lai_2d)
    # 12.2  Smooth salt-and-pepper artefacts
    veg_map = smooth_veg_map(veg_map_remapped)

    fig, (ax_lai, ax_veg) = plt.subplots(1, 2, figsize=(8, 6))

    # LAI plot
    im_lai = ax_lai.imshow(lai_2d, cmap="YlGn")
    ax_lai.set_title(f"LAI - {dt:%Y-%m}")
    fig.colorbar(im_lai, ax=ax_lai, label="LAI")

    # Vegetation map plot
    im_veg = ax_veg.imshow(veg_map)
    ax_veg.set_title(f"Vegetation map - {dt:%Y-%m}")
    fig.colorbar(im_veg, ax=ax_veg)

    fig.suptitle(f"{dt:%Y-%m}", fontsize=13, fontweight="bold")
    fig.tight_layout()

    outdir_veg = Path("figures_lai_temp")
    outdir_veg.mkdir(exist_ok=True)
    fig.savefig(outdir_veg / f"veg_map_{dt:%Y%m}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    exclude_veg, veg_valid = simu._check_outside_DEM(veg_map)
    simu.update_veg_map(veg_map, maxveg=N_VEG_CLASSES)

    inv_remap             = {v: k for k, v in remap.items()}
    veg_map_semantic      = np.vectorize(inv_remap.get)(veg_map).astype(int)
    return veg_map, remap, veg_valid, veg_map_semantic, inv_remap


def month_forcing_indices(valid_times, dt):
    """Indices of ``valid_times`` (ETp/rain time axis) that fall in month ``dt``."""
    month_mask = np.array([
        pd.Timestamp(t).year == dt.year and pd.Timestamp(t).month == dt.month
        for t in valid_times
    ])
    return np.where(month_mask)[0]


def plot_month_forcing(nv_month, dt, use_etp):
    """Optional monthly net-flux preview figure (``figures_ETp_temp/``)."""
    outdir_etp = Path("figures_ETp_temp")
    outdir_etp.mkdir(exist_ok=True)
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    ax2.plot(nv_month)
    ax2.set_title(f"Net flux ({'ETp' if use_etp else 'ERA5'}) - {dt:%Y-%m}")
    ax2.set_ylabel("m/s")
    fig2.savefig(outdir_etp / f"forcing_{dt:%Y%m}.png", dpi=150, bbox_inches="tight")
    plt.close(fig2)
