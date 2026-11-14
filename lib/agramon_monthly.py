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

import numpy as np
from pyCATHY import cathy_utils

from agramon_checks import add_datetime_col, apply_eta_artefact_floor
from agramon_config import TARGET_CRS


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
