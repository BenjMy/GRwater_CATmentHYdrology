"""
agramon_checks.py — Post-processing patches and sanity checks (monthly loop)
============================================================================
GRwater project — Agramon study site

  * ``apply_eta_artefact_floor`` — floors the ETa "net-forcing" artefact
    at ETp (temporary patch; see REFACTORING.md step 3.2 for the planned
    withIrr alternative that removes the cause).
  * ``add_datetime_col`` — time_sec → datetime helper.
  * ``check_monthly_sanity`` — physical-bounds / water-balance flags.

Moved verbatim out of Agramon_withLAI_withETp.py (REFACTORING.md, step 2.1).
Where the original file defined a function twice, the *second* definition
(the one Python actually used) was kept.
"""

import numpy as np
import pandas as pd

from agramon_config import ETA_TOL, PSI_SANITY_MIN, SW_TOL, WB_RESIDUAL_FRAC


# ── ETa "net-forcing" artefact patch ────────────────────────────────────────
# Known issue: netValue = rain - ETp is passed to CATHY as a single combined
# flux (HSPATM=1, IETO=1). CATHY has no way to separate "net infiltration"
# from "ET still happening at the same time" — so on any day/node where
# rain >= ETp, the reported ACT. ETRA collapses to 0, even though real ET was
# likely still occurring near the potential rate. This floors ONLY those
# artefact cells (ACT. ETRA ~ 0 AND rain >= ETp) at ETp; genuinely
# soil-moisture-limited zero-ETa days (dry soil, no rain) are left untouched.
def apply_eta_artefact_floor(
    df_et: pd.DataFrame,
    t_month: list,
    ETp_nodes: np.ndarray,
    rain_nodes: np.ndarray,
    idx_month: np.ndarray,
    eta_eps: float = 1e-9,
    net_neg_mode: str = "etp",
) -> tuple[pd.DataFrame, int, int]:
    """
    Add an 'ACT. ETRA_patched' column to df_et, flooring the net-forcing
    artefact at ETp. The original 'ACT. ETRA' column is left untouched.

    Parameters
    ----------
    df_et : DataFrame from simu.read_outputs("ET") for one month, with
        columns 'time_sec', 'SURFACE NODE', 'ACT. ETRA'.
    t_month : time axis (seconds from month start) passed to
        simu.update_atmbc() for this month.
    ETp_nodes, rain_nodes : ndarray (n_time_full, n_surf_nodes), m/s —
        the full-run forcing arrays built in step 9.
    idx_month : ndarray
        Indices into ETp_nodes/rain_nodes' time axis for this month
        (same mask used to build t_month/nv_month via
        subset_atmbc_for_month).
    eta_eps : float
        ACT. ETRA magnitude (m/s) at/below which a value is treated as
        CATHY's "collapsed to 0" artefact rather than a real reading.
    net_neg_mode : {"etp", "gross", "off"}
        Treatment of cells where net = rain - ETp < 0 (ETp exceeds rain).
        CATHY only receives the net flux, so there its ACT. ETRA is the
        *net* evaporation (~ ETp - rain at best), i.e. it under-reports the
        gross ET by the rain that evaporated first.
          "etp"   (default) patched ETa = ETp on every net < 0 cell.
          "gross" patched ETa = min(ETp, ACT. ETRA + rain): equals ETp
                  whenever the solver met the net demand, but keeps
                  genuinely soil-limited cells below ETp.
          "off"   leave net < 0 cells at the raw ACT. ETRA.
        Cells with net >= 0 are always handled by the artefact rule above.

    Returns
    -------
    (df_et, n_artefact, n_net_neg) : patched copy of df_et, the number of
    net >= 0 artefact cells patched, and the number of net < 0 cells whose
    value was changed by *net_neg_mode*.
    """
    if net_neg_mode not in ("etp", "gross", "off"):
        raise ValueError(f"net_neg_mode must be etp|gross|off, got {net_neg_mode!r}")
    if len(df_et) == 0 or len(t_month) == 0:
        df_et = df_et.copy()
        df_et["ACT. ETRA_patched"] = df_et.get("ACT. ETRA", pd.Series(dtype=float))
        return df_et, 0, 0

    t_month_arr = np.asarray(t_month, dtype=float)
    tsec        = df_et["time_sec"].to_numpy(dtype=float)

    # Nearest-time match (defensive against float round-trip through CATHY I/O)
    i_right  = np.clip(np.searchsorted(t_month_arr, tsec), 0, len(t_month_arr) - 1)
    i_left   = np.clip(i_right - 1, 0, len(t_month_arr) - 1)
    use_left = np.abs(t_month_arr[i_left] - tsec) < np.abs(t_month_arr[i_right] - tsec)
    i_time   = np.where(use_left, i_left, i_right)

    node_idx = df_et["SURFACE NODE"].to_numpy(dtype=int)
    n_surf   = ETp_nodes.shape[1]
    # CATHY node numbering is 1-based (Fortran) — auto-detect and correct
    if node_idx.max() == n_surf:
        node_idx = node_idx - 1
    assert node_idx.min() >= 0 and node_idx.max() < n_surf, (
        f"SURFACE NODE out of bounds for ETp_nodes (n_surf={n_surf}); "
        f"got range [{node_idx.min()}, {node_idx.max()}]"
    )

    rows      = idx_month[i_time]
    etp_vals  = ETp_nodes[rows, node_idx]    # m/s
    rain_vals = rain_nodes[rows, node_idx]   # m/s
    net_vals  = rain_vals - etp_vals

    eta_raw     = df_et["ACT. ETRA"].to_numpy(dtype=float)
    is_artefact = (np.abs(eta_raw) <= eta_eps) & (net_vals >= 0)
    eta_patched = np.where(is_artefact, etp_vals, eta_raw)

    # net < 0 (ETp > rain): patched ETa must equal ETp (see net_neg_mode).
    net_neg = net_vals < 0
    if net_neg_mode == "etp":
        eta_patched = np.where(net_neg, etp_vals, eta_patched)
    elif net_neg_mode == "gross":
        gross = np.minimum(np.abs(eta_raw) + rain_vals, etp_vals)
        eta_patched = np.where(net_neg, gross, eta_patched)
    n_net_neg = int((net_neg & ~np.isclose(eta_patched, eta_raw, rtol=0, atol=1e-15)).sum())

    df_et = df_et.copy()
    df_et["ACT. ETRA_patched"] = eta_patched
    return df_et, int(is_artefact.sum()), n_net_neg


def add_datetime_col(df: pd.DataFrame, base_dt) -> pd.DataFrame:
    """
    Add a 'datetime' column as base_dt + time_sec (returns df).

    Guards against a pyCATHY output-reading quirk where 'time_sec' is
    sometimes already a Timedelta dtype (instead of a plain float number
    of seconds) — calling pd.to_timedelta(..., unit="s") on data that's
    already Timedelta makes pandas try Timedelta+float arithmetic
    internally and raises "unsupported operand type(s) for +: 'Timedelta'
    and 'float'". When time_sec is already a Timedelta we use it as-is;
    otherwise we coerce to float first and apply unit="s" as before.
    """
    ts = df["time_sec"]
    if pd.api.types.is_timedelta64_dtype(ts):
        df["datetime"] = base_dt + ts
    else:
        df["datetime"] = base_dt + pd.to_timedelta(ts.astype(float), unit="s")
    return df


def check_monthly_sanity(
    dt,
    month_stats: dict,
    psi_month_df,
    sw_month_df,
    df_et,
    df_rec,
    t_month,
    idx_month,
    ETp_nodes,
    rain_nodes,
    sw_tol: float = SW_TOL,
    psi_sanity_min: float = PSI_SANITY_MIN,
    eta_tol: float = ETA_TOL,
    wb_residual_frac: float = WB_RESIDUAL_FRAC,
) -> list[str]:
    """
    Run cheap physical-bounds and (approximate) water-balance sanity
    checks on one month's outputs. Fills diagnostic fields into
    `month_stats` in place and returns a list of human-readable warning
    strings (empty if nothing looked off). Never raises — a check that
    can't be computed for this month (missing data, unexpected columns)
    is silently skipped rather than aborting the run.

    Water-balance note
    -------------------
    The residual is rain_in − ETa_out − recharge_out only: it does NOT
    account for the change in subsurface storage (ΔS), because that
    would require per-node porosity and control-volume information not
    readily available here. So this is not a strict mass-balance
    closure — treat a large residual as "worth checking the psi/sw
    trend for this month", not as proof of a bug on its own. It also
    assumes the 'recharge' output column(s) are in the same units
    (m/s) as ETa/rain; if your CATHY build reports recharge in
    different units, `recharge_out_mm` below will be wrong and should
    be ignored.
    """
    warnings: list[str] = []

    # ── Physical bounds: SW should stay within [0, 1] ──
    if sw_month_df is not None and len(sw_month_df):
        sw_vals = sw_month_df["sw"].to_numpy(dtype=float)
        n_bad = int(((sw_vals < -sw_tol) | (sw_vals > 1 + sw_tol)).sum())
        month_stats["sw_out_of_bounds_n"] = n_bad
        if n_bad:
            warnings.append(
                f"{dt:%Y-%m}: {n_bad} SW value(s) outside [0,1] "
                f"(range [{sw_vals.min():.3f}, {sw_vals.max():.3f}])."
            )

    # ── Physical bounds: PSI shouldn't be sitting at extreme drying values ──
    if psi_month_df is not None and len(psi_month_df):
        psi_vals = psi_month_df["psi"].to_numpy(dtype=float)
        n_extreme = int((psi_vals < psi_sanity_min).sum())
        month_stats["psi_extreme_n"] = n_extreme
        if n_extreme:
            warnings.append(
                f"{dt:%Y-%m}: {n_extreme} PSI value(s) below sanity threshold "
                f"{psi_sanity_min:.0f} m (min {psi_vals.min():.2f} m)."
            )

    # ── Physical bounds: ETa should not exceed ETp ──
    if (df_et is not None and len(df_et) and t_month is not None and len(t_month)
            and ETp_nodes is not None and idx_month is not None and len(idx_month)):
        try:
            t_month_arr = np.asarray(t_month, dtype=float)
            tsec = df_et["time_sec"].to_numpy(dtype=float)
            i_right  = np.clip(np.searchsorted(t_month_arr, tsec), 0, len(t_month_arr) - 1)
            i_left   = np.clip(i_right - 1, 0, len(t_month_arr) - 1)
            use_left = np.abs(t_month_arr[i_left] - tsec) < np.abs(t_month_arr[i_right] - tsec)
            i_time   = np.where(use_left, i_left, i_right)

            node_idx = df_et["SURFACE NODE"].to_numpy(dtype=int)
            n_surf   = ETp_nodes.shape[1]
            if node_idx.max() == n_surf:
                node_idx = node_idx - 1

            rows     = idx_month[i_time]
            etp_vals = ETp_nodes[rows, node_idx]
            eta_col  = "ACT. ETRA_patched" if "ACT. ETRA_patched" in df_et.columns else "ACT. ETRA"
            eta_vals = df_et[eta_col].to_numpy(dtype=float)

            n_exceed = int((eta_vals > etp_vals + eta_tol).sum())
            month_stats["eta_exceeds_etp_n"] = n_exceed
            if n_exceed:
                max_excess = float((eta_vals - etp_vals).max())
                warnings.append(
                    f"{dt:%Y-%m}: {n_exceed} ETa value(s) exceed ETp "
                    f"(max excess {max_excess:.3e} m/s)."
                )
        except Exception:
            pass  # sanity check itself must never break the run

    # ── Water balance: rain_in vs ETa_out vs recharge_out (flux-only) ──
    rain_in_mm = eta_out_mm = recharge_out_mm = None

    if rain_nodes is not None and idx_month is not None and len(idx_month):
        rain_in_mm = float((rain_nodes[idx_month].mean(axis=1) * 86400 * 1000).sum())
        month_stats["rain_in_mm"] = rain_in_mm

    if df_et is not None and len(df_et):
        eta_col = "ACT. ETRA_patched" if "ACT. ETRA_patched" in df_et.columns else "ACT. ETRA"
        if eta_col in df_et.columns:
            eta_by_time = df_et.groupby("time_sec")[eta_col].mean()
            eta_out_mm = float((eta_by_time * 86400 * 1000).sum())
            month_stats["eta_out_mm"] = eta_out_mm

    if df_rec is not None and len(df_rec):
        rec_value_cols = [c for c in df_rec.select_dtypes(include=[np.number]).columns
                           if c not in ("X", "Y", "time_sec")]
        if rec_value_cols:
            rec_by_time = df_rec.groupby("time_sec")[rec_value_cols].mean().mean(axis=1)
            recharge_out_mm = float((rec_by_time * 86400 * 1000).sum())
            month_stats["recharge_out_mm"] = recharge_out_mm

    if rain_in_mm is not None and eta_out_mm is not None and recharge_out_mm is not None:
        residual_mm = rain_in_mm - eta_out_mm - recharge_out_mm
        month_stats["water_balance_residual_mm"] = residual_mm
        if rain_in_mm > 1e-6 and abs(residual_mm) > wb_residual_frac * rain_in_mm:
            warnings.append(
                f"{dt:%Y-%m}: water balance residual {residual_mm:.1f} mm "
                f"(rain {rain_in_mm:.1f}, ETa {eta_out_mm:.1f}, "
                f"recharge {recharge_out_mm:.1f}) — large relative to rain "
                f"input. Excludes storage change (ΔS); also check the "
                f"psi/sw trend for this month before assuming a bug."
            )

    return warnings
