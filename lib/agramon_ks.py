"""
agramon_ks.py - post-fire Ks (PERMX) per month for the Agramon pipeline
=======================================================================
Thin layer between the main script and utils/ks_fire_utils.py (science and
equations: TO-DO.md 5.3). Opt-in: with --ks-fire-mode off (default) every
function returns the static --permx, so runs are unchanged.

Scope: the returned PERMX replaces --permx for EVERY soil zone (same as --permx);
surface-layer-only Ks (Atchley et al. 2018: top 2 cm) is not implemented here.
Before the fire date the UNBURNED Ks (ks_ul) is used in all modes: the time-based
util alone would hold Ks at the burned minimum for t <= 0.
"""
from __future__ import annotations

from datetime import date, datetime

import numpy as np

from utils import ks_fire_utils as KsFire

FIRE_DATE_DEFAULT = "2020-07-01"     # Agramon fire (TO-DO.md; FIRE_DATE in results_plotter)


def _as_date(d):
    return d.date() if isinstance(d, datetime) else d


def ks_config_from_args(args):
    """KsFireConfig from the --ks-fire-* args; validates the mode's required inputs now
    (not at the first month) by evaluating the config once."""
    cfg = KsFire.ks_fire_config_from_args(args)
    if cfg.mode != "off":
        if cfg.ks_ul is None:
            raise ValueError("--ks-fire-mode needs --ks-fire-ul (or --ks-fire-preset)")
        probe_date = cfg.fire_date or date(2000, 1, 1)
        KsFire.compute_permx_for_month(cfg, probe_date, cfg.lai_ul)   # raises if incomplete
    return cfg


def ks_for_month(cfg, dt, lai_2d, static_permx):
    """PERMX (m/s) for month `dt`. `lai_2d`: that month's LAI map (domain mean drives
    mode 'lai'). Returns `static_permx` when mode is off."""
    if cfg.mode == "off":
        return static_permx
    unburned = cfg.ks_ul
    if cfg.fire_date is not None and _as_date(dt) < _as_date(cfg.fire_date):
        return unburned
    lai = float(np.nanmean(lai_2d)) if cfg.mode == "lai" else None
    value = KsFire.compute_permx_for_month(cfg, dt, lai)
    return static_permx if value is None else value


def ks_tag(args) -> str:
    """Folder-name tag so Ks scenarios never share a folder / 'already done' key."""
    mode = getattr(args, "ks_fire_mode", "off")
    if mode == "off":
        return ""
    preset = getattr(args, "ks_fire_preset", None)
    return f"ks-{mode}" + (f"-{preset}" if preset else "")


def ks_banner_label(cfg) -> str:
    if cfg.mode == "off":
        return "off (static --permx)"
    mm = KsFire.m_s_to_mm_hr
    parts = [f"mode={cfg.mode}"]
    if cfg.ks_ll is not None:
        parts.append(f"Ks {mm(cfg.ks_ll):.3g}->{mm(cfg.ks_ul):.3g} mm/hr")
    if cfg.mode == "time" and cfg.r_ks is not None:
        parts.append(f"rate {cfg.r_ks:g}/yr")
    if cfg.mode == "lai":
        parts.append(f"LAI {cfg.lai_ll}->{cfg.lai_ul}")
    if cfg.fire_date is not None:
        parts.append(f"fire {cfg.fire_date}")
    return " | ".join(parts)
