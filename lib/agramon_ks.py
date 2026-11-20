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

FIRE_DATE_DEFAULT = "2020-07-01"     # Agramon fire (TO-DO.md; FIRE_DATE in plot_results_shared)

# CATHY vertical discretisation defaults (pyCATHY update_dem_parameters): 15 layers whose
# thicknesses are fractions (zratio, sum = 1) of the total column depth `base`.
CATHY_BASE_DEFAULT = 3.0
CATHY_ZRATIO_DEFAULT = (0.002, 0.004, 0.006, 0.008, 0.01, 0.01, 0.02, 0.02, 0.05, 0.05,
                        0.1, 0.1, 0.2, 0.2, 0.22)


def _as_date(d):
    return d.date() if isinstance(d, datetime) else d


def ks_config_from_args(args):
    """KsFireConfig from the --ks-fire-* args; validates the mode's required inputs now
    (not at the first month) by evaluating the config once."""
    cfg = KsFire.ks_fire_config_from_args(args)
    if cfg.mode != "off":
        if cfg.ks_ul is None:
            raise ValueError("--ks-fire-mode needs --ks-fire-ul (or --ks-fire-preset)")
        if cfg.mode == "lai" and (cfg.lai_ll is None or cfg.lai_ul is None):
            raise ValueError("--ks-fire-mode lai needs --ks-fire-lai-ll and --ks-fire-lai-ul "
                             "(LAI right after the fire / fully recovered LAI)")
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
    layers = (f"ks2L-d{args.ks_layer_depth:g}m-b{args.permx_bottom:.0e}".replace("-0", "-")
              if two_layer_enabled(args) else "")
    if mode == "off":
        return layers
    preset = getattr(args, "ks_fire_preset", None)
    return "_".join(t for t in (f"ks-{mode}" + (f"-{preset}" if preset else ""), layers) if t)


def ks_key(args) -> str:
    """Stable string identifying the post-fire Ks setup of a run ('off' if disabled).
    Logged in the CSV run log (column ks_fire) and compared by the shell runner's
    skip-if-done check, so Ks scenarios never satisfy each other."""
    mode = getattr(args, "ks_fire_mode", "off")
    layers = ks_layers_key(args)
    if mode == "off":
        return f"off|{layers}" if layers else "off"
    g = lambda k: getattr(args, k, None)
    num = lambda v: "" if v is None else f"{v:g}"
    return "|".join([mode, f"preset={g('ks_fire_preset') or ''}",
                     f"date={g('ks_fire_date') or ''}",
                     f"ll={num(g('ks_fire_ll'))}", f"ul={num(g('ks_fire_ul'))}",
                     f"rate={num(g('ks_fire_rate'))}",
                     f"lai_ll={num(g('ks_fire_lai_ll'))}", f"lai_ul={num(g('ks_fire_lai_ul'))}",
                     f"sev={num(g('ks_fire_severity'))}",
                     f"red={num(g('ks_fire_severity_reduction'))}"] + ([layers] if layers else []))


def ks_key_from_argv(argv) -> str:
    """ks_key of a list of --ks-fire-* flags (used by the shell runner's check)."""
    import argparse
    p = argparse.ArgumentParser()
    KsFire.add_ks_fire_args(p)
    p.set_defaults(ks_fire_date=FIRE_DATE_DEFAULT)
    return ks_key(p.parse_args(list(argv)))


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


# ── Two-layer Ks (top / bottom) ───────────────────────────────────────────────
# Opt-in with --permx-bottom. The top layer keeps --permx (or the post-fire Ks of the month),
# the layers below --ks-layer-depth get --permx-bottom. A CATHY layer belongs to the top
# soil when its CENTRE is shallower than the interface (layers are not split).

def layers_from_simu(simu):
    """(zratio, base) of the CATHY mesh from ``simu.dem_parameters`` (pyCATHY defaults if absent)."""
    dp = getattr(simu, "dem_parameters", None) or {}
    zr = dp.get("zratio(i),i=1,nstr", CATHY_ZRATIO_DEFAULT)
    base = dp.get("base", CATHY_BASE_DEFAULT)
    if isinstance(zr, str):                  # pyCATHY may store it as "0.1 0.2 ..." (or comma separated)
        zr = zr.replace(",", " ").split()
    elif np.isscalar(zr):
        zr = [zr]
    if isinstance(base, (list, tuple)):
        base = base[0]
    return [float(v) for v in zr], float(str(base).split()[0]) if isinstance(base, str) else float(base)


def layer_centres(zratio, base):
    """Depth (m below the surface) of the centre of each CATHY layer, top first."""
    z = np.asarray(zratio, float) * float(base)
    return np.cumsum(z) - z / 2.0


def layer_permx_profile(zratio, base, top_depth, ks_top, ks_bottom):
    """PERMX (m/s) of every CATHY layer for a two-layer soil: ks_top above ``top_depth``
    (layer centre rule), ks_bottom below. Returns an array of len(zratio), top layer first."""
    if top_depth is None or top_depth <= 0:
        raise ValueError("--ks-layer-depth must be > 0 m")
    if float(sum(zratio)) <= 0 or abs(float(sum(zratio)) - 1.0) > 1e-2:
        raise ValueError(f"zratio must sum to 1 (got {sum(zratio):.3f})")
    return np.where(layer_centres(zratio, base) < float(top_depth), float(ks_top), float(ks_bottom))


def two_layer_enabled(args) -> bool:
    return getattr(args, "permx_bottom", None) is not None


def validate_two_layer_args(args) -> None:
    """Fail at start-up (not at month 1) when the two-layer options are inconsistent."""
    if not two_layer_enabled(args):
        return
    if args.permx is None:
        raise ValueError("--permx-bottom needs --permx (Ks of the top layer)")
    if args.permx_bottom <= 0:
        raise ValueError("--permx-bottom must be > 0")
    if getattr(args, "ks_layer_depth", None) is None or args.ks_layer_depth <= 0:
        raise ValueError("--permx-bottom needs --ks-layer-depth > 0 (m, top-layer thickness)")


def ks_layers_key(args) -> str:
    """'' when single-layer, else a stable 'L2=<depth>:<bottom>' fragment."""
    if not two_layer_enabled(args):
        return ""
    return f"L2={args.ks_layer_depth:g}:{args.permx_bottom:g}"
