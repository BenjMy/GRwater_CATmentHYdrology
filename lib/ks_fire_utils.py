"""
ks_fire_utils.py — Post-fire saturated hydraulic conductivity (Ks) recovery
============================================================================
GRwater project — Agramon study site

Small, dependency-light module implementing the literature equations for
how fire-affected surface-soil Ks recovers with time and/or vegetation
regrowth, so ``Agramon_withLAI_withETp.py`` can drive its ``PERMX`` override
from a physically-grounded, published relationship instead of a single
static ``--permx`` value.

Three independent "modules" (equation families) are provided, matching the
three published approaches this was built from — pick whichever fits the
data you actually have for a given run:

1. ``time``          — Ks recovers as a logistic function of time-since-fire
                        alone (Ebel & Martin, 2017, *Hydrological Processes*
                        31(21):3682-3696; the same functional form is used
                        by Thomas et al., 2021, their Eq. 9).
2. ``lai``            — Ks is linearly rescaled from an *observed* LAI value
                        each month, anchored between a post-fire minimum and
                        a pre-fire/recovered maximum (Thomas et al., 2021,
                        *JGR: Earth Surface* 126:e2021JF006091, their Eq. 11,
                        re-derived below from first principles — see
                        ``ks_from_lai_linear`` docstring for the derivation).
                        This is the natural fit for this pipeline because it
                        already loads a real monthly LAI raster — no LAI
                        curve needs to be *fitted*, the observed value is
                        used directly.
3. ``burn_severity``  — A single, static Ks reduction scaled by a burn-
                        severity fraction (dNBR-style), in the spirit of
                        Atchley et al. (2018, *Vadose Zone Journal* 17(1):
                        1-13), who reduced Ks (and LAI) in the top 2 cm of
                        soil as a function of differenced Normalized Burn
                        Ratio for a single pre-fire vs. post-fire comparison
                        (no multi-year recovery trajectory). NOTE: Atchley
                        et al.'s own dNBR→Ks relation is imported from a
                        separate source they cite and is not reproduced
                        here — ``ks_burn_severity_step`` is a simple,
                        clearly-labelled linear placeholder for that slot;
                        replace ``severity_to_fraction`` with your own
                        calibrated relation if you have one.

All Ks values are handled internally in m/s (pyCATHY's PERMX unit); helper
constants below are given in mm/hr (as published) and converted for you.

Typical use from the main pipeline
-----------------------------------
    import ks_fire_utils as KsFire

    cfg = KsFire.KsFireConfig(
        mode="lai",
        ks_ll=KsFire.mm_hr_to_m_s(16.2),   # post-fire minimum Ks
        ks_ul=KsFire.mm_hr_to_m_s(43.8),   # recovered/unburned Ks
        lai_ll=0.05,                       # LAI observed right after the fire
        lai_ul=0.55,                       # pre-fire / fully-recovered LAI
    )

    # ... inside the monthly loop, once you have this month's `lai_2d` ...
    permx_month = KsFire.compute_permx_for_month(
        cfg, current_date=dt, lai_value=float(np.nanmean(lai_2d)),
    )
    update_soil_for_active_zones(..., permx=permx_month or args.permx, ...)

See ``add_ks_fire_args`` / ``ks_fire_config_from_args`` for wiring this up
as a CLI-configurable feature (matches the ``--permx``-group style already
used in ``Agramon_withLAI_withETp.py``).
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from datetime import date, datetime

# ── Unit conversion ──────────────────────────────────────────────────────
# Published Ks-recovery values are almost always given in mm/hr;
# pyCATHY's PERMX (SOIL SPP) expects m/s.

_MM_HR_TO_M_S = 1e-3 / 3600.0  # 1 mm/hr = 2.77778e-7 m/s


def mm_hr_to_m_s(value_mm_hr: float) -> float:
    """Convert a Ks value from mm/hr (as published) to m/s (PERMX unit)."""
    return value_mm_hr * _MM_HR_TO_M_S


def m_s_to_mm_hr(value_m_s: float) -> float:
    """Convert a Ks value from m/s (PERMX unit) back to mm/hr (for logging)."""
    return value_m_s / _MM_HR_TO_M_S


def years_since_fire(current_date: date | datetime, fire_date: date | datetime) -> float:
    """
    Fractional years elapsed between ``fire_date`` and ``current_date``
    (negative before the fire; callers generally clip time-based recovery
    at 0 — see ``ks_recovery_time_based``).
    """
    delta_days = (
        (current_date.date() if isinstance(current_date, datetime) else current_date)
        - (fire_date.date() if isinstance(fire_date, datetime) else fire_date)
    ).days
    return delta_days / 365.25


# ── Literature presets (Thomas et al., 2021, Table 1) ───────────────────
# Field-based fits of Eq. 9 to repeat tension-infiltrometer Ks data for the
# 2016 Fish Fire (San Gabriel Mountains, loamy sand). Values in mm/hr / yr-1.
# Keyed exactly as reported: upper/lower bound of the boxplot fits, the
# median fit, and the "overall average" of the three (their headline value).
THOMAS2021_FIELD_KS_PRESETS_MM_HR: dict[str, dict[str, float]] = {
    "upper_bound":     {"ks_ll": 25.5, "ks_ul": 87.5, "r": 4.86},
    "median":          {"ks_ll": 12.9, "ks_ul": 19.5, "r": 13.63},
    "lower_bound":     {"ks_ll": 0.5,  "ks_ul": 5.0,  "r": 3.21},
    "overall_average": {"ks_ll": 13.0, "ks_ul": 37.5, "r": 4.47},
}

# Satellite-LAI-based Ks anchor points Thomas et al. (2021) used for their
# LAI-driven Ks curve (Eq. 11) at the Arroyo Seco / Dunsmore sites: unburned
# Ks taken as the Rosetta loamy-sand default, post-fire minimum from a
# burned:unburned ratio of 0.37 (Ebel & Moody, 2020). No fitted LAI_LL/
# LAI_UL/r_LAI values are given in the source for general reuse — pull
# those directly from your own LAI record (see ``mode="lai"`` above).
THOMAS2021_SATELLITE_KS_ANCHORS_MM_HR: dict[str, float] = {
    "ks_ll": 16.2,  # post-fire minimum
    "ks_ul": 43.8,  # unburned / recovered (Rosetta loamy sand)
}


def get_field_preset_m_s(name: str = "overall_average") -> dict[str, float]:
    """
    Return one of ``THOMAS2021_FIELD_KS_PRESETS_MM_HR`` converted to m/s,
    ready to pass as ``ks_ll``/``ks_ul``/``r`` to ``ks_recovery_time_based``.
    """
    if name not in THOMAS2021_FIELD_KS_PRESETS_MM_HR:
        raise KeyError(
            f"Unknown preset '{name}'. Choices: "
            f"{sorted(THOMAS2021_FIELD_KS_PRESETS_MM_HR)}"
        )
    preset = THOMAS2021_FIELD_KS_PRESETS_MM_HR[name]
    return {
        "ks_ll": mm_hr_to_m_s(preset["ks_ll"]),
        "ks_ul": mm_hr_to_m_s(preset["ks_ul"]),
        "r": preset["r"],
    }


# ── Module 1: generic logistic recovery (Ebel & Martin, 2017) ───────────

def logistic_recovery(t_years: float, lower: float, upper: float, rate: float) -> float:
    """
    Generic logistic recovery curve shared by both the time-based Ks curve
    (Thomas et al., 2021, Eq. 9) and the LAI recovery curve (their Eq. 10):

        value(t) = lower + (upper - lower) * e^(rate * t) / (1 + e^(rate * t))

    ``t_years`` is time since fire (years); values before the fire
    (t_years < 0) are clipped to ``lower`` by the callers below rather than
    extrapolated, since the logistic form is only fit/valid post-fire.
    """
    # Guard against overflow for large positive rate*t (e^x with big x):
    # rewrite the logistic term as 1 / (1 + e^(-rate*t)) which is numerically
    # stable for both large positive and large negative exponents.
    frac = 1.0 / (1.0 + math.exp(-rate * t_years))
    return lower + (upper - lower) * frac


def ks_recovery_time_based(
    t_years: float,
    ks_ll: float,
    ks_ul: float,
    r_ks: float,
) -> float:
    """
    Thomas et al. (2021), Eq. 9 / Ebel & Martin (2017): Ks as a function of
    time since fire alone, all in m/s. Ks is held at ``ks_ll`` for
    ``t_years <= 0`` (i.e. no recovery before the fire happens).
    """
    if t_years <= 0:
        return ks_ll
    return logistic_recovery(t_years, ks_ll, ks_ul, r_ks)


def lai_recovery_time_based(
    t_years: float,
    lai_ll: float,
    lai_ul: float,
    r_lai: float,
) -> float:
    """
    Thomas et al. (2021), Eq. 10: LAI as a function of time since fire
    alone (their fit to annual postwildfire LAI maxima). Only useful if you
    don't have real LAI observations for the period you're simulating and
    need a forecast instead — this pipeline normally has real LAI, so
    prefer ``ks_from_lai_linear`` with the observed monthly value.
    """
    if t_years <= 0:
        return lai_ll
    return logistic_recovery(t_years, lai_ll, lai_ul, r_lai)


# ── Module 2: LAI-driven Ks (Thomas et al., 2021, Eq. 11) ───────────────

def ks_from_lai_linear(
    lai_value: float,
    lai_ll: float,
    lai_ul: float,
    ks_ll: float,
    ks_ul: float,
    clip: bool = True,
) -> float:
    """
    Map an *observed* LAI value onto Ks via the same linear rescaling
    Thomas et al. (2021) use in their Eq. 11.

    Derivation: Eq. 9 and Eq. 10 are both the same logistic shape function
    frac(t) = e^(r t)/(1+e^(r t)), just rescaled onto different [lower,
    upper] ranges. If LAI(t) = lai_ll + (lai_ul-lai_ll)*frac(t), then
    frac(t) = (LAI(t) - lai_ll) / (lai_ul - lai_ll), so substituting into
    Ks(t) = ks_ll + (ks_ul-ks_ll)*frac(t) gives a plain linear map from an
    LAI value to a Ks value — no curve-fitting required, only the four
    anchor points (lai_ll, lai_ul, ks_ll, ks_ul). This lets a *real*,
    already-loaded monthly LAI raster mean drive Ks directly, instead of
    needing to separately fit Eq. 10 to a time series first.

    ``clip=True`` (default) clamps the result to [ks_ll, ks_ul] so noisy or
    out-of-range LAI values (e.g. a wetter-than-reference spring) can't
    push Ks outside the physically-anchored recovery envelope.
    """
    if lai_ul == lai_ll:
        raise ValueError("lai_ul must differ from lai_ll (division by zero).")
    frac = (lai_value - lai_ll) / (lai_ul - lai_ll)
    if clip:
        frac = min(max(frac, 0.0), 1.0)
    return ks_ll + (ks_ul - ks_ll) * frac


# ── Module 3: static burn-severity step (Atchley et al., 2018-style) ────

def severity_to_fraction(dnbr: float, dnbr_max: float = 1300.0) -> float:
    """
    Simple linear normalisation of a dNBR value to a [0, 1] burn-severity
    fraction (0 = unburned, 1 = most severe in your dataset/domain).
    ``dnbr_max`` defaults to a generic high-severity dNBR ceiling; pass
    your own scene-specific max (or a calibrated dNBR→Ks relation via
    ``ks_burn_severity_step``'s ``severity_fraction`` argument directly)
    for anything beyond a rough placeholder.
    """
    return min(max(dnbr / dnbr_max, 0.0), 1.0)


def ks_burn_severity_step(
    ks_unburned: float,
    severity_fraction: float,
    max_reduction: float = 0.8,
) -> float:
    """
    Single, static Ks reduction scaled by burn severity — in the spirit of
    Atchley et al. (2018), who reduced Ks (and LAI) in the top 2 cm of soil
    as a function of dNBR for one pre-fire vs. post-fire comparison (no
    multi-year recovery trajectory; this is NOT the same concept as the
    time/LAI recovery curves above — see the module docstring).

    ``max_reduction`` is the fractional Ks reduction at severity_fraction
    == 1 (default 0.8, i.e. Ks drops to 20% of unburned at max severity —
    a generic placeholder, not Atchley et al.'s own fitted coefficients,
    which come from a separate source they cite and aren't reproduced
    here). Replace this function's body with your own calibrated dNBR→Ks
    relation if you have one; the signature is kept simple on purpose so
    it's a drop-in swap.
    """
    severity_fraction = min(max(severity_fraction, 0.0), 1.0)
    return ks_unburned * (1.0 - max_reduction * severity_fraction)


# ── Config + dispatcher ──────────────────────────────────────────────────

@dataclass
class KsFireConfig:
    """
    Everything needed to compute a fire-affected PERMX override for a
    given month, for whichever ``mode`` you choose. Unused fields for a
    given mode are simply ignored (kept optional so one config object
    covers all three modes without subclassing).
    """
    mode: str = "off"  # "off" | "time" | "lai" | "burn_severity"

    # mode="time"
    fire_date: date | None = None
    ks_ll: float | None = None       # m/s
    ks_ul: float | None = None       # m/s
    r_ks: float | None = None        # yr^-1

    # mode="lai" (reuses ks_ll / ks_ul above as the Ks anchors)
    lai_ll: float | None = None
    lai_ul: float | None = None

    # mode="burn_severity"
    ks_unburned: float | None = None  # m/s
    severity_fraction: float | None = None  # 0-1, pre-computed by caller
    max_reduction: float = 0.8


def compute_permx_for_month(
    cfg: KsFireConfig,
    current_date: date | datetime,
    lai_value: float | None = None,
) -> float | None:
    """
    Single entry point the pipeline calls once per month. Returns a PERMX
    value in m/s, or None if ``cfg.mode == "off"`` (caller should then fall
    back to its own static ``--permx``, if any).
    """
    if cfg.mode == "off":
        return None

    if cfg.mode == "time":
        if cfg.fire_date is None or cfg.ks_ll is None or cfg.ks_ul is None or cfg.r_ks is None:
            raise ValueError(
                "mode='time' requires fire_date, ks_ll, ks_ul, and r_ks to be set."
            )
        t_years = years_since_fire(current_date, cfg.fire_date)
        return ks_recovery_time_based(t_years, cfg.ks_ll, cfg.ks_ul, cfg.r_ks)

    if cfg.mode == "lai":
        if lai_value is None:
            raise ValueError("mode='lai' requires lai_value (this month's observed LAI).")
        if None in (cfg.lai_ll, cfg.lai_ul, cfg.ks_ll, cfg.ks_ul):
            raise ValueError(
                "mode='lai' requires lai_ll, lai_ul, ks_ll, and ks_ul to be set."
            )
        return ks_from_lai_linear(lai_value, cfg.lai_ll, cfg.lai_ul, cfg.ks_ll, cfg.ks_ul)

    if cfg.mode == "burn_severity":
        if cfg.ks_unburned is None or cfg.severity_fraction is None:
            raise ValueError(
                "mode='burn_severity' requires ks_unburned and severity_fraction."
            )
        return ks_burn_severity_step(
            cfg.ks_unburned, cfg.severity_fraction, cfg.max_reduction
        )

    raise ValueError(f"Unknown KsFireConfig.mode: '{cfg.mode}'")


# ── Optional CLI wiring (matches the --permx argument-group style) ──────

def add_ks_fire_args(parser: argparse.ArgumentParser) -> None:
    """
    Add a ``--ks-fire-*`` argument group to an existing argparse parser,
    mirroring the style of the ``--permx`` group in
    ``Agramon_withLAI_withETp.py``. Call this from ``_parse_args()`` and
    then build a ``KsFireConfig`` from the result with
    ``ks_fire_config_from_args``.
    """
    g = parser.add_argument_group(
        "Fire-affected Ks recovery (time/LAI/burn-severity, overrides --permx per month)"
    )
    g.add_argument(
        "--ks-fire-mode", default="off",
        choices=["off", "time", "lai", "burn_severity"],
        help="How PERMX is computed each month from post-fire recovery "
             "literature instead of a single static --permx. 'off' "
             "(default) leaves --permx as-is. See ks_fire_utils.py for "
             "the equations behind each mode."
    )
    g.add_argument(
        "--ks-fire-preset", default=None, metavar="NAME",
        choices=sorted(THOMAS2021_FIELD_KS_PRESETS_MM_HR),
        help="Shortcut: load ks_ll/ks_ul/r_ks from Thomas et al. (2021) "
             "Table 1 (field-based fits for the 2016 Fish Fire) instead of "
             "specifying --ks-fire-ll/--ks-fire-ul/--ks-fire-rate by hand. "
             "Only used when --ks-fire-mode=time."
    )
    g.add_argument(
        "--ks-fire-date", default=None, metavar="YYYY-MM-DD",
        help="Calendar date of the fire (ISO format). Required for "
             "--ks-fire-mode=time."
    )
    g.add_argument(
        "--ks-fire-ll", type=float, default=None, metavar="M/S",
        help="Post-fire minimum Ks (m/s). Required for --ks-fire-mode="
             "time/lai unless --ks-fire-preset is given."
    )
    g.add_argument(
        "--ks-fire-ul", type=float, default=None, metavar="M/S",
        help="Recovered/unburned Ks (m/s). Required for --ks-fire-mode="
             "time/lai unless --ks-fire-preset is given."
    )
    g.add_argument(
        "--ks-fire-rate", type=float, default=None, metavar="YR^-1",
        help="Logistic recovery rate (yr^-1). Required for --ks-fire-mode="
             "time unless --ks-fire-preset is given."
    )
    g.add_argument(
        "--ks-fire-lai-ll", type=float, default=None, metavar="LAI",
        help="LAI value anchoring ks-fire-ll (typically the LAI observed "
             "right after the fire). Required for --ks-fire-mode=lai."
    )
    g.add_argument(
        "--ks-fire-lai-ul", type=float, default=None, metavar="LAI",
        help="LAI value anchoring ks-fire-ul (typically pre-fire / fully-"
             "recovered LAI). Required for --ks-fire-mode=lai."
    )
    g.add_argument(
        "--ks-fire-severity", type=float, default=None, metavar="0-1",
        help="Static burn-severity fraction (0=unburned, 1=max severity). "
             "Required for --ks-fire-mode=burn_severity."
    )
    g.add_argument(
        "--ks-fire-severity-reduction", type=float, default=0.8, metavar="FRAC",
        help="Fractional Ks reduction at severity=1 for --ks-fire-mode="
             "burn_severity (default: 0.8, i.e. Ks drops to 20%% of "
             "unburned at max severity)."
    )


def ks_fire_config_from_args(args: argparse.Namespace) -> KsFireConfig:
    """
    Build a ``KsFireConfig`` from the ``--ks-fire-*`` args added by
    ``add_ks_fire_args``. Applies ``--ks-fire-preset`` (if given) before
    any explicit --ks-fire-ll/--ks-fire-ul/--ks-fire-rate, so a preset can
    still be fine-tuned by also passing an explicit override.
    """
    mode = args.ks_fire_mode
    if mode == "off":
        return KsFireConfig(mode="off")

    ks_ll, ks_ul, r_ks = args.ks_fire_ll, args.ks_fire_ul, args.ks_fire_rate
    if args.ks_fire_preset is not None:
        preset = get_field_preset_m_s(args.ks_fire_preset)
        ks_ll = ks_ll if ks_ll is not None else preset["ks_ll"]
        ks_ul = ks_ul if ks_ul is not None else preset["ks_ul"]
        r_ks = r_ks if r_ks is not None else preset["r"]

    fire_date = None
    if args.ks_fire_date is not None:
        fire_date = datetime.strptime(args.ks_fire_date, "%Y-%m-%d").date()

    return KsFireConfig(
        mode=mode,
        fire_date=fire_date,
        ks_ll=ks_ll,
        ks_ul=ks_ul,
        r_ks=r_ks,
        lai_ll=args.ks_fire_lai_ll,
        lai_ul=args.ks_fire_lai_ul,
        ks_unburned=ks_ul,  # burn_severity mode treats ks_ul as "unburned"
        severity_fraction=args.ks_fire_severity,
        max_reduction=args.ks_fire_severity_reduction,
    )
