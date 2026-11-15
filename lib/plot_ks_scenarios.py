"""
plot_ks_scenarios.py — what the post-fire Ks scenarios look like (no CATHY needed)
==================================================================================
Helps choose and interpret the Ks hypotheses before running anything:

* ``ks_trajectories``  tidy table (date, preset, mode, ks_mm_hr, ratio) of the monthly
                       Ks that each mode (time / lai / burn_severity) would feed to CATHY
* ``plot_ks_scenarios`` one panel per literature preset, log-y, with horizontal lines for
                       reference rain intensities. Ks far above the rain intensity means Ks
                       cannot limit infiltration for that forcing; Ks below it means
                       infiltration-excess runoff is possible.

Rain intensities: a daily forcing of R mm/day corresponds to R/24 mm/hr averaged over the
day, which understates sub-daily intensity (daily forcing hides the peaks that control
infiltration excess); the lines are therefore lower bounds on the true intensity.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from utils import ks_fire_utils as KF

MODES = ("time", "lai", "burn_severity")
COLORS = {"time": "#0072B2", "lai": "#009E73", "burn_severity": "#D55E00"}


def ks_trajectories(months, fire_date, presets=None, lai_ll: float = 0.3, lai_ul: float = 2.0,
                    severity: float = 0.5, max_reduction: float = 0.8) -> pd.DataFrame:
    """
    Monthly Ks (mm/hr) for each preset and mode. ``lai`` mode is driven by a
    LAI that follows the same logistic recovery (``lai_ll`` -> ``lai_ul``, rate of
    the preset) — i.e. the idealised case where LAI and Ks recover together;
    ``burn_severity`` is the constant step ``ks_ul * (1 - severity*max_reduction)``
    applied from the fire date. Before the fire every mode gives ``ks_ul``.
    """
    if not 0 <= severity <= 1:
        raise ValueError("severity must be in [0, 1]")
    fire = pd.Timestamp(fire_date).date()
    names = list(presets) if presets else sorted(KF.THOMAS2021_FIELD_KS_PRESETS_MM_HR)
    rows = []
    for name in names:
        p = KF.get_field_preset_m_s(name)
        base = dict(ks_ll=p["ks_ll"], ks_ul=p["ks_ul"])
        cfgs = {
            "time": KF.KsFireConfig(mode="time", fire_date=fire, r_ks=p["r"], **base),
            "lai": KF.KsFireConfig(mode="lai", fire_date=fire, lai_ll=lai_ll, lai_ul=lai_ul, **base),
            "burn_severity": KF.KsFireConfig(mode="burn_severity", ks_unburned=p["ks_ul"],
                                             severity_fraction=severity, max_reduction=max_reduction),
        }
        for m in pd.DatetimeIndex(months):
            d = m.date()
            t = KF.years_since_fire(d, fire)
            lai = KF.lai_recovery_time_based(max(t, 0.0), lai_ll, lai_ul, p["r"]) if t >= 0 else lai_ul
            for mode, cfg in cfgs.items():
                if t < 0:
                    ks = p["ks_ul"]
                else:
                    ks = KF.compute_permx_for_month(cfg, d, lai if mode == "lai" else None)
                rows.append(dict(date=m, preset=name, mode=mode, ks_mm_hr=KF.m_s_to_mm_hr(ks),
                                 ratio=ks / p["ks_ul"]))
    return pd.DataFrame(rows)


def plot_ks_scenarios(df: pd.DataFrame, fire_date, rain_mm_day=(10.0, 40.0),
                      out: Optional[str] = None):
    """One log-y panel per preset; horizontal dotted lines = rain intensity (mm/hr, daily mean)."""
    import matplotlib.pyplot as plt

    need = {"date", "preset", "mode", "ks_mm_hr"}
    if not need <= set(df.columns):
        raise ValueError(f"df needs columns {sorted(need)}")
    presets = list(dict.fromkeys(df["preset"]))
    nc = 2 if len(presets) > 1 else 1
    nr = int(np.ceil(len(presets) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(5.6 * nc, 3.2 * nr), sharex=True, squeeze=False)
    for ax, name in zip(axes.ravel(), presets):
        d = df[df["preset"] == name]
        for mode in MODES:
            s = d[d["mode"] == mode]
            if len(s):
                wide = mode == "time"          # 'time' and an idealised 'lai' coincide: keep both visible
                ax.plot(s["date"], s["ks_mm_hr"], "-" if wide else "-o", ms=2.5,
                        lw=4.5 if wide else 1.5, alpha=0.45 if wide else 1.0,
                        color=COLORS[mode], label=mode)
        for r in rain_mm_day:
            ax.axhline(r / 24.0, color="0.4", ls=":", lw=1)
            ax.text(d["date"].min(), r / 24.0, f" {r:g} mm/day", va="bottom", fontsize=7, color="0.3")
        ax.axvline(pd.Timestamp(fire_date), color="#CC0000", ls="--", lw=1)
        ax.set_yscale("log"); ax.set_title(name, fontsize=10); ax.set_ylabel("Ks (mm/hr)")
        ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7)
    for ax in axes.ravel()[len(presets):]:
        ax.axis("off")
    fig.autofmt_xdate(); fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
    return fig


if __name__ == "__main__":
    import argparse
    import matplotlib
    matplotlib.use("Agg")
    ap = argparse.ArgumentParser(description="Plot the post-fire Ks scenarios (all presets, all modes)")
    ap.add_argument("--fire-date", default="2020-07-01")
    ap.add_argument("--start", default="2019-01-01")
    ap.add_argument("--end", default="2023-12-31")
    ap.add_argument("--out", default="ks_scenarios.png")
    a = ap.parse_args()
    months = pd.date_range(a.start, a.end, freq="MS")
    plot_ks_scenarios(ks_trajectories(months, a.fire_date), a.fire_date, out=a.out)
    print("saved", a.out)
