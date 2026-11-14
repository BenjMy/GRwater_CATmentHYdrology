"""
recovery_metrics.py — post-fire recovery metrics (TO-DO §4)
===========================================================
Pure pandas/numpy functions, no CATHY needed. Inputs are monthly series
(``pd.Series`` indexed by month, e.g. mean LAI, ETa/ETp, soil moisture,
recharge) for the burned and the control catchment.

Definitions
-----------
``ratio``          burned / control, month by month (control == 0 -> NaN)
``recovered``      the ratio is within ``tol`` of 1 (|ratio - 1| <= tol) for at
                   least ``sustain`` consecutive months
``recovery_time``  months from the fire date to the *start* of the first such
                   run (NaN if it never recovers inside the series)
``deficit``        sum over months of (control - burned), i.e. the total
                   anomaly (water or LAI units * months) accumulated after the fire
"""
from __future__ import annotations

from typing import Optional

import pandas as pd


def _check_aligned(burned: pd.Series, control: pd.Series):
    if not isinstance(burned, pd.Series) or not isinstance(control, pd.Series):
        raise TypeError("burned and control must be pandas Series")
    idx = burned.index.intersection(control.index)
    if len(idx) == 0:
        raise ValueError("burned and control share no time index")
    return burned.loc[idx].astype(float), control.loc[idx].astype(float)


def relative_to_control(burned: pd.Series, control: pd.Series) -> pd.Series:
    """Monthly ratio burned / control on the common index (NaN where control == 0)."""
    b, c = _check_aligned(burned, control)
    return b / c.where(c != 0)


def recovery_time(ratio: pd.Series, fire_date, tol: float = 0.1, sustain: int = 3) -> float:
    """
    Months from ``fire_date`` to the start of the first run of ``sustain``
    consecutive months with ``|ratio - 1| <= tol``. NaN if none.
    Only months on/after ``fire_date`` are considered.
    """
    if sustain < 1 or tol < 0:
        raise ValueError("sustain must be >= 1 and tol >= 0")
    fire = pd.Timestamp(fire_date)
    r = ratio[ratio.index >= fire]
    ok = ((r - 1.0).abs() <= tol).to_numpy()          # NaN compares False
    run = 0
    for i, flag in enumerate(ok):
        run = run + 1 if flag else 0
        if run >= sustain:
            start = r.index[i - sustain + 1]
            return float((start.year - fire.year) * 12 + (start.month - fire.month))
    return float("nan")


def deficit(burned: pd.Series, control: pd.Series, fire_date) -> float:
    """Sum of (control - burned) over the months from ``fire_date`` (NaN months skipped)."""
    b, c = _check_aligned(burned, control)
    d = (c - b)[b.index >= pd.Timestamp(fire_date)]
    return float(d.sum(skipna=True))


def recovery_table(burned: pd.DataFrame, control: pd.DataFrame, fire_date,
                   tol: float = 0.1, sustain: int = 3,
                   variables: Optional[list] = None) -> pd.DataFrame:
    """
    One row per variable (columns shared by ``burned`` and ``control``):
    ``recovery_months``, ``min_ratio`` (worst month after the fire),
    ``final_ratio`` (last month) and ``deficit``.
    """
    cols = variables or [c for c in burned.columns if c in control.columns]
    if not cols:
        raise ValueError("burned and control share no variable")
    rows = {}
    for v in cols:
        r = relative_to_control(burned[v], control[v])
        post = r[r.index >= pd.Timestamp(fire_date)]
        rows[v] = {
            "recovery_months": recovery_time(r, fire_date, tol, sustain),
            "min_ratio": float(post.min()) if post.notna().any() else float("nan"),
            "final_ratio": float(post.dropna().iloc[-1]) if post.notna().any() else float("nan"),
            "deficit": deficit(burned[v], control[v], fire_date),
        }
    return pd.DataFrame(rows).T
