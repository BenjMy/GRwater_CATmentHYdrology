"""
plot_common.py - Shared constants and helpers (fire date, colour maps, time/node helpers, windows).
=================================================================================================
Split unchanged out of results_plotter.py (REFACTORING.md 2.11); results_plotter re-exports every name.
"""
from __future__ import annotations

from typing import Any
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr


CMAP_VEG = mcolors.ListedColormap(["#d4b483", "#78c679", "#006837"])


COLORS_HYDRO: dict[str, str] = {
    "uphill_surface":   "#0072B2",
    "uphill_1m":        "#56B4E9",
    "downhill_surface": "#D55E00",
    "downhill_1m":      "#E69F00",
    "ET_mean":          "#009E73",
    "LAI":              "#4daf4a",
    "ETp":              "#999999",
}


FIRE_DATE = pd.Timestamp("2020-07-01")


# ── ET variable choice ─────────────────────────────────────────────────────
# The upstream pipeline writes two candidate ET variables into ET_xr_all /
# ET_xr_mean: the raw solver output ("ACT. ETRA") and the artefact-floored
# version ("ACT. ETRA_patched", see apply_eta_artefact_floor in
# Agramon_withLAI_withETp.py). Which one is "the" ET series is a per-run
# choice made upstream (see ET_VAR / WB_ET_VAR in the calling scripts) — this
# module must not hardcode either name; it only falls back to this default
# when a caller doesn't pass et_var explicitly.
DEFAULT_ET_VAR = "ACT. ETRA_patched"


def add_fire_vline(
    ax: plt.Axes,
    d0: pd.Timestamp | None = None,
    d1: pd.Timestamp | None = None,
    *,
    label: str = "Fire (Jul 2020)",
    color: str = "#B22222",
) -> None:
    """Draw a dashed vertical line marking the Agramon fire date (2020-07).

    Intended for any axis whose x-values are real datetimes. If *d0*/*d1*
    (the plotted date window) are supplied and ``FIRE_DATE`` falls outside
    that window, nothing is drawn — this keeps scenarios that don't cover
    July 2020 from having their x-axis distorted by an out-of-range line.
    """
    if d0 is not None and FIRE_DATE < pd.Timestamp(d0):
        return
    if d1 is not None and FIRE_DATE > pd.Timestamp(d1):
        return
    ax.axvline(FIRE_DATE, color=color, linestyle="--", linewidth=1.3,
               alpha=0.85, zorder=4, label=label)


def save_fig(fig: plt.Figure, fig_path: Path, name: str, dpi: int = 150) -> Path:
    """Save *fig* to *fig_path / name* and print the destination."""
    out = fig_path / name
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    return out


def maybe_show(fig: plt.Figure, show: bool) -> None:
    """Optionally display *fig*, then close it."""
    if show:
        plt.show()
    plt.close(fig)


def time_label(t_val: Any) -> str:
    """Human-readable string for a time coordinate (timedelta64 or datetime64)."""
    if np.issubdtype(np.array(t_val).dtype, np.timedelta64):
        t_days = float(t_val / np.timedelta64(1, "s")) / 86400
        return f"t = {t_days:.1f} d"
    return pd.Timestamp(t_val).strftime("%Y-%m-%d")


def et_global_limits(et_dataset: xr.Dataset, var: str) -> tuple[float, float]:
    """Return (vmin, vmax) across the whole ET dataset for consistent colour axes."""
    vals = et_dataset[var].values
    return float(np.nanmin(vals)), float(np.nanmax(vals))


def veg_colorbar(
    ax: plt.Axes,
    cmap: mcolors.Colormap,
    labels: tuple[str, ...] = ("Bare", "Moderate veg.", "Dense veg."),
) -> plt.colorbar:
    """Attach a discrete 3-class colourbar to *ax*."""
    bounds = [0.5, 1.5, 2.5, 3.5]
    norm   = mcolors.BoundaryNorm(bounds, cmap.N)
    sm     = cm.ScalarMappable(cmap=cmap, norm=norm)
    cb     = plt.colorbar(sm, ax=ax, ticks=[1, 2, 3])
    cb.ax.set_yticklabels(labels)
    return cb


def date_window(
    date_start: str | None,
    date_end:   str | None,
    df_dates:   pd.Series | pd.DatetimeIndex,
) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Return (d0, d1) as Timestamps, falling back to series bounds."""
    d0 = pd.Timestamp(date_start) if date_start else pd.Timestamp(df_dates.min())
    d1 = pd.Timestamp(date_end)   if date_end   else pd.Timestamp(df_dates.max())
    return d0, d1


LAI_UNITS = "m² m⁻²"


LAI_VMAX  = 2.0


CMAP_LAI  = "YlGn"


def _xr_time_to_dates(da: xr.DataArray, time_dim: str) -> pd.DatetimeIndex:
    """Convert any xarray time coordinate to a pandas DatetimeIndex."""
    raw = da[time_dim].values
    # timedelta64 relative to first value → absolute if offset stored as attr
    if np.issubdtype(raw.dtype, np.timedelta64):
        raw = raw.astype("timedelta64[s]")
        # Fall back: treat t=0 as origin (caller should supply pev_series.index[0])
        origin = pd.Timestamp("1970-01-01")
        return origin + pd.to_timedelta(raw)
    return pd.DatetimeIndex(raw)


def unwrap_node_id(val) -> int:
    """Extract a single scalar node index from an arbitrarily-nested value
    (int, (int,), ((int,),), np.ndarray, …) — e.g. whatever
    ``simu.find_nearest_node(...)`` happens to hand back."""
    while hasattr(val, "__iter__") and not isinstance(val, (str, np.integer)):
        val = next(iter(val))
    return int(val)


def _node_id(nodes_dict: dict, key: str) -> int:
    """Extract a scalar node index from the nodes_dict regardless of nesting."""
    return unwrap_node_id(nodes_dict[key])


def _node_ids(nodes_dict: dict, key: str) -> int | list[int]:
    """Like ``_node_id``, but preserves a mid-elevation *band* instead of
    always collapsing to a single scalar.

    ``resolve_nodes()`` stores a single reference node (e.g. the outlet) as
    a possibly-nested scalar (``5``, ``(5,)``, ``((5,),)``, …), and a
    mid-elevation band (see ``find_representative_mid_node(...,
    return_band=True)``) as a flat list/array of plain ints. We tell the
    two apart by checking whether the first element is itself iterable:
    a flat band's first element is a bare int, a nested single-id
    container's first element is another container.

    Returns
    -------
    int for a single reference node, or list[int] for a band of nodes.
    """
    val = nodes_dict[key]
    if isinstance(val, (list, tuple, np.ndarray)) and len(val) > 0 \
            and not hasattr(val[0], "__iter__"):
        return [int(v) for v in val]
    return unwrap_node_id(val)


def _fmt_node_id(ids: int | list[int]) -> str:
    """Human-readable label for a reference node or a mid-elevation band."""
    if isinstance(ids, (list, tuple, np.ndarray)):
        ids = list(ids)
        if len(ids) > 1:
            return f"median of {len(ids)} nodes"
        return str(int(ids[0]))
    return str(ids)


def _values_at_nodes(
    da: "xr.DataArray",
    time_dim: str,
    time_mask: np.ndarray,
    ids: int | list[int],
) -> np.ndarray:
    """Extract a (masked-time,) series from an xarray DataArray at one or
    more node ids.

    A scalar id returns that single node's series, unchanged. A band
    (list/array of several ids — e.g. all pixels within ±dz_pct of the
    target mid-elevation) returns the per-timestep *median* across those
    nodes, so "mid elevation" reflects the whole band rather than one
    arbitrarily-chosen point.
    """
    sub = da.isel({time_dim: time_mask})
    if isinstance(ids, (list, tuple, np.ndarray)):
        ids = np.atleast_1d(np.asarray(ids, dtype=int))
        if ids.size > 1:
            return sub.sel(node=ids).median(dim="node").values
        ids = int(ids[0])
    return sub.sel(node=int(ids)).values


def _fire_year_window() -> tuple[pd.Timestamp, pd.Timestamp]:
    """(d0, d1) spanning the fire year and the following year, e.g. 2020-2021."""
    d0 = pd.Timestamp(year=FIRE_DATE.year, month=1, day=1)
    d1 = pd.Timestamp(year=FIRE_DATE.year + 1, month=12, day=31)
    return d0, d1


def _fire_month_window() -> tuple[pd.Timestamp, pd.Timestamp]:
    """(d0, d1) spanning only the calendar month the fire occurred in."""
    d0 = FIRE_DATE.replace(day=1)
    d1 = d0 + pd.offsets.MonthEnd(0)
    return d0, d1


CMAP_ET    = "YlGnBu"


MONTH_ABBR = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def _infer_time_dim(ds: xr.Dataset | xr.DataArray) -> str:
    """Return the time dimension name, regardless of capitalisation."""
    dims = list(ds.dims)
    for candidate in ("time", "datetime", "t", "date", "TIME", "Date"):
        if candidate in dims:
            return candidate
    for d in dims:
        if d in ds.coords and np.issubdtype(ds[d].dtype, np.datetime64):
            return d
    raise ValueError(f"Cannot detect time dimension in dims={dims}")


def _spatial_dims(da: xr.DataArray, time_dim: str) -> list[str]:
    """Return all non-time dimension names."""
    return [d for d in da.dims if d != time_dim]


def _da(ds: xr.Dataset, et_var: str) -> xr.DataArray:
    if et_var not in ds:
        raise KeyError(f"Variable '{et_var}' not found. Available: {list(ds.data_vars)}")
    return ds[et_var]


MONTHLY_RULE = "MS"


MONTHLY_BAR_WIDTH = 25.0   # days — keeps rain/ETp bars visible on a monthly axis


def _resample_monthly(obj, how: str = "mean"):
    """Monthly resample of a Dataset / DataArray / Series (None-safe)."""
    if obj is None:
        return None
    if isinstance(obj, pd.Series):
        return getattr(obj.resample(MONTHLY_RULE), how)() if len(obj) else obj
    if len(obj.dims) == 0:
        return obj
    td = _infer_time_dim(obj)
    return getattr(obj.resample({td: MONTHLY_RULE}), how)()
