"""
test_withLAI_results.py — Recharge & LAI Loop Results Viewer
=============================================================
GRwater project — Agramon study site

Loads pre-computed NetCDF / pickle artefacts saved by
``LT_SShydro_run_withLAI.py`` and reproduces all figures without
re-running the CATHY simulation.

The script mirrors the scenario-selection pattern of LT_SShydro_results.py:
  • It reads simulation_log.csv to locate the correct sim folder.
  • Artefacts are loaded from  <path2prj>/sim<N>/outputs/
  • Figures are written to     <path2prj>/figures/scenario_<N>/
    (or to --fig-dir if supplied).
  • --list-scenarios prints withLAI=1 rows from the log and exits.

Expected artefacts under <path2prj>/sim<N>/outputs/:
    et_output.nc          — xr.Dataset   (X, Y, datetime)
    psi_output.nc         — xr.Dataset   (datetime, node)
    sw_output.nc          — xr.Dataset   (datetime, node)
    veg_map_history.nc    — xr.DataArray (time, row, col)
    # recharge_monthly.csv  — pd.Series    monthly recharge [mm]  (optional)
    veg_map_history.pkl   — dict {datetime: lai_2d}              (optional)

Figures produced
----------------
  Fig 1  — Monthly recharge bar + cumulative line
  Fig 2  — ET spatial mean by month / year
  Fig 3  — ET time series per vegetation type
  Fig 4  — Spatial ET maps per month for SELECTED_YEAR (3×4 subplots)
  Fig 4b — Spatial surface soil-moisture saturation maps per month for
           SELECTED_YEAR (same 3×4 layout, surface-node values from
           sw_output.nc interpolated onto the ET raster grid)
  Fig 5  — Vegetation zone map samples (one per year)
  Fig 5b — Monthly vegetation maps for SELECTED_YEAR
  Fig 5c — Monthly LAI maps for SELECTED_YEAR
  LAI-A  — Spatial-mean vegetation index (unitless) per month, one line per year (→ LAI_monthly_mean_by_year_*.png)
  LAI-B  — Monthly vegetation index spatial maps for SELECTED_YEAR (4-col mosaic)(→ LAI_spatial_maps_<year>_*.png)
  LAI-C  — Histogram of per-pixel annual mean vegetation index (unitless)         (→ LAI_annual_mean_histogram_*.png)
  LAI-D  — Bar chart of domain-mean annual vegetation index (unitless)            (→ LAI_annual_mean_*.png)
  LAI-E  — Yearly time series of domain-median vegetation index (unitless)        (→ LAI_annual_median_*.png)
  Fig 6  — Hydro time-series: Rain/ETp + ψ evolution at reference nodes
  Fig 6b — 2D map of reference-node locations (where ψ / sw are extracted)
  Fig 7  — sw + Actual ETa time-series at reference nodes
  Fig 7b — sw + ETa + domain-mean LAI trend time-series (Fig 7 + LAI panel)
  Fig 8  — Catchment ET balance: monthly + annual bars  (→ ETa_balance_*.png)
  Fig 8b — ETa cumulative histogram: per-pixel annual sum distribution  (→ ETa_annual_sum_histogram_*.png)
  Fig 8c — ETa domain-total: monthly mean + annual total bar charts    (→ ETa_total_*.png)
  Fig 9  — Hydraulic graph                              (→ hgraph.png)
  Fig 10 — Recharge: monthly + annual bar charts        (→ recharge*.png)
  GIF A  — Animated vegetation zone maps
  GIF B  — Animated ET spatial maps

Usage
-----
    python test_withLAI_results.py                              # first withLAI=1 row
    python test_withLAI_results.py --scenario 3
    python test_withLAI_results.py --scenario 3 --selected-year 2020
    python test_withLAI_results.py --scenario 3 --date-start 2019-01-01 --date-end 2020-12-31
    python test_withLAI_results.py --list-scenarios
    python test_withLAI_results.py --no-plots --no-gifs --dpi 200
    python test_withLAI_results.py --outlet-xy 456200 4122500 --uphill-xy 456800 4123200
"""

import argparse
import pickle
import re
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")   # safe headless default; --backend can override

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.colors as mcolors
import rioxarray as rxr  # noqa: F401  (activates .rio accessor)
from Agramon_withLAI_withETp import load_dem_from_tif  # or move it to Agramon_utils

# ── Path setup ────────────────────────────────────────────────────────────────

MODULE_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
    #"/home/ben/Nextcloud/BenCSIC/Codes"
    "/Tech4agro_org/GRwater_geophy"
).resolve()

EO_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
    #"/home/ben/Nextcloud/BenCSIC/Codes"
    "/Tech4agro_org/GRwater_CATmentHYdrology"
).resolve()

LAI_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Training_Supervision"
    #"/home/ben/Nextcloud/BenCSIC/Training_Supervision"
    "/Supervision/Xela_Carracedo_Practicas_2026_data/301a-biophysical/agramon/input"
).resolve()

# Satellite ETp + rain forcing (replaces ERA5 for the Fig 6/7 rain/ETp
# panel — see Agramon_withLAI_withETp.py, which uses the same files to
# force the ATMBC boundary).
PathET    = Path("/home/z0272571a@campus.csic.es/Nextcloud/GRwater/data/satellite/Agramon")
#PathET    = Path("/home/ben/Nextcloud/GRwater/data/satellite/Agramon")
etp_path  = PathET / "20161001_20241231_ET0.nc"
rain_path = PathET / "20161001_20241231_TP.nc"

# LAI / ETp / rain source rasters are EPSG:32630 (UTM 30N, WGS84); DEM /
# shapefiles / mesh are TARGET_CRS = EPSG:25830 (UTM 30N, ETRS89) defined
# below. Same zone, different datum — always reproject before use.
NATIVE_CRS = "EPSG:32630"


if str(MODULE_PATH) not in sys.path:
    sys.path.append(str(MODULE_PATH))

import Agramon_utils as AgUtils
from pyCATHY import cathy_tools

# ── Shared plotter module ─────────────────────────────────────────────────────

import results_plotter as plotter
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=DeprecationWarning)

# ── Plot constants ────────────────────────────────────────────────────────────

TARGET_CRS = "EPSG:25830"
ET_VAR     = 'ACT. ETRA_patched' #"ACT. ETRA" ACT. ETRA_patched
#ET_VAR     = 'ACT. ETRA' #"ACT. ETRA" ACT. ETRA_patched
ET_SCALE   = 1e3 * 86400                           # m/s → mm/day
CMAP_ET    = "YlGnBu"
CMAP_SW    = "Blues"
CMAP_VEG   = mcolors.ListedColormap(["#d4b483", "#78c679", "#006837"])

LAI_VAR  = "LAI"
LAI_CMAP = "YlGn"
LAI_VMIN = 0.0
LAI_VMAX = 3

MONTH_NAMES = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
]


# =============================================================================
# Satellite ETp + rain forcing (Fig 6/7 rain panel — replaces ERA5)
# =============================================================================
def _load_forcing_nc(path: Path, label: str) -> xr.Dataset:
    """Open a satellite ETp/rain NetCDF, assign its native CRS (EPSG:32630),
    and reproject to TARGET_CRS (EPSG:25830). Mirrors ``_load_nc`` in
    Agramon_withLAI_withETp.py.
    """
    ds = (
        xr.open_dataset(path, decode_times=True, mask_and_scale=True)
          .drop_vars("spatial_ref", errors="ignore")
    )
    ds = ds.rio.write_crs(NATIVE_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  {label} native CRS : {ds.rio.crs}")

    if str(ds.rio.crs) != TARGET_CRS:
        print(f"  Reprojecting {label}  {ds.rio.crs} → {TARGET_CRS}")
        ds = ds.rio.reproject(TARGET_CRS)

    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    return ds


def load_etp(path: Path = etp_path) -> xr.Dataset:
    """Load satellite ETp (ET0) dataset and reproject to TARGET_CRS."""
    return _load_forcing_nc(path, "ETp")


def load_rain(path: Path = rain_path) -> xr.Dataset:
    """Load satellite rain (TP) dataset and reproject to TARGET_CRS."""
    return _load_forcing_nc(path, "rain")


def _domain_mean_daily_series(ds: xr.Dataset, var: str, gdf) -> pd.Series:
    """Clip a satellite ETp/rain dataset to the catchment polygon and
    return its daily domain-mean as a ``pd.Series`` (mm/day, DatetimeIndex).

    Uses rioxarray directly instead of ``AgUtils.extract_point_timeseries``:
    that helper was written for ERA5's lat/lon grid, whereas ``ds`` here is
    already projected to TARGET_CRS (EPSG:25830, x/y in metres) — calling
    it on this grid is what produced "satellite ETp/rain forcing not
    loaded" (it can't resolve ERA5-style coordinates/variable names on a
    plain x/y raster). A direct clip + spatial mean sidesteps that and
    also matches how the main run pipeline (Agramon_withLAI_withETp.py)
    already treats ETp/rain — as a mesh/domain average, not a point pick.
    """
    da = ds[var]

    gdf_proj = gdf
    if gdf.crs is not None and str(gdf.crs) != TARGET_CRS:
        gdf_proj = gdf.to_crs(TARGET_CRS)

    try:
        da = da.rio.clip(gdf_proj.geometry.values, gdf_proj.crs, drop=True)
    except Exception as clip_exc:
        print(f"  WARNING: clip to catchment failed for '{var}' — "
              f"using full raster extent instead.\n    ({clip_exc})")

    time_dim = "time" if "time" in da.dims else "datetime"
    spatial_dims = [d for d in da.dims if d != time_dim]
    series = da.mean(dim=spatial_dims, skipna=True).to_series()
    series.index = pd.DatetimeIndex(series.index)
    return series


# =============================================================================
# Spatial ET monthly subplot (Fig 4)
# =============================================================================
def plot_et_spatial_monthly(
    ET_xr_all: xr.Dataset,
    variable: str,
    scale: float,
    selected_year: int,
    out_dir: Path | None = None,
    cmap: str = "YlGnBu",
    dpi: int = 150,
    show: bool = False,
) -> None:
    """
    Plot monthly spatial ET maps for a selected year using xarray's imshow.

    Parameters
    ----------
    ET_xr_all : xr.Dataset
        Full ET dataset with dims (X, Y, datetime), as saved to et_output.nc.
    variable : str
        Variable to plot (e.g. 'ACT. ETRA').
    scale : float
        Multiplicative scale factor applied to the data (e.g. 1e3 * 86400 for m/s → mm/day).
    selected_year : int
        Calendar year to subset and display.
    out_dir : Path or None
        If given, saves the figure as a PNG. Otherwise shows interactively.
    cmap : str
        Matplotlib colormap name.
    dpi : int
        Resolution for the saved figure.
    show : bool
        If True, calls plt.show() after saving/building the figure.
    """
    da_all = ET_xr_all[variable] * scale
    _plot_spatial_monthly_generic(
        da_all, "datetime", selected_year, out_dir, cmap, dpi, show,
        title=f"Monthly ET ({variable})",
        cbar_label=f"{variable} (mm/day)",
        fname_prefix="et_spatial_monthly",
    )


def _plot_spatial_monthly_generic(
    da_all:        xr.DataArray,
    time_dim:      str,
    selected_year: int,
    out_dir:       Path | None,
    cmap:          str,
    dpi:           int,
    show:          bool,
    title:         str,
    cbar_label:    str,
    fname_prefix:  str,
    vmin:          float | None = None,
    vmax:          float | None = None,
) -> None:
    """
    Shared 3×4 (month) mosaic imshow, used by both plot_et_spatial_monthly
    and plot_sw_spatial_monthly so the two stay visually consistent.
    """
    da_year = da_all.sel({time_dim: da_all[time_dim].dt.year == selected_year})

    if da_year.sizes[time_dim] == 0:
        print(f"  {fname_prefix}: no data for {selected_year} — skipping.")
        return

    months   = np.unique(da_year[time_dim].dt.month.values)
    n_months = len(months)
    ncols    = min(4, n_months)
    nrows    = int(np.ceil(n_months / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5 * ncols, 4 * nrows),
        squeeze=False,
    )

    vmin_eff = float(da_year.min()) if vmin is None else vmin
    vmax_eff = float(da_year.max()) if vmax is None else vmax

    for idx, month in enumerate(months):
        ax = axes[idx // ncols][idx % ncols]

        da_month = da_year.sel(
            {time_dim: da_year[time_dim].dt.month == month}
        ).mean(time_dim)                          # collapse time steps → (X, Y)

        da_month.plot.imshow(
            ax=ax,
            cmap=cmap,
            x="X",
            y="Y",
            vmin=vmin_eff,
            vmax=vmax_eff,
            add_colorbar=(idx == 0),
            **( {"cbar_kwargs": {"label": cbar_label, "shrink": 0.8}}
                if idx == 0 else {} ),
        )

        ax.set_title(MONTH_NAMES[month - 1], fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_aspect("equal")

    # Hide unused axes
    for idx in range(n_months, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(f"{title} — {selected_year}", fontsize=12, y=1.01)
    plt.tight_layout()

    if out_dir is not None:
        out_path = Path(out_dir) / f"{fname_prefix}_{selected_year}.png"
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        print(f"  {title} map saved → {out_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


# =============================================================================
# Surface soil-moisture saturation → spatial grid (Fig 4b)
# =============================================================================
# sw_output.nc only carries values at 3-D mesh *nodes* (dims: datetime,
# node) — there's no ready-made (X, Y) raster the way there is for ET.
# build_sw_surface_grid() recovers one by interpolating the surface-layer
# nodes (grid3d['mesh3d_nodes'][:nnod]) onto the exact (X, Y) coordinates
# ET_xr_all already uses, so the result can be imshow'n with the same
# machinery as ET/LAI. nearest-neighbour is used deliberately (not
# linear/cubic): these are point values at the same physical locations
# the ET raster's own grid is built on, so no smoothing should be
# introduced by the regridding step itself.

def build_sw_surface_grid(
    sw_xr:     xr.Dataset,
    grid3d:    dict,
    ET_xr_all: xr.Dataset,
    sw_var:    str = "sw",
    x_dim:     str = "X",
    y_dim:     str = "Y",
    mask_to_dem: bool = True,
) -> xr.Dataset:
    """
    Interpolate surface-node soil-moisture saturation onto the same
    (X, Y) raster grid used by et_output.nc.

    Parameters
    ----------
    sw_xr : xr.Dataset
        Raw sw_output.nc — dims (datetime, node).
    grid3d : dict
        Return value of ``simu.read_outputs("grid3d")``; needs ``"nnod"``
        (surface node count) and ``"mesh3d_nodes"`` (all mesh nodes,
        columns x, y, z — the first ``nnod`` rows are the surface layer).
    ET_xr_all : xr.Dataset
        et_output.nc — supplies the target (X, Y) raster coordinates.
    sw_var : str
        Variable name inside sw_xr (default: "sw").
    mask_to_dem : bool
        If True (default), set pixels outside the DEM footprint to NaN.
        The footprint is taken directly from ``ET_xr_all`` itself (it is
        already NaN outside the catchment, since it's real CATHY model
        output on the DEM's own active cells) rather than from a
        separately-loaded DEM raster — that avoids any risk of shape or
        orientation mismatch between a standalone raster_DEM array and
        this (Y, X) grid, which is what silently skipped the masking
        before.

    Returns
    -------
    xr.Dataset with dims (datetime, Y, X) and variable *sw_var*,
    coordinates matching ET_xr_all exactly — ready for
    plot_sw_spatial_monthly() or any other ET-style spatial figure.
    """
    from scipy.interpolate import griddata

    if sw_var not in sw_xr:
        raise KeyError(
            f"'{sw_var}' not found in sw_output.nc. "
            f"Available: {list(sw_xr.data_vars)}"
        )

    nnod          = int(grid3d["nnod"])
    surface_nodes = grid3d["mesh3d_nodes"][:nnod]
    node_xy       = surface_nodes[:, :2]

    da_sw    = sw_xr[sw_var]
    time_dim = "datetime" if "datetime" in da_sw.dims else "time"
    sw_surf  = da_sw.isel(node=slice(0, nnod))

    X_vals = ET_xr_all[x_dim].values
    Y_vals = ET_xr_all[y_dim].values
    Xg, Yg = np.meshgrid(X_vals, Y_vals)

    n_t = sw_surf.sizes[time_dim]
    grid_vals = np.empty((n_t, len(Y_vals), len(X_vals)), dtype=float)
    for i in range(n_t):
        vals = sw_surf.isel({time_dim: i}).values
        grid_vals[i] = griddata(node_xy, vals, (Xg, Yg), method="nearest")

    # Crop to the DEM footprint — griddata's "nearest" fill extrapolates
    # a value into every pixel of the (Y, X) bounding box, including
    # outside the actual catchment. Reuse ET_xr_all's own NaN footprint
    # (same grid, guaranteed aligned) instead of a separate DEM raster.
    if mask_to_dem:
        et_var_any  = next(iter(ET_xr_all.data_vars))
        et_da       = ET_xr_all[et_var_any]
        et_time_dim = next((d for d in ("datetime", "time") if d in et_da.dims), None)
        if et_time_dim is not None:
            et_da = et_da.isel({et_time_dim: 0})
        et_da = et_da.transpose(y_dim, x_dim)
        in_dem = np.isfinite(et_da.values)

        if in_dem.shape == (len(Y_vals), len(X_vals)):
            grid_vals = np.where(in_dem[None, :, :], grid_vals, np.nan)
            print(f"  SW surface grid: masked {int((~in_dem).sum())} "
                  f"pixel(s) outside the DEM (via ET footprint).")
        else:
            print(f"  WARNING: ET footprint shape {in_dem.shape} != "
                  f"({len(Y_vals)}, {len(X_vals)}) — DEM masking skipped.")

    ds_sw_grid = xr.Dataset(
        {sw_var: ((time_dim, y_dim, x_dim), grid_vals)},
        coords={
            time_dim: sw_surf[time_dim].values,
            y_dim: Y_vals,
            x_dim: X_vals,
        },
    )
    print(f"  SW surface grid built: '{sw_var}' on ({y_dim}, {x_dim}) = "
          f"{ds_sw_grid.sizes[y_dim]}x{ds_sw_grid.sizes[x_dim]}, "
          f"{n_t} timesteps (from {nnod} surface nodes).")
    return ds_sw_grid


def plot_sw_spatial_monthly(
    ds_sw_grid:    xr.Dataset,
    variable:      str = "sw",
    selected_year: int = 2020,
    out_dir:       Path | None = None,
    cmap:          str = "Blues",
    dpi:           int = 150,
    show:          bool = False,
    vmin:          float | None = 0.0,
    vmax:          float | None = 1.0,
) -> None:
    """
    Monthly spatial maps of surface soil-moisture saturation for
    *selected_year*, in the same 3×4 mosaic layout as
    plot_et_spatial_monthly.

    Parameters
    ----------
    ds_sw_grid : xr.Dataset
        Output of build_sw_surface_grid() — dims (datetime|time, Y, X).
    variable : str
        Variable name inside ds_sw_grid (default: "sw").
    vmin, vmax : float or None
        Fixed colour-scale bounds (default 0-1, the usual saturation
        range). Pass None/None to auto-scale per year, like ET.
    """
    time_dim = "datetime" if "datetime" in ds_sw_grid[variable].dims else "time"
    _plot_spatial_monthly_generic(
        ds_sw_grid[variable], time_dim, selected_year, out_dir, cmap, dpi, show,
        title=f"Monthly surface saturation ({variable})",
        cbar_label=f"{variable} (-)",
        fname_prefix="sw_spatial_monthly",
        vmin=vmin, vmax=vmax,
    )


# =============================================================================
# Shared pre-processing helper
# =============================================================================

def _prepare_et_ds(
    ET_xr_all: xr.Dataset,
    variable: str,
    scale: float,
) -> xr.Dataset:
    """
    Return a scaled, plotter-compatible copy of *ET_xr_all*.

    ``results_plotter`` functions expect:
      • a ``"time"`` dimension (not ``"datetime"`` as in et_output.nc)
      • values already in mm/day (no further scaling needed)

    This helper applies *scale* in-place on *variable* and renames
    ``"datetime"`` → ``"time"`` if necessary, leaving all other variables
    and coordinates untouched.

    Parameters
    ----------
    ET_xr_all : xr.Dataset  — raw et_output.nc dataset (X, Y, datetime)
    variable  : str         — ET variable name (e.g. 'ACT. ETRA')
    scale     : float       — factor converting raw units to mm/day
    """
    if variable not in ET_xr_all:
        raise KeyError(
            f"Variable '{variable}' not found. "
            f"Available: {list(ET_xr_all.data_vars)}"
        )
    # Apply unit conversion on the target variable only
    ds = ET_xr_all.assign({variable: ET_xr_all[variable] * scale})

    # Rename 'datetime' → 'time' so results_plotter internals work unchanged.
    # et_output.nc often carries a scalar/auxiliary coord also called 'time',
    # which makes xarray raise "the new name 'time' conflicts" even though
    # 'time' is not a dimension.  Drop it first in that case.
    # if "datetime" in ds.dims and "time" not in ds.dims:
    #     if "time" in ds.coords:
    #         ds = ds.drop_vars("time")
    #     ds = ds.rename({"datetime": "time"})

    return ds


# =============================================================================
# ETa cumulated histogram (Fig 8b)  — delegates to results_plotter
# =============================================================================

def plot_eta_annual_hist(
    ET_xr_all:   xr.Dataset,
    variable:    str         = ET_VAR,
    scale:       float       = ET_SCALE,
    out_dir:     Path | None = None,
    fname:       str         = "ETa_annual_sum_histogram.png",
    scenario_id: int | str   = "",
    bins:        int         = 40,
    dpi:         int         = 150,
    show:        bool        = False,
) -> None:
    """
    Histogram of per-pixel annual ETa sums — thin wrapper around
    ``results_plotter.plot_eta_annual_hist``.

    *ET_xr_all* (dims X, Y, datetime; units m/s) is scaled to mm/day and its
    time dimension renamed to ``"time"`` before being forwarded, so the shared
    plotter receives exactly the interface it expects.

    Parameters
    ----------
    ET_xr_all   : xr.Dataset  — et_output.nc dataset (X, Y, datetime)
    variable    : str          — ET variable name (default: 'ACT. ETRA')
    scale       : float        — m/s → mm/day factor (default: 1e3 × 86400)
    out_dir     : Path | None  — output directory (created if absent)
    fname       : str          — output filename
    scenario_id : int | str    — appended to fname stem when non-empty
    bins        : int          — histogram bin count (default: 40)
    dpi         : int          — figure resolution (default: 150)
    show        : bool         — call plt.show() before closing
    """
    if out_dir is None:
        print("  plot_eta_annual_hist: out_dir is None — skipping.")
        return

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Embed scenario_id in the filename stem when provided
    if scenario_id != "":
        stem, suffix = fname.rsplit(".", 1) if "." in fname else (fname, "png")
        fname = f"{stem}_{scenario_id}.{suffix}"

    ds_ready = _prepare_et_ds(ET_xr_all, variable, scale)

    plotter.plot_eta_annual_hist(
        ds_ready,
        et_var   = variable,
        fig_path = out_dir,
        et_scale = 1.0,       # already scaled inside _prepare_et_ds
        bins     = bins,
        dpi      = dpi,
        show     = show,
        fname    = fname,
    )


# =============================================================================
# ETa domain-total bar charts (Fig 8c)  — delegates to results_plotter
# =============================================================================

def plot_eta_total(
    ET_xr_all:   xr.Dataset,
    variable:    str       = ET_VAR,
    scale:       float     = ET_SCALE,
    out_dir:     Path | None = None,
    label:       str       = "",
    scenario_id: int | str = "",
    dpi:         int       = 150,
    show:        bool      = False,
) -> None:
    """
    Monthly and annual domain-total ETa bar charts — thin wrapper around
    ``results_plotter.plot_eta_total``.

    Produces two PNGs in *out_dir*:
        ETa_total_monthly[_<label>].png
        ETa_total_annual[_<label>].png

    *ET_xr_all* (dims X, Y, datetime; units m/s) is scaled to mm/day and its
    time dimension renamed to ``"time"`` before being forwarded, so the shared
    plotter receives exactly the interface it expects.

    Parameters
    ----------
    ET_xr_all   : xr.Dataset  — et_output.nc dataset (X, Y, datetime)
    variable    : str          — ET variable name (default: 'ACT. ETRA')
    scale       : float        — m/s → mm/day factor (default: 1e3 × 86400)
    out_dir     : Path | None  — output directory (created if absent)
    label       : str          — scenario label appended to filenames
    scenario_id : int | str    — appended to label when non-empty
    dpi         : int          — figure resolution (default: 150)
    show        : bool         — call plt.show() before closing
    """
    if out_dir is None:
        print("  plot_eta_total: out_dir is None — skipping.")
        return

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build a combined label that includes the scenario id when available
    full_label = label
    if scenario_id != "":
        full_label = f"{label}_{scenario_id}" if label else str(scenario_id)

    # ds_ready = _prepare_et_ds(ET_xr_all, variable, scale)

    ds_ready = ET_xr_all.drop_vars("time")
    ds_ready = ds_ready.rename({"datetime": "time"})

    plotter.plot_eta_total(
        ds_ready,
        fig_path = out_dir,
        label    = full_label,
        dpi      = dpi,
        show     = show,
    )


# =============================================================================
# Vegetation index diagnostic wrappers  (Figs LAI-A … LAI-D)
# Accept the in-memory index_history dict {datetime: 2-D np.ndarray} that is
# already loaded from veg_map_history.pkl — no NetCDF file access needed.
# The dict is converted to an xr.Dataset via plotter.lai_history_to_dataset
# before being forwarded to the four plotter functions.
# =============================================================================

def _lai_ds_from_history(lai_history: dict, lai_var: str, caller: str) -> xr.Dataset | None:
    """Convert lai_history dict → xr.Dataset; return None and warn if empty."""
    if not lai_history:
        print(f"  {caller}: lai_history is empty — skipping.")
        return None
    return plotter.lai_history_to_dataset(lai_history, lai_var=lai_var)


def plot_index_monthly_by_year(
    index_history: dict,
    index_var:     str         = "LAI",
    out_dir:       Path | None = None,
    fname:         str | None  = None,
    scenario_id:   int | str   = "",
    dpi:           int         = 150,
    show:          bool        = False,
) -> None:
    """
    Spatial-mean vegetation index (unitless) per calendar month, one line per year.

    Thin wrapper around ``results_plotter.plot_lai_monthly_by_year``.

    Parameters
    ----------
    index_history : dict mapping ``datetime`` → 2-D index array (row, col),
                    as loaded from ``veg_map_history.pkl``
    index_var     : variable name used in the output dataset (default: 'LAI')
    out_dir       : output directory (created if absent)
    fname         : output filename base; auto-generated if None
    scenario_id   : appended to fname stem when non-empty
    dpi           : figure resolution
    show          : call plt.show() before closing
    """
    if out_dir is None:
        print("  plot_index_monthly_by_year: out_dir is None — skipping.")
        return

    ds = _lai_ds_from_history(index_history, index_var, "plot_index_monthly_by_year")
    if ds is None:
        return

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if fname is None:
        fname = f"veg_map_monthly_mean_by_year.png"

    if scenario_id != "":
        stem, suffix = fname.rsplit(".", 1) if "." in fname else (fname, "png")
        fname = f"{stem}_{scenario_id}.{suffix}"

    plotter.plot_lai_monthly_by_year(
        ds,
        lai_var  = index_var,
        fig_path = out_dir,
        dpi      = dpi,
        show     = show,
        fname    = fname,
    )


def plot_index_spatial_maps_year(
    index_history: dict,
    index_var:     str         = "LAI",
    selected_year: int         = 2021,
    out_dir:       Path | None = None,
    fname:         str | None  = None,
    scenario_id:   int | str   = "",
    dpi:           int         = 150,
    show:          bool        = False,
) -> None:
    """
    Monthly spatial maps (4-column mosaic) of a gridded index for *selected_year*.

    Thin wrapper around ``results_plotter.plot_index_spatial_maps_year``.

    Parameters
    ----------
    index_history : dict mapping ``datetime`` → 2-D array (row, col)
    index_var     : variable name used in the output dataset (default: 'LAI')
    selected_year : calendar year to display
    out_dir       : output directory (created if absent)
    fname         : explicit output filename; auto-generated if None
    scenario_id   : appended to fname stem when non-empty
    dpi           : figure resolution
    show          : call plt.show() before closing
    """
    if out_dir is None:
        print("  plot_index_spatial_maps_year: out_dir is None — skipping.")
        return

    ds = _lai_ds_from_history(index_history, index_var, "plot_index_spatial_maps_year")
    if ds is None:
        return

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if fname is None:
        suffix_id = f"_{scenario_id}" if scenario_id != "" else ""
        fname = f"{index_var}_spatial_maps_{selected_year}{suffix_id}.png"
    elif scenario_id != "":
        stem, ext = fname.rsplit(".", 1) if "." in fname else (fname, "png")
        fname = f"{stem}_{scenario_id}.{ext}"

    plotter.plot_index_spatial_maps_year(
        ds,
        index_var     = index_var,
        selected_year = selected_year,
        fig_path      = out_dir,
        dpi           = dpi,
        show          = show,
        fname         = fname,
    )


def plot_index_annual_hist(
    index_history: dict,
    index_var:     str         = "LAI",
    out_dir:       Path | None = None,
    fname:         str | None  = None,
    scenario_id:   int | str   = "",
    bins:          int         = 40,
    dpi:           int         = 150,
    show:          bool        = False,
) -> None:
    """
    Histogram of per-pixel annual mean values for a gridded index.

    Thin wrapper around ``results_plotter.plot_index_annual_hist``.

    Parameters
    ----------
    index_history : dict mapping ``datetime`` → 2-D array (row, col)
    index_var     : variable name used in the output dataset (default: 'LAI')
    out_dir       : output directory (created if absent)
    fname         : output filename base; auto-generated if None
    scenario_id   : appended to fname stem when non-empty
    bins          : histogram bin count
    dpi           : figure resolution
    show          : call plt.show() before closing
    """
    if out_dir is None:
        print("  plot_index_annual_hist: out_dir is None — skipping.")
        return

    ds = _lai_ds_from_history(index_history, index_var, "plot_index_annual_hist")
    if ds is None:
        return

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if fname is None:
        fname = f"{index_var}_annual_mean_histogram.png"

    if scenario_id != "":
        stem, suffix = fname.rsplit(".", 1) if "." in fname else (fname, "png")
        fname = f"{stem}_{scenario_id}.{suffix}"

    plotter.plot_index_annual_hist(
        ds,
        index_var = index_var,
        fig_path  = out_dir,
        bins      = bins,
        dpi       = dpi,
        show      = show,
        fname     = fname,
    )


def plot_index_annual_mean(
    index_history: dict,
    index_var:     str         = "LAI",
    out_dir:       Path | None = None,
    label:         str         = "",
    scenario_id:   int | str   = "",
    dpi:           int         = 150,
    show:          bool        = False,
) -> None:
    """
    Bar chart of domain-mean annual values for a gridded index.

    Thin wrapper around ``results_plotter.plot_index_annual_mean``.
    Produces one PNG:
        <index_var>_annual_mean[_<label>].png

    Parameters
    ----------
    index_history : dict mapping ``datetime`` → 2-D array (row, col)
    index_var     : variable name used in the output dataset (default: 'LAI')
    out_dir       : output directory (created if absent)
    label         : scenario label appended to the filename
    scenario_id   : appended to label when non-empty
    dpi           : figure resolution
    show          : call plt.show() before closing
    """
    if out_dir is None:
        print("  plot_index_annual_mean: out_dir is None — skipping.")
        return

    ds = _lai_ds_from_history(index_history, index_var, "plot_index_annual_mean")
    if ds is None:
        return

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    full_label = label
    if scenario_id != "":
        full_label = f"{label}_{scenario_id}" if label else str(scenario_id)

    plotter.plot_index_annual_mean(
        ds,
        index_var = index_var,
        fig_path  = out_dir,
        label     = full_label,
        dpi       = dpi,
        show      = show,
    )


def plot_index_annual_median(
    index_history: dict,
    index_var:     str         = "LAI",
    out_dir:       Path | None = None,
    label:         str         = "",
    scenario_id:   int | str   = "",
    dpi:           int         = 150,
    show:          bool        = False,
) -> None:
    """
    Yearly time series of the domain-median value for a gridded index.

    Thin wrapper around ``results_plotter.plot_lai_annual_median``.
    Produces one PNG:
        <index_var>_annual_median[_<label>].png

    Parameters
    ----------
    index_history : dict mapping ``datetime`` → 2-D array (row, col)
    index_var     : variable name used in the output dataset (default: 'LAI')
    out_dir       : output directory (created if absent)
    label         : scenario label appended to the filename
    scenario_id   : appended to label when non-empty
    dpi           : figure resolution
    show          : call plt.show() before closing
    """
    if out_dir is None:
        print("  plot_index_annual_median: out_dir is None — skipping.")
        return

    ds = _lai_ds_from_history(index_history, index_var, "plot_index_annual_median")
    if ds is None:
        return

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    full_label = label
    if scenario_id != "":
        full_label = f"{label}_{scenario_id}" if label else str(scenario_id)

    plotter.plot_lai_annual_median(
        ds,
        lai_var  = index_var,
        fig_path = out_dir,
        label    = full_label,
        dpi      = dpi,
        show     = show,
    )

# =============================================================================
# CLI
# =============================================================================

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="LAI loop post-processing — figure generation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Scenario selection — mirrors LT_SShydro_results.py
    p.add_argument(
        "--scenario", type=int, default=5, metavar="N",
        help=(
            "sim_index to visualise (from simulation_log.csv). "
            "Defaults to the first withLAI=1 row in the log."
        ),
    )
    p.add_argument(
        "--list-scenarios", action="store_true",
        help="Print available withLAI scenarios from the log file and exit.",
    )

    # Paths
    p.add_argument(
        "--path2prj", default="../SSHydro_withETp_withLAI", metavar="PATH",
        help=(
            "Parent directory containing sim<N>/ sub-folders. "
            "Artefacts are loaded from <path2prj>/sim<N>/outputs/."
        ),
    )
    p.add_argument(
        "--log-file", default="simulation_log_withETp_withLAI.csv", metavar="PATH",
        help=(
            "Shared CSV log produced by LT_SShydro_run_withLAI.py. "
            "Relative paths are resolved under <path2prj>."
        ),
    )
    p.add_argument(
        "--fig-dir", default=None, metavar="PATH",
        help=(
            "Directory where PNG figures are saved. "
            "Defaults to <path2prj>/figures/scenario_<N>/."
        ),
    )
    p.add_argument(
        "--dem-plot", type=int, default=3, choices=range(1, 9), metavar="N",
        help="DTM plot index (1-8) used to load the DEM raster.",
    )
    p.add_argument(
        "--dem-tif", default=str(MODULE_PATH / "DTMplots/20250618_AGRAMON100m_micasense_dtm.tif"),
        metavar="PATH",
        help=(
            "Full path to the GeoTIFF DEM used to build the mesh "
            "(must match --dem-tif used at generation time in "
            "Agramon_withLAI_withETp.py; default: "
            "20250618_AGRAMON100m_micasense_dtm.tif in MODULE_PATH)."
        ),
    )
    p.add_argument(
        "--et-var", default=ET_VAR, choices=["ACT. ETRA", "ACT. ETRA_patched"],
        metavar="VAR",
        help=(
            "Which ET variable to plot/aggregate throughout: 'ACT. ETRA' "
            "(raw solver output) or 'ACT. ETRA_patched' (artefact-floored, "
            "see apply_eta_artefact_floor in Agramon_withLAI_withETp.py). "
            "Overrides the ET_VAR module constant for this run."
        ),
    )

    # Content selection
    p.add_argument(
        "--selected-year", type=int, default=2021, metavar="YEAR",
        help="Year to show in the monthly veg-map and LAI-map panels.",
    )
    p.add_argument(
        "--lai-vmax", type=float, default=LAI_VMAX, metavar="VAL",
        help="Upper bound for LAI colour scale.",
    )
    p.add_argument(
        "--sw-var", default="sw", metavar="VAR",
        help="Variable name inside sw_output.nc to plot (default: 'sw').",
    )
    p.add_argument(
        "--sw-vmin", type=float, default=0.0, metavar="VAL",
        help="Lower bound for the surface-saturation colour scale (Fig 4b).",
    )
    p.add_argument(
        "--sw-vmax", type=float, default=1.0, metavar="VAL",
        help="Upper bound for the surface-saturation colour scale (Fig 4b).",
    )

    # Date filtering (Fig 6 / Fig 7 hydro plots)
    p.add_argument("--date-start", default=None, metavar="YYYY-MM-DD",
                   help="Start of the plotting window for hydro/sw/psi plots (inclusive).")
    p.add_argument("--date-end",   default=None, metavar="YYYY-MM-DD",
                   help="End of the plotting window for hydro/sw/psi plots (inclusive).")

    # Reference node positions (Fig 6 / Fig 7)
    p.add_argument("--outlet-xy", nargs=2, type=float, default=None, metavar=("X", "Y"),
                   help="X Y coordinates of the outlet reference node "
                        "(auto-selected as lowest surface node if omitted).")
    p.add_argument("--uphill-xy", nargs=2, type=float, default=None, metavar=("X", "Y"),
                   help="X Y coordinates of the mid-uphill reference node "
                        "(auto-selected as mid-elevation node if omitted).")
    p.add_argument("--depths", nargs="+", type=float, default=[0, 1, 2], metavar="M",
                   help="Depth offsets (m) below each reference surface node.")

    # Panel toggles for Fig 6
    p.add_argument("--no-rain", dest="show_rain", action="store_false")
    p.add_argument("--no-psi",  dest="show_psi",  action="store_false")
    p.add_argument("--no-sw",   dest="show_sw",   action="store_false")
    p.set_defaults(show_rain=True, show_psi=True, show_sw=True)

    # Series toggles for Fig 6 / Fig 7
    p.add_argument("--no-uphill-surf",   dest="show_uh_s", action="store_false")
    p.add_argument("--no-uphill-1m",     dest="show_uh_1", action="store_false")
    p.add_argument("--no-downhill-surf", dest="show_dh_s", action="store_false")
    p.add_argument("--no-downhill-1m",   dest="show_dh_1", action="store_false")
    p.set_defaults(show_uh_s=True, show_uh_1=True, show_dh_s=True, show_dh_1=True)
    p.add_argument(
        "--no-gifs", action="store_true",
        help="Skip GIF animation generation (can be slow).",
    )

    # Export
    p.add_argument("--dpi", type=int, default=150, metavar="DPI",
                   help="Resolution for all exported PNG figures.")
    p.add_argument("--no-plots", action="store_true",
                   help="Do not call plt.show(); save to disk only.")
    p.add_argument("--backend", default=None, metavar="BACKEND",
                   help="Matplotlib backend override (e.g. TkAgg).")

    return p.parse_args()


# =============================================================================
# Log helpers
# =============================================================================

def _resolve_log(path2prj: Path, log_file: str) -> Path:
    lp = Path(log_file)
    if not lp.is_absolute():
        # .resolve() elimina los '..' y convierte la ruta en absoluta
        lp = (path2prj / 'Codigos_HydroModel' / lp).resolve()
    else:
        lp = lp.resolve()

    if not lp.exists():
        raise FileNotFoundError(f"Log file not found: {lp}")
    return lp


def _find_scenario_dir(base_dir: Path, sim_index: int) -> Path:
    """
    Locate the folder for *sim_index* under *base_dir*.

    Simulation runs (Agramon_withLAI_withETp.py) name scenario folders
    ``scenario_<sim_index>`` (older runs, before this tagging existed) or
    ``scenario_<sim_index>_<tag>``, where <tag> describes the dem-plot /
    outlet / bottom / spin-up combination used for that run — see
    build_scenario_dirname() in Agramon_withLAI_withETp.py. Matching is
    anchored on sim_index as a whole underscore-delimited segment so e.g.
    sim_index=2 never matches a folder belonging to sim_index=20.

    Raises FileNotFoundError if nothing (or, after a warning, more than
    one candidate) matches.
    """
    if not base_dir.exists():
        raise FileNotFoundError(f"{base_dir} does not exist.")

    pattern = re.compile(rf"^scenario_{sim_index}(_.+)?$")
    matches = sorted(
        p for p in base_dir.iterdir() if p.is_dir() and pattern.match(p.name)
    )

    if not matches:
        available = sorted(p.name for p in base_dir.iterdir() if p.is_dir())
        raise FileNotFoundError(
            f"No folder matching 'scenario_{sim_index}' (optionally "
            f"'scenario_{sim_index}_<tag>') found under {base_dir}.\n"
            f"  Available: {available}"
        )
    if len(matches) > 1:
        print(f"  WARNING: multiple folders match scenario_{sim_index} under "
              f"{base_dir} — using '{matches[0].name}'. "
              f"All matches: {[m.name for m in matches]}")
    return matches[0]


def _load_log_row(path2prj: Path, log_file: str, scenario: int | None) -> tuple[int, pd.Series]:
    """
    Return (sim_index, log_row) for the requested withLAI scenario.

    If *scenario* is None, the first row where ``withLAI == 1`` is used.
    Raises ValueError when no matching row is found.
    """
    log_path = _resolve_log(path2prj, log_file)
    # log_path = _resolve_log(path2prj/args.log_file)

    df_log   = pd.read_csv(log_path)

    # Filter to withLAI=1 runs when the column exists
    lai_col = "withLAI" if "withLAI" in df_log.columns else None
    df_lai  = df_log[df_log[lai_col] == 1] if lai_col else df_log

    if df_lai.empty:
        raise ValueError(
            f"No withLAI=1 rows found in {log_path}. "
            "Run LT_SShydro_run_withLAI.py first."
        )

    sim_index = scenario if scenario is not None else int(df_lai["sim_index"].iloc[0])

    row = df_lai[df_lai["sim_index"] == sim_index]
    if row.empty:
        raise ValueError(
            f"Scenario {sim_index} not found (or not withLAI=1) in {log_path}. "
            f"Available withLAI scenarios: {list(df_lai['sim_index'])}"
        )
    return sim_index, row.iloc[0]


# =============================================================================
# Artefact loading
# =============================================================================

def _require(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Required artefact not found: {path}\n"
            f"  -> Run LT_SShydro_run_withLAI.py first to generate '{label}'."
        )


def load_artefacts(out_dir: Path, path2prj: Path | None = None, sim_index: int | None = None) -> dict:
    """Load all pre-computed NetCDF / pickle artefacts from *out_dir*.

    If *path2prj* and *sim_index* are provided the function also tries to
    instantiate a pyCATHY CATHY object so that ``simu.read_outputs("atmbc")``
    can be called for the hydro / atmbc plots (Fig 6 / 7).  If the project
    folder is not found the simu key is set to None and those figures are
    skipped gracefully.
    """
    print("Loading saved artefacts ...")

    # ET
    et_nc = out_dir / "et_output.nc"
    _require(et_nc, "et_output.nc")
    ET_xr_all = xr.open_dataset(et_nc)
    print(f"  ET      : {et_nc}  |  variables: {list(ET_xr_all.data_vars)}")

    # PSI
    psi_nc = out_dir / "psi_output.nc"
    _require(psi_nc, "psi_output.nc")
    psi_xr = xr.open_dataset(psi_nc)
    print(f"  PSI     : {psi_nc}")

    # SW
    sw_nc = out_dir / "sw_output.nc"
    _require(sw_nc, "sw_output.nc")
    sw_xr = xr.open_dataset(sw_nc)
    print(f"  SW      : {sw_nc}")

    # Vegetation map history
    veg_nc = out_dir / "veg_map_history.nc"
    _require(veg_nc, "veg_map_history.nc")
    veg_xr_full = xr.open_dataarray(veg_nc)
    sorted_dates = [pd.Timestamp(t).to_pydatetime() for t in veg_xr_full["time"].values]
    veg_map_history = {
        pd.Timestamp(t).to_pydatetime(): veg_xr_full.sel(time=t).values
        for t in veg_xr_full["time"].values
    }
    N, M = veg_xr_full.sizes["row"], veg_xr_full.sizes["col"]
    print(f"  VegMap  : {veg_nc}  |  {len(sorted_dates)} months, grid {N}x{M}")

    # Recharge (optional)
    # recharge_csv = out_dir / "recharge_monthly.csv"
    # if recharge_csv.exists():
    #     recharge_monthly = pd.read_csv(
    #         recharge_csv, index_col=0, parse_dates=True
    #     ).squeeze("columns")
    #     print(f"  Recharge: {recharge_csv}")
    # else:
    #     recharge_monthly = pd.Series(dtype=float)
    #     print("  Recharge: recharge_monthly.csv not found -- Fig 1 will be skipped.")

    # Recharge NetCDF (preferred over simu.read_outputs for Fig 10)
    recharge_nc = out_dir / "recharge_output.nc"
    if recharge_nc.exists():
        xr_recharge = xr.open_dataset(recharge_nc)
        print(f"  Recharge NC: {recharge_nc}  |  variables: {list(xr_recharge.data_vars)}")
    else:
        xr_recharge = None
        print("  Recharge NC: recharge_output.nc not found -- Fig 10 will fall back to simu.")

    # LAI history pickle (optional)
    lai_pkl = out_dir / "veg_map_history.pkl"
    lai_history: dict = {}
    if lai_pkl.exists():
        with open(lai_pkl, "rb") as fh:
            lai_history = pickle.load(fh)
        print(f"  LAI pkl : {lai_pkl}  |  {len(lai_history)} months")
    else:
        print("  LAI pkl : veg_map_history.pkl not found -- LAI scatter will be skipped.")

    print("Artefacts loaded.\n")

    # -- CATHY project object (optional — needed for grid3d / node resolution) --
    # psi, sw, and ET data come from the .nc artefacts above; we do NOT reload
    # them via simu.read_outputs to avoid a second (slow) read pass.
    simu = None

    if path2prj is not None and sim_index is not None:
        try:
            prj_folder = _find_scenario_dir(path2prj, sim_index)
        except FileNotFoundError:
            prj_folder = None
            print(f"  WARNING: No project folder matching 'scenario_{sim_index}' found under "
                  f"{path2prj} — Fig 6/7 node resolution will be skipped.")

        if prj_folder is not None:
            try:
                simu = cathy_tools.CATHY(dirName=path2prj, prj_name=prj_folder/'my_cathy_prj')
                print(f"  CATHY prj : {prj_folder.name}  (used for node resolution only)")
            except Exception as exc:
                print(f"  WARNING: Could not load CATHY project — "
                      f"Fig 6/7 node resolution will be skipped.\n    ({exc})")
                simu = None

    return dict(
        ET_xr_all=ET_xr_all,
        psi_xr=psi_xr,
        sw_xr=sw_xr,
        veg_map_history=veg_map_history,
        sorted_dates=sorted_dates,
        # recharge_monthly=recharge_monthly,
        xr_recharge=xr_recharge,
        lai_history=lai_history,
        simu=simu,
    )


# =============================================================================
# Node resolution (mirrors LT_SShydro_results.py)
# =============================================================================

def resolve_nodes(
    simu,
    grid3d: dict,
    args: argparse.Namespace,
    *,
    psi_xr: "xr.Dataset | None" = None,
    sw_xr:  "xr.Dataset | None" = None,
) -> dict:
    """
    Identify outlet and mid-uphill reference nodes at each requested depth.

    The outlet is always a single reference node. The mid-uphill site is,
    by default, a *band* of nodes within ±2% elevation of the domain's
    mid-elevation (see ``results_plotter.find_representative_mid_node(...,
    return_band=True)``) rather than a single arbitrarily-chosen point —
    downstream plotting takes the per-timestep median across the band, so
    "mid elevation" reflects the whole band. If the caller pins an exact
    XY location via ``args.uphill_xy``, that single point is used instead
    (a band doesn't make sense once the user has asked for one specific
    location).

    Parameters
    ----------
    simu     : pyCATHY CATHY object (already loaded).
    grid3d   : dict returned by ``simu.read_outputs("grid3d")``.
    args     : parsed CLI namespace; uses ``outlet_xy``, ``uphill_xy``,
               and ``depths``.
    psi_xr / sw_xr : optional pressure-head / soil-water xr.Datasets
               (already loaded, same ones passed to Fig 6/7). When given,
               the console summary also reports min/max/mean ψ and sw
               across every node in each band, over the full simulation —
               a quick sanity check that band nodes aren't just numerically
               identical to each other.

    Returns
    -------
    dict mapping label → node-index (or list[node-index] for a mid-uphill
    band), e.g.
        { "outlet_z-0m": 142, "outlet_z-1m": 998,
          "mid_z-0m": [57, 61, 64], "mid_z-1m": [1011, 1015, 1018], … }
    """
    nnod          = int(grid3d["nnod"])
    nodes         = grid3d["mesh3d_nodes"]
    surface_nodes = nodes[:nnod]

    # Outlet — lowest surface node or user-supplied XY
    if args.outlet_xy:
        ox, oy = args.outlet_xy
        dists  = np.hypot(surface_nodes[:, 0] - ox, surface_nodes[:, 1] - oy)
        outlet_idx = int(dists.argmin())
    else:
        outlet_idx = int(surface_nodes[:, 2].argmin())
    outlet_x, outlet_y, outlet_z = surface_nodes[outlet_idx]

    # Mid-uphill — a band of nodes around mid-elevation, or a single
    # user-supplied XY point.
    mid_band_idx: np.ndarray | None = None
    if args.uphill_xy:
        ux, uy = args.uphill_xy
        dists  = np.hypot(surface_nodes[:, 0] - ux, surface_nodes[:, 1] - uy)
        mid_idx = int(dists.argmin())
        mid_x, mid_y, mid_z = surface_nodes[mid_idx]
    else:
        mid_band_idx = plotter.find_representative_mid_node(
            surface_nodes, return_band=True
        )
        # Recentred single point too, purely for the console summary below
        # (the band itself is what actually gets used downstream).
        mid_idx = plotter.find_representative_mid_node(
            surface_nodes, return_band=False
        )
        mid_x, mid_y, mid_z = surface_nodes[mid_idx]

    print(f"  Outlet:     X={outlet_x:.2f}, Y={outlet_y:.2f}, Z={outlet_z:.2f}")
    print(f"  Mid-uphill: X={mid_x:.2f},   Y={mid_y:.2f},   Z={mid_z:.2f}", end="")
    if mid_band_idx is not None:
        band_z = surface_nodes[mid_band_idx, 2]
        print(
            f"  (band: {len(mid_band_idx)} node(s), "
            f"Z range {band_z.min():.2f}–{band_z.max():.2f})"
        )
    else:
        print("  (single user-supplied point)")

    nodes_dict: dict = {}
    for depth in args.depths:
        depth_label = int(depth) if float(depth).is_integer() else depth
        nodes_dict[f"outlet_z-{depth_label}m"] = plotter.unwrap_node_id(
            simu.find_nearest_node([outlet_x, outlet_y, outlet_z - depth])
        )
        if mid_band_idx is not None:
            if depth == 0:
                # Surface band indices already are node indices at depth 0.
                nodes_dict[f"mid_z-{depth_label}m"] = [int(i) for i in mid_band_idx]
            else:
                # Project every surface node in the band straight down by
                # `depth` and resolve each one to its own subsurface node,
                # so the band is preserved at every requested depth too.
                nodes_dict[f"mid_z-{depth_label}m"] = [
                    plotter.unwrap_node_id(
                        simu.find_nearest_node(
                            [surface_nodes[i, 0], surface_nodes[i, 1],
                             surface_nodes[i, 2] - depth]
                        )
                    )
                    for i in mid_band_idx
                ]
        else:
            nodes_dict[f"mid_z-{depth_label}m"] = plotter.unwrap_node_id(
                simu.find_nearest_node([mid_x, mid_y, mid_z - depth])
            )

    print("  Resolved nodes:")

    def _var_stats(ds: "xr.Dataset | None", ids) -> str:
        """Min/max/mean of ds's data variable across every id in `ids`
        (band or single node) and across the whole time series."""
        if ds is None:
            return ""
        var      = list(ds.data_vars)[0]
        time_dim = plotter._infer_time_dim(ds)
        da       = ds[var]
        sel_ids  = ids if isinstance(ids, list) else [ids]
        vals     = da.sel(node=sel_ids).values
        return (
            f" [{var}: min/max/mean = "
            f"{np.nanmin(vals):.3f}/{np.nanmax(vals):.3f}/{np.nanmean(vals):.3f}]"
        )

    for k, v in nodes_dict.items():
        if isinstance(v, list):
            v_arr = np.asarray(v)
            line = (
                f"    {k}: band of {len(v)} node(s) "
                f"(id min/max/mean = {v_arr.min()}/{v_arr.max()}/{v_arr.mean():.1f})"
            )
        else:
            line = f"    {k}: {v}"
        line += _var_stats(psi_xr, v) + _var_stats(sw_xr, v)
        print(line)
    return nodes_dict


# =============================================================================
# Pipeline
# =============================================================================

def run_pipeline(args: argparse.Namespace) -> None:
    global ET_VAR
    ET_VAR = args.et_var
    print(f"  ET variable in use : {ET_VAR}")

    # -- Backend ---------------------------------------------------------------
    if args.backend:
        matplotlib.use(args.backend)

    show     = not args.no_plots
    dpi      = args.dpi
    path2prj = Path(args.path2prj).resolve()

    # -- List-only mode --------------------------------------------------------
    if args.list_scenarios:
        # log_path = _resolve_log(path2prj/'..', args.log_file)
        log_path = _resolve_log(path2prj/args.log_file)
        df_log   = pd.read_csv(log_path)
        lai_col  = "withLAI" if "withLAI" in df_log.columns else None
        df_show  = df_log[df_log[lai_col] == 1] if lai_col else df_log
        print("Available withLAI simulations:")
        print(df_show.to_string(index=False))
        return

    # -- Resolve scenario from log ---------------------------------------------
    sim_index, log_row = _load_log_row(path2prj/'..', args.log_file, args.scenario)

    # -- Paths scoped to this simulation ---------------------------------------
    # Resolve the actual (possibly tagged) scenario folder the simulation
    # wrote to — e.g. scenario_7_plot2_outlet-free_bottom-closed_spinup-on3c
    # — rather than assuming the untagged 'scenario_<N>' name, so different
    # BC/spin-up combos for the same dem-plot don't collide and stay easy
    # to tell apart by folder name alone.
    out_dir = _find_scenario_dir(path2prj / "outputs", sim_index).resolve()
    scenario_dirname = out_dir.name

    # Figures reuse the same tagged name so a figures/ tree mirrors the
    # outputs/ tree one-to-one, and comparing scenarios visually is just a
    # matter of opening the two folders side by side.
    fig_dir = (
        Path(args.fig_dir).resolve()
        if args.fig_dir
        else path2prj / "figures" / scenario_dirname
    )
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'=' * 60}")
    print(f"  LAI Loop Results Viewer -- {scenario_dirname}")
    print(f"  Log row   : {log_row.to_dict()}")
    print(f"  Artefacts : {out_dir}")
    print(f"  Figures   : {fig_dir}")
    print(f"{'=' * 60}\n")

    # -- Load artefacts --------------------------------------------------------
    art = load_artefacts(out_dir, path2prj=path2prj, sim_index=sim_index)
    ET_xr_all        = art["ET_xr_all"]
    veg_map_history  = art["veg_map_history"]
    # recharge_monthly = art["recharge_monthly"]
    simu             = art["simu"]
    # (psi_output.nc, sw_output.nc, et_output.nc) — no pyCATHY simu reload.
    psi_xr = art["psi_xr"]
    sw_xr  = art["sw_xr"]

    # Cache grid3d early (needed below for the surface-saturation spatial
    # grid, and reused further down for node resolution / Fig 6b).
    grid3d = None
    if simu is not None:
        try:
            grid3d = simu.read_outputs('grid3d')
        except Exception as exc:
            print(f"  WARNING: grid3d read failed — surface-saturation "
                  f"spatial maps (Fig 4b) will be skipped.\n    ({exc})")

    # psi_xr_surf = art["psi_xr"].isel(node=slice(0,1060))
    # sw_xr_surf = art["sw_xr"].isel(node=slice(0,1060))

    # fig, ax = plt.subplots()
    # psi_xr_surf['psi'].isel(node=1).plot.scatter()
    # fig.savefig(fig_dir/'psi.png')


    # fig, ax = plt.subplots()
    # sw_xr_surf['sw'].isel(node=1).plot.scatter()
    # fig.savefig(fig_dir/'sw.png')



    #%%

    # %matplotlib inline
    # fig, ax = plt.subplots()
    # ET_xr_all['ACT. ETRA'].isel(datetime=0).plot.imshow(ax=ax)
    # plt.show()

    # -- DEM -------------------------------------------------------------------
    adf_folder = f"{MODULE_PATH}/DTMplots/dtmplot{log_row['watershed_nb']}/"
    
    gdf_Agramon = AgUtils.load_plot_shapefiles(['microcuencas_13'], MODULE_PATH / "shapefiles")
    dem_tif_path = Path(args.dem_tif)  # same CLI default/value used at generation time

    raster_DEM, raster_DEM_masked, xllcorner, yllcorner, res_x, res_y = load_dem_from_tif(
        dem_tif_path,
        gdf_clip=gdf_Agramon,
        fid=int(log_row["watershed_nb"]),
        resample_resolution=2,
    )


    # raster_DEM, raster_DEM_masked, xllcorner, yllcorner, res_x, res_y = AgUtils.load_dem(
    #    TARGET_CRS=TARGET_CRS,
    #    adf_folder=adf_folder,
    #)

    veg_map_history_masked = veg_map_history.copy()

    for vmh in veg_map_history:
        print(vmh)

        veg_map_history_masked[vmh] = np.where(
            raster_DEM_masked == -9999,
            np.nan,
            veg_map_history[vmh]
        )

    # -- Fig 1: Monthly recharge -----------------------------------------------
    print("-- Fig 1: Monthly recharge ------------------------------------------")
    # plotter.plot_recharge_bar(recharge_monthly, fig_dir, dpi=dpi, show=show)

    # -- Fig 2: Monthly mean ET by year ----------------------------------------
    print("-- Fig 2: Monthly mean ET -------------------------------------------")
    plotter.plot_et_monthly_mean(ET_xr_all, ET_VAR, ET_SCALE, fig_dir, dpi=dpi, show=show)

    # -- Fig 3: ET by vegetation type ------------------------------------------
    print("-- Fig 3: ET monthly trend per vegetation type ----------------------")
    # plotter.plot_et_by_veg_type(
    #     ET_xr_all, veg_map_history_masked, ET_VAR, ET_SCALE, fig_dir, dpi=dpi, show=show
    # )

    # -- Fig 4: Spatial ET maps per month for selected year --------------------
    print(f"-- Fig 4: Spatial ET monthly maps -- {args.selected_year} ----------")

    plot_et_spatial_monthly(
        ET_xr_all, ET_VAR, ET_SCALE,
        2019,
        fig_dir,
        cmap=CMAP_ET,
        dpi=dpi,
        show=show,
    )
    plot_et_spatial_monthly(
        ET_xr_all, ET_VAR, ET_SCALE,
        2021,
        fig_dir,
        cmap=CMAP_ET,
        dpi=dpi,
        show=show,
    )

    # -- Fig 4b: Spatial surface soil-moisture saturation maps per month ------
    print(f"-- Fig 4b: Spatial surface saturation monthly maps -- {args.selected_year} ----------")
    if grid3d is not None:
        try:
            ds_sw_grid = build_sw_surface_grid(
                sw_xr, grid3d, ET_xr_all, sw_var=args.sw_var,
            )
            for yr in (2019, args.selected_year):
                plot_sw_spatial_monthly(
                    ds_sw_grid,
                    variable=args.sw_var,
                    selected_year=yr,
                    out_dir=fig_dir,
                    cmap=CMAP_SW,
                    dpi=dpi,
                    show=show,
                    vmin=args.sw_vmin,
                    vmax=args.sw_vmax,
                )
        except Exception as exc:
            print(f"  WARNING: Surface saturation spatial maps skipped — {exc}")
    else:
        print("  Skipped — grid3d unavailable (CATHY project not loaded).")

    # -- Fig 5: Vegetation map samples (one per year) --------------------------
    print("-- Fig 5: Vegetation map samples -- by year -------------------------")
    plotter.plot_veg_map_samples(veg_map_history_masked, CMAP_VEG, fig_dir, dpi=dpi, show=show)

    # -- Fig 5b: Monthly vegetation maps for selected year ---------------------
    print(f"-- Fig 5b: Vegetation map -- {args.selected_year} (monthly) ----------")
    plotter.plot_veg_map_year(
        veg_map_history_masked, 2016, CMAP_VEG, fig_dir, dpi=dpi, show=show
    )
    plotter.plot_veg_map_year(
        veg_map_history_masked, args.selected_year, CMAP_VEG, fig_dir, dpi=dpi, show=show
    )

    # -- Fig 5c: Monthly LAI maps for selected year ----------------------------
    print(f"-- Fig 5c: LAI maps -- {args.selected_year} (monthly) ----------------")
    ds_lai = xr.open_dataset(LAI_PATH / "lai_monthly.nc")
    plotter.plot_lai_map_year(
        ds_lai, LAI_VAR, args.selected_year,
        LAI_CMAP, LAI_VMIN, args.lai_vmax,
        fig_dir, dpi=dpi, show=show,
    )
    ds_lai.close()


# def plot_lai_monthly_by_year(
#     lai_history: dict,
#     lai_var:     str         = "LAI",
#     out_dir:     Path | None = None,
#     fname:       str         = "LAI_monthly_mean_by_year.png",
#     scenario_id: int | str   = "",
#     dpi:         int         = 150,
#     show:        bool        = False,
# ) -> None:
#     """


    # -- LAI-A : monthly mean by year ------------------------------------------
    print("-- veg-A: veg indice monthly mean by year ----------------------------------")
    try:
        plot_index_monthly_by_year(
            art["lai_history"],
            index_var   = LAI_VAR,
            out_dir     = fig_dir,
            scenario_id = sim_index,
            dpi         = dpi,
            show        = show,
        )
    except Exception as exc:
        print(f"  WARNING: plot_index_monthly_by_year skipped — {exc}")

    # -- LAI-B : spatial maps for selected year --------------------------------
    print(f"-- LAI-B: LAI spatial maps -- {args.selected_year} -----------------")
    try:
        plot_index_spatial_maps_year(
            art["lai_history"],
            index_var     = LAI_VAR,
            selected_year = args.selected_year,
            out_dir       = fig_dir,
            scenario_id   = sim_index,
            dpi           = dpi,
            show          = show,
        )
    except Exception as exc:
        print(f"  WARNING: plot_index_spatial_maps_year skipped — {exc}")

    # -- LAI-C : per-pixel annual mean histogram --------------------------------
    print("-- LAI-C: LAI annual mean histogram ---------------------------------")
    try:
        plot_index_annual_hist(
            art["lai_history"],
            index_var   = LAI_VAR,
            out_dir     = fig_dir,
            scenario_id = sim_index,
            dpi         = dpi,
            show        = show,
        )
    except Exception as exc:
        print(f"  WARNING: plot_index_annual_hist skipped — {exc}")

    # -- LAI-D : domain-mean annual bar chart ----------------------------------
    print("-- LAI-D: LAI domain-mean annual bar chart --------------------------")
    try:
        plot_index_annual_mean(
            art["lai_history"],
            index_var   = LAI_VAR,
            out_dir     = fig_dir,
            scenario_id = sim_index,
            dpi         = dpi,
            show        = show,
        )
    except Exception as exc:
        print(f"  WARNING: plot_index_annual_mean skipped — {exc}")

    # -- LAI-E : domain-median annual time series ------------------------------
    print("-- LAI-E: LAI domain-median yearly time series ----------------------")
    try:
        plot_index_annual_median(
            art["lai_history"],
            index_var   = LAI_VAR,
            out_dir     = fig_dir,
            scenario_id = sim_index,
            dpi         = dpi,
            show        = show,
        )
    except Exception as exc:
        print(f"  WARNING: plot_index_annual_median skipped — {exc}")

    # -- GIFs ------------------------------------------------------------------
    if args.no_gifs:
        print("-- GIF A: Animated vegetation maps ----------------------------------")
        plotter.make_veg_gif(veg_map_history_masked, CMAP_VEG, fig_dir / "veg_map_animation.gif")

        print("-- GIF B: Animated ET maps ------------------------------------------")
        if ET_VAR in ET_xr_all:
            plotter.make_et_gif(
                ET_xr_all, ET_VAR, ET_SCALE,
                fig_dir / "ET_animation.gif",
                cmap=CMAP_ET,
            )
        else:
            print(f"  Skipped -- '{ET_VAR}' not in ET dataset.")
    else:
        print("-- GIFs skipped (--no-gifs) -----------------------------------------")

    # -- Fig 6: Hydro time-series (Rain/ETp + ψ + sw) -------------------------
    # -- Fig 7: sw + Actual ETa time-series ------------------------------------
    # Both figures are built directly from the xarray artefacts already loaded

    print("-- Resolving reference nodes for Fig 6 / 7 -------------------------")
    if simu is not None and grid3d is not None:
        # grid3d was already fetched above (for Fig 4b); simu is still
        # needed here for find_nearest_node() inside resolve_nodes().
        try:
            nodes_dict = resolve_nodes(simu, grid3d, args, psi_xr=psi_xr, sw_xr=sw_xr)
        except Exception as exc:
            print(f"  WARNING: Node resolution failed — Fig 6/7 skipped.\n    ({exc})")
            nodes_dict = None
    else:
        print("  WARNING: CATHY project folder not found — "
              "node resolution skipped, Fig 6/7 will be skipped.")
        nodes_dict = None

    if nodes_dict is not None:
        # Load satellite ETp + rain forcing for the rain/ETp panel in Fig 6
        # (replaces the ERA5 point time-series previously used here — see
        # Agramon_withLAI_withETp.py, step 1: "Load satellite ETp + rain
        # NetCDFs (EPSG:32630 → EPSG:25830)").
        pev_series: pd.Series | None = None
        tp_series:  pd.Series | None = None
        try:
            ds_etp  = load_etp()
            ds_rain = load_rain()
            gdf_Agramon = AgUtils.load_plot_shapefiles(['microcuencas_13'],
                                                       MODULE_PATH / "shapefiles")
            pev_series = _domain_mean_daily_series(ds_etp,  "ET_0-gf", gdf_Agramon)
            tp_series  = _domain_mean_daily_series(ds_rain, "TP-DD",   gdf_Agramon)

            # ET0 and TP come from separate NetCDFs and don't necessarily
            # share the same valid-date coverage (e.g. 3014 vs 3620 days)
            # — align to their common dates so every panel/mask downstream
            # sees two series of equal length.
            common_idx = pev_series.index.intersection(tp_series.index)
            n_dropped_pev = len(pev_series) - len(common_idx)
            n_dropped_tp  = len(tp_series)  - len(common_idx)
            if n_dropped_pev or n_dropped_tp:
                print(f"  Aligning ETp/rain to common dates: "
                      f"{len(common_idx)} shared days "
                      f"(dropped {n_dropped_pev} ETp-only, "
                      f"{n_dropped_tp} rain-only day(s)).")
            pev_series = pev_series.loc[common_idx]
            tp_series  = tp_series.loc[common_idx]
        except Exception as etp_exc:
            print(f"  WARNING: satellite ETp/rain forcing not loaded — "
                  f"rain/ETp panel will be empty.\n    ({etp_exc})")
#%%
        print("-- Fig 6: Hydro time-series (Rain/ETp + ψ + sw) -----------------")
        plotter.plot_hydro_xr(
            psi_xr, sw_xr, nodes_dict, fig_dir,
            pev_series=pev_series,
            tp_series=tp_series,
            date_start=args.date_start,
            date_end=args.date_end,
            show_rain=args.show_rain,
            show_psi=args.show_psi,
            show_sw=args.show_sw,
            show_uh_s=args.show_uh_s,
            show_uh_1=args.show_uh_1,
            show_dh_s=args.show_dh_s,
            show_dh_1=args.show_dh_1,
            scenario_id=sim_index,
            dpi=dpi, show=show,
        )

        print(f"-- Fig 6c: Hydro time-series -- fire year ({plotter.FIRE_DATE.year}) ---")
        try:
            plotter.plot_hydro_xr_fire_year(
                psi_xr, sw_xr, nodes_dict, fig_dir,
                pev_series=pev_series,
                tp_series=tp_series,
                show_rain=args.show_rain,
                show_psi=args.show_psi,
                show_sw=args.show_sw,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_hydro_xr_fire_year skipped — {exc}")

        print(f"-- Fig 6d: Hydro time-series -- fire month "
              f"({plotter.FIRE_DATE.strftime('%B %Y')}) -----------------")
        try:
            plotter.plot_hydro_xr_fire_month(
                psi_xr, sw_xr, nodes_dict, fig_dir,
                pev_series=pev_series,
                tp_series=tp_series,
                show_rain=args.show_rain,
                show_psi=args.show_psi,
                show_sw=args.show_sw,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_hydro_xr_fire_month skipped — {exc}")

        print("-- Fig 6e: Hydro time-series -- weekly ---------------------------")
        try:
            plotter.plot_hydro_xr_weekly(
                psi_xr, sw_xr, nodes_dict, fig_dir,
                pev_series=pev_series,
                tp_series=tp_series,
                date_start=args.date_start,
                date_end=args.date_end,
                show_rain=args.show_rain,
                show_psi=args.show_psi,
                show_sw=args.show_sw,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_hydro_xr_weekly skipped — {exc}")
#%%
        print("-- Fig 6b: Reference node locations (map) -----------------------")
        try:
            dem_da = None
            try:
                # Mask nodata (-9999) and reduce to a plain 2-D DataArray so
                # plot_node_locations can drop it straight into imshow().
                dem_da = raster_DEM.where(raster_DEM_masked != -9999)
                if hasattr(dem_da, "data_vars"):          # Dataset -> DataArray
                    dem_da = dem_da[list(dem_da.data_vars)[0]]
                if "band" in dem_da.dims:
                    dem_da = dem_da.squeeze("band", drop=True)
            except Exception as dem_exc:
                print(f"  WARNING: DEM basemap unavailable for Fig 6b ({dem_exc}).")
                dem_da = None

            plotter.plot_node_locations(
                grid3d, nodes_dict, fig_dir,
                dem_da=dem_da,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_node_locations skipped — {exc}")

        print("-- Fig 7: sw + Actual ETa time-series ---------------------------")
        plotter.plot_et_sw_xr(
            sw_xr, ET_xr_all, ET_VAR, ET_SCALE, nodes_dict, fig_dir,
            date_start=args.date_start,
            date_end=args.date_end,
            show_uh_s=args.show_uh_s,
            show_uh_1=args.show_uh_1,
            show_dh_s=args.show_dh_s,
            show_dh_1=args.show_dh_1,
            scenario_id=sim_index,
            dpi=dpi, show=show,
        )

        print("-- Fig 7b: sw + ETa + LAI trend time-series ----------------------")
        lai_ds = plotter.lai_history_to_dataset(art["lai_history"], lai_var=LAI_VAR)
        try:
            plotter.plot_et_sw_lai_xr(
                sw_xr, ET_xr_all, ET_VAR, ET_SCALE, lai_ds, nodes_dict, fig_dir,
                lai_var=LAI_VAR,
                date_start=args.date_start,
                date_end=args.date_end,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_et_sw_lai_xr skipped — {exc}")

        print(f"-- Fig 7bc: sw + ETa + LAI time-series -- fire years ({plotter.FIRE_DATE.year}-{plotter.FIRE_DATE.year + 1}) ---")
        try:
            plotter.plot_et_sw_lai_xr_fire_year(
                sw_xr, ET_xr_all, ET_VAR, ET_SCALE, lai_ds, nodes_dict, fig_dir,
                lai_var=LAI_VAR,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_et_sw_lai_xr_fire_year skipped — {exc}")

        print(f"-- Fig 7bd: sw + ETa + LAI time-series -- fire month "
              f"({plotter.FIRE_DATE.strftime('%B %Y')}) -----------------")
        try:
            plotter.plot_et_sw_lai_xr_fire_month(
                sw_xr, ET_xr_all, ET_VAR, ET_SCALE, lai_ds, nodes_dict, fig_dir,
                lai_var=LAI_VAR,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_et_sw_lai_xr_fire_month skipped — {exc}")

        print("-- Fig 7be: sw + ETa + LAI time-series -- weekly ------------------")
        try:
            plotter.plot_et_sw_lai_xr_weekly(
                sw_xr, ET_xr_all, ET_VAR, ET_SCALE, lai_ds, nodes_dict, fig_dir,
                lai_var=LAI_VAR,
                date_start=args.date_start,
                date_end=args.date_end,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_et_sw_lai_xr_weekly skipped — {exc}")

        print(f"-- Fig 7c: sw + ETa time-series -- fire years ({plotter.FIRE_DATE.year}-{plotter.FIRE_DATE.year + 1}) ---")
        try:
            plotter.plot_et_sw_xr_fire_year(
                sw_xr, ET_xr_all, ET_VAR, ET_SCALE, nodes_dict, fig_dir,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_et_sw_xr_fire_year skipped — {exc}")

        print(f"-- Fig 7d: sw + ETa time-series -- fire month "
              f"({plotter.FIRE_DATE.strftime('%B %Y')}) -----------------")
        try:
            plotter.plot_et_sw_xr_fire_month(
                sw_xr, ET_xr_all, ET_VAR, ET_SCALE, nodes_dict, fig_dir,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_et_sw_xr_fire_month skipped — {exc}")

        print("-- Fig 7e: sw + ETa time-series -- weekly ------------------------")
        try:
            plotter.plot_et_sw_xr_weekly(
                sw_xr, ET_xr_all, ET_VAR, ET_SCALE, nodes_dict, fig_dir,
                date_start=args.date_start,
                date_end=args.date_end,
                show_uh_s=args.show_uh_s,
                show_uh_1=args.show_uh_1,
                show_dh_s=args.show_dh_s,
                show_dh_1=args.show_dh_1,
                scenario_id=sim_index,
                dpi=dpi, show=show,
            )
        except Exception as exc:
            print(f"  WARNING: plot_et_sw_xr_weekly skipped — {exc}")

    # ── Fig 8: Catchment ET balance (monthly + annual) ─────────────────────
    # plot_catchment_et expects a data-dict with keys:
    #   ET_xr       — spatial xr.Dataset with a CRS (used for pixel-area)
    #   ET_xr_mean  — xr.Dataset/DataArray with spatial-mean ET_VAR values
    #   time_dates  — array-like of datetimes aligned with ET_xr_mean
    # (et_var=ET_VAR is passed explicitly below so this matches whichever
    #  variable name ET_xr_mean was actually built with, above.)
    # We build these from ET_xr_all which is already loaded.
    print("-- Fig 8: Catchment ET balance (monthly + annual) ----------------")
    # try:
    # Normalise: scale m/s → mm/day and rename 'datetime' → 'time'.
    # This also drops any conflicting 'time' coordinate that would
    # prevent the rename (the same fix applied to Figs 8b/8c).
    ET_xr_norm = _prepare_et_ds(ET_xr_all, ET_VAR, ET_SCALE)

    # Spatial-mean series (already mm/day after _prepare_et_ds)
    ET_xr_mean_da = ET_xr_norm[ET_VAR].mean(dim=["X", "Y"])
    ET_xr_mean_ds = ET_xr_mean_da.to_dataset(name=ET_VAR)

    # time_dates_arr = pd.DatetimeIndex(ET_xr_norm["time"].values).to_numpy()
    time_dates_arr = pd.to_datetime(ET_xr_all.datetime.values)

    # time_dates_arr = ET_xr_norm["time"].to_numpy()

    ET_xr_norm = ET_xr_norm.rename(
        {
            "X": "x",
            "Y": "y",
        }
    )
    data_catchment = dict(
        ET_xr      = ET_xr_norm,       # spatial dataset (X, Y, time)
        ET_xr_mean = ET_xr_mean_ds,
        time_dates = time_dates_arr,
    )
    plotter.plot_catchment_et(data_catchment,
                              fig_dir,
                              dpi=dpi,
                              show=show,
                              et_var=ET_VAR,
                              )


    # except Exception as exc:
    #     print(f"  WARNING: plot_catchment_et skipped — {exc}")

    # ── Fig 8b: ETa cumulated histogram (per-pixel annual sum) ────────────
    print("-- Fig 8b: ETa annual cumulative histogram -----------------------")
    # try:
    # plot_eta_annual_hist(
    #     ET_xr_all,
    #     variable=ET_VAR,
    #     scale=ET_SCALE,
    #     out_dir=fig_dir,
    #     scenario_id=sim_index,
    #     dpi=dpi,
    #     show=show,
    # )
    # except Exception as exc:
    #     print(f"  WARNING: plot_eta_annual_hist skipped — {exc}")

    # ── Fig 8c: ETa domain-total bar chart (monthly mean + annual total) ──
    print("-- Fig 8c: ETa domain-total monthly/annual bars ------------------")
    # try:
    # plot_eta_total(
    #     ET_xr_all,
    #     variable=ET_VAR,
    #     scale=ET_SCALE,
    #     out_dir=fig_dir,
    #     scenario_id=sim_index,
    #     dpi=dpi,
    #     show=show,
    # )
    # except Exception as exc:
    #     print(f"  WARNING: plot_eta_total skipped — {exc}")

    # ── Fig 9: Hydraulic graph ─────────────────────────────────────────────
    # Requires simu (pyCATHY object) to call simu.show(prop="hgraph").
    print("-- Fig 9: Hydraulic graph ----------------------------------------")
    if art["simu"] is not None:
        try:
            data_hgraph = dict(simu=art["simu"])
            plotter.plot_hgraph(data_hgraph, fig_dir, dpi=dpi, show=show)
        except Exception as exc:
            print(f"  WARNING: plot_hgraph skipped — {exc}")
    else:
        print("  Skipped — CATHY project not loaded (simu is None).")

    # ── Fig 10: Recharge (monthly + annual bar charts) ─────────────────────
    # Use recharge_output.nc saved by test_Agramon_withLAI.py instead of
    # simu.read_outputs("recharge"), so no live CATHY project is required.
    print("-- Fig 10: Recharge (recharge_output.nc) ----------------------------")
    xr_recharge = art.get("xr_recharge")
    if xr_recharge is not None:
        try:
            plotter.plot_recharge_xr(
                xr_recharge, fig_dir,
                dpi=dpi, show=show, scenario_id=sim_index,
            )
        except Exception as exc:
            print(f"  WARNING: plot_recharge_xr skipped — {exc}")
    elif art["simu"] is not None:
        # Fallback: live CATHY project if .nc is absent
        print("  recharge_output.nc not found — falling back to simu.read_outputs()")
        try:
            psi_time_dim = "datetime" if "datetime" in psi_xr.dims else "time"
            psi_var      = list(psi_xr.data_vars)[0]
            start_date_arr = psi_xr[psi_var][psi_time_dim].values[0]
            start_date_ts  = pd.Timestamp(start_date_arr)
            data_recharge = dict(
                simu=art["simu"],
                start_date=start_date_ts,
                scenario=sim_index,
            )
            plotter.plot_recharge(
                data_recharge, fig_dir,
                dpi=dpi, show=show, scenario_id=sim_index,
            )
        except Exception as exc:
            print(f"  WARNING: plot_recharge (fallback) skipped — {exc}")
    else:
        print("  Skipped — recharge_output.nc not found and CATHY project not loaded.")

    print(f"\n{'=' * 60}")
    print(f"  DONE -- scenario {sim_index} | all figures in {fig_dir}")
    print(f"{'=' * 60}")


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    args = _parse_args()
    run_pipeline(args)
