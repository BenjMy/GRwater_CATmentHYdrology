"""
analyse_ET_compare.py
=====================
GRwater project — Agramon study site

Compares two ETa sources side-by-side:
  A) Earth-Observation / Energy-Balance (EO-EB)
       — monthly NetCDF from the TSEB/biophysical pipeline
         (same file used by analyse_ET.py: 20161001_20241231_ET.nc)
  B) Water-Balance model (WB)
       — et_output.nc produced by LT_SShydro_run_withLAI.py
         (variable 'ACT. ETRA', units m/s, dims X/Y/datetime)

Figures produced
----------------
  compare_ET_monthly_mean_by_year_<label>.png
      Spatial-mean monthly ETa, one line per year — EO-EB vs WB on the same axes.

  compare_ET_annual_total_<label>.png
      Annual domain-mean ETa bar chart — EO-EB (blue) vs WB (orange).

  compare_ET_monthly_total_<label>.png
      Monthly domain-mean ETa bar chart — EO-EB vs WB grouped bars.

  compare_ET_scatter_<label>.png
      Scatter plot of coincident monthly spatial-mean values:
      EO-EB (x) vs WB (y), coloured by year.

  compare_ET_spatial_monthly_<year>_<label>.png  [one per source]
      Monthly spatial maps (mosaic) for a selected year — EO-EB and WB
      rendered with a shared colour scale for direct comparison.

Usage
-----
    python analyse_ET_compare.py                          # default plot 0, scenario 0
    python analyse_ET_compare.py --plots 0 2              # by index
    python analyse_ET_compare.py --plots Plot4Mulching    # by name
    python analyse_ET_compare.py --plots all              # all plots
    python analyse_ET_compare.py --scenario 3             # WB scenario index
    python analyse_ET_compare.py --year 2021              # spatial-map year
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray  # noqa: F401
import xarray as xr

# =============================================================================
# Paths & constants
# =============================================================================

EO_ET_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/GRwater/data/satellite/Agramon"
    "/20161001_20241231_ET.nc"
)

# Root of CATHY project — simulation artefacts live under <WB_PRJ_ROOT>/sim<N>/outputs/
# WB_PRJ_ROOT = Path(
#     "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
#     "/Tech4agro_org/GRwater_CATmentHYdrology"
# )
WB_PRJ_ROOT = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Training_Supervision/"
    "Supervision/Xela_Carracedo_Practicas_2026_data/SSHydro_withETp_withLAI"
)

# Shared CSV log produced by LT_SShydro_run_withLAI.py — same file used by
# withLAI_results.py to resolve a --scenario index into the actual (possibly
# tagged) scenario_<N>[_tag] output folder. Resolved under
# <WB_PRJ_ROOT>/../Codigos_HydroModel/<LOG_FILE> via _resolve_log().
LOG_FILE = "simulation_log_withETp_withLAI.csv"


MODULE_PATH = Path(
    "/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes"
    "/Tech4agro_org/GRwater_geophy"
).resolve()
if str(MODULE_PATH) not in sys.path:
    sys.path.append(str(MODULE_PATH))

import Agramon_utils as AgUtils

TARGET_CRS  = "EPSG:25830"
NATIVE_CRS  = "EPSG:32630"

#WB_ET_VAR   = "ACT. ETRA"
WB_ET_VAR   = "ACT. ETRA_patched"
WB_ET_SCALE = 1e3 * 86400        # m/s  → mm/day

CMAP_EO = "Blues"
CMAP_WB = "Oranges"
COLOR_EO = "#2166ac"
COLOR_WB = "#d6604d"

fig_path = Path().cwd() / "../Analysis/ET_compare"




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
# Shared utilities (mirrors analyse_ET.py helpers)
# =============================================================================

def load_shapefiles() -> gpd.GeoDataFrame:
    gdf = AgUtils.load_plot_shapefiles(['microcuencas_13'],
                                       MODULE_PATH / "shapefiles")
    if gdf.crs is None:
        raise ValueError("Shapefile has no CRS.")
    if str(gdf.crs) != TARGET_CRS:
        gdf = gdf.to_crs(TARGET_CRS)
    print(f"  Shapefile CRS : {gdf.crs}")
    for i, pid in enumerate(gdf["PlotID"]):
        print(f"    [{i}] {pid}")
    return gdf

def resolve_plot_selection(
    gdf: gpd.GeoDataFrame,
    selection: list[str | int] | None,
) -> gpd.GeoDataFrame:
    if selection is None or selection == ["all"]:
        return gdf.copy()

    # Ensure selection is a list
    if not isinstance(selection, list):
        selection = [selection]

    resolved_indices = []
    plot_ids = list(gdf["PlotID"])

    for item in selection:
        if isinstance(item, int) or (isinstance(item, str) and item.isdigit()):
            idx = int(item)
            if not (0 <= idx < len(gdf)):
                raise IndexError(f"Index {idx} out of range.")
            resolved_indices.append(idx)  # Store DataFrame index, not PlotID value
        else:
            if item not in plot_ids:
                raise ValueError(f"PlotID '{item}' not found. Available: {plot_ids}")
            matches = gdf[gdf["PlotID"] == item].index.tolist()
            resolved_indices.extend(matches)

    # Deduplicate indices
    seen: set[int] = set()
    unique_indices = [i for i in resolved_indices if not (i in seen or seen.add(i))]

    return gdf.iloc[unique_indices].copy().reset_index(drop=True)

def selection_label(gdf: gpd.GeoDataFrame, selection: list[str | int] | None) -> str:
    if selection is None or selection == ["all"]:
        return "all"

    if not isinstance(selection, list):
        selection = [selection]

    names: list[str] = []
    for item in selection:
        if isinstance(item, int) or (isinstance(item, str) and item.isdigit()):
            names.append(str(int(item)))  # Use the index directly in the label
        else:
            names.append(item)

    seen: set[str] = set()
    unique = [n for n in names if not (n in seen or seen.add(n))]
    return "_".join(unique)


def _clip_ds(ds: xr.Dataset, gdf: gpd.GeoDataFrame) -> xr.Dataset:
    clipped: dict[str, xr.DataArray] = {}
    for var in ds.data_vars:
        da = (
            ds[var]
            .rio.write_crs(TARGET_CRS, inplace=False)
            .rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
        )
        clipped[var] = da.rio.clip(gdf.geometry, gdf.crs, drop=True, all_touched=True)
    return xr.Dataset(clipped, attrs=ds.attrs)


# =============================================================================
# EO-EB dataset loader  (mirrors analyse_ET._load_nc)
# =============================================================================

def load_eo_et(path: Path = EO_ET_PATH) -> xr.Dataset:
    """Load EO ETa NetCDF, reproject to TARGET_CRS.  Returns dataset in mm/day."""
    ds = (
        xr.open_dataset(path, decode_times=True, mask_and_scale=True)
          .drop_vars("spatial_ref", errors="ignore")
    )
    ds = ds.rio.write_crs(NATIVE_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    if str(ds.rio.crs) != TARGET_CRS:
        ds = ds.rio.reproject(TARGET_CRS)
    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  EO-ET CRS : {ds.rio.crs}   vars: {list(ds.data_vars)}")
    return ds


# =============================================================================
# WB dataset loader  (reads et_output.nc from simulation artefacts)
# =============================================================================

def _find_sim_dir(scenario: int, log_file: str = LOG_FILE) -> Path:
    """
    Locate the WB simulation output folder for *scenario*.

    Mirrors withLAI_results.py's scenario-selection pattern instead of
    assuming a flat <WB_PRJ_ROOT>/sim<N>/outputs/ layout (which is not how
    LT_SShydro_run_withLAI.py actually writes its outputs):
      1. Resolve *scenario* to a sim_index via the shared log
         (<WB_PRJ_ROOT>/../Codigos_HydroModel/<log_file>).
      2. Find the matching (possibly tagged) scenario_<sim_index>[_tag]
         folder under <WB_PRJ_ROOT>/outputs/.
    """
    sim_index, log_row = _load_log_row(WB_PRJ_ROOT / "..", log_file, scenario)
    sim_dir = _find_scenario_dir(WB_PRJ_ROOT / "outputs", sim_index)
    print(f"  WB scenario  : sim_index={sim_index}  → {sim_dir.name}")
    return sim_dir


def load_wb_et(scenario: int = 0) -> xr.Dataset:
    """
    Load et_output.nc from a CATHY simulation, normalise to mm/day,
    and align spatial dims so rioxarray can clip it.

    The raw file has dims (X, Y, datetime) and values in m/s.
    Returns a Dataset with:
      • variable name preserved as WB_ET_VAR ('ACT. ETRA')
      • values scaled to mm/day
      • dims renamed to (x, y, time) for consistency with EO dataset
      • CRS written as TARGET_CRS (25830)
    """
    sim_dir = _find_sim_dir(scenario)
    nc_path = sim_dir / "et_output.nc"
    if not nc_path.exists():
        raise FileNotFoundError(f"et_output.nc not found in {sim_dir}")

    ds = xr.open_dataset(nc_path, decode_times=True, mask_and_scale=True)
    print(f"  WB-ET raw dims : {dict(ds.dims)}   vars: {list(ds.data_vars)}")

    if WB_ET_VAR not in ds:
        raise KeyError(
            f"Variable '{WB_ET_VAR}' not in WB dataset. "
            f"Available: {list(ds.data_vars)}"
        )

    # Scale m/s → mm/day
    ds = ds.assign({WB_ET_VAR: ds[WB_ET_VAR] * WB_ET_SCALE})

    # Normalise dimension names → (x, y, time)
    rename_map: dict[str, str] = {}
    if "X" in ds.dims:
        rename_map["X"] = "x"
    if "Y" in ds.dims:
        rename_map["Y"] = "y"
    if "datetime" in ds.dims and "time" not in ds.dims:
        # Drop any existing 'time' coord OR data variable that would
        # conflict with renaming the 'datetime' dim to 'time'. The raw
        # CATHY et_output.nc carries its own 'time' data var (seconds/
        # timestamp helper column) alongside the 'datetime' dim — using
        # ds.variables (not just ds.coords) catches that case too.
        if "time" in ds.variables:
            ds = ds.drop_vars("time")
        rename_map["datetime"] = "time"
    if rename_map:
        ds = ds.rename(rename_map)

    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  WB-ET final dims : {dict(ds.dims)}   CRS: {ds.rio.crs}")
    return ds


# =============================================================================
# Build common monthly spatial-mean time series
# =============================================================================

def _spatial_mean_monthly(ds: xr.Dataset, var: str, time_dim: str = "time") -> pd.Series:
    """Spatial mean (x,y) → resample to monthly mean → pd.Series."""
    sp_dims = [d for d in ds[var].dims if d != time_dim]
    da_mean = ds[var].mean(dim=sp_dims, skipna=True)
    df = pd.DataFrame(
        {var: da_mean.values},
        index=pd.DatetimeIndex(ds[time_dim].values),
    )
    return df[var].resample("ME").mean()


def build_common_series(
    ds_eo: xr.Dataset, var_eo: str,
    ds_wb: xr.Dataset,
) -> tuple[pd.Series, pd.Series]:
    """
    Return (s_eo, s_wb): coincident monthly spatial-mean ETa series (mm/day).
    Only months present in BOTH datasets are kept.
    """
    s_eo = _spatial_mean_monthly(ds_eo, var_eo, time_dim="time")
    s_wb = _spatial_mean_monthly(ds_wb, WB_ET_VAR, time_dim="time")

    # Align to the period common to both
    idx_common = s_eo.index.intersection(s_wb.index)
    if len(idx_common) == 0:
        raise ValueError(
            "EO-EB and WB datasets share no overlapping months.\n"
            f"EO range: {s_eo.index[0]} – {s_eo.index[-1]}\n"
            f"WB range: {s_wb.index[0]} – {s_wb.index[-1]}"
        )
    print(f"  Common period : {idx_common[0].date()} – {idx_common[-1].date()}  "
          f"({len(idx_common)} months)")
    return s_eo.loc[idx_common], s_wb.loc[idx_common]


# =============================================================================
# Figure A — monthly mean by year (EO vs WB on same axes)
# =============================================================================

def plot_compare_monthly_by_year(
    s_eo:     pd.Series,
    s_wb:     pd.Series,
    fig_path: Path,
    label:    str  = "",
    dpi:      int  = 150,
    show:     bool = False,
) -> None:
    """
    Spatial-mean monthly ETa (mm/day), one line per year.
    Solid lines = EO-EB, dashed lines = WB.  Matching colours per year.
    """
    years = sorted(set(s_eo.index.year) | set(s_wb.index.year))
    cmap_lines = plt.get_cmap("tab10", len(years))

    fig, ax = plt.subplots(figsize=(10, 5))
    months = np.arange(1, 13)
    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    for i, year in enumerate(years):
        color = cmap_lines(i)

        def _monthly_profile(s: pd.Series) -> pd.Series:
            sub = s[s.index.year == year]
            return sub.groupby(sub.index.month).mean()

        eo_m = _monthly_profile(s_eo)
        wb_m = _monthly_profile(s_wb)

        if not eo_m.empty:
            ax.plot(eo_m.index, eo_m.values,
                    marker="o", lw=2, color=color,
                    label=f"{year} EO-EB")
        if not wb_m.empty:
            ax.plot(wb_m.index, wb_m.values,
                    marker="s", lw=2, ls="--", color=color,
                    label=f"{year} WB")

    ax.set_xticks(months)
    ax.set_xticklabels(month_labels)
    ax.set_xlabel("Month")
    ax.set_ylabel("ETa (mm/day)")
    ax.set_title("Spatial-mean monthly ETa — EO-EB (solid) vs WB (dashed)",
                 fontweight="bold")
    ax.legend(ncol=2, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_monthly_mean_by_year{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure B — annual total bar chart (grouped EO vs WB)
# =============================================================================

def plot_compare_annual_total(
    s_eo:     pd.Series,
    s_wb:     pd.Series,
    fig_path: Path,
    label:    str  = "",
    dpi:      int  = 150,
    show:     bool = False,
) -> None:
    """Grouped annual bar chart: EO-EB (blue) and WB (orange)."""
    annual_eo = s_eo.resample("YE").sum()
    annual_wb = s_wb.resample("YE").sum()

    # Align to common years
    idx = annual_eo.index.intersection(annual_wb.index)
    annual_eo = annual_eo.loc[idx]
    annual_wb = annual_wb.loc[idx]
    years = annual_eo.index.year

    x   = np.arange(len(years))
    w   = 0.35
    fig, ax = plt.subplots(figsize=(max(8, len(years) * 1.1), 5))

    bars_eo = ax.bar(x - w / 2, annual_eo.values, w, label="EO-EB",
                     color=COLOR_EO, edgecolor="white")
    bars_wb = ax.bar(x + w / 2, annual_wb.values, w, label="WB",
                     color=COLOR_WB, edgecolor="white")

    for bars in (bars_eo, bars_wb):
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + max(annual_eo.max(), annual_wb.max()) * 0.01,
                    f"{h:.0f}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(years)
    ax.set_xlabel("Year")
    ax.set_ylabel("Annual ETa sum (mm/year)")
    ax.set_title("Annual domain-mean ETa — EO-EB vs Water-Balance model",
                 fontweight="bold")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_annual_total{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure C — monthly total bar chart (grouped EO vs WB)
# =============================================================================

def plot_compare_monthly_total(
    s_eo:     pd.Series,
    s_wb:     pd.Series,
    fig_path: Path,
    label:    str  = "",
    dpi:      int  = 150,
    show:     bool = False,
) -> None:
    """Grouped monthly bar chart for the full common record."""
    time_labels = s_eo.index.strftime("%Y-%m")
    x = np.arange(len(s_eo))
    w = 0.35

    fig, ax = plt.subplots(figsize=(max(12, len(s_eo) * 0.22), 5))
    ax.bar(x - w / 2, s_eo.values, w, label="EO-EB", color=COLOR_EO, edgecolor="none")
    ax.bar(x + w / 2, s_wb.values, w, label="WB",    color=COLOR_WB, edgecolor="none")

    # Year-boundary ticks only
    jan_pos    = [i for i, lbl in enumerate(time_labels) if lbl.endswith("-01")]
    jan_labels = [lbl[:4] for lbl in time_labels if lbl.endswith("-01")]
    ax.set_xticks(jan_pos)
    ax.set_xticklabels(jan_labels, rotation=0, ha="center", fontsize=9)
    for pos in jan_pos:
        ax.axvline(pos - 0.5 - w / 2, color="gray", lw=0.5, alpha=0.4, zorder=0)

    ax.set_ylabel("Monthly mean ETa (mm/day)")
    ax.set_title("Monthly domain-mean ETa — EO-EB vs Water-Balance model",
                 fontweight="bold")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_monthly_total{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure D — scatter EO-EB vs WB (coincident monthly means)
# =============================================================================

def plot_compare_scatter(
    s_eo:     pd.Series,
    s_wb:     pd.Series,
    fig_path: Path,
    label:    str  = "",
    dpi:      int  = 150,
    show:     bool = False,
) -> None:
    """Scatter of coincident monthly spatial-mean ETa, coloured by year."""
    years      = s_eo.index.year
    uniq_years = sorted(set(years))
    cmap_sc    = plt.get_cmap("tab10", len(uniq_years))
    year_color = {y: cmap_sc(i) for i, y in enumerate(uniq_years)}

    fig, ax = plt.subplots(figsize=(6, 6))
    for year in uniq_years:
        mask = years == year
        ax.scatter(s_eo[mask], s_wb[mask], color=year_color[year],
                   alpha=0.75, s=40, label=str(year))

    # 1:1 line
    all_vals = np.concatenate([s_eo.values, s_wb.values])
    lim = (0, np.nanmax(all_vals) * 1.05)
    ax.plot(lim, lim, "k--", lw=1, alpha=0.5, label="1:1")

    # Linear regression
    finite = np.isfinite(s_eo.values) & np.isfinite(s_wb.values)
    if finite.sum() > 2:
        m, b = np.polyfit(s_eo.values[finite], s_wb.values[finite], 1)
        corr = np.corrcoef(s_eo.values[finite], s_wb.values[finite])[0, 1]
        rmse = np.sqrt(np.mean((s_eo.values[finite] - s_wb.values[finite]) ** 2))
        xs   = np.array(lim)
        ax.plot(xs, m * xs + b, color="#555", lw=1.5,
                label=f"Reg.  R={corr:.2f}  RMSE={rmse:.2f}")
        ax.text(0.05, 0.93,
                f"R = {corr:.3f}\nRMSE = {rmse:.3f} mm/day\ny = {m:.2f}x + {b:.2f}",
                transform=ax.transAxes, fontsize=9,
                verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("EO-EB ETa (mm/day)")
    ax.set_ylabel("WB ETa (mm/day)")
    ax.set_title("Monthly spatial-mean ETa: EO-EB vs Water-Balance",
                 fontweight="bold")
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_aspect("equal")
    ax.grid(True, linestyle="--", alpha=0.35)
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    suffix = f"_{label}" if label else ""
    out = fig_path / f"compare_ET_scatter{suffix}.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


# =============================================================================
# Figure E — spatial monthly maps for a selected year (EO and WB, shared cscale)
# =============================================================================

def _spatial_maps_mosaic(
    ds:           xr.Dataset,
    var:          str,
    time_dim:     str,
    selected_year: int,
    vmin:         float,
    vmax:         float,
    cmap:         str,
    title_prefix: str,
    fig_path:     Path,
    fname:        str,
    dpi:          int  = 150,
    show:         bool = False,
) -> None:
    """Generic monthly spatial mosaic for one dataset."""
    da_year = ds[var].sel(
        {time_dim: ds[time_dim].dt.year == selected_year}
    )
    if da_year.sizes[time_dim] == 0:
        print(f"  {title_prefix}: no data for {selected_year} — skipped.")
        return

    months   = np.unique(da_year[time_dim].dt.month.values)
    n_months = len(months)
    ncols    = min(4, n_months)
    nrows    = int(np.ceil(n_months / ncols))

    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    sp_dims = [d for d in da_year.dims if d != time_dim]
    x_dim, y_dim = sp_dims[0], sp_dims[1]

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows),
                             squeeze=False)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap_obj = plt.get_cmap(cmap)

    for idx, month in enumerate(months):
        ax = axes[idx // ncols][idx % ncols]
        da_m = da_year.sel(
            {time_dim: da_year[time_dim].dt.month == month}
        ).mean(time_dim)

        im = ax.imshow(
            da_m.values,
            origin="upper",
            cmap=cmap_obj,
            norm=norm,
            extent=[float(da_m.x.min()), float(da_m.x.max()),
                    float(da_m.y.min()), float(da_m.y.max())],
            aspect="auto",
        )
        ax.set_title(month_names[month - 1], fontsize=9)
        ax.set_xlabel(x_dim, fontsize=7)
        ax.set_ylabel(y_dim, fontsize=7)
        ax.tick_params(labelsize=7)
        fig.colorbar(im, ax=ax, shrink=0.8).set_label("mm/day", fontsize=7)

    for idx in range(n_months, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(f"{title_prefix} — monthly ETa maps {selected_year}",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()

    fig_path.mkdir(parents=True, exist_ok=True)
    out = fig_path / fname
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    if show:
        plt.show()
    plt.close(fig)


def plot_compare_spatial_maps(
    ds_eo:         xr.Dataset,
    var_eo:        str,
    ds_wb:         xr.Dataset,
    selected_year: int,
    fig_path:      Path,
    label:         str  = "",
    dpi:           int  = 150,
    show:          bool = False,
) -> None:
    """
    Produce two spatial-map mosaics for *selected_year* (EO-EB and WB)
    using a shared colour scale so maps can be compared directly.
    """
    suffix = f"_{label}" if label else ""

    # Determine shared colour limits from both datasets for that year
    def _year_vals(ds: xr.Dataset, var: str, td: str) -> np.ndarray:
        mask = ds[td].dt.year == selected_year
        da   = ds[var].sel({td: mask})
        v    = da.values.ravel()
        return v[np.isfinite(v)]

    td_eo = "time"
    td_wb = "time"

    v_eo = _year_vals(ds_eo, var_eo, td_eo)
    v_wb = _year_vals(ds_wb, WB_ET_VAR, td_wb)

    if v_eo.size == 0 and v_wb.size == 0:
        print(f"  No data for {selected_year} in either dataset — spatial maps skipped.")
        return

    all_vals = np.concatenate([v for v in (v_eo, v_wb) if v.size > 0])
    vmin     = float(np.nanpercentile(all_vals, 2))
    vmax     = float(np.nanpercentile(all_vals, 98))

    _spatial_maps_mosaic(
        ds_eo, var_eo, td_eo, selected_year, vmin, vmax, CMAP_EO,
        "EO-EB ETa", fig_path,
        f"compare_ET_spatial_EO_{selected_year}{suffix}.png",
        dpi=dpi, show=show,
    )
    _spatial_maps_mosaic(
        ds_wb, WB_ET_VAR, td_wb, selected_year, vmin, vmax, CMAP_WB,
        "WB ETa", fig_path,
        f"compare_ET_spatial_WB_{selected_year}{suffix}.png",
        dpi=dpi, show=show,
    )


# =============================================================================
# Orchestrator
# =============================================================================

def run_comparison(
    ds_eo:         xr.Dataset,
    var_eo:        str,
    ds_wb:         xr.Dataset,
    label:         str,
    selected_year: int = 2020,
    dpi:           int = 150,
    show:          bool = False,
) -> None:
    """Run all five comparison figures."""
    fig_path.mkdir(parents=True, exist_ok=True)

    print("\n── Building common monthly series ───────────────────────────────")
    s_eo, s_wb = build_common_series(ds_eo, var_eo, ds_wb)

    print("\n── Fig A : monthly mean by year ─────────────────────────────────")
    plot_compare_monthly_by_year(s_eo, s_wb, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig B : annual total bar chart ───────────────────────────────")
    plot_compare_annual_total(s_eo, s_wb, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig C : monthly total bar chart ──────────────────────────────")
    plot_compare_monthly_total(s_eo, s_wb, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig D : scatter EO-EB vs WB ──────────────────────────────────")
    plot_compare_scatter(s_eo, s_wb, fig_path, label=label, dpi=dpi, show=show)

    print("\n── Fig E : spatial maps ─────────────────────────────────────────")
    plot_compare_spatial_maps(
        ds_eo, var_eo, ds_wb, selected_year,
        fig_path, label=label, dpi=dpi, show=show,
    )


# =============================================================================
# CLI argument parsing
# =============================================================================

def parse_args() -> tuple[list[str], int, int]:
    parser = argparse.ArgumentParser(
        description="Compare EO-EB vs Water-Balance ETa for Agramon"
    )
    parser.add_argument(
        "--plots", nargs="+", metavar="PLOT", default=["0"],
        help="Plots: row index, PlotID name, 'all', or mixed. Default: 0.",
    )
    parser.add_argument(
        "--year", type=int, default=2020,
        help="Year for spatial-map mosaics (default: 2020).",
    )
    parser.add_argument(
        "--scenario", type=int, default=3,
        help="sim_index to load from the simulation log (default: 0).",
    )

    known_flags = {"--plots", "--year", "--scenario"}
    clean: list[str] = []
    i = 1
    while i < len(sys.argv):
        tok = sys.argv[i]
        if tok in known_flags:
            clean.append(tok)
            i += 1
            while i < len(sys.argv) and not sys.argv[i].startswith("--"):
                clean.append(sys.argv[i])
                i += 1
        else:
            i += 1

    args = parser.parse_args(clean)
    return args.plots, args.year, args.scenario


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":
    selection, selected_year, scenario = parse_args()

    print("\n── Loading shapefiles ──────────────────────────────────────────")
    gdf_proj  = load_shapefiles()
    gdf_clip  = resolve_plot_selection(gdf_proj, selection)
    label     = selection_label(gdf_proj, selection)

    print("\n── Loading EO-EB ETa ───────────────────────────────────────────")
    ds_eo     = load_eo_et()
    ds_eo_c   = _clip_ds(ds_eo, gdf_clip)
    var_eo    = next(iter(ds_eo_c.data_vars))

    print(f"\n── Loading WB ETa (scenario {scenario}) ────────────────────────")
    ds_wb     = load_wb_et(scenario)
    ds_wb_c   = _clip_ds(ds_wb, gdf_clip)

    print(f"\n── Running comparison ──────────────────────────────────────────")
    print(f"  EO-EB var    : {var_eo}")
    print(f"  WB var       : {WB_ET_VAR}")
    print(f"  Plot label   : {label}")
    print(f"  Spatial year : {selected_year}")
    print(f"  Output dir   : {fig_path}\n")

    run_comparison(
        ds_eo_c, var_eo,
        ds_wb_c,
        label=label,
        selected_year=selected_year,
        show=False,
    )

    print("\nDone.")
