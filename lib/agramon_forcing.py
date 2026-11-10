"""
agramon_forcing.py — Input loaders for the Agramon pipeline
===========================================================
GRwater project — Agramon study site

Loaders for the atmospheric forcing (satellite ETp + rain rasters, the
Agramon ET CSV time series, and the legacy ERA5 series) and for the
monthly Earth-Observation LAI. Moved verbatim out of
Agramon_withLAI_withETp.py (REFACTORING.md, step 2.1).

NOTE: ``load_era5`` is legacy (kept for quick tests only) and is
deprecated; see REFACTORING.md, step 3.3.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import xarray as xr

import Agramon_utils as AgUtils
from agramon_config import (
    EO_PATH,
    ETP_TS_COLUMN_DEFAULT,
    NATIVE_CRS,
    TARGET_CRS,
    etp_path,
    rain_path,
)


def load_monthly_lai(lai_input_dir: Path) -> xr.Dataset:
    """
    Load monthly LAI dataset, assign its native CRS (EPSG:32630),
    and reproject to TARGET_CRS.
    """
    ds_lai = xr.open_dataset(lai_input_dir / "lai_monthly.nc")

    if ds_lai.rio.crs is None:
        ds_lai = ds_lai.rio.write_crs("EPSG:32630")
        print("  LAI CRS assigned : EPSG:32630")
    print(f"  LAI native CRS   : {ds_lai.rio.crs}")

    if str(ds_lai.rio.crs) != TARGET_CRS:
        print(f"  Reprojecting LAI {ds_lai.rio.crs} → {TARGET_CRS}")
        ds_lai = ds_lai.rio.reproject(TARGET_CRS)
    print(f"  LAI final CRS    : {ds_lai.rio.crs}")
    print(f"  LAI x range      : {float(ds_lai.x.min()):.1f} – {float(ds_lai.x.max()):.1f}")
    print(f"  LAI y range      : {float(ds_lai.y.min()):.1f} – {float(ds_lai.y.max()):.1f}")
    return ds_lai


def load_era5(start_year: int, end_year: int) -> xr.Dataset:
    """Load ERA5 climate forcing."""
    ds = AgUtils.load_era5_series(root_path=EO_PATH, start_year=start_year, end_year=end_year)
    if ds is not None:
        print(f"  Loaded ERA5 : {start_year} – {end_year - 1}")
    else:
        print("  No ERA5 datasets loaded — check EO_PATH.")
    return ds


def _finalise_ds(ds: xr.Dataset, label: str) -> xr.Dataset:
    """Write CRS + spatial dims after load or reproject."""
    ds = ds.rio.write_crs(TARGET_CRS, inplace=False)
    ds = ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    print(f"  {label} final CRS : {ds.rio.crs}")
    print(f"  {label} x range   : {float(ds.x.min()):.1f} – {float(ds.x.max()):.1f}")
    print(f"  {label} y range   : {float(ds.y.min()):.1f} – {float(ds.y.max()):.1f}")
    return ds


def _load_nc(path: Path, label: str) -> xr.Dataset:
    """Open NetCDF, assign native CRS (EPSG:32630), reproject to TARGET_CRS."""
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

    return _finalise_ds(ds, label)


def load_spatial_etp(path: Path = etp_path) -> xr.Dataset:
    """Load ETp dataset and reproject to TARGET_CRS."""
    return _load_nc(path, "ETp")


VALUE_COLS = ["ET", "T", "ET_0_dry", "T_0_dry", "ET_0_wet", "T_0_wet"]


def load_agramon_et(filename: Path = "agramon_ET.csv") -> pd.DataFrame:
    """Load the single-series Agramon ET file.

    Returns
    -------
    pd.DataFrame
        Indexed by a daily DatetimeIndex named ``t``, sorted ascending.
        Columns: year, month, ET, T, ET_0_dry, T_0_dry, ET_0_wet,
        T_0_wet, feature_id, ET/ET_dry, ET/ET_wet.
    """
    df = pd.read_csv(filename, parse_dates=["t"])
    df = df.sort_values("t").set_index("t")
    df.index.name = "t"
    return df


def load_treatments_et(
    filename: Path = "agramon_treatments_ET.csv",
    drop_empty: bool = True,
    add_ratios: bool = True,
) -> pd.DataFrame:
    """Load the per-treatment Agramon ET file (long format).

    Parameters
    ----------
    filename : str or Path
        Path to ``agramon_treatments_ET.csv``.
    drop_empty : bool, default True
        Drop treatments (feature_id) whose ET is entirely NaN
        (feature 9 in the current file).
    add_ratios : bool, default True
        Add ``ET/ET_dry`` and ``ET/ET_wet`` columns, as in
        ``agramon_ET.csv`` (NaN where the denominator is 0).

    Returns
    -------
    pd.DataFrame
        MultiIndex (t, feature_id), sorted. The constant ``band`` and
        ``spatial_ref`` columns are dropped. Columns: ET, T, ET_0_dry,
        T_0_dry, ET_0_wet, T_0_wet (+ ratios), year, month.
    """
    df = pd.read_csv(filename, parse_dates=["t"])
    df = df.drop(columns=["band", "spatial_ref"], errors="ignore")

    if drop_empty:
        has_et = df.groupby("feature_id")["ET"].transform(lambda s: s.notna().any())
        df = df[has_et]

    if add_ratios:
        df["ET/ET_dry"] = df["ET"] / df["ET_0_dry"].where(df["ET_0_dry"] != 0)
        df["ET/ET_wet"] = df["ET"] / df["ET_0_wet"].where(df["ET_0_wet"] != 0)

    df["year"] = df["t"].dt.year
    df["month"] = df["t"].dt.month

    return df.set_index(["t", "feature_id"]).sort_index()


def to_wide(df: pd.DataFrame, var: str = "ET") -> pd.DataFrame:
    """Pivot the long treatments frame to index=date, columns=feature_id."""
    return df[var].unstack("feature_id")


def monthly_mean(df: pd.DataFrame, cols: Optional[list] = None) -> pd.DataFrame:
    """Monthly mean per feature_id (treatments frame) or overall (single series)."""
    cols = cols or [c for c in VALUE_COLS if c in df.columns]
    if isinstance(df.index, pd.MultiIndex):
        t = df.index.get_level_values("t")
        fid = df.index.get_level_values("feature_id")
        return df[cols].groupby([fid, t.to_period("M")]).mean()
    return df[cols].groupby(df.index.to_period("M")).mean()


def load_all(folder: Path = Path("."), **treatment_kwargs) -> pd.DataFrame:
    """Load both files and stack them into a single flat DataFrame.

    Columns: t, source, feature_id, year, month, ET, T, ET_0_dry, T_0_dry,
    ET_0_wet, T_0_wet, ET/ET_dry, ET/ET_wet.

    ``source`` is ``"agramon_ET"`` for the single series and
    ``"treatments"`` for the per-treatment rows. Both use feature_id 0 for
    their first series, so always filter on ``source`` + ``feature_id``
    together (they are not the same series).
    """
    folder = Path(folder)
    single = load_agramon_et(folder / "agramon_ET.csv").reset_index()
    single["source"] = "agramon_ET"
    treat = load_treatments_et(
        folder / "agramon_treatments_ET.csv", **treatment_kwargs
    ).reset_index()
    treat["source"] = "treatments"

    out = pd.concat([single, treat], ignore_index=True)
    first = ["t", "source", "feature_id", "year", "month"]
    out = out[first + [c for c in out.columns if c not in first]]
    return out.sort_values(["source", "feature_id", "t"]).reset_index(drop=True)


def load_rain(path: Path = rain_path) -> xr.Dataset:
    """Load rain dataset and reproject to TARGET_CRS."""
    return _load_nc(path, "rain")


def build_etp_ds_from_timeseries(
    df: pd.DataFrame,
    ds_rain: xr.Dataset,
    plot_id: int,
    column: str = ETP_TS_COLUMN_DEFAULT,
    feature_id: int | None = None,
    source: str = "treatments",
    rain_var: str = "TP-DD",
    etp_var: str = "ET_0-gf",
) -> xr.Dataset:
    """
    Build an ETp dataset with the same (time, y, x) shape, coordinates and
    CRS as ``ds_rain`` from a single daily time series (``load_all`` output).

    The series is identical for every pixel, so the returned dataset is a
    drop-in replacement for ``load_spatial_etp()``: same variable name
    (``etp_var``), same time axis as the rain, units mm/day.

    Parameters
    ----------
    df : pd.DataFrame
        Output of ``load_all`` (columns t, source, feature_id, <column>, ...).
    ds_rain : xr.Dataset
        Rain dataset (already reprojected); provides time axis, grid, CRS.
    plot_id : int
        ``--dem-plot`` value. Used to pick the treatment when ``feature_id``
        is None (feature_id = plot_id - 1).
    column : str
        Column of ``df`` used as ETp (mm/day).
    feature_id : int, optional
        Explicit treatment feature_id (overrides the plot_id rule).
    source : str
        ``"treatments"`` or ``"agramon_ET"``.
    """
    fid = int(feature_id) if feature_id is not None else int(plot_id) - 1

    sel = df[(df["source"] == source) & (df["feature_id"] == fid)]
    if sel.empty:
        avail = sorted(int(v) for v in df.loc[df["source"] == source, "feature_id"].unique())
        raise ValueError(
            f"No rows for source={source!r}, feature_id={fid} "
            f"(plot_id={plot_id}). Available feature_ids: {avail}. "
            "Pass --etp-ts-feature-id to choose one explicitly."
        )
    if column not in sel.columns:
        raise ValueError(f"Column {column!r} not in df (columns: {list(sel.columns)})")

    series = sel.set_index("t")[column]
    if series.index.has_duplicates:
        raise ValueError(f"Duplicate dates for feature_id={fid} in column {column!r}")

    rain_times = pd.DatetimeIndex(ds_rain["time"].values)
    vals = series.reindex(rain_times.normalize()).to_numpy(dtype=float)
    n_missing = int(np.isnan(vals).sum())
    if n_missing:
        bad = rain_times[np.isnan(vals)]
        raise ValueError(
            f"ETp series (feature_id={fid}, column={column!r}) has {n_missing}/"
            f"{len(vals)} missing/NaN values over the rain time axis "
            f"({bad[0]:%Y-%m-%d} … {bad[-1]:%Y-%m-%d}). Check feature_id / column."
        )

    template = ds_rain[rain_var]
    ts = xr.DataArray(vals, dims="time", coords={"time": ds_rain["time"].values})
    _, ts_b = xr.broadcast(template, ts)
    etp_da = ts_b.transpose(*template.dims).astype("float64")
    # xr.broadcast keeps only the coords of `ts`/dims — restore the rest of the
    # template's coords (notably the scalar `spatial_ref` carrying the CRS, so
    # ds_etp.rio.crs / map_grid_to_mesh behave exactly as for the real ETp raster).
    etp_da = etp_da.assign_coords(
        {k: v for k, v in template.coords.items() if k not in etp_da.coords}
    )
    etp_da.attrs = {"units": "mm/day",
                    "long_name": f"ETp time series ({source}, feature_id={fid}, {column})"}

    ds_etp = etp_da.to_dataset(name=etp_var)
    print(f"  ETp time series : source={source}  feature_id={fid} "
          f"(plot {plot_id})  column={column}  "
          f"mean={vals.mean():.2f} mm/day  range=[{vals.min():.2f}, {vals.max():.2f}]")
    return ds_etp
