"""
lai_loop_utils.py - helpers of the CATHY "loop over LAI" model
================================================================
Pure helpers (no CATHY object needed) used by run_lai_loop_pycathy.py and compare_obs_pycathy.py.

Re-used from geoutils (nothing re-implemented when it already exists there):
    uCATHY.GridSpec / polygon_mask / cells_to_nodes / node_index_from_mesh / regrid_to_cells / eta_metrics
    uET.field_daily_series / norm_units / _grid_crs / _RATE_MM_DAY

New here: grid readers restricted to the model bounding box, pixel -> cell aggregation, monthly LAI cube,
LAI -> vegetation class / crop coefficient, water-table depth from psi, DEM sampling, observation loaders.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from geoutils import uCATHY as eo
from geoutils import uET

TIME_NAMES = ("time", "t", "valid_time", "date")
X_NAMES = ("x", "lon", "longitude")
Y_NAMES = ("y", "lat", "latitude")


# =============================================================================================
# 1. gridded NetCDF read through the model bounding box
# =============================================================================================
class GridSource:
    """Lazy reader of one variable of a (time, y, x) NetCDF, restricted to the model bounding box.

    bbox_model = (minx, miny, maxx, maxy) in ``model_crs``; it is transformed to the CRS of the file
    (``spatial_ref``, via uET._grid_crs) so only the needed pixels are ever read from disk
    (the *_s2-l2b.nc / *-ET-gf.nc files are 0.5-1.5 GB).
    """

    def __init__(self, path, var, bbox_model, model_crs, pad_px: int = 2):
        self.path, self.var = Path(path), var
        self.ds = xr.open_dataset(self.path)                       # lazy
        if var not in self.ds.data_vars:
            raise KeyError(f"'{var}' not in {self.path.name}; variables: {list(self.ds.data_vars)}")
        da = self.ds[var]
        self.tdim = next(d for d in TIME_NAMES if d in da.dims)
        self.xdim = next(d for d in X_NAMES if d in da.dims)
        self.ydim = next(d for d in Y_NAMES if d in da.dims)
        self.crs = uET._grid_crs(self.ds)                          # pyproj CRS or None
        x, y = self.ds[self.xdim].values.astype(float), self.ds[self.ydim].values.astype(float)
        minx, miny, maxx, maxy = bbox_model
        if self.crs is not None and not _same_crs(self.crs, model_crs):
            from pyproj import Transformer
            tr = Transformer.from_crs(model_crs, self.crs, always_xy=True)
            xs, ys = tr.transform([minx, maxx, minx, maxx], [miny, miny, maxy, maxy])
            minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
        resx, resy = abs(np.median(np.diff(x))), abs(np.median(np.diff(y)))
        ix = np.where((x >= minx - pad_px * resx) & (x <= maxx + pad_px * resx))[0]
        iy = np.where((y >= miny - pad_px * resy) & (y <= maxy + pad_px * resy))[0]
        if ix.size == 0 or iy.size == 0:
            raise ValueError(f"{self.path.name}: the model bounding box {bbox_model} (in {model_crs}) does not "
                             f"overlap the grid x [{x.min():.0f}, {x.max():.0f}] y [{y.min():.0f}, {y.max():.0f}] "
                             f"(grid CRS: {self.crs.to_string() if self.crs else 'unknown'})")
        self.xi, self.yi = slice(int(ix.min()), int(ix.max()) + 1), slice(int(iy.min()), int(iy.max()) + 1)
        self.x, self.y = x[self.xi], y[self.yi]
        self.times = pd.DatetimeIndex(self.ds[self.tdim].values)
        raw = str(da.attrs.get("units", ""))
        self.units = uET.norm_units(raw) if raw else ""
        self.to_mm_day = uET._RATE_MM_DAY.get(self.units, 1.0)    # only meaningful for ET / rain rates

    def read(self, t0, t1):
        """(times, array (T, ny, nx)) of the steps with t0 <= time < t1 (empty arrays if none)."""
        sel = np.where((self.times >= pd.Timestamp(t0)) & (self.times < pd.Timestamp(t1)))[0]
        if sel.size == 0:
            return self.times[:0], np.empty((0, self.y.size, self.x.size))
        a = self.ds[self.var].isel({self.tdim: slice(int(sel.min()), int(sel.max()) + 1),
                                    self.ydim: self.yi, self.xdim: self.xi})
        a = a.transpose(self.tdim, self.ydim, self.xdim).values.astype(float)
        return self.times[sel.min(): sel.max() + 1], a

    def close(self):
        self.ds.close()


def _same_crs(a, b) -> bool:
    from pyproj import CRS
    try:
        return CRS.from_user_input(a) == CRS.from_user_input(b)
    except Exception:       # noqa: BLE001
        return str(a) == str(b)


# =============================================================================================
# 2. pixels -> CATHY cells
# =============================================================================================
class CellAggregator:
    """Mean of the EO pixels falling inside each CATHY surface cell (M, N). Cells without pixel are NaN
    (see fill_nearest). Works for (ny, nx) and (T, ny, nx) arrays. Handles a CRS change pixel -> model."""

    def __init__(self, x, y, src_crs, grid: "eo.GridSpec"):
        X, Y = np.meshgrid(x, y)
        if src_crs is not None and grid.crs and not _same_crs(src_crs, grid.crs):
            from pyproj import Transformer
            X, Y = Transformer.from_crs(src_crs, grid.crs, always_xy=True).transform(X, Y)
        col = np.floor((X - grid.xll) / grid.dx).astype(int)
        top = grid.yll + grid.M * grid.dy
        row = np.floor((top - Y) / grid.dy).astype(int) if grid.north_up else np.floor((Y - grid.yll) / grid.dy).astype(int)
        ok = (col >= 0) & (col < grid.N) & (row >= 0) & (row < grid.M)
        self.M, self.N = grid.M, grid.N
        self.ok = ok.ravel()
        self.flat = (row * grid.N + col).ravel()[self.ok]
        if not self.ok.any():
            raise ValueError("no EO pixel falls inside the CATHY grid (check CRS / bounding box)")

    def _one(self, a2d):
        v = np.asarray(a2d, float).ravel()[self.ok]
        fin = np.isfinite(v)
        s = np.bincount(self.flat[fin], weights=v[fin], minlength=self.M * self.N)
        n = np.bincount(self.flat[fin], minlength=self.M * self.N)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(n > 0, s / n, np.nan).reshape(self.M, self.N)

    def __call__(self, arr):
        arr = np.asarray(arr, float)
        return self._one(arr) if arr.ndim == 2 else np.stack([self._one(a) for a in arr])


def fill_nearest(a2d: np.ndarray, fallback: float = 0.0) -> np.ndarray:
    """NaN cells of a 2-D field <- nearest valid cell (all NaN -> fallback)."""
    a = np.array(a2d, float)
    bad = ~np.isfinite(a)
    if not bad.any():
        return a
    if bad.all():
        return np.full_like(a, fallback)
    from scipy.ndimage import distance_transform_edt
    idx = distance_transform_edt(bad, return_distances=False, return_indices=True)
    return a[tuple(idx)]


def cells_to_node_flat(cells: np.ndarray, method: str = "mean") -> np.ndarray:
    """(T, M, N) cell fields -> (T, (M+1)*(N+1)) node fields, flattened like uCATHY.node_index_from_mesh expects."""
    cells = np.asarray(cells, float)
    if cells.ndim == 2:
        return eo.cells_to_nodes(cells, method).ravel()
    return np.stack([eo.cells_to_nodes(c, method).ravel() for c in cells])


# =============================================================================================
# 3. LAI cube (month x M x N) and its derivatives
# =============================================================================================
def month_starts(start, end) -> pd.DatetimeIndex:
    """first day of every month from the month of ``start`` to the month of ``end`` (inclusive)."""
    return pd.date_range(pd.Timestamp(start).to_period("M").to_timestamp(),
                         pd.Timestamp(end).to_period("M").to_timestamp(), freq="MS")


def build_lai_cube(src: GridSource, agg: CellAggregator, months: pd.DatetimeIndex, min_obs: int = 1,
                   log=print, stats: dict | None = None):
    """Returns (cube, n_obs): monthly LAI on the CATHY cells (len(months), M, N) and the number of acquisitions per month.

    Per month: median over the acquisitions of the month (pixel-wise, NaN = cloud), then mean of the pixels of
    each cell. Gaps: NaN cells <- previous month (carried forward), then next valid month (back-fill, start of the
    series), then nearest valid cell, then 0. Months with fewer than ``min_obs`` acquisitions are carried forward.
    Negative LAI is clipped to 0.

    ``stats`` (optional dict, filled in place) keeps what the input checks need: n_raw (acquisitions per month, valid or
    not), cloud_frac (fraction of NaN pixels per month), acq_times / acq_mean (bounding-box mean of every valid acquisition).
    """
    M, N = agg.M, agg.N
    cube = np.full((len(months), M, N), np.nan)
    n_obs = np.zeros(len(months), int)
    n_raw, cloud = np.zeros(len(months), int), np.full(len(months), np.nan)
    acq_t, acq_m = [], []
    for i, m in enumerate(months):
        times, a = src.read(m, m + pd.offsets.MonthBegin(1))
        n_raw[i] = len(times)
        if len(times) == 0:
            continue
        valid = np.isfinite(a).any(axis=(1, 2))
        n_obs[i] = int(valid.sum())
        cloud[i] = float(np.isnan(a).mean())
        if valid.any():
            acq_t += list(times[valid])
            acq_m += list(np.nanmean(a[valid], axis=(1, 2)))
        if n_obs[i] < min_obs:
            continue
        with np.errstate(all="ignore"):
            med = np.nanmedian(a[valid], axis=0)
        cube[i] = agg(med)
    log(f"LAI cube : {len(months)} months, {int((n_obs > 0).sum())} with >= {min_obs} acquisition(s) "
        f"({src.path.name}:{src.var})")
    for i in range(len(months)):                               # carry forward
        if i > 0:
            gap = ~np.isfinite(cube[i])
            cube[i][gap] = cube[i - 1][gap]
    for i in range(len(months) - 2, -1, -1):                   # back-fill the start of the series
        gap = ~np.isfinite(cube[i])
        cube[i][gap] = cube[i + 1][gap]
    cube = np.stack([fill_nearest(c, 0.0) for c in cube])
    if stats is not None:
        stats.update(n_raw=n_raw, cloud_frac=cloud, acq_times=pd.DatetimeIndex(acq_t), acq_mean=np.asarray(acq_m, float))
    return np.clip(cube, 0.0, None), n_obs


def lai_climatology(cube: np.ndarray, months: pd.DatetimeIndex, ref_years) -> np.ndarray:
    """Counterfactual LAI WITHOUT the LAI drop: for every calendar month the median over ``ref_years`` of the cube.
    Returns a cube of the same shape as ``cube`` (every year identical)."""
    ref = [y for y in ref_years if (months.year == y).any()]
    if not ref:
        raise ValueError(f"--ref-years {list(ref_years)} are outside the LAI months {months[0]:%Y-%m} .. {months[-1]:%Y-%m}")
    out = np.empty_like(cube)
    for mo in range(1, 13):
        sel = np.where((months.month == mo) & np.isin(months.year, ref))[0]
        tgt = np.where(months.month == mo)[0]
        if sel.size == 0:                                       # month not covered by the reference years
            out[tgt] = np.nanmedian(cube[months.month == mo], axis=0)
        else:
            out[tgt] = np.median(cube[sel], axis=0)
    return out


def lai_to_veg_map(lai_2d: np.ndarray, thresholds) -> np.ndarray:
    """Integer vegetation class map 1..K (K = len(thresholds)+1): class i+1 when thresholds[i-1] < LAI <= thresholds[i].
    The number of classes is FIXED for the whole run (the SOIL vegetation table needs a constant row count)."""
    veg = np.ones(lai_2d.shape, dtype=int)
    for i, th in enumerate(thresholds, start=2):
        veg[lai_2d > th] = i
    return veg


def kc_from_lai(lai: np.ndarray, kc_min: float, kc_max: float, k_ext: float = 0.5) -> np.ndarray:
    """Crop coefficient from LAI (Beer-Lambert cover fraction): Kc = kc_min + (kc_max-kc_min) * (1 - exp(-k*LAI))."""
    return kc_min + (kc_max - kc_min) * (1.0 - np.exp(-k_ext * np.clip(lai, 0.0, None)))


# =============================================================================================
# 4. CATHY outputs
# =============================================================================================
def water_table_depth(psi: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Water-table depth below the surface (m, >= 0) per column and output time.

    psi : (nt, nstr, nnod) pressure head (m); layer 0 = surface, last layer = mesh bottom.
    z   : (nstr, nnod) node elevations.
    The table is the elevation where psi = 0, searched from the bottom up: the deepest unsaturated node (psi < 0)
    and the saturated node below it are linearly interpolated. Columns saturated everywhere -> 0 (water at the surface),
    columns unsaturated down to the bottom -> depth of the mesh base (the table is below the mesh: lower bound).
    """
    nt, nstr, nnod = psi.shape
    unsat = psi < 0
    any_unsat = unsat.any(axis=1)                                                  # (nt, nnod)
    k_u = nstr - 1 - np.argmax(unsat[:, ::-1, :], axis=1)                          # deepest unsaturated layer
    k_s = np.minimum(k_u + 1, nstr - 1)
    ii = k_u[:, None, :]
    jj = k_s[:, None, :]
    p_u = np.take_along_axis(psi, ii, axis=1)[:, 0, :]
    p_s = np.take_along_axis(psi, jj, axis=1)[:, 0, :]
    zb = np.broadcast_to(z, psi.shape)
    z_u = np.take_along_axis(zb, ii, axis=1)[:, 0, :]
    z_s = np.take_along_axis(zb, jj, axis=1)[:, 0, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(p_s != p_u, (0.0 - p_u) / (p_s - p_u), 0.0)
    z_wt = z_u + np.clip(frac, 0.0, 1.0) * (z_s - z_u)
    depth = z[0][None, :] - z_wt
    depth = np.where(any_unsat, depth, 0.0)                                        # fully saturated column
    return np.clip(depth, 0.0, None)


def to_datetime(index, base) -> pd.DatetimeIndex:
    """CATHY output times (seconds, or Timedelta) -> datetimes from ``base``."""
    idx = pd.Index(index)
    if pd.api.types.is_timedelta64_dtype(idx):
        return pd.DatetimeIndex(pd.Timestamp(base) + idx)
    return pd.DatetimeIndex(pd.Timestamp(base) + pd.to_timedelta(np.asarray(idx, float), unit="s"))


# =============================================================================================
# 5. DEM
# =============================================================================================
def build_dem(bounds, crs, dx, margin, slope=0.01, dem_path=None):
    """DEM grid (row 0 = north) covering ``bounds`` (minx, miny, maxx, maxy in ``crs``) + margin.
    Returns dem (M, N), x_ll, y_ll. Without dem_path: gently sloping plane (descending to the south);
    with a GeoTIFF: bilinear sampling, reprojected on the fly if it is not in ``crs``."""
    minx, miny, maxx, maxy = bounds
    x_ll = np.floor((minx - margin) / dx) * dx
    y_ll = np.floor((miny - margin) / dx) * dx
    N = int(np.ceil((maxx + margin - x_ll) / dx))
    M = int(np.ceil((maxy + margin - y_ll) / dx))
    if dem_path is None:
        rows = np.arange(M)[:, None]
        return 100.0 + slope * dx * (M - 1 - rows) * np.ones((1, N)), float(x_ll), float(y_ll)
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.vrt import WarpedVRT
    from rasterio.windows import from_bounds
    with rasterio.open(dem_path) as src:
        same = src.crs is not None and _same_crs(src.crs, crs)
        ds = src if same or src.crs is None else WarpedVRT(src, crs=crs, resampling=Resampling.bilinear)
        win = from_bounds(x_ll, y_ll, x_ll + N * dx, y_ll + M * dx, ds.transform)
        dem = ds.read(1, window=win, out_shape=(M, N), boundless=True, resampling=Resampling.bilinear,
                      fill_value=np.nan).astype(float)
    if np.isnan(dem).all():
        raise ValueError("--dem does not overlap the model bounding box")
    return np.where(np.isnan(dem), np.nanmin(dem), dem), float(x_ll), float(y_ll)


def add_fake_outlet(dem, corner="sw", drop=1e-3):
    """Lower one corner cell by ``drop`` m so CATHY's preprocessor finds a unique outlet (row 0 = north)."""
    if corner == "none":
        return dem
    dem = dem.copy()
    dem[-1 if corner[0] == "s" else 0, 0 if corner[1] == "w" else -1] -= drop
    return dem


# =============================================================================================
# 6. observations (Lanzhot) - used by compare_obs_pycathy.py
# =============================================================================================
def read_obs_eta(path) -> pd.Series:
    """Eddy-covariance ETa, daily, mm/day. File: 'Date;LNZ' with dd.mm.yyyy dates (empty trailing lines tolerated)."""
    df = pd.read_csv(path, sep=";", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    d = pd.to_datetime(df.iloc[:, 0], format="%d.%m.%Y", errors="coerce")
    v = pd.to_numeric(df.iloc[:, 1].str.replace(",", ".", regex=False), errors="coerce")
    s = pd.Series(v.values, index=d).dropna()
    s = s[~s.index.isna()]
    s.index = pd.DatetimeIndex(s.index).normalize()
    s.name = "ETa_EC_mm_d"
    return s[~s.index.duplicated()].sort_index()


def read_obs_gwl(path) -> pd.DataFrame:
    """Piezometers: columns Date;Time;GWL1 (m);...;GWL4 (m), 30-min, dd.mm.yyyy. Values in m relative to the ground
    (negative = water table below the surface). NOTE: the file named *.xlsx of the project is a semicolon CSV in
    disguise; a real Excel file is read with pandas.read_excel."""
    path = Path(path)
    head = path.read_bytes()[:4]
    if head[:2] == b"PK":                                       # genuine xlsx (zip)
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep=";", decimal=".")
    df.columns = [str(c).strip() for c in df.columns]
    dt = pd.to_datetime(df["Date"].astype(str) + " " + df["Time"].astype(str), format="%d.%m.%Y %H:%M:%S", errors="coerce")
    out = df.drop(columns=["Date", "Time"]).apply(pd.to_numeric, errors="coerce")
    out.index = dt
    out = out[~out.index.isna()]
    out.columns = [c.split(" ")[0] for c in out.columns]       # 'GWL1 (m)' -> 'GWL1'
    return out.sort_index()


def metrics(model: pd.Series, obs: pd.Series, demean: bool = False) -> dict:
    """bias (model-obs), RMSE, Pearson r, NSE, KGE on the common dates. demean=True removes each series' mean first
    (compares the dynamics only; use it when the absolute datum of the model and the piezometers differ)."""
    c = model.dropna().index.intersection(obs.dropna().index)
    if len(c) < 3:
        return dict(n=len(c), bias=np.nan, rmse=np.nan, r=np.nan, nse=np.nan, kge=np.nan)
    m, o = model.loc[c].to_numpy(float), obs.loc[c].to_numpy(float)
    if demean:
        m, o = m - m.mean(), o - o.mean()
    rmse = float(np.sqrt(np.mean((m - o) ** 2)))
    r = float(np.corrcoef(m, o)[0, 1]) if m.std() > 0 and o.std() > 0 else np.nan
    nse = float(1 - np.sum((m - o) ** 2) / np.sum((o - o.mean()) ** 2)) if np.sum((o - o.mean()) ** 2) > 0 else np.nan
    kge = np.nan
    if np.isfinite(r) and o.mean() != 0 and o.std() > 0:
        kge = float(1 - np.sqrt((r - 1) ** 2 + (m.std() / o.std() - 1) ** 2 + (m.mean() / o.mean() - 1) ** 2))
    return dict(n=len(c), bias=float(m.mean() - o.mean()), rmse=rmse, r=r, nse=nse, kge=kge)
