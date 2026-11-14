"""
agramon_mesh.py - DEM loading, mesh set-up and boundary conditions for the Agramon pipeline
=========================================================================================
Moved unchanged from Agramon_withLAI_withETp.py (REFACTORING.md 2.3); the main script
re-exports the names, so withLAI_results.py and the synthetic Ks study keep working.

Planned (TO-DO.md 5.4): systematic soil / root properties in the mesh .vtk will be built on
top of this module.
"""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rioxarray as rxr  # noqa: F401  (activates .rio accessor)
from scipy.ndimage import gaussian_filter
from pyCATHY import cathy_tools

from agramon_config import TARGET_CRS


def load_dem_from_tif(
    tif_path: Path, target_crs: str = TARGET_CRS,
    gdf_clip: gpd.GeoDataFrame | None = None,
    fid: int | str | None = None,
    plot_id: str | None = None,
    gbuffer: float | None = None,
    show: bool = False,
    resample_resolution: float | None = None,
):
    """
    Load a DEM directly from a GeoTIFF, independent of ``AgUtils.load_dem``
    (which expects a folder of .adf grids). Mirrors that function's
    conventions exactly so it's a drop-in swap at the call site:

      - reprojects to ``target_crs`` (assumes EPSG:25830 if the GeoTIFF
        has no CRS embedded, same fallback as ``AgUtils.load_dem``);
      - by default, does NOT clip to any shapefile — it only crops to
        the DEM's own valid-data bounding box (trims all-NaN rows/cols
        at the edges), exactly like ``AgUtils.load_dem`` does (its
        shapefile-clip code path is present but disabled/dead there
        too). Pass ``gdf_clip`` to additionally clip to a shape's
        bounding box + buffer if the new DEM's extent turns out to be
        much larger than the catchment of interest;
      - fills nodata with ``-9999`` (CATHY's expected DEM no-data
        sentinel), NOT NaN — matching ``AgUtils.load_dem`` bit for bit,
        since the pyCATHY preprocessor/mesh writer expects this value;
      - sets ``xllcorner``/``yllcorner`` to the raw min(x)/min(y) pixel
        centre coordinates (no half-pixel edge offset), same as
        ``AgUtils.load_dem``.

    Parameters
    ----------
    tif_path : Path
        Full path to the GeoTIFF DEM.
    target_crs : str
        CRS every spatial dataset in the pipeline is aligned to.
    gdf_clip : gpd.GeoDataFrame, optional
        GeoDataFrame whose bounding box is used to clip the DEM.
        If provided, the DEM is clipped to this bounding box + buffer.
    fid : int or str, optional
        Feature ID to select a single polygon from ``gdf_clip`` by its
        ``fid_1`` column. Mutually exclusive with ``plot_id``.
    plot_id : str, optional
        Plot ID to select polygons from ``gdf_clip`` by its ``PlotID``
        column (e.g., ``args.dem_plot``). Mutually exclusive with ``fid``.
    gbuffer : float, optional
        Padding in metres around ``gdf_clip``'s bbox. Defaults to 1
        pixel of the current raster resolution.
    show : bool
        If True, plot the masked DEM for a quick visual sanity check.
    resample_resolution : float, optional
        Target resolution in metres for resampling (assumes square
        pixels). If provided, the DEM is resampled to this resolution
        after CRS transformation but before cropping to valid data.

    Returns
    -------
    raster_DEM : xr.DataArray
        DEM reprojected/cropped to its own valid-data extent, with
        ``.rio.crs`` set. Same x/y grid as ``raster_DEM_masked``.
    raster_DEM_masked : np.ndarray
        DEM elevations, with NaN → -9999 (CATHY no-data convention).
    xllcorner, yllcorner : float
        Lower-left corner of the grid (pixel-centre coordinate, not
        cell edge — matches ``AgUtils.load_dem``).
    res_x, res_y : float
        Pixel size (positive values).
    """
    tif_path = Path(tif_path)
    if not tif_path.is_file():
        raise FileNotFoundError(f"DEM GeoTIFF not found: {tif_path}")

    print(f"  DEM source (tif) : {tif_path}")
    raster_DEM = rxr.open_rasterio(tif_path, masked=True).isel(band=0)

    if raster_DEM.rio.crs is None:
        raster_DEM = raster_DEM.rio.write_crs("EPSG:25830")
        print("  DEM CRS was missing — assumed EPSG:25830")

    # Handle CRS transformation and resampling together
    needs_reproject = str(raster_DEM.rio.crs) != target_crs
    needs_resample = resample_resolution is not None

    if needs_reproject or needs_resample:
        dst_crs = target_crs if needs_reproject else raster_DEM.rio.crs
        kwargs = {"resolution": resample_resolution} if needs_resample else {}
        print(f"  Transforming DEM to {dst_crs} at {resample_resolution} m" if needs_resample
              else f"  Transforming DEM to {dst_crs}")
        raster_DEM = raster_DEM.rio.reproject(dst_crs, **kwargs)

    print(f"  DEM CRS : {raster_DEM.rio.crs}")

    res_x = abs(float(raster_DEM.rio.resolution()[0]))
    res_y = abs(float(raster_DEM.rio.resolution()[1]))
    print(f"  DEM resolution : {res_x:.2f} × {res_y:.2f} m")

    # Select and clip to gdf_clip bounding box + buffer if provided
    if gdf_clip is not None:
        # Select specific feature(s) from gdf_clip if fid or plot_id is provided
        if fid is not None and plot_id is not None:
            raise ValueError("Cannot specify both 'fid' and 'plot_id'. Use one or the other.")
        if fid is not None:
            gdf_clip = gdf_clip[gdf_clip["fid_1"] == fid].copy()
            print(f"  Selected fid={fid} from gdf_clip ({len(gdf_clip)} features)")
        elif plot_id is not None:
            gdf_clip = gdf_clip[gdf_clip["PlotID"] == plot_id].copy()
            print(f"  Selected PlotID={plot_id} from gdf_clip ({len(gdf_clip)} features)")

        if len(gdf_clip) == 0:
            raise ValueError(f"No features found in gdf_clip for {'fid' if fid is not None else 'PlotID'}={fid if fid is not None else plot_id}")

        if str(gdf_clip.crs) != target_crs:
            gdf_clip = gdf_clip.to_crs(target_crs)
            print(f"  Reprojected gdf_clip to {target_crs}")

        if gbuffer is None:
            gbuffer = res_x  # Default: 1 pixel of current resolution
        minx, miny, maxx, maxy = gdf_clip.total_bounds
        minx -= gbuffer
        miny -= gbuffer
        maxx += gbuffer
        maxy += gbuffer
        print(f"  Clipping DEM to gdf_clip bbox + {gbuffer:.2f} m buffer")
        #raster_DEM = raster_DEM.rio.clip_box(minx, miny, maxx, maxy)
        raster_DEM = raster_DEM.rio.clip(gdf_clip.geometry, all_touched=True)

    # Crop to valid-data bounding box (trim all-NaN rows/cols)
    valid_mask = ~np.isnan(raster_DEM)
    valid_x = raster_DEM["x"].where(valid_mask.any(dim="y"), drop=True)
    valid_y = raster_DEM["y"].where(valid_mask.any(dim="x"), drop=True)
    raster_DEM = (
        raster_DEM
        .where((raster_DEM["x"] >= float(valid_x.min())) &
               (raster_DEM["x"] <= float(valid_x.max())), drop=True)
        .where((raster_DEM["y"] >= float(valid_y.min())) &
               (raster_DEM["y"] <= float(valid_y.max())), drop=True)
    )

    xllcorner = float(raster_DEM["x"].min())
    yllcorner = float(raster_DEM["y"].min())

    raster_DEM_masked = np.where(np.isnan(raster_DEM), -9999, raster_DEM.values)

    print(f"  DEM shape      : {raster_DEM_masked.shape}")
    print(f"  xllcorner={xllcorner:.2f}, yllcorner={yllcorner:.2f}")

    if show:
        valid = raster_DEM_masked[raster_DEM_masked != -9999]
        fig, ax = plt.subplots(figsize=(8, 6))
        img = ax.imshow(raster_DEM_masked, cmap="terrain",
                        vmin=valid.min(), vmax=valid.max())
        plt.colorbar(img, ax=ax, label="Elevation (m)")
        ax.set_title(f"DEM — masked ({target_crs})")
        plt.tight_layout()

    return raster_DEM, raster_DEM_masked, xllcorner, yllcorner, res_x, res_y


def initialise_project(path2prj: str, start_year: int, end_year: int,
                       zroot: float, pmin: float, scenario_dirname: str):
    """Create and return the CATHY project object inside a sim-specific folder."""
    sim_dir = Path(path2prj) / scenario_dirname
    sim_dir.mkdir(parents=True, exist_ok=True)
    simu = cathy_tools.CATHY(dirName=str(sim_dir))
    print(f"  Working directory  : {sim_dir}")
    return simu


def smooth_dem(raster_DEM_masked: np.ndarray, sigma: float) -> np.ndarray:
    """
    Gaussian-smooth DEM elevations while leaving array shape and the
    NaN (masked/outside-catchment) footprint exactly unchanged.

    Plain ``scipy.ndimage.gaussian_filter`` on an array containing NaNs
    would leak NaN into every pixel within ~sigma of the mask edge (and
    vice-versa smear elevation into masked cells), silently growing or
    shrinking the effective valid-pixel footprint. That would change
    what "the DEM" is even though .shape stays (N, M) — the grid dims
    CATHY sees would still be right, but the terrain wouldn't be.

    Standard fix: zero-fill NaNs, smooth both the filled field and a
    0/1 validity mask with the same kernel, then divide the two
    (Gaussian-weighted normalisation) so only valid neighbours
    contribute to each valid pixel's new value. NaNs are then written
    back at the original locations, so the mask is bit-for-bit
    identical to the input — only masked-in elevations move.

    Parameters
    ----------
    raster_DEM_masked : np.ndarray
        DEM elevation grid, NaN outside the catchment/basin.
    sigma : float
        Gaussian sigma in pixels. <= 0 returns the input unchanged
        (no-op — smoothing was requested but effectively disabled).

    Returns
    -------
    np.ndarray
        Smoothed DEM, same shape and same NaN mask as the input.
    """
    if sigma is None or sigma <= 0:
        return raster_DEM_masked

    dem = np.asarray(raster_DEM_masked, dtype=float)
    nan_mask = ~np.isfinite(dem)

    filled = np.where(nan_mask, 0.0, dem)
    valid  = np.where(nan_mask, 0.0, 1.0)

    filled_smooth = gaussian_filter(filled, sigma=sigma, mode="nearest")
    valid_smooth  = gaussian_filter(valid,  sigma=sigma, mode="nearest")

    with np.errstate(invalid="ignore", divide="ignore"):
        dem_smooth = filled_smooth / valid_smooth

    # Re-impose the exact original mask — shape and valid-pixel count
    # (and therefore N, M and dem_valid_px used for area/logging) are
    # unchanged; only elevations at already-valid pixels move.
    dem_smooth[nan_mask] = np.nan
    return dem_smooth.astype(raster_DEM_masked.dtype, copy=False)


def setup_dem(simu, raster_DEM_masked: np.ndarray,
              xllcorner: float, yllcorner: float,
              delta_x: float, delta_y: float,
              show: bool = False) -> None:
    """Feed DEM into CATHY preprocessing inputs and optionally show 3-D terrain."""
    simu.update_prepo_inputs(
        xllcorner=xllcorner,
        yllcorner=yllcorner,
        DEM=raster_DEM_masked,
        delta_x=delta_x,
        delta_y=delta_y,
        ivert=1,
        #base=30,
    )
    if show:
        fig = plt.figure(figsize=(10, 6))
        ax  = plt.axes(projection="3d")
        simu.show_input(prop="dem", ax=ax)
        plt.tight_layout()
    simu.create_mesh_vtk(verbose=True)


def configure_boundary_conditions(
    simu, t_atmbc: list, outlet_flux: float = -1e-6, outlet_side: str = None,
    bottom_flux: float = 0.0, res_x: float = None, res_y: float = None,
) -> tuple[str, int]:
    """
    No-flow on the bottom and on 3 of the 4 lateral sides; a small
    prescribed outflow on the topographically lowest lateral side (the
    catchment outlet), so ponded/perched water can actually leave the
    domain instead of piling up against a fully closed box. Optionally
    also a free-bottom-drainage-style outflow on the mesh bottom.

    Thin orchestration layer only — the actual geometry (face areas),
    node identification, and file writing now live in pyCATHY itself
    (``CATHY.get_outlet_side``, ``CATHY.get_outlet_node``,
    ``CATHY.compute_nodal_face_areas``, ``CATHY.compute_bottom_face_area``,
    ``CATHY.update_nansfneubc``, ``CATHY.update_nansfdirbc``,
    ``CATHY.update_sfbc``, ``CATHY.check_neumann_flux_vs_ks``) rather than
    being reimplemented per-experiment here, so any other pyCATHY script
    needing an "open outlet" BC gets the same, already-reviewed logic
    instead of a fresh copy that can drift out of sync / re-introduce
    bugs already fixed once.

    ``outlet_flux`` / ``bottom_flux`` are in m/s, negative = water
    leaving the domain — tunable parameters, not physical constants. A
    reasonable starting point is a fraction of the saturated hydraulic
    conductivity (Ks) of the relevant soil zone; ``check_neumann_flux_vs_ks``
    below warns (does not block) if the magnitude exceeds Ks.

    UNITS NOTE: ``nansfneubc`` needs a volumetric flow rate Q [m^3/s]
    per node, not a flux density [m/s] like ``atmbc``. This is converted
    internally using each node's *own* face area
    (``compute_nodal_face_areas`` / ``compute_bottom_face_area``) —
    exact per node, not a side-wide average.

    CAVEAT (physical, not a code bug): this applies the same flux
    density to *every* node on the resolved side, from the DEM surface
    down to the mesh bottom — including nodes well above any realistic
    water table. That's an "always-open tap" approximation, not a true
    saturation-aware seepage face (pyCATHY's ``update_sfbc`` doesn't
    support non-trivial seepage faces yet — see its docstring). Treat
    behaviour near the shallow outlet-side nodes with more scepticism
    than behaviour near the water table / mesh bottom. If your
    catchment's real outlet is an interior pour point rather than a
    domain edge, treating a whole mesh side as the outlet is itself an
    approximation.

    Returns
    -------
    (outlet_side, outlet_node_id)
    """
    simu.create_mesh_bounds_df("nansfneubc", simu.grid3d["mesh3d_nodes"], t_atmbc)

    resolved_side = simu.get_outlet_side(outlet_side=outlet_side)
    outlet_node_id = simu.get_outlet_node()

    # Dirichlet & seepage-face channels stay fully closed for every time
    # block in t_atmbc (an explicit "no nodes imposed" record per block,
    # not an empty/0-byte file — see CATHY.update_nansfdirbc /
    # CATHY.update_sfbc docstrings for why time=[] used to be invalid).
    simu.update_nansfdirbc(time=t_atmbc, no_flow=True)
    simu.update_sfbc(time=t_atmbc, no_flow=True)

    # ── Neumann: per-node exact face-area conversion (m/s -> m^3/s) ──
    outlet_areas = simu.compute_nodal_face_areas(resolved_side)
    outlet_node_ids = np.array(sorted(outlet_areas.keys()), dtype=int)
    outlet_flux_m3s = np.array(
        [outlet_flux * outlet_areas[n] for n in outlet_node_ids]
    )
    simu.check_neumann_flux_vs_ks(outlet_flux, label="outlet_flux")

    bottom_active = bottom_flux != 0.0
    bottom_node_ids = np.array([], dtype=int)
    bottom_flux_m3s = np.array([])
    n_bottom_nodes = 0
    if bottom_active:
        if res_x is None or res_y is None:
            raise ValueError(
                "bottom_flux is non-zero but res_x/res_y (DEM pixel size) "
                "were not provided - needed to convert m/s -> m^3/s per "
                "bottom node."
            )
        df0 = simu.mesh_bound_cond_df[
            simu.mesh_bound_cond_df["time"] == simu.mesh_bound_cond_df["time"].iloc[0]
        ]
        bottom_node_ids = df0.loc[df0["bot_bound"], "id_node"].to_numpy(dtype=int)
        n_bottom_nodes = len(bottom_node_ids)
        bottom_nodal_area = simu.compute_bottom_face_area(res_x, res_y)
        bottom_flux_m3s = np.full(n_bottom_nodes, bottom_flux * bottom_nodal_area)
        simu.check_neumann_flux_vs_ks(bottom_flux, label="bottom_flux")

    # Combine both groups into ONE nansfneubc call. update_nansfneubc
    # rewrites the whole file on every call ("w+"), so calling it twice
    # (once per group) would silently lose the first group - merge
    # node/flux arrays here instead, exactly once.
    overlap_mask = np.isin(outlet_node_ids, bottom_node_ids)
    if overlap_mask.any():
        overlap_nodes = outlet_node_ids[overlap_mask]
        print(f"  WARNING: {len(overlap_nodes)} node(s) belong to both the "
              f"outlet side and the mesh bottom (e.g. a bottom corner also "
              f"on the lateral outlet edge): {sorted(overlap_nodes.tolist())}. "
              f"bottom_flux wins for these nodes (outlet contribution "
              f"dropped for them) since it's applied after the outlet "
              f"group below.")
    keep_outlet = ~overlap_mask
    all_nodes = np.concatenate([outlet_node_ids[keep_outlet], bottom_node_ids])
    all_flux = np.concatenate([outlet_flux_m3s[keep_outlet], bottom_flux_m3s])

    # Drop exactly-zero-flux nodes before writing. An explicit Neumann
    # node with flux=0 is physically identical to an unlisted node (the
    # solver's natural BC for any node absent from nansfneubc IS
    # zero-flux), but writing it explicitly makes the file noisier than
    # it needs to be, and needlessly inflates NQMAX via update_cathyH
    # below. When outlet_flux=0 and bottom_flux=0 (fully closed run) this
    # empties the list entirely, so we fall back to the same no_flow=True
    # convention update_nansfdirbc/update_sfbc already use above, instead
    # of writing a spurious "explicit list of N nodes, all flux=0" block.
    nonzero_mask = all_flux != 0.0
    all_nodes = all_nodes[nonzero_mask]
    all_flux = all_flux[nonzero_mask]

    if len(all_nodes) == 0:
        simu.update_nansfneubc(time=t_atmbc, no_flow=True)
    else:
        simu.update_nansfneubc(time=t_atmbc, nodes=all_nodes, flux=all_flux)

    outlet_total_area = sum(outlet_areas.values())
    if len(all_nodes) == 0:
        print("  BCs: Dirichlet + seepage-face + Neumann all closed "
              "(no_flow) - outlet_flux and bottom_flux both zero.")
    else:
        print(f"  BCs: Dirichlet + seepage-face closed; Neumann outflow on "
              f"outlet side '{resolved_side}': {outlet_flux:.1e} m/s x exact "
              f"per-node face area (side total {outlet_total_area:.3g} m^2 "
              f"over {len(outlet_node_ids)} nodes). Outlet diagnostic node: "
              f"{outlet_node_id}.")
    if bottom_active:
        print(f"  Free bottom drainage: {bottom_flux:.1e} m/s x nodal "
              f"area {bottom_nodal_area:.3g} m^2, over {n_bottom_nodes} "
              f"bottom nodes.")
    else:
        print("  Bottom: closed (no-flow).")
    return resolved_side, outlet_node_id
