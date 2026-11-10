"""
agramon_lai.py — LAI → vegetation class → soil/root-depth helpers
=================================================================
GRwater project — Agramon study site

Moved verbatim out of Agramon_withLAI_withETp.py (REFACTORING.md, step 2.1).
"""

import argparse
from typing import Optional

import numpy as np
import pandas as pd
import xarray as xr
from scipy.ndimage import binary_dilation, gaussian_filter, label as sp_label
from scipy.stats import mode

from agramon_config import LAI_THRESHOLDS, N_VEG_CLASSES, ZROOT_FACTORS, scf_label_for


def lai_to_veg_map(simu, lai_2d: np.ndarray) -> tuple[np.ndarray, dict[int, int]]:
    """
    Convert a 2-D LAI field to an integer vegetation-class map with
    *consecutive* classes starting at 1, as required by pyCATHY.

    Semantic classes (thresholds = ``agramon_config.LAI_THRESHOLDS``,
    calibrated to the observed Agramon LAI range ~1 – 3 m² m⁻²)
    ---------------------------------------------------------------------------
    1 : LAI ≤ 1.0            (very low / bare; burnt, post-fire)
    2 : 1.0 < LAI ≤ 1.5     (low / sparse)
    3 : 1.5 < LAI ≤ 2.0     (moderate)
    4 : 2.0 < LAI ≤ 2.5     (high / dense)
    5 : LAI > 2.5            (very high / very dense)

    Returns
    -------
    veg_remapped : (M, N) int array — consecutive classes 1 … K
    remap        : {original_semantic_class → remapped_class}
    """
    N = int(simu.hapin["N"])
    M = int(simu.hapin["M"])
    veg = np.ones((M, N), dtype=int)
    for i, thresh in enumerate(LAI_THRESHOLDS, start=2):
        veg[lai_2d > thresh] = i

    # NOTE: remap is intentionally the IDENTITY mapping over the fixed set of
    # semantic classes 1..N_VEG_CLASSES — NOT derived from np.unique(veg).
    # Re-deriving it per month (as before) renumbered whichever classes were
    # spatially present that month to consecutive 1..K, so the same veg-index
    # meant a different semantic class from one month to the next, and the
    # SOIL table size (and MAXVEG) mismatched whenever fewer than 5 classes
    # were present. Keeping it fixed means the soil table always has
    # N_VEG_CLASSES rows and row i always corresponds to semantic class i,
    # whether or not that class actually appears in the DEM this month.
    remap: dict[int, int] = {c: c for c in range(1, N_VEG_CLASSES + 1)}
    veg_remapped = veg  # already using the fixed semantic classes 1..N_VEG_CLASSES
    return veg_remapped, remap


def smooth_veg_map(
    veg: np.ndarray,
    min_patch_pixels: int = 4,
    n_classes: int | None = None,
) -> np.ndarray:
    """
    Remove isolated / salt-and-pepper pixels from a vegetation-class map
    using connected-component minimum-patch removal.
    """
    if n_classes is None:
        n_classes = int(veg.max())

    out = veg.copy()

    for cls in range(1, n_classes + 1):
        binary = out == cls
        labeled, n_features = sp_label(binary)

        for feat_id in range(1, n_features + 1):
            patch = labeled == feat_id
            if patch.sum() >= min_patch_pixels:
                continue

            border = binary_dilation(patch, iterations=1) & ~patch
            if not border.any():
                continue

            neighbour_vals = out[border]
            replacement = int(mode(neighbour_vals, keepdims=True).mode[0])
            out[patch] = replacement

    return out


def update_soil_for_active_zones(
    simu,
    active_zones: set[int],
    remap: dict[int, int],
    base_zroot: float,
    pmin: float,
    scf: float,
    zroot_factors: dict[int, float] | None = None,
    permx: float | None = None,
    permy_ratio: float = 1.0,
    permz_ratio: float = 1.0,
    vgncell: float | None = None,
    vgrmccell: float | None = None,
    vgpsatcell: float | None = None,
    pcana: float | None = None,
    pcref: float | None = None,
    pcwlt: float | None = None,
    pz: float | None = None,
    omgc: float | None = None,
    feddes_crop: str | None = None,
) -> None:
    """
    ``zroot_factors`` overrides the module-level ZROOT_FACTORS mapping
    (semantic veg class 1-5 → root-depth multiplier) when given —
    normally built from ``--zroot-factors`` — falling back to
    ZROOT_FACTORS itself when None, so existing callers that don't pass
    it keep the original hardcoded behaviour.

    ``permx`` (m/s), when given — normally from ``--permx`` — overrides
    PERMX uniformly across every row of the SPP_map returned by
    ``simu.set_SOIL_defaults(SPP_map_default=True)``; PERMY/PERMZ follow
    via ``permy_ratio``/``permz_ratio`` (anisotropy). Left untouched
    (pyCATHY's own SPP default) when permx is None.

    ``scf`` is the root water uptake stress compensation factor, a float
    in [0, 1]: 1 = Feddes approach (ET ~ all transpiration), 0 =
    PMIN-threshold approach (sparse/arid cover), values in between blend
    the two. Passed straight through to ``simu.update_soil(SCF=scf)``
    alongside PMIN — normally just ``args.scf`` at the call site.

    ``vgncell`` / ``vgrmccell`` / ``vgpsatcell``, when given, override the
    van Genuchten retention-curve columns (VGNCELL = n, VGRMCCELL =
    residual water content θr, VGPSATCELL = 1/alpha) uniformly across
    every row of df_spp, same broadcast pattern as permx above. Each is
    left at pyCATHY's own SPP default when None. Note VGRMCCELL must
    stay below POROS — simu.update_soil() itself warns (doesn't raise)
    if VGRMCCELL >= POROS for any row.

    ``pcana`` / ``pcref`` / ``pcwlt`` / ``pz`` / ``omgc``, when given,
    override the remaining Feddes (FP_map) columns — PCANA = anaerobiosis
    point h2, PCREF = field capacity h3, PCWLT = wilting point h4, PZ =
    empirical root-water-uptake exponent, OMGC = compensatory mechanism
    for root water uptake (0 < OMGC < 1) — uniformly across every row of
    df_fp, same broadcast pattern as the VG* overrides above. ZROOT (also
    part of FP_map) is handled separately by ``base_zroot``/``zroot_factors``,
    not here. Normally resolved once via ``feddes_params_from_args(args)``
    rather than being hardcoded at each call site — see that helper for how
    ``feddes_crop`` (pyCATHY's own literature lookup,
    ``pyCATHY.petro.Feddes_litterature``) feeds pcana/pcref/pcwlt when those
    are not explicitly overridden. ``feddes_crop`` itself is accepted here
    only for the run-log label below; the actual value resolution already
    happened in the caller.
    """
    inv_remap: dict[int, int] = {v: k for k, v in remap.items()}
    if zroot_factors is None:
        zroot_factors = ZROOT_FACTORS

    # Always build the SOIL vegetation table for the FULL fixed set of
    # N_VEG_CLASSES zones — matching MAXVEG — regardless of how many of
    # them are spatially present in the DEM this particular month. Zones
    # not currently active still get a physically-meaningful ZROOT (their
    # fixed zroot_factors value), they just aren't occupying any pixels
    # this month; they must not be dropped from the table, or the row
    # count stops matching MAXVEG and/or rows shift relative to the
    # semantic class each veg-index in the raster refers to.
    n_zones = N_VEG_CLASSES

    df_fp  = simu.set_SOIL_defaults(FP_map_default=True, nveg=n_zones)
    df_spp = simu.set_SOIL_defaults(SPP_map_default=True)

    template = df_fp.iloc[[0]]
    df_fp = pd.concat([template] * n_zones, ignore_index=True)
    # pd.concat(..., ignore_index=True) resets the index to a plain 0-based
    # RangeIndex(0, n_zones). But set_SOIL_defaults()/init_soil_FP_map_df()
    # index this table 1..n_zones ("Veg nb"), and cathy_tools.py's
    # _prepare_SOIL_vegetation_tb() reads it back with
    # FP_map[sfp].loc[iveg] for iveg in range(1, MAXVEG + 1) — so losing
    # the 1-based labels here means the last class (iveg == n_zones) is
    # never found, raising KeyError: n_zones. Re-attach the expected
    # labels; row order is already class 1..n_zones, so this is just
    # restoring the index, not reordering anything.
    df_fp.index = pd.RangeIndex(1, n_zones + 1)
    df_fp.index.name = "Veg nb"
    for key in df_fp.columns:
        if key == "ZROOT":
            df_fp[key] = [
                base_zroot * zroot_factors.get(inv_remap.get(z, z), 1.0)
                for z in range(1, n_zones + 1)
            ]

    permx_label = "pyCATHY default"
    if permx is not None:
        overrides = {
            "PERMX": permx,
            "PERMY": permx * permy_ratio,
            "PERMZ": permx * permz_ratio,
        }
        for col, val in overrides.items():
            if col in df_spp.columns:
                df_spp[col] = val
            else:
                print(f"    WARNING: SPP_map has no '{col}' column — "
                      f"--permx/--permy-ratio/--permz-ratio ignored for it.")
        permx_label = (f"PERMX={permx:.3e}, PERMY={overrides['PERMY']:.3e}, "
                       f"PERMZ={overrides['PERMZ']:.3e} m/s")

    # Retention curve (van Genuchten) overrides — VGNCELL/VGRMCCELL/
    # VGPSATCELL — same broadcast-to-every-row pattern as permx above.
    vg_overrides = {
        "VGNCELL": vgncell,
        "VGRMCCELL": vgrmccell,
        "VGPSATCELL": vgpsatcell,
    }
    vg_label_parts = []
    for col, val in vg_overrides.items():
        if val is None:
            continue
        if col in df_spp.columns:
            df_spp[col] = val
            vg_label_parts.append(f"{col}={val:.4g}")
        else:
            print(f"    WARNING: SPP_map has no '{col}' column — override ignored for it.")
    if (df_spp["VGRMCCELL"] >= df_spp["POROS"]).any():
        print("    WARNING: VGRMCCELL >= POROS for at least one zone — "
              "residual water content exceeds porosity, fix before running.")
    vg_label = ", ".join(vg_label_parts) if vg_label_parts else "pyCATHY default"

    # Feddes (FP_map) overrides — PCANA/PCREF/PCWLT/PZ/OMGC — same
    # broadcast-to-every-row pattern as the VG* overrides above. ZROOT is
    # already handled per-zone by the loop earlier in this function.
    fp_overrides = {
        "PCANA": pcana,
        "PCREF": pcref,
        "PCWLT": pcwlt,
        "PZ": pz,
        "OMGC": omgc,
    }
    fp_label_parts = []
    for col, val in fp_overrides.items():
        if val is None:
            continue
        if col in df_fp.columns:
            df_fp[col] = val
            fp_label_parts.append(f"{col}={val:.4g}")
        else:
            print(f"    WARNING: FP_map has no '{col}' column — override ignored for it.")
    if feddes_crop is not None:
        fp_label_parts.append(f"crop={feddes_crop} (literature)")
    fp_label = ", ".join(fp_label_parts) if fp_label_parts else "pyCATHY default"

    simu.update_soil(FP_map=df_fp, SPP_map=df_spp, PMIN=pmin, SCF=scf)
    semantic_active = sorted(inv_remap.get(z, z) for z in active_zones)
    scf_label = scf_label_for(scf)
    print(f"    → Soil updated | table rows: {n_zones} (all zones, fixed) "
          f"| spatially active this month: {sorted(active_zones)} "
          f"(semantic: {semantic_active}) | PMIN: {pmin:g} m | SCF: {scf_label} | "
          f"Ks: {permx_label} | retention curve: {vg_label} | "
          f"Feddes: {fp_label}")


def iter_monthly_lai(lai_ds: xr.Dataset):
    """Yield the first (datetime, lai_2d) of each calendar month."""
    seen = set()
    for t in lai_ds["time"].values:
        dt  = pd.Timestamp(t).to_pydatetime()
        key = (dt.year, dt.month)
        if key in seen:
            continue
        seen.add(key)
        yield dt, lai_ds["LAI"].sel(time=t).values
