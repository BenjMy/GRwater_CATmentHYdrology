"""
agramon_config.py — Paths and constants for the Agramon pyCATHY pipeline
========================================================================
GRwater project — Agramon study site

Single place for everything that is a *setting* rather than logic:
machine-dependent data paths, CRS codes, the LAI → vegetation-class
scheme, SCF default, and the monthly sanity-check tolerances.

Nothing in here runs code with side effects beyond resolving paths.
"""

import argparse
from pathlib import Path


# Column of the Agramon CSVs used as ETp (mm/day) in "etp_ts" mode: either
# "ET_0_wet" or "ET_0_dry" (selectable with --etp-ts-column).
ETP_TS_COLUMN_DEFAULT = "ET_0_wet"


# ── Data paths ────────────────────────────────────────────────────────────────
# TODO (REFACTORING.md, step 2.4): replace the machine auto-detection below
# with environment variables / a config file.

# ── Path setup ────────────────────────────────────────────────────────────────

# Auto-detect the machine from the home directory (no commenting/uncommenting).
HOME = Path.home()
if HOME.name == "z0272571a@campus.csic.es":
    ROOT = HOME / "Nextcloud"
else:
    ROOT = Path("/home/ben/Nextcloud")

CODES = ROOT / "BenCSIC" / "Codes" / "Tech4agro_org"

# MODULE_PATH = (CODES / "GRwater_geophy").resolve()
EO_PATH     = (CODES / "GRwater_CATmentHYdrology").resolve()

LAI_PATH = (
    ROOT / "BenCSIC" / "Training_Supervision" / "Supervision"
    / "Xela_Carracedo_Practicas_2026_data" / "301a-biophysical" / "agramon" / "input"
).resolve()

PathET    = ROOT / "GRwater" / "data" / "satellite" / "Agramon"
etp_path  = PathET / "20161001_20241231_ET0.nc"
rain_path = PathET / "20161001_20241231_TP.nc"


# ── Constants ─────────────────────────────────────────────────────────────────

# All spatial datasets are reprojected to this CRS before any operation.
# DEM / shapefiles : EPSG:25830 (UTM zone 30N, ETRS89)
# LAI / ETp / rain : EPSG:32630 (UTM zone 30N, WGS84) — same zone, different datum
TARGET_CRS = "EPSG:25830"
NATIVE_CRS = "EPSG:32630"

# LAI class boundaries  [very low | low | moderate | high | very high]
# Calibrated to the observed LAI range (~1 – 3 m² m⁻²), spaced evenly so each
# class is actually used: pre-fire canopy mostly in classes 4-5, burnt /
# post-fire pixels (LAI ~ 1) falling into class 1.
# (Previous values for a 0 – 0.6 range were [0.1, 0.25, 0.4, 0.55]; with
# LAI >= 1 those put every pixel in the top class.)
LAI_THRESHOLDS = [1.0, 1.5, 2.0, 2.5]

# Fixed number of semantic vegetation classes (very low / low / moderate /
# high / very high; the CLI calls them bare / sparse / moderate / dense / very
# dense). This must stay constant for the whole run and match MAXVEG passed
# to simu.update_cathyH(), since pyCATHY's SOIL vegetation table has a fixed
# row count for the whole project — it is NOT allowed to vary month to month
# with however many classes happen to be spatially present.
N_VEG_CLASSES = len(LAI_THRESHOLDS) + 1  # = 5

# Vegetation zone → root-depth multiplier
ZROOT_FACTORS = {
    1: 0.05,   # very low LAI  (LAI <= 1.0; burnt / post-fire)
    2: 0.25,   # low           (1.0 < LAI <= 1.5)
    3: 1.0,    # moderate      (1.5 < LAI <= 2.0)
    4: 1.5,    # high          (2.0 < LAI <= 2.5)
    5: 2.0,    # very high     (LAI > 2.5)
}


# CF (root water uptake stress compensation factor) ──────────────────────
# SCF is a float in [0, 1] blending pyCATHY's two ET-partitioning approaches:
#   SCF = 1  Feddes approach — assumes ET is essentially all transpiration.
#            Reasonable when vegetation cover is dense/continuous.
#   SCF = 0  PMIN-threshold approach — appropriate for sparse-cover, arid
#            ecosystems where bare soil evaporation is a major part of ET
#            (also depends on spatial resolution).
#   0 < SCF < 1  intermediate blend of the two.
# Set explicitly via --scf. No automatic inference from PMIN or anything
# else — SCF is just a plain CLI parameter like ZROOT/PMIN.
SCF_DEFAULT = 1.0


def scf_type(raw: str) -> float:
    """argparse type= for --scf: a float, restricted to [0, 1]."""
    val = float(raw)
    if not (0.0 <= val <= 1.0):
        raise argparse.ArgumentTypeError(
            f"--scf must be between 0 and 1 (got {raw})"
        )
    return val


def scf_label_for(scf: float) -> str:
    """Human-readable label for a given SCF value, used in console/markdown
    logs. SCF=0/1 are named after pyCATHY's two approaches; anything in
    between is reported as a blend of the two."""
    if scf == 1.0:
        return f"{scf:g} (Feddes)"
    if scf == 0.0:
        return f"{scf:g} (PMIN-threshold)"
    return f"{scf:g} (blend of PMIN-threshold/Feddes)"


# ── Sanity-check thresholds (monthly loop) ─────────────────────────────────
# Tolerances for the physical-bounds and water-balance checks run after
# every month (see `check_monthly_sanity`). Not physical hard limits —
# just flags for "this deserves a look", so kept generous on purpose.
SW_TOL           = 0.01     # SW allowed slightly outside [0, 1] before flagging
PSI_SANITY_MIN   = -50.0    # m — PSI below this is flagged as extreme drying
ETA_TOL          = 1e-8     # m/s — ETa allowed slightly above ETp before flagging
WB_RESIDUAL_FRAC = 0.5      # flag if |water balance residual| > this fraction of rain_in
