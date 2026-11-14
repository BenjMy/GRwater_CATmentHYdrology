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
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


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


def zroot_factors_from_args(args: argparse.Namespace) -> dict[int, float]:
    """
    Build the semantic-zone → root-depth-multiplier mapping consumed by
    ``update_soil_for_active_zones`` from ``--zroot-factors`` (5 values,
    positional order: bare, sparse, moderate, dense, very dense —
    matching ``ZROOT_FACTORS``' key order 1-5). Centralised here so the
    spin-up and main-loop call sites build the exact same mapping from
    the exact same CLI values instead of risking drift between the two.
    """
    return {i + 1: v for i, v in enumerate(args.zroot_factors)}


def feddes_params_from_args(args: argparse.Namespace) -> dict[str, float | None]:
    """
    Resolve PCANA/PCREF/PCWLT/PZ/OMGC (the 5 non-ZROOT Feddes parameters —
    ZROOT is handled separately via --zroot/--zroot-factors) from CLI args,
    consumed by ``update_soil_for_active_zones``.

    When ``--feddes-crop`` is given, ``pyCATHY.petro.Feddes_litterature()``
    is used to look up literature PCANA/PCREF/PCWLT values for that crop
    (pyCATHY's own "cathy tool" for inferring Feddes parameters — see
    https://benjmy.github.io/pycathy_wrapper/content/SSHydro/plot_2_pyCATHY_inputs.html
    for the underlying FP_map convention) as the base value for those three
    columns instead of pyCATHY's hardcoded FP_map_default. --pz/--omgc have
    no literature equivalent in Feddes_litterature() and always fall back
    to pyCATHY's own default when unset.

    An explicit --pcana/--pcref/--pcwlt always wins over the literature
    lookup, matching how every other override in this script (--permx,
    --vgncell, etc.) lets an explicit CLI value take precedence over
    whatever baseline it would otherwise fall back to. Centralised here so
    the spin-up and main-loop call sites resolve the exact same values from
    the exact same CLI inputs instead of risking drift between the two.
    """
    lit_pcana = lit_pcref = lit_pcwlt = None
    lit_ref = None
    if args.feddes_crop is not None:
        from pyCATHY.petro import Feddes_litterature   # lazy: keeps this module importable without pyCATHY
        lit = Feddes_litterature(crop=args.feddes_crop)
        lit_pcana = float(lit["PCANA"])
        lit_pcref = float(lit["PCREF"])
        lit_pcwlt = float(lit["PCWLT"])
        lit_ref = lit.get("ref")

    return {
        "pcana": args.pcana if args.pcana is not None else lit_pcana,
        "pcref": args.pcref if args.pcref is not None else lit_pcref,
        "pcwlt": args.pcwlt if args.pcwlt is not None else lit_pcwlt,
        "pz": args.pz,
        "omgc": args.omgc,
        "feddes_crop": args.feddes_crop,
        "feddes_crop_ref": lit_ref,
    }


# ── Forcing mode (replaces the former module globals, REFACTORING.md 2.13) ─────
#   "era5"   -> ERA5 point series (deprecated, quick-test only)
#   "etp"    -> satellite ETp raster + satellite rain raster (PathET files)
#   "etp_ts" -> ETp from the Agramon ET CSV time series (uniform over all pixels,
#               treatment selected by --dem-plot) + satellite rain raster
FORCING_MODE_DEFAULT = "etp"
FORCING_LABELS = {
    "era5": "ERA5",
    "etp": "satellite ETp + rain",
    "etp_ts": "ETp time series (CSV) + satellite rain",
}


@dataclass(frozen=True)
class ForcingMode:
    """Immutable forcing choice of one run (``--forcing``), carried in ``RunInputs``.

    "etp" and "etp_ts" share the whole per-node ETp/rain machinery (they differ only in
    how ``ds_etp`` is built), hence ``use_etp`` is True for both.
    """
    name: str

    def __post_init__(self):
        if self.name not in FORCING_LABELS:
            raise ValueError(f"forcing must be one of {list(FORCING_LABELS)}, got {self.name!r}")

    @classmethod
    def from_name(cls, name: str, warn: bool = True) -> "ForcingMode":
        fm = cls(name)
        if warn and name == "era5":
            warnings.warn(
                "--forcing era5 is deprecated (quick-test only) and will be removed; "
                "use 'etp' or 'etp_ts' (see REFACTORING.md, step 3.3).",
                DeprecationWarning, stacklevel=2)
        return fm

    @property
    def label(self) -> str:
        return FORCING_LABELS[self.name]

    @property
    def use_etp(self) -> bool:
        return self.name in ("etp", "etp_ts")

    @property
    def use_etp_ts(self) -> bool:
        return self.name == "etp_ts"



# ─────────────────────────────────────────────────────────────
# Boundary / initial conditions (REFACTORING.md 2.9)
# ─────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class BCConfig:
    """
    Boundary and initial conditions of one run (replaces seven loose
    ``args.*`` reads).

    Units: fluxes in m/s, negative = water leaving the domain (Neumann);
    ``ic_wtposition`` is the water-table depth below the surface (m, >= 0);
    ``ic_pressure_head`` is a uniform pressure head (m, negative = unsaturated).
    ``outlet_side`` is None for auto-detection from the DEM.
    """
    outlet_flux: float = 0.0
    outlet_side: Optional[str] = None
    bottom_flux: float = 0.0
    ic_mode: str = "wt"
    ic_wtposition: float = 1.0
    ic_pressure_head: float = -10.0

    def __post_init__(self):
        sides = ("xmin_bound", "xmax_bound", "ymin_bound", "ymax_bound")
        if self.ic_mode not in ("wt", "pressure"):
            raise ValueError(f"ic_mode must be 'wt' or 'pressure', got {self.ic_mode!r}")
        if self.outlet_side is not None and self.outlet_side not in sides:
            raise ValueError(f"outlet_side must be one of {sides} or None, "
                             f"got {self.outlet_side!r}")
        if self.ic_mode == "wt" and self.ic_wtposition < 0:
            raise ValueError("ic_wtposition is a depth below the surface (>= 0 m)")

    @classmethod
    def from_args(cls, args) -> "BCConfig":
        return cls(outlet_flux=args.outlet_flux, outlet_side=args.outlet_side,
                   bottom_flux=args.bottom_flux, ic_mode=args.ic_mode,
                   ic_wtposition=args.ic_wtposition,
                   ic_pressure_head=args.ic_pressure_head)

    @property
    def ic_label(self) -> str:
        if self.ic_mode == "wt":
            return f"WT depth {self.ic_wtposition} m (INDP=3)"
        return f"uniform pressure head {self.ic_pressure_head} m (INDP=0)"

    @property
    def bottom_label(self) -> str:
        return (f"{self.bottom_flux:.1e} m/s" if self.bottom_flux != 0.0
                else "disabled (closed)")

    def apply_ic(self, simu) -> None:
        """Set the initial condition on a pyCATHY ``simu``."""
        if self.ic_mode == "pressure":
            simu.update_ic(INDP=0, pressure_head_ini=self.ic_pressure_head)
        else:
            simu.update_ic(INDP=3, WTPOSITION=self.ic_wtposition)
