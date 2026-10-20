"""
plot_fig8_ETa_ETp_vs_psi.py
============================
GRwater project — Agramon study site

Reproduces Fig. 8 of Camporese et al. (2014, Adv. Water Resour. 69:95-105):
the ratio of daily actual to potential evapotranspiration (ETa/ETp) plotted
against pressure head averaged over the shallowest ~10 cm of soil, faceted
by psi_min (one row per scenario) and by location (two columns: outlet vs.
mid-uphill), with a total-least-squares (orthogonal regression) trend line
overlaid per panel — the same fitting procedure the paper describes.

Two reference locations, matching the paper's "close to the creek, at the
catchment outlet" vs. "upslope, near the water divide":
    outlet — lowest surface node
    mid    — a representative mid-elevation band, upslope
(same node-resolution machinery as Fig 6/7 in withLAI_results.py —
resolve_nodes() / find_representative_mid_node() — just evaluated at
shallow depths instead of [0, 1, 2] m.)

Requires a new artefact, etp_output.nc (per-surface-node ETp forcing,
dims datetime x node), written by Agramon_withLAI_withETp.py right after
out_dir is created — alongside psi_output.nc / sw_output.nc / et_output.nc.
If a scenario predates that patch, rerun it (USE_ETP_FORCING must be True)
before this script can use it.

Usage
-----
    python plot_fig8_ETa_ETp_vs_psi.py --scenarios 4 5 6
    python plot_fig8_ETa_ETp_vs_psi.py --scenarios 4 5 6 --path2prj ../SSHydro_withETp_withLAI
    python plot_fig8_ETa_ETp_vs_psi.py --scenarios 4 5 6 --etp-min-mmday 1.0
    python plot_fig8_ETa_ETp_vs_psi.py --scenarios 4 5 6 --shallow-depths 0 0.05 0.1

Notes / assumptions
--------------------
- "Shallowest 10 cm" is approximated by averaging psi across the requested
  --shallow-depths node set (default: 0 and 0.1 m). The exact node depths
  actually resolved depend on the mesh's own vertical discretisation
  (simu.find_nearest_node picks the nearest existing node to each
  requested depth), so this is a discretisation-dependent proxy for a
  true 0-10 cm average, not an exact integral.
- ETa/ETp ratios are computed from raw m/s fluxes (units cancel in the
  ratio), then clipped to [0, 1.05] to allow small numerical overshoot
  while keeping the physically-expected [0, 1] range visible.
- Days with ETp below --etp-min-mmday are dropped before computing the
  ratio, to avoid blow-up when potential demand is ~0 (e.g. on rainy
  days) — the paper doesn't detail this, but some such floor is implicit
  in any sane implementation of the ratio.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData

# Reuses the project's own node-resolution / artefact-loading machinery
# rather than re-implementing it — see resolve_nodes() / load_artefacts()
# in withLAI_results.py. Importing it (rather than duplicating everything
# the way analysis_ET_compare.py deliberately does for its lighter
# dependency set) is fine here because this script needs pyCATHY anyway
# for node resolution, same as withLAI_results.py itself.
import withLAI_results as wlr
import results_plotter as plotter

# =============================================================================
# Constants
# =============================================================================

fig_path = Path().cwd() / "../Analysis/Fig8_ETa_ETp_vs_psi"

LOCATIONS = ("outlet", "mid")
LOCATION_LABELS = {
    "outlet": "Near outlet (creek)",
    "mid":    "Upslope (mid-elevation)",
}
LOCATION_COLORS = {
    "outlet": "#2166ac",
    "mid":    "#d6604d",
}


# =============================================================================
# Artefact loading
# =============================================================================

def load_etp_forcing(out_dir: Path) -> xr.Dataset:
    """Load etp_output.nc (per-surface-node ETp forcing, dims datetime x node).

    Written by the ETp-persistence patch in Agramon_withLAI_withETp.py,
    right after out_dir is created — a static (time, node) artefact, not
    accumulated through the monthly loop like psi/sw/et.
    """
    etp_nc = out_dir / "etp_output.nc"
    if not etp_nc.exists():
        raise FileNotFoundError(
            f"'etp_output.nc' not found in {out_dir}.\n"
            "  -> This artefact is written by the ETp-persistence patch in "
            "Agramon_withLAI_withETp.py (right after out_dir is created). "
            "Rerun this scenario with USE_ETP_FORCING=True to generate it."
        )
    return xr.open_dataset(etp_nc)


# =============================================================================
# Reference locations — mirrors resolve_nodes() in withLAI_results.py, but
# also keeps the raw (X, Y) coordinates (needed to sample the gridded ETa
# raster) alongside the shallow node ids (needed to sample psi/ETp).
# =============================================================================

def resolve_shallow_nodes_and_xy(
    simu,
    grid3d: dict,
    *,
    outlet_xy: tuple[float, float] | None = None,
    uphill_xy: tuple[float, float] | None = None,
    shallow_depths: tuple[float, ...] = (0.0, 0.1),
) -> dict:
    """
    Resolve outlet / mid-uphill locations at each of *shallow_depths*.

    Returns a dict with:
      "outlet_xy", "mid_xy"       — (x, y) tuples, for sampling the ETa raster
      "outlet_nodes", "mid_nodes" — list (one entry per shallow_depths, in
                                     order) of node id (outlet) or list of
                                     node ids (mid — a band, mirroring
                                     resolve_nodes()'s mid-elevation handling)
    """
    nnod          = int(grid3d["nnod"])
    nodes         = grid3d["mesh3d_nodes"]
    surface_nodes = nodes[:nnod]

    if outlet_xy:
        ox, oy = outlet_xy
        dists  = np.hypot(surface_nodes[:, 0] - ox, surface_nodes[:, 1] - oy)
        outlet_idx = int(dists.argmin())
    else:
        outlet_idx = int(surface_nodes[:, 2].argmin())
    outlet_x, outlet_y, outlet_z = surface_nodes[outlet_idx]

    mid_band_idx = None
    if uphill_xy:
        ux, uy  = uphill_xy
        dists   = np.hypot(surface_nodes[:, 0] - ux, surface_nodes[:, 1] - uy)
        mid_idx = int(dists.argmin())
        mid_x, mid_y, mid_z = surface_nodes[mid_idx]
    else:
        mid_band_idx = plotter.find_representative_mid_node(surface_nodes, return_band=True)
        mid_idx      = plotter.find_representative_mid_node(surface_nodes, return_band=False)
        mid_x, mid_y, mid_z = surface_nodes[mid_idx]

    print(f"  Outlet: X={outlet_x:.2f}, Y={outlet_y:.2f}, Z={outlet_z:.2f}")
    print(f"  Mid   : X={mid_x:.2f}, Y={mid_y:.2f}, Z={mid_z:.2f}"
          + (f"  (band of {len(mid_band_idx)} nodes)" if mid_band_idx is not None else "  (single point)"))

    outlet_nodes_by_depth: list = []
    mid_nodes_by_depth:    list = []
    for depth in shallow_depths:
        outlet_nodes_by_depth.append(
            plotter.unwrap_node_id(simu.find_nearest_node([outlet_x, outlet_y, outlet_z - depth]))
        )
        if mid_band_idx is not None:
            if depth == 0:
                # Surface band indices are already node indices at depth 0.
                mid_nodes_by_depth.append([int(i) for i in mid_band_idx])
            else:
                mid_nodes_by_depth.append([
                    plotter.unwrap_node_id(
                        simu.find_nearest_node(
                            [surface_nodes[i, 0], surface_nodes[i, 1], surface_nodes[i, 2] - depth]
                        )
                    )
                    for i in mid_band_idx
                ])
        else:
            mid_nodes_by_depth.append(
                plotter.unwrap_node_id(simu.find_nearest_node([mid_x, mid_y, mid_z - depth]))
            )

    return dict(
        outlet_xy=(float(outlet_x), float(outlet_y)),
        mid_xy=(float(mid_x), float(mid_y)),
        outlet_nodes=outlet_nodes_by_depth,
        mid_nodes=mid_nodes_by_depth,
    )


# =============================================================================
# Series extraction
# =============================================================================

def _mean_over_nodes(da: xr.DataArray, node_ids) -> np.ndarray:
    """Per-timestep mean across one or more node ids (scalar id or a band)."""
    if isinstance(node_ids, (list, tuple, np.ndarray)):
        ids = np.atleast_1d(np.asarray(node_ids, dtype=int))
        return da.sel(node=ids).mean(dim="node").values
    return da.sel(node=int(node_ids)).values


def shallow_psi_series(psi_xr: xr.Dataset, nodes_by_depth: list) -> pd.Series:
    """Mean psi across the requested shallow depths — the paper's 'pressure
    head averaged over the shallowest 10 cm of soil' (see module docstring
    for the discretisation caveat)."""
    da       = psi_xr["psi"]
    time_dim = plotter._infer_time_dim(psi_xr)
    per_depth = [_mean_over_nodes(da, ids) for ids in nodes_by_depth]
    psi_mean  = np.mean(np.stack(per_depth, axis=0), axis=0)
    return pd.Series(psi_mean, index=pd.DatetimeIndex(da[time_dim].values))


def etp_series_at_surface(etp_xr: xr.Dataset, surface_node_ids) -> pd.Series:
    """ETp [m/s] at the surface node(s) for one location (the depth=0 entry
    of resolve_shallow_nodes_and_xy's node lists)."""
    da       = etp_xr["ETp"]
    time_dim = plotter._infer_time_dim(etp_xr)
    vals     = _mean_over_nodes(da, surface_node_ids)
    return pd.Series(vals, index=pd.DatetimeIndex(da[time_dim].values))


def eta_series_at_xy(ET_xr_all: xr.Dataset, et_var: str, xy: tuple[float, float]) -> pd.Series:
    """ETa [m/s] at the nearest et_output.nc raster pixel to *xy*."""
    x, y     = xy
    time_dim = plotter._infer_time_dim(ET_xr_all)
    da       = ET_xr_all[et_var].sel(X=x, Y=y, method="nearest")
    return pd.Series(da.values, index=pd.DatetimeIndex(da[time_dim].values))


def build_location_series(
    *,
    psi_xr: xr.Dataset,
    etp_xr: xr.Dataset,
    ET_xr_all: xr.Dataset,
    et_var: str,
    resolved: dict,
    location: str,
    etp_min_mmday: float,
) -> pd.DataFrame:
    """Coincident daily (psi, ETa/ETp) pairs for one location."""
    nodes_by_depth   = resolved[f"{location}_nodes"]
    surface_node_ids = nodes_by_depth[0]           # the depth=0 entry
    xy               = resolved[f"{location}_xy"]

    psi_s = shallow_psi_series(psi_xr, nodes_by_depth)
    etp_s = etp_series_at_surface(etp_xr, surface_node_ids)
    eta_s = eta_series_at_xy(ET_xr_all, et_var, xy)

    df = pd.concat({"psi": psi_s, "etp": etp_s, "eta": eta_s}, axis=1).dropna()

    etp_min_ms = etp_min_mmday * 1e-3 / 86400   # mm/day -> m/s
    df = df[df["etp"] > etp_min_ms].copy()
    df["ratio"] = (df["eta"] / df["etp"]).clip(lower=0.0, upper=1.05)
    return df[["psi", "ratio"]]


# =============================================================================
# Fitting — "standard total least squares", as in the paper's Fig. 8
# =============================================================================

def fit_tls(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Orthogonal-distance (total-least-squares) linear fit y = a*x + b."""
    def _linear(beta, x):
        return beta[0] * x + beta[1]

    p0  = np.polyfit(x, y, 1)   # OLS slope/intercept as the ODR initial guess
    odr = ODR(RealData(x, y), Model(_linear), beta0=p0)
    out = odr.run()
    return float(out.beta[0]), float(out.beta[1])


# =============================================================================
# Scenario loading
# =============================================================================

def load_scenario(
    sim_index: int,
    *,
    path2prj: Path,
    log_file: str,
    outlet_xy: tuple[float, float] | None,
    uphill_xy: tuple[float, float] | None,
    shallow_depths: tuple[float, ...],
    et_var: str,
    etp_min_mmday: float,
) -> dict:
    _, log_row = wlr._load_log_row(path2prj / "..", log_file, sim_index)
    out_dir    = wlr._find_scenario_dir(path2prj / "outputs", sim_index).resolve()
    pmin       = float(log_row["PMIN"]) if "PMIN" in log_row else float("nan")

    art = wlr.load_artefacts(out_dir, path2prj=path2prj, sim_index=sim_index)
    simu, psi_xr, ET_xr_all = art["simu"], art["psi_xr"], art["ET_xr_all"]
    if simu is None:
        raise RuntimeError(
            f"scenario {sim_index}: CATHY project could not be loaded — "
            "node resolution needs it (see load_artefacts warnings above)."
        )
    grid3d = simu.read_outputs("grid3d")
    etp_xr = load_etp_forcing(out_dir)

    resolved = resolve_shallow_nodes_and_xy(
        simu, grid3d,
        outlet_xy=outlet_xy, uphill_xy=uphill_xy,
        shallow_depths=shallow_depths,
    )

    result: dict = {"pmin": pmin, "out_dir": out_dir}
    for location in LOCATIONS:
        result[location] = build_location_series(
            psi_xr=psi_xr, etp_xr=etp_xr, ET_xr_all=ET_xr_all, et_var=et_var,
            resolved=resolved, location=location, etp_min_mmday=etp_min_mmday,
        )
    return result


# =============================================================================
# Figure
# =============================================================================

def plot_fig8(
    scenario_data: dict,
    fig_dir: Path,
    *,
    dpi: int = 150,
    show: bool = True,
) -> Path:
    sim_indices = sorted(scenario_data)
    n_rows = len(sim_indices)
    fig, axes = plt.subplots(
        n_rows, 2, figsize=(9, 3.2 * n_rows), squeeze=False, sharey=True
    )

    for row, sim_index in enumerate(sim_indices):
        entry = scenario_data[sim_index]
        pmin  = entry["pmin"]
        for col, location in enumerate(LOCATIONS):
            ax = axes[row][col]
            df = entry[location]

            if df.empty:
                ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
                ax.set_axis_off()
                continue

            ax.scatter(
                df["psi"], df["ratio"], s=10, alpha=0.35,
                color=LOCATION_COLORS[location], edgecolors="none",
            )

            if len(df) >= 5:
                slope, intercept = fit_tls(df["psi"].to_numpy(), df["ratio"].to_numpy())
                x_line = np.linspace(df["psi"].min(), df["psi"].max(), 50)
                ax.plot(x_line, slope * x_line + intercept, "k--", lw=1.3)
                ax.text(
                    0.03, 0.93, f"TLS slope = {slope:.3g}",
                    transform=ax.transAxes, fontsize=7.5, va="top",
                )

            ax.axhline(0, color="0.75", lw=0.6)
            ax.axhline(1, color="0.75", lw=0.6)
            ax.set_ylim(-0.05, 1.1)

            if row == 0:
                ax.set_title(LOCATION_LABELS[location], fontsize=10, fontweight="bold")
            if row == n_rows - 1:
                ax.set_xlabel(r"$\psi$ (m)")
            if col == 0:
                ax.set_ylabel(f"$\\psi_{{min}}$ = {pmin:g} m\n" + r"$ET_a / ET_p$")

    fig.suptitle(
        r"$ET_a/ET_p$ vs. shallow ($\leq$10 cm) pressure head"
        " — cf. Fig. 8, Camporese et al. (2014)",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    fig_dir.mkdir(parents=True, exist_ok=True)
    out = plotter.save_fig(fig, fig_dir, "fig8_ETa_ETp_vs_psi.png", dpi=dpi)
    plotter.maybe_show(fig, show)
    return out


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Reproduce Fig. 8 of Camporese et al. (2014): ETa/ETp vs. shallow psi",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--scenarios", nargs="+", type=int, required=True, metavar="N",
        help=(
            "sim_index values to compare, one row per scenario in the "
            "figure — each should have been run with a different --pmin "
            "(see Agramon_withLAI_withETp.py)."
        ),
    )
    p.add_argument("--path2prj", default="../SSHydro_withETp_withLAI", metavar="PATH")
    p.add_argument("--log-file", default="simulation_log_withETp_withLAI.csv", metavar="PATH")
    p.add_argument("--fig-dir", default=None, metavar="PATH")
    p.add_argument(
        "--et-var", default=wlr.ET_VAR, choices=["ACT. ETRA", "ACT. ETRA_patched"],
        help="Which actual-ET variable to read from et_output.nc.",
    )
    p.add_argument(
        "--shallow-depths", nargs="+", type=float, default=[0.0, 0.1], metavar="M",
        help="Depth offsets (m) averaged for the shallow-psi x-axis.",
    )
    p.add_argument("--outlet-xy", nargs=2, type=float, default=None, metavar=("X", "Y"))
    p.add_argument("--uphill-xy", nargs=2, type=float, default=None, metavar=("X", "Y"))
    p.add_argument(
        "--etp-min-mmday", type=float, default=0.5, metavar="MM/DAY",
        help="Drop days with ETp below this before computing the ratio.",
    )
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--no-plots", dest="show", action="store_false")
    p.set_defaults(show=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    path2prj = Path(args.path2prj).resolve()

    scenario_data: dict = {}
    for sim_index in args.scenarios:
        print(f"\n── Scenario {sim_index} ──────────────────────────────────")
        scenario_data[sim_index] = load_scenario(
            sim_index,
            path2prj=path2prj,
            log_file=args.log_file,
            outlet_xy=tuple(args.outlet_xy) if args.outlet_xy else None,
            uphill_xy=tuple(args.uphill_xy) if args.uphill_xy else None,
            shallow_depths=tuple(args.shallow_depths),
            et_var=args.et_var,
            etp_min_mmday=args.etp_min_mmday,
        )
        for location in LOCATIONS:
            n = len(scenario_data[sim_index][location])
            print(
                f"  {location:>6s}: {n} coincident daily points "
                f"(psi_min = {scenario_data[sim_index]['pmin']:g} m)"
            )

    fig_dir = (
        Path(args.fig_dir).resolve() if args.fig_dir
        else path2prj / "figures" / "fig8_ETa_ETp_vs_psi"
    )
    plot_fig8(scenario_data, fig_dir, dpi=args.dpi, show=args.show)
    print("\nDone.")


if __name__ == "__main__":
    main()
