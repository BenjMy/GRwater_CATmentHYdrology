"""
agramon_spinup.py - spin-up (pre-2016 equilibration) for the Agramon pipeline
============================================================================
Moved unchanged from Agramon_withLAI_withETp.py (REFACTORING.md 2.5); the main script
re-exports the names, so existing imports keep working.
"""
from __future__ import annotations

import argparse
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from agramon_config import N_VEG_CLASSES, feddes_params_from_args, zroot_factors_from_args
from agramon_forcing import subset_atmbc_for_month
from agramon_lai import (
    iter_monthly_lai,
    lai_to_veg_map,
    smooth_veg_map,
    update_soil_for_active_zones,
)

if TYPE_CHECKING:            # type hint only; agramon_runlog needs the private geoutils
    from agramon_runlog import MarkdownLog


def pick_default_spinup_year(valid_times: np.ndarray) -> int:
    """
    Pick a full calendar year of forcing to use as the repeating spin-up
    cycle. Skips the first calendar year present in the forcing record if
    it is incomplete (satellite ETp/rain start 2016-10-01, so 2016 only
    has Oct–Dec and cannot be cycled as a full year).
    """
    times      = pd.to_datetime(valid_times)
    first_year = int(times.year.min())
    months_in_first_year = times[times.year == first_year].month
    if len(np.unique(months_in_first_year)) < 12:
        return first_year + 1
    return first_year


def run_spinup(
    simu,
    args: argparse.Namespace,
    t_atmbc: list,
    net_flux: np.ndarray,
    valid_times: np.ndarray,
    ds_LAI_raster: xr.Dataset,
    first_iteration: bool,
    outlet_node_id: int = None,
    md_log: "MarkdownLog | None" = None,
) -> tuple[bool, bool]:
    """
    Equilibrate the pressure-head field before the real simulation period
    ("spin-up antes de 2016"): since satellite ETp/rain forcing only
    starts 2016-10-01, there is no real forcing to run before 2016, so
    the spin-up instead cycles one full calendar year of forcing
    (default: the first complete year available) ``args.spinup_cycles``
    times, warm-starting psi between runs, to bring the model to a
    quasi-equilibrium state before the real simulation loop begins.

    Vegetation is held static during spin-up (the vegetation state of the
    first month of the real record) — only ATMBC forcing and psi evolve.
    Spin-up runs are not written to the output NetCDFs; they only leave
    behind the final psi field on disk, which the main loop's existing
    "warm start from simu.read_outputs('psi')" logic then picks up.

    Parameters
    ----------
    first_iteration : bool
        Whether the CATHY preprocessor has not yet been run. If True,
        the first spin-up run triggers it, and the returned value is
        False so the main loop does not re-run it.
    outlet_node_id : int, optional
        Node id to track psi at across spin-up cycles (see
        ``get_outlet_node``). If given, a diagnostic CSV
        ("spinup_outlet_psi_diagnostic.csv") and plot
        ("spinup_outlet_psi_diagnostic.png") are written at the end of
        spin-up, one line per cycle per month. A monotonically rising
        trend toward 0 across cycles (rather than flattening out) means
        spin-up is drifting toward saturation because of closed
        boundaries, not equilibrating — the diagnostic this issue asked
        for.

    Returns
    -------
    (first_iteration, spinup_done) : tuple[bool, bool]
        Updated ``first_iteration`` flag and whether at least one
        spin-up run completed (used by the main loop to decide whether
        to warm-start ICs from disk even before its own first run).
    """
    spinup_year = args.spinup_year or pick_default_spinup_year(valid_times)
    print(f"  Spin-up cycle year : {spinup_year}  |  cycles : {args.spinup_cycles}")
    if md_log is not None:
        md_log.section("Spin-up (warm-up)")
        md_log.meta({
            "Spin-up year": spinup_year,
            "Cycles": args.spinup_cycles,
            "Outlet diagnostic node": outlet_node_id if outlet_node_id is not None else "n/a",
        })

    # Static vegetation for spin-up: first available month of the real record.
    try:
        dt0, lai0 = next(iter(
            iter_monthly_lai(ds_LAI_raster.isel(time=slice(1, None)))
        ))
    except StopIteration:
        print("  WARNING: no LAI data available — skipping spin-up.")
        if md_log is not None:
            md_log.warn("No LAI data available — spin-up skipped entirely.")
        return first_iteration, False

    veg_map_remapped, remap = lai_to_veg_map(simu, lai0)
    veg_map = smooth_veg_map(veg_map_remapped)
    _, veg_valid = simu._check_outside_DEM(veg_map)
    simu.update_veg_map(veg_map, maxveg=N_VEG_CLASSES)
    active_zones = set(int(z) for z in veg_valid)
    _feddes = feddes_params_from_args(args)
    update_soil_for_active_zones(
        simu, active_zones, remap=remap,
        base_zroot=args.zroot, pmin=args.pmin,
        scf=args.scf,
        zroot_factors=zroot_factors_from_args(args),
        permx=args.permx,
        permy_ratio=args.permy_ratio,
        permz_ratio=args.permz_ratio,
        pcana=_feddes["pcana"], pcref=_feddes["pcref"], pcwlt=_feddes["pcwlt"],
        pz=_feddes["pz"], omgc=_feddes["omgc"],
        feddes_crop=_feddes["feddes_crop"],
    )
    if md_log is not None:
        md_log.bullet(f"Static spin-up vegetation from {dt0:%Y-%m}: "
                       f"zones {sorted(active_zones)}")
        md_log.table_header(["Cycle", "Month", "LAI mean", "LAI min–max",
                              "ψ outlet (m)", "Status"])

    spinup_ran = False
    mesh_refreshed = False     # see note at the first update_ic below
    outlet_psi_by_cycle = []   # diagnostic: psi at outlet node, per cycle/month
    for cycle in range(1, args.spinup_cycles + 1):
        for month in range(1, 13):
            t_month, nv_month = subset_atmbc_for_month(
                t_atmbc, net_flux, valid_times, spinup_year, month
            )
            if t_month is None:
                continue

            simu.update_atmbc(HSPATM=1, IETO=1, time=t_month, netValue=nv_month)

            df_atmbc_month = simu.read_inputs("atmbc")
            simu.update_parm(TIMPRTi=list(df_atmbc_month.time), IPRT=4, VTKF=2)

            #if first_iteration:
                #simu.run_preprocessor(verbose=True)
                #first_iteration = False

            if spinup_ran:
                # FIX: same stale-mesh issue the main loop already handles
                # (see "this preprocessor call can finalize a DIFFERENT
                # (larger) mesh" note there). update_ic -> map_prop2mesh
                # validates against the mesh cached by create_mesh_vtk()
                # in setup_dem, which predates the processor finalizing the
                # 3D grid. psi has one value per real 3D node, so without
                # this refresh we get "'ic' has length N but the mesh has M
                # points". Rebuild once, after the first spin-up run.
                if not mesh_refreshed:
                    try:
                        simu.create_mesh_vtk(verbose=False)
                        mesh_refreshed = True
                    except Exception as exc:
                        print(f"  WARNING: could not refresh mesh before "
                              f"spin-up IC warm start ({exc}).")
                        if md_log is not None:
                            md_log.warn(f"Mesh refresh failed in spin-up: {exc}")

                df_psi   = simu.read_outputs("psi").copy()
                psi_ini  = df_psi.iloc[-1].values
                simu.update_ic(INDP=1, pressure_head_ini=psi_ini)

            print(f"  [spin-up {cycle}/{args.spinup_cycles}] "
                  f"{spinup_year}-{month:02d} …")
            # MAXVEG already fixed to N_VEG_CLASSES before the loop; no need
            # to re-set it here every cycle (it was previously re-set to a
            # mismatched value of 10 on every iteration).

            try:
                simu.run_processor(
                    IPRT1=2,
                    DTMIN=args.dtmin,
                    DTMAX=args.dtmax,
                    DELTAT=args.deltat,
                    TRAFLAG=0,
                    verbose=False,
                )
            except Exception as exc:
                print(f"  WARNING: spin-up run failed for "
                      f"{spinup_year}-{month:02d} ({exc}) — stopping spin-up early.")
                if md_log is not None:
                    md_log.table_row([cycle, f"{spinup_year}-{month:02d}",
                                       f"{np.nanmean(lai0):.3f}",
                                       f"{np.nanmin(lai0):.3f}–{np.nanmax(lai0):.3f}",
                                       "n/a", f"✗ FAILED: {exc}"])
                break
            spinup_ran = True

            # Diagnostic: track psi at the outlet node this cycle/month.
            psi_outlet_str = "n/a"
            if outlet_node_id is not None:
                try:
                    df_psi_diag = simu.read_outputs("psi")
                    psi_outlet = float(df_psi_diag.iloc[-1][outlet_node_id])
                    outlet_psi_by_cycle.append(
                        {"cycle": cycle, "month": month, "psi_outlet": psi_outlet}
                    )
                    psi_outlet_str = f"{psi_outlet:.4f}"
                except Exception as exc:
                    print(f"  WARNING: could not read outlet psi for "
                          f"{spinup_year}-{month:02d} ({exc})")
                    psi_outlet_str = "read failed"

            if md_log is not None:
                md_log.table_row([cycle, f"{spinup_year}-{month:02d}",
                                   f"{np.nanmean(lai0):.3f}",
                                   f"{np.nanmin(lai0):.3f}–{np.nanmax(lai0):.3f}",
                                   psi_outlet_str, "✓ OK"])

    if outlet_psi_by_cycle:
        df_diag = pd.DataFrame(outlet_psi_by_cycle)
        df_diag.to_csv("spinup_outlet_psi_diagnostic.csv", index=False)

        fig, ax = plt.subplots(figsize=(8, 4))
        for c, grp in df_diag.groupby("cycle"):
            ax.plot(grp["month"], grp["psi_outlet"], marker="o", label=f"cycle {c}")
        ax.axhline(0, color="k", ls="--", lw=1)
        ax.set_xlabel("Month (within spin-up year)")
        ax.set_ylabel("\u03c8 outlet (m)")
        ax.legend()
        ax.set_title(f"Spin-up convergence check — outlet node {outlet_node_id}")
        fig.savefig("spinup_outlet_psi_diagnostic.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        # Cheap monotonic-drift check: compare each cycle's December value
        # (or last available month) to the previous cycle's.
        last_by_cycle = df_diag.sort_values("month").groupby("cycle")["psi_outlet"].last()
        if len(last_by_cycle) > 1:
            diffs = last_by_cycle.diff().dropna()
            if (diffs > 0).all():
                verdict = ("end-of-cycle outlet psi rises monotonically toward 0 "
                            "across every spin-up cycle "
                            f"({list(last_by_cycle.round(4))}) — this looks like "
                            "accumulation from a closed boundary, not equilibration. "
                            "Re-check the outlet BC.")
                print(f"  DIAGNOSTIC WARNING: {verdict}")
                if md_log is not None:
                    md_log.warn(verdict)
            else:
                verdict = ("end-of-cycle outlet psi across cycles: "
                            f"{list(last_by_cycle.round(4))} — not monotonically "
                            "rising, consistent with equilibration rather than "
                            "closed-boundary drift.")
                print(f"  DIAGNOSTIC: {verdict}")
                if md_log is not None:
                    md_log.bullet(f"✓ {verdict}")
        print("  Spin-up outlet-psi diagnostic saved → "
              "spinup_outlet_psi_diagnostic.{csv,png}")
        if md_log is not None:
            md_log.bullet("Diagnostic detail saved → "
                           "`spinup_outlet_psi_diagnostic.csv` / `.png`")

    if not spinup_ran:
        print(f"  WARNING: no forcing found for spin-up year {spinup_year} "
              f"— spin-up skipped.")
        if md_log is not None:
            md_log.warn(f"No forcing found for spin-up year {spinup_year} "
                         f"— spin-up skipped.")
        return first_iteration, False

    print(f"  Spin-up complete ({args.spinup_cycles} cycle(s) of {spinup_year}).")
    if md_log is not None:
        md_log.bullet(f"Spin-up complete: {args.spinup_cycles} cycle(s) of "
                       f"{spinup_year}.")
    return first_iteration, True
