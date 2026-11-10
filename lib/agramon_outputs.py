"""
agramon_outputs.py — persistence of the running combined output files.

Moved verbatim from the nested ``save_combined_outputs`` closure of
``run_pipeline`` in ``Agramon_withLAI_withETp.py``. The closure read ten
names from the enclosing scope; they are now explicit arguments. Each of
them was bound exactly once in ``run_pipeline`` (the lists/dicts are only
mutated in place), so the behaviour is unchanged.
"""

import pickle

import numpy as np
import pandas as pd
import xarray as xr

from agramon_config import TARGET_CRS


def save_combined_outputs(
    out_dir,
    veg_map_history,
    lai_history,
    et_records,
    psi_records,
    sw_records,
    all_recharge_recs,
    all_hgraph_recs,
    M,
    N,
) -> None:
    """
    (Re)write the single running output files in out_dir from whatever
    has been accumulated in the *_records lists / veg_map_history so
    far. Safe to call after every month: each call overwrites the
    previous version of these files with the latest cumulative data,
    so — unlike the per-month files in monthly_dir — there's always
    exactly one psi_output.nc / sw_output.nc / et_output.nc /
    recharge_output.nc / veg_map_history.nc reflecting all months
    completed up to now, without waiting for the loop to finish.
    """
    if veg_map_history:
        veg_history_path = out_dir / "veg_map_history.pkl"
        with open(veg_history_path, "wb") as fh:
            pickle.dump(veg_map_history, fh)

        sorted_dates = sorted(veg_map_history)
        veg_arr = np.stack([veg_map_history[d] for d in sorted_dates], axis=0)
        veg_xr = xr.DataArray(
            veg_arr,
            dims=["time", "row", "col"],
            coords={"time": sorted_dates, "row": np.arange(M), "col": np.arange(N)},
            name="veg_map",
            attrs={
                "description": "Vegetation zone map (1=bare, 2=sparse, 3=moderate, 4=dense, 5=very dense)",
                "crs": TARGET_CRS,
            },
        )
        veg_xr.to_netcdf(out_dir / "veg_map_history.nc")

    if lai_history:
        sorted_lai_dates = sorted(lai_history)
        lai_arr = np.stack([lai_history[d] for d in sorted_lai_dates], axis=0)
        lai_xr = xr.DataArray(
            lai_arr,
            dims=["time", "row", "col"],
            coords={"time": sorted_lai_dates, "row": np.arange(M), "col": np.arange(N)},
            name="LAI",
            attrs={
                "description": "Monthly LAI raster used to derive the vegetation map that month.",
                "crs": TARGET_CRS,
            },
        )
        lai_xr.to_netcdf(out_dir / "lai_output.nc")

    if et_records:
        ET_xr_all = xr.concat(et_records, dim="datetime")
        ET_xr_all.time.attrs.pop("dtype", None)
        ET_xr_all.to_netcdf(out_dir / "et_output.nc",
                             encoding={"time": {"dtype": "float64"}})

    if psi_records:
        df_psi_all = pd.concat(psi_records, ignore_index=True)
        df_psi_all = (
            df_psi_all
            .groupby(["datetime", "node"], as_index=False)
            .last()
        )
        df_psi_all.set_index(["datetime", "node"]).to_xarray().to_netcdf(
            out_dir / "psi_output.nc"
        )

    if sw_records:
        df_sw_all = pd.concat(sw_records, ignore_index=True)
        df_sw_all = (
            df_sw_all
            .groupby(["datetime", "node"], as_index=False)
            .last()
        )
        df_sw_all.set_index(["datetime", "node"]).to_xarray().to_netcdf(
            out_dir / "sw_output.nc"
        )

    if all_recharge_recs:
        xr.concat(all_recharge_recs, dim="datetime").to_netcdf(
            out_dir / "recharge_output.nc")

    if all_hgraph_recs:
        pd.concat(all_hgraph_recs, ignore_index=True).to_csv(
            out_dir / "hgraph_output.csv", index=False)
