#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meshing from a Digital Elevation Model (DEM)
============================================

GRwater project 

*Estimated time to run the notebook = 5min*

"""
import os
import numpy as np
import pandas as pd
from pathlib import Path
import rioxarray as rxr
import xarray as xr
import matplotlib.pyplot as plt
import sys
# Add the path where Agramon_utils.py is located
module_path = Path('../../GRwater_geophy').resolve()  # replace with actual path
if str(module_path) not in sys.path:
    sys.path.append(str(module_path))
import Agramon_utils as AgUtils
from pyCATHY import cathy_tools
from pathlib import Path
import argparse
from datetime import datetime
import matplotlib.dates as mdates

import geoutils
from geoutils import uCATHY

#%%

# -----------------------------
# ARGPARSE
# -----------------------------
parser = argparse.ArgumentParser(description="Load ERA5 series and run simulation")
parser.add_argument("--start_year", type=int, default=2011, help="First year to include")
parser.add_argument("--end_year", type=int, default=2012, help="Last year (exclusive)")
parser.add_argument("--ZROOT", type=float, default=1.0, help="Root depth in meters")
parser.add_argument("--PMIN", type=float, default=-5, help="Minimum soil moisture or threshold")
parser.add_argument("--log_file", type=str, default="simulation_log.csv", help="Path to CSV log file")
args = parser.parse_args()

# -----------------------------
# Load ERA5
# -----------------------------
rootPath = Path('..').resolve()
ds = AgUtils.load_era5_series(
    root_path=rootPath,
    start_year=args.start_year,
    end_year=args.end_year
)
# ds = ds.sel(valid_time=slice(0, 60))
# ds = ds.sel(valid_time=slice(None, "2012-02-01"))
# ds.valid_time

if ds is not None:
    print(f"Loaded dataset from {args.start_year} to {args.end_year-1}")
else:
    print("No datasets loaded.")

# -----------------------------
# Log simulation
# -----------------------------
params = {
    "start_year": args.start_year,
    "end_year": args.end_year,
    "ZROOT": args.ZROOT,
    "PMIN": args.PMIN
}

sim_index = uCATHY.log_simulation(args.log_file, params)

# s

#%%
gdf_Agramon = AgUtils.load_plot_shapefiles(module_path/'shapefiles')
daily_ts = AgUtils.extract_point_timeseries(ds, gdf_Agramon)

#%%
# Convert to pandas Series
pev_series = daily_ts['pev'].to_series()
tp_series = daily_ts['tp'].to_series()

# Plot
fig, ax = plt.subplots(figsize=(15, 4))
ax.bar(pev_series.index, pev_series.values, width=1.0, color='skyblue')
ax.bar(tp_series.index, tp_series.values, width=1.0, color='orange')

# Format x-axis
locator = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)
ax.tick_params(axis='x', rotation=45)

ax.set_title('PEV over time')
ax.set_ylabel('mm/day or units')
ax.set_xlabel('Date')
ax.legend()
plt.tight_layout()
# plt.show()

#%% Init CATHY model
# ------------------------
prj_name = f"sc_{args.start_year}-{args.end_year}_rd{args.ZROOT}_pmin{args.PMIN}_sim{sim_index}"
path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
                         prj_name=prj_name
                         )
print(f"Project initialized with name: {prj_name}")


# Path to the directory containing the .adf file (not the file itself)
adf_folder = "../prepro/DTMsPlota/dtmplot1/"

# Open the raster (typically named 'hdr.adf', but you only need the folder)
raster_DEM = rxr.open_rasterio(adf_folder, masked=True).isel(band=0)

# Create a mask of valid (non-NaN) data
valid_mask = ~np.isnan(raster_DEM)

# Apply the mask to get valid coordinates
valid_x = raster_DEM['x'].where(valid_mask.any(dim='y'), drop=True)
valid_y = raster_DEM['y'].where(valid_mask.any(dim='x'), drop=True)

# Get min and max valid coordinates
min_lon, max_lon = float(valid_x.min()), float(valid_x.max())
min_lat, max_lat = float(valid_y.min()), float(valid_y.max())

raster_DEM_masked = raster_DEM.where(
    (raster_DEM['x'] >= min_lon) & (raster_DEM['x'] <= max_lon), drop=True
).where(
    (raster_DEM['y'] >= min_lat) & (raster_DEM['y'] <= max_lat), drop=True
)
       
raster_DEM_masked = np.where(np.isnan(raster_DEM_masked), -9999, raster_DEM_masked)
np.shape(raster_DEM_masked)
np.shape(raster_DEM)


#%% Fetch and show initial DEM

fig, ax = plt.subplots(1)
img = ax.imshow(raster_DEM_masked)
plt.colorbar(img)

simu.show_input(prop="dem")

simu.update_prepo_inputs(
    DEM=raster_DEM_masked,
    delta_x=5,
    delta_y=5,
    ivert=1,
)

fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)
simu.create_mesh_vtk(verbose=True)

#%%

simu.run_preprocessor(verbose=False)
# simu.run_processor(IPRT1=3,verbose=True)

# # simu.read_inputs('atmbc')
# simu.update_parm(TIMPRTi=[1800,7200])

grid3d = simu.read_outputs('grid3d')

# selec = [0,360]
t0 = daily_ts['pev'].valid_time[0].values  # first timestamp
t_seconds = (daily_ts['pev'].valid_time.values - t0) / np.timedelta64(1, "s")
t_atmbc = [float(ti) for ti in t_seconds]
netValue = np.array(daily_ts['tp']) + np.array(daily_ts['pev'])
# min(netValue)
# max(netValue)

#% Create an empty dataframe of SPP and set default SPP properties 
df_SPP_map = simu.init_soil_SPP_map_df(nzones=1,nstr=15)
SPP_map = simu.set_SOIL_defaults(SPP_map_default=True)

simu.update_veg_map()

df_FP_map = simu.init_soil_FP_map_df(nveg=1)
df_FP_map = simu.set_SOIL_defaults(FP_map_default=True)
df_FP_map['ZROOT']=args.ZROOT

# SPP_map['PERMX'] = 6.88e-4
# SPP_map['PERMY'] = 6.88e-4
# SPP_map['PERMZ'] = 6.88e-4

#% Update soil file
simu.update_soil(PMIN=args.PMIN,
                 SPP_map=SPP_map)


simu.update_ic(INDP=0,
               IPOND=0,
               pressure_head_ini=-15
                )


simu.update_atmbc(                                                                                                                                                                                                                                                                                                      
                    HSPATM=1,
                    IETO=1,
                    time=t_atmbc,
                    netValue=netValue
                  )

#%%
simu.update_parm(
                        IPRT=4,
                        VTKF=0, # dont write vtk files
                        # TIMPRTi=[t_atmbc[0],t_atmbc[-1]]
                        # TIMPRTi=[t_atmbc[0],t_atmbc[-1]]
                        TIMPRTi=t_atmbc
                        )

#%%
# ss
simu.create_mesh_bounds_df(
                            'nansfdirbc',
                            simu.grid3d["mesh3d_nodes"],
                            t_atmbc,
                            )
print(simu.mesh_bound_cond_df.head())
print(simu.mesh_bound_cond_df.columns)

simu.update_nansfneubc(no_flow=True)
simu.update_nansfdirbc(no_flow=True)
simu.update_sfbc(no_flow=True)
# simu.show_bc(time=0)
# meshbc = simu.mesh_bound_cond_df
# cplt.plot_mesh_bounds('nansfdirbc',meshbc, time=0)

#%%
# int(grid3d['nnod'])

simu.project_name
simu.workdir
simu.run_processor(IPRT1=2,
                    DTMIN=1e-2,
                    DTMAX=1e4,
                    DELTAT=1e3,
                    TRAFLAG=0,
                    verbose=True
                    )

