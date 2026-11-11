#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meshing from a Digital Elevation Model (DEM)
============================================

GRwater project 

*Estimated time to run the notebook = 5min*

"""
import os
import matplotlib.pyplot as plt
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

# Now import the module
import Agramon_utils as AgUtils

import pyCATHY.meshtools as mt
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as cplt
from pyCATHY import cathy_utils
import matplotlib.dates as mdates
import pyvista as pv
# %%
crs = 'EPGS:25830'
rootPath = Path('..').resolve()
scenario = 2
df_log = pd.read_csv('simulation_log.csv')
path2prj = Path("../SSHydro/").resolve()
prj_name = [p for p in path2prj.iterdir() if p.is_dir() and f"_sim{scenario}" in p.name][0]

ds = AgUtils.load_era5_series(
    root_path=rootPath,
    start_year=df_log[df_log['sim_index']==scenario]['start_year'].values[0],
    end_year=df_log[df_log['sim_index']==scenario]['end_year'].values[0],
)
gdf_Agramon = AgUtils.load_plot_shapefiles(module_path/'shapefiles')
daily_ts = AgUtils.extract_point_timeseries(ds, gdf_Agramon)
start_date = pd.to_datetime(daily_ts['tp'].valid_time[0].values)

#%% Init CATHY model
# ------------------------
path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
                         prj_name=prj_name
                         )
# --- Convert ERA5 to pandas Series ---
pev_series = daily_ts['pev'].to_series()
tp_series = daily_ts['tp'].to_series()

# --- Convert CATHY outputs ---
df_psi = simu.read_outputs('psi')
df_sw, _ = simu.read_outputs('sw')
df_dates = cathy_utils.change_x2date(df_psi.index, start_date)

df_ET = simu.read_outputs('ET')     
# --- Prepare ET data ---
ET_xr = df_ET.set_index(['X','Y','time']).to_xarray()
ET_xr = ET_xr.rio.set_spatial_dims('X','Y')
ET_xr_mean = ET_xr.mean(dim=['X','Y'])
ET_xr_mean = ET_xr_mean * (1e3 * 86400)  # to mm/day
ET_xr['time'] = start_date + ET_xr['time'].values
# ETtime_dates = start_date + ET_xr['time'].values

#%%
import centum 
from centum import plotting as pltC
from centum import accounting as acc

ET_xr_WA = acc.compute_water_accounting(ET_xr, variable='ACT. ETRA', freq='ME')

nrows = 4
ncols = 3
fig, axs = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), constrained_layout=True)


norm =  pltC.plot_monthly_volume_mm(ET_xr_WA, 
                                    variable='volume_mm',
                                    cmap='viridis',
                                    axs=axs
                                    )
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
fig.colorbar(sm, ax=axs[1], 
             orientation='vertical',
             fraction=0.02, pad=0.04, 
             label='volume_mm')
# fig.savefig(figurePath / 'map_net_irr_mm.png', dpi=300)

