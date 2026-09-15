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
from pyCATHY import cathy_utils

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
import matplotlib.dates as mdates

import pyvista as pv

crs = 'EPGS:25830'
# %%
rootPath = Path('..').resolve()
scenario = 2

#%%
df_log = pd.read_csv('simulation_log.csv')
# list file from path2prj, and find the one with patern '_sim{scenario}'

path2prj = Path("../SSHydro/").resolve()

# List folders containing the scenario pattern
prj_name = [p for p in path2prj.iterdir() if p.is_dir() and f"_sim{scenario}" in p.name][0]

ds = AgUtils.load_era5_series(
    root_path=rootPath,
    start_year=df_log[df_log['sim_index']==scenario]['start_year'].values[0],
    end_year=df_log[df_log['sim_index']==scenario]['end_year'].values[0],
)
gdf_Agramon = AgUtils.load_plot_shapefiles(module_path/'shapefiles')
daily_ts = AgUtils.extract_point_timeseries(ds, gdf_Agramon)

#%% Init CATHY model
# ------------------------

# --- Initialize CATHY project ---
simu = cathy_tools.CATHY(dirName=path2prj, prj_name=prj_name)

#%%
# --- Uphill nodes ---
uphill_surface_2d_node, uphill_surface_2d_pos = simu.find_nearest_node([250., 100.])
uphill_surface_node, uphill_surface_pos = simu.find_nearest_node([250., 100., uphill_surface_2d_pos[0][2]-0.2])
uphill_1m_node, uphill_1m_pos = simu.find_nearest_node([250., 100., uphill_surface_2d_pos[0][2]-1])
uphill_2m_node, uphill_2m_pos = simu.find_nearest_node([250., 100., uphill_surface_2d_pos[0][2]-2])

print("Uphill surface node position:", uphill_surface_pos[0])

# --- Downhill nodes ---
downhill_surface_2d_node, downhill_surface_2d_pos = simu.find_nearest_node([25., 50.])
downhill_surface_node, downhill_surface_pos = simu.find_nearest_node([25., 50., downhill_surface_2d_pos[0][2]-0.2])
downhill_1m_node, downhill_1m_pos = simu.find_nearest_node([25., 50., downhill_surface_2d_pos[0][2]-1])

# --- Access mesh3d nodes ---
mesh3d_nodes = simu.read_outputs('grid3d')['mesh3d_nodes']
print("First mesh3d node:", mesh3d_nodes[0])

# --- Plot using PyVista ---
pl = pv.Plotter(notebook=False)
mesh = pv.read(os.path.join(simu.workdir, simu.project_name, 'vtk') + '/' + prj_name.name + '.vtk')
pl.add_mesh(mesh,style='wireframe')
# cplt.show_vtk(
#     filename=prj_name.name + '.vtk',
#     unit="ic",
#     # timeStep=1,
#     path=os.path.join(simu.workdir, simu.project_name, 'vtk'),
#     style='wireframe',
#     opacity=0.1,
#     ax=pl,
# )
pl.add_bounding_box()
pl.add_axes()
# Add points for selected nodes (red)
for pos in [uphill_surface_pos[0], uphill_1m_pos[0], downhill_surface_pos[0], downhill_1m_pos[0]]:
    pl.add_points(pos, color='red')

pl.show()




#%%


# --- Convert ERA5 to pandas Series ---
pev_series = daily_ts['pev'].to_series()
tp_series = daily_ts['tp'].to_series()

# --- Convert CATHY outputs ---
df_psi = simu.read_outputs('psi')
df_sw, _ = simu.read_outputs('sw')

#%%
start_date = pd.to_datetime(daily_ts['tp'].valid_time[0].values)
df_dates = cathy_utils.change_x2date(df_psi.index, start_date)

# --- Setup subplots ---
fig, axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True,
                         gridspec_kw={"height_ratios": [1, 2, 2]})

# -------------------
# 1) Rain + ETp
# -------------------
ax0 = axes[0]
pev_dates = pev_series[pev_series.index<df_dates.max()].index
pev_values = pev_series[pev_series.index<df_dates.max()].values
tp_series_values = tp_series[pev_series.index<df_dates.max()].values

ax0.bar(pev_dates, 
        pev_values, width=1.0,
        color='grey', label="ETp", alpha=0.5)
ax0.bar(pev_dates, tp_series_values, width=1.0,
        color='skyblue', label="Rain")

# Invert y-axis so rain “falls down”
ax0.invert_yaxis()

ax0.set_ylabel("mm/day")
ax0.legend(loc="upper right")

# Colorblind-friendly palette
colors = {
    "uphill_surface": "#0072B2",  # blue
    "uphill_1m": "#56B4E9",       # light blue
    "downhill_surface": "#D55E00",# red-orange
    "downhill_1m": "#E69F00"      # orange
}

# -------------------
# 2) psi
# -------------------
ax1 = axes[1]

ax1.plot(df_dates, df_psi[uphill_surface_node].values,
         label="Uphill surface", color=colors["uphill_surface"],
         marker='+', linestyle='-')
ax1.plot(df_dates, df_psi[uphill_1m_node].values,
         label="Uphill -1m", color=colors["uphill_1m"], linestyle='--')

ax1.plot(df_dates, df_psi[downhill_surface_node].values,
         label="Downhill surface", color=colors["downhill_surface"],
         marker='o', linestyle='-')
ax1.plot(df_dates, df_psi[downhill_1m_node].values,
         label="Downhill -1m", color=colors["downhill_1m"], linestyle='--')

ax1.set_ylabel("ψ (m)")
ax1.legend(loc="upper right")
ax1.grid(True, linestyle='--', alpha=0.5)

# -------------------
# 3) sw
# -------------------
ax2 = axes[2]

ax2.plot(df_dates, df_sw[uphill_surface_node].values,
         label="Uphill surface", color=colors["uphill_surface"],
         marker='+', linestyle='-')
ax2.plot(df_dates, df_sw[uphill_1m_node].values,
         label="Uphill -1m", color=colors["uphill_1m"], linestyle='--')

ax2.plot(df_dates, df_sw[downhill_surface_node].values,
         label="Downhill surface", color=colors["downhill_surface"],
         marker='o', linestyle='-')
ax2.plot(df_dates, df_sw[downhill_1m_node].values,
         label="Downhill -1m", color=colors["downhill_1m"], linestyle='--')

ax2.set_ylabel("sw (-)")
ax2.set_xlabel("Date")
ax2.legend(loc="upper right")
ax2.grid(True, linestyle='--', alpha=0.5)


# -------------------
# Format x-axis (shared)
# -------------------
locator = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)

for ax in axes:
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.grid(True, which="both", linestyle="--", alpha=0.5)

fig.autofmt_xdate()
plt.tight_layout()
plt.show()



#%%

df_ET = simu.read_outputs('ET')     
# --- Prepare ET data ---
ET_xr = df_ET.set_index(['X','Y','time']).to_xarray()
ET_xr = ET_xr.rio.set_spatial_dims('X','Y')
ET_xr_mean = ET_xr.mean(dim=['X','Y'])
ET_xr_mean = ET_xr_mean * (1e3 * 86400)  # to mm/day


# --- Prepare time axis ---
start_date = pd.to_datetime(daily_ts['tp'].valid_time[0].values)
df_dates = cathy_utils.change_x2date(df_sw.index, start_date)
time_dates = start_date + ET_xr['time'].values


# Colorblind-friendly palette
colors = {
    "uphill_surface": "#0072B2",   # blue
    "uphill_1m": "#56B4E9",        # light blue
    "downhill_surface": "#D55E00", # red-orange
    "downhill_1m": "#E69F00",      # orange
    "ET_mean": "#009E73"           # green
}

# --- Create subplots ---
fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

# (1) Soil water content (sw)
axes[0].plot(df_dates, df_sw[uphill_surface_node].values,
             label="Uphill surface", color=colors["uphill_surface"], marker='+', linestyle='-')
axes[0].plot(df_dates, df_sw[uphill_1m_node].values,
             label="Uphill -1m", color=colors["uphill_1m"], linestyle='--')

axes[0].plot(df_dates, df_sw[downhill_surface_node].values,
             label="Downhill surface", color=colors["downhill_surface"], marker='o', linestyle='-')
axes[0].plot(df_dates, df_sw[downhill_1m_node].values,
             label="Downhill -1m", color=colors["downhill_1m"], linestyle='--')

axes[0].set_ylabel("sw (-)")
axes[0].legend(loc="upper right")
axes[0].grid(True, linestyle='--', alpha=0.5)

# (2) Actual ET
axes[1].scatter(time_dates, ET_xr_mean['ACT. ETRA'], c=colors["ET_mean"], s=12, label="ETa mean")
axes[1].set_ylabel("ETa (mm/day)")
axes[1].set_xlabel("Date")
axes[1].legend(loc="upper right")
axes[1].grid(True, linestyle='--', alpha=0.5)

# --- Format x-axis with auto dates ---
locator = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)
for ax in axes:
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

plt.tight_layout()
plt.show()


#%%
# df_hgraph = simu.read_outputs('hgraph')



# simu.show(prop="hgsfdet")

#%%

# simu.show(prop="dtcoupling", yprop="Atmpot-d")

#%%

# simu.show(prop="hgraph")

#%%

# simu.show(prop="cumflowvol")
