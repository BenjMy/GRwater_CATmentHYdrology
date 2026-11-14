"""
results_plotter.py — Shared Plotting Module
============================================
GRwater project — Agramon study site

Centralises every figure-generation function used by both:
  - test_withLAI_results.py   (LAI loop post-processing)
  - LT_SShydro_results.py     (long-term subsurface hydrology post-processing)

All functions are pure: they receive data, produce a matplotlib Figure,
save it to *fig_path*, and optionally display it.  No global state is
mutated.  Callers are responsible for loading artefacts and passing them in.

Public API
----------
Shared helpers
    save_fig(fig, fig_path, name, dpi)
    maybe_show(fig, show)
    time_label(t_val)
    et_global_limits(et_dataset, var)
    veg_colorbar(ax, cmap, labels)
    date_window(date_start, date_end, df_dates)
    add_fire_vline(ax, d0, d1, label, color)   # marks FIRE_DATE = 2020-07-01
    plot_catchment_et_freq(ET_xr, ET_xr_mean, time_dates, freq, fig_path, dpi, show, et_var)

LAI-loop figures (test_withLAI_results)
    plot_recharge_bar(recharge_monthly, fig_path, dpi, show)
    plot_et_monthly_mean(ET_xr_all, et_var, et_scale, fig_path, dpi, show)
    plot_et_by_veg_type(ET_xr_all, veg_map_history, et_var, et_scale, fig_path, dpi, show)
    plot_veg_map_samples(veg_map_history, cmap_veg, fig_path, dpi, show)
    plot_veg_map_year(veg_map_history, selected_year, cmap_veg, fig_path, dpi, show)
    plot_lai_map_year(ds_lai, var_name, selected_year, cmap, vmin, vmax, fig_path, dpi, show)
    make_veg_gif(veg_history, cmap_veg, out_path, duration)
    make_et_gif(et_dataset, var, scale, out_path, cmap, duration)

LAI diagnostic figures (analyse_LAI / withLAI_results)
    plot_lai_monthly_by_year(ds, lai_var, fig_path, dpi, show, fname)
    plot_lai_spatial_maps_year(ds, lai_var, selected_year, fig_path, cmap, vmin, vmax, dpi, show, fname)
    plot_lai_annual_hist(ds, lai_var, fig_path, bins, dpi, show, fname)
    plot_lai_annual_mean(ds, lai_var, fig_path, label, dpi, show, fname)
    plot_lai_annual_median(ds, lai_var, fig_path, label, dpi, show, fname)

SSHydro figures (LT_SShydro_results)
    plot_hydro(data, nodes_dict, args, fig_path)
    plot_et_sw(data, nodes_dict, args, fig_path, colors, et_var)
    plot_era5(data, args, fig_path)
    plot_catchment_et(data, fig_path, dpi, show, et_var)
    plot_hgraph(data, fig_path, dpi, show, hgraph_fname)
    plot_hgraph_csv(df_hgraph, fig_path, dpi, show, scenario_id, date_start, date_end, hgraph_fname)
    plot_node_locations(grid3d, nodes_dict, fig_path, dem_da, colors_dict, scenario_id, dpi, show)
    plot_recharge(data, fig_path, dpi, show, scenario_id)
    plot_recharge_xr(xr_recharge, fig_path, dpi, show, scenario_id)

LAI-loop xarray-based hydro figures (test_withLAI_results — Fig 6 / 7)
    plot_hydro_xr(psi_xr, sw_xr, pev_series, tp_series, nodes_dict,
                  fig_path, date_start, date_end, ...)
    plot_et_sw_xr(sw_xr, ET_xr_all, et_var, et_scale, nodes_dict,
                  fig_path, date_start, date_end, ...)
    plot_et_sw_lai_xr(sw_xr, ET_xr_all, et_var, et_scale, lai_ds, nodes_dict,
                  fig_path, lai_var, date_start, date_end, ...)   # Fig 7 + LAI panel

    # Fire-window / weekly variants (thin wrappers, same node/colour kwargs)
    plot_hydro_xr_fire_year(...)    # Fig 6c — fire year + following year
    plot_hydro_xr_fire_month(...)   # Fig 6d — fire month only
    plot_hydro_xr_weekly(...)       # Fig 6e — weekly-resampled
    plot_et_sw_xr_fire_year(...)    # Fig 7c — fire year + following year
    plot_et_sw_xr_fire_month(...)   # Fig 7d — fire month only
    plot_et_sw_xr_weekly(...)       # Fig 7e — weekly-resampled
    plot_et_sw_lai_xr_fire_year(...)    # Fig 7bc — fire year + following year
    plot_et_sw_lai_xr_fire_month(...)   # Fig 7bd — fire month only
    plot_et_sw_lai_xr_weekly(...)       # Fig 7be — weekly-resampled

LAI + ET + sw + psi + ETa/ETp ratio, at reference nodes (new)
    plot_lai_et_sw_psi_ratio_xr(lai_ds, ET_xr_all, et_var, et_scale,
                  ETp_xr, etp_var, etp_scale, sw_xr, psi_xr,
                  nodes_dict, grid3d, fig_path, ...)   # Fig 8 — 5-panel stack
"""
# Split by figure family (REFACTORING.md 2.11). Every name that used to live in this
# module is re-exported here, so `import results_plotter as plotter` keeps working.
from __future__ import annotations

from plot_common import (  # noqa: F401
    CMAP_ET,
    CMAP_LAI,
    CMAP_VEG,
    COLORS_HYDRO,
    DEFAULT_ET_VAR,
    FIRE_DATE,
    LAI_UNITS,
    LAI_VMAX,
    MONTHLY_BAR_WIDTH,
    MONTHLY_RULE,
    MONTH_ABBR,
    _da,
    _fire_month_window,
    _fire_year_window,
    _fmt_node_id,
    _infer_time_dim,
    _node_id,
    _node_ids,
    _resample_monthly,
    _spatial_dims,
    _values_at_nodes,
    _xr_time_to_dates,
    add_fire_vline,
    date_window,
    et_global_limits,
    maybe_show,
    save_fig,
    time_label,
    unwrap_node_id,
    veg_colorbar,
)
from plot_maps import (  # noqa: F401
    lai_history_to_dataset,
    make_et_gif,
    make_veg_gif,
    plot_et_by_veg_type,
    plot_et_monthly_mean,
    plot_lai_annual_hist,
    plot_lai_annual_mean,
    plot_lai_annual_median,
    plot_lai_map_year,
    plot_lai_monthly_by_year,
    plot_lai_spatial_maps_year,
    plot_recharge_bar,
    plot_veg_map_samples,
    plot_veg_map_year,
)
from plot_hydro_legacy import (  # noqa: F401
    plot_catchment_et,
    plot_catchment_et_freq,
    plot_era5,
    plot_et_sw,
    plot_hgraph,
    plot_hgraph_csv,
    plot_hydro,
    plot_node_locations,
    plot_recharge,
    plot_recharge_xr,
)
from plot_hydro_xr import (  # noqa: F401
    find_representative_mid_node,
    plot_et_psi_lai_xr,
    plot_et_psi_lai_xr_weekly,
    plot_et_psi_xr,
    plot_et_psi_xr_fire_year,
    plot_et_sw_lai_xr,
    plot_et_sw_lai_xr_fire_month,
    plot_et_sw_lai_xr_fire_year,
    plot_et_sw_lai_xr_weekly,
    plot_et_sw_xr,
    plot_et_sw_xr_fire_month,
    plot_et_sw_xr_fire_year,
    plot_et_sw_xr_weekly,
    plot_hydro_xr,
    plot_hydro_xr_fire_month,
    plot_hydro_xr_fire_year,
    plot_hydro_xr_weekly,
)
from plot_eta import (  # noqa: F401
    _plot_eta_total_freq,
    plot_catchment_eta_etp_ratio,
    plot_eta_annual_hist,
    plot_eta_monthly_by_year,
    plot_eta_spatial_maps_year,
    plot_eta_total,
    plot_lai_et_sw_psi_ratio_xr,
)
from plot_monthly import (  # noqa: F401
    plot_catchment_eta_etp_ratio_monthly,
    plot_et_psi_lai_xr_monthly,
    plot_et_psi_xr_monthly,
    plot_et_sw_lai_xr_monthly,
    plot_et_sw_xr_monthly,
    plot_hydro_xr_monthly,
    plot_lai_et_sw_psi_ratio_xr_monthly,
)
