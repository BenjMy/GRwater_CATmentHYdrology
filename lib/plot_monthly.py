"""
plot_monthly.py - Monthly-resampled variants of the xarray hydrology figures.
===========================================================================
Split unchanged out of results_plotter.py (REFACTORING.md 2.11); results_plotter re-exports every name.
"""
from __future__ import annotations


from plot_common import MONTHLY_BAR_WIDTH, MONTHLY_RULE, _resample_monthly
from plot_eta import plot_catchment_eta_etp_ratio, plot_lai_et_sw_psi_ratio_xr
from plot_hydro_xr import plot_et_psi_lai_xr, plot_et_psi_xr, plot_et_sw_lai_xr, plot_et_sw_xr, plot_hydro_xr


def plot_hydro_xr_monthly(
    psi_xr, sw_xr, nodes_dict, fig_path, *,
    pev_series=None, tp_series=None, date_start=None, date_end=None,
    show_rain=True, show_psi=True, show_sw=True,
    show_uh_s=True, show_uh_1=True, show_dh_s=True, show_dh_1=True,
    colors_dict=None, scenario_id="?", dpi=150, show=False,
    fname="fig6f_hydro_timeseries_monthly.png",
) -> None:
    """Fig 6f -- monthly version of Fig 6 (psi/sw mean; rain/ETp monthly sum)."""
    plot_hydro_xr(
        _resample_monthly(psi_xr), _resample_monthly(sw_xr), nodes_dict, fig_path,
        pev_series=_resample_monthly(pev_series, "sum"),
        tp_series=_resample_monthly(tp_series, "sum"),
        date_start=date_start, date_end=date_end,
        show_rain=show_rain, show_psi=show_psi, show_sw=show_sw,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        bar_width=MONTHLY_BAR_WIDTH, rain_unit="mm/month",
    )


def plot_et_sw_xr_monthly(
    sw_xr, ET_xr_all, et_var, et_scale, nodes_dict, fig_path, *,
    date_start=None, date_end=None,
    show_uh_s=True, show_uh_1=True, show_dh_s=True, show_dh_1=True,
    colors_dict=None, scenario_id="?", dpi=150, show=False,
    fname="fig7j_et_sw_timeseries_monthly.png",
) -> None:
    """Fig 7j -- monthly version of Fig 7 (sw mean; ETa monthly sum)."""
    plot_et_sw_xr(
        _resample_monthly(sw_xr), _resample_monthly(ET_xr_all, "sum"),
        et_var, et_scale, nodes_dict, fig_path,
        date_start=date_start, date_end=date_end,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        et_unit="mm/month", et_as_line=True,
    )


def plot_et_psi_xr_monthly(
    psi_xr, ET_xr_all, et_var, et_scale, nodes_dict, fig_path, *,
    date_start=None, date_end=None,
    show_uh_s=True, show_uh_1=True, show_dh_s=True, show_dh_1=True,
    colors_dict=None, scenario_id="?", dpi=150, show=False,
    fname="fig7gj_et_psi_timeseries_monthly.png",
) -> None:
    """Fig 7gj -- monthly version of Fig 7g (psi mean; ETa monthly sum)."""
    plot_et_psi_xr(
        _resample_monthly(psi_xr), _resample_monthly(ET_xr_all, "sum"),
        et_var, et_scale, nodes_dict, fig_path,
        date_start=date_start, date_end=date_end,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        et_unit="mm/month", et_as_line=True,
    )


def plot_et_sw_lai_xr_monthly(
    sw_xr, ET_xr_all, et_var, et_scale, lai_ds, nodes_dict, fig_path, *,
    lai_var="LAI", date_start=None, date_end=None,
    show_uh_s=True, show_uh_1=True, show_dh_s=True, show_dh_1=True,
    colors_dict=None, scenario_id="?", dpi=150, show=False,
    fname="fig7bj_et_sw_lai_timeseries_monthly.png",
) -> None:
    """Fig 7bj -- monthly version of Fig 7b (sw, LAI mean; ETa monthly sum)."""
    plot_et_sw_lai_xr(
        _resample_monthly(sw_xr), _resample_monthly(ET_xr_all, "sum"),
        et_var, et_scale, _resample_monthly(lai_ds), nodes_dict, fig_path,
        lai_var=lai_var, date_start=date_start, date_end=date_end,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        et_unit="mm/month", et_as_line=True,
    )


def plot_et_psi_lai_xr_monthly(
    psi_xr, ET_xr_all, et_var, et_scale, lai_ds, nodes_dict, fig_path, *,
    lai_var="LAI", date_start=None, date_end=None,
    show_uh_s=True, show_uh_1=True, show_dh_s=True, show_dh_1=True,
    colors_dict=None, scenario_id="?", dpi=150, show=False,
    fname="fig7fj_et_psi_lai_timeseries_monthly.png",
) -> None:
    """Fig 7fj -- monthly version of Fig 7f (psi, LAI mean; ETa monthly sum)."""
    plot_et_psi_lai_xr(
        _resample_monthly(psi_xr), _resample_monthly(ET_xr_all, "sum"),
        et_var, et_scale, _resample_monthly(lai_ds), nodes_dict, fig_path,
        lai_var=lai_var, date_start=date_start, date_end=date_end,
        show_uh_s=show_uh_s, show_uh_1=show_uh_1,
        show_dh_s=show_dh_s, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname,
        et_unit="mm/month", et_as_line=True,
    )


def plot_lai_et_sw_psi_ratio_xr_monthly(
    lai_ds, ET_xr_all, et_var, et_scale, ETp_xr, etp_var, etp_scale,
    sw_xr, psi_xr, nodes_dict, grid3d, fig_path, *,
    lai_var="LAI", date_start=None, date_end=None,
    show_uh=True, show_dh=True, show_uh_1=True, show_dh_1=True,
    colors_dict=None, scenario_id="?", dpi=150, show=False,
    fname="fig8j_lai_et_sw_psi_ratio_timeseries_monthly.png",
) -> None:
    """Fig 7hj -- monthly version of Fig 7h. LAI/sw/psi: monthly mean;
    ETa and ETp: monthly sum, so the ratio panel is sum(ETa)/sum(ETp)."""
    plot_lai_et_sw_psi_ratio_xr(
        _resample_monthly(lai_ds), _resample_monthly(ET_xr_all, "sum"),
        et_var, et_scale,
        _resample_monthly(ETp_xr, "sum"), etp_var, etp_scale,
        _resample_monthly(sw_xr), _resample_monthly(psi_xr),
        nodes_dict, grid3d, fig_path,
        lai_var=lai_var, date_start=date_start, date_end=date_end,
        show_uh=show_uh, show_dh=show_dh,
        show_uh_1=show_uh_1, show_dh_1=show_dh_1,
        colors_dict=colors_dict, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname, et_unit="mm/month",
    )


def plot_catchment_eta_etp_ratio_monthly(
    ET_xr_all, et_var, et_scale, etp_mean_series, etp_median_series, fig_path, *,
    date_start=None, date_end=None, scenario_id="?", dpi=150, show=False,
    fname="fig_catchment_eta_etp_ratio_monthly.png",
) -> None:
    """Fig 7ij -- monthly version of the catchment ETa/ETp ratio."""
    plot_catchment_eta_etp_ratio(
        ET_xr_all, et_var, et_scale, etp_mean_series, etp_median_series, fig_path,
        date_start=date_start, date_end=date_end, scenario_id=scenario_id,
        dpi=dpi, show=show, fname=fname, freq=MONTHLY_RULE,
    )
