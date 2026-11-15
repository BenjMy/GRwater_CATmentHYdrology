#!/usr/bin/env python3
"""
test_smoke_agramon.py - quick assertion-based smoke test of the pure Agramon helpers
====================================================================================
No pyCATHY, no data files, no simulation: runs in a few seconds and checks the
logic the paper results depend on (LAI -> vegetation class, SCF labels, ETa
artefact floor, scenario folder names, Ks conversions / recovery).

    python lib/test/test_smoke_agramon.py          # exits 0 and prints "smoke test OK"

`geoutils` (needed only to import agramon_runlog) is stubbed when not installed.
"""
import argparse
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd

LIB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LIB))
sys.path.insert(0, str(LIB / "utils"))

try:                                    # geoutils is a private dependency
    import geoutils  # noqa: F401
except ImportError:
    stub = types.ModuleType("geoutils")
    stub.uCATHY = types.SimpleNamespace()
    sys.modules["geoutils"] = stub

# heavy / private dependencies of the data-loading modules: stub whichever are missing so the
# pure helpers stay testable anywhere (a real install is used when present)
from unittest import mock  # noqa: E402

for _m in ("pyvista", "rioxarray", "geopandas", "scipy", "scipy.ndimage", "pyCATHY",
           "pyCATHY.meshtools", "pyCATHY.cathy_tools", "pyCATHY.cathy_utils", "pyCATHY.petro",
           "pyCATHY.plotters", "pyCATHY.importers", "rasterio", "shapely", "imageio"):
    try:
        __import__(_m)
    except ImportError:
        sys.modules[_m] = mock.MagicMock()

try:                                    # Agramon_utils pulls a long chain of optional packages
    import utils.Agramon_utils  # noqa: F401
except Exception:                       # noqa: BLE001  (any import-time failure)
    sys.modules["utils.Agramon_utils"] = mock.MagicMock()

import agramon_config as C                      # noqa: E402
from agramon_checks import apply_eta_artefact_floor   # noqa: E402
from agramon_lai import lai_to_veg_map          # noqa: E402
from agramon_runlog import build_scenario_dirname      # noqa: E402
import ks_fire_utils as K                       # noqa: E402
import run_paper_scenarios as RP                # noqa: E402
import agramon_ks as AK                         # noqa: E402
import agramon_spinup as SP                     # noqa: E402
import agramon_diagnostics as DG                # noqa: E402
import agramon_mesh as MS                       # noqa: E402
from agramon_forcing import subset_atmbc_for_month  # noqa: E402
from agramon_config import feddes_params_from_args, zroot_factors_from_args  # noqa: E402
from datetime import datetime                   # noqa: E402


def test_config():
    assert C.N_VEG_CLASSES == len(C.LAI_THRESHOLDS) + 1 == 5
    assert sorted(C.ZROOT_FACTORS) == list(range(1, C.N_VEG_CLASSES + 1))
    assert C.scf_type("0") == 0.0 and C.scf_type("1") == 1.0
    for bad in ("-0.1", "1.5"):
        try:
            C.scf_type(bad)
        except argparse.ArgumentTypeError:
            pass
        else:
            raise AssertionError(f"--scf {bad} must be rejected")
    assert "Feddes" in C.scf_label_for(1.0) and "PMIN" in C.scf_label_for(0.0)
    assert "blend" in C.scf_label_for(0.5)


def test_lai_to_veg_map():
    simu = types.SimpleNamespace(hapin={"N": 6, "M": 1})
    lai = np.array([[0.0, 1.0, 1.01, 1.5, 2.0, 2.6]])
    veg, remap = lai_to_veg_map(simu, lai)
    # thresholds are upper-inclusive: LAI == 1.0 stays in class 1
    assert veg.tolist() == [[1, 1, 2, 2, 3, 5]], veg
    assert remap == {c: c for c in range(1, 6)}, "remap must be the fixed identity"
    veg2, _ = lai_to_veg_map(simu, np.full((1, 6), 3.0))
    assert (veg2 == 5).all(), "all-dense field -> class 5, soil table size unchanged"


def _et_month(eta):
    return pd.DataFrame({"time_sec": [0.0, 86400.0], "SURFACE NODE": [1, 1],
                         "ACT. ETRA": eta})


def test_eta_artefact_floor():
    etp = np.full((2, 1), 4e-8)                 # m/s
    t = [0.0, 86400.0]
    idx = np.array([0, 1])
    # day 0: rain > ETp, ETa collapsed to 0 (artefact) -> floored at ETp
    rain = np.array([[8e-8], [0.0]])
    df, n_art, _ = apply_eta_artefact_floor(_et_month([0.0, 3e-8]), t, etp, rain, idx)
    assert n_art == 1
    assert np.isclose(df["ACT. ETRA_patched"].iloc[0], 4e-8), "artefact -> ETp"
    assert df["ACT. ETRA"].iloc[0] == 0.0, "raw column untouched"
    # day 1: no rain, mode 'etp' -> ETa = ETp; mode 'off' -> raw value kept
    assert np.isclose(df["ACT. ETRA_patched"].iloc[1], 4e-8)
    df_off, _, _ = apply_eta_artefact_floor(_et_month([0.0, 3e-8]), t, etp, rain, idx,
                                            net_neg_mode="off")
    assert np.isclose(df_off["ACT. ETRA_patched"].iloc[1], 3e-8)
    # 'gross': min(ETp, raw + rain); soil-limited cell stays below ETp
    df_g, _, _ = apply_eta_artefact_floor(_et_month([0.0, 1e-8]), t, etp, rain, idx,
                                          net_neg_mode="gross")
    assert np.isclose(df_g["ACT. ETRA_patched"].iloc[1], 1e-8)
    try:
        apply_eta_artefact_floor(_et_month([0.0, 0.0]), t, etp, rain, idx, net_neg_mode="bad")
    except ValueError:
        pass
    else:
        raise AssertionError("invalid net_neg_mode must raise")
    # empty month: returns a patched column without crashing
    none = pd.DataFrame({"time_sec": [], "SURFACE NODE": [], "ACT. ETRA": []})
    empty, a, b = apply_eta_artefact_floor(none, [], etp, rain, idx)
    assert len(empty) == 0 and (a, b) == (0, 0) and "ACT. ETRA_patched" in empty


def test_scenario_dirname():
    n = build_scenario_dirname(7, 2, 0.0, -1e-6, False, 3, -3.0, "wt")
    assert n == "scenario_7_plot2_outlet-closed_bottom-free_spinup-on3c_pmin-m3_ic-wt", n
    n2 = build_scenario_dirname(1, 3, -1e-6, 0.0, True, 0, -1e35, "pressure", "forcing-etpts-wet-fid3")
    assert n2.startswith("scenario_1_plot3_outlet-free_bottom-closed_spinup-off_")
    assert n2.endswith("_ic-pressure_forcing-etpts-wet-fid3"), n2
    assert build_scenario_dirname(1, 1, 0, 0, True, 0, -50, "wt") != \
        build_scenario_dirname(1, 1, 0, 0, True, 0, -50, "pressure"), "IC mode must differ"


def test_ks_utils():
    assert np.isclose(K.m_s_to_mm_hr(K.mm_hr_to_m_s(43.8)), 43.8)
    assert np.isclose(K.mm_hr_to_m_s(1.0), 1 / 3.6e6), "1 mm/hr = 2.778e-7 m/s"
    # logistic recovery is monotonic and bounded by the anchors
    v = [K.logistic_recovery(t, 16.2, 43.8, 1.0) for t in (-5, 0, 1, 5, 50)]
    assert all(a < b for a, b in zip(v, v[1:])) and 16.2 < v[0] and v[-1] <= 43.8 + 1e-9
    assert np.isclose(K.logistic_recovery(0, 16.2, 43.8, 1.0), (16.2 + 43.8) / 2)
    # LAI -> Ks linear map, clipped at the anchors
    f = lambda lai, clip=True: K.ks_from_lai_linear(lai, 0.0, 2.0, 16.2, 43.8, clip)
    assert np.isclose(f(0.0), 16.2) and np.isclose(f(2.0), 43.8) and np.isclose(f(1.0), 30.0)
    assert f(5.0) == 43.8 and f(-1.0) == 16.2 and f(5.0, clip=False) > 43.8
    try:
        K.ks_from_lai_linear(1.0, 1.0, 1.0, 1.0, 2.0)
    except ValueError:
        pass
    else:
        raise AssertionError("lai_ul == lai_ll must raise")


def test_paper_scenarios():
    sc = RP.load_config(RP.DEFAULT_CONFIG)
    assert len(sc) >= 10 and all(s.args for s in sc.values())
    good = RP.expand(sc["paper_good"])
    assert len(good) == 3, "paper_good must run all three Ks values (old shell recipe ran one)"
    assert [c[c.index("--permx") + 1] for c in good] == ["1e-4", "6e-4", "1e-3"]
    assert good[0][good[0].index("--pmin") + 1:good[0].index("--pmin") + 3] == ["-50", "-1e35"]
    z = RP.expand(sc["baseline_zroot_dem3"])
    assert len(z) == 2 and z[1][z[1].index("--zroot-factors") + 1:][:5] == ["2.5"] * 5
    assert RP._flag("force", True) == ["--force"] and RP._flag("force", False) == []
    assert [s.name for s in RP.select(sc, None, [])] == ["paper_good"], "default tag = paper"
    assert "paper_good" not in [s.name for s in RP.select(sc, ["all"], ["paper"])]
    try:
        RP.select(sc, ["nope"], [])
    except SystemExit:
        pass
    else:
        raise AssertionError("unknown scenario must stop")


def _ks_args(argv):
    p = argparse.ArgumentParser()
    K.add_ks_fire_args(p)
    p.set_defaults(ks_fire_date=AK.FIRE_DATE_DEFAULT)
    return p.parse_args(argv)


def test_agramon_ks():
    mm = K.mm_hr_to_m_s
    lai = np.array([[1.0, 2.0], [np.nan, 3.0]])
    off = AK.ks_config_from_args(_ks_args([]))
    assert AK.ks_for_month(off, datetime(2021, 1, 1), lai, 1e-4) == 1e-4, "off -> static permx"
    assert AK.ks_for_month(off, datetime(2021, 1, 1), lai, None) is None
    # time mode: unburned before the fire, burned minimum at the fire, recovering after
    a = _ks_args(["--ks-fire-mode", "time", "--ks-fire-preset", "overall_average"])
    cfg = AK.ks_config_from_args(a)
    pre = AK.ks_for_month(cfg, datetime(2019, 6, 1), lai, None)
    at = AK.ks_for_month(cfg, datetime(2020, 7, 1), lai, None)
    post = AK.ks_for_month(cfg, datetime(2023, 7, 1), lai, None)
    assert np.isclose(pre, mm(37.5)) and np.isclose(at, mm(13.0)), (pre, at)
    assert at < post <= pre + 1e-12, "Ks recovers after the fire, never above unburned"
    assert AK.ks_tag(a) == "ks-time-overall_average" and AK.ks_tag(_ks_args([])) == ""
    # lai mode: domain-mean LAI (NaN ignored) mapped linearly onto the anchors
    a = _ks_args(["--ks-fire-mode", "lai", "--ks-fire-ll", "1e-6", "--ks-fire-ul", "1e-5",
                  "--ks-fire-lai-ll", "0", "--ks-fire-lai-ul", "4"])
    cfg = AK.ks_config_from_args(a)
    assert np.isclose(AK.ks_for_month(cfg, datetime(2021, 1, 1), lai, None), 1e-6 + 9e-6 * 2.0 / 4)
    assert AK.ks_for_month(cfg, datetime(2019, 1, 1), lai, None) == 1e-5, "pre-fire -> unburned"
    assert "mode=lai" in AK.ks_banner_label(cfg)
    # run key: identifies the setup in the CSV log / shell skip-if-done check
    assert AK.ks_key(_ks_args([])) == "off" and AK.ks_key_from_argv([]) == "off"
    k1 = AK.ks_key_from_argv(["--ks-fire-mode=time", "--ks-fire-preset=median"])
    k2 = AK.ks_key_from_argv(["--ks-fire-mode=time", "--ks-fire-preset=overall_average"])
    k3 = AK.ks_key_from_argv(["--ks-fire-mode=time", "--ks-fire-preset=median",
                              "--ks-fire-date=2020-08-01"])
    assert len({k1, k2, k3, "off"}) == 4 and k1.startswith("time|preset=median|date=2020-07-01")
    assert AK.ks_key(_ks_args(["--ks-fire-mode", "time", "--ks-fire-preset", "median"])) == k1
    # incomplete configs fail at startup, not at the first month
    for bad in (["--ks-fire-mode", "time"], ["--ks-fire-mode", "lai", "--ks-fire-ul", "1e-5"]):
        try:
            AK.ks_config_from_args(_ks_args(bad))
        except ValueError:
            pass
        else:
            raise AssertionError(f"{bad} must raise")


def test_moved_helpers():
    """Functions moved out of the main script (REFACTORING 2.3-2.6) keep their behaviour."""
    # spin-up year: skip an incomplete first year (satellite forcing starts in October)
    v = pd.date_range("2016-10-01", "2019-12-31", freq="MS")
    assert SP.pick_default_spinup_year(v.values) == 2017
    assert callable(SP.run_spinup) and callable(DG.save_forcing_diagnostics)
    ns = argparse.Namespace(zroot_factors=[0.3, 0.6, 2.0, 3.0, 5.0])
    assert zroot_factors_from_args(ns) == {1: 0.3, 2: 0.6, 3: 2.0, 4: 3.0, 5: 5.0}
    fp = feddes_params_from_args(argparse.Namespace(
        feddes_crop=None, pcana=5.0, pcref=None, pcwlt=None, pz=None, omgc=None))
    assert fp["pcana"] == 5.0 and fp["pcref"] is None and fp["feddes_crop"] is None, fp
    assert callable(subset_atmbc_for_month)


def test_smooth_dem():
    if isinstance(sys.modules.get("scipy.ndimage"), mock.MagicMock):
        print("  (skipped: scipy not installed)")
        return
    dem = np.arange(36, dtype=float).reshape(6, 6)
    dem[0, :2] = np.nan
    dem[5, 4:] = np.nan
    assert MS.smooth_dem(dem, 0) is dem, "sigma <= 0 is a no-op"
    out = MS.smooth_dem(dem, 1.0)
    assert out.shape == dem.shape and np.array_equal(np.isnan(out), np.isnan(dem)), "mask unchanged"
    assert not np.array_equal(out[~np.isnan(out)], dem[~np.isnan(dem)]), "values do change"
    assert np.nanmin(dem) <= np.nanmin(out) and np.nanmax(out) <= np.nanmax(dem), "no overshoot"
    flat = np.full((4, 4), 7.0)
    flat[0, 0] = np.nan
    assert np.allclose(np.nan_to_num(MS.smooth_dem(flat, 2.0), nan=7.0), 7.0), "flat stays flat"


def test_forcing_mode():
    import warnings
    fm = C.ForcingMode.from_name("etp")
    assert (fm.use_etp, fm.use_etp_ts, fm.label) == (True, False, C.FORCING_LABELS["etp"])
    ts = C.ForcingMode.from_name("etp_ts")
    assert ts.use_etp and ts.use_etp_ts and ts.name == "etp_ts"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        era = C.ForcingMode.from_name("era5")
        assert not era.use_etp and not era.use_etp_ts
        assert any(issubclass(x.category, DeprecationWarning) for x in w), "era5 is deprecated"
        w.clear()
        C.ForcingMode.from_name("era5", warn=False)
        assert not w
    try:
        C.ForcingMode.from_name("nope")
    except ValueError:
        pass
    else:
        raise AssertionError("unknown forcing must raise")
    try:
        fm.name = "etp_ts"                       # immutable: no hidden global-style mutation
    except Exception:                            # noqa: BLE001 (FrozenInstanceError)
        pass
    else:
        raise AssertionError("ForcingMode must be frozen")
    assert C.FORCING_MODE_DEFAULT in C.FORCING_LABELS


def test_results_plotter_split():
    """results_plotter (REFACTORING 2.11) is a facade over plot_*.py: same names, same objects,
    and a figure that crosses two modules (plot_maps -> plot_common helpers) really renders."""
    import tempfile
    import matplotlib
    matplotlib.use("Agg")
    import results_plotter as RP
    import plot_common, plot_maps, plot_monthly, plot_hydro_xr, plot_eta
    for name, mod in (("plot_recharge_bar", plot_maps), ("plot_hydro_xr_monthly", plot_monthly),
                      ("plot_hydro_xr", plot_hydro_xr), ("plot_eta_total", plot_eta),
                      ("add_fire_vline", plot_common), ("FIRE_DATE", plot_common)):
        assert getattr(RP, name) is getattr(mod, name), name
    assert RP._infer_time_dim.__module__ == "plot_common", "last definition wins (as before the split)"
    assert RP.CMAP_ET == "YlGnBu" and RP.DEFAULT_ET_VAR == "ACT. ETRA_patched"
    idx = pd.date_range("2020-01-31", periods=12, freq="ME")   # month-end, as the callers pass
    with tempfile.TemporaryDirectory() as d:
        RP.plot_recharge_bar(pd.Series(np.arange(12.0), index=idx), Path(d), dpi=40)
        assert (Path(d) / "recharge_LAI_loop.png").stat().st_size > 1000


def test_dem_empty_polygon_regression():
    """Regression for the '--dem-plot 2 crashes' bug (REFACTORING.md): a polygon whose clip holds
    no valid DEM pixel used to die with 'zero-size array to reduction operation minimum';
    it must now raise a ValueError that names the cause. Needs real rasterio/rioxarray/geopandas."""
    import tempfile
    for m in ("rasterio", "rioxarray", "geopandas"):
        if isinstance(sys.modules.get(m), mock.MagicMock):
            print(f"  (skipped: {m} not installed)")
            return
    import geopandas as gpd
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import box

    x0, y0, res, n = 600000.0, 4400000.0, 5.0, 30
    dem = np.full((n, n), np.nan, dtype="float32")
    dem[:, :15] = 100 + np.arange(15)[None, :] + np.arange(n)[:, None] * 0.1   # data: left half
    with tempfile.TemporaryDirectory() as d:
        tif = Path(d) / "dem.tif"
        with rasterio.open(tif, "w", driver="GTiff", height=n, width=n, count=1, dtype="float32",
                           crs=C.TARGET_CRS, transform=from_origin(x0, y0, res, res),
                           nodata=np.nan) as dst:
            dst.write(dem, 1)
        xl = lambda a, b: box(x0 + a * res, y0 - 20 * res, x0 + b * res, y0 - 5 * res)
        gdf = gpd.GeoDataFrame({"fid_1": [1, 2, 3]}, crs=C.TARGET_CRS, geometry=[
            xl(2, 10),                                            # 1: over data
            xl(20, 28),                                           # 2: over nodata (the bug)
            box(x0 + 1e5, y0 + 1e5, x0 + 1e5 + 50, y0 + 1e5 + 50)])   # 3: outside the DEM
        r, masked, xll, yll, rx, ry = MS.load_dem_from_tif(tif, gdf_clip=gdf, fid=1)
        assert MS.count_valid_dem_pixels(masked) > 0 and (rx, ry) == (res, res)
        assert masked.shape[0] > 0 and masked.shape[1] > 0
        for bad in (2, 3):
            try:
                MS.load_dem_from_tif(tif, gdf_clip=gdf, fid=bad)
            except ValueError as e:
                assert "zero-size" not in str(e), f"fid {bad}: unclear error {e}"
            except Exception as e:                 # noqa: BLE001  (rioxarray NoDataInBounds)
                assert "zero-size" not in str(e), f"fid {bad}: unclear error {e}"
            else:
                raise AssertionError(f"fid {bad}: empty clip must raise")
        try:
            MS.load_dem_from_tif(tif, gdf_clip=gdf, fid=2)
        except ValueError as e:
            assert "no valid" in str(e) and "fid_1=2" in str(e), e
    assert MS.count_valid_dem_pixels(np.array([[-9999.0, 1.0], [np.nan, 2.0]])) == 2, "nodata counted"


def test_bc_and_forcing_containers():
    import types
    import numpy as np
    import agramon_forcing as F
    a = types.SimpleNamespace(outlet_flux=-1e-6, outlet_side=None, bottom_flux=0.0,
                              ic_mode="wt", ic_wtposition=2.5, ic_pressure_head=-10.0)
    bc = C.BCConfig.from_args(a)
    assert bc.ic_label == "WT depth 2.5 m (INDP=3)" and bc.bottom_label == "disabled (closed)"
    calls = []
    simu = types.SimpleNamespace(update_ic=lambda **k: calls.append(k))
    bc.apply_ic(simu)
    C.BCConfig(ic_mode="pressure", ic_pressure_head=-3.0).apply_ic(simu)
    assert calls == [dict(INDP=3, WTPOSITION=2.5), dict(INDP=0, pressure_head_ini=-3.0)], calls
    for bad in (dict(ic_mode="x"), dict(outlet_side="north"), dict(ic_wtposition=-1.0)):
        try:
            C.BCConfig(**bad)
        except ValueError:
            pass
        else:
            raise AssertionError(bad)
    T, S = 3, 4
    kw = dict(t_atmbc=[0, 86400, 172800], valid_times=np.arange(T), n_surf=S)
    f = F.Forcing(ETp_nodes=np.ones((T, S)), rain_nodes=np.zeros((T, S)),
                  net_flux=-np.ones(T), **kw)
    assert f.is_etp and f.n_steps == T
    assert not F.Forcing(None, None, np.zeros(T), t_atmbc=kw["t_atmbc"],
                         valid_times=kw["valid_times"]).is_etp
    for bad in (dict(ETp_nodes=np.ones((T, S + 1)), rain_nodes=np.zeros((T, S)), net_flux=np.zeros(T)),
                dict(ETp_nodes=np.ones((T, S)), rain_nodes=None, net_flux=np.zeros(T)),
                dict(ETp_nodes=np.ones((T, S)), rain_nodes=np.zeros((T, S)), net_flux=np.zeros(T + 1))):
        try:
            F.Forcing(**bad, **kw)
        except ValueError:
            pass
        else:
            raise AssertionError(bad)


def test_export_mesh_props(tmp_path=None):
    import os, tempfile, types, datetime as dtm
    import Agramon_withLAI_withETp as A
    d = tempfile.mkdtemp(); old = os.getcwd(); os.chdir(d)
    try:
        vtk = os.path.join(d, "proj", "vtk"); os.makedirs(vtk)
        open(os.path.join(vtk, "proj.vtk"), "w").write("x")
        calls = []
        ok = types.SimpleNamespace(workdir=d, project_name="proj",
                                   map_soil_root_to_mesh=lambda: calls.append(1))
        A.export_mesh_props(ok, dtm.datetime(2020, 7, 1))
        assert calls == [1] and os.path.exists("mesh_props/proj_202007.vtk")
        A._MESH_PROPS_WARNED = False
        old_py = types.SimpleNamespace(workdir=d, project_name="proj")
        A.export_mesh_props(old_py, dtm.datetime(2020, 8, 1))   # must not raise
        assert A._MESH_PROPS_WARNED and not os.path.exists("mesh_props/proj_202008.vtk")
    finally:
        os.chdir(old)


def test_run_pipeline_split():
    import inspect
    import Agramon_withLAI_withETp as A
    params = list(inspect.signature(A.run_monthly_loop).parameters)
    assert params == ["args", "inp", "ms", "atm", "ks_cfg", "ds_LAI_raster",
                      "first_iteration", "spinup_done"], params
    src = inspect.getsource(A.run_pipeline)
    assert [n for n in inspect.signature(A.configure_model).parameters] == ["args", "inp", "ms", "atm"]
    assert set(A.ModelSetup.__dataclass_fields__) == {"ds_LAI_raster", "outlet_node_id",
                                                      "first_iteration", "spinup_done"}
    assert "configure_model(" in src and "run_spinup(" not in src
    assert len(src.splitlines()) < 40, len(src.splitlines())
    assert "run_monthly_loop(" in src and "for dt, lai_2d in iter_monthly_lai" not in src
    assert "for dt, lai_2d in iter_monthly_lai" in inspect.getsource(A.run_monthly_loop)


def test_monthly_loop_helpers():
    import inspect, tempfile, types
    import numpy as np, pandas as pd
    from pathlib import Path
    import Agramon_withLAI_withETp as A
    assert list(inspect.signature(A.refresh_mesh_mapping).parameters) == [
        "simu", "fm", "ds_etp", "ds_rain", "raster_DEM_masked",
        "n_surf", "ds_mesh", "ETp_nodes", "rain_nodes"]
    src = inspect.getsource(A.run_monthly_loop)
    for call in ("write_etp_forcing(", "save_run_figures("):
        assert call in src, call
    assert "refresh_mesh_mapping(" in inspect.getsource(A.run_one_month)
    assert "simu.create_mesh_vtk(" not in src and "plot_run_summary(" not in src
    # write_etp_forcing: no-op without ETp, writes a (datetime, node) file otherwise
    with tempfile.TemporaryDirectory() as d:
        fm = types.SimpleNamespace(use_etp=True, use_etp_ts=True)
        A.write_etp_forcing(types.SimpleNamespace(use_etp=False, use_etp_ts=False), None, None, Path(d))
        assert not (Path(d) / "etp_output.nc").exists()
        times = pd.date_range("2020-01-01", periods=3)
        A.write_etp_forcing(fm, np.ones((3, 4)) * 1e-8, times, Path(d))
        assert (Path(d) / "etp_output.nc").exists()
        # save_run_figures must never raise, even with empty stats
        A.save_run_figures([], types.SimpleNamespace(), "x", Path(d))
        # refresh_mesh_mapping: when n_surf already matches, inputs are returned untouched
        simu = types.SimpleNamespace(create_mesh_vtk=lambda verbose=True: None,
                                     read_outputs=lambda k: {"nnod": 4})
        e, r = np.zeros((3, 4)), np.ones((3, 4))
        g, n, dm, e2, r2 = A.refresh_mesh_mapping(simu, fm, None, None, None, 4, "dsm", e, r)
        assert n == 4 and dm == "dsm" and e2 is e and r2 is r and g["nnod"] == 4


def test_config_env_paths():
    import importlib, os
    import agramon_config as C
    keys = ("AGRAMON_ROOT", "AGRAMON_PATH_ET", "AGRAMON_LAI_PATH", "AGRAMON_CODES")
    saved = {k: os.environ.pop(k, None) for k in keys}
    try:
        importlib.reload(C)
        default_et = C.PathET
        os.environ["AGRAMON_ROOT"] = "/data/nc"
        importlib.reload(C)
        assert str(C.PathET) == "/data/nc/GRwater/data/satellite/Agramon", C.PathET
        os.environ["AGRAMON_PATH_ET"] = "/data/et"
        importlib.reload(C)
        assert str(C.etp_path) == "/data/et/20161001_20241231_ET0.nc"
        assert str(C.rain_path).endswith("_TP.nc")
    finally:
        for k in keys:
            os.environ.pop(k, None)
            if saved[k] is not None:
                os.environ[k] = saved[k]
        importlib.reload(C)
    assert C.PathET == default_et


def test_run_one_month_split():
    import inspect
    import Agramon_withLAI_withETp as A
    assert list(inspect.signature(A.run_one_month).parameters) == [
        "args", "inp", "ms", "atm", "st", "ks_cfg", "dt", "lai_2d",
        "spinup_done", "out_dir", "monthly_dir"]
    a = A.LoopState(True, 4, None, None, None, None)
    b = A.LoopState(False, 5, None, None, None, None)
    a.psi_records.append(1)
    assert b.psi_records == [] and a.active_zones_prev is None   # no shared mutable defaults
    loop = inspect.getsource(A.run_monthly_loop)
    assert "run_one_month(" in loop and "run_processor" not in loop
    assert "run_processor" in inspect.getsource(A.run_one_month)
    assert len(loop.splitlines()) < 150, len(loop.splitlines())


def test_agramon_cli_parses():
    """The CLI must build (guards against duplicate option strings, which only fail at parse time)."""
    import sys
    import Agramon_withLAI_withETp as A
    old = sys.argv
    sys.argv = ["x", "--log-level", "DEBUG", "--log-to-file", "/tmp/x.log", "--log-file", "my.csv"]
    try:
        a = A._parse_args()
    finally:
        sys.argv = old
    assert a.log_file == "my.csv" and a.log_to_file == "/tmp/x.log" and a.log_level == "DEBUG"


def test_agramon_monthly_readers():
    import types, datetime as dtm
    import pandas as pd
    import agramon_monthly as MR
    dt = dtm.datetime(2020, 7, 1)

    class Boom:
        def read_outputs(self, k):
            raise RuntimeError("no file")

    recs, stats = [], {}
    assert MR.read_sw_month(Boom(), dt, recs, stats, None) is None and not recs
    assert MR.read_recharge_month(Boom(), dt, recs, stats, None) is None
    assert MR.read_et_month(Boom(), dt, recs, stats, None, None, use_etp=False,
                            t_month=None, ETp_nodes=None, rain_nodes=None,
                            idx_month=None, eta_net_neg_mode="x") is None
    MR.read_hgraph_month(Boom(), dt, recs, None)           # must not raise
    assert not recs and not stats
    ok = types.SimpleNamespace(
        read_outputs=lambda k: pd.DataFrame({"time": [0.0, 86400.0], "Q": [1.0, 2.0]}))
    hg = []
    MR.read_hgraph_month(ok, dt, hg, None)
    assert len(hg) == 1 and "time_sec" in hg[0].columns


def test_agramon_monthly_psi_and_stats():
    import tempfile, types, datetime as dtm
    from pathlib import Path
    import pandas as pd
    import agramon_monthly as MR
    dt = dtm.datetime(2020, 7, 1)

    class Log:
        def __init__(self): self.w, self.b = [], []
        def warn(self, m): self.w.append(m)
        def bullet(self, m): self.b.append(m)

    class Boom:
        def read_outputs(self, k): raise RuntimeError("no psi")

    log, recs, st = Log(), [], {"month": "2020-07"}
    df, status = MR.read_psi_month(Boom(), dt, recs, st, None, log, [86400.0])
    assert df is None and status.startswith("✗ PSI read failed") and log.w and not recs

    behind = types.SimpleNamespace(read_outputs=lambda k: pd.DataFrame(
        {0: [-1.0], 1: [-2.0]}, index=[100.0]))
    df, status = MR.read_psi_month(behind, dt, recs, st, None, log, [86400.0])
    assert df is None and status.startswith("⚠ not converged") and not recs

    out = Path(tempfile.mkdtemp()); stats = []
    MR.report_month_stats(dt, {"month": "2020-07", "psi_mean": -1.5}, "✓ OK", log, stats, out)
    assert len(stats) == 1 and (out / "monthly_stats.csv").exists()
    assert stats[0]["solver_status"] == "✓ OK" and any("psi mean" in b for b in log.b)


def test_agramon_monthly_prep_helpers():
    import os, tempfile, datetime as dtm
    import numpy as np
    import agramon_monthly as MR
    vt = np.array(["2020-06-30", "2020-07-01", "2020-07-02", "2020-08-01"], dtype="datetime64[D]")
    idx = MR.month_forcing_indices(vt, dtm.datetime(2020, 7, 1))
    assert list(idx) == [1, 2] and len(MR.month_forcing_indices(vt, dtm.datetime(2021, 1, 1))) == 0
    old = os.getcwd(); os.chdir(tempfile.mkdtemp())
    try:
        MR.plot_month_forcing(np.linspace(-1e-8, 0, 5), dtm.datetime(2020, 7, 1), True)
        assert os.path.exists("figures_ETp_temp/forcing_202007.png")
    finally:
        os.chdir(old)


def test_agramon_log():
    import contextlib, io, logging, tempfile
    from pathlib import Path
    import agramon_log as L
    log = L.get_logger("agramon_test")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        L.setup_logging("INFO")
        log.info("  hello")
        log.warning("  Warning: x")
        log.debug("hidden")
    assert buf.getvalue() == "  hello\n  Warning: x\n", repr(buf.getvalue())   # same as print
    buf = io.StringIO()
    f = Path(tempfile.mkdtemp()) / "sub" / "run.log"
    with contextlib.redirect_stdout(buf):
        L.setup_logging("WARNING", str(f)); L.setup_logging("WARNING", str(f))  # idempotent
        log.info("quiet"); log.error("boom")
    assert buf.getvalue() == "boom\n" and "ERROR" in f.read_text() and "quiet" not in f.read_text()
    assert len([h for h in logging.getLogger("agramon").handlers
                if isinstance(h, logging.FileHandler)]) == 1
    try:
        L.setup_logging("LOUD")
    except ValueError:
        pass
    else:
        raise AssertionError("bad level accepted")
    L.setup_logging("INFO")
    # no bare print() left in the converted modules except the 3 multi-arg debug dumps
    import ast
    for mod in ("agramon_forcing", "agramon_lai", "agramon_mesh", "agramon_monthly",
                "agramon_runlog", "agramon_spinup"):
        tree = ast.parse(open(Path(__file__).resolve().parents[1] / f"{mod}.py").read())
        assert not [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                    and getattr(n.func, "id", "") == "print"], mod


def test_recovery_metrics():
    import numpy as np
    import pandas as pd
    import recovery_metrics as R
    idx = pd.date_range("2020-01-31", periods=12, freq="ME")
    ctrl = pd.Series(1.0, index=idx)
    burn = pd.Series([1.0] * 6 + [0.4, 0.5, 0.7, 0.92, 0.95, 1.0], index=idx)   # fire 2020-07
    r = R.relative_to_control(burn, ctrl)
    assert np.isclose(r.iloc[6], 0.4)
    # recovered (|r-1|<=0.1) from Oct 2020 for >=3 months -> 3 months after July
    assert R.recovery_time(r, "2020-07-01", tol=0.1, sustain=3) == 3.0
    assert np.isnan(R.recovery_time(r, "2020-07-01", tol=0.01, sustain=3))   # never within 1 %
    assert np.isnan(R.recovery_time(r, "2020-07-01", tol=0.1, sustain=4))    # run too short
    assert np.isclose(R.deficit(burn, ctrl, "2020-07-01"), 0.6 + 0.5 + 0.3 + 0.08 + 0.05)
    df_b = pd.DataFrame({"LAI": burn, "ETa": burn * 2}); df_c = pd.DataFrame({"LAI": ctrl, "ETa": ctrl * 2})
    tab = R.recovery_table(df_b, df_c, "2020-07-01")
    assert list(tab.index) == ["LAI", "ETa"] and np.isclose(tab.loc["LAI", "min_ratio"], 0.4)
    assert np.isclose(tab.loc["ETa", "deficit"], 2 * tab.loc["LAI", "deficit"])
    z = R.relative_to_control(burn, ctrl * 0)
    assert z.isna().all()
    for bad in (lambda: R.recovery_time(r, "2020-07-01", sustain=0),
                lambda: R.relative_to_control(burn.values, ctrl)):
        try:
            bad()
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError("bad input accepted")


def test_plot_hypotheses():
    import os, tempfile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import plot_hypotheses as H
    d = tempfile.mkdtemp()
    # skill: perfect, biased, degenerate
    s = pd.Series(np.arange(10.0)); sk = H.skill_metrics(s, s)
    assert sk["rmse"] == 0 and all(np.isclose(sk[k], 1) for k in ("r", "nse", "kge")), sk
    sk = H.skill_metrics(s + 1, s)
    assert np.isclose(sk["bias"], 1) and np.isclose(sk["rmse"], 1) and sk["r"] > 0.999
    assert np.isnan(H.skill_metrics(pd.Series([1.0]), pd.Series([1.0]))["rmse"])
    assert np.isnan(H.skill_metrics(pd.Series([1.0, 1.0, 1.0]), pd.Series([1.0, 2.0, 3.0]))["r"])
    # run summary (+ missing columns tolerated, solver flag marker)
    st = pd.DataFrame({"month": ["2020-01", "2020-02", "2020-03"], "eta_mean": [1.0, 2.0, 1.5],
                       "psi_mean": [-3, -4, -5], "solver_status": ["✓ OK", "⚠ not converged (1 < 2)", "✓ OK"]})
    f = H.plot_run_summary(st, "2020-02-01", out=os.path.join(d, "rs.png"))
    assert os.path.getsize(os.path.join(d, "rs.png")) > 1000 and len(f.axes) == 2
    st["permx_mm_hr"] = [36.0, 20.0, 30.0]       # applied Ks -> extra log-scale panel on top
    f = H.plot_run_summary(st, "2020-02-01")
    assert len(f.axes) == 3 and f.axes[0].get_yscale() == "log" and f.axes[0].get_ylabel().startswith("Ks")
    for bad in (pd.DataFrame({"x": [1]}), pd.DataFrame({"month": ["2020-01"]})):
        try:
            H.plot_run_summary(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("bad stats accepted")
    # validation panel
    idx = pd.date_range("2020-01-01", periods=60, freq="D")
    obs = pd.Series(np.sin(np.arange(60) / 8) + 2, index=idx)
    f = H.plot_validation_panel(obs * 1.1 + 0.2, obs, out=os.path.join(d, "val.png"))
    assert len(f.axes) == 3
    try:
        H.plot_validation_panel(obs, obs.shift(1000).dropna())
    except ValueError:
        pass
    else:
        raise AssertionError("empty overlap accepted")
    # water balance
    mi = pd.date_range("2020-01-31", periods=6, freq="ME")
    wb = pd.DataFrame({"precip": [50, 40, 30, 10, 5, 0.0], "eta": [20, 25, 30, 20, 10, 5.0],
                       "recharge": [5, 5, 2, 0, 0, 0.0]}, index=mi)
    f = H.plot_water_balance(wb, out=os.path.join(d, "wb.png")); assert len(f.axes) == 1
    # recovery panels
    ic = pd.date_range("2020-01-31", periods=12, freq="ME")
    ctrl = pd.DataFrame({"LAI": 1.0, "ETa": 2.0}, index=ic)
    burn = ctrl * pd.Series([1.0] * 6 + [0.4, 0.5, 0.7, 0.92, 0.95, 1.0], index=ic).to_numpy()[:, None]
    f = H.plot_recovery_panels(burn, ctrl, "2020-07-01", out=os.path.join(d, "rec.png"))
    assert len(f.axes) == 4 and "recovered after 3 months" in f.axes[1].get_title()
    plt.close("all")


def test_check_monthly_sanity():
    import datetime as dtm
    import numpy as np
    import pandas as pd
    import agramon_checks as C
    dt = dtm.datetime(2020, 7, 1)
    t = [0.0, 86400.0, 2 * 86400.0]; idx = np.array([0, 1, 2])
    n_surf = 2
    etp = np.full((3, n_surf), 5e-8); rain = np.full((3, n_surf), 1e-8)
    sw = pd.DataFrame({"sw": [0.2, 0.5, 1.2]})
    psi = pd.DataFrame({"psi": [-1.0, -50.0, -2000.0]})
    et = pd.DataFrame({"time_sec": [0.0, 0.0, 86400.0, 86400.0],
                       "SURFACE NODE": [1, 2, 1, 2],
                       "ACT. ETRA": [1e-8, 1e-8, 9e-8, 1e-8]})            # one value > ETp
    rec = pd.DataFrame({"time_sec": [0.0, 86400.0], "X": [0.0, 0.0], "rech": [1e-9, 1e-9]})
    st = {}
    w = C.check_monthly_sanity(dt, st, psi, sw, et, rec, t, idx, etp, rain)
    assert st["sw_out_of_bounds_n"] == 1 and st["psi_extreme_n"] == 1 and st["eta_exceeds_etp_n"] == 1
    assert len([x for x in w if "SW value" in x]) == 1 and any("exceed ETp" in x for x in w)
    # rain 3 d * 1e-8 m/s * 86400 * 1000 = 2.592 mm; ETa mean per output time; residual is rain - ETa - recharge
    assert np.isclose(st["rain_in_mm"], 3 * 1e-8 * 86400 * 1e3)
    assert np.isclose(st["eta_out_mm"], (1e-8 + 5e-8) * 86400 * 1e3)
    assert np.isclose(st["water_balance_residual_mm"],
                      st["rain_in_mm"] - st["eta_out_mm"] - st["recharge_out_mm"])
    # clean month: no warnings, no keys for missing inputs, never raises on junk
    st2 = {}
    ok = C.check_monthly_sanity(dt, st2, pd.DataFrame({"psi": [-1.0]}), pd.DataFrame({"sw": [0.5]}),
                                None, None, None, None, None, None)
    assert ok == [] and st2["sw_out_of_bounds_n"] == 0 and "rain_in_mm" not in st2
    bad_et = pd.DataFrame({"time_sec": [0.0], "SURFACE NODE": [99], "ACT. ETRA": [1.0]})
    C.check_monthly_sanity(dt, {}, None, None, bad_et, None, t, idx, etp, rain)   # must not raise
    # big residual -> water-balance warning
    st3 = {}
    w3 = C.check_monthly_sanity(dt, st3, None, None, et.assign(**{"ACT. ETRA": 0.0}), rec, t, idx,
                                etp, rain * 100, wb_residual_frac=0.5)
    assert any("water balance residual" in x for x in w3)


def test_plot_ks_scenarios():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import plot_ks_scenarios as K
    months = pd.date_range("2019-07-01", "2022-07-01", freq="MS")
    df = K.ks_trajectories(months, "2020-07-01", presets=["median", "lower_bound"])
    assert set(df["mode"]) == set(K.MODES) and set(df["preset"]) == {"median", "lower_bound"}
    pre = df[df["date"] < "2020-07-01"]
    assert np.allclose(pre["ratio"], 1.0)                              # unburned before the fire
    med = df[(df["preset"] == "median")]
    t0 = med[(med["date"] == "2020-07-01")].set_index("mode")["ks_mm_hr"]
    assert np.isclose(t0["time"], 12.9) and np.isclose(t0["burn_severity"], 19.5 * (1 - 0.5 * 0.8))
    late = med[(med["date"] == "2022-07-01")].set_index("mode")["ratio"]
    assert late["time"] > 0.99 and late["lai"] > 0.99 and np.isclose(late["burn_severity"], 0.6)
    try:
        K.ks_trajectories(months, "2020-07-01", severity=2)
    except ValueError:
        pass
    else:
        raise AssertionError("bad severity accepted")
    fig = K.plot_ks_scenarios(df, "2020-07-01")
    assert len(fig.axes) == 2
    plt.close("all")


def test_cmp_eta_variants():
    import sys, tempfile
    from pathlib import Path
    import matplotlib
    matplotlib.use("Agg")
    import pandas as pd
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "synthetic_LAI_Ks"))
    import synthetic_ksfire_analysis as SA
    idx = pd.date_range("2020-01-01", periods=6, freq="MS")
    mon = pd.DataFrame({"etp_mm": 100.0, "eta_mm": 100.0, "eta_gross_mm": [60, 50, 40, 30, 30, 40.0],
                        "eta_raw_mm": [10.0] * 6}, index=idx)
    runs = {"a": {"mon": mon}, "b": {"mon": mon.drop(columns="eta_gross_mm")}}
    out = Path(tempfile.mkdtemp())
    SA.cmp_eta_variants(runs, {"a": {}, "b": {}}, "2020-03-01", out, 80)
    assert (out / "fig_cmp_eta_variants.png").stat().st_size > 1000
    try:
        SA.cmp_eta_variants({"b": runs["b"]}, {"b": {}}, None, out, 80)
    except ValueError:
        pass
    else:
        raise AssertionError("runs without eta_gross_mm accepted")
    assert SA.diagnose_eta_vs_etp({"a": {"mon": mon}}) == ["a"]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok  {name}")
    print("smoke test OK")
