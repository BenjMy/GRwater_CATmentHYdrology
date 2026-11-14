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


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok  {name}")
    print("smoke test OK")
