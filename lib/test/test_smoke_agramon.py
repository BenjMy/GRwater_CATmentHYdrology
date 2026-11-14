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

import agramon_config as C                      # noqa: E402
from agramon_checks import apply_eta_artefact_floor   # noqa: E402
from agramon_lai import lai_to_veg_map          # noqa: E402
from agramon_runlog import build_scenario_dirname      # noqa: E402
import ks_fire_utils as K                       # noqa: E402


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


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok  {name}")
    print("smoke test OK")
