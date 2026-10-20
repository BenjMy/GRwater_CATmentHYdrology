import numpy as np, pandas as pd
from ucathy_perf import *

def frames(n=50, artifact=False):
    dt = np.full(n, 3600.0); dt[0] = 0.01
    mbe = pd.DataFrame({"TIME": np.cumsum(dt), "NLIN": np.r_[np.full(n-1, 4), 12],
                        "REL. MBE": np.full(n, 0.01)})
    dtc = pd.DataFrame({"Deltat": dt, "Back": np.r_[1, np.zeros(n-1)],
                        "Atmact-vf": np.full(n, -1e-6)})
    t = np.repeat(np.arange(4) * 86400.0, 3)
    ev = np.full(12, 1e-9); ev[:3] = 0.0; tr = np.full(12, 2e-9)
    if artifact: ev[:3] = 1e-7
    et = pd.DataFrame({"time_sec": t, "ACT. EVAP": ev, "ACT. TRAN": tr, "ACT. ETRA": ev + tr})
    return dict(mbeconv=mbe, dtcoupling=dtc, ET=et)

cfg = PerfConfig(dtmin=1e-2, dtmax=3600.0)
r = evaluate(FrameSource(**frames()), "clean", cfg, {"cpu_s": 1.5})
m = r.metrics
assert m["n_steps"] == 50 and m["back_steps"] == 1 and m["nlin_max"] == 12
assert abs(m["frac_at_dtmin"] - 1/50) < 1e-12 and abs(m["frac_at_dtmax"] - 49/50) < 1e-12
assert m["et_t0_evap"] < 1e-12 and m["et_peak_at_t0"] is False and m["et_etra_sum_err"] < 1e-15
assert any("nlin_max" in w for w in r.warnings) and "cpu_s" in m
assert "mbe_cum_abs" not in m          # no MBE column -> silently skipped
a = evaluate(FrameSource(**frames(artifact=True)), "art", cfg)
assert any("t=0" in w for w in a.warnings) and any("et_t0_evap" in w for w in a.warnings)
class Broken:
    def read_outputs(self, n): raise IOError("nope")
b = evaluate(Broken(), "broken", cfg)
assert b.metrics == {} and len(b.notes) == 3
print("smoke tests OK; clean warnings:", r.warnings)

# --- aggregate_months ---
d = pd.DataFrame({"n_steps": [100, 300], "cpu_s": [1.0, 3.0], "back_steps": [2, 0],
                  "dt_min": [0.5, 0.1], "dt_max": [10., 20.], "dt_mean": [2., 6.],
                  "nlin_mean": [4., 8.], "nlin_max": [7, 11], "rel_mbe_max_pct": [.1, .3],
                  "cum_flux_m3_backward": [5., 7.], "warnings": ["", "nlin_max = 11 exceeds 10"],
                  "notes": ["", ""]})
g = aggregate_months(d)
assert g["n_steps"] == 400 and g["cpu_s"] == 4.0 and g["dt_min"] == 0.1 and g["nlin_max"] == 11
assert abs(g["nlin_mean"] - 7.0) < 1e-12 and abs(g["dt_mean"] - 5.0) < 1e-12   # step-weighted
assert abs(g["back_per_1000"] - 5.0) < 1e-12 and g["cum_flux_m3_backward"] == 12.0
assert g["n_warn_segments"] == 1 and g["n_note_segments"] == 0 and "lin_mean" not in g
assert aggregate_months(pd.DataFrame({"warnings": [""], "notes": ["x"]}))["n_note_segments"] == 1
print("aggregate tests OK")
