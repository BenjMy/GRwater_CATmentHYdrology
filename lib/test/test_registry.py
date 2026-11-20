"""Smoke test of agramon_registry on a toy run log (no CATHY needed).

    python -I lib/test/test_registry.py
"""
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import agramon_registry as reg  # noqa: E402

BASE = dict(start_year=2016, end_year=2026, ZROOT=1.0, PMIN=-50.0, SCF=1.0, watershed_nb=3,
            withLAI=1, outlet_flux=0.0, bottom_flux=0.0, spinup="off", ic_mode="pressure",
            ic_wtposition=1.0, ic_pressure_head=-10.0, zroot_factors="0.05,0.25,1,1.5,2",
            forcing_mode="etp_ts", etp_ts_column="ET_0_wet", etp_ts_feature_id=2.0)


def _log(path, idx, params):
    """Append a row the way log_simulation does (before the registry columns exist)."""
    df = pd.read_csv(path) if path.exists() else pd.DataFrame()
    pd.concat([df, pd.DataFrame([{"sim_index": idx, **params}])], ignore_index=True).to_csv(path, index=False)


def _finished_dir(prj, idx):
    d = prj / "outputs" / f"scenario_{idx}_tag"
    d.mkdir(parents=True)
    for f in reg.REQUIRED:
        (d / f).write_text("x")
    return d


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        prj, log = tmp / "prj", tmp / "simulation_log.csv"

        # hash: 2 == 2.0 == "2"; empty legacy fields read as defaults; etp_ts_* ignored off etp_ts
        assert reg.param_hash({**BASE, "etp_ts_feature_id": 2}) == reg.param_hash(BASE)
        assert reg.param_hash({**BASE, "PMIN": -15.0}) != reg.param_hash(BASE)
        etp = {**BASE, "forcing_mode": "etp"}
        assert reg.param_hash({**etp, "etp_ts_column": "x"}) == reg.param_hash({**etp, "forcing_mode": ""})
        assert reg.param_hash({**BASE, "ks_fire": "off"}) == reg.param_hash(BASE)

        # run 0: legacy row (no registry columns), finished folder on disk
        _log(log, 0, BASE)
        _finished_dir(prj, 0)
        assert reg.find_run(log, BASE, prj) == 0
        assert reg.find_run(log, BASE, prj, require_outputs=False) == 0
        assert reg.load_registry(log).loc[0, "status"] == "unknown"

        # run 1: different params, started then failed -> never matches
        p1 = {**BASE, "PMIN": -15.0}
        _log(log, 1, p1)
        assert reg.start_run(log, 1, p1)
        assert reg.find_run(log, p1, prj, require_outputs=False) is None       # running
        assert reg.finish_run(log, 1, "failed", seconds=1.26, error_tail="Boom: x\n  at y")
        df = reg.load_registry(log)
        assert df.loc[1, "status"] == "failed" and df.loc[1, "seconds"] == 1.3
        assert df.loc[1, "error_tail"] == "Boom: x | at y"
        assert reg.find_run(log, p1, prj, require_outputs=False) is None       # failed

        # run 2: same params as run 1, finished ok with a relative sim_dir -> reused
        _log(log, 2, p1)
        reg.start_run(log, 2, p1)
        d2 = _finished_dir(prj, 2)
        reg.finish_run(log, 2, "ok", seconds=10, sim_dir=d2, path2prj=prj)
        df = reg.load_registry(log)
        assert df.loc[2, "sim_dir"] == "outputs/scenario_2_tag", df.loc[2, "sim_dir"]
        assert reg.find_run(log, p1, prj) == 2
        assert reg.resolve_sim_dir(df.loc[2], prj) == d2
        # absolute (legacy) sim_dir that does not exist falls back to the scenario_<n>_* folder
        stale = df.loc[2].copy()
        stale["sim_dir"] = "/no/such/place/scenario_2_old"
        assert reg.resolve_sim_dir(stale, prj) == d2

        # an ok run whose outputs are gone is not reusable
        for f in d2.iterdir():
            f.unlink()
        assert reg.find_run(log, p1, prj) is None
        assert reg.find_run(log, p1, prj, require_outputs=False) == 2

        # partial runs: n_months far below the longest run of the same period are never reused
        p3 = {**BASE, "PMIN": -5.0}
        for idx, months in ((3, 107), (4, 24)):          # 4 stopped early but wrote all outputs
            _log(log, idx, p3)
            reg.start_run(log, idx, p3)
            d = _finished_dir(prj, idx)
            for f in (d / "monthly_stats.csv",):
                f.write_text("month\n" + "\n".join(f"2016-{m % 12 + 1:02d}" for m in range(months)))
            reg.finish_run(log, idx, "ok", seconds=5, sim_dir=d, path2prj=prj)
        df = reg.load_registry(log)
        assert list(df.loc[df.sim_index.isin([3, 4]), "n_months"]) == [107, 24]
        pm = reg.partial_rows(df)
        assert bool(pm[df.index[df.sim_index == 4][0]]) and not bool(pm[df.index[df.sim_index == 3][0]])
        assert reg.find_run(log, p3, prj) == 3            # latest run (4) is partial -> falls back to 3
        assert reg.count_months(prj / "no-such-dir") is None

        # unknown row id -> False, status validation, safe() swallows errors
        assert reg.finish_run(log, 99, "ok") is False
        try:
            reg.finish_run(log, 2, "done")
            raise AssertionError("bad status accepted")
        except ValueError:
            pass
        assert reg.safe(reg.finish_run, tmp / "missing.csv", 0, "ok") is None
    print("test_registry: all assertions passed")


if __name__ == "__main__":
    main()
