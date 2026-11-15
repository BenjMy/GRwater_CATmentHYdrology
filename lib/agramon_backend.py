"""
agramon_backend.py - CATHY build/run abstraction (REFACTORING.md 3.2)
======================================================================
A ``CathyBackend`` knows how to create the pyCATHY project object and how to run
one monthly processor call. Only the current behaviour is implemented
(``scf_lai``: the SCF/psi-min + LAI-vegetation CATHY build, rain and ETp combined
in ATMBC). A ``with_irr`` backend (rain through the irrigation file) is
deliberately NOT implemented: the A/B test of 2026-11-15 (TO-DO.md 5.2,
``lib/ab_withirr``) showed that the irrigation module adds to the same ATMPOT
array, so it reproduces ATMBC=rain-ETp and does not remove the ETa artefact.
"""
from __future__ import annotations

from agramon_log import get_logger

log = get_logger(__name__)

BACKENDS = ("scf_lai",)


class CathyBackend:
    name = "abstract"

    def create_simu(self, sim_dir):
        raise NotImplementedError

    def run_month(self, simu, dtmin, dtmax, deltat):
        raise NotImplementedError


class ScfLaiBackend(CathyBackend):
    """Current behaviour (moved verbatim from initialise_project / run_one_month)."""
    name = "scf_lai"

    def create_simu(self, sim_dir):
        from pyCATHY import cathy_tools
        return cathy_tools.CATHY(dirName=str(sim_dir))

    def run_month(self, simu, dtmin, dtmax, deltat):
        simu.run_processor(
            IPRT1=2,
            DTMIN=dtmin,
            DTMAX=dtmax,
            DELTAT=deltat,
            TRAFLAG=0,
            verbose=False,
        )


def get_backend(name: str = "scf_lai") -> CathyBackend:
    if name == "scf_lai":
        return ScfLaiBackend()
    if name == "with_irr":
        raise NotImplementedError(
            "with_irr backend not implemented: the A/B test showed it does not remove the "
            "ETa artefact (TO-DO.md 5.2). Use --eta-net-neg-mode gross instead.")
    raise ValueError(f"unknown CATHY backend {name!r}; choose from {BACKENDS}")
