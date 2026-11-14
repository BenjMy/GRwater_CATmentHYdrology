"""
agramon_log.py — logging for the Agramon pipeline (REFACTORING.md 2.10)
======================================================================
Console output is unchanged by default: records go to stdout with a bare
``%(message)s`` format, so what used to be ``print`` looks identical. What
is new: levels (warnings/errors are real ``WARNING``/``ERROR`` records),
``--log-level`` and an optional ``--log-file`` mirror.

Usage in a module::

    from agramon_log import get_logger
    log = get_logger(__name__)
    log.info("  something happened")

``get_logger`` installs the default stdout handler on first use, so modules
imported from other scripts still print. ``setup_logging`` (called once by
the CLI) changes the level and adds the file handler; it is idempotent.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

ROOT_NAME = "agramon"
_FORMAT = "%(message)s"
_FILE_FORMAT = "%(asctime)s %(levelname)-7s %(name)s: %(message)s"


class _StdoutHandler(logging.StreamHandler):
    """StreamHandler that always writes to the *current* ``sys.stdout``."""

    def __init__(self):
        super().__init__(sys.stdout)

    @property
    def stream(self):
        return sys.stdout

    @stream.setter
    def stream(self, _value):
        pass


def _root() -> logging.Logger:
    root = logging.getLogger(ROOT_NAME)
    if not any(isinstance(h, _StdoutHandler) for h in root.handlers):
        h = _StdoutHandler()
        h.setFormatter(logging.Formatter(_FORMAT))
        root.addHandler(h)
        root.setLevel(logging.INFO)
        root.propagate = False
    return root


def get_logger(name: str = "") -> logging.Logger:
    """Child logger of ``agramon`` (``__main__`` maps to ``agramon.main``)."""
    _root()
    name = "main" if name in ("", "__main__") else name
    return logging.getLogger(f"{ROOT_NAME}.{name}")


def setup_logging(level: str = "INFO", logfile: Optional[str] = None) -> logging.Logger:
    """Set the level (name or number) and optionally mirror records to ``logfile``."""
    root = _root()
    lv = level if isinstance(level, int) else logging.getLevelName(str(level).upper())
    if not isinstance(lv, int):
        raise ValueError(f"unknown log level {level!r}")
    root.setLevel(lv)
    if logfile:
        path = str(Path(logfile))
        if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(Path(path).resolve())
                   for h in root.handlers):
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(path, encoding="utf-8")
            fh.setFormatter(logging.Formatter(_FILE_FORMAT))
            root.addHandler(fh)
    return root
