"""Shared helpers for pipeline scripts: import path, thread limits, logging."""
from __future__ import annotations

import datetime as _dt
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

# one BLAS thread per process unless the caller says otherwise: results must not depend on
# thread scheduling, and the model runner parallelises over processes instead
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

# torch must load its OpenMP runtime before numpy loads another one: in this environment
# (torch 2.2.2 + conda-forge numpy on macOS) the reverse order segfaults in BatchNorm
import torch  # noqa: E402,F401


def log(msg: str) -> None:
    print(f"{_dt.datetime.now().strftime('%H:%M:%S')} {msg}", flush=True)
