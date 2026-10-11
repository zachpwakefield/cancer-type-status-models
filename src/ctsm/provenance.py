"""Provenance: file hashing and run records.

Every pipeline step writes a run record (config fingerprint, git commit, package versions,
inputs and outputs) so any output can be traced to the code, configuration and raw inputs.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

from .config import REPO_ROOT, Config

_CHUNK = 16 * 1024 * 1024


def sha256_file(path: Path | str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(_CHUNK):
            h.update(block)
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True,
                             text=True).stdout.strip() or None
        # generated results are outputs, so only code/config/docs changes count as dirty
        dirty = subprocess.run(["git", "status", "--porcelain", "--", ".", ":(exclude)results"],
                               cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()
        return f"{sha}+dirty" if sha and dirty else sha
    except (OSError, subprocess.SubprocessError):
        return None


def package_versions(names: Iterable[str] = ("numpy", "pandas", "pyarrow", "scipy", "sklearn", "torch",
                                             "optuna", "shap")) -> dict[str, str]:
    out = {"python": platform.python_version()}
    for name in names:
        try:
            out[name] = getattr(__import__(name), "__version__", "unknown")
        except ImportError:
            out[name] = "not installed"
    return out


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    tmp.replace(path)


def run_record(cfg: Config, step: str, inputs: list[Path], outputs: list[Path], extra: dict | None = None) -> dict:
    """Write <data_dir>/provenance/<step>.json and return it."""
    rec = {
        "step": step,
        "study_id": cfg.study_id,
        "config_fingerprint": cfg.fingerprint,
        "git_commit": git_commit(),
        "argv": sys.argv,
        "finished_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "platform": platform.platform(),
        "packages": package_versions(),
        "inputs": [str(p) for p in inputs],
        "outputs": [{"path": str(Path(p).relative_to(cfg.data_dir)) if str(p).startswith(str(cfg.data_dir)) else str(p),
                     "bytes": Path(p).stat().st_size} for p in outputs],
    }
    if extra:
        rec["extra"] = extra
    write_json(cfg.data_dir / "provenance" / f"{step}.json", rec)
    return rec
