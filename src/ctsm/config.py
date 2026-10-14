"""Configuration loading and path resolution.

Every rule is read from config/*.yaml; nothing is hard-coded in scripts.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "config"

TASKS = ("tumor_status", "cancer_type")
PAN = "PANCAN"                      # cohort name of the pan-cancer (cancer_type) task
RMATS_VIEWS = ("se", "mxe", "a3ss", "a5ss", "ri")
TERMINAL_VIEWS = ("afe", "ale")
INT16_VIEWS = RMATS_VIEWS + ("hit",)
LOCAL_ENV = CONFIG_DIR / "local.env"   # machine-specific paths, git-ignored; see config/local.env.example
_VAR = re.compile(r"\$\{([A-Z][A-Z0-9_]*)\}")


def load_yaml(name: str) -> dict[str, Any]:
    with open(CONFIG_DIR / name) as fh:
        return yaml.safe_load(fh)


def settings() -> dict[str, str]:
    """Machine-specific settings: KEY=VALUE lines of config/local.env, overridden by environment
    variables of the same name (CTSM_* and PYTHON)."""
    out: dict[str, str] = {}
    if LOCAL_ENV.exists():
        for line in LOCAL_ENV.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = os.path.expanduser(v.strip().strip("\"'"))
    out.update({k: v for k, v in os.environ.items() if k.startswith("CTSM_") or k == "PYTHON"})
    return out


def expand(value: Any, env: dict[str, str] | None = None) -> Any:
    """Fill ${VAR} in strings, recursively through dicts and lists. An unset variable stays as
    written, so a path that needs it fails with the variable's name in the error."""
    env = settings() if env is None else env
    if isinstance(value, str):
        return _VAR.sub(lambda m: env.get(m.group(1), m.group(0)), value)
    if isinstance(value, dict):
        return {k: expand(v, env) for k, v in value.items()}
    if isinstance(value, list):
        return [expand(v, env) for v in value]
    return value


def portable(path: Any, env: dict[str, str] | None = None) -> str:
    """A resolved path with the longest matching CTSM_* setting written back as ${VAR}, for
    records that are committed (no machine paths in the repository)."""
    env = settings() if env is None else env
    s = str(path)
    for k, v in sorted(env.items(), key=lambda kv: -len(kv[1])):
        if k.startswith("CTSM_") and v and (s == v or s.startswith(v.rstrip("/") + "/")):
            return "${" + k + "}" + s[len(v):]
    return s


def stable_hash(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()


@dataclass(frozen=True)
class Config:
    study: dict[str, Any]
    sources: dict[str, Any]

    @property
    def study_id(self) -> str:
        return self.study["study_id"]

    @property
    def views(self) -> list[str]:
        return list(self.study["views"])

    @property
    def fingerprint(self) -> str:
        return stable_hash({"study": self.study, "sources": self.sources})

    # ---- locations -------------------------------------------------------------------
    @property
    def data_dir(self) -> Path:
        """Processed data and run artefacts (git-ignored)."""
        return REPO_ROOT / self.study["data_root"] / self.study_id

    @property
    def matrix_dir(self) -> Path:
        return self.data_dir / "matrices"

    @property
    def runs_dir(self) -> Path:
        return self.data_dir / "runs"

    @property
    def results_dir(self) -> Path:
        """Small committed outputs."""
        return REPO_ROOT / "results"

    def run_dir(self, task: str, view: str, cohort: str) -> Path:
        return self.runs_dir / task / view / cohort

    def result_dir(self, task: str, view: str, cohort: str) -> Path:
        return self.results_dir / task / view / cohort

    # ---- raw sources -----------------------------------------------------------------
    def tcga_path(self, *keys: str) -> Path:
        node: Any = self.sources["tcga"]
        for k in keys:
            node = node[k]
        return Path(self.sources["tcga_root"]) / node

    @property
    def gencode_gtf(self) -> Path:
        return Path(self.sources["gencode_v45_gtf"])

    def raw_matrix(self, view: str) -> Path:
        return self.tcga_path("molecular", view)


def load_config() -> Config:
    cfg = Config(study=load_yaml("study.yaml"), sources=expand(load_yaml("sources.yaml")))
    unknown = set(cfg.views) - set(cfg.sources["tcga"]["molecular"])
    if unknown:
        raise ValueError(f"views without a raw source: {sorted(unknown)}")
    return cfg
