"""Configuration loading and path resolution.

Every rule is read from config/*.yaml; nothing is hard-coded in scripts.
"""
from __future__ import annotations

import hashlib
import json
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


def load_yaml(name: str) -> dict[str, Any]:
    with open(CONFIG_DIR / name) as fh:
        return yaml.safe_load(fh)


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
    cfg = Config(study=load_yaml("study.yaml"), sources=load_yaml("sources.yaml"))
    unknown = set(cfg.views) - set(cfg.sources["tcga"]["molecular"])
    if unknown:
        raise ValueError(f"views without a raw source: {sorted(unknown)}")
    return cfg
