"""Step 04: run one modelling job (nested CV + tuning + refit + SHAP) for (task, view, cohort).

Usage: python scripts/04_run_job.py <task> <view> <cohort>
       e.g. tumor_status se BRCA | cancer_type gex PANCAN
Threads: set OMP_NUM_THREADS in the environment BEFORE starting Python (the driver does this).
Checkpointed per outer fold; re-running resumes. Outputs: <data_dir>/runs/<task>/<view>/<cohort>/
(models, SHAP arrays, per-sample predictions; git-ignored) and results/<task>/<view>/<cohort>/
(metrics, fold table, SHAP summaries; committed).
"""
from __future__ import annotations

import os
import sys

from _common import log

import torch

from ctsm.config import TASKS, load_config
from ctsm.pipeline import run_job


def main(task: str, view: str, cohort: str) -> None:
    cfg = load_config()
    if task not in TASKS or view not in cfg.views:
        raise SystemExit(f"unknown task/view: {task} {view}")
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
    s = run_job(cfg, task, view, cohort, log=log)
    m = s["pooled_oof"]
    key = "auroc" if task == "tumor_status" else "accuracy"
    log(f"DONE {task}/{view}/{cohort}: pooled {key} {m[key]:.4f}, balanced accuracy {m['balanced_accuracy']:.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
