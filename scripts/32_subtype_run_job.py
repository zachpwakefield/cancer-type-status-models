"""Step 32: one subtype modelling job (nested CV + tuning + refit + SHAP) for (view, cohort).

Usage: python scripts/32_subtype_run_job.py <view> <cohort>       e.g. se BRCA | gex LGG_IDH
Threads: set OMP_NUM_THREADS before Python starts (step 33 does). Checkpointed per outer fold;
re-running resumes. Outputs: <data_dir>/subtype/runs/subtype/<view>/<cohort>/ (models, SHAP
arrays, per-patient predictions; git-ignored) and results/subtype/<view>/<cohort>/ (metrics.json,
folds.csv, per_class_metrics.csv, confusion_matrix.csv, shap_top_*.csv; committed).
"""
from __future__ import annotations

import os
import sys

from _common import log

import torch

from ctsm.subtype import cohort_specs, run_job, subtype_config


def main(view: str, cohort: str) -> None:
    scfg = subtype_config()
    if view not in scfg.views or cohort not in cohort_specs(scfg.spec):
        raise SystemExit(f"unknown view/cohort: {view} {cohort}")
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
    s = run_job(scfg, view, cohort, log=log)
    m = s["pooled_oof"]
    log(f"DONE subtype/{view}/{cohort}: balanced accuracy {m['balanced_accuracy']:.4f}, "
        f"macro F1 {m['macro_f1']:.4f}, macro AUROC {m['macro_auroc_ovr']:.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
