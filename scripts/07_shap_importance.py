"""Step 07: full SHAP importance tables for every job (all selected features, not only the top 200),
and per-class tables for the cancer-type models. Cached under <data_dir>/bio_cache/importance/
for the interpretation steps 08-13.

Usage: python scripts/07_shap_importance.py [--workers 6]
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from _common import log

from ctsm.biology import full_importance
from ctsm.config import PAN, load_config


def one(args) -> str:
    task, view, cohort, k = args
    cfg = load_config()
    full_importance(cfg, task, view, cohort, class_index=k)
    return f"{task}/{view}/{cohort}/{'all' if k is None else k}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    cfg = load_config()
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = list(scope.loc[scope.tumor_status_in_scope, "project"])
    n_cls = pd.read_parquet(cfg.data_dir / "tasks" / f"cancer_type__{PAN}.parquet")["y"].max() + 1
    jobs = [("tumor_status", v, c, None) for v in cfg.views for c in cancers]
    jobs += [("cancer_type", v, PAN, None) for v in cfg.views]
    jobs += [("cancer_type", v, PAN, k) for v in cfg.views for k in range(n_cls)]
    with ProcessPoolExecutor(a.workers) as ex:
        for i, r in enumerate(ex.map(one, jobs), 1):
            if i % 25 == 0 or i == len(jobs):
                log(f"{i}/{len(jobs)} {r}")


if __name__ == "__main__":
    main()
