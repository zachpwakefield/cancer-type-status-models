"""Step 02: raw TCGA tables -> compact sample x feature matrices (D-05..D-07).

Usage: python scripts/02_extract_matrices.py [view ...]       (default: every view)
Outputs: <data_dir>/matrices/<view>/ (see ctsm/matrix_io.py)
Summary (committed): results/data/extraction_<view>.csv
"""
from __future__ import annotations

import math
import sys

import numpy as np
import pandas as pd
from _common import log

from ctsm.annotation import gene_annotation
from ctsm.config import load_config
from ctsm.extract import extract_as_view, extract_gex
from ctsm.provenance import run_record
from ctsm.splits import min_inner_train


def cohort_rows(cfg, samples: pd.DataFrame) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    """Row positions of every task cohort in `samples`, and its smallest inner-train size."""
    sp = cfg.study["splits"]
    pos = pd.Series(np.arange(len(samples)), index=samples["file_id"])
    rows, n_ref = {}, {}
    for p in sorted((cfg.data_dir / "splits").glob("*.parquet")):
        cohort = p.stem.split("__", 1)[1]
        s = pd.read_parquet(p)
        rows[cohort] = np.sort(pos.loc[s["file_id"]].to_numpy())
        n_ref[cohort] = min_inner_train(s, sp["outer_folds"], sp["inner_folds"])
    return rows, n_ref


def main(views: list[str]) -> None:
    cfg = load_config()
    pp = cfg.study["preprocess"]
    samples = pd.read_parquet(cfg.data_dir / "cohort" / "samples.parquet")
    rows, n_ref = cohort_rows(cfg, samples)
    ann = gene_annotation(cfg)
    log(f"{len(samples):,} samples, {len(rows)} task cohorts, gene annotation {len(ann):,} genes")
    for view in views:
        # D-07: support threshold = the fold-level rule applied to the SMALLEST training set
        frac = pp["gex_min_expressed_frac"] if view == "gex" else pp["min_observed_frac"]
        k_min = {c: math.ceil(frac * n_ref[c]) for c in rows}
        log(f"=== {view} === k_min {k_min}")
        fn = extract_gex if view == "gex" else extract_as_view
        written, stats = fn(cfg, view, samples, rows, k_min, ann, log=log) if view != "gex" else \
            fn(cfg, samples, rows, k_min, ann, log=log)
        pd.DataFrame([stats]).to_csv(cfg.results_dir / "data" / f"extraction_{view}.csv", index=False)
        run_record(cfg, f"02_extract_{view}", [cfg.raw_matrix(view)], written, extra=stats)


if __name__ == "__main__":
    cfg0 = load_config()
    req = sys.argv[1:] or cfg0.views
    bad = set(req) - set(cfg0.views)
    if bad:
        raise SystemExit(f"unknown views {sorted(bad)}")
    main(req)
