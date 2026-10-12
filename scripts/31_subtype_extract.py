"""Step 31: raw TCGA tables -> matrices for the subtype cohorts (D-37).

Usage: python scripts/31_subtype_extract.py [view ...]        (default: every view)

The same extraction as step 02 (ctsm/extract.py: genome scope, host-gene biotype, AFE/ALE
structural zeros, lossless encoding), over the union of the subtype cohorts' samples, with the D-07
pre-filter recomputed for EACH subtype cohort from its own splits: support >= ceil(fraction x the
smallest inner-train set of that cohort) and non-constant over the cohort. The proof in
docs/LEAKAGE_POLICY.md then applies unchanged: no feature a fold-level filter could keep on any
training set of the cohort is dropped. The pre-filter reads values and the cohort's membership,
never the labels. Per-sample operations make every stored value identical to step 02's value for
the same sample and feature.

Outputs: <data_dir>/subtype/matrices/<view>/ (git-ignored); results/subtype/data/extraction_<view>.csv
"""
from __future__ import annotations

import math
import sys

import numpy as np
import pandas as pd
from _common import log

from ctsm.annotation import gene_annotation
from ctsm.extract import extract_as_view, extract_gex
from ctsm.provenance import run_record
from ctsm.splits import min_inner_train
from ctsm.subtype import TASK, load_tables, subtype_config


def main(views: list[str]) -> None:
    scfg = subtype_config()
    pp, sp = scfg.study["preprocess"], scfg.study["splits"]
    tables = load_tables(scfg)
    order = pd.concat(tables.values())[["file_id", "project", "case_id"]].drop_duplicates("file_id")
    samples = order.sort_values(["project", "case_id"]).reset_index(drop=True)[["file_id"]]
    pos = pd.Series(np.arange(len(samples)), index=samples["file_id"])
    rows, n_ref = {}, {}
    for cohort, tab in tables.items():
        s = pd.read_parquet(scfg.data_dir / "splits" / f"{TASK}__{cohort}.parquet")
        rows[cohort] = np.sort(pos.loc[s["file_id"]].to_numpy())
        n_ref[cohort] = min_inner_train(s, sp["outer_folds"], sp["inner_folds"])
    ann = gene_annotation(scfg)
    log(f"{len(samples):,} samples in {len(rows)} subtype cohorts")
    for view in views:
        frac = pp["gex_min_expressed_frac"] if view == "gex" else pp["min_observed_frac"]
        k_min = {c: math.ceil(frac * n_ref[c]) for c in rows}
        log(f"=== {view} === k_min {k_min}")
        if view == "gex":
            written, stats = extract_gex(scfg, samples, rows, k_min, ann, log=log)
        else:
            written, stats = extract_as_view(scfg, view, samples, rows, k_min, ann, log=log)
        out = scfg.results_dir / TASK / "data"
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([stats]).to_csv(out / f"extraction_{view}.csv", index=False)
        run_record(scfg, f"31_extract_{view}", [scfg.raw_matrix(view)], written, extra=stats)


if __name__ == "__main__":
    cfg0 = subtype_config()
    req = sys.argv[1:] or cfg0.views
    bad = set(req) - set(cfg0.views)
    if bad:
        raise SystemExit(f"unknown views {sorted(bad)}")
    main(req)
