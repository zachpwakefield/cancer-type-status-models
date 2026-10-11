"""Step 14 (B-07): is the splicing advantage for neighbouring cancer types a site effect?

In TCGA every tissue-source site (TSS, barcode field 2) contributes to exactly one project, so a
classifier of related projects (COAD vs READ, GBM vs LGG, ...) can in principle learn site
handling (ischaemia time, storage, extraction) rather than biology. RNA processing readouts
such as intron retention are plausibly more sensitive to handling than gene expression.

Probe (diagnostic only; it never touches the phase-1 models or their estimates): for each pair
of related classes and each view, the same simple, untuned classifier is cross-validated two ways
  case-grouped  StratifiedGroupKFold by case (as in phase 1)
  site-held-out StratifiedGroupKFold by TSS (every test site unseen in training)
Pipeline, fitted inside each training fold: top-2,000 variance (NaN-aware) -> median impute ->
standardise -> L2 logistic regression (C = 0.05, balanced classes). A view whose accuracy collapses when sites
are held out was relying on site-specific signal.

Outputs: results/biology/07_site_probe/ (tables; figures: scripts/15_bio_figures.py)
"""
from __future__ import annotations

import argparse
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from _common import log
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ctsm.config import PAN, load_config
from ctsm.matrix_io import load_view

warnings.filterwarnings("ignore")
PAIRS = [("COAD", "READ"), ("GBM", "LGG"), ("ESCA", "STAD"), ("LUAD", "LUSC"), ("KIRC", "KIRP"), ("CESC", "HNSC")]
SEEDS = [0, 1, 2]


class TopVariance(BaseEstimator, TransformerMixin):
    def __init__(self, k: int = 2000):
        self.k = k

    def fit(self, X, y=None):
        v = np.nanvar(X, axis=0)
        self.keep_ = np.argsort(-np.nan_to_num(v, nan=-1))[: self.k]
        return self

    def transform(self, X):
        return X[:, self.keep_]


def probe(args) -> list[dict]:
    a, b, view = args
    cfg = load_config()
    tab = pd.read_parquet(cfg.data_dir / "tasks" / f"cancer_type__{PAN}.parquet")
    sub = tab[tab["label"].isin([a, b])].reset_index(drop=True)
    M, _ = load_view(cfg.matrix_dir, view, sub["file_id"], PAN)
    X = M.get()
    keep = np.isfinite(X).mean(0) >= 0.5                      # label-free, same for both schemes
    X = X[:, keep]
    y = (sub["label"] == b).to_numpy().astype(int)
    rows = []
    for scheme, groups in [("case-grouped", sub["case_id"]), ("site-held-out", sub["tss"])]:
        for seed in SEEDS:
            p = np.full(len(y), np.nan)
            for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(X, y, groups):
                m = make_pipeline(TopVariance(2000), SimpleImputer(strategy="median"), StandardScaler(),
                                  LogisticRegression(C=0.05, class_weight="balanced", max_iter=3000))
                m.fit(X[tr], y[tr])
                p[te] = m.predict_proba(X[te])[:, 1]
            rows.append({"pair": f"{a} vs {b}", "view": view, "scheme": scheme, "seed": seed, "n": len(y),
                         "n_sites": sub["tss"].nunique(), "balanced_accuracy": balanced_accuracy_score(y, p >= 0.5),
                         "auroc": roc_auc_score(y, p)})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    cfg = load_config()
    out = cfg.results_dir / "biology" / "07_site_probe"
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(a, b, v) for a, b in PAIRS for v in cfg.views]
    with ProcessPoolExecutor(args.workers) as ex:
        res = pd.DataFrame([r for rr in ex.map(probe, jobs) for r in rr])
    res.to_csv(out / "site_probe_by_seed.csv", index=False)
    agg = res.groupby(["pair", "view", "scheme"])[["balanced_accuracy", "auroc"]].mean().reset_index()
    wide = agg.pivot_table(index=["pair", "view"], columns="scheme", values="balanced_accuracy").reset_index()
    wide["drop"] = wide["case-grouped"] - wide["site-held-out"]
    wide.to_csv(out / "site_probe_summary.csv", index=False)
    log("site probe:\n" + wide.to_string(index=False))

    log("site probe done")


if __name__ == "__main__":
    main()
