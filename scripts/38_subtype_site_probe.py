"""Step 38 (D-47): is the subtype signal a tissue-source-site signal?

The site probe of step 14 (D-20), applied to the subtype cohorts. For each cohort and view the
same untuned classifier is cross-validated two ways
  case-grouped   StratifiedGroupKFold(5) by case (one primary tumour per patient)
  site-held-out  StratifiedGroupKFold(5) by tissue-source site (every test site unseen in training)
Pipeline, fitted inside each training fold: top-2,000 variance (NaN-aware) -> median impute ->
standardise -> L2 logistic regression (C = 0.05, balanced classes; multinomial for > 2 classes).
A cohort is probed only if every class comes from >= min_sites_per_class sites. The site-class
confounding of each cohort is Cramer's V between TSS and class, from the labels alone.

Predictions (config/followup.yaml): P1, the median drop of the 8 AS views exceeds the GEX drop in
most cohorts (sign test); P2, the drop rises with Cramer's V (Spearman over cohorts).
Diagnostic only: it never touches a subtype model or its estimate.

Outputs: results/subtype/site_probe/ (tables; figure U9: scripts/36_subtype_figures.py)
Usage: python scripts/38_subtype_site_probe.py [--workers 6]
"""
from __future__ import annotations

import argparse
import importlib.util
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from _common import REPO, log
from scipy.stats import binomtest, chi2_contingency, spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ctsm.matrix_io import load_view
from ctsm.subtype import TASK, load_tables, subtype_config

warnings.filterwarnings("ignore")
_spec = importlib.util.spec_from_file_location("site14", Path(__file__).with_name("14_bio_site_probe.py"))
_site14 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_site14)
TopVariance, SEEDS = _site14.TopVariance, _site14.SEEDS
CFG = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["subtype_site_probe"]


def cramers_v(a: pd.Series, b: pd.Series) -> float:
    t = pd.crosstab(a, b).to_numpy()
    chi2 = chi2_contingency(t, correction=False)[0]
    return float(np.sqrt(chi2 / (t.sum() * (min(t.shape) - 1))))


def oof_proba(X: np.ndarray, y: np.ndarray, groups, seed: int, k: int) -> np.ndarray:
    """Pooled out-of-fold class probabilities of the probe pipeline, StratifiedGroupKFold(5) by
    `groups`. A class absent from a training fold gets probability 0 in that fold's test rows."""
    p = np.zeros((len(y), k))
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(X, y, groups):
        m = make_pipeline(TopVariance(2000), SimpleImputer(strategy="median"), StandardScaler(),
                          LogisticRegression(C=0.05, class_weight="balanced", max_iter=3000))
        m.fit(X[tr], y[tr])
        p[np.ix_(te, m.classes_)] = m.predict_proba(X[te])
    return p


def probe(args) -> list[dict]:
    cohort, view = args
    scfg = subtype_config()
    tab = load_tables(scfg)[cohort]
    M, _ = load_view(scfg.matrix_dir, view, tab["file_id"], cohort)
    X = M.get()
    X = X[:, np.isfinite(X).mean(0) >= 0.5]                   # label-free, same for both schemes
    y = tab["y"].to_numpy()
    k = int(y.max()) + 1
    rows = []
    for scheme, groups in [("case-grouped", tab["case_id"]), ("site-held-out", tab["tss"])]:
        for seed in SEEDS:
            p = oof_proba(X, y, groups, seed, k)
            rows.append({"cohort": cohort, "view": view, "scheme": scheme, "seed": seed, "n": len(y),
                         "n_classes": k, "n_sites": tab["tss"].nunique(),
                         "balanced_accuracy": balanced_accuracy_score(y, p.argmax(1))})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    scfg = subtype_config()
    out = scfg.results_dir / TASK / "site_probe"
    out.mkdir(parents=True, exist_ok=True)
    tables = load_tables(scfg)
    conf = []
    for cohort, t in tables.items():
        sites = t.groupby("label")["tss"].nunique()
        conf.append({"cohort": cohort, "n": len(t), "n_sites": t["tss"].nunique(),
                     "min_sites_per_class": int(sites.min()), "class_fewest_sites": sites.idxmin(),
                     "cramers_v_site_class": cramers_v(t["tss"], t["label"]),
                     "probed": bool(sites.min() >= CFG["min_sites_per_class"])})
    conf = pd.DataFrame(conf)
    conf.to_csv(out / "site_class_confounding.csv", index=False)
    jobs = [(c, v) for c in conf.loc[conf.probed, "cohort"] for v in scfg.views]
    with ProcessPoolExecutor(a.workers) as ex:
        res = pd.DataFrame([r for rr in ex.map(probe, jobs) for r in rr])
    res.to_csv(out / "site_probe_by_seed.csv", index=False)
    wide = (res.groupby(["cohort", "view", "scheme"])["balanced_accuracy"].mean().unstack("scheme").reset_index())
    wide["drop"] = wide["case-grouped"] - wide["site-held-out"]
    wide.to_csv(out / "site_probe_summary.csv", index=False)

    # P1 and P2
    as_med = wide[wide.view != "gex"].groupby("cohort")["drop"].median().rename("drop_as_median")
    per = (conf.set_index("cohort").join(wide[wide.view == "gex"].set_index("cohort")["drop"].rename("drop_gex"))
           .join(as_med).reset_index())
    per = per[per.probed]
    per["as_drop_exceeds_gex"] = per["drop_as_median"] > per["drop_gex"]
    per.to_csv(out / "site_probe_by_cohort.csv", index=False)
    d = per["drop_as_median"] - per["drop_gex"]
    tests = [{"id": "P1", "claim": "AS-median drop > GEX drop", "n": int((d != 0).sum()),
              "k": int((d > 0).sum()), "estimate": float(np.median(d)),
              "p": binomtest(int((d > 0).sum()), int((d != 0).sum())).pvalue}]
    for col, name in [("drop_gex", "GEX"), ("drop_as_median", "AS-view median")]:
        r = spearmanr(per["cramers_v_site_class"], per[col])
        tests.append({"id": "P2", "claim": f"drop rises with Cramer's V ({name})", "n": len(per),
                      "estimate": float(r.statistic), "p": float(r.pvalue)})
    pd.DataFrame(tests).to_csv(out / "site_probe_tests.csv", index=False)
    log("site probe:\n" + per.round(3).to_string(index=False) + "\n" + pd.DataFrame(tests).to_string(index=False))


if __name__ == "__main__":
    main()
