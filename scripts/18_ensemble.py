"""Step 18 (B-09): how much does combining the modalities buy, in each task separately?

Every model here is built from OUT-OF-FOLD probabilities only: each sample was predicted by a
model that never saw it. Three ways of combining them, in increasing order of what has to be
estimated:
    GEX + AS (1:1)   0.5 x GEX + 0.5 x the mean of the 8 AS views. Fixed in advance, nothing
                     fitted: the honest headline.
    all equal        the mean of the 9 views. Also fixed.
    stacked          a cross-fitted logistic regression on the 9 views' probabilities: for each
                     outer fold, fitted on the other folds' out-of-fold probabilities and applied
                     to this one. Nothing fitted on the samples it scores, but the base models
                     behind the stacker's training probabilities did see them, so a small
                     optimism remains (stated in BIOLOGY_METHODS, D-23). Reported as a ceiling,
                     never as the headline.

Per task the question differs, so the tables differ:
    tumour status    per cancer: AUROC, balanced accuracy, errors, bootstrap CI, McNemar vs GEX,
                     and which errors AS rescues.
    cancer type      overall accuracy and macro recall, per-class recall, the same tests.

Outputs: results/biology/09_ensemble/ (tables; figures: scripts/20_figures_status.py and scripts/21_figures_type.py)
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import binomtest, rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score, recall_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from ctsm.biology import bh, bio_config, full_importance
from ctsm.config import PAN, load_config

N_BOOT = 1000


def oof(cfg, task, view, cohort) -> pd.DataFrame:
    return pd.read_parquet(cfg.run_dir(task, view, cohort) / "oof_predictions.parquet")


def stack(P: np.ndarray, y: np.ndarray, folds: np.ndarray, n_classes: int, seed: int) -> np.ndarray:
    """Cross-fitted stacking: fold k is scored by a model fitted on the other folds only."""
    out = np.zeros((len(y), n_classes))
    for k in np.unique(folds):
        tr, te = folds != k, folds == k
        sc = StandardScaler().fit(P[tr])
        m = LogisticRegression(C=1.0, max_iter=1000, random_state=seed)
        m.fit(sc.transform(P[tr]), y[tr])
        pr = m.predict_proba(sc.transform(P[te]))
        out[np.ix_(te, m.classes_)] = pr
    return out


def auc(y: np.ndarray, score: np.ndarray) -> float:
    """AUROC from the rank sum of the positives (identical to roc_auc_score, ties averaged, ~20x faster:
    the bootstrap calls it a thousand times per model)."""
    npos = int(y.sum())
    nneg = len(y) - npos
    if npos == 0 or nneg == 0:
        raise ValueError("one class")
    return float((rankdata(score)[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def boot_ci(stat, *arrays, n_boot: int = N_BOOT, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(arrays[0])
    vals = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        try:
            vals.append(stat(*[a[i] for a in arrays]))
        except ValueError:                      # a resample can miss a class
            continue
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if vals else (np.nan, np.nan)


def mcnemar(a_correct: np.ndarray, b_correct: np.ndarray) -> dict:
    b01 = int((a_correct & ~b_correct).sum())
    b10 = int((~a_correct & b_correct).sum())
    p = binomtest(b01, b01 + b10, 0.5).pvalue if b01 + b10 else 1.0
    return {"only_a_correct": b01, "only_b_correct": b10, "p_mcnemar": float(p)}


def combine(views: list[str], probs: dict[str, np.ndarray], seed, y, folds, n_classes) -> dict[str, np.ndarray]:
    """The four combinations, plus every single view, as class-probability matrices."""
    as_views = [v for v in views if v != "gex"]
    out = {v.upper(): probs[v] for v in views}
    out["AS ensemble"] = np.mean([probs[v] for v in as_views], axis=0)
    out["all equal"] = np.mean([probs[v] for v in views], axis=0)
    out["GEX + AS"] = 0.5 * probs["gex"] + 0.5 * out["AS ensemble"]
    out["stacked"] = stack(np.hstack([probs[v] for v in views]), y, folds, n_classes, seed)
    return out


def tumour_status(cfg, cancers, seed) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metrics, tests, rescue = [], [], []
    for cancer in cancers:
        ref = oof(cfg, "tumor_status", "gex", cancer).sort_values("file_id").reset_index(drop=True)
        y, folds = ref["y"].to_numpy(), ref["outer_fold"].to_numpy()
        probs = {}
        for view in cfg.views:
            d = oof(cfg, "tumor_status", view, cancer).sort_values("file_id").reset_index(drop=True)
            assert (d["file_id"].to_numpy() == ref["file_id"].to_numpy()).all()
            probs[view] = d[["p_Normal", "p_Tumour"]].to_numpy()
        models = combine(cfg.views, probs, seed, y, folds, 2)
        correct = {}
        for name, P in models.items():
            p1 = P[:, 1]
            yhat = (p1 >= 0.5).astype(int)
            correct[name] = yhat == y
            lo, hi = boot_ci(auc, y, p1, seed=seed)
            metrics.append({"cancer": cancer, "model": name, "n": len(y), "n_normal": int((y == 0).sum()),
                            "auroc": roc_auc_score(y, p1), "auroc_lo": lo, "auroc_hi": hi,
                            "balanced_accuracy": balanced_accuracy_score(y, yhat),
                            "n_errors": int((yhat != y).sum())})
        for name in ("AS ensemble", "GEX + AS", "stacked"):
            tests.append({"cancer": cancer, "model": name, "reference": "GEX",
                          **mcnemar(correct["GEX"], correct[name])})
        wrong_gex = ~correct["GEX"]
        rescue.append({"cancer": cancer, "gex_errors": int(wrong_gex.sum()),
                       "rescued_by_as": int((wrong_gex & correct["AS ensemble"]).sum()),
                       "rescued_by_fused": int((wrong_gex & correct["GEX + AS"]).sum()),
                       "as_only_errors": int((correct["GEX"] & ~correct["AS ensemble"]).sum()),
                       "both_wrong": int((wrong_gex & ~correct["AS ensemble"]).sum())})
        log(f"  {cancer}: GEX {correct['GEX'].sum()}/{len(y)}, fused {correct['GEX + AS'].sum()}, "
            f"stacked {correct['stacked'].sum()}")
    t = pd.DataFrame(tests)
    t["q_mcnemar"] = bh(t["p_mcnemar"].to_numpy())
    return pd.DataFrame(metrics), t, pd.DataFrame(rescue)


def cancer_type(cfg, seed) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ref = oof(cfg, "cancer_type", "gex", PAN).sort_values("file_id").reset_index(drop=True)
    classes = [c[2:] for c in ref.columns if c.startswith("p_")]
    y, folds = ref["y"].to_numpy(), ref["outer_fold"].to_numpy()
    probs = {}
    for view in cfg.views:
        d = oof(cfg, "cancer_type", view, PAN).sort_values("file_id").reset_index(drop=True)
        assert (d["file_id"].to_numpy() == ref["file_id"].to_numpy()).all()
        probs[view] = d[[f"p_{c}" for c in classes]].to_numpy()
    models = combine(cfg.views, probs, seed, y, folds, len(classes))
    metrics, per_class, correct = [], [], {}
    for name, P in models.items():
        yhat = P.argmax(1)
        correct[name] = yhat == y
        lo, hi = boot_ci(lambda a, b: (a == b).mean(), y, yhat, seed=seed)
        metrics.append({"model": name, "n": len(y), "accuracy": float((yhat == y).mean()),
                        "accuracy_lo": lo, "accuracy_hi": hi,
                        "macro_recall": recall_score(y, yhat, average="macro"),
                        "macro_f1": f1_score(y, yhat, average="macro"), "n_errors": int((yhat != y).sum())})
        rec = recall_score(y, yhat, average=None, labels=np.arange(len(classes)))
        per_class.append(pd.DataFrame({"model": name, "label": classes, "recall": rec,
                                       "n": np.bincount(y, minlength=len(classes))}))
    tests = [{"model": name, "reference": "GEX", **mcnemar(correct["GEX"], correct[name])}
             for name in ("AS ensemble", "GEX + AS", "stacked")]
    t = pd.DataFrame(tests)
    t["q_mcnemar"] = bh(t["p_mcnemar"].to_numpy())
    return pd.DataFrame(metrics), t, pd.concat(per_class, ignore_index=True)


def task_overlap(cfg, cancers, k: int = 100) -> pd.DataFrame:
    """Do the two tasks rely on the same events? Top-k features of every tumour-status model
    against those of every cancer-type class, within the features both tasks could select."""
    classes = list(pd.read_parquet(cfg.run_dir("cancer_type", "gex", PAN) / "oof_predictions.parquet")
                   .columns.str.removeprefix("p_")[7:])
    rows = []
    for view in cfg.views:
        status, ctype = set(), set()
        for cancer in cancers:
            status |= set(full_importance(cfg, "tumor_status", view, cancer)["feature_id"].head(k))
        for i, _ in enumerate(classes):
            ctype |= set(full_importance(cfg, "cancer_type", view, PAN, class_index=i)["feature_id"].head(k))
        shared_universe = set(pd.read_parquet(cfg.matrix_dir / view / "features.parquet",
                                              columns=["feature_id"])["feature_id"])
        a, b = status & shared_universe, ctype & shared_universe
        both = a & b
        expected = len(a) * len(b) / len(shared_universe)
        rows.append({"view": view, "n_universe": len(shared_universe), "n_tumour_status": len(a),
                     "n_cancer_type": len(b), "n_both": len(both), "expected_by_chance": expected,
                     "fold_enrichment": len(both) / expected if expected else np.nan,
                     "jaccard": len(both) / len(a | b) if (a | b) else np.nan})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["tumor_status", "cancer_type", "overlap", "both"], default="both")
    args = ap.parse_args()
    cfg = load_config()
    seed = cfg.study["splits"]["seed"]
    out = cfg.results_dir / "biology" / "09_ensemble"
    out.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = sorted(scope.loc[scope["tumor_status_in_scope"], "project"])

    if args.task in ("tumor_status", "both"):
        log("tumour status")
        m, t, r = tumour_status(cfg, cancers, seed)
        m.to_csv(out / "tumor_status_model_metrics.csv", index=False)
        t.to_csv(out / "tumor_status_model_tests.csv", index=False)
        r.to_csv(out / "tumor_status_error_rescue.csv", index=False)
        piv = m.pivot_table(index="model", values=["auroc", "n_errors"], aggfunc={"auroc": "mean", "n_errors": "sum"})
        log("mean AUROC / total errors over 15 cancers:\n" + piv.sort_values("n_errors").to_string())

    if args.task in ("overlap", "both"):
        ov = task_overlap(cfg, cancers)
        ov.to_csv(out / "task_feature_overlap.csv", index=False)
        log("features shared by the two tasks:\n" + ov.round(3).to_string(index=False))

    if args.task in ("cancer_type", "both"):
        log("cancer type")
        m, t, pc = cancer_type(cfg, seed)
        m.to_csv(out / "cancer_type_model_metrics.csv", index=False)
        t.to_csv(out / "cancer_type_model_tests.csv", index=False)
        pc.to_csv(out / "cancer_type_recall_by_class.csv", index=False)
        log("\n" + m.round(4).to_string(index=False))
        log("\n" + t.round(6).to_string(index=False))


if __name__ == "__main__":
    main()
