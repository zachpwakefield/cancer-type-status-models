"""Classification metrics with case-level bootstrap confidence intervals (decision D-11).

All metrics are computed on out-of-fold predictions: every sample is predicted by the outer-fold
model that never saw it (or its case) during preprocessing, tuning or training.
"""
from __future__ import annotations

import numpy as np
from sklearn import metrics as skm


def binary_metrics(y: np.ndarray, p_tumour: np.ndarray) -> dict:
    pred = (p_tumour >= 0.5).astype(int)
    eps = 1e-12
    return {
        "auroc": float(skm.roc_auc_score(y, p_tumour)),
        "auprc_tumour": float(skm.average_precision_score(y, p_tumour)),
        "auprc_normal": float(skm.average_precision_score(1 - y, 1 - p_tumour)),
        "balanced_accuracy": float(skm.balanced_accuracy_score(y, pred)),
        "sensitivity_tumour": float(skm.recall_score(y, pred, pos_label=1, zero_division=0)),
        "specificity_normal": float(skm.recall_score(y, pred, pos_label=0, zero_division=0)),
        "f1_normal": float(skm.f1_score(y, pred, pos_label=0, zero_division=0)),
        "mcc": float(skm.matthews_corrcoef(y, pred)),
        "accuracy": float(skm.accuracy_score(y, pred)),
        "brier": float(np.mean((p_tumour - y) ** 2)),
        "log_loss": float(skm.log_loss(y, np.clip(p_tumour, eps, 1 - eps), labels=[0, 1])),
        "n": int(len(y)), "n_normal": int((y == 0).sum()), "n_tumour": int((y == 1).sum()),
    }


def multiclass_metrics(y: np.ndarray, proba: np.ndarray) -> dict:
    C = proba.shape[1]
    pred = proba.argmax(axis=1)
    top3 = np.argsort(-proba, axis=1)[:, :3]
    present = np.unique(y)
    return {
        "accuracy": float(skm.accuracy_score(y, pred)),
        "balanced_accuracy": float(skm.balanced_accuracy_score(y, pred)),
        "macro_f1": float(skm.f1_score(y, pred, average="macro", labels=present, zero_division=0)),
        "weighted_f1": float(skm.f1_score(y, pred, average="weighted", zero_division=0)),
        "mcc": float(skm.matthews_corrcoef(y, pred)),
        "top3_accuracy": float(np.mean((top3 == y[:, None]).any(axis=1))),
        "macro_auroc_ovr": float(skm.roc_auc_score(y, proba, multi_class="ovr", average="macro", labels=np.arange(C))),
        "log_loss": float(skm.log_loss(y, np.clip(proba, 1e-12, 1), labels=np.arange(C))),
        "n": int(len(y)), "n_classes": int(C),
    }


def per_class_table(y: np.ndarray, proba: np.ndarray, class_names: list[str]) -> list[dict]:
    pred = proba.argmax(axis=1)
    p, r, f, s = skm.precision_recall_fscore_support(y, pred, labels=np.arange(len(class_names)), zero_division=0)
    rows = []
    for c, name in enumerate(class_names):
        yc = (y == c).astype(int)
        auc = float(skm.roc_auc_score(yc, proba[:, c])) if 0 < yc.sum() < len(yc) else float("nan")
        rows.append({"class": name, "precision": float(p[c]), "recall": float(r[c]), "f1": float(f[c]),
                     "auroc_ovr": auc, "support": int(s[c])})
    return rows


def bootstrap_ci(y: np.ndarray, proba: np.ndarray, groups: np.ndarray, binary: bool, keys: list[str],
                 n_boot: int = 1000, seed: int = 0) -> dict:
    """Percentile 95 % CIs from resampling CASES with replacement (a case's samples move together)."""
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(groups, return_inverse=True)
    members = [np.flatnonzero(inv == g) for g in range(len(uniq))]
    vals = {k: [] for k in keys}
    for _ in range(n_boot):
        pick = rng.integers(0, len(uniq), len(uniq))
        idx = np.concatenate([members[g] for g in pick])
        yb, pb = y[idx], proba[idx]
        if len(np.unique(yb)) < 2:
            continue
        try:
            m = binary_metrics(yb, pb[:, 1]) if binary else multiclass_metrics(yb, pb)
        except ValueError:
            continue
        for k in keys:
            vals[k].append(m[k])
    return {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v else [None, None]
            for k, v in vals.items()}
