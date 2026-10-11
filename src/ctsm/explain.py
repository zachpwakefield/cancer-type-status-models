"""SHAP explanations (decision D-10).

Per outer fold, the refitted model is explained with shap.DeepExplainer (DeepLIFT-style SHAP
for deep networks) on EVERY held-out outer-test sample, against a background of random
outer-train samples. One scalar output is explained per sample, in logit units:
  tumor_status  the tumour log-odds (the model's single logit)
  cancer_type   the centred logit of the sample's TRUE class, z_k - mean_j z_j. Softmax is
                invariant to adding a constant to every logit; centring removes that free
                direction, so attributions measure evidence for class k relative to the others.
                It is a fixed linear map of the logits, so DeepExplainer stays exact.
Additivity (sum of SHAP values + expected output = output) is asserted for every sample.

Aggregation across folds is by feature ID. A feature not selected by a fold's model has an
attribution of exactly 0 for that fold's samples, so the pooled importance rewards features
that every fold relies on.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import shap
import torch
from scipy.stats import spearmanr
from torch import nn


def centred_logit_head(model: nn.Module, k: int, n_classes: int) -> nn.Module:
    """model.net followed by a fixed Linear(C, 1) computing z_k - mean_j z_j."""
    head = nn.Linear(n_classes, 1, bias=False)
    with torch.no_grad():
        w = torch.full((1, n_classes), -1.0 / n_classes)
        w[0, k] += 1.0
        head.weight.copy_(w)
    head.requires_grad_(False)
    return nn.Sequential(*model.net, head).eval()


def _explain(net: nn.Module, bg: torch.Tensor, X: np.ndarray, batch: int = 64) -> tuple[np.ndarray, float]:
    explainer = shap.DeepExplainer(net, bg)
    parts = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for a in range(0, len(X), batch):
            xb = torch.from_numpy(np.ascontiguousarray(X[a:a + batch], dtype=np.float32))
            sv = np.asarray(explainer.shap_values(xb, check_additivity=True))
            parts.append(sv.reshape(len(xb), X.shape[1]).astype(np.float32))
    return np.concatenate(parts, axis=0), float(np.ravel(explainer.expected_value)[0])


def explained_output(model: nn.Module, X: np.ndarray, y: np.ndarray, n_classes: int) -> np.ndarray:
    """The scalar that is explained, per sample (see module docstring)."""
    model.eval()
    with torch.no_grad():
        z = model(torch.from_numpy(np.ascontiguousarray(X, dtype=np.float32))).numpy().astype(np.float64)
    if n_classes == 2:
        return z[:, 0]
    return z[np.arange(len(y)), y] - z.mean(axis=1)


def deep_shap(model: nn.Module, X_background: np.ndarray, X_explain: np.ndarray, y_explain: np.ndarray,
              n_classes: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (SHAP values [n, p], expected output per sample [n])."""
    model.eval()
    bg = torch.from_numpy(np.ascontiguousarray(X_background, dtype=np.float32))
    sv = np.zeros(X_explain.shape, dtype=np.float32)
    ev = np.zeros(len(X_explain))
    if n_classes == 2:
        sv[:], ev[:] = _explain(model, bg, X_explain)
        return sv, ev
    for k in np.unique(y_explain):
        idx = np.flatnonzero(y_explain == k)
        sv[idx], ev[idx] = _explain(centred_logit_head(model, int(k), n_classes), bg, X_explain[idx])
    return sv, ev


def additivity_error(model: nn.Module, X: np.ndarray, y: np.ndarray, n_classes: int, sv: np.ndarray,
                     ev: np.ndarray) -> float:
    """Max |sum(SHAP) + E[f] - f(x)| over samples, in logit units."""
    return float(np.abs(sv.sum(axis=1) + ev - explained_output(model, X, y, n_classes)).max())


# ---- aggregation over folds -----------------------------------------------------------------

def _direction(values: list[np.ndarray], shaps: list[np.ndarray]) -> float:
    v = np.concatenate(values) if values else np.array([])
    s = np.concatenate(shaps) if shaps else np.array([])
    if len(v) < 3 or np.ptp(v) == 0 or np.ptp(s) == 0:
        return float("nan")
    return float(spearmanr(v, s).statistic)


def _accumulate(folds: list[dict], mask_of=lambda f: slice(None)) -> tuple[dict, dict, dict, int]:
    """Sum |SHAP|, sum SHAP and fold counts per feature over the selected samples of each fold."""
    abs_sum, sgn_sum, n_sel, n_total = {}, {}, {}, 0
    for f in folds:
        S = f["shap"][mask_of(f)].astype(np.float64)
        n_total += len(S)
        a, s = np.abs(S).sum(axis=0), S.sum(axis=0)
        for j, fid in enumerate(f["feature_ids"]):
            abs_sum[fid] = abs_sum.get(fid, 0.0) + float(a[j])
            sgn_sum[fid] = sgn_sum.get(fid, 0.0) + float(s[j])
            n_sel[fid] = n_sel.get(fid, 0) + 1
    return abs_sum, sgn_sum, n_sel, n_total


def importance_table(folds: list[dict], features: pd.DataFrame, mask_of=lambda f: slice(None),
                     n_direction: int | None = None) -> pd.DataFrame:
    """folds: [{feature_ids, index {feature_id: column}, X (standardised), shap [n, p], y}].
    One row per feature ever selected: mean |SHAP| and mean SHAP over ALL selected samples of
    all folds (0 where not selected), folds selected, and direction = Spearman(value, SHAP)
    over the samples where it was selected (> 0: higher values push the explained output up).
    Direction is computed for the top `n_direction` rows (all if None)."""
    abs_sum, sgn_sum, n_sel, n_total = _accumulate(folds, mask_of)
    fids = list(abs_sum)
    t = pd.DataFrame({"feature_id": fids, "mean_abs_shap": [abs_sum[x] / max(n_total, 1) for x in fids],
                      "mean_shap": [sgn_sum[x] / max(n_total, 1) for x in fids],
                      "n_folds_selected": [n_sel[x] for x in fids]})
    t = t.sort_values(["mean_abs_shap", "feature_id"], ascending=[False, True]).reset_index(drop=True)
    t.insert(0, "rank", np.arange(1, len(t) + 1))
    top = len(t) if n_direction is None else min(n_direction, len(t))
    dirs = np.full(len(t), np.nan)
    for i, fid in enumerate(t["feature_id"].iloc[:top]):
        vals, shs = [], []
        for f in folds:
            pos = f["index"].get(fid)
            if pos is not None:
                m = mask_of(f)
                vals.append(f["X"][m, pos])
                shs.append(f["shap"][m, pos].astype(np.float64))
        dirs[i] = _direction(vals, shs)
    t["direction_spearman"] = dirs
    return t.merge(features, on="feature_id", how="left")


def gene_table(feat_table: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Sum of feature importances per host gene (several AS events can share a gene)."""
    g = (feat_table.groupby(["gene_id", "symbol"], dropna=False)
         .agg(importance=(value_col, "sum"), n_features=("feature_id", "size"),
              top_feature=("label", "first"))
         .reset_index().sort_values("importance", ascending=False).reset_index(drop=True))
    g.insert(0, "rank", np.arange(1, len(g) + 1))
    return g


def class_importance_matrix(folds: list[dict], class_names: list[str]) -> pd.DataFrame:
    """features x classes: mean |SHAP| over explained samples of each class (0 where the feature
    was not selected by a fold). Multiclass runs explain each sample's own class, so column k
    is the importance of each feature for recognising class k."""
    C = len(class_names)
    n_class = np.zeros(C)
    sums: dict[str, np.ndarray] = {}
    for f in folds:
        n_class += np.bincount(f["y"], minlength=C)
        per = np.zeros((C, f["shap"].shape[1]))
        np.add.at(per, f["y"], np.abs(f["shap"].astype(np.float64)))
        for j, fid in enumerate(f["feature_ids"]):
            sums[fid] = sums.get(fid, np.zeros(C)) + per[:, j]
    fids = list(sums)
    return pd.DataFrame(np.vstack([sums[x] for x in fids]) / np.maximum(n_class, 1), index=fids, columns=class_names)


def load_folds(run_dir, y: np.ndarray, n_folds: int) -> list[dict]:
    """Read <run_dir>/fold*/shap.npz into the structure used by the aggregation functions."""
    folds = []
    for k in range(n_folds):
        z = np.load(run_dir / f"fold{k}" / "shap.npz")
        fids = z["feature_ids"].astype(str)
        folds.append({"feature_ids": fids, "index": {f: i for i, f in enumerate(fids)},
                      "X": z["X"].astype(np.float32), "shap": z["shap"].astype(np.float32), "y": y[z["rows"]]})
    return folds
