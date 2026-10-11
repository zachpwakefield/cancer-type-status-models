"""Fold-level preprocessing (decision D-07). Every statistic is computed on training rows only.

1. candidates   AS: observed in >= min_observed_frac of training samples
                GEX: log2(CPM+1) >= gex_min_log2cpm in >= gex_min_expressed_frac of them
                and nonzero training variance (observed values)
2. ranking      "variance": variance of observed training values (unsupervised)
                "anova":    one-way ANOVA F of observed training values across classes
3. top k_max    kept in rank order, so the top-K model for any K <= k_max is a column prefix
4. impute       training median of observed values (after ranking, so imputed values never
                influence which features are chosen)
5. scale        training mean / SD of the imputed values; clipped to +-clip_z
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .matrix_io import ViewMatrix

_CHUNK = 4096


def column_stats(M: ViewMatrix, rows: np.ndarray, y: np.ndarray, n_classes: int, view: str,
                 gex_min_log2cpm: float) -> dict[str, np.ndarray]:
    """Per-column observed count, expressed count (GEX), mean, variance and ANOVA F over `rows`.
    Computed in column chunks so temporaries stay small. y is aligned with rows."""
    n, p = len(rows), M.shape[1]
    out = {k: np.zeros(p) for k in ("n_obs", "n_expressed", "mean", "var", "anova_f")}
    Y = np.zeros((n, n_classes), dtype=np.float64)
    Y[np.arange(n), y] = 1.0
    for a in range(0, p, _CHUNK):
        b = min(p, a + _CHUNK)
        B = M.get(rows, slice(a, b)).astype(np.float64)
        obs = np.isfinite(B)
        B0 = np.where(obs, B, 0.0)
        n_obs = obs.sum(axis=0).astype(np.float64)
        s1 = B0.sum(axis=0)
        s2 = (B0 * B0).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = s1 / n_obs
            var = np.maximum(s2 / n_obs - mean * mean, 0.0)
            Nk = Y.T @ obs.astype(np.float64)                  # classes x cols
            mk = np.where(Nk > 0, (Y.T @ B0) / np.maximum(Nk, 1), 0.0)
            ssb = (Nk * (mk - mean[None, :]) ** 2).sum(axis=0)
            ssw = np.maximum(s2 - (Nk * mk * mk).sum(axis=0), 0.0)
            g = (Nk > 0).sum(axis=0)
            df_b, df_w = g - 1, n_obs - g
            f = (ssb / np.maximum(df_b, 1)) / (ssw / np.maximum(df_w, 1))
        f = np.where((df_b >= 1) & (df_w >= 1), f, 0.0)
        f = np.where(np.isnan(f), 0.0, f)                      # 0/0: no signal
        f = np.where(np.isposinf(f), np.finfo(np.float64).max, f)   # perfect separation
        out["n_obs"][a:b], out["mean"][a:b], out["var"][a:b], out["anova_f"][a:b] = n_obs, mean, var, f
        if view == "gex":
            out["n_expressed"][a:b] = (B0 >= gex_min_log2cpm).sum(axis=0)
    return out


def candidate_mask(stats: dict, n_rows: int, view: str, pp: dict) -> np.ndarray:
    if view == "gex":
        support = stats["n_expressed"] >= pp["gex_min_expressed_frac"] * n_rows
    else:
        support = stats["n_obs"] >= pp["min_observed_frac"] * n_rows
    return support & (stats["var"] > 0)


@dataclass
class FoldPreprocessor:
    cols: np.ndarray          # column indices into the job matrix, best-ranked first
    score: np.ndarray         # ranking score of each kept column
    median: np.ndarray
    mean: np.ndarray
    sd: np.ndarray
    clip_z: float
    method: str
    n_candidates: int

    def transform(self, M: ViewMatrix, rows: np.ndarray, k: int | None = None) -> np.ndarray:
        k = len(self.cols) if k is None else min(k, len(self.cols))
        Z = M.get(rows, self.cols[:k])
        Z = np.where(np.isfinite(Z), Z, self.median[:k][None, :])
        Z = (Z - self.mean[:k]) / self.sd[:k]
        return np.clip(Z, -self.clip_z, self.clip_z).astype(np.float32)

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in ("cols", "score", "median", "mean", "sd")} | {
            "clip_z": self.clip_z, "method": self.method, "n_candidates": self.n_candidates}


def fit_preprocessor(M: ViewMatrix, rows: np.ndarray, y: np.ndarray, view: str, method: str, k_max: int, pp: dict,
                     n_classes: int, stats: dict | None = None) -> FoldPreprocessor:
    """Fit on training `rows` only (y aligned with rows). `stats` from column_stats on the same
    rows may be passed in to share it between ranking methods."""
    if stats is None:
        stats = column_stats(M, rows, y, n_classes, view, pp["gex_min_log2cpm"])
    cand = np.flatnonzero(candidate_mask(stats, len(rows), view, pp))
    score = stats["var" if method == "variance" else "anova_f"][cand]
    order = np.lexsort((cand, -score))                   # descending score, ties by column index
    cols = cand[order[:k_max]]
    sub = M.get(rows, cols)
    med = np.nanmedian(sub, axis=0).astype(np.float32)
    filled = np.where(np.isfinite(sub), sub, med[None, :])
    mean = filled.mean(axis=0).astype(np.float32)
    sd = filled.std(axis=0).astype(np.float32)
    sd[sd <= 0] = 1.0
    return FoldPreprocessor(cols=cols, score=score[order[:k_max]], median=med, mean=mean, sd=sd,
                            clip_z=float(pp["clip_z"]), method=method, n_candidates=int(len(cand)))
