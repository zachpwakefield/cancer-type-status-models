"""Nested, case-grouped, label-stratified cross-validation (decision D-04).

Outer: StratifiedGroupKFold(outer_folds). Every sample is in exactly one outer test fold and
all samples of a case (its tumour and its normal) are in the same fold.
Inner: StratifiedGroupKFold(inner_folds) inside each outer-train set, used only for tuning.

Stored as one table per (task, cohort): file_id, case_id, y, outer, inner_o0..inner_o{K-1}
(inner fold of the sample when outer fold k is held out; -1 for the held-out samples).
"""
from __future__ import annotations

import zlib

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold


def _seed(base: int, *keys) -> int:
    return (int(base) + zlib.crc32("|".join(map(str, keys)).encode())) % (2**31 - 1)


def make_splits(table: pd.DataFrame, task: str, cohort: str, n_outer: int, n_inner: int, base_seed: int) -> pd.DataFrame:
    y = table["y"].to_numpy()
    groups = table["case_id"].to_numpy()
    out = table[["file_id", "case_id", "y"]].copy()
    out["outer"] = -1
    outer = StratifiedGroupKFold(n_splits=n_outer, shuffle=True, random_state=_seed(base_seed, task, cohort, "outer"))
    for k, (_, te) in enumerate(outer.split(np.zeros(len(y)), y, groups)):
        out.iloc[te, out.columns.get_loc("outer")] = k
    for k in range(n_outer):
        col = f"inner_o{k}"
        out[col] = -1
        tr = np.flatnonzero(out["outer"].to_numpy() != k)
        inner = StratifiedGroupKFold(n_splits=n_inner, shuffle=True,
                                     random_state=_seed(base_seed, task, cohort, "inner", k))
        for j, (_, va) in enumerate(inner.split(np.zeros(len(tr)), y[tr], groups[tr])):
            out.iloc[tr[va], out.columns.get_loc(col)] = j
    check_splits(out, n_outer, n_inner)
    return out


def check_splits(s: pd.DataFrame, n_outer: int, n_inner: int) -> None:
    if (s["outer"] < 0).any():
        raise AssertionError("sample without an outer fold")
    if s.groupby("case_id")["outer"].nunique().max() != 1:
        raise AssertionError("a case is split across outer folds")
    for k in range(n_outer):
        col = s[f"inner_o{k}"]
        if not ((col == -1) == (s["outer"] == k)).all():
            raise AssertionError(f"inner_o{k}: held-out samples must be exactly outer fold {k}")
        tr = s[s["outer"] != k]
        if tr.groupby("case_id")[f"inner_o{k}"].nunique().max() != 1:
            raise AssertionError(f"inner_o{k}: a case is split across inner folds")
        if set(tr[f"inner_o{k}"]) != set(range(n_inner)):
            raise AssertionError(f"inner_o{k}: empty inner fold")


def fold_indices(s: pd.DataFrame, outer: int, inner: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Row positions (train, held-out). inner=None -> (outer-train, outer-test);
    otherwise (inner-train, inner-validation) inside outer-train of `outer`."""
    if inner is None:
        o = s["outer"].to_numpy()
        return np.flatnonzero(o != outer), np.flatnonzero(o == outer)
    col = s[f"inner_o{outer}"].to_numpy()
    return np.flatnonzero((col >= 0) & (col != inner)), np.flatnonzero(col == inner)


def min_inner_train(s: pd.DataFrame, n_outer: int, n_inner: int) -> int:
    """Smallest training set any fold-level filter is fitted on (an inner-train set)."""
    return int(min(len(fold_indices(s, o, i)[0]) for o in range(n_outer) for i in range(n_inner)))
