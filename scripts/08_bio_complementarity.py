"""Step 08 (B-01): do alternative-RNA-processing views carry information that expression lacks?

Uses only out-of-fold probabilities (every sample predicted by a model that never saw it). The
fusion weights are fixed in config/biology.yaml (no tuning): 'all_equal' = mean of the 9
views; 'gex_as_half' = 0.5 GEX + 0.5 mean of the 8 AS views.

Outputs: results/biology/01_complementarity/ (tables; figures: scripts/15_bio_figures.py)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import binomtest
from sklearn.metrics import balanced_accuracy_score, roc_auc_score

from ctsm.config import PAN, load_config
from ctsm.figstyle import VIEW_NAMES


def oof(cfg, task, view, cohort) -> pd.DataFrame:
    return pd.read_parquet(cfg.run_dir(task, view, cohort) / "oof_predictions.parquet")


def mcnemar(a_correct: np.ndarray, b_correct: np.ndarray) -> tuple[int, int, float]:
    """Exact McNemar: b01 = only A correct, b10 = only B correct."""
    b01 = int((a_correct & ~b_correct).sum())
    b10 = int((~a_correct & b_correct).sum())
    p = binomtest(b01, b01 + b10, 0.5).pvalue if b01 + b10 else 1.0
    return b01, b10, float(p)


def boot_acc(y, pred, groups, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    u, inv = np.unique(groups, return_inverse=True)
    mem = [np.flatnonzero(inv == g) for g in range(len(u))]
    acc, bac = [], []
    for _ in range(n):
        idx = np.concatenate([mem[g] for g in rng.integers(0, len(u), len(u))])
        acc.append(np.mean(y[idx] == pred[idx]))
        bac.append(balanced_accuracy_score(y[idx], pred[idx]))
    return np.percentile(acc, [2.5, 97.5]), np.percentile(bac, [2.5, 97.5])


def main() -> None:
    cfg = load_config()
    out = cfg.results_dir / "biology" / "01_complementarity"
    out.mkdir(parents=True, exist_ok=True)
    views = cfg.views
    as_views = [v for v in views if v != "gex"]

    # ---------------- cancer type ----------------
    P = {v: oof(cfg, "cancer_type", v, PAN) for v in views}
    base = P["gex"]
    for v in views:
        assert (P[v]["file_id"].to_numpy() == base["file_id"].to_numpy()).all()
    pc = [c for c in base.columns if c.startswith("p_")]
    classes = [c[2:] for c in pc]
    y, groups = base["y"].to_numpy(), base["case_id"].to_numpy()
    M = {v: P[v][pc].to_numpy() for v in views}
    M["AS ensemble"] = np.mean([M[v] for v in as_views], axis=0)
    M["all 9 (equal)"] = np.mean([M[v] for v in views], axis=0)
    M["GEX + AS (1:1)"] = 0.5 * M["gex"] + 0.5 * M["AS ensemble"]
    rows = []
    for name, pr in M.items():
        pred = pr.argmax(1)
        (alo, ahi), (blo, bhi) = boot_acc(y, pred, groups)
        rows.append({"model": VIEW_NAMES.get(name, name), "accuracy": np.mean(pred == y), "acc_lo": alo, "acc_hi": ahi,
                     "balanced_accuracy": balanced_accuracy_score(y, pred), "bacc_lo": blo, "bacc_hi": bhi})
    perf = pd.DataFrame(rows)
    perf.to_csv(out / "cancer_type_fusion_metrics.csv", index=False)
    correct = {k: v.argmax(1) == y for k, v in M.items()}
    tests = []
    for a, b in [("gex", "AS ensemble"), ("gex", "GEX + AS (1:1)"), ("gex", "all 9 (equal)"), ("se", "gex")]:
        b01, b10, p = mcnemar(correct[a], correct[b])
        tests.append({"A": VIEW_NAMES.get(a, a), "B": b, "only_A_correct": b01, "only_B_correct": b10, "p_mcnemar": p})
    pd.DataFrame(tests).to_csv(out / "cancer_type_mcnemar.csv", index=False)
    gex_err = ~correct["gex"]
    rescue = {"gex_errors": int(gex_err.sum()), "rescued_by_AS_ensemble": int((correct["AS ensemble"] & gex_err).sum()),
              "AS_errors": int((~correct["AS ensemble"]).sum()),
              "rescued_by_GEX": int((correct["gex"] & ~correct["AS ensemble"]).sum())}
    pd.Series(rescue).to_csv(out / "cancer_type_error_rescue.csv")
    log(f"cancer type: {rescue}")

    # per-class recall: GEX vs AS ensemble vs fused
    rec = pd.DataFrame({name: [np.mean(M[key][y == k].argmax(1) == k) for k in range(len(classes))]
                        for name, key in [("GEX", "gex"), ("AS ensemble", "AS ensemble"), ("GEX + AS", "GEX + AS (1:1)")]},
                       index=classes)
    rec["n"] = np.bincount(y, minlength=len(classes))
    rec["AS_minus_GEX"] = rec["AS ensemble"] - rec["GEX"]
    rec.to_csv(out / "cancer_type_recall_gex_vs_as.csv")

    # error overlap between views: Jaccard of misclassified sample sets
    err = pd.DataFrame({VIEW_NAMES[v]: ~correct[v] for v in views})
    jac = pd.DataFrame(index=err.columns, columns=err.columns, dtype=float)
    for a in err.columns:
        for b in err.columns:
            jac.loc[a, b] = (err[a] & err[b]).sum() / max((err[a] | err[b]).sum(), 1)
    jac.to_csv(out / "cancer_type_error_jaccard.csv")

    # ---------------- tumour status ----------------
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = list(scope.loc[scope.tumor_status_in_scope, "project"])
    trows = []
    for c in cancers:
        Q = {v: oof(cfg, "tumor_status", v, c) for v in views}
        yy = Q["gex"]["y"].to_numpy()
        pt = {v: Q[v]["p_Tumour"].to_numpy() for v in views}
        pt["AS ensemble"] = np.mean([pt[v] for v in as_views], axis=0)
        pt["GEX + AS (1:1)"] = 0.5 * pt["gex"] + 0.5 * pt["AS ensemble"]
        for name, p in pt.items():
            trows.append({"cancer": c, "model": VIEW_NAMES.get(name, name), "auroc": roc_auc_score(yy, p),
                          "balanced_accuracy": balanced_accuracy_score(yy, (p >= 0.5).astype(int)),
                          "n_errors": int(((p >= 0.5).astype(int) != yy).sum())})
    ts = pd.DataFrame(trows)
    ts.to_csv(out / "tumor_status_fusion_metrics.csv", index=False)

    log("complementarity done")


if __name__ == "__main__":
    main()
