"""Step 12 (B-05): which splicing regulators could drive the events the models rely on?

Within the tumours of each cancer, the expression (log2 CPM) of ~70 RNA-binding proteins is
correlated with the PSI of (a) the top events of that cancer's tumour-status model and (b) a
random set of measured background events of the same view. Correlations are partial Spearman,
adjusting for a stromal and an immune expression score (mean log2 CPM of marker genes), because
tumour composition moves both splicing-factor expression and the PSI of tissue-specific events.

  association excess = share of top events with |partial rho| > threshold
                       - the same share among background events
  paired log2FC      = median over matched pairs of tumour - normal log2 CPM

A regulator is a candidate driver when its expression changes between tumour and normal AND it
tracks the top events beyond background. This is correlation, not causation: a factor can track
the events because it is co-regulated with them.

Outputs: results/biology/05_regulators/ (tables; figures: scripts/15_bio_figures.py)
"""
from __future__ import annotations

import argparse
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import rankdata, wilcoxon

from ctsm.biology import bh, bio_config, feature_values, full_importance, gex_by_symbol, matched_pairs, task_table
from ctsm.config import load_config

warnings.filterwarnings("ignore", category=RuntimeWarning)
RHO = 0.3            # |partial rho| threshold for "tracks"
MIN_FINITE = 0.7     # events need PSI in >= 70% of tumours (the rest median-imputed)


def ranks(A: np.ndarray) -> np.ndarray:
    return np.apply_along_axis(rankdata, 0, A)


def residualise(R: np.ndarray, C: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(C, R, rcond=None)
    return R - C @ beta


def corr_cols(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    A = (A - A.mean(0)) / (A.std(0) + 1e-12)
    B = (B - B.mean(0)) / (B.std(0) + 1e-12)
    return A.T @ B / len(A)


def one_cancer(args) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cancer, = args
    cfg = load_config()
    rc = bio_config()["regulators"]
    sfs = [g for grp in rc["splicing_factors"].values() for g in grp]
    tab = task_table(cfg, "tumor_status", cancer)
    tum = tab[tab["y"] == 1].reset_index(drop=True)
    G = gex_by_symbol(cfg, tab["file_id"], sfs + rc["stromal_genes"] + rc["immune_genes"])
    # paired differential expression of the regulators
    ti, ni, _ = matched_pairs(tab)
    de = []
    for g in [s for s in sfs if s in G]:
        d = G[g].to_numpy()[ti] - G[g].to_numpy()[ni]
        de.append({"cancer": cancer, "gene": g, "paired_log2fc": float(np.median(d)),
                   "p_paired": wilcoxon(d).pvalue if len(d) >= 5 and np.any(d != 0) else np.nan})
    Gt = G.loc[tum["file_id"]]
    cov = np.column_stack([np.ones(len(tum)),
                           rankdata(Gt[[g for g in rc["stromal_genes"] if g in Gt]].mean(1)),
                           rankdata(Gt[[g for g in rc["immune_genes"] if g in Gt]].mean(1))])
    sf_here = [s for s in sfs if s in Gt and Gt[s].std() > 0]
    Rs = residualise(ranks(Gt[sf_here].to_numpy()), cov)
    rng = np.random.default_rng(0)
    long, summ = [], []
    for view in rc["views"]:
        imp = full_importance(cfg, "tumor_status", view, cancer)
        top = imp["feature_id"].head(rc["top_events_per_view"]).tolist()
        pf = pd.read_parquet(cfg.matrix_dir / view / "prefilter.parquet", columns=[cancer])[cancer].to_numpy()
        fall = pd.read_parquet(cfg.matrix_dir / view / "features.parquet", columns=["feature_id", "symbol"])
        cand = np.setdiff1d(fall.loc[pf, "feature_id"].to_numpy(), top)
        bgr = rng.choice(cand, size=min(rc["n_background_events"], len(cand)), replace=False).tolist()
        fids = top + bgr
        X = feature_values(cfg, view, tum["file_id"], fids)
        ok = np.isfinite(X).mean(0) >= MIN_FINITE
        med = np.nanmedian(X, axis=0)
        X = np.where(np.isfinite(X), X, med)
        ok &= np.nanstd(X, axis=0) > 0
        Re = residualise(ranks(X[:, ok]), cov)
        rho = corr_cols(Rs, Re)                              # SF x events
        is_top = np.array([f in set(top) for f in np.array(fids)[ok]])
        sym = fall.set_index("feature_id").loc[np.array(fids)[ok], "symbol"].to_numpy()
        for i, s in enumerate(sf_here):
            r = np.abs(rho[i])
            summ.append({"cancer": cancer, "view": view, "gene": s, "n_top": int(is_top.sum()),
                         "n_background": int((~is_top).sum()),
                         "share_top_tracked": float(np.mean(r[is_top] > RHO)),
                         "share_background_tracked": float(np.mean(r[~is_top] > RHO))})
        for j in np.flatnonzero(is_top):
            i = int(np.argmax(np.abs(rho[:, j])))
            long.append({"cancer": cancer, "view": view, "feature_id": np.array(fids)[ok][j], "symbol": sym[j],
                         "best_regulator": sf_here[i], "partial_rho": float(rho[i, j])})
    return pd.DataFrame(de), pd.DataFrame(summ), pd.DataFrame(long)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    cfg = load_config()
    rc = bio_config()["regulators"]
    out = cfg.results_dir / "biology" / "05_regulators"
    out.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = list(scope.loc[scope.tumor_status_in_scope, "project"])
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(one_cancer, [(c,) for c in cancers]))
    de = pd.concat([r[0] for r in res], ignore_index=True)
    de["q_paired"] = bh(de["p_paired"].fillna(1).to_numpy())
    summ = pd.concat([r[1] for r in res], ignore_index=True)
    best = pd.concat([r[2] for r in res], ignore_index=True)
    summ["excess"] = summ["share_top_tracked"] - summ["share_background_tracked"]
    de.to_csv(out / "regulator_paired_expression.csv", index=False)
    summ.to_csv(out / "regulator_event_association.csv", index=False)
    best.to_csv(out / "top_event_best_regulator.csv", index=False)

    # pooled over views (weighted by number of top events)
    w = summ.assign(t=summ.share_top_tracked * summ.n_top, b=summ.share_background_tracked * summ.n_background)
    pooled = w.groupby(["cancer", "gene"]).agg(t=("t", "sum"), nt=("n_top", "sum"), b=("b", "sum"),
                                               nb=("n_background", "sum")).reset_index()
    pooled["excess"] = pooled.t / pooled.nt - pooled.b / pooled.nb
    pooled.to_csv(out / "regulator_association_excess_pooled.csv", index=False)

    log("regulators done")


if __name__ == "__main__":
    main()
