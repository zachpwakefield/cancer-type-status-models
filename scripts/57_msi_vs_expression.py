"""Step 57 (D-74): does the splicing model transfer across cancers better than expression models
trained the same way?

The repeat-exon logistic model (D-71 M4; A-19) was trained on COAD, STAD and UCEC and tested on the
30 other cancers (E1). Here expression models get the same tumours, the same label (MC3 >= 10
homopolymer deletions), the same split and the same tuning (C, and the number of genes, chosen by
leaving one training cancer out): X1 all genes as within-cancer percentiles (M4's input scale),
X2 all genes as log2(CPM + 1), X3 PreMSIm's 15 genes, X0 MLH1 alone (unfitted), and C1 the repeat
exons plus genes. EX1-EX6 compare held-out AUROC by a bootstrap over tumours within cancer x label.
config/followup.yaml: msi_vs_expression.

Inputs:  data.nosync/ctsm_v1/msi/models_{tumours,e1_predictions,e2_predictions,exon_matrix,exon_columns}
         (step 54), the GEX matrix, the D-69 caller table (all local)
Outputs: results/subtype/msi_cis/code/expression/ (aggregates); data.nosync/ctsm_v1/msi/
         expr_e1_predictions.parquet and expr_e2_predictions.parquet (per tumour)
Usage: python scripts/57_msi_vs_expression.py [--only fit|e2|report]
"""
from __future__ import annotations

import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import rankdata


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m54 = _load("msi54", "54_msi_models.py")
m55 = _load("msi55", "55_msi_benchmark.py")
m51 = m54.m51
LOCAL, OUT, TRAIN = m51.LOCAL, m51.OUT / "expression", m54.TRAIN
# PreMSIm's panel (Li et al. 2020; github.com/WangX-Lab/PreMSIm README, read 2026-10-13); TTC30A is now IFT70A
PREMSIM = ["DDX27", "EPM2AIP1", "HENMT1", "LYG1", "MLH1", "MSH4", "NHLRC1", "NOL4L", "RNLS", "RPL22L1",
           "RTF2", "SHROOM4", "SMAP1", "IFT70A", "ZSWIM3"]
C_GRID = [1e-4, 1e-3, 1e-2, 1e-1]
GRIDS = {"X1": [{"K": k, "C": c} for k in [500, 2000, None] for c in C_GRID],
         "X2": [{"K": k, "C": c} for k in [500, 2000, None] for c in C_GRID],
         "X3": [{"K": None, "C": c} for c in [1e-3, 1e-2, 1e-1, 1, 10]],
         "C1": [{"K": k, "C": c} for k in [500, 2000] for c in C_GRID]}
MLH1_LOW = 0.10                                     # EX5: MLH1 percentile within cancer below this
MIN_EX5 = 10
NAMES = {"M4": "splicing (repeat exons)", "X1": "expression, percentiles", "X2": "expression, log CPM",
         "X3": "PreMSIm 15 genes", "X3k": "PreMSIm recipe (k-NN)", "X0": "MLH1 alone",
         "C1": "repeat exons + genes"}


# ---------------------------------------------------------------- data

def tumours() -> pd.DataFrame:
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")[["file_id", "sample15", "project", "n_hp", "mmr"]]
    return t


def expression(t: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """log2(CPM + 1) and within-cancer percentiles (float32, tumours x genes) and the gene table."""
    from ctsm.config import load_config
    from ctsm.matrix_io import load_view
    M, F = load_view(load_config().matrix_dir, "gex", t.file_id)
    L = np.asarray(M.get(), dtype=np.float32)
    P = np.empty_like(L)
    for _, idx in t.groupby("project").indices.items():
        P[idx] = (rankdata(L[idx], axis=0) / len(idx)).astype(np.float32)
    return L, P, F.reset_index(drop=True)


# ---------------------------------------------------------------- models

def fit_predict(X, tr, te, y, params, order=None, always=None):
    """Standardise on tr, L2 logistic regression with balanced weights, decision function on te.
    `always` (columns of X) are kept in addition to the top-K of `order`."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    k = params.get("K")
    if always is None:
        cols = np.arange(X.shape[1]) if k is None else np.sort(order[:k])
    else:
        cols = np.concatenate([always, np.sort(order[:k])])
    A = np.asarray(X[np.ix_(tr, cols)], dtype=np.float64)
    sc = StandardScaler().fit(A)
    m = LogisticRegression(C=params["C"], class_weight="balanced", max_iter=5000).fit(sc.transform(A), y[tr])
    out = np.empty(len(te))
    for a in range(0, len(te), 1000):
        B = sc.transform(np.asarray(X[np.ix_(te[a:a + 1000], cols)], dtype=np.float64))
        out[a:a + 1000] = m.decision_function(B)
    return out, m, cols


def premsim_knn(L: np.ndarray, y: np.ndarray, tr: np.ndarray, te: np.ndarray) -> np.ndarray:
    """PreMSIm's own recipe on its genes (exploratory): log2 expression scaled to [0, 1] per gene over
    the training tumours, 5 nearest neighbours; score = share of MMR-deficient neighbours."""
    from sklearn.neighbors import KNeighborsClassifier
    lo, hi = L[tr].min(axis=0), L[tr].max(axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    m = KNeighborsClassifier(n_neighbors=5).fit((L[tr] - lo) / span, y[tr])
    return m.predict_proba((L[te] - lo) / span)[:, 1]


def run(name, X, y, groups, tr, te, offset=0, always=None):
    """Tune by leaving one training cancer out (pooled out-of-fold AUROC), refit on tr, score te.
    Top-K selection ranks columns offset.. of X (ANOVA F on the fitting rows); `always` are kept."""
    from sklearn.model_selection import LeaveOneGroupOut
    needs = any(p.get("K") for p in GRIDS[name])

    def order_on(rows):
        return offset + m54.anova_order(X[:, offset:], rows, y[rows]) if needs else None

    folds = list(LeaveOneGroupOut().split(tr, y[tr], groups[tr]))
    orders = [order_on(tr[a]) for a, _ in folds]
    table, best = [], None
    for params in GRIDS[name]:
        oof = np.full(len(tr), np.nan)
        for (a, b), od in zip(folds, orders):
            oof[b] = fit_predict(X, tr[a], tr[b], y, params, od, always)[0]
        au = m54.auc(oof, y[tr])
        table.append({"model": name, "K": params.get("K") or "all", "C": params["C"], "inner_auroc": au})
        if best is None or au > best[1]:
            best = (params, au)
    order = order_on(tr)
    score, m, cols = fit_predict(X, tr, te, y, best[0], order, always)
    for row in table:
        row["chosen"] = row["K"] == (best[0].get("K") or "all") and row["C"] == best[0]["C"]
    log(f"E1 {name}: chosen {best[0]}, inner AUROC {best[1]:.3f}, held-out AUROC {m54.auc(score, y[te]):.3f}")
    return {"score": score, "params": best[0], "inner": best[1], "table": table, "model": m, "cols": cols,
            "order": order}


def fit() -> None:
    t = tumours()
    y, groups = t.mmr.to_numpy(), t.project.to_numpy()
    held = ~t.project.isin(TRAIN).to_numpy()
    tr, te = np.flatnonzero(~held), np.flatnonzero(held)
    log(f"E1: train {len(tr):,} ({int(y[tr].sum())} MMR-deficient), test {len(te):,} ({int(y[te].sum())})")
    L, P, F = expression(t)
    sym = F.symbol.to_numpy()
    first = ~pd.Series(sym).duplicated().to_numpy()
    log(f"expression: {P.shape[1]:,} genes")
    pred = pd.DataFrame({"file_id": t.file_id.to_numpy()[te], "project": groups[te], "mmr": y[te]})
    e1 = pd.read_parquet(LOCAL / "models_e1_predictions.parquet")
    pred = pred.merge(e1[["file_id", "M4"]], on="file_id", how="left")
    mlh1 = int(np.flatnonzero((sym == "MLH1") & first)[0])
    pred["MLH1_pct"] = P[te, mlh1]
    pred["X0"] = 1 - P[te, mlh1]
    res, tables = {}, []
    res["X1"] = run("X1", P, y, groups, tr, te)
    res["X2"] = run("X2", L, y, groups, tr, te)
    pm = np.array([int(np.flatnonzero((sym == g) & first)[0]) for g in PREMSIM if ((sym == g) & first).any()])
    log(f"PreMSIm genes present: {len(pm)} of {len(PREMSIM)}")
    res["X3"] = run("X3", P[:, pm], y, groups, tr, te)
    pred["X3k"] = premsim_knn(L[:, pm], y, tr, te)
    del L
    cols = pd.read_parquet(LOCAL / "models_exon_columns.parquet")
    R = np.load(LOCAL / "models_exon_matrix.npy", mmap_mode="r")[:, (cols.acc_T >= 10).to_numpy()]
    nr = R.shape[1]
    C = np.hstack([np.asarray(R, dtype=np.float32), P])
    del P
    res["C1"] = run("C1", C, y, groups, tr, te, offset=nr, always=np.arange(nr))
    for k, r in res.items():
        pred[k] = r["score"]
        tables += r["table"]
    OUT.mkdir(parents=True, exist_ok=True)
    pred.to_parquet(LOCAL / "expr_e1_predictions.parquet", index=False)
    pd.DataFrame(tables).to_csv(OUT / "tuning_e1.csv", index=False, float_format="%.4g")
    # coefficients of the final expression fits (standardised inputs)
    co = []
    for k, names in [("X1", sym), ("X3", sym[pm])]:
        r = res[k]
        c = r["model"].coef_[0]
        co.append(pd.DataFrame({"model": k, "gene": names[r["cols"]] if k == "X1" else names, "coefficient": c}))
    r = res["C1"]
    lab = np.array(["repeat exon"] * nr + ["gene"] * (C.shape[1] - nr))
    co.append(pd.DataFrame({"model": "C1", "gene": lab[r["cols"]], "coefficient": r["model"].coef_[0]}))
    co = pd.concat(co, ignore_index=True)
    co.to_parquet(LOCAL / "expr_coefficients.parquet", index=False)
    pd.Series({k: str(r["params"]) for k, r in res.items()}).to_json(LOCAL / "expr_params.json")
    log("fit done")


def e2() -> None:
    """Exploratory: X1 leaving each cancer out, with E1's chosen parameters (no inner tuning)."""
    import ast
    import json
    t = tumours()
    y, groups = t.mmr.to_numpy(), t.project.to_numpy()
    _, P, _ = expression(t)
    params = ast.literal_eval(json.loads((LOCAL / "expr_params.json").read_text())["X1"])
    sc = np.full(len(t), np.nan)
    for c in sorted(set(groups)):
        te, tr = np.flatnonzero(groups == c), np.flatnonzero(groups != c)
        order = m54.anova_order(P, tr, y[tr]) if params.get("K") else None
        sc[te] = fit_predict(P, tr, te, y, params, order)[0]
        log(f"E2 X1: {c} done")
    e2m = pd.read_parquet(LOCAL / "models_e2_predictions.parquet")[["file_id", "M4"]]
    out = pd.DataFrame({"file_id": t.file_id, "project": groups, "mmr": y, "X1": sc}).merge(e2m, on="file_id")
    out.to_parquet(LOCAL / "expr_e2_predictions.parquet", index=False)
    log(f"E2 pooled AUROC: X1 {m54.auc(out.X1.to_numpy(), y):.3f}, M4 {m54.auc(out.M4.to_numpy(), y):.3f}")


# ---------------------------------------------------------------- report

def report() -> None:
    pred = pd.read_parquet(LOCAL / "expr_e1_predictions.parquet")
    models = ["M4", "X1", "X2", "X3", "X3k", "X0", "C1"]
    y = pred.mmr.to_numpy()
    sub = (y & (pred.MLH1_pct.to_numpy() >= MLH1_LOW)) | ~y          # EX5 rows
    n_sub = int((y & sub).sum())
    strata = [np.asarray(v) for v in pred.groupby(["project", "mmr"]).indices.values()]
    rng = np.random.default_rng(m51.SEED)
    S = {m: pred[m].to_numpy() for m in models}
    b_auc = {m: np.empty(m51.BOOT) for m in models}
    b_ap = {m: np.empty(m51.BOOT) for m in models}
    b_sub = {m: np.empty(m51.BOOT) for m in ["M4", "X1"]}
    for k in range(m51.BOOT):
        i = np.concatenate([rng.choice(g, len(g)) for g in strata])
        j = i[sub[i]]
        for m in models:
            b_auc[m][k] = m54.auc(S[m][i], y[i])
            b_ap[m][k] = m54.ap(S[m][i], y[i])
        for m in b_sub:
            b_sub[m][k] = m54.auc(S[m][j], y[j])
    o_auc = {m: m54.auc(S[m], y) for m in models}
    o_ap = {m: m54.ap(S[m], y) for m in models}
    o_sub = {m: m54.auc(S[m][sub], y[sub]) for m in b_sub}
    rows = []

    def add(pid, a, b, what, o, bs, n):
        d = bs[a] - bs[b]
        rows.append({"id": pid, "what": what, "effect": o[a] - o[b], "predicted_sign": 1, "p": m51.boot_p(d),
                     "n": n, "ci_low": float(np.quantile(d, 0.025)), "ci_high": float(np.quantile(d, 0.975))})

    pos = int(y.sum())
    for pid, b in [("EX1", "X1"), ("EX2", "X2"), ("EX3", "X3"), ("EX4", "X0")]:
        add(pid, "M4", b, f"AUROC(M4) - AUROC({b}), E1 held-out", o_auc, b_auc,
            f"M4 {o_auc['M4']:.3f}, {b} {o_auc[b]:.3f}; {pos} positives")
    if n_sub >= MIN_EX5:
        add("EX5", "M4", "X1", f"AUROC(M4) - AUROC(X1), MMR-deficient with MLH1 percentile >= {MLH1_LOW}", o_sub,
            b_sub, f"M4 {o_sub['M4']:.3f}, X1 {o_sub['X1']:.3f}; {n_sub} positives")
    else:
        rows.append({"id": "EX5", "what": "subgroup too small", "effect": np.nan, "predicted_sign": 1, "p": np.nan,
                     "n": f"{n_sub} positives"})
    add("EX6", "C1", "X1", "AUROC(C1) - AUROC(X1), E1 held-out", o_auc, b_auc,
        f"C1 {o_auc['C1']:.3f}, X1 {o_auc['X1']:.3f}")
    pr = m51.predictions(rows, "expression")
    # metrics
    met = []
    for m in models:
        met.append({"model": m, "name": NAMES[m], "n": len(y), "positives": pos, "auroc": o_auc[m],
                    "auroc_lo": np.quantile(b_auc[m], 0.025), "auroc_hi": np.quantile(b_auc[m], 0.975),
                    "auprc": o_ap[m], "auprc_lo": np.quantile(b_ap[m], 0.025), "auprc_hi": np.quantile(b_ap[m], 0.975),
                    "precision_at_recall_0.5": m54.precision_at_recall(S[m], y)})
    met = pd.DataFrame(met)
    tun = pd.read_csv(OUT / "tuning_e1.csv")
    ch = tun[tun.chosen.astype(bool)].set_index("model")
    met["chosen"] = met.model.map(lambda m: f"K={ch.loc[m, 'K']}, C={ch.loc[m, 'C']:g}" if m in ch.index else "")
    met["inner_auroc"] = met.model.map(lambda m: ch.loc[m, "inner_auroc"] if m in ch.index else np.nan)
    met.to_csv(OUT / "metrics.csv", index=False, float_format="%.4g")
    log(f"metrics:\n{met.round(3).to_string(index=False)}")
    # exploratory
    ex = []
    low = y & (pred.MLH1_pct.to_numpy() < MLH1_LOW)
    ex.append({"what": "held-out MMR-deficient tumours with MLH1 percentile < 0.10", "value": int(low.sum()),
               "n": pos})
    for proj, g in pred[pred.mmr].groupby("project"):
        ex.append({"what": "held-out MMR-deficient tumours by cancer: MLH1 percentile < 0.10", "project": proj,
                   "value": int((g.MLH1_pct < MLH1_LOW).sum()), "n": len(g)})
    for m in models:
        for lab, msk in [("MLH1 low (< 0.10)", low | ~y), ("MLH1 not low (>= 0.10)", sub)]:
            ex.append({"what": f"AUROC, MMR-deficient with {lab} vs all proficient", "model": m,
                       "value": m54.auc(S[m][msk], y[msk]), "n": int((y & msk).sum())})
    for proj, g in pred.groupby("project"):
        if g.mmr.sum() >= 3:
            for m in models:
                ex.append({"what": "per-cancer AUROC, E1", "project": proj, "model": m, "n": int(g.mmr.sum()),
                           "value": m54.auc(g[m].to_numpy(), g.mmr.to_numpy())})
    calls = pd.read_parquet(LOCAL / "pancancer_msicalls.parquet")[["file_id", "call_mantis"]]
    mt = pred.merge(calls, on="file_id", how="left")
    mt = mt[mt.call_mantis.isin(["MSI-H", "MSS"])]
    ym = mt.call_mantis.eq("MSI-H").to_numpy()
    for m in models:
        ex.append({"what": "E1 with MANTIS MSI-H vs MSS as the label: AUROC", "model": m, "n": int(ym.sum()),
                   "value": m54.auc(mt[m].to_numpy(), ym)})
    # rank agreement among held-out MMR-deficient tumours
    pw = pred[pred.mmr]
    for m in ["X1", "X3", "X0", "C1"]:
        ex.append({"what": "Spearman rho with M4 among held-out MMR-deficient tumours (within-cancer percentile)",
                   "model": m, "n": len(pw), "value": pd.Series(_pct(pred, m)[pred.mmr.to_numpy()]).corr(
                       pd.Series(_pct(pred, "M4")[pred.mmr.to_numpy()]), method="spearman")})
    co = pd.read_parquet(LOCAL / "expr_coefficients.parquet")
    for m in ["C1"]:
        c = co[co.model.eq(m)]
        for lab, g in c.groupby("gene"):
            ex.append({"what": "C1: share of the summed |coefficient|", "model": m, "project": lab,
                       "value": g.coefficient.abs().sum() / c.coefficient.abs().sum(), "n": len(g)})
    e2p = LOCAL / "expr_e2_predictions.parquet"
    if e2p.exists():
        d2 = pd.read_parquet(e2p)
        y2 = d2.mmr.to_numpy()
        h = ~d2.project.isin(TRAIN).to_numpy()
        for m in ["M4", "X1"]:
            ex.append({"what": "E2 pooled AUROC (X1 with E1's parameters)", "model": m, "n": int(y2.sum()),
                       "value": m54.auc(d2[m].to_numpy(), y2)})
            ex.append({"what": "E2 AUROC, the 30 held-out cancers", "model": m, "n": int(y2[h].sum()),
                       "value": m54.auc(d2[m].to_numpy()[h], y2[h])})
        b = m55.data().merge(d2[["file_id", "M4", "X1"]], on="file_id")
        b = b[b.mantis.notna() & b.msisensor.notna()]
        refs = m55.references(b)
        for judged, lab in refs.items():
            k = lab >= 0
            for m in ["M4", "X1"]:
                ex.append({"what": f"E2 within-cancer AUROC against the consensus excluding {judged}", "model": m,
                           "n": int((lab == 1).sum()), "value": m55.auc_within(b[m].to_numpy()[k], lab[k] == 1,
                                                                             b.project.to_numpy()[k])})
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "exploratory.csv", index=False, float_format="%.4g")
    log(f"exploratory:\n{ex.round(3).to_string(index=False)}")
    top = co[co.model.eq("X1")].assign(a=lambda d: d.coefficient.abs()).nlargest(30, "a").drop(columns="a")
    top.to_csv(OUT / "coefficients_top_x1.csv", index=False, float_format="%.4g")
    co[co.model.eq("X3")].to_csv(OUT / "coefficients_x3.csv", index=False, float_format="%.4g")
    log(f"X1 top coefficients:\n{top.round(3).to_string(index=False)}")
    # ROC curves (aggregate: thresholds dropped)
    from sklearn.metrics import roc_curve
    rc = []
    for m in models:
        fpr, tpr, _ = roc_curve(y, S[m])
        keep = np.unique(np.concatenate([np.searchsorted(fpr, np.linspace(0, 1, 201)).clip(0, len(fpr) - 1)]))
        rc.append(pd.DataFrame({"model": m, "fpr": fpr[keep], "tpr": tpr[keep]}))
    pd.concat(rc).to_csv(OUT / "roc_curves.csv", index=False, float_format="%.4g")
    log(f"predictions:\n{pr[['id', 'effect', 'ci_low', 'ci_high', 'p', 'q', 'verdict', 'n']].to_string(index=False)}")


def _pct(pred: pd.DataFrame, m: str) -> np.ndarray:
    return pred.groupby("project")[m].rank(pct=True).to_numpy()


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["fit", "e2", "report"])
    a = ap.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if a.only in (None, "fit"):
            fit()
        if a.only in (None, "e2"):
            e2()
        if a.only in (None, "report"):
            report()


if __name__ == "__main__":
    main()
