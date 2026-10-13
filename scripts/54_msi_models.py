"""Step 54 (D-71): fitted models on the MSI sequence code.

Does fitting weights on top of the code beat the fixed score G in cancers never used, and does
learning a weight per exon (as the subtype models do) transfer worse than learning a weight per
sequence class? config/followup.yaml: msi_code_models.

  features  label-blind per-tumour inputs for all 33 cancers: score G, the two arms, ten sequence
            bins (robust z within the cancer) and the per-exon matrix (background-adjusted
            percentiles of the 30,632 exons in the universe of >= 25 cancers; A-17)
  fit       M0-M6 on E1 (train COAD, STAD, UCEC; test the other 30) and M0-M5 on E2 (leave one
            cancer out); LR1-LR7 and the exploratory tables
  posthoc   A-18 (exploratory, chosen after the E1 AUROCs were read): what M4 and M5 weight, and
            three ablations (M5 without repeat exons, M4 on the arms only, M4 outside the arms)

Outputs: data.nosync/ctsm_v1/msi/models_* (per-tumour inputs and predictions, kept local)
         results/subtype/msi_cis/code/models/ (aggregates)
Usage: python scripts/54_msi_models.py --only features|fit|posthoc
"""
from __future__ import annotations

import argparse
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


m42 = _load("msi42", "42_msi_repeats.py")
m51 = _load("msi51", "51_msi_code.py")
m52 = _load("msi52", "52_msi_reversal.py")
from ctsm.matrix_io import load_view  # noqa: E402

DATA, LOCAL, PAN = m51.DATA, m51.LOCAL, m51.PAN
OUT = m51.OUT / "models"
TRAIN = m51.TRAIN                                  # COAD, STAD, UCEC
POOL_MIN_CANCERS = 25
SHORT_BINS = [(10, 11), (12, 13), (14, 15), (16, 19)]
LONG_BP_BINS = [(20, 29), (30, 39), (40, 49), (50, 69), (70, 100)]
BIN_NAMES = [f"run {a}-{b}" for a, b in SHORT_BINS] + [f"run >=20, bp {a}-{b}" for a, b in LONG_BP_BINS] + \
            ["run >=20, no bp"]
GRIDS = {
    "M1": [{"C": c} for c in [0.01, 0.1, 1, 10]],
    "M2": [{"C": c} for c in [0.01, 0.1, 1, 10]],
    "M3": [{"alpha": a, "hidden": 16} for a in [0.001, 0.01, 0.1]],
    "M4": [{"C": c} for c in [1e-4, 1e-3, 1e-2, 1e-1]],
    "M5": [{"K": k, "C": c} for k in [500, 2000] for c in [1e-4, 1e-3, 1e-2, 1e-1]],
    "M6": [{"K": k, "alpha": a, "hidden": 64} for k in [500, 2000] for a in [0.001, 0.01, 0.1]],
}
KIND = {"M1": "lr", "M2": "lr", "M3": "mlp", "M4": "lr", "M5": "lr", "M6": "mlp"}
INPUT = {"M1": "arms", "M2": "bins", "M3": "bins", "M4": "repeat", "M5": "all", "M6": "all"}
SEEDS, BOOT, E2_JOBS = (0, 1, 2, 3, 4), 2000, 4


# ---------------------------------------------------------------- features (label-blind)

def robust_z(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    if not np.isfinite(v).any():
        return np.zeros_like(v)
    med = np.nanmedian(v)
    s = 1.4826 * np.nanmedian(np.abs(v - med))
    if s == 0:
        q1, q3 = np.nanpercentile(v, [25, 75])
        s = (q3 - q1) / 1.349
    z = (v - med) / s if s > 0 else np.zeros_like(v)
    return np.nan_to_num(z, nan=0.0)


def bin_masks(at: np.ndarray, db: np.ndarray) -> list[np.ndarray]:
    out = [(at >= a) & (at <= b) for a, b in SHORT_BINS]
    long = at >= 20
    out += [long & (db >= a) & (db <= b) for a, b in LONG_BP_BINS]
    out.append(long & ~np.isfinite(db))
    return out


def safe_set_score(P, m, b, dec) -> np.ndarray:
    if m.sum() == 0:
        return np.full(P.shape[0], np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return m42.set_score(P, m, b, dec)


def tumours() -> pd.DataFrame:
    tab = pd.read_parquet(DATA / "tasks" / "cancer_type__PANCAN.parquet")
    t = tab[tab.role.eq("tumour")].sort_values("sample_id").drop_duplicates("case_id").copy()
    t["sample15"] = t.sample_id.str[:15]
    cnt = m51.mc3_counts(set(t.sample15))
    t = t[t.sample15.isin(cnt.index)].join(cnt, on="sample15")
    t["mmr"] = t.n_hp >= m51.MMR_CUT
    t["hyper"] = (t.n_hp < m51.MMR_CUT) & (t.n_snv >= m51.HYPER_SNV)
    return t[["file_id", "case_id", "sample15", "project", "n_hp", "n_snv", "mmr", "hyper"]].reset_index(drop=True)


def features() -> None:
    cls = m52.pan_classes()
    at, db = cls.acc_T.to_numpy(float), cls.d_bp.to_numpy(float)
    t = tumours()
    log(f"tumours with MC3 calls: {len(t):,} in {t.project.nunique()} cancers")
    nf = len(cls)
    count = np.zeros(nf, int)
    per, rows = {}, []
    for proj, s in t.groupby("project"):
        X = load_view(PAN, "se", s.file_id)[0].get()
        ok = cls.struct.to_numpy() & (np.isfinite(X).mean(axis=0) >= 0.8)
        count += ok
        Xo = X[:, ok]
        P = m52.percentiles(Xo)
        dec = np.asarray(pd.qcut(np.nanmean(Xo, axis=0), 10, labels=False, duplicates="drop"), dtype=int)
        a_, d_ = at[ok], db[ok]
        down, up, b = m51.arms(a_, d_)
        G = m51.score_g(P, a_, d_, dec)
        f = s[["file_id"]].copy()
        f["G"] = G
        f["pct_G"] = rankdata(G) / len(G)
        f["arm_down"] = robust_z(safe_set_score(P, down, b, dec))
        f["arm_up"] = robust_z(safe_set_score(P, up, b, dec))
        for name, m in zip(BIN_NAMES, bin_masks(a_, d_)):
            f[name] = robust_z(safe_set_score(P, m, b, dec))
        rows.append(f)
        # per-exon background adjustment: percentile minus the tumour's mean B percentile in its decile
        bmean = np.full((P.shape[0], int(np.nanmax(dec)) + 1), np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for d in range(bmean.shape[1]):
                k = b & (dec == d)
                if k.any():
                    bmean[:, d] = np.nanmean(P[:, k], axis=1)
        adj = (P - bmean[:, dec]).astype(np.float32)
        per[proj] = (s.index.to_numpy(), np.flatnonzero(ok), np.nan_to_num(adj, nan=0.0))
        log(f"{proj}: {len(s)} tumours, {int(s.mmr.sum())} MMR-deficient, {int(ok.sum()):,} exons")
    f = pd.concat(rows, ignore_index=True)
    t = t.merge(f, on="file_id", how="left")
    pool = np.flatnonzero(count >= POOL_MIN_CANCERS)
    pos = np.full(nf, -1)
    pos[pool] = np.arange(len(pool))
    M = np.zeros((len(t), len(pool)), dtype=np.float32)
    for proj, (ridx, oidx, adj) in per.items():
        keep = pos[oidx] >= 0
        M[np.ix_(ridx, pos[oidx[keep]])] = adj[:, keep]
    np.save(LOCAL / "models_exon_matrix.npy", M)
    cols = cls.iloc[pool][["feature_id", "symbol", "acc_T", "d_bp", "cls"]].reset_index(drop=True)
    cols.to_parquet(LOCAL / "models_exon_columns.parquet", index=False)
    t.to_parquet(LOCAL / "models_tumours.parquet", index=False)
    log(f"features: {len(t):,} tumours, pool {len(pool):,} exons ({int((cols.acc_T >= 10).sum()):,} with a run >= 10)")


# ---------------------------------------------------------------- models

def auc(score: np.ndarray, y: np.ndarray) -> float:
    """AUROC for positives scoring higher (ties half), by ranks."""
    r = rankdata(score)
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def ap(score: np.ndarray, y: np.ndarray) -> float:
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(y, score))


def precision_at_recall(score: np.ndarray, y: np.ndarray, r: float = 0.5) -> float:
    """Precision at the rank where the ceil(r * positives)-th positive is reached (score descending,
    ties broken by order, conservatively: tied negatives first)."""
    order = np.lexsort((y.astype(int), -score))      # within ties, negatives before positives
    k = int(np.ceil(r * y.sum()))
    hits = np.cumsum(y[order])
    rank = int(np.searchsorted(hits, k) + 1)
    return k / rank


def anova_order(X, rows: np.ndarray, y: np.ndarray, chunk: int = 4000) -> np.ndarray:
    """Columns by ANOVA F on the given rows (descending), computed in column chunks."""
    from sklearn.feature_selection import f_classif
    F = np.empty(X.shape[1])
    for c0 in range(0, X.shape[1], chunk):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f, _ = f_classif(np.asarray(X[rows, c0:c0 + chunk], dtype=np.float64), y)
        F[c0:c0 + chunk] = np.nan_to_num(f, nan=-np.inf)
    return np.argsort(-F, kind="stable")


def fit_predict(model: str, params: dict, X, tr: np.ndarray, te: np.ndarray, y: np.ndarray,
                order: np.ndarray | None = None, seeds=(0,)) -> np.ndarray:
    """Fit on rows `tr` of X (y over all rows) and score rows `te` (higher = MMR-deficient)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler
    if "K" in params:                                   # read only the selected columns
        cols = np.sort(order[:params["K"]])
        A = np.asarray(X[np.ix_(tr, cols)], dtype=np.float64)
        B = np.asarray(X[np.ix_(te, cols)], dtype=np.float64)
    else:
        A, B = np.asarray(X[tr], dtype=np.float64), np.asarray(X[te], dtype=np.float64)
    sc = StandardScaler().fit(A)
    A, B = sc.transform(A), sc.transform(B)
    ytr = y[tr]
    if KIND[model] == "lr":
        m = LogisticRegression(C=params["C"], class_weight="balanced", max_iter=5000)
        return m.fit(A, ytr).decision_function(B)
    out = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        pos, neg = np.flatnonzero(ytr), np.flatnonzero(~ytr)
        idx = np.concatenate([neg, rng.choice(pos, len(neg), replace=True)])
        m = MLPClassifier(hidden_layer_sizes=(params["hidden"],), activation="relu", solver="adam",
                          alpha=params["alpha"], max_iter=500, random_state=seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(A[idx], ytr[idx])
        out.append(m.predict_proba(B)[:, 1])
    return np.mean(out, axis=0)


def tune(model: str, X, rows: np.ndarray, y: np.ndarray, groups: np.ndarray, folds) -> tuple[dict, float, list]:
    """Grid search by pooled out-of-fold AUROC over `folds` (index pairs into `rows`)."""
    orders = [anova_order(X, rows[a], y[rows[a]]) if INPUT[model] == "all" else None for a, _ in folds]
    table, best = [], None
    for params in GRIDS[model]:
        oof = np.full(len(rows), np.nan)
        for (a, b), od in zip(folds, orders):
            oof[b] = fit_predict(model, params, X, rows[a], rows[b], y, order=od)
        au = auc(oof, y[rows])
        table.append({**params, "inner_auroc": au})
        if best is None or au > best[1]:
            best = (params, au)
    return best[0], best[1], table


def inputs(t: pd.DataFrame, M: np.ndarray, cols: pd.DataFrame) -> dict:
    return {"arms": t[["arm_down", "arm_up"]].to_numpy(float),
            "bins": t[BIN_NAMES].to_numpy(float),
            "repeat": M[:, (cols.acc_T >= 10).to_numpy()],
            "all": M}


def run_split(model: str, X, y, groups, tr_rows, te_rows, inner) -> dict:
    """Tune on tr_rows with the inner splitter, refit, score te_rows."""
    folds = list(inner.split(tr_rows, y[tr_rows], groups[tr_rows]))
    params, inner_au, table = tune(model, X, tr_rows, y, groups, folds)
    order = anova_order(X, tr_rows, y[tr_rows]) if INPUT[model] == "all" else None
    seeds = SEEDS if KIND[model] == "mlp" else (0,)
    score = fit_predict(model, params, X, tr_rows, te_rows, y, order=order, seeds=seeds)
    return {"params": params, "inner_auroc": inner_au, "table": table, "score": score}


def _e2_fold(model, X, y, groups, c):
    from sklearn.model_selection import GroupKFold
    te = np.flatnonzero(groups == c)
    tr = np.flatnonzero(groups != c)
    r = run_split(model, X, y, groups, tr, te, GroupKFold(n_splits=4))
    return c, te, r


def coefficients(model: str, X, y, rows, params, names) -> pd.DataFrame:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    A = StandardScaler().fit_transform(np.asarray(X[rows], dtype=np.float64))
    m = LogisticRegression(C=params["C"], class_weight="balanced", max_iter=5000).fit(A, y[rows])
    return pd.DataFrame({"model": model, "feature": names, "coefficient": m.coef_[0]})


def fit() -> None:
    from joblib import Parallel, delayed
    from sklearn.model_selection import LeaveOneGroupOut
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")
    cols = pd.read_parquet(LOCAL / "models_exon_columns.parquet")
    M = np.load(LOCAL / "models_exon_matrix.npy", mmap_mode="r")
    y, groups = t.mmr.to_numpy(), t.project.to_numpy()
    X = inputs(t, M, cols)
    held = ~t.project.isin(TRAIN).to_numpy()
    tr1, te1 = np.flatnonzero(~held), np.flatnonzero(held)
    log(f"E1: train {len(tr1):,} tumours ({int(y[tr1].sum())} MMR-deficient), test {len(te1):,} ({int(y[te1].sum())})")
    OUT.mkdir(parents=True, exist_ok=True)
    pred1 = pd.DataFrame({"file_id": t.file_id[te1].to_numpy(), "project": groups[te1], "mmr": y[te1],
                          "M0": 1 - t.pct_G.to_numpy()[te1]})
    tuning, e1 = [], {}
    for model in ["M1", "M2", "M3", "M4", "M5", "M6"]:
        r = run_split(model, X[INPUT[model]], y, groups, tr1, te1, LeaveOneGroupOut())
        pred1[model] = r["score"]
        e1[model] = r
        tuning += [{"split": "E1", "model": model, **row, "chosen": row_params == r["params"]}
                   for row in r["table"] for row_params in [{k: row[k] for k in r["params"]}]]
        log(f"E1 {model}: chosen {r['params']}, inner AUROC {r['inner_auroc']:.3f}, held-out AUROC "
            f"{auc(r['score'], y[te1]):.3f}, AUPRC {ap(r['score'], y[te1]):.3f}")
    pred1.to_parquet(LOCAL / "models_e1_predictions.parquet", index=False)
    # E2: leave one cancer out (M6 left out for compute, as pre-specified)
    pred2 = pd.DataFrame({"file_id": t.file_id, "project": groups, "mmr": y, "M0": 1 - t.pct_G.to_numpy()})
    e2_params = []
    for model in ["M1", "M2", "M3", "M4", "M5"]:
        res = Parallel(n_jobs=E2_JOBS)(delayed(_e2_fold)(model, X[INPUT[model]], y, groups, c)
                                       for c in sorted(set(groups)))
        sc_ = np.full(len(t), np.nan)
        for c, te, r in res:
            sc_[te] = r["score"]
            e2_params.append({"split": "E2", "model": model, "left_out": c, **r["params"], "inner_auroc": r["inner_auroc"]})
        pred2[model] = sc_
        log(f"E2 {model}: pooled AUROC {auc(sc_, y):.3f}, AUPRC {ap(sc_, y):.3f}")
    pred2.to_parquet(LOCAL / "models_e2_predictions.parquet", index=False)
    pd.DataFrame(tuning).to_csv(OUT / "tuning_e1.csv", index=False, float_format="%.4g")
    pd.DataFrame(e2_params).to_csv(OUT / "tuning_e2.csv", index=False, float_format="%.4g")
    report(t, pred1, pred2, e1, X, y, tr1)


def report(t, pred1, pred2, e1, X, y, tr1) -> None:
    models = ["M0", "M1", "M2", "M3", "M4", "M5", "M6"]
    y1 = pred1.mmr.to_numpy()
    strata = [np.asarray(v) for v in pred1.groupby(["project", "mmr"]).indices.values()]
    rng = np.random.default_rng(m51.SEED)
    boot_auc = {m: np.empty(BOOT) for m in models}
    boot_ap = {m: np.empty(BOOT) for m in models}
    for k in range(BOOT):
        i = np.concatenate([rng.choice(g, len(g)) for g in strata])
        for m in models:
            s_ = pred1[m].to_numpy()[i]
            boot_auc[m][k] = auc(s_, y1[i])
            boot_ap[m][k] = ap(s_, y1[i])
    obs_auc = {m: auc(pred1[m].to_numpy(), y1) for m in models}
    obs_ap = {m: ap(pred1[m].to_numpy(), y1) for m in models}
    pred = []

    def add(pid, a, b, sign, metric="AUROC"):
        o, bs = (obs_auc, boot_auc) if metric == "AUROC" else (obs_ap, boot_ap)
        d = bs[a] - bs[b]
        pred.append({"id": pid, "what": f"{metric}({a}) - {metric}({b}), E1 held-out", "effect": o[a] - o[b],
                     "predicted_sign": sign, "p": m51.boot_p(d), "n": f"{a} {o[a]:.3f}, {b} {o[b]:.3f}",
                     "ci_low": float(np.quantile(d, 0.025)), "ci_high": float(np.quantile(d, 0.975))})

    add("LR1", "M1", "M0", 1)
    add("LR2", "M2", "M1", 1)
    add("LR3", "M3", "M2", -1)
    add("LR4", "M4", "M2", -1)
    add("LR5", "M5", "M2", -1)
    add("LR6", "M6", "M2", -1)
    add("LR7", "M2", "M0", 1, metric="AUPRC")
    pr = m51.predictions(pred, "models")
    # metrics table
    rows = []
    prev1 = y1.mean()
    for m in models:
        rows.append({"split": "E1", "model": m, "n": len(y1), "positives": int(y1.sum()), "prevalence": prev1,
                     "auroc": obs_auc[m], "auroc_lo": np.quantile(boot_auc[m], 0.025), "auroc_hi": np.quantile(boot_auc[m], 0.975),
                     "auprc": obs_ap[m], "auprc_lo": np.quantile(boot_ap[m], 0.025), "auprc_hi": np.quantile(boot_ap[m], 0.975),
                     "precision_at_recall_0.5": precision_at_recall(pred1[m].to_numpy(), y1),
                     "inner_auroc": e1[m]["inner_auroc"] if m in e1 else np.nan,
                     "chosen": str(e1[m]["params"]) if m in e1 else "no fit"})
    y2 = pred2.mmr.to_numpy()
    for m in models[:-1]:
        s_ = pred2[m].to_numpy()
        rows.append({"split": "E2", "model": m, "n": len(y2), "positives": int(y2.sum()), "prevalence": y2.mean(),
                     "auroc": auc(s_, y2), "auprc": ap(s_, y2), "precision_at_recall_0.5": precision_at_recall(s_, y2)})
        h = ~pred2.project.isin(TRAIN).to_numpy()
        rows.append({"split": "E2, the 30 held-out cancers only", "model": m, "n": int(h.sum()),
                     "positives": int(y2[h].sum()), "prevalence": y2[h].mean(), "auroc": auc(s_[h], y2[h]),
                     "auprc": ap(s_[h], y2[h]), "precision_at_recall_0.5": precision_at_recall(s_[h], y2[h])})
    met = pd.DataFrame(rows)
    met.to_csv(OUT / "metrics.csv", index=False, float_format="%.4g")
    log(f"metrics:\n{met.round(3).to_string(index=False)}")
    # coefficients of the interpretable models (E1 final fit, standardised inputs)
    co = pd.concat([coefficients("M1", X["arms"], y, tr1, e1["M1"]["params"], ["falling arm (R_down)", "rising arm (R_up)"]),
                    coefficients("M2", X["bins"], y, tr1, e1["M2"]["params"], BIN_NAMES)], ignore_index=True)
    co.to_csv(OUT / "coefficients.csv", index=False, float_format="%.4g")
    log(f"coefficients:\n{co.round(3).to_string(index=False)}")
    # per cancer (E1, >= 3 positives), MANTIS labels, curves
    ex = []
    for proj, g in pred1.groupby("project"):
        if g.mmr.sum() >= 3:
            for m in models:
                ex.append({"what": "per-cancer AUROC, E1", "project": proj, "model": m, "positives": int(g.mmr.sum()),
                           "n": len(g), "value": auc(g[m].to_numpy(), g.mmr.to_numpy())})
    calls = pd.read_parquet(LOCAL / "pancancer_msicalls.parquet")[["file_id", "call_mantis"]]
    mt = pred1.merge(calls, on="file_id", how="left")
    mt = mt[mt.call_mantis.isin(["MSI-H", "MSS"])]
    ym = mt.call_mantis.eq("MSI-H").to_numpy()
    for m in models:
        ex.append({"what": "E1 with MANTIS MSI-H vs MSS as the test label: AUROC", "model": m,
                   "positives": int(ym.sum()), "n": len(mt), "value": auc(mt[m].to_numpy(), ym)})
        ex.append({"what": "E1 with MANTIS MSI-H vs MSS as the test label: AUPRC", "model": m,
                   "positives": int(ym.sum()), "n": len(mt), "value": ap(mt[m].to_numpy(), ym)})
    pd.DataFrame(ex).to_csv(OUT / "exploratory.csv", index=False, float_format="%.4g")
    curves = []
    for m in models:
        s_ = pred1[m].to_numpy()
        order = np.lexsort((y1.astype(int), -s_))
        hits = np.cumsum(y1[order])
        for r in np.linspace(0.02, 1.0, 50):
            k = int(np.ceil(r * y1.sum()))
            rank = int(np.searchsorted(hits, k) + 1)
            curves.append({"model": m, "recall": k / y1.sum(), "precision": k / rank, "rank": rank})
    pd.DataFrame(curves).drop_duplicates().to_csv(OUT / "pr_curves.csv", index=False, float_format="%.4g")
    log(f"predictions:\n{pr[['id', 'effect', 'p', 'q', 'verdict', 'n']].to_string(index=False)}")


# ---------------------------------------------------------------- A-18 (exploratory)

def _chosen(tun: pd.DataFrame, model: str) -> dict:
    r = tun[tun.model.eq(model) & tun.chosen.astype(bool)].iloc[0]
    return {k: (int(r[k]) if k in ("K", "hidden") else float(r[k])) for k in GRIDS[model][0]}


def posthoc() -> None:
    """A-18, chosen after the E1 AUROCs were read: does the per-exon gain come from the repeat sequence?"""
    from scipy.stats import spearmanr
    from sklearn.model_selection import LeaveOneGroupOut
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")
    cols = pd.read_parquet(LOCAL / "models_exon_columns.parquet")
    M = np.load(LOCAL / "models_exon_matrix.npy", mmap_mode="r")
    pred1 = pd.read_parquet(LOCAL / "models_e1_predictions.parquet")
    tun = pd.read_csv(OUT / "tuning_e1.csv")
    y, groups = t.mmr.to_numpy(), t.project.to_numpy()
    held = ~t.project.isin(TRAIN).to_numpy()
    tr1, te1 = np.flatnonzero(~held), np.flatnonzero(held)
    assert (pred1.file_id.to_numpy() == t.file_id.to_numpy()[te1]).all()
    at, db = cols.acc_T.to_numpy(float), cols.d_bp.to_numpy(float)
    down, up, _ = m51.arms(at, db)
    rep = at >= 10
    other_long = rep & ~down & ~up
    cls_ = np.select([down, up, other_long], ["run 10-19 (R_down)", "run >=20, bp >=36 nt (R_up)",
                                              "run >=20, bp <36 nt or none"], "no run >= 10")
    rows, ex = [], []

    def add(what, value, group="", n=None, lo=np.nan, hi=np.nan):
        rows.append({"what": what, "group": group, "n": n, "value": value, "ci_low": lo, "ci_high": hi})

    # (i) M4 coefficients by sequence class
    p4 = _chosen(tun, "M4")
    c4 = coefficients("M4", M[:, rep], y, tr1, p4, cols.feature_id[rep].tolist())
    c4 = c4.assign(symbol=cols.symbol[rep].to_numpy(), acc_T=at[rep], d_bp=db[rep], cls=cls_[rep])
    for k, g in c4.groupby("cls"):
        add("M4 standardised coefficient: median", g.coefficient.median(), k, len(g))
        add("M4 standardised coefficient: share > 0", (g.coefficient > 0).mean(), k, len(g))
    rho = spearmanr(c4.coefficient, c4.acc_T)[0]
    add("M4 coefficient against run length: Spearman", rho, "all repeat exons", len(c4))
    lg = c4[(c4.acc_T >= 20) & np.isfinite(c4.d_bp)]
    add("M4 coefficient against branch-point distance: Spearman", spearmanr(lg.coefficient, lg.d_bp)[0],
        "runs >= 20 with a branch point", len(lg))
    w = np.sort(np.abs(c4.coefficient.to_numpy()))[::-1]
    add("M4 share of the summed |coefficient| in the 20 largest", w[:20].sum() / w.sum(), "", len(w))
    add("M4 share of the summed |coefficient| in the 100 largest", w[:100].sum() / w.sum(), "", len(w))
    top4 = c4.reindex(c4.coefficient.abs().sort_values(ascending=False).index).head(20)
    for i, r in enumerate(top4.itertuples(), 1):
        ex.append({"list": "M4: 20 largest |standardised coefficient|", "rank": i, "feature_id": r.feature, "symbol": r.symbol,
                   "acc_T": r.acc_T, "d_bp": r.d_bp, "class": r.cls, "value": r.coefficient})
    # (ii) M5's top-500 ANOVA exons on the E1 training set
    order = anova_order(M, tr1, y[tr1])
    top = order[:500]
    A = np.asarray(M[np.ix_(tr1, top)], dtype=np.float64)
    diff = A[y[tr1]].mean(axis=0) - A[~y[tr1]].mean(axis=0)
    add("M5 top 500: share with a run >= 10", rep[top].mean(), "", 500)
    add("pool: share with a run >= 10", rep.mean(), "", len(rep))
    for k in ["run 10-19 (R_down)", "run >=20, bp >=36 nt (R_up)", "run >=20, bp <36 nt or none", "no run >= 10"]:
        m_ = cls_[top] == k
        add("M5 top 500: share in class", m_.mean(), k, int(m_.sum()))
        add("pool: share in class", (cls_ == k).mean(), k, int((cls_ == k).sum()))
        if m_.any():
            add("M5 top 500: share more included in MMR-deficient training tumours", (diff[m_] > 0).mean(), k, int(m_.sum()))
    mmr_genes = {"MLH1", "MSH2", "MSH6", "PMS2", "EPM2AIP1"}
    sym = cols.symbol.to_numpy()
    rank_of = np.empty(len(order), int)
    rank_of[order] = np.arange(1, len(order) + 1)
    for g_ in sorted(mmr_genes):
        k = np.flatnonzero(sym == g_)
        add("best ANOVA rank of the gene's exons (pool)", int(rank_of[k].min()) if len(k) else np.nan, g_, len(k))
    for i, j in enumerate(top[:20], 1):
        ex.append({"list": "M5: 20 largest ANOVA F (E1 training set)", "rank": i, "feature_id": cols.feature_id[j],
                   "symbol": sym[j], "acc_T": at[j], "d_bp": db[j], "class": cls_[j], "value": diff[i - 1]})

    # (ii b) M5 refitted at its chosen K and C: standardised coefficients by class
    p5 = _chosen(tun, "M5")
    sel = np.sort(order[:p5["K"]])
    c5 = coefficients("M5", M[:, sel], y, tr1, {"C": p5["C"]}, cols.feature_id[sel].tolist())
    c5 = c5.assign(symbol=sym[sel], acc_T=at[sel], d_bp=db[sel], cls=cls_[sel])
    for k, g in c5.groupby("cls"):
        add("M5 standardised coefficient: median", g.coefficient.median(), k, len(g))
        add("M5 standardised coefficient: share > 0", (g.coefficient > 0).mean(), k, len(g))
    w5 = np.sort(np.abs(c5.coefficient.to_numpy()))[::-1]
    add("M5 share of the summed |coefficient| in the 20 largest", w5[:20].sum() / w5.sum(), "", len(w5))
    top5 = c5.reindex(c5.coefficient.abs().sort_values(ascending=False).index).head(20)
    for i, r in enumerate(top5.itertuples(), 1):
        ex.append({"list": "M5: 20 largest |standardised coefficient|", "rank": i, "feature_id": r.feature,
                   "symbol": r.symbol, "acc_T": r.acc_T, "d_bp": r.d_bp, "class": r.cls, "value": r.coefficient})
    # (iii) ablations, E1, tuned as the model they modify
    variants = {"M5n": ("M5", ~rep, "all"), "M4a": ("M4", down | up, "repeat"), "M4x": ("M4", other_long, "repeat")}
    scores = {m: pred1[m].to_numpy() for m in ["M0", "M2", "M4", "M5"]}
    for name, (base, mask, inp) in variants.items():
        GRIDS[name], KIND[name], INPUT[name] = GRIDS[base], KIND[base], inp
        r = run_split(name, M[:, mask], y, groups, tr1, te1, LeaveOneGroupOut())
        scores[name] = r["score"]
        add("ablation: exons", int(mask.sum()), name)
        add("ablation: chosen", str(r["params"]), name)
        add("ablation: inner AUROC", r["inner_auroc"], name)
        log(f"A-18 {name}: {int(mask.sum()):,} exons, chosen {r['params']}, inner {r['inner_auroc']:.3f}, "
            f"held-out AUROC {auc(r['score'], y[te1]):.3f}")
    y1 = pred1.mmr.to_numpy()
    strata = [np.asarray(v) for v in pred1.groupby(["project", "mmr"]).indices.values()]
    rng = np.random.default_rng(m51.SEED)
    names = list(scores)
    ba = {m: np.empty(BOOT) for m in names}
    bp = {m: np.empty(BOOT) for m in names}
    for k in range(BOOT):
        i = np.concatenate([rng.choice(g, len(g)) for g in strata])
        for m in names:
            ba[m][k], bp[m][k] = auc(scores[m][i], y1[i]), ap(scores[m][i], y1[i])
    for name, (base, _, _) in variants.items():
        for metric, f_, b_ in [("AUROC", auc, ba), ("AUPRC", ap, bp)]:
            add(f"ablation: held-out {metric}", f_(scores[name], y1), name, len(y1),
                np.quantile(b_[name], 0.025), np.quantile(b_[name], 0.975))
            for ref in ["M2", base]:
                d = b_[name] - b_[ref]
                add(f"ablation: held-out {metric} minus {ref}", f_(scores[name], y1) - f_(scores[ref], y1), name,
                    len(y1), np.quantile(d, 0.025), np.quantile(d, 0.975))
    ph = pd.DataFrame(rows)
    ph.to_csv(OUT / "posthoc.csv", index=False, float_format="%.4g")
    pd.DataFrame(ex).to_csv(OUT / "posthoc_exons.csv", index=False, float_format="%.4g")
    log(f"posthoc:\n{ph.to_string(index=False)}")
    log(f"exons:\n{pd.DataFrame(ex).to_string(index=False)}")


def main() -> None:
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--only", choices=["features", "fit", "posthoc"], required=True)
    {"features": features, "fit": fit, "posthoc": posthoc}[ap_.parse_args().only]()


if __name__ == "__main__":
    main()
