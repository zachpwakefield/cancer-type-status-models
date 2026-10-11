"""Step 11 (B-04): which individual events carry the tumour signal, are they shared across
cancers, and what do they do to the encoded protein?

1. Recurrence. The top-k events (mean |SHAP|, held-out samples) of every tumour-status model;
   an event is 'recurrent' if it is in the top-k of >= min_cancers_recurrent cancers. The number
   expected by chance is computed from each event's measured cohorts (Poisson-binomial with
   p = k / events measured in that cohort). For each top event the median paired tumour - normal
   change is computed in EVERY cancer (not only where it ranked), so the cross-cancer pattern is
   not selected on.
2. Protein consequences (SE). Cassette exons in the top-k of any tumour-status model or any
   cancer-type class are compared with all measured cassette exons (PANCAN pre-filter): length
   divisible by 3, exact GENCODE v45 exon, coding, and overlap with Pfam/CDD domains, MobiDB-lite
   disordered regions and TM helices (residues of the best transcript: MANE > canonical > other).
   Fisher exact tests, BH across metrics and sets.
(Portraits of the most recurrent events are drawn by scripts/15_bio_figures.py from these tables.)

Outputs: results/biology/04_events/ (tables; figures: scripts/15_bio_figures.py)
"""
from __future__ import annotations

import argparse
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import fisher_exact, wilcoxon

from ctsm.biology import ProteinMap, bh, bio_config, feature_values, full_importance, matched_pairs, task_table
from ctsm.config import PAN, load_config

warnings.filterwarnings("ignore", category=RuntimeWarning)


def top_table(cfg, views, cancers, k) -> pd.DataFrame:
    rows = []
    for view in views:
        for c in cancers:
            imp = full_importance(cfg, "tumor_status", view, c).head(k)
            rows.append(imp.assign(view=view, cancer=c, norm_imp=imp["mean_abs_shap"] / imp["mean_abs_shap"].iloc[0]))
    return pd.concat(rows, ignore_index=True)


def paired_stats(args) -> pd.DataFrame:
    """Median paired tumour - normal change for the given events in one cancer."""
    view, cancer, fids = args
    cfg = load_config()
    tab = task_table(cfg, "tumor_status", cancer)
    ti, ni, _ = matched_pairs(tab)
    X = feature_values(cfg, view, tab["file_id"], fids)
    D = X[ti] - X[ni]
    n_obs = np.isfinite(D).sum(0)
    p = [wilcoxon(d[np.isfinite(d)]).pvalue if np.isfinite(d).sum() >= 5 and np.any(d[np.isfinite(d)] != 0)
         else np.nan for d in D.T]
    return pd.DataFrame({"view": view, "cancer": cancer, "feature_id": fids, "n_pairs": n_obs,
                         "median_paired_delta": np.nanmedian(D, axis=0), "p_paired": p,
                         "median_normal": np.nanmedian(X[ni], axis=0), "median_tumour": np.nanmedian(X[ti], axis=0)})


def chance_recurrence(cfg, views, cancers, k, m) -> dict:
    """Expected number of events in the top-k of >= m cancers if rankings were independent."""
    out = {}
    for view in views:
        pf = pd.read_parquet(cfg.matrix_dir / view / "prefilter.parquet", columns=cancers).to_numpy()
        p = k / pf.sum(0)                                   # per-cancer chance of being top-k
        P = pf * p                                          # events x cancers (0 where not measured)
        dist = np.zeros((P.shape[0], len(cancers) + 1))
        dist[:, 0] = 1
        for j in range(P.shape[1]):                          # Poisson-binomial by dynamic programming
            pj = P[:, j:j + 1]
            dist[:, 1:] = dist[:, 1:] * (1 - pj) + dist[:, :-1] * pj
            dist[:, 0] *= 1 - pj[:, 0]
        out[view] = float(dist[:, m:].sum())
    return out


def wilson(x, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = x / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def se_coords(region: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    se = region.str.split(":").str[1].str.split("-")
    return se.str[0].astype(int).to_numpy(), se.str[1].astype(int).to_numpy()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    cfg = load_config()
    ec = bio_config()["events"]
    k, m = ec["top_k"], ec["min_cancers_recurrent"]
    out = cfg.results_dir / "biology" / "04_events"
    out.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = list(scope.loc[scope.tumor_status_in_scope, "project"])
    as_views = [v for v in cfg.views if v != "gex"]

    # ---------------- 1. recurrence ----------------
    T = top_table(cfg, as_views, cancers, k)
    T.to_csv(out / "top_events_per_cancer.csv", index=False)
    rec = (T.groupby(["view", "feature_id"])
           .agg(n_cancers=("cancer", "nunique"), cancers=("cancer", lambda x: ",".join(sorted(x))),
                sum_norm_importance=("norm_imp", "sum"), best_cancer=("cancer", "first"),
                symbol=("symbol", "first"), gene_id=("gene_id", "first"), region=("region", "first"),
                label=("label", "first"), mean_direction=("direction_spearman", "mean"))
           .reset_index())
    best = T.sort_values("norm_imp", ascending=False).drop_duplicates(["view", "feature_id"])
    rec["best_cancer"] = rec.set_index(["view", "feature_id"]).index.map(
        best.set_index(["view", "feature_id"])["cancer"])
    expected = chance_recurrence(cfg, as_views, cancers, k, m)

    jobs = [(v, c, sorted(T.loc[T.view == v, "feature_id"].unique())) for v in as_views for c in cancers]
    with ProcessPoolExecutor(a.workers) as ex:
        P = pd.concat(list(ex.map(paired_stats, jobs)), ignore_index=True)
    P["q_paired"] = bh(P["p_paired"].fillna(1).to_numpy())
    P.to_csv(out / "paired_change_top_events_all_cancers.csv.gz", index=False)
    piv = P.pivot_table(index=["view", "feature_id"], columns="cancer", values="median_paired_delta")
    sgn = np.sign(piv)
    maj = np.sign(sgn.sum(axis=1))
    rec = rec.set_index(["view", "feature_id"])
    rec["n_cancers_measured"] = piv.notna().sum(axis=1)
    rec["n_cancers_same_sign_as_majority"] = (sgn.eq(maj, axis=0)).sum(axis=1)
    rec["majority_sign"] = maj.map({1.0: "up in tumour", -1.0: "down in tumour", 0.0: "mixed"})
    rec["median_delta_across_cancers"] = piv.median(axis=1)
    rec = rec.reset_index().sort_values(["n_cancers", "sum_norm_importance"], ascending=False)
    rec.to_csv(out / "event_recurrence_all_top_events.csv", index=False)
    recurrent = rec[rec.n_cancers >= m]
    recurrent.to_csv(out / "recurrent_events.csv", index=False)
    summ = pd.DataFrame({"view": as_views,
                         "n_top_events": [T[T.view == v].feature_id.nunique() for v in as_views],
                         "n_recurrent_observed": [int((recurrent.view == v).sum()) for v in as_views],
                         "n_recurrent_expected_by_chance": [expected[v] for v in as_views]})
    summ.to_csv(out / "recurrence_vs_chance.csv", index=False)
    log("recurrence:\n" + summ.to_string(index=False))

    # ---------------- 2. protein consequences of cassette exons ----------------
    pm = ProteinMap(cfg)
    feats = pd.read_parquet(cfg.matrix_dir / "se" / "features.parquet", columns=["feature_id", "gene_id", "region"])
    pf = pd.read_parquet(cfg.matrix_dir / "se" / "prefilter.parquet", columns=[PAN])[PAN].to_numpy()
    bg = feats[pf].reset_index(drop=True)
    ts_top = set(T.loc[T.view == "se", "feature_id"])
    n_cls = pd.read_parquet(cfg.data_dir / "tasks" / f"cancer_type__{PAN}.parquet")["y"].max() + 1
    ct_top = set()
    for c in range(n_cls):
        ct_top |= set(full_importance(cfg, "cancer_type", "se", PAN, class_index=c)["feature_id"].head(k))
    rng = np.random.default_rng(0)
    bg_s = bg.iloc[np.sort(rng.choice(len(bg), size=min(6000, len(bg)), replace=False))]
    ev = pd.concat([bg_s.assign(set="all measured exons (random 6,000)"),
                    feats[feats.feature_id.isin(ts_top)].assign(set="tumour vs normal, top events"),
                    feats[feats.feature_id.isin(ct_top)].assign(set="cancer type, top events")], ignore_index=True)
    s0, e0 = se_coords(ev["region"])
    ev["length"] = e0 - s0 + 1
    ev["len_mult3"] = ev["length"] % 3 == 0
    cons = pd.DataFrame([pm.exon_consequence(g, int(s), int(e)) for g, s, e in zip(ev.gene_id, s0, e0)])
    ev = pd.concat([ev, cons.drop(columns=["exon_length", "strand"], errors="ignore")], axis=1)
    # direction in tumour for the tumour-status top exons (median change in the cancers where top)
    tsd = (T[T.view == "se"][["feature_id", "cancer"]]
           .merge(P[P.view == "se"][["feature_id", "cancer", "median_paired_delta"]], on=["feature_id", "cancer"])
           .groupby("feature_id")["median_paired_delta"].median())
    ev["tumour_direction"] = ev["feature_id"].map(np.sign(tsd)).map({1.0: "included more in tumour",
                                                                     -1.0: "skipped more in tumour"})
    ev.to_csv(out / "cassette_exon_consequences.csv", index=False)
    metrics = [("length divisible by 3", "len_mult3", None),
               ("exact GENCODE exon", "annotated_exon", None),
               ("protein-coding exon", "coding", "annotated_exon"),
               ("overlaps a Pfam/CDD domain", "in_domain", "fully_coding"),
               ("overlaps a disordered region", "in_disordered", "fully_coding"),
               ("overlaps a TM helix", "in_transmembrane", "fully_coding")]
    sets = list(ev["set"].unique())
    rows = []
    for name, col, among in metrics:
        base = ev if among is None else ev[ev[among].fillna(False).astype(bool)]
        b = base[base.set == sets[0]][col].fillna(False).astype(bool)
        for st in sets:
            x = base[base.set == st][col].fillna(False).astype(bool)
            lo, hi = wilson(int(x.sum()), len(x))
            r = {"metric": name, "among": among or "all", "set": st, "n": len(x), "k": int(x.sum()),
                 "proportion": x.mean(), "ci_lo": lo, "ci_hi": hi}
            if st != sets[0]:
                odds, p = fisher_exact([[x.sum(), (~x).sum()], [b.sum(), (~b).sum()]])
                r.update(odds_ratio=odds, p_fisher=p)
            rows.append(r)
        for dname in ["included more in tumour", "skipped more in tumour"]:
            x = base[(base.set == sets[1]) & (base.tumour_direction == dname)][col].fillna(False).astype(bool)
            lo, hi = wilson(int(x.sum()), len(x))
            odds, p = fisher_exact([[x.sum(), (~x).sum()], [b.sum(), (~b).sum()]])
            rows.append({"metric": name, "among": among or "all", "set": f"tumour vs normal, {dname}", "n": len(x),
                         "k": int(x.sum()), "proportion": x.mean(), "ci_lo": lo, "ci_hi": hi,
                         "odds_ratio": odds, "p_fisher": p})
    pr = pd.DataFrame(rows)
    has_p = pr["p_fisher"].notna()
    pr.loc[has_p, "q_fisher"] = bh(pr.loc[has_p, "p_fisher"].to_numpy())
    pr.to_csv(out / "cassette_exon_protein_features.csv", index=False)
    log("protein consequences:\n" + pr[["metric", "set", "n", "proportion", "odds_ratio", "q_fisher"]].to_string())

    log("events done")


if __name__ == "__main__":
    main()
