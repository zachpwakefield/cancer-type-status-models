"""Step 13 (B-06): cell lineage in splicing and in expression.

1. Confusions. Which pairs of classes does each family mix up (out-of-fold, argmax)? GEX vs the
   AS ensemble (mean of the 8 AS views' probabilities).
2. Is the AS signal an expression proxy? (a) For each view, the per-class top events are
   correlated (Spearman, random 3,000 samples) with the expression of their host gene, compared
   with random measured events of the same view. (b) For each class, overlap of the top-100 GEX
   genes with the host genes of the top-100 events of every AS view, against the overlap
   expected by chance within genes measured by both.
3. A lineage atlas of cassette exons: the top SE events of each class (per-class SHAP), median
   PSI in every class.

Outputs: results/biology/06_lineage/ (tables; figures: scripts/15_bio_figures.py)
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import spearmanr

from ctsm.biology import feature_values, full_importance, gene_ranking, gex_by_symbol
from ctsm.config import PAN, load_config

warnings.filterwarnings("ignore", category=RuntimeWarning)
TOP_EVENTS = 20      # per class and view, for the expression-proxy test
ATLAS_PER_CLASS = 2  # SE events per class in the atlas


def main() -> None:
    cfg = load_config()
    out = cfg.results_dir / "biology" / "06_lineage"
    out.mkdir(parents=True, exist_ok=True)
    tab = pd.read_parquet(cfg.data_dir / "tasks" / f"cancer_type__{PAN}.parquet")
    classes = tab.drop_duplicates("y").sort_values("y")["label"].tolist()
    as_views = [v for v in cfg.views if v != "gex"]
    rng = np.random.default_rng(0)

    # ---------------- 1. confusions ----------------
    P = {v: pd.read_parquet(cfg.run_dir("cancer_type", v, PAN) / "oof_predictions.parquet") for v in cfg.views}
    pc = [f"p_{c}" for c in classes]
    y = P["gex"]["y"].to_numpy()
    preds = {"GEX": P["gex"][pc].to_numpy().argmax(1),
             "AS ensemble": np.mean([P[v][pc].to_numpy() for v in as_views], axis=0).argmax(1)}
    rows = []
    for i in range(len(classes)):
        for j in range(i + 1, len(classes)):
            r = {"pair": f"{classes[i]} / {classes[j]}"}
            for nm, pr in preds.items():
                r[nm] = int(((y == i) & (pr == j)).sum() + ((y == j) & (pr == i)).sum())
            rows.append(r)
    conf = pd.DataFrame(rows)
    conf["max"] = conf[["GEX", "AS ensemble"]].max(1)
    conf = conf.sort_values("max", ascending=False)
    conf.drop(columns="max").to_csv(out / "class_pair_confusions.csv", index=False)

    # ---------------- 2a. event value vs host-gene expression ----------------
    samp = np.sort(rng.choice(len(tab), size=min(3000, len(tab)), replace=False))
    fids_s = tab["file_id"].iloc[samp]
    rho_rows = []
    for view in as_views:
        feats = pd.read_parquet(cfg.matrix_dir / view / "features.parquet", columns=["feature_id", "symbol"])
        pf = pd.read_parquet(cfg.matrix_dir / view / "prefilter.parquet", columns=[PAN])[PAN].to_numpy()
        top = set()
        for k in range(len(classes)):
            top |= set(full_importance(cfg, "cancer_type", view, PAN, class_index=k)["feature_id"].head(TOP_EVENTS))
        cand = np.setdiff1d(feats.loc[pf, "feature_id"].to_numpy(), list(top))
        bg = rng.choice(cand, size=min(1000, len(cand)), replace=False)
        ev = pd.DataFrame({"feature_id": list(top) + list(bg), "set": ["top"] * len(top) + ["background"] * len(bg)})
        ev["symbol"] = ev["feature_id"].map(feats.set_index("feature_id")["symbol"])
        G = gex_by_symbol(cfg, fids_s, sorted(ev["symbol"].dropna().unique()))
        ev = ev[ev["symbol"].isin(G.columns)].reset_index(drop=True)
        X = feature_values(cfg, view, fids_s, ev["feature_id"])
        r = []
        for j, s in enumerate(ev["symbol"]):
            x, g = X[:, j], G[s].to_numpy()
            ok = np.isfinite(x) & np.isfinite(g)
            r.append(spearmanr(x[ok], g[ok]).statistic if ok.sum() >= 50 and np.nanstd(x[ok]) > 0 else np.nan)
        rho_rows.append(ev.assign(view=view, rho_host_expression=r))
    rho = pd.concat(rho_rows, ignore_index=True)
    rho["abs_rho"] = rho["rho_host_expression"].abs()
    rho.to_csv(out / "event_vs_host_expression.csv", index=False)
    rsum = rho.groupby(["view", "set"])["abs_rho"].agg(["median", lambda s: np.mean(s > 0.5), "size"])
    rsum.columns = ["median_abs_rho", "share_abs_rho_gt_0.5", "n"]
    rsum.to_csv(out / "event_vs_host_expression_summary.csv")
    log("event vs host expression:\n" + rsum.to_string())

    # ---------------- 2b. top-gene overlap, AS vs GEX, per class ----------------
    gex_genes = set(pd.read_parquet(cfg.matrix_dir / "gex" / "features.parquet", columns=["symbol"])["symbol"])
    ov = []
    for k, cls in enumerate(classes):
        gtop = set(gene_ranking(full_importance(cfg, "cancer_type", "gex", PAN, class_index=k))["symbol"].head(100))
        for view in as_views:
            feats = pd.read_parquet(cfg.matrix_dir / view / "features.parquet", columns=["symbol"])
            pf = pd.read_parquet(cfg.matrix_dir / view / "prefilter.parquet", columns=[PAN])[PAN].to_numpy()
            uni = set(feats.loc[pf, "symbol"].dropna()) & gex_genes
            atop = set(gene_ranking(full_importance(cfg, "cancer_type", view, PAN, class_index=k))["symbol"]
                       .dropna().head(100)) & uni
            g_in = gtop & uni
            ov.append({"class": cls, "view": view, "n_gex_top_in_universe": len(g_in), "n_as_top": len(atop),
                       "overlap": len(g_in & atop), "expected": len(g_in) * len(atop) / max(len(uni), 1),
                       "genes": ";".join(sorted(g_in & atop))})
    ov = pd.DataFrame(ov)
    ov.to_csv(out / "top_gene_overlap_gex_vs_as.csv", index=False)


    # ---------------- 3. lineage atlas of cassette exons ----------------
    feats = pd.read_parquet(cfg.matrix_dir / "se" / "features.parquet", columns=["feature_id", "symbol", "region"])
    picks, seen = [], set()
    for k, cls in enumerate(classes):
        imp = full_importance(cfg, "cancer_type", "se", PAN, class_index=k)
        n = 0
        for f, s in zip(imp["feature_id"], imp["symbol"]):
            if f in seen:
                continue
            picks.append({"class": cls, "feature_id": f, "symbol": s})
            seen.add(f)
            n += 1
            if n == ATLAS_PER_CLASS:
                break
    picks = pd.DataFrame(picks)
    X = feature_values(cfg, "se", tab["file_id"], picks["feature_id"])
    med = pd.DataFrame({c: np.nanmedian(X[y == k], axis=0) for k, c in enumerate(classes)},
                       index=picks["feature_id"]) * 100
    med.to_csv(out / "atlas_se_median_psi_by_class.csv")
    picks.to_csv(out / "atlas_se_events.csv", index=False)
    log("lineage done")


if __name__ == "__main__":
    main()
