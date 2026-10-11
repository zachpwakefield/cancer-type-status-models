"""Step 10 (B-03): which pathways do the models rely on?

For every (task, view, cohort) model, and for every cancer-type class, the top N genes by SHAP
importance (max feature mean|SHAP| per gene, held-out samples only) are tested for
over-representation in MSigDB Hallmark, Reactome and GO-BP gene sets. The universe is every
gene with >= 1 measurable feature in that view and cohort (extraction pre-filter). AS views
measure many events per gene, so long, highly expressed genes are more likely to be picked; the
GOseq correction (Wallenius test with a probability-weighting function of events per gene)
removes that bias. GEX has one feature per gene, where the test reduces to the hypergeometric.

Sensitivity (cancer type, AS views): the per-class AS results are dominated by translation terms
driven by alternative 5'/3' splice sites in ribosomal-protein genes (deeply covered, so precisely
measured in every tissue). The per-class AS tests are repeated with ribosomal-protein genes
(RPL*/RPS*, regex in config) removed from both the ranking and the universe.

Outputs: results/biology/03_pathways/ (per-model tables and summaries; figures: scripts/15_bio_figures.py)
Usage: python scripts/10_bio_pathways.py [--workers 6] [--rerun]
"""
from __future__ import annotations

import argparse
import re
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from _common import log

from ctsm.biology import add_q, bio_config, full_importance, gene_ranking, goseq_ora, n_testable, read_gmt, ref_path
from ctsm.config import PAN, load_config

_GS = {}


def gene_sets(cfg, coll):
    if coll not in _GS:
        _GS[coll] = read_gmt(ref_path(cfg, "gene_sets", coll))
    return _GS[coll]


def universe(cfg, view: str, cohort: str) -> pd.DataFrame:
    d = cfg.matrix_dir / view
    f = pd.read_parquet(d / "features.parquet", columns=["feature_id", "symbol"])
    keep = pd.read_parquet(d / "prefilter.parquet", columns=[cohort])[cohort].to_numpy()
    return f[keep].groupby("symbol").size().rename("n_events").reset_index()


def run_one(args) -> pd.DataFrame:
    task, view, cohort, k, label, exclude = args
    cfg = load_config()
    pc = bio_config()["pathways"]
    imp = full_importance(cfg, task, view, cohort, class_index=k)
    genes = gene_ranking(imp).dropna(subset=["symbol"])
    uni = universe(cfg, view, cohort)
    if exclude:
        rx = re.compile(exclude)
        genes = genes[~genes["symbol"].str.match(rx)]
        uni = uni[~uni["symbol"].str.match(rx)]
    top = set(genes["symbol"].head(pc["top_genes"]))
    usym = set(uni["symbol"])
    frames, ntest = [], {}
    for coll in pc["collections"]:
        gs = gene_sets(cfg, coll)
        ntest[coll] = n_testable(usym, gs, pc["min_set_size"], pc["max_set_size"])
        frames.append(goseq_ora(top, uni, gs, pc["min_set_size"], pc["max_set_size"], pc["min_overlap"], coll))
    res = pd.concat([f for f in frames if len(f)], ignore_index=True) if any(len(f) for f in frames) else pd.DataFrame()
    if len(res):
        res = add_q(res, ntest)
    res = res.assign(task=task, view=view, cohort=cohort, label=label, n_top=len(top & usym), universe=len(usym),
                     excluded="ribosomal" if exclude else "none")
    out = cfg.results_dir / "biology" / "03_pathways" / ("tables_noRP" if exclude else "tables")
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{task}__{view}__{cohort}" + ("" if k is None else f"__{label}")
    keep = res[(res.get("q", pd.Series(dtype=float)) < 0.25)] if len(res) else res
    keep.sort_values("p_wallenius").to_csv(out / f"{tag}.csv", index=False)
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--rerun", action="store_true", help="recompute even if the cached results exist")
    a = ap.parse_args()
    cfg = load_config()
    pc = bio_config()["pathways"]
    out = cfg.results_dir / "biology" / "03_pathways"
    out.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = list(scope.loc[scope.tumor_status_in_scope, "project"])
    tab = pd.read_parquet(cfg.data_dir / "tasks" / f"cancer_type__{PAN}.parquet")
    classes = tab.drop_duplicates("y").sort_values("y")["label"].tolist()
    jobs = [("tumor_status", v, c, None, c, None) for v in cfg.views for c in cancers]
    jobs += [("cancer_type", v, PAN, None, "all", None) for v in cfg.views]
    jobs += [("cancer_type", v, PAN, k, classes[k], None) for v in cfg.views for k in range(len(classes))]
    sens = [("cancer_type", v, PAN, k, classes[k], pc["ribosomal_regex"])
            for v in cfg.views if v != "gex" for k in range(len(classes))]
    parts = []
    for name, jj in [("pathways_all", jobs), ("pathways_noRP", sens)]:
        cache = cfg.data_dir / "bio_cache" / f"{name}.parquet"
        if cache.exists() and not a.rerun:
            parts.append(pd.read_parquet(cache))
            continue
        with ProcessPoolExecutor(a.workers) as ex:
            res = list(ex.map(run_one, jj))
        r = pd.concat([x for x in res if len(x) and "q" in x], ignore_index=True)
        r.to_parquet(cache, index=False)
        parts.append(r)
    allres, norp = parts
    sig = allres[allres["q"] < pc["fdr"]]
    sig.drop(columns="genes").to_csv(out / "significant_terms_all_models.csv", index=False)
    log(f"{len(jobs)} models tested; {len(sig)} significant (model, term) pairs")

    # ---------- summary 1: tumour status, Hallmark convergence across AS views and GEX ----------
    ts = allres[(allres.task == "tumor_status") & (allres.collection == "hallmark")]
    ts_sig = ts[ts.q < pc["fdr"]]
    as_count = (ts_sig[ts_sig.view != "gex"].groupby(["term", "cohort"]).size().unstack(fill_value=0)
                .reindex(columns=cancers, fill_value=0))
    gex_sig = (ts_sig[ts_sig.view == "gex"].assign(v=1).pivot_table(index="term", columns="cohort", values="v",
                                                                    aggfunc="max").reindex(columns=cancers).fillna(0))
    terms = sorted(set(as_count.index) | set(gex_sig.index))
    as_count = as_count.reindex(terms, fill_value=0)
    gex_sig = gex_sig.reindex(terms).fillna(0)
    score = as_count.sum(axis=1) + 2 * gex_sig.sum(axis=1)
    keep = score.sort_values(ascending=False).index[:22]
    as_count, gex_sig = as_count.loc[keep], gex_sig.loc[keep]
    pd.concat({"n_AS_views_significant": as_count, "GEX_significant": gex_sig}, axis=1).to_csv(
        out / "tumor_status_hallmark_convergence.csv")

    # ---------- summary 2: GO-BP recurrence across tumour-status models (table only) ----------
    go = allres[(allres.task == "tumor_status") & (allres.collection == "go_bp") & (allres.q < pc["fdr"])]
    rec_as = go[go.view != "gex"].groupby("term")["cohort"].nunique()
    rec_gex = go[go.view == "gex"].groupby("term")["cohort"].nunique()
    pd.DataFrame({"n_cancers_any_AS_view": rec_as, "n_cancers_GEX": rec_gex}).fillna(0).astype(int).sort_values(
        ["n_cancers_any_AS_view", "n_cancers_GEX"], ascending=False).to_csv(out / "tumor_status_gobp_recurrence.csv")

    # ---------- summary 3: pathway coherence = share of models whose top genes hit >= 1 term ----------
    def frac_with_term(df, task, units, col):
        hit = df[(df.task == task) & (df.q < pc["fdr"])].groupby("view")[col].nunique()
        return (hit.reindex(cfg.views).fillna(0) / len(units))
    ct_all = allres[allres.label != "all"]
    coh = pd.DataFrame({"tumour status (15 cancers)": frac_with_term(allres, "tumor_status", cancers, "cohort"),
                        "cancer type, per class (34)": frac_with_term(ct_all, "cancer_type", classes, "label"),
                        "cancer type, per class, ribosomal genes removed": frac_with_term(norp, "cancer_type", classes,
                                                                                            "label")})
    coh.loc["gex", "cancer type, per class, ribosomal genes removed"] = np.nan
    coh.to_csv(out / "pathway_coherence_by_view.csv")

    # ---------- summary 4: cancer type, best term per class for GEX and AS (with / without RP genes) ----------
    ct = ct_all[ct_all.q < pc["fdr"]]
    ctn = norp[norp.q < pc["fdr"]]
    rows = []
    for cls in classes:
        r = {"class": cls}
        for fam, g in [("GEX", ct[(ct.view == "gex") & (ct.label == cls)]), ("AS", ct[(ct.view != "gex") & (ct.label == cls)]),
                       ("AS_noRP", ctn[ctn.label == cls])]:
            g = g.sort_values("p_wallenius")
            r[f"{fam}_n_terms"] = len(g)
            r[f"{fam}_top_terms"] = "; ".join(f"{t} ({v})" for t, v in zip(g["term"].head(3), g["view"].head(3)))
        rows.append(r)
    pd.DataFrame(rows).to_csv(out / "cancer_type_top_terms_by_class.csv", index=False)
    transl = ct[(ct.view != "gex")]["term"].str.contains("TRANSLATION|RIBOSOM|EIF2AK4|SRP_DEPENDENT|NONSENSE|SELENO|"
                                                           "INFLUENZA|PEPTIDE_CHAIN|RRNA")
    log(f"cancer type AS significant terms: {len(transl)}, translation-type {transl.mean():.2f}; "
        f"after removing RP genes: {len(ctn)} terms")
    log("pathways done")


if __name__ == "__main__":
    main()
