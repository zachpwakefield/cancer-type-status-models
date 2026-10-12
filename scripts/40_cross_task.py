"""Step 40 (D-48..D-51): how do tumour status, cancer type and subtype relate?

Nothing is fitted. Everything comes from stored out-of-fold predictions, held-out SHAP
importance and the data; predictions and sources are in config/followup.yaml (cross_task).

1. gene_overlap   top-100 genes of the status model, the type model's own-class importance and
                  each subtype cohort's model, per view; pairwise overlap within the genes both
                  tasks could select; hypergeometric; BH over all (pair, cohort, view)    (D-48)
2. type_by_subtype   the cancer-type model's out-of-fold probabilities by subtype class, and the
                  predictions T1-T13 (summed probability over target projects, Mann-Whitney)  (D-49)
3. status_by_subtype the tumour-status model's out-of-fold P(tumour) by subtype class, the
                  predictions S1-S6, and C1 (P(tumour) against stromal and immune scores)     (D-50)
4. exon_character protein character of the subtype models' top cassette exons against the
                  background of D-17, next to the other two tasks' exons                       (D-51)

Models: GEX, AS ensemble (mean of the 8 AS views), GEX + AS (0.5 GEX + 0.5 AS ensemble).
Outputs: results/across_tasks/ (aggregates only; figures: scripts/41_figures_cross_task.py)
Usage: python scripts/40_cross_task.py [gene_overlap type_by_subtype status_by_subtype exon_character]
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from _common import REPO, log
from scipy.stats import fisher_exact, hypergeom, kruskal, mannwhitneyu, spearmanr, wilcoxon

from ctsm.biology import ProteinMap, bh, bio_config, full_importance, gene_ranking, gex_by_symbol
from ctsm.config import PAN, load_config
from ctsm.subtype import TASK, cohort_specs, load_tables, subtype_config

CFG = load_config()
SCFG = subtype_config()
XC = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["cross_task"]
OUT = CFG.results_dir / "across_tasks"
TABLES = load_tables(SCFG)
AS_VIEWS = [v for v in CFG.views if v != "gex"]
MODELS = ["GEX", "AS ensemble", "GEX + AS"]
SCOPE = pd.read_csv(CFG.results_dir / "data" / "task_scope.csv")
STATUS_CANCERS = set(SCOPE.loc[SCOPE.tumor_status_in_scope, "project"])
_spec = importlib.util.spec_from_file_location("ev11", Path(__file__).with_name("11_bio_events.py"))
_ev11 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ev11)


def label_type(cohort: str) -> str:
    perf = pd.read_csv(SCFG.results_dir / TASK / "summary" / "performance.csv")
    return perf.drop_duplicates("cohort").set_index("cohort")["label_type"][cohort]


def project_of(cohort: str) -> str:
    return TABLES[cohort]["project"].iloc[0]


# ============================================================ 1. gene overlap
def type_classes() -> list[str]:
    cols = pd.read_parquet(CFG.run_dir("cancer_type", "gex", PAN) / "oof_predictions.parquet").columns
    return [c[2:] for c in cols if c.startswith("p_")]


def universe(matrix_dir: Path, view: str, column: str) -> set[str]:
    F = pd.read_parquet(matrix_dir / view / "features.parquet", columns=["gene_id"])
    P = pd.read_parquet(matrix_dir / view / "prefilter.parquet", columns=[column])[column].to_numpy()
    return set(F.loc[P, "gene_id"])


def top_genes(imp: pd.DataFrame, uni: set[str], k: int) -> list[str]:
    g = gene_ranking(imp)
    g = g[g["gene_id"].isin(uni) & (g["importance"] > 0)]
    return g["gene_id"].head(k).tolist()


def compare(a: list[str], b: list[str], uni: set[str]) -> dict:
    A, B = set(a), set(b)
    n, k = len(uni), len(A & B)
    exp = len(A) * len(B) / n
    return {"n_universe": n, "n_a": len(A), "n_b": len(B), "n_both": k, "expected": exp,
            "fold_enrichment": k / exp if exp else np.nan, "jaccard": k / len(A | B) if A | B else np.nan,
            "p": float(hypergeom.sf(k - 1, n, len(A), len(B)))}


def step_gene_overlap() -> None:
    k = XC["gene_overlap"]["top_genes"]
    classes = type_classes()
    symbols = {}
    rows, triple = [], []
    for view in CFG.views:
        feats = pd.read_parquet(CFG.matrix_dir / view / "features.parquet", columns=["gene_id", "symbol"])
        symbols.update(dict(zip(feats.gene_id, feats.symbol)))
        u_type = universe(CFG.matrix_dir, view, PAN)
        status_done = set()
        for cohort in TABLES:
            proj = project_of(cohort)
            u_sub = universe(SCFG.matrix_dir, view, cohort)
            imp_sub = pd.read_parquet(SCFG.data_dir / "bio_cache" / "importance" / f"{view}__{cohort}__all.parquet")
            imp_type = full_importance(CFG, "cancer_type", view, PAN, class_index=classes.index(proj))
            meta = {"cohort": cohort, "project": proj, "label_type": label_type(cohort), "view": view}
            u = u_sub & u_type
            sub_t, typ_t = top_genes(imp_sub, u, k), top_genes(imp_type, u, k)
            rows.append({"pair": "subtype-type", **meta, **compare(sub_t, typ_t, u)})
            if proj not in STATUS_CANCERS:
                continue
            u_st = universe(CFG.matrix_dir, view, proj)
            imp_st = full_importance(CFG, "tumor_status", view, proj)
            u = u_sub & u_st
            rows.append({"pair": "subtype-status", **meta,
                         **compare(top_genes(imp_sub, u, k), top_genes(imp_st, u, k), u)})
            if proj not in status_done:
                u = u_st & u_type
                rows.append({"pair": "status-type", "cohort": proj, "project": proj, "label_type": "",
                             "view": view, **compare(top_genes(imp_st, u, k), top_genes(imp_type, u, k), u)})
                status_done.add(proj)
            u = u_sub & u_st & u_type
            common = set(top_genes(imp_sub, u, k)) & set(top_genes(imp_st, u, k)) & set(top_genes(imp_type, u, k))
            for g in sorted(common):
                triple.append({"cohort": cohort, "project": proj, "view": view, "gene_id": g, "symbol": symbols.get(g, "")})
    d = pd.DataFrame(rows)
    d["q"] = bh(d["p"].to_numpy())
    d.to_csv(OUT / "gene_overlap.csv", index=False)
    pd.DataFrame(triple).to_csv(OUT / "genes_in_all_three.csv", index=False)
    tests = []
    for pair, g in d.groupby("pair"):
        lf = np.log2(g["fold_enrichment"].replace(0, np.nan)).dropna()
        tests.append({"id": "G1", "pair": pair, "n": len(g), "median_fold_enrichment": g["fold_enrichment"].median(),
                      "share_q05": float((g["q"] < 0.05).mean()), "median_jaccard": g["jaccard"].median(),
                      "p_wilcoxon_log2fe_vs_0": float(wilcoxon(lf).pvalue) if len(lf) > 5 else np.nan})
    h = d[d.label_type == "histology"].pivot_table(index=["cohort", "view"], columns="pair", values="jaccard").dropna()
    diff = h["subtype-type"] - h["subtype-status"]
    tests.append({"id": "G2", "pair": "histology: subtype-type minus subtype-status Jaccard", "n": len(diff),
                  "median_difference": float(diff.median()), "share_positive": float((diff > 0).mean()),
                  "p_wilcoxon": float(wilcoxon(diff).pvalue)})
    pd.DataFrame(tests).to_csv(OUT / "gene_overlap_tests.csv", index=False)
    log("gene overlap:\n" + pd.DataFrame(tests).round(4).to_string(index=False))


# ============================================================ 2-3. probabilities by subtype
def oof_models(task: str, cohort: str, cols: list[str]) -> dict[str, pd.DataFrame]:
    P = {v: pd.read_parquet(CFG.run_dir(task, v, cohort) / "oof_predictions.parquet").set_index("file_id")
         for v in CFG.views}
    idx = P["gex"].index
    g = P["gex"].loc[idx, cols]
    a = sum(P[v].loc[idx, cols] for v in AS_VIEWS) / len(AS_VIEWS)
    return {"GEX": g, "AS ensemble": a, "GEX + AS": 0.5 * g + 0.5 * a}


def group_test(x: np.ndarray, lab: np.ndarray, hi: list[str], lo: list[str]) -> dict:
    a, b = x[np.isin(lab, hi)], x[np.isin(lab, lo)]
    r = mannwhitneyu(a, b, alternative="two-sided")
    return {"n_higher": len(a), "n_lower": len(b), "median_higher": float(np.median(a)),
            "median_lower": float(np.median(b)), "mean_higher": float(np.mean(a)), "mean_lower": float(np.mean(b)),
            "auroc": float(r.statistic / (len(a) * len(b))), "p": float(r.pvalue)}


def verdict(d: pd.DataFrame) -> pd.DataFrame:
    d["q"] = d.groupby("model")["p"].transform(lambda p: bh(p.to_numpy()))
    d["verdict"] = np.where(d["q"] >= 0.05, "not resolved", np.where(d["auroc"] > 0.5, "supported", "against"))
    return d


def step_type_by_subtype() -> None:
    classes = type_classes()
    M = oof_models("cancer_type", PAN, [f"p_{c}" for c in classes])
    mean_rows, call_rows, pred_rows, sens_rows = [], [], [], []
    for model, P in M.items():
        for cohort, tab in TABLES.items():
            proj = project_of(cohort)
            X = P.loc[tab["file_id"]].to_numpy()
            lab = tab["label"].to_numpy()
            called = np.array(classes)[X.argmax(1)]
            for c in sorted(set(lab)):
                m = lab == c
                mean_rows.append({"model": model, "cohort": cohort, "class": c, "n": int(m.sum()),
                                  **{p: float(v) for p, v in zip(classes, X[m].mean(0))}})
                other = called[m] != proj
                vc = pd.Series(called[m][other]).value_counts()
                call_rows.append({"model": model, "cohort": cohort, "project": proj, "class": c, "n": int(m.sum()),
                                  "n_called_other": int(other.sum()), "share_called_other": float(other.mean()),
                                  "most_common_other": vc.index[0] if len(vc) else "",
                                  "n_most_common_other": int(vc.iloc[0]) if len(vc) else 0})
        for pr in XC["type_probability_by_subtype"]["predictions"]:
            tab = TABLES[pr["cohort"]]
            X = P.loc[tab["file_id"]]
            lab = tab["label"].to_numpy()
            x = X[[f"p_{t}" for t in pr["targets"]]].sum(1).to_numpy()
            base = {"model": model, "id": pr["id"], "cohort": pr["cohort"], "targets": "+".join(pr["targets"]),
                    "higher_in": "+".join(pr["higher_in"]), "lower_in": "+".join(pr["lower_in"]), "source": pr["source"]}
            pred_rows.append({**base, **group_test(x, lab, pr["higher_in"], pr["lower_in"])})
            # sensitivity, not pre-specified: share of the probability NOT given to the own project
            own = X[f"p_{project_of(pr['cohort'])}"].to_numpy()
            sens_rows.append({**base, **group_test(x / np.clip(1 - own, 1e-9, None), lab, pr["higher_in"], pr["lower_in"])})
    pd.DataFrame(mean_rows).to_csv(OUT / "type_probability_by_subtype_class.csv", index=False)
    pd.DataFrame(call_rows).to_csv(OUT / "type_calls_by_subtype_class.csv", index=False)
    d = verdict(pd.DataFrame(pred_rows))
    d.to_csv(OUT / "type_predictions.csv", index=False)
    s = verdict(pd.DataFrame(sens_rows))
    s.to_csv(OUT / "type_predictions_offtarget_share.csv", index=False)
    log("type predictions:\n" + d[["model", "id", "cohort", "targets", "auroc", "p", "q", "verdict"]].round(4).to_string(index=False))


def step_status_by_subtype() -> None:
    rc = bio_config()["regulators"]
    pred_rows, desc_rows, comp_rows = [], [], []
    projects = sorted({project_of(c) for c in TABLES} & STATUS_CANCERS)
    for proj in projects:
        M = oof_models("tumor_status", proj, ["p_Tumour"])
        tab_st = pd.read_parquet(CFG.data_dir / "tasks" / f"tumor_status__{proj}.parquet")
        tum = tab_st.loc[tab_st["y"] == 1, "file_id"]
        G = gex_by_symbol(CFG, tum, rc["stromal_genes"] + rc["immune_genes"])
        scores = {"stromal": G[[g for g in rc["stromal_genes"] if g in G]].mean(1),
                  "immune": G[[g for g in rc["immune_genes"] if g in G]].mean(1)}
        for model, P in M.items():
            pt = P.loc[tum, "p_Tumour"].to_numpy()
            for name, sc in scores.items():
                r = spearmanr(pt, sc.loc[tum].to_numpy())
                comp_rows.append({"id": "C1", "model": model, "project": proj, "score": name, "n": len(tum),
                                  "rho": float(r.statistic), "p": float(r.pvalue)})
            for cohort in [c for c in TABLES if project_of(c) == proj]:
                tab = TABLES[cohort]
                x = P.loc[tab["file_id"], "p_Tumour"].to_numpy()
                lab = tab["label"].to_numpy()
                groups = [x[lab == c] for c in sorted(set(lab))]
                row = {"model": model, "cohort": cohort, "project": proj, "n": len(x),
                       "p_kruskal": float(kruskal(*groups).pvalue)}
                for c in sorted(set(lab)):
                    v = x[lab == c]
                    row.update({f"n|{c}": int(len(v)), f"median|{c}": float(np.median(v)),
                                f"share_below_0.5|{c}": float((v < 0.5).mean())})
                desc_rows.append(row)
        for pr in XC["status_probability_by_subtype"]["predictions"]:
            if project_of(pr["cohort"]) != proj:
                continue
            tab = TABLES[pr["cohort"]]
            for model, P in M.items():
                x = P.loc[tab["file_id"], "p_Tumour"].to_numpy()
                pred_rows.append({"model": model, "id": pr["id"], "cohort": pr["cohort"],
                                  "higher_in": "+".join(pr["higher_in"]), "lower_in": "+".join(pr["lower_in"]),
                                  "source": pr["source"],
                                  **group_test(x, tab["label"].to_numpy(), pr["higher_in"], pr["lower_in"])})
    d = verdict(pd.DataFrame(pred_rows)).sort_values(["model", "id"])
    d.to_csv(OUT / "status_predictions.csv", index=False)
    desc = pd.DataFrame(desc_rows)
    desc["q_kruskal"] = desc.groupby("model")["p_kruskal"].transform(lambda p: bh(p.to_numpy()))
    desc.to_csv(OUT / "status_probability_by_subtype_class.csv", index=False)
    comp = pd.DataFrame(comp_rows)
    comp["q"] = comp.groupby("model")["p"].transform(lambda p: bh(p.to_numpy()))
    comp.to_csv(OUT / "status_probability_vs_composition.csv", index=False)
    log("status predictions:\n" + d[["model", "id", "cohort", "auroc", "p", "q", "verdict"]].round(4).to_string(index=False))
    log("composition:\n" + comp.round(3).to_string(index=False))


# ============================================================ 4. exon protein character
def step_exon_character() -> None:
    stored = pd.read_csv(CFG.results_dir / "biology" / "04_events" / "cassette_exon_consequences.csv")
    bg_name = "all measured exons (random 6,000)"
    top = set()
    for cohort in TABLES:
        imp = pd.read_parquet(SCFG.data_dir / "bio_cache" / "importance" / f"se__{cohort}__all.parquet")
        top |= set(imp["feature_id"].head(100))
    F = pd.read_parquet(SCFG.matrix_dir / "se" / "features.parquet", columns=["feature_id", "gene_id", "region"])
    ev = F[F.feature_id.isin(top)].reset_index(drop=True).assign(set="subtype, top events")
    s0, e0 = _ev11.se_coords(ev["region"])
    ev["length"] = e0 - s0 + 1
    ev["len_mult3"] = ev["length"] % 3 == 0
    pm = ProteinMap(CFG)
    cons = pd.DataFrame([pm.exon_consequence(g, int(s), int(e)) for g, s, e in zip(ev.gene_id, s0, e0)])
    ev = pd.concat([ev, cons.drop(columns=["exon_length", "strand"], errors="ignore")], axis=1)
    keep = ["feature_id", "gene_id", "region", "set", "length", "len_mult3", "annotated_exon", "coding",
            "fully_coding", "frame_preserving", "in_domain", "in_disordered", "in_transmembrane"]
    allx = pd.concat([stored[[c for c in keep if c in stored]], ev[[c for c in keep if c in ev]]], ignore_index=True)
    metrics = [("length divisible by 3", "len_mult3", None), ("exact GENCODE exon", "annotated_exon", None),
               ("protein-coding exon", "coding", "annotated_exon"),
               ("overlaps a Pfam/CDD domain", "in_domain", "fully_coding"),
               ("overlaps a disordered region", "in_disordered", "fully_coding"),
               ("overlaps a TM helix", "in_transmembrane", "fully_coding")]
    sets = [bg_name, "tumour vs normal, top events", "cancer type, top events", "subtype, top events"]
    rows = []
    for name, col, among in metrics:
        base = allx if among is None else allx[allx[among].eq(True)]
        b = base[base.set == bg_name][col].eq(True)
        for st in sets:
            x = base[base.set == st][col].eq(True)
            lo, hi = _ev11.wilson(int(x.sum()), len(x))
            r = {"metric": name, "among": among or "all", "set": st, "n": len(x), "k": int(x.sum()),
                 "proportion": x.mean(), "ci_lo": lo, "ci_hi": hi}
            if st != bg_name:
                odds, p = fisher_exact([[x.sum(), (~x).sum()], [b.sum(), (~b).sum()]])
                r.update(odds_ratio=odds, p_fisher=p)
            rows.append(r)
    pr = pd.DataFrame(rows)
    has = pr["p_fisher"].notna()
    pr.loc[has, "q_fisher"] = bh(pr.loc[has, "p_fisher"].to_numpy())
    pr.to_csv(OUT / "exon_protein_character.csv", index=False)
    ev[keep].to_csv(OUT / "subtype_top_exons.csv.gz", index=False)
    log("exon character:\n" + pr[["metric", "set", "n", "proportion", "odds_ratio", "q_fisher"]].round(4).to_string(index=False))


STEPS = {"gene_overlap": step_gene_overlap, "type_by_subtype": step_type_by_subtype,
         "status_by_subtype": step_status_by_subtype, "exon_character": step_exon_character}

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for s in (sys.argv[1:] or list(STEPS)):
        STEPS[s]()
