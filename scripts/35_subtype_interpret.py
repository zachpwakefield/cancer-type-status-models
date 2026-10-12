"""Step 35: what the subtype models learned (D-41, D-45). Tables only; figures are step 36.

No model is refitted and nothing here feeds back into a model or a label. Everything is computed
from held-out SHAP values, the out-of-fold predictions and the data.

1. importance     every feature of every (cohort, view): mean |SHAP| over all held-out patients
                  and over the patients of each class (cached locally); top 30 per class committed
2. as_vs_gex      top 50 genes of each AS model against the top 50 of the GEX model of the same
                  cohort, within the genes both views measure in that cohort; hypergeometric
3. markers        the pre-specified literature markers: expression AUROC between the two named
                  classes (no model), and the gene's rank in the GEX and AS models; LGG 1p/19q dosage
4. event_audit    the top 5 events per class and AS view, and each literature marker's best event
                  (run after markers): overlap of the event's segment with
                  another GENCODE v45 gene's exon (>= 50 bp; antisense / same-strand / readthrough),
                  Spearman with that gene's and the host's log2 CPM over the cohort; per class the
                  median, 10th and 90th percentile and the share of patients with a value
5. cross_family   histology models' out-of-fold calls against the molecular class of the same patients
6. error_markers  the literature markers' expression in each (true class, called class) cell
7. shared_classes a class found in several cohorts (MSI, CIN, GS, squamous, adenocarcinoma): each
                  audited top event of it, measured in every cohort with that class (no model)
8. msi_confounds  each audited top event of the MSI class, in each MSI cohort: Spearman with immune
                  (PTPRC, CD8A, CIITA) and MLH1 log2 CPM, and MSI-vs-rest AUROC within the half of the
                  cohort with the lower PTPRC (is the event a readout of immune infiltrate?)

Outputs: results/subtype/interpret/*.csv (aggregates); <data_dir>/subtype/bio_cache/ (local).
Usage: python scripts/35_subtype_interpret.py [importance as_vs_gex markers event_audit cross_family
                                              error_markers shared_classes msi_confounds]
"""
from __future__ import annotations

import gzip
import re
import sys

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import hypergeom, mannwhitneyu, spearmanr

from ctsm.biology import bh, feature_values, gene_ranking, gex_by_symbol
from ctsm.explain import importance_table, load_folds
from ctsm.figstyle import VIEW_NAMES
from ctsm.matrix_io import load_view
from ctsm.subtype import TASK, cohort_specs, load_tables, subtype_config

SCFG = subtype_config()
SPECS = cohort_specs(SCFG.spec)
ITP = SCFG.spec["interpretation"]
OUT = SCFG.results_dir / TASK / "interpret"
CACHE = SCFG.data_dir / "bio_cache"
AS_VIEWS = [v for v in SCFG.views if v != "gex"]
TABLES = load_tables(SCFG)
CEN = {"chr1": 123_400_000, "chr19": 26_200_000}        # hg38 centromeres (p/q boundary)
MIN_CELL = 5                                              # smallest group a committed median is given for


def done(cohort: str, view: str) -> bool:
    return (SCFG.run_dir(TASK, view, cohort) / "job_done.json").exists()


def classes_of(cohort: str) -> list[str]:
    t = TABLES[cohort]
    return list(t.drop_duplicates("y").sort_values("y")["label"])


# ============================================================ 1. importance
def features(view: str) -> pd.DataFrame:
    return pd.read_parquet(SCFG.matrix_dir / view / "features.parquet")[
        ["feature_id", "label", "symbol", "gene_id", "chrom", "strand", "region"]]


def importance(cohort: str, view: str, k: int | None = None) -> pd.DataFrame:
    """Every feature selected by any fold; k = class index (patients of that class only). For two
    classes the direction is flipped for the first class, so > 0 always means evidence for it."""
    tag = "all" if k is None else f"class{k}"
    p = CACHE / "importance" / f"{view}__{cohort}__{tag}.parquet"
    if p.exists():
        return pd.read_parquet(p)
    tab = TABLES[cohort]
    folds = load_folds(SCFG.run_dir(TASK, view, cohort), tab["y"].to_numpy(), SCFG.study["splits"]["outer_folds"])
    mask = (lambda f: slice(None)) if k is None else (lambda f: f["y"] == k)
    t = importance_table(folds, features(view), mask_of=mask, n_direction=300)
    if k is not None and len(classes_of(cohort)) == 2 and k == 0:
        t["direction_spearman"] = -t["direction_spearman"]
    p.parent.mkdir(parents=True, exist_ok=True)
    t.to_parquet(p, index=False)
    return t


def step_importance() -> None:
    rows = []
    for cohort in TABLES:
        for view in SCFG.views:
            f = SCFG.result_dir(TASK, view, cohort) / "shap_top_features_by_class.csv"
            if not f.exists():
                continue
            rows.append(pd.read_csv(f).assign(cohort=cohort, view=view))
            importance(cohort, view)
            for k in range(len(classes_of(cohort))):
                importance(cohort, view, k)
    t = pd.concat(rows, ignore_index=True)
    t.to_csv(OUT / "top_features_by_class.csv", index=False)
    log(f"importance: {len(rows)} jobs")


# ============================================================ 2. AS vs GEX genes
def cohort_genes(cohort: str, view: str) -> set[str]:
    """Genes with at least one feature passing the cohort's pre-filter in this view."""
    F = pd.read_parquet(SCFG.matrix_dir / view / "features.parquet", columns=["gene_id"])
    P = pd.read_parquet(SCFG.matrix_dir / view / "prefilter.parquet", columns=[cohort])[cohort].to_numpy()
    return set(F.loc[P, "gene_id"])


def step_as_vs_gex() -> None:
    n_top = int(ITP["top_genes_overlap"])
    rows = []
    for cohort in TABLES:
        if not done(cohort, "gex"):
            continue
        g_gex = gene_ranking(importance(cohort, "gex"))
        U_gex = cohort_genes(cohort, "gex")
        for view in AS_VIEWS:
            if not done(cohort, view):
                continue
            U = U_gex & cohort_genes(cohort, view)
            a = gene_ranking(importance(cohort, view))
            top_as = list(a.loc[a.gene_id.isin(U), "gene_id"].head(n_top))
            top_gx = list(g_gex.loc[g_gex.gene_id.isin(U), "gene_id"].head(n_top))
            both = set(top_as) & set(top_gx)
            exp = len(top_as) * len(top_gx) / len(U)
            sym = a.set_index("gene_id")["symbol"].to_dict()
            rows.append({"cohort": cohort, "view": VIEW_NAMES[view], "n_universe": len(U), "n_top_as": len(top_as),
                         "n_top_gex": len(top_gx), "overlap": len(both), "expected": exp,
                         "fold_enrichment": len(both) / exp if exp else np.nan,
                         "p_hypergeom": float(hypergeom(len(U), len(top_gx), len(top_as)).sf(len(both) - 1)),
                         "shared_genes": ";".join(sorted(sym.get(g, g) for g in both))})
    pd.DataFrame(rows).to_csv(OUT / "as_vs_gex_genes.csv", index=False)
    log(f"as_vs_gex: {len(rows)} cohort x view pairs")


# ============================================================ 3. literature markers
def gene_positions() -> pd.DataFrame:
    ex = pd.read_parquet(SCFG.base_data_dir / "cache" / "gencode_v45_exons.parquet",
                         columns=["chrom", "start", "end", "strand", "gene_id", "gene_name"])
    ex["gene_id"] = ex["gene_id"].str.split(".").str[0]
    return ex.groupby("gene_id").agg(chrom=("chrom", "first"), start=("start", "min"), end=("end", "max"),
                                     strand=("strand", "first"), gene_name=("gene_name", "first"))


def arm_dosage(tab: pd.DataFrame) -> pd.DataFrame:
    """Per patient: median over the genes of 1p (19q) of log2 CPM minus the gene's cohort median."""
    M, F = load_view(SCFG.matrix_dir, "gex", tab["file_id"], tab["cohort_name"].iloc[0])
    X = M.get()
    pos = gene_positions().reindex(F["gene_id"])
    centred = X - np.nanmedian(X, axis=0)
    out = {}
    for arm, (ch, side) in {"1p": ("chr1", "p"), "19q": ("chr19", "q"), "1q": ("chr1", "q"),
                            "19p": ("chr19", "p")}.items():
        on = (pos["chrom"].to_numpy() == ch) & ((pos["start"].to_numpy() < CEN[ch]) if side == "p"
                                                 else (pos["start"].to_numpy() >= CEN[ch]))
        out[arm] = np.nanmedian(centred[:, on], axis=1)
        out[f"n_genes_{arm}"] = int(on.sum())
    return pd.DataFrame({k: v for k, v in out.items() if not k.startswith("n_")}, index=tab["file_id"]), \
        {k: v for k, v in out.items() if k.startswith("n_")}


def rank_of(imp: pd.DataFrame, symbol: str) -> float:
    g = gene_ranking(imp).reset_index(drop=True)
    hit = np.flatnonzero(g["symbol"].to_numpy() == symbol)
    return float(hit[0] + 1) if len(hit) else np.nan


def step_markers() -> None:
    top = int(ITP["marker_top_rank"])
    rows, dos_rows = [], []
    for cohort, markers in SCFG.spec["literature_markers"].items():
        if cohort not in TABLES or not done(cohort, "gex"):
            continue
        tab = TABLES[cohort].assign(cohort_name=cohort)
        cls = classes_of(cohort)
        for mk in markers:
            a, b = mk["higher_in"], mk["vs"]
            ia, ib = cls.index(a), cls.index(b)
            ya, yb = tab["y"].to_numpy() == ia, tab["y"].to_numpy() == ib
            if mk["gene"] == "1p/19q dosage":
                D, n_genes = arm_dosage(tab)
                for arm in ("1p", "19q", "1q", "19p"):
                    x = D[arm].to_numpy()
                    u = mannwhitneyu(x[ya], x[yb])
                    dos_rows.append({"cohort": cohort, "arm": arm, "n_genes": n_genes[f"n_genes_{arm}"],
                                     "higher_in": a, "vs": b, "n_a": int(ya.sum()), "n_b": int(yb.sum()),
                                     "median_a": float(np.median(x[ya])), "median_b": float(np.median(x[yb])),
                                     "auroc_a_over_b": u.statistic / (ya.sum() * yb.sum()), "p_mwu": float(u.pvalue),
                                     **{f"median_{c}": float(np.median(x[tab['y'].to_numpy() == i]))
                                        for i, c in enumerate(cls)}})
                # are the GEX model's top genes on 1p or 19q more than the measured genes are?
                g = gene_ranking(importance(cohort, "gex"))
                pos = gene_positions()
                arm_of = lambda gid: ("1p" if (pos.at[gid, "chrom"] == "chr1" and pos.at[gid, "start"] < CEN["chr1"])
                                      else "19q" if (pos.at[gid, "chrom"] == "chr19" and pos.at[gid, "start"] >= CEN["chr19"])
                                      else "") if gid in pos.index else ""
                U = sorted(cohort_genes(cohort, "gex"))
                in_arm = np.array([arm_of(x) != "" for x in U])
                tg = [x for x in g["gene_id"] if x in set(U)][:top]
                k = sum(arm_of(x) != "" for x in tg)
                rows.append({"cohort": cohort, "gene": "1p/19q dosage", "higher_in": a, "vs": b, "source": mk["source"],
                             "n_a": int(ya.sum()), "n_b": int(yb.sum()),
                             "auroc_expected": float(np.mean([r["auroc_a_over_b"] for r in dos_rows
                                                              if r["cohort"] == cohort and r["arm"] in ("1p", "19q")])),
                             "p_mwu": float(max(r["p_mwu"] for r in dos_rows
                                                if r["cohort"] == cohort and r["arm"] in ("1p", "19q"))),
                             "gex_top_on_1p19q": k, "gex_top_n": len(tg),
                             "gex_top_expected_on_1p19q": len(tg) * in_arm.mean(),
                             "p_hypergeom_1p19q": float(hypergeom(len(U), int(in_arm.sum()), len(tg)).sf(k - 1))})
                continue
            G = gex_by_symbol(SCFG, tab["file_id"], [mk["gene"]])
            row = {"cohort": cohort, "gene": mk["gene"], "higher_in": a, "vs": b, "source": mk["source"],
                   "n_a": int(ya.sum()), "n_b": int(yb.sum())}
            if mk["gene"] in G:
                x = G[mk["gene"]].to_numpy()
                u = mannwhitneyu(x[ya], x[yb])
                row.update({"auroc_expected": u.statistic / (ya.sum() * yb.sum()), "p_mwu": float(u.pvalue),
                            "median_log2cpm_a": float(np.median(x[ya])), "median_log2cpm_b": float(np.median(x[yb]))})
            imp_all = importance(cohort, "gex")
            row["gex_rank_all"] = rank_of(imp_all, mk["gene"])
            row["gex_rank_class_a"] = rank_of(importance(cohort, "gex", ia), mk["gene"])
            row["gex_rank_class_b"] = rank_of(importance(cohort, "gex", ib), mk["gene"])
            row["gex_in_top"] = bool(np.nanmin([row["gex_rank_all"], row["gex_rank_class_a"], row["gex_rank_class_b"],
                                                np.inf]) <= top)
            best = (np.inf, "")
            for v in AS_VIEWS:
                if done(cohort, v):
                    r = rank_of(importance(cohort, v), mk["gene"])
                    if np.isfinite(r) and r < best[0]:
                        best = (r, VIEW_NAMES[v])
            row["as_best_rank"], row["as_best_view"] = (best[0] if np.isfinite(best[0]) else np.nan), best[1]
            row["as_in_top"] = bool(best[0] <= top)
            rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "literature_markers.csv", index=False)
    pd.DataFrame(dos_rows).to_csv(OUT / "lgg_arm_dosage.csv", index=False)
    log(f"markers: {len(rows)}")


# ============================================================ 4. event audit
def readthrough_genes() -> set[str]:
    """GENCODE v45 genes with a transcript tagged readthrough_transcript (cached)."""
    p = CACHE / "readthrough_genes.txt"
    if p.exists():
        return set(p.read_text().split())
    rt = set()
    with gzip.open(SCFG.gencode_gtf, "rt") as fh:
        for line in fh:
            if "\ttranscript\t" in line and 'tag "readthrough_transcript"' in line:
                rt.add(re.search(r'gene_id "(ENSG\d+)', line).group(1))
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(sorted(rt)) + "\n")
    return rt


def segments(view: str, feature_id: str, region: str) -> list[tuple[str, int, int]]:
    """(role, start, end), 1-based inclusive: the sequence the event is about (D-41, as D-27).
    SE/AFE/ALE/HIT: the exon; MXE: both exons; RI: the intron; A3SS/A5SS: the extension."""
    if view in ("afe", "ale", "hit"):
        m = re.search(r"chr[\dXYM]+:(\d+)-(\d+)", feature_id)
        return [("exon", int(m.group(1)), int(m.group(2)))]
    head, *rest = feature_id.split(";")
    m = re.search(r"chr[\dXYM]+:(\d+)-(\d+)", head)
    parts = [(int(m.group(1)) + 1, int(m.group(2)))] + [(int(a) + 1, int(b)) for a, b in (p.split("-") for p in rest)]
    if view == "se":
        return [("exon", *parts[0])]
    if view == "mxe":
        return [("exon 1", *parts[0]), ("exon 2", *parts[1])]
    if view == "ri":
        return [("intron", parts[1][1] + 1, parts[2][0] - 1)]
    long_, short = parts[0], parts[1]
    return [("extension", short[1] + 1, long_[1])] if long_[0] == short[0] else [("extension", long_[0], short[0] - 1)]


class Exons:
    def __init__(self):
        ex = pd.read_parquet(SCFG.base_data_dir / "cache" / "gencode_v45_exons.parquet",
                             columns=["chrom", "start", "end", "strand", "gene_id", "gene_name", "gene_type"])
        ex["gene_id"] = ex["gene_id"].str.split(".").str[0]
        ex = ex.drop_duplicates(["chrom", "start", "end", "gene_id"])
        self.rt = readthrough_genes()
        self.by = {c: g.sort_values("start").reset_index(drop=True) for c, g in ex.groupby("chrom")}
        self.strand = ex.drop_duplicates("gene_id").set_index("gene_id")["strand"]
        self.maxlen = {c: int((g["end"] - g["start"]).max()) + 1 for c, g in self.by.items()}

    def overlaps(self, chrom: str, a: int, b: int, host: str, min_bp: int) -> pd.DataFrame:
        g = self.by.get(chrom)
        if g is None:
            return pd.DataFrame()
        st = g["start"].to_numpy()
        lo, hi = np.searchsorted(st, a - self.maxlen[chrom], "left"), np.searchsorted(st, b, "right")
        c = g.iloc[lo:hi]
        bp = np.minimum(c["end"].to_numpy(), b) - np.maximum(c["start"].to_numpy(), a) + 1
        c = c.assign(overlap_bp=bp)
        return c[(c["overlap_bp"] >= min_bp) & (c["gene_id"] != host)]


def rho(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 20 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return np.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def marker_events(cohort: str, cls: list[str]) -> list[pd.DataFrame]:
    """The splicing event a literature marker is represented by, so that it is audited before it is
    named (D-41): in the AS view where the marker gene ranks best (step 3), the gene's feature with
    the largest importance over the patients of the class the literature expects it higher in."""
    f = OUT / "literature_markers.csv"
    if not f.exists():
        return []
    mk = pd.read_csv(f)
    mk = mk[(mk.cohort == cohort) & mk["as_best_view"].notna() & mk["higher_in"].isin(cls)]
    out = []
    inv = {n: v for v, n in VIEW_NAMES.items()}
    for r in mk.itertuples(index=False):
        view, k = inv[r.as_best_view], cls.index(r.higher_in)
        t = importance(cohort, view, k).reset_index(drop=True)
        hit = np.flatnonzero(t["symbol"].to_numpy() == r.gene)
        if len(hit):
            out.append(t.iloc[[hit[0]]].assign(view=view, cls=r.higher_in, rank_in_class=int(hit[0]) + 1,
                                              selected_by="literature marker"))
    return out


def step_event_audit() -> None:
    n_top, min_bp, r_min = int(ITP["top_events_per_class"]), int(ITP["overlap_min_bp"]), float(ITP["overlap_tracking_rho"])
    EX = Exons()
    rows = []
    for cohort in TABLES:
        cls = classes_of(cohort)
        tab = TABLES[cohort]
        y = tab["y"].to_numpy()
        cand = []
        for view in AS_VIEWS:
            if not done(cohort, view):
                continue
            for k, c in enumerate(cls):
                t = importance(cohort, view, k).head(n_top)
                cand.append(t.assign(view=view, cls=c, rank_in_class=np.arange(1, len(t) + 1),
                                     selected_by=f"top {n_top} of the class"))
        cand += marker_events(cohort, cls)
        if not cand:
            continue
        cand = pd.concat(cand, ignore_index=True)
        G = None
        for view, g in cand.groupby("view", sort=False):
            X = feature_values(SCFG, view, tab["file_id"], list(g["feature_id"]))
            hosts = list(g["symbol"].dropna().unique())
            neigh = {}
            for r in g.itertuples(index=False):
                host = str(r.gene_id).split(".")[0]
                hs = r.strand if r.strand in ("+", "-") else EX.strand.get(host)
                best = None
                for role, a, b in segments(view, r.feature_id, r.region):
                    for x in EX.overlaps(r.chrom, a, b, host, min_bp).sort_values("overlap_bp", ascending=False) \
                            .drop_duplicates("gene_id").itertuples(index=False):
                        kind = ("antisense" if x.strand != hs else
                                "readthrough" if (x.gene_id in EX.rt or host in EX.rt) else "same-strand")
                        rk = {"antisense": 3, "same-strand": 2, "readthrough": 1}[kind]
                        if best is None or (rk, x.overlap_bp) > best[0]:
                            best = ((rk, x.overlap_bp), x.gene_name, kind, int(x.overlap_bp), role)
                neigh[r.feature_id] = best
            genes = sorted(set(hosts) | {b[1] for b in neigh.values() if b})
            G = gex_by_symbol(SCFG, tab["file_id"], genes)
            for j, r in enumerate(g.itertuples(index=False)):
                x = X[:, list(g["feature_id"]).index(r.feature_id)].astype(float)
                best = neigh[r.feature_id]
                other = best[1] if best else ""
                rho_o = rho(x, G[other].to_numpy()) if other in G else np.nan
                rho_h = rho(x, G[r.symbol].to_numpy()) if r.symbol in G else np.nan
                kind = best[2] if best else "none"
                tracks = bool(kind in ("antisense", "same-strand") and np.isfinite(rho_o) and abs(rho_o) >= r_min
                              and abs(rho_o) > (abs(rho_h) if np.isfinite(rho_h) else 0))
                med = {f"median_{c}": float(np.nanmedian(x[y == i])) if np.isfinite(x[y == i]).any() else np.nan
                       for i, c in enumerate(cls)}
                spread = {}
                for i, c in enumerate(cls):
                    xi = x[(y == i) & np.isfinite(x)]
                    spread[f"observed_{c}"] = float(np.isfinite(x[y == i]).mean())
                    spread[f"q10_{c}"] = float(np.quantile(xi, 0.1)) if len(xi) >= MIN_CELL else np.nan
                    spread[f"q90_{c}"] = float(np.quantile(xi, 0.9)) if len(xi) >= MIN_CELL else np.nan
                mk = cls.index(r.cls)
                others = [med[f"median_{c}"] for c in cls if c != r.cls]
                rows.append({"cohort": cohort, "view": VIEW_NAMES[view], "class": r.cls, "selected_by": r.selected_by,
                             "rank_in_class": int(r.rank_in_class), "feature_id": r.feature_id, "symbol": r.symbol,
                             "region": r.region, "mean_abs_shap_class": r.mean_abs_shap,
                             "direction_spearman": r.direction_spearman, "n_folds_selected": r.n_folds_selected,
                             "observed_frac": float(np.isfinite(x).mean()),
                             "overlap_class": kind, "other_gene": other, "overlap_bp": best[3] if best else 0,
                             "overlap_segment": best[4] if best else "", "rho_other_gene": rho_o, "rho_host": rho_h,
                             "tracks_other_gene": tracks, **med, **spread,
                             "gap_to_nearest_class": float(med[f"median_{r.cls}"] - min(others, key=lambda v: abs(
                                 v - med[f"median_{r.cls}"]) if np.isfinite(v) else np.inf)) if len(others) else np.nan,
                             "class_index": mk})
        log(f"event audit: {cohort} {len(cand)} events")
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "event_audit.csv", index=False)
    s = out.drop_duplicates(["cohort", "view", "feature_id"])
    log("overlap classes among audited events:\n" + s.groupby(["view", "overlap_class"]).size().unstack(fill_value=0)
        .to_string() + f"\ntracking another gene: {int(s.tracks_other_gene.sum())} of {len(s)}")


# ============================================================ 5. cross-family errors
PAIRS = [("LGG", "LGG_IDH"), ("UCEC", "UCEC_TCGA"), ("STAD", "STAD_TCGA"), ("BRCA", "BRCA_PAM50"),
         ("BRCA", "BRCA_ER"), ("PCPG", "PCPG_mRNA")]


def step_cross_family() -> None:
    rows = []
    for hist, mol in PAIRS:
        if hist not in TABLES or mol not in TABLES or not all(done(hist, v) for v in SCFG.views):
            continue
        cls = classes_of(hist)
        P = {v: pd.read_parquet(SCFG.run_dir(TASK, v, hist) / "oof_predictions.parquet").set_index("case_id")
             for v in SCFG.views}
        cols = [f"p_{c}" for c in cls]
        pred = {"GEX": P["gex"][cols].to_numpy().argmax(1),
                "AS ensemble": np.mean([P[v].loc[P["gex"].index, cols].to_numpy() for v in AS_VIEWS], axis=0).argmax(1)}
        lab = TABLES[mol].set_index("case_id")["label"]
        both = P["gex"].index.intersection(lab.index)
        true = P["gex"].loc[both, "label"]
        for model, p in pred.items():
            pp = pd.Series(np.array(cls)[p], index=P["gex"].index).loc[both]
            ct = pd.DataFrame({"true": true, "predicted": pp, "molecular": lab[both]}).value_counts().rename("patients")
            for (t_, p_, m_), n in ct.items():
                rows.append({"histology_cohort": hist, "molecular_cohort": mol, "model": model, "true_histology": t_,
                             "predicted_histology": p_, "molecular_class": m_, "patients": int(n)})
    pd.DataFrame(rows).to_csv(OUT / "cross_family_errors.csv", index=False)
    log(f"cross_family: {len(rows)} cells")


def step_error_markers() -> None:
    """For each cohort with literature markers: the median log2 CPM of each marker gene among the
    patients of every (true class, out-of-fold call) cell, for the GEX model and the AS ensemble.
    Says what the misclassified patients look like on the genes the literature names (descriptive).
    A cell of fewer than MIN_CELL patients gets no median: it would be close to one patient's value
    (D-12)."""
    rows = []
    for cohort, markers in SCFG.spec["literature_markers"].items():
        if cohort not in TABLES or not all(done(cohort, v) for v in SCFG.views):
            continue
        tab, cls = TABLES[cohort], classes_of(cohort)
        genes = [m["gene"] for m in markers if m["gene"] != "1p/19q dosage"]
        if not genes:
            continue
        G = gex_by_symbol(SCFG, tab["file_id"], genes)
        cols = [f"p_{c}" for c in cls]
        P = {v: pd.read_parquet(SCFG.run_dir(TASK, v, cohort) / "oof_predictions.parquet").set_index("file_id")
             .loc[tab["file_id"], cols].to_numpy() for v in SCFG.views}
        calls = {"GEX": P["gex"].argmax(1), "AS ensemble": np.mean([P[v] for v in AS_VIEWS], axis=0).argmax(1)}
        y = tab["y"].to_numpy()
        for model, c in calls.items():
            for t in range(len(cls)):
                for k in range(len(cls)):
                    m = (y == t) & (c == k)
                    for g in genes:
                        if g not in G:
                            continue
                        rows.append({"cohort": cohort, "model": model, "true_class": cls[t], "called_class": cls[k],
                                     "gene": g, "patients": int(m.sum()),
                                     "median_log2cpm": float(np.median(G[g].to_numpy()[m]))
                                     if m.sum() >= MIN_CELL else np.nan})
    pd.DataFrame(rows).to_csv(OUT / "error_markers.csv", index=False)
    log(f"error_markers: {len(rows)} cells")


# classes defined the same way in more than one cohort; "Classical" (THCA vs HNSC) and "Basal"
# (BRCA PAM50 vs HNSC) share a name only
SHARED_CLASSES = ("MSI", "CIN", "GS", "Squamous", "Adenocarcinoma")


def step_shared_classes() -> None:
    """Classes defined the same way in more than one cohort (MSI in three cancers, CIN and GS in two,
    squamous and adenocarcinoma in two): every audited top event of such a class, measured in every
    cohort that has the class, as that class against the rest of its cohort (AUROC; no model). Does
    an event a model chose in one cancer mark the same class in another? Exploratory (A-4); run
    after event_audit."""
    ev = pd.read_csv(OUT / "event_audit.csv")
    ev = ev[ev["selected_by"].str.startswith("top")]
    inv = {n: v for v, n in VIEW_NAMES.items()}
    holders = {}
    for c in TABLES:
        for cl in classes_of(c):
            if cl in SHARED_CLASSES:
                holders.setdefault(cl, []).append(c)
    rows = []
    for cl, cohorts in holders.items():
        cand = ev[(ev["class"] == cl) & ev["cohort"].isin(cohorts)].drop_duplicates(["view", "feature_id"])
        for view, g in cand.groupby("view"):
            for c in cohorts:
                tab = TABLES[c]
                ycl = tab["label"].to_numpy() == cl
                X = feature_values(SCFG, inv[view], tab["file_id"], list(g["feature_id"]))
                for j, r in enumerate(g.itertuples(index=False)):
                    x = X[:, j].astype(float)
                    a, b = x[ycl & np.isfinite(x)], x[~ycl & np.isfinite(x)]
                    ok = len(a) >= MIN_CELL and len(b) >= MIN_CELL
                    u = mannwhitneyu(a, b) if ok else None
                    rows.append({"class": cl, "view": view, "feature_id": r.feature_id, "symbol": r.symbol,
                                 "region": r.region, "chosen_in": r.cohort, "measured_in": c,
                                 "same_cohort": r.cohort == c, "n_class": int(len(a)), "n_rest": int(len(b)),
                                 "median_class": float(np.median(a)) if ok else np.nan,
                                 "median_rest": float(np.median(b)) if ok else np.nan,
                                 "auroc_class_over_rest": u.statistic / (len(a) * len(b)) if ok else np.nan,
                                 "p_mwu": float(u.pvalue) if ok else np.nan,
                                 "overlap_class": r.overlap_class, "tracks_other_gene": r.tracks_other_gene})
    out = pd.DataFrame(rows)
    rep = out[~out["same_cohort"]].copy()
    out["q_mwu_replication"] = np.nan
    if len(rep):
        out.loc[rep.index, "q_mwu_replication"] = bh(rep["p_mwu"].fillna(1).to_numpy())
    out["q_family"] = np.where(out["same_cohort"], "", "BH across every (event, other cohort) replication test")
    out.to_csv(OUT / "shared_class_events.csv", index=False)
    log(f"shared_classes: {len(out)} rows, {len(rep)} replication tests")


MSI_CONFOUND_GENES = ("PTPRC", "CD8A", "CIITA", "MLH1")


def step_msi_confounds() -> None:
    """MSI tumours carry dense immune infiltrate, and most sporadic ones have silenced MLH1. For each
    audited top event of the MSI class, in each cohort with an MSI class: Spearman with the log2 CPM of
    PTPRC (all leukocytes), CD8A (cytotoxic T cells), CIITA (interferon-gamma / MHC II) and MLH1, and the
    MSI-vs-rest AUROC within the patients at or below the cohort's median PTPRC. If an event read
    immune content, its correlation with PTPRC would be strong and its AUROC would fade in the
    low-immune half. Exploratory (A-4); run after event_audit. Aggregates only."""
    ev = pd.read_csv(OUT / "event_audit.csv")
    ev = ev[ev["selected_by"].str.startswith("top") & (ev["class"] == "MSI")]
    ev = ev.drop_duplicates(["view", "feature_id"])
    inv = {n: v for v, n in VIEW_NAMES.items()}
    rows = []
    for c in [c for c in TABLES if "MSI" in classes_of(c)]:
        tab = TABLES[c]
        y = tab["label"].to_numpy() == "MSI"
        G = gex_by_symbol(SCFG, tab["file_id"], list(MSI_CONFOUND_GENES))
        low = G["PTPRC"].to_numpy() <= np.median(G["PTPRC"].to_numpy())
        for s in G.columns:
            g = G[s].to_numpy()
            rows.append({"cohort": c, "kind": "gene", "symbol": s, "view": "GEX", "feature_id": s, "region": "",
                         "n_msi": int(y.sum()), "n_rest": int((~y).sum()),
                         "auroc_msi_over_rest": mannwhitneyu(g[y], g[~y]).statistic / (y.sum() * (~y).sum())})
        for view, g in ev.groupby("view"):
            X = feature_values(SCFG, inv[view], tab["file_id"], list(g["feature_id"]))
            for j, r in enumerate(g.itertuples(index=False)):
                x = X[:, j].astype(float)
                ok = np.isfinite(x)
                d = {"cohort": c, "kind": "event", "symbol": r.symbol, "view": view, "feature_id": r.feature_id,
                     "region": r.region, "chosen_in": r.cohort, "n_msi": int((ok & y).sum()),
                     "n_rest": int((ok & ~y).sum())}
                a, b = x[ok & y], x[ok & ~y]
                if min(len(a), len(b)) >= MIN_CELL:
                    d["auroc_msi_over_rest"] = mannwhitneyu(a, b).statistic / (len(a) * len(b))
                for s in G.columns:
                    d[f"rho_{s}"] = spearmanr(x[ok], G[s].to_numpy()[ok]).statistic
                m = ok & low
                d["n_msi_low_ptprc"], d["n_rest_low_ptprc"] = int((m & y).sum()), int((m & ~y).sum())
                a, b = x[m & y], x[m & ~y]
                if min(len(a), len(b)) >= MIN_CELL:
                    d["auroc_msi_over_rest_low_ptprc"] = mannwhitneyu(a, b).statistic / (len(a) * len(b))
                rows.append(d)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "msi_event_confounds.csv", index=False)
    log(f"msi_confounds: {len(out)} rows")


STEPS = {"importance": step_importance, "as_vs_gex": step_as_vs_gex, "markers": step_markers,
         "event_audit": step_event_audit, "cross_family": step_cross_family, "error_markers": step_error_markers,
         "shared_classes": step_shared_classes, "msi_confounds": step_msi_confounds}

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for name in sys.argv[1:] or list(STEPS):
        STEPS[name]()
