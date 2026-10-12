"""Step 23: do the models' events read another gene? (D-27)

An event is a statement about one gene's RNA. Two things can make it a statement about a different
gene instead, and both are expected in this data:

- TCGA mRNA libraries are mostly unstranded, so where an exon of the host gene overlaps an exon of
  a gene on the opposite strand, reads of either can be counted towards the event.
- GENCODE readthrough and nested gene models put exons of one gene inside another gene's model, so
  an "alternative" first, last or mutually exclusive exon of the host can be an exon of its
  neighbour, and its usage then measures how much of the neighbour is transcribed.

This step asks how often that happens, and whether it matters.

1. Annotation, every measured event of the 8 splicing views. The segment an event is about (the
   alternative exon; for MXE both exons of the pair; the retained intron; the splice-site
   extension) is intersected with the exons of every other GENCODE v45 gene. An overlap of
   >= 50 bp counts. Each event is classed as antisense (another gene's exon on the opposite
   strand), same-strand (another gene on the same strand), readthrough (the same-strand overlap
   is with a readthrough model of the host, or the host is itself a readthrough model and the
   exon belongs to one of its component genes: a naming question, not a second gene's
   expression), or none, in that order of precedence. "Overlapping" below means antisense or
   same-strand.
2. Enrichment. The share of overlapping events among the models' top events against all measured
   events of the view: tumour status, the top 100 per cancer and view; cancer type, the top 30
   per class and view. Fisher's exact test per view, BH over views within a task. Long segments
   overlap something more often by chance, so the same comparison is repeated within deciles of
   queried length (Mantel-Haenszel); this was added after the first run, as a robustness check.
3. Expression tracking, for every top event with an overlap. Spearman between the event's value
   and the overlapping gene's log2 CPM, and with the host gene's, over all samples of the task
   (within each cancer for tumour status, the median taken over the cancers whose top 100 holds
   the event; across all 9,729 samples for cancer type). An event *tracks* its neighbour when
   |rho_neighbour| >= 0.5 and |rho_neighbour| > |rho_host| — fixed before any result was seen.
   Whether the coupling is read misassignment or genuine co-transcription, the event is then a
   readout of the neighbour's expression and is not to be read as a splicing choice of the host.
4. Robustness. The survival enrichment of D-25 is recomputed with every overlapping event removed
   from both the candidate and the background set.

Writes results/biology/12_overlap_audit/. No model is refitted; nothing here feeds back into a
model or into the choice of an event set.

Usage: python scripts/23_overlap_audit.py
"""
from __future__ import annotations

import gzip
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
from _common import REPO, log
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.contingency_tables import StratifiedTable

from ctsm.biology import bh, feature_values, gex_by_symbol, task_table
from ctsm.config import PAN, load_config
from ctsm.matrix_io import view_dir

AS_VIEWS = ["afe", "ale", "hit", "se", "a5ss", "a3ss", "mxe", "ri"]
MIN_OVERLAP_BP = 50
TRACK_RHO = 0.5
OUT = REPO / "results" / "biology" / "12_overlap_audit"

_spec = importlib.util.spec_from_file_location("cards", Path(__file__).with_name("22_event_cards.py"))
cards = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cards)


# ------------------------------------------------------------------ GENCODE exons
def exon_table(cfg) -> pd.DataFrame:
    """Every GENCODE v45 exon with its gene, cached locally; `readthrough` marks genes whose
    models carry the readthrough_transcript tag."""
    cache = cfg.data_dir / "cache" / "gencode_v45_exons_tagged.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    gtf = cfg.gencode_gtf
    rows, rt = [], set()
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split("\t", 8)
            if t[2] not in ("exon", "transcript"):
                continue
            gid = re.search(r'gene_id "([^"]+)"', t[8]).group(1).split(".")[0]
            if t[2] == "transcript":
                if 'tag "readthrough_transcript"' in t[8]:
                    rt.add(gid)
                continue
            rows.append((t[0], int(t[3]), int(t[4]), t[6], gid,
                         re.search(r'gene_name "([^"]+)"', t[8]).group(1),
                         re.search(r'gene_type "([^"]+)"', t[8]).group(1)))
    ex = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene_id", "gene_name", "gene_type"])
    ex = ex.drop_duplicates(["chrom", "start", "end", "gene_id"])
    ex["readthrough"] = ex["gene_id"].isin(rt)
    cache.parent.mkdir(parents=True, exist_ok=True)
    ex.to_parquet(cache, index=False)
    return ex


class ExonIndex:
    """Per-chromosome sorted exon starts, for fast interval queries. Exons longer than LONG_EXON
    are kept apart and checked exhaustively, so one giant exon does not widen every window."""

    LONG_EXON = 20_000

    def __init__(self, ex: pd.DataFrame):
        self.readthrough = set(ex.loc[ex["readthrough"], "gene_id"])
        self.by_chrom = {}
        for c, g in ex.groupby("chrom"):
            g = g.sort_values("start").reset_index(drop=True)
            long = (g["end"] - g["start"]).to_numpy() > self.LONG_EXON
            short = np.flatnonzero(~long)
            st, en = g["start"].to_numpy(), g["end"].to_numpy()
            span = int((en[short] - st[short]).max()) + 1 if len(short) else 1
            self.by_chrom[c] = (g, st, en, g["gene_id"].to_numpy(), short, np.flatnonzero(long),
                                st[short], span)

    def query(self, chrom: str, a: int, b: int, exclude_gene: str = "", min_bp: int = 1) -> pd.DataFrame:
        """Exons overlapping [a, b] by at least `min_bp`, other than those of `exclude_gene`."""
        if chrom not in self.by_chrom:
            return pd.DataFrame()
        g, st, en, gid, short, long, s_st, span = self.by_chrom[chrom]
        lo, hi = np.searchsorted(s_st, a - span, "left"), np.searchsorted(s_st, b, "right")
        cand = np.concatenate([short[lo:hi], long])
        bp = np.minimum(en[cand], b) - np.maximum(st[cand], a) + 1
        keep = (bp >= min_bp) & (gid[cand] != exclude_gene)
        if not keep.any():
            return pd.DataFrame()
        return g.iloc[cand[keep]].assign(overlap_bp=bp[keep])


# ------------------------------------------------------------------ 1. annotation
def segments(view: str, feature_id: str, region) -> list[tuple[str, int, int]]:
    """(role, start, end) of the sequence an event is about, 1-based inclusive. For MXE both
    exons count: a swap between exons of two different genes is a readthrough artefact whichever
    of the two the PSI follows."""
    parts = cards.pieces(view, feature_id)
    if view == "mxe":
        tracked, other = cards.tracked_and_partner(parts, region)
        return [("tracked exon", *tracked), ("partner exon", *other)]
    a, b, kind = cards.alternative_segment(view, parts)
    return [(kind, a, b)]


def annotate_view(cfg, idx: ExonIndex, view: str) -> pd.DataFrame:
    F = pd.read_parquet(view_dir(cfg.matrix_dir, view) / "features.parquet")
    host_strand = None
    rows = []
    for r in F.itertuples(index=False):
        chrom = r.chrom
        host = str(r.gene_id).split(".")[0]
        host_rt = host in idx.readthrough
        best = None
        cls = "none"
        segs = segments(view, r.feature_id, r.region)
        query_nt = sum(b - a + 1 for _, a, b in segs)
        for role, a, b in segs:
            hit = idx.query(chrom, a, b, exclude_gene=host, min_bp=MIN_OVERLAP_BP)
            if hit.empty:
                continue
            g = hit.sort_values("overlap_bp", ascending=False).drop_duplicates("gene_id")
            hs = r.strand if r.strand in ("+", "-") else None
            if hs is None:
                hs = host_strand_of(idx, chrom, host)
            for x in g.itertuples(index=False):
                kind = ("antisense" if x.strand != hs else
                        "readthrough" if (x.readthrough or host_rt) else "same-strand")
                rank = {"antisense": 3, "same-strand": 2, "readthrough": 1}[kind]
                if best is None or (rank, x.overlap_bp) > best[0]:
                    best = ((rank, x.overlap_bp), x.gene_name, x.gene_id, x.gene_type, kind, int(x.overlap_bp),
                            role, b - a + 1)
        if best is not None:
            _, name, gid, gtype, cls, bp, role, seg_nt = best
            rows.append({"view": view, "feature_id": r.feature_id, "symbol": r.symbol, "overlap_class": cls,
                         "other_gene": name, "other_gene_id": gid, "other_gene_type": gtype,
                         "overlap_bp": bp, "segment_role": role, "segment_nt": seg_nt, "query_nt": query_nt})
        else:
            rows.append({"view": view, "feature_id": r.feature_id, "symbol": r.symbol, "overlap_class": cls,
                         "query_nt": query_nt})
    log(f"{view}: {len(rows):,} events annotated")
    return pd.DataFrame(rows)


_STRAND_CACHE: dict[str, str] = {}


def host_strand_of(idx: ExonIndex, chrom: str, gene_id: str) -> str | None:
    if gene_id not in _STRAND_CACHE:
        g = idx.by_chrom.get(chrom, (pd.DataFrame(columns=["gene_id", "strand"]),))[0]
        s = g.loc[g["gene_id"] == gene_id, "strand"]
        _STRAND_CACHE[gene_id] = s.iloc[0] if len(s) else None
    return _STRAND_CACHE[gene_id]


# ------------------------------------------------------------------ 2. enrichment
def top_sets(cfg) -> dict[str, pd.DataFrame]:
    bio = REPO / "results" / "biology"
    ts = pd.read_csv(bio / "04_events" / "top_events_per_cancer.csv", usecols=["view", "feature_id", "cancer"])
    ct = []
    for v in AS_VIEWS:
        p = cfg.results_dir / "cancer_type" / v / PAN / "shap_top_features_by_class.csv"
        ct.append(pd.read_csv(p, usecols=["feature_id", "class", "rank_in_class"]).assign(view=v))
    return {"tumour status": ts, "cancer type": pd.concat(ct, ignore_index=True)}


def length_matched_or(top_hit, top_bin, rest_hit, rest_bin) -> tuple[float, float]:
    """Mantel-Haenszel odds ratio over strata of queried length (deciles within the view): a long
    terminal exon overlaps something more often by chance alone."""
    tables = []
    for b in np.union1d(top_bin, rest_bin):
        t, r = top_hit[top_bin == b], rest_hit[rest_bin == b]
        tab = np.array([[t.sum(), len(t) - t.sum()], [r.sum(), len(r) - r.sum()]], dtype=float)
        if len(t) and len(r):
            tables.append(tab)
    if not tables:
        return np.nan, np.nan
    st = StratifiedTable(tables)
    return float(st.oddsratio_pooled), float(st.test_null_odds(correction=True).pvalue)


def enrichment(ann: pd.DataFrame, sets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for task, sel in sets.items():
        for v in AS_VIEWS:
            a = ann[ann["view"] == v].copy()
            a["len_bin"] = pd.qcut(a["query_nt"].rank(method="first"), 10, labels=False)
            chosen = set(sel.loc[sel["view"] == v, "feature_id"])
            top = a[a["feature_id"].isin(chosen)]
            rest = a[~a["feature_id"].isin(chosen)]
            for cls in ("any", "antisense", "same-strand"):
                hit = (lambda d: d["overlap_class"].isin(["antisense", "same-strand"])) if cls == "any" else \
                      (lambda d, c=cls: d["overlap_class"].eq(c))
                k1, n1, k0, n0 = int(hit(top).sum()), len(top), int(hit(rest).sum()), len(rest)
                odds, p = fisher_exact([[k1, n1 - k1], [k0, n0 - k0]]) if n1 and n0 else (np.nan, np.nan)
                lm_or, lm_p = length_matched_or(hit(top).to_numpy(), top["len_bin"].to_numpy(),
                                                hit(rest).to_numpy(), rest["len_bin"].to_numpy())
                rows.append({"task": task, "view": v, "overlap": cls, "top_events": n1, "top_overlapping": k1,
                             "top_share": k1 / n1 if n1 else np.nan, "all_events": n0, "all_overlapping": k0,
                             "all_share": k0 / n0 if n0 else np.nan, "odds_ratio": odds, "p_fisher": p,
                             "top_median_nt": float(top["query_nt"].median()),
                             "all_median_nt": float(rest["query_nt"].median()),
                             "odds_ratio_length_matched": lm_or, "p_length_matched": lm_p})
    out = pd.DataFrame(rows)
    for col in ("p_fisher", "p_length_matched"):
        q = col.replace("p_", "q_")
        out[q] = np.nan
        for (task, cls), g in out.groupby(["task", "overlap"]):
            m = g[col].notna()
            out.loc[g.index[m], q] = bh(g.loc[m, col].to_numpy())
    out["q_family"] = "BH over the 8 views, within task and overlap class"
    return out


# ------------------------------------------------------------------ 3. expression tracking
def rho(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 20 or np.nanstd(x[ok]) == 0 or np.nanstd(y[ok]) == 0:
        return np.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def tracking_status(cfg, flagged: pd.DataFrame, sel: pd.DataFrame) -> pd.DataFrame:
    """Within each cancer whose top 100 holds the event: rho with the neighbour and with the host."""
    pairs = sel.merge(flagged, on=["view", "feature_id"])
    rows = []
    for c, g in pairs.groupby("cancer"):
        tab = task_table(cfg, "tumor_status", c)
        genes = sorted(set(g["other_gene"]) | set(g["symbol"].dropna()))
        G = gex_by_symbol(cfg, tab["file_id"], genes)
        for v, gv in g.groupby("view"):
            X = feature_values(cfg, v, tab["file_id"], list(gv["feature_id"]))
            for j, r in enumerate(gv.itertuples(index=False)):
                x = X[:, j].astype(float)
                rows.append({"view": v, "feature_id": r.feature_id, "cancer": c, "symbol": r.symbol,
                             "other_gene": r.other_gene, "overlap_class": r.overlap_class,
                             "rho_neighbour": rho(x, G[r.other_gene].to_numpy()) if r.other_gene in G else np.nan,
                             "rho_host": rho(x, G[r.symbol].to_numpy()) if r.symbol in G else np.nan,
                             "neighbour_in_gex": r.other_gene in G})
    per = pd.DataFrame(rows)
    agg = (per.groupby(["view", "feature_id", "symbol", "other_gene", "overlap_class"], dropna=False)
           .agg(n_cancers=("cancer", "nunique"), rho_neighbour=("rho_neighbour", "median"),
                rho_host=("rho_host", "median"), neighbour_in_gex=("neighbour_in_gex", "all")).reset_index())
    return per, agg


def tracking_type(cfg, flagged: pd.DataFrame) -> pd.DataFrame:
    """Across all cancer-type samples, and over the 34 class medians."""
    tab = task_table(cfg, "cancer_type", PAN)
    genes = sorted(set(flagged["other_gene"]) | set(flagged["symbol"].dropna()))
    G = gex_by_symbol(cfg, tab["file_id"], genes)
    cls = tab["label"].to_numpy()
    rows = []
    for v, gv in flagged.groupby("view"):
        X = feature_values(cfg, v, tab["file_id"], list(gv["feature_id"]))
        for j, r in enumerate(gv.itertuples(index=False)):
            x = X[:, j].astype(float)
            o = G[r.other_gene].to_numpy() if r.other_gene in G else None
            med = (pd.DataFrame({"x": x, "o": o, "c": cls}).groupby("c").median() if o is not None else None)
            rows.append({"view": v, "feature_id": r.feature_id, "symbol": r.symbol, "other_gene": r.other_gene,
                         "overlap_class": r.overlap_class,
                         "rho_neighbour": rho(x, o) if o is not None else np.nan,
                         "rho_host": rho(x, G[r.symbol].to_numpy()) if r.symbol in G else np.nan,
                         "rho_neighbour_class_medians": (rho(med["x"].to_numpy(), med["o"].to_numpy())
                                                         if med is not None else np.nan),
                         "neighbour_top_class": med["o"].idxmax() if med is not None else None,
                         "neighbour_in_gex": o is not None})
    return pd.DataFrame(rows)


def tracks(df: pd.DataFrame) -> pd.Series:
    return (df["rho_neighbour"].abs() >= TRACK_RHO) & (df["rho_neighbour"].abs() > df["rho_host"].abs().fillna(0))


# ------------------------------------------------------------------ 4. survival, without them
def survival_without(ann: pd.DataFrame) -> pd.DataFrame:
    d = REPO / "results" / "biology" / "08_survival"
    idx = pd.read_csv(d / "event_index.csv")
    idx["feature_id"] = idx["event_id"].str.split(":", n=1).str[1]
    idx = idx.merge(ann[["view", "feature_id", "overlap_class"]], on=["view", "feature_id"], how="left")
    cells = pd.read_csv(d / "event_survival_cells.csv.gz",
                        usecols=["event_key", "endpoint", "arm", "cox_p", "cox_ok"])
    cells = cells[cells["cox_ok"]].merge(idx[["event_key", "role", "overlap_class"]], on="event_key")
    rows = []
    for (ep, arm), g in cells.groupby(["endpoint", "arm"]):
        for scope, keep in (("all events", g.index == g.index),
                            ("overlapping events removed",
                             ~g["overlap_class"].isin(["antisense", "same-strand"]).to_numpy())):
            x = g[keep]
            cand, bg = x[x["role"] == "candidate"], x[x["role"] == "background"]
            k1, n1 = int((cand["cox_p"] < 0.05).sum()), len(cand)
            k0, n0 = int((bg["cox_p"] < 0.05).sum()), len(bg)
            odds, p = fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])
            rows.append({"endpoint": ep, "arm": arm, "scope": scope, "candidate_cells": n1,
                         "candidate_rate": k1 / n1, "background_cells": n0, "background_rate": k0 / n0,
                         "odds_ratio": odds, "p_fisher": p})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ main
def main() -> None:
    cfg = load_config()
    OUT.mkdir(parents=True, exist_ok=True)
    idx = ExonIndex(exon_table(cfg))
    log("exon index built")

    ann = pd.concat([annotate_view(cfg, idx, v) for v in AS_VIEWS], ignore_index=True)
    ann.to_csv(OUT / "event_overlaps_all_measured.csv.gz", index=False)
    log(f"annotated {len(ann):,} events; overlapping another gene: "
        f"{ann['overlap_class'].isin(['antisense', 'same-strand']).mean():.1%}")

    sets = top_sets(cfg)
    enr = enrichment(ann, sets)
    enr.to_csv(OUT / "overlap_enrichment_by_view.csv", index=False)

    flag = ann[ann["overlap_class"].isin(["antisense", "same-strand"])]
    sel_ts = sets["tumour status"].drop_duplicates(["view", "feature_id", "cancer"])
    per, agg = tracking_status(cfg, flag[flag["feature_id"].isin(sel_ts["feature_id"])], sel_ts)
    per.to_csv(OUT / "tracking_tumour_status_by_cancer.csv", index=False)
    agg["tracks_neighbour"] = tracks(agg)
    agg.to_csv(OUT / "tracking_tumour_status.csv", index=False)

    sel_ct = sets["cancer type"].drop_duplicates(["view", "feature_id"])
    tt = tracking_type(cfg, flag[flag["feature_id"].isin(sel_ct["feature_id"])].drop_duplicates(["view", "feature_id"]))
    tt["tracks_neighbour"] = tracks(tt)
    tt.to_csv(OUT / "tracking_cancer_type.csv", index=False)

    surv = survival_without(ann)
    surv.to_csv(OUT / "survival_enrichment_without_overlaps.csv", index=False)

    # the named events of both narratives, and the recurrent set
    rec = pd.read_csv(REPO / "results" / "biology" / "04_events" / "recurrent_events.csv",
                      usecols=["view", "feature_id", "symbol", "n_cancers", "majority_sign"])
    rec = rec.merge(ann.drop(columns="symbol"), on=["view", "feature_id"], how="left").merge(
        agg[["view", "feature_id", "rho_neighbour", "rho_host", "tracks_neighbour"]], on=["view", "feature_id"],
        how="left")
    rec.to_csv(OUT / "recurrent_events_overlaps.csv", index=False)
    named = []
    for task, f, extra in (("tumour status", "tumour_status_cards.csv", agg),
                           ("cancer type", "cancer_type_cards.csv", tt)):
        c = pd.read_csv(REPO / "results" / "biology" / "11_event_cards" / f, usecols=["symbol", "view", "feature_id"])
        c = c.merge(ann.drop(columns="symbol"), on=["view", "feature_id"], how="left").merge(
            extra[["view", "feature_id", "rho_neighbour", "rho_host", "tracks_neighbour"]],
            on=["view", "feature_id"], how="left")
        named.append(c.assign(task=task))
    pd.concat(named, ignore_index=True).to_csv(OUT / "named_events_overlaps.csv", index=False)
    log("written " + str(OUT))


if __name__ == "__main__":
    main()
