"""Helpers for the biological interpretation layer (scripts 07-12, config/biology.yaml).

These analyses describe the fitted models (via their held-out SHAP values) and the data. They
never refit a model, and no result feeds back into a modelling choice.

Contents
  data access      cohort matrices, matched tumour/normal pairs, GEX by gene symbol
  SHAP             full per-feature importance tables (all features, cached)
  enrichment       GOseq-style over-representation that corrects for the number of
                   measured events per gene (Young et al. 2010, Genome Biol 11:R14)
  protein          map a cassette exon to a transcript, CDS, reading frame and domains
"""
from __future__ import annotations

import functools
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import hypergeom, nchypergeom_wallenius
from sklearn.isotonic import IsotonicRegression

from .config import CONFIG_DIR, Config
from .explain import importance_table, load_folds
from .matrix_io import decode, load_meta, load_view, view_dir


def bio_config() -> dict:
    with open(CONFIG_DIR / "biology.yaml") as fh:
        return yaml.safe_load(fh)


def ref_path(cfg: Config, *keys: str) -> Path:
    ref = cfg.sources["reference"]
    if keys[0] == "gene_sets":
        return Path(ref["root"]) / ref["msigdb_dir"] / ref["gene_sets"][keys[1]]
    if keys[0] == "protein":
        return Path(ref["root"]) / ref["protein_dir"] / keys[1]
    raise KeyError(keys)


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


# ---- data access -------------------------------------------------------------------------

def task_table(cfg: Config, task: str, cohort: str) -> pd.DataFrame:
    return pd.read_parquet(cfg.data_dir / "tasks" / f"{task}__{cohort}.parquet")


def cohort_matrix(cfg: Config, view: str, task: str, cohort: str):
    """(X float32 samples x features with NaN, features table, task table) for one cohort,
    over the features passing that cohort's extraction pre-filter."""
    tab = task_table(cfg, task, cohort)
    M, F = load_view(cfg.matrix_dir, view, tab["file_id"], cohort)
    return M.get(), F, tab


def matched_pairs(tab: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Row positions (tumour, normal) of cases with both, and the case IDs."""
    piv = tab.pivot_table(index="case_id", columns="role", values="file_id", aggfunc="first").dropna()
    pos = pd.Series(np.arange(len(tab)), index=tab["file_id"])
    return pos[piv["tumour"]].to_numpy(), pos[piv["normal"]].to_numpy(), piv.index.to_numpy()


def feature_values(cfg: Config, view: str, file_ids, feature_ids) -> np.ndarray:
    """Decoded values (len(file_ids) x len(feature_ids)) read straight from the stored matrix,
    in row chunks, without building the full cohort matrix."""
    d = view_dir(cfg.matrix_dir, view)
    samples = pd.read_parquet(d / "samples.parquet", columns=["file_id"])
    feats = pd.read_parquet(d / "features.parquet", columns=["feature_id"])
    rows = pd.Series(np.arange(len(samples)), index=samples["file_id"])[list(file_ids)].to_numpy()
    cols = pd.Series(np.arange(len(feats)), index=feats["feature_id"])[list(feature_ids)].to_numpy()
    mm = np.load(d / "values.npy", mmap_mode="r")
    order = np.argsort(rows)
    out = np.empty((len(rows), len(cols)), dtype=mm.dtype)
    for a in range(0, len(rows), 256):
        o = order[a:a + 256]
        out[o] = np.asarray(mm[rows[o]])[:, cols]
    return decode(out, load_meta(cfg.matrix_dir, view)["encoding"])


def gex_by_symbol(cfg: Config, file_ids, symbols: list[str]) -> pd.DataFrame:
    """log2(CPM+1) for the requested gene symbols (first gene per symbol), samples x symbols.
    Symbols absent from the stored GEX matrix are dropped."""
    M, F = load_view(cfg.matrix_dir, "gex", file_ids)
    F = F.drop_duplicates("symbol")
    col = F[F["symbol"].isin(symbols)]
    X = M.get(None, col.index.to_numpy())
    return pd.DataFrame(X, columns=col["symbol"].to_numpy(), index=list(file_ids))


@functools.lru_cache(maxsize=None)
def gene_strand_map(cfg_root: str) -> pd.Series:
    tx = pd.read_parquet(Path(cfg_root) / "transcripts.parquet", columns=["gene_id", "strand"])
    return tx.drop_duplicates("gene_id").set_index("gene_id")["strand"]


def gene_strand(cfg: Config) -> pd.Series:
    return gene_strand_map(str(ref_path(cfg, "protein", "x").parent))


def terminal_index(X: np.ndarray, F: pd.DataFrame, strand: pd.Series, min_usage: float) -> np.ndarray:
    """samples x genes position index (NaN where the kept siblings carry < min_usage)."""
    F = F.assign(start=F["region"].str.split(":").str[1].str.split("-").str[0].astype(int),
                 strand=F["gene_id"].map(strand)).reset_index(drop=True)
    F = F[F["strand"].notna()]
    F = F.assign(pos=np.where(F["strand"] == "+", F["start"], -F["start"]))
    F = F.assign(n=F.groupby("gene_id")["feature_id"].transform("size"))
    F = F[F["n"] >= 2]
    r = (F.groupby("gene_id")["pos"].rank(method="first") - 1) / (F["n"] - 1)
    g = pd.factorize(F["gene_id"])[0]
    U = np.nan_to_num(X[:, F.index.to_numpy()])
    G = g.max() + 1
    num = np.zeros((X.shape[0], G))
    den = np.zeros((X.shape[0], G))
    np.add.at(num.T, g, (U * r.to_numpy()).T)
    np.add.at(den.T, g, U.T)
    return np.where(den >= min_usage, num / np.where(den > 0, den, 1), np.nan)


# ---- SHAP --------------------------------------------------------------------------------

def full_importance(cfg: Config, task: str, view: str, cohort: str, n_direction: int = 300,
                    class_index: int | None = None) -> pd.DataFrame:
    """Importance of EVERY feature selected by any fold (mean |SHAP| over held-out samples,
    0 where not selected), with direction for the top n_direction. With class_index
    (cancer_type), only samples of that class are used. Cached under data_dir/bio_cache."""
    tag = "all" if class_index is None else f"class{class_index}"
    cache = cfg.data_dir / "bio_cache" / "importance" / f"{task}__{view}__{cohort}__{tag}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    tab = task_table(cfg, task, cohort)
    folds = load_folds(cfg.run_dir(task, view, cohort), tab["y"].to_numpy(), cfg.study["splits"]["outer_folds"])
    feats = pd.read_parquet(cfg.matrix_dir / view / "features.parquet")[["feature_id", "label", "symbol", "gene_id",
                                                                         "region", "chrom", "strand"]]
    mask = (lambda f: slice(None)) if class_index is None else (lambda f: f["y"] == class_index)
    t = importance_table(folds, feats, mask_of=mask, n_direction=n_direction)
    cache.parent.mkdir(parents=True, exist_ok=True)
    t.to_parquet(cache, index=False)
    return t


def gene_ranking(imp: pd.DataFrame) -> pd.DataFrame:
    """Gene-level ranking: max feature mean |SHAP| per gene (robust to many events per gene)."""
    g = (imp.groupby(["gene_id", "symbol"], dropna=False)
         .agg(importance=("mean_abs_shap", "max"), n_selected_features=("feature_id", "size"),
              top_feature=("label", "first"))
         .reset_index().sort_values(["importance", "gene_id"], ascending=[False, True]).reset_index(drop=True))
    return g


# ---- enrichment --------------------------------------------------------------------------

def read_gmt(path: Path) -> dict[str, set[str]]:
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = set(f[2:])
    return out


def pwf(n_events: np.ndarray, selected: np.ndarray) -> np.ndarray:
    """GOseq probability-weighting function: P(selected | number of measured events), binned on
    log2(events) and made monotone non-decreasing by isotonic regression."""
    b = np.floor(np.log2(np.maximum(n_events, 1))).astype(int)
    df = pd.DataFrame({"b": b, "s": selected.astype(float)}).groupby("b")["s"].agg(["mean", "size"])
    iso = IsotonicRegression(increasing=True, y_min=1e-6).fit(df.index.to_numpy(), df["mean"].to_numpy(),
                                                              sample_weight=df["size"].to_numpy())
    w = pd.Series(iso.predict(df.index.to_numpy()), index=df.index)
    return np.maximum(w.loc[b].to_numpy(), 1e-6)


def goseq_ora(selected: set[str], universe: pd.DataFrame, gene_sets: dict[str, set[str]], min_size: int,
              max_size: int, min_overlap: int, collection: str) -> pd.DataFrame:
    """Over-representation of `selected` symbols among `universe` (columns symbol, n_events).
    p_wallenius: Wallenius non-central hypergeometric with odds = mean PWF weight of set genes /
    mean weight of the other genes (GOseq). p_hyper: the uncorrected hypergeometric test."""
    u = universe.drop_duplicates("symbol").reset_index(drop=True)
    sym = u["symbol"].to_numpy()
    sel = np.isin(sym, list(selected))
    w = pwf(u["n_events"].to_numpy(), sel)
    M, N = len(u), int(sel.sum())
    idx = pd.Series(np.arange(M), index=sym)
    rows = []
    for name, genes in gene_sets.items():
        members = idx.reindex(list(genes)).dropna().astype(int).to_numpy()
        n = len(members)
        if n < min_size or n > max_size:
            continue
        k = int(sel[members].sum())
        if k < min_overlap:
            continue
        inside = np.zeros(M, bool)
        inside[members] = True
        odds = w[inside].mean() / w[~inside].mean()
        p_w = float(nchypergeom_wallenius(M, n, N, odds).sf(k - 1))
        p_h = float(hypergeom(M, n, N).sf(k - 1))
        expected = float(nchypergeom_wallenius(M, n, N, odds).mean())
        rows.append({"collection": collection, "term": name, "set_size": n, "overlap": k, "expected": expected,
                     "fold_enrichment": k / max(expected, 1e-9), "p_wallenius": p_w, "p_hyper": p_h,
                     "genes": ",".join(sorted(sym[members[sel[members]]]))})
    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values("p_wallenius").reset_index(drop=True)
    return out


def add_q(df: pd.DataFrame, n_tests_by_collection: dict[str, int]) -> pd.DataFrame:
    """BH within collection, counting every tested set (sets below min_overlap get p = 1)."""
    out = []
    for coll, g in df.groupby("collection"):
        p = np.concatenate([g["p_wallenius"].to_numpy(), np.ones(max(0, n_tests_by_collection[coll] - len(g)))])
        g = g.copy()
        g["q"] = bh(p)[:len(g)]
        out.append(g)
    return pd.concat(out, ignore_index=True) if out else df.assign(q=[])


def n_testable(universe_symbols: set[str], gene_sets: dict[str, set[str]], min_size: int, max_size: int) -> int:
    return sum(min_size <= len(g & universe_symbols) <= max_size for g in gene_sets.values())


# ---- protein consequences of cassette exons ----------------------------------------------

class ProteinMap:
    """GENCODE v45 transcripts/exons/CDS and protein features (1-based genomic coordinates)."""

    def __init__(self, cfg: Config):
        d = ref_path(cfg, "protein", "x").parent
        self.tx = pd.read_parquet(d / "transcripts.parquet")
        self.exons = pd.read_parquet(d / "exons.parquet")
        self.cds = pd.read_parquet(d / "cds.parquet")
        self.feat = pd.read_parquet(d / "protein_features.parquet")
        names = pd.read_csv(d / "feature_names.csv").set_index("acc")["name"]
        self.feat["name"] = self.feat["feature_id"].map(names).fillna(self.feat["source_db"])
        self.tx["prio"] = (~self.tx["mane_select"]).astype(int) * 2 + (~self.tx["canonical"]).astype(int)
        self.exon_index = self.exons.merge(self.tx[["transcript_id", "gene_id", "strand", "prio"]], on="transcript_id")
        self.exon_by_gene = dict(tuple(self.exon_index.groupby("gene_id")))
        self.has_cds = set(self.cds["transcript_id"])
        self.cds_by_tx = {t: g[["start", "end"]].to_numpy() for t, g in self.cds.groupby("transcript_id")}
        self.feat_by_tx = {t: g for t, g in self.feat.groupby("transcript_id")}

    def exon_consequence(self, gene_id: str, start: int, end: int, require_features: bool = False) -> dict:
        """Best transcript (MANE > canonical > other; coding first) with an exon exactly
        [start, end]; whether the exon is coding, frame-preserving, and which protein features
        its residues overlap. With require_features, a transcript that carries InterPro features
        is preferred over one that carries none -- the newest MANE transcripts are sometimes not
        in the feature table, and reading their empty annotation as "no domain" would be wrong."""
        g = self.exon_by_gene.get(gene_id)
        hit = g[(g["start"] == start) & (g["end"] == end)] if g is not None else None
        if hit is None or hit.empty:
            return {"annotated_exon": False}
        hit = hit.assign(coding=hit["transcript_id"].isin(self.has_cds),
                         featured=hit["transcript_id"].isin(self.feat_by_tx) if require_features else True)
        hit = hit.sort_values(["coding", "featured", "prio"], ascending=[False, False, True])
        t = hit.iloc[0]
        out = {"annotated_exon": True, "transcript_id": t["transcript_id"], "strand": t["strand"],
               "exon_length": int(end - start + 1)}
        cds = self.cds_by_tx.get(t["transcript_id"])
        if cds is None:
            return {**out, "coding": False}
        ov = np.clip(np.minimum(cds[:, 1], end) - np.maximum(cds[:, 0], start) + 1, 0, None)
        n_cds = int(ov.sum())
        if n_cds == 0:
            return {**out, "coding": False}
        # position of the exon's coding part along the CDS, 5'->3' of the transcript
        order = np.argsort(cds[:, 0]) if t["strand"] == "+" else np.argsort(-cds[:, 0])
        cum, first_nt = 0, None
        for i in order:
            s, e = cds[i]
            if ov[i] > 0 and first_nt is None:
                lo, hi = max(s, start), min(e, end)
                first_nt = cum + ((lo - s) if t["strand"] == "+" else (e - hi))
            cum += e - s + 1
        aa_lo, aa_hi = first_nt // 3 + 1, (first_nt + n_cds - 1) // 3 + 1
        f = self.feat_by_tx.get(t["transcript_id"])
        doms, regions = [], set()
        if f is not None:
            o = f[(f["start"] <= aa_hi) & (f["end"] >= aa_lo)]
            doms = sorted(set(o.loc[o["source_db"].isin(["pfam", "cdd"]), "name"]))
            regions = set(o["source_db"])
        return {**out, "coding": True, "fully_coding": n_cds == out["exon_length"],
                "frame_preserving": n_cds % 3 == 0, "aa_start": int(aa_lo), "aa_end": int(aa_hi),
                "protein_length_aa": int(cum // 3), "domains": "; ".join(doms), "in_domain": bool(doms),
                "in_disordered": "mobidblite" in regions, "in_transmembrane": "tmhmm" in regions}
