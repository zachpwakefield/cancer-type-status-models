"""Gene annotation, event-ID parsing and the genome-scope rule (decision D-06).

Adapted from the sister repository zachpwakefield/cancer-survival-predictions
(src/survpred/annotation.py, events.py), where the ID grammars were verified against every
TCGA event ID.

ID grammars:
  rMATS (se, mxe, a3ss, a5ss, ri):  ENSG...v#chr:start-end#strand;start-end;start-end[;...]
        0-based half-open, GENCODE v36. Every block is parsed:
          se   : cassette exon; genomic-upstream exon; genomic-downstream exon
          mxe  : 1st exon; 2nd exon; upstream exon; downstream exon
                 (PSI = inclusion of the TRANSCRIPT-upstream exon: 1st on +, 2nd on -)
          a3ss : long exon; short exon; flanking exon
          a5ss : long exon; short exon; flanking exon
          ri   : exon-intron-exon envelope; upstream exon; downstream exon
  HITindex (afe, ale, hit):         ENSG...v;chr:start-end;afe|ale|hit   (1-based metaexon)
  gex:                              ENSG... (unversioned)
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd

from .config import RMATS_VIEWS, Config


# ---- genes ----------------------------------------------------------------------------------

def _gtf_gene_chrom(path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            attr = f[8]
            i = attr.index('gene_id "') + 9
            gid = attr[i:attr.index('"', i)]
            rows.append((gid.split(".")[0], f[0], gid.endswith("_PAR_Y")))
    df = pd.DataFrame(rows, columns=["gene_id", "chrom", "par_y"])
    # PAR genes appear on X and Y with the same base ID; keep the X copy
    return df.sort_values(["gene_id", "par_y"]).drop_duplicates("gene_id")[["gene_id", "chrom"]]


def gene_annotation(cfg: Config) -> pd.DataFrame:
    """gene_id (unversioned), chrom, biotype, symbol. GENCODE v45 first, GDC v36 as fallback."""
    chrom = _gtf_gene_chrom(cfg.gencode_gtf)
    v45 = pd.read_csv(cfg.tcga_path("annotation", "gene_biotypes_v45")).rename(
        columns={"gene_id_unversioned": "gene_id", "gene_type": "biotype"})
    v36 = pd.read_csv(cfg.tcga_path("annotation", "gdc_genes_v36"))
    v36 = (v36[~v36.gene_id.str.endswith("_PAR_Y")]
           .drop(columns="gene_id")
           .rename(columns={"gene_id_unversioned": "gene_id", "gene_type": "biotype_v36", "gene_name": "symbol_v36"})
           [["gene_id", "biotype_v36", "symbol_v36"]].drop_duplicates("gene_id"))
    ann = v45.merge(chrom, on="gene_id", how="outer").merge(v36, on="gene_id", how="outer")
    ann["biotype"] = ann["biotype"].fillna(ann["biotype_v36"])
    ann["symbol"] = ann["symbol"].fillna(ann["symbol_v36"]).fillna(ann["gene_id"])
    return ann[["gene_id", "chrom", "biotype", "symbol"]].drop_duplicates("gene_id").reset_index(drop=True)


def scope_mask(genes: pd.Series, chroms: pd.Series, ann: pd.DataFrame, cfg: Config) -> np.ndarray:
    """D-06: feature chromosome in scope (chr1-22, X) and host-gene biotype allowed."""
    allowed_chr = set(cfg.study["features"]["chromosomes"])
    allowed_bt = set(cfg.study["features"]["host_biotypes"])
    bt = genes.map(ann.set_index("gene_id")["biotype"])
    return (chroms.isin(allowed_chr) & bt.isin(allowed_bt)).to_numpy()


# ---- event IDs ------------------------------------------------------------------------------

def parse_rmats_ids(ids: pd.Series, view: str) -> pd.DataFrame:
    parts = ids.str.split("#", n=2, expand=True)
    if parts.shape[1] != 3 or parts.isna().any().any():
        raise ValueError(f"malformed rMATS IDs in {view}")
    gene = parts[0].str.split(".", n=1).str[0]
    chrom = parts[1].str.split(":", n=1).str[0]
    first = parts[1].str.split(":", n=1).str[1]
    strand_rest = parts[2].str.split(";", n=1, expand=True)
    strand = strand_rest[0]
    blocks = (first + ";" + strand_rest[1]).str.split(";")
    n_expected = {"se": 3, "mxe": 4, "a3ss": 3, "a5ss": 3, "ri": 3}[view]
    if not (blocks.str.len() == n_expected).all():
        raise ValueError(f"{view}: expected {n_expected} intervals per event")
    coords = [blocks.str[i].str.split("-", expand=True).astype(np.int64) for i in range(n_expected)]
    starts = np.column_stack([c[0].to_numpy() for c in coords])
    ends = np.column_stack([c[1].to_numpy() for c in coords])
    if view in ("se", "a3ss", "a5ss"):
        ks, ke = starts[:, 0], ends[:, 0]
    elif view == "mxe":
        plus = (strand == "+").to_numpy()
        ks = np.where(plus, starts[:, 0], starts[:, 1])
        ke = np.where(plus, ends[:, 0], ends[:, 1])
    else:                                            # ri: the retained intron
        ks, ke = ends[:, 1], starts[:, 2]
    return pd.DataFrame({
        "feature_id": ids.to_numpy(),
        "gene_id": gene.to_numpy(),
        "chrom": chrom.to_numpy(),
        "strand": strand.to_numpy(),
        # the variable segment, 1-based inclusive for display (rMATS is 0-based half-open)
        "region": (chrom + ":" + pd.Series(ks + 1, index=ids.index).astype(str) + "-"
                   + pd.Series(ke, index=ids.index).astype(str)).to_numpy(),
    })


def parse_hit_ids(ids: pd.Series, view: str) -> pd.DataFrame:
    p = ids.str.extract(r"^(ENSG\d+)(?:\.\d+)?(?:_PAR_Y)?;(chr[^:]+):(\d+)-(\d+);(afe|ale|hit)$")
    if p.isna().any().any():
        bad = ids[p.isna().any(axis=1)].head(3).tolist()
        raise ValueError(f"malformed HITindex IDs in {view}: {bad}")
    if not (p[4] == view).all():
        raise ValueError(f"{view}: ID suffix mismatch")
    return pd.DataFrame({
        "feature_id": ids.to_numpy(),
        "gene_id": p[0].to_numpy(),
        "chrom": p[1].to_numpy(),
        "strand": None,
        "region": (p[1] + ":" + p[2] + "-" + p[3]).to_numpy(),   # metaexon, already 1-based
    })


def parse_ids(ids: pd.Series, view: str) -> pd.DataFrame:
    return parse_rmats_ids(ids, view) if view in RMATS_VIEWS else parse_hit_ids(ids, view)


def load_metaexon_terminal_support(path) -> pd.DataFrame:
    """hg38_gen_meta -> per-metaexon annotated first/last-exon transcript counts (1-based keys)."""
    meta = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "info"])
    info = meta["info"].str.extract(r";(ENSG\d+)\.\d+(?:_PAR_Y)?;TXPT:(\d+);FE:(\d+);internal:(\d+);LE:(\d+)")
    out = pd.DataFrame({
        "gene_id": info[0],
        "key": meta["chrom"] + ":" + meta["start"].astype(str) + "-" + meta["end"].astype(str),
        "FE": info[2].astype(float),
        "LE": info[4].astype(float),
    })
    return out.groupby(["gene_id", "key"], as_index=False)[["FE", "LE"]].max()
