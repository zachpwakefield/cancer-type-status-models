"""Step 42 (D-52): does genome sequence alone say which exons MSI tumours skip?

Mismatch-repair-deficient (MSI) tumours shorten mononucleotide repeats. A shortened poly-T tract
in the polypyrimidine tract of an acceptor weakens that 3' splice site, and the exon downstream
is skipped (HSP110 exon 9, MRE11 exon 5; genome-wide in colorectal cancer, Jonchere et al. 2024).
This step asks whether that cis route is visible across UCEC, STAD and COAD from the genome
sequence alone, with no model and no label used to choose the exons.

  annotate  label-blind: for every cassette exon of the subtype SE matrix, the longest sense-strand
            homopolymer run of each base overlapping three intronic windows (acceptor -1..-50,
            deep -201..-400 upstream, donor +1..+50 downstream), from the GRCh38 primary assembly
  test      the predictions of config/followup.yaml: msi_cis (M1-M8, K1-K2)

Coordinates: a feature_id reads gene#chrom:start-end#strand;flankA;flankB with 0-based half-open
intervals and the two flanking exons in genome order. Positions are counted in transcript
orientation from the exon: -1 is the last intron base before it, +1 the first after it.

Outputs: data.nosync/ctsm_v1/msi/ (per-sample scores, kept local)
         results/subtype/msi_cis/ (aggregates; figure U10: scripts/36_subtype_figures.py)
Usage: python scripts/42_msi_repeats.py [--only annotate|test]
"""
from __future__ import annotations

import argparse
import gzip
import re

import numpy as np
import pandas as pd
import yaml
from _common import REPO, log

FASTA = "/Users/zacharywakefield/Downloads/GRCh38.primary_assembly.genome.fa.gz"
DATA = REPO / "data.nosync" / "ctsm_v1"
LOCAL = DATA / "msi"
OUT = REPO / "results" / "subtype" / "msi_cis"
WINDOWS = {"acc": (1, 50), "deep": (201, 400), "don": (1, 50)}   # 1-based, inclusive
REACH = 30                     # bases read past each window edge so a run is measured in full
MIN_INTRON = {"acc": 60, "deep": 500, "don": 60}
COMP = str.maketrans("ACGTN", "TGCAN")


def parse_features(f: pd.DataFrame) -> pd.DataFrame:
    p = f["feature_id"].str.split("#", expand=True)
    ex = p[1].str.split(":", expand=True)
    se = ex[1].str.split("-", expand=True).astype(int)
    fl = p[2].str.split(";", expand=True)
    a = fl[1].str.split("-", expand=True).astype(int)
    b = fl[2].str.split("-", expand=True).astype(int)
    return pd.DataFrame({"feature_id": f["feature_id"], "symbol": f["symbol"], "chrom": ex[0],
                         "strand": fl[0], "es": se[0], "ee": se[1], "left_end": a[1],
                         "right_start": b[0]})


def read_genome(chroms: set[str]):
    """Yield (chrom, sequence) for the wanted chromosomes, streaming the gzip once."""
    name, parts = None, []
    with gzip.open(FASTA, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name in chroms:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif name in chroms:
                parts.append(line.rstrip())
    if name in chroms:
        yield name, "".join(parts).upper()


def max_runs(seq: str, lo: int, hi: int) -> dict[str, int]:
    """Longest run of each base in `seq` (index 0 = position 1) overlapping positions lo..hi."""
    out = {b: 0 for b in "ACGT"}
    for m in re.finditer(r"A+|C+|G+|T+", seq):
        s, e = m.start() + 1, m.end()                 # 1-based inclusive
        if e >= lo and s <= hi:
            b = m.group()[0]
            out[b] = max(out[b], e - s + 1)
    return out


def window_seq(chrom_seq: str, r, side: str) -> tuple[str | None, int]:
    """Sense-orientation intron sequence next to the exon, index 0 = position -1 (upstream) or
    +1 (downstream), with the intron length. Upstream intron for 'acc' and 'deep'."""
    plus = r.strand == "+"
    upstream = side in ("acc", "deep")
    if plus == upstream:                              # intron left of the exon in the genome
        lo, hi = r.left_end, r.es
        n = hi - lo
        far = WINDOWS[side][1] + REACH
        s = chrom_seq[max(lo, hi - far):hi][::-1]
    else:                                             # intron right of the exon
        lo, hi = r.ee, r.right_start
        n = hi - lo
        far = WINDOWS[side][1] + REACH
        s = chrom_seq[lo:min(hi, lo + far)]
    if not plus:
        s = s.translate(COMP)
    return s, n


def annotate() -> pd.DataFrame:
    f = parse_features(pd.read_parquet(DATA / "subtype" / "matrices" / "se" / "features.parquet"))
    rows = []
    for chrom, seq in read_genome(set(f["chrom"])):
        sub = f[f["chrom"].eq(chrom)]
        for r in sub.itertuples(index=False):
            rec = {"feature_id": r.feature_id}
            for side, (lo, hi) in WINDOWS.items():
                s, n = window_seq(seq, r, side)
                rec[f"{side}_intron_len"] = n
                ok = n >= MIN_INTRON[side]
                runs = max_runs(s, lo, hi) if ok else {b: np.nan for b in "ACGT"}
                for b, v in runs.items():
                    rec[f"{side}_{b}"] = v
                if side == "acc":
                    w = s[:hi] if ok else ""
                    rec["acc_pyr_frac"] = (w.count("C") + w.count("T")) / len(w) if w else np.nan
                    rec["acc_ag"] = s[:2][::-1] if ok else ""      # -2,-1 read 5'->3'
            rows.append(rec)
        log(f"{chrom}: {len(sub):,} exons")
    ann = f.merge(pd.DataFrame(rows), on="feature_id", how="left")
    LOCAL.mkdir(parents=True, exist_ok=True)
    ann.to_parquet(LOCAL / "se_repeats.parquet", index=False)
    log(f"annotated {len(ann):,} exons; canonical AG at -2,-1: {ann['acc_ag'].eq('AG').mean():.3f}")
    return ann


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["annotate", "test"])
    a = ap.parse_args()
    if a.only in (None, "annotate"):
        annotate()


if __name__ == "__main__":
    main()
