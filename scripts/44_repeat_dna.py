"""Step 44 (D-56), part 1: which loci and which exome BAMs, for reading the repeats themselves.

D-52..D-55 inferred repeat contraction from the reference genome and the MSI label. The decisive
check reads the acceptor T-runs in each tumour's exome sequencing. This script prepares that and
handles no credential and no controlled-access data:

  loci      label-blind: the acceptor T-run of SE repeat exons (all runs >= 20, a seeded random 250
            of 14-19 and of 10-13) and of A3SS repeat acceptors, with exact genome coordinates and
            8-nt flanks; BAT-25 (KIT) and BAT-26 (MSH2), the clinical MSI markers, as positive
            controls
  manifest  which exome BAMs (open GDC metadata), in fetch order: MSI tumour/normal pairs first,
            then POLE, then MSS; and the JSON region payload for the GDC slicing endpoint

The slices themselves are fetched by the data user, with their own token, using the GDC's
documented slicing request (docs.gdc.cancer.gov/API/Users_Guide/BAM_Slicing/), into
data.nosync/ctsm_v1/msi/slices/<file_id>.bam.

Outputs: data.nosync/ctsm_v1/msi/dna_loci.csv, dna_manifest.csv, dna_regions.json;
         results/subtype/msi_cis/dna/loci.csv (genome coordinates only)
Usage: python scripts/44_repeat_dna.py --only loci|manifest
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from _common import log

_spec = importlib.util.spec_from_file_location("msi42", Path(__file__).with_name("42_msi_repeats.py"))
m42 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m42)

DATA, LOCAL = m42.DATA, m42.LOCAL
OUT = m42.OUT / "dna"
META = DATA / "external" / "gdc_wxs" / "wxs_bam_metadata.csv"
COHORTS = ["UCEC_TCGA", "STAD_TCGA", "COAD_TCGA"]
FLANK, PAD = 8, 10
SEED = 0
# clinical MSI markers: the longest run of `base` in a window where the marker is documented -
# BAT-26 in MSH2 intron 5, just after exon 5 (MSH2-201 exon 5 ends at 47,414,418); BAT-25 in KIT intron 16
MARKERS = {"BAT-26": ("chr2", 47414418, 47414718, "A"), "BAT-25": ("chr4", 54731900, 54732200, "T")}


def run_coords(sense: str, lo: int, hi: int, base: str = "T") -> tuple[int, int] | None:
    """1-based (s, e), s <= e, of the longest run of `base` overlapping positions lo..hi of a
    sense-oriented intron string (index 0 = position 1, counting away from the exon)."""
    best, i, n = None, 0, len(sense)
    while i < n:
        if sense[i] == base:
            j = i
            while j < n and sense[j] == base:
                j += 1
            s, e = i + 1, j
            if e >= lo and s <= hi and (best is None or e - s > best[1] - best[0]):
                best = (s, e)
            i = j
        else:
            i += 1
    return best


def genome_run(strand: str, acc: int, s: int, e: int) -> tuple[int, int]:
    """0-based half-open genome interval of acceptor-relative positions -s..-e. `acc` is the
    exon start (0-based) on '+', the exon end (0-based, exclusive) on '-'."""
    return (acc - e, acc - s + 1) if strand == "+" else (acc + s - 1, acc + e)


def make_locus(seq: str, chrom: str, start: int, end: int, base: str, **info) -> dict | None:
    if seq[start:end] != base * (end - start) or seq[start - 1] == base or seq[end] == base:
        return None
    return {"chrom": chrom, "start": start, "end": end, "base": base, "ref_len": end - start,
            "left": seq[start - FLANK:start], "right": seq[end:end + FLANK],
            "region": f"{chrom}:{start + 1 - PAD}-{end + PAD}", **info}


def loci() -> pd.DataFrame:
    ex = pd.read_parquet(LOCAL / "exon_effects.parquet")
    se = ex[ex.acc_T >= 10].drop_duplicates("feature_id").copy()
    se["acc"] = np.where(se.strand.eq("+"), se.es, se.ee)
    se = se.drop_duplicates(["chrom", "acc", "strand"]).reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    pick = [se[se.acc_T >= 20]]
    for lo, hi in [(14, 19), (10, 13)]:
        d = se[(se.acc_T >= lo) & (se.acc_T <= hi)]
        pick.append(d.iloc[np.sort(rng.choice(len(d), size=min(250, len(d)), replace=False))])
    se = pd.concat(pick, ignore_index=True)
    a3 = pd.read_parquet(LOCAL / "a3ss_effects.parquet")
    a3 = a3[a3.side.ne(0) & (a3.run >= 10)].drop_duplicates("feature_id")
    a3ann = pd.read_parquet(LOCAL / "a3ss_repeats.parquet").set_index("feature_id")
    a3 = a3.join(a3ann[["chrom", "strand", "ls", "le", "ss", "se"]], on="feature_id")
    plus = a3["strand"].eq("+")                    # a3.le would be DataFrame.le, not the column
    a3["acc"] = np.where(a3["side"].eq(1), np.where(plus, a3["ls"], a3["le"]), np.where(plus, a3["ss"], a3["se"]))
    far = m42.WINDOWS["acc"][1] + m42.REACH
    rows = []
    for chrom, seq in m42.read_genome(set(se.chrom) | set(a3.chrom) | {m[0] for m in MARKERS.values()}):
        for kind, d in [("SE", se[se.chrom.eq(chrom)]), ("A3SS", a3[a3.chrom.eq(chrom)])]:
            for r in d.itertuples(index=False):
                acc = int(r.acc)
                w = seq[max(0, acc - far):acc][::-1] if r.strand == "+" else seq[acc:acc + far].translate(m42.COMP)
                rc = run_coords(w, 1, m42.WINDOWS["acc"][1])
                if rc is None:
                    continue
                g0, g1 = genome_run(r.strand, acc, *rc)
                loc = make_locus(seq, chrom, g0, g1, "T" if r.strand == "+" else "A", kind=kind,
                                 feature_id=r.feature_id, symbol=r.symbol, strand=r.strand,
                                 run_from_exon=f"-{rc[0]}..-{rc[1]}",
                                 acceptor_side=("distal" if r.side == 1 else "proximal") if kind == "A3SS" else "")
                if loc:
                    rows.append(loc)
        for name, (mc, a, b, base) in MARKERS.items():
            if mc == chrom:
                best = max(((x.start(), x.end()) for x in re.finditer(f"{base}+", seq[a:b])), key=lambda x: x[1] - x[0])
                rows.append(make_locus(seq, chrom, a + best[0], a + best[1], base, kind="marker", feature_id=name,
                                       symbol=name, strand="+", run_from_exon="", acceptor_side=""))
        log(f"loci {chrom}")
    L = pd.DataFrame(rows).drop_duplicates(["chrom", "start", "end"]).reset_index(drop=True)
    L.to_csv(LOCAL / "dna_loci.csv", index=False)
    OUT.mkdir(parents=True, exist_ok=True)
    L.to_csv(OUT / "loci.csv", index=False)
    (LOCAL / "dna_regions.json").write_text(json.dumps({"regions": L.region.tolist()}))
    log(f"{len(L)} loci: " + ", ".join(f"{k} {v}" for k, v in L.kind.value_counts().items())
        + f"; markers {L[L.kind.eq('marker')][['symbol', 'chrom', 'start', 'ref_len']].values.tolist()}")
    return L


def manifest() -> pd.DataFrame:
    meta = pd.read_csv(META)
    rng = np.random.default_rng(SEED)
    rows = []
    for c in COHORTS:
        t = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{c}.parquet")
        t["group"] = np.where(t.label.eq("MSI"), "MSI", np.where(t.label.eq("POLE"), "POLE", "rest"))
        mss = t[t.group.eq("rest")]
        mss = mss.iloc[rng.permutation(len(mss))].head(40)
        sel = pd.concat([t[t.group.isin(["MSI", "POLE"])], mss])
        with_normal = set(t.loc[t.group.isin(["MSI", "POLE"]), "case_id"]) | set(mss.case_id.head(20))
        for r in sel.itertuples(index=False):
            m = meta[meta.case_id.eq(r.case_id)]
            tum = m[m.sample_id.eq(r.sample_id)]
            if tum.empty:
                tum = m[m.sample_type.eq("Primary Tumor")]
            nor = m[m.sample_type.eq("Blood Derived Normal")]
            if nor.empty:
                nor = m[m.sample_type.eq("Solid Tissue Normal")]
            for role, f in [("tumour", tum)] + ([("normal", nor)] if r.case_id in with_normal else []):
                if f.empty:
                    continue
                f = f.sort_values("file_size").iloc[-1]
                rows.append({"cohort": c, "case_id": r.case_id, "group": r.group, "role": role,
                             "file_id": f.file_id, "sample_id": f.sample_id, "sample_type": f.sample_type})
    m = pd.DataFrame(rows)
    m["paired"] = m.groupby("case_id").role.transform("nunique").eq(2)
    m["tier"] = m.group.map({"MSI": 0, "POLE": 1, "rest": 2}) * 2 + (~m.paired).astype(int)
    case_order = pd.Series(rng.random(m.case_id.nunique()), index=m.case_id.unique())
    m["r"] = m.case_id.map(case_order)
    m = m.sort_values(["tier", "r", "role"], ascending=[True, True, False]).drop(columns="r").reset_index(drop=True)
    m["order"] = np.arange(len(m))
    m.to_csv(LOCAL / "dna_manifest.csv", index=False)
    log(f"manifest: {len(m)} BAMs, {m.case_id.nunique()} patients; "
        + ", ".join(f"{g}/{r} {n}" for (g, r), n in m.groupby(["group", "role"]).size().items()))
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["loci", "manifest"], required=True)
    a = ap.parse_args()
    loci() if a.only == "loci" else manifest()


if __name__ == "__main__":
    main()
