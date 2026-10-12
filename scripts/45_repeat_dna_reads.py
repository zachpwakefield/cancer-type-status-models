"""Step 45 (D-56), part 2: call the acceptor T-run length in each read and test G1-G6.

Reads the BAM slices of data.nosync/ctsm_v1/msi/slices/<file_id>.bam (fetched by the data user,
see step 44) with a small BAM reader (BGZF is multi-member gzip; samtools and pysam are not
installed), and calls each read's run length between the locus's exact 8-nt flanks.

  genotype  per sample and locus, the run length of every informative read (local only)
  test      the predictions of config/followup.yaml: msi_dna (G1-G6), family D

Outputs: data.nosync/ctsm_v1/msi/dna_reads.parquet, dna_tumour_locus.parquet (per sample, local);
         results/subtype/msi_cis/dna/ (aggregates only)
Usage: python scripts/45_repeat_dna_reads.py --only genotype|test
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import struct
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from _common import REPO, log
from scipy.stats import binomtest, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

from ctsm.matrix_io import load_view

_spec = importlib.util.spec_from_file_location("msi42", Path(__file__).with_name("42_msi_repeats.py"))
m42 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m42)

CFG = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["msi_dna"]
DATA, LOCAL = m42.DATA, m42.LOCAL
SLICES = LOCAL / "slices"
OUT = m42.OUT / "dna"
MIN_MAPQ, MIN_READS, MIN_TUMOURS = 20, 8, 10
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
_NT = np.frombuffer(b"=ACMGRSVTWYHKDBN", dtype=np.uint8)


# ============================================================ reading
def read_bam(path: Path):
    """Yield (chrom, pos0, mapq, flag, seq) for every record of a BAM file."""
    with gzip.open(path, "rb") as fh:
        data = fh.read()
    if data[:4] != b"BAM\x01":
        raise ValueError(f"{path}: not a BAM")
    off = 8 + struct.unpack_from("<i", data, 4)[0]
    n_ref = struct.unpack_from("<i", data, off)[0]
    off += 4
    refs = []
    for _ in range(n_ref):
        ln = struct.unpack_from("<i", data, off)[0]
        refs.append(data[off + 4:off + 4 + ln - 1].decode())
        off += 8 + ln
    n = len(data)
    while off < n:
        bs, ref, pos, lrn, mapq, _, ncig, flag, lseq = struct.unpack_from("<iiiBBHHHi", data, off)
        p = off + 36 + lrn + 4 * ncig
        packed = np.frombuffer(data, dtype=np.uint8, count=(lseq + 1) // 2, offset=p)
        codes = np.empty(2 * len(packed), dtype=np.uint8)
        codes[0::2], codes[1::2] = packed >> 4, packed & 15
        yield (refs[ref] if ref >= 0 else None), pos, mapq, flag, _NT[codes[:lseq]].tobytes().decode()
        off += 4 + bs


def run_length(seq: str, left: str, base: str, right: str) -> int | None:
    """Length of the `base` run between exact `left` and `right` flanks in a read, or None."""
    i = seq.find(left)
    while i != -1:
        j = k = i + len(left)
        while k < len(seq) and seq[k] == base:
            k += 1
        if k > j and seq.startswith(right, k):
            return k - j
        i = seq.find(left, i + 1)
    return None


def genotype_bam(path: Path, L: pd.DataFrame) -> dict[int, list[int]]:
    """{locus row of L: [run length of each informative read]}."""
    by_chrom = {}
    for c, d in L.groupby("chrom"):
        d = d.sort_values("start")
        by_chrom[c] = (d["start"].to_numpy(), list(zip(d.index, d["left"], d["base"], d["right"])))
    out: dict[int, list[int]] = {}
    for chrom, pos, mapq, flag, seq in read_bam(path):
        if chrom not in by_chrom or mapq < MIN_MAPQ or flag & SKIP_FLAGS:
            continue
        st, rows = by_chrom[chrom]
        for k in range(np.searchsorted(st, pos - len(seq)), np.searchsorted(st, pos + len(seq))):
            idx, left, base, right = rows[k]
            n = run_length(seq, left, base, right)
            if n is not None:
                out.setdefault(int(idx), []).append(n)
    return out


def genotype() -> None:
    L = pd.read_csv(LOCAL / "dna_loci.csv")
    m = pd.read_csv(LOCAL / "dna_manifest.csv")
    rows, per_bam = [], []
    for r in m.itertuples(index=False):
        path = SLICES / f"{r.file_id}.bam"
        if not path.exists():
            continue
        try:
            calls = genotype_bam(path, L)
        except (EOFError, OSError, ValueError, struct.error) as e:       # truncated or not a BAM
            log(f"skipped {r.order:>4} {path.name}: {type(e).__name__} {e}")
            continue
        for k, lens in calls.items():
            rows.append({"file_id": r.file_id, "case_id": r.case_id, "cohort": r.cohort, "group": r.group,
                         "role": r.role, "locus": k, "lengths": lens})
        n = np.array([len(v) for v in calls.values()] or [0])
        per_bam.append({"order": r.order, "cohort": r.cohort, "group": r.group, "role": r.role,
                        "mb": path.stat().st_size / 1e6, "loci_with_reads": len(calls),
                        "usable_loci": int((n >= MIN_READS).sum()), "median_reads": float(np.median(n))})
        log(f"genotyped {r.order:>4} {r.cohort} {r.group} {r.role}: {per_bam[-1]['mb']:.1f} MB, "
            f"{per_bam[-1]['usable_loci']} usable loci")
    g = pd.DataFrame(rows, columns=["file_id", "case_id", "cohort", "group", "role", "locus", "lengths"])
    g.to_parquet(LOCAL / "dna_reads.parquet", index=False)
    pb = pd.DataFrame(per_bam)
    OUT.mkdir(parents=True, exist_ok=True)
    pb.to_csv(OUT / "slices_by_bam.csv", index=False, float_format="%.4g")
    if len(pb):
        log(f"{len(pb)} BAMs, {pb.mb.sum():.0f} MB ({pb.mb.median():.1f} MB median); usable loci per BAM: "
            f"median {pb.usable_loci.median():.0f} of {len(L)}")


# ============================================================ contraction
def mode(x) -> int:
    v, c = np.unique(np.asarray(x), return_counts=True)
    return int(v[np.argmax(c)])


def contraction_table(g: pd.DataFrame) -> pd.DataFrame:
    """One row per usable tumour-locus cell: germline mode, c (stutter-corrected share shorter)."""
    g = g[g.lengths.map(len) >= MIN_READS]
    nor = g[g.role.eq("normal")]
    pooled = nor.groupby("locus").lengths.agg(lambda s: np.concatenate(s.to_list()))
    own = nor.set_index(["case_id", "locus"]).lengths
    rows = []
    for r in g[g.role.eq("tumour")].itertuples(index=False):
        if (r.case_id, r.locus) in own.index:
            ref, src = np.asarray(own.loc[(r.case_id, r.locus)]), "own normal"
        elif r.locus in pooled.index:
            ref, src = pooled.loc[r.locus], "pooled normals"
        else:
            continue
        gm = mode(ref)
        t = np.asarray(r.lengths)
        rows.append({"file_id": r.file_id, "case_id": r.case_id, "cohort": r.cohort, "group": r.group, "locus": r.locus,
                     "n_reads": len(t), "germline_mode": gm, "reference": src,
                     "c": float((t < gm).mean() - (ref < gm).mean()), "mean_shift": float(t.mean() - ref.mean())})
    return pd.DataFrame(rows)


def exon_percentiles(L: pd.DataFrame, tumours: pd.DataFrame) -> pd.DataFrame:
    """Each SE locus exon's PSI percentile within its cohort (all cohort patients), for the tumours."""
    out = []
    se = L[L.kind.eq("SE")]
    for c in tumours.cohort.unique():
        tab = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{c}.parquet")
        vm, feats = load_view(DATA / "subtype" / "matrices", "se", tab["file_id"])
        cols = feats.index[feats.feature_id.isin(se.feature_id)].to_numpy()
        P = pd.DataFrame(vm.get(cols=cols), columns=feats.feature_id.iloc[cols].to_numpy(),
                         index=tab["case_id"]).rank(pct=True)
        keep = tumours.loc[tumours.cohort.eq(c), "case_id"].unique()
        long = P.loc[P.index.isin(keep)].stack().rename("pct").reset_index()
        long.columns = ["case_id", "feature_id", "pct"]
        out.append(long.assign(cohort=c))
    return pd.concat(out, ignore_index=True)


def per_locus_rho(ct: pd.DataFrame, pct: pd.DataFrame, L: pd.DataFrame, loci: pd.Index) -> pd.DataFrame:
    d = ct[ct.group.eq("MSI") & ct.locus.isin(loci)].merge(L[["feature_id"]], left_on="locus", right_index=True)
    d = d.merge(pct, on=["case_id", "feature_id", "cohort"])
    rows = []
    for k, x in d.groupby("locus"):
        if len(x) >= MIN_TUMOURS and x.c.nunique() > 1 and x.pct.nunique() > 1:
            rows.append({"locus": k, "n_msi_tumours": len(x), "rho_c_pct": spearmanr(x.c, x.pct).statistic})
    return pd.DataFrame(rows)


def test() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L = pd.read_csv(LOCAL / "dna_loci.csv")
    g = pd.read_parquet(LOCAL / "dna_reads.parquet")
    ct = contraction_table(g)
    ct.to_parquet(LOCAL / "dna_tumour_locus.parquet", index=False)
    rep = L[L.kind.isin(["SE", "A3SS"])].index
    burden = ct[ct.locus.isin(rep)].groupby(["case_id", "cohort", "group"]).c.mean().rename("burden").reset_index()
    pred = []

    def add(pid, cohort, effect, sign, p, n, what):
        pred.append({"id": pid, "cohort": cohort, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    for c in ["UCEC_TCGA", "STAD_TCGA", "COAD_TCGA"]:
        b = burden[burden.cohort.eq(c)]
        x, y = b.loc[b.group.eq("MSI"), "burden"], b.loc[b.group.eq("rest"), "burden"]
        if len(x) and len(y):
            r = mannwhitneyu(x, y)
            add("G1", c, r.statistic / (len(x) * len(y)) - 0.5, 1, r.pvalue, f"MSI {len(x)}, rest {len(y)}",
                "AUROC(burden higher in MSI) - 0.5")
    b = burden[burden.cohort.eq("UCEC_TCGA")]
    x, y = b.loc[b.group.eq("MSI"), "burden"], b.loc[b.group.eq("POLE"), "burden"]
    if len(x) and len(y):
        r = mannwhitneyu(x, y)
        add("G2", "UCEC_TCGA", r.statistic / (len(x) * len(y)) - 0.5, 1, r.pvalue, f"MSI {len(x)}, POLE {len(y)}",
            "AUROC(burden higher in MSI than POLE) - 0.5")
    pct = exon_percentiles(L, ct[ct.group.eq("MSI")])
    se = L[L.kind.eq("SE")]
    for pid, lo, hi, sign in [("G3", 10, 19, -1), ("G5", 20, 999, 1)]:
        rl = per_locus_rho(ct, pct, L, se[(se.ref_len >= lo) & (se.ref_len <= hi)].index)
        rl.assign(set=pid).to_csv(OUT / f"per_locus_rho_{pid}.csv", index=False)
        if len(rl):
            pos, neg = int((rl.rho_c_pct > 0).sum()), int((rl.rho_c_pct < 0).sum())
            add(pid, "pooled", pos / max(1, pos + neg) - 0.5, sign, binomtest(pos, pos + neg).pvalue,
                f"{len(rl)} loci, {neg} with rho < 0; median rho {rl.rho_c_pct.median():.3f}",
                "share of loci with Spearman(c, exon PSI percentile) > 0, minus 0.5; MSI tumours")
    ex = pd.read_parquet(LOCAL / "exon_effects.parquet").groupby("feature_id").delta_adj.mean()
    msi_c = ct[ct.group.eq("MSI")].groupby("locus").agg(mean_c=("c", "mean"), n=("c", "size"))
    msi_c = msi_c[msi_c.n >= MIN_TUMOURS].join(L[["kind", "feature_id", "ref_len"]])
    d4 = msi_c[msi_c.kind.eq("SE") & msi_c.ref_len.between(10, 19)].join(ex, on="feature_id").dropna()
    if len(d4) > 3:
        r = spearmanr(d4.mean_c, d4.delta_adj)
        add("G4", "pooled", float(r.statistic), -1, float(r.pvalue), f"{len(d4)} loci",
            "Spearman(mean c in MSI, exon delta_adj), runs 10-19")
    d6 = msi_c[msi_c.kind.isin(["SE", "A3SS"])]
    if len(d6) > 3:
        r = spearmanr(d6.mean_c, d6.ref_len)
        add("G6", "pooled", float(r.statistic), 1, float(r.pvalue), f"{len(d6)} loci", "Spearman(mean c in MSI, run length)")
    pr = pd.DataFrame(pred)
    if len(pr):
        pr["q"] = multipletests(pr.p, method="fdr_bh")[1]
        pr["verdict"] = [m42.verdict(e, sg, q) for e, sg, q in zip(pr.effect, pr.predicted_sign, pr.q)]
        pr.to_csv(OUT / "predictions.csv", index=False)
        log("predictions (family D):\n" + pr[["id", "cohort", "effect", "p", "q", "verdict", "n"]].to_string(index=False))
    exploratory(ct, burden, pct, L, g)
    msi_c.reset_index().to_csv(OUT / "locus_mean_contraction_msi.csv", index=False)
    burden.groupby(["cohort", "group"]).burden.describe().reset_index().to_csv(OUT / "burden_by_group.csv", index=False)


def exploratory(ct, burden, pct, L, g) -> None:
    """BAT-25/26 by group; G3/G5 with c partialled on the tumour's burden; normals' share of reads
    shorter than the reference length, MSI patients against the rest."""
    mk = ct.merge(L[["symbol", "kind"]], left_on="locus", right_index=True)
    mk = mk[mk.kind.eq("marker")]
    mk.groupby(["symbol", "cohort", "group"]).c.describe().reset_index().to_csv(OUT / "markers_by_group.csv", index=False)
    log(f"markers:\n{mk.groupby(['symbol', 'group']).c.agg(['size', 'median']).round(3)}")
    se = L[L.kind.eq("SE")]
    d = ct[ct.group.eq("MSI")].merge(L[["feature_id"]], left_on="locus", right_index=True)
    d = d.merge(pct, on=["case_id", "feature_id", "cohort"]).merge(burden[["case_id", "burden"]], on="case_id")
    rows = []
    for name, lo, hi in [("G3", 10, 19), ("G5", 20, 999)]:
        loci = se[(se.ref_len >= lo) & (se.ref_len <= hi)].index
        rho = []
        for k, x in d[d.locus.isin(loci)].groupby("locus"):
            if len(x) >= MIN_TUMOURS and x.c.nunique() > 1 and x.pct.nunique() > 1:
                rho.append(m42.partial_spearman(x.c.to_numpy(), x.pct.to_numpy(), [x.burden.to_numpy()])[0])
        rho = np.asarray(rho)
        if len(rho):
            rows.append({"what": f"{name} with c partialled on burden", "loci": len(rho),
                         "median_rho": float(np.median(rho)), "share_negative": float((rho < 0).mean()),
                         "sign_p": binomtest(int((rho < 0).sum()), int((rho != 0).sum())).pvalue})
    nor = g[g.role.eq("normal") & g.lengths.map(len).ge(MIN_READS)].merge(L[["ref_len", "kind"]], left_on="locus",
                                                                          right_index=True)
    nor = nor[nor.kind.isin(["SE", "A3SS"])]
    if len(nor):
        nor["short"] = [float((np.asarray(v) < n).mean()) for v, n in zip(nor.lengths, nor.ref_len)]
        per = nor.groupby(["case_id", "group"]).short.mean().reset_index()
        x, y = per.loc[per.group.eq("MSI"), "short"], per.loc[per.group.eq("rest"), "short"]
        if len(x) and len(y):
            r = mannwhitneyu(x, y)
            rows.append({"what": "normals: share of reads shorter than the reference, MSI minus rest patients",
                         "loci": int(nor.locus.nunique()), "median_rho": float(x.median() - y.median()),
                         "share_negative": np.nan, "sign_p": float(r.pvalue)})
    pd.DataFrame(rows).to_csv(OUT / "exploratory.csv", index=False, float_format="%.4g")
    if rows:
        log("exploratory:\n" + pd.DataFrame(rows).to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["genotype", "test"], required=True)
    a = ap.parse_args()
    genotype() if a.only == "genotype" else test()


if __name__ == "__main__":
    main()
