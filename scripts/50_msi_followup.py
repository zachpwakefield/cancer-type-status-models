"""Step 50 (D-61..D-64): the MSI acceptor repeats followed up, from local data and GTEx sQTLs.

  geometry  label-blind: for every D-52 R exon (acceptor T-run >= 10), where the run lies, the
            branch-point A nearest the 3'ss upstream of it (yUnAy motif), the run's genomic
            interval and an antisense-Alu k-mer check, from the GRCh38 primary assembly
  load      D-61 (DL1-DL6): MC3 homopolymer-deletion load per tumour against the D-52 score
  nmd       D-62 (NM1-NM4): frameshifting coding repeat exons and their genes' mRNA in MSI
  branch    D-63 (BP1-BP4): the MSI shift against the branch-point distance, beyond run length
  germline  D-64 (GL1-GL3): GTEx v8 germline poly-T indels at the same exons

Positions count from the exon in transcript orientation: 1 is the last intron base before the
exon (-1 in D-52), larger numbers lie further upstream.

Outputs: data.nosync/ctsm_v1/msi/ (per-tumour and per-exon tables, kept local)
         results/subtype/msi_cis/followup/ (aggregates)
Usage: python scripts/50_msi_followup.py --only geometry|load|nmd|branch|germline
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
from _common import REPO, log
from scipy.stats import binomtest, mannwhitneyu, rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

from ctsm.biology import gex_by_symbol
from ctsm.config import load_config

_spec = importlib.util.spec_from_file_location("msi42", Path(__file__).with_name("42_msi_repeats.py"))
m42 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m42)

DATA, LOCAL = m42.DATA, m42.LOCAL
EXT = DATA / "external"
OUT = REPO / "results" / "subtype" / "msi_cis" / "followup"
COHORTS = m42.COHORTS
KEEP_FLAGS = {"PASS", "wga", "native_wga_mix"}
SEED, BOOT, PERM = 0, 2000, 1000
UP = 450                                   # intron bases read upstream of the 3'ss (run + 300 for the Alu check)
RUN_BINS = [(10, 11), (12, 13), (14, 15), (16, 19), (20, 24), (25, 999)]
ALU = ("TGAGCCACCGCGCCCGGCC", "GGCCGGGCGCGGTGGCTCA")
TISSUES = {"STAD_TCGA": "Stomach", "COAD_TCGA": "Colon_Transverse", "UCEC_TCGA": "Uterus"}
# members of https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_sQTL.tar,
# fetched 2026-10-13 into data.nosync/ctsm_v1/external/gtex_sqtl/ (bytes, md5)
GTEX_PINS = {"Colon_Transverse": (12261452, "e4b51c5561d8eab75e805216d93cd8af"),
             "Stomach": (8506496, "edea3e28b3e45ef5a63052f0b6a09d1b"),
             "Uterus": (3084115, "d51da0792652bbdee968060be13a177a")}
COMP = str.maketrans("ACGTN", "TGCAN")


def effects() -> pd.DataFrame:
    """D-52 per-exon MSI effects with the flanking-exon coordinates."""
    e = pd.read_parquet(LOCAL / "exon_effects.parquet")
    s = pd.read_parquet(LOCAL / "se_repeats.parquet", columns=["feature_id", "left_end", "right_start"])
    return e.merge(s, on="feature_id", how="left")


def add_verdicts(pred: list[dict]) -> pd.DataFrame:
    pr = pd.DataFrame(pred).reindex(columns=["id", "what", "effect", "predicted_sign", "p", "n"])
    ok = pr.p.notna()
    pr["q"] = np.nan
    if ok.any():
        pr.loc[ok, "q"] = multipletests(pr.loc[ok, "p"], method="fdr_bh")[1]
    pr["verdict"] = [m42.verdict(e, s, q) if ok_ else "not testable"
                     for e, s, q, ok_ in zip(pr.effect, pr.predicted_sign, pr.q, ok)]
    return pr


def save_predictions(pr: pd.DataFrame, sub: str) -> None:
    d = OUT / sub
    d.mkdir(parents=True, exist_ok=True)
    pr.to_csv(d / "predictions.csv", index=False)
    log(f"{sub} predictions:\n" + pr[["id", "effect", "p", "q", "verdict", "n"]].to_string(index=False))


# ---------------------------------------------------------------- geometry (label-blind)

def run_geometry(rev: str) -> dict:
    """rev: intron bases in position order (index 0 = position 1). Longest T-run overlapping 1..50
    (ties: nearest the 3'ss), and the nearest yUnAy branch-point A upstream of it within 100."""
    best = None
    for m in re.finditer(r"T+", rev):
        s, e = m.start() + 1, m.end()
        if s <= 50 and (best is None or e - s + 1 > best[1] - best[0] + 1):
            best = (s, e)
    if best is None:
        return {}
    s, e = best
    bp = np.nan
    for p in range(e + 1, min(100, len(rev) - 3) + 1):
        # position p is the A; p+1 lies 1 nt 5' of it, p-1 1 nt 3' of it (toward the 3'ss)
        if rev[p - 1] == "A" and rev[p + 1] == "T" and rev[p + 2] in "CT" and rev[p - 2] in "CT":
            bp = p
            break
    alu = any(k in rev[e:e + 300][::-1] for k in ALU)       # 5'->3' sense sequence upstream of the run
    return {"run_near": s, "run_far": e, "run_len": e - s + 1, "d_run": e, "d_bp": bp, "alu_kmer": alu}


def geometry() -> None:
    e = effects()
    ex = e[e.acc_T >= 10].drop_duplicates(["chrom", "strand", "es", "ee", "left_end"])
    rows = []
    for chrom, seq in m42.read_genome(set(ex.chrom)):
        for r in ex[ex.chrom.eq(chrom)].itertuples(index=False):
            if r.strand == "+":
                lo = max(int(r.left_end), int(r.es) - UP)
                up = seq[lo:int(r.es)]                       # 5'->3' sense, ends at position 1
            else:
                hi = min(int(r.right_start), int(r.ee) + UP)
                up = seq[int(r.ee):hi][::-1].translate(COMP)
            g = run_geometry(up[::-1])
            if not g:
                continue
            if r.strand == "+":                                # genomic 1-based interval of the run
                g["g_lo"], g["g_hi"] = int(r.es) - g["run_far"] + 1, int(r.es) - g["run_near"] + 1
            else:
                g["g_lo"], g["g_hi"] = int(r.ee) + g["run_near"], int(r.ee) + g["run_far"]
            rows.append({"chrom": chrom, "strand": r.strand, "es": r.es, "ee": r.ee, "left_end": r.left_end,
                         "acc_T": r.acc_T, **g})
    g = pd.DataFrame(rows)
    g.to_parquet(LOCAL / "run_geometry.parquet")
    log(f"geometry: {len(g)} exon-flank pairs, {len(g.drop_duplicates(['chrom', 'strand', 'es', 'ee']))} exons; run length agrees with D-52 acc_T for {(g.run_len == g.acc_T).mean():.3f}; "
        f"branch point found for {g.d_bp.notna().mean():.3f}; median d_bp {g.d_bp.median():.0f}")


# ---------------------------------------------------------------- D-61 deletion load

def purity_by_patient() -> pd.Series:
    t = pd.read_csv(EXT / "gdc_pancanatlas" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
    t = t[t["call status"].eq("called")].sort_values("array")
    t["case_id"] = t["array"].str[:12]
    return t.drop_duplicates("case_id").set_index("case_id")["purity"]


def homopolymer_deletions(samples: set[str]) -> tuple[pd.DataFrame, set[str]]:
    """MC3 homopolymer deletions of the RNA tumours (15-character sample barcodes; the first aliquot
    in sort order when a sample has several)."""
    cols = ["Tumor_Sample_Barcode", "Variant_Type", "Reference_Allele", "CONTEXT", "FILTER"]
    m = pd.read_csv(EXT / "gdc_pancanatlas" / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=cols, dtype=str)
    m["sample15"] = m.Tumor_Sample_Barcode.str[:15]
    m = m[m.sample15.isin(samples)]
    first = m.sort_values("Tumor_Sample_Barcode").drop_duplicates("sample15").Tumor_Sample_Barcode
    m = m[m.Tumor_Sample_Barcode.isin(first)]
    m["case_id"] = m.sample15.str[:12]
    present = set(m.case_id)
    m = m[m.FILTER.str.split(",").map(lambda f: set(f) <= KEEP_FLAGS)]
    d = m[m.Variant_Type.eq("DEL") & m.Reference_Allele.str.fullmatch(r"A+|T+")].copy()
    d["L"] = d.Reference_Allele.str.len()
    ctx_ok = [c[6:6 + n] == r for c, n, r in zip(d.CONTEXT, d.L, d.Reference_Allele)]
    log(f"MC3: {len(present)} cohort patients; {len(d):,} A/T deletions, CONTEXT agrees for {np.mean(ctx_ok):.4f}")
    d = d[ctx_ok]
    d["run"] = [context_run(c, r) for c, r in zip(d.CONTEXT, d.Reference_Allele)]
    return d[d.run >= 5], present


def context_run(context: str, deleted: str) -> int:
    """Run of the deleted base starting at a left-aligned MC3 deletion: the deleted bases plus the
    same base continuing into the right flank (CONTEXT = 6 bases, the deleted bases, 5 bases)."""
    right = context[6 + len(deleted):]
    return len(deleted) + len(right) - len(right.lstrip(deleted[0]))


def load() -> None:
    sc = pd.read_parquet(LOCAL / "sample_scores.parquet")
    sc = sc[sc.group.isin(["MSI", "rest"])].copy()
    sid = pd.concat([pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{c}.parquet", columns=["file_id", "sample_id"])
                     for c in COHORTS]).drop_duplicates("file_id").set_index("file_id").sample_id
    sc["sample15"] = sc.file_id.map(sid).str[:15]
    d, present = homopolymer_deletions(set(sc.sample15.dropna()))
    sc = sc[sc.case_id.isin(present)].copy()
    n = d.groupby("case_id").size()
    nl = d[d.L >= 2].groupby("case_id").size()
    sc["n_hp"] = sc.case_id.map(n).fillna(0)
    sc["n_long"] = sc.case_id.map(nl).fillna(0)
    sc["n_short"] = sc.n_hp - sc.n_long
    sc["load"], sc["long_load"], sc["short_load"] = np.log1p(sc.n_hp), np.log1p(sc.n_long), np.log1p(sc.n_short)
    sc["long_share"] = np.where(sc.n_hp >= 5, sc.n_long / sc.n_hp.clip(lower=1), np.nan)
    sc["purity"] = sc.case_id.map(purity_by_patient())
    sc.to_parquet(LOCAL / "deletion_load.parquet")
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    for pid, c in [("DL1", "STAD_TCGA"), ("DL2", "COAD_TCGA"), ("DL3", "UCEC_TCGA")]:
        s = sc[sc.cohort.eq(c) & sc.group.eq("MSI")]
        r, p, k = m42.partial_spearman(s.score_acc.to_numpy(), s.load.to_numpy(), [s.purity.to_numpy()])
        add(pid, r, -1, p, f"{k} MSI tumours", f"{c}: partial Spearman(score, homopolymer-deletion load | purity)")
    msi = sc[sc.group.eq("MSI")]
    u = msi[msi.cohort.eq("UCEC_TCGA")].long_share.dropna()
    for pid, c in [("DL4", "STAD_TCGA"), ("DL5", "COAD_TCGA")]:
        o = msi[msi.cohort.eq(c)].long_share.dropna()
        res = mannwhitneyu(o, u, alternative="two-sided")
        add(pid, res.statistic / (len(o) * len(u)) - 0.5, 1, float(res.pvalue),
            f"UCEC {len(u)}, {c.split('_')[0]} {len(o)}; median long share {u.median():.3f} vs {o.median():.3f}",
            f"AUROC(long share higher in {c} than UCEC_TCGA MSI) - 0.5")
    w = msi.copy()
    for col in ["score_acc", "long_load", "short_load", "purity"]:
        w[col + "_r"] = w.groupby("cohort")[col].rank(pct=True)
    r, p, k = m42.partial_spearman(w.score_acc_r.to_numpy(), w.long_load_r.to_numpy(),
                                   [w.short_load_r.to_numpy(), w.purity_r.to_numpy()])
    add("DL6", r, -1, p, f"{k} MSI tumours", "MSI, ranks within cohort: partial Spearman(score, long load | short load, purity)")
    pr = add_verdicts(pred)
    save_predictions(pr, "load")
    ex = []
    for c in COHORTS:
        for g in ["MSI", "rest"]:
            s = sc[sc.cohort.eq(c) & sc.group.eq(g)]
            r, p, k = m42.partial_spearman(s.score_acc.to_numpy(), s.load.to_numpy(), [s.purity.to_numpy()])
            r0 = spearmanr(s.score_acc, s.load).statistic
            ex.append({"cohort": c, "group": g, "n": len(s), "median_n_hp": s.n_hp.median(),
                       "q1_n_hp": s.n_hp.quantile(0.25), "q3_n_hp": s.n_hp.quantile(0.75),
                       "median_long_share": s.long_share.median(), "partial_rho_score_load": r, "partial_p": p,
                       "rho_unadjusted": r0})
    ex = pd.DataFrame(ex)
    ex.loc[ex.n < 5, ["median_n_hp", "q1_n_hp", "q3_n_hp", "median_long_share"]] = np.nan
    ex.to_csv(OUT / "load" / "by_cohort_group.csv", index=False, float_format="%.4g")
    lens = d.merge(sc[["case_id", "cohort", "group"]], on="case_id")
    lens["L_class"] = np.where(lens.L >= 4, ">=4", lens.L.astype(str))
    lt = pd.crosstab([lens.cohort, lens.group], lens.L_class, normalize="index").reset_index()
    lt.to_csv(OUT / "load" / "deletion_length_shares.csv", index=False, float_format="%.4g")
    log(f"by cohort and group:\n{ex.round(3).to_string(index=False)}\nlength shares:\n{lt.round(3).to_string(index=False)}")


# ---------------------------------------------------------------- D-62 frame and NMD

def canonical_coding_exons() -> pd.DataFrame:
    """Exons of Ensembl_canonical protein_coding transcripts: number, count, fully-CDS flag (1-based)."""
    gtf = load_config().gencode_gtf
    ex, cds = [], set()
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] not in ("exon", "CDS") or 'tag "Ensembl_canonical"' not in f[8] \
                    or 'transcript_type "protein_coding"' not in f[8]:
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            k = int(re.search(r"exon_number (\d+)", f[8]).group(1))
            if f[2] == "CDS":
                cds.add((t, k, int(f[3]), int(f[4])))
            else:
                ex.append((t, k, f[0], int(f[3]), int(f[4]), f[6]))
    e = pd.DataFrame(ex, columns=["transcript", "k", "chrom", "start", "end", "strand"])
    e["n_exons"] = e.groupby("transcript").k.transform("max")
    e["full_cds"] = [(t, k, s, x) in cds for t, k, s, x in zip(e.transcript, e.k, e.start, e.end)]
    e["internal"] = (e.k >= 2) & (e.k <= e.n_exons - 2)
    return e


def nmd() -> None:
    cfg = load_config()
    e = effects()
    cx = canonical_coding_exons()
    cx = cx[cx.full_cds & cx.internal].drop_duplicates(["chrom", "start", "end", "strand"])
    r = e[e.acc_T >= 10].copy()
    r["start"], r["end"] = r.es + 1, r.ee
    r = r.merge(cx[["chrom", "start", "end", "strand"]], on=["chrom", "start", "end", "strand"])
    r["frameshift"] = r.exon_len % 3 != 0
    r = r.sort_values("p").drop_duplicates(["cohort", "symbol"])
    log(f"coding internal R exons, one per gene and cohort: {len(r)} ({r.frameshift.mean():.2f} frameshifting)")
    fc = []
    for c in COHORTS:
        tab = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{c}.parquet")
        grp = np.where(tab.label.eq("MSI"), "MSI", np.where(tab.label.eq("POLE"), "POLE", "rest"))
        g = gex_by_symbol(cfg, tab.file_id, sorted(set(r.loc[r.cohort.eq(c), "symbol"])))
        d = g[grp == "MSI"].mean() - g[grp == "rest"].mean()
        fc.append(pd.DataFrame({"cohort": c, "symbol": d.index, "log2fc": d.to_numpy()}))
    r = r.merge(pd.concat(fc), on=["cohort", "symbol"], how="left")
    r.to_parquet(LOCAL / "frame_nmd_exons.parquet")
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    for pid, c in [("NM1", "STAD_TCGA"), ("NM2", "COAD_TCGA"), ("NM3", "UCEC_TCGA")]:
        s = r[r.cohort.eq(c) & (r.delta_adj < 0) & (r.q < 0.05) & r.log2fc.notna()]
        fs, inf = s[s.frameshift].log2fc, s[~s.frameshift].log2fc
        if min(len(fs), len(inf)) < 5:
            pred.append({"id": pid, "what": f"{c}: not testable", "n": f"frameshifting {len(fs)}, in-frame {len(inf)}"})
            continue
        res = mannwhitneyu(inf, fs, alternative="two-sided")
        add(pid, res.statistic / (len(fs) * len(inf)) - 0.5, 1, float(res.pvalue),
            f"frameshifting {len(fs)} (median {fs.median():+.3f}), in-frame {len(inf)} (median {inf.median():+.3f})",
            f"{c}, exons skipped in MSI: AUROC(gene log2FC lower when frameshifting) - 0.5")
    a = r[r.log2fc.notna() & r.delta_adj.notna()].reset_index(drop=True)
    rho = lambda s: spearmanr(s.delta_adj, s.log2fc).statistic
    obs = rho(a[a.frameshift]) - rho(a[~a.frameshift])
    rng = np.random.default_rng(SEED)
    bs = []
    for _ in range(BOOT):
        b = a.iloc[rng.integers(0, len(a), len(a))]
        bs.append(rho(b[b.frameshift]) - rho(b[~b.frameshift]))
    bs = np.asarray(bs)
    p = max(2 * min((bs <= 0).mean(), (bs >= 0).mean()), 1 / BOOT)
    add("NM4", float(obs), 1, float(p), f"frameshifting {int(a.frameshift.sum())} (rho {rho(a[a.frameshift]):.3f}), "
        f"in-frame {int((~a.frameshift).sum())} (rho {rho(a[~a.frameshift]):.3f})",
        "all cohorts: Spearman(delta_adj, gene log2FC), frameshifting minus in-frame")
    save_predictions(add_verdicts(pred), "nmd")
    sh = pd.read_csv(REPO / "results" / "subtype" / "msi_cis" / "shared_repeat_exons.csv")
    sh = sh[sh.skipped_all_three]
    p_ = sh.feature_id.str.split("#", expand=True)
    loc = p_[1].str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
    sh = sh.assign(chrom=loc[0], start=loc[1].astype(int) + 1, end=loc[2].astype(int), strand=p_[2].str[0])
    sh = sh.merge(cx[["chrom", "start", "end", "strand"]].assign(coding_internal=True),
                  on=["chrom", "start", "end", "strand"], how="left")
    sh["coding_internal"] = sh.coding_internal.notna()
    sh["frameshift"] = sh.exon_len % 3 != 0
    fcm = r.groupby("symbol").log2fc.mean()
    sh["mean_gene_log2fc"] = sh.symbol.map(fcm)
    sh = sh.drop_duplicates(["chrom", "start", "end"])
    cols = ["symbol", "chrom", "start", "end", "strand", "exon_len", "acc_T", "frameshift", "coding_internal",
            "mean_delta_adj", "mean_gene_log2fc"]
    sh[cols].sort_values("mean_delta_adj").to_csv(OUT / "nmd" / "shared_skipped_exons.csv", index=False,
                                                   float_format="%.4g")
    k = sh.coding_internal & sh.frameshift
    log(f"shared skipped repeat exons: {len(sh)}; coding internal {int(sh.coding_internal.sum())}; "
        f"frameshifting coding {int(k.sum())}: {', '.join(sh.loc[k, 'symbol'].head(25))}")


# ---------------------------------------------------------------- D-63 branch-point geometry

def adjusted_rho(y, x, bins, covs, rng, perm: int = PERM) -> tuple[float, float, int]:
    """Spearman of rank residuals (on run-length-bin dummies and covariate ranks); p by permuting x
    within bins."""
    ok = np.isfinite(y) & np.isfinite(x)
    for c in covs:
        ok &= np.isfinite(c)
    y, x, bins = y[ok], x[ok], bins[ok]
    covs = [c[ok] for c in covs]
    D = pd.get_dummies(bins).to_numpy(float)
    Z = np.column_stack([D] + [rankdata(c) for c in covs])
    Q, _ = np.linalg.qr(Z)
    res = lambda v: v - Q @ (Q.T @ v)
    ry, rx = res(rankdata(y)), rankdata(x)
    r = float(np.corrcoef(ry, res(rx))[0, 1])
    groups = [np.flatnonzero(bins == b) for b in np.unique(bins)]
    null = np.empty(perm)
    for k in range(perm):
        xp = rx.copy()
        for g in groups:
            xp[g] = rx[rng.permutation(g)]
        null[k] = np.corrcoef(ry, res(xp))[0, 1]
    p = (1 + (np.abs(null) >= abs(r)).sum()) / (1 + perm)
    return r, float(p), int(ok.sum())


def run_bin(t: np.ndarray) -> np.ndarray:
    out = np.full(len(t), "", dtype=object)
    for lo, hi in RUN_BINS:
        out[(t >= lo) & (t <= hi)] = f"{lo}-{hi}" if hi < 999 else f">={lo}"
    return out


def branch() -> None:
    g = pd.read_parquet(LOCAL / "run_geometry.parquet")
    e = effects()
    e = e[e.acc_T >= 10].merge(g[["chrom", "strand", "es", "ee", "left_end", "d_bp", "d_run", "alu_kmer"]],
                               on=["chrom", "strand", "es", "ee", "left_end"], how="inner")
    rng = np.random.default_rng(SEED)
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    for pid, c in [("BP1", "STAD_TCGA"), ("BP2", "COAD_TCGA"), ("BP3", "UCEC_TCGA")]:
        s = e[e.cohort.eq(c)]
        r, p, n = adjusted_rho(s.delta_adj.to_numpy(), s.d_bp.to_numpy(), run_bin(s.acc_T.to_numpy()),
                               [s.mean_psi.to_numpy()], rng)
        add(pid, r, 1, p, f"{n} exons", f"{c}: Spearman(delta_adj, d_bp | run-length bin, mean PSI)")
    m = e.groupby("feature_id").agg(delta_adj=("delta_adj", "mean"), mean_psi=("mean_psi", "mean"),
                                    acc_T=("acc_T", "first"), d_bp=("d_bp", "first"), d_run=("d_run", "first"),
                                    alu_kmer=("alu_kmer", "first"), n=("cohort", "size"))
    m2 = m[m.n >= 2]
    r, p, n = adjusted_rho(m2.delta_adj.to_numpy(), m2.d_bp.to_numpy(), run_bin(m2.acc_T.to_numpy()),
                           [m2.mean_psi.to_numpy()], rng)
    add("BP4", r, 1, p, f"{n} exons in >= 2 cohorts", "mean delta_adj over cohorts: Spearman(., d_bp | run-length bin, mean PSI)")
    save_predictions(add_verdicts(pred), "branch")
    ex = []
    for c in COHORTS + ["mean of >= 2 cohorts"]:
        s = e[e.cohort.eq(c)] if c in COHORTS else m2.reset_index()
        r, p, n = adjusted_rho(s.delta_adj.to_numpy(), s.d_run.to_numpy().astype(float), run_bin(s.acc_T.to_numpy()),
                               [s.mean_psi.to_numpy()], rng, perm=200)
        ex.append({"what": "d_run in place of d_bp", "cohort": c, "rho": r, "p": p, "n": n})
        for cls, k in [("10-19", s.acc_T < 20), (">=20", s.acc_T >= 20)]:
            t = s[k & s.d_bp.notna()]
            if len(t) >= 15:
                tert = pd.qcut(t.d_bp.rank(method="first"), 3, labels=["near", "middle", "far"])
                for lab, v in t.groupby(tert, observed=True):
                    ex.append({"what": "share moving up by d_bp tertile", "cohort": c, "run_class": cls, "tertile": lab,
                               "n": len(v), "share_up": float((v.delta_adj > 0).mean()),
                               "median_d_bp": float(v.d_bp.median())})
            a = s[k]
            for alu, v in a.groupby(a.alu_kmer):
                ex.append({"what": "antisense-Alu k-mer upstream of the run", "cohort": c, "run_class": cls,
                           "alu_kmer": bool(alu), "n": len(v), "share_up": float((v.delta_adj > 0).mean()),
                           "median_delta_adj": float(v.delta_adj.median())})
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "branch" / "exploratory.csv", index=False, float_format="%.4g")
    log(f"branch exploratory:\n{ex.round(3).to_string(index=False)}")


# ---------------------------------------------------------------- D-64 germline indels

def gtex_file(tissue: str) -> Path:
    """The pinned GTEx v8 significant-pairs file of a tissue, checked by size and md5."""
    f = EXT / "gtex_sqtl" / f"{tissue}.v8.sqtl_signifpairs.txt.gz"
    size, md5 = GTEX_PINS[tissue]
    if f.stat().st_size != size or hashlib.md5(f.read_bytes()).hexdigest() != md5:
        raise SystemExit(f"{f} does not match its pin ({size} bytes, md5 {md5})")
    return f


def poly_t_indels(g: pd.DataFrame) -> pd.DataFrame:
    """GTEx significant sQTL rows whose variant is an indel of the run's sense T inside an R-exon run."""
    rows = []
    for c, tissue in TISSUES.items():
        q = pd.read_csv(gtex_file(tissue), sep="\t",
                        usecols=["variant_id", "phenotype_id", "slope", "slope_se", "pval_nominal"])
        v = q.variant_id.str.split("_", expand=True)
        q = q.assign(chrom=v[0], pos=v[1].astype(int), ref=v[2], alt=v[3])
        q = q[q.ref.str.len() != q.alt.str.len()]
        q["tissue"] = tissue
        rows.append(q)
    q = pd.concat(rows, ignore_index=True)
    deleted = np.where(q.ref.str.len() > q.alt.str.len(),
                       [r[len(a):] if r.startswith(a) else "" for r, a in zip(q.ref, q.alt)],
                       [a[len(r):] if a.startswith(r) else "" for r, a in zip(q.ref, q.alt)])
    q["indel_bases"] = deleted
    q["delta"] = np.where(q.ref.str.len() > q.alt.str.len(), -1, 1)
    out = []
    for r in g.itertuples(index=False):
        b = "T" if r.strand == "+" else "A"
        k = q[q.chrom.eq(r.chrom) & (q.pos >= r.g_lo - 1) & (q.pos <= r.g_hi)]
        k = k[k.indel_bases.str.fullmatch(f"{b}+")]
        if len(k):
            out.append(k.assign(es=r.es, ee=r.ee, strand=r.strand, left_end=r.left_end, acc_T=r.acc_T))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def coverage(g: pd.DataFrame) -> None:
    """How many runs carry any GTEx significant variant (any type) inside them or near them."""
    rows = []
    for tissue in TISSUES.values():
        q = pd.read_csv(gtex_file(tissue), sep="\t", usecols=["variant_id"]).drop_duplicates()
        v = q.variant_id.str.split("_", expand=True)
        q = q.assign(chrom=v[0], pos=v[1].astype(int), indel=v[2].str.len() != v[3].str.len())
        by = {c: (d.sort_values("pos").pos.to_numpy(), d.sort_values("pos").indel.to_numpy())
              for c, d in q.groupby("chrom")}
        for pad, where in [(0, "inside the run (from start - 1)"), (50, "within 50 nt of the run"),
                           (300, "within 300 nt of the run")]:
            n_any = n_ind = 0
            for r in g.itertuples(index=False):
                if r.chrom not in by:
                    continue
                p, ind = by[r.chrom]
                i, j = np.searchsorted(p, r.g_lo - 1 - pad), np.searchsorted(p, r.g_hi + pad, side="right")
                n_any += j > i
                n_ind += bool(ind[i:j].any())
            rows.append({"tissue": tissue, "window": where, "runs": len(g), "runs_with_any_variant": n_any,
                         "runs_with_an_indel": n_ind})
    c = pd.DataFrame(rows)
    (OUT / "germline").mkdir(parents=True, exist_ok=True)
    c.to_csv(OUT / "germline" / "variant_coverage.csv", index=False)
    log(f"GTEx variant coverage of the runs:\n{c.to_string(index=False)}")


def germline() -> None:
    g = pd.read_parquet(LOCAL / "run_geometry.parquet")
    e = effects()
    e = e[e.acc_T >= 10]
    hits = poly_t_indels(g.drop_duplicates(["chrom", "strand", "es", "ee"]))
    log(f"GTEx significant pairs with a poly-T indel in an R-exon run: {len(hits)} rows, "
        f"{hits[['chrom', 'es', 'ee']].drop_duplicates().shape[0] if len(hits) else 0} exons")
    coverage(g.drop_duplicates(["chrom", "strand", "es", "ee"]))
    rows = []
    if len(hits):
        ph = hits.phenotype_id.str.split(":", expand=True)
        hits = hits.assign(ps=ph[1].astype(int), pe=ph[2].astype(int))
        feats = e.drop_duplicates("feature_id")
        for h in hits.itertuples(index=False):
            f = feats[feats.chrom.eq(h.chrom) & feats.es.eq(h.es) & feats.ee.eq(h.ee)]
            for fr in f.itertuples(index=False):
                junc = {(int(fr.left_end), int(fr.es) + 1): 1, (int(fr.ee), int(fr.right_start) + 1): 1,
                        (int(fr.left_end), int(fr.right_start) + 1): -1}
                s = junc.get((h.ps, h.pe))
                if s is not None:
                    rows.append({"chrom": h.chrom, "es": h.es, "ee": h.ee, "acc_T": h.acc_T, "tissue": h.tissue,
                                 "variant_id": h.variant_id, "phenotype_id": h.phenotype_id, "junction_sign": s,
                                 "delta": h.delta, "z": h.slope / h.slope_se, "signed_z": h.slope / h.slope_se * h.delta * s})
    j = pd.DataFrame(rows).drop_duplicates(["variant_id", "phenotype_id", "tissue"]) if rows else pd.DataFrame()
    d = OUT / "germline"
    d.mkdir(parents=True, exist_ok=True)
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    if len(j):
        j.to_csv(d / "variants_by_exon.csv", index=False, float_format="%.4g")
        ex = j.groupby(["chrom", "es", "ee"]).agg(acc_T=("acc_T", "first"), effect=("signed_z", "median"),
                                                  rows=("signed_z", "size"), tissues=("tissue", "nunique")).reset_index()
        md = e.groupby(["chrom", "es", "ee"]).delta_adj.mean().rename("mean_delta_adj").reset_index()
        ex = ex.merge(md, on=["chrom", "es", "ee"], how="left")
    else:
        ex = pd.DataFrame(columns=["acc_T", "effect", "mean_delta_adj"])
    for pid, cls, k, sign in [("GL1", "10-19", ex.acc_T < 20, 1), ("GL2", ">=20", ex.acc_T >= 20, -1)]:
        v = ex.loc[k, "effect"].dropna()
        if len(v) < 5:
            pred.append({"id": pid, "what": f"runs {cls}: not testable", "n": f"{len(v)} exons"})
            continue
        pos = int((v > 0).sum())
        add(pid, pos / len(v) - 0.5, sign, float(binomtest(pos, len(v)).pvalue), f"{pos} of {len(v)} exons positive",
            f"runs {cls}: share of exons where a longer germline tract raises inclusion - 0.5")
    v = ex.dropna(subset=["effect", "mean_delta_adj"])
    if len(v) >= 5:
        r = spearmanr(v.effect, v.mean_delta_adj)
        add("GL3", float(r.statistic), -1, float(r.pvalue), f"{len(v)} exons",
            "Spearman(germline lengthening effect, mean MSI delta_adj)")
    else:
        pred.append({"id": "GL3", "what": "not testable", "n": f"{len(v)} exons"})
    save_predictions(add_verdicts(pred), "germline")
    if len(ex):
        ex.to_csv(d / "exons.csv", index=False, float_format="%.4g")
        log(f"germline exons:\n{ex.round(3).to_string(index=False)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["geometry", "load", "nmd", "branch", "germline"], required=True)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {"geometry": geometry, "load": load, "nmd": nmd, "branch": branch, "germline": germline}[a.only]()


if __name__ == "__main__":
    main()
