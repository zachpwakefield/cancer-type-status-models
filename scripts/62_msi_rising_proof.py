"""Step 62 (D-77): a step-wise test of the MSI rising arm, in junction reads.

The claim under test (config/followup.yaml: msi_rising_proof): in mismatch-repair-deficient (MMRd)
tumours, cassette exons whose 3' splice site carries a long acceptor poly-T tract (longest T-run
overlapping intron -1..-50 >= 20 nt) and a distant predicted branch point (>= 36 nt from the AG)
are more included than in proficient tumours of the same cancer, and the rise needs both features.
Confirmatory tests use the 30 held-out cancers and junction reads only (rMATS JC); the repository's
SE matrix is rMATS JCEC, which also counts exon-body reads as inclusion.

  annotate  label-blind, for every SE exon of the pan-cancer matrix, from GRCh38 and GENCODE v45:
            the acceptor runs and branch point of steps 42/50/51 (recomputed and checked against
            the stored tables), the longest pyrimidine run at the acceptor and the branch point
            upstream of it, the pyrimidine share between that branch point and the AG, T-runs deep
            in the intron (-201..-400) and after the donor (+10..+100), a 3' splice-site log-odds
            PWM (intron -20..exon +3; every annotated acceptor against one matched intronic AG per
            intron) and its tract-blind core (-3..+3), exon length and GC, the canonical-transcript
            consequence of step 61 and the host gene's longest coding mononucleotide run, the Alu
            context of the acceptor run, the acceptor window and the exon (UCSC RepeatMasker), and
            the D-77 classes with their label-blind sizes
  branchpoints  label-blind: BPHunter's experimentally mapped branch points (hg19), lifted to
            GRCh38 through the UCSC chain (point liftover implemented here), checked on the genome,
            and placed on each exon's acceptor intron (positions 3..100)
  counts    label-blind: JC and JCEC read counts from the count store (${CTSM_SPLICE_COUNTS}) for
            D-71's 8,464 tumours and every Solid Tissue Normal sample of the 33 cancers, for the
            exons the tests need (every class member plus a fixed random sample of 12,000 B
            exons), stored per cancer; then each class's size in the JC universe of each cancer
  sizes     the class sizes of `counts` again, from the stored counts
  strata    label-blind: per stored exon, held-out and training eligibility, baseline band and
            coverage tertile; the pre-specified stratum coarsening of every class contrast with the
            exons kept and the effective control n; the collinearity of RP5a's covariates

  counts_extra  label-blind: the counts of the ppt_far >= 15 nt members `counts` did not store (a
            pre-specified sensitivity analysis needs every member), under proof/extra/
  compute   the tests of steps 1-8 and the mechanism tests, the sensitivity analyses and the
            descriptive items of the frozen block (commit 5471b27), with every label permutation,
            exon permutation and Freedman-Lane permutation it specifies
  report    tests.csv (q within family, test verdicts), sensitivity.csv, descriptive.csv and
            verdict.json - the claim-level verdict computed by claim_verdict from the pre-registered
            rule

compute and report need --labels: 'shuffled' is the dry run - every label (MC3 MMR status, the
hypermutated flag, the deletion and SNV counts, MANTIS and MSIsensor) is moved together among the
tumours of each cancer (seed 999; tumours with and without a normal pair separately) before anything
else, and every output goes under proof/dryrun/ (results and data.nosync); 'real' also needs
--real-labels-approved, given only after Gate 2, and never writes into a dryrun folder.

Positions count from the exon in transcript orientation: 1 is the last intron base before the exon
(-1 in D-52), larger numbers lie further upstream; after the donor, 1 is the first intron base.

Outputs: data.nosync/ctsm_v1/msi/proof/ (per-exon annotation, per-sample counts, adj matrices,
         per-tumour arms and the permutation engines' arrays; kept local)
         results/subtype/msi_cis/code/proof/ (aggregates only: class sizes, PWM, QC; tests.csv,
         sensitivity.csv, descriptive.csv, verdict.json and their dryrun/ versions)
Usage: python scripts/62_msi_rising_proof.py --only annotate|branchpoints|counts|counts_extra|sizes|strata
       python scripts/62_msi_rising_proof.py --only compute|report --labels shuffled
"""
from __future__ import annotations

import argparse
import gzip
import time
import warnings
import importlib.util
import json
import re
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd
from _common import log


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m61 = _load("msi61", "61_msi_rising_consequence.py")
m60 = m61.m60
m52, m51, m50, m42 = m60.m52, m60.m51, m60.m50, m60.m42
from ctsm.config import expand, load_config, portable  # noqa: E402

LOCAL, PAN, TRAIN = m51.LOCAL, m51.PAN, m51.TRAIN
PROOF = LOCAL / "proof"
OUT = m51.OUT / "proof"
COMP = m42.COMP
SEED = 0

# ---- class definitions (config/followup.yaml: msi_rising_proof, classes) ----
LONG, FAR = 20, m51.UP_CUT                 # rising arm: run >= 20 nt, branch point >= 36 nt
PPT_MIN, PPT_MAX_HOMO, PPT_MAX_ALT = 20, 7, 10   # PPT control: >= 20 pyrimidines, no T/C run >= 8,
CTRL_MAX_T = 7                                   # no CT/TC alternation >= 10; controls: acc_T <= 7
POLY_A, POLY_A_FEASIBLE, DEEP_T, DON_T = 20, 15, 20, 20   # poly-A >= 20 is too rare (37 exons); >= 15 kept
                                           # Alu-acceptor control (RP4e): acc_T <= 7 and 1..50 overlaps an Alu
DON2, DON2_MIN_INTRON = (10, 100), 250     # donor-side window (+10..+100), downstream intron length
REACH = m42.REACH                          # bases read past a window edge so a run is measured in full
PWM_UP, PWM_DOWN = 20, 3                   # 3'ss PWM window: intron -20..-1 and exon +1..+3
CORE = slice(PWM_UP - 3, PWM_UP + 3)       # tract-masked core 3'ss score: -3..-1 and +1..+3 (Gate 1)
DECOY_SPAN, DECOY_EDGE = (100, 600), 30    # decoy AG: 100-600 nt upstream of the true exon start,
                                           # its window >= 30 nt from the donor
BASES = "ACGT"
EVEN = {f"chr{i}" for i in range(2, 23, 2)}    # PWM check: fitted without these, scored on them
ALU_NEAR = 10                              # a run "in an Alu" overlaps it or lies within 10 nt of it

# ---- external files (git-ignored; pinned) ----
EXT = m51.DATA / "external"
# UCSC hg38 RepeatMasker table, downloaded 2026-10-23 with the data user's approval from
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz (Last-Modified 2022-10-19),
# with its schema rmsk.sql; genoStart is 0-based
RMSK = EXT / "ucsc_rmsk" / "rmsk.txt.gz"
RMSK_PIN = (155_633_856, "b2e108b535550ba9e3cf83c77417380f")
# BPHunter (Zhang et al. 2022 PNAS, PMID 36306325) Dataset S1, downloaded 2026-10-23 by the data
# user from the PNAS supplement: CHROM, BP_POS (1-based, hg19), STRAND, BP_NAME, SOURCE (eBP_* =
# experimentally mapped, cBP_* = predicted); lifted to GRCh38 here with the UCSC chain below
BPHUNTER = EXT / "bphunter" / "pnas.2211194119.sd01.txt"
BPHUNTER_PIN = (35_549_965, "47da2af8d85024e417ed3b7610d5fe1b")
# UCSC hg19 -> hg38 liftOver chain (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
# hg19ToHg38.over.chain.gz), copied 2026-10-23 from the data user's existing UCSC download
CHAIN = EXT / "ucsc_chain" / "hg19ToHg38.over.chain.gz"
CHAIN_PIN = (227_698, "35887f73fe5e2231656504d1f6430900")      # sha256 5c0598e500ceb5a7...b125a1d
BP_RANGE = (3, 100)                        # experimental branch points 3..100 nt from the exon
LIFT_MIN_A = 0.98                          # validation: lifted BPs named A must read A in GRCh38

# ---- counts ----
B_SAMPLE = 12_000                          # random background exons stored (seed 0)
MIN_READS, UNIVERSE_SHARE = 10, 0.5        # measured = JC IJC + SJC >= 10; universe = measured in
                                           # >= 50% of a cancer's tumours (label-blind)
BAND_EDGES = [0, .02, .05, .10, .20, .30, .40, .50, .60, .70, .80, .90, .95, 1.0001]   # 13 baseline bands
MIN_CTRL, MAX_SPAN = 3, 3                  # stratum coarsening: >= 3 control exons wherever the first
                                           # class has exons; a merged stratum spans <= 3 bands
VIEWS = ("jc", "jcec")
N_JOBS = 4


# ================================================================ sequence helpers (unit-tested)

def longest_run(rev: str, pattern: str, max_start: int = 50) -> tuple[int, int] | None:
    """(near, far) 1-based positions of the longest run matching `pattern` that starts within
    `max_start` (so it overlaps 1..max_start); ties go to the run nearest the 3'ss, as step 50."""
    best = None
    for m in re.finditer(pattern, rev):
        s, e = m.start() + 1, m.end()
        if s > max_start:
            break
        if best is None or e - s > best[1] - best[0]:
            best = (s, e)
    return best


def max_homopolymer(s: str, base: str) -> int:
    return max((len(m.group()) for m in re.finditer(f"{base}+", s)), default=0)


def max_alternation(s: str) -> int:
    """Longest perfect two-base alternation (CTCT..., TCTC...) in s, in nt (a dinucleotide
    microsatellite of a pyrimidine tract)."""
    best = 0
    for m in re.finditer(r"(?:CT)+C?|(?:TC)+T?", s):
        best = max(best, len(m.group()))
    return best


def ppt_geometry(rev: str) -> dict:
    """Longest pure-pyrimidine run overlapping positions 1..50, the yUnAy branch point upstream of
    it (step 50's find_bp, within 100 nt) and the pyrimidine share between that branch point and
    the AG (positions 3..d_bp - 1)."""
    r = longest_run(rev, r"[CT]+")
    if r is None:
        return {"ppt_near": np.nan, "ppt_far": np.nan, "ppt_len": 0, "ppt_max_T": 0, "ppt_max_C": 0,
                "ppt_alt": 0, "d_bp_ppt": m50.find_bp(rev, 3), "pyr_bp_ag": np.nan}
    s, e = r
    seq = rev[s - 1:e]
    d = m50.find_bp(rev, e + 1)
    pyr = np.nan
    if np.isfinite(d) and d > 3:
        w = rev[2:int(d) - 1]                       # positions 3..d-1
        pyr = (w.count("C") + w.count("T")) / len(w)
    return {"ppt_near": s, "ppt_far": e, "ppt_len": e - s + 1, "ppt_max_T": max_homopolymer(seq, "T"),
            "ppt_max_C": max_homopolymer(seq, "C"), "ppt_alt": max_alternation(seq), "d_bp_ppt": d,
            "pyr_bp_ag": pyr}


def downstream_seq(chrom_seq: str, r, far: int) -> tuple[str, int]:
    """Sense-orientation downstream intron, index 0 = position +1, at most `far` bases, with the
    intron length (the exon's right flank on +, its left flank on -)."""
    if r.strand == "+":
        lo, hi = int(r.ee), int(r.right_start)
        return chrom_seq[lo:min(hi, lo + far)], hi - lo
    lo, hi = int(r.left_end), int(r.es)
    return chrom_seq[max(lo, hi - far):hi][::-1].translate(COMP), hi - lo


def upstream_rev(chrom_seq: str, r, far: int) -> str:
    """Sense upstream intron in position order (index 0 = position 1), at most `far` bases (as
    step 51's geometry)."""
    if r.strand == "+":
        return chrom_seq[max(int(r.left_end), int(r.es) - far):int(r.es)][::-1]
    return chrom_seq[int(r.ee):min(int(r.right_start), int(r.ee) + far)].translate(COMP)


def acceptor_window(chrom_seq: str, strand: str, key: int) -> str:
    """Sense bases intron -20..-1 + exon +1..+3 of the acceptor whose key is the exon's 0-based
    first base (+) or its 0-based exclusive end (-)."""
    if strand == "+":
        return chrom_seq[key - PWM_UP:key + PWM_DOWN]
    return chrom_seq[key - PWM_DOWN:key + PWM_UP][::-1].translate(COMP)


def cds_max_run(chrom_seq: str, model: dict) -> float:
    """Longest single-base run in the coding sequence of a canonical transcript model (step 61's
    canonical_models: 1-based inclusive CDS parts in transcript order), in transcript orientation:
    the host gene's longest coding mononucleotide repeat, a frameshift target in MMR-deficient
    tumours. NaN without a CDS."""
    parts = [m61._seq(chrom_seq, a, b, model["strand"]) for c in model["cds"] if c is not None for a, b in [c]]
    s = "".join(parts)
    return float(max((len(m.group()) for m in re.finditer(r"A+|C+|G+|T+", s)), default=0)) if s else np.nan


def exon_gc(chrom_seq: str, es: int, ee: int) -> float:
    s = chrom_seq[es:ee]
    return (s.count("G") + s.count("C")) / max(len(s), 1)


def run_interval(strand: str, es: int, ee: int, near: float, far: float) -> tuple[int, int] | None:
    """Genomic 0-based half-open interval of upstream-intron positions near..far (position 1 = the
    last intron base before the exon)."""
    if not (np.isfinite(near) and np.isfinite(far)):
        return None
    near, far = int(near), int(far)
    return (es - far, es - near + 1) if strand == "+" else (ee + near - 1, ee + far)


# ================================================================ RepeatMasker (label-blind)

def md5_file(path: Path) -> str:
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_alu(path: Path = RMSK, pin=RMSK_PIN) -> dict:
    """{chrom: (starts, ends, strands, max_len)} of the Alu elements (repFamily 'Alu'), sorted by
    start; the file must match its pinned size and MD5."""
    if path.stat().st_size != pin[0] or md5_file(path) != pin[1]:
        raise SystemExit(f"{path}: size or MD5 differs from the pin {pin}")
    t = pd.read_csv(path, sep="\t", header=None, usecols=[5, 6, 7, 9, 12],
                    names=["chrom", "start", "end", "strand", "family"], dtype={"chrom": str, "family": str})
    t = t[t.family.eq("Alu")].sort_values(["chrom", "start"])
    return {c: (g.start.to_numpy(np.int64), g.end.to_numpy(np.int64), g.strand.to_numpy(), int((g.end - g.start).max()))
            for c, g in t.groupby("chrom")}


def alu_hit(alu: dict, chrom: str, lo: int, hi: int, pad: int = 0) -> tuple[bool, str]:
    """Whether any Alu overlaps [lo - pad, hi + pad), and the strand of the one overlapping most
    (ties: the first)."""
    if chrom not in alu:
        return False, ""
    s, e, st, mx = alu[chrom]
    a, b = lo - pad, hi + pad
    i0, i1 = np.searchsorted(s, a - mx), np.searchsorted(s, b)
    k = np.arange(i0, i1)
    k = k[(s[k] < b) & (e[k] > a)]
    if not len(k):
        return False, ""
    ov = np.minimum(e[k], b) - np.maximum(s[k], a)
    return True, str(st[k[int(np.argmax(ov))]])


# ================================================================ 3' splice-site PWM (label-blind)

def gtf_introns(gtf: Path) -> pd.DataFrame:
    """Unique introns of every GENCODE transcript: chrom, strand, i0, i1 (0-based half-open,
    genomic) and the acceptor key (0-based first exon base on +, i0 on -)."""
    tx: dict[str, tuple] = {}
    pat_t, pat_k = re.compile(r'transcript_id "([^"]+)"'), re.compile(r"exon_number (\d+)")
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "exon":
                continue
            t = pat_t.search(f[8]).group(1)
            tx.setdefault(t, (f[0], f[6], []))[2].append((int(pat_k.search(f[8]).group(1)), int(f[3]), int(f[4])))
    rows = []
    for chrom, strand, ex in tx.values():
        ex.sort()
        for (_, s1, e1), (_, s2, e2) in zip(ex, ex[1:]):
            if strand == "+":
                rows.append((chrom, strand, e1, s2 - 1))
            else:
                rows.append((chrom, strand, e2, s1 - 1))
    d = pd.DataFrame(rows, columns=["chrom", "strand", "i0", "i1"]).drop_duplicates()
    d = d[d.i1 > d.i0].sort_values(["chrom", "strand", "i0", "i1"]).reset_index(drop=True)
    d["key"] = np.where(d.strand.eq("+"), d.i1, d.i0)
    return d


def decoy(sense_tail: str, n: int, rng) -> int | None:
    """Index j (in the intron's sense coordinates, 0 = first intron base) of one AG chosen uniformly
    among those whose following base j lies DECOY_SPAN nt upstream of the true exon start (n) and
    whose window starts >= DECOY_EDGE nt into the intron. sense_tail holds the last len(sense_tail)
    intron bases."""
    off = n - len(sense_tail)
    cands = [off + m.start() + 2 for m in re.finditer(r"(?=AG)", sense_tail)]
    cands = [j for j in cands if DECOY_SPAN[0] <= n - j <= DECOY_SPAN[1] and j - PWM_UP >= DECOY_EDGE]
    return cands[int(rng.integers(len(cands)))] if cands else None


def pwm_counts(windows: list[str]) -> np.ndarray:
    ok = [w for w in windows if len(w) == PWM_UP + PWM_DOWN and set(w) <= set(BASES)]
    arr = np.frombuffer("".join(ok).encode(), dtype="S1").reshape(-1, PWM_UP + PWM_DOWN)
    return np.column_stack([(arr == b.encode()).sum(0) for b in BASES]).astype(float)


def pwm_fit(pos: list[str], neg: list[str]) -> np.ndarray:
    """Position x base log2-odds of acceptor against decoy, pseudocount 1 per base."""
    cp, cn = pwm_counts(pos), pwm_counts(neg)
    return np.log2((cp + 1) / (cp.sum(1, keepdims=True) + 4)) - np.log2((cn + 1) / (cn.sum(1, keepdims=True) + 4))


def pwm_score(windows, W: np.ndarray) -> np.ndarray:
    idx = {b: i for i, b in enumerate(BASES)}
    out = np.full(len(windows), np.nan)
    for k, w in enumerate(windows):
        if isinstance(w, str) and len(w) == W.shape[0] and set(w) <= set(BASES):
            out[k] = sum(W[i, idx[b]] for i, b in enumerate(w))
    return out


# ================================================================ annotate (label-blind)

def classes(a: pd.DataFrame) -> pd.DataFrame:
    """The D-77 classes. Main classes are exclusive (step 52's R_up/R_down/B_far/B_near, with
    long_near the other runs >= 20 nt); control flags are not."""
    at = a.acc_T.to_numpy(float)
    db = np.nan_to_num(a.d_bp.to_numpy(float), nan=-1)
    dp = np.nan_to_num(a.d_bp_ppt.to_numpy(float), nan=-1)
    st = a.struct.to_numpy(bool)
    a["cls"] = np.select([(at >= LONG) & (db >= FAR), (at >= 10) & (at < LONG), at >= LONG,
                          (at <= 5) & (db >= FAR), at <= 5], ["R_up", "R_down", "long_near", "B_far", "B_near"], "mid")
    no_t = at <= CTRL_MAX_T
    a["ctl_polyA"] = st & no_t & (a.acc_A.to_numpy(float) >= POLY_A)
    a["ctl_polyA15"] = st & no_t & (a.acc_A.to_numpy(float) >= POLY_A_FEASIBLE)
    a["ctl_ppt_far"] = (st & no_t & (a.acc_C.to_numpy(float) <= PPT_MAX_HOMO) & (a.ppt_len.to_numpy(float) >= PPT_MIN)
                        & (a.ppt_max_T.to_numpy(float) <= PPT_MAX_HOMO) & (a.ppt_max_C.to_numpy(float) <= PPT_MAX_HOMO)
                        & (a.ppt_alt.to_numpy(float) < PPT_MAX_ALT) & (dp >= FAR))
    a["ctl_deep"] = st & no_t & (a.deep_T.to_numpy(float) >= DEEP_T)
    a["ctl_don"] = st & no_t & (a.don2_T.to_numpy(float) >= DON_T)
    a["ctl_deep_far"] = a.ctl_deep & (dp >= FAR)
    a["ctl_don_far"] = a.ctl_don & (dp >= FAR)
    a["shield"] = st & a.sh_start.notna().to_numpy() & (at < 10)          # step 52's shield class
    a["ctl_alu"] = st & no_t & a.acc_alu.eq(True).to_numpy()              # RP4e (Gate 1)
    if "acc_alu_antisense" in a:
        a["ctl_alu_antisense"] = a.ctl_alu & a.acc_alu_antisense.eq(True).to_numpy()
    return a


def annotate() -> None:
    cfg = load_config()
    f = m42.parse_features(pd.read_parquet(PAN / "se" / "features.parquet"))
    introns = gtf_introns(cfg.gencode_gtf)
    acc_keys = set(zip(introns.chrom, introns.strand, introns.key))
    models = m61.canonical_models(cfg.gencode_gtf)
    log(f"GENCODE v45: {len(introns):,} unique introns, {len(acc_keys):,} acceptors, {len(models):,} canonical models")
    alu = read_alu()
    log(f"RepeatMasker: {sum(len(v[0]) for v in alu.values()):,} Alu elements")
    rg = pd.read_parquet(LOCAL / "reversal_geometry.parquet").set_index("feature_id")   # step 52: AGEZ, shield
    f = f.join(rg[["agez", "agez_censored", "sh_start", "sh_len", "sh_bp"]], on="feature_id")
    rng = np.random.default_rng(SEED)
    rows, pos_w, neg_w, cds_run = [], [], [], {}
    genes_by_chrom = {}
    for g, m in models.items():
        genes_by_chrom.setdefault(m["chrom"], []).append(g)
    for chrom, seq in m42.read_genome(set(f.chrom) | set(introns.chrom)):
        # PWM training windows: every annotated acceptor, one decoy AG per intron
        ic = introns[introns.chrom.eq(chrom)]
        for k in sorted(set(zip(ic.strand, ic.key))):
            pos_w.append((chrom, acceptor_window(seq, k[0], int(k[1]))))
        for r in ic.itertuples(index=False):
            n, far = r.i1 - r.i0, DECOY_SPAN[1] + PWM_UP + 10
            if r.strand == "+":
                tail = seq[max(r.i0, r.i1 - far):r.i1]
            else:
                tail = seq[r.i0:min(r.i1, r.i0 + far)][::-1].translate(COMP)
            j = decoy(tail, n, rng)
            if j is None:
                continue
            key = r.i0 + j if r.strand == "+" else r.i1 - j
            if (chrom, r.strand, key) in acc_keys:
                continue
            neg_w.append((chrom, acceptor_window(seq, r.strand, key)))
        # the SE exons of the pan-cancer matrix
        for r in f[f.chrom.eq(chrom)].itertuples(index=False):
            rec = {"feature_id": r.feature_id}
            for side, (lo, hi) in m42.WINDOWS.items():          # step 42, unchanged
                s, n = m42.window_seq(seq, r, side)
                rec[f"{side}_intron_len"] = n
                ok = n >= m42.MIN_INTRON[side]
                runs = m42.max_runs(s, lo, hi) if ok else {b: np.nan for b in "ACGT"}
                for b in "ACGT":
                    rec[f"{side}_{b}"] = runs[b]
                if side == "acc":
                    w = s[:hi] if ok else ""
                    rec["acc_pyr_frac"] = (w.count("C") + w.count("T")) / len(w) if w else np.nan
                    rec["acc_ag"] = s[:2][::-1] if ok else ""
            rev = upstream_rev(seq, r, m50.UP)                  # step 51's 450-nt read
            g = m50.run_geometry(rev)
            rec.update({"run_near": g.get("run_near", np.nan), "run_far": g.get("run_far", np.nan),
                        "run_len": g.get("run_len", np.nan), "d_bp_T": g.get("d_bp", np.nan),
                        "d_bp_11": m50.find_bp(rev, 11)})
            rec.update(ppt_geometry(rev))
            # Alu context of the acceptor T-run, the acceptor window -1..-50 and the exon body
            iv = run_interval(r.strand, int(r.es), int(r.ee), rec["run_near"], rec["run_far"])
            hit, ast = alu_hit(alu, chrom, *iv, pad=ALU_NEAR) if iv else (False, "")
            rec["run_alu"], rec["run_alu_antisense"] = hit, hit and ast != r.strand
            hit, ast = alu_hit(alu, chrom, *run_interval(r.strand, int(r.es), int(r.ee), 1, 50))
            rec["acc_alu"], rec["acc_alu_antisense"] = hit, hit and ast != r.strand
            rec["exon_alu"] = alu_hit(alu, chrom, int(r.es), int(r.ee))[0]
            iv = run_interval(r.strand, int(r.es), int(r.ee), r.sh_start, r.sh_start + r.sh_len - 1)
            hit, ast = alu_hit(alu, chrom, *iv, pad=ALU_NEAR) if iv else (False, "")
            rec["shield_alu"], rec["shield_alu_antisense"] = hit, hit and ast != r.strand
            dn, dlen = downstream_seq(seq, r, DON2[1] + REACH)
            rec["don2_intron_len"] = dlen
            rec["don2_T"] = m42.max_runs(dn, *DON2)["T"] if dlen >= DON2_MIN_INTRON else np.nan
            key = int(r.es) if r.strand == "+" else int(r.ee)
            rec["acc_window"] = acceptor_window(seq, r.strand, key)
            rec["exon_len"] = int(r.ee - r.es)
            rec["exon_gc"] = exon_gc(seq, int(r.es), int(r.ee))
            rec["acc_annotated"] = (chrom, r.strand, key) in acc_keys
            rows.append(rec)
        for g in genes_by_chrom.get(chrom, []):                 # host gene's coding mononucleotide runs
            cds_run[g] = cds_max_run(seq, models[g])
        log(f"{chrom}: {int(f.chrom.eq(chrom).sum()):,} exons, {len(ic):,} introns")
    a = f.merge(pd.DataFrame(rows), on="feature_id", how="left")
    a["struct"] = (a.acc_intron_len >= m42.MIN_INTRON["acc"]) & a.acc_ag.eq("AG")
    at = a.acc_T.to_numpy(float)
    a["d_bp"] = np.where(at >= 10, a.d_bp_T, np.where(at <= 5, a.d_bp_11, np.nan))   # step 51's d_bp
    # geometry of the acceptor T-run (runs >= 10): L, o (position of its 3' end), g (run to
    # branch point), d = o + L - 1 + g; no branch point within 100 nt is its own group
    a["geo_L"], a["geo_o"] = a.run_len, a.run_near
    a["geo_g"] = a.d_bp_T - a.run_far
    a["geo_d"] = a.d_bp_T
    a["bp_group"] = np.where(at < 10, "", np.where(a.d_bp_T.notna(), "branch point", "none within 100 nt"))
    # PWM: check on the even autosomes with a PWM fitted without them, then fit on everything
    pc = [c for c, _ in pos_w]
    nc = [c for c, _ in neg_w]
    W_odd = pwm_fit([w for c, w in pos_w if c not in EVEN], [w for c, w in neg_w if c not in EVEN])
    sp = pwm_score([w for c, w in pos_w if c in EVEN], W_odd)
    sn = pwm_score([w for c, w in neg_w if c in EVEN], W_odd)
    sp, sn = sp[np.isfinite(sp)], sn[np.isfinite(sn)]
    auc_even = float(m51.auc_lower(np.concatenate([sn, sp]), np.r_[np.ones(len(sn), bool), np.zeros(len(sp), bool)]))
    # the core score: the same fit restricted to -3..+3 (PWM columns are fitted independently, so
    # this is the refit on those positions); it cannot see the polypyrimidine tract upstream of -3
    cp = pwm_score([w[CORE] for c, w in pos_w if c in EVEN], W_odd[CORE])
    cn = pwm_score([w[CORE] for c, w in neg_w if c in EVEN], W_odd[CORE])
    cp, cn = cp[np.isfinite(cp)], cn[np.isfinite(cn)]
    auc_core = float(m51.auc_lower(np.concatenate([cn, cp]), np.r_[np.ones(len(cn), bool), np.zeros(len(cp), bool)]))
    W = pwm_fit([w for _, w in pos_w], [w for _, w in neg_w])
    a["pwm3"] = pwm_score(a.acc_window.tolist(), W)
    a["pwm_core"] = pwm_score([w[CORE] if isinstance(w, str) and len(w) == W.shape[0] else "" for w in a.acc_window], W[CORE])
    a["cds_max_run"] = a.feature_id.map(lambda x: cds_run.get(m61.parse_feature(x)["gene"], np.nan))
    # step 61's consequence for every exon with a valid acceptor
    a = a.join(consequence_all(a, models), on="feature_id")
    a = classes(a)
    PROOF.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    a.drop(columns=["acc_window"]).to_parquet(PROOF / "annotation.parquet", index=False)
    pd.DataFrame(W, columns=list(BASES), index=[f"{p}" for p in list(range(-PWM_UP, 0)) + [1, 2, 3]]).rename_axis(
        "position").to_csv(OUT / "pwm3_log2odds.csv", float_format="%.4f")
    qc = {"introns": len(introns), "acceptors": len(acc_keys), "pwm_positive_windows": len(pos_w),
          "pwm_decoy_windows": len(neg_w), "pwm_auroc_even_autosomes_fitted_without_them": auc_even,
          "pwm_n_even": [int(len(sp)), int(len(sn))], "pwm_core_auroc_even_autosomes": auc_core,
          "chromosomes_positive": len(set(pc)),
          "chromosomes_decoy": len(set(nc))}
    qc.update(reproduce(a))
    (OUT / "annotation_qc.json").write_text(json.dumps(qc, indent=1, default=float))
    log(f"annotate: {len(a):,} exons; PWM AUROC on even autosomes {auc_even:.3f} (core {auc_core:.3f}); "
        f"checks {json.dumps(qc, default=float)}")
    sizes_annotation(a)


def consequence_all(a: pd.DataFrame, models: dict) -> pd.DataFrame:
    """Step 61's annotation (place, then translate inserted exons) for every exon with a valid
    acceptor, indexed by feature_id."""
    rows = []
    sub = a[a.struct]
    for r in sub.itertuples(index=False):
        fe = m61.parse_feature(r.feature_id)
        where, i = m61.place(fe, models.get(fe["gene"]))
        rows.append({"feature_id": r.feature_id, "gene_id": fe["gene"], "where": where, "i": i, **fe})
    c = pd.DataFrame(rows)
    c["consequence"] = np.where(c["where"].eq("canonical"), "canonical", "other")
    c["dist"], c["peptide"] = np.nan, np.nan
    ins = c[c["where"].eq("inserted")]
    for chrom, cs in m42.read_genome(set(ins.chrom)):
        for k, r in ins[ins.chrom.eq(chrom)].iterrows():
            m = models[r.gene_id]
            ex, cds, i, st = m["exons"], m["cds"], int(r.i), r.strand
            coding = sum(c_[1] - c_[0] + 1 for c_ in cds[:i + 1] if c_ is not None) - m["phase"]
            p = coding % 3
            up = m61._seq(cs, *cds[i], st)
            tail = up[len(up) - p:] if p else ""
            res = m61.consequence(tail, m61._seq(cs, r.es + 1, r.ee, st), [m61._seq(cs, *e, st) for e in ex[i + 1:]])
            c.loc[k, ["consequence", "dist", "peptide"]] = [res["consequence"], res["dist"], res["peptide"]]
    return c.set_index("feature_id")[["gene_id", "where", "consequence", "dist", "peptide"]]


def reproduce(a: pd.DataFrame) -> dict:
    """Every recomputed column must equal the stored tables of steps 42, 51, 52 and 61."""
    out = {}
    old = pd.read_parquet(LOCAL / "se_repeats_pan.parquet").set_index("feature_id").loc[a.feature_id]
    for col in ["acc_A", "acc_C", "acc_G", "acc_T", "deep_T", "don_T", "acc_intron_len", "acc_ag", "acc_pyr_frac"]:
        x, y = a[col].to_numpy(), old[col].to_numpy()
        same = (x == y) | (pd.isna(x) & pd.isna(y)) if col == "acc_ag" else np.isclose(
            x.astype(float), y.astype(float), equal_nan=True)
        out[f"same_{col}"] = int(same.sum())
    geo = pd.read_parquet(LOCAL / "code_geometry.parquet").drop_duplicates("feature_id").set_index("feature_id")
    j = a.set_index("feature_id").join(geo[["d_bp"]].rename(columns={"d_bp": "old"}), how="inner")
    out["d_bp_rows_compared"] = len(j)
    out["same_d_bp"] = int(np.isclose(j.d_bp.astype(float), j.old.astype(float), equal_nan=True).sum())
    pc = m52.pan_classes().set_index("feature_id")
    j = a.set_index("feature_id").join(pc[["cls"]].rename(columns={"cls": "old"}))
    for c_ in ["R_up", "R_down", "B_far", "B_near"]:
        out[f"n_{c_}_step52"] = int(j.old.eq(c_).sum())
        out[f"same_{c_}_membership"] = int((j.old.eq(c_) == j.cls.eq(c_)).sum())
    ca = pd.read_parquet(LOCAL / "consequence_annotation.parquet")
    feats = pd.read_parquet(PAN / "se" / "features.parquet")["feature_id"]
    ca["feature_id"] = feats.to_numpy()[ca.feature_index.to_numpy()]
    j = ca.set_index("feature_id")[["consequence", "dist"]].join(a.set_index("feature_id")[["consequence", "dist"]],
                                                                  rsuffix="_new", how="inner")
    out["consequence_rows_compared"] = len(j)
    out["same_consequence"] = int((j.consequence == j.consequence_new).sum())
    out["same_dist"] = int(np.isclose(j.dist.astype(float), j.dist_new.astype(float), equal_nan=True).sum())
    out["n_exons"] = len(a)
    return out


CLASS_ROWS = [("R_up", lambda a: a.cls.eq("R_up")), ("R_down", lambda a: a.cls.eq("R_down")),
              ("long_near", lambda a: a.cls.eq("long_near")), ("B_far", lambda a: a.cls.eq("B_far")),
              ("B_near", lambda a: a.cls.eq("B_near")), ("mid (run 6-9)", lambda a: a.cls.eq("mid")),
              ("ctl_polyA", lambda a: a.ctl_polyA), ("ctl_polyA15", lambda a: a.ctl_polyA15),
              ("ctl_ppt_far", lambda a: a.ctl_ppt_far),
              ("ctl_deep", lambda a: a.ctl_deep), ("ctl_deep_far", lambda a: a.ctl_deep_far),
              ("ctl_don", lambda a: a.ctl_don), ("ctl_don_far", lambda a: a.ctl_don_far),
              ("shield", lambda a: a.shield),
              ("ctl_alu (acc_T <= 7, acceptor 1..50 in an Alu)", lambda a: a.ctl_alu),
              ("ctl_alu, the Alu antisense to the gene", lambda a: a.ctl_alu_antisense),
              ("R_up, acceptor 1..50 in an antisense Alu", lambda a: a.cls.eq("R_up") & a.acc_alu_antisense.eq(True)),
              ("R_up, run 20-24 nt", lambda a: a.cls.eq("R_up") & (a.acc_T < 25)),
              ("runs >= 20, branch point >= 50", lambda a: (a.acc_T >= LONG) & (a.d_bp >= 50)),
              ("runs >= 20, branch point > 60", lambda a: (a.acc_T >= LONG) & (a.d_bp > 60)),
              ("runs >= 20, no branch point within 100", lambda a: (a.acc_T >= LONG) & a.d_bp.isna()),
              ("runs >= 20, run in an Alu", lambda a: (a.acc_T >= LONG) & a.run_alu),
              ("runs >= 20, run in an antisense Alu", lambda a: (a.acc_T >= LONG) & a.run_alu_antisense),
              ("R_up, run in an antisense Alu", lambda a: a.cls.eq("R_up") & a.run_alu_antisense),
              ("R_down, run in an antisense Alu", lambda a: a.cls.eq("R_down") & a.run_alu_antisense),
              ("shield, run in an antisense Alu", lambda a: a.shield & a.shield_alu_antisense),
              ("B, exon overlaps an Alu", lambda a: a.cls.isin(["B_far", "B_near"]) & a.exon_alu)]


def sizes_annotation(a: pd.DataFrame) -> pd.DataFrame:
    """Label-blind class sizes over exons with a valid acceptor, with alternative control cuts."""
    s = a[a.struct]
    at, dp = s.acc_T.to_numpy(float), np.nan_to_num(s.d_bp_ppt.to_numpy(float), nan=-1)
    rows = [{"class": k, "rule": "pre-registered draft", "n": int(fn(s).sum())} for k, fn in CLASS_ROWS]
    alt = []
    for t in (7, 9):
        nt = at <= t
        for amin in (15, 20):
            alt.append(("ctl_polyA", f"acc_A >= {amin}, acc_T <= {t}", nt & (s.acc_A >= amin)))
        for pmin in (15, 20, 25):
            for homo in (7, 9):
                k = (nt & (s.ppt_len >= pmin) & (s.ppt_max_T <= homo) & (s.ppt_max_C <= homo) & (s.acc_C <= homo)
                     & (s.ppt_alt < PPT_MAX_ALT))
                alt.append(("ctl_ppt", f"PPT >= {pmin}, homopolymers <= {homo}, acc_T <= {t}, any bp", k))
                alt.append(("ctl_ppt_far", f"PPT >= {pmin}, homopolymers <= {homo}, acc_T <= {t}, bp >= 36", k & (dp >= FAR)))
        for c_, col in [("ctl_deep", "deep_T"), ("ctl_don", "don2_T")]:
            for rmin in (15, 20):
                alt.append((c_, f"{col} >= {rmin}, acc_T <= {t}", nt & (s[col] >= rmin)))
                alt.append((c_ + "_far", f"{col} >= {rmin}, acc_T <= {t}, bp >= 36", nt & (s[col] >= rmin) & (dp >= FAR)))
    rows += [{"class": c_, "rule": rule, "n": int(np.asarray(k).sum())} for c_, rule, k in alt]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "class_sizes_annotation.csv", index=False)
    log("class sizes (valid acceptor):\n" + t.to_string(index=False))
    return t


# ================================================================ experimental branch points (label-blind)

def read_chain(path: Path = CHAIN, pin=CHAIN_PIN) -> dict:
    """UCSC chain file -> {source chrom: DataFrame of ungapped blocks (t, size, q, q_chrom, q_size,
    q_minus, score)} sorted by t. Source (t) coordinates are 0-based on the + strand; a '-' query
    strand gives query coordinates on the reverse strand (forward = q_size - 1 - q)."""
    if pin is not None and (path.stat().st_size != pin[0] or md5_file(path) != pin[1]):
        raise SystemExit(f"{path}: size or MD5 differs from the pin {pin}")
    rows = []
    t = q = 0
    head = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            if f[0] == "chain":
                head = (f[2], f[7], int(f[8]), f[9] == "-", float(f[1]))
                if f[4] != "+":
                    raise ValueError("chain with a '-' reference strand")
                t, q = int(f[5]), int(f[10])
                continue
            size = int(f[0])
            rows.append((head[0], t, size, q, head[1], head[2], head[3], head[4]))
            if len(f) == 3:
                t += size + int(f[1])
                q += size + int(f[2])
    d = pd.DataFrame(rows, columns=["t_chrom", "t", "size", "q", "q_chrom", "q_size", "q_minus", "score"])
    return {c: g.sort_values(["t", "score"], ascending=[True, False]).reset_index(drop=True) for c, g in d.groupby("t_chrom")}


def lift_points(chain: dict, chrom: str, pos: np.ndarray, strand: np.ndarray) -> pd.DataFrame:
    """Lift 0-based source positions (with their strands) through the chain blocks. A position
    covered by several blocks takes the highest-scoring chain; unmapped positions get q_chrom ''."""
    pos = np.asarray(pos, np.int64)
    out = pd.DataFrame({"q_chrom": "", "q_pos": -1, "q_strand": ""}, index=range(len(pos)))
    if chrom not in chain or not len(pos):
        return out
    b = chain[chrom]
    t, size = b.t.to_numpy(), b["size"].to_numpy()
    end_max = np.maximum.accumulate(t + size)               # a block can start before and still cover
    best = np.full(len(pos), -1)
    best_score = np.full(len(pos), -np.inf)
    j = np.searchsorted(t, pos, side="right") - 1
    k = j.copy()
    active = k >= 0
    while active.any():                                     # walk back over blocks that may cover pos
        i = np.flatnonzero(active)
        kk = k[i]
        cover = (t[kk] <= pos[i]) & (pos[i] < t[kk] + size[kk])
        sc = b.score.to_numpy()[kk]
        better = cover & (sc > best_score[i])
        best[i[better]], best_score[i[better]] = kk[better], sc[better]
        k[i] -= 1
        active[i] = (k[i] >= 0) & (end_max[np.maximum(k[i], 0)] > pos[i])
    ok = best >= 0
    bb = b.iloc[best[ok]]
    q = bb.q.to_numpy() + (pos[ok] - bb.t.to_numpy())
    minus = bb.q_minus.to_numpy()
    out.loc[ok, "q_chrom"] = bb.q_chrom.to_numpy()
    out.loc[ok, "q_pos"] = np.where(minus, bb.q_size.to_numpy() - 1 - q, q)
    flip = {"+": "-", "-": "+"}
    st = np.asarray(strand)[ok]
    out.loc[ok, "q_strand"] = np.where(minus, [flip[s] for s in st], st)
    return out


def sense_base(chrom_seq: str, pos: np.ndarray, strand: np.ndarray) -> np.ndarray:
    b = np.array([chrom_seq[p] if 0 <= p < len(chrom_seq) else "N" for p in pos], dtype=object)
    return np.where(np.asarray(strand) == "-", [x.translate(COMP) for x in b], b)


def exon_bps(a: pd.DataFrame, bp: pd.DataFrame) -> pd.DataFrame:
    """Per exon, the experimental branch points on its strand at positions BP_RANGE (1 = the last
    intron base before the exon), within the intron: their number, the primary one (most eBP
    sources; ties: nearest the AG) and the nearest; whether any lies between the T-run's 3' end and
    the AG (p < run_near); g_exp = primary position - run_far."""
    out = []
    for (chrom, strand), sub in a.groupby(["chrom", "strand"]):
        q = bp[(bp.chrom == chrom) & (bp.strand == strand)].sort_values("pos")
        qp, qn = q.pos.to_numpy(), q.n_ebp.to_numpy()
        for r in sub.itertuples(index=False):
            top = min(BP_RANGE[1], int(r.acc_intron_len)) if np.isfinite(r.acc_intron_len) else BP_RANGE[1]
            if strand == "+":                                  # position p <-> genomic es - p
                lo, hi = int(r.es) - top, int(r.es) - BP_RANGE[0]
            else:                                              # position p <-> genomic ee + p - 1
                lo, hi = int(r.ee) + BP_RANGE[0] - 1, int(r.ee) + top - 1
            i0, i1 = np.searchsorted(qp, lo), np.searchsorted(qp, hi, side="right")
            if i1 <= i0:
                out.append({"feature_id": r.feature_id, "n_exp_bp": 0})
                continue
            p = (int(r.es) - qp[i0:i1]) if strand == "+" else (qp[i0:i1] - int(r.ee) + 1)
            n = qn[i0:i1]
            prim = np.lexsort((p, -n))[0]
            rec = {"feature_id": r.feature_id, "n_exp_bp": int(len(p)), "d_exp": int(p[prim]),
                   "d_exp_sources": int(n[prim]), "d_exp_near": int(p.min())}
            if np.isfinite(r.run_near):
                rec["exp_bp_run_to_ag"] = bool((p < r.run_near).any())
                rec["exp_primary_in_run"] = bool(r.run_near <= p[prim] <= r.run_far)
                rec["g_exp"] = int(p[prim]) - int(r.run_far)
            out.append(rec)
    return pd.DataFrame(out)


def branchpoints() -> None:
    """Lift BPHunter's experimental branch points to GRCh38, validate the lift on the genome, place
    them on the SE exons and report coverage and agreement with the sequence heuristic."""
    if not BPHUNTER.exists():
        raise SystemExit(f"{BPHUNTER} is not present; step 3's experimental branch-point test is not testable")
    if BPHUNTER.stat().st_size != BPHUNTER_PIN[0] or md5_file(BPHUNTER) != BPHUNTER_PIN[1]:
        raise SystemExit(f"{BPHUNTER}: size or MD5 differs from the pin {BPHUNTER_PIN}")
    raw = pd.read_csv(BPHUNTER, sep="\t", dtype={"CHROM": str, "STRAND": str, "BP_NAME": str, "SOURCE": str})
    raw["n_ebp"] = raw.SOURCE.str.split(",").apply(lambda s: sum(x.startswith("eBP_") for x in s))
    raw["named"] = raw.BP_NAME.str[-1]
    chain = read_chain()
    lifted, lifted0 = [], []
    for chrom, g in raw.groupby("CHROM"):
        lf = lift_points(chain, chrom, g.BP_POS.to_numpy() - 1, g.STRAND.to_numpy())   # BP_POS read as 1-based
        lf.index = g.index
        lifted.append(lf)
        l0 = lift_points(chain, chrom, g.BP_POS.to_numpy(), g.STRAND.to_numpy())       # check: read as 0-based
        l0.index = g.index
        lifted0.append(l0.add_suffix("_0"))
    raw = raw.join(pd.concat(lifted)).join(pd.concat(lifted0))
    raw["same_chrom"] = raw.q_chrom.eq(raw.CHROM)
    raw["grch38_base"] = ""
    raw["grch38_base_0"] = ""
    for chrom, seq in m42.read_genome(set(raw.q_chrom[raw.same_chrom]) | set(raw.q_chrom_0[raw.q_chrom_0.eq(raw.CHROM)])):
        for sfx in ["", "_0"]:
            k = raw.index[raw["q_chrom" + sfx].eq(chrom) & raw["q_chrom" + sfx].eq(raw.CHROM)]
            raw.loc[k, "grch38_base" + sfx] = sense_base(seq, raw.loc[k, "q_pos" + sfx].to_numpy(),
                                                         raw.loc[k, "q_strand" + sfx].to_numpy())
    exp = raw.n_ebp >= 1
    lifted_ok = raw.same_chrom
    named_a = exp & lifted_ok & raw.named.eq("A")
    qc = {"rows": len(raw), "experimental_rows": int(exp.sum()),
          "experimental_unliftable": int((exp & raw.q_chrom.eq("")).sum()),
          "experimental_to_another_chrom": int((exp & ~raw.q_chrom.eq("") & ~raw.same_chrom).sum()),
          "experimental_named_A_read_A_grch38": float((raw.grch38_base[named_a] == "A").mean()),
          "experimental_named_A_read_A_if_0based": float((raw.grch38_base_0[named_a] == "A").mean()),
          "experimental_named_base_equals_grch38": float((raw.grch38_base[exp & lifted_ok] == raw.named[exp & lifted_ok]).mean()),
          "experimental_read_A_grch38_any_name": float((raw.grch38_base[exp & lifted_ok] == "A").mean()),
          "all_named_base_equals_grch38": float((raw.grch38_base[lifted_ok] == raw.named[lifted_ok]).mean())}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "branchpoints_qc.json").write_text(json.dumps(qc, indent=1, default=float))
    log(f"lift: {json.dumps(qc, default=float)}")
    if qc["experimental_named_A_read_A_grch38"] < LIFT_MIN_A:
        raise SystemExit(f"liftover validation failed: {qc['experimental_named_A_read_A_grch38']:.3f} < {LIFT_MIN_A}")
    bp = raw[exp & lifted_ok].rename(columns={"q_chrom": "chrom", "q_pos": "pos", "q_strand": "strand"})[
        ["chrom", "pos", "strand", "n_ebp", "named", "grch38_base", "SOURCE"]]
    bp.to_parquet(PROOF / "bphunter_experimental_grch38.parquet", index=False)
    a = pd.read_parquet(PROOF / "annotation.parquet")
    e = exon_bps(a, bp)
    e.to_parquet(PROOF / "bp_experimental.parquet", index=False)
    coverage(a.merge(e, on="feature_id", how="left"))


def coverage(a: pd.DataFrame) -> pd.DataFrame:
    """Label-blind: share of each class with an experimental branch point, and agreement of the
    heuristic d_bp with the experimental one."""
    from scipy.stats import spearmanr
    a = a[a.struct]
    rows = []
    groups = [("R_up", a.cls.eq("R_up")), ("long_near", a.cls.eq("long_near")), ("runs >= 20", a.acc_T >= LONG),
              ("R_down", a.cls.eq("R_down")), ("B_far", a.cls.eq("B_far")), ("B_near", a.cls.eq("B_near")),
              ("B", a.cls.isin(["B_far", "B_near"])), ("shield", a.shield)]
    for name, m in groups:
        s = a[m]
        has = s.n_exp_bp.fillna(0) > 0
        both = s[has & s.d_bp.notna()]
        rho = spearmanr(both.d_bp, both.d_exp).statistic if len(both) > 2 else np.nan
        row = {"class": name, "n": len(s), "with_exp_bp": int(has.sum()), "share_with_exp_bp": float(has.mean()),
               "median_exp_bps": float(s.n_exp_bp[has].median()) if has.any() else np.nan,
               "heuristic_and_exp_bp": len(both), "spearman_d_bp_vs_d_exp": rho,
               "share_within_3nt": float(((both.d_bp - both.d_exp).abs() <= 3).mean()) if len(both) else np.nan,
               "share_within_3nt_nearest": float(((both.d_bp - both.d_exp_near).abs() <= 3).mean()) if len(both) else np.nan}
        if name in ("R_up", "long_near", "runs >= 20", "R_down", "shield"):
            row["exp_bp_between_run_and_ag"] = int(s.exp_bp_run_to_ag.eq(True).sum())
            up = s[has & (s.g_exp > 0)]
            row["primary_exp_bp_upstream_of_run"] = len(up)
            row["primary_exp_bp_in_run"] = int(s.exp_primary_in_run.eq(True).sum())
            row["upstream_d_exp_ge_36"] = int((up.d_exp >= FAR).sum())
            row["upstream_d_exp_lt_36"] = int((up.d_exp < FAR).sum())
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "branchpoints_coverage.csv", index=False, float_format="%.4g")
    log("experimental branch points by class:\n" + t.to_string(index=False))
    return t


# ================================================================ counts (label-blind)

def csr_take(indptr, indices, chans: dict, rows, cols) -> dict:
    """Dense [len(rows) x len(cols)] arrays of each channel for the given CSR rows (samples) and
    columns (event indices, any order; indices must ascend within a row). Absent entries are 0 in
    integer channels and NaN in float ones."""
    cols = np.asarray(cols, np.int64)
    order = np.argsort(cols, kind="stable")
    cs = cols[order]
    out = {k: np.full((len(rows), len(cols)), np.nan if v.dtype.kind == "f" else 0,
                      dtype=np.float32 if v.dtype.kind == "f" else np.int32) for k, v in chans.items()}
    for i, r in enumerate(rows):
        lo, hi = int(indptr[r]), int(indptr[r + 1])
        idx = np.asarray(indices[lo:hi])
        pos = np.searchsorted(idx, cs)
        hit = pos < len(idx)
        hit[hit] = idx[pos[hit]] == cs[hit]
        src = lo + pos[hit]
        dst = order[hit]
        for k, v in chans.items():
            out[k][i, dst] = np.asarray(v[src])
    return out


def psi_from_counts(ijc, sjc, ratio):
    """rMATS PSI = (IJC/inc_len) / (IJC/inc_len + SJC/skip_len) = IJC / (IJC + ratio * SJC), with
    ratio = inc_len / skip_len; NaN where there are no reads."""
    ijc, sjc = np.asarray(ijc, float), np.asarray(sjc, float)
    den = ijc + np.asarray(ratio, float) * sjc
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((ijc + sjc > 0) & (den > 0), ijc / den, np.nan)


def store_samples(store: Path) -> pd.DataFrame:
    s = pd.read_csv(store / "samples.tsv", sep="\t", dtype=str, keep_default_na=False)
    s["sample_idx"] = s.sample_idx.astype(int)
    return s


def count_exons(a: pd.DataFrame) -> pd.DataFrame:
    """The exons stored: every member of a class or control with a valid acceptor (B itself
    excepted), and the B sample - a fixed random sample of 12,000 of all B exons (acc_T <= 5, valid
    acceptor; seed 0), drawn whatever their other memberships, so that B-sample members of any B
    subset (B_far, an Alu acceptor) are a random sample of it. b_sample flags the B sample."""
    s = a[a.struct].copy()
    ctl = s[[c for c in s.columns if c.startswith("ctl_")] + ["shield"]].any(axis=1)
    keep = s.cls.isin(["R_up", "R_down", "long_near"]) | ctl
    b = s[s.cls.isin(["B_far", "B_near"])]
    rng = np.random.default_rng(SEED)
    pick = b.index[np.sort(rng.choice(len(b), min(B_SAMPLE, len(b)), replace=False))]
    s["b_sample"] = s.index.isin(pick)
    s["stored_as"] = np.where(keep, "class member", np.where(s.b_sample, "background sample", ""))
    return s[s.stored_as.ne("")].reset_index(drop=True)


def _count_cancer(proj: str, samp: pd.DataFrame, cols: np.ndarray, store: str, dest: str | None = None) -> dict:
    store, dest = Path(store), Path(dest) if dest else PROOF
    qc = {"project": proj, "n_samples": len(samp)}
    for v in VIEWS:
        d = store / "sparse" / f"se_{v}"
        mm = {k: np.load(d / f"{k}.npy", mmap_mode="r") for k in ["indptr", "indices", "ijc", "sjc", "inc_len", "skip_len", "psi"]}
        chans = {k: mm[k] for k in ["ijc", "sjc", "inc_len", "skip_len", "psi"]}
        got = csr_take(mm["indptr"], mm["indices"], chans, samp.sample_idx.to_numpy(), cols)
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(got["skip_len"] > 0, got["inc_len"] / np.maximum(got["skip_len"], 1), np.nan).astype(np.float32)
        psi = psi_from_counts(got["ijc"], got["sjc"], ratio)
        diff = np.abs(psi - got["psi"])
        qc[f"{v}_entries"] = int(np.isfinite(got["psi"]).sum())
        qc[f"{v}_max_abs_psi_diff"] = float(np.nanmax(diff)) if np.isfinite(diff).any() else 0.0
        qc[f"{v}_psi_nan_mismatch"] = int((np.isfinite(psi) != np.isfinite(got["psi"])).sum())
        for k, x in [("ijc", got["ijc"]), ("sjc", got["sjc"]), ("ratio", ratio)]:
            np.save(dest / f"{proj}__{v}__{k}.npy", x)
    samp.to_parquet(dest / f"{proj}__samples.parquet", index=False)
    return qc


def counts() -> None:
    from joblib import Parallel, delayed
    store = Path(expand("${CTSM_SPLICE_COUNTS}"))
    if not (store / "samples.tsv").exists():
        raise SystemExit("CTSM_SPLICE_COUNTS is not set or the store is missing (config/local.env)")
    a = load_annotation()
    ex = count_exons(a)
    ev = pd.read_csv(store / "events" / "se.events.tsv.gz", sep="\t", usecols=["event_idx", "id"], dtype={"id": str})
    ev = ev.set_index("id").event_idx
    miss = ~ex.feature_id.isin(ev.index)
    if miss.any():
        raise SystemExit(f"{int(miss.sum())} exons are not in the count store")
    ex["event_idx"] = ev.loc[ex.feature_id].to_numpy()
    ex.to_parquet(PROOF / "exons.parquet", index=False)
    cols = ex.event_idx.to_numpy(np.int64)
    s = store_samples(store)
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")[["file_id", "case_id", "project"]]
    pos = s.set_index("File.ID")
    if not t.file_id.isin(pos.index).all():
        raise SystemExit("some D-71 tumours are not in the count store")
    t = t.assign(role="tumour", sample_idx=pos.loc[t.file_id, "sample_idx"].to_numpy(),
                 sample_type=pos.loc[t.file_id, "Sample.Type"].to_numpy())
    projects = set(t.project)
    nm = s[s["Sample.Type"].eq("Solid Tissue Normal") & s["Project.ID"].str.replace("TCGA-", "", regex=False).isin(projects)]
    nm = pd.DataFrame({"file_id": nm["File.ID"], "case_id": nm["Case.ID"], "project": nm["Project.ID"].str[5:],
                       "role": "normal", "sample_idx": nm.sample_idx, "sample_type": nm["Sample.Type"]})
    samples = pd.concat([t, nm], ignore_index=True)
    pan = set(pd.read_parquet(PAN / "se" / "samples.parquet").file_id)
    samples["in_pan_matrix"] = samples.file_id.isin(pan)          # pair rule: prefer the matrix's normal
    log(f"counts: {len(ex):,} exons ({ex.stored_as.value_counts().to_dict()}), {len(t):,} tumours, {len(nm):,} normals")
    PROOF.mkdir(parents=True, exist_ok=True)
    jobs = [(p, g.sort_values(["role", "sample_idx"], ascending=[False, True]).reset_index(drop=True))
            for p, g in samples.groupby("project")]
    jobs.sort(key=lambda x: -len(x[1]))
    qc = Parallel(n_jobs=N_JOBS)(delayed(_count_cancer)(p, g, cols, str(store)) for p, g in jobs)
    meta = json.loads((store / "config.json").read_text())
    sm = json.loads((store / "sparse" / "se_jc" / "meta.json").read_text())
    info = {"store": portable(store), "layout_sha256": meta.get("layout_sha256"),
            "manifest_sha256": meta.get("manifest_sha256"), "se_jc_nnz": sm.get("nnz"),
            "n_exons": len(ex), "n_tumours": len(t), "n_normals": len(nm),
            "max_abs_psi_diff": {v: max(q[f"{v}_max_abs_psi_diff"] for q in qc) for v in VIEWS},
            "psi_nan_mismatch": {v: sum(q[f"{v}_psi_nan_mismatch"] for q in qc) for v in VIEWS}}
    pd.DataFrame(qc).to_csv(OUT / "counts_qc_by_cancer.csv", index=False)
    (OUT / "counts_qc.json").write_text(json.dumps(info, indent=1, default=float))
    log(f"counts written: {json.dumps(info, default=float)}")
    check_jcec_matches_matrix(ex, samples)
    sizes_counts(ex, samples)


def load_counts(proj: str, view: str = "jc", mmap: bool = True):
    """(samples, IJC, SJC, ratio) of one cancer; columns follow proof/exons.parquet."""
    mode = "r" if mmap else None
    samp = pd.read_parquet(PROOF / f"{proj}__samples.parquet")
    return (samp, *(np.load(PROOF / f"{proj}__{view}__{k}.npy", mmap_mode=mode) for k in ["ijc", "sjc", "ratio"]))


def check_jcec_matches_matrix(ex: pd.DataFrame, samples: pd.DataFrame, n: int = 3) -> None:
    """The stored JCEC PSI of a few tumours equals the repository's SE matrix (step 02)."""
    from ctsm.matrix_io import load_view
    feats = pd.read_parquet(PAN / "se" / "features.parquet")["feature_id"]
    col = pd.Series(np.arange(len(feats)), index=feats).loc[ex.feature_id].to_numpy()
    out = []
    for proj in sorted(samples.project.unique())[:n]:
        samp, ijc, sjc, ratio = load_counts(proj, "jcec")
        k = int(np.flatnonzero(samp.role.eq("tumour"))[0])
        mine = psi_from_counts(ijc[k], sjc[k], ratio[k])
        repo = load_view(PAN, "se", [samp.file_id[k]])[0].get()[0, col]
        both = np.isfinite(mine) & np.isfinite(repo)
        out.append({"project": proj, "compared": int(both.sum()),
                    "max_abs_diff": float(np.max(np.abs(mine[both] - repo[both]))) if both.any() else np.nan,
                    "finite_mismatch": int((np.isfinite(mine) != np.isfinite(repo)).sum())})
    pd.DataFrame(out).to_csv(OUT / "counts_jcec_vs_matrix.csv", index=False)
    log(f"JCEC vs SE matrix:\n{pd.DataFrame(out).to_string(index=False)}")


BP_ROWS = [("runs >= 20, experimental bp upstream of the run", lambda a: (a.acc_T >= LONG) & (a.g_exp > 0)),
           ("  of which d_exp >= 36", lambda a: (a.acc_T >= LONG) & (a.g_exp > 0) & (a.d_exp >= FAR)),
           ("  of which d_exp < 36", lambda a: (a.acc_T >= LONG) & (a.g_exp > 0) & (a.d_exp < FAR)),
           ("runs >= 20, experimental bp between the run and the AG", lambda a: (a.acc_T >= LONG) & a.exp_bp_run_to_ag.eq(True)),
           ("R_up with an experimental bp", lambda a: a.cls.eq("R_up") & (a.n_exp_bp > 0)),
           ("B, experimental d_exp >= 36", lambda a: a.cls.isin(["B_far", "B_near"]) & (a.d_exp >= FAR)),
           ("B, experimental d_exp < 36", lambda a: a.cls.isin(["B_far", "B_near"]) & (a.d_exp < FAR))]
MIN_MMR_ELIGIBLE = 20          # an exon enters exon-level tests if the cancers whose universe holds it
                               # carry >= 20 MMRd tumours (a per-cancer label count, fixed under
                               # within-cancer permutation)


def load_annotation() -> pd.DataFrame:
    """The step's annotation with the classes (re)applied and the experimental branch points joined."""
    a = classes(pd.read_parquet(PROOF / "annotation.parquet"))
    f = PROOF / "bp_experimental.parquet"
    if f.exists():
        a = a.merge(pd.read_parquet(f), on="feature_id", how="left")
        a["n_exp_bp"] = a.n_exp_bp.fillna(0).astype(int)
    return a


def universe_shares(ex: pd.DataFrame, projects) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Label-blind, per exon x cancer: the share of the cancer's tumours measured in JC, and the
    median JC reads (IJC + SJC) over its tumours."""
    share, cov = {}, {}
    for p in projects:
        samp, ijc, sjc, _ = load_counts(p, "jc")
        tum = samp.role.eq("tumour").to_numpy()
        n = np.asarray(ijc[tum], np.int64) + np.asarray(sjc[tum], np.int64)
        share[p], cov[p] = (n >= MIN_READS).mean(0), np.median(n, axis=0)
    return pd.DataFrame(share, index=ex.feature_id), pd.DataFrame(cov, index=ex.feature_id)


def sizes_counts(ex: pd.DataFrame | None = None, samples: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per class: exons in the JC universe (measured, JC IJC + SJC >= 10, in >= 50% of a cancer's
    tumours; also 80%) of held-out and training cancers, those eligible for exon-level tests (the
    cancers holding them carry >= 20 MMRd tumours: per-cancer label counts only, no PSI compared),
    and their median JC reads where in the universe. Also tumour-normal pairs per cancer."""
    a = load_annotation().set_index("feature_id")
    if ex is None:
        ex = pd.read_parquet(PROOF / "exons.parquet")
    ex = ex[["feature_id", "stored_as", "b_sample"]].join(a.drop(columns=["stored_as", "b_sample"], errors="ignore"),
                                                         on="feature_id")
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")
    n_mmr = t.groupby("project").mmr.sum()
    projects = sorted(t.project.unique())
    held = [p for p in projects if p not in TRAIN]
    S, C = universe_shares(ex, projects)
    rows = []
    member = {k: fn(ex).to_numpy() for k, fn in CLASS_ROWS}
    if "d_exp" in ex:
        member.update({k: np.asarray(fn(ex), bool) for k, fn in BP_ROWS})
    member["B sample (all)"] = ex.b_sample.to_numpy()
    member["B_far in the B sample (RP3a)"] = (ex.b_sample & ex.cls.eq("B_far")).to_numpy()
    member["B in the B sample, acceptor 1..50 in an Alu"] = (ex.b_sample & ex.acc_alu.eq(True)).to_numpy()
    for k, m in member.items():
        for thr in (UNIVERSE_SHARE, 0.8):
            U = S.loc[m] >= thr
            mm_h = (U[held] * n_mmr[held]).sum(1)
            mm_t = (U[TRAIN] * n_mmr[TRAIN]).sum(1)
            cov = C.loc[m].where(U).stack()
            rows.append({"class": k, "universe_share": thr, "n_annotated": int(m.sum()),
                         "in_universe_ge1_held_out": int((U[held].sum(1) >= 1).sum()),
                         "eligible_held_out": int((mm_h >= MIN_MMR_ELIGIBLE).sum()),
                         "in_universe_ge1_training": int((U[TRAIN].sum(1) >= 1).sum()),
                         "eligible_training": int((mm_t >= MIN_MMR_ELIGIBLE).sum()),
                         "median_jc_reads_in_universe": float(cov.median()) if len(cov) else np.nan})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "class_sizes_jc_universe.csv", index=False)
    if samples is None:
        samples = pd.concat([pd.read_parquet(f) for f in sorted(PROOF.glob("*__samples.parquet"))], ignore_index=True)
    tum = samples[samples.role.eq("tumour")].merge(t[["file_id", "mmr"]], on="file_id")
    nm = samples[samples.role.eq("normal")]
    tum = tum[tum.case_id.isin(nm.case_id)]
    pairs = nm.groupby("project").size().rename("normal_files").to_frame().join(
        tum.groupby("project").agg(pairs=("case_id", "size"), pairs_mmr_deficient=("mmr", "sum")), how="outer").fillna(0)
    pairs.astype(int).reset_index().to_csv(OUT / "normals_by_cancer.csv", index=False)
    pd.set_option("display.width", 250)
    log("class sizes in the JC universe:\n" + out.to_string(index=False))
    log("tumour-normal pairs:\n" + pairs.astype(int).to_string())
    return out


# ================================================================ strata (label-blind)

def coarsen(band, terc, first, ctrl, min_ctrl: int = MIN_CTRL, per_first: int = 0, max_span: int = MAX_SPAN):
    """The pre-specified stratum coarsening of a class contrast (Gate 1). Start from band x coverage
    tertile. In a band where a cell holding first-class exons has fewer than the required control
    exons, merge the band's tertiles. Then, lowest band first, merge a band group still short of
    controls with its neighbouring group holding more control exons (ties - the lower), until every
    group holding first-class exons has enough (or a group spans max_span bands). Required =
    max(min_ctrl, per_first x first-class exons in the stratum); per_first > 0 is for drawing
    matched controls without replacement. Returns (stratum label per exon, '' outside the two
    classes; kept = the exon's stratum holds both classes)."""
    band, terc = np.asarray(band, int), np.asarray(terc, int)
    first, ctrl = np.asarray(first, bool), np.asarray(ctrl, bool)

    def need(nf):
        return max(min_ctrl, per_first * nf)

    def cnt(mask, lo, hi, t=None):
        m = mask & (band >= lo) & (band <= hi)
        return int((m & (terc == t)).sum()) if t is not None else int(m.sum())

    groups = [[b, b, True] for b in range(len(BAND_EDGES) - 1)]
    for g in groups:
        if any(cnt(first, g[0], g[1], t) and cnt(ctrl, g[0], g[1], t) < need(cnt(first, g[0], g[1], t))
               for t in range(3)):
            g[2] = False
    stuck = set()
    while len(groups) > 1:
        short = [i for i, g in enumerate(groups) if not g[2] and i not in stuck and cnt(first, g[0], g[1])
                 and cnt(ctrl, g[0], g[1]) < need(cnt(first, g[0], g[1]))]
        if not short:
            break
        i = short[0]
        nb = [j for j in (i - 1, i + 1) if 0 <= j < len(groups)
              and max(groups[i][1], groups[j][1]) - min(groups[i][0], groups[j][0]) + 1 <= max_span]
        if not nb:                                         # no merge within the span limit
            stuck.add(i)
            continue
        j = max(nb, key=lambda k: (cnt(ctrl, groups[k][0], groups[k][1]), -k))
        lo, hi = min(groups[i][0], groups[j][0]), max(groups[i][1], groups[j][1])
        groups[min(i, j)] = [lo, hi, False]
        del groups[max(i, j)]
        stuck = set()
    lab = np.full(len(band), "", dtype=object)
    for lo, hi, split in groups:
        m = (first | ctrl) & (band >= lo) & (band <= hi)
        lab[m] = [f"{lo}-{hi}:{t}" if split else f"{lo}-{hi}" for t in terc[m]]
    kept = np.zeros(len(band), bool)
    for s_ in set(lab[lab != ""]):
        m = lab == s_
        if (m & first).any() and (m & ctrl).any():
            kept |= m
    return lab, kept


def effective_n(lab, first, ctrl, kept) -> dict:
    """First-class exons kept, control exons kept, and the first-class-weighted effective number
    of controls, 1 / sum_s (w_s^2 / c_s) with w_s = the stratum's share of kept first-class exons."""
    f = pd.Series(lab[kept & first]).value_counts()
    c = pd.Series(lab[kept & ctrl]).value_counts().reindex(f.index)
    w = f / f.sum() if len(f) else f
    core = [s_.split("|")[-1].split(":")[0].split("-") for s_ in f.index]     # "[L bin|]lo-hi[:tertile]"
    span = [int(b) - int(a) + 1 for a, b in core]
    return {"first_kept": int(f.sum()), "control_kept": int((kept & ctrl).sum()),
            "control_effective_n": float(1 / (w ** 2 / c).sum()) if len(f) else 0.0,
            "strata_with_first": len(f), "max_bands_per_stratum": max(span, default=0)}


def exon_strata() -> pd.DataFrame:
    """Per stored exon and per side (held out, training), label-blind: eligibility, the baseline
    (the all-tumour median JC PSI of each universe cancer, weighted by its MMRd count), its band,
    and the median JC reads over all tumours of those cancers with its tertile over the side's
    eligible exons. Only per-cancer MMRd counts are used (the eligibility and baseline weights)."""
    a = load_annotation().set_index("feature_id")
    ex = pd.read_parquet(PROOF / "exons.parquet")[["feature_id", "stored_as", "b_sample"]]
    ex = ex.join(a, on="feature_id")
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")
    n_mmr = t.groupby("project").mmr.sum()
    projects = sorted(t.project.unique())
    sides = {"held_out": [p for p in projects if p not in TRAIN], "training": list(TRAIN)}
    E = len(ex)
    U, M = {}, {}
    for p in projects:
        samp, ijc, sjc, ratio = load_counts(p, "jc")
        tum = samp.role.eq("tumour").to_numpy()
        n = np.asarray(ijc[tum], np.int64) + np.asarray(sjc[tum], np.int64)
        meas = n >= MIN_READS
        U[p] = meas.mean(0) >= UNIVERSE_SHARE
        x = np.where(meas, psi_from_counts(ijc[tum], sjc[tum], ratio[tum]), np.nan)
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")                    # all-NaN columns outside the universe
            M[p] = np.where(U[p], np.nanmedian(np.where(U[p], x, np.nan), axis=0), np.nan)
    for side, ps in sides.items():
        w = np.array([[n_mmr[p]] for p in ps], float) * np.vstack([U[p] for p in ps])
        elig = w.sum(0) >= (MIN_MMR_ELIGIBLE if side == "held_out" else 1)
        with np.errstate(all="ignore"):
            base = np.nansum(np.vstack([np.nan_to_num(M[p]) for p in ps]) * w, 0) / w.sum(0)
        cov = np.full(E, np.nan)
        for lo in range(0, E, 2000):                       # pooled median reads, in exon chunks
            cols = np.arange(lo, min(E, lo + 2000))
            pool = [[] for _ in cols]
            for p in ps:
                samp, ijc, sjc, _ = load_counts(p, "jc")
                tum = samp.role.eq("tumour").to_numpy()
                n = (np.asarray(ijc[:, cols[0]:cols[-1] + 1], np.int64)
                     + np.asarray(sjc[:, cols[0]:cols[-1] + 1], np.int64))[tum]
                for k in np.flatnonzero(U[p][cols]):
                    pool[k].append(n[:, k])
            for k, v in enumerate(pool):
                if v and elig[cols[k]]:
                    cov[cols[k]] = float(np.median(np.concatenate(v)))
        q = np.nanquantile(cov[elig], [1 / 3, 2 / 3])
        ex[f"{side}_eligible"] = elig
        ex[f"{side}_baseline"] = np.where(elig, base, np.nan)
        ex[f"{side}_band"] = np.where(elig, np.digitize(np.nan_to_num(base), BAND_EDGES) - 1, -1)
        ex[f"{side}_reads"] = cov
        ex[f"{side}_tertile"] = np.where(elig, np.digitize(np.nan_to_num(cov), q), -1)
        log(f"strata, {side}: {int(elig.sum()):,} eligible stored exons; reads tertile cuts {q.tolist()}")
    ex.to_parquet(PROOF / "exon_strata.parquet", index=False)
    return ex


def contrast_masks(ex: pd.DataFrame) -> list[tuple]:
    """(id, side, first class, control class, per_first) of every pre-registered class contrast."""
    up, bs = ex.cls.eq("R_up"), ex.b_sample.eq(True)
    long_ = ex.acc_T >= LONG
    exp_up = long_ & (ex.g_exp > 0)
    return [("RP3a", "held_out", up, bs & ex.cls.eq("B_far"), 0),
            ("RP3b", "held_out", up, ex.cls.eq("long_near") & ex.d_bp.notna(), 0),
            ("RP3d", "held_out", exp_up & (ex.d_exp >= FAR), exp_up & (ex.d_exp < FAR), 0),
            ("RP4a", "held_out", up, ex.ctl_polyA15, 0),
            ("RP4b", "held_out", up, ex.ctl_ppt_far, 0),
            ("RP4c", "held_out", up, ex.ctl_deep, 0),
            ("RP4d", "held_out", up, ex.ctl_don, 0),
            ("RP4e", "held_out", up, ex.ctl_alu, 0),
            ("RP4e, antisense Alu (sensitivity)", "held_out", up & ex.acc_alu_antisense.eq(True), ex.ctl_alu_antisense, 0),
            ("RP4a, R_up run 20-24 (sensitivity)", "held_out", up & (ex.acc_T < 25), ex.ctl_polyA15, 0),
            ("RPM1", "held_out", up, ex.shield, 0),
            ("RPM3, L 20-24", "held_out", long_ & (ex.acc_T < 25) & ex.run_alu_antisense.eq(True),
             long_ & (ex.acc_T < 25) & ~ex.run_alu_antisense.eq(True), 0),
            ("RPM3, L >= 25", "held_out", long_ & (ex.acc_T >= 25) & ex.run_alu_antisense.eq(True),
             long_ & (ex.acc_T >= 25) & ~ex.run_alu_antisense.eq(True), 0),
            ("RP2c (5 matched B-sample exons per R_up exon)", "held_out", up, bs, 5),
            ("RPD3a", "training", up, bs & ex.cls.eq("B_far"), 0)]


def strata() -> None:
    ex = exon_strata()
    rows = []
    for cid, side, first, ctrl, per in contrast_masks(ex):
        el = ex[f"{side}_eligible"].to_numpy()
        f, c = first.to_numpy(bool) & el, ctrl.to_numpy(bool) & el
        base = {"contrast": cid, "side": side, "first_eligible": int(f.sum()), "control_eligible": int(c.sum())}
        cell = ex[f"{side}_band"].astype(str).to_numpy() + ":" + ex[f"{side}_tertile"].astype(str).to_numpy()
        both = set(cell[f]) & set(cell[c])
        base["first_kept_band_x_tertile_uncoarsened"] = int((f & np.isin(cell, list(both))).sum())
        lab, kept = coarsen(ex[f"{side}_band"], ex[f"{side}_tertile"], f, c, per_first=per)
        rows.append({**base, **effective_n(lab, f, c, kept)})
    # instrument check: R_down exons drawn to R_up's MIRRORED baseline (the band of 1 - baseline,
    # room to move down as R_up has room to move up) x coverage tertile, one R_down per R_up exon
    el = ex.held_out_eligible.to_numpy()
    f, c = ex.cls.eq("R_up").to_numpy() & el, ex.cls.eq("R_down").to_numpy() & el
    band = ex.held_out_band.to_numpy().copy()
    band[f] = np.digitize(1 - ex.held_out_baseline.to_numpy()[f], BAND_EDGES) - 1
    cell = np.array([f"{b_}:{t_}" for b_, t_ in zip(band, ex.held_out_tertile)])
    lab, kept = coarsen(band, ex.held_out_tertile, f, c, per_first=1)
    rows.append({"contrast": "instrument check (R_down drawn to R_up's mirrored baseline)", "side": "held_out",
                 "first_eligible": int(f.sum()), "control_eligible": int(c.sum()),
                 "first_kept_band_x_tertile_uncoarsened": int(sum(   # cells with enough R_down to draw from
                     (cell[f] == k).sum() for k in set(cell[f]) if (cell[c] == k).sum() >= (cell[f] == k).sum())),
                 **effective_n(lab, f, c, kept)})
    need = pd.Series(band[f]).value_counts().sort_index()
    have = pd.Series(ex.held_out_band.to_numpy()[c]).value_counts().reindex(need.index).fillna(0).astype(int)
    log("instrument check, R_up exons per mirrored band and eligible R_down exons in that band:\n"
        + pd.DataFrame({"R_up_mirrored": need, "R_down_eligible": have}).to_string())
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "strata_contrasts.csv", index=False, float_format="%.1f")
    pd.set_option("display.width", 250)
    log("class contrasts after the coarsening rule:\n" + out.to_string(index=False))
    rp5a_vif(ex)


def vif(X: pd.DataFrame) -> pd.Series:
    """Variance inflation factor of each column (OLS of the column on the others, with intercept)."""
    out = {}
    A = X.to_numpy(float)
    for j, col in enumerate(X.columns):
        y, Z = A[:, j], np.column_stack([np.ones(len(A)), np.delete(A, j, axis=1)])
        r = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
        r2 = 1 - r.var() / y.var() if y.var() > 0 else 1.0
        out[col] = 1 / max(1 - r2, 1e-12)
    return pd.Series(out)


def rp5a_vif(ex: pd.DataFrame) -> None:
    """Label-blind collinearity of RP5a's exon covariates (the host-gene change, which uses the
    label, left out): the primary model (core 3'ss score, exon in an Alu) and the full model (also
    the -20..+3 PWM and the acceptor Alu), over eligible held-out R_up and B-sample exons."""
    el = ex.held_out_eligible & (ex.cls.eq("R_up") | (ex.b_sample.eq(True) & ex.cls.isin(["B_far", "B_near"])))
    d = ex[el].copy()
    X = pd.DataFrame({"R_up": d.cls.eq("R_up").astype(float),
                      "log_exon_length": np.log(d.exon_len.astype(float)),
                      "log_median_reads": np.log(d.held_out_reads.astype(float) + 1),
                      "poison": d.consequence.eq("ptc_nmd").astype(float),
                      "core_3ss": d.pwm_core.astype(float), "exon_gc": d.exon_gc.astype(float),
                      "exon_alu": d.exon_alu.eq(True).astype(float)}, index=d.index)
    bands = pd.get_dummies(d.held_out_band, prefix="band", drop_first=True).astype(float)
    full = X.assign(pwm_20_3=d.pwm3.astype(float), acceptor_alu=d.acc_alu.eq(True).astype(float))
    rows = []
    for name, M_ in [("primary", X), ("full (sensitivity)", full)]:
        M_ = pd.concat([M_, bands], axis=1).dropna()
        v = vif(M_.loc[:, M_.std() > 0])
        for col in [c for c in M_.columns if not c.startswith("band_")]:
            rows.append({"model": name, "covariate": col, "vif": float(v.get(col, np.nan)), "n_exons": len(M_)})
    r = np.corrcoef(X.R_up, X.exon_alu)[0, 1]
    rows.append({"model": "pairwise", "covariate": "exon_alu with R_up", "vif": float(1 / (1 - r ** 2)), "n_exons": len(X)})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "rp5a_vif.csv", index=False, float_format="%.3f")
    log("RP5a covariate VIFs (label-blind):\n" + out.to_string(index=False))


# ================================================================ Phase 2: the tests (D-77)
# Everything below uses the MMR label. `--labels shuffled` (the dry run) moves every label among the
# tumours of each cancer (seed 999) before anything else and writes only under .../proof/dryrun/;
# the real labels need `--labels real --real-labels-approved`, given only after Gate 2.

N_PERM = 2000                  # inference.label_perm
N_CLASS_PERM = 10_000          # class_perm, and every exon, tumour or pair permutation set at 10,000
N_FL_PERM = 10_000             # Freedman-Lane permutations: not fixed by the block (proposed at Gate 2)
N_SETS = 1000                  # RP1's matched sets
N_DRAWS = 1000                 # instrument check
N_COUNT_PERM = 200             # count model, label permutations (models refitted)
RPD1_PERM = 1000               # A-22: RPD1's permutation null (X rows permuted within fold)
SEP_LIMIT = 10.0               # A-23: a count-model fit with |coefficient| > 10 or SE > 10 (logit) is separated
DRY_SEED = 999
MIN_DELTA_MMR, MIN_ARM, MIN_NORM = 5, 5, 10   # delta, arm and adj minimums (quantities)
MIN_CLASS, MIN_TUMOUR, MIN_RP3D, MIN_RP3E, MIN_SHAP = 20, 10, 15, 40, 10   # inference.min_n
SESOI = 0.003
INSTR_PASS = 0.8
RP6C_DRAWS, RP6C_GROUP, RP6C_MIN_CLASS, RP6C_MIN_MEAS, RP6C_NULL = 50, 20, 20, 10, 57
M4_AUROC, M4_TOL = 0.900, 0.001
M4_RANDOM_DRAWS = 200          # A-25: random removals of as many non-R_up repeat inputs (seed 0)
ALPHA = 0.05
LABEL_COLS = ["n_hp", "n_snv", "mmr", "hyper", "mantis", "msisensor", "call_mantis", "call_msisensor"]
L_BINS = {"RP3c": [20, 25, 30], "RP3e": [20, 25], "RPM3": [20, 25]}     # left edges; last bin open
O_BINS = {"RP3c": [3, 5, 10], "RP3e": [3, 5]}
D_BINS = [20, 30, 40, 50, 70]                                            # A-16's branch-point bins (descriptive)
EXTRA = PROOF / "extra"                                                  # ppt_far >= 15 members not stored by `counts`


class Run:
    """Where a compute or report run writes. The dry run (shuffled labels) may write only under
    results/.../proof/dryrun/ and data.nosync/.../proof/dryrun/; the real-label run needs an explicit
    approval flag and never writes into a dryrun folder."""

    def __init__(self, labels: str, approved: bool = False):
        if labels not in ("shuffled", "real"):
            raise SystemExit("--labels must be 'shuffled' (the dry run) or 'real'")
        if labels == "real" and not approved:
            raise SystemExit("the real MMR labels need --labels real --real-labels-approved (after Gate 2); "
                             "the dry run is --labels shuffled")
        self.labels, self.dry = labels, labels == "shuffled"
        self.out = OUT / "dryrun" if self.dry else OUT
        self.local = PROOF / "dryrun" if self.dry else PROOF / "tests"

    def allowed(self, p: Path) -> bool:
        p = Path(p).resolve()
        if self.dry:
            return any(p.is_relative_to(b.resolve()) for b in (OUT / "dryrun", PROOF / "dryrun"))
        return "dryrun" not in p.parts and any(p.is_relative_to(b.resolve()) for b in (OUT, PROOF / "tests"))

    def path(self, name: str, local: bool = False) -> Path:
        p = (self.local if local else self.out) / name
        if not self.allowed(p):
            raise RuntimeError(f"a {self.labels}-label run may not write {p}")
        p.parent.mkdir(parents=True, exist_ok=True)
        return p


# ---------------------------------------------------------------- labels

def pair_rows(s: pd.DataFrame) -> pd.DataFrame:
    """The pair rule on one cancer's samples (rows = the count arrays' rows): each tumour with a Solid
    Tissue Normal file of the same case, the first normal file id in sort order if several."""
    tum = s[s.role.eq("tumour")].reset_index()
    nm = s[s.role.eq("normal")].sort_values("file_id").drop_duplicates("case_id").reset_index()
    m = tum.merge(nm, on="case_id", suffixes=("_t", "_n"))
    return pd.DataFrame({"case_id": m.case_id, "tumour_file_id": m.file_id_t, "normal_file_id": m.file_id_n,
                         "tumour_row": m.index_t, "normal_row": m.index_n})


def pairs_table() -> pd.DataFrame:
    """Label-blind tumour-normal pairs of every cancer (pair_rows), with their rows in the count arrays."""
    return pd.concat([pair_rows(pd.read_parquet(f)).assign(project=f.name.split("__")[0])
                      for f in sorted(PROOF.glob("*__samples.parquet"))], ignore_index=True)


def shuffle_labels(t: pd.DataFrame, seed: int = DRY_SEED) -> pd.DataFrame:
    """The dry run's labels: all label columns moved together among the tumours of each cancer, within
    tumours with and without a normal pair separately (so each cancer keeps its MMRd pair count)."""
    rng = np.random.default_rng(seed)
    out = t.copy()
    for _, idx in sorted(t.groupby(["project", "has_pair"]).indices.items()):
        src = idx[rng.permutation(len(idx))]
        for c in LABEL_COLS:
            col = out[c].to_numpy().copy()
            col[idx] = t[c].to_numpy()[src]
            out[c] = col
    return out


def tumour_labels(run: Run) -> pd.DataFrame:
    """D-71's tumours with every label the tests use (MC3 status, hypermutation, MANTIS and MSIsensor);
    shuffled within cancer first in the dry run."""
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")[["file_id", "case_id", "sample15", "project",
                                                          "n_hp", "n_snv", "mmr", "hyper"]]
    c = pd.read_parquet(LOCAL / "pancancer_msicalls.parquet")[["file_id", "mantis", "msisensor",
                                                              "call_mantis", "call_msisensor"]]
    t = t.merge(c, on="file_id", how="left")
    t["has_pair"] = t.file_id.isin(pairs_table().tumour_file_id)
    t["held_out"] = ~t.project.isin(TRAIN)
    if run.dry:
        t = shuffle_labels(t, DRY_SEED)
    t["mmr"] = t.mmr.astype(bool)
    t["hyper"] = t.hyper.astype(bool)
    return t


def label_perms(y, groups, n_perm: int = N_PERM, seed: int = SEED) -> np.ndarray:
    """[n x (n_perm + 1)]: column 0 the observed labels, then n_perm permutations of them among the
    members of each group (counts kept), drawn permutation by permutation over the groups in sorted
    order - so the first k columns of a longer run equal a run of k (the count model's 200)."""
    y, groups = np.asarray(y), np.asarray(groups)
    out = np.empty((len(y), n_perm + 1), dtype=y.dtype)
    out[:, 0] = y
    idx = [np.flatnonzero(groups == g) for g in sorted(set(groups.tolist()))]
    rng = np.random.default_rng(seed)
    for k in range(1, n_perm + 1):
        for ii in idx:
            out[ii, k] = y[ii[rng.permutation(len(ii))]]
    return out


def perm_p(obs: float, null) -> float:
    """(1 + number of |null| >= |observed|) / (1 + permutations); undefined permutations dropped."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs):
        return np.nan
    return float((1 + np.sum(np.abs(null) >= abs(obs) - 1e-12)) / (len(null) + 1))


def iut(*ps) -> float:
    """inference.intersection_union: the larger of the available p values (exon or class permutation,
    label permutation)."""
    ok = [float(x) for x in ps if x is not None and np.isfinite(x)]
    return max(ok) if ok else np.nan


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q over the finite p (NaN stays NaN)."""
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        o = np.argsort(pv)
        r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        out = np.empty(len(pv))
        out[o] = np.minimum(r, 1)
        q[ok] = out
    return q


# ---------------------------------------------------------------- per-cancer quantities (label-blind)

class Spec(NamedTuple):
    """One read-out: the view, the measured-cell minimum, the universe share and the units."""
    name: str
    view: str = "jc"
    min_reads: int = MIN_READS
    share: float = UNIVERSE_SHARE
    units: str = "psi"


SPECS = {"primary": Spec("primary"), "jcec": Spec("jcec", view="jcec"), "reads30": Spec("reads30", min_reads=30),
         "univ80": Spec("univ80", share=0.8), "pct": Spec("pct", units="pct")}


def ppt_far_at(a: pd.DataFrame, length: int) -> pd.Series:
    """ctl_ppt_far with the pyrimidine-run cut at `length` (sensitivity: 15 and 25 nt)."""
    dp = np.nan_to_num(a.d_bp_ppt.to_numpy(float), nan=-1)
    return pd.Series(a.struct.to_numpy(bool) & (a.acc_T.to_numpy(float) <= CTRL_MAX_T)
                     & (a.acc_C.to_numpy(float) <= PPT_MAX_HOMO) & (a.ppt_len.to_numpy(float) >= length)
                     & (a.ppt_max_T.to_numpy(float) <= PPT_MAX_HOMO) & (a.ppt_max_C.to_numpy(float) <= PPT_MAX_HOMO)
                     & (a.ppt_alt.to_numpy(float) < PPT_MAX_ALT) & (dp >= FAR), index=a.index)


def counts_extra() -> None:
    """Label-blind: JC and JCEC counts of the ppt_far >= 15 nt members `counts` did not store (the
    sensitivity analysis needs every member), in the same sample order, under proof/extra/."""
    from joblib import Parallel, delayed
    store = Path(expand("${CTSM_SPLICE_COUNTS}"))
    a = load_annotation()
    stored = set(pd.read_parquet(PROOF / "exons.parquet").feature_id)
    x = a[ppt_far_at(a, 15) & ~a.feature_id.isin(stored)].copy()
    ev = pd.read_csv(store / "events" / "se.events.tsv.gz", sep="\t", usecols=["event_idx", "id"], dtype={"id": str})
    x["event_idx"] = ev.set_index("id").event_idx.loc[x.feature_id].to_numpy()
    x["b_sample"], x["stored_as"] = False, "sensitivity extra"
    EXTRA.mkdir(parents=True, exist_ok=True)
    x.to_parquet(EXTRA / "exons.parquet", index=False)
    jobs = [(f.name.split("__")[0], pd.read_parquet(f)) for f in sorted(PROOF.glob("*__samples.parquet"))]
    qc = Parallel(n_jobs=N_JOBS)(delayed(_count_cancer)(p, s, x.event_idx.to_numpy(np.int64), str(store), str(EXTRA))
                                 for p, s in jobs)
    for p, s in jobs:                                   # same rows as the stored arrays
        assert pd.read_parquet(EXTRA / f"{p}__samples.parquet").file_id.equals(s.file_id)
    log(f"extra counts: {len(x)} ppt_far >= 15 exons; max |PSI - store| "
        f"{max(q['jc_max_abs_psi_diff'] for q in qc):.2g} (JC), {max(q['jcec_max_abs_psi_diff'] for q in qc):.2g} (JCEC)")


def exon_table() -> pd.DataFrame:
    """Every stored exon (proof/exons.parquet, in column order) and, after them, the extra ppt_far >= 15
    exons, with the annotation, the classes and the sensitivity class flags."""
    a = load_annotation().set_index("feature_id")
    parts = [pd.read_parquet(PROOF / "exons.parquet")[["feature_id", "stored_as", "b_sample"]]]
    if (EXTRA / "exons.parquet").exists():
        parts.append(pd.read_parquet(EXTRA / "exons.parquet")[["feature_id", "stored_as", "b_sample"]])
    ex = pd.concat(parts, ignore_index=True)
    ex = ex.join(a.drop(columns=["stored_as", "b_sample"], errors="ignore"), on="feature_id")
    ex["b_sample"] = ex.b_sample.astype(bool)
    ex["ctl_ppt_far15"], ex["ctl_ppt_far25"] = ppt_far_at(ex, 15).to_numpy(), ppt_far_at(ex, 25).to_numpy()
    ex["g_exp_near"] = ex.d_exp_near - ex.run_far
    return ex


def load_counts_all(proj: str, view: str):
    """(samples, IJC, SJC, ratio) over the stored and the extra exons (extra columns last)."""
    samp, ijc, sjc, ratio = load_counts(proj, view, mmap=True)
    f = EXTRA / f"{proj}__{view}__ijc.npy"
    if not f.exists():
        return samp, ijc, sjc, ratio
    ext = [np.load(EXTRA / f"{proj}__{view}__{k}.npy", mmap_mode="r") for k in ["ijc", "sjc", "ratio"]]
    return (samp, *(np.concatenate([np.asarray(a_), np.asarray(b_)], axis=1) for a_, b_ in zip([ijc, sjc, ratio], ext)))


def band_of(x) -> np.ndarray:
    return np.digitize(np.nan_to_num(np.asarray(x, float)), BAND_EDGES) - 1


def band_adjust(dev: np.ndarray, band: np.ndarray, bmask: np.ndarray, min_n: int = MIN_NORM):
    """adj = dev minus the tumour's mean dev over its measured B-sample exons of the exon's band (the
    band and its two neighbours when fewer than min_n are measured; the cell is dropped when fewer
    still); also each tumour's mean dev over all its measured B-sample exons (the global axis)."""
    nb = len(BAND_EDGES) - 1
    bi = np.flatnonzero(bmask)
    Bd = dev[:, bi].astype(np.float64)
    fin = np.isfinite(Bd)
    B0 = np.where(fin, Bd, 0.0)
    oh = np.zeros((len(bi), nb))
    oh[np.arange(len(bi)), band[bi]] = 1
    S, N = B0 @ oh, fin.astype(float) @ oh
    S3, N3 = S.copy(), N.copy()
    S3[:, 1:] += S[:, :-1]
    S3[:, :-1] += S[:, 1:]
    N3[:, 1:] += N[:, :-1]
    N3[:, :-1] += N[:, 1:]
    with np.errstate(invalid="ignore", divide="ignore"):
        ref = np.where(N >= min_n, S / np.maximum(N, 1), np.where(N3 >= min_n, S3 / np.maximum(N3, 1), np.nan))
        ga = np.where(fin.sum(1) > 0, B0.sum(1) / np.maximum(fin.sum(1), 1), np.nan)
    adj = np.full(dev.shape, np.nan, np.float32)
    cols = np.flatnonzero(band >= 0)
    adj[:, cols] = dev[:, cols] - ref[:, band[cols]]
    return adj, ga


def class_arm(x: np.ndarray, mask: np.ndarray, min_n: int = MIN_ARM) -> np.ndarray:
    """Per row, the mean over the measured exons of `mask` (missing below min_n)."""
    v = x[:, np.flatnonzero(mask)].astype(np.float64)
    fin = np.isfinite(v)
    n = fin.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n >= min_n, np.where(fin, v, 0).sum(1) / np.maximum(n, 1), np.nan)


def pct_within(x: np.ndarray) -> np.ndarray:
    """Mid-rank percentile (rank / number finite) of each column's finite values; NaN stays NaN."""
    from scipy.stats import rankdata
    x = np.asarray(x, float)
    r = rankdata(x, axis=0, nan_policy="omit")
    n = np.isfinite(x).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(np.isfinite(x), r / np.maximum(n, 1), np.nan)


def pct_vector(v: np.ndarray) -> np.ndarray:
    return pct_within(np.asarray(v, float)[:, None])[:, 0]


def matched_draws(fb, ft, pb, pt, n_sets: int, rng, log: dict | None = None) -> np.ndarray:
    """[n_sets x first exons] indices into the pool: each first exon replaced by a pool exon of the same
    band and tertile, without replacement within a set where the cell holds enough; otherwise from the
    band's tertiles merged (with replacement if still too few), then from neighbouring bands."""
    fb, ft, pb, pt = (np.asarray(v, int) for v in (fb, ft, pb, pt))
    out = np.empty((n_sets, len(fb)), int)
    for key in sorted(set(zip(fb.tolist(), ft.tolist()))):
        fi = np.flatnonzero((fb == key[0]) & (ft == key[1]))
        cand = np.flatnonzero((pb == key[0]) & (pt == key[1]))
        rule = "band x tertile"
        if len(cand) < len(fi):
            cand = np.flatnonzero(pb == key[0])
            rule = "tertiles merged"
        w = 1
        while len(cand) == 0 and w < len(BAND_EDGES):
            cand = np.flatnonzero(np.abs(pb - key[0]) <= w)
            rule = "neighbouring bands"
            w += 1
        if len(cand) >= len(fi):
            pick = rng.random((n_sets, len(cand))).argsort(axis=1)[:, :len(fi)]
        else:
            pick = rng.integers(0, len(cand), (n_sets, len(fi)))
            rule += ", with replacement"
        if log is not None:
            log[rule] = log.get(rule, 0) + len(fi)
        out[:, fi] = cand[pick]
    return out


def set_arms(dev: np.ndarray, first: np.ndarray, pool: np.ndarray, band, terc, n_sets: int, rng,
             log: dict | None = None) -> np.ndarray:
    """RP1: [tumours x (1 + n_sets)] mid-rank percentiles within the cancer of the dev-arm over the
    first-class exons (column 0) and over each matched set of pool exons."""
    draws = matched_draws(band[first], terc[first], band[pool], terc[pool], n_sets, rng, log)
    cnt = np.zeros((len(pool), n_sets), np.float64)
    np.add.at(cnt, (draws.ravel(), np.repeat(np.arange(n_sets), len(first))), 1)
    D = dev[:, pool].astype(np.float64)
    fin = np.isfinite(D)
    s, c = np.where(fin, D, 0) @ cnt, fin.astype(float) @ cnt
    with np.errstate(invalid="ignore", divide="ignore"):
        sets = np.where(c >= MIN_ARM, s / np.maximum(c, 1), np.nan)
    obs = class_arm(dev, np.isin(np.arange(dev.shape[1]), first))
    return pct_within(np.column_stack([obs, sets])).astype(np.float32)


def build_cancer(proj: str, spec: Spec, flags: dict, workdir: str, k: int, rp1: bool, pairs: pd.DataFrame | None) -> dict:
    """Label-blind, one cancer, one read-out: universe, m and its band, adj (saved as float32 to
    workdir), the R_up and R_down arms, the global axis, RP1's matched-set percentiles and the
    tumour-normal pair arms."""
    samp, ijc, sjc, ratio = load_counts_all(proj, spec.view)
    tum = samp.role.eq("tumour").to_numpy()
    I, S_ = np.asarray(ijc[tum], np.float64), np.asarray(sjc[tum], np.float64)
    R = np.asarray(ratio[tum], np.float64)
    n = I + S_
    meas = n >= spec.min_reads
    psi = np.where(meas, psi_from_counts(I, S_, R), np.nan)          # float64, as Phase 1 (band edges)
    U = meas.mean(0) >= spec.share
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        m = np.where(U, np.nanmedian(np.where(U, psi, np.nan), axis=0), np.nan)
        band = np.where(U, band_of(m), -1)
        X = psi
        if spec.units == "pct":
            X = np.full(psi.shape, np.nan)
            X[:, U] = pct_within(psi[:, U])
        mx = np.where(U, np.nanmedian(np.where(U, X, np.nan), axis=0), np.nan) if spec.units == "pct" else m
        dev = np.where(U, X - mx, np.nan)
    adj, ga = band_adjust(dev, band, flags["b_sample"] & U)
    np.save(Path(workdir) / f"{proj}.npy", adj)
    reads = np.median(n, axis=0).astype(np.float32)
    res = {"project": proj, "file_id": samp.file_id[tum].to_numpy(), "U": U, "m": m, "band": band.astype(np.int8),
           "arm_up": class_arm(adj, flags["R_up"] & U), "arm_down": class_arm(adj, flags["R_down"] & U),
           "ga": ga, "ga_pct": pct_vector(ga), "n_cols": adj.shape[1]}
    if rp1:
        uni = np.flatnonzero(U)
        terc = np.digitize(reads, np.quantile(reads[uni], [1 / 3, 2 / 3])) if len(uni) else np.zeros(len(U), int)
        pool = np.flatnonzero(flags["b_sample"] & U)
        for arm, key in [("R_up", 1), ("R_down", 2)]:
            first = np.flatnonzero(flags[arm] & U)
            if len(first) and len(pool):
                rng = np.random.default_rng([SEED, k, key])
                res[f"rp1_{arm}_rules"] = {}
                res[f"rp1_{arm}"] = set_arms(dev, first, pool, band, terc, N_SETS, rng, res[f"rp1_{arm}_rules"])
    if pairs is not None and len(pairs):
        rows = []
        Ia, Sa, Ra = np.asarray(ijc, np.float64), np.asarray(sjc, np.float64), np.asarray(ratio, np.float64)
        P_all = np.where(Ia + Sa >= spec.min_reads, psi_from_counts(Ia, Sa, Ra), np.nan)
        for r in pairs.itertuples():
            d = np.where(U, P_all[r.tumour_row] - P_all[r.normal_row], np.nan)[None, :].astype(np.float32)
            pa, _ = band_adjust(d, band, flags["b_sample"] & U)
            rows.append({"case_id": r.case_id, "pair_up": class_arm(pa, flags["R_up"] & U)[0],
                         "pair_down": class_arm(pa, flags["R_down"] & U)[0]})
        res["pairs"] = pd.DataFrame(rows)
    return res


def pooled_reads(spec: Spec, projects: list, U: dict, cols: np.ndarray, chunk: int = 2000) -> np.ndarray:
    """Per column in `cols`: the median reads (IJC + SJC, the spec's view) over all tumours of the
    listed cancers whose universe holds the exon."""
    out = np.full(len(cols), np.nan)
    arrays = {}
    for p in projects:
        samp, ijc, sjc, _ = load_counts_all(p, spec.view)
        arrays[p] = (samp.role.eq("tumour").to_numpy(), ijc, sjc)
    for lo in range(0, len(cols), chunk):
        cc = cols[lo:lo + chunk]
        pool = [[] for _ in cc]
        for p in projects:
            tum, ijc, sjc = arrays[p]
            sel = np.flatnonzero(U[p][cc])
            if not len(sel):
                continue
            nn = (np.asarray(ijc[:, cc[sel]], np.int64) + np.asarray(sjc[:, cc[sel]], np.int64))[tum]
            for j, k in enumerate(sel):
                pool[k].append(nn[:, j])
        for k, v in enumerate(pool):
            if v:
                out[lo + k] = float(np.median(np.concatenate(v)))
    return out


def side_tables(spec: Spec, res: dict, n_mmr: pd.Series, thresholds=(10, 20, 30), n_stored: int | None = None) -> dict:
    """Label-blind, per side: eligibility (held out at each MMRd threshold; training = in a training
    universe), baseline, its band, the pooled median reads and their tertile, with the tertile cuts over
    the side's eligible stored exons (the first n_stored columns; extra exons take the same cuts).
    Returns {(side, threshold): DataFrame over all columns}."""
    projects = sorted(res)
    E = res[projects[0]]["n_cols"]
    U = {p: res[p]["U"] for p in projects}
    out = {}
    for side, ps, thr_list in [("held_out", [p for p in projects if p not in TRAIN], thresholds),
                               ("training", list(TRAIN), (1,))]:
        if not ps or not set(ps) <= set(projects):          # a held-out-only read-out has no training side
            continue
        w = np.array([[n_mmr[p]] for p in ps], float) * np.vstack([U[p] for p in ps])
        mm = w.sum(0)
        with np.errstate(all="ignore"):
            base = np.nansum(np.vstack([np.nan_to_num(res[p]["m"]) for p in ps]) * w, 0) / mm
        lo = min(thr_list)
        cand = np.flatnonzero(mm >= lo)
        reads = np.full(E, np.nan)
        reads[cand] = pooled_reads(spec, ps, U, cand)
        stored = np.arange(E) < (n_stored if n_stored is not None else E)
        for thr in thr_list:
            el = mm >= thr
            q = np.nanquantile(reads[el & stored], [1 / 3, 2 / 3])
            out[(side, thr)] = pd.DataFrame({"eligible": el, "mmr_in_universe": mm,
                                             "baseline": np.where(el, base, np.nan),
                                             "band": np.where(el, band_of(base), -1),
                                             "reads": np.where(el, reads, np.nan),
                                             "tertile": np.where(el, np.digitize(np.nan_to_num(reads), q), -1)})
    return out


# ---------------------------------------------------------------- statistics (unit-tested on synthetic data)

def mmrd_index(Y: np.ndarray) -> np.ndarray:
    """[columns x MMRd count]: the rows of the MMRd tumours in each label column (the count is fixed)."""
    nd = int(Y[:, 0].sum())
    return np.argsort(~Y, axis=0, kind="stable")[:nd].T


def exon_engine(adjs: list, Ys: list, min_mmr: int = MIN_DELTA_MMR, block_elems: float = 1.5e7):
    """Per exon and label column: delta, its studentized value and the measured MMRd count N.
    adjs - per cancer [tumours x exons] adj (NaN = not measured); Ys - per cancer [tumours x columns]
    MMRd labels. delta = the mean over measured MMRd tumours (of cancers holding a measured MMRp
    tumour) of adj minus the cancer's measured-MMRp mean, NaN when N < min_mmr. Studentized as in
    inference.studentized, each measured MMRd tumour contributing adj / N and each measured MMRp
    tumour adj x n_D,c / (N n_P,c): SE^2 = [SSD x N / (N - K) + VP] / N^2 with SSD the within-cancer
    sum of squares of MMRd adj and VP = sum_c n_D,c^2 / n_P,c x var_P,c(adj)."""
    P, E = Ys[0].shape[1], adjs[0].shape[1]
    num, SSD, VP = np.zeros((P, E)), np.zeros((P, E)), np.zeros((P, E))
    N, K = np.zeros((P, E), np.int32), np.zeros((P, E), np.int16)
    for A, Y in zip(adjs, Ys):
        nd = int(Y[:, 0].sum())
        if nd == 0:
            continue
        M = np.isfinite(A)
        A0 = np.where(M, A, 0).astype(np.float32)
        tot_s, tot_n = A0.sum(0, dtype=np.float64), M.sum(0)
        tot_q = np.square(A0.astype(np.float64)).sum(0)
        idx = mmrd_index(Y)
        b = max(1, int(block_elems // (nd * E)))
        for p0 in range(0, P, b):
            ii = idx[p0:p0 + b]
            G = A0[ii].astype(np.float64)
            sD, qD, nD = G.sum(1), np.square(G).sum(1), M[ii].sum(1)
            sP, nP, qP = tot_s - sD, tot_n - nD, tot_q - qD
            ok = (nD > 0) & (nP > 0)
            sl = slice(p0, p0 + len(ii))
            with np.errstate(all="ignore"):
                num[sl] += np.where(ok, sD - nD * sP / np.maximum(nP, 1), 0)
                N[sl] += np.where(ok, nD, 0).astype(np.int32)
                K[sl] += ok
                SSD[sl] += np.where(ok, np.maximum(qD - sD ** 2 / np.maximum(nD, 1), 0), 0)
                varP = np.where(nP > 1, np.maximum(qP - sP ** 2 / np.maximum(nP, 1), 0) / np.maximum(nP - 1, 1), 0)
                VP[sl] += np.where(ok, nD ** 2 / np.maximum(nP, 1) * varP, 0)
    with np.errstate(all="ignore"):
        delta = np.where(N >= min_mmr, num / np.maximum(N, 1), np.nan)
        se2 = (SSD * N / np.maximum(N - K, 1) + VP) / np.maximum(N, 1).astype(float) ** 2
        t = np.where((N >= min_mmr) & (N > K) & (se2 > 0), delta / np.sqrt(se2), np.nan)
    return delta.astype(np.float32), t.astype(np.float32), N


def studentized_sum(Us: list, Vs: list, Yds: list, Yps: list):
    """inference.studentized from per-tumour contributions. Per cancer: U, V [tumours x columns] the
    contribution of each tumour if MMRd (u) or MMRp (v); Yd, Yp the MMRd and MMRp tumours contributing
    in each column. T = sum_D u - sum_P v; SE^2 = N s2_D + sum_c n_P,c s2_P,c, s2_D = sum_c sum_D (u -
    mean_c u)^2 / (N - K) (K = cancers holding an MMRd tumour), s2_P,c the variance (ddof 1) of v over
    cancer c's MMRp tumours (cancers holding an MMRd tumour). Returns T, SE2, N, K and column 0's
    variance components with their df for the Welch-Satterthwaite interval."""
    C = Us[0].shape[1]
    T, ssd, ND, K, VP = (np.zeros(C) for _ in range(5))
    comp = []
    for U, V, Yd, Yp in zip(Us, Vs, Yds, Yps):
        nd, npp = Yd.sum(0), Yp.sum(0)
        su, su2 = (U * Yd).sum(0), (U * U * Yd).sum(0)
        sv, sv2 = (V * Yp).sum(0), (V * V * Yp).sum(0)
        has = nd > 0
        with np.errstate(all="ignore"):
            ssd += np.where(has, np.maximum(su2 - su ** 2 / np.maximum(nd, 1), 0), 0)
            vp = np.where(has & (npp > 1), npp * np.maximum(sv2 - sv ** 2 / np.maximum(npp, 1), 0)
                          / np.maximum(npp - 1, 1), 0)
        ND += nd
        K += has
        VP += vp
        T += np.where(has, su - sv, 0)
        if has[0] and npp[0] > 1:
            comp.append((float(vp[0]), int(npp[0]) - 1))
    with np.errstate(all="ignore"):
        s2d = np.where(ND > K, np.maximum(ssd, 0) / np.maximum(ND - K, 1), np.nan)
    SE2 = ND * s2d + VP
    comps = {"V_mmrd": float(ND[0] * s2d[0]), "df_mmrd": int(ND[0] - K[0]), "mmrp": comp}
    return T, SE2, ND, K, comps


def ws_df(se2: float, comps: dict) -> float:
    """Welch-Satterthwaite df, SE^4 / sum_k V_k^2 / df_k, over the MMRd component (df N - K) and each
    cancer's MMRp component (df n_P,c - 1)."""
    den = (comps["V_mmrd"] ** 2 / comps["df_mmrd"] if comps["df_mmrd"] > 0 else 0.0)
    den += sum(v ** 2 / d for v, d in comps["mmrp"] if d > 0)
    return float(se2 ** 2 / den) if den > 0 else np.nan


def interval(est: float, se2: float, comps: dict, level: float = 0.95) -> tuple[float, float, float]:
    """inference.interval: estimate +/- t x SE at the Welch-Satterthwaite df; returns (lo, hi, df)."""
    from scipy.stats import t as tdist
    df = ws_df(se2, comps)
    q = tdist.ppf(0.5 + level / 2, df) if np.isfinite(df) else np.nan
    se = float(np.sqrt(se2))
    return est - q * se, est + q * se, df


def equal_weights(valid: np.ndarray, sl: slice | None = None) -> np.ndarray:
    return valid / np.maximum(valid.sum(1, keepdims=True), 1)


def gene_weights(genes) -> callable:
    """Weights for 'delta averaged within gene before the class mean': each gene with a valid exon
    gets 1 / (genes), split equally over its valid exons."""
    g = pd.factorize(pd.Series(genes).fillna("none").to_numpy())[0]
    G = np.zeros((len(g), g.max() + 1))
    G[np.arange(len(g)), g] = 1

    def fn(valid, sl=None):
        per_gene = valid.astype(float) @ G                       # valid exons per gene
        n_genes = (per_gene > 0).sum(1, keepdims=True)
        with np.errstate(all="ignore"):
            w = valid / np.maximum(per_gene @ G.T, 1) / np.maximum(n_genes, 1)
        return np.where(valid, w, 0)
    return fn


def class_stat(As: list, Ys: list, weight_fn=None, min_mmr: int = MIN_DELTA_MMR, block: int = 400) -> dict:
    """A class statistic sum_e w_e delta_e (default the mean delta over the class exons with a delta in
    that column) written as a sum of per-tumour contributions and studentized (inference.studentized,
    N = every MMRd tumour of the cancers given). As - per cancer [tumours x class exons] adj; Ys - per
    cancer [tumours x columns] MMRd labels (label permutations, or the observed labels repeated across
    instrument draws, with weight_fn selecting each draw's exons). weight_fn(valid [columns x exons],
    column slice) returns weights summing to 1 over the valid exons of each column."""
    weight_fn = weight_fn or equal_weights
    C = Ys[0].shape[1]
    out = {k: np.zeros(C) for k in ("T", "SE2", "N", "K", "n_exons")}
    comps0 = None
    Ms = [np.isfinite(A) for A in As]
    A0s = [np.where(M, A, 0.0) for A, M in zip(As, Ms)]
    for c0 in range(0, C, block):
        sl = slice(c0, min(C, c0 + block))
        nDs, nPs, N = [], [], 0
        for M, Y in zip(Ms, Ys):
            nD = Y[:, sl].T.astype(np.float64) @ M
            nP = M.sum(0)[None, :] - nD
            nD = np.where((nD > 0) & (nP > 0), nD, 0)
            nDs.append(nD)
            nPs.append(nP)
            N = N + nD
        valid = N >= min_mmr
        w = np.asarray(weight_fn(valid, sl), float)
        with np.errstate(all="ignore"):
            w = np.where(valid, w / np.maximum(N, 1), 0)
        Us, Vs, Yds, Yps = [], [], [], []
        for A0, Y, nD, nP in zip(A0s, Ys, nDs, nPs):
            Wd = np.where(nD > 0, w, 0)
            with np.errstate(all="ignore"):
                Wp = np.where(nD > 0, w * nD / np.maximum(nP, 1), 0)
            Us.append(A0 @ Wd.T)
            Vs.append(A0 @ Wp.T)
            y = Y[:, sl].astype(bool)
            Yds.append(y)
            Yps.append(~y)
        T, SE2, ND, K, comps = studentized_sum(Us, Vs, Yds, Yps)
        for k, v in zip(("T", "SE2", "N", "K", "n_exons"), (T, SE2, ND, K, valid.sum(1))):
            out[k][sl] = v
        if c0 == 0:
            comps0 = comps
    with np.errstate(all="ignore"):
        out["z"] = np.where(out["SE2"] > 0, out["T"] / np.sqrt(out["SE2"]), np.nan)
    out["comps"] = comps0
    return out


def arm_stat(values, Y: np.ndarray, groups, block: int = 500) -> dict:
    """arm_contrast - the mean over MMRd tumours of (arm minus the mean arm of the MMRp tumours of the
    same cancer) - for each label column, studentized; tumours without an arm, and MMRd tumours of a
    cancer with no MMRp arm, leave the column."""
    v = np.asarray(values, float)
    ok = np.isfinite(v)
    v0 = np.where(ok, v, 0)
    groups = np.asarray(groups)
    idx = [np.flatnonzero(groups == g) for g in sorted(set(groups.tolist()))]
    C = Y.shape[1]
    out = {k: np.zeros(C) for k in ("T", "SE2", "N", "K")}
    comps0 = None
    for c0 in range(0, C, block):
        sl = slice(c0, min(C, c0 + block))
        Yb = Y[:, sl].astype(bool)
        nds = [(Yb[ii] & ok[ii, None]).sum(0) for ii in idx]
        nps = [(~Yb[ii] & ok[ii, None]).sum(0) for ii in idx]
        va = [(a > 0) & (b > 0) for a, b in zip(nds, nps)]
        ND = sum(np.where(x, a, 0) for x, a in zip(va, nds))
        Us, Vs, Yds, Yps = [], [], [], []
        for ii, nd, npp, x in zip(idx, nds, nps, va):
            with np.errstate(all="ignore"):
                Us.append(v0[ii, None] / np.maximum(ND, 1)[None, :])
                Vs.append(v0[ii, None] * np.where(x, nd / np.maximum(ND * npp, 1), 0)[None, :])
            Yds.append(Yb[ii] & ok[ii, None] & x[None, :])
            Yps.append(~Yb[ii] & ok[ii, None] & x[None, :])
        T, SE2, N_, K_, comps = studentized_sum(Us, Vs, Yds, Yps)
        for k, val in zip(("T", "SE2", "N", "K"), (T, SE2, N_, K_)):
            out[k][sl] = val
        if c0 == 0:
            comps0 = comps
    with np.errstate(all="ignore"):
        out["z"] = np.where(out["SE2"] > 0, out["T"] / np.sqrt(out["SE2"]), np.nan)
    out["comps"] = comps0
    return out


def auroc_matrix(S: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """[score columns x label columns] AUROC of positives (Y) scoring higher, by pooled mid-ranks over
    the tumours with a finite score (ties half)."""
    from scipy.stats import rankdata
    S = np.asarray(S, float)
    fin = np.isfinite(S)
    R = np.where(fin, rankdata(np.where(fin, S, np.inf), axis=0), 0.0)
    Yf = Y.astype(np.float64)
    n1 = fin.T.astype(np.float64) @ Yf
    n0 = fin.sum(0)[:, None] - n1
    s1 = R.T @ Yf
    with np.errstate(all="ignore"):
        return (s1 - n1 * (n1 + 1) / 2) / (n1 * n0)


def mean_by_strata(D: np.ndarray, first: np.ndarray, ctrl: np.ndarray, lab: np.ndarray) -> np.ndarray:
    """Per row of D [columns x exons]: the stratum-weighted difference of class means over the exons with
    a finite value (weights = the first class's exons per stratum, strata lacking a class dropped)."""
    D = np.atleast_2d(np.asarray(D, float))
    strata = sorted(set(lab[(first | ctrl) & (lab != "")].tolist()))
    F = np.column_stack([(lab == s_) & first for s_ in strata]).astype(float)
    Cm = np.column_stack([(lab == s_) & ctrl for s_ in strata]).astype(float)
    fin = np.isfinite(D)
    D0 = np.where(fin, D, 0)
    nf, nc = fin.astype(float) @ F, fin.astype(float) @ Cm
    with np.errstate(all="ignore"):
        diff = (D0 @ F) / nf - (D0 @ Cm) / nc
    use = (nf > 0) & (nc > 0)
    w = np.where(use, nf, 0)
    with np.errstate(all="ignore"):
        return np.where(w.sum(1) > 0, (np.where(use, diff, 0) * w).sum(1) / w.sum(1), np.nan)


def class_perm(v: np.ndarray, first: np.ndarray, lab: np.ndarray, n_perm: int = N_CLASS_PERM,
               seed: int = SEED) -> tuple[float, np.ndarray]:
    """class_perm: the observed stratum-weighted difference of class means and its null, the class
    labels permuted among the exons of each stratum (weights = first-class exons per stratum).
    v - values of the exons in the comparison (finite); first - first-class flag; lab - stratum."""
    rng = np.random.default_rng(seed)
    strata = sorted(set(lab.tolist()))
    obs, null, wsum = 0.0, np.zeros(n_perm), 0
    for s_ in strata:
        ii = np.flatnonzero(lab == s_)
        f = int(first[ii].sum())
        if f == 0 or f == len(ii):
            continue
        x = v[ii]
        tot = x.sum()
        o = x[first[ii]].sum()
        obs += f * (o / f - (tot - o) / (len(ii) - f))
        pick = rng.random((n_perm, len(ii))).argsort(axis=1)[:, :f]
        sp = x[pick].sum(1)
        null += f * (sp / f - (tot - sp) / (len(ii) - f))
        wsum += f
    return (obs / wsum, null / wsum) if wsum else (np.nan, np.full(n_perm, np.nan))


def residualize(Z: np.ndarray | None, n: int) -> np.ndarray:
    """[n x n] residual-maker of an OLS on an intercept and Z."""
    X = np.ones((n, 1)) if Z is None else np.column_stack([np.ones(n), Z])
    return np.eye(n) - X @ np.linalg.pinv(X)


def partial_spearman(x, y, Z=None) -> float:
    """Correlation of the residuals of rank(x) and rank(y), each regressed on an intercept and Z (Z
    already in the units wanted, e.g. ranks or indicators)."""
    from scipy.stats import rankdata
    H = residualize(Z, len(x))
    rx, ry = H @ rankdata(x), H @ rankdata(y)
    den = np.sqrt((rx @ rx) * (ry @ ry))
    return float(rx @ ry / den) if den > 0 else np.nan


def within_perm_index(cells, n_perm: int, seed: int = SEED) -> np.ndarray:
    """[n_perm x n] permutation indices that move items only within their cell."""
    cells = np.asarray(cells)
    rng = np.random.default_rng(seed)
    out = np.tile(np.arange(len(cells)), (n_perm, 1))
    for c in sorted(set(cells.tolist())):
        ii = np.flatnonzero(cells == c)
        if len(ii) > 1:
            out[:, ii] = ii[rng.random((n_perm, len(ii))).argsort(axis=1)]
    return out


def partial_spearman_perm(x, y, Z, cells, n_perm: int = N_CLASS_PERM, seed: int = SEED) -> tuple[float, np.ndarray]:
    """Partial Spearman of x and y given Z, and its null with y permuted within cells (the ranks of a
    within-cell permutation of y are the same permutation of y's ranks)."""
    from scipy.stats import rankdata
    H = residualize(Z, len(x))
    rx, ry = H @ rankdata(x), rankdata(y)
    obs = partial_spearman(x, y, Z)
    pidx = within_perm_index(cells, n_perm, seed)
    R = H @ ry[pidx].T                                   # n x n_perm
    den = np.sqrt((rx @ rx) * (R * R).sum(0))
    with np.errstate(all="ignore"):
        return obs, (rx @ R) / den


def ols_hc3(X: np.ndarray, y: np.ndarray, j: int) -> tuple[float, float]:
    """OLS coefficient j and its HC3 t."""
    XtXi = np.linalg.pinv(X.T @ X)
    Pm = XtXi @ X.T
    b = Pm @ y
    e = y - X @ b
    h = np.einsum("ij,ji->i", X, Pm)
    var = (Pm[j] ** 2 * e ** 2 / (1 - h) ** 2).sum()
    return float(b[j]), float(b[j] / np.sqrt(var)) if var > 0 else np.nan


def freedman_lane(X: np.ndarray, y: np.ndarray, j: int, cells, n_perm: int = N_FL_PERM, seed: int = SEED,
                  chunk: int = 1000) -> tuple[float, float, np.ndarray]:
    """Freedman-Lane: residuals of the model without column j permuted within cells and added to its
    fit; the full model refitted; the HC3 t of column j. Returns (coefficient, t, null t)."""
    b, t = ols_hc3(X, y, j)
    Xr = np.delete(X, j, axis=1)
    fit_r = Xr @ (np.linalg.pinv(Xr) @ y)
    res_r = y - fit_r
    XtXi = np.linalg.pinv(X.T @ X)
    Pm = XtXi @ X.T
    h = np.einsum("ij,ji->i", X, Pm)
    c = Pm[j] ** 2 / (1 - h) ** 2
    pidx = within_perm_index(cells, n_perm, seed)
    null = np.empty(n_perm)
    for a in range(0, n_perm, chunk):
        Ys = fit_r[:, None] + res_r[pidx[a:a + chunk]].T
        B = Pm @ Ys
        E = Ys - X @ B
        with np.errstate(all="ignore"):
            null[a:a + chunk] = B[j] / np.sqrt(c @ E ** 2)
    return b, t, null


def quasi_binomial(y, n, off, X, max_iter: int = 50, tol: float = 1e-8, beta0=None):
    """Binomial GLM (logit) of y successes out of n with an offset, by IRLS; quasi-binomial scale =
    Pearson X2 / (rows - columns). Returns (beta, scaled covariance, scale, converged)."""
    from scipy.special import expit
    y, n, off = (np.asarray(v, float) for v in (y, n, off))
    p = X.shape[1]
    if beta0 is None:
        mu = (y + 0.5) / (n + 1)
        eta = np.log(mu / (1 - mu))
    else:
        eta = X @ beta0 + off
        mu = expit(eta)
    dev_old, conv, beta = np.inf, False, beta0
    for _ in range(max_iter):
        w = n * mu * (1 - mu)
        z = (eta - off) + (y - n * mu) / np.maximum(w, 1e-300)
        XtW = X.T * w
        beta = np.linalg.solve(XtW @ X + 1e-12 * np.eye(p), XtW @ z)
        eta = X @ beta + off
        mu = np.clip(expit(eta), 1e-12, 1 - 1e-12)
        with np.errstate(all="ignore"):
            dev = 2 * np.sum(np.where(y > 0, y * np.log(y / (n * mu)), 0)
                             + np.where(n - y > 0, (n - y) * np.log((n - y) / (n * (1 - mu))), 0))
        if abs(dev - dev_old) < tol * (abs(dev) + 0.1):
            conv = True
            break
        dev_old = dev
    w = n * mu * (1 - mu)
    x2 = float(np.sum((y - n * mu) ** 2 / np.maximum(w, 1e-300)))
    scale = x2 / max(len(y) - p, 1)
    cov = scale * np.linalg.pinv((X.T * w) @ X)
    return beta, cov, scale, conv


def coarsened_draws(lab: np.ndarray, first: np.ndarray, ctrl: np.ndarray, per: int, n_draws: int, rng) -> np.ndarray:
    """[n_draws x first exons x per] control exons drawn without replacement within each coarsened
    stratum (strata in sorted order, first-class exons in index order); -1 where a stratum is short."""
    fi = np.flatnonzero(first)
    out = np.full((n_draws, len(fi), per), -1, int)
    pos = {e: k for k, e in enumerate(fi)}
    for s_ in sorted(set(lab[first].tolist())):
        f = np.flatnonzero(first & (lab == s_))
        c = np.flatnonzero(ctrl & (lab == s_))
        need = len(f) * per
        if len(c) < need:
            continue
        pick = c[rng.random((n_draws, len(c))).argsort(axis=1)[:, :need]].reshape(n_draws, len(f), per)
        out[:, [pos[e] for e in f], :] = pick
    return out


def draw_weights(draws: np.ndarray, cols: np.ndarray):
    """class_stat weights for draws: column d of the statistic averages delta over the exons drawn in
    draw d (draws [draws x k] exon indices, -1 = none; cols the class_stat exon columns)."""
    pos = {e: k for k, e in enumerate(cols)}
    W = np.zeros((len(draws), len(cols)), bool)
    for d in range(len(draws)):
        W[d, [pos[e] for e in draws[d] if e >= 0]] = True

    def fn(valid, sl):
        x = valid & W[sl]
        return x / np.maximum(x.sum(1, keepdims=True), 1)
    return fn


def rp6c_draws(mm: np.ndarray, project: np.ndarray, prof: dict, rng, draws: int = RP6C_DRAWS,
               group: int = RP6C_GROUP) -> tuple[list, list]:
    """RP6c's held-out draws (D-70's RV1): 20 of the MMRd tumours, each with a random MMRp tumour of its
    cancer, no MMRp tumour twice within a draw."""
    ia, ib = [], []
    for _ in range(draws):
        a = rng.choice(mm, group, replace=False)
        b, used = [], set()
        for i in a:
            pool = [j for j in prof[project[i]] if j not in used]
            j = int(rng.choice(pool))
            used.add(j)
            b.append(j)
        ia.append(a)
        ib.append(np.array(b))
    return ia, ib


def instrument_verdict(draw_z: np.ndarray, null_z: np.ndarray) -> dict:
    """instrument_check.pass: the median draw statistic < 0 and >= 80% of draws below the 2.5th
    percentile of RPC2a's studentized permutation null."""
    crit = float(np.nanpercentile(null_z, 2.5))
    dz = draw_z[np.isfinite(draw_z)]
    share = float((dz < crit).mean()) if len(dz) else np.nan
    med = float(np.median(dz)) if len(dz) else np.nan
    return {"critical_value": crit, "median_draw": med, "share_below": share,
            "passes": bool(len(dz) and med < 0 and share >= INSTR_PASS)}


def rise_share(rise: np.ndarray, sets: np.ndarray) -> np.ndarray:
    """RP2c: per row of rise [columns x exons], the share of first exons that rise minus the share of
    their matched controls. sets - [sets x 6] exon indices, the first exon in column 0."""
    r = rise[:, sets]
    return r[:, :, 0].mean(1) - r[:, :, 1:].mean((1, 2))


def rise_class_perm(rise_obs: np.ndarray, sets: np.ndarray, n_perm: int = N_CLASS_PERM, seed: int = SEED) -> np.ndarray:
    """RP2c's p_exon null: the 'first' label moved to a random member of each matched set."""
    rng = np.random.default_rng(seed)
    r = rise_obs[sets].astype(float)                      # sets x 6
    k = rng.integers(0, sets.shape[1], (n_perm, len(sets)))
    pick = r[np.arange(len(sets))[None, :], k]
    tot = r.sum(1)
    return pick.mean(1) - ((tot[None, :] - pick) / (sets.shape[1] - 1)).mean(1)


def pair_stat(arm: np.ndarray, Y: np.ndarray, groups) -> np.ndarray:
    """RP7: per label column, the mean over MMRd pairs of (pair arm minus the mean pair arm of the MMRp
    pairs of the same cancer); pairs without an arm leave."""
    ok = np.isfinite(arm)
    a, Y, groups = arm[ok], Y[ok].astype(bool), np.asarray(groups)[ok]
    num, den = np.zeros(Y.shape[1]), np.zeros(Y.shape[1])
    for g in sorted(set(groups.tolist())):
        ii = np.flatnonzero(groups == g)
        y = Y[ii]
        nd, npp = y.sum(0), (~y).sum(0)
        sd, sp = (a[ii, None] * y).sum(0), (a[ii, None] * ~y).sum(0)
        use = (nd > 0) & (npp > 0)
        with np.errstate(all="ignore"):
            num += np.where(use, sd - nd * sp / np.maximum(npp, 1), 0)
        den += np.where(use, nd, 0)
    with np.errstate(all="ignore"):
        return np.where(den > 0, num / den, np.nan)


def contrast_draw_stat(A: np.ndarray, ia: list, ib: list, min_meas: int = RP6C_MIN_MEAS) -> float:
    """RP6c (D-70's RV1 in JC): the mean over draws of the mean over exons (>= min_meas measured
    tumours per group) of the difference in mean adj, group a minus group b."""
    out = []
    for a, b in zip(ia, ib):
        X, Z = A[a], A[b]
        fa, fb = np.isfinite(X).sum(0), np.isfinite(Z).sum(0)
        ok = (fa >= min_meas) & (fb >= min_meas)
        if ok.any():
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out.append(float(np.mean(np.nanmean(X[:, ok], 0) - np.nanmean(Z[:, ok], 0))))
    return float(np.mean(out)) if out else np.nan


def null_contrasts(tasks: dict, mmr_samples: set, d71: set, min_class: int = RP6C_MIN_CLASS) -> list:
    """RP6c's null: every pair of classes with >= min_class tumours in each subtype cohort, D-71 tumours
    only, MSI- and POLE-labelled and MMRd tumours removed. tasks - {cohort: table with file_id,
    sample15, label}. Returns [(cohort, a, b, file_ids a, file_ids b)] in sorted order."""
    out = []
    for cohort in sorted(tasks):
        u = tasks[cohort]
        u = u[~(u.label.isin(["MSI", "POLE"]) | u.sample15.isin(mmr_samples)) & u.file_id.isin(d71)]
        c = u.label.value_counts()
        labs = sorted(c[c >= min_class].index)
        for i, a in enumerate(labs):
            for b in labs[i + 1:]:
                out.append((cohort, a, b, u.file_id[u.label.eq(a)].tolist(), u.file_id[u.label.eq(b)].tolist()))
    return out


def msi_class_index(job_done: dict, task: pd.DataFrame) -> int:
    """The MSI class of a subtype run: class 1 in its job_done.json, and y == 1 is labelled MSI."""
    assert job_done["classes"][1] == "MSI", job_done["classes"]
    assert set(task.label[task.y.eq(1)]) == {"MSI"}, set(task.label[task.y.eq(1)])
    return 1


def shap_directions(X: np.ndarray, S: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Per input: Spearman of its standardised value with its SHAP over the given rows (NaN with fewer
    than 3 rows or no spread) - explain._direction per fold."""
    from scipy.stats import rankdata
    X, S = np.asarray(X, float)[rows], np.asarray(S, float)[rows]
    out = np.full(X.shape[1], np.nan)
    if len(rows) < 3:
        return out
    rx, rs = rankdata(X, axis=0), rankdata(S, axis=0)
    rx, rs = rx - rx.mean(0), rs - rs.mean(0)
    den = np.sqrt((rx ** 2).sum(0) * (rs ** 2).sum(0))
    ok = (np.ptp(X, axis=0) > 0) & (np.ptp(S, axis=0) > 0) & (den > 0)
    out[ok] = (rx * rs).sum(0)[ok] / den[ok]
    return out


def m4_guard(auroc: float, real: bool) -> None:
    """Stop when the refitted M4 does not reproduce step 54's held-out AUROC (0.900 within 0.001); the
    dry run's shuffled labels cannot reproduce it, so there it is reported, not enforced."""
    if real and not abs(auroc - M4_AUROC) <= M4_TOL + 1e-12:
        raise SystemExit(f"M4 refit: held-out AUROC {auroc:.4f}, not {M4_AUROC} +/- {M4_TOL}; stopping")


# ---------------------------------------------------------------- verdicts (code, never by hand)

def test_verdict(est: float, q: float, sign: int, testable: bool = True) -> str:
    """supported = the predicted sign and q < 0.05; against = the opposite sign and q < 0.05."""
    if not testable:
        return "not testable"
    if not (np.isfinite(est) and np.isfinite(q)):
        return "not resolved"
    if q < ALPHA and np.sign(est) == sign:
        return "supported"
    if q < ALPHA and np.sign(est) == -sign:
        return "against"
    return "not resolved"


def claim_verdict(tests: pd.DataFrame, instrument: dict, sesoi: float = SESOI) -> dict:
    """The claim-level verdict of msi_rising_proof.verdict from the tests table (columns id, verdict,
    estimate, ci_high, family) and the instrument check; precedence uninterpretable > refuted >
    established > partial."""
    t = tests.set_index("id")

    def v(i):
        return t.verdict.get(i, "not testable")

    def est(i):
        return float(t.estimate.get(i, np.nan))

    rp2a, ub = v("RP2a"), float(t.ci_high.get("RP2a", np.nan))
    e3c, e3e = est("RP3c"), est("RP3e")
    rp3e_testable = v("RP3e") != "not testable"
    branch_id = "RP3e" if rp3e_testable else "RP3c"
    clause_i = rp2a == "not resolved" and (est("RP2a") <= 0 or ub < sesoi)
    refuted = {"(0) RP2a against": rp2a == "against",
               "(i) RP2a not resolved, instrument check passes, no rise of interest": clause_i and bool(instrument.get("passes")),
               "(ii) RP2a supported and RP3a's estimate <= 0": rp2a == "supported" and est("RP3a") <= 0,
               "(iii) RP2a supported and the branch-point estimates <= 0":
                   rp2a == "supported" and ((e3e <= 0 and e3c <= 0) if rp3e_testable else e3c <= 0)}
    rp_ids = ["RP2a", "RP3a", branch_id, "RP4a", "RP4b", "RP4c", "RP4d", "RP4e", "RP5a", "RP5b", "RP6a", "RP6c"]
    unmet = []
    if v("RPC2a") != "supported":
        unmet.append("RPC2a supported")
    unmet += [f"{i} supported ({v(i)})" for i in rp_ids if v(i) != "supported"]
    if not ub >= sesoi:
        unmet.append(f"RP2a's 95% upper bound reaches the SESOI ({ub:.4g} < {sesoi})")
    if not (est("RP7") > 0 and v("RP7") != "against"):
        unmet.append(f"RP7's estimate > 0 and not against ({est('RP7'):.4g}, {v('RP7')})")
    reasons = []
    if v("RPC2a") != "supported":
        verdict = "uninterpretable"
        reasons = [f"RPC2a {v('RPC2a')}: junction reads held out cannot see the known falling arm"]
    elif any(refuted.values()):
        verdict = "refuted"
        reasons = [k for k, x in refuted.items() if x]
    elif not unmet:
        verdict = "established"
        reasons = ["all conditions met"]
    else:
        verdict = "partial"
        reasons = list(unmet)
        if rp2a == "supported" and ub < sesoi:
            reasons.append("rise below the smallest effect of interest")
        if clause_i and not instrument.get("passes"):
            reasons.append("instrument cannot see a near-silent falling class")
    return {"verdict": verdict, "reasons": reasons, "branch_point_test": branch_id,
            "refuted_clauses": refuted, "instrument_passes": bool(instrument.get("passes")),
            "rp2a": {"verdict": rp2a, "estimate": est("RP2a"), "ci_high": ub, "sesoi": sesoi},
            "rp7": {"verdict": v("RP7"), "estimate": est("RP7")}}


# ---------------------------------------------------------------- compute: shared state

N_WORKERS = 8                  # the machine is shared: at most 8 processes

TEST_META = {   # id: (predicted sign, short claim); families are read from the frozen block
    "RP1": (1, "R_up dev-arm AUROC minus the median of 1,000 matched B sets"),
    "RPC1": (-1, "R_down dev-arm AUROC minus the median of 1,000 matched B sets"),
    "RPD1a": (1, "share of R_up SHAP directions > 0, minus 0.5"),
    "RPD1b": (-1, "share of R_down SHAP directions > 0, minus 0.5"),
    "RP2a": (1, "mean delta over eligible R_up exons (JC, PSI)"),
    "RP2b": (1, "R_up arm_contrast"),
    "RP2c": (1, "share of R_up exons rising minus the share of matched B-sample exons"),
    "RPC2a": (-1, "mean delta over eligible R_down exons"),
    "RPC2b": (-1, "R_down arm_contrast"),
    "RPD2a": (1, "training: mean delta over R_up exons"),
    "RPD2b": (-1, "training: mean delta over R_down exons"),
    "RP3a": (1, "R_up minus B_far, stratum-weighted mean delta"),
    "RP3b": (1, "R_up minus near_bp"),
    "RP3c": (1, "partial Spearman of delta with g given L and o (heuristic branch point)"),
    "RP3d": (1, "experimental branch point d_exp >= 36 minus < 36"),
    "RP3e": (1, "partial Spearman of delta with g_exp given L and o"),
    "RPD3a": (1, "training: R_up minus B_far"),
    "RPD3b": (1, "training: partial Spearman of delta with g_exp given L and o"),
    "RP4a": (1, "R_up minus polyA (>= 15 nt)"),
    "RP4b": (1, "R_up minus ppt_far"),
    "RP4c": (1, "R_up minus deep"),
    "RP4d": (1, "R_up minus donor"),
    "RP4e": (1, "R_up minus alu"),
    "RP5a": (1, "OLS coefficient of R_up on delta, exon covariates (vs B sample)"),
    "RP5b": (1, "OLS coefficient of MMRd on the R_up arm, tumour covariates"),
    "RPC5a": (-1, "OLS coefficient of R_down on delta, exon covariates (vs B sample)"),
    "RPC5b": (-1, "OLS coefficient of MMRd on the R_down arm, tumour covariates"),
    "RP6a": (1, "AUROC of the R_up arm z, MMRd vs hypermutated MMRp, minus 0.5"),
    "RP6b": (1, "R_up arm_contrast, MANTIS MSI-H vs MSS"),
    "RP6c": (1, "MMRd contrast s (JC adj) against 57 subtype contrasts"),
    "RPC6": (-1, "AUROC of the R_down arm z, MMRd vs hypermutated MMRp, minus 0.5"),
    "RP7": (1, "tumour-normal pair arm, MMRd minus MMRp pairs"),
    "RPC7": (-1, "tumour-normal pair arm (R_down), MMRd minus MMRp pairs"),
    "RP8a": (1, "partial Spearman of the R_up arm with MANTIS given cancer, MMRd"),
    "RP8b": (1, "partial Spearman of the R_up arm with log(1 + deletions) given cancer, MMRd"),
    "RP8c": (1, "Spearman of delta with run length over R_up"),
    "RPC8": (-1, "partial Spearman of the R_down arm with MANTIS given cancer, MMRd"),
    "RPM1": (1, "R_up minus shield"),
    "RPM2": (1, "mean delta over shield exons"),
    "RPM3": (1, "runs >= 20 in an antisense Alu minus not"),
}


class XCtx:
    """Exon-level read-out of one spec and side: the engine's columns (exons eligible at the lowest
    threshold), their table, and delta / studentized delta / N over the label columns; elig, band,
    terc, base and reads at the threshold in use."""

    def __init__(self, spec, side, cols, tab, D, T, N, projects, elig=None, band=None, terc=None, base=None, reads=None):
        self.spec, self.side, self.cols, self.tab, self.D, self.T, self.N = spec, side, cols, tab, D, T, N
        self.projects, self.elig, self.band, self.terc, self.base, self.reads = projects, elig, band, terc, base, reads

    def at(self, side_tab: pd.DataFrame) -> "XCtx":
        """The same columns at another eligibility threshold (its band, tertile and baseline)."""
        s = side_tab.iloc[self.cols]
        return XCtx(self.spec, self.side, self.cols, self.tab, self.D, self.T, self.N, self.projects,
                    s.eligible.to_numpy(bool), s.band.to_numpy(int), s.tertile.to_numpy(int),
                    s.baseline.to_numpy(float), s.reads.to_numpy(float))

    def m(self, col: str) -> np.ndarray:
        return self.tab[col].to_numpy()


class State:
    """What one compute run shares."""

    def __init__(self, run: Run, fresh: bool = False):
        import yaml
        self.run, self.fresh, self.t0 = run, fresh, time.time()
        self.timing = {}
        self.pre = yaml.safe_load(open(Path(__file__).resolve().parents[1] / "config" / "followup.yaml"))["msi_rising_proof"]
        self.family = {t["id"]: t["family"] for t in self.pre["tests"]}
        assert set(self.family) == set(TEST_META), "TEST_META does not match the frozen block"
        self.lab = tumour_labels(run)
        self.pos = pd.Series(np.arange(len(self.lab)), index=self.lab.file_id)
        self.ex = exon_table()
        self.pairs = pairs_table()
        self.n_mmr = self.lab.groupby("project").mmr.sum()
        self.projects = sorted(self.lab.project.unique())
        self.held = [p for p in self.projects if p not in TRAIN]
        self.Y = label_perms(self.lab.mmr.to_numpy(), self.lab.project.to_numpy(), N_PERM, SEED)
        self.res, self.sides, self.ctx = {}, {}, {}
        self.rows, self.sens, self.desc = [], [], []
        log(f"labels ({run.labels}): {int(self.lab.mmr.sum())} MMRd tumours "
            f"({int(self.lab[self.lab.held_out].mmr.sum())} held out), {len(self.lab):,} tumours; "
            f"{len(self.ex):,} exon columns")

    def tick(self, what: str) -> None:
        self.timing[what] = round(time.time() - self.t0, 1)
        log(f"[{self.timing[what]:.0f} s] {what}")

    def workdir(self, spec: Spec) -> Path:
        return self.run.path(f"adj_{spec.name}/README", local=True).parent

    def adj(self, spec: Spec, proj: str, cols) -> np.ndarray:
        return np.asarray(np.load(self.workdir(spec) / f"{proj}.npy", mmap_mode="r")[:, cols], np.float64)

    def rows_of(self, spec: Spec, proj: str) -> np.ndarray:
        return self.pos.loc[self.res[spec.name][proj]["file_id"]].to_numpy()

    def Yp(self, spec: Spec, proj: str, Y=None) -> np.ndarray:
        return (self.Y if Y is None else Y)[self.rows_of(spec, proj)]

    def with_mmr(self, ps) -> list:
        return [p for p in ps if self.n_mmr.get(p, 0) > 0]

    # ---- builds
    def build(self, spec: Spec, projects: list) -> dict:
        import pickle
        from joblib import Parallel, delayed
        f = self.run.path(f"build_{spec.name}.pkl", local=True)
        if f.exists() and not self.fresh:
            self.res[spec.name] = pickle.loads(f.read_bytes())
            return self.res[spec.name]
        ex = self.ex
        flags = {"b_sample": ex.b_sample.to_numpy(bool), "R_up": ex.cls.eq("R_up").to_numpy(),
                 "R_down": ex.cls.eq("R_down").to_numpy()}
        wd = self.workdir(spec)
        k_of = {p: k for k, p in enumerate(self.projects)}
        prim = spec.name == "primary"
        jobs = sorted(projects, key=lambda p: -int((self.lab.project == p).sum()))
        out = Parallel(n_jobs=N_WORKERS)(delayed(build_cancer)(
            p, spec, flags, str(wd), k_of[p], prim and p not in TRAIN,
            self.pairs[self.pairs.project.eq(p)] if prim else None) for p in jobs)
        self.res[spec.name] = {r["project"]: r for r in out}
        f.write_bytes(pickle.dumps(self.res[spec.name]))
        return self.res[spec.name]

    def side(self, spec: Spec) -> dict:
        if spec.name in self.sides:
            return self.sides[spec.name]
        f = self.run.path(f"sides_{spec.name}.parquet", local=True)
        if f.exists() and not self.fresh:
            d = pd.read_parquet(f)
            out = {(s_, int(t_)): g.drop(columns=["side", "thr"]).reset_index(drop=True) for (s_, t_), g in d.groupby(["side", "thr"])}
        else:
            src = spec if spec.units == "psi" else SPECS["primary"]          # percentile units: same universe
            res = self.res[src.name]
            n_stored = len(pd.read_parquet(PROOF / "exons.parquet", columns=["feature_id"]))
            out = side_tables(src, res, self.n_mmr, thresholds=(10, 20, 30), n_stored=n_stored)
            pd.concat([v.assign(side=k[0], thr=k[1]) for k, v in out.items()]).to_parquet(f, index=False)
        self.sides[spec.name] = out
        return out

    # ---- exon engine
    def xctx(self, spec: Spec, side: str) -> XCtx:
        key = (spec.name, side)
        if key in self.ctx:
            return self.ctx[key]
        sides = self.side(spec)
        base = sides[(side, 10 if side == "held_out" else 1)]
        cols = np.flatnonzero(base.eligible.to_numpy())
        ps = self.with_mmr(self.held if side == "held_out" else list(TRAIN))
        f = self.run.path(f"engine_{spec.name}_{side}.npz", local=True)
        if f.exists() and not self.fresh:
            z = np.load(f)
            D, T, N = z["D"], z["T"], z["N"]
        else:
            adjs = [np.load(self.workdir(spec) / f"{p}.npy", mmap_mode="r")[:, cols] for p in ps]
            D, T, N = exon_engine([np.asarray(a) for a in adjs], [self.Yp(spec, p) for p in ps])
            np.savez(f, D=D, T=T, N=N)
        tab = self.ex.iloc[cols].reset_index(drop=True)
        ctx = XCtx(spec, side, cols, tab, D, T, N, ps).at(sides[(side, 20 if side == "held_out" else 1)])
        self.ctx[key] = ctx
        return ctx

    # ---- result rows
    def add(self, tid: str, estimate: float, p_class=np.nan, p_label=np.nan, p=None, n: str = "",
            ci=(np.nan, np.nan), statistic=np.nan, value=np.nan, testable: bool = True, inference: str = "",
            note: str = "") -> None:
        if p is None:
            p = iut(p_class, p_label)
        self.rows.append({"id": tid, "family": self.family[tid], "claim": TEST_META[tid][1], "n": n,
                          "estimate": estimate, "ci_low": ci[0], "ci_high": ci[1], "statistic": statistic,
                          "value": value, "p_class": p_class, "p_label": p_label, "p": p,
                          "predicted_sign": TEST_META[tid][0], "testable": bool(testable),
                          "inference": inference, "note": note})
        log(f"{tid}: estimate {estimate:+.4g}, p {p:.4g}{'' if testable else ' (not testable)'} {n}")

    def add_sens(self, item: str, test: str, estimate: float, p=np.nan, n: str = "", note: str = "", **kw) -> None:
        self.sens.append({"item": item, "test": test, "estimate": estimate, "p": p, "n": n, "note": note, **kw})

    def add_desc(self, item: str, group: str, value: float, n="", note: str = "", **kw) -> None:
        self.desc.append({"item": item, "group": group, "value": value, "n": n, "note": note, **kw})


# ---------------------------------------------------------------- exon-level tests

def class_mean(S: State, ctx: XCtx, sel: np.ndarray, weight_fn=None, drop: tuple = ()) -> dict:
    """A studentized class statistic (RP2a's design) over the ctx columns in `sel`."""
    ps = [p for p in ctx.projects if p not in drop]
    cols = ctx.cols[sel]
    out = class_stat([S.adj(ctx.spec, p, cols) for p in ps], [S.Yp(ctx.spec, p) for p in ps], weight_fn)
    out["p"] = perm_p(out["z"][0], out["z"][1:])
    out["ci"] = interval(out["T"][0], out["SE2"][0], out["comps"])
    out["n_exons_obs"] = int(out["n_exons"][0])
    return out


def contrast(ctx: XCtx, first: np.ndarray, ctrl: np.ndarray, per_first: int = 0, lab_kept=None,
             min_n: int = MIN_CLASS, seed: int = SEED, D=None) -> dict:
    """A class contrast: stratum-weighted difference of mean delta within the coarsened strata, with
    class_perm (exons) and label_perm (tumours); p = max (intersection-union)."""
    D = ctx.D if D is None else D
    f, c = first & ctx.elig, ctrl & ctx.elig
    lab, kept = lab_kept if lab_kept is not None else coarsen(ctx.band, ctx.terc, f, c, per_first=per_first)
    use = kept & (f | c)                              # label-blind; each column drops exons without a delta
    stat = mean_by_strata(D[:, use], f[use], c[use], lab[use])
    fin = np.isfinite(D[0, use])
    obs, null = class_perm(D[0, use][fin].astype(float), f[use][fin], lab[use][fin], N_CLASS_PERM, seed)
    if np.isfinite(obs):
        assert abs(obs - stat[0]) < 1e-6 * max(1, abs(obs)), (obs, stat[0])
    obs_use = use.copy()
    obs_use[use] = fin
    good = {s_ for s_ in set(lab[obs_use & f].tolist()) if (obs_use & c & (lab == s_)).any()}
    obs_use &= np.isin(lab, list(good))
    nf, nc = int((obs_use & f).sum()), int((obs_use & c).sum())
    eff = effective_n(lab, f, c, obs_use)
    return {"estimate": float(stat[0]), "p_class": perm_p(stat[0], null), "p_label": perm_p(stat[0], stat[1:]),
            "n_first": nf, "n_ctrl": nc, "eff_ctrl": eff["control_effective_n"],
            "testable": nf >= min_n and nc >= min_n,
            "n": f"{nf} vs {nc} exons (effective controls {eff['control_effective_n']:.0f})"}


def geometry(ctx: XCtx, sel: np.ndarray, g: np.ndarray, L: np.ndarray, o: np.ndarray, lb: list, ob: list,
             min_n: int = 0) -> dict:
    """RP3c / RP3e: partial Spearman of delta with g given rank(L) and rank(o); p_exon by permuting g
    within L bin x o bin cells, p_label by recomputing delta; p = max."""
    from scipy.stats import rankdata
    rows = np.flatnonzero(sel & ctx.elig & np.isfinite(g) & np.isfinite(L) & np.isfinite(o))
    X = ctx.D[:, rows].astype(float)
    g_, L_, o_ = g[rows], L[rows], o[rows]
    f0 = np.isfinite(X[0])
    cells = np.array([f"{a}:{b}" for a, b in zip(np.digitize(L_, lb), np.digitize(o_, ob))])

    def ps(k):
        f = np.isfinite(X[k])
        return partial_spearman(X[k][f], g_[f], np.column_stack([rankdata(L_[f]), rankdata(o_[f])])) if f.sum() > 4 else np.nan
    Z0 = np.column_stack([rankdata(L_[f0]), rankdata(o_[f0])])
    obs, null = partial_spearman_perm(X[0][f0], g_[f0], Z0, cells[f0])
    lab_null = np.array([ps(k) for k in range(1, X.shape[0])])
    n = int(f0.sum())
    return {"estimate": obs, "p_class": perm_p(obs, null), "p_label": perm_p(obs, lab_null), "n_exons": n,
            "testable": n >= max(min_n, 5), "n": f"{n} exons"}


def spearman_within(ctx: XCtx, sel: np.ndarray, x: np.ndarray, cells: np.ndarray) -> dict:
    """RP8c-type: Spearman of delta with x over the selected exons; p_exon by permuting x within cells,
    p_label by recomputing delta; p = max."""
    rows = np.flatnonzero(sel & ctx.elig & np.isfinite(x))
    X = ctx.D[:, rows].astype(float)
    f0 = np.isfinite(X[0])
    obs, null = partial_spearman_perm(X[0][f0], x[rows][f0], None, cells[rows][f0])
    lab_null = []
    for k in range(1, X.shape[0]):
        f = np.isfinite(X[k])
        lab_null.append(partial_spearman(X[k][f], x[rows][f]) if f.sum() > 4 else np.nan)
    return {"estimate": obs, "p_class": perm_p(obs, null), "p_label": perm_p(obs, np.array(lab_null)),
            "n_exons": int(f0.sum()), "testable": int(f0.sum()) >= 5, "n": f"{int(f0.sum())} exons"}


def rp2c(ctx: XCtx) -> dict:
    """RP2c: the share of eligible R_up exons that rise (studentized delta > c) minus the share among
    5 matched B-sample exons each; c = the 97.5th percentile of the studentized delta of all eligible
    B-sample exons over the label permutations."""
    up, bs = ctx.m("cls") == "R_up", ctx.m("b_sample").astype(bool)
    f, c = up & ctx.elig, bs & ctx.elig
    cut = float(np.nanpercentile(ctx.T[1:, c], 97.5))
    lab, kept = coarsen(ctx.band, ctx.terc, f, c, per_first=5)
    dr = coarsened_draws(lab, f & kept, c & kept, 5, 1, np.random.default_rng(SEED))[0]
    fi = np.flatnonzero(f & kept)
    okr = (dr >= 0).all(1)
    sets = np.column_stack([fi[okr], dr[okr]])
    rise = np.nan_to_num(ctx.T, nan=-np.inf) > cut
    stat = rise_share(rise, sets)
    null = rise_class_perm(rise[0], sets)
    return {"estimate": float(stat[0]), "p_class": perm_p(stat[0], null), "p_label": perm_p(stat[0], stat[1:]),
            "cut": cut, "n_sets": len(sets), "share_up": float(rise[0, sets[:, 0]].mean()),
            "share_ctrl": float(rise[0, sets[:, 1:]].mean()), "testable": len(sets) >= MIN_CLASS,
            "n": f"{len(sets)} R_up exons x 5 matched B-sample exons; c = {cut:.3f}"}


def instrument_check(S: State, ctx: XCtx, rpc2a_null_z: np.ndarray) -> dict:
    """The instrument check: 1,000 draws of one eligible R_down exon per R_up exon from the coarsened
    stratum of R_up's mirrored baseline band (band of 1 - baseline) x coverage tertile, without
    replacement within a draw; per draw RP2a's studentized statistic with the observed labels."""
    up, dn = (ctx.m("cls") == "R_up") & ctx.elig, (ctx.m("cls") == "R_down") & ctx.elig
    band = ctx.band.copy()
    band[up] = band_of(1 - ctx.base[up])
    lab, kept = coarsen(band, ctx.terc, up, dn, per_first=1)
    draws = coarsened_draws(lab, up & kept, dn & kept, 1, N_DRAWS, np.random.default_rng(SEED))[:, :, 0]
    dcols = np.flatnonzero(dn)
    wfn = draw_weights(draws, dcols)
    ps = ctx.projects
    Ys = [np.repeat(S.Yp(ctx.spec, p)[:, :1], N_DRAWS, axis=1) for p in ps]
    out = class_stat([S.adj(ctx.spec, p, ctx.cols[dcols]) for p in ps], Ys, wfn)
    v = instrument_verdict(out["z"], rpc2a_null_z)
    hi = ctx.base >= BAND_EDGES[-2]                       # baseline in the top band [0.95, 1]
    v.update({"draw_z": out["z"], "draw_estimate": out["T"], "n_up_drawn": int((draws[0] >= 0).sum()),
              "r_down_pool": int((dn & kept).sum()), "r_up_mirrored_top_band": int((up & (band == len(BAND_EDGES) - 2)).sum()),
              "drawn_r_down_baseline_ge_095_mean": float(np.mean([hi[d[d >= 0]].sum() for d in draws]))})
    return v


def gene_change(Gs: list, Ys: list, Us: list, gi: np.ndarray, block: int = 200) -> np.ndarray:
    """delta's definition applied to expression, for every label column: per exon, sum over its universe
    cancers of [sum over MMRd tumours of (expression - the cancer's MMRp mean)] / their MMRd count.
    Gs - per cancer [tumours x genes]; Ys - per cancer [tumours x columns] labels; Us - per cancer
    [exons] universe membership (and gene measured); gi - each exon's gene column. NaN where no MMRd."""
    P, E = Ys[0].shape[1], len(gi)
    num, den = np.zeros((P, E)), np.zeros((P, E))
    for G, Y, U in zip(Gs, Ys, Us):
        nd = int(Y[:, 0].sum())
        if nd == 0 or not U.any() or nd == len(G):
            continue
        idx, tot = mmrd_index(Y), G.sum(0)
        for p0 in range(0, P, block):
            sD = G[idx[p0:p0 + block]].sum(1)
            num[p0:p0 + block] += (sD - nd * (tot - sD) / (len(G) - nd))[:, gi] * U
        den += nd * U
    with np.errstate(all="ignore"):
        return np.where(den > 0, num / den, np.nan)


def host_change(S: State, ctx: XCtx, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """RP5a's host-gene change (D-76): gene_change on log2(CPM + 1) of the exons' host genes over the
    held-out cancers with MMRd tumours, recomputed in every label column; 0 with an indicator where the
    gene is not measured."""
    from ctsm.matrix_io import load_view
    cfg = load_config()
    genes = ctx.tab.gene_id.to_numpy()[rows]
    feats = pd.read_parquet(Path(cfg.matrix_dir) / "gex" / "features.parquet")
    have = pd.Series(np.arange(len(feats)), index=feats.gene_id)
    measured = pd.Series(genes).isin(have.index).to_numpy()
    ug = sorted(set(genes[measured].tolist()))
    gi = np.where(measured, pd.Series(np.arange(len(ug)), index=ug).reindex(genes).fillna(0).to_numpy(int), 0)
    res = S.res["primary"]
    Gs, Ys, Us = [], [], []
    for p in ctx.projects:
        Gs.append(np.asarray(load_view(cfg.matrix_dir, "gex", res[p]["file_id"])[0].get()[:, have.loc[ug].to_numpy()], np.float64))
        Ys.append(S.Yp(ctx.spec, p))
        Us.append(res[p]["U"][ctx.cols[rows]] & measured)
    hc = gene_change(Gs, Ys, Us, gi)
    return np.where(measured, np.nan_to_num(hc), 0.0), ~measured


def rp5a(S: State, ctx: XCtx, cls: str, extra: tuple = (), hc=None) -> dict:
    """RP5a / RPC5a: OLS delta ~ class + C(band of baseline) + log exon length + log median JC reads +
    host-gene change (+ not-measured indicator) + poison + core_3ss + exon GC + exon Alu (+ extra);
    the class coefficient's HC3 t; p = max(Freedman-Lane, exon residuals within band x tertile;
    label_perm, delta and the host-gene change recomputed)."""
    t = ctx.tab
    first = (t.cls.eq(cls)).to_numpy()
    rows = np.flatnonzero(ctx.elig & (first | t.b_sample.to_numpy(bool)))
    def num(col):
        v = t[col]
        return (v.map({True: 1.0, False: 0.0}).astype(float) if v.dtype == object or v.dtype == bool
                else v.astype(float)).to_numpy()[rows]
    cov = pd.DataFrame({"cls": first[rows].astype(float), "log_len": np.log(num("exon_len")),
                        "log_reads": np.log(ctx.reads[rows] + 1), "poison": t.consequence.eq("ptc_nmd").to_numpy(float)[rows],
                        "core_3ss": num("pwm_core"), "exon_gc": num("exon_gc"),
                        "exon_alu": t.exon_alu.eq(True).to_numpy(float)[rows]})
    for e in extra:
        cov[e] = num(e)
    if hc is None:
        hc = host_change(S, ctx, rows)
    HC, miss = hc
    cov["host_missing"] = miss.astype(float)
    bands = pd.get_dummies(pd.Series(ctx.band[rows]), prefix="band", drop_first=True).astype(float)
    ok = cov.notna().all(1).to_numpy()
    cells = np.array([f"{a}:{b}" for a, b in zip(ctx.band[rows], ctx.terc[rows])])

    def design(hcv, m):
        Xd = pd.concat([cov.reset_index(drop=True), bands.reset_index(drop=True)], axis=1).assign(host_change=hcv)
        Xd = Xd.loc[m]
        Xd = Xd.loc[:, Xd.std() > 0]
        return np.column_stack([np.ones(len(Xd)), Xd.to_numpy(float)]), list(Xd.columns)
    y0 = ctx.D[0, rows].astype(float)
    m0 = ok & np.isfinite(y0)
    X0, names = design(HC[0], m0)
    j = 1 + names.index("cls")
    b, tval, null = freedman_lane(X0, y0[m0], j, cells[m0], N_FL_PERM, SEED)
    lab_null = []
    for k in range(1, ctx.D.shape[0]):
        yk = ctx.D[k, rows].astype(float)
        mk = ok & np.isfinite(yk)
        Xk, nk = design(HC[k], mk)
        lab_null.append(ols_hc3(Xk, yk[mk], 1 + nk.index("cls"))[1] if "cls" in nk else np.nan)
    n1 = int((first[rows] & m0).sum())
    return {"estimate": b, "statistic": tval, "p_class": perm_p(tval, null), "p_label": perm_p(tval, np.array(lab_null)),
            "testable": n1 >= MIN_CLASS, "n": f"{n1} {cls} + {int(m0.sum()) - n1} B-sample exons",
            "hc": hc, "names": names}


# ---------------------------------------------------------------- tumour-level tests

def tumour_values(S: State) -> pd.DataFrame:
    """Per D-71 tumour (label order): the JC arms, the global axis and its percentile (primary read-out)."""
    out = pd.DataFrame(index=np.arange(len(S.lab)), columns=["arm_up", "arm_down", "ga", "ga_pct"], dtype=float)
    for p, r in S.res["primary"].items():
        i = S.rows_of(SPECS["primary"], p)
        for k in out.columns:
            out.loc[i, k] = r[k]
    return out


def rp1(S: State, arm: str) -> dict:
    """RP1 / RPC1: pooled held-out AUROC of the class dev-arm percentile minus the median AUROC of 1,000
    matched B-sample sets; p_exon over the sets, p_label with every AUROC recomputed; p = max."""
    held = S.lab.held_out.to_numpy()
    Sc = np.full((len(S.lab), N_SETS + 1), np.nan, np.float32)
    for p in S.held:
        r = S.res["primary"][p]
        if f"rp1_{arm}" in r:
            Sc[S.rows_of(SPECS["primary"], p)] = r[f"rp1_{arm}"]
    A = auroc_matrix(Sc[held], S.Y[held])
    med = np.median(A[1:], axis=0)
    stat = A[0] - med
    p_exon = float((1 + np.sum(np.abs(A[1:, 0] - med[0]) >= abs(stat[0]) - 1e-12)) / (N_SETS + 1))
    fin = np.isfinite(Sc[held, 0])
    nd = int((fin & S.lab.mmr.to_numpy()[held]).sum())
    rules = {}
    for p in S.held:
        for k, v in S.res["primary"][p].get(f"rp1_{arm}_rules", {}).items():
            rules[k] = rules.get(k, 0) + v
    for k, v in sorted(rules.items()):
        S.add_desc(f"RP1 matched sets ({arm}): exon x cancer slots by draw rule", k, float(v), n=int(sum(rules.values())))
    return {"estimate": float(stat[0]), "value": float(A[0, 0]), "p_class": p_exon, "p_label": perm_p(stat[0], stat[1:]),
            "testable": nd >= MIN_TUMOUR, "n": f"{int(fin.sum())} tumours ({nd} MMRd); null median {med[0]:.4f}"}


def arm_test(S: State, arm: np.ndarray, sel: np.ndarray, Y=None) -> dict:
    """RP2b-type: studentized arm_contrast over the selected tumours, label_perm."""
    Y = S.Y if Y is None else Y
    out = arm_stat(arm[sel], Y[sel], S.lab.project.to_numpy()[sel])
    nd = int(out["N"][0])
    lo, hi, df = interval(out["T"][0], out["SE2"][0], out["comps"])
    return {"estimate": float(out["T"][0]), "statistic": float(out["z"][0]), "p_label": perm_p(out["z"][0], out["z"][1:]),
            "ci": (lo, hi), "testable": nd >= MIN_TUMOUR, "n": f"{nd} MMRd tumours with an arm", "out": out}


def covariates(S: State) -> pd.DataFrame:
    c = pd.read_parquet(LOCAL / "pancancer_covariates.parquet")[["file_id", "stromal_p", "immune_p", "ri_mean_p", "purity_p"]]
    return S.lab[["file_id"]].merge(c, on="file_id", how="left")


def rp5b(S: State, arm: np.ndarray, purity: bool = True, extra: np.ndarray | None = None) -> dict:
    """RP5b / RPC5b: OLS arm ~ MMRd + C(project) + global_axis + stromal_p + immune_p + ri_mean_p
    (+ purity_p); the MMRd coefficient's HC3 t; Freedman-Lane, tumour residuals within cancer."""
    tv, cv = S.tv, covariates(S)
    cols = ["stromal_p", "immune_p", "ri_mean_p"] + (["purity_p"] if purity else [])
    m = S.lab.held_out.to_numpy() & np.isfinite(arm) & np.isfinite(tv.ga_pct.to_numpy(float)) & cv[cols].notna().all(1).to_numpy()
    if extra is not None:
        m &= np.isfinite(extra)
    proj = S.lab.project.to_numpy()[m]
    dums = pd.get_dummies(pd.Series(proj), drop_first=True).astype(float).to_numpy()
    parts = [np.ones(m.sum()), dums, tv.ga_pct.to_numpy(float)[m], cv[cols].to_numpy(float)[m]]
    if extra is not None:
        parts.append(extra[m])
    X = np.column_stack(parts + [S.lab.mmr.to_numpy(float)[m]])
    b, t, null = freedman_lane(X, arm[m], X.shape[1] - 1, proj, N_FL_PERM, SEED)
    nd = int(S.lab.mmr.to_numpy()[m].sum())
    return {"estimate": b, "statistic": t, "p_class": perm_p(t, null), "testable": nd >= MIN_TUMOUR,
            "n": f"{int(m.sum())} tumours ({nd} MMRd)"}


def rp6a(S: State, arm: np.ndarray, within: bool = False) -> dict:
    """RP6a / RPC6: each held-out tumour's arm as a z within its cancer against the cancer's
    non-hypermutated MMRp arms; AUROC of MMRd against hypermutated MMRp; Mann-Whitney."""
    from scipy.stats import mannwhitneyu
    lab = S.lab
    z = np.full(len(lab), np.nan)
    for p in S.held:
        i = np.flatnonzero(lab.project.eq(p).to_numpy())
        ref = i[(~lab.mmr.to_numpy()[i]) & (~lab.hyper.to_numpy()[i]) & np.isfinite(arm[i])]
        if len(ref) > 2:
            z[i] = (arm[i] - arm[ref].mean()) / arm[ref].std(ddof=1)
    pos = lab.held_out.to_numpy() & lab.mmr.to_numpy() & np.isfinite(z)
    neg = lab.held_out.to_numpy() & lab.hyper.to_numpy() & ~lab.mmr.to_numpy() & np.isfinite(z)
    if within:
        both = set(lab.project[pos]) & set(lab.project[neg])
        keep = lab.project.isin(both).to_numpy()
        pos, neg = pos & keep, neg & keep
    if not pos.any() or not neg.any():
        return {"estimate": np.nan, "p_class": np.nan, "testable": False, "n": "no tumours"}
    u = mannwhitneyu(z[pos], z[neg], alternative="two-sided")
    au = float(u.statistic / (pos.sum() * neg.sum()))
    return {"estimate": au - 0.5, "value": au, "p_class": float(u.pvalue), "testable": pos.sum() >= MIN_TUMOUR,
            "n": f"{int(pos.sum())} MMRd vs {int(neg.sum())} hypermutated MMRp tumours"}


def msi_call_test(S: State, arm: np.ndarray, call: str) -> dict:
    """RP6b: as RP2b with a DNA MSI call as the label (MSI-H against MSS; between / no score left out),
    MC3 status ignored; label_perm with this label."""
    lab = S.lab
    inc = lab.held_out.to_numpy() & lab[call].isin(["MSI-H", "MSS"]).to_numpy()
    y = lab[call].eq("MSI-H").to_numpy()
    Y = np.zeros((len(lab), N_PERM + 1), bool)
    Y[inc] = label_perms(y[inc], lab.project.to_numpy()[inc], N_PERM, SEED)
    r = arm_test(S, arm, inc, Y)
    r["n"] = f"{int(r['out']['N'][0])} {call} MSI-H tumours with an arm against {int((inc & ~y & np.isfinite(arm)).sum())} MSS"
    return r


def dose(S: State, arm: np.ndarray, y: np.ndarray, sel: np.ndarray, purity: np.ndarray | None = None) -> dict:
    """RP8a/b, RPC8: partial Spearman of the arm with y given cancer (ranks residualised on cancer
    indicators, and on rank(purity) if given) over the selected tumours; y permuted within cancer."""
    from scipy.stats import rankdata
    m = sel & np.isfinite(arm) & np.isfinite(y)
    if purity is not None:
        m &= np.isfinite(purity)
    proj = S.lab.project.to_numpy()[m]
    Z = pd.get_dummies(pd.Series(proj)).astype(float).to_numpy()
    if purity is not None:
        Z = np.column_stack([Z, rankdata(purity[m])])
    obs, null = partial_spearman_perm(arm[m], y[m], Z, proj, N_CLASS_PERM, SEED)
    return {"estimate": obs, "p_class": perm_p(obs, null), "testable": m.sum() >= MIN_TUMOUR, "n": f"{int(m.sum())} MMRd tumours"}


def pair_frame(S: State) -> pd.DataFrame:
    pr = pd.concat([r["pairs"].assign(project=p) for p, r in S.res["primary"].items() if "pairs" in r], ignore_index=True)
    lab = S.lab.set_index("case_id")
    pr["mmr"] = lab.mmr.reindex(pr.case_id).to_numpy().astype(bool)
    return pr


def rp7(S: State, col: str) -> dict:
    """RP7 / RPC7: mean over MMRd pairs of (pair arm minus the mean pair arm of the MMRp pairs of the same
    cancer), in the cancers with an MMRd pair; MMR labels permuted among each cancer's pairs."""
    pr = S.pairs_arm
    pr = pr[pr.project.isin(set(pr.project[pr.mmr]))]
    pr = pr[np.isfinite(pr[col].to_numpy(float))]
    Y = label_perms(pr.mmr.to_numpy(), pr.project.to_numpy(), N_CLASS_PERM, SEED)
    st = pair_stat(pr[col].to_numpy(float), Y, pr.project.to_numpy())
    nd = int(pr.mmr.sum())
    return {"estimate": float(st[0]), "p_class": perm_p(st[0], st[1:]), "testable": nd >= MIN_TUMOUR,
            "n": f"{nd} MMRd vs {int((~pr.mmr).sum())} MMRp pairs ({', '.join(sorted(set(pr.project)))})"}


def rp6c(S: State) -> dict:
    """RP6c: D-70's RV1 in JC adj - the held-out MMRd statistic s against every non-MSI subtype contrast."""
    lab, prim = S.lab, SPECS["primary"]
    up = np.flatnonzero(S.ex.cls.eq("R_up").to_numpy())
    A = np.full((len(lab), len(up)), np.nan)
    for p in S.projects:
        A[S.rows_of(prim, p)] = S.adj(prim, p, up)
    rng = np.random.default_rng(SEED)
    mm = np.flatnonzero((lab.mmr & lab.held_out).to_numpy())
    prof = {p: np.flatnonzero((lab.project.eq(p) & ~lab.mmr).to_numpy()) for p in S.projects}
    ia, ib = rp6c_draws(mm, lab.project.to_numpy(), prof, rng)
    s = contrast_draw_stat(A, ia, ib)
    pan = pd.read_parquet(m51.DATA / "tasks" / "cancer_type__PANCAN.parquet").set_index("file_id")
    tasks = {}
    for f in sorted((m51.DATA / "subtype" / "tasks").glob("subtype__*.parquet")):
        t = pd.read_parquet(f)
        t["sample15"] = t.file_id.map(pan.sample_id).str[:15]
        tasks[f.stem.split("__")[1]] = t
    nulls = null_contrasts(tasks, set(lab.sample15[lab.mmr]), set(lab.file_id))
    if not S.run.dry and len(nulls) != RP6C_NULL:
        raise SystemExit(f"RP6c: {len(nulls)} null contrasts, not the pre-registered {RP6C_NULL}")
    null = []
    for cohort, a_, b_, fa, fb in nulls:
        pa, pb = S.pos.loc[fa].to_numpy(), S.pos.loc[fb].to_numpy()
        xa = [rng.choice(pa, RP6C_GROUP, replace=False) for _ in range(RP6C_DRAWS)]
        xb = [rng.choice(pb, RP6C_GROUP, replace=False) for _ in range(RP6C_DRAWS)]
        null.append({"cohort": cohort, "a": a_, "b": b_, "n_a": len(fa), "n_b": len(fb), "s": contrast_draw_stat(A, xa, xb)})
    null = pd.DataFrame(null)
    p = float((1 + np.sum(np.abs(null.s.to_numpy()) >= s - 1e-12)) / (len(null) + 1))
    return {"estimate": s, "p_class": p, "testable": True, "null": null,
            "n": f"{len(null)} null contrasts; max |s| {np.abs(null.s).max():.4f}"}


def direction_shares(folds: list, rng=None) -> tuple[dict, pd.DataFrame]:
    """RPD1's statistic: per input x cohort, the mean over the folds that selected it of the Spearman
    of its standardised value with its SHAP over the fold's MSI test tumours; per arm, the number of
    units with a positive direction and the number with a direction. folds - [{cohort, X, S, fids,
    arm}] restricted to the MSI test tumours and the arm inputs. With rng, the rows of X (all columns
    together) are permuted among each fold's tumours first, independently per fold (A-22's null)."""
    acc = {}
    for f in folds:
        X = f["X"] if rng is None else f["X"][rng.permutation(len(f["X"]))]
        d = shap_directions(X, f["S"], np.arange(len(X)))
        for fid, a, v in zip(f["fids"], f["arm"], d):
            r = acc.setdefault((f["cohort"], fid), [a, 0.0, 0])
            if np.isfinite(v):
                r[1] += v
                r[2] += 1
    t = pd.DataFrame([{"cohort": c, "feature_id": fid, "cls": a, "n_folds_with_direction": n,
                       "direction": s_ / n if n else np.nan} for (c, fid), (a, s_, n) in acc.items()])
    out = {}
    for arm in ["R_up", "R_down"]:
        x = t[t.cls.eq(arm) & t.direction.notna()] if len(t) else t
        out[arm] = (int((x.direction > 0).sum()) if len(x) else 0, len(x))
    return out, t


def rpd1_test(folds: list, n_perm: int = RPD1_PERM, seed: int = SEED) -> dict:
    """RPD1a/b (A-22): the share of positive directions per arm and its permutation p - X rows permuted
    among each fold's MSI test tumours (independently per fold and cohort), directions and shares
    recomputed; p = perm_p(share - 0.5, null share - 0.5). The binomial test the block named treats
    the inputs as independent; their directions are correlated within a fold."""
    obs, table = direction_shares(folds)
    rng = np.random.default_rng(seed)
    null = {a: np.empty(n_perm) for a in obs}
    for b in range(n_perm):
        sh, _ = direction_shares(folds, rng)
        for a, (k, n) in sh.items():
            null[a][b] = k / n if n else np.nan
    out = {"table": table}
    for a, (k, n) in obs.items():
        share = k / n if n else np.nan
        out[a] = {"estimate": share - 0.5, "value": share, "p_class": perm_p(share - 0.5, null[a] - 0.5),
                  "testable": n >= MIN_SHAP, "null_sd": float(np.nanstd(null[a])),
                  "n": f"{n} input x cohort units ({k} positive)"}
    return out


def rpd1(S: State) -> dict:
    """RPD1a/b on the subtype MLPs (SE view; STAD, COAD, UCEC): SHAP explains each tumour's own class,
    so the MSI class's SHAP exists for the MSI test tumours only, and directions are taken over those
    (per fold, from shap.npz for every input). Dry run: the observed X rows are first permuted among
    each fold's MSI tumours (seed 999), so no real direction is computed."""
    cls = load_annotation().set_index("feature_id").cls
    rng = np.random.default_rng(DRY_SEED)
    folds = []
    for cohort in ["COAD_TCGA", "STAD_TCGA", "UCEC_TCGA"]:
        rd = m51.DATA / "subtype" / "runs" / "subtype" / "se" / cohort
        tab = pd.read_parquet(m51.DATA / "subtype" / "tasks" / f"subtype__{cohort}.parquet")
        k_msi = msi_class_index(json.loads((rd / "job_done.json").read_text()), tab)
        y = tab.y.to_numpy()
        for k in range(5):
            z = np.load(rd / f"fold{k}" / "shap.npz")
            msi = np.flatnonzero(y[z["rows"]] == k_msi)
            fids = z["feature_ids"].astype(str)
            arm = cls.reindex(fids).fillna("").to_numpy()
            j = np.flatnonzero(np.isin(arm, ["R_up", "R_down"]))
            X = z["X"].astype(np.float64)[msi][:, j]
            if S.run.dry:
                X = X[rng.permutation(len(X))]
            folds.append({"cohort": cohort, "X": X, "S": z["shap"].astype(np.float64)[msi][:, j], "fids": fids[j], "arm": arm[j]})
    return rpd1_test(folds)


def m4_refit(S: State) -> list:
    """Descriptive: M4 refitted with step 54's code at its chosen C (stop if the held-out AUROC is not
    0.900 +/- 0.001 with the real labels); standardised coefficients by class; M4 without its R_up
    inputs and the held-out AUROC change, paired bootstrap within cancer x MMR status (2,000, seed 0)."""
    m54 = _load("msi54", "54_msi_models.py")
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")
    cols = pd.read_parquet(LOCAL / "models_exon_columns.parquet")
    M = np.load(LOCAL / "models_exon_matrix.npy", mmap_mode="r")
    y = S.lab.set_index("file_id").mmr.reindex(t.file_id).to_numpy().astype(bool)
    held = ~t.project.isin(TRAIN).to_numpy()
    tr1, te1 = np.flatnonzero(~held), np.flatnonzero(held)
    at, db = cols.acc_T.to_numpy(float), cols.d_bp.to_numpy(float)
    _, up, _ = m51.arms(at, db)
    rep = at >= 10
    p4 = m54._chosen(pd.read_csv(m54.OUT / "tuning_e1.csv"), "M4")
    X = np.asarray(M[:, np.flatnonzero(rep)], np.float64)
    s1 = m54.fit_predict("M4", p4, X, tr1, te1, y)
    au1 = m54.auc(s1, y[te1])
    m4_guard(au1, not S.run.dry)
    out = [{"item": "M4 refit: held-out AUROC", "group": "all repeat exons (acc_T >= 10)", "value": au1, "n": len(te1),
            "note": f"C = {p4['C']}; guard {M4_AUROC} +/- {M4_TOL} " + ("not enforced (shuffled labels)" if S.run.dry else "passed")}]
    c4 = m54.coefficients("M4", X, y, tr1, p4, cols.feature_id[rep].tolist()).assign(cls=cols.cls.to_numpy()[rep])
    for k, g in c4.groupby("cls"):
        out.append({"item": "M4 standardised coefficient: median", "group": k, "value": g.coefficient.median(), "n": len(g)})
        out.append({"item": "M4 standardised coefficient: share > 0", "group": k, "value": (g.coefficient > 0).mean(), "n": len(g)})
    keep = np.flatnonzero(~up[np.flatnonzero(rep)])
    s2 = m54.fit_predict("M4", p4, X[:, keep], tr1, te1, y)
    au2 = m54.auc(s2, y[te1])
    rng = np.random.default_rng(SEED)
    pr = pd.DataFrame({"project": t.project.to_numpy()[te1], "mmr": y[te1]})
    strata = [np.asarray(v) for v in pr.groupby(["project", "mmr"]).indices.values()]
    boot = np.empty(m51.BOOT)
    for b in range(m51.BOOT):
        i = np.concatenate([rng.choice(g, len(g)) for g in strata])
        boot[b] = m54.auc(s2[i], y[te1][i]) - m54.auc(s1[i], y[te1][i])
    n_up = int(up[rep].sum())
    out.append({"item": "M4 without its R_up inputs: held-out AUROC", "group": "", "value": au2, "n": len(te1),
                "note": f"{n_up} R_up inputs removed"})
    out.append({"item": "M4 without its R_up inputs: AUROC change", "group": "", "value": au2 - au1, "n": len(te1),
                "ci_low": float(np.percentile(boot, 2.5)), "ci_high": float(np.percentile(boot, 97.5)),
                "note": "paired bootstrap within cancer x MMR status, 2,000, seed 0 - test-set resampling only; "
                        "fits held fixed; not evidence about R_up (A-25); no q"})
    # A-25, exploratory reference: the same number of inputs removed at random from the non-R_up repeat inputs
    from joblib import Parallel, delayed
    others = np.flatnonzero(~up[np.flatnonzero(rep)])
    rng_d = np.random.default_rng(SEED)
    drops = [rng_d.choice(others, n_up, replace=False) for _ in range(M4_RANDOM_DRAWS)]
    chunks = np.array_split(np.arange(M4_RANDOM_DRAWS), N_WORKERS)
    parts = Parallel(n_jobs=N_WORKERS)(delayed(_m4_auc_without)(X, tr1, te1, y, p4, [drops[i] for i in ch])
                                       for ch in chunks)
    au_r = np.concatenate([np.asarray(a, float) for a in parts])
    ch_r = au_r - au1
    pct = float(100 * (np.sum(ch_r < au2 - au1) + 0.5 * np.sum(ch_r == au2 - au1)) / len(ch_r))
    out.append({"item": f"M4 without {n_up} random non-R_up repeat inputs: AUROC change (exploratory, A-25)", "group": "",
                "value": float(np.median(ch_r)), "n": M4_RANDOM_DRAWS,
                "ci_low": float(np.percentile(ch_r, 2.5)), "ci_high": float(np.percentile(ch_r, 97.5)),
                "note": f"median and 2.5-97.5th percentiles over {M4_RANDOM_DRAWS} draws of {n_up} inputs (seed 0), "
                        f"M4 refitted at C = {p4['C']}; no q"})
    out.append({"item": "M4 without its R_up inputs: percentile of the AUROC change among random removals (A-25)",
                "group": "", "value": pct, "n": M4_RANDOM_DRAWS,
                "note": "share of random-removal changes below the R_up removal's change (ties half), x 100; exploratory; no q; "
                        "not a calibrated p - R_up is a correlated block, the reference removes scattered inputs, and on "
                        "shuffled labels the R_up removal already sat at the 2.0th percentile (A-25)"})
    return out


def _m4_auc_without(X, tr1, te1, y, p4, drops: list) -> list:
    """Held-out AUROC of M4 (step 54's fit_predict at the same C) refitted without each draw's inputs
    (A-25's random-removal reference)."""
    m54 = _load("msi54", "54_msi_models.py")
    out = []
    for d in drops:
        keep = np.setdiff1d(np.arange(X.shape[1]), d)
        out.append(m54.auc(m54.fit_predict("M4", p4, X[:, keep], tr1, te1, y), y[te1]))
    return out


# ---------------------------------------------------------------- count model (sensitivity)

def _glm_chunk(I, Sj, R, U, cidx, ga, Y, n_cancers):
    """Per exon column: the quasi-binomial MMRd coefficient and SE for every label column, and whether
    each fit converged (1), did not (0) or was not fitted (-1)."""
    out_b = np.full((I.shape[1], Y.shape[1]), np.nan)
    out_s = np.full((I.shape[1], Y.shape[1]), np.nan)
    conv_m = np.full((I.shape[1], Y.shape[1]), -1, np.int8)
    for e in range(I.shape[1]):
        n = I[:, e] + Sj[:, e]
        keep_c = [c for c in range(n_cancers) if U[c, e] and I[cidx == c, e].sum() > 0]
        m = np.isin(cidx, keep_c) & (n > 0) & np.isfinite(R[:, e]) & (R[:, e] > 0)
        if m.sum() < 10 or not Y[m, 0].any():
            continue
        cc = pd.factorize(cidx[m])[0]
        D_ = np.zeros((m.sum(), cc.max() + 1))
        D_[np.arange(m.sum()), cc] = 1
        base = np.column_stack([D_, ga[m]])
        off = np.log(R[m, e])
        b0 = None
        for k in range(Y.shape[1]):
            X = np.column_stack([base, Y[m, k].astype(float)])
            try:
                beta, cov, _, conv = quasi_binomial(I[m, e], n[m], off, X, beta0=b0)
            except np.linalg.LinAlgError:
                continue
            if k == 0:
                b0 = beta
            conv_m[e, k] = int(conv)
            out_b[e, k], out_s[e, k] = beta[-1], np.sqrt(max(cov[-1, -1], 0))
    return out_b, out_s, conv_m


def count_model(S: State, ctx: XCtx, sel: np.ndarray) -> dict:
    """Per selected eligible held-out exon: quasi-binomial GLM of IJC out of IJC + SJC over the tumours
    of its universe cancers (cancers without any inclusion read dropped), offset log(IncFormLen /
    SkipFormLen), cancer fixed effects, global axis (percentile) and MMRd, for the observed labels
    and the first 200 label permutations."""
    from joblib import Parallel, delayed
    cols = ctx.cols[sel]
    I, Sj, R, cidx, ga, Ys = [], [], [], [], [], []
    held = [p for p in S.held]
    Umat = np.vstack([S.res["primary"][p]["U"][cols] for p in held])
    for c, p in enumerate(held):
        samp, ijc, sjc, ratio = load_counts_all(p, "jc")
        tum = samp.role.eq("tumour").to_numpy()
        I.append(np.asarray(ijc[tum][:, cols], np.float64))
        Sj.append(np.asarray(sjc[tum][:, cols], np.float64))
        R.append(np.asarray(ratio[tum][:, cols], np.float64))
        cidx.append(np.full(tum.sum(), c))
        ga.append(S.tv.ga_pct.to_numpy(float)[S.rows_of(SPECS["primary"], p)])
        Ys.append(S.Yp(SPECS["primary"], p)[:, :N_COUNT_PERM + 1])
    I, Sj, R, cidx, ga, Ys = (np.concatenate(v) for v in (I, Sj, R, cidx, ga, Ys))
    n_imputed = int(np.isnan(ga).sum())                  # no measured B-sample exon: global axis set to 0.5
    ga = np.nan_to_num(ga, nan=0.5)
    chunks = np.array_split(np.arange(len(cols)), max(1, len(cols) // 40))
    res = Parallel(n_jobs=N_WORKERS)(delayed(_glm_chunk)(I[:, ch], Sj[:, ch], R[:, ch], Umat[:, ch], cidx, ga, Ys, len(held))
                                     for ch in chunks)
    B = np.vstack([r[0] for r in res])
    SE = np.vstack([r[1] for r in res])
    conv = np.vstack([r[2] for r in res])
    return {"beta": B, "se": SE, "conv": conv, "converged": int((conv[:, 0] == 1).sum()), "n": len(cols),
            "ga_imputed": n_imputed, "n_tumours": len(ga)}


def separated(beta: np.ndarray, se: np.ndarray, limit: float = SEP_LIMIT) -> np.ndarray:
    """A-23: count-model fits on a separation plateau (|coefficient| > 10 or SE > 10, logit units)."""
    with np.errstate(invalid="ignore"):
        return (np.abs(beta) > limit) | (se > limit)


def ivw_contrast(beta: np.ndarray, se: np.ndarray, first: np.ndarray, ctrl: np.ndarray, lab: np.ndarray) -> np.ndarray:
    """Per column of beta [exons x columns]: the stratum-weighted difference of inverse-variance-weighted
    class means (weights = first-class exons per stratum); ctrl None -> the first class's IVW mean."""
    w = np.where(np.isfinite(beta) & np.isfinite(se) & (se > 0), 1 / se ** 2, 0)
    b0 = np.where(w > 0, beta, 0)
    if ctrl is None:
        f = first[:, None] * w
        return (f * b0).sum(0) / f.sum(0)
    num, den = 0, 0
    for s_ in sorted(set(lab[first].tolist())):
        F, C = (first & (lab == s_))[:, None] * w, (ctrl & (lab == s_))[:, None] * w
        with np.errstate(all="ignore"):
            d = (F * b0).sum(0) / F.sum(0) - (C * b0).sum(0) / C.sum(0)
        nf = int((first & (lab == s_)).sum())
        ok = np.isfinite(d)
        num = num + np.where(ok, d, 0) * nf
        den = den + ok * nf
    with np.errstate(all="ignore"):
        return num / den


def ivw_class_perm(beta, se, first, ctrl, lab, n_perm: int = N_CLASS_PERM, seed: int = SEED) -> np.ndarray:
    """Null of ivw_contrast (observed labels) with class labels permuted within each stratum."""
    rng = np.random.default_rng(seed)
    w = np.where(np.isfinite(beta) & np.isfinite(se) & (se > 0), 1 / se ** 2, 0)
    b0 = np.where(w > 0, beta, 0)
    num, den = np.zeros(n_perm), 0
    for s_ in sorted(set(lab[first].tolist())):
        ii = np.flatnonzero((first | ctrl) & (lab == s_))
        f = int(first[ii].sum())
        if f == 0 or f == len(ii):
            continue
        pick = rng.random((n_perm, len(ii))).argsort(axis=1)[:, :f]
        wb, ww = (w * b0)[ii], w[ii]
        sf, wf = wb[pick].sum(1), ww[pick].sum(1)
        with np.errstate(all="ignore"):
            d = sf / wf - (wb.sum() - sf) / (ww.sum() - wf)
        num += np.where(np.isfinite(d), d, 0) * f
        den += f
    return num / den if den else np.full(n_perm, np.nan)


# ---------------------------------------------------------------- compute driver

def check_strata(S: State) -> None:
    """The generalised side tables and coarsening reproduce the label-blind Phase 1 tables."""
    st = pd.read_parquet(PROOF / "exon_strata.parquet")
    n = len(st)
    for side, thr in [("held_out", 20), ("training", 1)]:
        mine = S.side(SPECS["primary"])[(side, thr)].iloc[:n]
        for c in ["eligible", "band", "tertile"]:
            a, b = mine[c].to_numpy(), st[f"{side}_{c}"].to_numpy()
            if not np.array_equal(a.astype(float), b.astype(float)):
                raise SystemExit(f"side table {side} {c} differs from exon_strata.parquet ({int((a != b).sum())} exons)")
    ref = pd.read_csv(OUT / "strata_contrasts.csv").set_index("contrast")
    ctx = S.xctx(SPECS["primary"], "held_out")
    for cid, side, first, ctrl, per in contrast_masks(ctx.tab):
        if side != "held_out" or cid not in ref.index:
            continue
        f, c = first.to_numpy(bool) & ctx.elig, ctrl.to_numpy(bool) & ctx.elig
        lab, kept = coarsen(ctx.band, ctx.terc, f, c, per_first=per)
        e = effective_n(lab, f, c, kept)
        if e["first_kept"] != ref.loc[cid, "first_kept"] or abs(e["control_effective_n"] - ref.loc[cid, "control_effective_n"]) > 0.051:
            raise SystemExit(f"coarsening of {cid} does not reproduce strata_contrasts.csv: {e}")
    log("side tables and coarsened strata reproduce Phase 1 (exon_strata.parquet, strata_contrasts.csv)")


def exon_tests(S: State, ctx: XCtx, prefix: str = "", only: set | None = None) -> dict:
    """The exon-level tests on one read-out; returns their result dicts by id."""
    t = ctx.tab
    up, dn = (t.cls.eq("R_up")).to_numpy(), (t.cls.eq("R_down")).to_numpy()
    bs = t.b_sample.to_numpy(bool)
    long_ = (t.acc_T >= LONG).to_numpy()
    exp_up = long_ & (t.g_exp.to_numpy(float) > 0)
    near = (t.cls.eq("long_near") & t.d_bp.notna()).to_numpy()
    R = {}

    def want(i):
        return only is None or i in only
    if want("RP2a"):
        R["RP2a"] = class_mean(S, ctx, up & ctx.elig)
    if want("RPC2a"):
        R["RPC2a"] = class_mean(S, ctx, dn & ctx.elig)
    if want("RP2c"):
        R["RP2c"] = rp2c(ctx)
    pairs = {"RP3a": (up, bs & t.cls.eq("B_far").to_numpy()), "RP3b": (up, near),
             "RP3d": (exp_up & (t.d_exp.to_numpy(float) >= FAR), exp_up & (t.d_exp.to_numpy(float) < FAR)),
             "RP4a": (up, t.ctl_polyA15.to_numpy(bool)), "RP4b": (up, t.ctl_ppt_far.to_numpy(bool)),
             "RP4c": (up, t.ctl_deep.to_numpy(bool)), "RP4d": (up, t.ctl_don.to_numpy(bool)),
             "RP4e": (up, t.ctl_alu.to_numpy(bool)), "RPM1": (up, t.shield.to_numpy(bool))}
    for i, (f, c) in pairs.items():
        if want(i):
            R[i] = contrast(ctx, f, c, min_n=MIN_RP3D if i == "RP3d" else MIN_CLASS)
    L, o = t.geo_L.fillna(t.acc_T).to_numpy(float), t.geo_o.fillna(t.run_near).to_numpy(float)
    if want("RP3c"):
        R["RP3c"] = geometry(ctx, long_ & t.geo_g.notna().to_numpy(), t.geo_g.to_numpy(float), L, o, L_BINS["RP3c"], O_BINS["RP3c"])
    if want("RP3e"):
        R["RP3e"] = geometry(ctx, exp_up, t.g_exp.to_numpy(float), L, o, L_BINS["RP3e"], O_BINS["RP3e"], min_n=MIN_RP3E)
    if want("RP8c"):
        R["RP8c"] = spearman_within(ctx, up, t.acc_T.to_numpy(float), ctx.band.astype(str))
    if want("RPM2"):
        sh = t.shield.to_numpy(bool) & ctx.elig
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            st = np.nanmean(ctx.D[:, sh], axis=1)
        R["RPM2"] = {"estimate": float(st[0]), "p_label": perm_p(st[0], st[1:]),
                     "testable": int(np.isfinite(ctx.D[0, sh]).sum()) >= MIN_CLASS, "n": f"{int(np.isfinite(ctx.D[0, sh]).sum())} exons"}
    if want("RPM3"):
        ant = t.run_alu_antisense.eq(True).to_numpy()
        labs, keeps = np.full(len(t), "", object), np.zeros(len(t), bool)
        for lo, hi in [(20, 25), (25, 10 ** 6)]:
            b_ = long_ & (t.acc_T.to_numpy(float) >= lo) & (t.acc_T.to_numpy(float) < hi)
            lab, kept = coarsen(ctx.band, ctx.terc, b_ & ant & ctx.elig, b_ & ~ant & ctx.elig)
            sel = b_ & (lab != "")
            labs[sel] = [f"L{lo}|{x}" for x in lab[sel]]
            keeps |= kept & b_
        R["RPM3"] = contrast(ctx, long_ & ant, long_ & ~ant, lab_kept=(labs, keeps))
    return R


def compute(run: Run, fresh: bool = False, quick: bool = False) -> None:
    S = State(run, fresh)
    prim = SPECS["primary"]
    S.build(prim, S.projects)
    S.tick("build: primary (33 cancers)")
    S.tv = tumour_values(S)
    S.pairs_arm = pair_frame(S)
    for k in ["jcec", "reads30", "univ80", "pct"]:
        S.build(SPECS[k], S.held)
        S.tick(f"build: {k} (30 held-out cancers)")
    check_strata(S)
    S.tick("side tables")
    H = S.xctx(prim, "held_out")
    TR = S.xctx(prim, "training")
    S.tick("engine: primary held out and training")
    lab = S.lab
    held = lab.held_out.to_numpy()
    arm_up, arm_dn = S.tv.arm_up.to_numpy(float), S.tv.arm_down.to_numpy(float)

    # ---- confirmatory, secondary, controls and mechanism on the primary read-out
    R = exon_tests(S, H)
    S.tick("exon-level tests (primary)")
    instr = instrument_check(S, H, R["RPC2a"]["z"][1:])
    S.tick("instrument check")
    for i in ["RP2a", "RPC2a"]:
        r = R[i]
        S.add(i, r["T"][0], p_label=r["p"], ci=r["ci"][:2], statistic=r["z"][0], testable=r["n_exons_obs"] >= MIN_CLASS,
              n=f"{r['n_exons_obs']} exons, {int(r['N'][0])} MMRd tumours, {int(r['K'][0])} cancers; "
                f"Welch-Satterthwaite df {r['ci'][2]:.1f}", inference="label_perm, studentized")
        D0 = np.nanmean(H.D[0, (H.tab.cls.eq('R_up' if i == 'RP2a' else 'R_down')).to_numpy() & H.elig])
        assert abs(D0 - r["T"][0]) < 1e-5, (i, D0, r["T"][0])
    for i in ["RP2c", "RP3a", "RP3b", "RP3c", "RP3d", "RP3e", "RP4a", "RP4b", "RP4c", "RP4d", "RP4e", "RP8c", "RPM1", "RPM3"]:
        r = R[i]
        S.add(i, r["estimate"], p_class=r["p_class"], p_label=r["p_label"], testable=r["testable"], n=r["n"],
              inference="class_perm / exon permutation and label_perm, p = max")
    S.add("RPM2", R["RPM2"]["estimate"], p_label=R["RPM2"]["p_label"], testable=R["RPM2"]["testable"],
          n=R["RPM2"]["n"], inference="label_perm")
    for i, cls in [("RP5a", "R_up"), ("RPC5a", "R_down")]:
        r = rp5a(S, H, cls)
        R[i] = r
        S.add(i, r["estimate"], p_class=r["p_class"], p_label=r["p_label"], statistic=r["statistic"], testable=r["testable"],
              n=r["n"], inference="Freedman-Lane (exon residuals within band x tertile) and label_perm on the HC3 t, p = max")
    S.tick("RP5a / RPC5a")
    for i, arm in [("RP1", "R_up"), ("RPC1", "R_down")]:
        r = rp1(S, arm)
        S.add(i, r["estimate"], p_class=r["p_class"], p_label=r["p_label"], value=r["value"], testable=r["testable"],
              n=r["n"], inference="matched sets and label_perm, p = max")
    S.tick("RP1 / RPC1")
    for i, a in [("RP2b", arm_up), ("RPC2b", arm_dn)]:
        r = arm_test(S, a, held)
        R[i] = r
        S.add(i, r["estimate"], p_label=r["p_label"], statistic=r["statistic"], testable=r["testable"], n=r["n"],
              inference="label_perm, studentized")
    for i, a in [("RP5b", arm_up), ("RPC5b", arm_dn)]:
        r = rp5b(S, a)
        S.add(i, r["estimate"], p_class=r["p_class"], statistic=r["statistic"], testable=r["testable"], n=r["n"],
              inference="Freedman-Lane, tumour residuals within cancer, HC3 t")
    for i, a in [("RP6a", arm_up), ("RPC6", arm_dn)]:
        r = rp6a(S, a)
        S.add(i, r["estimate"], p_class=r["p_class"], value=r.get("value", np.nan), testable=r["testable"], n=r["n"],
              inference="Mann-Whitney")
    r = msi_call_test(S, arm_up, "call_mantis")
    S.add("RP6b", r["estimate"], p_label=r["p_label"], statistic=r["statistic"], testable=r["testable"], n=r["n"],
          inference="label_perm with the MANTIS call, studentized")
    r6c = rp6c(S)
    S.add("RP6c", r6c["estimate"], p_class=r6c["p_class"], testable=True, n=r6c["n"],
          inference="rank against the null contrasts (p floor 1/58)")
    S.tick("tumour-level tests")
    for i, col in [("RP7", "pair_up"), ("RPC7", "pair_down")]:
        r = rp7(S, col)
        S.add(i, r["estimate"], p_class=r["p_class"], testable=r["testable"], n=r["n"],
              inference="MMR labels permuted among each cancer's pairs (10,000)")
    mantis = lab.mantis.to_numpy(float)
    mm_h = lab.mmr.to_numpy() & held
    for i, a in [("RP8a", arm_up), ("RPC8", arm_dn)]:
        r = dose(S, a, mantis, mm_h)
        S.add(i, r["estimate"], p_class=r["p_class"], testable=r["testable"], n=r["n"], inference="MANTIS permuted within cancer")
    r = dose(S, arm_up, np.log1p(lab.n_hp.to_numpy(float)), lab.mmr.to_numpy())
    S.add("RP8b", r["estimate"], p_class=r["p_class"], testable=r["testable"], n=r["n"],
          inference="log(1 + deletions) permuted within cancer")
    S.tick("pairs and dose")

    # ---- discovery (training)
    for i, cls in [("RPD2a", "R_up"), ("RPD2b", "R_down")]:
        r = class_mean(S, TR, TR.tab.cls.eq(cls).to_numpy() & TR.elig)
        S.add(i, r["T"][0], p_label=r["p"], ci=r["ci"][:2], statistic=r["z"][0], testable=r["n_exons_obs"] >= MIN_CLASS,
              n=f"{r['n_exons_obs']} exons, {int(r['N'][0])} MMRd tumours", inference="label_perm, studentized")
    RT = exon_tests(S, TR, only={"RP3a", "RP3e"})
    S.add("RPD3a", RT["RP3a"]["estimate"], p_class=RT["RP3a"]["p_class"], p_label=RT["RP3a"]["p_label"],
          testable=RT["RP3a"]["testable"], n=RT["RP3a"]["n"], inference="class_perm and label_perm, p = max")
    S.add("RPD3b", RT["RP3e"]["estimate"], p_class=RT["RP3e"]["p_class"], p_label=RT["RP3e"]["p_label"],
          testable=RT["RP3e"]["testable"], n=RT["RP3e"]["n"], inference="exon permutation in cells and label_perm, p = max")
    d1 = rpd1(S)
    for i, arm in [("RPD1a", "R_up"), ("RPD1b", "R_down")]:
        r = d1[arm]
        S.add(i, r["estimate"], p_class=r["p_class"], value=r["value"], testable=r["testable"], n=r["n"],
              inference="permutation of X rows among each fold's MSI test tumours (1,000; A-22)",
              note=f"null SD of the share {r['null_sd']:.4f}")
    S.tick("discovery tests")

    # ---- descriptive
    S.add_desc("instrument check: critical value (2.5th percentile of RPC2a's studentized null)", "", instr["critical_value"])
    S.add_desc("instrument check: median draw statistic", "", instr["median_draw"], n=N_DRAWS)
    S.add_desc("instrument check: share of draws below the critical value", "", instr["share_below"], n=N_DRAWS)
    S.add_desc("instrument check: passes", "", float(instr["passes"]),
               note=f"{instr['n_up_drawn']} R_up exons matched; pool {instr['r_down_pool']} R_down exons; "
                    f"{instr['r_up_mirrored_top_band']} R_up exons have their mirrored baseline in [0.95, 1], and a draw "
                    f"holds on average {instr['drawn_r_down_baseline_ge_095_mean']:.1f} R_down exons with baseline >= 0.95")
    S.add_desc("RP2c: the rise cut c (97.5th percentile of eligible B-sample studentized deltas over the label permutations)",
               "", R["RP2c"]["cut"], n=R["RP2c"]["n_sets"],
               note="per-exon studentized deltas are heavy-tailed (few MMRd tumours per exon), so c lies well above 1.96")
    lo, hi, df = R["RP2a"]["ci"]
    S.add_desc("RP2a 95% interval against the SESOI", f"SESOI {SESOI}", R["RP2a"]["T"][0], ci_low=lo, ci_high=hi,
               note=f"Welch-Satterthwaite df {df:.1f}; upper bound {'>=' if hi >= SESOI else '<'} SESOI")
    try:
        for d in m4_refit(S):
            S.add_desc(d.pop("item"), d.pop("group"), d.pop("value"), n=d.pop("n", ""), **d)
    except SystemExit:
        raise
    S.tick("M4 refit")
    descriptives(S, H)
    S.tick("descriptives")

    # ---- sensitivity
    sensitivity(S, H, R, quick=quick)
    S.tick("sensitivity")

    # ---- save
    pd.DataFrame(S.rows).to_csv(run.path("raw_tests.csv", local=True), index=False)
    pd.DataFrame(S.sens).to_csv(run.path("raw_sensitivity.csv", local=True), index=False)
    pd.DataFrame(S.desc).to_csv(run.path("raw_descriptive.csv", local=True), index=False)
    pd.DataFrame({"draw": np.arange(N_DRAWS), "estimate": instr["draw_estimate"], "studentized": instr["draw_z"]}).to_csv(
        run.path("instrument_draws.csv"), index=False, float_format="%.6g")
    inst = {k: v for k, v in instr.items() if k not in ("draw_z", "draw_estimate")}
    run.path("instrument.json", local=True).write_text(json.dumps(inst, indent=1, default=float))
    r6c["null"].to_csv(run.path("rp6c_null_contrasts.csv"), index=False, float_format="%.6g")
    d1["table"].to_parquet(run.path("rpd1_directions.parquet", local=True), index=False)
    tv = pd.concat([lab[["file_id", "project", "mmr", "hyper"]].reset_index(drop=True), S.tv], axis=1)
    tv.to_parquet(run.path("tumours.parquet", local=True), index=False)
    S.pairs_arm.to_parquet(run.path("pairs.parquet", local=True), index=False)
    run.path("timing.json").write_text(json.dumps({"labels": run.labels, **S.timing}, indent=1))
    log(f"compute done in {time.time() - S.t0:.0f} s")


def descriptives(S: State, H: XCtx) -> None:
    """Descriptive items: mean delta by L bin x d bin (heuristic and experimental), each control class's
    own mean delta (label_perm p), per-cancer RP2b (>= 3 MMRd), the largest R_up movers, mean delta by
    baseline band for every class."""
    t = H.tab
    long_ = (t.acc_T >= LONG).to_numpy() & H.elig
    L = t.geo_L.fillna(t.acc_T).to_numpy(float)
    Lb = np.digitize(L, L_BINS["RP3c"])
    for kind, d in [("heuristic d", t.d_bp.to_numpy(float)), ("experimental d_exp", np.where(t.g_exp > 0, t.d_exp, np.nan).astype(float))]:
        db = np.digitize(d, D_BINS)
        for li in range(1, len(L_BINS["RP3c"]) + 1):
            for di in range(1, len(D_BINS) + 1):
                m = long_ & (Lb == li) & (db == di) & np.isfinite(H.D[0])
                if m.any():
                    S.add_desc(f"mean delta, runs >= 20 by L bin x {kind} bin", f"L bin {li} (from {L_BINS['RP3c'][li - 1]}), "
                               f"d bin {di} (from {D_BINS[di - 1]})", float(H.D[0, m].mean()), n=int(m.sum()))
    ctl = {"B_far (B sample)": t.b_sample.to_numpy(bool) & t.cls.eq("B_far").to_numpy(), "near_bp": (t.cls.eq("long_near") & t.d_bp.notna()).to_numpy(),
           "polyA15": t.ctl_polyA15.to_numpy(bool), "ppt_far": t.ctl_ppt_far.to_numpy(bool), "deep": t.ctl_deep.to_numpy(bool),
           "donor": t.ctl_don.to_numpy(bool), "alu": t.ctl_alu.to_numpy(bool), "shield": t.shield.to_numpy(bool),
           "B sample (all)": t.b_sample.to_numpy(bool)}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for k, m in ctl.items():
            st = np.nanmean(H.D[:, m & H.elig], axis=1)
            S.add_desc("control class: own mean delta against 0 (label_perm p)", k, float(st[0]),
                       n=int(np.isfinite(H.D[0, m & H.elig]).sum()), p=perm_p(st[0], st[1:]))
        classes = {"R_up": t.cls.eq("R_up").to_numpy(), "R_down": t.cls.eq("R_down").to_numpy(), **ctl}
        for k, m in classes.items():
            for b in range(len(BAND_EDGES) - 1):
                mm = m & H.elig & (H.band == b) & np.isfinite(H.D[0])
                if mm.any():
                    S.add_desc("mean delta by baseline band", f"{k}, band {b} [{BAND_EDGES[b]}, {min(BAND_EDGES[b + 1], 1)})",
                               float(H.D[0, mm].mean()), n=int(mm.sum()))
    lab, held = S.lab, S.lab.held_out.to_numpy()
    for p in S.held:
        sel = held & lab.project.eq(p).to_numpy()
        if lab.mmr.to_numpy()[sel].sum() >= 3:
            r = arm_test(S, S.tv.arm_up.to_numpy(float), sel)
            S.add_desc("per-cancer RP2b (R_up arm_contrast, >= 3 MMRd)", p, r["estimate"], n=r["n"], p=r["p_label"],
                       n_mmrd=int(r["out"]["N"][0]))
    up = np.flatnonzero(t.cls.eq("R_up").to_numpy() & H.elig & np.isfinite(H.D[0]))
    for k, i in enumerate(up[np.argsort(-H.D[0, up])][:20], 1):
        S.add_desc("largest R_up movers (delta, PSI units)", f"{k}. {t.symbol.iat[i]} {t.feature_id.iat[i]}", float(H.D[0, i]),
                   n=int(H.N[0, i]), note=f"baseline {H.base[i]:.3f}; acc_T {t.acc_T.iat[i]:.0f}; d_bp {t.d_bp.iat[i]}")


def sensitivity(S: State, H: XCtx, R: dict, quick: bool = False) -> None:
    """The pre-specified sensitivity analyses (no q)."""
    t = H.tab
    up, dn = t.cls.eq("R_up").to_numpy(), t.cls.eq("R_down").to_numpy()
    bs = t.b_sample.to_numpy(bool)
    tv, lab = S.tv, S.lab
    held = lab.held_out.to_numpy()

    def cm(item, ctx, sel, label, **kw):
        r = class_mean(S, ctx, sel, **kw)
        S.add_sens(item, label, r["T"][0], r["p"], n=f"{r['n_exons_obs']} exons", ci_low=r["ci"][0], ci_high=r["ci"][1])

    def ct(item, label, r):
        S.add_sens(item, label, r["estimate"], max(r["p_class"], r["p_label"]), n=r["n"], p_class=r["p_class"], p_label=r["p_label"])

    # JCEC for RP2a-RP4e and RPC2a-b
    J = S.xctx(SPECS["jcec"], "held_out")
    RJ = exon_tests(S, J, only={"RP2a", "RPC2a", "RP2c", "RP3a", "RP3b", "RP3c", "RP3d", "RP3e", "RP4a", "RP4b", "RP4c", "RP4d", "RP4e"})
    for i in ["RP2a", "RPC2a"]:
        S.add_sens("JCEC counts", i, RJ[i]["T"][0], RJ[i]["p"], n=f"{RJ[i]['n_exons_obs']} exons")
    for i in ["RP2c", "RP3a", "RP3b", "RP3c", "RP3d", "RP3e", "RP4a", "RP4b", "RP4c", "RP4d", "RP4e"]:
        ct("JCEC counts", i, RJ[i])
    jt = tumour_values_spec(S, SPECS["jcec"])
    for i, a in [("RP2b", jt[0]), ("RPC2b", jt[1])]:
        r = arm_test(S, a, held)
        S.add_sens("JCEC counts", i, r["estimate"], r["p_label"], n=r["n"])
    S.tick("sensitivity: JCEC")
    # cells >= 30 reads; universe at 80%
    for k, item in [("reads30", "cells with JC reads >= 30"), ("univ80", "universe at 80%")]:
        C_ = S.xctx(SPECS[k], "held_out")
        RC = exon_tests(S, C_, only={"RP2a", "RPC2a", "RP3a"})
        for i in ["RP2a", "RPC2a"]:
            S.add_sens(item, i, RC[i]["T"][0], RC[i]["p"], n=f"{RC[i]['n_exons_obs']} exons")
        ct(item, "RP3a", RC["RP3a"])
        a_up, a_dn = tumour_values_spec(S, SPECS[k])
        for i, a in [("RP2b", a_up), ("RPC2b", a_dn)]:
            r = arm_test(S, a, held)
            S.add_sens(item, i, r["estimate"], r["p_label"], n=r["n"])
    # eligibility at >= 10 and >= 30 MMRd tumours
    for thr in (10, 30):
        C_ = H.at(S.side(SPECS["primary"])[("held_out", thr)])
        RC = exon_tests(S, C_, only={"RP2a", "RPC2a", "RP3a"})
        for i in ["RP2a", "RPC2a"]:
            S.add_sens(f"eligibility at >= {thr} MMRd tumours", i, RC[i]["T"][0], RC[i]["p"], n=f"{RC[i]['n_exons_obs']} exons")
        ct(f"eligibility at >= {thr} MMRd tumours", "RP3a", RC["RP3a"])
    S.tick("sensitivity: reads, universe, eligibility")
    # baseline >= 0.90 left out
    lo_base = H.base < 0.90
    cm("exons with baseline >= 0.90 left out", H, up & H.elig & lo_base, "RP2a")
    cm("exons with baseline >= 0.90 left out", H, dn & H.elig & lo_base, "RPC2a")
    Hb = XCtx(H.spec, H.side, H.cols, H.tab, H.D, H.T, H.N, H.projects, H.elig & lo_base, H.band, H.terc, H.base, H.reads)
    ct("exons with baseline >= 0.90 left out", "RP3a", contrast(Hb, up, bs & t.cls.eq("B_far").to_numpy()))
    # gene level; non-poison; without BRCA / CESC; median; percentile units
    cm("gene level (delta averaged within gene first)", H, up & H.elig, "RP2a",
       weight_fn=gene_weights(t.gene_id.to_numpy()[up & H.elig]))
    cm("R_up exons that are not D-76 ptc_nmd (poison)", H, up & H.elig & ~t.consequence.eq("ptc_nmd").to_numpy(), "RP2a")
    for c in ["BRCA", "CESC"]:
        cm(f"without {c}", H, up & H.elig, "RP2a", drop=(c,))
        r = arm_test(S, tv.arm_up.to_numpy(float), held & ~lab.project.eq(c).to_numpy())
        S.add_sens(f"without {c}", "RP2b", r["estimate"], r["p_label"], n=r["n"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(H.D[:, up & H.elig], axis=1)
    S.add_sens("median instead of mean delta", "RP2a", float(med[0]), perm_p(med[0], med[1:]), n=f"{int((up & H.elig).sum())} exons",
               note="label_perm on the median (not studentized)")
    P_ = S.xctx(SPECS["pct"], "held_out")
    cm("D-66 percentile units on JC", P_, (P_.tab.cls.eq("R_up")).to_numpy() & P_.elig, "RP2a")
    S.tick("sensitivity: RP2a variants")
    # RP3b with R_up cut at >= 50 and > 60; Spearman of delta with d; no_bp against B_far
    near = (t.cls.eq("long_near") & t.d_bp.notna()).to_numpy()
    long_ = (t.acc_T >= LONG).to_numpy()
    d = t.d_bp.to_numpy(float)
    ct("R_up cut at d_bp >= 50", "RP3b", contrast(H, long_ & (d >= 50), near))
    ct("R_up cut at d_bp > 60", "RP3b", contrast(H, long_ & (d > 60), near))
    r = spearman_within(H, long_ & np.isfinite(d), d, H.band.astype(str))
    ct("Spearman of delta with d, runs >= 20 with a heuristic branch point (d permuted within baseline band)", "RP3b", r)
    ct("no_bp against B_far", "RP3b", contrast(H, long_ & t.d_bp.isna().to_numpy(), bs & t.cls.eq("B_far").to_numpy()))
    # nearest experimental branch point
    gn, dn_ = t.g_exp_near.to_numpy(float), t.d_exp_near.to_numpy(float)
    en = long_ & (gn > 0)
    ct("nearest experimental branch point", "RP3d", contrast(H, en & (dn_ >= FAR), en & (dn_ < FAR), min_n=MIN_RP3D))
    L, o = t.geo_L.fillna(t.acc_T).to_numpy(float), t.geo_o.fillna(t.run_near).to_numpy(float)
    ct("nearest experimental branch point", "RP3e", geometry(H, en, gn, L, o, L_BINS["RP3e"], O_BINS["RP3e"], min_n=MIN_RP3E))
    # control variants
    ct("polyA >= 20 nt", "RP4a", contrast(H, up, t.ctl_polyA.to_numpy(bool)))
    ct("R_up runs of 20-24 nt against polyA", "RP4a", contrast(H, up & (t.acc_T < 25).to_numpy(), t.ctl_polyA15.to_numpy(bool)))
    ct("deep restricted to a branch point >= 36 nt", "RP4c", contrast(H, up, t.ctl_deep_far.to_numpy(bool)))
    ct("donor restricted to a branch point >= 36 nt", "RP4d", contrast(H, up, t.ctl_don_far.to_numpy(bool)))
    for L_ in (15, 25):
        ct(f"ppt_far at >= {L_} nt", "RP4b", contrast(H, up, t[f"ctl_ppt_far{L_}"].to_numpy(bool)))
    ct("within the antisense-Alu subset", "RP4e", contrast(H, up & t.acc_alu_antisense.eq(True).to_numpy(), t.ctl_alu_antisense.eq(True).to_numpy()))
    na = ~t.acc_alu.eq(True).to_numpy()
    for i, c in [("RP3a", bs & t.cls.eq("B_far").to_numpy()), ("RP4a", t.ctl_polyA15.to_numpy(bool)), ("RP4b", t.ctl_ppt_far.to_numpy(bool)),
                 ("RP4c", t.ctl_deep.to_numpy(bool)), ("RP4d", t.ctl_don.to_numpy(bool))]:
        ct("acceptor not in an Alu", i, contrast(H, up & na, c & na))
    S.tick("sensitivity: class variants")
    # RP5a full model; RP5a with cds_max_run
    hc = R["RP5a"]["hc"]
    r = rp5a(S, H, "R_up", extra=("pwm3", "acc_alu"), hc=hc)
    ct("RP5a's full model (+ -20..+3 PWM, acceptor Alu)", "RP5a", r)
    r = rp5a(S, H, "R_up", extra=("cds_max_run",), hc=hc)
    ct("RP5a with cds_max_run added", "RP5a", r)
    # RP5b without purity; with MKI67
    r = rp5b(S, tv.arm_up.to_numpy(float), purity=False)
    S.add_sens("RP5b without purity (all held-out tumours)", "RP5b", r["estimate"], r["p_class"], n=r["n"])
    r = rp5b(S, tv.arm_up.to_numpy(float), extra=mki67_pct(S))
    S.add_sens("RP5b with MKI67 log2(CPM + 1) percentile within cancer", "RP5b", r["estimate"], r["p_class"], n=r["n"])
    # RP6a within cancer; RP6b with MSIsensor; RP8b with purity
    r = rp6a(S, tv.arm_up.to_numpy(float), within=True)
    S.add_sens("RP6a within cancer", "RP6a", r["estimate"], r["p_class"], n=r["n"])
    r = msi_call_test(S, tv.arm_up.to_numpy(float), "call_msisensor")
    S.add_sens("RP6b with MSIsensor MSI-H", "RP6b", r["estimate"], r["p_label"], n=r["n"])
    pur = lab.case_id.map(m50.purity_by_patient()).to_numpy(float)
    r = dose(S, tv.arm_up.to_numpy(float), np.log1p(lab.n_hp.to_numpy(float)), lab.mmr.to_numpy(), purity=pur)
    S.add_sens("RP8b with purity", "RP8b", r["estimate"], r["p_class"], n=r["n"])
    S.tick("sensitivity: tumour-level variants")
    # count model
    if not quick:
        cmodel(S, H)
        S.tick("sensitivity: count model")


def tumour_values_spec(S: State, spec: Spec) -> tuple[np.ndarray, np.ndarray]:
    up, dn = np.full(len(S.lab), np.nan), np.full(len(S.lab), np.nan)
    for p, r in S.res[spec.name].items():
        i = S.rows_of(spec, p)
        up[i], dn[i] = r["arm_up"], r["arm_down"]
    return up, dn


def mki67_pct(S: State) -> np.ndarray:
    """MKI67 log2(CPM + 1) as a percentile within cancer, held-out tumours."""
    from ctsm.matrix_io import load_view
    cfg = load_config()
    out = np.full(len(S.lab), np.nan)
    for p in S.held:
        i = np.flatnonzero(S.lab.project.eq(p).to_numpy())
        G, F = load_view(cfg.matrix_dir, "gex", S.lab.file_id.to_numpy()[i])
        j = int(np.flatnonzero(F.symbol.eq("MKI67").to_numpy())[0])
        out[i] = pct_vector(np.asarray(G.get()[:, j], float))
    return out


def cmodel(S: State, H: XCtx) -> None:
    """The count-level sensitivity: RP2a, RP3a, RP3b and RP4a-RP4e on the quasi-binomial MMRd coefficient."""
    t = H.tab
    up = t.cls.eq("R_up").to_numpy()
    bs = t.b_sample.to_numpy(bool)
    ctrls = {"RP3a": bs & t.cls.eq("B_far").to_numpy(), "RP3b": (t.cls.eq("long_near") & t.d_bp.notna()).to_numpy(),
             "RP4a": t.ctl_polyA15.to_numpy(bool), "RP4b": t.ctl_ppt_far.to_numpy(bool), "RP4c": t.ctl_deep.to_numpy(bool),
             "RP4d": t.ctl_don.to_numpy(bool), "RP4e": t.ctl_alu.to_numpy(bool)}
    sel = H.elig & (up | np.logical_or.reduce(list(ctrls.values())))
    f = S.run.path("count_model.npz", local=True)
    z = dict(np.load(f)) if f.exists() and not S.fresh else {}
    if "conv" not in z:                                  # older cache without the per-fit convergence flags
        cmr = count_model(S, H, sel)
        z = {k: np.asarray(cmr[k]) for k in ("beta", "se", "conv", "converged", "ga_imputed", "n_tumours")}
        np.savez(f, **z)
    B, SE, conv = z["beta"], z["se"], z["conv"]
    fitted = conv >= 0
    bad = separated(B, SE)                               # A-23: separated fits dropped, column by column
    B, SE = np.where(bad, np.nan, B), np.where(bad, np.nan, SE)
    idx = np.flatnonzero(sel)
    pos = pd.Series(np.arange(len(idx)), index=idx)

    def sub(m):
        return pos.reindex(np.flatnonzero(m & sel)).to_numpy()
    st = ivw_contrast(B[sub(up)], SE[sub(up)], np.ones(len(sub(up)), bool), None, None)
    nc_perm = (conv[:, 1:] == 0).sum() / max(fitted[:, 1:].sum(), 1)
    S.add_sens("count model (quasi-binomial MMRd coefficient, logit)", "RP2a", float(st[0]), perm_p(st[0], st[1:]),
               n=f"{len(sub(up))} exons; 200 label permutations",
               note=(f"{int(z['converged'])} of {int(fitted[:, 0].sum())} observed fits converged; non-converged share over "
                     f"the permutation columns {nc_perm:.4f}; separated fits (|coefficient| or SE > {SEP_LIMIT:g}) dropped "
                     f"(A-23): {int(bad[:, 0].sum())} of {len(idx)} exons in the observed column ({int(bad[sub(up), 0].sum())} R_up), "
                     f"mean {bad[:, 1:].sum(0).mean():.1f} per permutation column; global axis imputed at 0.5 for "
                     f"{int(z['ga_imputed'])} of {int(z['n_tumours'])} held-out tumours"))
    for i, c in ctrls.items():
        fm, cm_ = up & H.elig, c & H.elig
        lab_, kept = coarsen(H.band, H.terc, fm, cm_)
        m = kept & (fm | cm_) & sel
        rows = sub(m)
        fl, lb = fm[np.flatnonzero(m & sel)], lab_[np.flatnonzero(m & sel)]
        stat = ivw_contrast(B[rows], SE[rows], fl, ~fl, lb)
        null = ivw_class_perm(B[rows, 0], SE[rows, 0], fl, ~fl, lb)
        pc, pl = perm_p(stat[0], null), perm_p(stat[0], stat[1:])
        S.add_sens("count model (quasi-binomial MMRd coefficient, logit)", i, float(stat[0]), max(pc, pl),
                   n=f"{int(fl.sum())} vs {int((~fl).sum())} exons", p_class=pc, p_label=pl)


# ---------------------------------------------------------------- report

def blank_small(desc: pd.DataFrame, min_mmrd: int = 5) -> pd.DataFrame:
    """A-24: in the committed descriptive table, value, p and interval are blanked where fewer than 5
    MMRd tumours (with an arm) enter; the row stays with its n and a note. Raw values stay local."""
    d = desc.copy()
    if "n_mmrd" not in d:
        return d
    m = pd.to_numeric(d.n_mmrd, errors="coerce").lt(min_mmrd).to_numpy()
    for c in ["value", "p", "ci_low", "ci_high"]:
        if c in d:
            d[c] = d[c].astype(float)
            d.loc[m, c] = np.nan
    note = d.note.fillna("").astype(str) if "note" in d else pd.Series("", index=d.index)
    tag = f"blanked (< {min_mmrd} MMRd tumours)"
    d["note"] = [(f"{x}; {tag}" if x else tag) if b else x for x, b in zip(note, m)]
    return d


def report(run: Run) -> None:
    """tests.csv (q within family, verdicts), sensitivity.csv, descriptive.csv and verdict.json - the
    claim-level verdict computed by claim_verdict from the pre-registered rule."""
    t = pd.read_csv(run.path("raw_tests.csv", local=True))
    t["testable"] = t.testable.astype(bool)
    # RP3c takes RP3e's place in RP if RP3e is not testable (and leaves RPS)
    if not t.set_index("id").testable.get("RP3e", True):
        t.loc[t.id.eq("RP3c"), "family"] = "RP"
    t["q"] = np.nan
    for fam, g in t.groupby("family"):
        m = g.testable & g.p.notna()
        t.loc[g.index[m], "q"] = bh(g.p[m].to_numpy())
    t["verdict"] = [test_verdict(e, q, s, ok) for e, q, s, ok in zip(t.estimate, t.q, t.predicted_sign, t.testable)]
    order = ["RP", "RPS", "RPP", "RPC", "RPD", "RPM"]
    t = t.sort_values(["family", "id"], key=lambda c: c.map({f: i for i, f in enumerate(order)}) if c.name == "family" else c)
    cols = ["id", "family", "claim", "n", "estimate", "ci_low", "ci_high", "statistic", "value", "p_class", "p_label", "p", "q",
            "predicted_sign", "verdict", "inference", "note"]
    t[cols].to_csv(run.path("tests.csv"), index=False, float_format="%.6g")
    pd.read_csv(run.path("raw_sensitivity.csv", local=True)).to_csv(run.path("sensitivity.csv"), index=False, float_format="%.6g")
    blank_small(pd.read_csv(run.path("raw_descriptive.csv", local=True))).to_csv(run.path("descriptive.csv"), index=False,
                                                                               float_format="%.6g")
    inst = json.loads(run.path("instrument.json", local=True).read_text())
    v = claim_verdict(t, inst, SESOI)
    v["labels"] = run.labels
    v["families"] = {f: {"tests": int((t.family == f).sum()), "testable": int((t.family.eq(f) & t.testable).sum()),
                         "supported": int((t.family.eq(f) & t.verdict.eq("supported")).sum()),
                         "against": int((t.family.eq(f) & t.verdict.eq("against")).sum())} for f in order}
    run.path("verdict.json").write_text(json.dumps(v, indent=1, default=lambda x: x if not isinstance(x, np.generic) else x.item()))
    pd.set_option("display.width", 250)
    log("tests:\n" + t[["id", "family", "estimate", "p_class", "p_label", "p", "q", "verdict"]].to_string(index=False))
    log(f"claim verdict: {v['verdict']} - {v['reasons']}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", required=True, choices=["annotate", "branchpoints", "counts", "counts_extra", "sizes",
                                                      "strata", "compute", "report"])
    ap.add_argument("--labels", choices=["shuffled", "real"],
                    help="compute/report: 'shuffled' is the dry run (labels permuted within cancer, seed 999, "
                         "outputs under proof/dryrun/); 'real' also needs --real-labels-approved")
    ap.add_argument("--real-labels-approved", action="store_true", help="given only after Gate 2")
    ap.add_argument("--fresh", action="store_true", help="recompute cached builds and engines")
    ap.add_argument("--quick", action="store_true", help="development: leave out the count model")
    a = ap.parse_args()
    if a.only in ("compute", "report"):
        if a.labels is None:
            raise SystemExit("--labels is required for compute and report ('shuffled' = the dry run)")
        run = Run(a.labels, a.real_labels_approved)
        if a.only == "compute":
            compute(run, a.fresh, a.quick)
        else:
            report(run)
        return
    {"annotate": annotate, "branchpoints": branchpoints, "counts": counts, "counts_extra": counts_extra,
     "sizes": sizes_counts, "strata": strata}[a.only]()


if __name__ == "__main__":
    main()
