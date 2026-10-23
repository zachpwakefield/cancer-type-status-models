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

The tests themselves (steps 1-8 of D-77 and the mechanism tests) are added after the
pre-registration is reviewed and committed.

Positions count from the exon in transcript orientation: 1 is the last intron base before the exon
(-1 in D-52), larger numbers lie further upstream; after the donor, 1 is the first intron base.

Outputs: data.nosync/ctsm_v1/msi/proof/ (per-exon annotation, per-sample counts; kept local)
         results/subtype/msi_cis/code/proof/ (label-blind aggregates: class sizes, PWM, QC)
Usage: python scripts/62_msi_rising_proof.py --only annotate|branchpoints|counts|sizes|strata
"""
from __future__ import annotations

import argparse
import gzip
import warnings
import importlib.util
import json
import re
from pathlib import Path

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


def _count_cancer(proj: str, samp: pd.DataFrame, cols: np.ndarray, store: str) -> dict:
    store = Path(store)
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
            np.save(PROOF / f"{proj}__{v}__{k}.npy", x)
    samp.to_parquet(PROOF / f"{proj}__samples.parquet", index=False)
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
    span = [int(s_.split(":")[0].split("-")[1]) - int(s_.split(":")[0].split("-")[0]) + 1 for s_ in f.index]
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["annotate", "branchpoints", "counts", "sizes", "strata"], required=True)
    a = ap.parse_args()
    {"annotate": annotate, "branchpoints": branchpoints, "counts": counts, "sizes": sizes_counts,
     "strata": strata}[a.only]()


if __name__ == "__main__":
    main()
