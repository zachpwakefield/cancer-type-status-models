"""Step 42 (D-52): does genome sequence alone say which exons MSI tumours skip?

Mismatch-repair-deficient (MSI) tumours shorten mononucleotide repeats. A shortened poly-T tract
in the polypyrimidine tract of an acceptor weakens that 3' splice site, and the exon downstream
is skipped (HSP110 exon 9, MRE11 exon 5; genome-wide in colorectal cancer, Jonchere et al. 2024).
This step asks whether that cis route is visible across UCEC, STAD and COAD from the genome
sequence alone, with no model and no label used to choose the exons.

  annotate  label-blind: for every cassette exon of the subtype SE matrix, the longest sense-strand
            homopolymer run of each base overlapping three intronic windows (acceptor -1..-50,
            deep -201..-400 upstream, donor +1..+50 downstream), from the GRCh38 primary assembly
  test      the predictions of config/followup.yaml: msi_cis (M1-M7, K1) and exploratory tables
  tracts    msi_long_tracts (D-53, N1-N3): does the MSI shift of a repeat exon depend on its
            baseline inclusion, as hnRNP C repression of long U-tracts would make it?
  a3ss      msi_a3ss_replication (D-54, Q1-Q2): the same repeats at the two acceptors of A3SS events

Coordinates: a feature_id reads gene#chrom:start-end#strand;flankA;flankB with 0-based half-open
intervals and the two flanking exons in genome order. Positions are counted in transcript
orientation from the exon: -1 is the last intron base before it, +1 the first after it.

Outputs: data.nosync/ctsm_v1/msi/ (per-sample scores, kept local)
         results/subtype/msi_cis/ (aggregates; figure U10: scripts/36_subtype_figures.py)
Usage: python scripts/42_msi_repeats.py [--only annotate|test|tracts|a3ss]
"""
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from _common import REPO, log
from scipy.stats import fisher_exact, mannwhitneyu, rankdata, spearmanr
from scipy.stats import t as t_dist
from statsmodels.stats.multitest import multipletests

from ctsm.matrix_io import load_view

FASTA = "/Users/zacharywakefield/Downloads/GRCh38.primary_assembly.genome.fa.gz"
DATA = REPO / "data.nosync" / "ctsm_v1"
LOCAL = DATA / "msi"
OUT = REPO / "results" / "subtype" / "msi_cis"
WINDOWS = {"acc": (1, 50), "deep": (201, 400), "don": (1, 50)}   # 1-based, inclusive
REACH = 30                     # bases read past each window edge so a run is measured in full
MIN_INTRON = {"acc": 60, "deep": 500, "don": 60}
COMP = str.maketrans("ACGTN", "TGCAN")
CFG = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["msi_cis"]
COHORTS = CFG["cohorts"]
MRE11_E5 = "ENSG00000020922.13#chr11:94479673-94479761#-;94478734-94478876;94485923-94486084"
TOP_K = 50
AS_VIEWS = ["afe", "ale", "hit", "se", "a5ss", "a3ss", "mxe", "ri"]


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


def annotate(features: Path = DATA / "subtype" / "matrices" / "se" / "features.parquet",
             out: Path = LOCAL / "se_repeats.parquet") -> pd.DataFrame:
    f = parse_features(pd.read_parquet(features))
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
    ann.to_parquet(out, index=False)
    log(f"annotated {len(ann):,} exons; canonical AG at -2,-1: {ann['acc_ag'].eq('AG').mean():.3f}")
    return ann



# ============================================================ test (needs the labels)
def cohort_data(cohort: str):
    """Task table with a `group` column, decoded SE matrix (patients x cohort pre-filter features)
    and the features with their repeat annotation."""
    tab = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{cohort}.parquet").reset_index(drop=True)
    grp = np.where(tab["label"].eq("MSI"), "MSI", np.where(tab["label"].eq("POLE"), "POLE", "rest")).astype(object)
    ts = pd.read_csv(DATA / "external" / "gdc_pancanatlas" / "TCGASubtype.20170308.tsv", sep="\t")
    pcr = tab["case_id"].map(ts.set_index("pan.samplesID")["Subtype_other"])
    grp[(grp == "POLE") & pcr.eq("MSI-H").to_numpy()] = "POLE_MSI-H"
    tab["group"] = grp
    vm, feats = load_view(DATA / "subtype" / "matrices", "se", tab["file_id"], cohort)
    ann = pd.read_parquet(LOCAL / "se_repeats.parquet").set_index("feature_id")
    feats = feats[["feature_id", "symbol"]].join(ann.drop(columns="symbol"), on="feature_id")
    return tab, vm.get(), feats


def universe(X: np.ndarray, feats: pd.DataFrame, grp: np.ndarray) -> np.ndarray:
    fin = np.isfinite(X)
    ok = fin.mean(axis=0) >= 0.8
    for g in ("MSI", "rest"):
        ok &= fin[grp == g].sum(axis=0) >= 10
    ok &= (feats["acc_intron_len"].to_numpy() >= MIN_INTRON["acc"]) & feats["acc_ag"].eq("AG").to_numpy()
    return ok


def set_score(P: np.ndarray, in_r: np.ndarray, in_b: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """Per patient: mean percentile over R minus the mean over B, B weighted to R's distribution
    over mean-PSI deciles (missing values left out of both means)."""
    r, b = np.flatnonzero(in_r), np.flatnonzero(in_b)
    share_r = pd.Series(dec[r]).value_counts(normalize=True)
    count_b = pd.Series(dec[b]).value_counts()
    w = pd.Series(dec[b]).map(share_r / count_b).fillna(0.0).to_numpy()
    s_r = np.nanmean(P[:, r], axis=1)
    pb = P[:, b]
    m = np.isfinite(pb)
    s_b = (np.where(m, pb, 0.0) * w).sum(axis=1) / (m * w).sum(axis=1)
    return s_r - s_b


def auroc_lower(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """P(a < b) (ties half) and the two-sided Mann-Whitney p."""
    res = mannwhitneyu(b, a, alternative="two-sided")
    return float(res.statistic / (len(a) * len(b))), float(res.pvalue)


def partial_spearman(x, y, covs) -> tuple[float, float, int]:
    ok = np.isfinite(x) & np.isfinite(y)
    for c in covs:
        ok &= np.isfinite(c)
    rx, ry = rankdata(x[ok]), rankdata(y[ok])
    Z = np.column_stack([np.ones(ok.sum())] + [rankdata(c[ok]) for c in covs])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    r = float(np.corrcoef(ex, ey)[0, 1])
    df = ok.sum() - 2 - len(covs)
    tt = r * np.sqrt(df / max(1e-12, 1 - r * r))
    return r, float(2 * t_dist.sf(abs(tt), df)), int(ok.sum())


def boot_auroc_diff(s1, s2, grp, seed: int, reps: int) -> tuple[float, float]:
    """Paired bootstrap (patients resampled within MSI and rest) of AUROC(s1) - AUROC(s2)."""
    rng = np.random.default_rng(seed)
    i_m, i_r = np.flatnonzero(grp == "MSI"), np.flatnonzero(grp == "rest")
    d = np.empty(reps)
    for k in range(reps):
        m, r = rng.choice(i_m, len(i_m)), rng.choice(i_r, len(i_r))
        d[k] = auroc_lower(s1[m], s1[r])[0] - auroc_lower(s2[m], s2[r])[0]
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975)), max(p, 1 / reps)


def verdict(effect: float, sign: int, q: float) -> str:
    if not np.isfinite(q) or q >= 0.05:
        return "not resolved"
    return "supported" if np.sign(effect) == sign else "against"


def analyse_cohort(cohort: str) -> dict:
    tab, X, feats = cohort_data(cohort)
    grp = tab["group"].to_numpy()
    ok = universe(X, feats, grp)
    X, f = X[:, ok], feats[ok].reset_index(drop=True)
    P = pd.DataFrame(X).rank(pct=True, axis=0).to_numpy()
    mean_psi = np.nanmean(X, axis=0)
    dec = pd.qcut(mean_psi, 10, labels=False, duplicates="drop")
    accT, deepT = f["acc_T"].to_numpy(), f["deep_T"].to_numpy()
    R, B = accT >= 10, accT <= 5
    R_deep = (accT <= 6) & (deepT >= 10)
    B_deep = (accT <= 5) & (deepT <= 5)
    scores = {"acc": set_score(P, R, B, dec), "deep": set_score(P, R_deep, B_deep, dec)}
    for lo, hi in [(6, 7), (8, 9), (10, 11), (12, 13), (14, 999)]:
        scores[f"bin_{lo}"] = set_score(P, (accT >= lo) & (accT <= hi), B, dec)
    scores["A"] = set_score(P, (f["acc_A"].to_numpy() >= 10) & (accT <= 6), B, dec)
    scores["C"] = set_score(P, (f["acc_C"].to_numpy() >= 8) & (accT <= 6), B, dec)
    scores["G"] = set_score(P, (f["acc_G"].to_numpy() >= 8) & (accT <= 6), B, dec)
    donT = f["don_T"].to_numpy()
    scores["donor"] = set_score(P, (donT >= 10) & (accT <= 6), B & (donT <= 5), dec)
    msi, rest = grp == "MSI", grp == "rest"
    # per-exon shift and test
    delta = np.nanmean(P[msi], axis=0) - np.nanmean(P[rest], axis=0)
    pex = np.array([mannwhitneyu(a[np.isfinite(a)], b[np.isfinite(b)]).pvalue
                    for a, b in zip(X[msi].T, X[rest].T)])
    qex = multipletests(pex, method="fdr_bh")[1]
    bmed = pd.Series(delta[B]).groupby(dec[B]).median()
    delta_adj = delta - pd.Series(dec).map(bmed).to_numpy()
    ex = f[["feature_id", "symbol", "chrom", "strand", "es", "ee", "acc_T", "acc_A", "acc_pyr_frac", "deep_T",
            "don_T"]].copy()
    ex["exon_len"] = ex["ee"] - ex["es"]
    ex["mean_psi"], ex["psi_decile"] = mean_psi, dec
    ex["delta"], ex["delta_adj"], ex["p"], ex["q"] = delta, delta_adj, pex, qex
    ex["cohort"] = cohort
    log(f"{cohort}: {ok.sum():,} exons, R {R.sum():,}, B {B.sum():,}, R_deep {R_deep.sum():,}; "
        f"MSI {msi.sum()}, rest {rest.sum()}, POLE {(grp == 'POLE').sum()}")
    sc = tab[["file_id", "case_id", "label", "group"]].copy()
    for k, v in scores.items():
        sc[f"score_{k}"] = v
    sc["cohort"] = cohort
    return {"exons": ex, "scores": sc, "sizes": {"exons": int(ok.sum()), "R": int(R.sum()), "B": int(B.sum()),
            "R_deep": int(R_deep.sum()), "B_deep": int(B_deep.sum())}}


def model_top(cohort: str, k: int) -> pd.Series:
    t = pd.read_parquet(DATA / "subtype" / "bio_cache" / "importance" / f"se__{cohort}__class1.parquet")
    return t.sort_values("rank")["feature_id"].head(k)


def oof_auroc(cohort: str) -> dict:
    run = DATA / "subtype" / "runs" / "subtype"
    def load(v):
        return pd.read_parquet(run / v / cohort / "oof_predictions.parquet").set_index("file_id")
    g = load("gex")
    y = g["label"].eq("MSI").to_numpy()
    as_p = np.mean([load(v).loc[g.index, "p_MSI"].to_numpy() for v in AS_VIEWS], axis=0)
    out = {}
    for name, p in [("GEX", g["p_MSI"].to_numpy()), ("AS ensemble", as_p), ("GEX + AS", 0.5 * g["p_MSI"].to_numpy() + 0.5 * as_p)]:
        out[name] = auroc_lower(p[~y], p[y])[0]     # P(p_MSI higher in MSI)
    return out


def test() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reps, seed = int(CFG["bootstrap"]["reps"]), int(CFG["bootstrap"]["seed"])
    res = {c: analyse_cohort(c) for c in COHORTS}
    ex = pd.concat([r["exons"] for r in res.values()], ignore_index=True)
    sc = pd.concat([r["scores"] for r in res.values()], ignore_index=True)
    sc.to_parquet(LOCAL / "sample_scores.parquet", index=False)
    ex.to_parquet(LOCAL / "exon_effects.parquet", index=False)
    pred = []

    def add(pid, cohort, effect, sign, p, n, what, lo=np.nan, hi=np.nan, family="M"):
        pred.append({"id": pid, "family": family, "cohort": cohort, "what": what, "effect": effect,
                     "predicted_sign": sign, "ci_low": lo, "ci_high": hi, "p": p, "n": n})

    for c in COHORTS:
        s, e = res[c]["scores"], res[c]["exons"]
        g = s["group"].to_numpy()
        a, p = auroc_lower(s.loc[g == "MSI", "score_acc"].to_numpy(), s.loc[g == "rest", "score_acc"].to_numpy())
        add("M1", c, a - 0.5, 1, p, f"MSI {int((g == 'MSI').sum())}, rest {int((g == 'rest').sum())}",
            "AUROC(score lower in MSI) - 0.5")
        sub = e[e["acc_T"] >= 6]
        r, p, n = partial_spearman(sub["acc_T"].to_numpy(), sub["delta"].to_numpy(), [sub["mean_psi"].to_numpy()])
        add("M2", c, r, -1, p, n, "partial Spearman(acceptor T-run, exon delta | mean PSI)")
        r, p, n = partial_spearman(sub["acc_T"].to_numpy(), sub["delta"].to_numpy(),
                                   [sub["mean_psi"].to_numpy(), sub["acc_pyr_frac"].to_numpy()])
        add("M3", c, r, -1, p, n, "as M2, also | pyrimidine fraction -1..-50")
        s1, s2 = s["score_acc"].to_numpy(), s["score_deep"].to_numpy()
        a1 = auroc_lower(s1[g == "MSI"], s1[g == "rest"])[0]
        a2 = auroc_lower(s2[g == "MSI"], s2[g == "rest"])[0]
        lo, hi, p = boot_auroc_diff(s1, s2, g, seed, reps)
        add("M4", c, a1 - a2, 1, p, f"R_deep {res[c]['sizes']['R_deep']}", "AUROC acceptor - AUROC deep", lo, hi)
    s = res["UCEC_TCGA"]["scores"]
    a, p = auroc_lower(s.loc[s.group.eq("MSI"), "score_acc"].to_numpy(), s.loc[s.group.eq("POLE"), "score_acc"].to_numpy())
    add("M5", "UCEC_TCGA", a - 0.5, 1, p, f"MSI {int(s.group.eq('MSI').sum())}, POLE {int(s.group.eq('POLE').sum())}",
        "AUROC(score lower in MSI than POLE) - 0.5")
    ann = pd.read_parquet(LOCAL / "se_repeats.parquet").set_index("feature_id")
    pf = pd.read_parquet(DATA / "subtype" / "matrices" / "se" / "prefilter.parquet")
    feats = pd.read_parquet(DATA / "subtype" / "matrices" / "se" / "features.parquet")["feature_id"]
    m6 = []
    for c in COHORTS:
        u = ann.loc[feats[pf[c].to_numpy()]]
        u = u[(u["acc_intron_len"] >= MIN_INTRON["acc"]) & u["acc_ag"].eq("AG")]
        top = model_top(c, TOP_K)
        in_top = u.index.isin(top)
        rr = (u["acc_T"] >= 10).to_numpy()
        tbl = [[int((in_top & rr).sum()), int((in_top & ~rr).sum())], [int((~in_top & rr).sum()), int((~in_top & ~rr).sum())]]
        orr, p = fisher_exact(tbl)
        add("M6", c, float(np.log(orr)) if orr > 0 else -np.inf, 1, p, f"top {in_top.sum()}: {tbl[0][0]} with T>=10",
            "log odds ratio (top-50 MSI events vs other exons, acceptor T-run >= 10)")
        m6.append({"cohort": c, "top_n": int(in_top.sum()), "top_R": tbl[0][0], "other_n": int((~in_top).sum()),
                   "other_R": tbl[1][0], "odds_ratio": orr, "p": p})
    pairs = [("UCEC_TCGA", "STAD_TCGA"), ("UCEC_TCGA", "COAD_TCGA"), ("STAD_TCGA", "COAD_TCGA")]
    m7 = []
    for a_, b_ in pairs:
        for setname, cond in [("R", lambda d: d["acc_T"] >= 10), ("B", lambda d: d["acc_T"] <= 5)]:
            da = ex[ex.cohort.eq(a_) & cond(ex)].set_index("feature_id")["delta_adj"]
            db = ex[ex.cohort.eq(b_) & cond(ex)].set_index("feature_id")["delta_adj"]
            j = da.index.intersection(db.index)
            r = spearmanr(da.loc[j], db.loc[j])
            m7.append({"pair": f"{a_} ~ {b_}", "set": setname, "n": len(j), "rho": r.statistic, "p": r.pvalue})
            if setname == "R":
                add("M7", f"{a_} ~ {b_}", float(r.statistic), 1, float(r.pvalue), len(j),
                    "Spearman of background-adjusted exon delta between cohorts (R exons)")
    # K1: MRE11 exon 5 from the full matrix, whatever the pre-filter
    for c in COHORTS:
        tab = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{c}.parquet")
        g = np.where(tab["label"].eq("MSI"), "MSI", np.where(tab["label"].eq("POLE"), "POLE", "rest"))
        vm, fe = load_view(DATA / "subtype" / "matrices", "se", tab["file_id"])
        j = int(np.flatnonzero(fe["feature_id"].eq(MRE11_E5))[0])
        v = vm.get(cols=[j])[:, 0]
        a = v[(g == "MSI") & np.isfinite(v)]
        b = v[(g == "rest") & np.isfinite(v)]
        au, p = auroc_lower(a, b)
        add("K1", c, au - 0.5, 1, p, f"MSI {len(a)}, rest {len(b)}; median {np.median(a):.3f} vs {np.median(b):.3f}",
            "AUROC(MRE11 exon 5 PSI lower in MSI) - 0.5", family="K")
    pr = pd.DataFrame(pred)
    pr["q"] = np.nan
    for fam, idx in pr.groupby("family").groups.items():
        pr.loc[idx, "q"] = multipletests(pr.loc[idx, "p"], method="fdr_bh")[1]
    pr["verdict"] = [verdict(e, sg, q) for e, sg, q in zip(pr["effect"], pr["predicted_sign"], pr["q"])]
    pr.to_csv(OUT / "predictions.csv", index=False)
    pd.DataFrame(m6).to_csv(OUT / "model_top_events.csv", index=False)
    pd.DataFrame(m7).to_csv(OUT / "cross_cohort_concordance.csv", index=False)
    log("predictions:\n" + pr[["id", "cohort", "effect", "p", "q", "verdict"]].to_string(index=False))
    exploratory(res, ex, sc)


def exploratory(res: dict, ex: pd.DataFrame, sc: pd.DataFrame) -> None:
    rows = []
    for c in COHORTS:
        s = sc[sc.cohort.eq(c)]
        g = s["group"].to_numpy()
        for k in [col for col in s.columns if col.startswith("score_")]:
            v = s[k].to_numpy()
            a, p = auroc_lower(v[g == "MSI"], v[g == "rest"])
            row = {"cohort": c, "score": k[6:], "auroc_msi_vs_rest": a, "p": p}
            if c == "UCEC_TCGA":
                row["auroc_pole_vs_rest"] = auroc_lower(v[g == "POLE"], v[g == "rest"])[0]
                row["auroc_msi_vs_cn_low"] = auroc_lower(v[g == "MSI"], v[s["label"].eq("CN-low").to_numpy()])[0]
            rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "score_aurocs.csv", index=False)
    # group summaries (aggregates only)
    summ = (sc.groupby(["cohort", "group"])["score_acc"]
            .agg(n="size", median="median", q1=lambda x: x.quantile(0.25), q3=lambda x: x.quantile(0.75))
            .reset_index())
    summ.to_csv(OUT / "score_by_group.csv", index=False)
    # detector: MSI vs every other patient, next to the models' out-of-fold AUROCs
    det = []
    for c in COHORTS:
        s = sc[sc.cohort.eq(c)]
        y = s["group"].eq("MSI").to_numpy()
        det.append({"cohort": c, "method": "repeat score (no training)",
                    "auroc": auroc_lower(s["score_acc"].to_numpy()[y], s["score_acc"].to_numpy()[~y])[0]})
        for k, v in oof_auroc(c).items():
            det.append({"cohort": c, "method": f"{k} model (out-of-fold)", "auroc": v})
    pd.DataFrame(det).to_csv(OUT / "detector.csv", index=False)
    # share of MSI-shifted exons that carry an acceptor repeat
    sh = []
    for c in COHORTS:
        e = ex[ex.cohort.eq(c)]
        rr = (e["acc_T"] >= 10).to_numpy()
        for d, cond in [("less included in MSI", (e["q"] < 0.05) & (e["delta"] < 0)),
                        ("more included in MSI", (e["q"] < 0.05) & (e["delta"] > 0))]:
            cond = cond.to_numpy()
            orr, p = fisher_exact([[int((cond & rr).sum()), int((cond & ~rr).sum())],
                                   [int((~cond & rr).sum()), int((~cond & ~rr).sum())]])
            sh.append({"cohort": c, "direction": d, "n_exons": int(cond.sum()), "n_with_T10": int((cond & rr).sum()),
                       "share_with_T10": (cond & rr).sum() / max(1, cond.sum()), "share_all": rr.mean(),
                       "odds_ratio": orr, "p": p})
        skipped_r = (e["q"] < 0.05).to_numpy() & (e["delta"] < 0).to_numpy() & rr
        sh[-2]["frame_shift_share_skipped_R"] = float((e.loc[skipped_r, "exon_len"] % 3 != 0).mean())
        sh[-2]["frame_shift_share_all_R"] = float((e.loc[rr, "exon_len"] % 3 != 0).mean())
    pd.DataFrame(sh).to_csv(OUT / "share_explained.csv", index=False)
    # repeat exons per cohort (aggregates per exon) and the shared list
    keep = ["feature_id", "symbol", "chrom", "strand", "es", "ee", "exon_len", "acc_T", "acc_pyr_frac", "mean_psi",
            "delta", "delta_adj", "p", "q", "cohort"]
    rex = ex[ex["acc_T"] >= 10][keep]
    rex.to_csv(OUT / "repeat_exon_effects.csv", index=False, float_format="%.4g")
    w = rex.pivot_table(index=["feature_id", "symbol", "acc_T", "exon_len"], columns="cohort",
                        values=["delta_adj", "q"]).dropna()
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w["mean_delta_adj"] = w[[f"delta_adj_{c}" for c in COHORTS]].mean(axis=1)
    w["skipped_all_three"] = np.all([(w[f"delta_adj_{c}"] < 0) & (w[f"q_{c}"] < 0.05) for c in COHORTS], axis=0)
    w.sort_values("mean_delta_adj").reset_index().to_csv(OUT / "shared_repeat_exons.csv", index=False, float_format="%.4g")
    # dose-response summary of exon deltas by bin (aggregates)
    bins = pd.cut(ex["acc_T"], [-1, 5, 7, 9, 11, 13, 999], labels=["<=5", "6-7", "8-9", "10-11", "12-13", ">=14"])
    dose = (ex.assign(bin=bins).groupby(["cohort", "bin"], observed=True)["delta_adj"]
            .agg(n="size", median="median", q1=lambda x: x.quantile(0.25), q3=lambda x: x.quantile(0.75),
                 share_negative=lambda x: (x < 0).mean())
            .reset_index())
    dose.to_csv(OUT / "dose_exon_delta.csv", index=False)
    # the ten U7 events (SE ones carry an annotation)
    u7 = pd.read_csv(REPO / "results" / "subtype" / "interpret" / "shared_class_events.csv")
    u7 = u7[u7["class"].eq("MSI")]
    ann = pd.read_parquet(LOCAL / "se_repeats.parquet").set_index("feature_id")
    u7se = u7[u7["view"].eq("SE")].drop_duplicates("feature_id")[["feature_id", "symbol"]]
    u7se = u7se.join(ann[["acc_T", "acc_A", "deep_T", "don_T"]], on="feature_id")
    u7se.to_csv(OUT / "u7_se_events_repeats.csv", index=False)
    s = sc[sc.cohort.eq("UCEC_TCGA")]
    ph = s[s.group.eq("POLE_MSI-H")]
    if len(ph):
        v = ph["score_acc"].iloc[0]
        log(f"POLE MSI-H tumour: score {v:.4f}; percentile among MSI "
            f"{(s.loc[s.group.eq('MSI'), 'score_acc'] < v).mean():.2f}, among rest {(s.loc[s.group.eq('rest'), 'score_acc'] < v).mean():.2f}")
    sizes = pd.DataFrame({c: r["sizes"] for c, r in res.items()}).T
    sizes.to_csv(OUT / "set_sizes.csv")
    log("written " + str(OUT))


# ============================================================ tracts (D-53)
def fisher_z(r1: float, n1: int, r2: float, n2: int) -> tuple[float, float]:
    from scipy.stats import norm
    z = (np.arctanh(r1) - np.arctanh(r2)) / np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    return float(z), float(2 * norm.sf(abs(z)))


def tracts() -> None:
    ex = pd.read_parquet(LOCAL / "exon_effects.parquet")
    u7 = pd.read_csv(REPO / "results" / "subtype" / "interpret" / "shared_class_events.csv")
    drop = set(u7.loc[u7["class"].eq("MSI") & u7["view"].eq("SE"), "feature_id"])
    ex = ex[~ex["feature_id"].isin(drop)]
    log(f"left out {len(drop)} MSI-class SE events of SUBTYPE.md section 8")
    rows, expl = [], []
    for c in COHORTS:
        e = ex[ex.cohort.eq(c)]
        accT, deepT = e["acc_T"], e["deep_T"]
        sets = {"R": e[accT >= 10], "R_deep": e[(accT <= 6) & (deepT >= 10)],
                "T10-14": e[(accT >= 10) & (accT <= 14)], "T>=20": e[accT >= 20]}
        rho = {k: spearmanr(v["delta_adj"], v["mean_psi"]) for k, v in sets.items()}
        n = {k: len(v) for k, v in sets.items()}
        for k in sets:
            expl.append({"cohort": c, "set": k, "n": n[k], "rho_delta_adj_vs_mean_psi": rho[k].statistic, "p": rho[k].pvalue})
        rows.append({"id": "N1", "cohort": c, "effect": rho["R"].statistic, "predicted_sign": -1, "p": rho["R"].pvalue,
                     "n": f"R {n['R']}", "what": "Spearman(delta_adj, mean PSI), R exons"})
        z, p = fisher_z(rho["R"].statistic, n["R"], rho["R_deep"].statistic, n["R_deep"])
        rows.append({"id": "N2", "cohort": c, "effect": rho["R"].statistic - rho["R_deep"].statistic, "predicted_sign": -1,
                     "p": p, "n": f"R {n['R']}, R_deep {n['R_deep']}", "what": "rho(R) - rho(R_deep), Fisher z"})
        z, p = fisher_z(rho["T>=20"].statistic, n["T>=20"], rho["T10-14"].statistic, n["T10-14"])
        rows.append({"id": "N3", "cohort": c, "effect": rho["T>=20"].statistic - rho["T10-14"].statistic, "predicted_sign": -1,
                     "p": p, "n": f"T>=20 {n['T>=20']}, T10-14 {n['T10-14']}", "what": "rho(T>=20) - rho(T10-14), Fisher z"})
    pr = pd.DataFrame(rows)
    pr["q"] = multipletests(pr["p"], method="fdr_bh")[1]
    pr["verdict"] = [verdict(e, sg, q) for e, sg, q in zip(pr["effect"], pr["predicted_sign"], pr["q"])]
    pr.to_csv(OUT / "long_tract_predictions.csv", index=False)
    pd.DataFrame(expl).to_csv(OUT / "long_tract_correlations.csv", index=False)
    # exploratory: share up by baseline tertile x run length; |delta_adj| by run length
    tab = []
    for c in COHORTS:
        e = ex[ex.cohort.eq(c) & (ex["acc_T"] >= 10)].copy()
        e["baseline"] = pd.qcut(e["mean_psi"], 3, labels=["low", "middle", "high"])
        e["run"] = pd.cut(e["acc_T"], [9, 14, 19, 999], labels=["10-14", "15-19", ">=20"])
        g = e.groupby(["run", "baseline"], observed=True)["delta_adj"]
        t = g.agg(n="size", share_up=lambda x: (x > 0).mean(), median="median").reset_index()
        t["cohort"] = c
        tab.append(t)
    pd.concat(tab).to_csv(OUT / "long_tract_share_up.csv", index=False)
    ab = []
    for c in COHORTS:
        e = ex[ex.cohort.eq(c)]
        for win, col, cond in [("acceptor", "acc_T", np.ones(len(e), bool)), ("deep", "deep_T", (e["acc_T"] <= 6).to_numpy())]:
            sub = e[cond & np.isfinite(e[col]).to_numpy()]
            r = spearmanr(sub[col], sub["delta_adj"].abs())
            ab.append({"cohort": c, "window": win, "n": len(sub), "rho_abs_delta_adj_vs_run": r.statistic, "p": r.pvalue})
            for lo, hi in [(0, 5), (6, 9), (10, 14), (15, 19), (20, 999)]:
                b = sub[(sub[col] >= lo) & (sub[col] <= hi)]["delta_adj"].abs()
                ab.append({"cohort": c, "window": win, "bin": f"{lo}-{hi}" if hi < 999 else f">={lo}", "n": len(b),
                           "median_abs_delta_adj": b.median()})
    pd.DataFrame(ab).to_csv(OUT / "abs_shift_by_run.csv", index=False)
    log("long tracts:\n" + pr[["id", "cohort", "effect", "p", "q", "verdict", "n"]].to_string(index=False))
    log("retained-intron check:\n" + ri_check().round(3).to_string(index=False))


def ri_check() -> pd.DataFrame:
    """Exploratory: for repeat exons whose upstream intron is also an RI event of the subtype matrix,
    the RI event's MSI shift beside the exon's - would retention explain an apparent inclusion?"""
    ex = pd.read_parquet(LOCAL / "exon_effects.parquet")
    ann = pd.read_parquet(LOCAL / "se_repeats.parquet").set_index("feature_id")[["left_end", "right_start"]]
    ex = ex[ex["acc_T"] >= 10].join(ann, on="feature_id")
    ex["intron"] = np.where(ex["strand"].eq("+"),
                            ex["chrom"] + ":" + ex["left_end"].astype(str) + "-" + ex["es"].astype(str),
                            ex["chrom"] + ":" + ex["ee"].astype(str) + "-" + ex["right_start"].astype(str))
    ri = pd.read_parquet(DATA / "subtype" / "matrices" / "ri" / "features.parquet")
    q = ri["feature_id"].str.split("#", expand=True)
    fl = q[2].str.split(";", expand=True)
    ri["intron"] = (q[1].str.split(":", expand=True)[0] + ":" + fl[1].str.split("-", expand=True)[1] + "-"
                    + fl[2].str.split("-", expand=True)[0])
    out = []
    for c in COHORTS:
        tab = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{c}.parquet")
        msi, rest = tab["label"].eq("MSI").to_numpy(), ~tab["label"].isin(["MSI", "POLE"]).to_numpy()
        e = ex[ex["cohort"].eq(c)]
        hit = ri[ri["intron"].isin(e["intron"])].drop_duplicates("intron")
        vm, fe = load_view(DATA / "subtype" / "matrices", "ri", tab["file_id"])
        cols = np.flatnonzero(fe["feature_id"].isin(hit["feature_id"]).to_numpy())
        X = vm.get(cols=cols)
        P = pd.DataFrame(X).rank(pct=True, axis=0).to_numpy()
        m = pd.DataFrame({"intron": hit.set_index("feature_id").loc[fe["feature_id"].iloc[cols], "intron"].to_numpy(),
                          "ri_delta": np.nanmean(P[msi], axis=0) - np.nanmean(P[rest], axis=0),
                          "ri_finite": np.isfinite(X).mean(axis=0)})
        j = e.merge(m, on="intron")
        j = j[j["ri_finite"] >= 0.5]
        for lab, sub in [("10-19", j[j["acc_T"] < 20]), (">=20", j[j["acc_T"] >= 20])]:
            out.append({"cohort": c, "run": lab, "n": len(sub),
                        "rho_exon_vs_ri": sub[["delta_adj", "ri_delta"]].corr("spearman").iloc[0, 1] if len(sub) > 2 else np.nan,
                        "median_ri_delta": sub["ri_delta"].median(), "median_exon_delta_adj": sub["delta_adj"].median()})
    r = pd.DataFrame(out)
    r.to_csv(OUT / "ri_check.csv", index=False)
    return r


# ============================================================ a3ss (D-54)
def annotate_a3ss() -> pd.DataFrame:
    """Label-blind: T-run and AG at the distal (long-form) and proximal (short-form) acceptor of
    every A3SS event. feature_id = gene#chrom:long#strand;short;upstream flank (0-based)."""
    f = pd.read_parquet(DATA / "subtype" / "matrices" / "a3ss" / "features.parquet")
    q = f["feature_id"].str.split("#", expand=True)
    ch = q[1].str.split(":", expand=True)
    lg = ch[1].str.split("-", expand=True).astype(int)
    fl = q[2].str.split(";", expand=True)
    sh = fl[1].str.split("-", expand=True).astype(int)
    fk = fl[2].str.split("-", expand=True).astype(int)
    t = pd.DataFrame({"feature_id": f["feature_id"], "symbol": f["symbol"], "chrom": ch[0], "strand": fl[0],
                      "ls": lg[0], "le": lg[1], "ss": sh[0], "se": sh[1], "fs": fk[0], "fe": fk[1]})
    far = WINDOWS["acc"][1] + REACH
    rows = []
    for chrom, seq in read_genome(set(t["chrom"])):
        for r in t[t["chrom"].eq(chrom)].itertuples(index=False):
            rec = {"feature_id": r.feature_id}
            if r.strand == "+":
                acc = {"distal": r.ls, "proximal": r.ss}
                intron = r.ls - r.fe
                get = lambda pos: seq[max(0, pos - far):pos][::-1]
            else:
                acc = {"distal": r.le, "proximal": r.se}
                intron = r.fs - r.le
                get = lambda pos: seq[pos:pos + far].translate(COMP)
            rec["intron_len"], rec["segment_len"] = intron, abs(acc["proximal"] - acc["distal"])
            for k, pos in acc.items():
                w = get(pos)
                rec[f"{k}_T"] = max_runs(w, 1, WINDOWS["acc"][1])["T"]
                rec[f"{k}_ag"] = w[:2][::-1]
            rows.append(rec)
        log(f"a3ss {chrom}")
    ann = t.merge(pd.DataFrame(rows), on="feature_id")
    ann.to_parquet(LOCAL / "a3ss_repeats.parquet", index=False)
    log(f"a3ss annotated {len(ann):,}; AG distal {ann['distal_ag'].eq('AG').mean():.3f}, "
        f"proximal {ann['proximal_ag'].eq('AG').mean():.3f}")
    return ann


def a3ss() -> None:
    from scipy.stats import wilcoxon
    ann = annotate_a3ss().set_index("feature_id")
    rows, keep = [], []
    for c in COHORTS:
        tab = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{c}.parquet")
        grp = np.where(tab["label"].eq("MSI"), "MSI", np.where(tab["label"].eq("POLE"), "POLE", "rest"))
        vm, feats = load_view(DATA / "subtype" / "matrices", "a3ss", tab["file_id"], c)
        X = vm.get()
        a = ann.loc[feats["feature_id"]].reset_index()
        fin = np.isfinite(X)
        ok = (fin.mean(axis=0) >= 0.8) & (fin[grp == "MSI"].sum(axis=0) >= 10) & (fin[grp == "rest"].sum(axis=0) >= 10)
        ok &= a["distal_ag"].eq("AG").to_numpy() & a["proximal_ag"].eq("AG").to_numpy()
        ok &= (a["intron_len"] >= MIN_INTRON["acc"]).to_numpy()
        ok &= ~a.duplicated(["chrom", "ls", "le", "ss", "se"]).to_numpy()
        X, a = X[:, ok], a[ok].reset_index(drop=True)
        P = pd.DataFrame(X).rank(pct=True, axis=0).to_numpy()
        mean_psi = np.nanmean(X, axis=0)
        dec = pd.qcut(mean_psi, 10, labels=False, duplicates="drop")
        delta = np.nanmean(P[grp == "MSI"], axis=0) - np.nanmean(P[grp == "rest"], axis=0)
        dT, pT = a["distal_T"].to_numpy(), a["proximal_T"].to_numpy()
        B = (dT <= 5) & (pT <= 5)
        bmed = pd.Series(delta[B]).groupby(dec[B]).median()
        adj = delta - pd.Series(dec).map(bmed).to_numpy()
        run = np.where(pT <= 6, dT, np.where(dT <= 6, pT, np.nan))
        side = np.where((pT <= 6) & (dT > 6), 1, np.where((dT <= 6) & (pT > 6), -1, 0))
        oriented = adj * side
        e = a[["feature_id", "symbol", "distal_T", "proximal_T", "segment_len"]].copy()
        e["mean_psi"], e["delta"], e["delta_adj"], e["run"], e["side"], e["oriented"] = mean_psi, delta, adj, run, side, oriented
        e["cohort"] = c
        keep.append(e)
        log(f"{c}: {ok.sum():,} A3SS events, B {B.sum():,}")
        for pid, lo, hi, sign in [("Q1", 10, 19, -1), ("Q2", 20, 999, 1)]:
            m = (run >= lo) & (run <= hi) & (side != 0)
            v = oriented[m]
            w = wilcoxon(v)
            rows.append({"id": pid, "cohort": c, "effect": float(np.median(v)), "predicted_sign": sign,
                         "p": float(w.pvalue), "n": int(m.sum()), "n_distal": int((m & (side == 1)).sum()),
                         "share_used_less": float((v < 0).mean()),
                         "what": f"median oriented delta_adj, runs {lo}-{hi if hi < 999 else ''}"})
    pr = pd.DataFrame(rows)
    pr["q"] = multipletests(pr["p"], method="fdr_bh")[1]
    pr["verdict"] = [verdict(e, sg, q) for e, sg, q in zip(pr["effect"], pr["predicted_sign"], pr["q"])]
    pr.to_csv(OUT / "a3ss_predictions.csv", index=False)
    ev = pd.concat(keep, ignore_index=True)
    ev.to_parquet(LOCAL / "a3ss_effects.parquet", index=False)
    ev[ev["side"].ne(0) & (ev["run"] >= 10)].to_csv(OUT / "a3ss_repeat_events.csv", index=False, float_format="%.4g")
    log("a3ss:\n" + pr[["id", "cohort", "effect", "share_used_less", "p", "q", "verdict", "n", "n_distal"]].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["annotate", "test", "tracts", "a3ss"])
    a = ap.parse_args()
    if a.only in (None, "annotate"):
        annotate()
    if a.only in (None, "test"):
        test()
    if a.only in (None, "tracts"):
        tracts()
    if a.only in (None, "a3ss"):
        a3ss()

if __name__ == "__main__":
    main()
