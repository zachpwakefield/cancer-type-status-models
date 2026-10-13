"""Step 49 (D-60): D-58's two HER2s again, with ERBB2 amplification that has no ceiling.

ABSOLUTE never reports more than 7 copies at ERBB2, so D-58's amplified group held only the
highest-level amplifications. Here the copy number comes from the PanCanAtlas SNP6 copy ratio at
ERBB2's anchor (the same GRCh37 position as D-58), made absolute with ABSOLUTE purity a and ploidy t:
    q = (R x (a t + 2 (1 - a)) - 2 (1 - a)) / a,   R = 2^Segment_Mean
and amplified = q / t >= 2, D-58's rule. Probabilities, AUROCs and the knockout are step 47's.

Outputs: results/subtype/brca/her2_uncapped/ (aggregates only); the per-tumour table stays in
data.nosync (D-12).
Usage: python scripts/49_brca_her2_uncapped.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from _common import REPO, log
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

_spec = importlib.util.spec_from_file_location("brca47", Path(__file__).with_name("47_brca_her2.py"))
m47 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m47)
m42, auc = m47.m42, m47.auc

SEG = m47.EXT / "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg"
OUT = REPO / "results" / "subtype" / "brca" / "her2_uncapped"
LOCAL = m47.DATA / "brca"
SEED, BOOT, MIN_NEW = 0, 2000, 10
LUM = ["LumA", "LumB"]


def weighted_median(v: np.ndarray, w: np.ndarray) -> float:
    o = np.argsort(v)
    c = np.cumsum(w[o])
    return float(v[o][np.searchsorted(c, c[-1] / 2)])


def snp6_at_erbb2(samples: set[str]) -> pd.DataFrame:
    """Segment_Mean at ERBB2's anchor per 15-character sample (first aliquot in barcode order), and each
    sample's length-weighted median Segment_Mean over the genome (the convention check)."""
    anchor = m47.erbb2_anchor()
    seg = pd.read_csv(SEG, sep="\t", dtype={"Chromosome": str})
    seg["s15"] = seg.Sample.str[:15]
    seg = seg[seg.s15.isin(samples)]
    first = seg.sort_values("Sample").drop_duplicates("s15").set_index("s15").Sample
    seg = seg[seg.Sample.isin(first)]
    at = seg[seg.Chromosome.eq("17") & (seg.Start <= anchor) & (seg.End >= anchor)].drop_duplicates("s15")
    wmed = seg.assign(len=seg.End - seg.Start).groupby("s15").apply(
        lambda g: weighted_median(g.Segment_Mean.to_numpy(), g["len"].to_numpy()))
    d = pd.DataFrame({"seg_mean": at.set_index("s15").Segment_Mean, "genome_median": wmed})
    log(f"SNP6: {len(first)} cohort samples, {d.seg_mean.notna().sum()} with a segment at chr17:{anchor:,.0f}")
    return d


def boot_auc_diff(y: np.ndarray, a: np.ndarray, b: np.ndarray, rng) -> tuple[float, float, float]:
    """AUROC(a) - AUROC(b) for y; paired bootstrap over tumours: (diff, two-sided p floored at 1/(B+1), CI half)."""
    obs = auc(a[y], a[~y]) - auc(b[y], b[~y])
    bs = []
    for _ in range(BOOT):
        i = rng.integers(0, len(y), len(y))
        yy = y[i]
        if yy.all() or not yy.any():
            continue
        bs.append(auc(a[i][yy], a[i][~yy]) - auc(b[i][yy], b[i][~yy]))
    bs = np.asarray(bs)
    p = min(1.0, 2 * (1 + min((bs <= 0).sum(), (bs >= 0).sum())) / (1 + len(bs)))
    return obs, float(p), float(np.percentile(bs, 97.5) - np.percentile(bs, 2.5)) / 2


def boot_rho_diff(x: np.ndarray, a: np.ndarray, b: np.ndarray, rng) -> tuple[float, float, float]:
    """Spearman(a, x) - Spearman(b, x); paired bootstrap as boot_auc_diff."""
    rho = lambda u, i: spearmanr(u[i], x[i]).statistic
    idx = np.arange(len(x))
    obs = rho(a, idx) - rho(b, idx)
    bs = np.array([rho(a, i) - rho(b, i) for i in (rng.integers(0, len(x), len(x)) for _ in range(BOOT))])
    bs = bs[np.isfinite(bs)]
    p = min(1.0, 2 * (1 + min((bs <= 0).sum(), (bs >= 0).sum())) / (1 + len(bs)))
    return float(obs), float(p), float(np.percentile(bs, 97.5) - np.percentile(bs, 2.5)) / 2


def tumours() -> pd.DataFrame:
    d = m47.oof_her2()
    d["s15"] = d.sample_id.str[:15]
    s = snp6_at_erbb2(set(d.s15))
    ab = m47.erbb2_copy_number().rename(columns={"cn": "cn_absolute", "rel": "rel_absolute",
                                                 "amplified": "amp_absolute"})
    d = d.join(s, on="s15").join(ab, on="s15")
    a, t = d.purity, d.ploidy
    d["R"] = 2.0 ** d.seg_mean
    d["q"] = (d.R * (a * t + 2 * (1 - a)) - 2 * (1 - a)) / a
    d["rel"] = d.q / t
    d["amp"] = np.where(d.rel.notna(), d.rel >= 2, np.nan)
    d["amp_raw"] = np.where(d.R.notna(), d.R >= 2, np.nan)
    d["clin"] = d.case_id.map(m47.clinical_her2())
    LOCAL.mkdir(parents=True, exist_ok=True)
    d.to_parquet(LOCAL / "her2_uncapped.parquet")
    return d


def tests(d: pd.DataFrame, rule: str, rng, verdicts: bool) -> pd.DataFrame:
    """HU1-HU7 under one amplification rule ('amp' = D-60, 'amp_raw' = the exploratory raw rule)."""
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    k = d[rule].notna()
    lum = d[k & d.label.isin(LUM)].assign(y=lambda x: x[rule].astype(bool))
    her = d[k & d.label.eq("Her2")].assign(y=lambda x: x[rule].astype(bool))
    a, n = lum[lum.y], lum[~lum.y]
    add("HU1", auc(a.p_gex, n.p_gex) - 0.5, 1, float(mannwhitneyu(a.p_gex, n.p_gex).pvalue),
        f"amplified {len(a)}, not {len(n)}", "LumA/LumB: AUROC(GEX P(Her2) higher when amplified) - 0.5")
    diff, p, half = boot_auc_diff(lum.y.to_numpy(), lum.p_as.to_numpy(), lum.p_gex.to_numpy(), rng)
    add("HU2", diff, 1, p, f"{len(lum)} tumours; 95% CI half-width {half:.3f}; AS {auc(a.p_as, n.p_as):.3f}, "
        f"GEX {auc(a.p_gex, n.p_gex):.3f}", "LumA/LumB: AUROC(AS ensemble) - AUROC(GEX), amplified against not")
    new = lum[~lum.amp_absolute.astype("boolean").fillna(False).astype(bool)]
    n_new = int(new.y.sum())
    if n_new >= MIN_NEW:
        diff, p, half = boot_auc_diff(new.y.to_numpy(), new.p_as.to_numpy(), new.p_gex.to_numpy(), rng)
        na, nn = new[new.y], new[~new.y]
        add("HU3", diff, 1, p, f"newly amplified {n_new}, not {len(nn)}; 95% CI half-width {half:.3f}; "
            f"AS {auc(na.p_as, nn.p_as):.3f}, GEX {auc(na.p_gex, nn.p_gex):.3f}",
            "LumA/LumB without D-58's amplified: AUROC(AS ensemble) - AUROC(GEX), newly amplified against not")
    else:
        pred.append({"id": "HU3", "what": "not testable", "n": f"newly amplified {n_new} (< {MIN_NEW})"})
    ha, hn = her[her.y], her[~her.y]
    diff, p, half = boot_auc_diff(her.y.to_numpy(), her.p_gex.to_numpy(), her.p_as.to_numpy(), rng)
    add("HU4", diff, 1, p, f"amplified {len(ha)}, not {len(hn)}; 95% CI half-width {half:.3f}; "
        f"GEX {auc(ha.p_gex, hn.p_gex):.3f}, AS {auc(ha.p_as, hn.p_as):.3f}",
        "PAM50 Her2: AUROC(GEX) - AUROC(AS ensemble), amplified against not")
    add("HU5", auc(hn.p_as, a.p_as) - 0.5, 1, float(mannwhitneyu(hn.p_as, a.p_as).pvalue),
        f"non-amplified Her2 {len(hn)}, amplified luminal {len(a)}",
        "AUROC(AS P(Her2) higher in non-amplified Her2 than amplified LumA/LumB) - 0.5")
    level = a.rel if rule == "amp" else a.R
    r = spearmanr(a.p_gex, level)
    add("HU6", float(r.statistic), 1, float(r.pvalue), f"{len(a)} amplified luminal tumours",
        "amplified LumA/LumB: Spearman(GEX P(Her2), amplification level)")
    diff, p, half = boot_rho_diff(level.to_numpy(), a.p_gex.to_numpy(), a.p_as.to_numpy(), rng)
    add("HU7", diff, 1, p, f"{len(a)} tumours; 95% CI half-width {half:.3f}; AS rho "
        f"{spearmanr(a.p_as, level).statistic:.3f}", "amplified LumA/LumB: Spearman(GEX) - Spearman(AS ensemble)")
    pr = pd.DataFrame(pred)
    ok = pr.p.notna()
    pr.loc[ok, "q"] = multipletests(pr.loc[ok, "p"], method="fdr_bh")[1]
    if verdicts:
        pr["verdict"] = [m42.verdict(e, s, q) if np.isfinite(q) else "not testable"
                         for e, s, q in zip(pr.effect, pr.predicted_sign, pr.q.fillna(np.nan))]
    return pr


def exploratory(d: pd.DataFrame, rng) -> pd.DataFrame:
    ex = []
    raw = tests(d, "amp_raw", rng, verdicts=False)
    for r in raw.itertuples(index=False):
        ex.append({"what": f"raw rule (R >= 2): {r.id}, {r.what}", "value": r.effect, "p": r.p, "n": r.n})
    for rule, name in [("amp", "D-60 rule"), ("amp_raw", "raw rule"), ("amp_absolute", "D-58 ABSOLUTE rule")]:
        b = d[d.clin.notna() & d[rule].notna()]
        for c in ["positive", "negative"]:
            s = b[b.clin.eq(c)]
            ex.append({"what": f"{name}: amplified among clinical HER2 {c}", "value": s[rule].astype(bool).mean(),
                       "n": len(s)})
    u = d[d.cn_absolute.notna() & d.q.notna() & (d.cn_absolute < 7)]
    ex.append({"what": "Spearman(q, ABSOLUTE copy number) where ABSOLUTE is below 7",
               "value": spearmanr(u.q, u.cn_absolute).statistic, "n": len(u)})
    am = d[d.amp.eq(True)]
    ex.append({"what": "D-60 amplified: share with ABSOLUTE copy number at its ceiling (7)",
               "value": float((am.cn_absolute >= 7).mean()), "n": len(am)})
    ex.append({"what": "D-60 amplified: median q / ploidy", "value": float(am.rel.median()), "n": len(am)})
    g = d.genome_median.dropna()
    ex.append({"what": "convention check: median of per-sample genome-median Segment_Mean", "value": float(g.median()),
               "n": len(g)})
    ex.append({"what": "convention check: share of samples with |genome median| <= 0.1",
               "value": float((g.abs() <= 0.1).mean()), "n": len(g)})
    return pd.DataFrame(ex)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    d = tumours()
    pr = tests(d, "amp", rng, verdicts=True)
    pr.to_csv(OUT / "predictions.csv", index=False)
    log("predictions:\n" + pr[["id", "effect", "p", "q", "verdict", "n"]].to_string(index=False))
    ex = exploratory(d, rng)
    ex.to_csv(OUT / "exploratory.csv", index=False, float_format="%.4g")
    log("exploratory:\n" + ex.to_string(index=False))
    k = d[d.amp.notna()].assign(amplified=lambda x: x.amp.astype(bool))
    s = k.groupby(["label", "amplified"]).agg(n=("p_gex", "size"), median_p_gex=("p_gex", "median"),
                                              median_p_as=("p_as", "median"), median_rel=("rel", "median")).reset_index()
    s.loc[s.n < 5, ["median_p_gex", "median_p_as", "median_rel"]] = np.nan
    s.to_csv(OUT / "her2_probability_by_class_and_amplification.csv", index=False, float_format="%.4g")
    log(f"by class and amplification:\n{s.to_string(index=False)}")
    # luminal tumours by amplification group: who the ABSOLUTE ceiling left out (exploratory)
    lum = k[k.label.isin(LUM)]
    ab = lum.amp_absolute.astype("boolean").fillna(False).astype(bool)
    grp = np.where(ab, "D-58 amplified", np.where(lum.amplified, "newly amplified", "not amplified"))
    gt = lum.groupby(grp).agg(n=("ploidy", "size"), median_ploidy=("ploidy", "median"), median_rel=("rel", "median"),
                              median_absolute_cn=("cn_absolute", "median"),
                              lumb_share=("label", lambda s: (s == "LumB").mean()),
                              auroc_gex=("p_gex", lambda s: np.nan), auroc_as=("p_as", lambda s: np.nan))
    rest = lum[grp == "not amplified"]
    for g_ in ["D-58 amplified", "newly amplified"]:
        x = lum[grp == g_]
        gt.loc[g_, "auroc_gex"] = auc(x.p_gex, rest.p_gex)
        gt.loc[g_, "auroc_as"] = auc(x.p_as, rest.p_as)
    gt.loc[gt.n < 5, ["median_ploidy", "median_rel", "median_absolute_cn"]] = np.nan
    gt.rename_axis("group").reset_index().to_csv(OUT / "luminal_amplification_groups.csv", index=False,
                                                 float_format="%.4g")
    log(f"luminal groups:\n{gt.round(3).to_string()}")
    cn = k.drop_duplicates("s15").set_index("s15")[["amplified"]]
    m47.OUT = OUT                                     # the D-58 knockout, with this amplified group
    m47.knockout(cn)


if __name__ == "__main__":
    main()
