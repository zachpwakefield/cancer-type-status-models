"""Step 46 (D-57): inherited splicing variation in the BRCA subtypes - LOH and ancestry.

  scan      label-blind: genotype-like events in the 8 AS views (BRCA tumours and normals); each
            gene anchored to GRCh37 through its MC3 call span; arms split at each chromosome's
            largest MC3 gap
  test      L1-L6 (LOH, ABSOLUTE) and A1-A3 (ancestry) of config/followup.yaml: brca_germline,
            and the exploratory tables
  check     A-13, after the results: a germline check from tumour-normal pairs (a homozygous normal
            cannot flip to the opposite homozygote in its own tumour); the class profile of the
            models' top genotype-like events
  confirmed A-13: L1-L6 and A3 again, on the events that pass the germline check

Autosomes only: a tumour is a clone, so its X is inactivated on one homologue throughout and an X
heterozygote would look like LOH without any loss (A-11).

Outputs: data.nosync/ctsm_v1/brca/ (per sample, local); results/subtype/brca/germline/ (aggregates)
Usage: python scripts/46_brca_germline.py --only scan|test|check|confirmed
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from _common import REPO, log
from scipy.special import gammaln
from scipy.stats import binomtest, chi2, fisher_exact, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

from ctsm.biology import cohort_matrix, gex_by_symbol
from ctsm.config import load_config

_spec = importlib.util.spec_from_file_location("msi42", Path(__file__).with_name("42_msi_repeats.py"))
m42 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m42)

CFG = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["brca_germline"]
DATA = m42.DATA
EXT = DATA / "external" / "gdc_pancanatlas"
LOCAL = DATA / "brca"
OUT = REPO / "results" / "subtype" / "brca" / "germline"
VIEWS = ["se", "a3ss", "a5ss", "mxe", "ri", "afe", "ale", "hit"]
GRID = np.round(np.arange(0.02, 0.99, 0.02), 2)
MIN_STEP, MIN_GROUP = 4, 5            # thresholds at least 4 grid steps (0.08) apart; >= 5 per group
BAND, MIN_ARM = 0.15, 5
ACROCENTRIC = {"13", "14", "15", "21", "22"}
SEED, BOOT, PERM = 0, 2000, 1000
MAX_FLIP, MIN_HOM_PAIRS = 0.05, 20       # A-13 germline check
EXAMPLES = ["CAST", "TMEM101"]           # U12: a germline event and the models' top genotype-like event


# ============================================================ scan (label-blind)
def three_groups(v: np.ndarray):
    """Thresholds (a, b) on GRID minimising the within-group sum of squares of sorted v."""
    n = len(v)
    s1, s2 = np.concatenate([[0], np.cumsum(v)]), np.concatenate([[0], np.cumsum(v * v)])
    cut = np.searchsorted(v, GRID)                       # number of values < each grid point
    i, k = np.triu_indices(len(GRID), MIN_STEP)
    ci, ck = cut[i], cut[k]
    n0, n1, n2 = ci, ck - ci, n - ck
    ok = (n0 >= MIN_GROUP) & (n1 >= MIN_GROUP) & (n2 >= MIN_GROUP)
    if not ok.any():
        return None
    i, k, ci, ck, n0, n1, n2 = i[ok], k[ok], ci[ok], ck[ok], n0[ok], n1[ok], n2[ok]

    def sse(a, b, m):
        return (s2[b] - s2[a]) - (s1[b] - s1[a]) ** 2 / m

    tot = sse(0, ci, n0) + sse(ci, ck, n1) + sse(ck, n, n2)
    j = int(np.argmin(tot))
    a, b = ci[j], ck[j]
    mu = [v[:a].mean(), v[a:b].mean(), v[b:].mean()]
    return GRID[i[j]], GRID[k[j]], mu, float(np.sqrt(max(tot[j], 0.0) / n))   # rounding can dip below 0


def scan_view(cfg, view: str) -> pd.DataFrame:
    X, F, tab = cohort_matrix(cfg, view, "tumor_status", "BRCA")
    nor = tab["role"].eq("normal").to_numpy()
    ok = (np.isfinite(X).mean(0) >= 0.9) & (np.nanstd(X, 0) > 0.12)
    rows = []
    for j in np.flatnonzero(ok):
        x = X[:, j]
        fin = np.isfinite(x)
        g = three_groups(np.sort(x[fin]))
        if g is None:
            continue
        a, b, mu, sd = g
        rng = mu[2] - mu[0]
        xn = x[fin & nor]
        c = np.array([(xn < a).sum(), ((xn >= a) & (xn < b)).sum(), (xn >= b).sum()])
        n = c.sum()
        p = (2 * c[0] + c[1]) / (2 * n)
        e = n * np.array([p * p, 2 * p * (1 - p), (1 - p) ** 2])
        hwe = float(chi2.sf(((c - e) ** 2 / np.maximum(e, 1e-9)).sum(), 1))
        xpos = (xn - mu[0]) / rng
        rows.append({"view": view, "feature_id": F["feature_id"].iat[j], "gene_id": F["gene_id"].iat[j],
                     "symbol": F["symbol"].iat[j], "chrom": F["chrom"].iat[j], "region": F["region"].iat[j],
                     "a": a, "b": b, "c_low": mu[0], "c_mid": mu[1], "c_high": mu[2], "sd": sd,
                     "spacing": abs((mu[1] - mu[0]) - (mu[2] - mu[1])) / rng, "gap": min(mu[1] - mu[0], mu[2] - mu[1]) / sd,
                     "n_low": c[0], "n_mid": c[1], "n_high": c[2], "maf": min(p, 1 - p), "hwe_p": hwe,
                     "band_rate_normals": float((np.abs(xpos - 0.5) < BAND).mean())})
    r = pd.DataFrame(rows)
    r["genotype_like"] = (r.spacing < 0.25) & (r.gap > 3) & (r.maf > 0.1) & (r.hwe_p > 0.05)
    log(f"{view}: {ok.sum():,} eligible, {int(r.genotype_like.sum())} genotype-like")
    return r


def mc3_anchors() -> tuple[pd.DataFrame, pd.Series]:
    """Per Ensembl gene: GRCh37 chromosome and the midpoint of its MC3 call span; per chromosome:
    the centromere, taken as the largest gap between consecutive MC3 positions."""
    m = pd.read_csv(EXT / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=["Gene", "Chromosome", "Start_Position"],
                    dtype={"Gene": str, "Chromosome": str, "Start_Position": np.int64})
    m = m[m.Chromosome.isin([str(i) for i in range(1, 23)])]
    cen = {}
    for c, d in m.groupby("Chromosome"):
        p = np.unique(d.Start_Position.to_numpy())
        k = int(np.argmax(np.diff(p)))
        cen[c] = (p[k] + p[k + 1]) / 2
    g = m.groupby("Gene").agg(chrom19=("Chromosome", lambda s: s.mode().iat[0]), lo=("Start_Position", "min"),
                              hi=("Start_Position", "max"))
    g["anchor"] = (g.lo + g.hi) / 2
    g["arm"] = [f"{c}{'q' if (c in ACROCENTRIC or a > cen[c]) else 'p'}" for c, a in zip(g.chrom19, g.anchor)]
    return g, pd.Series(cen)


def scan() -> None:
    cfg = load_config()
    LOCAL.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    r = pd.concat([scan_view(cfg, v) for v in VIEWS], ignore_index=True)
    anc, cen = mc3_anchors()
    r = r.join(anc[["chrom19", "anchor", "arm"]], on="gene_id")
    same = ("chr" + r.chrom19.astype(str)).eq(r.chrom)
    r.loc[~same, ["chrom19", "anchor", "arm"]] = np.nan
    r.to_parquet(LOCAL / "scan.parquet", index=False)
    gl = r[r.genotype_like].copy()
    gl["width"] = gl.c_high - gl.c_low
    gl = gl.sort_values("width", ascending=False)
    gl["gene_pick"] = ~gl.duplicated("gene_id")
    gl.to_csv(OUT / "genotype_like_events.csv", index=False, float_format="%.4g")
    cen.rename("centromere_grch37").to_csv(OUT / "centromeres_from_mc3.csv")
    one = gl[gl.gene_pick]
    log(f"genotype-like: {len(gl)} events, {len(one)} genes; anchored autosomal genes "
        f"{int(one.anchor.notna().sum())}; by view " + ", ".join(f"{k} {v}" for k, v in gl.view.value_counts().items()))


# ============================================================ test
def sample_positions(cfg, ev: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Allelic position x (samples x genes) of the chosen event of each gene, and the task table."""
    cols = {}
    for v, d in ev.groupby("view"):
        X, F, tab = cohort_matrix(cfg, v, "tumor_status", "BRCA")
        pos = pd.Series(np.arange(len(F)), index=F.feature_id)
        for r in d.itertuples(index=False):
            cols[r.gene_id] = (X[:, pos[r.feature_id]] - r.c_low) / (r.c_high - r.c_low)
    return pd.DataFrame(cols, index=tab["file_id"]), tab


def absolute_at(samples: list[str], genes: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Per tumour sample (15-char barcode): LOH, HSCN1, HSCN2 at each gene's anchor."""
    seg = pd.read_csv(EXT / "TCGA_mastercalls.abs_segtabs.fixed.txt", sep="\t",
                      usecols=["Sample", "Chromosome", "Start", "End", "Modal_HSCN_1", "Modal_HSCN_2", "LOH"])
    seg = seg[seg.Sample.isin(samples) & seg.Chromosome.between(1, 22)]
    seg["chrom"] = seg.Chromosome.astype(int).astype(str)
    out = {}
    for s, d in seg.groupby("Sample"):
        res = pd.DataFrame(np.nan, index=genes.index, columns=["loh", "c1", "c2"])
        for c, dc in d.groupby("chrom"):
            dc = dc.sort_values("Start")
            gi = genes.index[genes.chrom19.eq(c)]
            a = genes.loc[gi, "anchor"].to_numpy()
            k = np.searchsorted(dc.Start.to_numpy(), a, side="right") - 1
            hit = (k >= 0) & (a <= dc.End.to_numpy()[np.clip(k, 0, None)])
            kk = k[hit]
            res.loc[gi[hit], "loh"] = dc.LOH.to_numpy()[kk]
            res.loc[gi[hit], "c1"] = dc.Modal_HSCN_1.to_numpy()[kk]
            res.loc[gi[hit], "c2"] = dc.Modal_HSCN_2.to_numpy()[kk]
        out[s] = res
    return out


def d_pred(p: float, c1, c2):
    big, tot = np.maximum(c1, c2), c1 + c2
    f = (p * big + 1 - p) / (p * tot + 2 * (1 - p))
    return np.where(tot > 0, f - 0.5, np.nan)


def _lchoose(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def fisher_p(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p of [[a, b], [c, d]], as scipy.stats.fisher_exact (tables with probability
    up to 1 + 1e-7 times the observed one), on log-gamma for speed (A3 runs it a million times)."""
    M, n, N = a + b + c + d, a + b, a + c
    k = np.arange(max(0, N - (M - n)), min(n, N) + 1)
    lp = _lchoose(n, k) + _lchoose(M - n, N - k) - _lchoose(M, N)
    lp0 = _lchoose(n, a) + _lchoose(M - n, N - a) - _lchoose(M, N)
    return float(min(1.0, np.exp(lp[lp <= lp0 + np.log1p(1e-7)]).sum()))


def otsu(v: np.ndarray, bins: int = 256) -> float:
    h, e = np.histogram(v, bins=bins)
    c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h)
    w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1)
    m1 = ((h * c).sum() - np.cumsum(h * c)) / np.maximum(w1, 1)
    s = w0 * w1 * (m0 - m1) ** 2
    return float(c[s >= s.max() * (1 - 1e-9)].mean())         # middle of the plateau when the gap is empty


def log_or(t: np.ndarray) -> float:
    a, b, c, d = t.ravel() + 0.5
    return float(np.log(a * d / (b * c)))


def mh_log_or(tables: list[np.ndarray]) -> float:
    num = den = 0.0
    for t in tables:
        a, b, c, d = t.ravel() + 0.5
        n = a + b + c + d
        num += a * d / n
        den += b * c / n
    return float(np.log(num / den))


def auroc(pos: np.ndarray, neg: np.ndarray) -> tuple[float, float]:
    r = mannwhitneyu(pos, neg)
    return float(r.statistic / (len(pos) * len(neg))), float(r.pvalue)


def test(confirmed: bool = False) -> None:
    cfg = load_config()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    tag = "_germline_confirmed" if confirmed else ""
    gl = pd.read_csv(OUT / "genotype_like_events.csv")
    ev_all = gl[gl.gene_pick].set_index("gene_id", drop=False)
    if confirmed:
        chk = pd.read_csv(OUT / "germline_check.csv").set_index("gene_id")
        ev_all = ev_all[ev_all.index.isin(chk.index[chk.confirmed])]
    X_all, tab = sample_positions(cfg, ev_all)
    ev_all = ev_all.loc[X_all.columns]
    ev = ev_all[ev_all.anchor.notna()].copy()                 # autosomal genes with a GRCh37 anchor
    ev["chrom19"] = pd.to_numeric(ev.chrom19).astype(int).astype(str)
    X = X_all[ev.index]
    tab = tab.set_index("file_id")
    tab["s15"] = tab.sample_id.str[:15]
    tum = tab[tab.role.eq("tumour")].sort_values("sample_id")
    tum = tum[~tum.case_id.duplicated()]                       # one tumour per patient
    nor = tab[tab.role.eq("normal")]
    by_case = lambda f: pd.read_parquet(DATA / "subtype" / "tasks" / f).drop_duplicates("case_id").set_index("case_id")["label"]
    pur = pd.read_csv(EXT / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
    pur = pur[pur["call status"].eq("called")].set_index("array")["purity"]
    called = tum[tum.s15.isin(pur.index)]
    ab = absolute_at(called.s15.tolist(), ev[["chrom19", "anchor"]])
    log(f"{len(ev)} genes, {len(tum)} tumours ({len(ab)} with ABSOLUTE), {len(nor)} normals")
    pred, expl = [], []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    # ---- L1-L3: tumour-normal pairs
    lo, mid_hi = ev.a, ev.b                      # thresholds in PSI units
    rows_pair, per_t = [], []
    for fid, r in called[called.case_id.isin(nor.case_id)].iterrows():
        if r.s15 not in ab:
            continue
        nf = nor.index[nor.case_id.eq(r.case_id)][0]
        psi_n = X.loc[nf] * (ev.c_high - ev.c_low) + ev.c_low
        het = (psi_n >= lo) & (psi_n < mid_hi)
        d = (X.loc[fid] - 0.5).abs()
        a = ab[r.s15]
        p = float(pur[r.s15])
        dp = pd.Series(d_pred(p, a.c1.to_numpy(), a.c2.to_numpy()), index=a.index)
        k = het & d.notna() & a.loh.notna()
        if k.sum() == 0:
            continue
        sub = pd.DataFrame({"d": d[k], "loh": a.loh[k], "d_pred": dp[k]})
        rows_pair.append(sub.assign(case_id=r.case_id, purity=p))
        n1, n0 = int(sub.loh.eq(1).sum()), int(sub.loh.eq(0).sum())
        rec = {"case_id": r.case_id, "purity": p, "het_events": len(sub), "het_loh": n1, "het_not": n0,
               "diff": sub.d[sub.loh.eq(1)].mean() - sub.d[sub.loh.eq(0)].mean() if min(n1, n0) >= 3 else np.nan}
        s3 = sub.dropna(subset=["d_pred"])
        rec["rho_d_dpred"] = (spearmanr(s3.d, s3.d_pred).statistic
                              if len(s3) >= 10 and s3.d_pred.nunique() >= 2 else np.nan)
        per_t.append(rec)
    pt = pd.DataFrame(per_t)
    pairs = pd.concat(rows_pair, ignore_index=True)
    pairs.to_parquet(LOCAL / f"pair_events{tag}.parquet", index=False)
    x1 = pt["diff"].dropna()
    k1 = int((x1 > 0).sum())
    add("L1", k1 / len(x1) - 0.5, 1, binomtest(k1, len(x1)).pvalue,
        f"{k1} of {len(x1)} tumours; median difference {x1.median():.3f}",
        "share of tumours whose heterozygous events sit further from the middle inside LOH, minus 0.5")
    l2 = pt.dropna(subset=["diff"])
    r2 = spearmanr(l2["diff"], l2.purity)
    add("L2", float(r2.statistic), 1, float(r2.pvalue), f"{len(l2)} tumours", "Spearman(L1 difference, purity)")
    x3 = pt.rho_d_dpred.dropna()
    k3 = int((x3 > 0).sum())
    add("L3", k3 / len(x3) - 0.5, 1, binomtest(k3, len(x3)).pvalue,
        f"{k3} of {len(x3)} tumours; median rho {x3.median():.3f}", "share of tumours with Spearman(d, d_pred) > 0, minus 0.5")

    # ---- depletion and arm LOH (tumour-only)
    band = ((X - 0.5).abs() < BAND).astype(float).where(X.notna())
    exp_rate = ev.band_rate_normals
    arms = ev.groupby("arm").gene_id.apply(list)
    arms = arms[arms.map(len) >= MIN_ARM]
    dep = pd.DataFrame({arm: 1 - band.loc[tum.index, g].sum(1) / X.loc[tum.index, g].notna().mul(exp_rate[g]).sum(1)
                        for arm, g in arms.items()})
    aloh = pd.DataFrame(np.nan, index=tum.index, columns=arms.index)
    for fid, r in called.iterrows():
        if r.s15 in ab:
            a = ab[r.s15].loh
            aloh.loc[fid] = [a[g].mean() for g in arms]
    dep.to_parquet(LOCAL / f"arm_depletion{tag}.parquet")
    aloh.to_parquet(LOCAL / f"arm_loh{tag}.parquet")
    lab = tum.case_id.map(by_case("subtype__BRCA_PAM50.parquet"))
    h = tum.case_id.map(by_case("subtype__BRCA.parquet"))
    rows4 = []
    for arm in arms.index:
        la, ba = lab.eq("LumA"), lab.eq("Basal")
        hasA = aloh[arm].notna()
        rows4.append({"arm": arm, "genes": len(arms[arm]),
                      "dep_diff": dep.loc[la, arm].mean() - dep.loc[ba, arm].mean(),
                      "loh_diff": (aloh.loc[la & hasA, arm] >= 0.5).mean() - (aloh.loc[ba & hasA, arm] >= 0.5).mean()})
    a4 = pd.DataFrame(rows4)
    r4 = spearmanr(a4.dep_diff, a4.loh_diff)
    add("L4", float(r4.statistic), 1, float(r4.pvalue), f"{len(a4)} arms",
        "Spearman over arms of LumA-minus-Basal depletion and LumA-minus-Basal LOH share")
    if "16q" in dep:
        au, p5 = auroc(dep.loc[h.eq("Lobular"), "16q"].dropna(), dep.loc[h.eq("Ductal"), "16q"].dropna())
        add("L5", au - 0.5, 1, p5, f"lobular {int(h.eq('Lobular').sum())}, ductal {int(h.eq('Ductal').sum())}",
            "AUROC(16q depletion higher in lobular than ductal) - 0.5")
    rho6 = []
    for arm in arms.index:
        k = dep[arm].notna() & aloh[arm].notna()
        rho6.append({"arm": arm, "rho": spearmanr(dep.loc[k, arm], aloh.loc[k, arm]).statistic, "n": int(k.sum())})
    a6 = pd.DataFrame(rho6).merge(a4, on="arm")
    k6 = int((a6.rho > 0).sum())
    add("L6", k6 / len(a6) - 0.5, 1, binomtest(k6, len(a6)).pvalue,
        f"{k6} of {len(a6)} arms; median rho {a6.rho.median():.3f}", "share of arms with Spearman(depletion, arm LOH) > 0, minus 0.5")
    a6.to_csv(OUT / f"arms{tag}.csv", index=False, float_format="%.4g")

    # ---- ancestry
    cl = pd.read_csv(EXT / "clinical_PANCAN_patient_with_followup.tsv", sep="\t", usecols=["bcr_patient_barcode", "race"],
                     encoding="latin-1").set_index("bcr_patient_barcode")["race"]
    race_case = cl.map({"BLACK OR AFRICAN AMERICAN": "Black", "WHITE": "white"})
    race = tab.case_id.map(race_case)
    gst = gex_by_symbol(cfg, tab.index, ["GSTM1"])["GSTM1"]
    thr = otsu(gst.dropna().to_numpy())
    null = gst < thr
    tt = tum.index
    bw = race.reindex(tt).notna()
    t1 = pd.crosstab(race[tt][bw], null[tt][bw]).reindex(index=["Black", "white"], columns=[True, False]).fillna(0).to_numpy()
    r1 = fisher_exact(t1)
    add("A1", log_or(t1), -1, float(r1.pvalue), f"Black null {int(t1[0, 0])}/{int(t1[0].sum())}, white null "
        f"{int(t1[1, 0])}/{int(t1[1].sum())}; threshold log2CPM {thr:.2f}", "log odds ratio of GSTM1 null, Black vs white")
    er = tum.case_id.map(by_case("subtype__BRCA_ER.parquet"))
    d2 = pd.DataFrame({"present": ~null[tt], "erneg": er.eq("ER-negative"), "race": race[tt], "has_er": er.notna()})
    d2 = d2[d2.has_er & d2.race.notna()].reset_index(drop=True)

    def tabs(d):
        t = lambda x: np.array([[(x.present & x.erneg).sum(), (x.present & ~x.erneg).sum()],
                                [(~x.present & x.erneg).sum(), (~x.present & ~x.erneg).sum()]], float)
        return log_or(t(d)), mh_log_or([t(d[d.race.eq(g)]) for g in ("Black", "white")])

    crude, mh = tabs(d2)
    bs = []
    for _ in range(BOOT):
        bs.append(np.subtract(*tabs(d2.iloc[rng.integers(0, len(d2), len(d2))])))
    bs = np.asarray(bs)
    add("A2", crude - mh, 1, float(min(1.0, 2 * (1 + min((bs <= 0).sum(), (bs >= 0).sum())) / (1 + len(bs)))),
        f"{len(d2)} patients; crude log OR {crude:.3f}, race-stratified {mh:.3f}",
        "crude minus Mantel-Haenszel log OR of ER-negative for GSTM1 present")
    # A3: allele frequency by race, tumour groups, every genotype-like gene
    psi = X_all.loc[tt] * (ev_all.c_high - ev_all.c_low) + ev_all.c_low
    geno = np.where(psi.isna(), np.nan, (psi >= ev_all.a).astype(float) + (psi >= ev_all.b).astype(float))
    geno = pd.DataFrame(geno, index=tt, columns=ev_all.index)
    rb = race[tt][bw]
    G = geno.loc[rb.index]

    def n_sig(lbl: np.ndarray) -> int:
        ps = []
        for g in G.columns:
            v = G[g].to_numpy()
            ok = ~np.isnan(v)
            t = []
            for grp in ("Black", "white"):
                w = v[ok & (lbl == grp)]
                t += [int(w.sum()), int(2 * len(w) - w.sum())]
            ps.append(fisher_p(*t))
        return int((multipletests(ps, method="fdr_bh")[1] < 0.05).sum())

    lbl = rb.to_numpy()
    obs = n_sig(lbl)
    null_counts = np.array([n_sig(rng.permutation(lbl)) for _ in range(PERM)])
    add("A3", float(obs - np.median(null_counts)), 1, float((1 + (null_counts >= obs).sum()) / (1 + PERM)),
        f"{obs} of {G.shape[1]} genes; permutation median {np.median(null_counts):.0f}, max {null_counts.max()}",
        "genes with a race allele-frequency difference at q < 0.05, minus the permutation median")

    pr = pd.DataFrame(pred)
    pr["family"] = pr.id.str[0]
    pr["q"] = pr.groupby("family").p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
    pr["verdict"] = [m42.verdict(e, s, q) for e, s, q in zip(pr.effect, pr.predicted_sign, pr.q)]
    pr["reading"] = ""
    if "L5" not in set(pr.id):                               # A-13 sensitivity: fewer anchored events
        pr = pd.concat([pr, pd.DataFrame([{"id": "L5", "family": "L", "verdict": "not testable",
                                           "n": f"16q has fewer than {MIN_ARM} anchored events"}])], ignore_index=True)
        pr = pr.sort_values("id", key=lambda s: s.map({k: i for i, k in enumerate(
            ["L1", "L2", "L3", "L4", "L5", "L6", "A1", "A2", "A3"])}), ignore_index=True)
    if crude < 0:                                            # A-13 (2)
        pr.loc[pr.id.eq("A2"), "reading"] = ("against the claim as worded: the crude log OR is negative, so a positive "
                                             "difference means the association grew under race stratification (A-13)")
    pr.to_csv(OUT / f"predictions{tag}.csv", index=False)
    log("predictions:\n" + pr[["id", "effect", "p", "q", "verdict", "n"]].to_string(index=False))
    pt.drop(columns="case_id").describe().T.to_csv(OUT / f"pairs_summary{tag}.csv", float_format="%.4g")
    if not confirmed:
        exploratory(lab, dep, aloh, arms, race_case)


def pair_kappa(a: np.ndarray, b: np.ndarray) -> float:
    """Cohen's kappa between the normal's and the tumour's group (0/1/2) over complete pairs."""
    k = ~np.isnan(a) & ~np.isnan(b)
    if k.sum() < MIN_HOM_PAIRS:
        return np.nan
    a, b = a[k].astype(int), b[k].astype(int)
    po = float((a == b).mean())
    pe = float(sum((a == g).mean() * (b == g).mean() for g in range(3)))
    return (po - pe) / (1 - pe) if pe < 1 else np.nan


def check() -> None:
    """A-13: germline check over tumour-normal pairs (role only), and the class profile of the subtype
    models' top genotype-like events."""
    cfg = load_config()
    gl = pd.read_csv(OUT / "genotype_like_events.csv")
    ev = gl[gl.gene_pick].set_index("gene_id", drop=False)
    X, tab = sample_positions(cfg, ev)
    ev = ev.loc[X.columns]
    tab = tab.set_index("file_id")
    psi = (X * (ev.c_high - ev.c_low) + ev.c_low).to_numpy()
    G = pd.DataFrame(np.where(np.isnan(psi), np.nan, (psi >= ev.a.to_numpy()).astype(float)
                              + (psi >= ev.b.to_numpy()).astype(float)), index=X.index, columns=X.columns)
    tum = tab[tab.role.eq("tumour")].sort_values("sample_id")
    tum = tum[~tum.case_id.duplicated()]
    nor = tab[tab.role.eq("normal")]
    nor = nor[~nor.case_id.duplicated()]
    pr = tum.reset_index().merge(nor.reset_index(), on="case_id", suffixes=("_t", "_n"))
    gt, gn = G.loc[pr.file_id_t].to_numpy(), G.loc[pr.file_id_n].to_numpy()
    hom = np.isin(gn, [0, 2]) & ~np.isnan(gt)
    flip = hom & (np.abs(gt - gn) == 2)
    chk = pd.DataFrame({"gene_id": X.columns, "symbol": ev.symbol.to_numpy(), "view": ev["view"].to_numpy(),
                        "region": ev.region.to_numpy(), "arm": ev.arm.to_numpy(), "homozygous_pairs": hom.sum(0),
                        "flips": flip.sum(0)})
    chk["flip_rate"] = chk.flips / chk.homozygous_pairs.clip(lower=1)
    chk["af_tumours"] = (G.loc[tum.index].mean() / 2).to_numpy()
    chk["af_normals"] = (G.loc[nor.index].mean() / 2).to_numpy()
    chk["confirmed"] = (chk.homozygous_pairs >= MIN_HOM_PAIRS) & (chk.flip_rate <= MAX_FLIP)
    chk["pair_kappa"] = [pair_kappa(gn[:, j], gt[:, j]) for j in range(gt.shape[1])]   # A-14, exploratory
    by_case = lambda f: pd.read_parquet(DATA / "subtype" / "tasks" / f).drop_duplicates("case_id").set_index("case_id")["label"]
    pam = tum.case_id.map(by_case("subtype__BRCA_PAM50.parquet"))
    for c in ["LumA", "LumB", "Her2", "Basal"]:
        chk[f"af_{c}"] = (G.loc[pam.index[pam.eq(c)]].mean() / 2).to_numpy()
    chk.to_csv(OUT / "germline_check.csv", index=False, float_format="%.4g")
    # pair tables (normal group x own tumour group) for the two examples of U12
    ex = []
    for sym in EXAMPLES:
        j = int(np.flatnonzero(chk.symbol.eq(sym).to_numpy())[0])
        k = ~np.isnan(gn[:, j]) & ~np.isnan(gt[:, j])
        tab2 = pd.crosstab(gn[k, j], gt[k, j]).reindex(index=[0, 1, 2], columns=[0, 1, 2], fill_value=0)
        ex += [{"symbol": sym, "view": chk.view.iat[j], "region": chk.region.iat[j], "normal_group": gi,
                "tumour_group": gj, "pairs": int(tab2.loc[gi, gj])} for gi in range(3) for gj in range(3)]
    pd.DataFrame(ex).to_csv(OUT / "pair_tables_examples.csv", index=False)
    log(f"germline check: {len(chk)} genes, median flip rate {chk.flip_rate.median():.3f}; confirmed "
        f"{int(chk.confirmed.sum())} (anchored {int((chk.confirmed & chk.arm.notna()).sum())})")
    # class profile of the models' top genotype-like events, and the confirmed share among top events
    top = pd.read_csv(OUT / "top_events_genotype_like_list.csv")
    best = top.groupby(["view", "region"]).mean_abs_shap.max().rename("best_shap").reset_index()
    prof = best.merge(chk, on=["view", "region"]).sort_values("best_shap", ascending=False)
    rows = []
    for r in prof.itertuples(index=False):
        g = G.loc[tum.index, r.gene_id]
        rec = r._asdict()
        rec["het_normals"] = float((G.loc[nor.index, r.gene_id] == 1).mean())
        for c in ["LumA", "LumB", "Her2", "Basal"]:
            k = pam.eq(c).to_numpy()
            rec[f"het_{c}"] = float((g[k] == 1).mean())
            rec[f"high_allele_{c}"] = float(g[k].mean() / 2)
        rows.append(rec)
    pd.DataFrame(rows).to_csv(OUT / "top_events_genotype_like_profile.csv", index=False, float_format="%.4g")
    scan = pd.read_parquet(LOCAL / "scan.parquet")
    conf = set(chk.loc[chk.confirmed, "gene_id"])
    top_all = pd.read_csv(REPO / "results" / "subtype" / "interpret" / "top_features_by_class.csv")
    top_all = top_all[top_all.cohort.isin(["BRCA", "BRCA_ER", "BRCA_PAM50"]) & top_all.view.ne("gex")]
    t = top_all.drop(columns="gene_id").merge(scan[["view", "feature_id", "gene_id", "genotype_like"]],
                                              on=["view", "feature_id"])
    pick = gl.loc[gl.gene_pick, ["view", "feature_id"]].assign(pick=True)
    t = t.merge(pick, on=["view", "feature_id"], how="left")
    scan = scan.merge(pick, on=["view", "feature_id"], how="left")
    rows = []
    for name, d in [("top events (3 BRCA cohorts)", t), ("all eligible events", scan)]:
        gk = d.genotype_like.astype(bool)
        pk = gk & d.pick.fillna(False).astype(bool)          # the check covers one event per gene
        ck = pk & d.gene_id.isin(conf)
        rows.append({"set": name, "events": len(d), "genotype_like": int(gk.sum()), "checked": int(pk.sum()),
                     "confirmed_germline": int(ck.sum())})
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "top_events_germline_share.csv", index=False)
    head = pd.read_csv(OUT / "top_events_genotype_like_profile.csv").head(12)
    cols = ["symbol", "view", "flip_rate", "high_allele_LumA", "high_allele_Basal", "het_normals", "het_LumA", "het_Basal"]
    log(f"top-event shares:\n{s.to_string(index=False)}\nprofile head:\n{head[cols].round(2).to_string(index=False)}")


def balanced_accuracy(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean([(p[y == c] == c).mean() for c in np.unique(y)]))


def exploratory(lab, dep, aloh, arms, race_case) -> None:
    # depletion and LOH maps by PAM50 class
    rows = []
    for c in ["LumA", "LumB", "Her2", "Basal", "Normal-like"]:
        k = lab.eq(c)
        for arm in arms.index:
            ha = k & aloh[arm].notna()
            rows.append({"class": c, "arm": arm, "n": int(k.sum()), "mean_depletion": dep.loc[k, arm].mean(),
                         "n_absolute": int(ha.sum()), "loh_share": (aloh.loc[ha, arm] >= 0.5).mean()})
    pd.DataFrame(rows).to_csv(OUT / "arm_maps_by_pam50.csv", index=False, float_format="%.4g")
    # genotype-like events among the models' top events
    top = pd.read_csv(REPO / "results" / "subtype" / "interpret" / "top_features_by_class.csv")
    top = top[top.cohort.isin(["BRCA", "BRCA_ER", "BRCA_PAM50"]) & top.view.ne("gex")]
    scan = pd.read_parquet(LOCAL / "scan.parquet")
    elig = scan[["view", "feature_id", "genotype_like"]]
    t = top.merge(elig, on=["view", "feature_id"], how="left")
    t["eligible"] = t.genotype_like.notna()
    summ = t.groupby("cohort").agg(top_events=("feature_id", "size"), eligible=("eligible", "sum"),
                                   genotype_like=("genotype_like", lambda s: int(s.fillna(False).sum()))).reset_index()
    summ["share_all_eligible"] = scan.genotype_like.mean()
    summ.to_csv(OUT / "top_events_genotype_like.csv", index=False, float_format="%.4g")
    t[t.genotype_like.fillna(False).astype(bool)][["cohort", "class", "view", "rank_in_class", "symbol", "region",
                                                   "mean_abs_shap"]].to_csv(OUT / "top_events_genotype_like_list.csv",
                                                                            index=False, float_format="%.4g")
    log(f"top events genotype-like:\n{summ.to_string(index=False)}")
    # out-of-fold balanced accuracy by race
    rr = []
    for coh in ["BRCA", "BRCA_ER", "BRCA_PAM50"]:
        P = {}
        for v in VIEWS + ["gex"]:
            o = pd.read_parquet(DATA / "subtype" / "runs" / "subtype" / v / coh / "oof_predictions.parquet")
            o = o.sort_values("file_id").set_index("file_id")
            P[v] = o
        y = P["gex"].y.to_numpy()
        pc = [c for c in P["gex"].columns if c.startswith("p_")]
        mods = {"GEX": P["gex"][pc].to_numpy(), "AS ensemble": np.mean([P[v][pc].to_numpy() for v in VIEWS], axis=0)}
        rc = P["gex"].case_id.map(race_case).to_numpy()
        for m, prob in mods.items():
            hat = prob.argmax(1)
            for g in ["Black", "white"]:
                k = (rc == g)
                rr.append({"cohort": coh, "model": m, "race": g, "n": int(k.sum()),
                           "balanced_accuracy": balanced_accuracy(y[k], hat[k]) if k.sum() >= 5 else np.nan})
    rr = pd.DataFrame(rr)
    rr.to_csv(OUT / "oof_by_race.csv", index=False, float_format="%.4g")
    log(f"OOF balanced accuracy by race:\n{rr.pivot_table(index=['cohort', 'model'], columns='race', values='balanced_accuracy').round(3)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["scan", "test", "check", "confirmed"], required=True)
    a = ap.parse_args()
    {"scan": scan, "test": test, "check": check, "confirmed": lambda: test(confirmed=True)}[a.only]()


if __name__ == "__main__":
    main()
