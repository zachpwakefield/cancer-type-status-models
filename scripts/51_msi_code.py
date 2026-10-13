"""Step 51 (D-65..D-67): the MSI sequence code: specificity, thirty held-out cancers, endometrium.

The code (config/followup.yaml, shared definitions above msi_code_specificity): repeat exons whose
acceptor carries a T-run of 10-19 nt are skipped in mismatch-repair-deficient tumours (R_down);
those with a run >= 20 nt and a branch point >= 36 nt from the AG gain inclusion (R_up). Score G
reads both arms; score A is D-52's, which reads all repeat exons as falling.

  annotate     label-blind: D-52's acceptor annotation for every exon of the pan-cancer SE matrix
  geometry     label-blind: run and branch point (D-63) for R exons, branch point for B exons, of
               the subtype and pan-cancer matrices
  specificity  D-65 (BS1-BS4)
  pancancer    D-66 (PC1-PC8): the 30 TCGA projects other than COAD, STAD and UCEC
  endometrium  D-67 (EG1-EG2)
  infiltration D-68 (IF1-IF6): score G against stromal, immune, retained-intron and purity covariates
  msicalls     D-69 (MC1-MC3): score G against MANTIS and MSIsensor calls (step 00c fetches them)

Outputs: data.nosync/ctsm_v1/msi/ (per-tumour scores, kept local)
         results/subtype/msi_cis/code/ (aggregates)
Usage: python scripts/51_msi_code.py --only annotate|geometry|specificity|pancancer|endometrium|
                                             infiltration|msicalls
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from _common import REPO, log
from scipy.stats import rankdata, spearmanr


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m42 = _load("msi42", "42_msi_repeats.py")
m50 = _load("msi50", "50_msi_followup.py")
from ctsm.matrix_io import load_view  # noqa: E402

DATA, LOCAL = m42.DATA, m42.LOCAL
PAN = DATA / "matrices"
OUT = REPO / "results" / "subtype" / "msi_cis" / "code"
TRAIN = ["COAD", "STAD", "UCEC"]
RELATIVES = ["READ", "ESCA", "UCS"]
UP_CUT, MMR_CUT, HYPER_SNV = 36, 10, 500
SEED, BOOT = 0, 2000


# ---------------------------------------------------------------- label-blind annotation

def annotate() -> None:
    m42.annotate(features=PAN / "se" / "features.parquet", out=LOCAL / "se_repeats_pan.parquet")


def geometry() -> None:
    a = pd.concat([pd.read_parquet(LOCAL / "se_repeats.parquet"), pd.read_parquet(LOCAL / "se_repeats_pan.parquet")])
    a = a.drop_duplicates("feature_id")
    a = a[(a.acc_T >= 10) | (a.acc_T <= 5)]
    rows = []
    for chrom, seq in m42.read_genome(set(a.chrom)):
        for r in a[a.chrom.eq(chrom)].itertuples(index=False):
            if r.strand == "+":
                up = seq[max(int(r.left_end), int(r.es) - m50.UP):int(r.es)]
            else:
                up = seq[int(r.ee):min(int(r.right_start), int(r.ee) + m50.UP)][::-1].translate(m50.COMP)
            rev = up[::-1]
            if r.acc_T >= 10:
                g = m50.run_geometry(rev)
                rec = {"cls": "R", "d_run": g.get("d_run", np.nan), "d_bp": g.get("d_bp", np.nan),
                       "run_len": g.get("run_len", np.nan)}
            else:
                rec = {"cls": "B", "d_run": np.nan, "d_bp": m50.find_bp(rev, 11), "run_len": np.nan}
            rows.append({"feature_id": r.feature_id, "acc_T": r.acc_T, **rec})
    g = pd.DataFrame(rows)
    g.to_parquet(LOCAL / "code_geometry.parquet", index=False)
    for c, v in g.groupby("cls"):
        log(f"geometry {c}: {len(v):,} exons, branch point found for {v.d_bp.notna().mean():.3f}, "
            f"median d_bp {v.d_bp.median():.0f}")


def arms(acc_t: np.ndarray, d_bp: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(R_down, R_up, B) masks of the rule."""
    down = (acc_t >= 10) & (acc_t <= 19)
    up = (acc_t >= 20) & (np.nan_to_num(d_bp, nan=-1) >= UP_CUT)
    return down, up, acc_t <= 5


def score_g(P, acc_t, d_bp, dec) -> np.ndarray:
    down, up, b = arms(acc_t, d_bp)
    return m42.set_score(P, down, b, dec) - m42.set_score(P, up, b, dec)


# ---------------------------------------------------------------- shared statistics

def resid_rho(y, x, bins, covs) -> float:
    """Spearman of rank residuals on bin dummies and covariate ranks (D-63's statistic)."""
    D = pd.get_dummies(bins).to_numpy(float)
    Z = np.column_stack([D] + [rankdata(c) for c in covs])
    Q, _ = np.linalg.qr(Z)
    ry, rx = rankdata(y), rankdata(x)
    return float(np.corrcoef(ry - Q @ (Q.T @ ry), rx - Q @ (Q.T @ rx))[0, 1])


def boot_p(d: np.ndarray) -> float:
    return float(max(2 * min((d <= 0).mean(), (d >= 0).mean()), 1 / len(d)))


def auc_lower(s: np.ndarray, pos: np.ndarray) -> float:
    """P(score of a positive < score of a negative), ties half, by ranks."""
    r = rankdata(s)
    n1, n0 = pos.sum(), (~pos).sum()
    return float(1 - (r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def predictions(pred: list[dict], sub: str) -> pd.DataFrame:
    pr = m50.add_verdicts(pred)
    (OUT / sub).mkdir(parents=True, exist_ok=True)
    pr.to_csv(OUT / sub / "predictions.csv", index=False)
    log(f"{sub} predictions:\n" + pr[["id", "effect", "p", "q", "verdict", "n"]].to_string(index=False))
    return pr


def delta_adj(P, a, b, dec, is_b) -> np.ndarray:
    """D-52's per-exon shift: mean percentile of group a minus group b, minus the median B-exon shift
    of the same mean-PSI decile."""
    d = np.nanmean(P[a], axis=0) - np.nanmean(P[b], axis=0)
    bmed = pd.Series(d[is_b]).groupby(dec[is_b]).median()
    return d - pd.Series(dec).map(bmed).to_numpy()


def cohort_matrix(cohort: str):
    """D-52's universe, percentiles and deciles for a cohort, with the geometry joined."""
    tab, X, feats = m42.cohort_data(cohort)
    grp = tab["group"].to_numpy()
    ok = m42.universe(X, feats, grp)
    X, f = X[:, ok], feats[ok].reset_index(drop=True)
    P = pd.DataFrame(X).rank(pct=True, axis=0).to_numpy()
    mean_psi = np.nanmean(X, axis=0)
    dec = pd.qcut(mean_psi, 10, labels=False, duplicates="drop")
    geo = pd.read_parquet(LOCAL / "code_geometry.parquet").set_index("feature_id")
    f["d_bp"] = f.feature_id.map(geo.d_bp).to_numpy()
    f["d_run"] = f.feature_id.map(geo.d_run).to_numpy()
    f["mean_psi"] = mean_psi
    return tab, grp, X, P, f, dec


# ---------------------------------------------------------------- D-65 specificity

def specificity() -> None:
    geo = pd.read_parquet(LOCAL / "code_geometry.parquet")
    e = m50.effects().merge(geo[["feature_id", "d_bp", "d_run"]], on="feature_id", how="left")
    mean = e.groupby("feature_id").agg(delta_adj=("delta_adj", "mean"), mean_psi=("mean_psi", "mean"),
                                       acc_T=("acc_T", "first"), d_bp=("d_bp", "first"), d_run=("d_run", "first"),
                                       n=("cohort", "size"))
    m2 = mean[mean.n >= 2]
    R = m2[(m2.acc_T >= 10) & m2.d_bp.notna()]
    B = m2[(m2.acc_T <= 5) & m2.d_bp.notna()]
    rng = np.random.default_rng(SEED)
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    gap = (R.d_bp - R.d_run).to_numpy(float)
    r, p, n = m50.adjusted_rho(R.delta_adj.to_numpy(), gap, m50.run_bin(R.acc_T.to_numpy()),
                               [R.mean_psi.to_numpy(), R.d_run.to_numpy(float)], rng)
    add("BS1", r, 1, p, f"{n} R exons in >= 2 cohorts", "Spearman(mean delta_adj, d_bp - d_run | run bins, mean PSI, d_run)")
    # BS2: R against B exons
    stat_r = lambda d: resid_rho(d.delta_adj.to_numpy(), d.d_bp.to_numpy(), m50.run_bin(d.acc_T.to_numpy()),
                                 [d.mean_psi.to_numpy()])
    stat_b = lambda d: resid_rho(d.delta_adj.to_numpy(), d.d_bp.to_numpy(), np.zeros(len(d)), [d.mean_psi.to_numpy()])
    rr, rb = stat_r(R), stat_b(B)
    bs = np.array([stat_r(R.iloc[rng.integers(0, len(R), len(R))]) - stat_b(B.iloc[rng.integers(0, len(B), len(B))])
                   for _ in range(BOOT)])
    add("BS2", rr - rb, 1, boot_p(bs), f"R {len(R)} (rho {rr:.3f}), B {len(B):,} (rho {rb:.3f})",
        "adjusted rho(mean delta_adj, d_bp): R exons minus B exons")
    # BS3, BS4: mismatch-repair-proficient contrasts
    ex = []
    for pid, cohort, (ga, gb), lab in [("BS3", "UCEC_TCGA", ("POLE", "rest"), "POLE vs rest"),
                                       ("BS4", "STAD_TCGA", ("CIN", "GS"), "CIN vs GS")]:
        tab, grp, X, P, f, dec = cohort_matrix(cohort)
        lab_col = tab["label"].to_numpy() if cohort == "STAD_TCGA" else grp
        is_b = f.acc_T.to_numpy() <= 5
        d_msi = delta_adj(P, grp == "MSI", grp == "rest", dec, is_b)
        d_ctl = delta_adj(P, lab_col == ga, lab_col == gb, dec, is_b)
        ref = e[e.cohort.eq(cohort)].set_index("feature_id").delta_adj.reindex(f.feature_id).to_numpy()
        log(f"{cohort}: MSI delta_adj reproduces D-52 to {np.nanmax(np.abs(ref - d_msi)):.2e}")
        fin = np.isfinite(X)
        keep = ((f.acc_T >= 10) & f.d_bp.notna()).to_numpy() & (fin[lab_col == ga].sum(0) >= 10) & \
            (fin[lab_col == gb].sum(0) >= 10)
        s = f[keep].assign(d_msi=d_msi[keep], d_ctl=d_ctl[keep]).reset_index(drop=True)
        rho = lambda d, col: resid_rho(d[col].to_numpy(), d.d_bp.to_numpy(), m50.run_bin(d.acc_T.to_numpy()),
                                       [d.mean_psi.to_numpy()])
        r_m, r_c = rho(s, "d_msi"), rho(s, "d_ctl")
        bs = []
        for _ in range(BOOT):
            b = s.iloc[rng.integers(0, len(s), len(s))]
            bs.append(rho(b, "d_msi") - rho(b, "d_ctl"))
        add(pid, r_m - r_c, 1, boot_p(np.asarray(bs)),
            f"{len(s)} R exons; MSI rho {r_m:.3f}, {lab} rho {r_c:.3f}; groups {ga} {(lab_col == ga).sum()}, "
            f"{gb} {(lab_col == gb).sum()}", f"{cohort}: adjusted rho(MSI shift, d_bp) minus adjusted rho({lab} shift, d_bp)")
        for col, name in [("d_msi", "MSI vs rest"), ("d_ctl", lab)]:
            for cls, k in [("10-19", s.acc_T < 20), (">=20", s.acc_T >= 20)]:
                t = s[k]
                tert = pd.qcut(t.d_bp.rank(method="first"), 3, labels=["near", "middle", "far"])
                for tl, v in t.groupby(tert, observed=True):
                    ex.append({"what": "share up by d_bp tertile", "cohort": cohort, "contrast": name, "run_class": cls,
                               "bin": tl, "n": len(v), "share_up": float((v[col] > 0).mean()),
                               "median_d_bp": float(v.d_bp.median())})
    predictions(pred, "specificity")
    for cls, k in [("10-19", R.acc_T < 20), (">=20", R.acc_T >= 20)]:
        t = R[k]
        for lo in range(10, 100, 10):
            v = t[(t.d_bp >= lo) & (t.d_bp < lo + 10)]
            if len(v) >= 10:
                ex.append({"what": "share up by 10-nt d_bp bin, mean over >= 2 cohorts", "run_class": cls,
                           "bin": f"{lo}-{lo + 9}", "n": len(v), "share_up": float((v.delta_adj > 0).mean()),
                           "median_d_bp": float(v.d_bp.median())})
    gb = geo[geo.cls.eq("B")]
    ex.append({"what": "B-exon branch points", "n": len(gb), "share_up": np.nan,
               "median_d_bp": float(gb.d_bp.median()), "bin": f"found for {gb.d_bp.notna().mean():.3f}"})
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "specificity" / "exploratory.csv", index=False, float_format="%.4g")
    log(f"specificity exploratory:\n{ex.round(3).to_string(index=False)}")


# ---------------------------------------------------------------- D-66 thirty held-out cancers

def mc3_counts(samples: set[str]) -> pd.DataFrame:
    """Per 15-character tumour barcode (first aliquot in sort order): calls passing the D-61 flags,
    homopolymer deletions (D-61) and SNVs."""
    cols = ["Tumor_Sample_Barcode", "Variant_Type", "Reference_Allele", "CONTEXT", "FILTER"]
    m = pd.read_csv(m50.EXT / "gdc_pancanatlas" / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=cols, dtype=str)
    m["sample15"] = m.Tumor_Sample_Barcode.str[:15]
    m = m[m.sample15.isin(samples)]
    first = m.sort_values("Tumor_Sample_Barcode").drop_duplicates("sample15").Tumor_Sample_Barcode
    m = m[m.Tumor_Sample_Barcode.isin(first)]
    m = m[m.FILTER.str.split(",").map(lambda f: set(f) <= m50.KEEP_FLAGS)]
    d = m[m.Variant_Type.eq("DEL") & m.Reference_Allele.str.fullmatch(r"A+|T+")]
    ok = [c[6:6 + len(r)] == r and m50.context_run(c, r) >= 5 for c, r in zip(d.CONTEXT, d.Reference_Allele)]
    out = pd.DataFrame({"n_calls": m.groupby("sample15").size()})
    out["n_hp"] = d[ok].groupby("sample15").size().reindex(out.index).fillna(0).astype(int)
    out["n_snv"] = m[m.Variant_Type.eq("SNP")].groupby("sample15").size().reindex(out.index).fillna(0).astype(int)
    return out


def pancancer() -> None:
    tab = pd.read_parquet(DATA / "tasks" / "cancer_type__PANCAN.parquet")
    t = tab[tab.role.eq("tumour") & ~tab.project.isin(TRAIN)].sort_values("sample_id").drop_duplicates("case_id").copy()
    t["sample15"] = t.sample_id.str[:15]
    cnt = mc3_counts(set(t.sample15))
    t = t[t.sample15.isin(cnt.index)].join(cnt, on="sample15")
    t["mmr"] = t.n_hp >= MMR_CUT
    t["hyper"] = (t.n_hp < MMR_CUT) & (t.n_snv >= HYPER_SNV)
    log(f"held-out tumours with MC3 calls: {len(t):,} in {t.project.nunique()} projects; MMR-deficient "
        f"{int(t.mmr.sum())}; hypermutated proficient {int(t.hyper.sum())}")
    ann = pd.read_parquet(LOCAL / "se_repeats_pan.parquet").set_index("feature_id")
    geo = pd.read_parquet(LOCAL / "code_geometry.parquet").set_index("feature_id")
    feats = pd.read_parquet(PAN / "se" / "features.parquet")
    a = ann.reindex(feats.feature_id)
    acc_t = a.acc_T.to_numpy(float)
    d_bp = geo.d_bp.reindex(feats.feature_id).to_numpy(float)
    struct = (a.acc_intron_len.to_numpy(float) >= m42.MIN_INTRON["acc"]) & a.acc_ag.eq("AG").to_numpy()
    nf = len(feats)
    acc = {k: np.zeros(nf) for k in ["s_m", "c_m", "s_o", "c_o", "s_psi", "c_psi"]}
    mre = int(np.flatnonzero(feats.feature_id.eq(m42.MRE11_E5))[0]) if feats.feature_id.eq(m42.MRE11_E5).any() else None
    scores, per = [], []
    for proj, s in t.groupby("project"):
        vm, _ = load_view(PAN, "se", s.file_id)
        X = vm.get()
        ok = struct & (np.isfinite(X).mean(axis=0) >= 0.8)
        Xo = X[:, ok]
        P = pd.DataFrame(Xo).rank(pct=True, axis=0).to_numpy()
        mpsi = np.nanmean(Xo, axis=0)
        dec = pd.qcut(mpsi, 10, labels=False, duplicates="drop")
        at, db = acc_t[ok], d_bp[ok]
        A = m42.set_score(P, at >= 10, at <= 5, dec)
        G = score_g(P, at, db, dec)
        mm = s.mmr.to_numpy()
        down, up_, b_ = arms(at, db)
        sc = s[["file_id", "case_id", "project", "n_hp", "n_snv", "mmr", "hyper"]].assign(
            score_A=A, score_G=G, score_down=m42.set_score(P, down, b_, dec), score_up=m42.set_score(P, up_, b_, dec))
        for col in ["A", "G", "down", "up"]:
            sc[f"pct_{col}"] = rankdata(sc[f"score_{col}"]) / len(sc)
        scores.append(sc)
        idx = np.flatnonzero(ok)
        fin = np.isfinite(P)
        acc["s_m"][idx] += np.nansum(P[mm], axis=0)
        acc["c_m"][idx] += fin[mm].sum(axis=0)
        acc["s_o"][idx] += np.nansum(P[~mm], axis=0)
        acc["c_o"][idx] += fin[~mm].sum(axis=0)
        acc["s_psi"][idx] += np.nansum(Xo, axis=0)
        acc["c_psi"][idx] += fin.sum(axis=0)
        row = {"project": proj, "tumours": len(s), "mmr_deficient": int(mm.sum()),
               "hypermutated_proficient": int(s.hyper.sum()), "exons": int(ok.sum()),
               "R_down": int(arms(at, db)[0].sum()), "R_up": int(arms(at, db)[1].sum())}
        if mm.sum() >= 3:
            row["auroc_A"], row["auroc_G"] = auc_lower(A, mm), auc_lower(G, mm)
        if mre is not None:
            v = X[:, mre]
            for g_, k in [("mmr", mm), ("other", ~mm)]:
                w = v[k & np.isfinite(v)]
                row[f"mre11_e5_median_{g_}"] = float(np.median(w)) if len(w) >= 5 else np.nan
        per.append(row)
        log(f"{proj}: {len(s)} tumours, {int(mm.sum())} MMR-deficient, {int(ok.sum()):,} exons")
    sc = pd.concat(scores, ignore_index=True)
    sc["purity"] = sc.case_id.map(m50.purity_by_patient())
    sc.to_parquet(LOCAL / "pancancer_scores.parquet", index=False)
    per = pd.DataFrame(per)
    (OUT / "pancancer").mkdir(parents=True, exist_ok=True)
    per.to_csv(OUT / "pancancer" / "by_cancer.csv", index=False, float_format="%.4g")
    log(f"by cancer:\n{per.round(3).to_string(index=False)}")
    # tumour-level tests
    rng = np.random.default_rng(SEED)
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    mm = sc.mmr.to_numpy()
    for pid, col in [("PC1", "pct_A"), ("PC2", "pct_G")]:
        au, p = m42.auroc_lower(sc.loc[mm, col].to_numpy(), sc.loc[~mm, col].to_numpy())
        add(pid, au - 0.5, 1, p, f"MMR-deficient {int(mm.sum())}, other {int((~mm).sum())}",
            f"pooled AUROC(score {col[-1]} lower in MMR-deficient) - 0.5")
    strata = [np.asarray(v) for v in sc.groupby(["project", "mmr"]).indices.values()]
    gA, gG = sc.pct_A.to_numpy(), sc.pct_G.to_numpy()
    bs = []
    for _ in range(BOOT):
        i = np.concatenate([rng.choice(g, len(g)) for g in strata])
        bs.append(auc_lower(gG[i], mm[i]) - auc_lower(gA[i], mm[i]))
    obs = auc_lower(gG, mm) - auc_lower(gA, mm)
    add("PC3", obs, 1, boot_p(np.asarray(bs)), f"AUROC G {auc_lower(gG, mm):.3f}, A {auc_lower(gA, mm):.3f}",
        "AUROC(G) - AUROC(A), pooled")
    k = ~sc.project.isin(RELATIVES).to_numpy()
    au, p = m42.auroc_lower(gG[k & mm], gG[k & ~mm])
    add("PC4", au - 0.5, 1, p, f"MMR-deficient {int((k & mm).sum())}, other {int((k & ~mm).sum())}",
        "as PC2 without READ, ESCA and UCS")
    hy = sc.hyper.to_numpy()
    if min(mm.sum(), hy.sum()) >= 10:
        au, p = m42.auroc_lower(gG[mm], gG[hy])
        add("PC5", au - 0.5, 1, p, f"MMR-deficient {int(mm.sum())}, hypermutated proficient {int(hy.sum())}",
            "AUROC(G lower in MMR-deficient than hypermutated proficient) - 0.5")
    else:
        pred.append({"id": "PC5", "what": "not testable", "n": f"MMR-deficient {int(mm.sum())}, hyper {int(hy.sum())}"})
    # exon-level tests
    el = feats[["feature_id", "symbol"]].assign(acc_T=acc_t, d_bp=d_bp, c_m=acc["c_m"], c_o=acc["c_o"])
    keep = (acc["c_m"] >= 10) & (acc["c_o"] >= 100)
    with np.errstate(invalid="ignore", divide="ignore"):
        el["delta"] = acc["s_m"] / acc["c_m"] - acc["s_o"] / acc["c_o"]
        el["mean_psi"] = acc["s_psi"] / acc["c_psi"]
    el = el[keep].reset_index(drop=True)
    el["dec"] = pd.qcut(el.mean_psi, 10, labels=False, duplicates="drop")
    bmed = el[el.acc_T <= 5].groupby("dec").delta.median()
    el["delta_adj"] = el.delta - el.dec.map(bmed)
    tri = m50.effects().groupby("feature_id").delta_adj.mean().rename("train_delta_adj")
    el = el.join(tri, on="feature_id")
    R = el[el.acc_T >= 10]
    v = R.dropna(subset=["train_delta_adj"])
    r = spearmanr(v.delta_adj, v.train_delta_adj)
    add("PC6", float(r.statistic), 1, float(r.pvalue), f"{len(v)} R exons",
        "Spearman(held-out delta_adj, mean delta_adj of the three cohorts)")
    down, up, _ = arms(R.acc_T.to_numpy(), R.d_bp.to_numpy())
    au, p = m42.auroc_lower(R.delta_adj.to_numpy()[down], R.delta_adj.to_numpy()[up])
    add("PC7", au - 0.5, 1, p, f"R_down {int(down.sum())} (median {np.median(R.delta_adj.to_numpy()[down]):+.4f}), "
        f"R_up {int(up.sum())} (median {np.median(R.delta_adj.to_numpy()[up]):+.4f})",
        "AUROC(held-out delta_adj higher for R_up than R_down) - 0.5")
    Rb = R[R.d_bp.notna()]
    r, p, n = m50.adjusted_rho(Rb.delta_adj.to_numpy(), Rb.d_bp.to_numpy(), m50.run_bin(Rb.acc_T.to_numpy()),
                               [Rb.mean_psi.to_numpy()], rng)
    add("PC8", r, 1, p, f"{n} R exons", "adjusted rho(held-out delta_adj, d_bp)")
    predictions(pred, "pancancer")
    R.drop(columns=["dec"]).to_csv(OUT / "pancancer" / "exon_shifts_R.csv", index=False, float_format="%.5g")
    # exploratory
    ex = []
    m_ = sc[sc.mmr]
    r, p, n = m42.partial_spearman(m_.pct_G.to_numpy(), np.log1p(m_.n_hp.to_numpy(float)), [m_.purity.to_numpy()])
    ex.append({"what": "dose: partial Spearman(G percentile, log(1 + homopolymer deletions) | purity), MMR-deficient",
               "value": r, "p": p, "n": n})
    for name, k in [("all held-out", np.ones(len(sc), bool)), ("without READ, ESCA, UCS", k)]:
        for col in ["pct_A", "pct_G"]:
            ex.append({"what": f"AUROC {col[-1]}, {name}", "value": auc_lower(sc[col].to_numpy()[k], mm[k]),
                       "n": int((mm & k).sum())})
    for col, lab in [("pct_down", "R_down arm against B, lower in MMR-deficient"),
                     ("pct_up", "R_up arm against B, lower in MMR-deficient (below 0.5 = higher)")]:
        ex.append({"what": f"A-15: AUROC of the {lab}, all held-out", "value": auc_lower(sc[col].to_numpy(), mm),
                   "n": int(mm.sum())})
    for lab, s_ in [("R_down", down), ("R_up", up)]:
        ex.append({"what": f"share of {lab} exons more included in held-out MMR-deficient tumours",
                   "value": float((R.delta_adj.to_numpy()[s_] > 0).mean()), "n": int(s_.sum())})
    nm = tab[tab.role.eq("normal")]
    ex.append({"what": "MMR-deficient patients with an adjacent normal in the matrix",
               "value": int(nm.case_id.isin(sc.loc[sc.mmr, "case_id"]).sum()), "n": int(sc.mmr.sum())})
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "pancancer" / "exploratory.csv", index=False, float_format="%.4g")
    log(f"pancancer exploratory:\n{ex.round(4).to_string(index=False)}")


# ---------------------------------------------------------------- D-67 endometrium

def endometrium() -> None:
    pred, ex, local = [], [], []
    for cohort in ["UCEC_TCGA", "STAD_TCGA", "COAD_TCGA"]:
        tab, grp, X, P, f, dec = cohort_matrix(cohort)
        at, db = f.acc_T.to_numpy(float), f.d_bp.to_numpy(float)
        A = m42.set_score(P, at >= 10, at <= 5, dec)
        G = score_g(P, at, db, dec)
        ref = pd.read_parquet(LOCAL / "sample_scores.parquet").set_index("file_id").score_acc.reindex(tab.file_id)
        log(f"{cohort}: score A reproduces D-52 to {np.nanmax(np.abs(ref.to_numpy() - A)):.2e}; "
            f"R_down {int(arms(at, db)[0].sum())}, R_up {int(arms(at, db)[1].sum())}")
        msi, rest = grp == "MSI", grp == "rest"
        au_g, p_g = m42.auroc_lower(G[msi], G[rest])
        au_a, _ = m42.auroc_lower(A[msi], A[rest])
        if cohort == "UCEC_TCGA":
            pred.append({"id": "EG1", "what": "UCEC_TCGA: AUROC(G lower in MSI than the rest) - 0.5", "effect": au_g - 0.5,
                         "predicted_sign": 1, "p": p_g, "n": f"MSI {msi.sum()}, rest {rest.sum()}"})
            lo, hi, p = m42.boot_auroc_diff(G, A, grp, SEED, BOOT)
            pred.append({"id": "EG2", "what": "UCEC_TCGA: AUROC(G) - AUROC(A)", "effect": au_g - au_a, "predicted_sign": 1,
                         "p": p, "n": f"AUROC G {au_g:.3f}, A {au_a:.3f}; 95% CI of the difference {lo:.3f} to {hi:.3f}"})
            pole = grp == "POLE"
            ex.append({"cohort": cohort, "what": "AUROC(G lower in MSI than POLE)",
                       "value": m42.auroc_lower(G[msi], G[pole])[0], "n": f"MSI {msi.sum()}, POLE {pole.sum()}"})
        ex.append({"cohort": cohort, "what": "AUROC(score lower in MSI than the rest), A", "value": au_a,
                   "n": f"MSI {msi.sum()}, rest {rest.sum()}"})
        ex.append({"cohort": cohort, "what": "AUROC(score lower in MSI than the rest), G", "value": au_g, "p": p_g,
                   "n": f"MSI {msi.sum()}, rest {rest.sum()}"})
        for arm, k in zip(["R_down", "R_up"], arms(at, db)[:2]):
            s_ = m42.set_score(P, k, at <= 5, dec)
            ex.append({"cohort": cohort, "what": f"AUROC(arm {arm} vs B lower in MSI)", "value": m42.auroc_lower(s_[msi], s_[rest])[0],
                       "n": f"{int(k.sum())} exons"})
        local.append(tab[["file_id", "case_id", "group"]].assign(cohort=cohort, score_A=A, score_G=G))
    predictions(pred, "endometrium")
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "endometrium" / "exploratory.csv", index=False, float_format="%.4g")
    pd.concat(local, ignore_index=True).to_parquet(LOCAL / "code_scores_three.parquet", index=False)
    log(f"endometrium exploratory:\n{ex.round(3).to_string(index=False)}")


# ---------------------------------------------------------------- D-68 infiltration

def infiltration() -> None:
    import statsmodels.formula.api as smf

    from ctsm.biology import gex_by_symbol
    from ctsm.config import load_config
    m43 = _load("msi43", "43_msi_robustness.py")
    st, im = m43.BIO["stromal_genes"], m43.BIO["immune_genes"]
    sc = pd.read_parquet(LOCAL / "pancancer_scores.parquet")
    g = gex_by_symbol(load_config(), sc.file_id, st + im)
    sc["stromal"] = g[[x for x in st if x in g]].mean(axis=1).to_numpy()
    sc["immune"] = g[[x for x in im if x in g]].mean(axis=1).to_numpy()
    ri = np.full(len(sc), np.nan)
    for _, s in sc.groupby("project"):
        X = load_view(PAN, "ri", s.file_id)[0].get()
        ri[s.index.to_numpy()] = np.nanmean(X[:, np.isfinite(X).mean(axis=0) >= 0.8], axis=1)
    sc["ri_mean"] = ri
    for c in ["stromal", "immune", "ri_mean", "purity"]:
        sc[f"{c}_p"] = sc.groupby("project")[c].rank(pct=True)
    sc["mmr_i"] = sc.mmr.astype(int)
    sc.to_parquet(LOCAL / "pancancer_covariates.parquet", index=False)
    pred, ex = [], []

    def add(pid, effect, sign, p, n, what, lo=np.nan, hi=np.nan):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n,
                     "ci_low": lo, "ci_high": hi})

    base = "pct_G ~ mmr_i + C(project)"
    adj = base + " + stromal_p + immune_p + ri_mean_p"
    for pid, formula, d, what in [("IF1", adj, sc, "stromal, immune, ri_mean"),
                                  ("IF2", adj + " + purity_p", sc[sc.purity.notna()], "the same and ABSOLUTE purity")]:
        fit = smf.ols(formula, data=d).fit(cov_type="HC3")
        raw = smf.ols(base, data=d).fit(cov_type="HC3")
        ci = fit.conf_int().loc["mmr_i"]
        add(pid, float(fit.params["mmr_i"]), -1, float(fit.pvalues["mmr_i"]),
            f"{len(d):,} tumours, {int(d.mmr.sum())} MMR-deficient; unadjusted {raw.params['mmr_i']:+.3f}",
            f"held-out OLS coefficient of MMR deficiency on the G percentile, cancer fixed effects, adjusted for {what}",
            ci[0], ci[1])
        ex.append({"what": f"{pid}: unadjusted coefficient (cancer fixed effects only)", "value": float(raw.params["mmr_i"]),
                   "n": len(d)})
        ex.append({"what": f"{pid}: share of the unadjusted coefficient kept", "value": float(fit.params["mmr_i"] / raw.params["mmr_i"]),
                   "n": len(d)})
        for cov in [c for c in ["stromal_p", "immune_p", "ri_mean_p", "purity_p"] if c in formula]:
            ex.append({"what": f"{pid}: coefficient of {cov}", "value": float(fit.params[cov]), "p": float(fit.pvalues[cov]),
                       "n": len(d)})
    low = sc[sc.immune <= sc.groupby("project").immune.transform("median")]
    au, p = m42.auroc_lower(low.loc[low.mmr, "pct_G"].to_numpy(), low.loc[~low.mmr, "pct_G"].to_numpy())
    add("IF3", au - 0.5, 1, p, f"MMR-deficient {int(low.mmr.sum())}, other {int((~low.mmr).sum())}",
        "AUROC(G lower in MMR-deficient) - 0.5, lower-immune half of each cancer")
    for proj, v in low.groupby("project"):
        if v.mmr.sum() >= 3:
            ex.append({"what": f"IF3 by cancer: {proj}", "value": auc_lower(v.pct_G.to_numpy(), v.mmr.to_numpy()),
                       "n": int(v.mmr.sum())})
    three = pd.read_parquet(LOCAL / "code_scores_three.parquet")
    for pid, cohort in [("IF4", "UCEC_TCGA"), ("IF5", "STAD_TCGA"), ("IF6", "COAD_TCGA")]:
        tab = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{cohort}.parquet")
        s = three[three.cohort.eq(cohort)].merge(tab[["file_id", "tss"]], on="file_id")
        s = s.join(m43.covariates(cohort, s), on="file_id")
        s = s[s.group.isin(["MSI", "rest"])].reset_index(drop=True)
        s["msi"] = s.group.eq("MSI").astype(int)
        keep = s.groupby("tss")["tss"].transform("size") > 1
        fit = smf.ols("score_G ~ msi + C(tss) + stromal + immune + ri_mean", data=s[keep]).fit(cov_type="HC3")
        raw = smf.ols("score_G ~ msi + C(tss)", data=s[keep]).fit(cov_type="HC3")
        ci = fit.conf_int().loc["msi"]
        add(pid, float(fit.params["msi"]), -1, float(fit.pvalues["msi"]),
            f"{int(keep.sum())} patients ({int((~keep).sum())} from single-patient sites left out); site only {raw.params['msi']:+.4f}",
            f"{cohort}: OLS coefficient of MSI on score G, adjusted for site, stromal, immune, ri_mean (HC3)", ci[0], ci[1])
        ex.append({"what": f"{pid}: share of the site-only coefficient kept", "value": float(fit.params["msi"] / raw.params["msi"]),
                   "n": int(keep.sum())})
        lw = s[s.immune <= s.immune.median()]
        ex.append({"what": f"{pid}: AUROC(G lower in MSI), lower-immune half", "n": int(lw.msi.sum()),
                   "value": m42.auroc_lower(lw.loc[lw.msi.eq(1), "score_G"].to_numpy(), lw.loc[lw.msi.eq(0), "score_G"].to_numpy())[0]})
    predictions(pred, "infiltration")
    ex.append({"what": "AUROC of the immune score itself, higher in MMR-deficient (held-out)",
               "value": 1 - auc_lower(sc.immune_p.to_numpy(), sc.mmr.to_numpy()), "n": int(sc.mmr.sum())})
    o = sc[~sc.mmr]
    for col in ["pct_G", "pct_down", "pct_up"]:
        r = spearmanr(o[col], o.immune_p)
        ex.append({"what": f"Spearman({col}, immune percentile), non-MMR-deficient held-out tumours",
                   "value": float(r.statistic), "p": float(r.pvalue), "n": len(o)})
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "infiltration" / "exploratory.csv", index=False, float_format="%.4g")
    log(f"infiltration exploratory:\n{ex.round(4).to_string(index=False)}")


# ---------------------------------------------------------------- D-69 MSI calls not from MC3

MSI_TSV = DATA / "external" / "cbioportal_msi" / "tcga_pancan_msi_scores.tsv"
CALLS = {"mantis": (0.4, 0.6), "msisensor": (4.0, 10.0)}   # MSS below the first, MSI-H above the second


def msi_call(x: pd.Series, col: str) -> pd.Series:
    lo, hi = CALLS[col]
    return pd.Series(np.select([x > hi, x < lo], ["MSI-H", "MSS"], np.where(x.notna(), "between", "no score")),
                     index=x.index)


def msicalls() -> None:
    scores = pd.read_csv(MSI_TSV, sep="\t").drop_duplicates("sample_id").set_index("sample_id")[list(CALLS)]
    tab = pd.read_parquet(DATA / "tasks" / "cancer_type__PANCAN.parquet")[["file_id", "sample_id"]]
    sc = pd.read_parquet(LOCAL / "pancancer_scores.parquet").merge(tab, on="file_id")
    sc = sc.join(scores, on=sc.sample_id.str[:15])
    for col in CALLS:
        sc[f"call_{col}"] = msi_call(sc[col], col)
    log(f"held-out tumours {len(sc):,}; with a MANTIS score {int(sc.mantis.notna().sum()):,}, "
        f"MSIsensor {int(sc.msisensor.notna().sum()):,}")
    pred, ex = [], []
    for pid, col in [("MC1", "mantis"), ("MC2", "msisensor")]:
        c = sc[f"call_{col}"]
        hi, ms = sc.loc[c.eq("MSI-H"), "pct_G"].to_numpy(), sc.loc[c.eq("MSS"), "pct_G"].to_numpy()
        au, p = m42.auroc_lower(hi, ms) if len(hi) and len(ms) else (np.nan, np.nan)
        pred.append({"id": pid, "what": f"pooled held-out AUROC(G percentile lower in {col} MSI-H than MSS) - 0.5",
                     "effect": au - 0.5, "predicted_sign": 1, "p": p if len(hi) >= 10 else np.nan,
                     "n": f"MSI-H {len(hi)}, MSS {len(ms)}" + ("" if len(hi) >= 10 else "; fewer than 10 MSI-H")})
    k = sc.mantis.notna()
    mp = sc[k].groupby("project").mantis.rank(pct=True)
    r = spearmanr(sc.loc[k, "pct_G"], mp)
    pred.append({"id": "MC3", "what": "Spearman(G percentile, MANTIS score as a percentile within cancer), held-out",
                 "effect": float(r.statistic), "predicted_sign": -1, "p": float(r.pvalue), "n": f"{int(k.sum()):,} tumours"})
    predictions(pred, "msicalls")
    # exploratory: agreement with D-66's MC3 call, and score G where the calls disagree
    for col in CALLS:
        for mc3, d in sc.groupby("mmr"):
            for call, v in d.groupby(f"call_{col}"):
                ex.append({"what": f"{col} {call}, MC3 {'MMR-deficient' if mc3 else 'proficient'}", "n": len(v),
                           "value": float(v.pct_G.median()) if len(v) >= 5 else np.nan,
                           "value_is": "median G percentile (blank below 5 tumours)"})
    three = pd.read_parquet(LOCAL / "code_scores_three.parquet")
    for cohort in ["UCEC_TCGA", "STAD_TCGA", "COAD_TCGA"]:
        t = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{cohort}.parquet")[["file_id", "sample_id"]]
        s = three[three.cohort.eq(cohort)].merge(t, on="file_id")
        s = s.join(scores, on=s.sample_id.str[:15])
        for col in CALLS:
            c = msi_call(s[col], col)
            hi, ms = s.loc[c.eq("MSI-H"), "score_G"].to_numpy(), s.loc[c.eq("MSS"), "score_G"].to_numpy()
            ex.append({"what": f"{cohort} (in-sample for the rule): AUROC(G lower in {col} MSI-H than MSS)",
                       "value": m42.auroc_lower(hi, ms)[0], "n": len(hi), "value_is": f"AUROC; MSS {len(ms)}"})
            ex.append({"what": f"{cohort}: {col} MSI-H among the subtype's MSI class", "n": int(s.group.eq("MSI").sum()),
                       "value": int((c.eq("MSI-H") & s.group.eq("MSI")).sum()), "value_is": "count"})
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "msicalls" / "exploratory.csv", index=False, float_format="%.4g")
    sc[["file_id", "project", "mmr", "pct_G", "mantis", "msisensor", "call_mantis", "call_msisensor"]].to_parquet(
        LOCAL / "pancancer_msicalls.parquet", index=False)
    log(f"msicalls exploratory:\n{ex.round(4).to_string(index=False)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", required=True,
                    choices=["annotate", "geometry", "specificity", "pancancer", "endometrium", "infiltration",
                             "msicalls"])
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {"annotate": annotate, "geometry": geometry, "specificity": specificity, "pancancer": pancancer,
     "endometrium": endometrium, "infiltration": infiltration, "msicalls": msicalls}[a.only]()


if __name__ == "__main__":
    main()
