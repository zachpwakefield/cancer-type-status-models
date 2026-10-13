"""Step 52 (D-70): is the rising arm of the MSI sequence code real?

The code (D-65): exons after a T-run of 10-19 nt are skipped in mismatch-repair-deficient tumours
(R_down); exons after a run >= 20 nt whose branch point lies >= 36 nt from the AG gain inclusion
(R_up). D-70 asks whether the rise is specific to mismatch repair, survives a tumour-wide inclusion
axis, needs a tract between the branch point and the AG, follows the AG exclusion zone and the
deletion load, and whether PTB relief explains it instead (config/followup.yaml: msi_code_reversal).

  annotate  label-blind: AG exclusion zone and the 'shield' class (a T-run >= 20 nt upstream of a
            branch point) for every exon of the pan-cancer SE matrix
  test      RV1-RV8 and the exploratory tables
  posthoc   A-16, exploratory and written after RV1-RV8 were read: branch-point distance within long
            runs, run length, exon length, shield exons against 0, and regulator levels

Outputs: data.nosync/ctsm_v1/msi/reversal_*.parquet (per-tumour and per-exon values, kept local)
         results/subtype/msi_cis/code/reversal/ (aggregates)
Usage: python scripts/52_msi_reversal.py --only annotate|test|posthoc
"""
from __future__ import annotations

import argparse
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import binomtest, rankdata, spearmanr


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m42 = _load("msi42", "42_msi_repeats.py")
m50 = _load("msi50", "50_msi_followup.py")
m51 = _load("msi51", "51_msi_code.py")
from ctsm.matrix_io import load_view  # noqa: E402

DATA, LOCAL, PAN = m51.DATA, m51.LOCAL, m51.PAN
OUT = m51.OUT / "reversal"
READ, AGEZ_MIN = 300, 12                 # upstream bases read; AGs closer than 12 nt are ignored
SHIELD_RUN, SHIELD_START, SHIELD_END, BP_MIN = 20, 25, 200, 18
DRAWS, PER_GROUP, MIN_FINITE, MIN_GROUP = 50, 20, 10, 20
TRAIN_COHORTS = ["COAD_TCGA", "STAD_TCGA", "UCEC_TCGA"]


# ---------------------------------------------------------------- label-blind annotation

def agez(rev: str) -> tuple[int, bool]:
    """Distance from the exon to the G of the first upstream AG at >= AGEZ_MIN nt (rev: index 0 =
    position 1); (length read, True) when there is none."""
    for p in range(AGEZ_MIN, len(rev)):
        if rev[p - 1] == "G" and rev[p] == "A":
            return p, False
    return len(rev), True


def shield(rev: str) -> dict:
    """First T-run >= 20 nt starting >= 25 and ending <= 200 nt from the exon with a yUnAy branch-
    point A between it and the AG (at 18 nt or more, >= 2 nt downstream of the run)."""
    for m in re.finditer(rf"T{{{SHIELD_RUN},}}", rev):
        s, e = m.start() + 1, m.end()
        if s >= SHIELD_START and e <= SHIELD_END:
            bp = m50.find_bp(rev, BP_MIN, s - 2)
            if np.isfinite(bp):
                return {"sh_start": s, "sh_len": e - s + 1, "sh_bp": bp}
    return {"sh_start": np.nan, "sh_len": np.nan, "sh_bp": np.nan}


def annotate() -> None:
    f = m42.parse_features(pd.read_parquet(PAN / "se" / "features.parquet"))
    rows = []
    for chrom, seq in m42.read_genome(set(f.chrom)):
        for r in f[f.chrom.eq(chrom)].itertuples(index=False):
            if r.strand == "+":
                up = seq[max(int(r.left_end), int(r.es) - READ):int(r.es)]
            else:
                up = seq[int(r.ee):min(int(r.right_start), int(r.ee) + READ)][::-1].translate(m50.COMP)
            rev = up[::-1]
            z, cens = agez(rev)
            rows.append({"feature_id": r.feature_id, "exon_len": int(r.ee - r.es), "agez": z,
                         "agez_censored": cens, **shield(rev)})
    g = pd.DataFrame(rows)
    g.to_parquet(LOCAL / "reversal_geometry.parquet", index=False)
    log(f"annotate: {len(g):,} exons; AGEZ censored {g.agez_censored.mean():.3f}; shield runs {g.sh_start.notna().sum()}")


# ---------------------------------------------------------------- shared pieces

def pan_classes() -> pd.DataFrame:
    """Per pan-cancer SE feature: run, branch point, AGEZ, shield, structure and the D-70 classes."""
    feats = pd.read_parquet(PAN / "se" / "features.parquet")[["feature_id", "symbol"]]
    ann = pd.read_parquet(LOCAL / "se_repeats_pan.parquet").set_index("feature_id")
    geo = pd.read_parquet(LOCAL / "code_geometry.parquet").set_index("feature_id")
    rg = pd.read_parquet(LOCAL / "reversal_geometry.parquet").set_index("feature_id")
    a = ann.reindex(feats.feature_id)
    c = feats.assign(acc_T=a.acc_T.to_numpy(float), d_bp=geo.d_bp.reindex(feats.feature_id).to_numpy(float),
                     struct=((a.acc_intron_len.to_numpy(float) >= m42.MIN_INTRON["acc"]) & a.acc_ag.eq("AG").to_numpy()))
    c = c.join(rg, on="feature_id")
    down, up, b = m51.arms(c.acc_T.to_numpy(), c.d_bp.to_numpy())
    c["cls"] = np.select([up, down, b & (c.d_bp >= m51.UP_CUT).to_numpy(), b, c.sh_start.notna() & (c.acc_T < 10)],
                         ["R_up", "R_down", "B_far", "B_near", "shield"], "other")
    c["is_b"] = b
    return c.reset_index(drop=True)


def delta_adj(s_a, c_a, s_b, c_b, dec, is_b, min_a=MIN_FINITE, min_b=MIN_FINITE) -> np.ndarray:
    """Per exon: mean percentile of group a minus group b (NaN under the finite minimum), minus the
    median B-exon shift of the same mean-PSI decile."""
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.where((c_a >= min_a) & (c_b >= min_b), s_a / np.maximum(c_a, 1) - s_b / np.maximum(c_b, 1), np.nan)
    ok = np.isfinite(d) & is_b & np.isfinite(dec)
    bmed = pd.Series(d[ok]).groupby(dec[ok]).median()
    return d - pd.Series(dec).map(bmed).to_numpy()


def percentiles(X: np.ndarray) -> np.ndarray:
    return pd.DataFrame(X).rank(pct=True, axis=0).to_numpy()


def draw_stat(P, ia_list, ib_list, dec, is_b, r_up) -> float:
    """RV1 statistic: mean over draws of the mean delta_adj of R_up exons."""
    out = []
    for ia, ib in zip(ia_list, ib_list):
        A, B = P[ia], P[ib]
        fa, fb = np.isfinite(A), np.isfinite(B)
        d = delta_adj(np.nansum(A, 0), fa.sum(0), np.nansum(B, 0), fb.sum(0), dec, is_b)
        out.append(np.nanmean(d[r_up]))
    return float(np.mean(out))


# ---------------------------------------------------------------- RV1 null contrasts

def contrasts(cls: pd.DataFrame, mmr_samples: set[str], rng) -> pd.DataFrame:
    """Every pair of classes of the 22 subtype cohorts (MSI, POLE and MC3 MMR-deficient tumours
    removed), and the three MSI contrasts, each at 20 against 20 tumours over 50 draws."""
    pan = pd.read_parquet(DATA / "tasks" / "cancer_type__PANCAN.parquet").set_index("file_id")
    rows = []
    for f in sorted((DATA / "subtype" / "tasks").glob("subtype__*.parquet")):
        cohort = f.stem.split("__")[1]
        t = pd.read_parquet(f)
        t["sample15"] = t.file_id.map(pan.sample_id).str[:15]
        msi = t.label.isin(["MSI", "POLE"]) | t.sample15.isin(mmr_samples)
        sets = [("null", t[~msi].reset_index(drop=True))]
        if cohort in TRAIN_COHORTS:
            keep = t.label.eq("MSI") | (~t.label.eq("POLE") & ~t.sample15.isin(mmr_samples))
            u = t[keep].reset_index(drop=True)
            u["label"] = np.where(u.label.eq("MSI"), "MSI", "rest")
            sets.append(("msi", u))
        for kind, u in sets:
            counts = u.label.value_counts()
            labs = sorted(counts[counts >= MIN_GROUP].index) if kind == "null" else ["MSI", "rest"]
            if len(labs) < 2:
                continue
            X = load_view(PAN, "se", u.file_id)[0].get()
            ok = cls.struct.to_numpy() & (np.isfinite(X).mean(axis=0) >= 0.8)
            Xo = X[:, ok]
            P = percentiles(Xo)
            dec = pd.qcut(np.nanmean(Xo, axis=0), 10, labels=False, duplicates="drop").astype(float)
            c = cls[ok]
            is_b, r_up = c.is_b.to_numpy(), c.cls.eq("R_up").to_numpy()
            pairs = [(a, b) for i, a in enumerate(labs) for b in labs[i + 1:]] if kind == "null" else [("MSI", "rest")]
            for a, b in pairs:
                ia_all, ib_all = np.flatnonzero(u.label.eq(a)), np.flatnonzero(u.label.eq(b))
                ia = [rng.choice(ia_all, PER_GROUP, replace=False) for _ in range(DRAWS)]
                ib = [rng.choice(ib_all, PER_GROUP, replace=False) for _ in range(DRAWS)]
                s = draw_stat(P, ia, ib, dec, is_b, r_up)
                rows.append({"kind": kind, "cohort": cohort, "a": a, "b": b, "n_a": len(ia_all), "n_b": len(ib_all),
                             "r_up_exons": int(r_up.sum()), "stat": s})
                log(f"contrast {kind} {cohort} {a} vs {b}: {s:+.4f}")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- held-out pass

def held_out(cls: pd.DataFrame, rng) -> tuple[pd.DataFrame, pd.DataFrame, float, np.ndarray]:
    """One pass over the 30 held-out cancers: pooled exon shifts, the global axis, the RV1 held-out
    draws and PTB sensitivity of the long-run exons."""
    from ctsm.biology import gex_by_symbol
    from ctsm.config import load_config
    sc = pd.read_parquet(LOCAL / "pancancer_covariates.parquet").reset_index(drop=True)
    nf = len(cls)
    acc = {k: np.zeros(nf) for k in ["s_m", "c_m", "s_o", "c_o", "s_psi", "c_psi"]}
    # RV1 held-out draws: 20 MMR-deficient tumours, each paired with a proficient one of its cancer
    mm_idx = np.flatnonzero(sc.mmr.to_numpy())
    prof = {p: np.flatnonzero((sc.project.eq(p) & ~sc.mmr).to_numpy()) for p in sc.project.unique()}
    draws = []
    for _ in range(DRAWS):
        a = rng.choice(mm_idx, PER_GROUP, replace=False)
        b, used = [], set()
        for i in a:
            pool = [j for j in prof[sc.project[i]] if j not in used]
            j = int(rng.choice(pool))
            used.add(j)
            b.append(j)
        draws.append((set(a.tolist()), set(b)))
    dk = {k: np.zeros((DRAWS, nf)) for k in ["s_a", "c_a", "s_b", "c_b"]}
    long_run = cls.acc_T.to_numpy() >= 20
    ptb_sum, ptb_w = np.zeros(nf), np.zeros(nf)
    g = np.full(len(sc), np.nan)
    cfg = load_config()
    for proj, s in sc.groupby("project"):
        X = load_view(PAN, "se", s.file_id)[0].get()
        ok = cls.struct.to_numpy() & (np.isfinite(X).mean(axis=0) >= 0.8)
        idx = np.flatnonzero(ok)
        P = percentiles(X[:, ok])
        rows = s.index.to_numpy()
        g[rows] = np.nanmean(P, axis=1)
        mm = s.mmr.to_numpy()
        fin = np.isfinite(P)
        acc["s_m"][idx] += np.nansum(P[mm], axis=0)
        acc["c_m"][idx] += fin[mm].sum(axis=0)
        acc["s_o"][idx] += np.nansum(P[~mm], axis=0)
        acc["c_o"][idx] += fin[~mm].sum(axis=0)
        acc["s_psi"][idx] += np.nansum(X[:, ok], axis=0)
        acc["c_psi"][idx] += np.isfinite(X[:, ok]).sum(axis=0)
        for k, (a, b) in enumerate(draws):
            for side, members in [("a", a), ("b", b)]:
                sel = np.array([r in members for r in rows])
                if sel.any():
                    dk[f"s_{side}"][k, idx] += np.nansum(P[sel], axis=0)
                    dk[f"c_{side}"][k, idx] += fin[sel].sum(axis=0)
        # PTB sensitivity: proficient tumours, partial Spearman(PSI, PTBP1 | global axis)
        ptb = gex_by_symbol(cfg, s.file_id, ["PTBP1"])["PTBP1"].to_numpy()
        gz = g[rows]
        for j in np.flatnonzero(long_run[idx]):
            k_ = fin[:, j] & ~mm & np.isfinite(ptb)
            if k_.sum() < 20:
                continue
            x, y, z = rankdata(P[k_, j]), rankdata(ptb[k_]), rankdata(gz[k_])
            if min(np.ptp(x), np.ptp(y), np.ptp(z)) == 0:      # a constant column has no correlation
                continue
            rxy, rxz, ryz = np.corrcoef(x, y)[0, 1], np.corrcoef(x, z)[0, 1], np.corrcoef(y, z)[0, 1]
            den = np.sqrt(max(1e-12, (1 - rxz ** 2) * (1 - ryz ** 2)))
            ptb_sum[idx[j]] += k_.sum() * (rxy - rxz * ryz) / den
            ptb_w[idx[j]] += k_.sum()
        log(f"{proj}: {len(s)} tumours, {int(mm.sum())} MMR-deficient, {len(idx):,} exons")
    sc["g"] = g
    sc["g_p"] = sc.groupby("project").g.rank(pct=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_psi = acc["s_psi"] / acc["c_psi"]
    dec_all = np.full(nf, np.nan)
    fm = np.isfinite(mean_psi)
    dec_all[fm] = pd.qcut(mean_psi[fm], 10, labels=False, duplicates="drop")
    # pooled exon shift (D-66's rule: >= 10 MMR-deficient and >= 100 other finite values)
    el = cls.assign(c_m=acc["c_m"], c_o=acc["c_o"], mean_psi=mean_psi)
    keep = (acc["c_m"] >= 10) & (acc["c_o"] >= 100)
    el = el[keep].reset_index(drop=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        el["delta"] = acc["s_m"][keep] / acc["c_m"][keep] - acc["s_o"][keep] / acc["c_o"][keep]
    el["dec"] = pd.qcut(el.mean_psi, 10, labels=False, duplicates="drop")
    bmed = el[el.is_b].groupby("dec").delta.median()
    el["delta_adj"] = el.delta - el.dec.map(bmed)
    with np.errstate(invalid="ignore", divide="ignore"):
        el["ptb_sens"] = np.where(ptb_w[keep] > 0, ptb_sum[keep] / ptb_w[keep], np.nan)
    stats = []
    for k in range(DRAWS):
        d = delta_adj(dk["s_a"][k], dk["c_a"][k], dk["s_b"][k], dk["c_b"][k], dec_all, cls.is_b.to_numpy())
        stats.append(np.nanmean(d[cls.cls.eq("R_up").to_numpy()]))
    return sc, el, float(np.mean(stats)), np.asarray(stats)


# ---------------------------------------------------------------- RV6: the three cohorts

def train_up_arm() -> pd.DataFrame:
    """R_up arm percentile within each training cohort (D-52 universe), with MC3 counts and purity."""
    dl = pd.read_parquet(LOCAL / "deletion_load.parquet")[["file_id", "n_hp", "purity"]]
    out = []
    for cohort in TRAIN_COHORTS:
        tab, grp, X, P, f, dec = m51.cohort_matrix(cohort)
        down, up, b = m51.arms(f.acc_T.to_numpy(float), f.d_bp.to_numpy(float))
        s_up = m42.set_score(P, up, b, dec)
        t = tab[["file_id"]].assign(cohort=cohort, group=grp, pct_up=rankdata(s_up) / len(s_up))
        out.append(t.merge(dl, on="file_id", how="left"))
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------- tests

def test() -> None:
    import statsmodels.formula.api as smf
    rng = np.random.default_rng(m51.SEED)
    cls = pan_classes()
    log("classes: " + ", ".join(f"{k} {v:,}" for k, v in cls.cls.value_counts().items()))
    sc, el, s_held, s_draws = held_out(cls, rng)
    pan = pd.read_parquet(DATA / "tasks" / "cancer_type__PANCAN.parquet")
    samples = set(pan[pan.role.eq("tumour")].sample_id.str[:15])
    cnt = m51.mc3_counts(samples)
    mmr_samples = set(cnt.index[cnt.n_hp >= m51.MMR_CUT])
    con = contrasts(cls, mmr_samples, np.random.default_rng(m51.SEED))
    OUT.mkdir(parents=True, exist_ok=True)
    sc[["file_id", "project", "mmr", "hyper", "n_hp", "pct_up", "pct_down", "g", "g_p", "purity"]].to_parquet(
        LOCAL / "reversal_tumours.parquet", index=False)
    el.to_parquet(LOCAL / "reversal_exons.parquet", index=False)
    pred, ex = [], []

    def add(pid, effect, sign, p, n, what, lo=np.nan, hi=np.nan):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n,
                     "ci_low": lo, "ci_high": hi})

    # RV1
    null = con[con.kind.eq("null")].stat.abs().to_numpy()
    p1 = (1 + (null >= s_held).sum()) / (1 + len(null))
    add("RV1", s_held, 1, p1, f"{len(null)} null contrasts, largest |s| {null.max():.4f}, 95th pct {np.quantile(null, 0.95):.4f}",
        "held-out R_up statistic (20 v 20, 50 draws) against every-pair non-MSI subtype contrasts' |statistic|")
    # RV2
    sc["mmr_i"] = sc.mmr.astype(int)
    adj = "pct_up ~ mmr_i + C(project) + g_p + stromal_p + immune_p + ri_mean_p"
    fit = smf.ols(adj, data=sc).fit(cov_type="HC3")
    raw = smf.ols("pct_up ~ mmr_i + C(project)", data=sc).fit(cov_type="HC3")
    ci = fit.conf_int().loc["mmr_i"]
    add("RV2", float(fit.params["mmr_i"]), 1, float(fit.pvalues["mmr_i"]),
        f"{len(sc):,} tumours, {int(sc.mmr.sum())} MMR-deficient; cancer only {raw.params['mmr_i']:+.3f}",
        "OLS coefficient of MMR deficiency on the R_up arm percentile, cancer fixed effects, global axis, stromal, immune, ri_mean (HC3)",
        ci[0], ci[1])
    # RV3 / RV5
    d = el.delta_adj.to_numpy()
    for pid, other in [("RV3", "B_far"), ("RV5", "shield")]:
        a_, b_ = d[el.cls.eq("R_up")], d[el.cls.eq(other)]
        au, p = m42.auroc_lower(b_, a_)
        add(pid, au - 0.5, 1, p, f"R_up {len(a_)} (median {np.median(a_):+.4f}), {other} {len(b_)} (median {np.median(b_):+.4f})",
            f"AUROC(held-out delta_adj higher for R_up than {other}) - 0.5")
    # RV4
    lr = el[(el.acc_T >= 20) & el.agez.notna()]
    bins = np.where(lr.acc_T >= 25, ">=25", "20-24")
    r, p, n = m50.adjusted_rho(lr.delta_adj.to_numpy(), lr.agez.to_numpy(float), bins, [lr.mean_psi.to_numpy()],
                               np.random.default_rng(m51.SEED), perm=1000)
    add("RV4", r, 1, p, f"{n} exons with a run >= 20", "adjusted rho(held-out delta_adj, AGEZ | run bins, mean PSI)")
    # RV6
    tr = train_up_arm()
    tr = tr[tr.group.eq("MSI") & tr.n_hp.notna()]
    pool = pd.concat([sc.loc[sc.mmr, ["pct_up", "n_hp", "purity"]].assign(source="held-out"),
                      tr[["pct_up", "n_hp", "purity"]].assign(source="three cohorts")], ignore_index=True)
    r, p, n = m42.partial_spearman(pool.pct_up.to_numpy(), np.log1p(pool.n_hp.to_numpy(float)), [pool.purity.to_numpy(float)])
    add("RV6", r, 1, p, f"{n} tumours ({int(sc.mmr.sum())} held-out, {len(tr)} in the three cohorts)",
        "partial Spearman(R_up arm percentile, log(1 + homopolymer deletions) | purity)")
    pool.to_parquet(LOCAL / "reversal_dose.parquet", index=False)
    # RV7
    mm, hy = sc.mmr.to_numpy(), sc.hyper.to_numpy()
    au, p = m42.auroc_lower(sc.pct_up.to_numpy()[hy], sc.pct_up.to_numpy()[mm])
    add("RV7", au - 0.5, 1, p, f"MMR-deficient {int(mm.sum())}, hypermutated proficient {int(hy.sum())}",
        "AUROC(R_up arm higher in MMR-deficient than hypermutated proficient) - 0.5")
    # RV8
    v = el[(el.acc_T >= 20) & el.ptb_sens.notna()]
    rr = spearmanr(v.delta_adj, v.ptb_sens)
    add("RV8", float(rr.statistic), -1, float(rr.pvalue), f"{len(v)} exons with a run >= 20",
        "Spearman(held-out delta_adj, PTB sensitivity); negative = PTB relief")
    pr = m51.predictions(pred, "reversal")

    # ---- exploratory
    nullall = np.concatenate([null, -null])
    ex.append({"what": "RV1 z: held-out statistic / SD of the null (both orientations)", "value": s_held / nullall.std(ddof=1),
               "n": len(null)})
    ex.append({"what": "RV1 held-out statistic, 2.5-97.5% over draws", "value": s_held,
               "lo": float(np.quantile(s_draws, 0.025)), "hi": float(np.quantile(s_draws, 0.975)), "n": DRAWS})
    for r_ in con[con.kind.eq("msi")].itertuples():
        ex.append({"what": f"training MSI contrast {r_.cohort} (in-sample): statistic and share of null |s| below it",
                   "value": r_.stat, "p": float((null < r_.stat).mean()), "n": len(null)})
    ex.append({"what": "RV2 without the global axis: coefficient of mmr", "value": float(
        smf.ols("pct_up ~ mmr_i + C(project) + stromal_p + immune_p + ri_mean_p", data=sc).fit(cov_type="HC3").params["mmr_i"]),
        "n": len(sc)})
    ex.append({"what": "RV2: coefficient of the global axis", "value": float(fit.params["g_p"]), "p": float(fit.pvalues["g_p"]),
               "n": len(sc)})
    calls = pd.read_parquet(LOCAL / "pancancer_msicalls.parquet")[["file_id", "call_mantis"]]
    sm = sc.merge(calls, on="file_id", how="left")
    sm = sm[sm.call_mantis.isin(["MSI-H", "MSS"])].copy()
    sm["msih"] = sm.call_mantis.eq("MSI-H").astype(int)
    fm = smf.ols(adj.replace("mmr_i", "msih"), data=sm).fit(cov_type="HC3")
    ex.append({"what": "RV2 with MANTIS MSI-H vs MSS in place of MMR deficiency: coefficient", "value": float(fm.params["msih"]),
               "p": float(fm.pvalues["msih"]), "n": int(sm.msih.sum())})
    au, p = m42.auroc_lower(sc.pct_up.to_numpy()[hy], sm.loc[sm.msih.eq(1), "pct_up"].to_numpy())
    ex.append({"what": "RV7 with MANTIS MSI-H: AUROC(R_up arm higher than hypermutated proficient)", "value": au, "p": p,
               "n": int(sm.msih.sum())})
    for c_, v_ in el.groupby("cls"):
        ex.append({"what": f"held-out share more included, {c_}", "value": float((v_.delta_adj > 0).mean()),
                   "median": float(v_.delta_adj.median()), "n": len(v_)})
    fr = el[(el.exon_len % 3 == 0)]
    au, p = m42.auroc_lower(fr.delta_adj[fr.cls.eq("B_far")].to_numpy(), fr.delta_adj[fr.cls.eq("R_up")].to_numpy())
    ex.append({"what": "RV3 within exons whose length is a multiple of 3: AUROC", "value": au, "p": p,
               "n": int(fr.cls.eq("R_up").sum())})
    for lab, k in [("MMR-deficient", sc.mmr), ("proficient", ~sc.mmr)]:
        r_ = spearmanr(sc.pct_up[k], 1 - sc.pct_down[k])
        ex.append({"what": f"arm coupling: Spearman(R_up arm, minus R_down arm), {lab}", "value": float(r_.statistic),
                   "p": float(r_.pvalue), "n": int(k.sum())})
    # the training cohorts' mean shift (D-52 universes, in-sample for the rule)
    tm = m50.effects().groupby("feature_id").agg(delta_adj=("delta_adj", "mean"), mean_psi=("mean_psi", "mean"),
                                                 n=("cohort", "size"))
    tm = tm[tm.n >= 2].join(cls.set_index("feature_id")[["cls", "acc_T", "agez"]], how="inner")
    for lab, other in [("RV3", "B_far"), ("RV5", "shield")]:
        a_, b_ = tm.delta_adj[tm.cls.eq("R_up")].to_numpy(), tm.delta_adj[tm.cls.eq(other)].to_numpy()
        if min(len(a_), len(b_)) >= 5:
            au, p = m42.auroc_lower(b_, a_)
            ex.append({"what": f"{lab} on the three cohorts' mean shift: AUROC (R_up {len(a_)}, {other} {len(b_)})",
                       "value": au, "p": p, "n": len(b_)})
    lt = tm[(tm.acc_T >= 20) & tm.agez.notna()]
    r, p, n = m50.adjusted_rho(lt.delta_adj.to_numpy(), lt.agez.to_numpy(float), np.where(lt.acc_T >= 25, ">=25", "20-24"),
                               [lt.mean_psi.to_numpy()], np.random.default_rng(m51.SEED), perm=1000)
    ex.append({"what": "RV4 on the three cohorts' mean shift: adjusted rho", "value": r, "p": p, "n": n})
    lr2 = el[(el.acc_T >= 20) & el.d_bp.notna() & el.agez.notna()]
    ex.append({"what": "Spearman(AGEZ, d_bp), held-out exons with a run >= 20", "value": float(spearmanr(lr2.agez, lr2.d_bp).statistic),
               "n": len(lr2)})
    for lab, k in [("R_up", el.cls.eq("R_up")), ("runs >= 20", el.acc_T >= 20)]:
        v_ = el[k & el.ptb_sens.notna()]
        ex.append({"what": f"median PTB sensitivity, {lab}", "value": float(v_.ptb_sens.median()), "n": len(v_)})
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "exploratory.csv", index=False, float_format="%.4g")
    log(f"reversal exploratory:\n{ex.round(4).to_string(index=False)}")
    # aggregates for the figure
    con.assign(kind=con.kind.replace({"null": "other"})).to_csv(OUT / "contrasts.csv", index=False,
                                                                  float_format="%.5g")
    pd.DataFrame({"draw": range(DRAWS), "stat": s_draws}).to_csv(OUT / "held_out_draws.csv", index=False, float_format="%.5g")
    el[el.cls.isin(["R_up", "R_down", "B_far", "shield"]) | (el.acc_T >= 20)][
        ["feature_id", "symbol", "cls", "acc_T", "d_bp", "agez", "agez_censored", "sh_start", "exon_len", "c_m", "c_o",
         "mean_psi", "delta_adj", "ptb_sens"]].to_csv(OUT / "exon_shifts.csv", index=False, float_format="%.5g")
    log(f"reversal predictions done:\n{pr[['id', 'effect', 'p', 'q', 'verdict']].to_string(index=False)}")


# ---------------------------------------------------------------- A-16 (exploratory, after RV1-RV8)

BP_BINS = [(18, 29), (30, 35), (36, 49), (50, 69), (70, 101)]
RUN_BINS = [(10, 11), (12, 13), (14, 15), (16, 19), (20, 24), (25, 29), (30, 99)]
REGULATORS = ["PTBP1", "HNRNPC", "U2AF2", "MKI67"]


def posthoc() -> None:
    """Exploratory checks added after RV1-RV8 were read (amendment A-16). Nothing here is a test of
    D-70; every row is descriptive, with a nominal p where one is given."""
    import statsmodels.formula.api as smf

    from ctsm.biology import gex_by_symbol
    from ctsm.config import load_config
    el = pd.read_parquet(LOCAL / "reversal_exons.parquet")
    rows = []

    def share(what, v, **kw):
        k = int((v.delta_adj > 0).sum())
        rows.append({"what": what, "n": len(v), "share_up": k / len(v) if len(v) else np.nan,
                     "median": float(v.delta_adj.median()) if len(v) else np.nan,
                     "p": binomtest(k, len(v)).pvalue if len(v) else np.nan, **kw})

    lr = el[el.acc_T >= 20]
    for lo, hi in BP_BINS:
        share("runs >= 20 by branch-point distance", lr[(lr.d_bp >= lo) & (lr.d_bp <= hi)], bin=f"{lo}-{hi}")
    share("runs >= 20, no branch point found within 100 nt", lr[lr.d_bp.isna()], bin="none")
    au, p = m42.auroc_lower(lr.delta_adj[lr.d_bp < m51.UP_CUT].to_numpy(), lr.delta_adj[lr.d_bp >= m51.UP_CUT].to_numpy())
    rows.append({"what": "runs >= 20: AUROC(shift higher with the branch point >= 36 nt than < 36 nt)", "value": au, "p": p,
                 "n": int((lr.d_bp < m51.UP_CUT).sum())})
    for lo, hi in RUN_BINS:
        share("all repeat exons by run length", el[(el.acc_T >= lo) & (el.acc_T <= hi)], bin=f"{lo}-{hi}")
    r = spearmanr(lr.acc_T, lr.delta_adj)
    rows.append({"what": "runs >= 20: Spearman(shift, run length)", "value": float(r.statistic), "p": float(r.pvalue), "n": len(lr)})
    share("shield exons against 0", el[el.cls.eq("shield")], bin="shield")
    share("R_up exons against 0", el[el.cls.eq("R_up")], bin="R_up")
    for c in ["R_up", "B_far"]:
        v = el[el.cls.eq(c)]
        r = spearmanr(v.exon_len, v.delta_adj)
        rows.append({"what": f"{c}: Spearman(shift, exon length)", "value": float(r.statistic), "p": float(r.pvalue), "n": len(v)})
    up = el[el.cls.eq("R_up") & el.ptb_sens.notna()]
    share("R_up exons with PTB sensitivity <= 0", up[up.ptb_sens <= 0], bin="ptb<=0")
    share("R_up exons with PTB sensitivity > 0", up[up.ptb_sens > 0], bin="ptb>0")
    # regulator levels: higher in MMR-deficient tumours, and do they take the R_up arm's MMR effect?
    t = pd.read_parquet(LOCAL / "pancancer_covariates.parquet").merge(
        pd.read_parquet(LOCAL / "reversal_tumours.parquet")[["file_id", "g_p"]], on="file_id")
    g = gex_by_symbol(load_config(), t.file_id, REGULATORS)
    for c in REGULATORS:
        t[f"{c}_p"] = pd.Series(g[c].to_numpy(), index=t.index).groupby(t.project).rank(pct=True)
        au, p = m42.auroc_lower(t.loc[~t.mmr, f"{c}_p"].to_numpy(), t.loc[t.mmr, f"{c}_p"].to_numpy())
        rows.append({"what": f"{c} expression (percentile in cancer): AUROC higher in MMR-deficient", "value": au, "p": p,
                     "n": int(t.mmr.sum())})
    base = "pct_up ~ mmr_i + C(project) + g_p + stromal_p + immune_p + ri_mean_p"
    for extra in ["", " + PTBP1_p + HNRNPC_p", " + " + " + ".join(f"{c}_p" for c in REGULATORS)]:
        f = smf.ols(base + extra, data=t).fit(cov_type="HC3")
        rows.append({"what": f"RV2 model{extra.replace('_p', '') or ''}: coefficient of MMR deficiency",
                     "value": float(f.params["mmr_i"]), "p": float(f.pvalues["mmr_i"]), "n": len(t)})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "posthoc.csv", index=False, float_format="%.4g")
    log(f"A-16 posthoc:\n{out.round(4).to_string(index=False)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["annotate", "test", "posthoc"], required=True)
    {"annotate": annotate, "test": test, "posthoc": posthoc}[ap.parse_args().only]()


if __name__ == "__main__":
    main()
