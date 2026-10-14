"""Step 60 (D-75): how large is the MSI splicing shift, and how is it spread over tumours?

Every earlier MSI result was in percentiles within a cancer. Here the shifts are measured in PSI
units (0-1), within cancer, for the code's classes of exons: R_down (acceptor run 10-19 nt), R_up
(run >= 20 nt with the branch point >= 36 nt), other long runs, and B (run <= 5 nt, background).
MG1-MG3 ask whether the shift scales with tumour purity among MMR-deficient tumours (a change in
tumour cells is diluted by normal cells); MG4-MG5 whether each exon's change sits in a subset of
MMR-deficient tumours rather than shifting all of them a little. config/followup.yaml: msi_magnitude.

Inputs:  data.nosync/ctsm_v1/msi/models_tumours.parquet (step 54), the pan-cancer SE matrix, the
         exon classes of step 52, ABSOLUTE purity (all local)
Outputs: results/subtype/msi_cis/code/magnitude/ (aggregates; per-exon shifts are medians over
         >= 20 tumours); data.nosync/ctsm_v1/msi/magnitude_{tumours,exons}.parquet (local)
Usage: python scripts/60_msi_magnitude.py [--only compute|report]
"""
from __future__ import annotations

import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from _common import log


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m52 = _load("msi52", "52_msi_reversal.py")
m51 = m52.m51
m50, m42 = m51.m50, m51.m42
LOCAL, PAN, OUT, TRAIN = m51.LOCAL, m51.PAN, m51.OUT / "magnitude", m51.TRAIN
MIN_FINITE, SHIFT, N_PERM, MIN_HELD, MIN_MEDIAN = 0.8, 0.10, 10_000, 20, 5
CLASSES = ["R_down", "R_up", "long, other", "B"]
BINS = [("10–11", 10, 11, None), ("12–13", 12, 13, None), ("14–15", 14, 15, None), ("16–19", 16, 19, None),
        ("≥ 20, bp < 36 or none", 20, 999, "near"), ("≥ 20, bp ≥ 36", 20, 999, "far")]


def exon_classes() -> pd.DataFrame:
    """Step 52's classes, regrouped: R_down, R_up, other runs >= 20 nt, B; and run-length bins."""
    c = m52.pan_classes()
    at = c.acc_T.to_numpy(float)
    down, up, b = m51.arms(at, c.d_bp.to_numpy(float))
    c["mclass"] = np.select([down, up, at >= 20, b], CLASSES, "")
    far = np.nan_to_num(c.d_bp.to_numpy(float), nan=-1) >= m51.UP_CUT
    c["bin"] = ""
    for name, lo, hi, bp in BINS:
        k = (at >= lo) & (at <= hi)
        k &= ~far if bp == "near" else far if bp == "far" else True
        c.loc[k, "bin"] = name
    c.loc[b, "bin"] = "B (run ≤ 5)"
    return c


# ---------------------------------------------------------------- statistics (unit-tested)

def deviations(X: np.ndarray, prof: np.ndarray, dec: np.ndarray, is_b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """dev = PSI minus the proficient median of each exon; adj = dev minus the tumour's mean dev over
    the B exons of the exon's mean-PSI decile."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        dev = X - np.nanmedian(X[prof], axis=0)
        bmean = np.full((X.shape[0], int(dec.max()) + 1), np.nan)
        for d in range(bmean.shape[1]):
            k = is_b & (dec == d)
            if k.any():
                bmean[:, d] = np.nanmean(dev[:, k], axis=1)
    return dev, dev - bmean[:, dec]


def tail_excess(P: np.ndarray, M: np.ndarray, shift: np.ndarray, lower: bool) -> tuple[np.ndarray, np.ndarray]:
    """Per column: the share of M (MMR-deficient) beyond the 5th (lower) or 95th percentile of P
    (proficient), and the share expected if every M value were a P value moved by `shift`."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q = np.nanpercentile(P, 5 if lower else 95, axis=0)
        fm, fp = np.isfinite(M), np.isfinite(P)
        beyond = (M < q) if lower else (M > q)
        moved = ((P + shift) < q) if lower else ((P + shift) > q)
        obs = (beyond & fm).sum(0) / np.maximum(fm.sum(0), 1)
        exp = (moved & fp).sum(0) / np.maximum(fp.sum(0), 1)
    return np.where(fm.any(0), obs, np.nan), np.where(fp.any(0), exp, np.nan)


def strat_perm(x: np.ndarray, is_a: np.ndarray, strata: np.ndarray, rng, n: int,
               chunk: int = 500) -> tuple[float, float, float]:
    """mean(x[a]) - mean(x[~a]); its null from permuting `is_a` within strata. Returns the observed
    difference, the stratum-matched difference (within-stratum a minus not-a, weighted by the number
    of a) and a two-sided permutation p."""
    x, is_a = np.asarray(x, float), np.asarray(is_a, bool)
    na, nb, tot = int(is_a.sum()), int((~is_a).sum()), float(x.sum())
    obs = float(x[is_a].mean() - x[~is_a].mean())
    s_a = np.zeros(n)
    within, weight = 0.0, 0
    for s in np.unique(strata):
        m = strata == s
        v, k = x[m], int(is_a[m].sum())
        if k == 0:
            continue
        if k < len(v):
            within += k * (v[is_a[m]].mean() - v[~is_a[m]].mean())
            weight += k
        if k == len(v):
            s_a += v.sum()
            continue
        for i in range(0, n, chunk):
            r = rng.random((min(chunk, n - i), len(v)))
            s_a[i:i + len(r)] += v[np.argpartition(r, k - 1, axis=1)[:, :k]].sum(1)
    null = s_a / na - (tot - s_a) / nb
    p = 2 * min((null >= obs).sum() + 1, (null <= obs).sum() + 1) / (n + 1)
    return obs, within / max(weight, 1) if weight else np.nan, float(min(1.0, p))


# ---------------------------------------------------------------- compute (per-sample data stays local)

def compute() -> None:
    from ctsm.matrix_io import load_view
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")[["file_id", "case_id", "project", "n_hp", "mmr"]]
    t["purity"] = t.case_id.map(m50.purity_by_patient())
    cls = exon_classes()
    mclass_all = cls.mclass.to_numpy()
    per_tumour, per_exon, held = [], [], []
    for proj, s in t.groupby("project"):
        X = load_view(PAN, "se", s.file_id)[0].get()
        if np.nanmax(X) > 1.0001:
            raise ValueError("PSI expected on 0-1")
        ok = cls.struct.to_numpy() & (np.isfinite(X).mean(axis=0) >= MIN_FINITE)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mean_psi = np.nanmean(X[:, ok], axis=0)
        dec_ok = np.asarray(pd.qcut(mean_psi, 10, labels=False, duplicates="drop"), dtype=int)
        use = mclass_all[ok] != ""                          # deciles over all measured exons, as steps 52-54
        idx = np.flatnonzero(ok)[use]
        Xo, dec, mc, mean_psi = X[:, idx], dec_ok[use], mclass_all[idx], mean_psi[use]
        y = s.mmr.to_numpy()
        dev, adj = deviations(Xo, ~y, dec, mc == "B")
        fin = np.isfinite(dev)
        f = s[["file_id", "project", "mmr", "purity", "n_hp"]].reset_index(drop=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for c_, col in [("R_down", "arm_down"), ("R_up", "arm_up")]:
                f[col] = np.nanmean(adj[:, mc == c_], axis=1)
        for c_, sign, col in [("R_down", -1, "down10_R_down"), ("R_up", 1, "up10_R_up"),
                              ("B", -1, "down10_B"), ("B", 1, "up10_B")]:
            k = mc == c_
            hit = (dev[:, k] * sign >= SHIFT) & fin[:, k]
            f["n_" + col] = hit.sum(1)
            f[col] = hit.sum(1) / np.maximum(fin[:, k].sum(1), 1)
        per_tumour.append(f)
        P, M = Xo[~y], Xo[y]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med_p = np.nanmedian(P, axis=0)
        if proj in TRAIN:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                shift = np.nanmedian(M, axis=0) - med_p
                grand = np.nanmean(Xo, axis=0)
                tot = np.nansum((Xo - grand) ** 2, axis=0)
                between = (np.isfinite(M).sum(0) * (np.nanmean(M, 0) - grand) ** 2 +
                           np.isfinite(P).sum(0) * (np.nanmean(P, 0) - grand) ** 2)
            lo_o, lo_e = tail_excess(P, M, shift, lower=True)
            hi_o, hi_e = tail_excess(P, M, shift, lower=False)
            per_exon.append(pd.DataFrame({
                "group": proj, "feature_index": idx, "decile": dec, "mclass": mc, "shift": shift,
                "low_obs": lo_o, "low_exp": lo_e, "high_obs": hi_o, "high_exp": hi_e,
                "eta2": between / np.where(tot > 0, tot, np.nan), "n_mmr": np.isfinite(M).sum(0),
                "mean_psi": mean_psi, "median_proficient": med_p}))
        elif y.any():
            # each held-out MMR-deficient tumour: dev, and its percentile among its cancer's proficient tumours
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                u = (P[None, :, :] < M[:, None, :]).sum(1) / np.maximum(np.isfinite(P).sum(0), 1)
            r, j = np.nonzero(np.isfinite(M))
            held.append(pd.DataFrame({"feature_index": idx[j], "dev": dev[y][r, j], "u": u[r, j]}))
        log(f"{proj}: {len(s)} tumours, {int(y.sum())} MMR-deficient, {len(idx):,} exons")
    tum = pd.concat(per_tumour, ignore_index=True)
    tum.to_parquet(LOCAL / "magnitude_tumours.parquet", index=False)
    h = pd.concat(held, ignore_index=True).groupby("feature_index").agg(
        shift=("dev", "median"), n_mmr=("dev", "size"),
        low_obs=("u", lambda v: float((v < 0.05).mean())), high_obs=("u", lambda v: float((v > 0.95).mean())))
    h = h.reset_index().assign(group="held out")
    h["mclass"] = mclass_all[h.feature_index]
    ex = pd.concat([*per_exon, h], ignore_index=True)
    ex = ex.join(cls[["feature_id", "symbol", "acc_T", "d_bp", "bin"]], on="feature_index")
    ex.to_parquet(LOCAL / "magnitude_exons.parquet", index=False)
    log(f"compute done: {len(tum):,} tumours, {len(ex):,} exon x group rows")


# ---------------------------------------------------------------- report

def _q(v: pd.Series, q: float) -> float:
    return float(v.quantile(q)) if len(v) >= MIN_MEDIAN else np.nan


def tests(tum: pd.DataFrame, ex: pd.DataFrame) -> pd.DataFrame:
    import statsmodels.formula.api as smf
    rng = np.random.default_rng(m51.SEED)
    rows = []
    md = tum[tum.mmr & tum.purity.notna()].reset_index(drop=True)
    dums = [md.project.eq(c).to_numpy(float) for c in sorted(md.project.unique())[1:]]
    for pid, col, sign in [("MG1", "arm_down", -1), ("MG2", "arm_up", 1)]:
        r, p, n = m42.partial_spearman(md[col].to_numpy(float), md.purity.to_numpy(float), dums)
        rows.append({"id": pid, "what": f"partial Spearman of {col} with purity given cancer, MMR-deficient tumours",
                     "effect": r, "predicted_sign": sign, "p": p, "n": f"{n} tumours"})
    d3 = tum[tum.purity.notna() & tum.arm_down.notna()].assign(mmr_i=lambda d: d.mmr.astype(int))
    fit = smf.ols("arm_down ~ purity * mmr_i + C(project)", data=d3).fit(cov_type="HC3")
    ci = fit.conf_int().loc["purity:mmr_i"]
    rows.append({"id": "MG3", "what": "OLS arm_down ~ purity * MMR + cancer (HC3): purity x MMR coefficient",
                 "effect": float(fit.params["purity:mmr_i"]), "predicted_sign": -1,
                 "p": float(fit.pvalues["purity:mmr_i"]), "ci_low": float(ci[0]), "ci_high": float(ci[1]),
                 "n": f"{len(d3):,} tumours ({int(d3.mmr.sum())} MMR-deficient)"})
    tr = ex[ex.group.isin(TRAIN)].assign(low_x=lambda d: d.low_obs - d.low_exp, high_x=lambda d: d.high_obs - d.high_exp)
    for pid, c_, col in [("MG4", "R_down", "low_x"), ("MG5", "R_up", "high_x")]:
        sub = tr[tr.mclass.isin([c_, "B"]) & tr[col].notna()]
        strata = (sub.group + "_" + sub.decile.astype(str)).to_numpy()
        obs, matched, p = strat_perm(sub[col].to_numpy(float), sub.mclass.eq(c_).to_numpy(), strata, rng, N_PERM)
        rows.append({"id": pid, "what": f"mean {col} of {c_} minus B exons, training cancers", "effect": obs,
                     "predicted_sign": 1, "p": p, "decile_matched_difference": matched,
                     "mean_class": float(sub.loc[sub.mclass.eq(c_), col].mean()),
                     "mean_B": float(sub.loc[sub.mclass.eq("B"), col].mean()),
                     "n": f"{int(sub.mclass.eq(c_).sum())} {c_} and {int(sub.mclass.eq('B').sum())} B exon x cancer rows"})
    return m51.predictions(rows, "magnitude")


def exploratory(tum: pd.DataFrame, ex: pd.DataFrame) -> pd.DataFrame:
    import statsmodels.formula.api as smf
    out = []

    def add(what, value, group, cls_, n, **kw):
        out.append({"what": what, "group": group, "class": cls_, "value": value, "n": n, **kw})

    ex = ex[ex.group.isin(TRAIN) | ex.n_mmr.ge(MIN_HELD)]
    for g, sub in ex.groupby("group"):
        for c_ in CLASSES:
            s2 = sub[sub.mclass.eq(c_)]
            v = s2["shift"].dropna()
            for q in [0.05, 0.25, 0.5, 0.75, 0.95]:
                add(f"exon shift (PSI), quantile {q}", _q(v, q), g, c_, len(v))
            for thr in [0.05, 0.10]:
                add(f"share of exons with shift <= -{thr:.2f}", float((v <= -thr).mean()), g, c_, len(v))
                add(f"share of exons with shift >= +{thr:.2f}", float((v >= thr).mean()), g, c_, len(v))
            if g in TRAIN:
                e2 = s2.eta2.dropna()
                add("variance of PSI explained by MMR status (eta2), median", _q(e2, 0.5), g, c_, len(e2))
                add("variance explained, 95th percentile", _q(e2, 0.95), g, c_, len(e2))
                for side in ["low", "high"]:
                    add(f"{side}-tail share of MMR-deficient tumours, observed (mean over exons)",
                        float(s2[f"{side}_obs"].mean()), g, c_, len(s2))
                    add(f"{side}-tail share expected from a uniform shift (mean over exons)",
                        float(s2[f"{side}_exp"].mean()), g, c_, len(s2))
            else:
                add("low-tail share (below the proficient 5th percentile)", float(s2.low_obs.mean()), g, c_, len(s2))
                add("high-tail share (above the proficient 95th percentile)", float(s2.high_obs.mean()), g, c_, len(s2))
    pooled = ex.group.isin(TRAIN).map({True: "training (3 cancers)", False: "held out"})
    for g, sub in ex.groupby(pooled):
        for b, s2 in sub.groupby("bin"):
            if not b:
                continue
            v = s2["shift"].dropna()
            add("exon shift (PSI) by run-length bin: median", _q(v, 0.5), g, b, len(v),
                q25=_q(v, 0.25), q75=_q(v, 0.75), share_le_m005=float((v <= -0.05).mean()),
                share_ge_p005=float((v >= 0.05).mean()))
    tum = tum.assign(grp=np.where(tum.project.isin(TRAIN), tum.project, "held out"))
    cols = ["arm_down", "arm_up", "down10_R_down", "up10_R_up", "down10_B", "up10_B",
            "n_down10_R_down", "n_up10_R_up", "n_down10_B", "n_up10_B"]
    for (g, m_), sub in tum.groupby(["grp", "mmr"]):
        for col in cols:
            v = sub[col].dropna()
            add(f"per tumour {col}: median", _q(v, 0.5), g, "MMR-deficient" if m_ else "proficient", len(v),
                q25=_q(v, 0.25), q75=_q(v, 0.75))
    for col in ["arm_down", "arm_up"]:
        for lab, sub in [("MMR-deficient", tum[tum.mmr]), ("proficient", tum[~tum.mmr])]:
            sub = sub[sub.purity.notna() & sub[col].notna()]
            f2 = smf.ols(f"{col} ~ purity + C(project)", data=sub).fit()
            # average over cancers of the fitted value at purity 1 and at each cancer's median purity
            cancers = sub.project.unique()
            new = pd.DataFrame({"project": cancers, "purity": 1.0})
            med = pd.DataFrame({"project": cancers, "purity": [sub.purity[sub.project.eq(c)].median() for c in cancers]})
            w = sub.project.value_counts().reindex(cancers).to_numpy(float)
            add(f"{col} ~ purity + cancer: fitted at purity 1 (tumour-weighted over cancers)",
                float(np.average(f2.predict(new), weights=w)), "all", lab, len(sub), slope=float(f2.params.purity),
                at_median_purity=float(np.average(f2.predict(med), weights=w)),
                median_purity=float(sub.purity.median()))
    return pd.DataFrame(out)


def movers(ex: pd.DataFrame) -> pd.DataFrame:
    """The largest R_down and R_up shifts measured in all three training cancers, and MRE11 exon 5."""
    tr = ex[ex.group.isin(TRAIN)]
    w = tr.pivot_table(index="feature_index", columns="group", values="shift").dropna()
    w["mean_shift"] = w[TRAIN].mean(axis=1)
    tails = tr.groupby("feature_index")[["low_obs", "high_obs", "eta2"]].mean()
    meta = ex.drop_duplicates("feature_index").set_index("feature_index")[["feature_id", "symbol", "acc_T", "d_bp", "mclass"]]
    ho = ex[ex.group.eq("held out") & ex.n_mmr.ge(MIN_HELD)].set_index("feature_index")["shift"]
    w = w.join(meta).join(tails).assign(held_out_shift=ho.reindex(w.index))
    mre = w[w.symbol.eq("MRE11") & w.acc_T.eq(11)]
    out = pd.concat([w[w.mclass.eq("R_down")].nsmallest(12, "mean_shift"),
                     w[w.mclass.eq("R_up")].nlargest(8, "mean_shift"), mre])
    out = out[~out.index.duplicated()].reset_index()
    out["role"] = np.where(out.symbol.eq("MRE11") & out.acc_T.eq(11), "positive control (MRE11 exon 5)",
                           np.where(out.mclass.eq("R_down"), "largest R_down decrease", "largest R_up increase"))
    return out


def report() -> None:
    tum = pd.read_parquet(LOCAL / "magnitude_tumours.parquet")
    ex = pd.read_parquet(LOCAL / "magnitude_exons.parquet")
    OUT.mkdir(parents=True, exist_ok=True)
    tests(tum, ex)
    e = exploratory(tum, ex)
    e.to_csv(OUT / "exploratory.csv", index=False, float_format="%.5g")
    mv = movers(ex)
    mv.drop(columns=["feature_index"]).to_csv(OUT / "largest_movers.csv", index=False, float_format="%.4g")
    keep = ex[ex.group.isin(TRAIN) | ex.n_mmr.ge(MIN_HELD)]
    keep.drop(columns=["feature_index"]).to_csv(OUT / "exon_shifts.csv.gz", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    log("exploratory:\n" + e[~e.what.str.contains("quantile")].round(4).to_string(index=False))
    log("largest movers:\n" + mv.drop(columns=["feature_index", "feature_id"]).round(3).to_string(index=False))


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["compute", "report"])
    a = ap.parse_args()
    if a.only in (None, "compute"):
        compute()
    if a.only in (None, "report"):
        report()


if __name__ == "__main__":
    main()
