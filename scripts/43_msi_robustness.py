"""Step 43 (D-55): is the MSI repeat score of D-52 tumour-intrinsic?

U9 showed that splicing models lean on the hospital that sent the sample in some cohorts, STAD_TCGA
among them. The repeat score of step 42 involves no model, but MSI tumours could still cluster in
hospitals, compositions or RNA-handling conditions that move repeat exons. Four checks, fixed in
config/followup.yaml (msi_robustness) before they were computed:
  S1  within hospital: AUROC over MSI-rest pairs from the same tissue-source site; labels permuted
      within site
  S2  OLS score ~ MSI + site + stromal + immune + ri_mean (HC3); samples from a site with a single
      patient are left out, because their leverage is 1 and HC3 is undefined for them (A-10)
  S3  tumour minus own adjacent normal, patients with an MSI tumour against the other patients,
      scores built in the main matrix with each project's tumours and normals ranked together
  S4  the score rebuilt in the main matrix (its own feature filter), MSI against the rest
Family S = STAD_TCGA and COAD_TCGA (S3 pooled over UCEC, STAD, COAD); UCEC_TCGA reported alongside.

Outputs: results/subtype/msi_cis/robustness/ (aggregates); data.nosync/ctsm_v1/msi/ (per sample)
Usage: python scripts/43_msi_robustness.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import yaml
from _common import REPO, log
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

from ctsm.biology import gex_by_symbol, task_table
from ctsm.config import load_config
from ctsm.matrix_io import load_view

_spec = importlib.util.spec_from_file_location("msi42", Path(__file__).with_name("42_msi_repeats.py"))
m42 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m42)

CFG = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["msi_robustness"]
BIO = yaml.safe_load(open(REPO / "config" / "biology.yaml"))["regulators"]
DATA, LOCAL = m42.DATA, m42.LOCAL
OUT = m42.OUT / "robustness"
COHORTS = {"UCEC_TCGA": "UCEC", "STAD_TCGA": "STAD", "COAD_TCGA": "COAD"}
FAMILY = ["STAD_TCGA", "COAD_TCGA"]
N_PERM, SEED = 10_000, 0


def within_site_auroc(score: np.ndarray, msi: np.ndarray, site: np.ndarray) -> float:
    """P(score_MSI < score_rest) over MSI-rest pairs from the same site (ties half)."""
    num = den = 0.0
    for s in np.unique(site):
        m = site == s
        a, b = score[m & msi], score[m & ~msi]
        if len(a) and len(b):
            num += (a[:, None] < b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()
            den += len(a) * len(b)
    return num / den


def perm_within_site(score, msi, site, n: int, seed: int) -> tuple[float, float]:
    obs = within_site_auroc(score, msi, site)
    rng = np.random.default_rng(seed)
    idx = [np.flatnonzero(site == s) for s in np.unique(site)]
    null = np.empty(n)
    lab = msi.copy()
    for k in range(n):
        for ix in idx:
            lab[ix] = rng.permutation(msi[ix])
        null[k] = within_site_auroc(score, lab, site)
    p = (1 + (np.abs(null - 0.5) >= abs(obs - 0.5) - 1e-12).sum()) / (1 + n)
    return obs, float(p)


def covariates(cohort: str, tab: pd.DataFrame) -> pd.DataFrame:
    cfg = load_config()
    g = gex_by_symbol(cfg, tab["file_id"], BIO["stromal_genes"] + BIO["immune_genes"])
    out = pd.DataFrame(index=tab["file_id"])
    out["stromal"] = g[[x for x in BIO["stromal_genes"] if x in g]].mean(axis=1).to_numpy()
    out["immune"] = g[[x for x in BIO["immune_genes"] if x in g]].mean(axis=1).to_numpy()
    vm, _ = load_view(DATA / "subtype" / "matrices", "ri", tab["file_id"], cohort)
    out["ri_mean"] = np.nanmean(vm.get(), axis=1)
    return out


def main_matrix_scores(project: str, cohort: str) -> pd.DataFrame:
    """Repeat score per sample in the main SE matrix: the project's tumours and normals ranked
    together (S3) and the cohort's labelled tumours alone (S4), exons chosen by the msi_cis rules
    on the main matrix's own pre-filter."""
    cfg = load_config()
    ann_path = LOCAL / "se_repeats_main.parquet"
    if not ann_path.exists():
        m42.annotate(DATA / "matrices" / "se" / "features.parquet", ann_path)
    ann = pd.read_parquet(ann_path).set_index("feature_id")
    tab = task_table(cfg, "tumor_status", project).reset_index(drop=True)
    sub = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{cohort}.parquet").set_index("case_id")["label"]
    lab = tab["case_id"].map(sub)
    grp = np.where(tab["role"].eq("normal"), "normal",
                   np.where(lab.eq("MSI"), "MSI", np.where(lab.eq("POLE"), "POLE",
                            np.where(lab.isna(), "unlabelled", "rest")))).astype(object)
    tab["group"], tab["tumour_label"] = grp, lab.to_numpy()
    vm, feats = load_view(cfg.matrix_dir, "se", tab["file_id"], project)
    X = vm.get()
    a = ann.loc[feats["feature_id"]].reset_index()
    out = []
    for setting, rows in [("tumours+normals", np.ones(len(tab), bool)),
                          ("labelled tumours", np.isin(grp, ["MSI", "rest", "POLE"]))]:
        Xs, g = X[rows], grp[rows]
        fin = np.isfinite(Xs)
        ok = (fin.mean(axis=0) >= 0.8) & (fin[g == "MSI"].sum(axis=0) >= 10) & (fin[g == "rest"].sum(axis=0) >= 10)
        ok &= (a["acc_intron_len"] >= m42.MIN_INTRON["acc"]).to_numpy() & a["acc_ag"].eq("AG").to_numpy()
        P = pd.DataFrame(Xs[:, ok]).rank(pct=True, axis=0).to_numpy()
        dec = pd.qcut(np.nanmean(Xs[:, ok], axis=0), 10, labels=False, duplicates="drop")
        accT = a.loc[ok, "acc_T"].to_numpy()
        sc = m42.set_score(P, accT >= 10, accT <= 5, dec)
        d = tab.loc[rows, ["file_id", "case_id", "role", "tss", "group", "tumour_label"]].copy()
        d["score"], d["setting"], d["cohort"], d["n_exons"], d["n_R"] = sc, setting, cohort, int(ok.sum()), int((accT >= 10).sum())
        out.append(d)
    log(f"{cohort} main matrix: {len(tab)} samples ({(grp == 'normal').sum()} normals)")
    return pd.concat(out, ignore_index=True)


def pair_summary(pairs: pd.DataFrame) -> pd.DataFrame:
    """n and median tumour-minus-normal difference per cohort and group; medians of fewer than
    5 patients are left blank (D-12)."""
    s = pairs.groupby(["cohort", "group"])["diff"].agg(n="size", median="median").reset_index()
    s.loc[s.n < 5, "median"] = np.nan
    return s


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sc = pd.read_parquet(LOCAL / "sample_scores.parquet")
    pred, expl, local = [], [], []

    def add(pid, cohort, effect, sign, p, n, what, lo=np.nan, hi=np.nan):
        pred.append({"id": pid, "cohort": cohort, "in_family": cohort in FAMILY or cohort == "pooled",
                     "what": what, "effect": effect, "predicted_sign": sign, "ci_low": lo, "ci_high": hi,
                     "p": p, "n": n})

    for cohort in COHORTS:
        tab = pd.read_parquet(DATA / "subtype" / "tasks" / f"subtype__{cohort}.parquet")
        s = sc[sc.cohort.eq(cohort)].merge(tab[["file_id", "tss"]], on="file_id")
        s = s.join(covariates(cohort, s), on="file_id")
        s = s[s.group.isin(["MSI", "rest"])].reset_index(drop=True)
        s["msi"] = s.group.eq("MSI").astype(int)
        local.append(s.assign(cohort=cohort))
        score, msi, site = s["score_acc"].to_numpy(), s["msi"].to_numpy().astype(bool), s["tss"].to_numpy()
        # S1
        obs, p = perm_within_site(score, msi, site, N_PERM, SEED)
        mixed = s.groupby("tss").msi.transform(lambda x: 0 < x.sum() < len(x))
        add("S1", cohort, obs - 0.5, 1, p, f"MSI {int(msi[mixed].sum())} in {s.loc[mixed, 'tss'].nunique()} mixed sites",
            "within-site AUROC(score lower in MSI) - 0.5")
        # S2
        keep = s.groupby("tss")["tss"].transform("size") > 1
        fit = smf.ols("score_acc ~ msi + C(tss) + stromal + immune + ri_mean", data=s[keep]).fit(cov_type="HC3")
        ci = fit.conf_int().loc["msi"]
        add("S2", cohort, float(fit.params["msi"]), -1, float(fit.pvalues["msi"]),
            f"{int(keep.sum())} patients ({int((~keep).sum())} from single-patient sites left out)",
            "OLS coefficient of MSI, adjusted for site, stromal, immune, ri_mean (HC3)", ci[0], ci[1])
        # exploratory
        lo_site = []
        for t in s.loc[s.msi.eq(1), "tss"].unique():
            r = s[s.tss.ne(t)]
            lo_site.append(m42.auroc_lower(r.loc[r.msi.eq(1), "score_acc"].to_numpy(), r.loc[r.msi.eq(0), "score_acc"].to_numpy())[0])
        low = s[s.immune <= s.immune.median()]
        a_low = m42.auroc_lower(low.loc[low.msi.eq(1), "score_acc"].to_numpy(), low.loc[low.msi.eq(0), "score_acc"].to_numpy())
        row = {"cohort": cohort, "leave_one_site_out_auroc_min": min(lo_site), "leave_one_site_out_auroc_max": max(lo_site),
               "auroc_lower_immune_half": a_low[0], "p_lower_immune_half": a_low[1],
               "n_msi_lower_immune_half": int(low.msi.sum()),
               "fe_only_msi_coef": float(smf.ols("score_acc ~ msi + C(tss)", data=s[keep]).fit(cov_type="HC3").params["msi"])}
        for cov in ["stromal", "immune", "ri_mean"]:
            for g, gg in [("MSI", 1), ("rest", 0)]:
                d = s[s.msi.eq(gg)]
                row[f"rho_{cov}_{g}"] = spearmanr(d["score_acc"], d[cov]).statistic
            row[f"coef_{cov}"] = float(fit.params[cov])
            row[f"p_{cov}"] = float(fit.pvalues[cov])
        expl.append(row)
        log(f"{cohort}: S1 {obs:.3f} (p {p:.4f}); S2 {fit.params['msi']:+.4f} (p {fit.pvalues['msi']:.2g})")
    # S3, S4 in the main matrix
    mm = pd.concat([main_matrix_scores(p, c) for c, p in COHORTS.items()], ignore_index=True)
    mm.to_parquet(LOCAL / "main_matrix_scores.parquet", index=False)
    tn = mm[mm.setting.eq("tumours+normals")]
    pairs = []
    for c in COHORTS:
        d = tn[tn.cohort.eq(c)]
        t = d[d.role.eq("tumour")].groupby("case_id").agg(score_t=("score", "mean"), group=("group", "first"))
        n = d[d.role.eq("normal")].groupby("case_id")["score"].mean().rename("score_n")
        j = t.join(n, how="inner")
        j["diff"], j["cohort"] = j.score_t - j.score_n, c
        pairs.append(j.reset_index())
    pairs = pd.concat(pairs, ignore_index=True)
    a = pairs.loc[pairs.group.eq("MSI"), "diff"].to_numpy()
    b = pairs.loc[pairs.group.ne("MSI"), "diff"].to_numpy()
    au, p = m42.auroc_lower(a, b)
    add("S3", "pooled", au - 0.5, 1, p, f"MSI patients {len(a)}, others {len(b)}",
        "AUROC(tumour - normal difference lower in MSI patients) - 0.5")
    b2 = pairs.loc[pairs.group.eq("rest"), "diff"].to_numpy()
    nm = tn[tn.role.eq("normal")].copy()
    nm["patient_msi"] = nm.case_id.isin(pairs.loc[pairs.group.eq("MSI"), "case_id"])
    expl.append({"cohort": "pooled", "s3_auroc_vs_labelled_rest_only": m42.auroc_lower(a, b2)[0],
                 "s3_p_vs_labelled_rest_only": m42.auroc_lower(a, b2)[1], "n_labelled_rest_pairs": len(b2),
                 "normals_auroc_msi_patient_lower": m42.auroc_lower(nm.loc[nm.patient_msi, "score"].to_numpy(),
                                                                   nm.loc[~nm.patient_msi, "score"].to_numpy())[0],
                 "msi_pairs_with_tumour_below_normal": int((a < 0).sum()), "msi_pairs": len(a),
                 "other_pairs_with_tumour_below_normal": int((b < 0).sum()), "other_pairs": len(b)})
    pairs.drop(columns=[]).to_parquet(LOCAL / "tumour_normal_pairs.parquet", index=False)
    lt = mm[mm.setting.eq("labelled tumours")]
    for c in COHORTS:
        d = lt[lt.cohort.eq(c)]
        au, p = m42.auroc_lower(d.loc[d.group.eq("MSI"), "score"].to_numpy(), d.loc[d.group.eq("rest"), "score"].to_numpy())
        add("S4", c, au - 0.5, 1, p, f"{int(d.n_exons.iloc[0])} exons, R {int(d.n_R.iloc[0])}",
            "AUROC(score lower in MSI) - 0.5, main matrix")
        r = sc[sc.cohort.eq(c)].merge(d[["file_id", "score"]], on="file_id")
        expl.append({"cohort": c, "s4_rho_with_subtype_matrix_score": spearmanr(r.score_acc, r.score).statistic})
    pr = pd.DataFrame(pred)
    fam = pr.in_family
    pr["q"] = np.nan
    pr.loc[fam, "q"] = multipletests(pr.loc[fam, "p"], method="fdr_bh")[1]
    pr["verdict"] = [m42.verdict(e, sg, q) if f else ("reported, not in family: " + ("direction as predicted" if np.sign(e) == sg else "direction opposite"))
                     for e, sg, q, f in zip(pr.effect, pr.predicted_sign, pr.q, pr.in_family)]
    pr.to_csv(OUT / "predictions.csv", index=False)
    pd.DataFrame(expl).to_csv(OUT / "exploratory.csv", index=False)
    pd.concat(local, ignore_index=True).to_parquet(LOCAL / "robustness_samples.parquet", index=False)
    pair_summary(pairs).to_csv(OUT / "tumour_normal_pairs_summary.csv", index=False)
    log("predictions:\n" + pr[["id", "cohort", "effect", "p", "q", "verdict", "n"]].to_string(index=False))
    log("exploratory:\n" + pd.DataFrame(expl).round(3).T.to_string())


if __name__ == "__main__":
    main()
