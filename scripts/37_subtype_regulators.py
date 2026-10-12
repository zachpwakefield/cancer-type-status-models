"""Step 37 (D-46): which splicing regulators differ between subtypes, and do they track the
events the subtype models use?

On the tumours of each subtype cohort (no model refitted):
1. scan         each of the 70 splicing factors of D-18 (log2 CPM) between the classes:
                Mann-Whitney (2 classes) or Kruskal-Wallis (> 2); BH within the cohort; each
                class's AUROC against the rest
2. association  as step 12 (D-18): the top 50 events of the cohort's model per view (SE, A3SS,
                A5SS, MXE, RI) against 500 random measured events, partial Spearman with every
                factor adjusting for the stromal and immune scores; excess = share of top
                events with |rho| > 0.3 minus the share among the background
3. predictions  R1-R6 of config/followup.yaml (ESRP1/2 in MESO and HNSC_mRNA; PTBP1 and the
                ACTN1/TPM1 smooth-muscle exons in SARC); BH over R1-R6
4. sensitivity  (A-6, exploratory, after R1-R6 were read) for R5/R6: the same partial Spearman
                also adjusting for the sample's mean MXE PSI (the global axis of D-24), and the
                rank of PTBP1 among all 70 factors for each event
5. aggregates   per-class quartiles of each predicted factor and, for R5/R6, the event's median
                PSI by PTBP1 quintile within LMS and within the other sarcomas (figure U8)

Outputs: results/subtype/regulators/ (tables; figure U8: scripts/36_subtype_figures.py)
Usage: python scripts/37_subtype_regulators.py [--workers 6]
"""
from __future__ import annotations

import argparse
import importlib.util
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from _common import REPO, log
from scipy.stats import kruskal, mannwhitneyu, rankdata, spearmanr

from ctsm.biology import bh, bio_config, feature_values, gex_by_symbol
from ctsm.matrix_io import load_view
from ctsm.subtype import TASK, load_tables, subtype_config

warnings.filterwarnings("ignore", category=RuntimeWarning)
_spec = importlib.util.spec_from_file_location("reg12", Path(__file__).with_name("12_bio_regulators.py"))
_reg12 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_reg12)
ranks, residualise, corr_cols = _reg12.ranks, _reg12.residualise, _reg12.corr_cols
FCFG = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["subtype_regulators"]
RC = bio_config()["regulators"]
SFS = [g for grp in RC["splicing_factors"].values() for g in grp]
FAMILY = {g: fam for fam, grp in RC["splicing_factors"].items() for g in grp}


def auroc(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """P(a > b) + 0.5 P(a = b) and the two-sided Mann-Whitney p."""
    r = mannwhitneyu(a, b, alternative="two-sided")
    return float(r.statistic / (len(a) * len(b))), float(r.pvalue)


def composition(G: pd.DataFrame) -> np.ndarray:
    return np.column_stack([rankdata(G[[g for g in RC["stromal_genes"] if g in G]].mean(1)),
                            rankdata(G[[g for g in RC["immune_genes"] if g in G]].mean(1))])


def one_cohort(args) -> tuple[pd.DataFrame, pd.DataFrame]:
    cohort, = args
    scfg = subtype_config()
    tab = load_tables(scfg)[cohort]
    G = gex_by_symbol(scfg, tab["file_id"], SFS + RC["stromal_genes"] + RC["immune_genes"])
    lab = tab["label"].to_numpy()
    classes = list(tab.drop_duplicates("y").sort_values("y")["label"])
    scan = []
    for g in [s for s in SFS if s in G and G[s].std() > 0]:
        x = G[g].to_numpy()
        groups = [x[lab == c] for c in classes]
        p = (mannwhitneyu(*groups).pvalue if len(classes) == 2 else kruskal(*groups).pvalue)
        row = {"cohort": cohort, "gene": g, "family": FAMILY[g], "n": len(x), "p": float(p)}
        for c in classes:
            row[f"auroc_{c}"] = auroc(x[lab == c], x[lab != c])[0]
            row[f"median_{c}"] = float(np.median(x[lab == c]))
        aur = np.array([row[f"auroc_{c}"] for c in classes])
        row["max_abs_auroc_dev"] = float(np.max(np.abs(aur - 0.5)))
        row["class_highest"] = classes[int(np.argmax(aur))]
        scan.append(row)
    scan = pd.DataFrame(scan)
    scan["q_within_cohort"] = bh(scan["p"].to_numpy())

    cov = np.column_stack([np.ones(len(tab)), composition(G)])
    sf_here = [s for s in SFS if s in G and G[s].std() > 0]
    Rs = residualise(ranks(G[sf_here].to_numpy()), cov)
    rng = np.random.default_rng(0)
    assoc = []
    for view in FCFG["association"]["views"]:
        imp = pd.read_parquet(scfg.data_dir / "bio_cache" / "importance" / f"{view}__{cohort}__all.parquet")
        top = imp["feature_id"].head(FCFG["association"]["top_events_per_view"]).tolist()
        pf = pd.read_parquet(scfg.matrix_dir / view / "prefilter.parquet", columns=[cohort])[cohort].to_numpy()
        fall = pd.read_parquet(scfg.matrix_dir / view / "features.parquet", columns=["feature_id"])
        cand = np.setdiff1d(fall.loc[pf, "feature_id"].to_numpy(), top)
        bgr = rng.choice(cand, size=min(FCFG["association"]["n_background_events"], len(cand)), replace=False).tolist()
        fids = np.array(top + bgr)
        X = feature_values(scfg, view, tab["file_id"], fids)
        ok = np.isfinite(X).mean(0) >= FCFG["association"]["min_finite"]
        X = np.where(np.isfinite(X), X, np.nanmedian(X, axis=0))
        ok &= np.nanstd(X, axis=0) > 0
        rho = corr_cols(Rs, residualise(ranks(X[:, ok]), cov))
        is_top = np.isin(fids[ok], top)
        thr = FCFG["association"]["rho"]
        for i, s in enumerate(sf_here):
            r = np.abs(rho[i])
            assoc.append({"cohort": cohort, "view": view, "gene": s, "n_top": int(is_top.sum()),
                          "n_background": int((~is_top).sum()),
                          "share_top_tracked": float(np.mean(r[is_top] > thr)),
                          "share_background_tracked": float(np.mean(r[~is_top] > thr))})
    assoc = pd.DataFrame(assoc)
    assoc["excess"] = assoc["share_top_tracked"] - assoc["share_background_tracked"]
    return scan, assoc


def predictions(scfg) -> pd.DataFrame:
    tables = load_tables(scfg)
    rows = []
    for pr in FCFG["predictions"]:
        tab = tables[pr["cohort"]]
        G = gex_by_symbol(scfg, tab["file_id"], [pr["gene"]] + RC["stromal_genes"] + RC["immune_genes"])
        x = G[pr["gene"]].to_numpy()
        lab = tab["label"].to_numpy()
        row = {"id": pr["id"], "cohort": pr["cohort"], "gene": pr["gene"], "source": pr["source"]}
        if "event" not in pr:
            hi, lo = np.isin(lab, pr["higher_in"]), np.isin(lab, pr["lower_in"])
            a, p = auroc(x[hi], x[lo])
            row.update(test="AUROC (higher_in vs lower_in)", n_a=int(hi.sum()), n_b=int(lo.sum()),
                       median_a=float(np.median(x[hi])), median_b=float(np.median(x[lo])),
                       estimate=a, predicted_sign=1, observed_sign=int(np.sign(a - 0.5)), p=p,
                       claim=f"{pr['gene']} higher in {'+'.join(pr['higher_in'])} than {'+'.join(pr['lower_in'])}")
        else:
            psi = feature_values(scfg, "mxe", tab["file_id"], [pr["event"]])[:, 0].astype(float)
            ok = np.isfinite(psi)
            lms = lab == "LMS"
            shift = float(np.median(psi[ok & lms]) - np.median(psi[ok & ~lms]))
            dummies = np.column_stack([(lab == c).astype(float) for c in sorted(set(lab))[1:]])
            cov = np.column_stack([np.ones(ok.sum()), dummies[ok], composition(G.loc[tab["file_id"][ok]])])
            rx = residualise(rankdata(x[ok])[:, None], cov)[:, 0]
            ry = residualise(rankdata(psi[ok])[:, None], cov)[:, 0]
            r = spearmanr(rx, ry)
            raw = spearmanr(x[ok], psi[ok]).statistic
            within = spearmanr(x[ok & lms], psi[ok & lms])
            row.update(test="partial Spearman (class, stromal, immune)", n_a=int(ok.sum()), n_b=np.nan,
                       event=pr["event"], event_name=pr["event_name"], lms_minus_rest_psi=shift,
                       estimate=float(r.statistic), predicted_sign=-int(np.sign(shift)),
                       observed_sign=int(np.sign(r.statistic)), p=float(r.pvalue), rho_unadjusted=float(raw),
                       rho_within_lms=float(within.statistic), p_within_lms=float(within.pvalue),
                       claim=f"PTBP1 moves the {pr['event_name']} against the LMS direction")
        rows.append(row)
    d = pd.DataFrame(rows)
    d["q"] = bh(d["p"].to_numpy())
    d["verdict"] = np.where(d["q"] >= 0.05, "not resolved",
                            np.where(d["predicted_sign"] == d["observed_sign"], "supported", "against"))
    return d


def _event_inputs(scfg, tab, G, pr):
    lab = tab["label"].to_numpy()
    psi = feature_values(scfg, "mxe", tab["file_id"], [pr["event"]])[:, 0].astype(float)
    ok = np.isfinite(psi)
    dummies = np.column_stack([(lab == c).astype(float) for c in sorted(set(lab))[1:]])
    cov = np.column_stack([np.ones(ok.sum()), dummies[ok], composition(G.loc[tab["file_id"][ok]])])
    return psi, ok, cov


def sensitivity(scfg) -> pd.DataFrame:
    tab = load_tables(scfg)["SARC"]
    G = gex_by_symbol(scfg, tab["file_id"], SFS + RC["stromal_genes"] + RC["immune_genes"])
    M, _ = load_view(scfg.matrix_dir, "mxe", tab["file_id"], "SARC")
    glob = np.nanmean(M.get(), axis=1)
    sfs = [s for s in SFS if s in G and G[s].std() > 0]
    rows = []
    for pr in [p for p in FCFG["predictions"] if "event" in p]:
        psi, ok, cov = _event_inputs(scfg, tab, G, pr)
        for adj, C in [("class, stromal, immune", cov),
                       ("class, stromal, immune, mean MXE PSI", np.column_stack([cov, rankdata(glob[ok])]))]:
            ry = residualise(rankdata(psi[ok])[:, None], C)[:, 0]
            rr = []
            for g in sfs:
                r = spearmanr(residualise(rankdata(G[g].to_numpy()[ok])[:, None], C)[:, 0], ry)
                rr.append((g, float(r.statistic), float(r.pvalue)))
            d = pd.DataFrame(rr, columns=["gene", "rho", "p"]).sort_values("rho", key=np.abs, ascending=False)
            d = d.reset_index(drop=True)
            pt = d[d.gene == "PTBP1"].iloc[0]
            rows.append({"id": pr["id"], "event_name": pr["event_name"], "adjusted_for": adj,
                         "rho_ptbp1": pt.rho, "p_ptbp1": pt.p, "rank_ptbp1_of": f"{int(pt.name) + 1} of {len(d)}",
                         "top_factors": "; ".join(f"{r.gene} {r.rho:+.2f}" for r in d.head(5).itertuples())})
    return pd.DataFrame(rows)


def aggregates(scfg) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables = load_tables(scfg)
    q = []
    for pr in [p for p in FCFG["predictions"] if "event" not in p]:
        tab = tables[pr["cohort"]]
        x = gex_by_symbol(scfg, tab["file_id"], [pr["gene"]])[pr["gene"]].to_numpy()
        for c in sorted(set(tab["label"])):
            v = x[tab["label"].to_numpy() == c]
            q.append({"id": pr["id"], "cohort": pr["cohort"], "gene": pr["gene"], "class": c, "n": len(v),
                      "predicted": "higher" if c in pr["higher_in"] else "lower" if c in pr["lower_in"] else "",
                      "q25": np.quantile(v, 0.25), "median": np.median(v), "q75": np.quantile(v, 0.75)})
    tab = tables["SARC"]
    G = gex_by_symbol(scfg, tab["file_id"], ["PTBP1"])
    lab = tab["label"].to_numpy()
    b = []
    for pr in [p for p in FCFG["predictions"] if "event" in p]:
        psi = feature_values(scfg, "mxe", tab["file_id"], [pr["event"]])[:, 0].astype(float)
        for grp, m in [("LMS", lab == "LMS"), ("other sarcomas", lab != "LMS")]:
            m = m & np.isfinite(psi)
            x, y = G["PTBP1"].to_numpy()[m], psi[m]
            bins = pd.qcut(rankdata(x), 5, labels=False)
            for k in range(5):
                s = bins == k
                b.append({"id": pr["id"], "event_name": pr["event_name"], "group": grp, "ptbp1_quintile": k + 1,
                          "n": int(s.sum()), "median_ptbp1_log2cpm": float(np.median(x[s])),
                          "q25_psi": float(np.quantile(y[s], 0.25)), "median_psi": float(np.median(y[s])),
                          "q75_psi": float(np.quantile(y[s], 0.75))})
    return pd.DataFrame(q), pd.DataFrame(b)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", choices=["all", "post"], default="all", help="post: steps 4-5 only")
    a = ap.parse_args()
    scfg = subtype_config()
    out = scfg.results_dir / TASK / "regulators"
    out.mkdir(parents=True, exist_ok=True)
    if a.only == "post":
        post(scfg, out)
        return
    cohorts = list(load_tables(scfg))
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(one_cohort, [(c,) for c in cohorts]))
    scan = pd.concat([r[0] for r in res], ignore_index=True)
    assoc = pd.concat([r[1] for r in res], ignore_index=True)
    scan.to_csv(out / "factor_by_class.csv", index=False)
    assoc.to_csv(out / "factor_event_association.csv", index=False)
    w = assoc.assign(t=assoc.share_top_tracked * assoc.n_top, b=assoc.share_background_tracked * assoc.n_background)
    pooled = w.groupby(["cohort", "gene"]).agg(t=("t", "sum"), nt=("n_top", "sum"), b=("b", "sum"),
                                               nb=("n_background", "sum")).reset_index()
    pooled["excess"] = pooled.t / pooled.nt - pooled.b / pooled.nb
    pooled.drop(columns=["t", "b"]).to_csv(out / "factor_association_excess_pooled.csv", index=False)
    summ = (scan.groupby("cohort").agg(n_factors=("gene", "size"), n_q05=("q_within_cohort", lambda q: int((q < 0.05).sum())),
                                       n_dev_0p2=("max_abs_auroc_dev", lambda v: int((v >= 0.2).sum())))
            .reset_index())
    summ.to_csv(out / "factor_scan_summary.csv", index=False)
    pr = predictions(scfg)
    pr.to_csv(out / "predictions.csv", index=False)
    log("scan:\n" + summ.to_string(index=False))
    log("predictions:\n" + pr[["id", "cohort", "gene", "estimate", "p", "q", "verdict"]].to_string(index=False))
    post(scfg, out)


def post(scfg, out) -> None:
    sens = sensitivity(scfg)
    sens.to_csv(out / "predictions_sensitivity.csv", index=False)
    q, b = aggregates(scfg)
    q.to_csv(out / "predicted_factor_by_class.csv", index=False)
    b.to_csv(out / "sarc_event_by_ptbp1_quintile.csv", index=False)
    log("sensitivity:\n" + sens.to_string(index=False))


if __name__ == "__main__":
    main()
