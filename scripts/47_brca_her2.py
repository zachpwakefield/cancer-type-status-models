"""Step 47 (D-58): which HER2 did each modality learn - ERBB2 amplification or the HER2-enriched programme?

Amplification from ABSOLUTE: total copy number of the segment holding ERBB2's GRCh37 anchor (the
midpoint of its MC3 call span) >= 2 x the sample's ploidy. Probabilities are the stored
out-of-fold P(Her2) of BRCA_PAM50: GEX, and the AS ensemble (mean of the 8 AS views, as step 34).

After H2 came out against its prediction, an exploratory knockout was added (A-12): in every stored
BRCA_PAM50 fold model, the inputs from the 17q12 amplicon are set to their training mean (z = 0)
and P(Her2) is recomputed for the held-out tumours; 20 random input sets of the same size are the
control.

Outputs: results/subtype/brca/her2/ (aggregates only)
Usage: python scripts/47_brca_her2.py
"""
from __future__ import annotations

import functools
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from _common import REPO, log
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

from ctsm.biology import cohort_matrix
from ctsm.config import load_config
from ctsm.model import MLP, logits_to_proba, predict_logits

_spec = importlib.util.spec_from_file_location("msi42", Path(__file__).with_name("42_msi_repeats.py"))
m42 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m42)

DATA = m42.DATA
EXT = DATA / "external" / "gdc_pancanatlas"
OUT = REPO / "results" / "subtype" / "brca" / "her2"
AS_VIEWS = ["afe", "ale", "hit", "se", "a5ss", "a3ss", "mxe", "ri"]
ERBB2 = "ENSG00000141736"
CORE = (39_600_000, 39_800_000)          # GRCh38 chr17: PPP1R1B-STARD3-TCAP-PNMT-PGAP3-ERBB2-MIEN1-GRB7
WIDE = (39_400_000, 40_200_000)
CORE_GENES = ["PPP1R1B", "STARD3", "TCAP", "PNMT", "PGAP3", "ERBB2", "MIEN1", "GRB7"]
WIDE_GENES = CORE_GENES + ["NEUROD2", "IKZF3", "ZPBP2", "GSDMB", "ORMDL3", "LRRC3C", "GSDMA", "PSMD3", "MED1", "CDK12"]
SEED, BOOT = 0, 2000


@functools.lru_cache(maxsize=1)
def erbb2_anchor() -> float:
    """GRCh37 position of ERBB2: the midpoint of its MC3 call span."""
    m = pd.read_csv(EXT / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=["Gene", "Chromosome", "Start_Position"],
                    dtype={"Gene": str, "Chromosome": str})
    m = m[m.Gene.eq(ERBB2) & m.Chromosome.eq("17")]
    return float((m.Start_Position.min() + m.Start_Position.max()) / 2)


def erbb2_copy_number() -> pd.DataFrame:
    anchor = erbb2_anchor()
    seg = pd.read_csv(EXT / "TCGA_mastercalls.abs_segtabs.fixed.txt", sep="\t",
                      usecols=["Sample", "Chromosome", "Start", "End", "Modal_Total_CN"])
    seg = seg[seg.Chromosome.eq(17) & (seg.Start <= anchor) & (seg.End >= anchor)]
    pur = pd.read_csv(EXT / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
    pur = pur[pur["call status"].eq("called")].set_index("array")
    cn = seg.drop_duplicates("Sample").set_index("Sample")["Modal_Total_CN"]
    d = pd.DataFrame({"cn": cn}).join(pur[["ploidy", "purity"]], how="inner")
    d["rel"] = d.cn / d.ploidy
    d["amplified"] = d.rel >= 2
    log(f"ERBB2 anchor chr17:{anchor:,.0f} (GRCh37); {len(d)} samples with a called copy number")
    return d


def clinical_her2() -> pd.Series:
    cl = pd.read_csv(EXT / "clinical_PANCAN_patient_with_followup.tsv", sep="\t", encoding="latin-1",
                     usecols=["bcr_patient_barcode", "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                              "lab_procedure_her2_neu_in_situ_hybrid_outcome_type"]).set_index("bcr_patient_barcode")
    ihc, ish = cl.iloc[:, 0], cl.iloc[:, 1]
    pos = ihc.eq("Positive") | ish.eq("Positive")
    neg = (ihc.eq("Negative") & ~ish.eq("Positive")) | (ish.eq("Negative") & ~ihc.eq("Positive"))
    return pd.Series(np.where(pos, "positive", np.where(neg, "negative", None)), index=cl.index)


def oof_her2(cohort: str = "BRCA_PAM50") -> pd.DataFrame:
    P = {}
    for v in AS_VIEWS + ["gex"]:
        o = pd.read_parquet(DATA / "subtype" / "runs" / "subtype" / v / cohort / "oof_predictions.parquet")
        P[v] = o.sort_values("file_id").set_index("file_id")
    d = P["gex"][["case_id", "sample_id", "label"]].copy()
    d["p_gex"] = P["gex"]["p_Her2"]
    d["p_as"] = np.mean([P[v]["p_Her2"].to_numpy() for v in AS_VIEWS], axis=0)
    return d


def auc(pos, neg) -> float:
    return float(mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg)))


def boot_diff(d: pd.DataFrame, grp: str, rng) -> tuple[float, float, float]:
    """AUROC(GEX) - AUROC(AS) for grp; paired bootstrap over tumours. (diff, two-sided p, CI half)"""
    y = d[grp].to_numpy()
    a, b = d.p_gex.to_numpy(), d.p_as.to_numpy()
    obs = auc(a[y], a[~y]) - auc(b[y], b[~y])
    bs = []
    for _ in range(BOOT):
        i = rng.integers(0, len(d), len(d))
        yy = y[i]
        if yy.all() or not yy.any():
            continue
        bs.append(auc(a[i][yy], a[i][~yy]) - auc(b[i][yy], b[i][~yy]))
    bs = np.asarray(bs)
    return obs, float(min(1.0, 2 * min((bs <= 0).mean(), (bs >= 0).mean()))), float(np.percentile(bs, 97.5) - np.percentile(bs, 2.5)) / 2


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    cn = erbb2_copy_number()
    d = oof_her2()
    d["s15"] = d.sample_id.str[:15]
    d = d.join(cn, on="s15")
    d["clin"] = d.case_id.map(clinical_her2())
    lum = d[d.label.isin(["LumA", "LumB"]) & d.amplified.notna()].copy()
    lum["amplified"] = lum.amplified.astype(bool)
    her = d[d.label.eq("Her2") & d.amplified.notna()].copy()
    her["amplified"] = her.amplified.astype(bool)
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    a, n1 = lum[lum.amplified], lum[~lum.amplified]
    add("H1", auc(a.p_gex, n1.p_gex) - 0.5, 1, float(mannwhitneyu(a.p_gex, n1.p_gex).pvalue),
        f"amplified {len(a)}, not {len(n1)}", "LumA/LumB: AUROC(GEX P(Her2) higher when amplified) - 0.5")
    diff, p, half = boot_diff(lum, "amplified", rng)
    add("H2", diff, 1, p, f"{len(lum)} tumours; 95% CI half-width {half:.3f}; AS AUROC {auc(a.p_as, n1.p_as):.3f}",
        "LumA/LumB: AUROC(GEX) - AUROC(AS ensemble), amplified against not")
    diff, p, half = boot_diff(her, "amplified", rng)
    ha, hn = her[her.amplified], her[~her.amplified]
    add("H3", diff, 1, p, f"amplified {len(ha)}, not {len(hn)}; 95% CI half-width {half:.3f}; "
        f"GEX {auc(ha.p_gex, hn.p_gex):.3f}, AS {auc(ha.p_as, hn.p_as):.3f}",
        "PAM50 Her2: AUROC(GEX) - AUROC(AS ensemble), amplified against not")
    add("H4", auc(hn.p_as, a.p_as) - 0.5, 1, float(mannwhitneyu(hn.p_as, a.p_as).pvalue),
        f"non-amplified Her2 {len(hn)}, amplified luminal {len(a)}",
        "AUROC(AS P(Her2) higher in non-amplified Her2 than amplified LumA/LumB) - 0.5")
    pr = pd.DataFrame(pred)
    pr["q"] = multipletests(pr.p, method="fdr_bh")[1]
    pr["verdict"] = [m42.verdict(e, s, q) for e, s, q in zip(pr.effect, pr.predicted_sign, pr.q)]
    pr.to_csv(OUT / "predictions.csv", index=False)
    log("predictions:\n" + pr[["id", "effect", "p", "q", "verdict", "n"]].to_string(index=False))

    # ---- exploratory
    ex = []
    gex_h4 = auc(hn.p_gex, a.p_gex)
    ex.append({"what": "H4 counterpart, GEX: AUROC(P(Her2) higher in non-amplified Her2 than amplified luminal)",
               "value": gex_h4, "n": f"{len(hn)}, {len(a)}"})
    for m in ["p_gex", "p_as"]:
        r = spearmanr(lum[m], lum.rel)
        ex.append({"what": f"LumA/LumB: Spearman({m}, ERBB2 copy number / ploidy)", "value": r.statistic, "n": len(lum)})
    cl = d[d.clin.notna()]
    for name, sub in [("LumA/LumB", cl[cl.label.isin(["LumA", "LumB"])]), ("Her2", cl[cl.label.eq("Her2")])]:
        y = sub.clin.eq("positive")
        if y.sum() >= 5 and (~y).sum() >= 5:
            for m in ["p_gex", "p_as"]:
                ex.append({"what": f"{name}: AUROC({m}, clinical HER2 positive against negative)",
                           "value": auc(sub.loc[y, m], sub.loc[~y, m]), "n": f"{int(y.sum())}, {int((~y).sum())}"})
    hn_c, a_c = cl[cl.label.eq("Her2") & cl.clin.eq("negative")], cl[cl.label.isin(["LumA", "LumB"]) & cl.clin.eq("positive")]
    if len(hn_c) >= 5 and len(a_c) >= 5:
        for m in ["p_gex", "p_as"]:
            ex.append({"what": f"H4 with clinical HER2 ({m}): HER2-negative Her2 against HER2-positive luminal",
                       "value": auc(hn_c[m], a_c[m]), "n": f"{len(hn_c)}, {len(a_c)}"})
    both = d[d.clin.notna() & d.amplified.notna()]
    ct = pd.crosstab(both.amplified.astype(bool), both.clin)
    ex.append({"what": "amplified among clinical HER2 positive", "value": ct.loc[True, "positive"] / ct["positive"].sum(),
               "n": int(ct["positive"].sum())})
    ex.append({"what": "amplified among clinical HER2 negative", "value": ct.loc[True, "negative"] / ct["negative"].sum(),
               "n": int(ct["negative"].sum())})
    ex = pd.DataFrame(ex)
    ex.to_csv(OUT / "exploratory.csv", index=False, float_format="%.4g")
    log("exploratory:\n" + ex.to_string(index=False))
    # amplification by PAM50 class; group medians of the two probabilities (cells under 5 left blank)
    g = d[d.amplified.notna()].assign(amplified=lambda x: x.amplified.astype(bool)).groupby(["label", "amplified"])
    s = g.agg(n=("p_gex", "size"), median_p_gex=("p_gex", "median"), median_p_as=("p_as", "median")).reset_index()
    s.loc[s.n < 5, ["median_p_gex", "median_p_as"]] = np.nan
    s.to_csv(OUT / "her2_probability_by_class_and_amplification.csv", index=False, float_format="%.4g")
    log(f"by class and amplification:\n{s.to_string(index=False)}")
    erbb2_events(d)
    knockout(cn)


def knockout(cn: pd.DataFrame) -> None:
    D = DATA / "subtype"
    tab = pd.read_parquet(D / "tasks" / "subtype__BRCA_PAM50.parquet")
    a = tab.sample_id.str[:15].map(cn.amplified)
    lum = tab.label.isin(["LumA", "LumB"]).to_numpy() & a.notna().to_numpy()
    amp = a.fillna(False).astype(bool).to_numpy()[lum]
    rng = np.random.default_rng(SEED)
    rows = []
    for v in AS_VIEWS + ["gex"]:
        F = pd.read_parquet(D / "matrices" / v / "features.parquet").set_index("feature_id")
        P = {k: np.full(len(tab), np.nan) for k in ["base", "core", "wide", "erbb2", "random"]}
        n_core = []
        for k in range(5):
            fd = D / "runs" / "subtype" / v / "BRCA_PAM50" / f"fold{k}"
            prm = json.load(open(fd / "done.json"))
            z = np.load(fd / "shap.npz", allow_pickle=True)
            X, te, f = z["X"].astype(np.float32), z["rows"], F.reindex(z["feature_ids"])
            net = MLP(X.shape[1], 5, prm["param_n_layers"], prm["param_width"], prm["param_dropout"], prm["param_input_dropout"])
            net.load_state_dict(torch.load(fd / "model.pt", map_location="cpu"))
            if v == "gex":
                sets = {"core": f.symbol.isin(CORE_GENES), "wide": f.symbol.isin(WIDE_GENES), "erbb2": f.symbol.eq("ERBB2")}
            else:
                r = f.region.str.extract(r"^chr17:(\d+)-(\d+)$").astype(float)
                inwin = lambda w: f.chrom.eq("chr17") & (r[1] >= w[0]) & (r[0] <= w[1])
                sets = {"core": inwin(CORE), "wide": inwin(WIDE), "erbb2": f.symbol.eq("ERBB2")}
            sets = {s: m.to_numpy() for s, m in sets.items()}
            n_core.append(int(sets["core"].sum()))

            def her2(mask):
                Z = X.copy()
                Z[:, mask] = 0.0
                return logits_to_proba(predict_logits(net, Z), False)[:, 2]

            P["base"][te] = her2(np.zeros(X.shape[1], bool))
            for s, m in sets.items():
                P[s][te] = her2(m)
            ctrl = []
            for _ in range(20):
                m = np.zeros(X.shape[1], bool)
                m[rng.choice(X.shape[1], max(1, n_core[-1]), replace=False)] = True
                ctrl.append(her2(m))
            P["random"][te] = np.median(ctrl, axis=0)
        rec = {"view": v, "amplicon_inputs_per_fold": "/".join(map(str, n_core))}
        for s, p in P.items():
            q = p[lum]
            rec[f"auroc_{s}"] = auc(q[amp], q[~amp])
            rec[f"median_p_amplified_{s}"] = float(np.median(q[amp]))
        rows.append(rec)
    r = pd.DataFrame(rows)
    r.to_csv(OUT / "knockout.csv", index=False, float_format="%.4g")
    log(f"knockout (LumA/LumB, amplified against not):\n{r[['view', 'amplicon_inputs_per_fold', 'auroc_base', 'auroc_core', 'auroc_wide', 'auroc_random']].round(3).to_string(index=False)}")


def erbb2_events(d: pd.DataFrame) -> None:
    """Every ERBB2 event of the main AS matrices: AUROC of PSI, amplified against not, all BRCA tumours."""
    cfg = load_config()
    rows = []
    amp = d.dropna(subset=["amplified"]).set_index("case_id")["amplified"].astype(bool)
    for v in AS_VIEWS:
        X, F, tab = cohort_matrix(cfg, v, "tumor_status", "BRCA")
        keep = tab.role.eq("tumour").to_numpy() & tab.case_id.isin(amp.index).to_numpy()
        y = tab.case_id[keep].map(amp).to_numpy()
        for j in np.flatnonzero(F.symbol.eq("ERBB2").to_numpy()):
            x = X[keep, j]
            ok = np.isfinite(x)
            if ok.sum() < 50 or y[ok].sum() < 10:
                continue
            rows.append({"view": v, "region": F.region.iat[j], "n_amplified": int(y[ok].sum()), "n_not": int((~y[ok]).sum()),
                         "auroc_amplified_higher": auc(x[ok][y[ok]], x[ok][~y[ok]]),
                         "median_amplified": float(np.median(x[ok][y[ok]])), "median_not": float(np.median(x[ok][~y[ok]]))})
    r = pd.DataFrame(rows).sort_values("auroc_amplified_higher")
    r.to_csv(OUT / "erbb2_events_by_amplification.csv", index=False, float_format="%.4g")
    log(f"ERBB2 events:\n{r.to_string(index=False)}")


if __name__ == "__main__":
    main()
