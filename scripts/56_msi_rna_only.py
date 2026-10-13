"""Step 56 (D-73): tumours the repeat-exon model calls MSI that no exome measure calls.

The repeat-exon logistic model (D-71 M4; A-19) calls a tumour MSI when its leave-one-cancer-out
decision function is above 0. Among tumours with both caller scores, 'RNA-only' calls are those
that neither MC3 (>= 10 homopolymer deletions), MANTIS (> 0.6) nor MSIsensor (> 10) calls. RO1-RO5
ask whether they look like missed mismatch-repair deficiency (low MMR-gene expression, sub-threshold
instability) or like false calls. config/followup.yaml: msi_rna_only.

Inputs:  data.nosync/ctsm_v1/msi/models_{tumours,e2_predictions}.parquet (step 54), the D-69 caller
         table, MC3, ABSOLUTE purity and gene expression (all local)
Outputs: results/subtype/msi_cis/code/rna_only/ (aggregates); data.nosync/ctsm_v1/msi/rna_only_tumours.parquet
         and rna_dna_calls.parquet (every tumour, for figure U25)
Usage: python scripts/56_msi_rna_only.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import rankdata


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m51 = _load("msi51", "51_msi_code.py")
m50 = m51.m50
LOCAL, OUT = m51.LOCAL, m51.OUT / "rna_only"
MMR_GENES = ["MLH1", "MSH2", "MSH6", "PMS2"]
MUT_GENES = MMR_GENES + ["POLE"]
NON_SILENT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
              "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation"}
N_PERM, MIN_N = 10_000, 5


def mmr_mutations(samples: set[str]) -> pd.DataFrame:
    """Per 15-character barcode (first aliquot), whether MC3 holds a non-silent call passing the D-61
    flags in each of MUT_GENES."""
    cols = ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "FILTER"]
    m = pd.read_csv(m50.EXT / "gdc_pancanatlas" / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=cols, dtype=str)
    m["sample15"] = m.Tumor_Sample_Barcode.str[:15]
    m = m[m.sample15.isin(samples)]
    first = m.sort_values("Tumor_Sample_Barcode").drop_duplicates("sample15").Tumor_Sample_Barcode
    m = m[m.Tumor_Sample_Barcode.isin(first) & m.FILTER.str.split(",").map(lambda f: set(f) <= m50.KEEP_FLAGS)]
    m = m[m.Hugo_Symbol.isin(MUT_GENES) & m.Variant_Classification.isin(NON_SILENT)]
    out = pd.DataFrame(index=sorted(samples))
    for g in MUT_GENES:
        out[f"mut_{g}"] = out.index.isin(m[m.Hugo_Symbol.eq(g)].sample15)
    return out


def within_pct(v: pd.Series, grp: pd.Series) -> pd.Series:
    """Percentile within the tumour's cancer (missing values stay missing)."""
    return v.groupby(grp).rank(pct=True)


def auroc_perm(x: np.ndarray, case: np.ndarray, grp: np.ndarray, rng, lower: bool) -> tuple[float, float, int, int]:
    """AUROC of case against non-case on x (higher, or lower if `lower`), and a two-sided permutation p
    with case labels shuffled within each group. Rows with missing x are dropped."""
    ok = np.isfinite(x)
    x, case, grp = (-x if lower else x)[ok], case[ok], grp[ok]
    n1, n0 = int(case.sum()), int((~case).sum())
    if n1 < MIN_N:
        return np.nan, np.nan, n1, n0

    r = rankdata(x)

    def auc(c):
        return (r[c].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    obs = auc(case)
    idx = [np.flatnonzero(grp == g) for g in np.unique(grp)]
    hits = 0
    for _ in range(N_PERM):
        c = case.copy()
        for i in idx:
            c[i] = rng.permutation(c[i])
        hits += abs(auc(c) - 0.5) >= abs(obs - 0.5) - 1e-12
    return float(obs), (hits + 1) / (N_PERM + 1), n1, n0


def main() -> None:
    from ctsm.biology import gex_by_symbol
    from ctsm.config import load_config
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")[["file_id", "case_id", "sample15", "project", "n_hp",
                                                             "n_snv", "pct_G"]]
    p = pd.read_parquet(LOCAL / "models_e2_predictions.parquet")[["file_id", "M4"]]
    t = t.merge(p, on="file_id")
    sc = pd.read_csv(m51.MSI_TSV, sep="\t").drop_duplicates("sample_id").set_index("sample_id")[["mantis", "msisensor"]]
    t = t.join(sc, on="sample15")
    t["rna"] = t.M4 > 0
    t["dna_mc3"], t["dna_mantis"], t["dna_msisensor"] = t.n_hp >= 10, t.mantis > 0.6, t.msisensor > 10
    t["n_dna"] = t[["dna_mc3", "dna_mantis", "dna_msisensor"]].sum(axis=1)
    t["both_scored"] = t.mantis.notna() & t.msisensor.notna()
    g = gex_by_symbol(load_config(), t.file_id, MMR_GENES)
    for gene in MMR_GENES:
        t[f"{gene}_pct"] = within_pct(t.file_id.map(g[gene]), t.project)
    t["min_mmr_pct"] = t[[f"{gene}_pct" for gene in MMR_GENES]].min(axis=1)
    t["purity"] = t.case_id.map(m50.purity_by_patient())
    t = t.join(mmr_mutations(set(t.sample15)), on="sample15")
    for c in ["mantis", "msisensor", "n_hp", "purity", "n_snv"]:
        t[f"{c}_pct"] = within_pct(t[c].astype(float), t.project)
    b = t[t.both_scored]
    t["group"] = np.select(
        [t.both_scored & t.rna & t.n_dna.eq(0), t.both_scored & ~t.rna & t.n_dna.eq(0),
         t.both_scored & t.rna & t.n_dna.ge(1), t.both_scored & ~t.rna & t.n_dna.ge(1), ~t.both_scored],
        ["rna_only", "reference", "rna_and_dna", "dna_only", "caller score missing"], "")
    log(f"{len(t):,} tumours; RNA calls {int(t.rna.sum())}; both caller scores {len(b):,}")
    OUT.mkdir(parents=True, exist_ok=True)
    # counts
    counts = (t.assign(rna_call=np.where(t.rna, "RNA call", "no RNA call"),
                       dna=np.where(t.both_scored, t.n_dna.astype(str) + " DNA calls", "caller score missing"))
              .groupby(["rna_call", "dna"]).size().rename("n").reset_index())
    counts.to_csv(OUT / "counts.csv", index=False)
    log(f"counts:\n{counts.to_string(index=False)}")
    ro = t[t.group.eq("rna_only")]
    by = ro.groupby("project").size().rename("rna_only").to_frame()
    by["rna_calls"] = t[t.rna].groupby("project").size()
    by["tumours"] = t.groupby("project").size()
    by = by.reset_index().sort_values("rna_only", ascending=False)
    by.to_csv(OUT / "by_cancer.csv", index=False)
    log(f"RNA-only by cancer:\n{by.to_string(index=False)}")
    # tests: rna_only against reference, cancers holding >= 1 rna_only tumour
    cmp_ = t[t.group.isin(["rna_only", "reference"]) & t.project.isin(ro.project.unique())]
    case, grp = cmp_.group.eq("rna_only").to_numpy(), cmp_.project.to_numpy()
    rng = np.random.default_rng(m51.SEED)
    pred, ex = [], []
    for pid, col, lower, what in [("RO1", "MLH1_pct", True, "MLH1 expression lower"),
                                  ("RO2", "min_mmr_pct", True, "lowest MMR-gene expression percentile lower"),
                                  ("RO3", "mantis_pct", False, "MANTIS score higher"),
                                  ("RO4", "msisensor_pct", False, "MSIsensor score higher"),
                                  ("RO5", "n_hp_pct", False, "MC3 homopolymer deletions higher")]:
        a, pv, n1, n0 = auroc_perm(cmp_[col].to_numpy(float), case, grp, rng, lower)
        pred.append({"id": pid, "what": f"AUROC({what} in RNA-only than reference, within-cancer percentile) - 0.5",
                     "effect": a - 0.5 if np.isfinite(a) else np.nan, "predicted_sign": 1, "p": pv,
                     "n": f"RNA-only {n1}, reference {n0}" + ("" if n1 >= MIN_N else f"; fewer than {MIN_N}")})
    pr = m50.add_verdicts(pred)
    pr.to_csv(OUT / "predictions.csv", index=False)
    log("rna_only predictions:\n" + pr[["id", "effect", "p", "q", "verdict", "n"]].to_string(index=False))
    for col, lower, what in [("purity_pct", True, "purity lower"), ("n_snv_pct", False, "SNV count higher"),
                             ("MSH2_pct", True, "MSH2 lower"), ("MSH6_pct", True, "MSH6 lower"),
                             ("PMS2_pct", True, "PMS2 lower")]:
        a, pv, n1, n0 = auroc_perm(cmp_[col].to_numpy(float), case, grp, rng, lower)
        ex.append({"what": f"AUROC({what} in RNA-only than reference) [exploratory, no q]", "group": "rna_only",
                   "n": n1, "value": a, "p_perm": pv})
    # group summaries (medians blanked below 5)
    for name, d in t.groupby("group"):
        n = len(d)
        row = {"group": name, "n": n}
        for c in ["M4", "pct_G", "MLH1_pct", "min_mmr_pct", "mantis", "msisensor", "n_hp", "purity", "n_snv"]:
            row[f"median_{c}"] = float(d[c].median()) if n >= MIN_N else np.nan
        for gene in MUT_GENES:
            row[f"n_mut_{gene}"] = int(d[f"mut_{gene}"].sum())
        row["n_MLH1_below_10th_pct"] = int((d.MLH1_pct < 0.1).sum())
        ex.append({"what": "group summary", **row})
    exd = pd.DataFrame(ex)
    exd.to_csv(OUT / "exploratory.csv", index=False, float_format="%.4g")
    log(f"exploratory:\n{exd.round(3).to_string(index=False)}")
    keep = ["file_id", "sample15", "project", "group", "M4", "pct_G", "n_hp", "n_snv", "mantis", "msisensor",
            "purity", *[f"{g_}_pct" for g_ in MMR_GENES], "min_mmr_pct", *[f"mut_{g_}" for g_ in MUT_GENES]]
    t[t.rna | t.n_dna.ge(1)][keep].sort_values(["group", "M4"], ascending=[True, False]).to_parquet(
        LOCAL / "rna_only_tumours.parquet", index=False)
    allc = ["dna_mc3", "dna_mantis", "dna_msisensor", "rna", "both_scored", "n_dna"]
    t[[*keep, *allc]].to_parquet(LOCAL / "rna_dna_calls.parquet", index=False)      # all tumours, for U25
    log(f"per-tumour rows written locally: {LOCAL / 'rna_only_tumours.parquet'}, {LOCAL / 'rna_dna_calls.parquet'}")


if __name__ == "__main__":
    main()
