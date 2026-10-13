"""Step 55 (D-72): the splicing model against the exome MSI callers.

D-71's all-exon logistic regression (M5) reads only splicing; MANTIS and MSIsensor (D-69) and the
MC3 homopolymer-deletion count that defines our label read the exome. None is a gold standard, so
each DNA method is judged against a reference built from the other two, and M5 against the same
reference on the same tumours, by within-cancer AUROC. config/followup.yaml: msi_tool_benchmark.

Inputs:  data.nosync/ctsm_v1/msi/models_{tumours,e1_predictions,e2_predictions}.parquet (step 54)
         data.nosync/ctsm_v1/external/cbioportal_msi/tcga_pancan_msi_scores.tsv (D-69)
Outputs: results/subtype/msi_cis/code/benchmark/ (aggregates only; per-tumour scores stay local)
Usage: python scripts/55_msi_benchmark.py [--e2 FILE]
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
LOCAL, OUT = m51.LOCAL, m51.OUT / "benchmark"
THRESH = {"MC3": ("n_hp", 10, ">="), "MANTIS": ("mantis", 0.6, ">"), "MSIsensor": ("msisensor", 10.0, ">")}
LOW = {"MANTIS": 0.4, "MSIsensor": 4.0}              # MSS below these (D-69)
DNA = ["MC3", "MANTIS", "MSIsensor"]
RNA = {"E2": ["M0", "M2", "M4", "M5"], "E1": ["M0", "M2", "M4", "M5", "M6"]}
TESTS = [("BM1", "E2", "MANTIS"), ("BM2", "E2", "MSIsensor"), ("BM3", "E2", "MC3"),
         ("BM4", "E1", "MANTIS"), ("BM5", "E1", "MSIsensor"), ("BM6", "E1", "MC3")]
MARGIN = 0.05


def calls(d: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """(positive, negative) masks of each DNA method's own call."""
    return {"MC3": (d.n_hp.to_numpy() >= 10, d.n_hp.to_numpy() < 10),
            "MANTIS": (d.mantis.to_numpy() > 0.6, d.mantis.to_numpy() < LOW["MANTIS"]),
            "MSIsensor": (d.msisensor.to_numpy() > 10, d.msisensor.to_numpy() < LOW["MSIsensor"])}


def references(d: pd.DataFrame) -> dict[str, np.ndarray]:
    """Reference label per tumour: 1 positive, 0 negative, -1 left out. R_<tool> excludes <tool>."""
    c = calls(d)
    out = {}
    for judged in DNA:
        a, b = [k for k in DNA if k != judged]
        lab = np.full(len(d), -1)
        lab[c[a][0] & c[b][0]] = 1
        lab[c[a][1] & c[b][1]] = 0
        out[judged] = lab
    lab = np.full(len(d), -1)
    lab[c["MC3"][0] & c["MANTIS"][0] & c["MSIsensor"][0]] = 1
    lab[c["MC3"][1] & c["MANTIS"][1] & c["MSIsensor"][1]] = 0
    out["all"] = lab
    return out


def auc_within(score: np.ndarray, y: np.ndarray, grp: np.ndarray) -> float:
    """P(positive outranks a negative of the same group), ties half, over all within-group pairs."""
    won = pairs = 0.0
    for g in np.unique(grp):
        k = grp == g
        yk = y[k]
        n1 = yk.sum()
        n0 = len(yk) - n1
        if n1 == 0 or n0 == 0:
            continue
        r = rankdata(score[k])
        won += r[yk].sum() - n1 * (n1 + 1) / 2
        pairs += n1 * n0
    return float(won / pairs) if pairs else np.nan


def auc_pooled(score, y) -> float:
    r = rankdata(score)
    n1 = y.sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * (len(y) - n1)))


def paired_bootstrap(Sk: dict, y: np.ndarray, g: np.ndarray, rng) -> dict:
    """Within-cancer AUROC of each method over BOOT resamples of tumours within cancer x label strata
    (the same draws for every method). Equal to auc_within on each resample, computed per cancer."""
    cancers = np.unique(g)
    strata = [(a, b, np.flatnonzero((g == a) & (y == b))) for a in cancers for b in (False, True)]
    strata = [x for x in strata if len(x[2])]
    both = [a for a in cancers if (y[g == a].any() and (~y[g == a]).any())]
    out = {m: np.empty(m51.BOOT) for m in Sk}
    for k in range(m51.BOOT):
        draw = {(a, b): rng.choice(s, len(s)) for a, b, s in strata}
        for m, sm in Sk.items():
            won = pairs = 0.0
            for a in both:
                pi, ni = draw[(a, True)], draw[(a, False)]
                r = rankdata(np.concatenate([sm[pi], sm[ni]]))
                n1 = len(pi)
                won += r[:n1].sum() - n1 * (n1 + 1) / 2
                pairs += n1 * len(ni)
            out[m][k] = won / pairs
    return out


def ap(score, y) -> float:
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(y, score))


def matched_volume(score: np.ndarray, y: np.ndarray, k: int) -> tuple[float, float]:
    """Sensitivity and PPV of flagging the top k by score (tied negatives ranked first)."""
    if k == 0:
        return np.nan, np.nan
    order = np.lexsort((y.astype(int), -score))
    tp = int(y[order[:k]].sum())
    return tp / y.sum(), tp / k


def data() -> pd.DataFrame:
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")[["file_id", "sample15", "project", "n_hp", "mmr"]]
    sc = pd.read_csv(m51.MSI_TSV, sep="\t").drop_duplicates("sample_id").set_index("sample_id")[["mantis", "msisensor"]]
    t = t.join(sc, on="sample15")
    return t


def main() -> None:
    import argparse
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--e2", default="models_e2_predictions.parquet",
                     help="E2 predictions file in the local MSI folder (default: step 54's)")
    e2_file = ap_.parse_args().e2
    t = data()
    n_all = len(t)
    cov = [{"what": "coverage: share of D-71's tumours scored", "method": m, "n": int(v.sum()), "value": v.mean()}
           for m, v in [("MC3 count", t.n_hp.notna()), ("MANTIS", t.mantis.notna()), ("MSIsensor", t.msisensor.notna()),
                        ("both exome callers", t.mantis.notna() & t.msisensor.notna()), ("RNA (M0-M5)", pd.Series(True, index=t.index))]]
    common = t.mantis.notna() & t.msisensor.notna()
    log(f"{n_all:,} tumours; common set {int(common.sum()):,}")
    pred, met, op, pc, disc = [], [], [], [], []
    rng = np.random.default_rng(m51.SEED)
    for split in ["E2", "E1"]:
        p = pd.read_parquet(LOCAL / (e2_file if split == "E2" else "models_e1_predictions.parquet"))
        d = t[common].merge(p.drop(columns=["project", "mmr"]), on="file_id", how="inner").reset_index(drop=True)
        rna = [m for m in RNA[split] if m in d.columns]
        S = {"MC3": d.n_hp.to_numpy(float), "MANTIS": d.mantis.to_numpy(float), "MSIsensor": d.msisensor.to_numpy(float),
             **{m: d[m].to_numpy(float) for m in rna}}
        grp = d.project.to_numpy()
        refs = references(d)
        for ref, lab in refs.items():
            keep = lab >= 0
            y, g = lab[keep] == 1, grp[keep]
            methods = ([ref] if ref in DNA else []) + rna      # the judged DNA method, then the RNA models
            Sk = {m: S[m][keep] for m in methods}
            obs = {m: auc_within(Sk[m], y, g) for m in methods}
            boot = paired_bootstrap(Sk, y, g, rng)
            for m in methods:
                met.append({"split": split, "reference": f"R_{ref}", "method": m, "positives": int(y.sum()),
                            "negatives": int((~y).sum()), "cancers_with_positives": int(len(np.unique(g[y]))),
                            "auroc_within": obs[m], "lo": np.quantile(boot[m], 0.025), "hi": np.quantile(boot[m], 0.975),
                            "auroc_pooled": auc_pooled(Sk[m], y), "auprc_pooled": ap(Sk[m], y)})
            if ref in DNA:
                for m in rna:
                    dd = boot[m] - boot[ref]
                    lo, hi = np.quantile(dd, [0.025, 0.975])
                    test = [tid for tid, sp, j in TESTS if sp == split and j == ref and m == "M5"]
                    row = {"what": f"{split}, R_{ref}: within-cancer AUROC({m}) - AUROC({ref})", "effect": obs[m] - obs[ref],
                           "ci_low": lo, "ci_high": hi, "non_inferior": bool(lo > -MARGIN),
                           "n": f"{int(y.sum())} positives, {int((~y).sum())} negatives; {m} {obs[m]:.3f}, {ref} {obs[ref]:.3f}"}
                    if test:
                        pred.append({"id": test[0], **row, "predicted_sign": -1, "p": m51.boot_p(dd)})
                    else:
                        met.append({"split": split, "reference": f"R_{ref}", "method": f"{m} minus {ref}",
                                    "auroc_within": obs[m] - obs[ref], "lo": lo, "hi": hi})
                # matched-volume operating point
                col, thr, op_ = THRESH[ref]
                flag = Sk[ref] >= thr if op_ == ">=" else Sk[ref] > thr
                kf = int(flag.sum())
                tp = int((flag & y).sum())
                op.append({"split": split, "reference": f"R_{ref}", "method": ref, "flagged": kf,
                           "sensitivity": tp / y.sum(), "ppv": tp / kf if kf else np.nan})
                for m in rna:
                    se, pp = matched_volume(Sk[m], y, kf)
                    op.append({"split": split, "reference": f"R_{ref}", "method": m, "flagged": kf, "sensitivity": se, "ppv": pp})
                # per cancer
                for c in np.unique(g):
                    k = g == c
                    if y[k].sum() >= 3:
                        for m in methods:
                            pc.append({"split": split, "reference": f"R_{ref}", "project": c, "method": m,
                                       "positives": int(y[k].sum()), "n": int(k.sum()),
                                       "auroc": auc_pooled(Sk[m][k], y[k])})
            log(f"{split} R_{ref}: {int(y.sum())} positives, {int((~y).sum())} negatives; "
                + ", ".join(f"{m} {obs[m]:.3f}" for m in methods))
        # discordance: M4's and M5's within-cancer percentiles over all of the split's tumours (M4: A-19)
        full = t.merge(p.drop(columns=["project", "mmr"]), on="file_id", how="inner")
        for m in ["M4", "M5"]:
            full[f"pct_{m}"] = full.groupby("project")[m].rank(pct=True)
        cm = full.mantis > 0.6
        cs = full.msisensor > 10
        both_scored = full.mantis.notna() & full.msisensor.notna()
        groups = {"MC3 MMR-deficient, both callers MSI-H": full.mmr & cm & cs,
                  "MC3 MMR-deficient, one caller MSI-H": full.mmr & both_scored & (cm ^ cs),
                  "MC3 MMR-deficient, neither caller MSI-H": full.mmr & both_scored & ~cm & ~cs,
                  "MC3 MMR-deficient, a caller score missing": full.mmr & ~both_scored,
                  "MC3 proficient, both callers MSI-H": ~full.mmr & cm & cs,
                  "MC3 proficient, neither caller MSI-H": ~full.mmr & both_scored & ~cm & ~cs}
        for name, k in groups.items():
            n = int(k.sum())
            row = {"split": split, "group": name, "n": n}
            for m in ["M4", "M5"]:
                row[f"median_pct_{m}"] = float(full.loc[k, f"pct_{m}"].median()) if n >= 5 else np.nan
                row[f"share_top_decile_{m}"] = float((full.loc[k, f"pct_{m}"] > 0.9).mean()) if n >= 5 else np.nan
            disc.append(row)
    OUT.mkdir(parents=True, exist_ok=True)
    pr = m51.predictions(pred, "benchmark")
    pd.DataFrame(met).to_csv(OUT / "metrics.csv", index=False, float_format="%.4g")
    pd.DataFrame(op).to_csv(OUT / "operating_points.csv", index=False, float_format="%.4g")
    pd.DataFrame(pc).to_csv(OUT / "per_cancer.csv", index=False, float_format="%.4g")
    pd.DataFrame(disc).to_csv(OUT / "discordance.csv", index=False, float_format="%.4g")
    pd.DataFrame(cov).to_csv(OUT / "coverage.csv", index=False, float_format="%.4g")
    log(f"metrics:\n{pd.DataFrame(met).round(3).to_string(index=False)}")
    log(f"operating points:\n{pd.DataFrame(op).round(3).to_string(index=False)}")
    log(f"discordance:\n{pd.DataFrame(disc).round(3).to_string(index=False)}")
    log(f"coverage:\n{pd.DataFrame(cov).round(4).to_string(index=False)}")
    log(f"non-inferiority (margin {MARGIN}):\n{pr[['id', 'what', 'ci_low', 'non_inferior']].to_string(index=False)}")


if __name__ == "__main__":
    main()
