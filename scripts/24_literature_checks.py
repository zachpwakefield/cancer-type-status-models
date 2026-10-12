"""Step 24 (D-28): test two published mechanisms against the events the tumour-status models use.

Each quantity is computed inside every tumour-status cancer and summarised by a sign test across
the 15 cancers. Where a paper fixes the direction the test is one-sided in that direction; where
it does not (ENOSF1) the test is two-sided and asks only whether the sign is consistent. These
are concordance checks, not blind replications: the numbers were first seen in an exploratory
probe while the event cards were written. Tumours only, unless stated.

1. macroH2A1 (MXE; the tracked exon is the one in macroH2A1.2, see step 22).
   QKI promotes the macroH2A1.1 exon (Novikov et al. 2011, Mol Cell Biol 31:4244), so the
   macroH2A1.2 exon should FALL as QKI rises (rho < 0), also once proliferation is held fixed.
   Proliferating cells favour macroH2A1.2 (Novikov 2011; Sporn et al. 2009, Oncogene 28:3423),
   so it should RISE with a proliferation score (rho > 0), and tumours should sit above their
   matched normals (paired tumour - normal > 0). One exploratory addition, made after seeing which
   cancers do not shift: across the 15 cancers, does the shift follow the tumour - normal
   change in proliferation? It is labelled as such in the output.
2. ENOSF1 (HIT). ENOSF1 was first cloned as rTS, an antisense RNA to thymidylate synthase
   (Dolnick 1993, Nucleic Acids Res 21:1747); its rTS-alpha 3' end is complementary to the TYMS
   3' UTR and last intron. If the event partly reads TYMS, it should follow TYMS beyond what the
   shared proliferation signal explains, and still carry some proliferation signal once TYMS is
   held fixed if part of it is the cell-cycle programme itself.
The proliferation score is the mean log2(CPM+1) of the ten genes step 09 uses.

Partial correlations are Spearman: ranks residualised on the ranked covariate by least squares.
Written: results/biology/13_literature_checks/{per_cancer,summary}.csv (aggregates only, D-12).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from _common import REPO, log
from scipy.stats import binomtest, rankdata, spearmanr

from ctsm.biology import feature_values, gex_by_symbol, matched_pairs, task_table
from ctsm.config import load_config

PROLIF = ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "CENPF", "TPX2", "NUSAP1", "UBE2C"]
OUT = REPO / "results" / "biology" / "13_literature_checks"
# (column, expected sign, statement): the direction a source predicts, 0 where none does
CHECKS = [
    ("mH2A1.2_tumour_minus_normal", +1, "macroH2A1.2 exon higher in tumour than matched normal"),
    ("mH2A1.2_vs_QKI", -1, "macroH2A1.2 exon falls as QKI rises"),
    ("mH2A1.2_vs_QKI_given_prolif", -1, "... with proliferation held fixed"),
    ("mH2A1.2_vs_prolif", +1, "macroH2A1.2 exon rises with proliferation"),
    ("ENOSF1_vs_TYMS_given_prolif", 0, "ENOSF1 event follows antisense TYMS beyond proliferation"),
    ("ENOSF1_vs_prolif_given_TYMS", 0, "ENOSF1 event still tracks proliferation with TYMS fixed"),
]


def rho(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    return float(spearmanr(x[ok], y[ok]).statistic) if ok.sum() >= 10 else np.nan


def partial(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < 10:
        return np.nan
    rx, ry, rz = (rankdata(v[ok]) for v in (x, y, z))
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return float(np.corrcoef(ex, ey)[0, 1])


def main() -> None:
    cfg = load_config()
    OUT.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = list(scope.loc[scope.tumor_status_in_scope, "project"])
    cards = pd.read_csv(REPO / "results" / "biology" / "11_event_cards" / "tumour_status_cards.csv")
    fid = cards.set_index("symbol")["feature_id"]
    rows = []
    for c in cancers:
        tab = task_table(cfg, "tumor_status", c)
        files = tab["file_id"]
        mh = feature_values(cfg, "mxe", files, [fid["MACROH2A1"]])[:, 0].astype(float)
        en = feature_values(cfg, "hit", files, [fid["ENOSF1"]])[:, 0].astype(float)
        G = gex_by_symbol(cfg, files, PROLIF + ["QKI", "TYMS"])
        prol = G[[g for g in PROLIF if g in G]].mean(axis=1).to_numpy()
        q, ty = G["QKI"].to_numpy(), G["TYMS"].to_numpy()
        ti, ni, _ = matched_pairs(tab)
        d = mh[ti] - mh[ni]
        t = tab["role"].eq("tumour").to_numpy()
        rows.append({
            "cancer": c, "n_tumours": int(t.sum()), "n_pairs": int(np.isfinite(d).sum()),
            "mH2A1.2_tumour_minus_normal": float(np.nanmedian(d)) if np.isfinite(d).any() else np.nan,
            "prolif_tumour_minus_normal": float(np.nanmedian(prol[ti] - prol[ni])),
            "mH2A1.2_vs_QKI": rho(mh[t], q[t]),
            "mH2A1.2_vs_QKI_given_prolif": partial(mh[t], q[t], prol[t]),
            "mH2A1.2_vs_prolif": rho(mh[t], prol[t]),
            "QKI_vs_prolif": rho(q[t], prol[t]),
            "ENOSF1_vs_TYMS": rho(en[t], ty[t]),
            "ENOSF1_vs_TYMS_given_prolif": partial(en[t], ty[t], prol[t]),
            "ENOSF1_vs_prolif": rho(en[t], prol[t]),
            "ENOSF1_vs_prolif_given_TYMS": partial(en[t], prol[t], ty[t]),
            "TYMS_vs_prolif": rho(ty[t], prol[t]),
        })
        log(f"{c} done")
    per = pd.DataFrame(rows)
    per.round(4).to_csv(OUT / "per_cancer.csv", index=False)
    summ = []
    for col, sign, text in CHECKS:
        v = per[col].dropna().to_numpy()
        if sign:
            k = int((np.sign(v) == sign).sum())
            p = binomtest(k, len(v), 0.5, alternative="greater").pvalue
        else:
            k = int(max((v > 0).sum(), (v < 0).sum()))
            p = binomtest(k, len(v), 0.5, alternative="two-sided").pvalue
        summ.append({"quantity": col, "expected_sign": {1: "+", -1: "-", 0: "none"}[sign],
                     "prediction": text, "n_cancers": len(v),
                     "n_in_expected_or_majority_direction": k,
                     "median_over_cancers": float(np.median(v)), "p_sign_test": float(p)})
    # exploratory, added after seeing which cancers do not shift: is the tumour - normal move
    # toward macroH2A1.2 larger where tumours out-proliferate their normals by more?
    across = spearmanr(per["mH2A1.2_tumour_minus_normal"], per["prolif_tumour_minus_normal"])
    summ.append({"quantity": "exploratory: across cancers, mH2A1.2 shift vs proliferation shift",
                 "expected_sign": "+", "prediction": "Spearman over the 15 cancers (post hoc)",
                 "n_cancers": len(per), "n_in_expected_or_majority_direction": np.nan,
                 "median_over_cancers": float(across.statistic), "p_sign_test": float(across.pvalue)})
    summ = pd.DataFrame(summ)
    summ.round(4).to_csv(OUT / "summary.csv", index=False)
    log("\n" + summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
