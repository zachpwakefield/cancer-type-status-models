"""Step 09 (B-02): global, directional RNA-processing shifts between tumour and matched normal.

Uses ALL measured events of each cohort (no SHAP selection, so the tests are not circular).
  PSI views   per event: delta = median(tumour) - median(normal); fraction of events shifting
              up / down by > min_delta. Per patient (matched pairs): mean paired dPSI over the
              events observed in both samples; Wilcoxon signed-rank across patients.
              Meaning of "up": RI = more intron retention; A5SS/A3SS = more use of the long exon
              form (intron-proximal splice site, extra intronic sequence); SE / MXE = more
              inclusion of the cassette / transcript-upstream exon.
  AFE / ALE   gene-level terminal-exon position index per sample = usage-weighted rank of the
              used first / last exons along the direction of transcription (0 = most upstream,
              1 = most downstream); paired tumour - normal change per patient.
Also: does per-sample intron retention track a proliferation score within tumours?

Outputs: results/biology/02_global_shifts/ (tables; figures: scripts/15_bio_figures.py; per-patient values stay in <data_dir>/bio_cache, D-12)
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import spearmanr, wilcoxon

from ctsm.biology import bh, bio_config, cohort_matrix, gene_strand, gex_by_symbol, matched_pairs, terminal_index
from ctsm.config import load_config

warnings.filterwarnings("ignore", category=RuntimeWarning)
PSI_VIEWS = ["ri", "a5ss", "a3ss", "se", "mxe"]
PROLIF = ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "CENPF", "TPX2", "NUSAP1", "UBE2C"]


def main() -> None:
    cfg = load_config()
    bc = bio_config()["shifts"]
    out = cfg.results_dir / "biology" / "02_global_shifts"
    out.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = list(scope.loc[scope.tumor_status_in_scope, "project"])
    strand = gene_strand(cfg)
    rows, patient_rows, prolif_rows = [], [], []
    for c in cancers:
        for view in PSI_VIEWS + ["ale", "afe"]:
            X, F, tab = cohort_matrix(cfg, view, "tumor_status", c)
            if view in ("ale", "afe"):
                X = terminal_index(X, F, strand, bc["terminal_min_usage"])
            ti, ni, cases = matched_pairs(tab)
            yv = tab["y"].to_numpy()
            delta = np.nanmedian(X[yv == 1], axis=0) - np.nanmedian(X[yv == 0], axis=0)
            ok = np.isfinite(delta)
            d_pat = np.nanmean(X[ti] - X[ni], axis=1)
            p = wilcoxon(d_pat).pvalue if len(d_pat) >= 5 else np.nan
            rows.append({"cancer": c, "view": view, "n_events": int(ok.sum()), "n_pairs": len(ti),
                         "frac_up": float(np.mean(delta[ok] > bc["min_delta"])),
                         "frac_down": float(np.mean(delta[ok] < -bc["min_delta"])),
                         "median_paired_shift": float(np.median(d_pat)), "p_wilcoxon": float(p),
                         "frac_patients_up": float(np.mean(d_pat > 0))})
            for cs, dp in zip(cases, d_pat):
                patient_rows.append({"cancer": c, "view": view, "case_id": cs, "paired_shift": float(dp)})
            if view == "ri":   # exploratory: retention vs proliferation within tumours
                tum = np.flatnonzero(yv == 1)
                ir = np.nanmean(X[tum], axis=1)
                gx = gex_by_symbol(cfg, tab["file_id"].iloc[tum], PROLIF)
                pscore = gx.mean(axis=1).to_numpy()
                rho = spearmanr(ir, pscore).statistic
                prolif_rows.append({"cancer": c, "n_tumours": len(tum), "spearman_IR_vs_proliferation": float(rho)})
        log(f"{c} done")
    res = pd.DataFrame(rows)
    res["q_wilcoxon"] = bh(res["p_wilcoxon"].fillna(1).to_numpy())
    res["log2_up_down"] = np.log2((res["frac_up"] + 1e-3) / (res["frac_down"] + 1e-3))
    res.to_csv(out / "shift_summary.csv", index=False)
    pats = pd.DataFrame(patient_rows)
    pats.to_csv(cfg.data_dir / "bio_cache" / "paired_shift_per_patient.csv", index=False)   # individual-level: local (D-12)
    pr = pd.DataFrame(prolif_rows)
    pr.to_csv(out / "intron_retention_vs_proliferation.csv", index=False)

    log("global shifts done")


if __name__ == "__main__":
    main()
