**Figure U10. Genome sequence alone says which exons MSI tumours skip.** Pre-specified before any PSI was read by
repeat (config/followup.yaml: msi_cis, D-52). For every cassette exon, the longest sense-strand T-run overlapping
intron positions −1..−50 before it was read from GRCh38; no label or model chose the exons. MSI against the other
classes (POLE apart) of UCEC_TCGA, STAD_TCGA and COAD_TCGA.

(a) Share of exons with a run ≥ 10 among all exons measured in the cohort and among the top 50 SE events of the
cohort's MSI model (mean |SHAP| over MSI patients); Fisher exact (M6). (b) The repeat score of each patient: each
exon's PSI as a percentile within the cohort, mean over exons with a run ≥ 10 minus the mean over exons with ≤ 5
(weighted to the same mean-PSI deciles, which cancels the global PSI axis of D-24). Box = interquartile range,
diamond = median, whiskers = 1.5 × IQR, dots = patients (n under each box); Mann–Whitney (M1, M5). (c) AUROC of the
same score built from each run-length bin, and (hollow) from T-runs ≥ 10 at −201..−400, deep in the intron, with the
acceptor ≤ 6.

- M6, the trained MSI models chose repeat exons: UCEC 10 of 50 against 4.9% (OR 5, p = 1.4e-4); STAD 40 of 50
  against 4.5% (OR 84, p = 1.8e-44); COAD 34 of 50 against 4.7% (OR 43, p = 2.9e-33).
- M1, the score is lower in MSI: AUROC UCEC 0.51 (p = 0.85); STAD 0.71 (p = 1.1e-6); COAD 0.78 (p = 3.1e-8).
  Endometrium does not separate on this score.
- M4, acceptor window against deep window: AUROC difference UCEC +0.17 (95% CI +0.04 to +0.31); STAD +0.28 (95% CI
  +0.19 to +0.37); COAD +0.41 (95% CI +0.29 to +0.52); supported in all three, endometrium included.
- M5, MSI against POLE (UCEC): AUROC 0.53, p = 0.73, not resolved; POLE against the rest 0.49. With MSI itself at
  0.51 in this cohort, the POLE control cannot discriminate here.
- Not drawn: K1, MRE11 exon 5 (T11) is less included in MSI in all three (AUROC UCEC 0.63, STAD 0.81, COAD 0.76);
  M7, the repeat exons' shifts agree between cohorts (Spearman 0.44, 0.43, 0.65).

Stars: * p < 0.05, ** p < 0.01, *** p < 0.001, ns otherwise. q over family M (19 tests) is in the CSV with every
verdict in `results/subtype/msi_cis/predictions.csv`. Per-patient scores are kept locally (D-12).
