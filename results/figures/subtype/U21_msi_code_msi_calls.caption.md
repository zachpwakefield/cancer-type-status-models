**Figure U21. The code finds the tumours that two exome MSI callers call MSI-high.** Pre-specified
(config/followup.yaml: msi_code_msi_calls, D-69) before the scores were downloaded. MANTIS (Bonneville et al. 2017)
and MSIsensor scores of the TCGA PanCancer Atlas, fetched from cBioPortal (step 00c) and matched by sample barcode.
Calls: MANTIS MSI-H > 0.6, MSS < 0.4; MSIsensor MSI-H > 10, MSS < 4; scores in between (mid) are kept apart.
Tumours: D-66's 7,685 held-out tumours (30 cancers, none used to define the code); score G as there, a percentile
within each cancer, not re-fitted.

(a) AUROC for positive tumours scoring lower on G (positives in brackets): D-66's MC3 call (≥ 10 homopolymer
deletions) against the rest, MANTIS and MSIsensor MSI-H against MSS; then, in-sample for the rule, MSI-H against MSS
in UCEC, STAD and COAD. (b, c) Score G percentile of the held-out tumours by D-66's MC3 call and the caller's call
(n in brackets); box = interquartile range, diamond = median, whiskers = 1.5 × IQR, dots = tumours (not drawn for
the large MSS group); dashed line = the middle of the cancer.

- MC1, MANTIS MSI-H against MSS: AUROC 0.97 (17 MSI-H), q = 6.3e-11, supported.
- MC2, MSIsensor: AUROC 0.89 (23 MSI-H), q = 1.4e-10, supported.
- MC3, Spearman of score G with the MANTIS score (percentile within cancer), all 7,388 tumours: ρ = +0.011, q =
  0.32, not resolved. Over all tumours, nearly all of them MSS, the MANTIS score shows no relation to score G.

Agreement (exploratory): every MANTIS MSI-H tumour (17) is MC3-deficient; of the 50 MC3-deficient tumours, MANTIS
calls 27 mid and 6 MSS, and these still score low on G. MSIsensor calls 4 MC3-proficient tumours MSI-H. q over
family MC (3 tests). Per-tumour values are kept locally (D-12).
