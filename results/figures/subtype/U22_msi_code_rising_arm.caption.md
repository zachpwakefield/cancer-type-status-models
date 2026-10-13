**Figure U22. The rising arm is specific to mismatch repair and needs the tract, but its mechanism is not settled.**
Pre-specified (config/followup.yaml: msi_code_reversal, D-70) before any test was computed; panel (c) is exploratory
(A-16), added after the tests were read. Exon shift = D-66's background-adjusted difference in inclusion percentile,
MMR-deficient minus other tumours, pooled over the 30 held-out cancers (×100 = percentile points). R_up = acceptor
T-run ≥ 20 nt with the sequence-predicted branch point ≥ 36 nt from the AG; R_down = runs of 10–19 nt; B_far = no
run (≤ 5) with the branch point ≥ 36 nt; shield = a run ≥ 20 nt lying upstream of a branch point, so that
contracting it cannot move the branch point. R_up arm = each tumour's mean R_up percentile against matched exons
without a run, as a percentile within its cancer.

(a) Mean R_up shift in 20 against 20 tumours (50 draws): every pair of subtypes in the 22 cohorts with MSI, POLE and
MC3-deficient tumours removed (grey, absolute value; dashed = their 95th percentile), the held-out MMR-deficient
tumours each paired with a proficient tumour of the same cancer (blue; bar = 95% of draws) and, in-sample for the
rule, MSI in the three training cohorts (hollow). (b) Held-out exon shifts by class (n in brackets); box =
interquartile range, diamond = median, whiskers = 1.5 × IQR, dots = exons (not drawn for the large B_far set); the
axis is clipped to −16 to 21. (c) Exploratory: exons with a run ≥ 20 nt, share more included in MMR-deficient
tumours by branch-point distance (n in brackets); bars = Wilson 95% interval. (d) R_up arm of held-out tumours by
DNA class; dots = tumours (not drawn for the proficient group). (e) R_up arm against homopolymer deletions in
MMR-deficient tumours (blue, held-out) and MSI tumours of the three cohorts (grey).

- RV1, specificity: held-out +3.49 points, above 55 of the 57 other contrasts (largest 4.2, sarcoma LMS against
  MFS); q = 0.083, not resolved. The three in-sample MSI contrasts lie above all or all but one.
- RV2, a tumour-wide inclusion axis: the MMR coefficient on the R_up arm is +0.236 adjusted for it, stromal, immune
  and retained-intron levels (cancer only +0.237); q = 5.9e-8, supported.
- RV3, the tract is needed: R_up above B_far, AUROC 0.69; q = 1.9e-13, supported.
- RV4, AG exclusion zone in place of the branch point, runs ≥ 20: adjusted ρ = −0.041; q = 0.51, not resolved.
- RV5, R_up above shield exons: AUROC 0.56; q = 0.34, not resolved. Shield exons move up in 26 of 42 (exploratory).
- RV6, dose (panel e): partial ρ = +0.077 given purity; q = 0.36, not resolved.
- RV7, mismatch repair, not mutation burden (panel d): AUROC 0.74 against hypermutated proficient tumours; q =
  1.2e-7, supported.
- RV8, PTB relief: Spearman of the shift with PTB sensitivity +0.200, opposite to the sign PTB relief predicts; q =
  0.003, against.

Panel (c), exploratory: with the branch point 18–29 nt from the AG, 0.20 of long-run exons move up (n = 25); from
30–35 nt 0.73, the shape seen in the training cohorts (Figure U19a). q over family RV (8 tests). Per-tumour values
are kept locally (D-12).
