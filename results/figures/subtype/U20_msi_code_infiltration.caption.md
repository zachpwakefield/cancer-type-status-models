**Figure U20. The code is not an infiltration readout.** Pre-specified (config/followup.yaml: msi_code_infiltration,
D-68). Covariates as in D-55: stromal and immune = mean log2 CPM of eight stromal and eight immune genes
(config/biology.yaml); RI = the tumour's mean retained-intron PSI, an RNA-handling proxy; purity = ABSOLUTE, where
called. In the 30 held-out cancers (D-66) each covariate and score G are percentiles within the tumour's cancer.

(a) AUROC for MMR-deficient tumours (MSI in the three cohorts) scoring lower on G: all held-out tumours, then the
lower and upper halves of each cancer's immune score; grey = the immune score itself, read as higher in
MMR-deficient tumours; then MSI against the rest in the lower-immune half of UCEC, STAD and COAD (MMR-deficient or
MSI tumours in brackets). (b) OLS coefficient of MMR deficiency (or MSI) on score G after adjustment, as a share of
the coefficient with cancer (or site) fixed effects only; bar = 95% HC3 interval; p of the adjusted coefficient.

- IF1, IF2, held-out cancers: adjusted coefficient −0.339 and −0.334 (unadjusted −0.335); q = 4.0e-26, supported.
- IF3, lower-immune half: AUROC 0.84 (21 MMR-deficient), supported; upper half 0.83 (exploratory).
- IF4–IF6, the three cohorts with site fixed effects: UCEC −0.043, STAD −0.079, COAD −0.080; all supported, each
  keeping 1.01–1.07 of the site-only effect.

The immune score marks MMR-deficient tumours only weakly here (AUROC 0.59), and in the adjusted model more immune
expression goes with a higher G, the opposite of the MMR effect. q over family IF (6 tests). Per-tumour values are
kept locally (D-12).
