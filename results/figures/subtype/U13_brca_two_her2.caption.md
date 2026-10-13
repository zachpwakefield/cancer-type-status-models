**Figure U13. Splicing reads HER2 through the HER2-enriched programme, not through the 17q12 amplicon.**
Pre-specified before any probability was compared with copy number (config/followup.yaml: brca_her2, D-58). P(Her2)
is the stored out-of-fold probability of the BRCA_PAM50 models: GEX, and the AS ensemble (mean of the 8 splicing
views). ERBB2 amplified = ABSOLUTE total copy number at ERBB2 ≥ 2 × ploidy (31 of 714 LumA/LumB and 30 of 54 PAM50
Her2 tumours). ABSOLUTE's copy number at ERBB2 does not exceed 7 in any tumour, so the rule catches the
highest-level amplifications only: 34% of clinically HER2-positive tumours meet it. Clinical HER2 (IHC or ISH) is
the pre-specified second reading.

(a) P(Her2) by PAM50 group and amplification; box = interquartile range, diamond = median, whiskers = 1.5 × IQR,
dots = tumours (n under each group). (b) AUROC of P(Her2) for each comparison, GEX (grey) and AS ensemble (blue).
(c) In-silico knockout (exploratory, A-12): in each stored fold model the inputs from chr17:39.6–39.8 Mb (GRCh38;
PPP1R1B to GRB7) are set to their training mean, and the LumA/LumB AUROC is recomputed (hollow); tick = 20 random
input sets of the same size.

- H1, GEX P(Her2) is higher in amplified luminal tumours: AUROC 0.72, q = 9.3e-5, supported.
- H2, splicing would see amplification less well than GEX: AUROC difference −0.09 (q = 0.016), **against**: the AS
  ensemble separates amplified luminal tumours better. Without the ABSOLUTE ceiling (D-60, U15) this lead is not
  resolved, and absent on the newly amplified tumours.
- H3, within PAM50 Her2, neither model separates amplified from non-amplified tumours (difference −0.003, not
  resolved).
- H4, AS P(Her2) is higher in non-amplified Her2 than in amplified luminal tumours: AUROC 0.87, supported; GEX 0.97.
- c, removing the amplicon inputs lowers the GEX AUROC from 0.716 to 0.686; in the splicing views the change is at
  most 0.007. The splicing models do not reach these tumours through the amplicon; they give them part of the
  HER2-enriched splicing programme.

q over family H (4 tests). Per-tumour probabilities are kept locally (D-12).
