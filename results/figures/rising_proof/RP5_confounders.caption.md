**Figure RP5. Exon and tumour covariates leave the rise unchanged.** Held-out cancers, junction reads. Exon level:
OLS of delta on an R_up indicator against the B sample of background exons (92 R_up + 4,756 exons), with exon
covariates; the core splice-site score (−3..+3) is blind to the tract, while the full −20..+3 score and the
acceptor-Alu flag measure the tract itself and enter only a sensitivity model. Tumour level: OLS of each held-out
tumour's R_up arm on MMR status with cancer indicators and the tumour's global inclusion axis, stroma, immune
infiltrate, retained-intron level and purity. Blue = the pre-registered models; dark grey = no covariates (the
difference of the two class means at exon level; RP2b at tumour level); light grey = sensitivity models. No
intervals were computed for these coefficients; p is shown instead.

- RP5a, R_up coefficient, Freedman–Lane within band × tertile and MMR-label permutation, p = max: +0.0084, p =
  0.004, supported. RPC5a, R_down: −0.0022, p = 0.034, supported.
- RP5b, MMRd coefficient, Freedman–Lane with tumour residuals permuted within cancer: +0.0082, p = 1.0e-4,
  supported. RPC5b, R_down arm: −0.0033, p = 2.0e-4, supported.
- No covariates: R_up mean minus B-sample mean +0.0070; RP2b, arm contrast: +0.0084, p = 5.0e-4, supported.
