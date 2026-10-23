**Figure RP5. The rise remains after exon and tumour covariates.** Held-out cancers, junction reads. MMRd =
MMR-deficient (≥ 10 MC3 homopolymer deletions); MMRp = MMR-proficient. (a) OLS of each exon's change in inclusion on
a rising-arm indicator against the B sample of background exons (92 rising-arm + 4,756 exons). The pre-registered
model adds eight exon covariates: baseline, length, coverage, host-gene change, poison status, the core splice-site
score (−3..+3, blind to the tract), GC and Alu exon. The full −20..+3 score and the acceptor-Alu flag measure the
tract itself and enter only a sensitivity model. (b) OLS of each held-out tumour's rising-arm level on MMR status
with cancer indicators. The pre-registered model adds the tumour's global inclusion axis, stroma, immune infiltrate,
retained-intron level and purity. Blue = the pre-registered models; dark grey = no covariates; light grey =
sensitivity models. No intervals were computed for these coefficients; their p values are below.

- RP5a, rising-arm coefficient, Freedman–Lane within band × tertile and MMR-label permutation, p = max: +0.0084, p =
  0.004, supported. RPC5a, falling arm: −0.0022, p = 0.034, supported.
- RP5b, MMRd coefficient, Freedman–Lane with tumour residuals permuted within cancer: +0.0082, p = 1.0e-4,
  supported. RPC5b, falling-arm level: −0.0033, p = 2.0e-4, supported.
- No covariates: rising-arm mean minus background mean +0.0070; RP2b, arm contrast: +0.0084, p = 5.0e-4, supported.
