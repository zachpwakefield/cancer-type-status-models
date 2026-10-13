**Is the code about the branch point, and about mismatch repair?**

Pre-specified (config/followup.yaml: msi_code_specificity, D-65; msi_code_endometrium, D-67). Exon shift = D-52's
background-adjusted difference in inclusion percentile. Adjusted ρ = D-63's: shift and branch-point distance
rank-residualised on run-length bins and mean inclusion.

a, exploratory: share of repeat exons more included in MSI (mean over ≥ 2 of the three cancers), by 10-nt bins of
branch-point distance, bins with ≥ 10 exons; grey line = the rule's 36-nt cut for runs ≥ 20. Exons whose branch
point lies within 30 nt of the AG move up rarely (0.12 of runs ≥ 20 at 20–29 nt, n = 24; 0.28 of runs 10–19); from
30–39 nt 0.73 and 0.51. b, adjusted ρ for MSI (blue) and the matching control (grey): repeat exons (R, acceptor
T-run ≥ 10) and exons without a tract (B, ≤ 5) measured in ≥ 2 cancers; UCEC MSI and POLE, each against the rest;
STAD MSI against the rest and CIN against GS. p = the test of the difference. c, AUROC of score A and score G for
MSI against the rest.

- BS1, the distance from the run's far end to the branch point adds beyond the run: ρ = 0.116, q = 0.001, supported.
- BS2, exons without a tract (acceptor T-run ≤ 5) show no such relation: R minus B +0.156, q = 0.001, supported.
- BS3, POLE (hypermutated, mismatch-repair proficient) against MSI in UCEC: difference +0.075, q = 0.012, supported.
  BS4, CIN against GS in STAD: +0.249, q = 0.001, supported.
- EG1, score G separates UCEC MSI: AUROC 0.84, q = 2.8e-9; EG2, ahead of score A (+0.33), q = 5.0e-4. Both
  supported.

The rule behind score G was read from these three cohorts (D-53 and U17b), so c is in-sample for the rule, and its
STAD and COAD values are exploratory; U18 is the out-of-sample test. q over families BS (4 tests) and EG (2 tests).
