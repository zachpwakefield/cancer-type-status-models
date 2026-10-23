**Figure RP0. The rising arm passes 9 of the 12 pre-registered tests; the claim is partial.** Pre-registered as D-77
(config/followup.yaml: msi_rising_proof; commit 5471b27) before any junction-read value was compared by MMR status.
Confirmatory tests use the 30 held-out cancers (50 MMR-deficient tumours in 16 of them) and junction reads only
(rMATS JC), in PSI units. delta = an exon's mean adjusted PSI in MMR-deficient (MMRd) tumours minus the proficient
(MMRp) mean of the same cancer. One row per test shown, in the order of the argument; the id is bold for the 12
tests the verdict needs (family RP). RPD3b is the gap test of RP3e repeated in the training cancers. Not shown:
RP2b, RPC2b, RPD2a, RPD2b, RPD3a, RPM2 (all in tests.csv). Effects in different units sit in different columns: PSI
differences (class contrasts, regression coefficients, arm contrasts), Spearman ρ (rank associations and dose), and
AUROC or share differences (RP1: AUROC minus the median of matched exon sets; RPD1a: share of positive SHAP
directions minus 0.5; RP2c: share of R_up exons rising minus that of matched exons; RP6a: AUROC minus 0.5). Filled =
supported (predicted sign and q < 0.05 within the family), hollow = not resolved; bars = 95% intervals where the
test has one (RP2a and RPC2a: studentized, Welch–Satterthwaite). Grey diamonds = the same test on the falling arm
(runs of 10–19 nt), predicted below 0. p = raw p of the test (the larger of the exon and the label permutation p
where both exist); q and families in RP0_ladder.csv. The line under the panel gives the claim's verdict by the
pre-registered rule.

- Family RP (verdict): 9 of 12 supported; RPS (secondary held out): 5 of 9; RPP (not held out): 1 of 2; RPC (falling
  arm): 7 of 8; RPD (training, discovery): 5 of 6; RPM (mechanism): 0 of 3.
- RP2a: +0.0070 (95% interval 0.0033 to 0.0107), p = 1.0e-3; the whole interval lies above the SESOI of 0.003.
  RPC2a: −0.0030 (−0.0047 to −0.0013), p = 0.002.
- Not resolved in family RP: RP4a (poly-A) +0.0042, p = 0.25; RP4d (donor-side T-run) +0.0047, p = 0.088; RP6c
  +0.0087, p = 0.086.
- Claim verdict: partial; not met: RP4a, RP4d and RP6c not resolved (needed: supported); no refuting clause holds.
