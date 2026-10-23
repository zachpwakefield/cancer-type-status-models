**Figure RP0. In junction reads the rising arm rises in held-out cancers, but the claim is only partly supported.**
Pre-registered as D-77 (config/followup.yaml: msi_rising_proof; commit 5471b27) before any junction-read value was
compared by MMR status. The claim: in MMR-deficient tumours, exons after an acceptor T-run of ≥ 20 nt with the
sequence-predicted branch point ≥ 36 nt from the AG (the rising arm) are more included than in MMR-proficient
tumours of the same cancer, and the rise needs both features. Confirmatory tests use the 30 held-out cancers (50
MMR-deficient tumours in 16 of them) and junction reads only (rMATS JC), in PSI units. The eight steps are read in
order; each box gives the step's question, a plain answer and a marker: filled = every test of the step supported
(predicted sign and Benjamini–Hochberg q < 0.05 within its family), half = some supported, hollow = none resolved.
The mechanism box (dashed) is outside the claim's verdict. The verdict box gives the claim's verdict by the
pre-registered rule, under which "established" needs all 12 tests of family RP. Every test, with its effect, p and
q, is in section 5 of results/MSI_RISING_ARM_PROOF.md and in tests.csv.

- Step 1 (Do models use the rising arm?): RP1 supported; RPD1a supported.
- Step 2 (Does it rise in junction reads?): RP2a supported; RP2b supported; RP2c supported.
- Step 3 (Does it need both features?): RP3a supported; RP3b not resolved; RP3c supported; RP3d not resolved; RP3e
  supported.
- Step 4 (Is it the microsatellite?): RP4a not resolved; RP4b supported; RP4c supported; RP4d not resolved; RP4e
  supported.
- Step 5 (Does it survive confounders?): RP5a supported; RP5b supported.
- Step 6 (Is it mismatch repair?): RP6a supported; RP6b supported; RP6c not resolved.
- Step 7 (Is it there within the same patient?): RP7 supported.
- Step 8 (Does it follow dose?): RP8a not resolved; RP8b not resolved; RP8c not resolved.
- Step M (Mechanism: branch-point geometry or hnRNP C?): RPM1 not resolved; RPM2 not resolved; RPM3 not resolved.
- Families: RP (verdict) 9 of 12 supported; RPS (secondary, held out) 5 of 9; RPP (not held out) 1 of 2; RPC
  (falling arm) 7 of 8; RPD (training) 5 of 6; RPM (mechanism) 0 of 3.
- Claim verdict: partial; not met: RP4a, RP4d and RP6c not resolved (needed: supported); no refuting clause holds.
