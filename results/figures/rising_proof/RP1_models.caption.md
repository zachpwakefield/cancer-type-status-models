**Figure RP1. Held out, the rising arm separates MMR-deficient tumours better than matched exon sets, but removing
it from the fitted model costs no detectable accuracy.** Rising arm = acceptor T-run ≥ 20 nt with the
sequence-predicted branch point ≥ 36 nt from the AG; falling arm = acceptor T-run of 10–19 nt. MMRd = MMR-deficient
(≥ 10 MC3 homopolymer deletions); MMRp = MMR-proficient. (a) Held-out cancers (7,685 tumours, 50 MMRd): AUROC for
MMRd against MMRp tumours of each arm's level (a tumour's mean PSI deviation from its cancer's median over the arm's
exons, as a percentile within the cancer); dot = the arm; grey bar = middle 95% of 1,000 sets of background exons
matched to the arm on baseline band and coverage tertile within each cancer, which carry any tumour-wide shift as
well; tick = their median; dotted line = 0.5. The falling arm is predicted below 0.5. (b) Training cancers (COAD,
STAD, UCEC): the subtype networks (SE view, D-38) explained with SHAP. For each rising- or falling-arm input and
cohort, the direction is the Spearman correlation of the input's value with its SHAP over a fold's MSI test tumours,
averaged over the folds; bar = share of input × cohort units whose higher inclusion pushes towards MSI; grey band
behind each bar = the range expected by chance, the middle 95% (2.5th–97.5th percentiles) of that share over 1,000
permutations of the inputs' rows among each fold's MSI test tumours (A-22); dotted line = 0.5. (c) M4, step 54's
logistic regression on all repeat exons (fitted in the training cancers; held-out AUROC 0.900), refitted at the same
C without its 119 rising-arm inputs (line) and without 119 inputs drawn at random from its other repeat inputs, 200
times (grey; exploratory, A-25).

- RP1, rising-arm AUROC minus the matched median: +0.213 (AUROC 0.719 against 0.506), p = 1.0e-3, supported; RPC1,
  falling arm: −0.155 (0.337 against 0.492), p = 1.0e-3, supported.
- RPD1a, share of rising-arm directions minus 0.5: +0.378 (36 of 41), p = 1.0e-3, supported; RPD1b, falling arm:
  −0.228 (108 of 397), p = 1.0e-3, supported.
- M4 without its rising-arm inputs: AUROC change +0.0001 (paired bootstrap over test tumours with the fits held
  fixed, −0.0145 to 0.0146); random removals: median −0.0031, 2.5th–97.5th percentiles −0.0153 to 0.0071; the
  rising-arm removal lies at the 68.5th percentile. No test; no loss was detected (test-set interval −0.015 to
  +0.015). M4's standardised coefficients are positive for 66% of its rising-arm inputs and 46% of its falling-arm
  inputs.
