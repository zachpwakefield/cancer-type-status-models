**Figure RP1. The rising arm separates MMR-deficient tumours better than matched exon sets, but the fitted model
does not need it.** (a) Training cancers (COAD, STAD, UCEC): the subtype MLPs (SE view, D-38) explain each tumour's
own class, so the MSI class's SHAP exists for the MSI test tumours. For each R_up or R_down input and cohort, the
direction is the Spearman correlation of the input's value with its SHAP over a fold's MSI test tumours, averaged
over the folds; bars = share of input × cohort units with a positive direction (higher inclusion pushes towards
MSI); grey band = 2.5th–97.5th percentiles of 1,000 permutations of the inputs' rows among each fold's MSI test
tumours (A-22); dashed line = 0.5. (b, c) Held-out cancers (7,685 tumours, 50 MMRd): AUROC for MMRd against MMRp of
the R_up (b) and R_down (c) arm (the tumour's mean PSI deviation from its cancer's median over the class's exons, as
a percentile within the cancer); line = the class; grey = 1,000 sets of background exons matched to the class on
baseline band and coverage tertile within each cancer, which carry any tumour-wide shift as well; dotted line =
AUROC 0.5. (d) M4, step 54's logistic regression on all repeat exons (fitted in the training cancers; held-out AUROC
0.900), refitted at the same C without its 119 R_up inputs (line) and without 119 inputs drawn at random from its
other repeat inputs, 200 times (grey; exploratory, A-25).

- RPD1a, share of R_up directions above 0.5: +0.378 (36 of 41), p = 1.0e-3, supported; RPD1b, R_down: −0.228 (108 of
  397), p = 1.0e-3, supported.
- RP1, R_up AUROC minus the matched median: +0.213 (AUROC 0.719 against 0.506), p = 1.0e-3, supported; RPC1, R_down:
  −0.155 (0.337 against 0.492), p = 1.0e-3, supported.
- M4 without its R_up inputs: AUROC change +0.0001 (paired bootstrap over test tumours with the fits held fixed,
  −0.0145 to 0.0146); random removals: median −0.0031, 2.5th–97.5th percentiles −0.0153 to 0.0071; the R_up removal
  lies at the 68.5th percentile. No test; the fitted model does not need the rising arm. M4's standardised
  coefficients are positive for 66% of its R_up inputs and 46% of its R_down inputs.
