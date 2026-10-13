**Figure U23. Ten sequence bins beat the two arms held out, and per-exon models did not transfer worse.**
Pre-specified (config/followup.yaml: msi_code_models, D-71) before any model was fitted. Label: ≥ 10 MC3 homopolymer
deletions. E1: trained on all COAD, STAD and UCEC tumours (779; 150 MMR-deficient), tested on the 30 other cancers
(7,685 tumours, 50 MMR-deficient, prevalence 0.0065). Inputs, all within the tumour's cancer: score G; the two arms
and ten sequence bins (run-length bins for runs of 10–19 nt; runs ≥ 20 nt by 10-nt branch-point bins), each against
matched exons without a run, as robust z; and each exon's background-adjusted percentile (repeat exons only, or all
30,632 exons measured in ≥ 25 cancers, top K by ANOVA). Hyperparameters chosen by leaving one training cancer out.

(a) AUROC held out (filled; bar = 95% bootstrap interval over tumours within cancer and label), in the inner
leave-one-training-cancer-out loop (hollow) and on the same 30 cancers when each is left out of a model trained on
the other 32 (cross; not run for the all-exon MLP); dashed = score G held out. (b) AUPRC held out; dashed =
prevalence. (c) Precision against recall held out for score G, the 10-bin logistic model and the all-exon logistic
model; dashed = prevalence. (d) Standardised coefficients of the 10-bin logistic model (blue = more inclusion of the
bin's exons raises the odds of MMR deficiency).

- LR1, two weighted arms against G: AUROC +0.003 (95% CI −0.026 to +0.029), q = 0.87, not resolved.
- LR2, ten bins against two arms: +0.034 (95% CI +0.006 to +0.062), q = 0.045, supported.
- LR3, MLP against logistic on the ten bins: +0.005 (95% CI −0.015 to +0.027), q = 0.75, not resolved.
- LR4, per-exon weights on repeat exons against the bins: +0.032 (95% CI −0.019 to +0.089), q = 0.39, not resolved.
- LR5, all exons against the bins: +0.034 (95% CI −0.010 to +0.080), q = 0.28, not resolved.
- LR6, all-exon MLP against the bins: +0.017 (95% CI −0.028 to +0.065), q = 0.65, not resolved.
- LR7, AUPRC of the bins against G: +0.173 (95% CI +0.108 to +0.237), q = 0.004, supported.

Held-out AUROC: M0 0.832, M1 0.835, M2 0.868, M3 0.873, M4 0.900, M5 0.902, M6 0.885. AUPRC: M0 0.083, M1 0.117, M2
0.257, M3 0.265, M4 0.439, M5 0.470, M6 0.492. Leaving one cancer out (exploratory), AUROC over all 33 cancers (200
MMR-deficient): M0 0.843, M1 0.852, M2 0.895, M3 0.851, M4 0.967, M5 0.957; over the 30: M0 0.832, M1 0.833, M2
0.870, M3 0.817, M4 0.914, M5 0.905. Exploratory (A-18, chosen after these results): the all-exon model without
repeat exons fell to AUROC 0.777 (0.720–0.837); per-exon weights on the two arms' exons alone gave AUROC 0.888
(0.848–0.927), and on the 107 other repeat exons AUROC 0.776 (0.713–0.839). q over family LR (7 tests). Per-tumour
predictions are kept locally (D-12).
