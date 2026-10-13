**Figure U23. Logistic regression on repeat exons, fitted in three cancers, finds MMR deficiency in 30 others.**
Pre-specified (config/followup.yaml: msi_code_models, D-71) before any model was fitted; A-19 restricts this figure
to one of the seven pre-registered models, chosen after the results. Label: ≥ 10 MC3 homopolymer deletions. The
model is L2 logistic regression (balanced class weights) on the 1,562 exons with an acceptor T-run ≥ 10 nt that are
measured in ≥ 25 cancers; each input is the exon's inclusion percentile within the tumour's cancer minus the
tumour's mean percentile of matched exons without a run. It was trained on all COAD, STAD and UCEC tumours (779; 150
MMR-deficient), with C chosen by leaving one training cancer out (C = 0.1), and tested on the 30 other cancers
(7,685 tumours, 50 MMR-deficient, prevalence 0.0065). Score G is the unfitted code.

(a) ROC held out; grey line = chance. (b) Precision against recall held out; dashed = prevalence. (c) AUROC per
held-out cancer with ≥ 3 MMR-deficient tumours (in brackets). (d) The 20 largest standardised coefficients, refitted
on the training cancers (gene, acceptor T-run length, branch-point distance; positive = more inclusion raises the
odds of MMR deficiency).

- Held-out AUROC 0.900 (95% CI 0.864–0.936) against 0.832 for score G; AUPRC 0.439 (0.348–0.543) against 0.083;
  precision 0.30 at recall 0.5. In the leave-one-training-cancer-out loop, 0.960.
- Leaving each cancer out of a model trained on the other 32 (exploratory): 0.914 on the same 30 cancers, 0.967 over
  all 33 (200 MMR-deficient).
- With MANTIS MSI-H against MSS as the held-out label (17 against 7,321; exploratory): AUROC 0.9997, AUPRC 0.897.
- LR4, the pre-registered test of this model (AUROC against the 10-bin logistic model, predicted lower): +0.032 (95%
  CI −0.019 to +0.089), q = 0.39, not resolved.

Coefficients by class (A-18, exploratory): median −0.0033 for runs of 10–19 nt (46% positive), +0.0144 for runs ≥ 20
nt with the branch point ≥ 36 nt (66% positive). The same model on the two arms' exons only reached 0.888
(0.848–0.927); on the 107 other repeat exons, 0.776 (0.713–0.839). The other pre-registered models, held-out AUROC:
two-arm logistic 0.835, 10-bin logistic 0.868, 10-bin MLP 0.873, all-exon logistic 0.902, all-exon MLP 0.885 (tables
in code/models/; family LR, 7 tests). Per-tumour predictions are kept locally (D-12).
