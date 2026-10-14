**Figure U26. Expression models trained the same way transfer almost as well; splicing adds to them.** Pre-specified
(config/followup.yaml: msi_vs_expression, D-74) before any expression model was scored on a held-out cancer. Every
model was trained on COAD, STAD and UCEC with the same label (MC3 ≥ 10 homopolymer deletions) and tuned by leaving
one training cancer out; each was tested on the 30 other cancers (7,685 tumours, 50 MMR-deficient). All are L2
logistic regressions with balanced class weights except the PreMSIm recipe (5 nearest neighbours on its 15 genes
scaled to [0, 1]) and MLH1 alone (unfitted). The two whole-transcriptome models keep the top K of 23,640 genes by
ANOVA F on the training tumours (chosen for both: K = 500, C = 0.1).

(a) Held-out AUROC and (b) AUPRC; bars = 95% bootstrap intervals over tumours within cancer and label; dashed = the
splicing model; grey = chance and prevalence. (c) AUROC of each model for the held-out MMR-deficient tumours whose
MLH1 is low (filled; below the 10th percentile of their cancer) or not (hollow), each against all held-out
proficient tumours. (d) The held-out MMR-deficient tumours: rank of each within its cancer by the splicing score
against the rank by the expression model on percentile inputs (fraction of the cancer at or above it; log axes, best
at the top right); grey line = equality. (e) Leaving each cancer out of training (E2), AUROC within cancer against
the tumours on which the two named exome measures agree (as Figure U24), for the splicing model and the expression
model on percentile inputs (refitted with its E1 parameters, without inner tuning).

- EX1, splicing against expression on percentile inputs: +0.034 (−0.026 to +0.099), q = 0.44, not resolved.
- EX2, against expression on log CPM: +0.008 (−0.042 to +0.057), q = 0.93, not resolved.
- EX3, against PreMSIm's genes refitted: +0.115 (+0.047 to +0.186), q = 0.001, supported.
- EX4, against MLH1 alone: +0.142 (+0.066 to +0.229), q = 0.001, supported.
- EX5, MMR-deficient tumours with MLH1 not low: +0.005 (−0.084 to +0.093), q = 0.95, not resolved.
- EX6, splicing + expression against expression: +0.060 (+0.023 to +0.100), q = 0.001, supported.

3 of 6 supported, none against (family EX). Held-out AUROC: splicing 0.900 (0.864–0.936); expression, log CPM 0.893
(0.849–0.930); percentiles 0.866 (0.809–0.916); PreMSIm genes 0.785 (0.722–0.845); PreMSIm recipe 0.702
(0.634–0.769); MLH1 0.758 (0.694–0.818); splicing + expression 0.926 (0.898–0.951). AUPRC 0.44, 0.40, 0.43 and 0.53
for splicing, the two expression models and the combination; precision at recall 0.5: 0.30, 0.26, 0.39 and 0.68.
Against MANTIS MSI-H (17 held out): splicing 1.000, expression 0.969 (log CPM) and 0.924 (percentiles). Panel (d),
Spearman ρ = 0.36. Panel (e): splicing 0.998, expression 0.975–0.979; the exome measures score 0.999–1.000 against
each other (Figure U24). E2 pooled AUROC, all 33 cancers: splicing 0.967, expression 0.941. Per-tumour points are
kept locally (D-12).
