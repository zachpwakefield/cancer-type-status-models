# Methods

A manuscript-style description of the pipeline. Decision IDs (D-xx) point to
[DECISIONS.md](DECISIONS.md); leakage rules (R-x) point to [LEAKAGE_POLICY.md](LEAKAGE_POLICY.md).

## Data

The data are TCGA RNA-seq files with GDC sample metadata. Alternative RNA processing was
quantified per file:
- rMATS percent-spliced-in (PSI) for skipped exons (SE), mutually exclusive exons (MXE),
  alternative 3′ and 5′ splice sites (A3SS, A5SS) and retained introns (RI);
- HITindex usage of alternative first and last exons (AFE, ALE) and the HIT index (HIT), which
  scores whether a metaexon behaves as a first (−1) or last (+1) exon.

Gene expression (GEX) is the STAR gene count of the same files, converted to log2(CPM + 1)
(D-05). Raw tables are read in place, and their SHA-256 hashes are recorded
(`results/data/source_manifest.json`).

**Sample selection (D-01).** We keep one file per case and role. Tumour means primary tumour,
or peripheral blood for LAML; normal means solid tissue normal. The result is 9,729 samples
from 9,075 cases across 33 projects (9,034 tumours, 695 normals).

**Feature preparation (D-06).**
- Features are restricted to chr1–22 and chrX, in protein-coding or lncRNA host genes.
- Missing AFE/ALE usage of a gene's unexpressed terminal exons is set to 0 when that sample's
  observed siblings sum to 1. After that, only annotated first/last exons are kept.
- Values are stored losslessly. Missing values stay missing until the fold-level imputation.

## Tasks

1. **Tumour status (D-02).** A binary classifier, tumour vs adjacent normal, trained
   separately for each of the 15 projects with ≥ 10 normal cases.
2. **Cancer type (D-03).** A 34-class pan-cancer classifier: the 33 TCGA projects (primary
   tumours) plus one pooled Normal class.

Each task is fitted independently for each of the nine modalities (D-05). That makes 15 × 9 +
9 = 144 modelling jobs.

## Validation design (D-04)

We use nested, case-grouped, label-stratified cross-validation:
- **Outer folds:** 5. Every sample is predicted once, by a model whose preprocessing, tuning and
  training never saw the sample's case (R1, R2).
- **Inner folds:** 3, inside each outer-train set, for tuning (R3). The pan-cancer task uses
  inner fold 0 only.
- **Shared splits:** splits depend only on the (task, cohort), so all nine modalities are
  evaluated on identical folds.

## Model fitting, per outer fold (D-07 to D-09)

1. **Preprocessing.** Fitted on training rows only.
   - Candidates are AS features observed in ≥ 50 % of samples, or genes with CPM ≥ 1 in ≥ 10 %
     of samples, with nonzero variance.
   - Candidates are ranked by variance or by ANOVA F, and the top K are kept.
   - Missing values are imputed with the training median after ranking, then features are
     standardised.
2. **Model.** A multilayer perceptron:
   input dropout → [Linear → BatchNorm → ReLU → Dropout] × L → Linear.
   It is trained with class-balanced cross-entropy (AdamW) and early-stopped on the validation
   balanced log loss.
3. **Tuning.** Optuna TPE (30 trials for task 1, 20 for task 2) with median pruning. The search
   covers the ranking method, K ∈ {500, 1000, 2000, 4000}, depth 1–3, width 64–1024, dropout,
   input dropout, learning rate, weight decay and batch size. The objective is the mean
   inner-validation balanced log loss.
4. **Refit.** The best configuration is refitted on the whole outer-train set for the median
   best epoch of its inner folds. It then predicts the outer-test fold.

## Evaluation (D-11)

- Metrics are computed on pooled out-of-fold predictions, with 95 % case-level bootstrap CIs
  and per-fold mean ± SD.
- **Task 1:** AUROC, AUPRC (normal as positive), balanced accuracy, specificity, MCC.
- **Task 2:** accuracy, balanced accuracy, macro F1, macro one-vs-rest AUROC, top-3 accuracy,
  per-class metrics and the confusion matrix.

## Interpretability (D-10)

- **Method.** Every outer-fold model is explained with `shap.DeepExplainer`. This is DeepLIFT
  with Shapley-style averaging over a background of 100 random outer-train samples.
- **Samples.** Explanations cover all of that fold's held-out samples (R7).
- **Explained output.**
  - Task 1: the tumour log-odds.
  - Task 2: for each sample, the logit of its true class centred on the mean logit, zₖ − z̄.
    This measures evidence for that class relative to the others. It is invariant to the
    softmax shift and stays an exact linear target for DeepExplainer.
- **Additivity.** SHAP values sum to the explained output minus its background expectation;
  this is asserted for every sample.
- **Feature importance.**
  - Importance is the mean |SHAP| over all explained samples of all folds. A fold that did not
    select a feature contributes 0, so consistently selected features rank first.
  - Direction is the Spearman correlation between the standardised feature value and its SHAP
    value. A positive value means higher PSI, usage, HIT index or expression pushes towards
    tumour (task 1) or towards the sample's class (task 2).
  - Per-gene importance sums its events.

## Software

Python 3.10, numpy 1.26, pandas 2.3, scikit-learn 1.7, PyTorch 2.2 (CPU), Optuna 4.6 and
SHAP 0.48 (`environment.yml`). Each job records the git commit, package versions and
configuration in its `metrics.json`.
