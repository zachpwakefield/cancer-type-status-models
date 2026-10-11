# Decisions

Every rule that shapes the data or the models, with its rationale and where it is implemented.
Values live in `config/study.yaml`; this file explains them. All decisions were fixed before any
model was evaluated, and none was changed after seeing outer-test results.

### D-01 Samples · `cohort` · `ctsm/cohort.py`
- **Unit:** one RNA-seq file per (case, role). Role is *tumour* (GDC `Primary Tumor`, or
  `Primary Blood Derived Cancer - Peripheral Blood` for LAML, which has no solid primary) or
  *normal* (`Solid Tissue Normal`, barcode type 11).
- **Tie-break:** 46 cases have more than one primary-tumour file; we keep the earliest vial,
  then the lowest barcode, then the lowest UUID. Selection never looks at molecular values.
- **Not used:** metastatic (243 files, mostly SKCM), recurrent (47) and "additional new
  primary" (10) samples. They are a different biological question (and SKCM metastases carry
  the tissue of the metastatic site). They remain available for an external test later.
- **Required:** the file must be in all nine molecular tables. All 10,101 clinical files are;
  358 AS-table files have no clinical row (no sample type or project) and cannot be labelled.
- **Result:** 9,729 samples from 9,075 cases: 9,034 tumours (33 projects) and 695 normals.
  Summary: `results/data/sample_accounting.csv`.

### D-02 Task 1, tumour status · `tasks.tumor_status` · `ctsm/cohort.py`
- **Question:** within one cancer type, is a sample tumour (1) or adjacent normal tissue (0)?
- **Scope rule:** a TCGA project is modelled iff it has ≥ 10 normal cases. That puts at least
  2 normals in every outer test fold. The rule selects 15 projects: BLCA, BRCA, COAD, ESCA,
  HNSC, KICH, KIRC, KIRP, LIHC, LUAD, LUSC, PRAD, STAD, THCA and UCEC. It excludes CHOL and
  READ (9 normals each). `results/data/task_scope.csv`.
- **Caveat:** ESCA (11), UCEC (13) and BLCA (18) have few normals; their CIs are wide.
- **Pairing:** 654 of the 695 normals have a matched tumour from the same patient. Splits are by
  case (D-04), so a patient's tumour and normal are never on opposite sides of a split.

### D-03 Task 2, cancer type · `tasks.cancer_type` · `ctsm/cohort.py`
- **Question:** which of the 33 TCGA projects is a tumour from, or is it normal tissue?
- **Classes:** 34 = 33 projects (primary tumours) + one pooled `Normal` class holding every
  solid-tissue normal, whatever its tissue. The smallest class is CHOL (35 cases); every class
  must have ≥ 30 cases (asserted).
- **Known hard pairs:** COAD/READ, LUAD/LUSC, KIRC/KIRP/KICH and GBM/LGG share tissue of origin.
  Confusions among them are expected and reported.

### D-04 Nested, case-grouped, stratified CV · `splits` · `ctsm/splits.py`
- **Outer:** `StratifiedGroupKFold(5)`, with groups = case and strata = label. Every sample is
  predicted exactly once, by a model that never saw its case. Performance is pooled over these
  out-of-fold predictions.
- **Inner:** `StratifiedGroupKFold(3)` inside each outer-train set, used only for
  hyper-parameter tuning.
- **Inner folds used:** all three for task 1 (small cohorts, 88–1,185 samples). Task 2 uses
  only inner fold 0: one ⅔ : ⅓ train/validation split of about 7,800 outer-train samples, so
  validation has about 2,600 samples. This keeps 9 × 5 pan-cancer tunings tractable on a shared
  machine. It is a split of training data only, so it costs efficiency, not validity.
- **Seeds:** derived from `splits.seed` and the (task, cohort) name, so splits are identical
  for every modality of a cohort. Modalities are therefore compared on the same samples.

### D-05 Modalities, modelled independently · `views`
- Nine views, each with its own model, preprocessing and tuning; there is no fusion:
  - AS (alternative RNA processing): SE, MXE, A3SS, A5SS and RI (rMATS PSI); AFE and ALE
    (HITindex first/last-exon usage); HIT (the HITindex HIT index, −1…1);
  - GEX (gene expression).
- **GEX source:** `gex.csv` holds log2(raw STAR count + 1) for **every** file, normals included.
  We recover the integer counts (exactness is asserted) and store log2(CPM + 1), where library
  size is the total count over all 57,656 genes of the sample.
  - **Why not TPM:** the official GDC TPM matrix used by the sister survival project contains
    only primary tumours. Mixing TPM tumours with CPM normals would give the model a trivial
    batch signal. CPM uses no information from other samples.
  - **Check:** V4 in `results/data/VALIDATION.md` (median per-sample Spearman vs TPM 0.96).

### D-06 Feature scope and AS handling · `features` · `ctsm/annotation.py`, `ctsm/extract.py`
- **Chromosomes:** chr1–22 and chrX. chrY would re-encode sex, which identifies the sex-specific
  cancers trivially.
- **Host-gene biotype:** `protein_coding` or `lncRNA` (GENCODE v45; GDC v36 fallback).
- **AFE/ALE structural zeros:** HITindex reports only terminal exons with reads. Within one
  sample and gene, a blank sibling becomes 0 iff ≥ 1 observed sibling is nonzero, the observed
  siblings sum to 1 ± 1e-6, and ≥ 1 sibling is blank. The rule is per sample, so it uses no
  other sample's data.
- **AFE/ALE annotation:** after the zero rule, keep only metaexons annotated as first (AFE) or
  last (ALE) exons in `hg38_gen_meta`.
- **Lineage:** same rules as `cancer-survival-predictions` D-050, D-051 and D-053.
- **Storage:** rMATS PSI and HIT are exactly 3-decimal and stored as int16 × 1000 (asserted
  lossless); AFE, ALE and GEX are stored as float32. Missing stays missing: no imputation at
  the data level.

### D-07 Fold-level preprocessing and the extraction pre-filter · `preprocess` · `ctsm/preprocess.py`
Every step below is fitted on the training rows of the current fold only:
1. **Candidates.**
   - AS: value observed in ≥ 50 % of training samples.
   - GEX: CPM ≥ 1 in ≥ 10 % of training samples.
   - Both: nonzero training variance.
2. **Ranking:** either variance (unsupervised) or the one-way ANOVA F across classes, both on
   observed values. The method is a tuned hyper-parameter.
3. **Selection:** the top K features, K ∈ {500, 1000, 2000, 4000}, also tuned.
4. **Imputation:** training median, applied after ranking, so imputed values never influence
   which features are chosen.
5. **Scaling:** standardise with training mean/SD; clip to ±10.

**Extraction pre-filter (storage only).** Unfiltered matrices would be about 30 GB. Extraction
therefore keeps a feature only if it is non-constant and has support ≥ ⌈fraction × smallest
inner-train size⌉ in at least one task cohort. Because support counts and ranges can only shrink
on a subset, this keeps every feature that any fold-level rule could keep. It cannot change
which features a fold selects (proof in `LEAKAGE_POLICY.md`; tested in `tests/test_core.py`;
thresholds checked by V5).

### D-08 Model · `model` · `ctsm/model.py`
- **Architecture:** a multilayer perceptron:
  input dropout → [Linear → BatchNorm1d → ReLU → Dropout] × L → Linear.
  - widths taper by half per layer;
  - output: one logit (tumour log-odds) or 34 logits.
- **Why an MLP:** it is a deep model on tabular features of a few thousand dimensions, and
  every layer is handled exactly by `shap.DeepExplainer` (D-10).
- **Loss:** class-balanced cross-entropy, so each class contributes equally. The 0.5 threshold
  (binary) and the argmax (multiclass) therefore correspond to balanced decisions.
- **Optimiser and early stopping:** AdamW; early stopping on the validation balanced log loss
  (patience 20 for task 1, 12 for task 2; at most 300 / 100 epochs).
- **Final refit:** the whole outer-train set, for the median best epoch of the winning trial's
  inner folds.
- **Hardware:** CPU, single-threaded processes, seeded, so results are deterministic per job.

### D-09 Hyper-parameter tuning · `tuning` · `ctsm/pipeline.py`
- **Search:** Optuna TPE (seeded), 30 trials for task 1 and 20 for task 2, run per outer fold.
- **Objective:** mean inner-validation balanced log loss. It is a proper scoring rule that is
  insensitive to class imbalance and still discriminates between near-perfect classifiers,
  where AUROC saturates.
- **Pruning:** `MedianPruner` after each inner fold, once 6 trials are complete.
- **Search space:**
  - K and the ranking method (D-07);
  - layers 1–3; width 64–1024 (log scale);
  - dropout 0–0.6; input dropout 0–0.5;
  - learning rate 1e-4–3e-3 (log); weight decay 1e-6–1e-2 (log);
  - batch size {32, 64, 128, 256} for task 1 and {128, 256, 512} for task 2.
- **Records:** the chosen values per fold are in `folds.csv` of each job.

### D-10 SHAP · `shap` · `ctsm/explain.py`
- **Method:** `shap.DeepExplainer` on every outer-fold model, against 100 random outer-train
  samples as background.
- **Which samples:** every held-out (outer-test) sample is explained, so each sample is
  explained by a model that never saw it.
- **Explained output, task 1:** the tumour log-odds.
- **Explained output, task 2:** each sample's own true-class centred logit, zₖ − mean(z).
  - Softmax is unchanged by adding a constant to every logit; centring removes that free
    direction.
  - It is a fixed linear map, so DeepExplainer stays exact.
  - It needs one backward pass per sample instead of 34.
- **Additivity:** sum of SHAP values + expected value = output, asserted ≤ 1e-3 per fold.
- **Aggregation:** by feature ID across folds. A feature not selected by a fold contributes 0
  for that fold's samples, so features used by every fold rank highest. `n_folds_selected` is
  reported. Direction is the Spearman correlation of feature value vs SHAP.

### D-11 Metrics · `ctsm/metrics.py`
- **Task 1:** AUROC; AUPRC with tumour and with normal as positive; balanced accuracy;
  sensitivity; specificity; MCC; Brier score; log loss.
- **Task 2:** accuracy; balanced accuracy; macro and weighted F1; MCC; top-3 accuracy; macro
  one-vs-rest AUROC; log loss; per-class precision, recall, F1 and AUROC; confusion matrix.
- **Uncertainty:** pooled out-of-fold values with 95 % case-level bootstrap CIs (1,000
  resamples for task 1, 200 for task 2), plus the per-fold mean ± SD.

### D-12 What is published · `.gitignore`
- **Committed:** code, configuration, documentation, and aggregate results (metrics, fold
  tables, SHAP summaries and figures).
- **Local only (`data.nosync/`):** processed matrices, model weights, SHAP arrays and
  per-sample predictions.
- **Why:** they are individual-level data derived from controlled-access TCGA BAM files.
