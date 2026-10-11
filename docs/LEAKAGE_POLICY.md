# Leakage policy

Binding rules for every model in this repository, and how each one is enforced.

| # | Rule | Enforcement |
|---|---|---|
| R1 | A held-out sample never influences anything fitted for the model that predicts it: feature filters, rankings, imputation values, scaling, hyper-parameters, epochs or weights. | `ctsm/pipeline.py`: the outer-test rows are touched only after refit, for prediction and SHAP. `tests/test_core.py::TestPreprocessorUsesTrainingRowsOnly` |
| R2 | All samples of one patient sit on the same side of every split (outer and inner). A matched tumour and normal are never split. | `StratifiedGroupKFold` with group = case; re-checked by `splits.check_splits` and V2 |
| R3 | Hyper-parameters are chosen on inner folds of the outer-train set only. | Optuna runs on `fold_indices(spl, k, j)`, rows inside outer-train k |
| R4 | No test-set result feeds back into any choice. All rules in `config/study.yaml` were fixed before the first outer-test evaluation. | `docs/DECISIONS.md` |
| R5 | Data-level operations are per-sample only: CPM normalisation and AFE/ALE structural zeros. No cross-sample statistic is computed before splitting, except the pre-filter below, which provably cannot affect selection. | `ctsm/extract.py` |
| R6 | Labels come only from GDC sample metadata (sample type, project), never from molecular values. | `ctsm/cohort.py` |
| R7 | Interpretation uses held-out samples only. | SHAP is computed on outer-test samples of each fold's model |

## Why the extraction pre-filter cannot leak (R5)

**Setup.** Fix a task cohort C (one cancer for task 1, all samples for task 2).

- The fold-level candidate rule (D-07) keeps a feature on a training set S ⊆ C iff:
  - `support_S ≥ f · |S|`, where support is the observed count (AS) or the CPM ≥ 1 count
    (GEX), and f is 0.5 (AS) or 0.1 (GEX);
  - and `var_S > 0`.
- Let `m_C` be the smallest training set the rule is ever applied to: the smallest inner-train
  set over all outer/inner folds of C.
- At extraction, a feature is kept if, for at least one cohort C:
  - `support_C ≥ ⌈f · m_C⌉`;
  - and the feature is non-constant over C.

**Claim.** Any feature that the fold rule keeps on some training set S ⊆ C also passes the
pre-filter.

**Proof.** Suppose the fold rule keeps the feature on S.
- Support counts are monotone in the sample set, so `support_C ≥ support_S ≥ f · |S| ≥ f · m_C`.
  Since `support_C` is an integer, `support_C ≥ ⌈f · m_C⌉`.
- If the feature has nonzero variance on S, it is non-constant on S ⊆ C, and hence on C.

So the feature passes the pre-filter for C. ∎

**Consequence.** Dropping features that fail the pre-filter changes no fold's candidate set,
and so changes no ranking, selection or model. The thresholds are recomputed from the splits
and checked by V5 (`scripts/03_validate_data.py`). The property itself is tested on random
data in `tests/test_core.py::TestPrefilterSuperset`.

## Known non-leakage caveats

- **Tissue-source-site and batch effects.** TCGA normals and tumours of the same patient share
  a tissue source site and often a plate. Case-grouped splitting does not remove site-level
  batch effects shared by different patients. A tumour-vs-normal or cancer-type model may
  partly use them, as in all TCGA classifiers. SHAP tables make it possible to check that top
  features are biologically plausible.
- **Read length.** Read length differs by project (about 48–76 bp), which changes how many AS
  events are quantifiable. In the pan-cancer task, missingness patterns therefore carry some
  project information. Median imputation after ranking prevents imputed values from driving
  selection. But a model can still learn that "value = training median" is common in some
  projects.
- **GEX normalisation.** It is CPM (depth only), not TPM. Gene-length effects cancel within a
  gene across samples, so cross-sample comparisons of a gene are unaffected.
