# Subtype task: methods and decisions

The third task of the repository: **within one cancer, which subtype is a primary tumour?** It is
modelled, like the other two, independently in each of the nine views (AFE, ALE, HIT, SE, A5SS,
A3SS, MXE, RI, GEX). Decisions are numbered from D-30 (D-48..D-51, which span the three tasks, are in
[ACROSS_TASKS_METHODS.md](ACROSS_TASKS_METHODS.md), as is amendment A-8; D-52..D-56, the MSI-repeat
follow-up, and D-57..D-60, the BRCA follow-up, D-61..D-64, the MSI follow-up, and D-65..D-76, the MSI sequence code, are here; D-01..D-12 are in
[DECISIONS.md](DECISIONS.md), D-13..D-29 in [BIOLOGY_METHODS.md](BIOLOGY_METHODS.md)). Values live
in [config/subtype.yaml](../config/subtype.yaml); this file explains them.

**Order of work.** Every label map, scope rule, threshold and literature expectation below was
written into `config/subtype.yaml` before the first subtype model was fitted. The molecular and
IHC label set (D-43..D-45) was added on 2026-10-12, also before any model was fitted. Nothing
about model performance has changed any of them; the one amendment made after a model number
was printed is logged under [Amendments and deviations](#amendments-and-deviations) with its
reason.

| Step | Script | What it does |
|---|---|---|
| 30 | `scripts/30_subtype_labels.py` | labels, task tables, nested splits, label accounting |
| 31 | `scripts/31_subtype_extract.py` | matrices for the subtype cohorts, pre-filter recomputed per cohort |
| 32 | `scripts/32_subtype_run_job.py` | one (view, cohort) job: nested CV, tuning, refit, SHAP |
| 33 | `scripts/33_subtype_run_all.py` | every job, 6 single-threaded workers |
| 34 | `scripts/34_subtype_summarise.py` | performance tables, baselines, AS ensemble and fusion, McNemar |
| 35 | `scripts/35_subtype_interpret.py` | SHAP by class, AS vs GEX genes, literature markers, event audit, cross-family, shared-class and MSI-confound tables |
| 36 | `scripts/36_subtype_figures.py` | figures U1-U9 |
| 37 | `scripts/37_subtype_regulators.py` | splicing-factor scan, event association and predictions R1-R6 (D-46) |
| 38 | `scripts/38_subtype_site_probe.py` | patient-grouped vs site-held-out probe (D-47) |

Library: [`src/ctsm/subtype.py`](../src/ctsm/subtype.py). Tests: `tests/test_subtype.py`.

---

### D-30 The question and the label families · `label_sets` · `ctsm/subtype.py`
- **Question:** within one TCGA project, which subtype is a primary tumour? One model per
  (cohort, view); a cohort is one project and one label family.
- **Two label families**, kept apart throughout because they test different things:
  - *Histology* (`label_set: histology`): the pathologist's morphological type, assessed without
    the RNA. It is the only label family for the cancers with no published molecular subtype
    (SARC, TGCT, MESO, THYM, CESC).
  - *Molecular and IHC* (`label_set: molecular`, D-43..D-45): published subtype calls and
    immunohistochemistry/ISH status. Each carries a `label_type` saying how it was made
    (D-44), because a label derived from mRNA expression is a different test from one derived
    from DNA or protein.
- **Histology sources:** the class comes from the TCGA Clinical Data Resource (TCGA-CDR; Liu et
  al. 2018, *Cell* 173:400), Table S1, column `histological_type`, one row per patient. The GDC
  clinical table's `primary_diagnosis` (ICD-O-3 text of the patient's selected primary-tumour
  file) must map to the same class (D-33).
- **The CDR file** is registered in `config/sources.yaml` (`subtype_labels.tcga_cdr`), hashed by
  step 00, and read with the standard library (`read_xlsx_sheet`: an xlsx is a zip of XML; shared
  strings are resolved, `#N/A` cells kept as text) because no environment here has openpyxl.
  The parser reproduces an independently parsed copy of the sheet cell for cell (11,160 × 33,
  0 mismatches, checked once by hand on 2026-10-12).
- **Leakage rule R6 holds:** labels come only from clinical and published annotation tables,
  never from a molecular value of this repository's matrices.
- **Extensibility:** a further label family is a new `label_sets` entry plus a loader in
  `LABEL_SOURCES`; the modelling code does not change.

### D-31 Unit: one primary tumour per patient · `ctsm/subtype.py: label_cohort`
- The samples are the primary tumours of the D-01 cohort for the project: exactly one RNA-seq
  file per patient (D-01's tie-break), asserted. Normals, metastases, recurrences and additional
  primaries are never used.
- Every label is patient-level in its source and is attached to that patient's selected primary
  tumour. For the PanCanAtlas subtype table, rows of other sample types are dropped first (D-43).

### D-32 Label maps: exhaustive, with explicit exclusions · `classes`, `exclude`, `ignore`
- **Exhaustive.** For each cohort and source, every raw value met in the cohort is listed either
  under a class or as excluded. An unlisted value stops step 30 with an error, so a typo in the
  map cannot silently drop or relabel patients. `results/subtype/data/label_accounting.csv`
  lists every raw value with its patients, class and how many were labelled.
- **Excluded categories** are those that do not name one subtype: NOS, "Other, specify",
  "Mixed Histology", "Not Available", mixed serous-endometrioid, adenosquamous, untyped
  adenocarcinoma. They never appear in a task table (tested).
- **Multi-component values** (TGCT and THYM in the CDR, and joined columns, D-44) are split at
  `|`. Components listed in `ignore` (a test not done) are dropped; any `exclude` component
  excludes the patient; the remaining components must all map to one class, otherwise the
  patient is excluded as mixed.
- **Choices inside the maps**, each with its reason:
  - *SARC:* "Pleomorphic 'MFH'", "Undifferentiated Pleomorphic Sarcoma (UPS)" and the giant-cell
    MFH are one class, UPS: MFH was renamed UPS in the WHO 2002 classification. Desmoid tumour
    is excluded: it is an intermediate-grade fibromatosis, not a sarcoma.
  - *STAD (Lauren):* signet-ring cell carcinoma joins the diffuse type, as in Lauren's
    definition; mucinous, tubular and papillary adenocarcinoma join the intestinal type because
    TCGA's own form files them under "Intestinal Adenocarcinoma". "Adenocarcinoma, NOS" is not
    Lauren-typed and is excluded, never merged into intestinal.
  - *THCA:* GDC renders TCGA's tall-cell category as "Papillary carcinoma, columnar cell"; the
    35 concordant patients confirm the mapping is one-to-one.
  - *TGCT:* seminoma = every component seminoma; non-seminoma = every component non-seminomatous.
    Mixed tumours with a seminoma component are excluded (17): clinically they are managed as
    non-seminomas, but the sampled piece of an RNA-seq aliquot could be either. A component
    recorded as "[Not Available]" excludes the patient, because it could be seminoma.
  - *THYM:* the WHO types are grouped A/AB against B1-B3; a value within one group (e.g.
    B2|B3) keeps its group. Thymic carcinoma (type C) is its own candidate class (D-34).
  - *LGG:* oligoastrocytoma is excluded: it is not an entity in WHO 2016/2021.
  - *BRCA:* the special types (mucinous, metaplastic, medullary) are candidate classes; mixed
    ductal-lobular and other mixed types are excluded.
  - *MESO:* sarcomatoid is a candidate class (1 patient); "Diffuse malignant mesothelioma - NOS"
    is excluded.

### D-33 Histology: two sources must agree · `sources.concordance`
- A histology label is kept only if the GDC `primary_diagnosis` maps to the same class as the
  CDR `histological_type`. Disagreements and GDC values outside the map exclude the patient.
- **Agreement:** of the 3,433 patients with a class in both sources, 3,416 agree (99.5%); per
  cohort in `results/subtype/data/source_concordance.csv`. The rule removes 17 discordant
  patients and 31 whose GDC value is excluded.
- **Why:** the two fields are entered separately; agreement removes transcription errors and
  the mixed tumours that one field records as mixed and the other as pure.

### D-34 Minimum class size · `rules`
- **≥ 20 patients per class, ≥ 2 classes per cohort.** With 5 outer folds a class of 20 puts
  about 4 patients in every test fold and about 5 in every inner-validation fold — the least at
  which a class's inner-fold loss is not decided by one or two patients, and at which its pooled
  recall has a 95% interval narrower than about ±0.2.
- A class below the minimum is **dropped, not merged**: its patients leave the cohort. The only
  merges are the biological groupings fixed in D-32 and D-44 (THYM A/AB and B1-B3; LGG_IDH).
- **What the rule removed:** BRCA mucinous (17), metaplastic (9), medullary (5); SARC MPNST (9),
  synovial (10); MESO sarcomatoid (1); THYM thymic carcinoma (11); STAD_TCGA HM-SNV (7);
  COAD_TCGA HM-SNV (2); PRAD_ETS ETV4 (14), FOXA1 (9), FLI1 (4), IDH1 (3). The threshold was 20
  before these counts were tabulated against it; at 15, BRCA mucinous (17) would have entered.

### D-35 Histology cohorts in scope, and why the others are not · `cancers`
Twelve projects, in three tiers. Counts are patients after D-32..D-34
(`results/subtype/data/cohort_scope.csv`).

**Primary: histology that names a biologically distinct entity.**

| Cohort | Classes (patients) | Why it is a fair question |
|---|---|---|
| BRCA | ductal 759, lobular 193 | Invasive lobular carcinoma is defined by E-cadherin loss and discohesive growth; CDH1 is mutated in most ILC (Ciriello et al. 2015, *Cell* 163:506). |
| THCA | classical 351, follicular 100, tall cell 35 | The follicular variant is RAS-like; classical and tall-cell tumours are BRAF-like and less differentiated (TCGA 2014, *Cell* 159:676). Tall cell vs classical is expected to be hard: both are BRAF-like. |
| UCEC | endometrioid 104, serous 58 | Serous carcinomas fall almost entirely into the copy-number-high, TP53-mutant group (TCGA 2013, *Nature* 497:67). |
| SARC | LMS 102, DDLPS 57, UPS 47, MFS 24 | Distinct lineages; DDLPS carries 12q13-15 (MDM2/CDK4) amplification; UPS and myxofibrosarcoma are genomically similar and form a spectrum (TCGA 2017, *Cell* 171:950). UPS vs MFS is expected to be hard. |
| TGCT | seminoma 64, non-seminoma 50 | Developmentally distinct germ-cell tumours: seminoma retains a primordial-germ-cell programme (KIT, SOX17), embryonal carcinoma a pluripotent one (SOX2) (Shen et al. 2018, *Cell Rep* 23:3392). |
| MESO | epithelioid 55, biphasic 21 | Biphasic tumours carry a sarcomatoid component; TCGA describes an epithelioid-to-sarcomatoid gradient (Hmeljak et al. 2018, *Cancer Discov* 8:1548). |
| STAD | intestinal 157, diffuse 70 | Lauren types; diffuse tumours are enriched in the genomically stable group with CDH1 and RHOA alterations (TCGA 2014, *Nature* 513:202). |

**Weak label: kept, and flagged as weak wherever it is reported.**

| Cohort | Classes | Why it is weak |
|---|---|---|
| LGG | astrocytoma 191, oligodendroglioma 187 | Histological typing of diffuse glioma is poorly reproducible between pathologists, and since WHO 2016 the types are defined by IDH mutation and 1p/19q codeletion, not by morphology (Louis et al. 2016, *Acta Neuropathol* 131:803; TCGA, Brat et al. 2015, *NEJM* 372:2481). In this cohort only 128 of 186 histological oligodendrogliomas are IDH-mutant codeleted and 6 of 189 astrocytomas are (`label_crosstab.csv`). The molecular question is asked directly in LGG_IDH (D-44). |
| THYM | A/AB 51, B1-B3 57 | The WHO scheme is ordinal (A, AB, B1, B2, B3) and the cohort is small. A/AB against B1-B3 follows the medullary/cortical divide and the GTF2I mutation, which is frequent in A and AB and rare in B thymomas (Petrini et al. 2014, *Nat Genet* 46:844); but AB thymomas contain lymphocyte-rich B-like areas, so the groups overlap by construction. |

**Controls: near-trivial contrasts, to show the pipeline finds what it should.**

| Cohort | Classes | What it controls for |
|---|---|---|
| ESCA | squamous 80, adenocarcinoma 81 | Two epithelial lineages: squamous oesophageal cancer resembles head-and-neck squamous cancer, adenocarcinoma resembles CIN gastric cancer (TCGA 2017, *Nature* 541:169). |
| CESC | squamous 240, adenocarcinoma 45 | The same lineage contrast in the cervix (TCGA 2017, *Nature* 543:378). |
| PCPG | pheochromocytoma 148, paraganglioma 30 | A site label (adrenal vs extra-adrenal), not a morphological one. Adrenaline synthesis (PNMT) is adrenal-medullary, and adrenal tumours carry admixed adrenal cortex, which TCGA saw as a "cortical admixture" expression subtype (Fishbein et al. 2017, *Cancer Cell* 31:181). It tests whether a model reads the tumour or its surroundings. |

**Not in scope, and why.**
- *No second class of ≥ 20 in the histology field:* ACC, BLCA, CHOL, DLBC, GBM, HNSC, KICH,
  KIRC, KIRP, LIHC, LUSC, OV, PAAD (the second class is "Other subtype"), PRAD, READ, UCS
  (heterologous 19, homologous 13).
- *LUAD:* the histology field is a single class; the GDC growth patterns are mostly "NOS" or
  "mixed subtypes".
- *UVM:* the Callender cell type is a graded mixture (37 of 78 are mixed spindle-epithelioid),
  pure epithelioid is 13, and the molecular axis that matters (monosomy 3/BAP1) is not local.
- *COAD mucinous vs conventional:* the label is a threshold (≥ 50% extracellular mucin) with a
  residual comparator, and its known molecular correlate (MSI) is asked directly in COAD_TCGA.
- *LAML FAB morphology:* superseded by the WHO genetic classification; the three classes ≥ 20
  (M1, M2, M4) are adjacent maturation stages.
- *LGG grade (G2/G3):* an ordinal grade, poorly reproducible, and since WHO 2021 assigned within
  a molecular type.
- *SKCM:* 88 primary tumours and no histology field.

### D-36 Validation and model settings · `modelling` · `ctsm/subtype.py: subtype_config`
- The D-04..D-10 machinery is used unchanged, through `pipeline.run_fold`: 5 outer folds
  (`StratifiedGroupKFold`, groups = patient, strata = label), 3 inner folds inside each
  outer-train set, Optuna TPE with median pruning on the mean inner-validation balanced log
  loss, fold-level preprocessing (D-07), class-balanced MLP (D-08), refit on outer-train for the
  median best epoch, SHAP on every held-out patient.
- **Settings:** the tumour-status values, because the cohorts (76-1,073 patients) are the size
  of the tumour-status cohorts: 3 inner folds used, 30 trials, at most 300 epochs, patience 20,
  batch size {32, 64, 128, 256}.
- **Seeds** derive from `splits.seed` and ("subtype", cohort), so all nine views of a cohort see
  the same folds; the subtype splits are independent of the tumour-status splits.
- **Folds:** `StratifiedGroupKFold` with shuffling is greedy, so the smallest class can be thin
  in one test fold (2 patients in the worst case: MESO, UCEC_TCGA, COAD_TCGA;
  `split_summary.csv`). All metrics are pooled over the out-of-fold predictions, which every
  patient contributes exactly once.
- **Configuration isolation:** `config/study.yaml` is not modified. `subtype_config()` deep-copies
  the study block and adds the `subtype` entries, and writes under `<data_dir>/subtype/`, so the
  two existing tasks run with exactly the configuration and the files they were run with.
- **Code change to the shared pipeline, behaviour-preserving:** `run_fold` decided "one-logit
  model" by `task == "tumor_status"`; it now decides by `C == 2` (identical for the two existing
  tasks, which always have 2 and 34 classes) and accepts an optional `fold_metrics` hook for the
  per-fold test metrics. A re-run of an existing tumour-status fold with the changed code
  reproduces its stored predictions and SHAP values bit for bit (end of this file).
- **Run order,** fixed in advance: histology cohorts, then the genomic/IHC cohorts, then the
  expression-defined ones; within a group, the largest first. 6 workers at most on the shared
  machine, `OMP_NUM_THREADS=1`.

### D-37 Feature universe: a dedicated extraction with the pre-filter recomputed · `scripts/31`
- **Question asked:** does any cohort-level feature pre-filter see the subtype labels?
- **Answer:** the only cohort-level step is the D-07 extraction pre-filter. It keeps a feature if
  its support (observed count for AS, CPM ≥ 1 count for GEX) reaches ⌈fraction × the smallest
  inner-train set of the cohort⌉ and it is non-constant over the cohort. It reads values and
  cohort membership, never a label: it is label-blind.
- **But the existing store is the wrong one.** Its per-cohort columns were computed for the
  tumour-status and pan-cancer cohorts. A subtype cohort is a different sample set (BRCA subtype
  ⊂ BRCA tumour-status; LGG, SARC, TGCT, THYM, MESO, CESC and PCPG are not tumour-status cohorts
  at all), and its training sets are smaller, so those columns are not guaranteed
  supersets of what a subtype fold could select. Using them would not leak, but it would restrict
  the feature universe in a way the LEAKAGE_POLICY proof does not cover.
- **So step 31 re-extracts** the raw tables over the 4,618 patients of the subtype cohorts, with
  one pre-filter column per subtype cohort computed from that cohort's own splits. The proof of
  LEAKAGE_POLICY.md then applies verbatim: every training set S of a subtype fold is a subset of
  the cohort C, support and range are monotone, so no feature a fold could keep is dropped.
- **Cohort membership** is "has an admissible label". That uses the labels only to decide who is
  in the cohort, which is the cohort's definition; every training and test set is inside it, and
  no class information about a held-out patient enters any fitted quantity. Stratification of
  the folds uses the labels of all cohort patients, as in D-04.
- **Independence of cohorts:** a job loads only its own cohort's pre-filter column, so adding or
  removing a cohort changes no other cohort's features or values.
- **Values are identical** to step 02's for the same sample and feature (per-sample operations
  only; GEX library size over all 57,656 genes). Checked on 300 × 300 random sample-feature
  pairs of RI, GEX, A3SS and A5SS: equal to the last bit, NaN pattern included. The subtype
  store has 1,712 genes and 125-509 events per view that only the subtype cohorts' pre-filters
  admit (e.g. testis- and thymus-expressed genes).

### D-38 SHAP: what is explained · `ctsm/subtype.py: summarise_shap`
- **Two classes:** the model's single logit, the log-odds of the second class in YAML order
  (e.g. lobular in BRCA; each `metrics.json` lists the classes in order).
  `shap_top_features.csv` reports `direction_towards_<second class>`.
- **More than two:** each held-out patient's centred true-class logit zₖ − z̄ (as D-10 for
  cancer type); no pooled direction is reported.
- **By class** (`shap_top_features_by_class.csv`, top 30 per class): mean |SHAP| over the
  patients of that class; `direction_spearman` > 0 always means a higher value raises the
  evidence for that class (for the first of two classes the sign of the log-odds is flipped).
- Additivity is asserted per fold (≤ 1e-3), as D-10.

### D-39 Metrics, uncertainty and baselines · `ctsm/subtype.py: subtype_metrics, bootstrap, baselines`
- **Pooled out-of-fold:** balanced accuracy (mean per-class recall), macro F1 over all classes,
  macro one-vs-rest AUROC (the AUROC itself for two classes), accuracy, MCC, log loss and
  per-class recall, precision and F1; confusion matrix.
- **95% CIs:** 1,000 patient-level bootstrap resamples (percentile); a resample missing a class
  is skipped and the number used is recorded.
- **Baselines, per cohort:** chance balanced accuracy 1/C, chance macro AUROC 0.5, chance macro
  F1 1/C (guessing each class at its prevalence); the majority-class predictor's accuracy,
  balanced accuracy (1/C) and macro F1.

### D-40 AS ensemble, fusion and tests · `scripts/34`
- As D-14/D-23, per cohort, from out-of-fold probabilities only:
  - **AS ensemble:** the mean of the 8 splicing views' class probabilities;
  - **GEX + AS (1:1):** 0.5 × GEX + 0.5 × AS ensemble. The weights were fixed before any result.
  No learned stacking is reported as a headline.
- **Tests against GEX:** exact McNemar on paired per-patient correctness (argmax); and, because
  McNemar weighs classes by their size, a paired patient-bootstrap 95% CI of the difference in
  balanced accuracy (1,000 resamples).
- **FDR family:** BH across the cohorts of the subtype task, separately for each comparison (AS
  ensemble vs GEX; GEX + AS vs GEX). Figures show raw p; q is in the CSVs with its family named.

### D-41 Interpretation, fixed before any SHAP value was read · `interpretation`, `literature_markers` · `scripts/35`
- **Top features per class and view:** `shap_top_features_by_class.csv` of each job; a gene's
  importance is its best event (`biology.gene_ranking`).
- **Do splicing and expression read the same genes?** For each cohort and AS view, the top 50
  genes of the AS model against the top 50 of the GEX model, within the genes both could select;
  overlap against its hypergeometric expectation.
- **Literature markers** (`literature_markers`, with sources): for each, (i) is the expression
  difference between the two named classes in the expected direction (AUROC of the gene's log2
  CPM, Mann-Whitney p) — a property of the data, no model involved; (ii) does the gene rank in
  the GEX model's top 100 for either class; (iii) does any event of the gene rank in an AS
  model's top 100. LGG's "1p/19q dosage" is, per patient, the median over the genes of 1p (and
  separately 19q) of the gene's log2 CPM minus its median over the cohort; codeletion predicts it
  is lower in codeleted tumours. Arms are split at the hg38 centromeres (chr1 123.4 Mb, chr19
  26.2 Mb).
- **Single events: the neighbour-gene audit** (as D-27). Before an event is named as a splicing
  finding, its segment (SE/AFE/ALE/HIT: the exon; MXE: both exons; RI: the intron; A3SS/A5SS:
  the extension) is intersected with every exon of every other GENCODE v45 gene (≥ 50 bp), and
  classed antisense, same-strand, readthrough or none. An overlapping event *tracks* the other
  gene if |ρ(event, other gene's log2 CPM)| ≥ 0.5 and exceeds |ρ(event, host log2 CPM)|, over
  the cohort's patients. A tracking event is reported as a readout of the other gene, not as a
  splicing choice of its host. The audit is reported for every named event, whatever it finds.
- **Which events are named:** the audit covers the top 5 events per class and AS view of every
  cohort. The figure of single events (U6) shows, for each histology and genomic/IHC cohort, the
  one AS event with the largest class-specific mean |SHAP| across the 8 views; the rule uses the
  importance only, never the audit's outcome. Any other event named in the text comes from the
  audit table and is given with its audit status; the table also holds each literature marker's
  best event (A-2).
- **Cross-family errors** (descriptive): where a histology cohort has a molecular counterpart in
  the same project (LGG and LGG_IDH, UCEC and UCEC_TCGA, STAD and STAD_TCGA, BRCA and BRCA_PAM50
  / BRCA_ER, PCPG and PCPG_mRNA), the histology models' out-of-fold calls (GEX, AS ensemble) are
  cross-tabulated against the molecular class of the same patients: are the "errors" of a
  histology model molecularly coherent? Counts only; no test.

### D-42 What is published · D-12
- **Committed:** `config/subtype.yaml`, code, this file, `results/subtype/` (per job: metrics,
  folds, per-class metrics, confusion matrix, SHAP summaries; per task: label counts, scope,
  concordance, cross-tabulations, extraction summaries, performance and interpretation tables),
  `results/figures/subtype/` and `results/SUBTYPE.md`.
- **Local only** (`<data_dir>/subtype/`): the per-patient label table, task tables, splits,
  matrices, models, SHAP arrays and per-patient predictions, and the PanCanAtlas tables
  (`<data_dir>/external/`). Committed outputs hold counts and aggregates only.

### D-43 Molecular and IHC label sources · `pancan_subtype`, `pancan_clinical`
- **Files** (open access, GDC publication pages; fetched 2026-10-12 and MD5-checked by
  `scripts/00b_fetch_gdc_labels.py` into `<data_dir>/external/gdc_pancanatlas/`, with a
  `MANIFEST.json` that the loaders re-check before reading):

  | File | URL | Bytes | MD5 |
  |---|---|---|---|
  | `TCGASubtype.20170308.tsv` (PanImmune page) | https://api.gdc.cancer.gov/data/0f31b768-7f67-4fc4-abc3-06ac5bd90bf0 | 603,836 | d0a6054b799a8a61488636ad400c3af9 |
  | `clinical_PANCAN_patient_with_followup.tsv` (PanCanAtlas page) | https://api.gdc.cancer.gov/data/0fc78496-818b-4896-bd83-52db1f533c5c | 18,633,685 | ffcb35edda305dd8d615497f9214eb92 |

- **Handling:** downloaded data are untrusted: read with pandas only (`dtype=str`; the clinical
  table with `encoding="latin-1"`), nothing in that directory is executed, per-patient rows are
  never committed.
- **Sample types:** `pan.samplesID` holds patient barcodes (12 characters) for some cancers and
  sample or aliquot barcodes for others. Only patient rows and primary-tumour rows (sample type
  01; 03 for LAML) of the cohort's cancer type (`AML` = LAML, `OVCA` = OV) are used. This
  matters: 114 BRCA rows are solid-tissue normals, all but one called "BRCA.Normal", and 266
  SKCM rows are metastases. A patient with two usable rows that disagree is left unlabelled
  (none occurred).
- **Published labels were assigned once over the whole TCGA cohort** (clusterings and centroid
  calls fitted on all patients, including those that land in our test folds). That is a property
  of the label, not leakage in this pipeline: our feature selection, tuning and training still
  see training folds only, and the labels are targets, not inputs.

### D-44 Molecular and IHC cohorts · `label_sets.molecular`
`label_type` says how each label was made. **Genomic** (DNA: mutation, copy number,
methylation, fusion) and **IHC** (protein or in-situ hybridisation) labels are tests of the
transcriptome that do not depend on it. **Expression-defined** labels were derived from mRNA
expression: a GEX model predicting them is close to circular, and the question there is how
close the splicing views get.

| Cohort | Type | Column | Classes (patients) | Notes |
|---|---|---|---|---|
| LGG_IDH | genomic | `Subtype_Selected` | IDHmut-codel 171, IDHmut-non-codel 242, IDHwt 94 | Ceccarelli et al. 2016 (*Cell* 164:550) groups merged to IDH/codel status: Codel; G-CIMP-high + G-CIMP-low; Mesenchymal-like + Classic-like + PA-like + LGm6-GBM. |
| UCEC_TCGA | genomic | `Subtype_Selected` | POLE 24, MSI 34, CN-low 28, CN-high 83 | TCGA 2013. |
| STAD_TCGA | genomic | `Subtype_Selected` | CIN 196, MSI 55, GS 43, EBV 26 | Liu et al. 2018 (*Cancer Cell* 33:721); HM-indel = MSI. EBV was called from viral sequence reads; GS is the remainder. |
| COAD_TCGA | genomic | `Subtype_Selected` | CIN 162, MSI 40, GS 33 | As STAD. |
| PRAD_ETS | genomic | `Subtype_Selected` | ERG 152, ETV1 28, SPOP 37, Other 86 | TCGA 2015 (*Cell* 163:1011). Fusion calls used RNA-seq reads as well as DNA, and an ETS fusion raises the partner's expression; "Other" is a residual class. |
| BRCA_ER | IHC | `breast_carcinoma_estrogen_receptor_status` | ER-negative 235, ER-positive 790 | Indeterminate and not evaluated excluded. The assay measures the ESR1 product, so ESR1 expression is expected to separate it. |
| HNSC_HPV | IHC | `hpv_status_by_p16_testing` + `hpv_status_by_ish_testing` | HPV-negative 79, HPV-positive 33 | A test not done is ignored; every test done must agree (1 discordant excluded); 388 patients have neither test. HPV-positive tumours are mostly oropharyngeal, so site is confounded with the label. |
| BRCA_PAM50 | expression-defined | `Subtype_Selected` | LumA 558, LumB 212, Her2 77, Basal 187, Normal-like 39 | Parker et al. 2009 (*J Clin Oncol* 27:1160); TCGA 2012 (*Nature* 490:61). Normal-like largely reflects normal-tissue admixture. |
| HNSC_mRNA | expression-defined | `Subtype_Selected` | atypical 68, basal 86, classical 49, mesenchymal 75 | TCGA 2015 (*Nature* 517:576). |
| PCPG_mRNA | expression-defined | `Subtype_Selected` | kinase signalling 68, pseudohypoxia 61, Wnt-altered 22, cortical admixture 22 | Fishbein et al. 2017. |

- **Column choice:** `Subtype_Selected` throughout, except where it does not hold the label
  wanted: for LGG it mixes IDH/codel status with methylation clusters, so the classes are merged
  to IDH/codel (above); BRCA ER and HNSC HPV come from the clinical table.
- **Not run, and why:** expression-defined subtypes of LUSC, LUAD, KIRC, OV, BLCA, THCA, LAML and
  UCS (three expression-defined cohorts answer the question "how close does splicing get";
  more would repeat it); integrative clusters of KIRP, LIHC, LUAD and ACC (mixed label types);
  GBM (methylation subgroups of 97 IDH-wildtype tumours); ESCA (the GI.ESCC/GI.CIN split is the
  histology contrast again); SKCM (61 primary tumours, one class ≥ 20); KICH eosinophilic (19).
- **Cross-family tables** (`label_crosstab.csv`, patients labelled in both, no model involved)
  show how the label families relate: e.g. 48 of 53 UCEC serous tumours are CN-high; every
  paraganglioma is pseudohypoxic or kinase-signalling and every cortical-admixture tumour is a
  pheochromocytoma.

### D-45 Literature markers of the molecular cohorts · `literature_markers`
- As D-41, keyed by cohort: MLH1 lower in MSI (UCEC, STAD, COAD), CD274 higher in EBV gastric
  cancer, ERG and ETV1 higher in their fusion classes, ESR1 in ER-positive, CDKN2A (p16) in
  HPV-positive, ESR1/ERBB2/KRT5 for the PAM50 classes, CYP11B1 for cortical admixture and PNMT
  against pseudohypoxia, 1p/19q dosage for codeletion. Several are near-tautological (ESR1 for
  ER IHC, CDKN2A for p16, ERG for ERG fusions) and are there as positive controls of the check.

### D-46 Splicing regulators of the subtypes · `config/followup.yaml: subtype_regulators` · `scripts/37`
Added 2026-10-12, after every subtype model was summarised and before any regulator was looked
at by class. Three parts, all on the cohort's tumours, no model refitted:
- **Scan.** Each of the 70 splicing factors of D-18 (log2 CPM), between the cohort's classes:
  Mann–Whitney for two classes, Kruskal–Wallis for more; BH within the cohort. Effect: each
  class's AUROC against the rest.
- **Association** (as D-18). The top 50 events of the cohort's model in SE, A3SS, A5SS, MXE and
  RI against 500 random measured events of the view: partial Spearman with every factor,
  adjusting for the stromal and immune scores; excess = share of top events with |ρ| > 0.3 minus
  the share among the background.
- **Predictions** R1–R6, fixed in the config with their sources: ESRP1 and ESRP2 lower in
  biphasic than epithelioid mesothelioma, ESRP1 lower in the mesenchymal HNSC class (Warzecha et
  al. 2009, 2010); PTBP1 lower in leiomyosarcoma than in the other sarcomas, and, adjusting for
  class and composition, PTBP1 moves the ACTN1 and TPM1 mutually exclusive exons against the
  leiomyosarcoma direction (Southby et al. 1999; Gooding et al. 1998). BH over R1–R6.

### D-47 Is the subtype signal a tissue-source-site signal? · `config/followup.yaml: subtype_site_probe` · `scripts/38`
The site probe of D-20, applied to every subtype cohort: the same untuned pipeline
cross-validated with patients grouped (5 folds) and with tissue-source sites held out (5 folds
of sites), three seeds, all nine views. A cohort is probed only if every class comes from ≥ 5
sites (SARC is not: myxofibrosarcoma comes from 4). The site–class confounding of each cohort
is Cramér's V between TSS and class, from the labels alone. Predictions: P1, the AS views lose
more than GEX when sites are held out in most cohorts (sign test), as they did between
neighbouring projects in D-20; P2, the loss rises with Cramér's V (Spearman over cohorts).
The probe is a diagnostic: it never touches a subtype model or its estimate.

### D-52 Does genome sequence say which exons MSI tumours skip? · `config/followup.yaml: msi_cis` · `scripts/42`
Added 2026-10-12, after §8 of `results/SUBTYPE.md` and before any PSI value was read by repeat
class. Mismatch-repair-deficient tumours shorten mononucleotide repeats; a shortened poly-T tract
in an acceptor's polypyrimidine tract weakens the 3′ splice site and the exon is skipped (MRE11
exon 5, Giannini et al. 2002; HSP110 exon 9, Dorard et al. 2011; genome-wide in colorectal
cancer, Jonchère et al. 2024). §8 found MSI splicing marks shared by three organs; this asks
whether a cis route explains part of them, choosing exons from the genome alone.
- **Label-blind annotation** (`--only annotate`, GRCh38 primary assembly). For every cassette exon
  of the subtype SE matrix, the longest sense-strand run of each base overlapping three intronic
  windows: acceptor −1..−50, deep −201..−400 (intron ≥ 500 nt) and donor +1..+50. 99.9% of
  exons read AG at −2,−1, which checks the coordinates; MRE11 exon 5 reads its T11.
- **Groups.** UCEC_TCGA, STAD_TCGA and COAD_TCGA as in the subtype task: MSI against every other
  class except POLE; the UCEC POLE tumours (hypermutated, mismatch-repair proficient) as a
  negative control, without the one called MSI-H by PCR.
- **Score** (per patient, no model): each exon's PSI as a percentile within the cohort; the mean
  over repeat exons (acceptor T-run ≥ 10) minus the mean over background exons (≤ 5), the
  background weighted to the repeat exons' mean-PSI deciles so the global PSI axis (D-24) cancels.
- **Predictions** M1–M7 and K1 (config): the score is lower in MSI (M1); per-exon MSI shift falls
  with repeat length (M2), also after adjusting for the acceptor's pyrimidine content (M3); the
  acceptor window beats the deep window (M4); MSI scores below POLE (M5); the trained MSI models'
  top SE events are enriched for acceptor T-runs (M6); the repeat exons' shifts agree between
  organs (M7); MRE11 exon 5 is skipped in MSI (K1). FDR families M (19 tests) and K (3).
  HSPH1 exon 9 is not a cassette event of either SE matrix and cannot be tested.

### D-53 Long U-tracts: can MSI raise inclusion? · `config/followup.yaml: msi_long_tracts` · `scripts/42 --only tracts`
Added 2026-10-12 after the D-52 results were read, before this relation was computed. The two
U7 loci more included in MSI (CCNDBP1, ANXA7) carry 30- and 29-nt acceptor T-runs, the longest
of the U7 events, and the spread of the MSI shift widens for runs ≥ 14 in both directions.
hnRNP C binds long U-tracts and competes with U2AF65 at the acceptor, repressing the exon
(Zarnack et al. 2013). If MSI shortens such a tract it may relieve repression, while a moderate
tract that U2AF65 needs is weakened. Repressed exons sit at low baseline inclusion, so the
predictions are about baseline: N1, among repeat exons the background-adjusted MSI shift falls
with the exon's mean PSI (Spearman < 0); N2, more so than for deep-intron T-runs; N3, more so for
runs ≥ 20 than for runs 10–14. The nine MSI-class SE events of §8 are left out. BH over N1–N3
and the three cohorts (9 tests).

### D-54 Replication in alternative 3′ splice sites · `config/followup.yaml: msi_a3ss_replication` · `scripts/42 --only a3ss`
Added 2026-10-12 after the D-52 and D-53 results were read, before any A3SS value was read in
this analysis. Exploratory tables of D-53 showed SE repeat exons with acceptor T-runs ≥ 20 moving
toward inclusion in MSI in all three cohorts, while runs of 10–19 move toward skipping. An A3SS
event has two acceptors, and the long form's PSI is the use of the distal one, so a run at the
distal acceptor and a run at the proximal one predict opposite signs. Oriented shift = the
background-adjusted shift × (+1 distal, −1 proximal); negative means the acceptor carrying the
run is used less in MSI. Q1: runs 10–19 at one acceptor (other ≤ 6), oriented shift < 0. Q2:
runs ≥ 20, oriented shift > 0. Wilcoxon signed-rank per cohort, BH over 6 tests.

### D-55 Is the MSI repeat score tumour-intrinsic? · `config/followup.yaml: msi_robustness` · `scripts/43`
Added 2026-10-12 after D-52..D-54 were reported, before any of this was computed. U9 showed that
splicing models lose accuracy at unseen hospitals, STAD_TCGA among them, so the repeat score is
checked against the hospital, tumour composition, RNA handling and the feature filter. S1: within
hospital (MSI–rest pairs from the same tissue-source site; permutation within site). S2: adjusted
for site, stromal and immune expression and the sample's mean retained-intron PSI. S3: MSI is
somatic, so the patient's own adjacent normal should not carry it; the tumour-minus-normal
difference is compared between patients with and without an MSI tumour (9 against the others,
three projects pooled, scores built in the main matrix). S4: the score rebuilt in the main matrix
with its own feature filter. Family S: STAD and COAD (7 tests); UCEC reported alongside.

### D-56 Reading the repeats in the tumours' DNA · `config/followup.yaml: msi_dna` · `scripts/44`, `scripts/45`
Added 2026-10-12 after D-52..D-55, before any read was fetched. D-52..D-55 inferred repeat
shortening from the reference genome and the MSI label; this reads it. 785 loci (label-blind:
every SE acceptor T-run ≥ 20, seeded samples of 250 runs of 14–19 and of 10–13, the A3SS repeat
acceptors, and the clinical MSI markers BAT-25 and BAT-26) are sliced from harmonised exome BAMs
of every MSI and POLE tumour with its normal and 120 MSS tumours (60 with normals), in that
order, stopping at 6 GB. The data user fetches the slices with their own GDC token; the
repository never handles it. Each read's run length is called between exact 8-nt flanks, with
no trust in the aligner's indel placement. Contraction per tumour and locus is the share of reads
shorter than the germline mode, minus the same share in the normal. Predictions G1–G6: MSI
tumours carry more contraction (G1), POLE less than MSI (G2); within MSI tumours, a shorter allele
goes with less inclusion of its exon at runs of 10–19 (G3) and with more at runs ≥ 20 (G5); loci
that contract more skip more (G4); contraction rises with run length (G6). Family D (8 tests).


### D-57 Inherited splicing variation in the BRCA subtypes · `config/followup.yaml: brca_germline` · `scripts/46`
Added 2026-10-12, before any test was computed, with three open GDC tables and the MC3 calls
fetched for it (step 00b; about 1 GB, approved by the data user). Many splicing events are set by
a common germline variant: their PSI falls in three evenly spaced groups, the two homozygotes
and the heterozygote between them. A label-blind scan finds these events in the BRCA matrices.
Two routes could carry them into a subtype signal. (1) **LOH**: when a tumour loses one
homologue, a heterozygote's PSI moves toward the retained allele, more so the purer the tumour
(L1–L3, in the 108 tumour–normal pairs with ABSOLUTE calls). If splicing reads LOH, the subtypes'
different LOH maps show in their splicing (L4, LumA against Basal; L5, 16q in lobular tumours;
L6, arm by arm against the DNA call). (2) **Ancestry**: GSTM1, a deletion polymorphism, is a top
ER feature of the GEX model; ER-negative disease is more common in Black patients, who less often
carry the deletion (A1, A2); genotype-like splicing events differ in allele frequency by race
(A3). GRCh37 tables are joined to GRCh38 events through each gene's MC3 call span (no liftover is
installed). Families L (6 tests) and A (3).

### D-58 Two meanings of HER2 · `config/followup.yaml: brca_her2` · `scripts/47`
Added with D-57. The PAM50 Her2 class is an expression programme; clinical HER2 is ERBB2
amplification. The GEX model's top Her2 genes are the 17q12 amplicon; the AS models' are not, and
a PSI does not move with copy number. H1–H4 ask which HER2 each modality learned, with
amplification from ABSOLUTE copy number at ERBB2 (≥ 2 × ploidy). Family H (4 tests). Re-tested
without ABSOLUTE's ceiling in D-60.

### D-59 Splice-site mutations and SF3B1 in BRCA · `config/followup.yaml: brca_mutations` · `scripts/48`
Added with D-57. MC3 splice-site SNVs are mapped to the cassette exon they hit (C1, every gene; C2,
CDH1, the open item of the lobular event card); SF3B1 hotspot tumours should use cryptic
acceptors 10–30 nt upstream (C3). Family C (3 tests).
### D-60 HER2 amplification without the ABSOLUTE ceiling · `config/followup.yaml: brca_her2_uncapped` · `scripts/49`
Written 2026-10-13, after D-58's results and before the file was downloaded. ABSOLUTE never reports
more than 7 copies at ERBB2, so D-58's amplified group was a third of the clinically HER2-positive
tumours. D-60 re-reads amplification from the PanCanAtlas SNP6 copy ratio, made absolute with
ABSOLUTE purity and ploidy and given no ceiling, and keeps D-58's rule (≥ 2 × ploidy). HU1–HU5
repeat H1–H4 (HU2 now predicts D-58's result, splicing ahead); HU3 is its independent part, on the
tumours only the new call makes amplified. HU6–HU7 ask whether P(Her2) follows the amplification
level, which the amplicon-reading GEX model should and a programme-reading splicing model should
not. Family HU (7 tests).

### D-61 to D-64 The MSI acceptor repeats, followed up · `config/followup.yaml` · `scripts/50`
Written 2026-10-13, before any of these tests was computed. They follow the questions of
[MSI_REPEATS_LITERATURE.md](MSI_REPEATS_LITERATURE.md).
- **D-61, deletion load (DL1–DL6).** From MC3: does each MSI tumour's homopolymer-deletion load
  set its repeat score, adjusted for purity? Are endometrial MSI deletions shorter, as published
  marker shifts suggest?
- **D-62, frame and NMD (NM1–NM4).** When a coding repeat exon that MSI skips is frameshifting,
  does its gene's mRNA fall?
- **D-63, branch-point geometry (BP1–BP4).** A long T-run pushes the branch point away from the
  AG. Does the MSI shift follow the branch-point distance, beyond run length? This is the
  candidate explanation of the reversal at ≥ 20 nt.
- **D-64, germline indels (GL1–GL3).** Do GTEx v8 germline poly-T indels at the same exons move
  splicing the way MSI contraction does, with the same reversal? The Hamanaka et al. catalogue
  genotypes no homopolymer, so it could not be used.
- **Data.** GTEx v8 significant sQTL pairs for Colon_Transverse, Stomach and Uterus: three members
  of `GTEx_Analysis_v8_sQTL.tar` (storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl),
  24 MB, fetched 2026-10-13 and pinned by size and md5 in `scripts/50`. The six Hamanaka et al.
  files (Zenodo record 7086007, the three tissues' `allpairs` and `genes` tables, 146 MB) were
  fetched and md5-checked the same day, then not used.
- **Outcome, recorded after running.** GL1–GL3 are not testable: the GTEx pairs hold no variant
  inside any of the 1,580 runs (`followup/germline/variant_coverage.csv`). No amendment was needed;
  every other test ran as written.

### D-65 to D-67 The MSI sequence code · `config/followup.yaml` · `scripts/51`
Written 2026-10-13, after D-61..D-64 and before any of these tests was computed. Nothing outside
UCEC_TCGA, STAD_TCGA and COAD_TCGA had been read for MSI: no splicing value, MC3 count or mismatch-
repair call of any other tumour. The claim they build toward: mismatch-repair deficiency moves 3′
splice-site use by a code in the reference genome, tract length and branch-point distance, that
says exon by exon whether inclusion falls or rises. Read in RNA, that code should find
mismatch-repair-deficient tumours in organs never used to define it.
- **The rule.** Fixed from D-53 and the U17b table; STAD and COAD show the pattern on their own.
  Runs of 10–19 nt fall; runs ≥ 20 nt with the branch point ≥ 36 nt from the AG rise. 36 is the
  first tertile of that distance among runs ≥ 20, a label-blind cut. Score G subtracts the rising
  arm from the falling one.
- **D-65, specificity (BS1–BS4).**
  - Does the branch point add beyond where the run ends?
  - Is the effect absent from exons without a tract (B exons)?
  - Is it absent from mismatch-repair-proficient contrasts (POLE against the rest in UCEC; CIN
    against GS in STAD)?
- **D-66, thirty held-out cancers (PC1–PC8).**
  - DNA truth: mismatch-repair deficiency is called from MC3 as ≥ 10 homopolymer deletions, the
    D-61 count. In the three cohorts it puts 129 of 129 MSI tumours and 2 of 571 others above
    the cut.
  - Scores A and G are computed with no re-fitting, within each cancer.
  - Tested on these tumours: tumour-level detection, the advantage of G over A, organ independence
    (READ, ESCA and UCS left out), and mismatch repair against mutation burden (hypermutated
    proficient tumours).
  - Tested on the exons: replication of the exon shifts and of the direction rule.
- **D-67, endometrium (EG1–EG2).** Does score G separate UCEC MSI, where score A could not?

### D-68 and D-69 The code against infiltration and against independent MSI calls · `config/followup.yaml` · `scripts/51`
Written 2026-10-13, after D-65..D-67 and before either was computed.
- **D-68, infiltration (IF1–IF6).** MSI tumours are heavily infiltrated, and D-55 ruled out
  infiltration for score A only. Does score G keep its separation in the 30 held-out cancers, and
  in the three cohorts, after stromal and immune expression, the retained-intron level (an
  RNA-handling proxy) and, where called, ABSOLUTE purity? Is it kept within the less-infiltrated
  half of each cancer?
- **D-69, MSI calls not from MC3 (MC1–MC3).** D-66's DNA truth was a count of MC3 homopolymer
  deletions. MANTIS (Bonneville et al. 2017) and MSIsensor call MSI from the exome reads by
  different means. Does score G follow those calls in the 30 held-out cancers? D-69 runs only if
  the user approves the download: per-sample attributes from the cBioPortal PanCancer Atlas
  studies, about 20,000 records.
- **Outcome, D-68 (computed 2026-10-13; IF1–IF6 supported).** Held out, the MMR coefficient on
  the G percentile is −0.339 adjusted for stromal, immune and retained-intron levels (−0.335
  unadjusted; IF1) and −0.334 with purity (IF2). In the lower-immune half of each cancer the
  AUROC is 0.84 (21 MMR-deficient tumours; IF3). In UCEC, STAD and COAD, with site fixed effects
  and HC3 errors, the MSI coefficient keeps 1.01–1.07 of its site-only value (IF4–IF6). The
  immune score alone gives AUROC 0.59, and its adjusted coefficient on G is positive (+0.055),
  against the MMR effect. Tables: `results/subtype/msi_cis/code/infiltration/`; figure U20.
- **Outcome, D-69 (computed 2026-10-13; MC1–MC2 supported, MC3 not resolved).** Download
  approved by the data user: 64 cBioPortal REST requests, 6.2 MB of JSON (above the 5 MB first
  sized), kept as a 320 KB TSV pinned by MD5 (step 00c). Of the 7,685 held-out tumours, 7,388
  have a MANTIS and 7,470 an MSIsensor score. Score G separates MANTIS MSI-H from MSS at AUROC
  0.97 (17 MSI-H; MC1) and MSIsensor MSI-H at 0.89 (23; MC2). Its Spearman with the MANTIS
  percentile over all scored tumours is +0.011 (MC3, not resolved). All 17 MANTIS MSI-H tumours
  are MC3-deficient; MANTIS calls 27 of the 50 MC3-deficient tumours mid and 6 MSS. Tables:
  `results/subtype/msi_cis/code/msicalls/`; figure U21.

### D-70 Is the rising arm real? · `config/followup.yaml` · `scripts/52`
Written 2026-10-13, after D-68 and D-69 and before any of these tests was computed. What had been
seen is listed in the config header. The claim under test: in mismatch-repair-deficient tumours,
exons after a T-run of 20 nt or more, with the branch point at least 36 nt from the AG (R_up),
become *more* included, and this happens because contraction shortens a branch-point-to-AG distance
that is too long.
- **Specificity (RV1, RV7).** Do R_up exons rise in MSI more than in any of about 60 other subtype
  contrasts? Do they rise more than in hypermutated tumours whose repair is intact? Every contrast is
  drawn at 20 against 20 tumours, so noise is equal.
- **A tumour-wide axis (RV2).** Is the rise left once each tumour's average inclusion is in the
  model?
- **The tract is needed (RV3), and its side matters (RV5).** Exons with a far branch point but no
  tract should not rise. Long runs that lie *upstream* of a branch point should not rise either,
  because contracting them cannot move the branch point closer to the AG.
- **Geometry without a predictor (RV4).** The AG exclusion zone bounds where a branch point can lie
  (Gooding et al. 2006) and needs no branch-point model. Its length should predict the rise among
  long runs.
- **Dose (RV6).** The rise should grow with the tumour's homopolymer deletion count.
- **The alternative (RV8).** PTB represses many exons with long tracts and distant branch points. If
  contraction relieved PTB, the exons that PTBP1 represses most should rise most.
- **Not testable here.** The public MC3 file holds almost no intronic deletions (458 A/T deletions
  within 50 nt of any exon), so the tumour's own allele at each locus cannot be read. Only 2 R_up
  introns have a retained-intron event, so a read-out from intronic reads is not possible either.
- **Outcome (computed 2026-10-13; RV2, RV3 and RV7 supported, RV8 against, RV1 and RV4–RV6 not
  resolved).**
  - **The rise is real.** Held out, R_up exons move up in 98 of 135; the median shift is +0.028.
  - **Not a tumour-wide axis or infiltration.** The MMR coefficient on the R_up arm is +0.236 with
    the global inclusion axis, stromal, immune and RI levels in the model, against +0.237 with
    cancer alone (RV2).
  - **It needs the tract.** Exons with a far branch point but no run do not move (median −0.0004;
    AUROC 0.69; RV3).
  - **Mismatch repair, not burden.** Against hypermutated proficient tumours the AUROC is 0.74
    (RV7); with MANTIS MSI-H it is 0.87 (exploratory).
  - **Not unique to MSI.** The held-out rise (+3.5 percentile points at 20 against 20) beats 55 of
    57 other subtype contrasts (RV1, q = 0.083). Sarcoma LMS against MFS (4.2) and LMS against
    UPS (3.6) move these exons as much.
  - **Geometry not confirmed.** The AG exclusion zone does not predict the rise (ρ −0.04; RV4), and
    it does not track the predicted branch point either (ρ −0.04). Runs upstream of a branch point
    rise about as often (26 of 42; RV5). Within MMR-deficient tumours the rise does not follow
    deletion load (ρ +0.08; RV6).
  - **PTB relief is ruled out.** The more PTBP1 represses an exon, the *less* it rises (ρ +0.20
    against the predicted negative; RV8).
  - Tables: `results/subtype/msi_cis/code/reversal/`; figure U22.

### D-71 Fitted models on the code · `config/followup.yaml` · `scripts/54`
Written 2026-10-13, after D-70 and A-16 and before any model was fitted. Score G has no fitted
weight: each arm counts once, and every exon in an arm counts equally. Two questions:
- **Do weights help?** Can a few weights fitted on top of the code beat it in cancers never used?
- **Does the sequence prior help transfer?** Does learning a weight per exon, as the subtype models
  do, transfer worse than learning a weight per sequence class?

**The ladder.** All models share one label (≥ 10 MC3 homopolymer deletions) and one split:
- **M0.** Score G, unfitted.
- **M1.** Logistic regression on the two arms.
- **M2.** Logistic regression on ten sequence bins: four run-length bins for runs of 10–19 nt, and
  runs ≥ 20 by 10-nt branch-point bins.
- **M3.** An MLP on the same ten bins.
- **M4.** Penalised logistic regression on each repeat exon. It keeps the sequence prior but
  learns one weight per exon.
- **M5 and M6.** Penalised logistic regression and an MLP on all exons, with no sequence prior: the
  pan-cancer cousin of the subtype models.

**Split.** Training uses the three cancers the code came from; testing uses the 30 held out
(primary). Hyperparameters are chosen by leaving one training cancer out, so tuning rewards
transfer, not fit. Leave-one-cancer-out over all 33 is a second, exploratory estimate.

**Tests.** LR1–LR7 compare AUROC up and across the ladder, and AUPRC of M2 against G, by paired
bootstrap of the held-out tumours.

**Seen when written.** The 30 held-out cancers have been read many times (D-66 to D-70, A-15,
A-16), and the 10-nt branch-point bins were seen to split held out (A-16). The bins themselves
were set on the training cancers (U19a), and no weight is fitted on a held-out tumour.

**Outcome (computed 2026-10-13; LR2 and LR7 supported, LR1 and LR3–LR6 not resolved).**
- **Finer bins help; weighting the arms does not.** Logistic regression on the two arms matched
  score G held out (AUROC 0.835 against 0.832; LR1). On the ten bins it reached 0.868 (LR2,
  q = 0.046), and its AUPRC rose from 0.083 to 0.257 (LR7). The bin weights reproduce the code
  without being told it: every run of 10–19 nt, and runs of 20 or more with the branch point
  20–29 nt away, get negative weights; runs of 20 or more with the branch point 30 nt or further
  away, or none found, get positive weights (U23d).
- **An MLP adds nothing** on the ten bins (0.873; LR3).
- **Per-exon models did not transfer worse.** This is against the prediction's sign, although not
  resolved. Repeat-exon logistic regression reached 0.900, all-exon 0.902 and the all-exon MLP
  0.885, against 0.868 for the bins (LR4–LR6, q ≥ 0.28). Their AUPRC was 0.44–0.49.
- **Leaving one cancer out (exploratory).** AUROC 0.957 (all-exon) and 0.967 (repeat-exon) over
  all 33 cancers, and 0.905 and 0.914 on the 30. The MLP on the bins fell to 0.851.
- **The per-exon gain comes from the code (A-18).**
  - The all-exon model's 500 exons hold 6 times the pool's share of repeat exons (32% against 5%).
    82% of its run-10–19 exons fall in MMR-deficient training tumours, and 13 of its 15 rising-arm
    exons rise. MLH1's best exon ranks 732nd, outside the 500.
  - Without repeat exons the model falls to 0.777. Per-exon weights on the two arms' exons alone
    give 0.888; on the 107 other repeat exons, 0.776.
- **MANTIS labels (exploratory).** With MANTIS MSI-H against MSS as the held-out label (17 against
  7,321), the AUROC is 0.997 for the all-exon model and 0.970 for score G.
- Tables: `results/subtype/msi_cis/code/models/`; figure U23.

---

### D-72 Benchmark against the exome MSI callers · `config/followup.yaml` · `scripts/55`
Written 2026-10-13, after D-71's E1 results and before any RNA model other than score G was scored
against MANTIS or MSIsensor. D-71's M5 reads only splicing; MANTIS and MSIsensor (D-69) read
microsatellites in the exome, as does the MC3 deletion count that defines our label. None is a gold
standard.

**Design.**
- **Triangulated references.** Each DNA method is judged against a reference built from the other
  two: R_MANTIS needs MC3 and MSIsensor to agree, R_MSIsensor needs MC3 and MANTIS, R_MC3 needs
  MANTIS and MSIsensor. M5 is judged against the same reference on the same tumours.
- **Common set.** Every method is scored on the tumours that have both exome scores.
- **Metric.** Within-cancer AUROC: a positive is only ranked against negatives of its own cancer.
  This needs no calibration across cancers, or across the 33 models of leave-one-cancer-out.
- **Split.** Leave-one-cancer-out (E2) is primary here. It is the only split with enough reference
  positives and the only one that scores COAD, STAD and UCEC, where MSI testing is routine. The
  30 held-out cancers (E1) are secondary.

**Tests.** BM1–BM6 predict that M5 ranks below each DNA tool, in E2 and in E1. The DNA callers
measure the instability itself and share the exome, with its coverage artefacts, with the
reference, so correlated errors favour them. A non-inferiority margin of 0.05 AUROC is set here
and reported descriptively.

**Seen when written.** D-69's call concordance and score G against both callers; D-71's E1 results,
including M5's per-cancer AUROC, curves and coefficients (A-18). No DNA tool had been scored
against another, and M5 had not been scored in E2.

**Outcome (computed 2026-10-13; BM1, BM2, BM4 and BM5 supported, BM3 and BM6 not resolved).**
- **The exome measures agree almost perfectly.** Each ranks the other two's consensus at
  within-cancer AUROC 0.999–1.000. No MC3-proficient tumour is called MSI by both callers.
- **Within each cancer, the splicing model is close behind.** The all-exon model reaches 0.998–0.999
  against every reference, in both splits. It ranks below MANTIS and MSIsensor by 0.001–0.002
  (BM1, BM2, BM4, BM5; q ≤ 0.013) and does not differ from the MC3 count (BM3, BM6). It is
  non-inferior everywhere at the 0.05 margin (lowest bound −0.006). Per cancer it reaches 1.000 in
  COAD, STAD and BRCA, 0.995 in UCEC and 0.987 in CESC. Score G reaches 0.955–0.959.
- **Across cancers it is not calibrated.** Flagging as many tumours as each caller calls MSI, pooled
  over cancers, the all-exon model finds 80–84% of the reference positives, at PPV 0.80–0.84. The
  callers find 94–100%, at PPV 0.95–1.00. Its inputs are percentiles within a cancer, so a threshold
  must be set per cancer.
- **Where the DNA measures disagree (exploratory).** Of 197 tumours that MC3 calls MMR-deficient
  with both caller scores, 50 are called MSI by neither caller. The all-exon model still places 42%
  of them in its top decile within their cancer (median percentile 0.86), against 69% when both
  callers agree and 9% of proficient tumours.
- **Coverage.** Both caller scores exist for 96.3% of the 8,464 tumours, against 100% for RNA.
- Tables: `results/subtype/msi_cis/code/benchmark/`; figure U24.

---

### D-73 What the splicing model calls that the exome does not · `config/followup.yaml` · `scripts/56`
Written 2026-10-13, after D-72 and A-19 and before the DNA status of any tumour the model calls was
read. The repeat-exon model (D-71 M4) calls a tumour MSI when its decision function is above 0. Leaving
each cancer out of training, that is 168 tumours. Some may be called by no exome measure (MC3,
MANTIS, MSIsensor). Two readings predict different things:
- **Missed MMR deficiency.** MLH1 silenced by promoter methylation is read in RNA but not in the
  exome, so MLH1 expression should be low (RO1, RO2). Instability should sit just under the
  callers' thresholds (RO3–RO5).
- **False calls.** Neither should differ from tumours that no one calls.

Each covariate is a percentile within the tumour's cancer, compared by AUROC against tumours that
no one calls, with permutations within cancer. **Seen when written.** D-72's discordance table and
the label-blind count of calls (168); nothing about the DNA, expression or purity of a called tumour.

**Outcome (computed 2026-10-13; RO4 supported, RO1–RO3 and RO5 not resolved: mixed, mostly false
calls).**
- **Most calls are confirmed.** Of the 168 calls, 115 are called MSI by all three exome measures, 6
  by two and 12 by one. 32 (19%) are called by none; 3 lack a caller score.
- **The 32 RNA-only calls look mostly like false calls.**
  - They are scattered over 18 cancers where MSI is rare (1–4 each), and only one is in STAD.
  - They sit just above the model's threshold (median decision function 1.0, against 8.9 for
    confirmed calls).
  - Their MLH1 expression is not low (RO1, AUROC 0.47), and neither is any MMR gene's (RO2).
  - Their MANTIS scores and MC3 deletion counts are no higher than those of tumours no one calls
    (RO3, RO5).
- **Their MSIsensor scores are higher** (RO4, AUROC 0.65, q = 0.019), but far below the stable cut
  (median 0.16 against 0.04; the cut is 4). They also carry more SNVs for their cancer (AUROC 0.64,
  exploratory).
- **A few are plausible missed cases (exploratory).** These three are listed locally:
  - a STAD tumour with the highest score (6.0): 9 homopolymer deletions, one short of the cut,
    MANTIS 0.45 (the between zone) and one MMR gene at the 1st percentile of expression;
  - a hypermutated PAAD tumour (21,427 SNVs; MSIsensor 4.5, between zone) with MC3 calls in MLH1,
    MSH2, MSH6, PMS2 and POLE;
  - a BLCA tumour with MLH1 at the 1st percentile, but no DNA instability.
- **The reverse.** The model misses 68 tumours that an exome measure calls. 20 of them are called by
  all three, and 19 of those 20 are UCEC, where the code was always weakest (D-67).
- Tables: `results/subtype/msi_cis/code/rna_only/`; per-tumour rows in
  `data.nosync/ctsm_v1/msi/rna_only_tumours.parquet`.

---

### D-74 Splicing against expression, trained the same way · `config/followup.yaml` · `scripts/57`
Written 2026-10-13, after D-73 and before any expression model was scored on a held-out cancer.
MSI can already be called from expression (MSIsensor-RNA, PreMSIm and others). The repeat-exon model
is only worth having if it does something those do not. One candidate is transfer to cancers it
never saw:
- **Expression models** learn how MMR deficiency looks in their training cancers. In COAD, STAD and
  UCEC that is mostly MLH1 promoter silencing.
- **The splicing code** reads the consequence (contracted repeats), whatever the cause.

**Design.** Expression models get M4's tumours, label, split (train on COAD, STAD and UCEC; test on
the 30 others) and tuning (leaving one training cancer out):
- X1: all genes as within-cancer percentiles (M4's input scale);
- X2: all genes as log2(CPM + 1);
- X3: PreMSIm's 15 genes;
- X0: MLH1 alone, unfitted;
- C1: the repeat exons plus genes.

**Tests (family EX).**
- EX1–EX4: M4 beats X1, X2, X3 and X0 on held-out AUROC.
- EX5: M4 beats X1 on the held-out MMR-deficient tumours whose MLH1 is not low.
- EX6: C1 beats X1.

**Seen when written.** M4's scores, benchmark and discordance (D-71–D-73). Within COAD, STAD and UCEC,
the GEX subtype models already classify MSI (U1) and MLH1 is lower there (D-46). No expression model
had been scored on a held-out cancer, and no held-out MMR-deficient tumour's MLH1 had been read.

**Outcome (computed 2026-10-13; EX3, EX4, EX6 supported, EX1, EX2, EX5 not resolved, none against).**
- **Transfer is not special to splicing.** Retrained on all genes, expression came close: 0.893 on
  log CPM and 0.866 on percentiles, against M4's 0.900 (EX2 +0.008, EX1 +0.034; not resolved). Its
  largest weights went to MLH1, EPM2AIP1 and five of PreMSIm's genes.
- **Small fixed panels do not transfer.** PreMSIm's genes refitted reached 0.785 (EX3 supported);
  PreMSIm's own k-NN recipe 0.702; MLH1 alone 0.758 (EX4 supported).
- **The cause of deficiency mattered less than predicted.** MLH1 was low (bottom tenth of the
  cancer) in 28 of the 50 held-out MMR-deficient tumours. On the other 22, M4 and X1 tied (0.854
  and 0.849; EX5 not resolved), and MLH1 alone fell to 0.47.
- **The two are complementary.** They rank the deficient tumours differently (Spearman ρ 0.36).
  Together (C1) they reached 0.926 (EX6 +0.060, q = 0.001), AUPRC 0.53, and precision 0.68 at
  recall 0.5, against 0.30 for M4 and 0.39 for X1. Repeat exons carry 71% of C1's summed |weights|.
- **Exploratory: splicing stays closer to the exome callers.** Leaving each cancer out (X1 with
  E1's parameters), within-cancer AUROC against the callers' consensus was 0.998 for M4 and
  0.975–0.979 for X1. Against MANTIS MSI-H held out: 1.000 against 0.969 (X2) and 0.924 (X1).
- Tables: `results/subtype/msi_cis/code/expression/`; Figure U26.

---

### D-75 How large is the MSI splicing shift? · `config/followup.yaml` · `scripts/60`
Written 2026-10-13, after D-74. Every MSI result so far was in percentiles within a cancer, which
rank tumours but say nothing about how far an exon moves. D-75 measures the shifts in PSI units
and asks two questions about their nature:
- **Where do they happen?** Deletions in tumour cells are diluted by normal cells, so the shift
  should grow with purity among MMR-deficient tumours (MG1–MG3). A signal from the immune infiltrate
  would grow at lower purity instead.
- **How are they spread?** A deletion at one tract occurs in some MMR-deficient tumours and not
  others. Each exon's change should then sit in a subset of tumours, beyond what a uniform small
  shift would put in the tail of the proficient distribution (MG4, MG5). Background exons absorb any
  general widening of MMR-deficient distributions.

**Seen when written.** Every earlier MSI result in percentile units, and D-68's purity adjustment of
score G over all tumours. Nothing in PSI units, nothing on purity within MMR-deficient tumours, and
no tail analysis.

**Outcome (computed 2026-10-13; none of MG1–MG5 supported: MG4 against, the others not resolved).**
- **Purity (MG1–MG3).** Among MMR-deficient tumours the R_down arm leans the predicted way with
  purity (MG1 ρ = −0.15, p = 0.047, q = 0.12), but the R_up arm does not (MG2 ρ = −0.05). The
  purity × MMR interaction is zero (MG3 −0.001, 95% CI −0.010 to +0.007). Across purity fifths the
  R_down arm stays at −0.006 to −0.007 PSI, where dilution predicts a 2.5-fold change. Purity
  neither places the shift in tumour cells nor in the infiltrate.
- **Tails (MG4 against, MG5 not resolved).** R_down exons have a smaller lower-tail excess than
  background exons (−0.027, q = 0.001). The statistic was mis-specified (A-20): its location-shift
  expectation moves proficient values past 0 and 1, which inflates the expected tail of any exon
  that shifts near a bound. Background exons score −0.072 below a proficient median of 0.1 and
  +0.001 above 0.98. Where neither bound is reached, R_down exons exceed background by
  0.002–0.010. The verdict stands as computed; the reading is that the test does not show
  concentration in a subset of tumours, not that it shows the opposite.
- **Size** (exploratory). Typical shifts are a few hundredths of PSI: the median R_down exon whose
  proficient median lies within 0.1–0.9 falls by 0.021, 0.028 and 0.015 in COAD, STAD and UCEC.
  16–18% of these exons fall by ≥ 0.10 in COAD and STAD, against 4% of background exons. An
  MMR-deficient COAD or STAD tumour has about 33 more R_down exons down by ≥ 0.10 than the median
  proficient tumour of its cancer. MMR status explains a median 0.4–0.5% of an R_down exon's PSI
  variance. The largest movers: METTL8 T14 0.58 → 0.22, BOD1L1 T24 0.13 → 0.43, CCNDBP1 T30
  0.04 → 0.25. Figure U27.


### D-76 What does a rising exon do to its transcript? · `config/followup.yaml` · `scripts/61`
Written 2026-10-13, after D-75. Rising-arm exons are mostly skipped in proficient tumours (median
PSI 0.06), so MSI raises the use of exons that are normally near-silent. D-76 asks what that
inclusion does:
- **Annotation** (label-blind). Each exon is placed in its gene's canonical protein-coding
  transcript. An exon whose flanking exons are joined directly inside the coding sequence is an
  *inserted* exon. Translating through it classes it as an in-frame peptide insertion, a premature
  stop that NMD should see (a "poison" exon), or a premature stop that escapes NMD.
- **Positive control** (CQ1). Across tumours, poison background exons should be included less
  where their gene's mRNA is higher. If this fails, the read-out is not validated.
- **The rising arm** (CQ2–CQ4). The same coupling for R_up exons, and whether MSI's rise in poison
  R_up exons lowers their genes' mRNA, as D-62 asked of the skipped exons.

**Seen when written.** D-62's null NMD result for skipped exons and D-75's sizes and examples. No
consequence annotation of rising or background exons, no expression of their genes, no
normal-tissue PSI.

---

## Amendments and deviations

| # | When | What | Why | Informed by results? |
|---|---|---|---|---|
| A-1 | 2026-10-12, after the pilot's first fold metric was printed in a log, before any SHAP value or event was read | Event audit: overlap ≥ 50 bp and tracking = \|ρ\| ≥ 0.5 and > \|ρ_host\| (was: any overlap, \|ρ\| ≥ 0.4) | To use the same definition as BIOLOGY_METHODS.md D-27, written in parallel the same day | No: it concerns the event audit only and was set from D-27 |
| A-2 | 2026-10-12, after the BRCA, THCA and SARC importance tables had been read | Event audit also covers each literature marker's best event: in the AS view where the marker gene ranks best (step 3), the gene's event with the largest importance over the patients of the class the literature expects it higher in (`selected_by` = "literature marker") | The rule "audit every named event" (D-41) would otherwise not cover marker events outside a class's top 5, such as CDH1's in BRCA, which the text names | It widens the audit; it changes no model, label, threshold or ranking, and the audit's result cannot change which events are audited |
| A-3 | 2026-10-12, after the histology cohorts' cross-tabulations were read | Figures renumbered to follow the order of `results/SUBTYPE.md` (histology errors U6 → U3, AS vs GEX genes U3 → U5, single events U5 → U6) and the cross-family figure widened from LGG to the four histology cohorts with a molecular counterpart (LGG, UCEC, STAD, PCPG) | Presentation; all four tables were pre-specified in D-41 | Presentation only |
| A-4 | 2026-10-12, after the event audit of the genomic cohorts had been read | Two descriptive steps added to `scripts/35`: `error_markers` (the literature markers' median expression in each true × called cell, cells under 5 patients left blank) and `shared_classes` (every audited top event of a class defined the same way in several cohorts — MSI, CIN, GS, squamous, adenocarcinoma — measured as that class against the rest in every cohort with the class; BH across the cross-cohort tests); figure U7. Then `msi_confounds`: each audited MSI top event's Spearman with PTPRC, CD8A, CIITA and MLH1 expression, and its MSI-vs-rest AUROC within the half of the cohort with the lower PTPRC (added after U7 was drawn, to ask whether the MSI events read immune infiltrate) | The same CCNDBP1 exon appeared in the MSI top 5 of three cohorts' models; the step asks the question for every event of every shared class, not only the one noticed | Yes: **exploratory**, motivated by an observation. The cross-cohort test is a replication of events chosen in another cohort, never in the cohort that chose them; no model, label or threshold changed |
| A-5 | 2026-10-12, after R1–R6 were computed | R6's predicted sign (PTBP1 moves the TPM1 exon 2/3 event *against* the LMS direction) rested on a misreading of Gooding et al. 1998: in TPM1, PTB, with MBNL1 (Gooding et al. 2013), represses exon 3, the *non*-smooth-muscle exon, so the published mechanism predicts PTBP1 *with* the LMS direction. The test and the verdict computed under the rule are kept as written ("supported"); the text reports R6 as running against the published mechanism. R4's rationale ("PTBP1 represses the smooth-muscle exons of ACTN1 and TPM1") holds for ACTN1 only | Found when the mechanism was re-read to write the result | Yes: the reading of R6 changed after its result was seen; no test, threshold or number changed |
| A-6 | 2026-10-12, after R1–R6 were read | Two exploratory checks of R5/R6 (`scripts/37 --only post`): the partial Spearman also adjusting for each sample's mean MXE PSI (the global axis of D-24), and PTBP1's rank among the 70 factors for each event | D-24 had shown that factor–PSI correlations need that covariate; D-46 mirrored D-18, which predates it | Yes: exploratory; reported beside the pre-specified test, never instead of it |
| A-7 | 2026-10-12, after the site probe's p-values were seen | D-47 did not name the FDR family for P1 and P2; it is taken as D-47's three tests (P1, P2 for GEX, P2 for the splicing median), the stricter of the two readings: P1 q = 0.080 ("not resolved") rather than 0.027 alone | The config's header rule needs a family | Yes: chosen after the p-values were seen, in the direction that weakens the result |
| A-9 | 2026-10-12, after the D-52 code was written, before it was run | K1 reads MRE11 exon 5 from the full SE matrix whatever the cohort pre-filter (`msi_cis` names the event, not the matrix subset), and `test` also writes exploratory tables (`share_explained`, `detector`, `shared_repeat_exons`, `dose_exon_delta`) that the config lists as exploratory. After D-53, an exploratory retained-intron check (`ri_check.csv`) was added to ask whether unspliced RNA could explain the increases | The config did not say which matrix subset a named event is read from | K1 rule: no. The retained-intron check: yes, it was added after the reversal at ≥ 20 nt was seen; it is exploratory and too thin to decide (4 long-run exons per cohort) |
| A-10 | 2026-10-12, S2: when the D-55 code was written, before it was run; S4: after it was run | S2 leaves out patients from tissue-source sites that sent only one patient of the cohort (UCEC 10, STAD 4, COAD 8): with a site fixed effect their leverage is 1, so their HC3 error is undefined, and they carry no within-site information. S4's score correlates ρ = 1.000 with the subtype-matrix score in all three cohorts (the same RNA files and nearly the same exons), so S4 checks the feature filter only, as D-55 says, and is not a second sample | `msi_robustness` did not say how S2 treats single-patient sites | S2: no. S4: a reading, added after the result; no number changed |
| A-11 | 2026-10-12, when the D-57 code was written, before any test was computed | D-57 uses autosomes only for L1–L6 (the GRCh37 anchors come from autosomal MC3 calls); a patient with two tumour samples contributes the first; A3 counts every genotype-like gene, anchored or not | A tumour is a clone, so one X is inactive in all its cells and an X heterozygote looks like LOH without any loss; `brca_germline` did not say how repeated samples or unanchored genes are treated | No: fixed before any PSI was compared with LOH, race or a label |
| A-12 | 2026-10-12, after H1–H4 were computed | Two exploratory additions to step 47: an in-silico knockout (in every stored BRCA_PAM50 fold model, the inputs from the 17q12 amplicon are set to their training mean and P(Her2) is recomputed for the held-out tumours, against 20 random input sets of the same size), and a check of the ERBB2 terminal events against the neighbouring genes in GENCODE v45 (MIEN1 and PGAP3 lie on the opposite strand over ERBB2's alternative last and first exons) | H2 came out against its prediction: the splicing models separate amplified luminal tumours better than GEX; the knockout asks whether they do it through amplicon events | Yes: exploratory, motivated by a result; no test or verdict changed |
| A-13 | 2026-10-12, after L1–L6 and A1–A3 were computed and the models' top genotype-like events were profiled by class | (1) A germline check of the genotype-like events, using tumour–normal pairs only: a homozygous normal cannot flip to the opposite homozygote in its own tumour, so an event is kept as confirmed germline when such flips occur in ≤ 5% of at least 20 homozygous pairs; L1–L6 and A3 are re-run on the confirmed events as a sensitivity analysis. (2) A2's test (crude minus race-stratified log odds ratio > 0) assumed the crude log odds ratio was positive, as the SHAP direction suggested; it is negative (GSTM1 present goes with ER-positive status), so a positive difference means the association grew, not shrank; A2 is read by its claim, as against it. (3) Bootstrap p-values are floored at 1/(B + 1) instead of reading 0 | The class profile showed allele frequencies that no inherited allele could have (TMEM101 high group: LumA 0.44, Basal 0.99); in pairs, the median flip rate of the 1,110 events is 0.20 and 82 pass the check | Yes: added after the results; the pre-specified verdicts are kept and reported with these readings |
| A-14 | 2026-10-12, after the A-13 germline check was computed | Exploratory, beside the A-13 flip rate: agreement between each normal's group and its own tumour's group (Cohen's kappa over complete pairs, at least 20). It is not a pre-specified filter, and L1–L6 and A3 are not re-run on it | The A-13 rule passed 82 events, but their minor-allele frequency in normals was about half that of all events (median 0.18 vs 0.33). With a rare allele, almost every pair is homozygous for the common one, so a flip can hardly happen, and 42 of the 82 had kappa ≤ 0.2. Kappa does not depend on allele frequency in the same way | Yes: chosen after the A-13 set was seen; it is reported as exploratory, and no verdict uses it |
| A-15 | 2026-10-13, after PC1–PC8 were computed | Exploratory, beside D-66: the AUROC of each arm of score G alone (R_down against B, R_up against B) in the 30 held-out cancers; the per-tumour arm scores are kept locally. No test, cut or verdict changed (the predictions were re-run and are identical) | PC2 supported score G and PC1 did not support score A; the arm split asks which part carries the held-out signal | Yes: added after the results; reported as exploratory |
| A-16 | 2026-10-13, after RV1–RV8 were computed | Exploratory checks of the rising arm (`scripts/52 --only posthoc`, `posthoc.csv`): among held-out exons with runs ≥ 20, the share more included by branch-point bin (18–29, 30–35, 36–49, 50–69, 70–101 nt); the shift against run length and exon length; shield exons tested against 0; R_up exons split by the sign of their PTB sensitivity; PTBP1, HNRNPC, U2AF2 and MKI67 expression in MMR-deficient tumours, and the RV2 model with them added. Also, the RV1 contrast table labels its non-MSI rows `other` (pandas reads `null` as missing). No test, cut or verdict changed | RV4 and RV5 did not support the geometry reading, and RV8 came out against its sign; these checks ask what the branch point does inside long runs, and whether the RV8 sign reflects trans regulator levels | Yes: chosen after the results; reported as exploratory |
| A-17 | 2026-10-13, after the D-71 features were built, before any model was fitted | The per-exon pool is computed over D-71's own tumours (the 8,464 with MC3 calls), as its tumour definition requires: 30,632 exons, 1,562 with an acceptor run ≥ 10. The pre-specification quoted 30,568 and 1,555, a count made beforehand over all 9,034 tumours | The pre-registration's count and its rule disagreed; the rule was followed | No: no model had been fitted and no label was used |
| A-18 | 2026-10-13, after the E1 AUROCs of M0–M6 were read, before the LR family was computed | Exploratory checks of what the per-exon models learn (`scripts/54 --only posthoc`, `posthoc.csv`, `posthoc_exons.csv`): M4 refitted on the E1 training set at its chosen C, standardised coefficients summarised by sequence class (runs 10–19; runs ≥ 20 with the branch point ≥ 36 nt; other runs ≥ 20) and the 20 largest; the share of M5's top-500 ANOVA exons with a run ≥ 10 against the pool's share, and their genes; M5 refitted on the E1 training set at its chosen K and C, its standardised coefficients by sequence class and the 20 largest; three ablations fitted on E1 with the grid and tuning of the model they modify, described with bootstrap intervals against M2 and no q: M5 on the pool exons without a run ≥ 10 (M5n); M4 on the exons of the 2 arms only (M4a, per-exon weights on the exons score G uses); and M4 on the repeat exons outside the arms (M4x). No test, cut or verdict changed | M4 and M5 beat the 10 bins held out, against LR4 and LR5's predicted sign; these checks ask whether the gain comes from the repeat sequence or from other exons | Yes: chosen after the results; reported as exploratory |
| A-19 | 2026-10-13, after D-71, D-72 and A-18 were computed and read | At the user's request, figures U23 and U24 and the write-ups show one fitted model: logistic regression on the repeat exons (M4), with score G as the unfitted reference in U23. The pre-registered tests are unchanged and still reported: the LR family compares M0–M6, and BM1–BM6 test the all-exon model (M5). M4's benchmark differences are exploratory (same paired bootstrap, no q). Step 55's discordance table is also computed for M4; nothing else is recomputed | Presentation: M4 keeps the sequence prior, so it is the fitted model that reads the code | Yes: chosen after the results, when M4 had the highest leave-one-cancer-out AUROC (0.967); its numbers are labelled exploratory where the tests were on M5 |
| A-20 | 2026-10-13, after MG1–MG5 and the exploratory tables of D-75 were computed and read | Four exploratory additions to `scripts/60` (`figure_data`): (1) shift summaries and figure panels restricted to exons whose proficient median PSI lies within 0.1–0.9; (2) mid-rank percentiles of MMR-deficient tumours among proficient ones, binned by class (ties count half; the held-out tail shares now use them too, where strict ranks had counted a tumour tied with every proficient value as a low-tail case); (3) MG4 and MG5's statistic by proficient-median band and on the 0.1–0.9 exons; (4) PSI histograms of MRE11 exon 5 and three exons chosen after the largest-movers table was read (METTL8 and BOD1L1, the largest mean shifts of each arm; CCNDBP1, the event of Figure U7) | The median exon shift was exactly 0 in every class (most exons sit at PSI 0 or 1), and MG4's location-shift expectation was found to move proficient values past 0 and 1: the positive control MRE11 exon 5 scored 57% observed against 100% expected in COAD | Yes: all four were added after the results were read; they are exploratory and reported beside the pre-specified tests, whose verdicts are unchanged |

## Verification that the existing tasks are unchanged
- `tests/test_subtype.py::TestExistingTasksUnchanged` rebuilds every existing task table from
  `samples.parquet` and compares it with the stored one, and checks that `subtype_config()`
  leaves `config/study.yaml`'s study block untouched.
- The 25 pre-existing tests pass unchanged.
- One existing fold (`tumor_status/ri/KICH`, outer fold 0) was re-run on 2026-10-12 with the
  changed `pipeline.run_fold`, into a scratch directory, and compared with its stored outputs:
  probabilities, logits, SHAP values, standardised inputs, expected values, row and feature
  order are bit-identical, and so are the chosen trial, inner loss, refit epochs, selected K and
  method, test AUROC, balanced accuracy, MCC and additivity error.
