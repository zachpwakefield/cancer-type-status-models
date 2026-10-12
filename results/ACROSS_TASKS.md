# Across the three tasks: tumour status, cancer type and subtype

The three tasks were modelled separately, on the same samples and the same nine views: is a
sample tumour ([TUMOUR_STATUS.md](TUMOUR_STATUS.md)), which tissue it came from
([CANCER_TYPE.md](CANCER_TYPE.md)) and, within a cancer, which subtype it is
([SUBTYPE.md](SUBTYPE.md)). This section asks how their results relate. Nothing here is refitted:
every number comes from the stored out-of-fold predictions, held-out SHAP values and the data.

- **Methods and decisions:** [docs/ACROSS_TASKS_METHODS.md](../docs/ACROSS_TASKS_METHODS.md)
  (D-48 to D-51, amendment A-8). The 13 + 6 + 4 predictions, their sources and the rule for
  "supported" were committed in [config/followup.yaml](../config/followup.yaml) before anything was
  computed, together with a list of what had already been seen.
- **Scripts:** `40` computes the tables, `41` draws the figures.
- **Figures:** X1–X4 in [figures/across_tasks/](figures/across_tasks/), each with a caption file and
  a CSV of the plotted values.
- **Tables:** `results/across_tasks/`.
- **Models compared:** GEX, the AS ensemble (mean of the 8 splicing views) and GEX + AS (1:1).
- **Supported** = the predicted direction with BH q < 0.05 in its family; **against** = the
  opposite direction with q < 0.05; otherwise **not resolved**.

## The argument

The three tasks look like three questions, but their models share much of what they read, and
their uncertainty lines up with biology the literature already describes.

**What the models read.** Any two tasks' top-100 gene lists overlap 5–8 times more than chance
(X1). The genes in all three lists of a cancer are mostly splicing switches of smooth-muscle and
stromal genes — TPM1's mutually exclusive exons in 10 cohort-views, CHN1, ACTN1 — rather than
expression markers. That points to tissue composition as common ground.

**Which other cancers a subtype resembles** (X2). The pan-cancer type model was never told about
subtypes, yet where it spreads probability onto other projects it does so as the literature
predicts:
- squamous oesophageal and cervical tumours lean towards head-and-neck and lung squamous
  carcinoma, and oesophageal adenocarcinoma towards stomach (TCGA 2017);
- cervical adenocarcinoma leans towards endometrium (TCGA 2017);
- IDH-wildtype lower-grade glioma towards glioblastoma (Brat et al. 2015);
- copy-number-high and serous endometrial tumours towards ovary and uterine carcinosarcoma (TCGA
  2013);
- biphasic mesothelioma towards sarcoma;
- HPV-positive head-and-neck tumours towards cervix (Pyeon et al. 2007).

Twelve of thirteen predictions hold in GEX + AS. The thirteenth, intestinal-type gastric cancer
leaning towards colon, is not resolved. One supported result is weaker than it looks: basal-like
breast tumours carry more ovarian probability, but only because they carry more off-target
probability of every kind.

**How tumour-like the status model finds each subtype** (X3). All six subtype contrasts predicted
in advance appear in the splicing models; GEX resolves two of them:
- normal-like breast cancer looks less tumour-like than the other PAM50 classes;
- in papillary thyroid cancer, follicular variant less than classical, and tall cell more;
- diffuse and genomically stable gastric cancer look less tumour-like;
- lobular breast cancer looks less tumour-like than ductal.

GEX's probabilities sit at the ceiling, so this says as much about calibration as about biology.
Stromal or immune content does not explain it: the correlation runs both ways across cancers.

**What the exons do to the protein** (X4). The subtype models' cassette exons share the other two
tasks' profile: more often frame-preserving, less often in folded domains, more often in
disordered segments.

## Summary

1. **Shared genes** (X1, G1 supported): median 7.8-fold the chance overlap between tumour status
   and cancer type, 5.0-fold between subtype and status, 6.1-fold between subtype and type; 87–94%
   of comparisons at q < 0.05. **G2 not supported**: histology subtypes are not closer to the type
   model than to the status model (median Jaccard difference 0.000).
2. **Type probabilities follow the literature** (X2): T1–T12 supported in GEX + AS (AUROC
   0.65–0.95), T13 not resolved. AS ensemble 13 of 13, GEX 11 of 13. Off-target-share sensitivity
   (A-8): T6 reverses.
3. **Status probabilities follow differentiation, in splicing** (X3): S1–S6 supported in the AS
   ensemble and in GEX + AS (AUROC 0.57–0.84), 2 of 6 in GEX. C1 (composition explains it) not
   supported.
4. **Exon character** (X4, E1 supported): subtype top exons divisible by 3 49% against 41%, in a
   domain 39% against 51%, in a disordered region 18% against 14%.

---

## 1. The three tasks read overlapping genes

![X1](figures/across_tasks/X1_gene_overlap.png)

For every subtype cohort (22), its cancer's tumour-status model (in 8 cancers, covering 12
cohorts) and the type model's importance over that cancer's own samples, the top 100 genes of
each model per view were compared within the genes both tasks could select (D-48).

- **Every pair overlaps far above chance** (result). Status–type 7.8-fold (72 cohort-views,
  94% q < 0.05), subtype–status 5.0-fold (117, 87%), subtype–type 6.1-fold (198, 88%). Six of the
  387 comparisons share no gene (all subtype comparisons). The fold enrichment is largest in HIT
  (18–32-fold), whose 170,000 measured events make chance overlap small, and smallest in MXE and
  RI (3–5-fold).
- **No lineage preference** (result, G2 not supported). For the five histology cohorts with a status
  model, the subtype list overlaps the status list as much as the type list (median Jaccard 0.053
  for both; paired difference 0.000, p = 0.24). The prediction was that histology, being lineage,
  would sit closer to the type model.
- **What all three share** (result, descriptive; `genes_in_all_three.csv`). 413 gene × cohort-view
  entries are in all three top-100 lists, 402 of them from splicing views. The most frequent are:
  - TPM1 (10 cohort-views, mostly its exon 2/3 MXE);
  - SNHG29 (9, A3SS and RI);
  - CHN1 (5, MXE);
  - DNM2 (4) and GCOM1 (4).

  TPM1 and ACTN1 carry the smooth-muscle exons of §9 of SUBTYPE.md. GCOM1's MXE is the MYZAP exon
  that the D-27 audit of TUMOUR_STATUS.md found inside POLR2M readthrough transcripts.

**Reading.** The shared core is not a set of expression markers: GEX contributes 11 of the 413
entries. It is a small set of splicing switches in contractile and stromal genes, which move with
how much smooth muscle and stroma a sample contains. That is a property of every task's samples
— tumours differ from normal tissue, tissues differ from each other, and subtypes differ in
stroma — so each task's models pick it up (interpretation).

## 2. Which other cancers a subtype resembles

![X2](figures/across_tasks/X2_type_probability_by_subtype.png)

The cancer-type model assigns each tumour probabilities over 34 classes (33 projects and normal
tissue). The pre-specified test is whether a subtype carries more probability for the projects the
literature says it resembles (D-49).

| | Prediction (source) | Mean target probability, GEX + AS | AUROC GEX / AS / GEX + AS |
|---|---|---|---|
| T1 | ESCA squamous → HNSC + LUSC + CESC (TCGA 2017) | 0.009 vs 0.001 | 0.79 / 0.78 / 0.81 |
| T2 | ESCA adenocarcinoma → STAD (TCGA 2017) | 0.287 vs 0.025 | 0.90 / 0.93 / 0.95 |
| T3 | CESC squamous → HNSC + LUSC (TCGA 2017) | 0.059 vs 0.009 | 0.68 / 0.85 / 0.85 |
| T4 | CESC adenocarcinoma → UCEC (TCGA 2017) | 0.104 vs 0.019 | 0.73 / 0.86 / 0.86 |
| T5 | LGG IDH-wildtype → GBM (Brat 2015; Ceccarelli 2016) | 0.068 vs 0.004 | 0.85 / 0.88 / 0.87 |
| T6 | BRCA basal-like → OV (TCGA 2012) | 0.0007 vs 0.0003 | 0.54 / 0.78 / 0.72 |
| T7 | UCEC CN-high → OV + UCS (TCGA 2013) | 0.047 vs 0.025 | 0.73 / 0.74 / 0.75 |
| T8 | UCEC serous → OV + UCS (TCGA 2013) | 0.056 vs 0.024 | 0.71 / 0.71 / 0.73 |
| T9 | MESO biphasic → SARC (Hmeljak 2018) | 0.105 vs 0.020 | 0.70 / 0.80 / 0.80 |
| T10 | HNSC HPV-positive → CESC (Pyeon 2007) | 0.242 vs 0.021 | 0.90 / 0.92 / 0.92 |
| T11 | PCPG cortical admixture → ACC (control) | 0.039 vs 0.0004 | 0.73 / 0.80 / 0.78 |
| T12 | STAD CIN → ESCA (TCGA 2017) | 0.114 vs 0.044 | 0.64 / 0.64 / 0.65 |
| T13 | STAD intestinal → COAD + READ (Laurén 1965) | 0.005 vs 0.002 | 0.54 / 0.62 / 0.57 |

All AUROCs are q < 0.05 except T6 in GEX and T13 in GEX and GEX + AS (CSV). In words:
- **Squamous and glandular carcinomas lean across organs by histology.** An oesophageal
  adenocarcinoma's second choice is the stomach, not oesophageal squamous carcinoma: 18 of 81 are
  called STAD outright, against 0 of 80 squamous tumours. Cervical adenocarcinomas lean towards
  endometrium (8 of 45 called another project, 5 of them UCEC). The TCGA oesophageal and cervical
  papers drew the same lines from multi-omic clustering (context).
- **Molecular classes lean the same way.** IDH-wildtype lower-grade gliomas lean to GBM (3 of 94
  called another project, 2 of them GBM). Copy-number-high endometrial tumours lean to ovary and
  carcinosarcoma. HPV-positive head-and-neck tumours lean to cervix (3 of 33 called CESC).
- **The positive control works:** PCPG tumours with adrenal-cortex admixture carry on average
  about 100 times the ACC probability of the rest.

**Two cautions** (sensitivity A-8, not pre-specified). A higher summed probability can also mean
that a tumour is less confidently its own project, so the test was repeated on the targets'
*share* of the off-target probability.
- **T6 reverses.** Basal-like breast tumours give ovary a smaller share of their off-target
  probability than other breast tumours do (AUROC 0.29). They carry more off-target probability of
  every kind; 4 of the 6 called another project are called SARC.
- **T5 weakens in the AS ensemble** (0.57).

The other eleven hold on the share.

**Not pre-specified, descriptive.** No MSI colon tumour is called READ (0 of 40), against 26 of
162 CIN tumours. MSI colon cancers sit mostly in the proximal colon, far from the rectum (TCGA
2012; context).

## 3. How tumour-like the status model finds each subtype

![X3](figures/across_tasks/X3_status_probability_by_subtype.png)

The tumour-status model of each cancer scores every tumour out of fold. Six predictions name a
subtype expected to look less tumour-like, because it is closer to the normal tissue or carries
more of it (D-50).

- **Normal-like breast cancer** looks least tumour-like (AUROC 0.84 in the AS ensemble, 0.69 in
  GEX). In each model 13% of normal-like tumours are scored below 0.5, against ≤ 0.5% of any
  other PAM50 class. That fits the view that this class reflects normal-tissue content (Parker et
  al. 2009; context).
- **Thyroid differentiation orders the histologies** in the splicing models:
  - follicular variant less tumour-like than classical (AUROC 0.67 in the AS ensemble; GEX 0.49,
    not resolved);
  - tall cell more tumour-like than classical (0.62; GEX 0.65).

  RAS-like follicular tumours are the most differentiated papillary carcinomas and the tall-cell
  cluster the least (TCGA 2014; context).
- **Gastric and breast histology:**
  - diffuse less tumour-like than intestinal (0.61 AS ensemble; GEX not resolved);
  - genomically stable less than the other molecular classes (0.70; GEX not resolved);
  - lobular less than ductal (0.58; GEX 0.55, q = 0.060).

**Why GEX resolves only two** (result). GEX's P(tumour) for these primary tumours sits at the
ceiling: the comparison groups' medians are 0.998–1.000. Ranks among them carry little. The AS
ensemble averages eight models and spreads out (medians 0.95–0.999), so it can order the
subtypes. This is partly calibration, not only biology: the splicing models are less certain
overall, and their uncertainty concentrates on the subtypes closest to normal tissue.

**Not composition** (result, C1 not supported). If stromal and immune content explained it,
P(tumour) would fall as those scores rise in every cancer. Over the 8 cancers' tumours in GEX + AS:
- falls with the scores in BRCA, COAD and PRAD (4 of 16 correlations negative at q < 0.05);
- rises with them in HNSC and THCA (3 positive; THCA stromal ρ +0.25).

So the status models' doubt about these subtypes is not simply their stroma.

## 4. The exons each task uses do the same kind of thing to the protein

![X4](figures/across_tasks/X4_exon_protein_character.png)

The union of the subtype models' top-100 cassette exons (1,857), against the 6,000 random measured
exons of D-17, alongside the tumour-status (1,281) and cancer-type (1,914) sets of S6 (D-51):

| | Background | Status | Type | Subtype (q) |
|---|---|---|---|---|
| length divisible by 3 | 41% | 51% | 54% | 49% (1.6e-9) |
| overlaps a Pfam/CDD domain (fully coding) | 51% | 37% | 35% | 39% (1.3e-10) |
| overlaps a disordered region (fully coding) | 14% | 21% | 20% | 18% (0.005) |
| overlaps a TM helix (fully coding) | 6.7% | 3.7% | 4.6% | 4.1% (0.005) |

The subtype exons carry the same signature, slightly weaker on every metric (result, E1
supported). Exons that preserve the frame and sit in disordered, interaction-prone segments are
the kind tissues regulate (Buljan et al. 2012; Ellis et al. 2012; context). Whatever separates a
tumour from normal tissue, one tissue from another and one subtype from another, the splicing
models reach for the same class of exon.

## What would change the reading

- **Calibrated status probabilities.** If recalibrating GEX's tumour-status scores (e.g. on the
  logit scale) recovered the S2–S6 orderings, §3's contrast between GEX and splicing would be
  calibration only. If it did not, the splicing models carry differentiation information that
  the expression model does not.
- **Purity estimates.** C1 used marker-gene scores. ABSOLUTE or ESTIMATE purity (not local) could
  still show that the less tumour-like subtypes are simply less pure.
- **A different type model.** X2 reads one pan-cancer model per view. A model trained without the
  target projects, or a different classifier, would test whether the "resemblances" belong to the
  tissues or to this model's decision boundaries.

## Limitations

- **One resource, out-of-fold only.** Every comparison is within TCGA. The type and status models
  saw the same patients' other samples in training folds, though never the scored sample itself.
- **Predictions were written after the main results.** The pre-specification lists what had been
  seen, including project-level type confusions (ESCA ↔ STAD, CESC ↔ HNSC). T2 and T3 were
  therefore expected at the project level, and only their subtype split was untested.
- **A-8 is not pre-specified.** It is reported as a sensitivity, and it changes the reading of T6.
- **Gene identity is coarse.** X1 counts a gene once, whatever event or view carried it. The
  overlap is between genes, not between the same exon.

## References

Literature used as *context* (never to choose a feature, a model or an event):

- Brat DJ et al. (2015) Comprehensive, integrative genomic analysis of diffuse lower-grade gliomas. *N Engl J Med* 372:2481.
- Buljan M et al. (2012) Tissue-specific splicing of disordered segments that embed binding motifs rewires protein interaction networks. *Mol Cell* 46:871.
- Cancer Genome Atlas Network (2012) Comprehensive molecular characterization of human colon and rectal cancer. *Nature* 487:330.
- Cancer Genome Atlas Network (2012) Comprehensive molecular portraits of human breast tumours. *Nature* 490:61.
- Cancer Genome Atlas Research Network (2013) Integrated genomic characterization of endometrial carcinoma. *Nature* 497:67.
- Cancer Genome Atlas Research Network (2014) Integrated genomic characterization of papillary thyroid carcinoma. *Cell* 159:676.
- Cancer Genome Atlas Research Network (2014) Comprehensive molecular characterization of gastric adenocarcinoma. *Nature* 513:202.
- Cancer Genome Atlas Research Network (2017) Integrated genomic characterization of oesophageal carcinoma. *Nature* 541:169.
- Cancer Genome Atlas Research Network (2017) Integrated genomic and molecular characterization of cervical cancer. *Nature* 543:378.
- Ceccarelli M et al. (2016) Molecular profiling reveals biologically discrete subsets and pathways of progression in diffuse glioma. *Cell* 164:550.
- Ellis JD et al. (2012) Tissue-specific alternative splicing remodels protein-protein interaction networks. *Mol Cell* 46:884.
- Fishbein L et al. (2017) Comprehensive molecular characterization of pheochromocytoma and paraganglioma. *Cancer Cell* 31:181.
- Hmeljak J et al. (2018) Integrative molecular characterization of malignant pleural mesothelioma. *Cancer Discov* 8:1548.
- Laurén P (1965) The two histological main types of gastric carcinoma: diffuse and so-called intestinal-type carcinoma. *Acta Pathol Microbiol Scand* 64:31.
- Parker JS et al. (2009) Supervised risk predictor of breast cancer based on intrinsic subtypes. *J Clin Oncol* 27:1160.
- Pyeon D et al. (2007) Fundamental differences in cell cycle deregulation in human papillomavirus-positive and human papillomavirus-negative head/neck and cervical cancers. *Cancer Res* 67:4605.

Pyeon 2007, the cervical (2017) and thyroid (2014) TCGA papers were checked against PubMed, PMC or
journal records on 2026-10-12. The thyroid paper's full text was read for which histology is most
differentiated.
