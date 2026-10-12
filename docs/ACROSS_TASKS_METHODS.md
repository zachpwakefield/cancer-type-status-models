# Across the three tasks: methods and decisions

Tumour status, cancer type and subtype were modelled separately, on the same samples and the
same nine views. This layer asks how the three relate. It fits nothing: every number comes from
stored out-of-fold predictions, stored held-out SHAP values and the data. Decisions D-48..D-51
follow D-46 and D-47 in [SUBTYPE_METHODS.md](SUBTYPE_METHODS.md); the predictions, with their
sources, are in [config/followup.yaml](../config/followup.yaml) (`cross_task`), committed on
2026-10-12 before any of these analyses was run. That file also lists what had already been seen
when the predictions were written.

| Step | Script | What it does |
|---|---|---|
| 40 | `scripts/40_cross_task.py` | gene overlap, type and status probabilities by subtype, exon protein character |
| 41 | `scripts/41_figures_cross_task.py` | figures X1-X4 |

Results: [results/ACROSS_TASKS.md](../results/ACROSS_TASKS.md).

**Models compared.** Three out-of-fold scores per task: GEX, the AS ensemble (mean of the 8 AS
views' probabilities) and GEX + AS (0.5 GEX + 0.5 AS ensemble), as in D-23 and D-40. The
subtype patients are all in the cancer-type cohort (one primary tumour each); 12 of the 22
subtype cohorts belong to a cancer with a tumour-status model.

**Reading a prediction.** Two-sided tests. "Supported": the stated direction with BH q < 0.05
in its family; "against": the opposite direction with q < 0.05; otherwise "not resolved".

---

### D-48 Do the three tasks read the same genes? · `gene_overlap`
For each subtype cohort, its cancer's tumour-status model (if any) and the cancer-type model's
class for that cancer, the top 100 genes of each model per view (a gene ranks by its best
feature's mean |SHAP|; for cancer type, over the samples of the cancer's own class). Each pair is
compared within the genes both tasks could select (a feature passing both pre-filters):
overlap, its hypergeometric expectation and p, Jaccard; BH over every (pair, cohort, view).
G1: every pair overlaps above chance. G2: for histology cohorts, subtype overlaps type more
than it overlaps status (paired Wilcoxon on Jaccard). The genes in all three top-100 lists of a
cancer are listed.

### D-49 Which other cancers does a subtype resemble? · `type_probability_by_subtype`
The cancer-type model gives each tumour a probability for each of 34 projects. Per subtype class,
the mean of each project's probability, and the share of patients called another project.
Thirteen literature predictions (T1–T13) name a subtype, a set of target projects and a
direction: e.g. squamous oesophageal tumours carry more HNSC/LUSC/CESC probability than
oesophageal adenocarcinomas; IDH-wildtype lower-grade gliomas more GBM probability; basal-like
breast tumours more ovarian. Measure: the patient's summed probability over the targets;
Mann–Whitney between the named groups, effect as AUROC; BH over T1–T13 within each model.
T11 (cortical admixture in PCPG → ACC) is a positive control.

### D-50 Does the tumour-status model see some subtypes as less tumour-like? · `status_probability_by_subtype`
For the 12 subtype cohorts with a tumour-status model, the out-of-fold P(tumour) of each primary
tumour, by class: Kruskal–Wallis per cohort (descriptive, BH over cohorts) and six predictions
(S1–S6): normal-like breast, follicular-variant thyroid, diffuse and genomically stable gastric
and lobular breast tumours look less tumour-like; tall-cell thyroid more. C1: within each such
cancer's tumours, P(tumour) falls as the stromal and immune scores rise (Spearman), which would
make composition the common reason.

### D-51 Protein character of the exons each task uses · `exon_protein_character`
D-17 compared the top cassette exons of tumour status and cancer type with 6,000 random measured
exons. The same metrics for the union of every subtype cohort's top 100 SE events, against the
same background; Fisher exact tests, BH over metrics and sets. E1: subtype exons share the
other tasks' bias (frame-preserving, outside domains, in disordered regions).

### Amendments and deviations
Numbered with the subtype amendments (A-1..A-7 in [SUBTYPE_METHODS.md](SUBTYPE_METHODS.md)).

| # | When | What | Why | Informed by results? |
|---|---|---|---|---|
| A-8 | 2026-10-12, written into `scripts/40` before it was first run, but not into the config | D-49 sensitivity: the targets' share of the probability *not* given to the patient's own project, P(targets) / (1 − P(own)), same test | A summed target probability also rises when a tumour is simply less confidently its own project; the share separates "resembles the target" from "resembles nothing in particular" | No, but it is not pre-specified and is reported as a sensitivity |

### What is published (D-12)
Aggregates only: per class, cohort and view. Per-patient probabilities stay in
`<data_dir>`.
