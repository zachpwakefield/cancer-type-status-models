# Cancer type, tumour status and subtype from TCGA RNA processing

These are deep-learning classifiers that predict, from a single TCGA RNA-seq sample:

1. **Tumour status, per cancer.** Is the sample tumour or adjacent normal tissue? This is a
   separate binary model for each of 15 cancer types with ≥ 10 normal samples.
2. **Cancer type, pan-cancer.** Which of the 33 TCGA cancer types is the sample from, or is it
   normal tissue? This is a 34-class model.
3. **Subtype, within a cancer.** Which histological, molecular or IHC-defined subtype is a
   primary tumour? This is a separate model for each of 22 subtype cohorts in 15 cancers.

A fourth layer fits no new model. It compares the three tasks' predictions and explanations on
the same patients.

Each task is modelled **independently in each of nine modalities**:

| View | Feature | Source |
|---|---|---|
| AFE | alternative first-exon usage | HITindex |
| ALE | alternative last-exon usage | HITindex |
| HIT | HIT index (first ↔ last exon behaviour, −1…1) | HITindex |
| SE | skipped-exon PSI | rMATS |
| A5SS | alternative 5′ splice-site PSI | rMATS |
| A3SS | alternative 3′ splice-site PSI | rMATS |
| MXE | mutually-exclusive-exon PSI | rMATS |
| RI | retained-intron PSI | rMATS |
| GEX | gene expression, log2(CPM + 1) | STAR counts |

Every model is a hyper-parameter-tuned multilayer perceptron, evaluated by nested,
patient-grouped cross-validation. It is explained with SHAP on held-out samples only.

## Results

Three questions, each answered on its own and then compared. Every task section below has the
same four parts: what was asked, performance, what the models learned (figure by figure), and the
checks that bound the claims. Each also has a full write-up. A fourth section asks how the three
tasks relate.

| | [Tumour status](#tumour-status-tumour-vs-normal-tissue-within-a-cancer) | [Cancer type](#cancer-type-34-classes-pan-cancer) | [Subtype](#subtype-within-a-cancer-22-cohorts-in-15-cancers) |
|---|---|---|---|
| Question | tumour or adjacent normal tissue? | which of 33 projects, or normal tissue? | which subtype of this cancer? |
| Scope | 15 cancers, 6,815 samples, 135 models | 34 classes, 9,729 samples, 9 models | 22 cohorts in 15 cancers, 7,552 patient–cohort pairs, 198 models |
| GEX | mean AUROC 0.992; 99 errors | accuracy 0.954 | balanced accuracy 0.84 / 0.88 / 0.90 (histology / genomic or IHC / expression-defined) |
| AS ensemble (8 splicing views) | 0.994; 91 errors | 0.950 | 0.80 / 0.82 / 0.81 |
| GEX + AS (1:1) | 77 errors (McNemar p = 1.1e-4) | 0.961 (p = 3.9e-9) | −0.04 to +0.04 per cohort; no McNemar p < 0.05 |
| Headline | splicing is a second, partly independent view of the tumour state | splicing is a lineage code, and fails where a lineage code should | splicing matches expression on lineage and genome-state subtypes and loses where one gene defines the class |
| Write-up · figures | [TUMOUR_STATUS.md](results/TUMOUR_STATUS.md) · S1–S10 | [CANCER_TYPE.md](results/CANCER_TYPE.md) · T1–T9 | [SUBTYPE.md](results/SUBTYPE.md) · U1–U26 |

[**Across the three tasks**](#across-the-three-tasks) ([ACROSS_TASKS.md](results/ACROSS_TASKS.md),
X1–X4) compares them:
- the tasks' models share genes 5–8 times more often than chance;
- 12 of 13 literature predictions about which other cancers a subtype resembles hold in the
  cancer-type model's probabilities;
- the splicing status models find all six predicted subtypes less tumour-like;
- the exons all three tasks use have the same protein character.

One subtype result goes furthest. The MSI splicing marks are largely written in the genome: the
exons the MSI models chose sit after long poly-T tracts that mismatch-repair deficiency destabilises.
Tract length and the branch point's distance from the splice site say which exons lose and which
gain inclusion. That two-part code, read in RNA without re-fitting, finds mismatch-repair-deficient
tumours in 30 cancers never used to define it (AUROC 0.83)
([U10–U11, U16–U26](#8-msi-splicing-follows-acceptor-poly-t-repeats-in-the-genome)). The MSI line of work
is written up as a journal-style draft in [manuscript/MSI_SPLICING_CODE.md](manuscript/MSI_SPLICING_CODE.md).
In breast cancer, splicing reads HER2 through a programme rather than through the amplified genes
([U12–U15](#9-breast-cancer-with-dna-two-her2s-splice-site-mutations-and-what-the-pair-check-removed)).

Full tables, CIs and the top SHAP features of every cohort and view are in
[results/RESULTS.md](results/RESULTS.md). Methods are in three files:
- [docs/BIOLOGY_METHODS.md](docs/BIOLOGY_METHODS.md), D-13 to D-29: interpretation of status and
  type;
- [docs/SUBTYPE_METHODS.md](docs/SUBTYPE_METHODS.md), D-30 to D-47 and D-52 to D-74: the subtype task
  and its MSI-repeat and BRCA follow-ups;
- [docs/ACROSS_TASKS_METHODS.md](docs/ACROSS_TASKS_METHODS.md), D-48 to D-51: the cross-task
  layer.

## Tumour status: tumour vs normal tissue within a cancer

**Full write-up: [results/TUMOUR_STATUS.md](results/TUMOUR_STATUS.md)** (figures S1–S10) — *splicing
is a second, partly independent view of the tumour state.* One binary model per cancer and per
view, for the 15 TCGA projects with at least 10 normal samples (6,815 samples, 135 models).

### Performance

<!-- RESULTS-STATUS:START (generated by scripts/06_summarise_results.py) -->
**Out-of-fold, mean over the 15 per-cancer models:**

| View | AUROC | Balanced accuracy | MCC |
|---|---|---|---|
| AFE | 0.985 | 0.956 | 0.859 |
| ALE | 0.973 | 0.929 | 0.806 |
| HIT | 0.976 | 0.954 | 0.860 |
| SE | 0.973 | 0.946 | 0.866 |
| A5SS | 0.982 | 0.939 | 0.850 |
| A3SS | 0.979 | 0.938 | 0.835 |
| MXE | 0.982 | 0.938 | 0.831 |
| RI | 0.983 | 0.952 | 0.882 |
| GEX | 0.992 | 0.977 | 0.906 |

<!-- RESULTS-STATUS:END -->

<img src="results/summary/fig_tumor_status_auroc.png" width="511" alt="Out-of-fold AUROC of every tumour-vs-normal model, 15 cancers by 9 views">

**Figure S0. Every view separates tumour from normal tissue almost perfectly in most cancers.**
Pooled out-of-fold AUROC of each per-cancer tumour-vs-normal model: rows are the 15 cancers with ≥
10 normal samples, columns the nine views. Each sample is scored by the outer-fold model that never
saw its patient.

The weak spot is ESCA, the cohort with the fewest normal samples (11): every splicing view drops
there (SE to 0.711) while expression holds at 0.976. PRAD is the next hardest. CIs and balanced
accuracy are in [results/RESULTS.md](results/RESULTS.md).

[caption](results/summary/fig_tumor_status_auroc.caption.md)

### What the models learned, in figures

Six of the nine figures, in the order the argument runs. Each has a caption file and a CSV of the
plotted values next to it.

#### 1. Expression and splicing fail on different samples

<img src="results/figures/tumour_status/S1_modality_gain.png" width="866" alt="S1: out-of-fold errors per cancer and pooled, for GEX, the AS ensemble and their fusion">

**Figure S1. Expression and splicing err on different samples, so their fusion beats expression
alone.** (a) Out-of-fold errors as a share of each cancer's samples (n in brackets) for expression
(black), the AS ensemble (the mean of the 8 splicing views, light blue) and their equal-weight
fusion (dark blue); stars = exact McNemar, fusion against expression. (b) Errors pooled over all
6,815 samples; 'stacked' is a cross-fitted logistic regression that learns the weights, shown as a
ceiling.

Expression makes 99 errors and splicing 91, but not the same ones, so the fusion makes 77: 27 of
expression's errors corrected against 5 introduced (exact McNemar p = 1.1e-4). The weights were
fixed before any result was seen. The stacker, allowed to learn them, does no better (76).

[caption](results/figures/tumour_status/S1_modality_gain.caption.md) ·
[§1](results/TUMOUR_STATUS.md#1-combining-the-modalities-beats-expression-alone)

#### 2. A transcriptome-wide shift that no model selection could have biased

<img src="results/figures/tumour_status/S2_global_shifts.png" width="467" alt="S2: median paired tumour minus normal shift, 15 cancers by 7 event views">

**Figure S2. Tumours retain more intronic sequence and move their last exons upstream, relative to
their own normal tissue.** Every measured event, with no SHAP ranking anywhere in the calculation:
the median over matched patients of the patient's own tumour − normal change, by cancer (rows, pairs
in brackets) and event class (columns). Blue = higher in the tumour, orange = higher in the normal
tissue (the house tumour/normal colours); stars = raw Wilcoxon p across patients.

Tumours **retain more intronic sequence** (7 cancers up, BRCA down), **use more long splice-site
forms** (A5SS up in 10, down in none) and **move their last exons upstream** (7 cancers, none the
other way). Caveat: RNA integrity also raises intronic reads, and TCGA has no RIN for these samples.

[caption](results/figures/tumour_status/S2_global_shifts.caption.md) ·
[§2](results/TUMOUR_STATUS.md#2-tumours-retain-more-intronic-sequence-and-choose-upstream-last-exons)

#### 3. The splicing signal organises by regulator, not by pathway

<img src="results/figures/tumour_status/S4_regulons.png" width="900" alt="S4: regulon enrichment of model events, and regulon overlap">

**Figure S4. The models' events sit in small, specific splicing-factor programmes, not in the large
shared one.** A factor's *core regulon* is the events whose PSI tracks its expression within
tumours, reproducibly in ≥ 5 of the 15 cancers, after residualising out stromal and immune scores
and the sample's own mean PSI. (a) The odds that a model's top event sits in a factor's regulon,
against a background matched on coverage and spread; only q < 0.05 shown (BH over 32 factors),
regulon size in brackets. (b) Each factor's regulon size against its largest overlap with another
factor's regulon (Jaccard), coloured by the odds of (a).

Without the mean-PSI covariate every nuclear factor appears to regulate a thousand events, because
factor expression and overall PSI share one sample-level axis. The large regulons are not separate
programmes: 1,074 events carry 7,529 memberships, about 7 factors each. The models' events
concentrate in the small, distinct programmes of PTBP1, RBM47, HNRNPU, HNRNPF and QKI and are
depleted from the big shared one.

[caption](results/figures/tumour_status/S4_regulons.caption.md) ·
[§4](results/TUMOUR_STATUS.md#4-the-models-events-live-in-small-specific-splicing-factor-programmes)

The complement of this result is that the splicing models' genes form almost no pathways: their
top genes reach any enriched Hallmark, Reactome or GO-BP term in at most 9 of the 15 cancers, and
in none for retained introns, against 14 of 15 for expression
([S3](results/figures/tumour_status/S3_pathways.png),
[§3](results/TUMOUR_STATUS.md#3-expression-reads-the-microenvironment-splicing-does-not-form-pathways)).
The splicing signal is not the expression signal re-encoded, which is the mechanical reason the
fusion in S1 works.

#### 4. Different cancers rely on the same events

<img src="results/figures/tumour_status/S5_recurrent_events.png" width="820" alt="S5: the 40 most recurrent events, their paired change in every cancer">

**Figure S5. Different cancers rely on the same events, and those events move the same way across
cancers.** The 40 events that rank in the top 100 of the most per-cancer models. Colour = the median
paired tumour − normal change in **every** cancer, including those where the event did not rank, so
the pattern across a row is not selected on (blue higher in the tumour, orange higher in the normal
tissue); dots = the cancers whose top 100 contain the event.

404 events rank in the top 100 of at least 4 of the 15 independent per-cancer models; about 3 would
be expected if the rankings were independent. Recurrent events tend to move the same way everywhere
(median 73% of cancers), and several genes contribute switch *pairs* moving in opposite directions:
two EXOSC7 last exons, two NDE1 first exons, two TPM1 mutually exclusive exons. Two of those pairs
are partly another gene (Figure S10): the falling EXOSC7 exon follows the stromal gene CLEC3B, and
the falling NDE1 exon follows smooth-muscle MYH11.

[caption](results/figures/tumour_status/S5_recurrent_events.caption.md) ·
[§5](results/TUMOUR_STATUS.md#5-different-cancers-rely-on-the-same-events)

#### 5. What those events do to the protein

<img src="results/figures/tumour_status/S9_event_protein.png" width="882" alt="S9: the alternative segment of nine recurrent events on their proteins">

**Figure S9. The most shared events exchange protein segments through mutually exclusive and
terminal exons.** Nine of the ten most shared events, each on its own protein scaled N to C: grey
blocks = Pfam/CDD domains, the line underneath = disordered regions, coloured block = the segment
the event adds or removes (blue higher in tumours, orange lower). The tenth belongs to no coding
transcript and cannot be drawn.

Not one of the ten is a cassette exon: six are mutually exclusive pairs and four are terminal exons,
including the macroH2A1.1/1.2 switch, the chimaerin SH2 switch and a length-neutral 23-amino-acid
exchange in the collagen V C-terminal domain. Three transcripts carry no InterPro features at all
and say so, because an empty annotation is not evidence of no domain. The macroH2A1 exon is the one
in macroH2A1.2, and it behaves as the literature predicts: it falls as QKI rises within tumours (12
of 15 cancers; QKI promotes macroH2A1.1, Novikov et al. 2011) and rises with proliferation (13 of
15).

[caption](results/figures/tumour_status/S9_event_protein.caption.md) ·
[§9](results/TUMOUR_STATUS.md#9-the-events-one-by-one)

#### 6. Does any of it predict outcome?

<img src="results/figures/tumour_status/S8_survival.png" width="833" alt="S8: share of event-by-cancer pairs with Cox p below 0.05, and the strongest single events">

**Figure S8. The models' events are modestly enriched for survival associations, more so the more
models agree.** Every model event was tested against a background event matched on view, coverage
and spread, through the same tool with the same gates; nothing about outcome entered either
selection. (a) The share of event × cancer pairs reaching a nominal Cox p < 0.05 for two
pre-specified arms (PSI + host gene; + age, sex and stage), with 95% intervals; dashed line = the
same pipeline on outcomes shuffled between patients. (b) The twelve strongest single events: hazard
ratio per SD with 95% CI (cancers with p < 0.05 / cancers tested).

Model events beat the background (11.1% against 9.0%), and the rate **rises with how many models
agree** on an event (odds 1.27 → 1.56 → 1.62), which is what one expects if the models select
biology rather than cohort noise. Shuffled outcomes give 5.6%, so the 9% background is real signal
too.

The careful reading is in the narrative: the effect is modest, it halves once age, sex and stage
enter, the events that clear FDR *inside* a cancer are 94% renal (though the set-level gradient
holds without the kidney cohorts), and directions flip between cancers. Separating tumour from
normal tissue and predicting outcome are distinct axes.

[caption](results/figures/tumour_status/S8_survival.caption.md) ·
[§8](results/TUMOUR_STATUS.md#8-do-these-events-carry-information-about-outcome)

#### 7. Which events measure another gene

<img src="results/figures/tumour_status/S10_overlap_audit.png" width="863" alt="S10: share of events overlapping another gene's exon, all measured against the models' top events, and how many track that gene">

**Figure S10. Some of the models' events, terminal exons above all, read a neighbouring gene's
expression.** Most TCGA libraries are unstranded and GENCODE nests genes inside one another, so an
event named after one gene can count another's reads. (a) The share of events whose segment overlaps
another gene's exon, for every measured event (grey) and the models' top events (dark), by view. (b)
Of those overlapping top events, how many follow the other gene's expression more than their own (a
rule fixed in advance).

Terminal (HIT) exons stand out: 44% against 7%, odds 6.6 with exon length held fixed. 163 of 1,105
overlapping top events track the other gene. Three of the ten events of Figure S9 do: EXOSC7 →
CLEC3B, ENOSF1 → TYMS, to which it was first described as an antisense RNA, and PLCG1 → ZHX3. The
models do not change, and the survival enrichment of Figure S8 is the same without any overlapping
event (odds 1.25 against 1.27); what changes is how these events are named.

[caption](results/figures/tumour_status/S10_overlap_audit.caption.md) ·
[§10](results/TUMOUR_STATUS.md#10-which-events-measure-another-gene)

Not shown here, both in the write-up:
- the protein character of the models' cassette exons against all measured exons (S6, extended to
  all three tasks in X4);
- the SHAP dependence of the six most shared events (S7).

### Checks and caveats

- **Neighbouring genes.** 163 of 1,105 overlapping top events read another gene more than their
  own (S10). They are named as such throughout.
- **RNA quality.** Tumours retain more intronic sequence (S2), but RNA integrity raises intronic
  reads too, and TCGA has no RIN for these samples.
- **Few normal samples.** ESCA has 11, and every splicing view drops there.
- **Outcome.** The survival signal of S8 is modest, halves with age, sex and stage, and is mostly
  renal at the level of single events.
- **Site.** Not a confound for this task: a patient's tumour and normal sample come from the same
  site.

## Cancer type: 34 classes, pan-cancer

**Full write-up: [results/CANCER_TYPE.md](results/CANCER_TYPE.md)** (figures T1–T9) — *splicing is
a lineage code, and it fails where a lineage code should.* One 34-class model per view — the 33
TCGA projects and a pooled normal class — over 9,729 samples.

### Performance

<!-- RESULTS-TYPE:START (generated by scripts/06_summarise_results.py) -->
**Out-of-fold, one 34-class model per view:**

| View | Accuracy | Balanced accuracy | Macro F1 |
|---|---|---|---|
| AFE | 0.918 | 0.914 | 0.904 |
| ALE | 0.898 | 0.890 | 0.882 |
| HIT | 0.920 | 0.920 | 0.909 |
| SE | 0.941 | 0.923 | 0.918 |
| A5SS | 0.917 | 0.907 | 0.901 |
| A3SS | 0.920 | 0.912 | 0.905 |
| MXE | 0.909 | 0.903 | 0.893 |
| RI | 0.933 | 0.921 | 0.918 |
| GEX | 0.954 | 0.936 | 0.930 |

<!-- RESULTS-TYPE:END -->

<img src="results/summary/fig_cancer_type_metrics.png" width="417" alt="Out-of-fold accuracy, balanced accuracy and macro F1 of the 34-class model, by view">

**Figure T0. Every view identifies the cancer type, expression best and the splicing views close
behind.** Out-of-fold accuracy, balanced accuracy and macro F1 of the 34-class model, by view, with
95% case-level bootstrap CIs.

Expression is the best single view (accuracy 0.954); the best splicing view is SE (0.941) and the
weakest ALE (0.898). Balanced accuracy and macro F1 sit below accuracy for every view because the
small classes are the hard ones; per-class recall is in
[fig_cancer_type_recall](results/summary/fig_cancer_type_recall.png) and T2.

[caption](results/summary/fig_cancer_type_metrics.caption.md)

### What the models learned, in figures

Five of the eight figures, in the order the argument runs. Each has a caption file and a CSV of
the plotted values next to it.

#### 1. No splicing view matches expression alone, but together they do, and fusing beats both

<img src="results/figures/cancer_type/T1_modality_gain.png" width="442" alt="T1: 34-class accuracy by modality, with the ensembles">

**Figure T1. No splicing view matches expression alone, but their ensemble does, and fusing the two
beats both.** Out-of-fold accuracy over all 9,729 samples and 34 classes for each single view (grey;
GEX black), the AS ensemble (the mean of the 8 splicing views, blue), the equal-weight fusion with
expression, and a cross-fitted stacker (diamond). Whiskers = 95% bootstrap CIs; brackets = exact
McNemar, raw p.

Every splicing view alone is below expression (0.898–0.941), yet their average reaches 0.950 against
expression's 0.954, with a *higher* macro recall (0.940 against 0.936): the views' errors are
largely independent of each other. Fusing the two gives 0.961: 109 of expression's errors corrected,
38 introduced, p = 3.9e-9.

[caption](results/figures/cancer_type/T1_modality_gain.caption.md) ·
[§1](results/CANCER_TYPE.md#1-combining-the-modalities-identifies-the-tissue-better-than-expression-alone)

#### 2. The two modalities fail in ways that say what each is reading

<img src="results/figures/cancer_type/T3_class_confusions.png" width="376" alt="T3: the 16 most-confused class pairs, GEX against the AS ensemble">

**Figure T3. Expression confuses anatomical neighbours; splicing confuses cell types across
organs.** The 16 class pairs with the most out-of-fold confusions, either direction, for expression
and the AS ensemble.

Expression mixes up **anatomical neighbours**: colon with rectum (96 errors), oesophagus with
stomach (40), lung adenocarcinoma with squamous carcinoma (37), glioma grades (20). Splicing makes
almost none of those (GBM against LGG: zero, printed as a number because a zero draws no bar) and
instead confuses **cell types across organs**: squamous carcinomas of cervix, head and neck, lung
and bladder, and several tumours with the normal tissue they arose in. That is what a readout of the
cell of origin should do.

[caption](results/figures/cancer_type/T3_class_confusions.caption.md) ·
[§3](results/CANCER_TYPE.md#3-the-two-modalities-confuse-different-classes)

Two checks stand behind that reading. The class-defining splice-site events correlate with their
host gene's expression no more than random events do, and each class's top splicing genes and top
expression genes share about one gene ([T5](results/figures/cancer_type/T5_not_expression.png)).
And an apparent pathway richness in the splicing explanations turned out to come from deeply
covered ribosomal-protein splice sites; what survives their removal is lineage biology and the
splicing machinery's own transcripts ([T4](results/figures/cancer_type/T4_pathway_coherence.png),
[§4](results/CANCER_TYPE.md#4-do-the-views-agree-on-the-processes-behind-a-class)).

#### 3. One exon can mark one lineage

<img src="results/figures/cancer_type/T6_lineage_atlas.png" width="880" alt="T6: two cassette exons per class, median inclusion across all 34 classes">

**Figure T6. Single cassette exons mark single lineages.** The two cassette exons with the largest
per-class importance for each of the 34 classes (rows), against every class (columns, clustered).
Colour = the exon's median inclusion in a class minus its median over all classes; dots = the class
that ranks it.

The coloured diagonal is the finding: MLANA is skipped in the melanocytic tumours, RAG2 in thymoma,
SLC37A2 included in the squamous carcinomas, TAGLN and DENND2B in LAML. Off-diagonal blocks are
exons shared by related classes.

[caption](results/figures/cancer_type/T6_lineage_atlas.caption.md) ·
[§6](results/CANCER_TYPE.md#6-single-exons-mark-single-lineages)

#### 4. One event can name a class on its own

<img src="results/figures/cancer_type/T8_class_events.png" width="683" alt="T8: one event per class, its class median against the nearest class and the other 32">

**Figure T8. A single splicing event can name a lineage on its own.** For each class, of the five
events the models rank highest across the eight splicing views, the one whose class median is
furthest from the nearest other class; the ten clearest are drawn. Filled dot = the class the event
marks, open dot = its nearest class, grey dots = the other 32. No model is involved in these
numbers.

TOMM40L's terminal exon is used as a last exon in every class except liver, where it reverses: a gap
of 1.63. KRT19 and BSCL2 are the two graded cases, and they are where one event stops being enough,
which is the argument for the ensembles of Figure T1. Two of the ten read a neighbouring gene rather
than their own (Figure T9): PLEKHG3's exon follows SPTB, and BSCL2's follows LRRN4CL.

[caption](results/figures/cancer_type/T8_class_events.caption.md) ·
[§8](results/CANCER_TYPE.md#8-one-event-per-lineage-in-detail)

#### 5. One headline gain that is not biology

<img src="results/figures/cancer_type/T7_site_probe.png" width="826" alt="T7: a fixed probe on six class pairs, grouped by patient against whole sites held out">

**Figure T7. Only the colon–rectum separation, and part of oesophagus–stomach, depends on the
tissue-source site.** In TCGA every tissue-source site contributes to exactly one project, so
separating neighbouring projects can exploit site handling rather than biology. One plot per pair of
related classes: a fixed, untuned probe cross-validated grouped by patient (open) and with whole
sites held out (filled); GEX grey, splicing views blue.

Splicing readouts such as intron retention are plausibly more sensitive to handling than expression
is. Glioma grade, kidney subtype and lung histology **hold** (drop ≤ 0.015). Colon against rectum
**fails**: no splicing view separates them at unseen sites, while expression keeps 0.66. Oesophagus
against stomach loses 0.08–0.14 in the splicing views (expression 0.06). The colon/rectum advantage
in T2 and T3 is not biology and is marked as such throughout. The confound is specific to this task:
the tumour-status models compare tissue within patient and within site.

[caption](results/figures/cancer_type/T7_site_probe.caption.md) ·
[§7](results/CANCER_TYPE.md#7-how-much-of-this-is-the-tissue-source-site)

#### 6. Which class-defining events measure another gene

<img src="results/figures/cancer_type/T9_overlap_audit.png" width="860" alt="T9: share of events overlapping another gene's exon, all measured against the class-defining events, and how many track that gene">

**Figure T9. Many of the cancer-type models' terminal exons read a neighbouring lineage gene.** The
audit of Figure S10 on the cancer-type models' top events. (a) The share of events whose segment
overlaps another gene's exon, for every measured event (grey) and the class-defining top events
(dark), by view. (b) Of those overlapping top events, how many follow the neighbour's expression
more than their own gene's.

61% of the top terminal (HIT) exons overlap another gene's exon, against 7% of all measured ones
(odds 12.6 with length held fixed), and 101 of 695 overlapping top events follow the neighbour. Here
the neighbours are mostly lineage genes (ceruloplasmin in liver, HNF1A in rectum, PAX6 in glioma),
so these events carry lineage information the expression model already reads. The comparison of T5
is with each event's *own* gene and stands; this is the other gene.

[caption](results/figures/cancer_type/T9_overlap_audit.caption.md) ·
[§9](results/CANCER_TYPE.md#9-which-class-defining-events-measure-another-gene)

Not shown here, all in the write-up: per-class recall against expression (T2), pathway coherence
(T4), and whether class-defining events are expression in disguise (T5).

### Checks and caveats

- **Site.** The splicing views' colon/rectum advantage is a tissue-source-site effect (T7) and is
  not counted as biology. Glioma grade, kidney subtype and lung histology hold at unseen sites.
- **Neighbouring genes.** 101 of 695 overlapping top events read another gene (T9). That gene is
  mostly a lineage gene the expression model already reads.
- **Small classes.** Balanced accuracy and macro F1 sit below accuracy for every view because the
  small classes are the hard ones (T2).

## Subtype within a cancer: 22 cohorts in 15 cancers

**Full write-up: [results/SUBTYPE.md](results/SUBTYPE.md)** (Figures U1–U26) — *within a cancer,
splicing recovers lineage and genome state about as well as expression, and loses where a subtype
is a single gene switched on.* One model per (cohort, view) over 22 subtype cohorts from 15 TCGA
projects:
- 12 histology cohorts, 7 genomic or IHC, and 3 expression-defined;
- 198 models over 7,552 patient–cohort pairs, one primary tumour per patient.

Every label map, exclusion and minimum class size (20 patients) was fixed in
[config/subtype.yaml](config/subtype.yaml) before the first model was fitted
([docs/SUBTYPE_METHODS.md](docs/SUBTYPE_METHODS.md), D-30 to D-45); the regulator, site,
MSI-repeat and BRCA analyses of Figures U8–U26 were fixed in [config/followup.yaml](config/followup.yaml)
before they were run (D-46, D-47, D-52 to D-74).

### Performance

**Out-of-fold balanced accuracy, mean over the cohorts of each label family:**

| View | Histology (12) | Genomic / IHC (7) | Expression-defined (3) |
|---|---|---|---|
| AFE | 0.765 | 0.759 | 0.721 |
| ALE | 0.740 | 0.729 | 0.659 |
| HIT | 0.760 | 0.800 | 0.751 |
| SE | 0.802 | 0.821 | 0.794 |
| A5SS | 0.778 | 0.741 | 0.693 |
| A3SS | 0.753 | 0.747 | 0.710 |
| MXE | 0.745 | 0.750 | 0.705 |
| RI | 0.784 | 0.781 | 0.739 |
| GEX | 0.836 | 0.882 | 0.897 |
| AS ensemble (mean of 8) | 0.804 | 0.822 | 0.812 |
| GEX + AS (1:1) | 0.834 | 0.883 | 0.909 |

Expression-defined labels (PAM50, HNSC and PCPG mRNA classes) were computed from mRNA, so GEX is
near-circular on them; for the splicing views they are a fair test. Per-cohort numbers with
chance and majority baselines, CIs and fusion tests are in
[SUBTYPE.md §2](results/SUBTYPE.md#2-how-well-each-modality-recovers-the-subtype) and
`results/subtype/summary/`.

<img src="results/figures/subtype/U1_performance.png" width="885" alt="U1: balanced accuracy per cohort for the eight splicing views, GEX, the AS ensemble and the fusion">

**Figure U1. Splicing matches expression on lineage and genome-state subtypes, and fusing the two
adds little.** Out-of-fold balanced accuracy of every cohort, grouped by how its label was made. The
grey bar spans the eight splicing views; the other marks are GEX (black), the AS ensemble (blue) and
GEX + AS fused 1:1 (dark blue). The tick is chance.
- **Close to expression.** Outside the expression-defined labels, the AS ensemble is within 0.03
  of GEX in 11 of 19 cohorts. It is ahead in 4: LGG histology +0.026, CESC, STAD and ESCA.
- **Well behind.** It falls more than 0.10 short in 3: PCPG site −0.17, PRAD ETS class −0.15 and
  UCEC molecular class −0.12.
- **Fusion.** Fixed fusion moves balanced accuracy by −0.04 to +0.04. Its CI clears zero only in
  LGG and UCEC_TCGA, and neither McNemar test reaches p < 0.05.

[caption](results/figures/subtype/U1_performance.caption.md) ·
[§2](results/SUBTYPE.md#2-how-well-each-modality-recovers-the-subtype)

### What the models learned, in figures

Each figure has a caption file and a CSV of the plotted values next to it. Every biological
statement in the write-up is marked as a result of this data or as context from the literature.

#### 1. Splicing loses the small classes that one programme defines

<img src="results/figures/subtype/U2_recall_by_class.png" width="900" alt="U2: recall of every class for GEX, the AS ensemble and the fusion">

**Figure U2. Splicing loses most on small classes defined by one programme.** (a) Histology cohorts.
(b) Genomic and IHC cohorts, then the expression-defined ones. Out-of-fold recall of every class for
GEX (black), the AS ensemble (blue) and GEX + AS fused 1:1 (dark blue); patients in brackets.
- **Losses.** Splicing's five largest losses against expression are all in classes of fewer than
  40 patients (GEX → AS ensemble):
  - POLE-ultramutated endometrial cancer, 0.54 → 0.08;
  - myxofibrosarcoma, 0.67 → 0.29;
  - ETV1-fused prostate cancer, 0.93 → 0.57;
  - paraganglioma, 0.83 → 0.57;
  - tall-cell thyroid carcinoma, 0.71 → 0.51.
- **Ensemble shrinkage.** In all five, at least one single splicing view does better than the
  ensemble (tall cell: SE 0.77; the best view is picked after the fact). Averaging eight views that
  are weak on a small class pulls its probability down.
- **Gains.** The ensemble gains on large classes, such as copy-number-high endometrial cancer
  (0.84 → 0.94).

[caption](results/figures/subtype/U2_recall_by_class.caption.md) ·
[§3](results/SUBTYPE.md#3-which-subtypes-splicing-loses)

#### 2. A histology model's errors follow the molecular class

<img src="results/figures/subtype/U3_histology_errors.png" width="816" alt="U3: true histology against the GEX model's call, broken down by the published molecular class">

**Figure U3. Histology models err in the molecular direction.** Four histology cohorts with a
published molecular class for the same patients. Counts only; no molecular label entered any
histology model.
- **Glioma.** The model calls 41 of 58 non-codeleted oligodendrogliomas astrocytoma, which is the
  reclassification WHO made in 2016, and all 6 codeleted astrocytomas oligodendroglioma.
- **Endometrial.** 9 of the 12 endometrioid tumours it calls serous are copy-number high, against
  28 of 100 overall.
- **PCPG.** Every adrenal tumour it calls extra-adrenal (11 of 11) is pseudohypoxic. Their PNMT,
  the adrenal enzyme the pseudohypoxic tumours lack, is at a median of 1.6 log2 CPM, against 7.9
  in the adrenal tumours called correctly.

The "wrong" calls are mostly the patients whose histology and molecular state disagree.
[caption](results/figures/subtype/U3_histology_errors.caption.md) ·
[§4](results/SUBTYPE.md#4-when-a-histology-model-is-wrong-its-error-is-molecularly-coherent)

#### 3. The literature's markers are in the data and the expression models, rarely in splicing

<img src="results/figures/subtype/U4_literature_markers.png" width="900" alt="U4: AUROC of each literature marker and its rank in the GEX and splicing models">

**Figure U4. The literature's markers are in the data and in the expression models, far less in the
splicing models.** The 40 gene markers below were fixed in advance with their sources; one row per
marker.
- **In the data.** 39 of 40 differ in the expected direction (p < 0.05); SLC5A5 in THCA does not.
- **In the models.** 27 rank in the top 100 of their GEX model; only 8 have an event in a splicing
  model's top 100. Among those 8 are CDH1's first exon in BRCA (rank 2), ERG's cassette exon in
  PRAD (rank 2) and CDKN2A's p16 first exon in HPV-positive HNSC (rank 9).
- **Glioma 1p/19q.** The 1p/19q arm dosage separates the classes in both glioma cohorts (AUROC
  0.85 for 1p in histology, 0.996 against codeletion). The GEX model's top genes are enriched on
  1p/19q for the molecular label (33 of 100, against 9.4 expected) but not for histology (7 of
  100).

[caption](results/figures/subtype/U4_literature_markers.caption.md) ·
[§5](results/SUBTYPE.md#5-the-literatures-biology-in-the-data-and-in-the-models)

<img src="results/figures/subtype/U5_as_vs_gex_genes.png" width="496" alt="U5: genes shared between each splicing model's top 50 and the GEX model's top 50">

**Figure U5. The splicing and expression models of a cohort read different genes.** How many genes a
cohort's splicing and expression models share in their top 50.
- **Overall.** The median is 1, and 55 of the 176 cohort × view pairs share none.
- **Retained introns** overlap most. Where a gene is barely transcribed, its intron looks
  retained: ESCA's RI model shares HNF1A, LGALS4, MYO7B and KIF12 with expression.
- **MXE in sarcoma** shares the smooth-muscle machinery (ACTN1, TPM1, TPM2, PDLIM7, LDB3).
- **Errors.** The two modalities still repeat most of each other's errors: the AS ensemble repeats
  68% of expression's.

[caption](results/figures/subtype/U5_as_vs_gex_genes.caption.md) ·
[§6](results/SUBTYPE.md#6-splicing-and-expression-read-different-genes)

#### 4. The strongest single events, audited for neighbouring genes

<img src="results/figures/subtype/U6_event_audit.png" width="900" alt="U6: the strongest splicing event per cohort, its class medians and its neighbour-gene audit">

**Figure U6. Four of the nineteen strongest splicing events read a neighbouring gene.** For each
histology and genomic/IHC cohort, the one splicing event with the largest class-specific importance,
with its median in each class. On the right is the D-27 audit: overlap with another gene's exon
(most TCGA libraries are unstranded) and Spearman with that gene's and the host's expression.
- **Readouts of a neighbour (4 of 19):** BRCA ZNF276 → FANCA, UCEC and UCEC_TCGA ERCC2 → KLC3,
  and HNSC_HPV VRK2 → FANCL.
- **Not cleared (2):** they overlap a gene the expression matrix does not hold.
- **Across all top events:** 97 of 2,248 top-5 events track a neighbour, 62 of them in HIT.
- **Named events that pass the audit** (SUBTYPE.md §7):
  - the smooth-muscle mutually exclusive exons of ACTN1 and TPM1 in leiomyosarcoma;
  - two frameshifting CDH1 exons whose lowest decile is skipped in lobular breast cancer;
  - the p16INK4a first exon of CDKN2A in HPV-positive head-and-neck cancer.
- **A named event that fails it:** the leiomyosarcoma "NDE1 promoter" is MYH11 read on the other
  strand.

[caption](results/figures/subtype/U6_event_audit.caption.md) ·
[§7](results/SUBTYPE.md#7-the-strongest-single-events-and-which-of-them-measure-a-neighbouring-gene)

#### 5. MSI leaves the same splicing marks in three cancers (exploratory)

<img src="results/figures/subtype/U7_msi_events.png" width="656" alt="U7: AUROC of ten MSI events in endometrial, gastric and colon cancer">

**Figure U7. Ten MSI splicing events replicate across endometrial, gastric and colorectal cancer.**
Splicing events from the MSI top 5 of the UCEC, STAD or COAD models, measured as MSI against the
rest in all three cancers, with no model involved. Ten events separate MSI in all three (AUROC ≥
0.75 or ≤ 0.25 in each).
- **Chosen independently:** 7 of the 10 were chosen by at least two of the three models. A 57-bp
  CCNDBP1 exon was chosen by all three (AUROC 0.83, 0.93 and 0.97).
- **Audit:** none overlaps another gene.
- **Immune infiltrate:** none follows it. With PTPRC, CD8A and CIITA, |ρ| ≤ 0.37, and the
  separation holds in the low-infiltrate half of each cohort.
- **RPL22L1:** the RPL22L1 event fits the RPL22 frameshifts common in MSI. This is context; RPL22
  mutation status was not tested.
- **Status:** the analysis was added after the CCNDBP1 recurrence was noticed (amendment A-4), so
  it is exploratory.

[caption](results/figures/subtype/U7_msi_events.caption.md) ·
[§8](results/SUBTYPE.md#8-msi-leaves-the-same-splicing-marks-in-three-cancers-exploratory)

#### 6. Splicing regulators of the subtypes

<img src="results/figures/subtype/U8_regulators.png" width="900" alt="U8: predicted splicing-factor differences between subtypes, and the ACTN1 and TPM1 events by PTBP1 level in sarcoma">

**Figure U8. Named splicing regulators differ between subtypes as predicted in four of six tests.**
(a) Each factor's expression by class: dot = median, bar = interquartile range, blue = the class
predicted lower. (b) The ACTN1 and TPM1 events by PTBP1 quintile, cut within leiomyosarcoma (blue)
and within the other sarcomas (grey). Six predictions, fixed before any splicing factor was looked
at by subtype (D-46).
- **Supported:**
  - ESRP2 is lower in biphasic than in epithelioid mesothelioma (AUROC 0.82);
  - ESRP1 is lower in the mesenchymal HNSC class (0.76);
  - PTBP1 is lower in leiomyosarcoma (0.66);
  - adjusting for class and composition, higher PTBP1 moves the ACTN1 smooth-muscle exon away
    from the leiomyosarcoma pattern (partial ρ +0.43; +0.57 within leiomyosarcoma), as PTB's
    repression of that exon predicts.
- **Not resolved:** ESRP1 in mesothelioma, where it is barely expressed.
- **Against the mechanism:** the TPM1 exon 2/3 event also moves away from the leiomyosarcoma
  pattern as PTBP1 rises (ρ −0.22). That matches the sign written down in advance, but the sign
  misread its source: in TPM1, PTB and MBNL1 repress the *non*-smooth-muscle exon (A-5). MBNL1
  follows the event most closely, in the published direction.
- **Specificity** (exploratory): PTBP1 ranks 2nd of 70 factors for ACTN1 and 12th for TPM1, so
  this is a PTB–MBNL axis rather than one factor.

[caption](results/figures/subtype/U8_regulators.caption.md) ·
[§9](results/SUBTYPE.md#9-splicing-regulators-of-the-subtypes)

#### 7. Part of the splicing accuracy is specific to the site

<img src="results/figures/subtype/U9_site_probe.png" width="833" alt="U9: loss of balanced accuracy when tissue-source sites are held out, per cohort, and against site–class association">

**Figure U9. Part of the splicing models' subtype signal is tied to the site that sent the sample.**
The site probe of T7 on 21 subtype cohorts (D-47): a fixed classifier is cross-validated grouped by
patient, then with whole tissue-source sites held out. (a) Loss of balanced accuracy per cohort:
single splicing views (grey), their median (blue), GEX (black). (b) The loss against the site–class
association (Cramér's V).
- **Expression transfers.** GEX loses ≤ 0.035 balanced accuracy in 20 of 21 cohorts; THCA, where
  tall-cell tumours come from 5 sites, loses 0.116.
- **Splicing loses more, in particular cohorts** (median of the 8 views): HNSC_HPV 0.098, CESC
  0.077, THCA 0.077, HNSC_mRNA 0.073. Outside THCA, GEX loses −0.001 to 0.014 in these cohorts.
- **Pre-specified tests:**
  - "splicing loses more than GEX in most cohorts" holds in 16 of 21 (sign-test q = 0.080): not
    resolved;
  - "the loss grows with site–class association" is not supported. ESCA, the most confounded
    cohort, loses almost nothing.

[caption](results/figures/subtype/U9_site_probe.caption.md) ·
[§10](results/SUBTYPE.md#10-is-the-subtype-signal-tied-to-the-site-that-sent-the-sample)

#### 8. MSI splicing follows acceptor poly-T repeats in the genome

*A journal-style account of this subsection, with its own legends and a model schematic (Figure M1):
[manuscript/MSI_SPLICING_CODE.md](manuscript/MSI_SPLICING_CODE.md).*

<img src="results/figures/subtype/U10_msi_repeats.png" width="900" alt="U10: share of MSI-model top events with an acceptor poly-T run, a genome-only repeat score by MSI status, and its AUROC by run length">

**Figure U10. Genome sequence alone says which exons MSI tumours skip.** (a) Share of exons with an
acceptor T-run ≥ 10: all exons against the top 50 SE events of the cohort's MSI model. (b) A
genome-only repeat score per patient, MSI against the rest. (c) Its AUROC by run length, and for
runs deep in the intron (hollow).

Why do the Figure U7 events change in MSI? Mismatch-repair-deficient tumours cannot repair
slippage at mononucleotide repeats, and a shortened poly-T just before an exon weakens its 3′
splice site. That much is known for MRE11 exon 5, for HSP110 exon 9, and genome-wide in
colorectal cancer. Exons were chosen from the reference genome alone: the longest T-run in the 50
intron bases before each cassette exon. The predictions were fixed before any PSI was read this
way (D-52).
- **The models had found the repeats.** 40 of the STAD and 34 of the COAD MSI model's top 50 SE
  events carry a T-run ≥ 10, and 10 of UCEC's, against 5% of all exons. So do 7 of the 9 SE
  events behind Figure U7, among them CCNDBP1 (T30).
- **The genome alone ranks MSI tumours.** A per-patient score with no model and no training
  reaches AUROC 0.71 in STAD and 0.78 in COAD. Endometrium does not separate (0.51).
- **Only at the acceptor.** The same score built from T-runs deep in the intron is worse by
  0.17–0.41 AUROC in all three cancers, endometrium included.
- **Beyond pyrimidine content.** The MSI shift grows with run length after adjusting for the
  acceptor's pyrimidine content (STAD, COAD).
- **The same exons in three organs.** Their shifts correlate between cancers (ρ 0.43–0.65).
- **Positive control.** MRE11 exon 5 replicates in all three.

<img src="results/figures/subtype/U11_msi_long_tracts.png" width="870" alt="U11: share of repeat exons more included in MSI by run length, the A3SS replication, and the exons with the largest shared shifts">

**Figure U11. The longest acceptor repeats move toward inclusion in MSI.** (a) Share of repeat exons
more included in MSI, by run length. (b) The A3SS replication: use of the acceptor that carries the
run, MSI minus the rest. (c) The exons with the largest shared shifts, labelled with run and exon
length.

The longest tracts go the other way. Exons after runs ≥ 20 nt are *more* included in MSI
in all three cancers (56–67% move up, against 38–49% for runs of 10–19). The largest increases
(ANXA7, CCNDBP1, TBK1, …) all carry runs of 21–32. This was seen after the fact.
- **Replication on unread A3SS events** (D-54): 2 of 6 tests supported, all 6 in the predicted
  direction.
- **No mechanism yet:** the hnRNP C mechanism written down for it was not supported (D-53).
- **Not the hospital** (D-55): the score holds within hospital (AUROC 0.68 STAD, 0.81 COAD), after
  adjusting for site, infiltrate and RNA handling, and in the less-infiltrated half of each cohort.
- **Caveat:** the intronic contractions are inferred, not seen. Reading them in the exome reads is
  pre-specified (D-56).

<img src="results/figures/subtype/U16_msi_deletion_load.png" width="820" alt="U16: repeat score against MC3 homopolymer deletions per tumour in STAD, COAD and UCEC, and MSI tumours' deletion counts by cancer">

**Figure U16. Mismatch-repair damage tracks the repeat score in stomach cancer only.** (a) Repeat
score against MC3 homopolymer deletions per tumour in each cancer; the grey column is the rest of
the cohort. (b) Deletion counts of the MSI tumours by cancer, with each cancer's D-52 AUROC.

Dose shows in one cancer (D-61). Within STAD MSI tumours, more MC3 homopolymer deletions go
with a lower repeat score (partial ρ −0.45, q = 0.005); COAD and UCEC are not resolved. UCEC MSI
tumours carry the most deletions (median 147) yet have the weakest score, and their deletions are
not shorter, so deletion size does not explain the endometrial gap.
- **No NMD signature in gene expression** (D-62, not resolved).
- **GTEx genotypes no variant inside any run**, so the germline test (D-64) could not be run.

<img src="results/figures/subtype/U17_msi_branchpoint.png" width="820" alt="U17: branch-point distance from the 3' splice site by run length, the share of repeat exons more included in MSI by branch-point tertile, and the adjusted correlations">

**Figure U17. The branch point's distance from the 3′ splice site sets the direction of the MSI
shift.** (a) Branch-point distance from the AG by run length; grey band = the 19–23 nt optimum. (b)
Share of repeat exons more included in MSI by branch-point tertile, within each run class. (c)
Adjusted correlation of the MSI shift with branch-point distance (filled) and with the run's far end
(hollow).

The branch point sets the direction (D-63, fixed in advance; supported in all three
cancers). An exon whose branch point, read from sequence, lies further from the AG gains inclusion
in MSI, beyond its run length (adjusted ρ 0.11–0.16, permutation p ≤ 0.001). The long-run reversal
of Figure U11 sits in exons whose branch point is far: among runs ≥ 20, 48% move up when it is near and
73% when it is far. A contraction that brings a distant branch point toward the AG would fit
(Chua and Reed 2001). The run's far end carries nearly the same information, so the two are not
separated.

<img src="results/figures/subtype/U18_msi_code_pancancer.png" width="900" alt="U18: AUROC of score A and score G for MMR-deficient tumours in 30 held-out cancers, per-exon shifts in held-out cancers against the training cohorts, and score G against deletion load">

**Figure U18. A two-arm sequence code finds mismatch-repair-deficient tumours in 30 cancers never
used to define it.** (a) AUROC of score A (hollow) and score G (filled) for MMR-deficient tumours:
pooled, without the tissue relatives, against hypermutated proficient tumours, and per cancer. (b)
Per-exon shift in the 30 cancers against the training cancers (blue = R_up). (c) Score G against the
deletion count within the 50 MMR-deficient tumours.

The sequence code finds mismatch-repair-deficient tumours in 30 cancers never used to
define it (D-66, fixed in advance; 7 of 8 tests supported).
- **The code** (score G): runs of 10–19 nt fall; runs ≥ 20 nt with the branch point ≥ 36 nt from
  the AG rise. Score A, D-52's, reads every repeat exon as falling.
- **The test.** Among 7,685 tumours, 50 carry ≥ 10 MC3 homopolymer deletions. Score G finds them
  at AUROC 0.83, against 0.55 for score A, with no re-fitting.
- **Organ independence and specificity.** 0.81 without READ, ESCA and UCS; 0.85 against
  hypermutated proficient tumours.
- **Exon by exon.** The per-exon shifts replicate (ρ 0.51), and so do the direction rule and the
  branch-point relation.

<img src="results/figures/subtype/U19_msi_code_specificity.png" width="820" alt="U19: share of repeat exons more included in MSI by branch-point distance, the adjusted correlations for MSI and for control contrasts, and AUROC of score A and score G in the three cohorts">

**Figure U19. The code is specific to the branch point and to mismatch repair.** (a) Share of repeat
exons more included in MSI by branch-point distance, for runs of 10–19 and ≥ 20 nt. (b) Adjusted
branch-point correlation for MSI (blue) and for controls (grey): exons without a tract, POLE in
UCEC, CIN against GS in STAD. (c) AUROC of score A and score G in the three training cancers.

The code is about the branch point and about mismatch repair (D-65, D-67).
- **Beyond the run.** The branch point adds beyond the run's own extent.
- **Controls.** Exons without a tract show no relation, and neither do POLE or CIN-vs-GS
  contrasts.
- **The switch** sits near 30 nt from the AG (exploratory).
- **Endometrium.** Score G separates UCEC MSI at 0.84, where score A gave 0.51; in-sample for the
  rule.

<img src="results/figures/subtype/U20_msi_code_infiltration.png" width="820" alt="U20: AUROC of score G by immune half and of the immune score alone, and the share of the MMR effect on score G kept after adjusting for infiltration, RNA handling, purity and site">

**Figure U20. The code is not an infiltration readout.** (a) AUROC of score G in all held-out
tumours and in each immune half, of the immune score alone (grey), and in the less-infiltrated half
of each training cancer. (b) The effect of MMR deficiency (or MSI) on score G kept after adjusting
for stromal, immune and retained-intron levels, with purity or site; 95% CI.

MSI tumours are heavily infiltrated, so the code could in principle read the infiltrate (D-68,
fixed in advance; 6 of 6 tests supported).
- **Adjusted, the effect stays whole.** Adjusting for stromal and immune expression and the
  retained-intron level keeps all of it: 1.01 in the 30 held-out cancers, 1.01–1.07 in the
  training cancers.
- **Low infiltrate, same separation.** AUROC 0.84 in the less-infiltrated half of each held-out
  cancer (21 MMR-deficient tumours), and 0.85–0.98 in the training cancers' lower halves.
- **Infiltration points the other way.** The immune score alone marks MMR deficiency weakly
  (0.59), and in the adjusted model more immune expression goes with a *higher* score G.

<img src="results/figures/subtype/U21_msi_code_msi_calls.png" width="900" alt="U21: AUROC of score G against MC3, MANTIS and MSIsensor calls, and score G of held-out tumours by MC3 call and by each caller's call">

**Figure U21. The code finds the tumours that two exome MSI callers call MSI-high.** (a) AUROC of
score G against D-66's MC3 call and against MANTIS and MSIsensor MSI-H, held out, then in-sample in
the three training cancers. (b, c) Score G percentile of the held-out tumours by MC3 call (grey
proficient, blue deficient) and by the caller's call (mid = between its two cut-offs).

Every earlier held-out test took its truth from MC3 deletions. MANTIS and MSIsensor call MSI from
the exome reads by other means (D-69, fixed before the 6.2 MB cBioPortal download).
- **Both callers agree with the code.** Score G separates MSI-H from MSS at 0.97 (MANTIS, 17 MSI-H)
  and 0.89 (MSIsensor, 23), in 30 cancers it never saw.
- **It marks the MSI-high state, not the score below it.** Over all 7,388 scored tumours, nearly all
  MSS, the MANTIS score shows no relation to score G (ρ +0.01; MC3 not resolved).
- **Where the calls disagree.** All 17 MANTIS MSI-H tumours are MC3-deficient; the 33 MC3-deficient
  tumours MANTIS calls mid or MSS still score low on G.

<img src="results/figures/subtype/U22_msi_code_rising_arm.png" width="900" alt="U22: the held-out rise of R_up exons against other subtype contrasts, exon shifts by class, long runs by branch-point distance, the R_up arm by DNA class, and against deletion count">

**Figure U22. The rising arm is specific to mismatch repair and needs the tract, but its mechanism
is not settled.** (a) Mean shift of R_up exons in 20 against 20 tumours: 57 other subtype contrasts
(grey, absolute value; dashed = their 95th percentile), held-out MMR-deficient tumours (blue; bar =
95% of draws) and, in-sample, MSI in the three training cohorts (hollow). (b) Held-out exon shifts
by class. (c) Exploratory: long runs (≥ 20 nt), share more included by branch-point distance. (d)
The R_up arm of held-out tumours by DNA class. (e) The R_up arm against homopolymer deletions in
MMR-deficient and MSI tumours.

Long tracts should weaken a splice site, yet the R_up exons are included *more* in MMR-deficient
tumours. D-70 tested that eight ways, fixed in advance (3 supported, 1 against, 4 not resolved).
- **The rise is real.** In 30 cancers never used, 98 of 135 R_up exons rise (median +0.028).
- **Not a tumour-wide shift.** Adjusting for each tumour's average inclusion, infiltration and
  retained-intron level leaves the effect unchanged (+0.236 against +0.237).
- **It needs the tract.** Exons with an equally far branch point but no run do not move.
- **Mismatch repair, not mutation burden.** Hypermutated proficient tumours do not show it
  (AUROC 0.74).
- **Not unique to MSI.** The rise beats 55 of 57 other subtype contrasts, but two sarcoma contrasts
  move the same exons as much (not resolved).
- **Mechanism open.** The rise follows the branch point's distance (panel c, exploratory: within 30
  nt of the AG only 5 of 25 long-run exons rise, beyond it 60–84%). But the branch-point-free AG
  exclusion zone does not predict it, runs upstream of the branch point rise too, it does not grow
  with deletion count, and PTB relief runs the wrong way.

<img src="results/figures/subtype/U23_msi_code_models.png" width="900" alt="U23: held-out ROC and precision-recall of the repeat-exon logistic model against score G, AUROC per held-out cancer, and its largest coefficients">

**Figure U23. Logistic regression on repeat exons, fitted in three cancers, finds MMR deficiency in
30 others.** (a) ROC in the 30 held-out cancers, against score G. (b) Precision against recall
held out. (c) AUROC per held-out cancer with three or more MMR-deficient tumours. (d) The 20
largest standardised coefficients (gene, acceptor T-run, branch-point distance).

A weight per exon, on the 1,562 exons with an acceptor run of 10 nt or more (D-71; shown alone
by A-19):
- **It transfers.** Trained on COAD, STAD and UCEC, it reaches AUROC 0.900 in the other 30 cancers
  (score G 0.832), AUPRC 0.44 (0.08), and 0.9997 against MANTIS MSI-H.
- **Its weights are the code.** Long runs with a far branch point carry the largest positive
  weights (CCNDBP1 T30); runs of 10–19 nt carry negative ones.
- **Pre-registered.** It was predicted to transfer worse than ten sequence bins; it did not (LR4,
  not resolved).

<img src="results/figures/subtype/U24_msi_callers_benchmark.png" width="900" alt="U24: within-cancer AUROC of each exome MSI measure and the repeat-exon model against the other two measures' consensus, per cancer, at matched call volume, and where the DNA measures disagree">

**Figure U24. Within each cancer, the repeat-exon model ranks tumours almost as well as the exome
MSI callers.** (a) Within-cancer AUROC of each exome measure judged by the other two, and of the
model on the same tumours (filled: each cancer left out of training; hollow: the 30 held-out
cancers). (b) Per cancer. (c) Sensitivity and PPV at each caller's call volume, pooled over
cancers. (d) MC3-deficient tumours by the callers' calls: share in the model's top decile.

Each of MANTIS, MSIsensor and the MC3 count is judged against the tumours on which the other two
agree (D-72):
- **Within a cancer, the model reaches 0.998–1.000.** That is 0.002 or less below the exome
  measures, from RNA alone.
- **Across cancers it is not calibrated.** At each caller's call volume it finds 82–84% of
  positives, against the callers' 94–100%.
- **Pre-registered tests.** These were on the all-exon model (4 of 6 supported, small deficits);
  the repeat-exon numbers are exploratory.

<img src="results/figures/subtype/U25_msi_model_vs_callers.png" width="900" alt="U25: tumours called MSI by each combination of the repeat-exon model, MC3, MANTIS and MSIsensor; calls per cancer; and the model's score against each exome measure">

**Figure U25. The model agrees with the exome measures on its confident calls; its own extra calls
are borderline.** (a) Tumours called MSI by each combination of the model, MC3, MANTIS and
MSIsensor. (b) Calls per cancer by agreement. (c–e) The model's score against each exome measure;
red = model only, black = exome only, rings = three possible missed cases.

Of the model's 168 calls, 115 are confirmed by all three exome measures (D-73):
- **Extra calls.** The 32 it makes alone sit just above its threshold, show no MMR-gene silencing,
  and mostly look like false calls.
- **Clear misses.** Of the tumours all three exome measures call, it misses 20; 19 are endometrial.

<img src="results/figures/subtype/U26_msi_splicing_vs_expression.png" width="900" alt="U26: held-out AUROC and AUPRC of the splicing model against expression models trained the same way; performance by MLH1 status; the two models' ranks of the held-out MMR-deficient tumours; within-cancer agreement with the exome callers">

**Figure U26. Expression models trained the same way transfer almost as well; splicing adds to
them.** (a, b) Held-out AUROC and AUPRC in the 30 cancers, for models trained on COAD, STAD and
UCEC. (c) MMR-deficient tumours with and without MLH1 silencing. (d) Each deficient tumour's rank
under splicing and expression. (e) Within-cancer agreement with the exome callers.

Expression retrained on all genes came close to splicing held out: 0.89 against 0.90 (D-74).
- **Fixed panels.** PreMSIm's 15 genes did not transfer (0.70–0.79).
- **Complementary.** The two models ranked deficient tumours differently. Together they reached
  0.93 and doubled the precision at half recall.

[U10 caption](results/figures/subtype/U10_msi_repeats.caption.md) ·
[U11 caption](results/figures/subtype/U11_msi_long_tracts.caption.md) ·
[U16 caption](results/figures/subtype/U16_msi_deletion_load.caption.md) ·
[U17 caption](results/figures/subtype/U17_msi_branchpoint.caption.md) ·
[U18 caption](results/figures/subtype/U18_msi_code_pancancer.caption.md) ·
[U19 caption](results/figures/subtype/U19_msi_code_specificity.caption.md) ·
[U20 caption](results/figures/subtype/U20_msi_code_infiltration.caption.md) ·
[U21 caption](results/figures/subtype/U21_msi_code_msi_calls.caption.md) ·
[U22 caption](results/figures/subtype/U22_msi_code_rising_arm.caption.md) ·
[U23 caption](results/figures/subtype/U23_msi_code_models.caption.md) ·
[U24 caption](results/figures/subtype/U24_msi_callers_benchmark.caption.md) ·
[U25 caption](results/figures/subtype/U25_msi_model_vs_callers.caption.md) ·
[U26 caption](results/figures/subtype/U26_msi_splicing_vs_expression.caption.md) ·
[§11](results/SUBTYPE.md#11-the-msi-splicing-marks-are-written-in-the-genome-acceptor-poly-t-repeats)

#### 9. Breast cancer with DNA: two HER2s, splice-site mutations, and what the pair check removed

<img src="results/figures/subtype/U13_brca_two_her2.png" width="900" alt="U13: PAM50 Her2 probability by ERBB2 amplification for GEX and the AS ensemble, AUROCs for each comparison, and the amplicon knockout">

**Figure U13. Splicing reads HER2 through the HER2-enriched programme, not through the 17q12
amplicon.** (a) P(Her2) by PAM50 group and ERBB2 amplification, for GEX and the AS ensemble. (b)
AUROC for each comparison. (c) The amplicon knockout: the luminal AUROC with the 17q12 inputs set to
their mean (hollow), against random input sets of the same size (tick).

Which HER2 did the PAM50 models learn? ERBB2 amplification was first taken from ABSOLUTE
copy number (D-58).
- **Splicing seemed to see amplification better.** In luminal tumours, the AS ensemble separated
  amplified from not at AUROC 0.81, against 0.72 for GEX (H2, *against* our prediction). Figure U15 does
  not hold this up.
- **It reads a programme, not the amplicon.** Setting the 17q12 amplicon's inputs to their mean
  costs GEX 0.03 and no splicing view more than 0.007 (A-12, exploratory).
- **Within PAM50 Her2, neither model tells amplified from not** (0.48).

<img src="results/figures/subtype/U15_brca_her2_uncapped.png" width="660" alt="U15: AUROC of P(Her2) for three amplified groups, GEX against the AS ensemble, and P(Her2) against ERBB2 copy number per ploidy in amplified luminal tumours">

**Figure U15. Without the ABSOLUTE ceiling, splicing matches expression on ERBB2 amplification but
does not beat it.** (a) AUROC of P(Her2), GEX against the AS ensemble, for D-58's amplified tumours,
all tumours amplified by the uncapped call, and the newly amplified ones. (b) P(Her2) against copy
number / ploidy in amplified luminal tumours (hollow = newly amplified).

ABSOLUTE never reports more than 7 copies at ERBB2, so D-58 could not call any tumour of
ploidy above 3.5. D-60, fixed before the file was downloaded, takes the SNP6 copy ratio without a
ceiling.
- **The ceiling had left out genome-doubled tumours.** 15 of 45 amplified luminal tumours are new,
  with median ploidy 3.85.
- **Splicing's lead does not survive.** AS 0.81 against GEX 0.75 (HU2, not resolved). On the newly
  amplified tumours, AS 0.82 and GEX 0.83.
- **Both models follow the amplification level** (ρ 0.47 and 0.45). Splicing still loses nothing
  when the amplicon's inputs are removed. It reads a downstream programme whose strength follows
  the dose.

<img src="results/figures/subtype/U14_brca_mutations.png" width="620" alt="U14: PSI percentile of exons next to splice-site mutations, and a cryptic A3SS score in SF3B1-hotspot tumours">

**Figure U14. Splice-site mutations skip their exons, and SF3B1 hotspot tumours use cryptic
acceptors.** (a) PSI percentile of the exon next to each splice-site mutation, in the mutated
tumour; dashed = no effect. (b) A cryptic-acceptor score in SF3B1-hotspot tumours against the rest.

Mutations leave the expected marks (D-59).
- **Splice-site mutations skip their exon.** 110 of 157 put it below the cohort median (median
  percentile 0.09).
- **SF3B1-hotspot tumours use cryptic acceptors.** The score separates them at AUROC 0.90, and all
  14 are luminal.
- **CDH1 explains part of the lobular tail.** Its 9 splice-mutant lobular tumours are the extreme
  of the lobular CDH1 skipping tail, not all of it.

<img src="results/figures/subtype/U12_brca_germline_check.png" width="900" alt="U12: tumour-normal pair tables for CAST and TMEM101, tumour-normal agreement across events, and high-group frequency in LumA against Basal tumours">

**Figure U12. Most of the BRCA models' genotype-like events are subtype programmes, not inherited
alleles.** (a, b) Each normal's group against its own tumour's group, for CAST (a germline event)
and TMEM101 (the models' top genotype-like event). (c) Tumour–normal agreement (Cohen's κ) for all
events, those passing the check and the models' top events. (d) High-group frequency in LumA against
Basal tumours, per event.

Are the models' "genotype-like" splicing events inherited? A label-blind scan found 1,594
events whose PSI falls into three genotype-like groups, and the BRCA models favour them (D-57). The
pre-specified LOH and ancestry tests mostly held. Then a check from 112 tumour–normal pairs was
added (A-13): an inherited allele cannot turn a homozygous normal into the opposite homozygote.
- **Most events change group in their own tumour.** Only 82 of 1,110 pass the check, and 20
  agree at κ > 0.6. CAST is one of them (a).
- **The models' top events are not inherited.** 3 of 45 pass, and none agrees above κ 0.24
  (TMEM101, b). Their frequency differs between LumA and Basal tumours (d): they are subtype
  programmes.
- **On the confirmed events, LOH still shows** (L1–L3), and ancestry differences are commoner
  (55% of genes, against 32%).

[U12 caption](results/figures/subtype/U12_brca_germline_check.caption.md) ·
[U13 caption](results/figures/subtype/U13_brca_two_her2.caption.md) ·
[U14 caption](results/figures/subtype/U14_brca_mutations.caption.md) ·
[U15 caption](results/figures/subtype/U15_brca_her2_uncapped.caption.md) ·
[§12](results/SUBTYPE.md#12-breast-cancer-what-the-subtype-models-read-beyond-the-labels)

### Checks and caveats

- **Site** (Figure U9). Part of the splicing models' accuracy on HPV status, the HNSC expression classes,
  cervical histology and gastric molecular class is specific to the sites that sent the samples;
  expression's is not.
- **Neighbouring genes** (Figure U6). 97 of 2,248 top events read a neighbour, and every named event is
  given with its audit status.
- **Label provenance.**
  - Expression-defined labels make GEX near-circular.
  - The PRAD ETS-fusion calls used RNA-seq.
  - HPV status is known for 112 of 501 HNSC patients and is confounded with tumour site.
- **Small classes.** Classes of 20–40 patients have recall intervals up to ±0.20.

## Across the three tasks

**Full write-up: [results/ACROSS_TASKS.md](results/ACROSS_TASKS.md)** (figures X1–X4). *The three
tasks' models share genes and exon types, and where the type and status models hesitate, they
hesitate as the literature predicts.*

Nothing is refitted: every number comes from the three tasks' out-of-fold predictions and
held-out SHAP values. The predictions, their sources and a list of what had already been seen
were committed in [config/followup.yaml](config/followup.yaml) before anything was computed
([docs/ACROSS_TASKS_METHODS.md](docs/ACROSS_TASKS_METHODS.md), D-48 to D-51).

| Predictions | GEX | AS ensemble | GEX + AS |
|---|---|---|---|
| T1–T13: a subtype carries more type probability for the projects it is said to resemble | 11 of 13 | 13 of 13 | 12 of 13 |
| S1–S6: a subtype looks less tumour-like to its cancer's status model | 2 of 6 | 6 of 6 | 6 of 6 |

The model-level predictions:
- **Supported:** G1 (every pair of tasks shares genes above chance) and E1 (subtype exons have the
  other tasks' protein character).
- **Not supported:** G2 (histology subtypes sit closer to the type model than to the status model)
  and C1 (stromal and immune content explains the status probabilities).

#### 1. The three tasks read overlapping genes

<img src="results/figures/across_tasks/X1_gene_overlap.png" width="744" alt="X1: top-100 gene overlap between each pair of tasks, and for histology cohorts subtype against status and type">

**Figure X1. The three tasks' models read overlapping genes, several-fold more than chance.** For
each subtype cohort, the top 100 genes per view of three models: the subtype model, its cancer's
status model, and the type model's importance over that cancer's own samples. (a) Overlap of each
pair of tasks against chance, one dot per cohort and view; bar = median. (b) The five histology
cohorts with a status model: Jaccard of the subtype list with the status list and with the type
list.

- **Every pair overlaps** 5–8 times more than chance (87–94% of 387 comparisons at q < 0.05).
- **No lineage preference** (b). Histology subtypes are no closer to the type model than to the status
  model.
- **The shared core is splicing.** 402 of the 413 gene entries found in all three lists of a cancer
  come from splicing views. Most are switches in smooth-muscle and stromal genes: TPM1 in 10
  cohort-views, CHN1, ACTN1. They track tissue composition, which every task's samples vary in.

[caption](results/figures/across_tasks/X1_gene_overlap.caption.md) ·
[§1](results/ACROSS_TASKS.md#1-the-three-tasks-read-overlapping-genes)

#### 2. Which other cancers a subtype resembles

<img src="results/figures/across_tasks/X2_type_probability_by_subtype.png" width="900" alt="X2: cancer-type probability of the predicted target projects in each subtype group, and AUROC by model">

**Figure X2. In the cancer-type model, subtypes resemble the cancers the literature predicts.**
Thirteen pre-specified predictions, one per row: a subtype, the projects it should resemble (→) and
its comparison group. (a) The mean summed out-of-fold probability of the target projects in each
group, GEX + AS (log scale). (b) AUROC per model, higher group against lower; filled = BH q < 0.05
within the model.

The pan-cancer type model never saw a subtype label. Where it spreads probability beyond a tumour's
own project, it does so as the literature predicts:

- squamous oesophageal and cervical tumours lean towards head-and-neck and lung squamous
  carcinoma;
- oesophageal adenocarcinoma leans towards stomach (18 of 81 called STAD, against 0 of 80 squamous
  tumours);
- cervical adenocarcinoma leans towards endometrium;
- IDH-wildtype lower-grade glioma leans towards glioblastoma;
- copy-number-high and serous endometrial cancer lean towards ovary and carcinosarcoma;
- biphasic mesothelioma leans towards sarcoma;
- HPV-positive head-and-neck cancer leans towards cervix (AUROC 0.92).

Intestinal-type gastric cancer leaning towards colon is not resolved. One caution (A-8, not
pre-specified): basal-like breast tumours carry more ovarian probability only because they carry
more off-target probability of every kind.

[caption](results/figures/across_tasks/X2_type_probability_by_subtype.caption.md) ·
[§2](results/ACROSS_TASKS.md#2-which-other-cancers-a-subtype-resembles)

#### 3. How tumour-like the status model finds each subtype

<img src="results/figures/across_tasks/X3_status_probability_by_subtype.png" width="900" alt="X3: mean one minus P(tumour) by subtype group, AUROC by model, and P(tumour) against stromal and immune scores">

**Figure X3. The splicing status models see the six predicted subtypes as less tumour-like.** Six
pre-specified predictions, one per row: the subtype expected to look less tumour-like against its
comparison. (a) Mean 1 − P(tumour) in the predicted group (filled) and the comparison (open), for
GEX and the AS ensemble (log scale). (b) AUROC per model; filled = BH q < 0.05 within the model. (c)
Spearman ρ of P(tumour) with the stromal and the immune score over each cancer's tumours.

The splicing status models find all six predicted subtypes less tumour-like:

- normal-like breast cancer (AUROC 0.84; 13% scored below 0.5, against ≤ 0.5% of other PAM50
  classes);
- follicular-variant thyroid cancer (tall cell is more tumour-like than classical);
- diffuse and genomically stable gastric cancer;
- lobular breast cancer.

GEX resolves two of the six, because its probabilities sit at the ceiling (medians 0.998–1.000), so
this is partly calibration. Stromal and immune content does not explain it (c): P(tumour) falls with
those scores in three cancers and rises in two.

[caption](results/figures/across_tasks/X3_status_probability_by_subtype.caption.md) ·
[§3](results/ACROSS_TASKS.md#3-how-tumour-like-the-status-model-finds-each-subtype)

#### 4. The same kind of exon in every task

<img src="results/figures/across_tasks/X4_exon_protein_character.png" width="546" alt="X4: protein character of the top cassette exons of the three tasks against all measured exons">

**Figure X4. The subtype models' cassette exons share the protein character of the other two tasks'
exons.** Cassette exons in the top 100 of any tumour-status model (blue), any cancer-type class
(dark blue) or any subtype cohort's SE model (amber), against 6,000 random measured exons (grey):
the share divisible by 3, in a Pfam/CDD domain, in a disordered region and in a transmembrane helix;
Fisher exact tests against the background.

The subtype models' top cassette exons (1,857) carry the signature first seen in the other two tasks
(S6), slightly weaker:

- more often divisible by 3 (49% against 41% of all measured exons);
- less often in a Pfam/CDD domain (39% against 51%);
- more often in a disordered region (18% against 14%).

Whatever separates tumour from normal, tissue from tissue or subtype from subtype, the splicing
models reach for exons that rewire interaction segments rather than break folded domains.

[caption](results/figures/across_tasks/X4_exon_protein_character.caption.md) ·
[§4](results/ACROSS_TASKS.md#4-the-exons-each-task-uses-do-the-same-kind-of-thing-to-the-protein)

### Checks and caveats

- **When the predictions were written.** They were written after the main results. Some
  project-level confusions (ESCA ↔ STAD, CESC ↔ HNSC) had already been seen, as the config
  records; only the subtype split was new.
- **Calibration.** X3 reads probabilities that GEX pushes to the ceiling.
- **A sensitivity that is not pre-specified.** The off-target share (A-8) changes the reading of
  T6.

## How the models are built

```
9,729 TCGA samples ─► 5 outer folds (grouped by patient, stratified by label)
                         │
                         ├─ outer-train ─► 3 inner folds ─► Optuna (30 / 20 trials)
                         │                   preprocessing fitted on inner-train only:
                         │                   missingness filter → rank (variance | ANOVA) → top K
                         │                   → median-impute → standardise → MLP
                         │                 best config refitted on outer-train
                         │
                         └─ outer-test ──► predictions (pooled → metrics + bootstrap CIs)
                                           SHAP DeepExplainer on every held-out sample
```

- **No leakage.** Every statistic, from feature filtering to epochs, is fitted on training rows
  only. A patient's tumour and normal always share a fold. Hyper-parameters never see the
  outer-test fold. See [docs/LEAKAGE_POLICY.md](docs/LEAKAGE_POLICY.md); a proof covers the one
  cohort-level pre-filter.
- **Same samples, every modality.** Splits depend only on the cohort, so modalities are
  compared on identical folds.
- **SHAP-interpretable by construction.** The MLP uses only layers that `shap.DeepExplainer`
  handles exactly (Linear, BatchNorm, ReLU, Dropout). Additivity is asserted for every
  explained sample.
- **Every choice is recorded** in [docs/DECISIONS.md](docs/DECISIONS.md) (D-01…D-12), with values in
  [config/study.yaml](config/study.yaml).

## Reading the outputs

For each job, `results/<task>/<view>/<cohort>/` contains:

| File | Content |
|---|---|
| `metrics.json` | pooled out-of-fold metrics, 95 % case-bootstrap CIs, per-fold mean ± SD, config, git commit |
| `folds.csv` | per outer fold: tuned hyper-parameters, number of features, inner loss, refit epochs, test metrics, SHAP additivity error |
| `shap_top_features.csv` | top 200 features by mean \|SHAP\|: `label` (gene + genomic region of the variable segment), `n_folds_selected` (0–5), `direction_spearman` |
| `shap_top_genes.csv` | importance summed per host gene |
| `shap_top_features_by_class.csv` | *cancer type only*: top 30 features for each of the 34 classes, with per-class direction |
| `per_class_metrics.csv`, `confusion_matrix.csv` | *cancer type only* |
| `fig_shap_beeswarm.png` | *tumour status*: each dot is one held-out sample, x = SHAP (tumour log-odds), colour = feature value |
| `fig_shap_class_heatmap.png`, `fig_confusion_matrix.png` | *cancer type* |
| `*.caption.md` | what each figure shows (figures carry no in-plot titles; [docs/FIGURE_GUIDELINES.md](docs/FIGURE_GUIDELINES.md)) |

**How to read `direction_spearman`:**
- Tumour status: > 0 means higher PSI, usage, HIT or expression pushes a sample towards
  *tumour*; < 0 means towards *normal*.
- Cancer type: per class, > 0 means a higher value raises the evidence for that class.

**Units:** SHAP values are in logit units. Task 1 uses log-odds of tumour. Task 2 uses the
centred logit zₖ − z̄ of the sample's true class (see [docs/METHODS.md](docs/METHODS.md)).

## Repository layout

```
config/       study.yaml (every rule and the search space), sources.yaml (raw data and reference paths),
              biology.yaml (interpretation thresholds), subtype.yaml (label maps, fixed before any fit),
              followup.yaml (predictions of the subtype regulator/site and cross-task analyses)
src/ctsm/     library
  cohort.py splits.py annotation.py extract.py matrix_io.py   data
  preprocess.py model.py pipeline.py metrics.py explain.py    modelling
  biology.py                                                  interpretation helpers (GOseq, protein mapping, ...)
  subtype.py                                                  subtype labels, configuration and summaries
  figstyle.py plots.py                                        house figure style and phase-1 figures
scripts/      00_register_sources → 01_build_cohort → 02_extract_matrices → 03_validate_data
              → 04_run_job (one job) / 05_run_all (all jobs, parallel) → 06_summarise_results
              → 07_shap_importance → 08–19, 22–25 interpretation analyses (tables only)
              → 20_figures_status / 21_figures_type (every figure)
              subtype task: 00b_fetch_gdc_labels, 00c_fetch_msi_scores → 30_subtype_labels → 31_subtype_extract
              → 32_subtype_run_job / 33_subtype_run_all → 34_subtype_summarise
              → 35_subtype_interpret → 37_subtype_regulators, 38_subtype_site_probe → 36_subtype_figures
              across tasks: 40_cross_task → 41_figures_cross_task
              run_pipeline.sh
tests/        unit tests (leakage guards, encodings, preprocessing, SHAP additivity, biology helpers)
docs/         METHODS, DECISIONS, LEAKAGE_POLICY, BIOLOGY_METHODS, SUBTYPE_METHODS, ACROSS_TASKS_METHODS,
              FIGURE_GUIDELINES
results/      data/ (cohort, splits, extraction, validation summaries)
              tumor_status/<view>/<cancer>/, cancer_type/<view>/PANCAN/, summary/, RESULTS.md
              biology/01_complementarity … 14_survival_robustness/ (tables)
              subtype/<view>/<cohort>/, subtype/{data,summary,interpret,regulators,site_probe}/
              across_tasks/ (tables)
              figures/tumour_status/ (S1–S10), figures/cancer_type/ (T1–T9), figures/subtype/ (U1–U9),
              figures/across_tasks/ (X1–X4)
              TUMOUR_STATUS.md, CANCER_TYPE.md, SUBTYPE.md, ACROSS_TASKS.md
```

## Reproducing

```bash
conda env create -f environment.yml
```

```bash
bash scripts/run_pipeline.sh
```

```bash
OMP_NUM_THREADS=1 python scripts/04_run_job.py tumor_status se BRCA
```

- The second command runs everything end to end. The third runs a single job.
- **Raw inputs:** paths are in `config/sources.yaml`. These are TCGA rMATS/HITindex tables
  and STAR counts produced from GDC BAM files, plus GENCODE annotation. They are not
  redistributed.
- **Local-only outputs:** processed matrices, model weights, SHAP arrays and per-sample
  predictions are written to the git-ignored `data.nosync/`.
- **Provenance of the published results:** the 144 jobs were run on 2026-10-10/11 with the
  code of commit `b249c3c`. Each `metrics.json` records `git_commit: null` because the runs
  finished before that first commit. Its `config` block records the exact settings. The 198
  subtype jobs ran on 2026-10-12 from a working tree whose subtype code was not yet committed,
  so their `git_commit` is the commit of the day plus `+dirty`; that code was committed
  afterwards, with the results.
- **Follow-up analyses:** the predictions of `config/followup.yaml` were committed in `741a1a1`
  before steps 37, 38 and 40 were first run (2026-10-12). Those steps refit no task model.
- **Subtype labels:** the TCGA-CDR histology table (path in `config/sources.yaml`) and two
  open-access GDC PanCanAtlas tables (published subtype calls and clinical fields), fetched and
  MD5-checked by `scripts/00b_fetch_gdc_labels.py` into `data.nosync/`. The MANTIS and MSIsensor
  scores for D-69 come from the public cBioPortal API (`scripts/00c_fetch_msi_scores.py`,
  6.2 MB, MD5-pinned).
- **Environment note:** set `OMP_NUM_THREADS` before Python starts. In this environment
  (torch 2.2.2 + conda-forge numpy, macOS), importing numpy before torch with multiple OpenMP
  threads segfaults; `scripts/_common.py` handles it.

## Related

The same TCGA processing and rules are used by
[cancer-survival-predictions](https://github.com/zachpwakefield/cancer-survival-predictions)
(survival modelling with alternative RNA processing).
