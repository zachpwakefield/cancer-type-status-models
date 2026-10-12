# Biology layer: methods and decisions

Scripts `07`–`22` interpret the fitted models. None of them changes or refits a phase-1 model;
steps 14 and 18 fit separate, untuned diagnostic models. No result here feeds back into a
modelling choice, so the outer-test estimates in [`RESULTS.md`](../results/RESULTS.md) stay
unbiased. Every threshold is fixed in [`config/biology.yaml`](../config/biology.yaml) or in the
script's header. Steps 07–19 and 22 write tables only; steps 20 and 21 draw the figures from those
tables in the house style ([FIGURE_GUIDELINES.md](FIGURE_GUIDELINES.md)), each with a caption
file and a CSV of the plotted values.

The two tasks answer different questions and are written up separately:
[`results/TUMOUR_STATUS.md`](../results/TUMOUR_STATUS.md) (figures S1–S9) and
[`results/CANCER_TYPE.md`](../results/CANCER_TYPE.md) (figures T1–T8).

Decisions continue the numbering of [`DECISIONS.md`](DECISIONS.md).

### D-13 SHAP importance for interpretation · `scripts/07` · `ctsm/biology.full_importance`
- **What:** every feature that any outer fold selected, ranked by mean |SHAP| over all
  held-out samples (0 for samples whose fold did not select it). Phase 1 kept only the top 200.
- **Per class (cancer type):** the same, restricted to held-out samples of that class. The
  explained output is the sample's own centred class logit (D-10), so a class table shows what
  makes the model call that class.
- **Direction:** Spearman of the (standardised) feature value against its SHAP value, over the
  samples where it was selected.
- **Gene level:** a gene's importance is the maximum over its events. A sum would favour genes
  with many measured events.
- **Cache:** `<data_dir>/bio_cache/importance/` (local; derived from per-sample SHAP arrays).

### D-14 Complementarity of views · `scripts/08`
- **Inputs:** out-of-fold probabilities only. Every sample was predicted by a model that never
  saw it.
- **Fusion:** the weights are fixed in advance, with no tuning.
  - `all_equal`: the mean of the 9 views.
  - `GEX + AS (1:1)`: 0.5 × GEX + 0.5 × the mean of the 8 AS views.
- **Tests:** exact McNemar on paired correctness. CIs come from a case-level bootstrap
  (1,000 resamples).
- **Caveat:** fused out-of-fold predictions are an honest estimate for fixed weights. They
  are not a nested estimate of a fusion method that chooses its weights.
- **Superseded for the narratives.** Step 18 (D-23) repeats this per cancer and per class, adds
  bootstrap CIs and a stacked comparator, and is the source of every number in
  `TUMOUR_STATUS.md` and `CANCER_TYPE.md`. Step 08's tables are kept for the error-overlap and
  per-class recall views that step 18 does not reproduce.

### D-15 Global shifts in matched pairs · `scripts/09`
- **Events tested:** ALL measured events of each cohort (extraction pre-filter), with no SHAP
  selection. This avoids circularity, because the models were trained to separate exactly these
  samples.
- **Patient-level statistic:** for each patient with a tumour and a matched normal sample, the
  mean paired difference over events observed in both. A Wilcoxon signed-rank test runs across
  patients, with BH correction over all cancer × view tests.
- **What "up" means in each view:**
  - RI: more retention.
  - A5SS/A3SS: more use of the longer exon form.
  - SE: more inclusion.
  - MXE: more of the upstream exon.
- **Terminal exons (AFE/ALE):** a gene-level position index, the usage-weighted rank of the
  used first or last exons along the direction of transcription (0 = most upstream,
  1 = most downstream). It needs ≥ 0.5 total usage among kept siblings.
- **Proliferation control:** the within-tumour Spearman of per-sample mean RI PSI against a
  10-gene proliferation score.
- **Caveat:** RNA integrity and the nuclear RNA fraction also raise intron reads. TCGA has no
  RIN for most samples, so a technical part of the RI shift cannot be excluded.

### D-16 Pathway over-representation · `scripts/10` · `ctsm/biology.goseq_ora`
- **Selection:** the top 200 genes by SHAP importance (D-13) of each model, and of each class
  for cancer type.
- **Universe:** genes with ≥ 1 measurable feature in that view and cohort.
- **Bias correction:** AS views measure several events per gene, so genes with many events
  are more likely to be picked. A probability weighting function of events per gene (monotone
  fit) gives the Wallenius non-central hypergeometric test (GOseq). For GEX it reduces to the
  hypergeometric test.
- **Gene sets and multiplicity:** MSigDB 2026.1 Hallmark, Reactome and GO-BP, with 15–500
  genes in the universe and ≥ 3 overlapping genes. BH runs within each collection and model;
  the main threshold is q < 0.05.
- **Sensitivity:** the per-class AS results are dominated by translation terms. These come from
  alternative 5'/3' splice sites in cytosolic ribosomal-protein genes, which are deeply covered
  and so measured precisely in every tissue. The per-class AS tests were repeated with
  RPL*/RPS* genes removed from both the ranking and the universe.

### D-17 Events and protein consequences · `scripts/11` · `ctsm/biology.ProteinMap`
- **Recurrence:** an event is recurrent if it is in the top 100 of ≥ 4 tumour-status models
  of the same view.
- **Chance expectation:** a Poisson-binomial model with p = 100 / events measured in each
  cohort, over the cohorts where the event is measured.
- **Cross-cancer change:** the median paired tumour − normal change is shown in every cancer,
  not only where the event ranked.
- **Protein mapping:** cassette exons are matched exactly to GENCODE v45 exons. The transcript
  priority is coding, then MANE Select, canonical and other. The exon's codons are located in
  the CDS and overlapped with InterPro features:
  - Pfam/CDD domains;
  - MobiDB-lite disordered regions;
  - TMHMM helices.
- **Comparison:** top exons vs a random 6,000 of all measured exons. Fisher exact tests, with
  BH correction across metrics and sets.

### D-18 Splicing regulators · `scripts/12`
- **Scope:** within the tumours of each cancer.
- **Statistic:** partial Spearman between RNA-binding-protein expression (log2 CPM) and event
  PSI. Both are rank-transformed and residualised on a stromal and an immune marker score.
- **Event sets:** the 50 top events of each view, and 500 random measured events of the same
  view as background.
- **Excess:** the share of top events with |ρ| > 0.3 minus the background share.
- **Paired expression change:** median tumour − normal log2 CPM over matched pairs.
- **Caveat:** these are correlations; a factor may be co-regulated with its "targets" rather
  than drive them.

### D-19 Lineage · `scripts/13`
- **Confusions:** pairwise out-of-fold confusions, GEX vs the AS ensemble.
- **Expression-proxy test:** Spearman of each top per-class event with its host gene's log2 CPM
  (3,000 random samples), against random measured events.
- **Gene overlap:** top-100 GEX genes vs top-100 AS host genes, per class. The expected
  overlap is computed within genes measured by both.
- **Atlas:** the top 2 SE events per class (per-class SHAP), with median PSI in every class.

### D-20 Site-confounding probe · `scripts/14`
- **The confound:** every TCGA tissue-source site (TSS) contributes to exactly one project. For
  example, COAD and READ share 0 of 38 sites, and GBM and LGG share 0 of 40. A model separating
  related projects could use site handling rather than biology.
- **Probe:** for 6 pairs of related classes and each view, a fixed, untuned pipeline is
  cross-validated by case and by TSS (every test site unseen). The pipeline, fitted inside each
  training fold, is top-2,000 variance (NaN-aware) → median impute → standardise → L2 logistic
  regression (C = 0.05).
- **Reading:** a view whose accuracy falls sharply when sites are held out relied on
  site-specific signal.
- **Scope:** diagnostic only. It does not re-estimate the phase-1 models.

### D-21 What is published
- **Committed:** aggregate tables and figures, as under D-12.
- **Local only:** per-patient and per-sample values, such as per-patient shifts and importance
  caches, stay in `<data_dir>/bio_cache/`.
- **Figures:** they may plot per-sample points, without identifiers. A figure's published CSV
  holds aggregates. The exact per-sample plotted values go to `<data_dir>/figure_data/`
  ([FIGURE_GUIDELINES.md](FIGURE_GUIDELINES.md)).

### D-22 Outcome-blind candidates and a matched background · `scripts/16`
- **Candidates:** the top 10 features by mean |SHAP| of every tumour-status model (8 AS views ×
  15 cancers), pooled: 1,014 distinct events.
- **Background:** for each candidate, one event of the same view drawn from the same quintile of
  observed fraction and of SD (both over that view's tumours in the 15 cancers), never in the
  top 100 of any model: 910 distinct events.
- **Why:** candidates are well measured and variable, and those two properties alone make an
  event easier to test for anything. The background holds them fixed, so a difference in hit
  rate is not a difference in measurability.
- **Blindness:** no survival time, event indicator, stage, age or sex enters either set. The
  models themselves never saw an outcome.
- The same split is reused by the regulon enrichment (D-24) and the survival layer (D-25).

### D-23 Stacked ensemble · `scripts/18`
- **Fixed-weight fusion is the headline** (D-14): nothing is estimated, so it cannot overfit.
- **Stacked:** a logistic regression on the 9 views' out-of-fold probabilities (standardised,
  C = 1, no tuning), cross-fitted over the same outer folds — fold *k* is scored by a stacker
  fitted on the other folds only.
- **Residual optimism:** the base models that produced the stacker's *training* probabilities
  were themselves trained on folds that include the samples the stacker later scores. The effect
  is second order but not zero, so the stacked number is reported as a ceiling and never as the
  headline.
- Bootstrap CIs resample cases (1,000 resamples). AUROC in the bootstrap is computed from the
  rank sum of the positives, which equals `roc_auc_score` with ties averaged.

### D-24 Splicing-factor regulons · `scripts/19`
- **Scope:** within the tumours of each cancer (normals would turn any tumour-normal difference
  into a correlation).
- **Universe:** 2,000 random events per view (one seeded list, the same in every cancer) plus
  every D-22 candidate and background event; in a given cancer an event needs ≥ 80% observed
  values and non-zero variance, and the remaining missing values are median-imputed.
- **Statistic:** partial Spearman of factor expression against event PSI, both rank-transformed
  and residualised on a stromal and an immune marker score (D-18). A pair is a member in that
  cancer if |ρ| ≥ 0.3 and p < 0.05.
- **Core regulon:** the pairs that recur with the same sign in ≥ 5 of the 15 cancers.
- **Enrichment:** candidates against the matched background (D-22), Fisher exact.
- **Caveat:** these are correlations inside tumours. A factor may be co-regulated with its
  apparent targets rather than drive them, and a 0.3 threshold on ρ is a convention, not a
  biological boundary.

### D-25 Survival association of the model's events · `scripts/16`, `scripts/17`
- **Tool:** splice-assay 0.4.0 (an external, separately tested package), run over a local
  dataset built from the same matrices the models used.
- **Cohorts:** the 15 tumour-status cancers, 6,152 tumours, 1,620 deaths (OS) and 960 (DSS).
  Survival and clinical values come from the TCGA clinical table; `is_missing` is read as
  missing.
- **Models per event × cancer:** Cox on PSI with the host gene's own expression (base), and the
  same plus age, sex and stage, against stage I and female (adjusted). Also a median-split
  Kaplan-Meier with a log-rank test, and tumour-vs-normal group tests.
- **Gates (splice-assay defaults):** Cox needs 30 patients and 10 events, and is marked
  low-power below 20 events; the log-rank test needs 10 patients per arm and 10 events; an
  event needs 10 patients off its modal value. Cells that fail a gate are not counted in either
  arm of the enrichment test.
- **Multiplicity, two families, do not mix them.** `cox_q_cohort` / `km_q_cohort` are BH over
  every event tested in one cancer × endpoint × model (about 1,600 Cox tests) — the family that
  answers "does this one event stand out in this cancer", and the one the narrative quotes.
  `cox_q` / `km_q` are splice-assay's own, BH within one event over the cancers it was run in,
  which is what its assay pages star; the tool leaves them empty for families under 10, so they
  are mostly missing here. Figures show raw p. The enrichment test itself is one Fisher test per
  endpoint × model × candidate set, BH over those.
- **What is claimed:** association in observational cohorts, never cause and never a
  therapeutic target. PSI can reflect tumour purity or composition as well as splicing.
- **What is published:** the per-event × cohort statistics, as a small `event_index.csv` (one
  row per event) and `event_survival_cells.csv.gz` keyed by it — repeating the long rMATS ids
  over 115,000 rows would be most of the file. Numbers are written to 4 significant digits, and
  the untrimmed table stays in `<data_dir>/survival/cells_full.parquet`. The PSI, expression and
  clinical tables behind them stay in `<data_dir>/survival/` (D-12).
- **A permuted-outcome replicate** (survival shuffled between patients within each cancer, OS,
  base model, one replicate) measures what this pipeline reports when nothing is associated:
  `survival_null.csv`. It is the reference line in figure S8, in place of the nominal 5%.

### D-26 Event cards · `scripts/22`
- **What they are:** for the handful of events each narrative names, one row holding everything
  the pipeline already knows — geometry, protein consequence, per-cancer or per-class movement,
  how each model reads the event, regulon membership and survival. The narratives quote these
  rows; nothing here is a new statistic.
- **Which events, tumour status:** the recurrent events of D-17 (top 100 of >= 4 models), most
  models first, one per gene. No new selection.
- **Which events, cancer type:** per class, the five events with the largest mean |SHAP| in that
  class's own explanation across the 8 splicing views (the models' own ranking), then the one of
  those whose class median is furthest from the nearest other class, then the ten clearest
  classes. The second step decides **which examples to show and nothing else** — a feature can
  rank first for a class and still not separate it alone, because the model reads it together
  with the others, and an example that does not separate cannot be described concretely.
- **The alternative segment** is the piece the two forms differ by: the exon for SE and for the
  terminal views, the **difference between the two exons** for MXE (the swap replaces one with
  the other, so neither length alone sets the frame), the retained intron for RI, and the
  extension the long form adds for A3SS/A5SS. The reading frame is recorded only for internal
  events; a first or last exon replaces one end of the protein, so the field is left empty.
- **Protein consequence:** GENCODE v45, coding transcripts first, then MANE > canonical, and
  within that a transcript carrying InterPro features in preference to one carrying none — some
  of the newest MANE transcripts are absent from the feature table and reading their empty
  annotation as "no domain" would be wrong. An exon in no coding transcript of its gene is
  recorded as such rather than as "no consequence"; the protein map holds coding transcripts
  only, so that is what the field means.
- **Published:** the cards, the InterPro features of each card's transcript (for S9), and the
  per-class medians of each cancer-type card (for T8). Per-sample values stay local (D-12).
