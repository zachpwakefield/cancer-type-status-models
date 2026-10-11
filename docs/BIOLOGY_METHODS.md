# Biology layer: methods and decisions

Scripts `07`–`14` interpret the fitted models. None of them changes or refits a phase-1 model;
step 14 fits separate, untuned diagnostic probes. No result here feeds back into a modelling
choice, so the outer-test estimates in [`RESULTS.md`](../results/RESULTS.md) stay unbiased.
Every threshold is fixed in [`config/biology.yaml`](../config/biology.yaml) or in the
script's header. Steps 08–14 write tables only. Step 15 draws every figure from those tables in
the house style ([FIGURE_GUIDELINES.md](FIGURE_GUIDELINES.md)), each with a caption file and a
CSV of the plotted values.
The narrative is in [`results/BIOLOGY.md`](../results/BIOLOGY.md).

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
