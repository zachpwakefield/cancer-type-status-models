# Tumour versus normal tissue: what splicing adds

One classifier per cancer, per modality: 15 TCGA projects with at least 10 normal cases, 9
modalities each, 135 models. Every number here is out-of-fold — each sample was scored by a
model that never saw it ([RESULTS.md](RESULTS.md)).

This document asks what the splicing models learned about becoming a tumour. The other task,
telling the cancers apart, is a different question with a different answer and is written up
separately in [CANCER_TYPE.md](CANCER_TYPE.md).

- **Methods and decisions:** [docs/BIOLOGY_METHODS.md](../docs/BIOLOGY_METHODS.md) (D-13 to D-25).
- **Scripts:** `07`–`19` write the tables, `20` draws these figures.
- **Figures:** S1–S8 in [figures/tumour_status/](figures/tumour_status/). Each has a caption file
  stating the question, the method, n and the test, and a CSV of the plotted values.
- **Statistics in figures:** raw p. Any q is in the CSV with its FDR family named.

## Summary

1. **Splicing adds to expression.** Over 6,815 samples, the expression model makes 99 errors and
   the splicing ensemble 91 — but not the same ones. Fused at fixed weights they make 77:
   27 of expression's errors corrected, 5 introduced (McNemar p = 1.1e-4).
2. **Tumours move RNA processing in one direction.** They retain more intronic sequence (RI up in
   7 cancers, A5SS long forms in 10, each down in at most 1) and choose more upstream last exons
   (7 cancers, none the other way).
3. **The two modalities see different things.** Expression's tumour signal runs through
   microenvironment programmes (EMT in 13 of 15 cancers); the splicing models' genes mostly do
   not form pathways at all.
4. **The models' events sit in specific splicing-factor programmes.** They are enriched in the
   small regulons of PTBP1, RBM47, QKI and HNRNPU, and depleted from the large shared programme
   of nuclear RNA-processing factors. The events in those small regulons are also the ones most
   often associated with survival (§8).
5. **Different cancers rely on the same events.** 404 events rank in the top 100 of at least 4
   models, against 3 expected by chance.
6. **The exons look tissue-regulated**: more often frame-preserving and disordered, less often
   inside folded domains, than measured exons in general.
7. **The events carry some outcome information, but less than the discrimination suggests.**
   Model-selected events reach a nominal Cox p < 0.05 in 10.5% of event × cancer pairs against
   9.1% for matched events (odds 1.18, p = 3.6e-4), with the host gene's expression in the
   model, and the rate rises with how many models agree on an event (odds 1.28 in the same arm
   for events recurrent in four or more, 1.56 without the clinical covariates). The events that
   survive FDR inside a cancer are almost all renal, and directions flip between cancers.
   Separating tumour from normal tissue and predicting outcome are distinct axes.

---

## 1. Combining the modalities beats expression alone

![S1](figures/tumour_status/S1_modality_gain.png)

**Pooled over the 15 cancers** (6,815 out-of-fold samples):

| Model | Errors | Mean AUROC |
|---|---|---|
| expression (GEX) | 99 | 0.992 |
| AS ensemble (mean of 8 views) | 91 | 0.994 |
| GEX + AS, 1:1 | **77** | **0.995** |
| stacked (cross-fitted) | 76 | 0.983 |

- **The splicing ensemble is not better than expression; it is different.** On its own it is
  level with GEX (54 samples right where GEX is wrong, 46 the other way, p = 0.48).
- **Fusion is where the gain is**: 27 of expression's errors corrected against 5 introduced,
  exact McNemar p = 1.1e-4. The weights were fixed before any result was seen (D-14).
- **The gain is not clever weighting.** A cross-fitted stacker, which is allowed to learn the
  weights and is mildly optimistic by construction (D-23), reaches the same error count and a
  *worse* mean AUROC. Equal weights are enough.
- **Where it helps most:** ESCA (8.7% → 7.0% of samples misclassified), KICH (1.1% → 0), LUAD
  (1.1% → 0.4%) and PRAD (4.4% → 3.6%). ESCA is both the hardest cohort and the one that gains
  most; it has the fewest normal samples of the fifteen (11).
- No single cancer has enough errors for its own McNemar test to reach significance; the pooled
  test is the one to read.

**The two tasks read mostly different events.** Taking each model's top 100 features, the events
the 15 tumour-status models use overlap those of the 34 cancer-type classes well above chance
(2–78 fold, `09_ensemble/task_feature_overlap.csv`) but still only partly: 28–36% of the
tumour-status events for the cassette-exon, splice-site and intron views, against 48% (AFE) and
63% (ALE) for the terminal-exon views, which read promoter and poly(A) choice and are the most
lineage-bound of the modalities. That is why the two tasks are written up separately.

## 2. Tumours retain more intronic sequence and choose upstream last exons

![S2](figures/tumour_status/S2_global_shifts.png)

This uses **every measured event**, not the models' selections, so the models cannot have biased
it: in matched tumour-normal pairs from the same patient, does RNA processing move in a common
direction?

**Retained introns.**
- 7 cancers retain more (BLCA, COAD, KICH, KIRC, KIRP, LUSC, PRAD); BRCA retains less.
- Retained introns that gain ≥ 5 PSI points in tumours outnumber those that lose them about 3:1
  (median 15% vs 4.5% of events), in 13 of 15 cancers.
- The shift is carried by most patients: in KICH, BLCA, PRAD and KIRC, 70–91% of patients retain
  more in their tumour than in their own normal tissue.

**Alternative splice sites.** A5SS long forms rise in 10 cancers and fall in none; A3SS long
forms rise in 5.

**Last exons.** They move upstream in 7 cancers and downstream in none — the direction of the
proximal polyadenylation widely reported in tumours.

**Alternative explanations.**
- *Proliferation.* Within tumours, per-sample retention does not track a 10-gene proliferation
  score (Spearman −0.25 to 0.19, no consistent sign).
- *3′-end factors.* How much a cancer shifts is not explained by how much its cleavage and
  polyadenylation factors rise (|ρ| ≤ 0.39 across 15 cancers).
- *RNA quality.* RNA integrity and nuclear-RNA content also raise intronic reads. TCGA has no RIN
  for these samples, so a technical component cannot be excluded.

## 3. Expression reads the microenvironment; splicing does not form pathways

![S3](figures/tumour_status/S3_pathways.png)

For each model, the 200 genes with the largest SHAP importance were tested against Hallmark,
Reactome and GO-BP, correcting for how many events each gene has (GOseq, D-16).

- **Expression** draws on stromal and microenvironment programmes: EMT in 13 of 15 cancers,
  coagulation in 10, KRAS signalling in 9.
- **Splicing** reaches the same programmes far less often (EMT in 7, myogenesis in 4), and its
  top genes form *any* enriched pathway in at most 9 of 15 cancers (SE) and in none for RI.

The same holds for GO-BP, where no term is reached by any AS view in more than 2 of the 15
cancers (the most recurrent being actin-filament-based movement and morphogenesis), against 3–4
for expression (`03_pathways/tumor_status_gobp_recurrence.csv`).

This is the clearest statement of why fusion helps: the splicing signal is spread over genes that
do not share a programme, so it is not a restatement of the expression signal. What the splicing
models use is better described event by event (§5–§7) than pathway by pathway.

## 4. The models' events live in small, specific splicing-factor programmes

![S4](figures/tumour_status/S4_regulons.png)

A factor's **core regulon** is the set of events whose PSI tracks its expression within tumours —
partial Spearman |ρ| ≥ 0.3 at p < 0.05, after residualising both on a stromal score, an immune
score and the sample's own mean PSI in that view — reproducibly, with the same sign, in at least
5 of the 15 cancers (D-24).

**That third covariate changes the answer.** Without it, every nuclear RNA-processing factor
appears to "regulate" more than a thousand events, 94% of them in the same direction, because
factor expression and overall PSI level share one sample-level axis. Residualising it out cuts
those regulons by about 60% (PABPN1 1,434 → 566 events) and is the difference between measuring
regulation and re-measuring that axis.

**Most large regulons are one programme, not many.** Among the 47 factors with at least 20 core
events, the median factor shares 42% of its regulon (Jaccard) with its closest neighbour, and
the 1,074 distinct events carry 7,529 memberships — about 7 factors each. PABPN1, RBM39, SRSF11,
HNRNPH1 and RBM5 are largely the same set of events read five times.

**The models' events avoid that programme.** Against the matched background (D-22), the top
events of the tumour-status models are:
- **enriched** in the small regulons of PTBP1 (15 candidate events against 0 background;
  odds 29 with a 0.5 correction, infinite without), RBM47 (5.0), HNRNPU (4.2), HNRNPF (3.6),
  QKI (3.4) and PRPF8 (2.7);
- **depleted** from the large shared ones: RBM5, HNRNPH1, PABPN1, SRSF11 (odds 0.49–0.57),
  RBM39, SRSF10, RBM4, SRSF4, SRSF6 (0.26–0.38).

All of these pass BH at q < 0.05 over the 32 factors with enough marked events to test. PTBP1,
QKI and RBM47 are differentiation regulators whose loss or gain is a recurring theme in
epithelial tumours, which is the kind of programme §2 and §6 point at. PCBP1 is the most
distinct regulon of all (Jaccard 0.10 with its closest neighbour) but is not enriched for model
events.

**Which factors change.** Of about 70 pre-specified RNA-binding proteins, ESRP1 rises in 10
cancers and the cleavage and polyadenylation factors in up to 8, while CELF2 (11 cancers),
RBM20 and RBM24 (10 each), NOVA2 (9), QKI and RBM47 (7 each) fall. Part of that fall is
composition: normal tissue carries more stroma and smooth muscle.

These are correlations inside tumours. A factor may be co-regulated with its apparent targets
rather than drive them.

## 5. Different cancers rely on the same events

![S5](figures/tumour_status/S5_recurrent_events.png)

Across the 8 AS views, 404 events rank in the top 100 (mean |SHAP|) of at least 4 of the 15
models. About 3 would be expected if the rankings were independent (Poisson-binomial over each
event's measured cohorts).

- By view: ALE 74, MXE 73, HIT 62, A3SS 51, RI 51, A5SS 47, AFE 39 and SE 7.
- Measured in every cancer, not only where they ranked, recurrent events share their majority
  direction in a median 73% of cancers.
- Recurrent RI and A5SS events mostly rise in tumours (73% and 77%), in line with §2.

**Examples.** A CENPU HIT event rises in all 15 cancers; an ENOSF1 HIT event and a GCOM1 MXE
event fall in 14; PDLIM7 and PLCG1 mutually exclusive exon events rise in all 15. Several genes
contribute pairs of events that move in opposite directions: two EXOSC7 last exons, two NDE1
first exons, two TPM1 mutually exclusive exon events.

**A caution.** The NR1D1 and THRA HIT events overlap at chr17:40.09 Mb, where the two genes lie on
opposite strands, and move in opposite directions in 14 of 15 cancers. Most TCGA RNA-seq is
unstranded, so these are better read as one finding than two.

## 6. The exons the models use look tissue-regulated

![S6](figures/tumour_status/S6_exon_consequences.png)

The top cassette exons of the tumour-status models (n = 1,281) and of the cancer-type classes
(n = 1,914), against 6,000 random measured exons:

- more often of a length divisible by 3: 51% and 54%, against 41% (p = 1.3e-10 and 3.1e-24);
- among fully coding exons, less often overlapping a Pfam/CDD domain: 37% and 35%, against 51%;
- among the same exons, more often overlapping a disordered region: 21% and 20%, against 14%;
- less often overlapping a transmembrane helix: 3.7% and 4.6%, against 6.7%.

So the signal sits in exons that change protein segments in frame, mostly in disordered,
interaction-prone regions, rather than breaking folded domains or shifting the reading frame.
These are the properties of exons that tissues regulate.

**Direction.** Exons skipped more in tumours are the most frame-preserving (55% divisible by 3)
and the most domain-depleted (30%, OR 0.42); exons included more in tumours are the most enriched
for disorder (24%, OR 1.98).

## 7. What the most shared events look like

![S7](figures/tumour_status/S7_event_portraits.png)

Six events are shown: the most recurrent, one per gene, whose paired change has the same sign in
at least 80% of cancers. Left, each cancer's matched normal and tumour medians; right, the SHAP
dependence in the cancer where the event ranks highest. The model reads all six monotonically —
the further a sample's value lies in the tumour direction, the more it is pushed towards
"tumour".

## 8. Do these events carry information about outcome?

![S8](figures/tumour_status/S8_survival.png)

A classifier that separates tumour from normal tissue is not, by construction, a prognostic
model. This section asks whether the events it relies on are nonetheless associated with how
long patients live.

**The design is what makes the question answerable.** Each of the models' top events is matched
to a background event of the same view, the same coverage quintile and the same spread quintile
(D-22), and both sets go through the same tool with the same gates: splice-assay 0.4.0, Cox on
the event's value with the host gene's own expression in the model, with and without age, sex
and stage (D-25). Nothing about outcome entered the choice of either set. The comparison is
therefore not about which events are easier to test.

**Model events are associated with survival more often than matched events, and the more
models agree on an event the more often it is.** Each cell below is one event in one cancer;
the rate is the share with Cox p < 0.05, against the matched background in the same arm.

| Set | Pairs | Base Cox | Adjusted Cox | Odds ratio, base / adjusted |
|---|---|---|---|---|
| every candidate | 11,702 | 11.1% | 10.5% | 1.27 / 1.18 |
| recurrent in ≥ 4 models (D-17) | 2,687 | 13.3% | 11.4% | 1.56 / 1.28 |
| top 10 for ≥ 2 cancers | 1,553 | 13.8% | 10.8% | 1.62 / 1.21 |
| matched background | 10,185 | 9.0% | 9.1% | — |

Every row is significant against the background (base arm p = 1.6e-7, 8.7e-11, 1.2e-8). The
gradient is the part worth keeping: an event that several independent per-cancer models pick out
is more likely to track outcome than a candidate picked by one, which is what we would expect if
the models are selecting biology rather than cohort-specific noise. Adding age, sex and stage
shrinks it by about half (1.18 → 1.28 → 1.21), so much of it runs through stage. The same
pattern is stronger on disease-specific survival throughout — base Cox 11.6% → 14.7% against
9.3% background, OR 1.29 → 1.69 — as a splicing signal specific to the tumour should be
(`08_survival/survival_enrichment.csv`, 24 arm × test × set combinations, all q < 0.05).

**Both sets are above chance.** A nominal p < 0.05 would fire on 5% of tests if nothing were
there, so a background at 9% needs explaining before the 11% means anything. Running the whole
pipeline again with survival shuffled between patients inside each cancer (OS, base model, one
replicate, D-25) gives **5.6%** — 1,148 hits in 20,578 event × cancer pairs. The gates, the
host-gene covariate and the test battery introduce no excess of their own. So the background's
9% is real signal too: splicing events matched to the candidates on coverage and spread are
themselves associated with outcome more often than chance, and the model's events are
associated more often still.

**The candidates inside a model-enriched regulon carry more of it.** A third cut of the same
cells, and the one that ties this section to §4: candidates that sit in one of the six
splicing-factor programmes the models are enriched for — PTBP1, RBM47, HNRNPU, HNRNPF, QKI,
PRPF8 — reach Cox p < 0.05 in 13.3% of their pairs against 10.3% for the candidates that sit in
none (934 against 10,762 pairs, OR 1.33, p = 0.005; `08_survival/regulon_vs_survival.csv`).
Three independent ways of asking which candidates matter — how many models agree, whether the
event recurs, whether it belongs to a coherent factor programme — pick out the same events.

**Some single events do survive multiple testing, and they are almost all kidney.** Applying BH
inside each cancer's own family of ~1,600 Cox tests (a different family from the one
splice-assay prints on its pages, which is BH within one event across cancers), 229 of 11,702
candidate cells clear q < 0.05 in the base arm against 158 of 10,185 background cells — 2.0%
against 1.6%, OR 1.27, p = 0.024; on disease-specific survival 2.6% against 1.6%, OR 1.63,
p = 1.7e-6. But 94% of those cells are KIRC or KIRP, and once age, sex and stage enter the model
the candidate advantage at this threshold disappears entirely (1.17% against 1.15%, OR 1.02,
p = 0.9). The honest reading: the candidates' edge is a broad shift of many small associations,
which survives adjustment at the nominal level, and not an excess of very strong ones, which
does not. Outside the two kidney cohorts only six candidate events clear FDR in the adjusted
model at all — SLC12A4 HIT, SLC39A14 MXE and CYP4F3 MXE in LIHC, NDRG4 SE in LUAD, PHKA1 SE in
STAD and CP A5SS in UCEC.

**Direction is mostly cancer-specific.** Of the 135 model events significant in at least three
cancers, only 30% point the same way in all of them (median consistency 0.67; the background
behaves the same, 22%). An event that shortens survival in one cancer often lengthens it in
another, so a pan-cancer prognostic reading of these events is not supported.

**The rate depends strongly on the cancer.** In the adjusted model it ranges from 25.4% of
tested events in KIRC and 15.2% in KIRP down to 5.0–5.7% in UCEC, THCA and PRAD
(`08_survival/survival_rate_by_cancer.csv`). KICH has too few deaths (9) to fit a Cox model at
all under the tool's gates.

**The events associated in the most cancers.** "Top 10" counts the per-cancer models whose top
ten the event entered — the selection rule of D-22 — and "top 100" the models whose top hundred
it entered, which §5 only tracks for events reaching four or more. "Hits" is the cancers with
Cox p < 0.05 out of those tested, with the number that also clear BH inside their cancer in
brackets, and "worse" the hits with HR > 1. The hazard ratio is the event's strongest cancer,
adjusted model, overall survival.

| Event | Top 10 | Top 100 | In tumour | Hits (FDR) | Worse | Strongest cancer, HR per SD (95% CI) | Enriched regulon |
|---|---|---|---|---|---|---|---|
| RGS1 RI | 1 | < 4 | up | 6/14 (0) | 3 | HNSC 0.77 (0.64–0.92) | PRPF8 |
| ENOSF1 HIT | 6 | 10 | down | 5/14 (0) | 2 | LIHC 0.73 (0.61–0.88) | PTBP1, HNRNPF, HNRNPU |
| COQ6 RI | 1 | < 4 | up | 5/14 (1) | 3 | KIRC 1.26 (1.09–1.45) | – |
| MACROH2A1 MXE | 3 | 9 | up | 5/14 (0) | 1 | KIRC 0.79 (0.67–0.94) | HNRNPF, HNRNPU |
| HSPA12A MXE | 1 | 6 | up | 5/14 (0) | 1 | KIRC 1.25 (1.06–1.46) | – |
| GLIPR1 HIT | 1 | 8 | down | 5/14 (0) | 3 | BRCA 1.37 (1.07–1.74) | – |

The "worse" column makes the direction problem concrete: no event in this list points the same
way in every cancer where it reaches p < 0.05. Of the six, only COQ6 RI clears FDR inside a
cancer (KIRC); the candidates that clear it in two are EPB41 MXE and ADAM15 MXE, in KIRC and
KIRP, both with longer survival at a higher value.

`08_survival/event_dossiers.csv` joins, for every candidate, how many models rank it, how it
moves between tumour and normal, what it does to the protein, whose regulon it belongs to and
its strongest survival association. Assay pages for the leading events — the
tumour-versus-normal comparison, the Kaplan-Meier split and the full Cox model in each of their
three most promising cancers — are in
[figures/tumour_status/assay/](figures/tumour_status/assay/).

**What this does and does not say.** Splicing carries outcome information in these cohorts, and
the events the tumour-status models select carry somewhat more of it than matched events. But
the effect is modest, the individually significant events concentrate in two kidney cohorts, and
directions flip between cancers. Tumour-versus-normal discrimination and prognosis are related but distinct axes, and these
numbers are associations in observational cohorts, not causes and not targets.

## Limitations

- **SHAP explains the models, not causation.** A feature can rank highly because it is a precise,
  deeply covered proxy.
- **RNA quality.** Intronic and long-form signals overlap with RNA quality and nuclear-RNA
  content; TCGA lacks RIN for these samples.
- **Composition.** Tumour versus normal is confounded with stromal, immune and muscle content.
  The regulon and survival analyses adjust for marker scores or host expression, but no
  deconvolution was done.
- **Fusion.** Fused numbers are honest for fixed weights; they are not a nested estimate of a
  method that chooses its weights. The stacked model is a ceiling (D-23).
- **Observational cohorts.** Every survival number is an association in TCGA, with nominal p
  values. None of it establishes cause, and none of it is a therapeutic claim.
- **Reference data.** Protein annotation is GENCODE v45 with InterPro features of one transcript
  per exon (MANE > canonical); gene sets are MSigDB 2026.1.
