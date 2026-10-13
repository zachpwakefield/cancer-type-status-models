# Tumour versus normal tissue: what splicing adds

One classifier per cancer, per modality: 15 TCGA projects with at least 10 normal cases, 9
modalities each, 135 models. Every number here is out-of-fold — each sample was scored by a
model that never saw it ([RESULTS.md](RESULTS.md)).

This document asks what the splicing models learned about becoming a tumour. The other task,
telling the cancers apart, is a different question with a different answer and is written up
separately in [CANCER_TYPE.md](CANCER_TYPE.md).

- **Methods and decisions:** [docs/BIOLOGY_METHODS.md](../docs/BIOLOGY_METHODS.md) (D-13 to D-29).
- **Scripts:** `07`–`19` and `22` write the tables, `20` draws these figures.
- **Figures:** S1–S10 in [figures/tumour_status/](figures/tumour_status/). Each has a caption file
  stating the question, the method, n and the test, and a CSV of the plotted values.
- **Statistics in figures:** raw p. Any q is in the CSV with its FDR family named.

## The argument

Expression already separates tumour from normal tissue almost perfectly, so the question is not
whether splicing can do the same — it is whether splicing is *looking at something else*. Three
findings say it is. The two modalities fail on different samples, so fusing them at fixed
weights removes a fifth of the errors (§1). Expression's signal runs through microenvironment
pathways; the splicing models' genes mostly form no pathway at all (§3), so the splicing signal
is not the expression signal re-encoded. And where the splicing signal does organise, it
organises by *regulator* rather than by pathway: the models' events are concentrated in the
small, distinct programmes of PTBP1, QKI, RBM47 and HNRNPU and depleted from the large shared
nuclear-processing programme (§4).

One result anchors the rest because no model touched it. Using every measured event, with no
SHAP ranking anywhere in the calculation, tumours retain more intronic sequence, use more long
splice-site forms and move their last exons upstream relative to each patient's own normal
tissue (§2). Whatever the models are selecting, they are selecting it from a transcriptome that
really has moved.

What those events look like is consistent with that reading. They recur across cancers far more
than chance allows (§5), they have the protein character of exons a tissue regulates rather than
of damage — in frame, disordered, outside folded domains (§6) — and the models read them
monotonically, as markers rather than as interactions (§7).

The last two sections are where the story is most careful. The events do carry outcome
information, more than matched events do, and more the more models agree on them (§8) — but the
effect is modest, concentrates in two kidney cohorts once multiplicity is enforced inside a
cancer, and reverses direction between cancers. Separating tumour from normal tissue and
predicting outcome are related but distinct axes. §9 then spends ten events' worth of detail
showing that the headline set is not what §6 would lead one to expect: six mutually exclusive
pairs and four terminal exons, not one cassette exon.

§10 is an audit that cuts against the rest, and it changes how some events should be named.
Most TCGA libraries are unstranded and GENCODE nests genes inside one another, so an event named
after one gene can count another gene's reads. The models' top events overlap another gene's
exon more often than measured events do — 44% against 7% for hybrid/terminal (HIT) exons — and
163 of the 1,105 overlapping ones follow the other gene's expression more closely than their
own. Three of §9's ten are among them: EXOSC7's last exon reads CLEC3B, ENOSF1's terminal exon
partly reads its antisense partner TYMS, and PLCG1's exon pair reads ZHX3. Several of the
neighbours are composition markers — smooth-muscle MYH11, stromal CLEC3B — so part of what looks
like splicing is the microenvironment signal of §3 arriving by another route. No model changes,
and the survival result of §8 holds with every overlapping event removed.

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
   for events recurrent in four or more, 1.56 without the clinical covariates), and the gradient
   holds without the two kidney cohorts. The events that survive FDR inside a cancer are almost
   all renal, and directions flip between cancers.
   Separating tumour from normal tissue and predicting outcome are distinct axes.
8. **Ten of those events, in detail.** Six mutually exclusive pairs and four terminal exons,
   not one cassette exon: the macroH2A1.1/1.2 switch, the chimaerin SH2 switch, a 23-amino-acid
   exchange in the collagen V C-terminal domain, and a 4-kb ENOSF1 exon that moves from last to
   first in 14 of 15 cancers and is the most prognostic event in the set. Six of the nine with
   protein coordinates land *inside* an annotated domain — the opposite of the cassette-exon
   picture in point 6, because the most shared events are not cassette exons at all.
9. **Some of the events measure a neighbouring gene.** Overlap with another gene's exon is
   enriched among the models' top events in seven of eight views (44% against 7% for HIT, odds
   6.6 with exon length held fixed), and 163 of 1,105 overlapping top events follow the other
   gene's expression more than their own — EXOSC7 → CLEC3B, ENOSF1 → TYMS, PLCG1 → ZHX3,
   NDE1 → MYH11. macroH2A1, CHN1, COL5A1 and PDLIM7 do not. The survival enrichment of point 7
   is unchanged without any overlapping event (odds 1.25 base, 1.17 adjusted).

---

## 1. Combining the modalities beats expression alone

![S1](figures/tumour_status/S1_modality_gain.png)

<!-- legend S1: copied from figures/tumour_status/S1_modality_gain.caption.md by the figure script's 'legends' step -->
> **Figure S1. Expression and splicing err on different samples, so their fusion beats expression alone.** (a)
> Out-of-fold errors as a share of each cancer's samples (n in brackets) for the expression model (black), the AS
> ensemble (the mean of the 8 AS views' probabilities, blue) and an equal-weight fusion of the two (dark blue). The
> fusion weights were fixed before any result was seen (D-14). Stars: exact McNemar, fusion against GEX, on paired
> per-sample correctness (* p < 0.05, ** p < 0.01, *** p < 0.001); q (BH over the 45 tests) is in the CSV.
>
> (b) Errors pooled over all cancers. 'stacked' is a cross-fitted logistic regression on the 9 views' probabilities
> (D-23): it is shown as a ceiling, not as the headline, because the base models behind its training probabilities
> did see the samples it scores.
>
> Pooled over the cohort, fusion corrects 27 of the expression model's 99 errors and introduces 5 (exact McNemar p =
> 1.1e-4). No single cancer has enough errors to reach significance on its own; the per-cancer tests are in the CSV.
<!-- /legend -->

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

<!-- legend S2: copied from figures/tumour_status/S2_global_shifts.caption.md by the figure script's 'legends' step -->
> **Figure S2. Tumours retain more intronic sequence and move their last exons upstream, relative to their own
> normal tissue.** Each cell is the median over matched patients (n in brackets) of the patient's mean paired tumour
> − normal change, averaged over every event measured in that cancer, so there is no SHAP selection and no
> circularity. The columns say what an increase means:
> - RI: more intron retention.
> - A5SS/A3SS: more use of the long exon form, which includes extra intronic sequence.
> - SE: more cassette inclusion.
> - MXE: more of the upstream exon.
> - ALE/AFE: the used last or first exons move downstream (a 0-1 position index along transcription).
>
> Blue is higher in the tumour and orange higher in the normal tissue, the house tumour/normal colours. Stars:
> two-sided Wilcoxon signed-rank test across patients, raw p (* < 0.05, ** < 0.01, *** < 0.001). BH q (family: 15
> cancers × 7 views) is in the CSV.
>
> Retention and long splice-site forms rise in most cancers. Last exons move upstream and never downstream. BRCA is
> the exception for RI. Caveat: RNA integrity and nuclear-RNA content also raise intronic reads, and RIN is not
> available (D-15).
<!-- /legend -->

This uses **every measured event**, not the models' selections, so the models cannot have biased
it: in matched tumour-normal pairs from the same patient, does RNA processing move in a common
direction?

**Retained introns.**
- 7 cancers retain more (BLCA, COAD, KICH, KIRC, KIRP, LUSC, PRAD); BRCA retains less.
- Retained introns that gain ≥ 5 PSI points in tumours outnumber those that lose them about 3:1
  (median 15% vs 4.5% of events), in 13 of 15 cancers.
- The shift is carried by most patients: in KICH, BLCA, PRAD and KIRC, 70–91% of patients retain
  more in their tumour than in their own normal tissue.
- This is the pattern Dvinge and Bradley (2015) reported across 16 TCGA cancer types with matched
  normals, including its exception: retention rose in tumours in almost every type, while in
  breast it was enriched in the *normal* tissue. BRCA is the one cancer that retains less here.

**Alternative splice sites.** A5SS long forms rise in 10 cancers and fall in none; A3SS long
forms rise in 5.

**Last exons.** They move upstream in 7 cancers and downstream in none — the direction of the
proximal polyadenylation widely reported in tumours (Mayr and Bartel 2009; across TCGA, Xia et
al. 2014). Those studies measured tandem 3′ UTRs, which this one does not; an upstream
alternative last exon is the same choice made between exons rather than within one.

**Alternative explanations.**
- *Proliferation.* Within tumours, per-sample retention does not track a 10-gene proliferation
  score (Spearman −0.25 to 0.19, no consistent sign).
- *3′-end factors.* How much a cancer shifts is not explained by how much its cleavage and
  polyadenylation factors rise (|ρ| ≤ 0.39 across 15 cancers).
- *RNA quality.* RNA integrity and nuclear-RNA content also raise intronic reads. TCGA has no RIN
  for these samples, so a technical component cannot be excluded.

## 3. Expression reads the microenvironment; splicing does not form pathways

![S3](figures/tumour_status/S3_pathways.png)

<!-- legend S3: copied from figures/tumour_status/S3_pathways.caption.md by the figure script's 'legends' step -->
> **Figure S3. Expression models draw on stromal programmes; the splicing models' genes rarely form pathways.** For
> every cancer (columns) and view, the 200 genes with the largest SHAP importance were tested for enrichment of the
> MSigDB Hallmark sets against all genes measurable in that view and cohort. The test is GOseq/Wallenius, correcting
> for the number of events per gene (D-16).
> - A dot marks a significant set in at least one AS view; its size is the number of AS views (of 8).
> - A ring marks significance in the GEX model.
> - Significance is q < 0.05, BH within the Hallmark collection of each model.
>
> GEX models consistently use stromal and microenvironment programmes (EMT, coagulation, KRAS, myogenesis). The AS
> models reach the same programmes far less often, mostly through EMT and myogenesis genes. Their top genes form any
> enriched set in at most 9 of 15 cancers (SE) and in none for RI, so the tumour splicing signal is largely not
> organised by pathway membership.
<!-- /legend -->

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

<!-- legend S4: copied from figures/tumour_status/S4_regulons.caption.md by the figure script's 'legends' step -->
> **Figure S4. The models' events sit in small, specific splicing-factor programmes, not in the large shared one.**
> A factor's *core regulon* is the events whose PSI tracks its expression within tumours — partial Spearman |rho| >=
> 0.3 at p < 0.05 after rank-residualising both on a stromal score, an immune score and the sample's own mean PSI in
> that view — reproducibly, with the same sign, in at least 5 of the 15 cancers (D-24).
>
> (a) For every factor with at least 15 marked events, the odds that one of the models' top events is in its core
> regulon, against the coverage- and spread-matched background of D-22. Only factors with q < 0.05 (BH over the 32
> factors tested) are shown; the regulon's size is in brackets. Odds ratios use a 0.5 correction so that an empty
> cell stays finite — PTBP1's 15 candidate events against 0 background events is otherwise infinite.
>
> (b) Every factor with at least 20 core events. Large regulons are not separate programmes: they overlap each other
> heavily (Jaccard up to 0.7), and the 1,074 distinct events carry 7,529 memberships, about 7 factors each.
>
> The two panels say the same thing. The models' events avoid the big shared programme of nuclear RNA-processing
> factors (SRSF6, SRSF4, RBM4, SRSF10, odds 0.26-0.50) and concentrate in the small, distinct regulons of tissue and
> differentiation regulators (PRPF8, QKI, HNRNPF, HNRNPU, RBM47, PTBP1).
>
> These are correlations inside tumours: a factor may be co-regulated with its apparent targets rather than drive
> them. `regulon_core_events.csv` lists the members, `global_axis_factors.csv` the factors most tied to the
> per-sample PSI level that was residualised out.
<!-- /legend -->

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
composition: normal tissue carries more stroma and smooth muscle. QKI's fall matches the
literature, which finds it among the splicing factors most often lost in lung adenocarcinoma
(Zong et al. 2014) and reduced across the cancers that lose macroH2A1.1 (Novikov et al. 2011);
§9 tests the second of those links on the models' own macroH2A1 event.

These are correlations inside tumours. A factor may be co-regulated with its apparent targets
rather than drive them.

## 5. Different cancers rely on the same events

![S5](figures/tumour_status/S5_recurrent_events.png)

<!-- legend S5: copied from figures/tumour_status/S5_recurrent_events.caption.md by the figure script's 'legends' step -->
> **Figure S5. Different cancers rely on the same events, and those events move the same way across cancers.**
> Across the 8 AS views, 404 events rank in the top 100 (mean |SHAP|) of at least 4 of the 15 tumour-vs-normal
> models. 3.1 would be expected if the rankings were independent (Poisson-binomial over each event's measured
> cohorts; `recurrence_vs_chance.csv`).
>
> The figure shows the 40 most recurrent events. Rows give the host gene (HGNC), view and coordinates (hg38,
> 1-based). Colour is the median paired tumour − normal change in every cancer, including cancers where the event
> did not rank (so the pattern across cancers is not selected on): blue where the event is higher in the tumour,
> orange where it is higher in the normal tissue, clipped at ±30. Dots mark the cancers whose top 100 contain the
> event.
>
> Recurrent events tend to move the same way across cancers. For the median event, 73% of the cancers where it is
> measured share its majority direction; 46% of events reach 80%.
>
> Several genes appear as switch pairs, with two exons of the same gene moving in opposite directions (EXOSC7 last
> exons, NDE1 first exons, TPM1 mutually exclusive exons).
>
> Recurrent RI and A5SS events mostly rise in tumours (73% and 77%), matching the global shift in S2. The full list
> with direction consistency is in `recurrent_events.csv`.
<!-- /legend -->

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

**A caution, now measured (§10).** 107 of the 404 recurrent events overlap an exon of another
gene and 30 of them follow that gene's expression more than their own host's. Two of the pairs
above are among them. The EXOSC7 last exon that falls in tumours follows CLEC3B, a stromal gene
nested in EXOSC7's 3′ end on the same strand (ρ +0.76, against −0.28 with EXOSC7), and the NDE1
first exon that falls follows MYH11, smooth-muscle myosin (+0.80 against −0.34): both report how
much stroma or muscle a sample holds, which normal tissue has more of. The TPM1 pair overlaps
nothing. The NR1D1 and THRA HIT events overlap at chr17:40.09 Mb on opposite strands and move in
opposite directions in 14 of 15 cancers; the THRA event follows NR1D1's expression (−0.62,
against +0.42 with THRA), while the NR1D1 event follows its own gene (+0.59) more than THRA
(−0.41). The pair is one NR1D1 signal read twice, not two events.

## 6. The exons the models use look tissue-regulated

![S6](figures/tumour_status/S6_exon_consequences.png)

<!-- legend S6: copied from figures/tumour_status/S6_exon_consequences.caption.md by the figure script's 'legends' step -->
> **Figure S6. The models' cassette exons keep the reading frame and avoid folded domains.** Cassette exons in the
> top 100 of any tumour-vs-normal model (blue) or any class of the cancer-type model (dark blue), compared with a
> random 6,000 of all measured exons (grey).
> - Length divisible by 3 is computed from the coordinates.
> - The remaining features need an exact GENCODE v45 exon. Codons are taken from the best transcript (MANE >
>   canonical > other) and overlapped with InterPro Pfam/CDD domains, MobiDB-lite disordered regions and TMHMM
>   helices (D-17).
>
> p: two-sided Fisher exact test against all measured exons, raw p. BH q is in the CSV, with its family. Exons
> skipped or included more in tumours are in `cassette_exon_protein_features.csv`.
<!-- /legend -->

The top cassette exons of the tumour-status models (n = 1,281) and of the cancer-type classes
(n = 1,914), against 6,000 random measured exons:

- more often of a length divisible by 3: 51% and 54%, against 41% (p = 1.3e-10 and 3.1e-24);
- among fully coding exons, less often overlapping a Pfam/CDD domain: 37% and 35%, against 51%;
- among the same exons, more often overlapping a disordered region: 21% and 20%, against 14%;
- less often overlapping a transmembrane helix: 3.7% and 4.6%, against 6.7%.

So the signal sits in exons that change protein segments in frame, mostly in disordered,
interaction-prone regions, rather than breaking folded domains or shifting the reading frame.
These are the properties of exons that tissues regulate: tissue-specific exons are enriched in
disordered segments carrying binding motifs, and switching them rewires protein interactions
(Buljan et al. 2012; Ellis et al. 2012). The frame bias has also been seen in tumours before.
Tsai et al. (2015) found 53% of cancer-specific cassette exons across TCGA to be
frame-preserving, against 42% of controls, close to the 51% against 41% here.

**Direction.** Exons skipped more in tumours are the most frame-preserving (55% divisible by 3)
and the most domain-depleted (30%, OR 0.42); exons included more in tumours are the most enriched
for disorder (24%, OR 1.98).

## 7. The models read these events as markers, not as interactions

![S7](figures/tumour_status/S7_event_portraits.png)

<!-- legend S7: copied from figures/tumour_status/S7_event_portraits.caption.md by the figure script's 'legends' step -->
> **Figure S7. The most shared events shift the typical tumour, and the models read them as markers.** The six most
> recurrent events (S5) whose paired change has the same sign in at least 80% of the cancers where measured, one per
> gene.
>
> (a) The median value in matched normal (orange) and tumour (blue) samples of each cancer. Cancers whose top 100
> contains the event are in black, the others in grey. The heading gives the host gene, view, coordinates and, for
> cassette exons, the protein consequence (GENCODE v45; D-17).
>
> (b) SHAP dependence in the cancer where the event ranks highest. Each held-out sample's event value is plotted
> against its SHAP value on the tumour log-odds. The value is the PSI; missing values were median-imputed by the
> model.
>
> The CSV gives the per-cancer medians and paired Wilcoxon p (raw). Per-sample values stay local (D-12).
<!-- /legend -->

Six events are shown: the most recurrent, one per gene, whose paired change has the same sign in
at least 80% of cancers. (a) Each cancer's matched normal and tumour medians; (b) the SHAP
dependence in the cancer where the event ranks highest.

**The medians say the shift moves the typical tumour, not a subset of them.** Each of the six
moves the same way in 14 or 15 of the 15 cancers, and the median shift is large for four of them
(|ΔPSI| 0.17–0.46). There is one substantial reversal in the whole panel: GCOM1 rises by 0.34 in
KICH against a median fall of 0.24 elsewhere — the single off-colour cell of its row in S5.
These are not events that pick out an aggressive subset; they are events whose typical value
differs between a patient's tumour and that patient's own normal tissue.

**The dependence plots say the models use them as levels.** Value and SHAP are monotone in all
six panels (Spearman 0.75–0.93 in five; 0.52 in EXOSC7, where 61% of the samples sit at exactly
one value and the correlation is bounded by the ties). Within any fifth of the value range,
tumour and normal samples get almost the same push — the two groups' mean SHAP differs by at
most 0.05–0.18 against total SHAP ranges of 0.43–0.70 — so the model is reading the event's
level, not the label it happens to come with.

That is also the mechanical reason §1 came out the way it did. If each view contributes a nearly
additive, nearly monotone piece of evidence, a fixed 1:1 average is close to the right
combination rule and a stacker with freedom to learn weights has little left to find, which is
what the error counts show (77 against 76, with a worse mean AUROC for the stacker).

## 8. Do these events carry information about outcome?

![S8](figures/tumour_status/S8_survival.png)

<!-- legend S8: copied from figures/tumour_status/S8_survival.caption.md by the figure script's 'legends' step -->
> **Figure S8. The models' events are modestly enriched for survival associations, more so the more models agree.**
> Every event was tested with splice-assay 0.4.0 in each of the 15 cancers (6,152 tumours, 1,620 deaths): Cox on the
> event's value with the host gene's own expression in the model, in the pre-specified base form and again with age,
> sex and stage added (D-25). The models' top events are compared with the coverage- and spread-matched background
> of D-22, so this is not a comparison of which events are easier to test.
>
> (a) The share of testable event × cancer pairs with a nominal Cox p < 0.05, with 95% Clopper-Pearson intervals;
> the three candidate sets are nested and ordered by how many models agree on the event. Odds ratios against the
> background of the same arm, base / adjusted: every candidate 1.27 / 1.18, recurrent 1.56 / 1.28, top 10 for two or
> more cancers 1.62 / 1.21 (Fisher exact, every set q < 0.05). The rate rises along that order.
>
> The dashed line is the same pipeline run over outcomes shuffled between patients inside each cancer: 5.6% (1,148
> of 20,578 pairs), so both sets sit above what the gates, the host-gene covariate and the test battery produce on
> their own.
>
> Candidates inside one of the six regulons the models are enriched for (S4) reach p < 0.05 in 13.3% of their pairs
> against 10.3% for the rest (OR 1.33, p = 0.005), so the three ways of ranking a candidate -- how many models
> agree, whether it recurs, whether it sits in a coherent factor programme -- pick out the same events.
>
> (b) The 12 candidate events associated in the most cancers, each in its strongest cancer: hazard ratio per SD of
> PSI (per SD of the HIT index for HIT events) with its 95% CI, and (cancers with p < 0.05 / cancers tested). Red: a
> higher value goes with worse survival.
>
> Nominal p values from observational cohorts. BH inside each cancer's own family of ~1,600 tests (`cox_q_cohort`)
> leaves 1.2% of candidate cells significant in the adjusted arm, no more than the background's 1.2%, and 94% of
> them are KIRC or KIRP: the set-level excess is a broad shift of small associations, not an excess of strong ones.
> `event_survival_cells.csv.gz` holds every statistic with both BH families (keyed by `event_index.csv`), and
> `event_dossiers.csv` joins each event's recurrence, protein consequence, regulon membership and survival. PSI can
> reflect tumour purity or composition, so these are associations, not causes, and not therapeutic claims.
<!-- /legend -->

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

**Not a neighbour's expression.** Some of the candidates read another gene's expression (§10).
Removing every event whose segment overlaps another gene's exon — from the candidates and the
background alike — leaves the result where it was: 10.9% against 8.9%, OR 1.25 (p = 4.5e-6) in
the base arm and 1.17 (p = 0.001) adjusted, against 1.27 and 1.18 with them
(`12_overlap_audit/survival_enrichment_without_overlaps.csv`).

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

**The set-level result is not a kidney result.** Pooling again without KIRC and KIRP (D-29), the
rates fall — 8.6% of candidate pairs against 7.3% of background in the base arm — but the odds
and the gradient stay: 1.19 for every candidate, 1.48 for recurrent events and 1.58 for the top
10 of two or more cancers (adjusted 1.19, 1.30, 1.24; every row p < 0.05). Leaving out any one
cancer instead keeps the base-arm odds for every candidate between 1.23 and 1.32
(`14_survival_robustness/`). So the kidney cohorts carry the strong single associations, and
the broad shift is spread across the cancers.

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
way in every cancer where it reaches p < 0.05. Of the six, only ENOSF1 follows a neighbouring
gene by the rule of §10 (TYMS); COQ6 RI sits just under it (ENTPD5, ρ 0.499 against 0.38). Of
the six, only COQ6 RI clears FDR inside a cancer (KIRC); the candidates that clear it in two are
EPB41 MXE and ADAM15 MXE, in KIRC and KIRP, both with longer survival at a higher value.

`08_survival/event_dossiers.csv` joins, for every candidate, how many models rank it, how it
moves between tumour and normal, what it does to the protein, whose regulon it belongs to and
its strongest survival association. Assay pages for the leading events — the
tumour-versus-normal comparison, the Kaplan-Meier split and the full Cox model in each of their
three most promising cancers — are in
[figures/tumour_status/assay/](figures/tumour_status/assay/).

**What this does and does not say.** Splicing carries outcome information in these cohorts, and
the events the tumour-status models select carry somewhat more of it than matched events. But
the effect is modest, the individually significant events concentrate in two kidney cohorts, and
directions flip between cancers. Tumour-versus-normal discrimination and prognosis are related
but distinct axes, and these numbers are associations in observational cohorts, not causes and
not targets.

## 9. The events one by one

![S9](figures/tumour_status/S9_event_protein.png)

<!-- legend S9: copied from figures/tumour_status/S9_event_protein.caption.md by the figure script's 'legends' step -->
> **Figure S9. The most shared events exchange protein segments through mutually exclusive and terminal exons.** The
> 9 recurrent events of S5 whose exon belongs to a protein-coding transcript (of ten, by number of models;
> `11_event_cards/tumour_status_cards.csv` holds all ten with the rest of what is known about them). Each bar is one
> protein, scaled to its own length (printed at the right), from its N terminus to its C terminus: grey blocks are
> Pfam/CDD domains of the chosen transcript, the line underneath marks disordered regions (MobiDB-lite), and the
> coloured block is the segment the event adds or removes -- blue where the segment is higher in tumours, orange
> where it is lower (D-17, GENCODE v45). Three of the transcripts have no InterPro annotation at all, which is said
> so in the margin: a bare bar there means the reference is silent, not that the segment misses every domain.
>
> The annotation is what the swap does to the reading frame. For a mutually exclusive pair that is the difference
> between the two exons, not either length, because the swap replaces one with the other. A first or last exon
> replaces one end of the protein rather than shifting a frame, so it is annotated by the length of the segment it
> contributes instead. A retained intron has no protein coordinates and is not drawn.
>
> These are annotations of one transcript per gene, not isoform quantification: the models see the event, not the
> protein, and no protein was measured here.
<!-- /legend -->

Everything above is about sets of events. This section takes ten of them — the recurrent events
of §5, most models first, one per gene — and says what the pipeline knows about each: where the
alternative segment falls in the protein, how it moves between tumour and normal, how the models
read it, whose programme it belongs to, and whether it tracks outcome (D-26). Every number is
from `11_event_cards/tumour_status_cards.csv`, which holds the rest of each card. Where a
sentence names a known isoform, that is context from the literature and not a result here.

Two conventions matter. For a mutually exclusive pair the length that counts is the *difference*
between the two exons, because the swap replaces one with the other: COL5A1's two exons are both
69 nt, so that swap exchanges 23 amino acids and changes nothing else. And for a first or last
exon the reading frame is not the question — the event replaces one end of the protein.

**GCOM1 · mutually exclusive exons · chr15:57,674,968-57,675,068 · down in tumour.** The most
widely shared event in the study: in the top 100 of 14 of the 15 models, ranked first in one of
them, median rank 8. The tracked exon is 101 nt against a 205-nt partner, so the two forms
differ by 104 nt and by the frame downstream, and the segment lands inside the RNA polymerase II
subunit domain (aa 402–435 of 765). It falls in tumours in 14 of 15 cancers (median ΔPSI
−0.24), most steeply in LUSC (0.73 in matched normal tissue to 0.10 in tumour), and in 13 of
the 14 models a higher value argues for normal (median Spearman between value and SHAP −0.87). It
belongs to the HNRNPU regulon (ρ −0.37) and reaches Cox p < 0.05 in 2 of its 14 testable cancers
(BLCA, HR 0.71 per SD). The name needs stating (§10): GCOM1 is a GENCODE readthrough model
spanning MYZAP and POLR2M. The tracked exon is a MYZAP exon and its partner a POLR2M exon, both
joined to a MYZAP exon upstream and a POLR2M exon downstream, so the event weighs one gene's
exon against the other's inside readthrough transcripts, and how much MYZAP–POLR2M readthrough
a sample makes is part of what it measures. The domain the segment lands in is POLR2M's: Gdown1,
an RNA polymerase II-associated protein.

**PLCG1 · mutually exclusive exons · chr20:41,174,124-41,174,311 · up in tumour.** 188 nt against
161 nt: a 27-nt, 9-amino-acid, in-frame difference at aa 1,216–1,278 of 1,291, inside the last
6% of phospholipase C-γ1. It is one of two events here that move the same way in all 15 cancers
(median +0.11, KIRC 0.50 → 0.80), and all 13 models that rank it read a higher value as tumour
(+0.91). Two of 14 cancers reach p < 0.05 (BLCA, HR 0.75). The reference carries no InterPro
features for this transcript, so the pipeline cannot say what the segment sits in. And it may
not be a PLCG1 choice at all (§10): the partner exon overlaps an exon of ZHX3, which runs the
other way across PLCG1's 3′ end, and within cancers the event follows ZHX3's expression
(ρ −0.60) more than PLCG1's (+0.20). More ZHX3 reads on the partner exon lower the PSI, so the
rise in tumours is at least partly ZHX3 falling.

**PDLIM7 · mutually exclusive exons · chr5:177,491,807-177,491,976 · up in tumour.** The largest
in-frame change here: a 170-nt exon against a 17-nt partner, 153 nt and 51 amino acids, at
aa 94–150 of a 195-residue isoform (ENST00000393546). Up in all 15 cancers, with the single
biggest shift of any event in this list — COAD moves from 0.25 to 0.88 — and all 12 models read
it the same way. It is also the one event here whose higher values go with *worse* survival
where it is strongest: KIRC HR 1.49 per SD, significant in 3 of 14 cancers.

**EXOSC7 · alternative last exon · chr3:45,035,524-45,036,066 · down in tumour.** A 543-nt last
exon that belongs to no protein-coding transcript of the gene in GENCODE v45, so the pipeline
can place it but not translate it; it is the one card of the ten that S9 cannot draw. It is
still among the cleanest signals in the study: down in 14 of 15 cancers, median ΔPSI −0.44, and
in BLCA it goes from 0.93 in normal tissue to 0.00 in tumour — switched off. All 10 models that
rank it read a higher value as normal (−0.89), and 2 of 14 cancers reach p < 0.05 (HNSC, HR
0.77). §5 notes that EXOSC7 also contributes a second last exon that moves the other way. The
reason it belongs to no coding transcript of EXOSC7 is that it is CLEC3B's last exon: CLEC3B,
which encodes the stromal and adipose protein tetranectin, is nested in EXOSC7's 3′ end on the
same strand, and the event follows its expression (ρ +0.76, against −0.28 with EXOSC7). What
switches off in tumours is CLEC3B, and the opposite-moving EXOSC7 exon is the other side of the
same ratio.

**ENOSF1 · terminal exon · chr18:670,318-674,406 · down in tumour.** The HIT index of this 4.1-kb
exon measures whether it is used as a first or a last exon, and in tumours it moves towards
first-exon use in 14 of 15 cancers (median −0.46, UCEC 0.68 → −0.62). Used as a last exon it
contributes the final 33 residues of a 443-aa protein. It is the most connected event here — in
13 factor regulons, among them PTBP1 (ρ −0.41), HNRNPU (−0.41) and HNRNPF (−0.38), all three
programmes the models' events are enriched in (§4) — and the most prognostic, with Cox p < 0.05
in 5 of 14 cancers, strongest in LIHC (HR 0.73, 0.61–0.88). Its assay page is in
[figures/tumour_status/assay/](figures/tumour_status/assay/). ENOSF1 was first cloned as rTS,
an RNA antisense to thymidylate synthase (Dolnick 1993), and this exon overlaps TYMS on the other
strand by 719 bp. Within tumours the event follows TYMS (ρ −0.62, against +0.34 with ENOSF1).
TYMS is a proliferation gene, so the step that matters is the partial correlation (D-28): with a
ten-gene proliferation score held fixed the event still follows TYMS in all 15 cancers (median
ρ −0.35), and with TYMS held fixed it still follows proliferation in all 15 (−0.30). It reads
two things at once — its antisense partner and the cell cycle — and the survival association is
best read as theirs.

**NOP16 · alternative last exon · chr5:176,383,946-176,384,374 · up in tumour.** 429 nt
contributing aa 132–178 of 178 — the whole C-terminal quarter, inside the Nop16
ribosome-biogenesis domain. The median shift is small (+0.05) because the event is near
all-or-nothing in only some cancers: LUAD goes from 0.00 in normal tissue to 1.00 in tumour, and
11 of 15 cancers share the direction. Two of 14 cancers reach p < 0.05 (HNSC, HR 1.29).

**NR1D1 · terminal exon · chr17:40,092,793-40,093,282 · up in tumour.** The last 66 residues of a
614-aa protein, inside the nuclear-hormone-receptor ligand-binding domain. It shifts towards
last-exon use in 14 of 15 cancers (median +0.29, THCA −0.82 → −0.09). It is the one card of
the ten that was never tested for survival, because it entered no model's top ten and so was not a
candidate under D-22. It is also the cautionary one: as §5 says, the NR1D1 and THRA terminal
exons overlap at chr17:40.09 Mb on opposite strands and most TCGA libraries are unstranded, so
the pair is one observation, not two. §10 says which: the THRA event follows NR1D1's
expression, while this one follows NR1D1 itself (ρ +0.59) more than THRA (−0.41), so the signal
is NR1D1's.

**COL5A1 · mutually exclusive exons · chr9:134,829,976-134,830,044 · up in tumour.** Two 69-nt
exons, so the swap is exactly length-neutral and in frame: 23 amino acids exchanged at
aa 1,690–1,712 of 1,838, inside the fibrillar-collagen C-terminal domain. Up in 11 of 15 cancers
(median +0.11, ESCA 0.63 → 0.87), and all nine models that rank it read a higher value as tumour
(+0.90). One of 14 cancers reaches p < 0.05 (LUSC, HR 0.85).

**CHN1 · mutually exclusive exons · chr2:174,877,840-174,878,128 · down in tumour.** A 289-nt
exon against a 78-nt partner — 211 nt, frame-shifting — inside the SH2 domain at aa 87–183 of
459. This is the switch that distinguishes the two chimaerin isoforms, only one of which carries
the SH2 domain. Tumours use less of the SH2-containing form in 13 of 15 cancers (median −0.11,
HNSC 0.38 → 0.15), 8 of 9 models read a higher value as normal (−0.75), and 2 of 14 cancers
reach p < 0.05 (KIRC, HR 0.83).

**MACROH2A1 · mutually exclusive exons · chr5:135,352,946-135,353,045 · up in tumour.** 100 nt
against 91 nt: a 9-nt, 3-amino-acid, in-frame difference inside the macro domain (aa 197–230 of
372). These are the two exons that distinguish the macroH2A1.1 and macroH2A1.2 histone variants.
Tumours shift towards the tracked exon in 10 of 15 cancers (median +0.09, UCEC 0.56 → 0.90) —
the least consistent direction of the ten, which the S5 heat map shows row by row. It sits in
the HNRNPF (ρ +0.37) and HNRNPU (+0.32) regulons, both enriched for model events, and reaches p
< 0.05 in 5 of 14 cancers (KIRC, HR 0.79), with an assay page. The tracked exon is
macroH2A1.2's; the 369-aa partner form is macroH2A1.1, whose macro domain binds ADP-ribose
metabolites, which macroH2A1.2's cannot (Kustatscher et al. 2005). That makes two published
predictions testable here (D-28). QKI promotes the macroH2A1.1 exon (Novikov et al. 2011), and
within tumours the macroH2A1.2 exon falls as QKI rises in 12 of 15 cancers (median ρ −0.23; sign
test p = 0.018), 13 of 15 with proliferation held fixed (−0.32, p = 0.004). Proliferating
tumours lose macroH2A1.1 (Novikov 2011; Sporn et al. 2009), and the exon rises with a
proliferation score in 13 of 15 (+0.26, p = 0.004). The five cancers that do not shift towards
it — KIRC, KIRP, KICH, THCA, HNSC — are mostly the slow-growing ones, and across the 15 the size
of the shift follows how much more tumours proliferate than their normals (ρ 0.71; exploratory,
added after seeing which cancers did not shift). An exon of PITX1-AS1 overlaps the partner exon
on the other strand, but the event does not follow it (ρ +0.10). Of the ten, this is the one
whose movement the literature predicts in both direction and regulator.

**What the ten have in common.** Not one is a cassette exon: six are mutually exclusive pairs
and four are terminal or last exons. That is §5's view-level pattern at its sharpest — SE
contributes 7 of the 404 recurrent events — and it means the protein statistics of §6, which are
about cassette exons, do not transfer here. Six of the nine with protein coordinates fall inside
an annotated domain, and the other three sit on transcripts the feature reference does not
annotate at all — so not one of the nine is known to miss every domain, which is the opposite of
the §6 picture. Only two shift the reading frame, both by swapping an exon for a much longer
partner. All nine that were tested reach nominal Cox significance in at least one of their 14
cancers, where the permuted-outcome rate of §8 would put about five of nine at that mark; the
two that reach it in five cancers, ENOSF1 and MACROH2A1, are both in regulons the models' events
are enriched for. After §10, the ten sort into three groups: five are what they say they are and
overlap nothing that moves them (PDLIM7, NOP16, COL5A1, CHN1, MACROH2A1); two are real splicing
events whose name or neighbour needs stating (GCOM1, a MYZAP–POLR2M readthrough choice, and
NR1D1, which is NR1D1's but is read twice); and three follow another gene's expression (EXOSC7 →
CLEC3B, PLCG1 → ZHX3, ENOSF1 → TYMS, the last also following proliferation).

## 10. Which events measure another gene

![S10](figures/tumour_status/S10_overlap_audit.png)

<!-- legend S10: copied from figures/tumour_status/S10_overlap_audit.caption.md by the figure script's 'legends' step -->
> **Figure S10. Some of the models' events, terminal exons above all, read a neighbouring gene's expression.** TCGA
> mRNA libraries are mostly unstranded, and GENCODE readthrough and nested models put exons of one gene inside
> another, so an event named after one gene can count reads of a neighbour (D-27). (a) The share of events whose
> segment (the alternative exon, both exons of an MXE pair, the retained intron or the splice-site extension)
> overlaps an exon of another GENCODE v45 gene by >= 50 bp, on either strand, for every measured event (grey; 6.8%
> over all 8 views) and for the models' top events (top 100 of each cancer's model, per view; dark). The gap is
> widest for HIT: 44% against 7%. Odds ratios with the queried length held fixed (Mantel-Haenszel over length
> deciles, because long terminal exons overlap something more often by chance) are 1.2 to 6.6 across views; every
> statistic is in `S10_overlap_audit.csv`.
>
> (b) Of the top events with such an overlap, how many *track* the other gene — Spearman with its log2 CPM of at
> least 0.5 in size and larger than with the host gene's own (within each cancer, median over the cancers whose top
> 100 holds the event); the criterion was fixed before any result was seen. 163 of 1,105 do. Those are better read
> as a readout of the neighbour's expression than as a splicing choice of the host, whether the coupling is
> misassigned reads or genuine co-transcription; the two cannot be told apart without stranded libraries. Among the
> events the narrative names: PLCG1 MXE follows ZHX3 (ρ −0.60, own gene +0.20); EXOSC7 ALE follows CLEC3B (ρ +0.76,
> own gene −0.28); ENOSF1 HIT follows TYMS (ρ −0.62, own gene +0.34).
>
> No model was refitted: the audit changes how events are described, not what the models are.
<!-- /legend -->

An event is meant to be a statement about one gene's RNA. Two features of this data can make it
a statement about another gene instead (D-27). Most TCGA mRNA libraries are unstranded, so where
an exon of the host overlaps an exon of a gene on the other strand, reads of either count towards
the event. And GENCODE's readthrough and nested models put exons of one gene inside another's, so
an "alternative" first, last or mutually exclusive exon of the host can be its neighbour's exon,
whose use then measures how much of the neighbour is made.

**How often.** Of the 357,832 measured events, 6.8% have a segment that overlaps an exon of
another gene by at least 50 bp (5.6% antisense, 1.2% same strand); another 5.5% overlap only a
readthrough model of their own gene, which is a naming matter. The models' top events overlap
more often in seven of the eight views, and in HIT by far the most: 44% of the top HIT events,
against 7% of all measured HIT events (odds 10.1). Long terminal exons overlap something more
often by chance, and the top HIT exons are longer (median 585 against 161 nt), but with length
held fixed the odds are still 6.6. The other views sit between 1.2 and 1.8 once length is
matched, and MXE is not enriched at all (q = 0.28).

**How often it matters.** Overlap is not the same as contamination: an overlapping event is a
problem only if it follows the other gene. The rule, fixed before any result was seen, is that an
event *tracks* its neighbour when, within the cancers whose model ranks it, its Spearman
correlation with the neighbour's expression is at least 0.5 in size and larger than with its own
gene's. 163 of the 1,105 overlapping top events do (15%); in HIT, 124 of 450 (28%). A quarter
of the overlapping events cannot be checked because the neighbour is not in the stored
expression matrix (mostly lncRNAs and small RNAs — snoRNAs, miRNAs — which it does not hold).

**What they read.** The neighbours that do get tracked are often not splicing biology at all
but the composition of the sample: smooth-muscle MYH11 (NDE1's first exon), stromal CLEC3B
(EXOSC7's last exon), the liver and plasma protein CP (an HPS3 terminal exon), proliferation's
TYMS (ENOSF1). §3 found that expression separates tumour from normal tissue through the
microenvironment; these events carry part of that same signal, in splicing's clothes.

**What does not change.** No model is refitted, and nothing above removes an event from one; the
models were trained on these measurements, and they classify as well as they do. The global
shifts of §2 use every measured event and cannot be driven by a few percent of them. The
survival enrichment of §8 is the same without any overlapping event (OR 1.25 against 1.27, base
arm). And whether a tracked event is misassigned reads or genuine co-regulation of two
neighbouring genes — convergent transcription is itself regulated — cannot be told apart without
stranded libraries. What changes is how the events are named, which is what §9 now does.

## What would change the reading

The limitations below are the standing caveats. These are the specific observations that would
overturn the argument, and how each could be checked:

- **An RIN or nuclear-RNA measurement.** §2's intron-retention and long-form shifts are the one
  global result with a plausible technical twin. If retention tracked RNA integrity within
  tumours in a cohort that has RIN, §2 would become a quality finding rather than a biological
  one. §4–§9 would survive it: they are about which events the models pick, against a background
  matched on coverage and spread.
- **A deconvolution.** Tumour and matched normal tissue differ in stromal, immune and muscle
  content, and §4 already shows that one global axis — the sample's own mean PSI — manufactures
  regulons when it is not residualised out. If the enrichments of §4 vanished under cell-type
  deconvolution rather than marker scores, the "small, specific programmes" reading would go
  with them.
- **A second cohort.** The recurrence of §5 is across cancers inside TCGA, with one processing
  pipeline. Recovering the same events in an independent cohort processed differently would
  separate biology from pipeline; failing to would not be decisive, but would bound the claim.
- **The kidney cohorts.** Drop KIRC and KIRP and §8's FDR-level result is gone. The set-level
  gradient was recomputed without them and holds (§8, D-29), so what remains open is the
  strong single associations, not the broad shift; a renal cohort outside TCGA would test them.
- **Strandedness.** §10 identifies the events whose overlap with another gene could be read
  misassignment, and tracks which follow the neighbour, but in unstranded data it cannot tell
  misassigned reads from co-regulated neighbours. A stranded cohort would: an event that kept
  its tumour shift with only sense-strand reads counted would be the host's.

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
- **Unstranded libraries.** An event can count reads of an overlapping gene on the other strand,
  or of a nested gene on the same one (§10). The audit flags these but does not remove them from
  the models.

## References

Literature cited for context and for the concordance checks of §2, §4, §6 and §9 (D-28). None
of it was used to choose a feature, a model or an event.

- Buljan M, et al. Tissue-specific splicing of disordered segments that embed binding motifs
  rewires protein interaction networks. *Mol Cell* 2012;46:871–883.
- Dolnick BJ. Cloning and characterization of a naturally occurring antisense RNA to human
  thymidylate synthase mRNA. *Nucleic Acids Res* 1993;21:1747–1752.
- Dvinge H, Bradley RK. Widespread intron retention diversifies most cancer transcriptomes.
  *Genome Med* 2015;7:45.
- Ellis JD, et al. Tissue-specific alternative splicing remodels protein-protein interaction
  networks. *Mol Cell* 2012;46:884–892.
- Kustatscher G, Hothorn M, Pugieux C, Scheffzek K, Ladurner AG. Splicing regulates NAD
  metabolite binding to histone macroH2A. *Nat Struct Mol Biol* 2005;12:624–625.
- Mayr C, Bartel DP. Widespread shortening of 3′UTRs by alternative cleavage and
  polyadenylation activates oncogenes in cancer cells. *Cell* 2009;138:673–684.
- Novikov L, et al. QKI-mediated alternative splicing of the histone variant macroH2A1 regulates
  cancer cell proliferation. *Mol Cell Biol* 2011;31:4244–4255.
- Sporn JC, et al. Histone macroH2A isoforms predict the risk of lung cancer recurrence.
  *Oncogene* 2009;28:3423–3428.
- Tsai YS, Dominguez D, Gomez SM, Wang Z. Transcriptome-wide identification and study of
  cancer-specific splicing events across multiple tumors. *Oncotarget* 2015;6:6825–6839.
- Xia Z, et al. Dynamic analyses of alternative polyadenylation from RNA-seq reveal a 3′-UTR
  landscape across seven tumour types. *Nat Commun* 2014;5:5274.
- Zong FY, et al. The RNA-binding protein QKI suppresses cancer-associated aberrant splicing.
  *PLoS Genet* 2014;10:e1004289.
