# Cancer type: splicing as a lineage code

One classifier per modality over 34 classes — 33 TCGA projects and a pooled normal class —
9,729 samples, 9 modalities. Every number here is out-of-fold ([RESULTS.md](RESULTS.md)).

This document asks what the splicing models learned about *which tissue a tumour came from*.
The other task, telling tumour from normal tissue within a cancer, is a different question with
a different answer and is written up separately in [TUMOUR_STATUS.md](TUMOUR_STATUS.md).

- **Methods and decisions:** [docs/BIOLOGY_METHODS.md](../docs/BIOLOGY_METHODS.md) (D-13 to D-29).
- **Scripts:** `07`–`19` and `22` write the tables, `21` draws these figures.
- **Figures:** T1–T9 in [figures/cancer_type/](figures/cancer_type/), each with a caption file
  and a CSV of the plotted values.
- **Statistics in figures:** raw p. Any q is in the CSV with its FDR family named.

## The argument

A tumour's tissue of origin is the easiest thing about it to predict from RNA: expression gets
34 classes right 95.4% of the time. So the interesting question is the same as in the other
document — is splicing looking at something else? — and here the answer has a sharper edge,
because splicing is not merely complementary. It is a *lineage code*, and it fails where a
lineage code should fail.

No single splicing view matches expression, but their average does, with a better macro recall
(§1); fused with expression at fixed weights it removes 71 of 450 errors. The two modalities
then fail on different classes, and the pattern is informative rather than random. Expression
confuses anatomical neighbours — colon with rectum, oesophagus with stomach, glioma grades.
Splicing makes almost none of those errors (GBM against LGG: 20 for expression, zero for
splicing) and instead confuses cell types *across* organs: squamous carcinomas of head and neck,
lung, cervix and bladder look alike to it, and several tumours look like the normal tissue they
arose in (§2, §3). That is what a readout of the cell of origin should do.

The evidence that this is lineage rather than expression re-encoded is direct. The class-defining
splice-site events correlate with their host gene's expression no more than random events do,
and each class's top splicing genes and top expression genes overlap by about one gene (§5).
Single exons mark single lineages (§6), and across all eight views a single event is often
enough on its own: ten of them separate a class from the nearest of the 33 others by 0.58 to
1.63 — in PSI, or in the HIT index for terminal exons — with no model involved in that number
(§8).

Three sections keep the claim honest. §4 shows that the pathway enrichments of the splicing
explanations were largely an artefact of deeply covered ribosomal-protein splice sites, which
are measured precisely everywhere; what survives their removal is lineage biology and the
splicing machinery itself. §7 shows that one of the headline gains is not biology at all:
because each TCGA tissue-source site contributes to exactly one project, the colon/rectum
advantage disappears when whole sites are held out, while glioma, kidney and lung survive. And
§9 qualifies §5: the class-defining events are not their *own* gene's expression in disguise,
but in unstranded libraries a terminal exon that overlaps another gene can count that gene's
reads, and 101 of the 695 overlapping top events follow a neighbour's expression — often a
lineage gene such as CP in liver, HNF1A in rectum or PAX6 in glioma. Two of the ten events of
§8 are among them.

## Summary

1. **No single splicing view matches expression, but their ensemble does**: 0.950 against 0.954
   accuracy, and a better macro recall (0.940 vs 0.936).
2. **Fusing the two beats either alone**: 0.961, correcting 109 of expression's errors and
   introducing 38 (McNemar p = 3.9e-9). A learned stacker does worse than equal weights.
3. **The two modalities fail differently.** Expression confuses anatomical neighbours
   (colon/rectum, oesophagus/stomach, glioma grades); splicing confuses cell types across organs
   (squamous carcinomas) and some tumours with their own normal tissue.
4. **The class-defining events are not their own gene's expression in disguise.** For the
   splice-site views they correlate with their host gene's expression no more than random events
   do, and each class's top splicing genes and top expression genes overlap by about one gene.
   Some terminal exons read a *neighbouring* gene instead (point 8).
5. **Single exons mark single lineages** — MLANA in melanocytic tumours, RAG2 in thymoma,
   SLC37A2 in squamous carcinomas. Across all eight views, one event is often enough on its own:
   TOMM40L's terminal exon separates liver from its nearest class by 1.63 of HIT index, and
   MEF2B, EVL and NPM1 fire in DLBC, BRCA and GBM and essentially nowhere else.
6. **One advantage is a site effect, not biology.** With whole tissue-source sites held out, the
   glioma, kidney and lung advantages survive; the colon/rectum advantage disappears.
7. **Lineage and tumour state are the same kind of exon, but not the same events.** The
   class-defining exons have the protein character of the tumour-status ones — in frame,
   disordered, outside folded domains — yet the two tasks' top-100 sets overlap only 28–36% for
   the internal views, and none of the ten events of §8 sits in a splicing-factor regulon the
   tumour-status models' events are enriched for.
8. **Some class-defining events measure a neighbouring gene.** 61% of the top HIT events overlap
   another gene's exon, against 7% of measured HIT events (odds 12.6 with exon length held
   fixed), and 101 of 695 overlapping top events follow the neighbour's expression more than
   their own gene's — PLEKHG3 → SPTB in PCPG, BSCL2 → LRRN4CL in UVM, HPS3 → CP in liver.

---

## 1. Combining the modalities identifies the tissue better than expression alone

![T1](figures/cancer_type/T1_modality_gain.png)

| Model | Accuracy | Macro recall | Errors |
|---|---|---|---|
| best single AS view (SE) | 0.941 | 0.923 | 577 |
| expression (GEX) | 0.954 | 0.936 | 450 |
| AS ensemble (mean of 8 views) | 0.950 | 0.940 | 491 |
| all 9 views, equal | 0.955 | 0.944 | 443 |
| GEX + AS, 1:1 | **0.961** | **0.948** | **379** |
| stacked (cross-fitted) | 0.957 | 0.942 | 414 |

- **The ensemble recovers what single views lose.** Every AS view alone is below expression
  (0.898–0.941), yet averaging their probabilities reaches 0.950 and a *higher* macro recall
  than expression, so the views' errors are largely independent of each other.
- **Fusion is the real gain**: 109 of expression's errors corrected against 38 introduced,
  exact McNemar p = 3.9e-9.
- **Learned weights do not help.** A cross-fitted stacker — which is mildly optimistic by
  construction (D-23) — makes 414 errors against the fixed fusion's 379. The gain comes from the
  modalities being complementary, not from how they are weighted.

**The two tasks read mostly different events.** Taking each model's top 100 features, the events
the 15 tumour-status models use overlap those of the 34 cancer-type classes well above chance
(2–78 fold, `09_ensemble/task_feature_overlap.csv`) but still only partly: 28–36% of the
tumour-status events for the cassette-exon, splice-site and intron views, against 48% (AFE) and
63% (ALE) for the terminal-exon views, which read promoter and poly(A) choice and are the most
lineage-bound of the modalities. That is why the two tasks are written up separately.

## 2. Which classes splicing resolves better, and which worse

![T2](figures/cancer_type/T2_recall_by_class.png)

**Splicing gains** for classes expression confuses with an anatomical neighbour: COAD (+0.114
recall), ESCA (+0.093), READ, GBM/LGG, and for CHOL. §7 shows that only some of these gains
survive held-out sites.

**Expression is better** for the pooled normal class (−0.102), uterine carcinosarcoma, and
squamous and urothelial tumours (BLCA −0.062, HNSC −0.050, CESC).

**Fusion is asymmetric, in a useful way.** Across the 11 classes where splicing alone is worse
than expression, the median shortfall shrinks from −0.034 to −0.002 when the two are fused;
across the 17 where splicing is better, most of the gain survives (median +0.016 to +0.012).
The pooled normal class makes the point: −0.102 for splicing alone, −0.007 fused. The gain in
overall accuracy is not bought at the cost of particular classes.

## 3. The two modalities confuse different classes

![T3](figures/cancer_type/T3_class_confusions.png)

**Expression** confuses anatomical neighbours: colon and rectum (96 errors), oesophagus and
stomach (40), lung adenocarcinoma and squamous carcinoma (37), the kidney subtypes, and GBM with
LGG (20).

**Splicing** makes fewer of those errors — none at all between GBM and LGG — and instead
confuses cell types across organs: CESC with HNSC (28 against 7 for expression), HNSC with
LUSC, BLCA with CESC and HNSC. It also confuses several tumours with the normal tissue they
arose in (PRAD, THCA, LIHC, STAD), which is consistent with splicing tracking the cell of
origin more than the tumour state.

The squamous pattern has been seen without a classifier. In an embedding of exon-skipping PSI
across 32 TCGA cancers, Kahles et al. (2018) found LUSC, CESC and HNSC grouping together,
in their words "even overpowering" the tissue of origin, and found the lineage structure less
pronounced in expression. That is the confusion matrix here, from the other direction.

## 4. Do the views agree on the processes behind a class?

![T4](figures/cancer_type/T4_pathway_coherence.png)

The per-class splicing explanations looked pathway-rich at first: 2,390 significant terms, a
third of them translation-related. That came largely from alternative 5′/3′ splice sites in
cytosolic ribosomal-protein genes (RPS24, RPL13A, RPL3), which are deeply covered and therefore
measured precisely in every tissue — a measurement property, not a lineage programme. With
RPL*/RPS* genes removed from both the ranking and the universe, 950 terms remain (D-16).

What remains is two kinds of biology:
- **lineage programmes** — synapse and cell-projection genes in PCPG (SE), smooth-muscle
  contraction in PRAD, the RHO GTPase cycle in DLBC, extracellular matrix in THCA, adhesion and
  calcium transport in LGG (MXE);
- **the splicing machinery's own transcripts** — RNA-splicing terms from A3SS in ACC, KICH, LUAD
  and UVM.

## 5. The class-defining events are not expression in disguise

![T5](figures/cancer_type/T5_not_expression.png)

- For the splice-site views, top events correlate with their host gene's expression no more than
  random measured events do (median |ρ| 0.13–0.18 against 0.10–0.15).
- The terminal-exon views track expression more, which is expected: an alternative promoter or
  poly(A) site also changes how much of a gene is made.
- Each class's top-100 expression genes share about **one** gene with each splicing view's
  top-100 host genes. The shared ones are lineage markers that are both expressed and spliced
  lineage-specifically: GFAP, MBP, DUOXA2, MYH11, LGALS4, FXYD2.
- The comparison is with the **host** gene. An event can also follow a *different* gene whose
  exon overlaps it, which this section does not test; §9 does, and finds it for 101 of the
  top events, mostly terminal exons.

## 6. Single exons mark single lineages

![T6](figures/cancer_type/T6_lineage_atlas.png)

Across the 34 classes, each class's top cassette exons are switched specifically in that class:
MLANA is skipped in the melanocytic tumours (SKCM, UVM), RAG2 in thymoma; SLC37A2 is included in
the squamous carcinomas (HNSC, LUSC, BLCA, CESC), TAGLN and DENND2B in LAML, PEX10 in PRAD.
Exons shared by related classes form off-diagonal blocks.

These class-defining exons have the same protein character as the tumour-status ones: more often
frame-preserving (54% divisible by 3 against 41% of measured exons), more often in disordered
regions and less often inside folded domains (figure S6 in
[TUMOUR_STATUS.md](TUMOUR_STATUS.md), which tests both selections against the same background).

**The same kind of exon, but not the same exons.** Those two statements sit together and are
worth keeping apart. A tissue switching a lineage exon and a tumour switching a state exon are
both doing the thing cells do with alternative splicing — exchanging a segment in frame, usually
in a disordered, interaction-prone region, rarely through a folded domain — which is why the two
selections look alike against a random background. But the *events* differ: the two tasks' top
100 features overlap only 28–36% for the cassette-exon, splice-site and intron views (§1), none
of the ten lineage events of §8 belongs to a splicing-factor regulon the tumour-status models
are enriched for, and nine of those ten are terminal or first exons, against one of ten on the
tumour-status side. The grammar is shared; the sentences are not.

## 7. How much of this is the tissue-source site?

![T7](figures/cancer_type/T7_site_probe.png)

In TCGA every tissue-source site contributes to exactly one project, so separating neighbouring
projects can exploit site handling rather than biology — and splicing readouts such as intron
retention are plausibly more sensitive to handling than expression is. A fixed, untuned probe was
cross-validated two ways: grouped by patient, and with whole sites held out (D-20).

- **Holds:** glioma grade, kidney subtype, lung histology. Separation survives unseen sites in
  every view (drop ≤ 0.015). The near-perfect splicing separation of GBM from LGG is robust.
- **Fails:** colon versus rectum. No splicing view separates the two at unseen sites (0.48–0.57),
  while expression keeps 0.66. The colon/rectum gain in §2 and §3 should not be read as biology.
- **Partly site-specific:** oesophagus versus stomach, where splicing views lose 0.08–0.14 at
  unseen sites against 0.06 for expression.

This confound is specific to this task. Tumour-versus-normal comparisons are within patient and
within site, so [TUMOUR_STATUS.md](TUMOUR_STATUS.md) is not exposed to it.

## 8. One event per lineage, in detail

![T8](figures/cancer_type/T8_class_events.png)

§6 showed that single cassette exons mark single lineages. This section asks how far one event
alone can go, across all eight splicing views, and names ten of them. The rule has two steps
(D-26): of the five events with the largest mean |SHAP| in a class's own explanation — the
models' own ranking — keep the one whose median in that class is furthest from the nearest other
class, then take the ten clearest classes. The second step decides which examples to show and
nothing else; a feature can rank first for a class and still not separate it by itself, because
the model reads it together with 3,999 others.

The margin below is the gap between the class and its nearest neighbour among the other 33, in
PSI, or in the HIT index for terminal exons. It is a description of the data, not a prediction:
no model is involved in that number.

| Event | Class | Median there | Nearest class | Margin | Over all classes | What the segment is |
|---|---|---|---|---|---|---|
| TOMM40L HIT | LIHC | −0.66 | CHOL 0.97 | −1.63 | 1.00 | last 46 aa of 308, in a domain |
| PLEKHG3 HIT | PCPG | −0.70 | LUAD 0.61 | −1.32 | 1.00 | last 240 aa of 1,219, disordered |
| WDR55 ALE | TGCT | 0.00 | ACC 1.00 | −1.00 | 1.00 | not in a coding transcript |
| MARK3 MXE | LAML | 0.00 | LUAD 0.97 | −0.97 | 1.00 | −48 nt, in frame, aa 615–624 of 753 |
| MEF2B AFE | DLBC | 0.86 | ACC 0.00 | +0.86 | 0.00 | first exon, outside the CDS |
| EVL AFE | BRCA | 0.82 | ACC 0.00 | +0.82 | 0.00 | first exon, aa 1–8 of 390 |
| NPM1 ALE | GBM | 0.77 | ACC 0.00 | +0.77 | 0.00 | not in a coding transcript |
| CTSA A5SS | THYM | 0.81 | LAML 0.11 | +0.69 | 0.01 | 24-nt extension, outside the CDS |
| KRT19 HIT | HNSC | 0.38 | SKCM 0.98 | −0.60 | −1.00 | first 140 aa of 400, in a domain |
| BSCL2 HIT | UVM | 0.14 | SKCM 0.73 | −0.58 | 0.99 | not in a coding transcript |

**Four of them are near-perfect switches.** TOMM40L's terminal exon is used as a last exon in
every class (median HIT index 1.00) except liver, where it reverses to −0.66 — a 1.63 gap to the
nearest class, and the largest of any event here. PLEKHG3 does the same in PCPG, WDR55's last
exon is used everywhere but testis, and MARK3's mutually exclusive pair switches only in LAML.
These are binary, tissue-level decisions: one event, read alone, separates one class from all 33
others.

**Four mark a class against a silent background.** MEF2B, EVL and NPM1 sit at 0.00 in the median
class and at 0.77–0.86 in DLBC, BRCA and GBM; CTSA's 24-nt extension sits at 0.01 everywhere and
0.81 in thymoma. The alternative form exists essentially nowhere else, which for a classifier is
the cheapest possible evidence, and it is part of why §5 found so little overlap with
expression: a promoter or poly(A) choice that fires in one lineage need not change how much RNA
the gene makes.

**Two are graded rather than binary, and both are the interesting cases.** KRT19's terminal exon
runs at −1.00 across most classes, rises to 0.98 in melanoma, and sits between the two at 0.38
in HNSC; BSCL2 separates uveal from cutaneous melanoma by 0.58 with the rest of the classes
beyond SKCM at 0.99. These are the only two rows whose nearest class is a biologically close
one, and they are where a single event stops being enough — which is the whole argument for the
ensembles of §1.

**Two of the ten read a neighbouring gene (§9).** PLEKHG3's terminal exon overlaps an exon of
SPTB, which runs the other way across PLEKHG3's 3′ end, by 3.2 kb, and across all samples the
event follows SPTB's expression (ρ −0.61, against +0.16 with PLEKHG3); SPTB's own expression
peaks in PCPG, the class the event marks. BSCL2's terminal exon overlaps the same-strand
neighbour LRRN4CL and follows it (ρ −0.78, against −0.03 with BSCL2), with LRRN4CL highest in
UVM. Both separate their class, but what they separate it by is how much of a neighbouring gene
the tumour makes. Two more cannot be cleared. TOMM40L's exon overlaps NR1I3, the nuclear
receptor CAR, on the other strand; NR1I3 is a liver gene and the event's liver reversal fits
it, but the correlation (−0.45, against −0.05 with TOMM40L) falls short of the pre-specified
0.5, so by the rule the event stays TOMM40L's. WDR55's last exon overlaps DND1, a germ-cell
gene, on the other strand, and the event switches off in testis — but DND1 is not in the
expression matrix, so the check cannot be run. The other six overlap nothing.

**What the segments do.** Only five of the ten change protein sequence at all: TOMM40L and
KRT19 inside an annotated domain, PLEKHG3 in a disordered region, MARK3 a 16-amino-acid in-frame
swap, EVL the first 8 residues. Two more (MEF2B's first exon, CTSA's 24-nt extension) sit in a
coding transcript but outside its coding sequence, and three belong to no coding transcript of
their gene in GENCODE v45. The lineage signal is carried about as often by untranslated
structure as by protein sequence.

Nine of the ten are terminal exons, first exons or a splice-site shift, and only MARK3 is a
mutually exclusive swap — the opposite balance to the tumour-status events of
[TUMOUR_STATUS.md](TUMOUR_STATUS.md) §9, where six of ten are mutually exclusive pairs. It is
also what §5 predicts: a class is named by where transcription starts and ends, which is the
same decision that makes the terminal-exon views track host-gene expression. It is also what
GTEx shows for normal tissues: alternative start and termination sites, not splicing, account
for most tissue-dependent isoform usage, and most tissue-dependent splicing involves untranslated
exons (Reyes and Huber 2018) — the two features of this list. None of the ten is in a
splicing-factor regulon the tumour-status models' events are enriched for. Lineage decisions and
the tumour programme are different sets of events.

## 9. Which class-defining events measure another gene

![T9](figures/cancer_type/T9_overlap_audit.png)

The same audit as [TUMOUR_STATUS.md](TUMOUR_STATUS.md) §10 (D-27), on the cancer-type models'
top events (top 30 of each class's explanation, per view). Most TCGA mRNA libraries are
unstranded, so where an event's exon overlaps an exon of another gene on the other strand, reads
of either can count towards it; and GENCODE's nested and readthrough models put one gene's exons
inside another's.

**How often.** 6.8% of the 357,832 measured events overlap another gene's exon by at least
50 bp. The cancer-type models' top events overlap more often in every view, and in HIT by far
the most: 61% of the top HIT events against 7% of all measured ones (odds 19.6; 12.6 with the
exon's length held fixed, since long terminal exons overlap something more often by chance). The
other views sit between 1.8 and 2.6 once length is matched.

**How often it matters.** An event *tracks* its neighbour when, over all 9,729 samples, its
Spearman correlation with the neighbour's expression is at least 0.5 in size and larger than
with its own gene's — a rule fixed before any result was seen. 101 of the 695 overlapping top
events do (15%); in HIT, 63 of 224 (28%). Another 22% cannot be checked because the neighbour is
not in the expression matrix.

**What they read.** For this task the neighbours are mostly lineage genes, and the coupling
holds at the level of classes as well as samples. HPS3's terminal exon follows ceruloplasmin
(CP; ρ −0.93), which peaks in liver; a C12orf43 exon follows HNF1A in rectum; an ELP4 exon
follows PAX6 in lower-grade glioma; PKIG's last exons follow ADA in thymoma. Of the classes in
which the tracked neighbours peak, LAML is the most common (20 of 101). These events carry real
lineage information, but it is the neighbour's expression, which is the signal the expression
model already reads. That is a reason to expect less complementarity from them, not a reason to
doubt the classification.

**What does not change.** No model is refitted and no event removed: the models were trained on
these measurements and classify as well as §1 says. §5's comparison with the host gene stands,
and §6's cassette-exon atlas is the view least affected (7 of 62 overlapping top SE events
track). Whether a tracked event is misassigned reads or two neighbouring genes regulated
together cannot be told apart without stranded libraries.

## What would change the reading

The limitations below are the standing caveats. These are the specific observations that would
overturn the argument:

- **More held-out-site failures.** §7 probed six class pairs with a fixed linear probe. If the
  same test were run on every pair and the advantages collapsed broadly rather than only for
  colon against rectum, the lineage-code reading would become a site-handling reading. The probe
  bounds this; it does not settle it.
- **A batch or platform correlate.** The events of §8 are near-binary switches. A switch that
  tracked library preparation or sequencing centre rather than tissue would look exactly the
  same in this data. The site probe is the available check; a cohort processed elsewhere is the
  real one.
- **The pooled normal class.** It mixes tissues, and it is the class splicing does worst on
  (−0.102 recall against expression). Splitting it by tissue would say whether that is a real
  weakness or an artefact of pooling.
- **A stranded cohort.** §9 marks the events that follow an overlapping gene, but in unstranded
  data misassigned reads and co-regulated neighbours look the same. Counting only sense-strand
  reads would settle PLEKHG3/SPTB, BSCL2/LRRN4CL and TOMM40L/NR1I3, and a stranded expression
  matrix that holds DND1 would settle WDR55.
- **The small classes.** CHOL (35 samples), DLBC (47), UCS (55), UVM (78) and ACC (79) are the
  smallest of the 34, and they appear all over §8 — DLBC and UVM as the marked class, CHOL and
  ACC as the nearest one in five of the ten rows. Those medians rest on tens of samples, so the
  margins are the least stable numbers in the document. The ten events themselves were chosen by
  SHAP, which does not see class size; the ranking of the margins does.

## Limitations

- **SHAP explains the models, not causation.** A feature can rank highly because it is precisely
  measured, as the ribosomal-protein splice sites show.
- **Class imbalance.** Classes range from 36 to 1,073 samples; macro recall is the metric to read
  for the small ones.
- **The normal class is pooled** across tissues, so its errors mix tissues.
- **Fusion.** Fused numbers are honest for fixed weights; they are not a nested estimate of a
  method that chooses its weights (D-23).
- **The site probe is a linear probe** with fixed settings. It bounds, but does not measure, how
  site-dependent the tuned MLPs are.
- **Reference data.** Gene sets are MSigDB 2026.1; protein annotation is GENCODE v45.
- **Unstranded libraries.** An event can count reads of an overlapping gene (§9). The audit
  flags these but does not remove them from the models.

## References

Literature cited for context (§3, §8). None of it was used to choose a feature, a model or an
event.

- Kahles A, et al. Comprehensive analysis of alternative splicing across tumors from 8,705
  patients. *Cancer Cell* 2018;34:211–224.
- Reyes A, Huber W. Alternative start and termination sites of transcription drive most
  transcript isoform differences across human tissues. *Nucleic Acids Res* 2018;46:582–592.
