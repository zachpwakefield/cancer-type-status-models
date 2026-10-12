# Cancer type: splicing as a lineage code

One classifier per modality over 34 classes — 33 TCGA projects and a pooled normal class —
9,729 samples, 9 modalities. Every number here is out-of-fold ([RESULTS.md](RESULTS.md)).

This document asks what the splicing models learned about *which tissue a tumour came from*.
The other task, telling tumour from normal tissue within a cancer, is a different question with
a different answer and is written up separately in [TUMOUR_STATUS.md](TUMOUR_STATUS.md).

- **Methods and decisions:** [docs/BIOLOGY_METHODS.md](../docs/BIOLOGY_METHODS.md) (D-13 to D-25).
- **Scripts:** `07`–`19` write the tables, `21` draws these figures.
- **Figures:** T1–T7 in [figures/cancer_type/](figures/cancer_type/), each with a caption file
  and a CSV of the plotted values.
- **Statistics in figures:** raw p. Any q is in the CSV with its FDR family named.

## Summary

1. **No single splicing view matches expression, but their ensemble does**: 0.950 against 0.954
   accuracy, and a better macro recall (0.940 vs 0.936).
2. **Fusing the two beats either alone**: 0.961, correcting 109 of expression's errors and
   introducing 38 (McNemar p = 3.9e-9). A learned stacker does worse than equal weights.
3. **The two modalities fail differently.** Expression confuses anatomical neighbours
   (colon/rectum, oesophagus/stomach, glioma grades); splicing confuses cell types across organs
   (squamous carcinomas) and some tumours with their own normal tissue.
4. **The class-defining events are not expression in disguise.** For the splice-site views they
   correlate with their host gene's expression no more than random events do, and each class's
   top splicing genes and top expression genes overlap by about one gene.
5. **Single exons mark single lineages** — MLANA in melanocytic tumours, RAG2 in thymoma,
   SLC37A2 in squamous carcinomas.
6. **One advantage is a site effect, not biology.** With whole tissue-source sites held out, the
   glioma, kidney and lung advantages survive; the colon/rectum advantage disappears.

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
Lineage-switched exons and tumour-switched exons are drawn from the same kind of exon.

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
