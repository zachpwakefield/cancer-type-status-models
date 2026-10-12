# Subtype within a cancer: what splicing reads

One classifier per cohort and modality: 22 subtype cohorts in 15 TCGA projects, 9 modalities each,
198 models over 7,552 patient–cohort pairs (4,618 patients, one primary tumour each). Every number
here is out-of-fold — each patient was scored by models that never saw it.

The other two tasks ask whether a sample is tumour ([TUMOUR_STATUS.md](TUMOUR_STATUS.md)) and which
tissue it came from ([CANCER_TYPE.md](CANCER_TYPE.md)). This one asks, inside one cancer, which
kind of that cancer a tumour is — and whether splicing knows. How the three tasks' results relate
is in [ACROSS_TASKS.md](ACROSS_TASKS.md).

- **Methods and decisions:** [docs/SUBTYPE_METHODS.md](../docs/SUBTYPE_METHODS.md) (D-30 to D-47,
  amendments A-1 to A-7). Every label map, cohort and threshold was fixed in
  [config/subtype.yaml](../config/subtype.yaml) before the first model was fitted; the regulator
  and site analyses of §9–§10 were fixed in [config/followup.yaml](../config/followup.yaml) before
  they were run.
- **Scripts:** `30`–`35` build the labels, matrices, models and tables; `37` and `38` the
  regulator and site analyses; `36` draws these figures.
- **Figures:** U1–U9 in [figures/subtype/](figures/subtype/). Each has a caption file stating the
  question, the method, n and the test, and a CSV of the plotted values.
- **Tables:** `results/subtype/summary/` (performance, baselines, recall, fusion tests),
  `results/subtype/interpret/` (markers, gene overlap, event audit, cross-family tables),
  `results/subtype/regulators/`, `results/subtype/site_probe/` and one folder per (view, cohort)
  under `results/subtype/<view>/<cohort>/`.
- **Statistics in figures:** raw p. Any q is in the CSV with its FDR family named.
- **Result or context.** Each biological statement is marked as a *result* of this data or as
  *context* from the literature; the literature was never used to choose a feature, a model or an
  event.

## The argument

Inside one cancer, the classes are close relatives, and several were defined from the expression
data a GEX model reads. So the question is not how high the numbers go but where splicing keeps
up with expression and where it does not — and the answer follows what the subtype *is*.

Where a subtype is a cell lineage or a state of the genome, splicing alone matches expression.
Squamous against glandular carcinoma in oesophagus and cervix, seminoma against non-seminoma,
IDH and 1p/19q status in glioma, the MSI and EBV classes of gastric cancer, oestrogen-receptor
status: the AS ensemble is within 0.015 of expression in all of them, and on glioma histology it
is ahead (§2). Where a subtype is one gene switched on — ERG or ETV1 driven by a fusion, PNMT in
adrenal chromaffin cells — splicing falls 0.15–0.17 behind. The AS ensemble also loses much of
several small classes: POLE-ultramutated endometrial cancer (recall 0.08 against 0.54),
myxofibrosarcoma and ETV1-fused prostate cancer (§3). Equal-weight fusion, which helped the pan-cancer task, barely
moves anything here: −0.04 to +0.04 in balanced accuracy, with a CI clear of zero only in glioma
histology and endometrial molecular class.

What the models learned is biology rather than label noise, and their errors show it best. A
histology model that is "wrong" is wrong in the molecular direction: in glioma it calls 41 of 58
non-codeleted oligodendrogliomas astrocytoma, the reclassification WHO made in 2016; its "serous"
endometrioid tumours are copy-number high; every adrenal tumour it calls extra-adrenal is
pseudohypoxic and has lost PNMT (§4). The literature's markers are in the data (39 of 40 in the
expected direction) and in the expression models (27 in a top 100), far less in the splicing models
(8; §5), and a cohort's splicing and expression models share about one gene in fifty (§6). Where
splicing carries biology of its own it is recognisable: leiomyosarcomas switch the mutually
exclusive exons of α-actinin-1 and tropomyosin-1; lobular breast tumours skip two frameshifting
CDH1 exons in their tail; HPV-positive head-and-neck tumours use the p16INK4a first exon of CDKN2A
(§7). One result goes beyond what was asked: a 57-bp CCNDBP1 exon is in the MSI top 5 of three
independently trained models, in three cancers, and nine more MSI events separate MSI from the
rest in all three cancers, including those whose model did not pick them (§8).

Two checks added afterwards, with predictions written down first, bound the claims. The splicing
regulators the literature names behave as predicted in four of six tests (ESRP2 in biphasic
mesothelioma, ESRP1 in mesenchymal head-and-neck cancer, PTBP1 and the ACTN1 smooth-muscle exon
in leiomyosarcoma). The TPM1 exon goes against the published mechanism, and PTBP1 is one of
several factors that move with the leiomyosarcoma pattern, not the only one (§9). And part of
the splicing models' accuracy on HPV status, the HNSC expression classes, cervical histology and
gastric molecular class is specific to the sites that sent the samples, where expression's is
not (§10).

The neighbour-gene audit keeps this honest. Of 2,248 top events, 97 follow an overlapping gene
more than their own — 62 of them in HIT — and four of the nineteen strongest single events are
such readouts; the leiomyosarcoma first exon that looks like an NDE1 promoter is MYH11
transcription read on the other strand (§7).

## Summary

1. **Splicing matches expression on lineage and genome-state subtypes** (U1): ESCA, CESC, TGCT,
   LGG_IDH, STAD_TCGA and BRCA_ER within 0.015 balanced accuracy; LGG histology +0.026 for
   splicing (CI −0.002 to +0.054). In AUROC the AS ensemble equals or beats expression in 7 of 22
   cohorts.
2. **It falls behind where one gene defines the class**: PCPG site −0.174 (CI −0.267 to −0.089;
   McNemar p = 1.1e-5, 21 against 1 discordant patients), PRAD ETS class −0.145 (p = 2.5e-5),
   UCEC molecular class −0.117 (CI −0.205 to −0.029). The losses concentrate in small classes
   (recall: POLE 0.54 → 0.08, myxofibrosarcoma 0.67 → 0.29, ETV1 0.93 → 0.57; U2).
3. **Fixed 1:1 fusion adds little**: between −0.042 and +0.036 in every cohort; CI above zero only
   for LGG (+0.023) and UCEC_TCGA (+0.036), with McNemar p 0.064 and 0.070.
4. **Histology errors are molecularly coherent** (U3): codeletion in glioma, copy-number class in
   endometrial cancer, genomic stability in gastric cancer, pseudohypoxia and PNMT loss in PCPG.
5. **The literature's markers are in the data and in the expression models, rarely in the
   splicing models** (U4), and the two modalities' top-50 genes share a median of one gene (U5).
6. **Named splicing biology, audited** (U6): smooth-muscle ACTN1/TPM1 exons in leiomyosarcoma; a
   CDH1 skipping tail in lobular breast cancer; the p16INK4a first exon in HPV-positive HNSC; the
   NDE1 "promoter" of leiomyosarcoma is MYH11; 97 of 2,248 top events read a neighbour.
7. **MSI leaves the same splicing marks in three cancers** (U7, exploratory): 10 events separate
   MSI from the rest in UCEC, STAD and COAD alike (AUROC ≥ 0.75 or ≤ 0.25 in all three), 7 of
   them chosen by at least two of the three models independently.
8. **Regulators** (U8; R1–R6 fixed in advance): ESRP2 lower in biphasic mesothelioma (AUROC
   0.82), ESRP1 lower in mesenchymal HNSC (0.76), PTBP1 lower in leiomyosarcoma (0.66) and moving
   the ACTN1 smooth-muscle exon against the leiomyosarcoma pattern (partial ρ +0.43). ESRP1 in
   mesothelioma is not resolved, and TPM1 runs against the published mechanism (A-5).
9. **Site** (U9): holding tissue-source sites out costs GEX ≤ 0.035 in 20 of 21 cohorts, but the
   splicing views 0.056–0.098 in HNSC_HPV, CESC, THCA, HNSC_mRNA, TGCT and STAD_TCGA. "Splicing
   loses more than GEX in most cohorts" holds in 16 of 21 (q = 0.080): not resolved.

---

## 1. The labels, and what they already say about each other

No model is involved in this section. Twenty-two cohorts were fixed in `config/subtype.yaml`
before any subtype model was fitted (D-30..D-35, D-43..D-45); counts are patients, one primary
tumour each (`results/subtype/data/`).

| Cohort | Label | Tier | Patients | Classes (patients) | Dropped: class below 20 |
|---|---|---|---|---|---|
| BRCA | histology | primary | 952 of 1,073 | ductal 759, lobular 193 | mucinous 17, metaplastic 9, medullary 5 |
| THCA | histology | primary | 486 of 503 | classical 351, follicular 100, tall cell 35 | – |
| UCEC | histology | primary | 162 of 178 | endometrioid 104, serous 58 | – |
| SARC | histology | primary | 230 of 256 | LMS 102, DDLPS 57, UPS 47, MFS 24 | MPNST 9, synovial 10 |
| TGCT | histology | primary | 114 of 134 | seminoma 64, non-seminoma 50 | – |
| MESO | histology | primary | 76 of 85 | epithelioid 55, biphasic 21 | sarcomatoid 1 |
| STAD | histology | primary | 227 of 352 | intestinal 157, diffuse 70 | – |
| LGG | histology | weak | 378 of 510 | astrocytoma 191, oligodendroglioma 187 | – |
| THYM | histology | weak | 108 of 119 | A/AB 51, B1–B3 57 | thymic carcinoma 11 |
| ESCA | histology | control | 161 of 161 | squamous 80, adenocarcinoma 81 | – |
| CESC | histology | control | 285 of 291 | squamous 240, adenocarcinoma 45 | – |
| PCPG | histology (site) | control | 178 of 179 | pheochromocytoma 148, paraganglioma 30 | – |
| LGG_IDH | genomic | primary | 507 of 510 | IDH-mut codel 171, IDH-mut non-codel 242, IDH-wt 94 | – |
| UCEC_TCGA | genomic | primary | 169 of 178 | POLE 24, MSI 34, CN-low 28, CN-high 83 | – |
| STAD_TCGA | genomic | primary | 320 of 352 | CIN 196, MSI 55, GS 43, EBV 26 | HM-SNV 7 |
| COAD_TCGA | genomic | primary | 235 of 263 | CIN 162, MSI 40, GS 33 | HM-SNV 2 |
| PRAD_ETS | genomic | primary | 303 of 497 | ERG 152, ETV1 28, SPOP 37, other 86 | ETV4 14, FOXA1 9, FLI1 4, IDH1 3 |
| BRCA_ER | IHC | primary | 1,025 of 1,073 | ER-negative 235, ER-positive 790 | – |
| HNSC_HPV | IHC / ISH | primary | 112 of 501 | HPV-negative 79, HPV-positive 33 | – |
| BRCA_PAM50 | expression-defined | – | 1,073 of 1,073 | LumA 558, LumB 212, Her2 77, basal 187, normal-like 39 | – |
| HNSC_mRNA | expression-defined | – | 278 of 501 | atypical 68, basal 86, classical 49, mesenchymal 75 | – |
| PCPG_mRNA | expression-defined | – | 173 of 179 | kinase 68, pseudohypoxia 61, Wnt-altered 22, cortical admixture 22 | – |

Dropped by the size rule (D-34, fewer than 20 patients): seven histology classes and six
molecular classes, listed in the last column; they are dropped, never merged into a neighbour.
Every excluded patient has a recorded reason (`exclusions.csv`): an unmapped "NOS", "other" or
"mixed" value (e.g. 119 STAD adenocarcinomas NOS, 129 LGG oligoastrocytomas), the two
histology sources disagreeing (17 patients over all cohorts; concordance 99.5%), or no value
in the source (e.g. 388 HNSC patients without a p16 or ISH result).

**How each label was made matters for what a model can claim** (D-44):

- *Histology* — a pathologist's reading of the slide, made without RNA.
- *Genomic* — DNA: IDH mutation and 1p/19q copy number (LGG_IDH); mutation burden, MSI, copy
  number and EBV (UCEC_TCGA, STAD_TCGA, COAD_TCGA); fusions and SPOP mutation (PRAD_ETS, whose
  fusion calls also used RNA-seq).
- *IHC* — protein: oestrogen-receptor staining (BRCA_ER); p16 staining or HPV in-situ
  hybridisation (HNSC_HPV).
- *Expression-defined* — clusters of mRNA expression (PAM50, HNSC and PCPG mRNA subtypes). An
  expression model predicting them is reading the data the label was computed from, so its score
  is near-circular and is reported apart from the others. Splicing was not used to make these
  labels, so for the AS views they are a fair test — with the caveat that some splicing
  features partly encode expression (§6, §7).

Published cluster labels were assigned over the full TCGA cohort before this work; that is a
property of the label, not leakage, and it applies equally to every view.

**What the labels already say about each other** (`label_crosstab.csv`, no model):

| Pair | What the counts show | Literature |
|---|---|---|
| LGG histology × IDH/codeletion | 128 of 186 oligodendrogliomas are IDH-mutant codeleted, 38 IDH-mutant non-codeleted, 20 IDH-wildtype; 6 of 189 astrocytomas are codeleted | WHO 2016 redefined both types by IDH and 1p/19q (Louis et al. 2016) — *consistent* |
| UCEC histology × TCGA class | 48 of 53 serous tumours are copy-number high; endometrioid tumours spread over all four classes (28 of 100 CN-high) | serous tumours cluster with CN-high (TCGA 2013) — *consistent* |
| BRCA histology × ER, PAM50 | 180 of 188 lobular tumours are ER-positive and 161 of 193 luminal A | ILC is mostly ER+ luminal A (Ciriello et al. 2015) — *consistent* |
| STAD histology × TCGA class | genomically stable: 17 of 66 diffuse vs 12 of 138 intestinal tumours | GS is enriched for diffuse histology (TCGA 2014) — *consistent* |
| PCPG site × mRNA class | all 27 extra-adrenal tumours with an mRNA class are pseudohypoxia (25) or kinase (2); cortical admixture only in adrenal tumours | SDHx paragangliomas are pseudohypoxic; cortex is adrenal (Fishbein et al. 2017) — *consistent* |
| HNSC HPV × mRNA class | all 8 HPV-positive tumours with an mRNA class are atypical | HPV+ tumours concentrate in the atypical subtype (TCGA 2015) — *consistent* |

## 2. How well each modality recovers the subtype

![U1](figures/subtype/U1_performance.png)

Balanced accuracy (the mean of per-class recall) is the headline metric because most cohorts are
imbalanced; chance and the majority-class predictor both score 1/classes on it. "Majority acc."
is what always guessing the largest class would score on plain accuracy. CIs are 1,000 patient
bootstraps; the McNemar test is fused against GEX on paired per-patient correctness.

| Cohort | n | Classes | Chance BA | Majority acc. | GEX BA [95% CI] | Best AS view | AS ensemble | GEX + AS | Δ fused − GEX [95% CI] | McNemar p |
|---|---|---|---|---|---|---|---|---|---|---|
| BRCA | 952 | 2 | 0.500 | 0.797 | 0.860 [0.829–0.887] | 0.842 (RI) | 0.835 | 0.866 | +0.006 [−0.013, +0.025] | 0.73 |
| THCA | 486 | 3 | 0.333 | 0.722 | 0.686 [0.628–0.744] | 0.675 (SE) | 0.660 | 0.705 | +0.019 [−0.007, +0.042] | 1.2e-4 |
| UCEC | 162 | 2 | 0.500 | 0.642 | 0.899 [0.849–0.945] | 0.848 (RI) | 0.855 | 0.913 | +0.013 [−0.010, +0.042] | 0.69 |
| SARC | 230 | 4 | 0.250 | 0.443 | 0.749 [0.679–0.814] | 0.699 (SE) | 0.661 | 0.738 | −0.011 [−0.056, +0.027] | 0.77 |
| TGCT | 114 | 2 | 0.500 | 0.561 | 1.000 [1.000–1.000] | 1.000 (AFE) | 1.000 | 1.000 | +0.000 [+0.000, +0.000] | 1.00 |
| MESO | 76 | 2 | 0.500 | 0.724 | 0.787 [0.681–0.881] | 0.708 (SE) | 0.759 | 0.745 | −0.042 [−0.106, +0.000] | 0.25 |
| STAD | 227 | 2 | 0.500 | 0.692 | 0.647 [0.577–0.711] | 0.661 (SE) | 0.654 | 0.641 | −0.006 [−0.043, +0.027] | 0.61 |
| LGG (weak) | 378 | 2 | 0.500 | 0.505 | 0.812 [0.772–0.848] | 0.843 (A3SS) | 0.838 | 0.835 | +0.023 [+0.002, +0.045] | 0.064 |
| THYM (weak) | 108 | 2 | 0.500 | 0.528 | 0.798 [0.714–0.871] | 0.794 (HIT) | 0.748 | 0.796 | −0.002 [−0.046, +0.043] | 1.00 |
| ESCA (control) | 161 | 2 | 0.500 | 0.503 | 0.988 [0.969–1.000] | 0.994 (SE) | 0.994 | 0.988 | +0.000 [+0.000, +0.000] | 1.00 |
| CESC (control) | 285 | 2 | 0.500 | 0.842 | 0.924 [0.873–0.966] | 0.924 (SE) | 0.935 | 0.922 | −0.002 [−0.006, +0.000] | 1.00 |
| PCPG (control) | 178 | 2 | 0.500 | 0.831 | 0.880 [0.804–0.943] | 0.725 (AFE) | 0.706 | 0.856 | −0.023 [−0.068, +0.003] | 0.38 |
| LGG_IDH | 507 | 3 | 0.333 | 0.477 | 0.985 [0.973–0.995] | 0.976 (SE) | 0.972 | 0.982 | −0.004 [−0.012, +0.000] | 1.00 |
| UCEC_TCGA | 169 | 4 | 0.250 | 0.491 | 0.678 [0.598–0.755] | 0.555 (HIT) | 0.560 | 0.713 | +0.036 [+0.006, +0.074] | 0.070 |
| STAD_TCGA | 320 | 4 | 0.250 | 0.613 | 0.879 [0.837–0.919] | 0.905 (SE) | 0.877 | 0.878 | −0.001 [−0.019, +0.016] | 1.00 |
| COAD_TCGA | 235 | 3 | 0.333 | 0.689 | 0.835 [0.771–0.894] | 0.749 (A3SS) | 0.772 | 0.813 | −0.022 [−0.061, +0.010] | 1.00 |
| PRAD_ETS | 303 | 4 | 0.250 | 0.502 | 0.920 [0.878–0.957] | 0.807 (SE) | 0.775 | 0.928 | +0.008 [−0.004, +0.027] | 0.62 |
| BRCA_ER | 1,025 | 2 | 0.500 | 0.771 | 0.936 [0.918–0.953] | 0.920 (HIT) | 0.928 | 0.933 | −0.003 [−0.011, +0.004] | 0.55 |
| HNSC_HPV | 112 | 2 | 0.500 | 0.705 | 0.942 [0.887–0.988] | 0.893 (RI) | 0.871 | 0.936 | −0.006 [−0.019, +0.000] | 1.00 |
| BRCA_PAM50 | 1,073 | 5 | 0.200 | 0.520 | 0.862 [0.827–0.896] | 0.818 (HIT) | 0.818 | 0.882 | +0.020 [−0.001, +0.042] | 0.36 |
| HNSC_mRNA | 278 | 4 | 0.250 | 0.309 | 0.894 [0.854–0.930] | 0.807 (HIT) | 0.798 | 0.903 | +0.009 [−0.010, +0.031] | 0.51 |
| PCPG_mRNA | 173 | 4 | 0.250 | 0.393 | 0.935 [0.892–0.976] | 0.816 (SE) | 0.820 | 0.943 | +0.008 [+0.000, +0.019] | 0.50 |

Macro F1, macro AUROC, accuracy and MCC for every model, with CIs, are in
`summary/performance.csv`; every single view's per-class recall is in `summary/recall_by_class.csv`.

**Mean balanced accuracy by label family** (unweighted over cohorts):

| Model | Histology (12) | Genomic / IHC (7) | Expression-defined (3) |
|---|---|---|---|
| best single view (SE) | 0.802 | 0.821 | 0.794 |
| AS ensemble | 0.804 | 0.822 | 0.812 |
| GEX | 0.836 | 0.882 | 0.897 |
| GEX + AS (1:1) | 0.834 | 0.883 | 0.909 |

- **Where splicing matches expression** (result). The controls do what controls should: oesophageal
  squamous against adenocarcinoma is 0.96–0.99 in every view, seminoma against non-seminoma is
  perfect in eight of nine views (ALE 0.992), cervical squamous against adenocarcinoma 0.92–0.94
  for GEX and both ensembles. Beyond the controls, the AS ensemble is within 0.015 of GEX for IDH/codeletion in
  glioma (0.972 vs 0.985; SE alone 0.976), the four TCGA gastric classes (0.877 vs 0.879; SE alone
  0.905, above GEX) and ER status (0.928 vs 0.936). On glioma histology splicing is ahead: 0.838
  against 0.812, and five of the eight single views beat GEX (A3SS 0.843, SE 0.835).
- **Where it falls behind** (result). PCPG site (−0.174), PRAD ETS class (−0.145) and UCEC
  molecular class (−0.117) have CIs well below zero; SARC (−0.089), HNSC HPV (−0.071) and COAD
  molecular class (−0.063) too. §3 shows where the losses sit.
- **Ranking against thresholding.** In macro AUROC the AS ensemble equals or beats GEX in 7 of 22
  cohorts (ESCA, TGCT, LGG, THCA, STAD_TCGA, BRCA_ER, HNSC_HPV): part of its balanced-accuracy
  deficit is where the argmax falls for small classes, not whether the classes are ordered
  correctly.
- **Fusion barely moves.** Fixed 1:1 fusion changes balanced accuracy by −0.042 (MESO) to +0.036
  (UCEC_TCGA). Its paired CI clears zero only for LGG (+0.023, 0.002 to 0.045) and UCEC_TCGA
  (+0.036, 0.006 to 0.074), and neither McNemar test is below 0.05 (0.064, 0.070; q 0.52). THCA's
  McNemar p of 1.2e-4 comes from 31 against 7 discordant patients, mostly classical tumours —
  plain accuracy rises, balanced accuracy only +0.019 (CI −0.007 to +0.042). The pan-cancer task's
  fusion gain (CANCER_TYPE.md §1) does not carry over: within a cancer, the two modalities' errors
  overlap too much.
- **Expression-defined labels** (BRCA PAM50, HNSC and PCPG mRNA classes). GEX reads the data
  these labels were computed from (0.862–0.935). The AS ensemble, for which they are an
  independent test, reaches 0.798–0.820 against chance 0.20–0.25, and on PAM50 basal-like tumours
  (recall 0.979) and normal-like (0.744) it matches or beats GEX.

## 3. Which subtypes splicing loses

![U2](figures/subtype/U2_recall_by_class.png)

Recall of every class, GEX against the AS ensemble (result):

| Class (patients) | GEX | AS ensemble | GEX + AS | What defines the class |
|---|---|---|---|---|
| UCEC_TCGA POLE (24) | 0.54 | 0.08 | 0.58 | ultramutation from POLE exonuclease mutations |
| SARC myxofibrosarcoma (24) | 0.67 | 0.29 | 0.54 | myxoid stroma; genomically close to UPS |
| PRAD_ETS ETV1 (28) | 0.93 | 0.57 | 0.93 | ETV1 switched on by a fusion |
| PCPG extra-adrenal (30) | 0.83 | 0.57 | 0.80 | site; mostly SDHx-pseudohypoxic |
| THCA tall cell (35) | 0.71 | 0.51 | 0.69 | BRAF-like, least differentiated |
| COAD_TCGA GS (33) | 0.73 | 0.55 | 0.64 | no hypermutation, no chromosomal instability |
| LGG astrocytoma (191) | 0.86 | 0.92 | 0.93 | histology (mostly IDH-mutant non-codeleted or IDH-wildtype) |
| THCA classical (351) | 0.56 | 0.66 | 0.63 | BRAF-like |
| UCEC_TCGA CN-high (83) | 0.84 | 0.94 | 0.88 | copy-number high, TP53-mutant, serous-like |

- **The losses are in small classes.** Splicing's five largest recall losses are all in classes
  of fewer than 40 patients, and over the 12 classes that small the AS ensemble is below GEX in 9.
  The gains are in large classes. Single views lose the small classes too, but the ensemble loses
  more: in five of the six small classes above its recall is below that of its best single view
  (POLE 0.08 against MXE 0.38; tall cell 0.51 against SE 0.77, which is above GEX). Averaging
  eight views whose evidence for a small class is weak and disagrees pulls that class's
  probability below the others' (an interpretation; it fits the better AUROC of §2). Picking the
  best view after the fact is optimistic, so this is a reason to look at the views, not a
  result.
- **POLE tumours are called MSI.** Of the 24 POLE-ultramutated endometrial tumours the AS ensemble
  calls 12 MSI and GEX calls 8 (result; `summary/confusion_long.csv`). Both classes are
  hypermutated (context). The MSI events of §8 do not explain the confusion: on them the POLE
  tumours sit with the non-MSI classes (CCNDBP1 exon median 0.08 in POLE, 0.20 in MSI, 0.06 in
  CN-low and CN-high; `interpret/event_audit.csv`). What the splicing models use to call a POLE
  tumour MSI is not identified here.
- **ETV1 against ERG.** GEX recognises ETV1-fused prostate tumours by ETV1 itself (rank 1 in the
  ETV1 class; AUROC 1.00, §5). The splicing models have ERG's own exon (SE rank 2, §7) but no
  ETV1 event in the top 100, and call 6 of 28 ETV1 tumours ERG and 5 "other".

## 4. When a histology model is wrong, its error is molecularly coherent

![U3](figures/subtype/U3_histology_errors.png)

Four histology cohorts have a published molecular classification of the same patients.
Cross-tabulating each histology model's out-of-fold calls against it (counts only, no test; D-41)
shows that the "errors" are not noise (result; the AS ensemble's tables, in the CSV, show the same
pattern):

- **LGG.** The expression model calls 126 of the 128 codeleted oligodendrogliomas correctly, but
  41 of the 58 oligodendrogliomas *without* a codeletion astrocytoma, and all 6 codeleted
  astrocytomas oligodendroglioma. Where pathology and 1p/19q status agree (311 patients) it is right
  for 93%; where they disagree (64) it sides with the molecular class for 47. Since WHO 2016
  an oligodendroglioma is defined by IDH mutation with 1p/19q codeletion (Louis et al. 2016;
  context), so most of the model's errors fall where the morphological label is the superseded
  one. The AS ensemble does the same: 127 of 128 codeleted and 14 of 58 non-codeleted
  oligodendrogliomas called oligodendroglioma.
- **UCEC.** 9 of the 12 endometrioid tumours called serous are copy-number high, against 28 of all
  100 endometrioid tumours with a molecular class. TCGA found serous carcinomas and a minority of
  high-grade endometrioid tumours in one serous-like, copy-number-high cluster (TCGA 2013;
  context): the model's "serous" is that cluster.
- **STAD.** All 17 genomically stable diffuse tumours are called diffuse, and so are 7 of the 12
  genomically stable intestinal ones; the GS class is enriched for diffuse histology (TCGA 2014;
  context).
- **PCPG.** All 11 adrenal tumours the expression model calls extra-adrenal are in the
  pseudohypoxia cluster, and all 22 with cortical admixture are called adrenal.

**What the misclassified tumours look like on the literature's genes** (median log2 CPM by true
and called class, GEX model; `interpret/error_markers.csv`, cells under 5 patients blank; result):

| Cohort | True → called (patients) | Marker | Correct cell | Error cell |
|---|---|---|---|---|
| PCPG | adrenal → extra-adrenal (11) | PNMT | 7.91 (adrenal → adrenal, 137) | 1.62 |
| SARC | LMS → UPS (11) | MYOCD / CNN1 | 7.62 / 11.00 (LMS → LMS, 85) | 0.24 / 4.71 |
| SARC | DDLPS → MFS (6) | MDM2 / CDK4 | 9.83 / 10.04 (DDLPS → DDLPS, 49) | 6.64 / 5.97 |
| BRCA | lobular → ductal (38) | CDH1 | 5.99 (lobular → lobular, 155) | 7.76 |
| UCEC | endometrioid → serous (12) | PGR / L1CAM | 7.04 / 1.17 (correct, 92) | 2.55 / 5.04 |
| THCA | classical → follicular (71) | TPO | 6.57 (classical → classical, 198) | 10.67 |

Each error cell carries the molecular profile of the class the model chose: adrenal tumours
without PNMT (the noradrenergic pseudohypoxic phenotype; Eisenhofer et al. 2004, context),
leiomyosarcomas without the smooth-muscle programme, dedifferentiated liposarcomas without the
12q13–15 amplicon's expression, lobular tumours that kept E-cadherin, endometrioid tumours with a
serous hormone-receptor and L1CAM profile, classical papillary tumours with follicular-variant
iodide-organification levels. These are the cases where the morphological label and the tumour's
molecular state disagree, and the models follow the molecular state.

The labels in these comparisons were assigned by other groups over the full TCGA cohorts; that is
a property of the label, not a leak, and no molecular label entered any histology model.

## 5. The literature's biology, in the data and in the models

![U4](figures/subtype/U4_literature_markers.png)

The 42 markers of U4 were written into `config/subtype.yaml`, with their sources, before any model
was fitted (D-41, D-45). For each the question is twofold: is the expected difference in the data
(AUROC of the gene's expression between the two named classes; no model), and does a model use it?

| Check (cohort) | In the data | GEX model | Splicing models | Reading |
|---|---|---|---|---|
| CDH1 lost in lobular (BRCA) | AUROC 0.90 | rank 1 | AFE rank 2 | *result consistent with* E-cadherin loss defining ILC (Ciriello et al. 2015) |
| RAS-like follicular variant more differentiated: TPO, TG (THCA) | 0.80, 0.77 | TPO rank 63 | TG SE rank 63 | *consistent* with TCGA 2014; SLC5A5 is not (0.55, p = 0.14) |
| Tall cell less differentiated than classical: TPO (THCA) | 0.68 | rank 46 | none in top 100 | *consistent*, weakly |
| Serous: low ESR1/PGR, high CDKN2A and L1CAM (UCEC) | 0.79–0.89 | CDKN2A rank 4, ESR1 18 | none in top 100 | *consistent* with serous being p16-high, hormone-receptor-low |
| Serous falls in CN-high (UCEC × UCEC_TCGA) | 48 of 53 | – | – | *consistent* with TCGA 2013 (§1, §4) |
| Seminoma KIT, SOX17; non-seminoma SOX2 (TGCT) | 0.96, 0.93, 0.99 | SOX2 19, SOX17 31 | none in top 100 | *consistent*; perfect separation by every view |
| LMS smooth muscle MYOCD, CNN1; DDLPS MDM2, CDK4 (SARC) | 0.91–0.94 | MDM2 rank 1, CDK4 5, CNN1 15 | none in top 100; the smooth-muscle splicing of §7 instead | *consistent* with TCGA 2017 |
| Epithelioid MSLN (MESO) | 0.81 | rank 42 | none in top 100 | *consistent* |
| Intestinal CDX2, CDH1 (STAD) | 0.65, 0.69 | ranks > 300 | none in top 100 | weak; Lauren typing is the hardest histology here (0.65) |
| 1p/19q dosage lower in oligodendroglioma (LGG) | 1p 0.85, 19q 0.81 (1q 0.45, 19p 0.50) | 7 of top 100 genes on 1p/19q (9.4 expected) | – | the *data* carry the codeletion; the histology model does not read it gene by gene |
| 1p/19q dosage lower in codeleted (LGG_IDH) | 0.96 | 33 of top 100 on 1p/19q (9.4 expected, p = 5.3e-11) | – | the molecular model does |
| Thymocyte DNTT in B thymoma (THYM) | 0.67 | rank 616 | none in top 100 | weak, as expected for a grouped ordinal scheme |
| Squamous TP63, KRT5 (ESCA, CESC) | 0.95–1.00 | TP63 rank 5 (CESC) | TP63 SE rank 26 (ESCA) | *consistent* |
| Adrenal PNMT, cortex CYP11B1 (PCPG) | 0.84, 0.72 | ~110 | PNMT AFE rank 48 | *consistent*; §4 |
| MLH1 silenced in MSI (UCEC, STAD, COAD) | 0.88–0.96 | rank 1, 1, 1 | A3SS/A5SS ranks 144–299 | *consistent*: sporadic MSI is MLH1 promoter methylation |
| PD-L1 (CD274) high in EBV gastric (STAD_TCGA) | 0.86 | rank 462 | no event selected | *consistent* with TCGA 2014; the model does not need it |
| ERG, ETV1 over-expressed by fusion (PRAD_ETS) | 1.00, 1.00 | rank 1, 1 | ERG SE rank 2; ETV1 none | *consistent*; §3 |
| ESR1 (BRCA_ER) | 0.96 | rank 1 | AFE rank 136 | *consistent* |
| CDKN2A (p16) high in HPV+ (HNSC_HPV) | 0.93 | rank 19 | AFE rank 9 | *consistent*: E7 inactivates RB, p16 rises; §7 |
| Expression-defined: PAM50 ESR1, ERBB2, KRT5; HNSC atypical CDKN2A; PCPG cortical CYP11B1, kinase PNMT | 0.69–0.99 | ERBB2 2, KRT5 4, CYP11B1 10 | KRT5 AFE 95, PNMT AFE 30 | *consistent*; near-circular for GEX |

(Ranks are a gene's best rank in the model, over all patients or within either named class;
splicing ranks are the best of the gene's events across eight views.)

- **In the data, nearly always** (result): 39 of 40 gene markers differ in the expected direction
  at raw p < 0.05; SLC5A5 in THCA is the exception.
- **In the expression models, mostly** (result): 27 of the 40 genes rank in their GEX model's top
  100, several at rank 1 (CDH1, MDM2, MLH1 three times, ERG, ETV1, ESR1).
- **In the splicing models, rarely** (result): 8 of 40 have an event in a splicing top 100 — CDH1,
  TG, TP63, PNMT twice, ERG, CDKN2A and KRT5. The splicing models reach the same subtypes through
  other genes (§6).
- **LGG is the instructive case.** The arm-level dosage signal is plainly in the expression data
  (AUROC 0.85 for 1p and 0.81 for 19q, against 0.45 and 0.50 for the control arms), and the
  molecular model concentrates on it (33 of its top 100 genes on 1p/19q). The histology model does
  not (7), yet its errors follow codeletion (§4): it reads the oligodendroglial state that
  codeletion produces, not the dosage itself.

## 6. Splicing and expression read different genes

![U5](figures/subtype/U5_as_vs_gex_genes.png)

For each cohort and splicing view, the top 50 genes of the splicing model (a gene ranked by its
best event) against the top 50 of the cohort's expression model, within the genes both could
select (D-41). Result:

- **About one gene in fifty.** The median overlap is 1, 55 of 176 cohort–view pairs share none,
  and the total, 236 shared genes, is about twice the 109 expected by chance. Only 7 pairs exceed
  chance at q < 0.05 (BH over the 176 pairs).
- **Retained introns overlap most** (median 2, 53 genes in all), and that is expected: an
  intron's PSI can be high where its gene is barely transcribed, so RI partly measures expression.
  ESCA's RI model shares 8 genes with its GEX model, most of them intestinal or glandular (HNF1A,
  LGALS4, MYO7B, KIF12, ERN2); plausibly their silence in squamous carcinoma is what makes their
  introns look "retained" (interpretation, not tested).
- **MXE overlaps next** (median 2). In SARC it shares six genes with expression, five of them the
  contractile machinery (ACTN1, TPM1, TPM2, PDLIM7, LDB3) plus SNAP25; ACTN1 and TPM1 are the genes
  whose mutually exclusive exons switch in leiomyosarcoma (§7).
- **HIT overlaps least** (median 0): terminal-exon choice is the view most independent of the
  expression model's genes, and also the view most exposed to neighbouring genes (§7).

So the near-equal performance of §2 is not the same model twice: where the two modalities reach
the same answer they mostly get there through different genes. Different genes do not mean
different mistakes, though. Over the 22 cohorts the AS ensemble repeats 657 of expression's 963
out-of-fold errors (68%; from `performance.csv` and the discordant counts in `fusion_tests.csv`),
and §4 shows why: many of those patients are ones whose label and molecular state disagree. That
is why fixed fusion gains little here.

## 7. The strongest single events, and which of them measure a neighbouring gene

![U6](figures/subtype/U6_event_audit.png)

Before an event is named here it is audited as in D-27 (D-41). Its segment is intersected with
every other GENCODE v45 gene's exons (≥ 50 bp); the segment is the exon, both exons of an MXE, the
intron of an RI, or the extension of an A3SS/A5SS. An overlapping event *tracks* the other gene when,
across the cohort's patients, its Spearman correlation with that gene's expression is at least 0.5
in size and larger than its correlation with its own gene. Most TCGA libraries are unstranded, so
an exon that overlaps another gene on the opposite strand can count that gene's reads.

**How often** (result; `interpret/event_audit.csv`). The top 5 events of a class, across every view
and cohort, give 2,248 distinct (cohort, view, event) entries:
- 371 overlap an exon of another gene (323 antisense, 48 same-strand) and 133 sit inside a
  readthrough model;
- 97 track the other gene;
- 65 overlap a gene that the expression matrix does not hold, so they cannot be cleared.

HIT is the outlier, as it is in the cancer-type task. 64% of its top events overlap another gene
antisense or same-strand, because its long terminal metaexons meet the 3′ ends of convergent
genes, and 62 of 272 track the neighbour. The internal views sit at 4–10% (A3SS, A5SS, MXE, RI, SE),
with 1–5 tracking events each; AFE and ALE are at 14% and 21%.

**U6** applies the pre-specified rule: for each histology and genomic/IHC cohort, take the one event
with the largest class-specific importance across the eight views, and print every audit result.
- **Readouts of a neighbour (4 of 19):**
  - the BRCA HIT event of ZNF276 follows FANCA (ρ −0.66, host +0.16);
  - the ERCC2 HIT event chosen in both UCEC cohorts follows KLC3 (ρ −0.88 and −0.89);
  - the HNSC_HPV HIT event of VRK2 follows FANCL (ρ −0.60).
- **Not cleared (2):** THCA YY1 and PCPG RINL overlap an antisense gene that is not in the
  expression matrix.
- **Follow their own gene (7, with |ρ| ≥ 0.5):** SARC SLMAP −0.76, PRAD_ETS HLA-DMB −0.81, ESCA
  NHERF4 intron −0.72, LGG_IDH KHNYN +0.80, CESC MNX1 −0.60, BRCA_ER SEC16A +0.52, and THCA YY1.
  Such an event re-encodes its gene's expression at least in part, so its rank alone does not make
  it splicing biology.

**Named events, with their audit.** The table lists every event named in this document; each row
is in the audit table (D-41, A-2). "Result" is what this data shows; "context" is the literature.

| Cohort · class | Event | Result | Audit | Context |
|---|---|---|---|---|
| SARC · LMS | ACTN1 MXE, 81-bp exon (chr14:68.88 Mb); LMS MXE rank 1 | median PSI 0.28 in LMS vs 0.87–0.95 in DDLPS, UPS and MFS | no overlap; host ρ −0.44 | ACTN1's mutually exclusive exon pair is a classic smooth-muscle-regulated switch (Southby et al. 1999) |
| SARC · LMS | TPM1 MXE, 126-bp exon (chr15:63.04 Mb); rank 2 | 0.63 in LMS vs 0.06–0.14 | no overlap; host ρ +0.53 | tropomyosin-1 has smooth-muscle-specific exon choices (Gooding et al. 1994) |
| SARC · LMS | two TPM4 first exons (chr19:16.08 and 16.07 Mb); LMS AFE ranks 3 and 5 | one used at 0.35 in LMS vs 0.86–1.00 elsewhere; the other at 0.65 vs 0–0.15 | first: antisense to a gene not in the matrix (not cleared); second: no overlap; host ρ −0.52 / +0.52 | TPM4 has alternative promoters; a promoter switch in LMS is *consistent* with a smooth-muscle isoform but is not shown here |
| SARC · LMS | NDE1 first exon (chr16:15.72 Mb); LMS AFE rank 2 | 0.90 in LMS, 0 elsewhere | **antisense MYH11, ρ +0.75 vs host +0.30: a readout of MYH11**, the smooth-muscle myosin heavy chain | not an NDE1 promoter |
| BRCA · lobular | CDH1 exons 6 (145 bp) and 11 (146 bp), SE; lobular SE ranks 1 and 3 | median 1.00 in both classes, but the lowest decile of lobular tumours reaches 0.95 and 0.96, vs 0.996 and 0.992 in ductal | no overlap; host ρ −0.14 and −0.00 | neither exon is a multiple of 3, so skipping shifts the frame; CDH1 alterations define ILC (Ciriello et al. 2015). Mutation calls are not local, so whether the skipping tumours carry splice-site mutations is *not tested* |
| BRCA · ductal | CDH1 first exon (AFE); ductal AFE rank 2 | observed in 76% of ductal and 33% of lobular tumours | no overlap; host ρ +0.09 | the AFE value mostly records whether CDH1 is transcribed at all, which is E-cadherin loss again |
| THCA · follicular | DUOXA1 HIT; follicular HIT rank 2 | median −0.02 in follicular vs 0.44 in classical and 0.65 in tall cell | **antisense DUOXA2, ρ −0.82 vs host −0.36: a readout of DUOXA2**; the same event is the HNSC_mRNA atypical class's HIT rank 2 (ρ −0.75) | DUOXA2 belongs to the thyroid-hormone synthesis machinery; higher expression in the more differentiated, RAS-like follicular variant is *consistent* with TCGA 2014 |
| HNSC_HPV · HPV-positive | CDKN2A first exon 1α, the p16INK4a promoter (AFE event rank 11; gene rank 9) | usage 0.36 in HPV-positive vs 0 in HPV-negative; observed in 94% vs 47% | no overlap; host ρ +0.67 | E7 inactivates RB and p16INK4a rises, which is why p16 is the clinical HPV surrogate. The model reads p16's exon, not ARF's exon 1β. *Consistent* |
| PRAD_ETS · ERG | ERG 95-bp cassette exon (chr21:38.47 Mb); ERG SE rank 2 | 0.38 in ERG-fused tumours, 0 elsewhere | no overlap; host ρ +0.64 | the exon is seen only where ERG is transcribed, so this is the fusion read through splicing |
| UCEC, STAD, COAD · MSI | CCNDBP1 57-bp exon and nine more (§8) | see U7 | no overlap; host ρ −0.72 to +0.25 | see §8 |

**Leiomyosarcoma** is the clearest case here of splicing carrying biology of its own.
- **Result.** The class's top two MXE events are exon switches in the contractile apparatus, which
  the literature describes as smooth-muscle-regulated (context).
- **Same genes, different readout.** The expression model ranks the same genes in its top list
  (§6), but as expression, not exon choice.
- **Breadth (result).** The switches are not specific to the tumour. The same ACTN1 exon is the top
  MXE event of the PAM50 Normal-like class, at a smaller shift (median 0.85 vs 0.91–0.95), and the
  NDE1/MYH11 readout is that class's top AFE event (0.48 vs 0). So the switches record smooth-muscle
  content wherever it is, and in leiomyosarcoma that content *is* the lineage.
- **Caveat.** One of the five top LMS first-exon events is not splicing at all.

## 8. MSI leaves the same splicing marks in three cancers (exploratory)

![U7](figures/subtype/U7_msi_events.png)

**Status.** This section is exploratory (A-4). The same 57-bp CCNDBP1 exon turned up in the MSI
top 5 of three models trained separately on three cancers. The analysis was added after that was
seen, and it was then widened to every audited top event of every class that more than one cohort
defines the same way. No model, label or threshold changed.

**What was tested** (result; `interpret/shared_class_events.csv`). Every top-5 MSI event of the
UCEC, STAD or COAD models (102 distinct events) was measured as MSI against the rest of the cohort
in the other two cancers. No model is involved, only the event's values and the labels.
- Of the 204 cross-cancer tests, 109 have BH q < 0.05, and 37 also lie ≥ 0.25 from AUROC 0.5.
- Other shared classes replicate less: CIN 30 of 77 (3 strong), genomically stable 13 of 80 (0).
- Squamous and adenocarcinoma, between ESCA and CESC, replicate as much as MSI: 55 of 74 (31 strong)
  and 61 of 74 (39 strong).

**U7** keeps the events that are ≥ 0.25 from 0.5 in the same direction in all three cancers. There
are 10 events at 9 loci; CCNDBP1 counts twice, because its exon is read both as a cassette (SE)
and as an alternative 3′ splice site (A3SS).
- **Chosen independently:** 7 of the 10 were chosen by at least two of the three models, and three
  by all three: the CCNDBP1 exon, a RECQL4 alternative 5′ site and a PTP4A2 cassette exon.
- **CCNDBP1:** its exon is included at a median of 0.20 against 0.06 in UCEC, 0.28 against 0.04 in
  STAD and 0.28 against 0.04 in COAD (MSI against the rest). The AUROCs are 0.83, 0.93 and 0.97.
- **Direction:** two loci are more included in MSI (CCNDBP1, ANXA7). Seven are less included or
  shifted away from the usual site: RPL22L1, KLC1, EHMT1, RECQL4, CDC16, TMX3 and PTP4A2.
- **Audit:** none of the ten overlaps another gene. All but one barely follow their own gene's
  expression (|ρ| ≤ 0.30). The exception is RPL22L1, whose A3SS value falls as RPL22L1 expression
  rises (ρ −0.72 in UCEC, −0.70 in STAD).
- **Splicing as a detector:** the AS ensemble recalls MSI about as well as expression does. UCEC
  0.65 against 0.65, STAD 0.95 against 0.91, COAD 0.93 against 0.95 (U2).
- **MSI, not mutation load:** in UCEC_TCGA the POLE-ultramutated tumours are hypermutated but
  mismatch-repair proficient, and on these events they sit with the non-MSI classes (CCNDBP1 0.08
  against 0.20 in MSI and 0.06 in CN-low and CN-high; TMX3 1.00 against 0.89; RPL22L1 0.16 against
  0.02; `event_audit.csv`, 24 POLE patients). That fits a mismatch-repair route to these events
  rather than mutation load alone (interpretation).

**Is it MSI, or the immune infiltrate that comes with MSI?** (Result; `interpret/msi_event_confounds.csv`,
added under A-4 after U7 was drawn.) MSI tumours are immune-rich, and in COAD the immune genes
themselves separate MSI from the rest (AUROC 0.66 for PTPRC, 0.74 for CD8A and 0.83 for CIITA).
The ten U7 events do not follow the infiltrate:
- **Correlation with immune genes:** in all three cancers, |ρ| ≤ 0.28 with PTPRC, ≤ 0.32 with CD8A
  and ≤ 0.37 with CIITA.
- **Low-infiltrate half:** within the half of each cohort with the lower PTPRC (11 to 24 MSI
  patients), every one of the 30 event × cancer cells still lies ≥ 0.15 from AUROC 0.5. Sixteen
  move further from 0.5.
- **MLH1:** the events correlate with MLH1 expression at |ρ| 0.14–0.40, in the direction of MLH1
  loss. That is expected, because silencing of MLH1 underlies most sporadic MSI. It does not
  separate cause from co-occurrence.
- **Where the infiltrate does show:** across all 306 MSI event × cancer cells, five correlate with
  PTPRC at |ρ| ≥ 0.5. They include the COAD model's HIT event in CIITA (ρ 0.88 with CIITA
  expression), which is the interferon-γ response read as a terminal exon.

**Context, and what is not shown.**
- **RPL22L1 (consistent, not tested).** RPL22 carries a frameshift hotspot in MSI tumours (Novetsky
  et al. 2013). RPL22 represses its paralogue RPL22L1 (O'Leary et al. 2013) and regulates its
  splicing (Zhang et al. 2017). An RPL22L1 event that is lower in MSI, and that follows RPL22L1
  expression, fits RPL22 loss. RPL22 mutation calls are not local, so this is not tested here.
- **Mechanism for the others (candidate, not tested).** A known route from MSI to splicing is the
  deletion of intronic mononucleotide repeats near splice sites; HSP110 exon 9 skipping in MSI
  colorectal cancer is the textbook case (Dorard et al. 2011). No genome sequence is local, so
  whether the nine loci sit next to such repeats was not checked.
- **CCNDBP1.** I know of no prior report linking CCNDBP1 splicing to MSI.
- **A3SS caveat.** Only the RPL22L1 event has a host-expression correlation large enough that it
  may partly encode expression rather than exon choice.

## 9. Splicing regulators of the subtypes

![U8](figures/subtype/U8_regulators.png)

Six predictions were written into `config/followup.yaml` (D-46) before any splicing factor was
looked at by subtype. Then all 70 factors of D-18 were scanned in every cohort
(`results/subtype/regulators/`).

**The predictions** (result; BH over R1–R6):
- **ESRP2 is lower in biphasic mesothelioma** than in epithelioid (AUROC 0.82, q = 3.5e-5). ESRP1
  is not (0.55, q = 0.52): it is barely expressed in either class (median 1.2 log2 CPM), so R1 is
  not resolved. ESRP1 and ESRP2 drive the epithelial splicing programme and are lost in
  epithelial–mesenchymal transition (Warzecha et al. 2009, 2010), and a biphasic tumour carries a
  sarcomatoid component (context).
- **ESRP1 is lower in the mesenchymal HNSC class** than in the other three (AUROC 0.76,
  q = 1.6e-10).
- **PTBP1 is lower in leiomyosarcoma** than in the other sarcomas (AUROC 0.66, q = 3.5e-5).
- **PTBP1 and the ACTN1 smooth-muscle exon.** Adjusting for class, stromal and immune scores,
  PTBP1 moves the ACTN1 event against the leiomyosarcoma direction (partial ρ +0.43 against a
  leiomyosarcoma shift of −0.62; p = 1.4e-11). The same holds within leiomyosarcoma alone
  (ρ +0.57). That is what PTB's repression of the α-actinin smooth-muscle exon predicts (Southby
  et al. 1999).
- **PTBP1 and TPM1 run against the published mechanism.** The TPM1 exon 2/3 event also moves
  against the leiomyosarcoma direction as PTBP1 rises (ρ −0.22, shift +0.53; p = 7.1e-4). That
  matches the sign written down in advance, but the sign was a misreading of the source (A-5). In
  TPM1, PTB works with MBNL1 to *repress exon 3*, the non-smooth-muscle exon (Gooding et
  al. 1994, 1998, 2013), so the mechanism predicts the opposite. MBNL1, not PTBP1, is the factor
  that follows the TPM1 event most closely, and in the published direction (ρ +0.35).

**How specific** (exploratory, A-6). The ACTN1 and TPM1 correlations survive adjustment for each
sample's mean MXE PSI, the global axis of D-24 (ρ +0.43 and −0.22). But PTBP1 is not unique.
- **ACTN1:** PTBP1 ranks 2nd of 70 factors, behind HNRNPC (+0.45). CSTF2 and HNRNPA2B1 follow,
  and MBNL1 runs the other way.
- **TPM1:** PTBP1 ranks 12th.

Read together, the sarcoma result is a PTB–MBNL axis rather than PTBP1 alone. The ACTN1 half
agrees with the literature and the TPM1 half is only half explained by it.

**The scan** (result). How many factors differ between a cohort's classes (q < 0.05) varies from
4 of 70 (PCPG) to 68 (BRCA_PAM50). The highest counts are in subtypes that differ in
proliferation or lineage: PAM50, the HNSC expression classes, gastric molecular classes, sarcoma
histology and glioma IDH status.
- **RBM47** separates classes more than any other factor. It is higher in adenocarcinoma than in
  squamous carcinoma of the oesophagus (AUROC 0.97) and cervix (0.89), in luminal A breast
  tumours, and in cortical admixture in PCPG. Loss of RBM47 marks breast cancer progression
  (Vanharanta et al. 2014; context).
- **Three glioma factors are a gene-dosage result.** SRSF4, SRSF11 and PTBP2 are lowest in
  codeleted tumours (median 5.9 against 6.7 log2 CPM for SRSF4), and all three sit on 1p. SRRM3
  (chromosome 7) runs the other way, highest in codeleted tumours (5.8 against 3.0 in
  IDH-wildtype).

**Event association** (descriptive). The excess of factors tracking a cohort's top events, over
background events, is largest where a whole lineage differs: TGCT (30 factors with excess ≥ 0.1),
ESCA (13; RBM47 +0.66, QKI +0.59) and THYM (10). In sarcoma, RBM24, PTBP1, MBNL1 and PTBP2
lead (+0.21 to +0.24). Unlike R5/R6, this association runs across classes. So a factor that
differs between classes tracks every event that does: the excess measures how much a factor moves
with the subtype, not that it regulates the subtype's exons.

## 10. Is the subtype signal tied to the site that sent the sample?

![U9](figures/subtype/U9_site_probe.png)

In TCGA each tissue-source site sends one project, and within a project sites differ in which
subtypes they send: Cramér's V between site and class is 0.21 (BRCA_ER) to 0.80 (ESCA). The probe
of D-20 was rerun on 21 cohorts (D-47; SARC excluded because myxofibrosarcoma comes from 4 sites):
one untuned classifier per view, cross-validated once with patients grouped and once with every
test site unseen.

- **Expression barely depends on site.** GEX loses ≤ 0.035 balanced accuracy in 20 of 21 cohorts;
  the exception is THCA (0.116), where tall-cell tumours come from 5 sites.
- **Splicing depends on it more, and in particular cohorts** (result): HNSC_HPV loses 0.098, CESC
  0.077, THCA 0.077, HNSC_mRNA 0.073, TGCT 0.057 and STAD_TCGA 0.056 (median of the 8 views).
  GEX loses −0.002 to 0.014 in the same cohorts, THCA apart.
- **P1** (splicing loses more than GEX in most cohorts): 16 of 21, median difference +0.010, sign
  test p = 0.027. Under the family rule decided after the p-values were seen (A-7), q = 0.080, so
  it is not resolved.
- **P2** (the loss grows with site–class association): ρ +0.12 (GEX) and +0.18 (splicing),
  not supported. The cohorts that lose most are not the most site-confounded. ESCA, with the
  strongest confounding (V = 0.80), loses almost nothing (0.006), because squamous and glandular
  carcinoma differ in every view.

**Reading.** Where the classes differ in lineage, splicing transfers across sites. Where they
differ less sharply (HPV status, expression subtypes, gastric molecular class), part of the
splicing models' accuracy is site-specific, plausibly how each site handles its samples
(interpretation). It is not confined to the readouts most sensitive to handling: in HNSC_HPV, RI
loses most (0.142) but SE and A5SS are close behind (0.125, 0.112). The probe's accuracies are
those of a simple classifier, not of the subtype models (D-20).

## What would change the reading

- **Mutation and copy-number calls joined to the events.** Three named results in §7–§8 rest on a
  mechanism this repository cannot see:
  - CDH1 exon skipping in lobular tumours, if those tumours turn out to carry CDH1 splice-site
    mutations;
  - the RPL22L1 event, if RPL22 frameshift status explains it;
  - the POLE class that splicing misses (§3), if its POLE-mutant tumours share no splicing
    programme at all.

  A finding that the skipping tumours carry no CDH1 lesion, or that the RPL22L1 event is unrelated
  to RPL22, would remove those mechanisms from the reading.
- **The MSI events in another cohort.** U7 was found and tested within TCGA (three cancers, one
  pipeline). The claim would weaken considerably if the events fail in an independent MSI cohort
  (e.g. CPTAC colon or endometrial, or a cell-line panel with known MMR status) processed through
  the same quantification. It would also weaken if they follow tumour purity once purity is
  estimated properly; only PTPRC-based halves were tried here.
- **Feature budget rather than biology.** Every view gets the same pre-filter and K search.
  Suppose GEX's lead in PCPG and PRAD_ETS disappears when a splicing model is allowed the single
  host-gene events of PNMT, ERG or ETV1 by construction. Then §3's "one gene switched on" reading
  is a statement about the pre-filter, not about splicing.
- **Labels re-reviewed.** If centrally re-reviewed histology moved the LGG and THYM results, the two
  weak cohorts (and the "errors are molecularly coherent" reading of §4 for LGG) would need
  revisiting. The cross-family tables suggest such a review would mostly agree with the models.
- **Different folds.** Small classes (20–40 patients) have 95% recall CIs with half-widths up to 0.20, and the
  greedy splitter leaves 2 patients of a class in one test fold in the worst case. A repeat with
  other fold seeds that changed the small-class ranking of GEX against the AS ensemble would change
  §3.

- **A site-balanced cohort.** §10's site-specific share of the splicing accuracy would be
  confirmed if a cohort sequenced at one centre, or with handling recorded, showed the same
  losses; it would be refuted if HPV status or the HNSC classes were recovered as well across
  sites once RNA quality (not recorded in TCGA) was matched.
- **Perturbation of PTBP1 and MBNL1 in leiomyosarcoma cells.** §9's PTB–MBNL reading is
  correlational. Knock-down data showing the ACTN1 and TPM1 events move with the factors, in
  the directions seen here, would make it a mechanism; no movement would make it co-regulation.

## Limitations

- **One resource, no external test.** Every cohort is a TCGA project, scored out-of-fold. There is
  no independent cohort, and the expression-defined labels (PAM50, HNSC mRNA, PCPG mRNA) were made
  from expression data overlapping the GEX input, so GEX on them is near-circular (D-44).
- **Small classes.** Classes of 20–40 patients set the bar for inclusion (D-34), so their recall
  intervals are wide. MESO, UCEC_TCGA and COAD_TCGA each have an outer test fold with only 2
  patients of their smallest class (`split_summary.csv`). Metrics are pooled out-of-fold, so no
  per-fold figure is reported.
- **Two weak histology cohorts.** LGG histology (oligoastrocytoma excluded, but the
  astrocytoma/oligodendroglioma split is pre-2016) and THYM (an ordinal WHO scheme, grouped) are
  labelled weak in advance (D-30) and kept apart in every summary.
- **Label provenance.**
  - The PRAD ETS-fusion calls used RNA-seq as well as DNA, so a model reading ERG or ETV1 is
    partly reading the evidence the label was built from.
  - HNSC HPV status is known for a minority of patients (112 of 501) and is confounded with tumour
    site, since oropharynx tumours are mostly HPV-positive. Site was not modelled.
- **No mutation, methylation or copy-number data locally** for CDH1, POLE, RPL22 or MLH1. Every
  mechanism in §7–§8 is therefore context, not result.
- **Unstranded libraries.** The audit can flag but not resolve antisense overlap. 65 top events
  overlap a gene the expression matrix does not hold and stay uncleared, including two of U6's
  nineteen.
- **Combination rules.** The AS ensemble averages eight views with equal weight, which shrinks
  minority-class probabilities when most views are weak on that class (§3). Fusion is fixed at 1:1
  with no learned weights (D-40). The "best AS view" columns are picked after the fact, are
  descriptive only, and are not used in any test.
- **Regulators are read from expression.** Factor activity is approximated by mRNA level; the
  event association of §9 runs across classes, so it cannot separate regulation from co-variation
  with the subtype. Only R5/R6 adjust for class.
- **Exploratory additions.** The error-marker, shared-class and MSI-confound tables and U7 were
  added after results had been read (A-4), as were the R5/R6 sensitivity checks (A-6). They are
  descriptive, and none feeds a model. A-1 to A-3 changed presentation or widened an audit only.

## References

Literature used as *context* (never to choose a feature, a model or an event):

- Brat DJ et al. (2015) Comprehensive, integrative genomic analysis of diffuse lower-grade gliomas. *N Engl J Med* 372:2481.
- Cancer Genome Atlas Research Network (2013) Integrated genomic characterization of endometrial carcinoma. *Nature* 497:67.
- Cancer Genome Atlas Research Network (2014) Comprehensive molecular characterization of gastric adenocarcinoma. *Nature* 513:202.
- Cancer Genome Atlas Research Network (2014) Integrated genomic characterization of papillary thyroid carcinoma. *Cell* 159:676.
- Cancer Genome Atlas Research Network (2015) The molecular taxonomy of primary prostate cancer. *Cell* 163:1011.
- Cancer Genome Atlas Network (2015) Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature* 517:576.
- Cancer Genome Atlas Research Network (2017) Comprehensive and integrated genomic characterization of adult soft tissue sarcomas. *Cell* 171:950.
- Cancer Genome Atlas Research Network (2017) Integrated genomic characterization of oesophageal carcinoma. *Nature* 541:169.
- Ceccarelli M et al. (2016) Molecular profiling reveals biologically discrete subsets and pathways of progression in diffuse glioma. *Cell* 164:550.
- Ciriello G et al. (2015) Comprehensive molecular portraits of invasive lobular breast cancer. *Cell* 163:506.
- Dorard C et al. (2011) Expression of a mutant HSP110 sensitizes colorectal cancer cells to chemotherapy and improves disease prognosis. *Nat Med* 17:1283.
- Eisenhofer G et al. (2004) Distinct gene expression profiles in norepinephrine- and epinephrine-producing hereditary and sporadic pheochromocytomas: activation of hypoxia-driven angiogenic pathways in von Hippel-Lindau syndrome. *Endocr Relat Cancer* 11:897.
- Fishbein L et al. (2017) Comprehensive molecular characterization of pheochromocytoma and paraganglioma. *Cancer Cell* 31:181.
- Gooding C et al. (1994) Smooth muscle-specific switching of α-tropomyosin mutually exclusive exon selection by specific inhibition of the strong default exon. *EMBO J* 13:3861.
- Gooding C, Roberts GC, Smith CW (1998) Role of an inhibitory pyrimidine element and polypyrimidine tract binding protein in repression of a regulated α-tropomyosin exon. *RNA* 4:85.
- Gooding C et al. (2013) MBNL1 and PTB cooperate to repress splicing of Tpm1 exon 3. *Nucleic Acids Res* 41:4765.
- Liu J et al. (2018) An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell* 173:400.
- Liu Y et al. (2018) Comparative molecular analysis of gastrointestinal adenocarcinomas. *Cancer Cell* 33:721.
- Louis DN et al. (2016) The 2016 World Health Organization classification of tumors of the central nervous system. *Acta Neuropathol* 131:803.
- Novetsky AP et al. (2013) Frequent mutations in the RPL22 gene and its clinical and functional implications. *Gynecol Oncol* 128:470.
- O'Leary MN et al. (2013) The ribosomal protein Rpl22 controls ribosome composition by directly repressing expression of its own paralog, Rpl22l1. *PLoS Genet* 9:e1003708.
- Parker JS et al. (2009) Supervised risk predictor of breast cancer based on intrinsic subtypes. *J Clin Oncol* 27:1160.
- Southby J, Gooding C, Smith CWJ (1999) Polypyrimidine tract binding protein functions as a repressor to regulate alternative splicing of α-actinin mutually exclusive exons. *Mol Cell Biol* 19:2699.
- Vanharanta S et al. (2014) Loss of the multifunctional RNA-binding protein RBM47 as a source of selectable metastatic traits in breast cancer. *eLife* 3:e02734.
- Warzecha CC et al. (2009) ESRP1 and ESRP2 are epithelial cell-type-specific regulators of FGFR2 splicing. *Mol Cell* 33:591.
- Warzecha CC et al. (2010) An ESRP-regulated splicing programme is abrogated during the epithelial–mesenchymal transition. *EMBO J* 29:3286.
- Zhang Y et al. (2017) Ribosomal proteins Rpl22 and Rpl22l1 control morphogenesis by regulating pre-mRNA splicing. *Cell Rep* 18:545.

The ten citations outside the TCGA marker papers were checked against PubMed, PMC or publisher
records on 2026-10-12; two were corrected (the title of Eisenhofer 2004, and the year and journal
of Gooding 1994). The five added with §9 (Gooding 1998 and 2013, Vanharanta 2014, Warzecha 2009
and 2010) were checked against the same kinds of record; re-reading the two Gooding abstracts is
what exposed the R6 misreading (A-5).
