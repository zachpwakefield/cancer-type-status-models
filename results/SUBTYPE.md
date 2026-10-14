# Subtype within a cancer: what splicing reads

One classifier per cohort and modality: 22 subtype cohorts in 15 TCGA projects, 9 modalities each,
198 models over 7,552 patient–cohort pairs (4,618 patients, one primary tumour each). Every number
here is out-of-fold — each patient was scored by models that never saw it.

The other two tasks ask whether a sample is tumour ([TUMOUR_STATUS.md](TUMOUR_STATUS.md)) and which
tissue it came from ([CANCER_TYPE.md](CANCER_TYPE.md)). This one asks, inside one cancer, which
kind of that cancer a tumour is — and whether splicing knows. How the three tasks' results relate
is in [ACROSS_TASKS.md](ACROSS_TASKS.md).

- **Methods and decisions:** [docs/SUBTYPE_METHODS.md](../docs/SUBTYPE_METHODS.md) (D-30 to D-47
  and D-52 to D-69, amendments A-1 to A-7 and A-9 to A-15). Every label map, cohort and threshold was
  fixed in [config/subtype.yaml](../config/subtype.yaml) before the first model was fitted; the
  regulator, site, MSI-repeat and BRCA analyses of §9–§12 were fixed in
  [config/followup.yaml](../config/followup.yaml) before they were run.
- **Scripts:** `30`–`35` build the labels, matrices, models and tables; `37` and `38` the
  regulator and site analyses; `42`, `43` and `45` the MSI-repeat analysis, its robustness checks and
  the DNA-read genotyping; `50` and `51` its follow-ups and the sequence code; `46`–`49` the BRCA
  follow-ups; `36` draws these figures.
- **Figures:** U1–U22 in [figures/subtype/](figures/subtype/). Each has a caption file stating the
  question, the method, n and the test, and a CSV of the plotted values. Each figure below comes
  with a short guide to reading it, then its legend (copied from the caption file), then the
  findings.
- **Tables:** `results/subtype/summary/` (performance, baselines, recall, fusion tests),
  `results/subtype/interpret/` (markers, gene overlap, event audit, cross-family tables),
  `results/subtype/regulators/`, `results/subtype/site_probe/`, `results/subtype/msi_cis/`,
  `results/subtype/brca/` and one
  folder per (view, cohort) under `results/subtype/<view>/<cohort>/`.
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

A third follow-up explains §8. MSI tumours cannot repair slippage at mononucleotide repeats, and
the exons their models chose sit, far more often than chance, right after a long poly-T tract in
the intron: 40 and 34 of the top 50 in stomach and colon, against 5% of exons. A score built from
the reference genome alone, with no model, ranks MSI tumours in both (AUROC 0.71, 0.78). It works
only for tracts at the acceptor, and the same exons move in all three organs. The longest tracts
(≥ 20 nt) move the other way, toward inclusion. That reversal was seen after the fact and
replicated only in part on alternative 3′ splice sites, and the mechanism we wrote down for it
did not hold (§11).

Followed further, with each step written down before it was computed, the reversal has a cause in
the sequence, and the cause makes a test. Where the branch point lies far from the AG, a repeat
exon gains inclusion in MSI. Where it lies close, the exon is skipped. The relation holds beyond
the run's own extent, is absent from exons without a tract, and is absent from
mismatch-repair-proficient contrasts. A two-part code follows: runs of 10–19 nt fall; runs ≥ 20 nt
with a distant branch point rise.
- **Endometrium.** The code separates endometrial MSI (AUROC 0.84, against 0.51 for the one-sided
  score). The rule was partly read from these cohorts, so that figure is in-sample.
- **Thirty cancers never used.** Unchanged and unfitted, the code finds the 50 tumours that
  carry a DNA signature of mismatch-repair deficiency among 7,685: AUROC 0.83. It does so without
  the gastrointestinal and uterine relatives (0.81), and against hypermutated proficient tumours
  (0.85). There, the one-sided score fails (0.55).
- **The exons and the dose.** The same exons move (ρ 0.51), and within those 50 tumours the
  signal deepens with the deletion count (§11).
- **Not infiltration.** MSI tumours are heavily infiltrated, yet adjusting for immune and stromal
  expression and an RNA-handling proxy keeps all of the effect. The code separates as well in
  the less-infiltrated half of each cancer (0.84).
- **Other MSI callers.** Against MANTIS and MSIsensor, which call MSI from the exome reads by other
  means, the code separates MSI-high from stable tumours at 0.97 and 0.89.

A fourth follow-up brings DNA into breast cancer (§12). Splicing reads HER2 without the amplicon.
It separates ERBB2-amplified luminal tumours about as well as expression does, and setting the
amplicon's inputs to their mean leaves it unchanged, while expression loses 0.03. An apparent lead
over expression did not survive a copy-number call without ABSOLUTE's ceiling. Splice-site mutations skip their exons, and
SF3B1-hotspot tumours use cryptic acceptors, as they should. The events that looked inherited
mostly are not. In tumour–normal pairs, the models' top "genotype-like" events often change group
between a patient's normal tissue and tumour: they are subtype programmes.

The neighbour-gene audit keeps this honest. Of 2,248 top events, 97 follow an overlapping gene
more than their own — 62 of them in HIT — and four of the nineteen strongest single events are
such readouts; the leiomyosarcoma first exon that looks like an NDE1 promoter is MYH11
transcription read on the other strand (§7).

## Summary

1. **Splicing matches expression on lineage and genome-state subtypes** (Figure U1): ESCA, CESC, TGCT,
   LGG_IDH, STAD_TCGA and BRCA_ER within 0.015 balanced accuracy; LGG histology +0.026 for
   splicing (CI −0.002 to +0.054). In AUROC the AS ensemble equals or beats expression in 7 of 22
   cohorts.
2. **It falls behind where one gene defines the class**: PCPG site −0.174 (CI −0.267 to −0.089;
   McNemar p = 1.1e-5, 21 against 1 discordant patients), PRAD ETS class −0.145 (p = 2.5e-5),
   UCEC molecular class −0.117 (CI −0.205 to −0.029). The losses concentrate in small classes
   (recall: POLE 0.54 → 0.08, myxofibrosarcoma 0.67 → 0.29, ETV1 0.93 → 0.57; Figure U2).
3. **Fixed 1:1 fusion adds little**: between −0.042 and +0.036 in every cohort; CI above zero only
   for LGG (+0.023) and UCEC_TCGA (+0.036), with McNemar p 0.064 and 0.070.
4. **Histology errors are molecularly coherent** (Figure U3): codeletion in glioma, copy-number class in
   endometrial cancer, genomic stability in gastric cancer, pseudohypoxia and PNMT loss in PCPG.
5. **The literature's markers are in the data and in the expression models, rarely in the
   splicing models** (Figure U4), and the two modalities' top-50 genes share a median of one gene (Figure U5).
6. **Named splicing biology, audited** (Figure U6): smooth-muscle ACTN1/TPM1 exons in leiomyosarcoma; a
   CDH1 skipping tail in lobular breast cancer; the p16INK4a first exon in HPV-positive HNSC; the
   NDE1 "promoter" of leiomyosarcoma is MYH11; 97 of 2,248 top events read a neighbour.
7. **MSI leaves the same splicing marks in three cancers** (Figure U7, exploratory): 10 events separate
   MSI from the rest in UCEC, STAD and COAD alike (AUROC ≥ 0.75 or ≤ 0.25 in all three), 7 of
   them chosen by at least two of the three models independently.
8. **Regulators** (Figure U8; R1–R6 fixed in advance): ESRP2 lower in biphasic mesothelioma (AUROC
   0.82), ESRP1 lower in mesenchymal HNSC (0.76), PTBP1 lower in leiomyosarcoma (0.66) and moving
   the ACTN1 smooth-muscle exon against the leiomyosarcoma pattern (partial ρ +0.43). ESRP1 in
   mesothelioma is not resolved, and TPM1 runs against the published mechanism (A-5).
9. **Site** (Figure U9): holding tissue-source sites out costs GEX ≤ 0.035 in 20 of 21 cohorts, but the
   splicing views 0.056–0.098 in HNSC_HPV, CESC, THCA, HNSC_mRNA, TGCT and STAD_TCGA. "Splicing
   loses more than GEX in most cohorts" holds in 16 of 21 (q = 0.080): not resolved.
10. **MSI splicing follows acceptor poly-T repeats** (Figures U10–U11; M1–M7, N1–N3, Q1–Q2 fixed in
    advance): 80% and 68% of the STAD and COAD MSI models' top 50 SE events carry an acceptor
    T-run ≥ 10 (5% of exons); a genome-only score ranks MSI at AUROC 0.71 and 0.78 (UCEC 0.51), at
    the acceptor only; MRE11 exon 5 replicates in all three. Runs ≥ 20 move toward inclusion
    (exploratory; 2 of 6 replication tests supported, 6 of 6 in direction). Followed up (Figures U16–U17;
    D-61 to D-64 fixed in advance): an exon whose branch point lies further from the AG gains
    inclusion in MSI, beyond its run length, in all three cancers (adjusted ρ 0.11–0.16). Deletion
    dose tracks the score in STAD only, and does not explain the weak endometrial score. No NMD
    signature shows in gene expression. GTEx genotypes no variant inside the runs. The resulting
    two-arm code (Figures U18–U19; D-65 to D-67 fixed in advance) is specific to the branch point and to
    mismatch repair. It finds MMR-deficient tumours in 30 cancers never used (AUROC 0.83 against
    0.55 for the one-sided score; 0.85 against hypermutated proficient tumours), and it separates
    endometrial MSI (0.84, in-sample for the rule). It is not an infiltration readout (Figure U20;
    D-68 fixed in advance): adjusting for immune, stromal and RNA-handling covariates keeps the
    whole effect, and it separates at 0.84 in the less-infiltrated half. Against two exome MSI
    callers it separates MSI-H from MSS at 0.97 (MANTIS) and 0.89 (MSIsensor) (Figure U21; D-69).
    The rising arm (Figure U22; D-70) is held out, needs the tract and is specific to mismatch repair
    against mutation burden. It is not unique to MSI among subtype contrasts, and its cause is open:
    it follows the branch point's distance but not the AG exclusion zone, the run's side of the
    branch point or deletion dose, and PTB relief runs the wrong way.
11. **BRCA with DNA** (Figures U12–U15; D-57 to D-60 fixed in advance): the AS ensemble separates
    ERBB2-amplified luminal tumours about as well as GEX (0.81 against 0.75 without the ABSOLUTE
    ceiling; HU2 not resolved, and no gap on the 15 newly amplified tumours), and not through the
    amplicon (knockout ≤ 0.007). 110 of 157 splice-site mutations skip their exon;
    SF3B1-hotspot tumours use cryptic acceptors (AUROC 0.90). Only 3 of 45 of the models' top
    genotype-like events pass a tumour–normal germline check; on the 82 that pass overall, LOH still
    shows (L1–L3).

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

**How to read Figure U1.**
- **Rows are cohorts**, grouped by how the label was made; the tick on each row is chance.
- **Compare the blue dot with the black one.** Blue is the AS ensemble, black is GEX. A blue dot
  at or beyond the black one means splicing recovers the subtype as well as expression. The grey
  bar shows how far the eight single splicing views spread.
- **The dark-blue dot is the 1:1 fusion.** The columns at right say whether it beats GEX on the
  same patients.
- **Read the bottom block apart.** Its labels were made from mRNA clusters, so GEX is near-circular
  there (D-44) and splicing is the independent test.

![U1](figures/subtype/U1_performance.png)

<!-- legend U1: copied from figures/subtype/U1_performance.caption.md by the figure script's 'legends' step -->
> **Figure U1. Splicing matches expression on lineage and genome-state subtypes, and fusing the two adds little.**
> Out-of-fold balanced accuracy (mean per-class recall) of one model per cohort and view, 22 cohorts and 7,552
> patient-cohort pairs: the eight splicing views (grey; the bar spans them), expression (GEX, black), the AS
> ensemble (the mean of the eight splicing views' class probabilities, blue) and GEX + AS fused 1:1 (dark blue).
> Fusion weights were fixed before any result (D-40). Tick: chance, 1/classes, which is also the majority-class
> predictor's balanced accuracy. Rows: cohorts, patients in brackets, grouped by how the label was made: histology
> (weak = a poorly reproducible label, control = a lineage or site contrast), genomic or IHC (DNA, protein), and
> expression-defined (mRNA clusters, so GEX is near-circular there; D-44).
>
> Columns at right: fused minus GEX balanced accuracy, and the exact McNemar test of fused against GEX on paired
> per-patient correctness, raw p. McNemar weighs classes by their size and balanced accuracy weighs them equally, so
> the two can disagree (THCA). The CSV holds every point's 95% patient-bootstrap CI, macro F1, macro AUROC, the
> paired-bootstrap CI of the fused − GEX difference and BH q (family: the subtype cohorts, per comparison).
>
> Outside the expression-defined labels the AS ensemble is within 0.03 of GEX in 11 of 19 cohorts and above it in 4:
> LGG (+0.026), CESC (+0.011), STAD (+0.006), ESCA (+0.006). It falls more than 0.10 short in 3: PCPG (−0.17),
> PRAD_ETS (−0.15), UCEC_TCGA (−0.12). Fusing helps balanced accuracy with a CI above zero only in LGG (+0.023,
> McNemar p = 0.064) and UCEC_TCGA (+0.036, McNemar p = 0.070).
<!-- /legend -->

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

**How to read Figure U2.**
- **One row per class**, patients in brackets; the x axis is the class's out-of-fold recall.
- **Look for rows where the blue dot sits far left of the black one.** Those are the classes
  splicing loses. Then read the class size: the large losses sit in small classes.
- **(a)** holds the histology cohorts; **(b)** the genomic, IHC and expression-defined ones, in the
  order of Figure U1.

![U2](figures/subtype/U2_recall_by_class.png)

<!-- legend U2: copied from figures/subtype/U2_recall_by_class.caption.md by the figure script's 'legends' step -->
> **Figure U2. Splicing loses most on small classes defined by one programme.** Out-of-fold recall of every class
> (patients in brackets) for GEX (black), the AS ensemble (blue) and GEX + AS fused 1:1 (dark blue). (a) Histology
> cohorts. (b) Genomic and IHC cohorts, then the expression-defined ones, in the order of U1. 95% patient-bootstrap
> CIs and every single view are in the CSV and in `recall_by_class.csv`.
>
> Splicing's five largest losses against expression (recall, GEX → AS ensemble) are all in classes of fewer than 40
> patients: UCEC_TCGA POLE (24): 0.54 → 0.08; SARC MFS (24): 0.67 → 0.29; PRAD_ETS ETV1 (28): 0.93 → 0.57; PCPG
> Paraganglioma (30): 0.83 → 0.57; THCA Tall cell (35): 0.71 → 0.51. Its largest gains: UCEC_TCGA CN-high (83): 0.84
> → 0.94; THCA Classical (351): 0.56 → 0.66; SARC UPS (47): 0.64 → 0.72; STAD Intestinal (157): 0.64 → 0.71. Over
> the 12 classes with fewer than 40 patients, the AS ensemble's recall is lower than GEX's in 9.
<!-- /legend -->

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

**How to read Figure U3.**
- **Each panel is one cohort's GEX histology model.** Rows pair the true histology with the
  model's out-of-fold call; columns are the published molecular class of the same patients.
- **Compare an error row with the correct rows above and below it.** If errors were random, a
  misclassified tumour would look molecularly like its true histology. Instead it looks like the
  class it was called.
  - LGG: oligodendrogliomas called astrocytoma are mostly IDH-mutant without codeletion (28 of
    43), like true astrocytomas.
  - UCEC: endometrioid tumours called serous are mostly copy-number high (9 of 12), like true
    serous tumours (45 of 48).
- **Counts only.** No test was pre-specified (D-41).

![U3](figures/subtype/U3_histology_errors.png)

<!-- legend U3: copied from figures/subtype/U3_histology_errors.caption.md by the figure script's 'legends' step -->
> **Figure U3. Histology models err in the molecular direction.** Each panel is one histology cohort's expression
> (GEX) model. Rows: true histology → the model's out-of-fold call, with patients in brackets; columns: the
> published molecular class of the same patients (PanCanAtlas, D-43); cells: patients, coloured by the share of the
> row. Only patients with both labels are counted. Counts only, no test (D-41). The AS ensemble's table is in the
> CSV (`plotted` = False); it shows the same pattern.
>
> The labels here were assigned by other groups over the full TCGA cohort; that is a property of the label, not a
> leak, and no molecular label entered any histology model.
<!-- /legend -->

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

**How to read Figure U4.**
- **One row per marker** fixed before any model was fitted, with the class it should be higher in.
- **The dot asks whether the marker is in the data.** It is the AUROC of the gene's expression for
  its class against the other, with no model involved: 0.5 = no difference, filled = p < 0.05.
- **The two columns at right ask whether the models use it.** They give the gene's best rank in
  the GEX model, and the best rank of any of its events in a splicing model; a dash means never
  selected.
- **The pattern to look for** is a dot far from 0.5, a GEX rank, and a dash or a low rank in the
  splicing column.

![U4](figures/subtype/U4_literature_markers.png)

<!-- legend U4: copied from figures/subtype/U4_literature_markers.caption.md by the figure script's 'legends' step -->
> **Figure U4. The literature's markers are in the data and in the expression models, far less in the splicing
> models.** Each row is a marker fixed in `config/subtype.yaml` before any model was fitted (D-41, D-45; sources in
> the CSV), with the class the literature expects it to be higher in. Dot: AUROC of the gene's log2 CPM for that
> class over the other (0.5 = no difference), computed on all patients of the two classes with no model involved;
> filled = Mann-Whitney raw p < 0.05, open = not. Columns: the gene's best rank in the GEX model (over all patients
> or within either class) and the best rank of any of its events across the eight splicing models, with the view (–
> = no event of the gene was selected). For LGG and LGG_IDH, "1p/19q dosage" is the median over the arm's genes of
> each gene's log2 CPM minus its cohort median, 1p and 19q averaged; the GEX column counts the model's top 100 genes
> that lie on 1p or 19q.
>
> 39 of 40 gene markers differ in the expected direction (p < 0.05). 27 rank in the top 100 of their GEX model; 8
> have an event in the top 100 of a splicing model. LGG: 7 of the top 100 GEX genes on 1p or 19q against 9.4
> expected (hypergeometric p = 0.84); LGG_IDH: 33 of the top 100 GEX genes on 1p or 19q against 9.4 expected
> (hypergeometric p = 5.3e-11).
<!-- /legend -->

The 42 markers of Figure U4 were written into `config/subtype.yaml`, with their sources, before any model
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

**How to read Figure U5.**
- **Each cell** counts the genes in both the top 50 of one splicing view (columns) and the top 50
  of the same cohort's GEX model (rows; patients in brackets).
- **Chance is 0.2–1.6 genes per cell.** Most cells hold 0–2: within a cohort, the two modalities
  lean on different genes. The few dark cells are named in the legend.

![U5](figures/subtype/U5_as_vs_gex_genes.png)

<!-- legend U5: copied from figures/subtype/U5_as_vs_gex_genes.caption.md by the figure script's 'legends' step -->
> **Figure U5. The splicing and expression models of a cohort read different genes.** For each cohort (patients in
> brackets) and splicing view: how many genes are in both the top 50 genes of the splicing model (a gene ranked by
> its best event) and the top 50 of the cohort's GEX model, counted within the genes both views could select in that
> cohort (D-41). Under independent draws the expected count is 0.2–1.6 (median 0.5); the CSV holds it with the fold
> enrichment, the hypergeometric raw p and BH q (family named in the CSV).
>
> The median overlap is 1 gene of 50 and 55 of 176 pairs share none; 7 pairs exceed chance at q < 0.05. The largest
> overlaps are ESCA RI: 8 (ATP7B, ERN2, HNF1A, IKBKB-DT, IL17RB, KIF12, LGALS4, MYO7B); HNSC_HPV RI: 6 (CBARP,
> KLHL35, MIR9-3HG, MSL3, STAG3, TMEM161B-DT); SARC MXE: 6 (ACTN1, LDB3, PDLIM7, SNAP25, TPM1, TPM2).
<!-- /legend -->

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

**How to read Figure U6.**
- **One row per cohort**: the single splicing event with the largest class-specific importance,
  labelled by host gene and view.
- **The dots are class medians** on the event's own scale; the filled, named dot is the class the
  event ranks for. The further it sits from the hollow ones, the better the event separates that
  class.
- **The text at right is the overlap audit.** It names any other gene the event's segment
  overlaps, and correlates the event with that gene and with its host. Red means the event follows
  the neighbour rather than its host: it is a readout of another gene's expression.

![U6](figures/subtype/U6_event_audit.png)

<!-- legend U6: copied from figures/subtype/U6_event_audit.caption.md by the figure script's 'legends' step -->
> **Figure U6. Four of the nineteen strongest splicing events read a neighbouring gene.** For each histology and
> genomic/IHC cohort, the one splicing event with the largest class-specific mean |SHAP| across the eight splicing
> models (rule fixed before the audit was run, D-41), labelled by cohort, host gene and view. Dots: the event's
> median over the patients of each class, PSI or usage (0 to 1) or HIT index (−1 to 1); filled and named = the class
> the event ranks for. Columns at right: the D-27 audit. The event's segment is intersected with every other GENCODE
> v45 gene's exons (≥ 50 bp): antisense, same-strand, inside a readthrough model, or no overlap. ρ: Spearman
> correlation over the cohort's patients with the other gene's log2 CPM and with the host's. An overlapping event
> that follows the other gene (|ρ| ≥ 0.5 and above |ρ host|) is a readout of it (red). "Not measured": the other
> gene is not in the expression matrix, so the audit cannot clear the event.
>
> 4 of 19 are readouts of a neighbour (BRCA ZNF276 → FANCA, UCEC ERCC2 → KLC3, UCEC_TCGA ERCC2 → KLC3, HNSC_HPV VRK2
> → FANCL); 2 overlap a gene the audit cannot measure; 7 follow their own gene's expression with |ρ| ≥ 0.5. The full
> audit of the top 5 events per class and view, and of every literature marker's best event, is
> `results/subtype/interpret/event_audit.csv`.
<!-- /legend -->

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

**Figure U6** applies the pre-specified rule: for each histology and genomic/IHC cohort, take the one event
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
| UCEC, STAD, COAD · MSI | CCNDBP1 57-bp exon and nine more (§8) | see Figure U7 | no overlap; host ρ −0.72 to +0.25 | see §8 |

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

**How to read Figure U7.**
- **Rows are the 10 events; columns are the three cancers.** Each cell is the event's AUROC for MSI
  against the rest of its cohort, with no model involved; red = higher in MSI.
- **A dot (•) marks that the cohort's own model picked the event.** A cell without one is a
  replication in a cohort whose model did not need it.
- **The columns at right** are the overlap audit of Figure U6.

![U7](figures/subtype/U7_msi_events.png)

<!-- legend U7: copied from figures/subtype/U7_msi_events.caption.md by the figure script's 'legends' step -->
> **Figure U7. Ten MSI splicing events replicate across endometrial, gastric and colorectal cancer.** Every splicing
> event in the top 5 of the MSI class of the UCEC, STAD or COAD model (any view) whose AUROC, MSI against the rest
> of the cohort, is at least 0.25 from 0.5 in the same direction in all three cohorts: 10 events (AUROC in each
> cell; red = higher in MSI). • = in that cohort's MSI top 5; a cell without • is a replication in a cohort whose
> model did not pick the event. No model is involved in the AUROC; the patients with a value in each cell (MSI and
> rest) are in the CSV. Columns at right: the D-27 audit, as in U6. Exploratory (A-4): the step was added after the
> same CCNDBP1 exon was seen in the MSI top 5 of all three models; it then tests every event of every shared class.
>
> 7 of the 10 events were chosen independently by at least two of the three models. CCNDBP1's 57-bp exon is included
> at a median UCEC 0.20 vs 0.06; STAD 0.28 vs 0.04; COAD 0.28 vs 0.04 (MSI vs rest). All 10 replication cells have
> Mann-Whitney raw p < 7.3e-05; BH q over every cross-cohort test of every shared class is in the CSV.
<!-- /legend -->

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

**Figure U7** keeps the events that are ≥ 0.25 from 0.5 in the same direction in all three cancers. There
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
  0.65 against 0.65, STAD 0.95 against 0.91, COAD 0.93 against 0.95 (Figure U2).
- **MSI, not mutation load:** in UCEC_TCGA the POLE-ultramutated tumours are hypermutated but
  mismatch-repair proficient, and on these events they sit with the non-MSI classes (CCNDBP1 0.08
  against 0.20 in MSI and 0.06 in CN-low and CN-high; TMX3 1.00 against 0.89; RPL22L1 0.16 against
  0.02; `event_audit.csv`, 24 POLE patients). That fits a mismatch-repair route to these events
  rather than mutation load alone (interpretation).

**Is it MSI, or the immune infiltrate that comes with MSI?** (Result; `interpret/msi_event_confounds.csv`,
added under A-4 after Figure U7 was drawn.) MSI tumours are immune-rich, and in COAD the immune genes
themselves separate MSI from the rest (AUROC 0.66 for PTPRC, 0.74 for CD8A and 0.83 for CIITA).
The ten Figure U7 events do not follow the infiltrate:
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

**How to read Figure U8.**
- **Panel (a)** tests whether a factor's expression differs between subtypes as predicted. One row
  per class: dot = median, bar = interquartile range, blue = the class predicted lower.
- **Panel (b)** asks whether PTBP1 tracks its target exons within a class, not only between
  classes. Quintiles of PTBP1 are cut within leiomyosarcoma (blue) and within the other sarcomas
  (grey) separately. A slope inside a group is what a regulator would produce; a step between
  the groups is not enough.

![U8](figures/subtype/U8_regulators.png)

<!-- legend U8: copied from figures/subtype/U8_regulators.caption.md by the figure script's 'legends' step -->
> **Figure U8. Named splicing regulators differ between subtypes as predicted in four of six tests.** Six
> predictions fixed before the data were looked at by class (config/followup.yaml, D-46). (a) A factor's log2 CPM in
> each class: dot = median, bar = interquartile range, patients in brackets; blue = the class predicted lower. AUROC
> = P(predicted-higher patient > predicted-lower patient), two-sided Mann–Whitney.
>
> (b) The event's PSI by PTBP1 quintile, the quintiles cut within LMS and within the other sarcomas separately (dot
> = median, bar = interquartile range, about 20–26 patients per point). Test: partial Spearman of PTBP1 with the PSI
> over all 230 sarcomas, adjusting for class, stromal and immune scores; predicted sign = opposite to the LMS shift.
>
> Sensitivity (A-6, exploratory): adding each sample's mean MXE PSI (D-24) leaves ρ at +0.43 and −0.22; PTBP1 ranks
> 2 of 70 factors for ACTN1 and 12 of 70 for TPM1, where MBNL1 leads (MBNL1 +0.35). BH q over R1–R6 is in the CSV.
> The scan of all 70 factors in every cohort is `results/subtype/regulators/factor_by_class.csv`.
<!-- /legend -->

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

**How to read Figure U9.**
- **Panel (a)** gives each cohort's loss of balanced accuracy when the tissue-source sites of the
  test patients are kept out of training (patients and sites in brackets). 0 means the classifier
  does not lean on the site; grey dots are single splicing views, blue their median, black GEX.
- **Panel (b)** asks whether the loss is larger where site and class are confounded in the labels
  themselves (Cramér's V). It is not: the four largest splicing losses sit at modest V.

![U9](figures/subtype/U9_site_probe.png)

<!-- legend U9: copied from figures/subtype/U9_site_probe.caption.md by the figure script's 'legends' step -->
> **Figure U9. Part of the splicing models' subtype signal is tied to the site that sent the sample.** The probe of
> D-20 on the subtype cohorts (D-47): one untuned classifier per cohort and view (top-2,000 variance, median impute,
> standardise, L2 logistic regression C = 0.05, balanced), cross-validated with patients grouped and again with
> tissue-source sites held out (5 folds, 3 seeds, pooled out-of-fold calls). (a) The loss of balanced accuracy,
> grouped minus site-held-out, per cohort (patients; sites): grey = each splicing view, blue = their median, black =
> GEX. Ordered by the splicing median. SARC is not probed: myxofibrosarcoma comes from 4 sites (rule ≥ 5). (b) The
> loss against Cramér's V between site and class, computed from the labels alone.
>
> Largest splicing losses: HNSC_HPV +0.098 (GEX +0.014), CESC +0.077 (GEX −0.001), THCA +0.077 (GEX +0.116),
> HNSC_mRNA +0.073 (GEX +0.008). The probe's accuracies are those of the untuned classifier, not of the subtype
> models; per view, scheme and seed in `results/subtype/site_probe/`.
<!-- /legend -->

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

## 11. The MSI splicing marks are written in the genome: acceptor poly-T repeats

**How to read Figure U10.** The exons were chosen by genome sequence alone, before any PSI was
read by repeat.
- **Panel (a)**: how much more often the MSI model's top exons carry an acceptor T-run ≥ 10 than
  exons at large (5% in every cohort).
- **Panel (b)**: each patient's repeat score, the inclusion of exons with a run ≥ 10 relative to
  exons without one. Lower means the repeat exons are skipped; MSI patients sit lower in STAD and
  COAD, not in UCEC.
- **Panel (c)**: the same score built from each run-length bin. If the tract causes the skipping,
  the AUROC should rise with run length near the acceptor and vanish for runs deep in the intron
  (hollow).

![U10](figures/subtype/U10_msi_repeats.png)

<!-- legend U10: copied from figures/subtype/U10_msi_repeats.caption.md by the figure script's 'legends' step -->
> **Figure U10. Genome sequence alone says which exons MSI tumours skip.** Pre-specified before any PSI was read by
> repeat (config/followup.yaml: msi_cis, D-52). For every cassette exon, the longest sense-strand T-run overlapping
> intron positions −1..−50 before it was read from GRCh38; no label or model chose the exons. MSI against the other
> classes (POLE apart) of UCEC_TCGA, STAD_TCGA and COAD_TCGA.
>
> (a) Share of exons with a run ≥ 10 among all exons measured in the cohort and among the top 50 SE events of the
> cohort's MSI model (mean |SHAP| over MSI patients); Fisher exact (M6). (b) The repeat score of each patient: each
> exon's PSI as a percentile within the cohort, mean over exons with a run ≥ 10 minus the mean over exons with ≤ 5
> (weighted to the same mean-PSI deciles, which cancels the global PSI axis of D-24). Box = interquartile range,
> diamond = median, whiskers = 1.5 × IQR, dots = patients (n under each box); Mann–Whitney (M1, M5). (c) AUROC of
> the same score built from each run-length bin, and (hollow) from T-runs ≥ 10 at −201..−400, deep in the intron,
> with the acceptor ≤ 6.
>
> Stars: * p < 0.05, ** p < 0.01, *** p < 0.001, ns otherwise. q over family M (19 tests) is in the CSV with every
> verdict in `results/subtype/msi_cis/predictions.csv`. Per-patient scores are kept locally (D-12).
<!-- /legend -->

**Why.** §8 found MSI splicing marks shared by three cancers but could not say where they come
from. Mismatch-repair-deficient cells cannot correct slippage at mononucleotide repeats. A
shortened poly-T tract in an acceptor's polypyrimidine tract weakens the 3′ splice site, and the
exon after it is skipped (context): MRE11 exon 5 (Giannini et al. 2002, 2004), HSP110 exon 9
(Dorard et al. 2011), and genome-wide in colorectal cancer (Jonchère et al. 2024). D-52 asks
whether that route shows across UCEC, STAD and COAD when the exons are chosen from the reference
genome alone. Its predictions were committed before any PSI value was read by repeat class.

**What was done.** For every cassette exon in the subtype SE matrix, the longest sense-strand
T-run overlapping intron positions −1..−50 was read from GRCh38. 99.9% of exons read AG at
−2,−1, a check on the coordinates; 5% carry a run of ≥ 10. Each patient gets a repeat score with
no model and no label: each exon's PSI as a percentile within the cohort, averaged over the
repeat exons (≥ 10) minus the average over background exons (≤ 5). The background is weighted to
the repeat exons' mean-PSI deciles, which cancels the global PSI axis of D-24. MSI is compared
with every other class except POLE (result; `results/subtype/msi_cis/`).

- **The models had found the repeats** (M6, supported in all three). Of the top 50 SE events of
  each cohort's MSI model, 40 (STAD), 34 (COAD) and 10 (UCEC) carry an acceptor T-run ≥ 10,
  against 4.5–4.9% of the exons measured (odds ratio 84, 43 and 5). Of the nine SE events among
  §8's MSI top events, seven carry runs of 13–30: PTP4A2 14, CDC16 13, TMX3 15, ECT2 16, BOD1L1 24,
  ANXA7 29 and CCNDBP1 30. The other two, GEMIN2 and COA8, carry 6.
- **The genome alone ranks MSI tumours** (M1). The score is lower in MSI with AUROC 0.71 in STAD
  and 0.78 in COAD (p ≤ 1.1e-6). In UCEC it is 0.51, not resolved. The trained models reach
  0.84–0.99 (`detector.csv`), so the score is a test of mechanism, not a detector.
- **Only at the acceptor** (M4, supported in all three, UCEC included). The same score built
  from T-runs ≥ 10 deep in the intron (−201..−400) separates MSI by 0.17, 0.28 and 0.41 AUROC
  less (UCEC, STAD, COAD; bootstrap CIs exclude 0). It even runs the other way (0.34–0.43).
- **Longer runs, larger shifts** (M2, M3). Across exons with runs ≥ 6, the per-exon MSI shift
  falls as the run lengthens. That holds in STAD (partial ρ −0.055). In COAD it misses by a hair
  (q = 0.050) until the acceptor's pyrimidine content is also held fixed; then it holds in both
  (M3: −0.062 and −0.048). So the effect follows repeat length beyond the tract's pyrimidine
  richness, which is what slippage predicts. A trans change in 3′-site strength would follow
  pyrimidine richness instead (interpretation). UCEC is not resolved on either.
- **The same exons in three organs** (M7). For repeat exons measured in two cohorts, the MSI
  shifts correlate between cancers (Spearman 0.44 UCEC–STAD, 0.43 UCEC–COAD, 0.65 STAD–COAD). For
  background exons the figures are 0.25, 0.19 and 0.35.
- **Positive control** (K1, supported in all three). MRE11 exon 5, an 88-nt exon after a T11,
  whose skipping makes the known 484del88 transcript, is less included in MSI: AUROC 0.63, 0.81
  and 0.76. Median PSI is 0.93 against 1.00 in STAD and COAD.
- **POLE** (M5, not resolved). Because the UCEC score does not separate MSI itself, MSI against
  POLE (0.53) cannot test anything here.
- **How much it explains** (exploratory). Among exons significantly less included in MSI, 22%
  (UCEC), 17% (STAD) and 16% (COAD) carry a run ≥ 10, against 5% of all exons (odds ratio
  3.8–5.5). Most MSI-associated skipping is therefore not explained by a long acceptor T-run. That
  agrees with Jonchère et al., who explained only part of it in colorectal cancer.

**Endometrium.** UCEC is the weak cohort on the sample score. Yet its repeat exons shift with
those of STAD and COAD (ρ 0.44), and the acceptor window beats the deep one. Endometrial MSI
tumours show smaller allele shifts at mononucleotide repeats than colorectal and gastric ones
(Chung et al. 2023; context). Smaller shifts would mean a weaker cis effect, and UCEC also has
the fewest MSI patients (34). This reading was formed after the result and was not tested.

**How to read Figure U11.**
- **Panel (a)**: the share of repeat exons more included in MSI, by run length; 0.5 = no
  direction. Below 0.5 means skipping, above means inclusion.
- **Panel (b)**: the replication on A3SS events, where the run sits at one of two acceptors. Below
  0 means that acceptor is used less in MSI.
- **Panel (c)**: the exons with the largest shared shifts in all three cancers, each labelled with
  its run length and exon length. The increases carry the longest runs.

![U11](figures/subtype/U11_msi_long_tracts.png)

<!-- legend U11: copied from figures/subtype/U11_msi_long_tracts.caption.md by the figure script's 'legends' step -->
> **Figure U11. The longest acceptor repeats move toward inclusion in MSI.** (a) Exploratory (read after D-52): SE
> exons with an acceptor T-run ≥ 10 (the nine MSI events of §8 left out), share whose background-adjusted MSI shift
> is positive, by run length; bars = Wilson 95% interval. Runs of 10–19 lean toward skipping (share 0.38–0.49), runs
> ≥ 20 toward inclusion (0.56–0.67) in all three cancers. (b) The replication fixed after that (D-54) on A3SS
> events, which this analysis had not read: one acceptor with the run, the other ≤ 6; shift in the use of the
> acceptor carrying the run (negative = used less in MSI); marker = median, bar = interquartile range; Wilcoxon
> signed-rank against 0.
>
> (c) Exons measured in all three cohorts with the largest mean shift, one per gene: the 7 largest increases (top
> block) and the 9 largest decreases, then MRE11 exon 5, the known MSI target. Label = gene · acceptor T-run · exon
> length. The increases carry runs of 21–32 nt; the decreases 14–21 nt.
>
> The pre-specified mechanism was not supported: D-53 predicted that the shift would depend on the exon's baseline
> inclusion, as hnRNP C repression of long U-tracts would make it (Zarnack et al. 2013), and N1–N3 are not resolved
> in any cohort. What the data show is a length rule, not its cause. Whether the rMATS PSI here counts exon-body
> reads is not recorded; if it does, unspliced pre-mRNA could add to an apparent increase. Stars: * p < 0.05, ** p <
> 0.01, *** p < 0.001, ns otherwise; q in the CSV.
<!-- /legend -->

**The longest runs go the other way.** The two §8 loci that are *more* included in MSI, CCNDBP1
and ANXA7, carry the two longest runs of the set (30 and 29 nt).
- **Exploratory, after D-52.** Among all repeat exons (§8's nine left out), runs of 10–19 lean
  toward skipping in MSI: 38–49% move up. Runs ≥ 20 lean toward inclusion: 56–67% move up, in all
  three cancers (Figure U11a).
- **The extremes split by length.** The exons with the largest shared increase carry runs of
  21–32 nt (ANXA7, CCNDBP1, TBK1, CDC42EP1, ELOF1, EPN2, BOD1L1). The largest decreases carry
  14–21 nt (Figure U11c).
- **A mechanism, pre-specified and not supported** (D-53). hnRNP C binds long U-tracts and
  competes with U2AF65, repressing the acceptor (Zarnack et al. 2013; context). That predicts an
  increase mainly for exons at low baseline inclusion. N1–N3 test that dependence on baseline,
  and all nine tests are not resolved. The direction follows run length, not baseline inclusion.
- **A replication on events not read before** (D-54, Figure U11b). An A3SS event has two acceptors, so
  a run at the distal and at the proximal one predicts opposite signs. Runs of 10–19 make their
  acceptor used less in MSI; supported in STAD (median −0.017, p = 0.009, n = 85). Runs ≥ 20 make
  it used more; supported in UCEC (+0.044, p = 0.016, n = 16). The other four tests are not
  resolved, but all six medians fall in the predicted direction. The long-run sets hold only 16–19
  events per cohort.

**Reading.** In stomach and colon, most of the MSI splicing signal the models learned (§8) is a
DNA-sequence signal; in endometrium, a smaller part (interpretation). A fixed set of acceptors with
long poly-T tracts sits in every genome, mismatch-repair deficiency makes them unstable, and their
exons shift together in stomach, colon and endometrium.
- **Length sets the shift.** The longer the tract, the larger the shift.
- **The longest tracts reverse.** From about 20 nt, an MSI tumour includes the exon more rather
  than less. That reversal is a length rule we can describe (exploratory, partly replicated in
  A3SS), not one we can explain. The hnRNP C version we wrote down did not hold.

**Not the hospital, the infiltrate or RNA handling** (D-55; family S, 6 of 7 supported). Figure U9
showed that splicing models can learn the hospital, so the score was checked against the sending
site, tumour composition and RNA handling (`results/subtype/msi_cis/robustness/`).
- **Within hospital** (S1). Using only MSI and other tumours from the same tissue-source site, the
  score still separates MSI: AUROC 0.68 in STAD (13 sites, permutation p = 0.0008) and 0.81 in
  COAD (15 sites, p = 0.0001). Leaving out one site at a time keeps it within 0.69–0.73 and
  0.76–0.79.
- **Adjusted** (S2). With site, stromal and immune expression and the mean retained-intron PSI in
  the model, the MSI coefficient is −0.017 in STAD (95% CI −0.026 to −0.009) and −0.025 in COAD
  (−0.034 to −0.017). With site alone it is −0.019 and −0.028.
- **Not the infiltrate.** MSI tumours are heavily infiltrated (context), so this was the confound
  most worth ruling out. Immune expression lowers the score on its own (p ≤ 0.002 in STAD and
  COAD) but takes little of the MSI effect. In the less-infiltrated half of each cohort the AUROC
  is 0.73 (STAD, 21 MSI) and 0.75 (COAD, 12 MSI).
- **The normal tissue** (S3, not resolved). If the cause were the patient, the hospital or the
  handling, the adjacent normal would carry it too. 6 of 9 MSI patients' tumours score below their
  normals, against 25 of 47 other patients (AUROC 0.61, p = 0.29). Nine pairs cannot decide this.
- **The feature filter** (S4). Rebuilt in the main SE matrix, the score separates MSI with AUROC
  0.71 (STAD) and 0.78 (COAD). It correlates ρ = 1.000 with the original, so this rules out the
  filter and nothing more (A-10).
- **Endometrium** stays at 0.54 within hospital (outside the family).

**What this does not show.**
- **The deletions themselves.** The mechanism is inferred from the reference sequence and the MSI
  label, not from contraction alleles seen at these tracts in each tumour. MC3 adds only coding
  homopolymer deletions (D-61), a measure of mismatch-repair damage, not of the intronic tracts.
  D-56 pre-specifies reading the 785 repeats in each tumour's exome reads (G1–G6); it waits on the
  BAM slices.
- **Exon-body reads.** The repository does not record whether the rMATS PSI counts exon-body
  reads. If it does, unspliced pre-mRNA could add to an apparent increase.
- **The retained-intron check.** It was too thin to settle that: only 4 repeat exons with runs
  ≥ 20 per cohort have a matching RI event (`ri_check.csv`).

### Followed up: dose, branch point, frame and germline (D-61 to D-64, Figures U16–U17)

Four questions from [MSI_REPEATS_LITERATURE.md](../docs/MSI_REPEATS_LITERATURE.md) (Q1, Q2, Q5,
Q3), each fixed in `config/followup.yaml` before it was computed; results are in
`results/subtype/msi_cis/followup/`. Dose (Figure U16) held in one cancer and did not explain the
weak endometrial score. The branch point (Figure U17) held in all three cancers and became the
basis of the next round.

#### Is the shift a dose response? (D-61, Figure U16)

If contracting tracts cause the skipping, a tumour with more mismatch-repair damage should skip
more. MC3 has no reads at the intronic tracts, but it does count each tumour's homopolymer
deletions in the exome: deletions of A or T bases in runs of ≥ 5. That count measures how much
damage the tumour has accumulated.

**How to read Figure U16.**
- **Panel (a)**, one plot per cancer: each tumour's D-52 repeat score against its deletion count,
  on a log scale. Lower scores mean the repeat exons are less included.
  - The grey column on the left is the rest of each cohort, with almost no homopolymer deletions.
  - The coloured cloud is the MSI tumours. A dose response would make it slope down.
- **Panel (b)** compares the MSI tumours of the three cancers: their deletion counts, with each
  cancer's D-52 AUROC under its name.

![U16](figures/subtype/U16_msi_deletion_load.png)

<!-- legend U16: copied from figures/subtype/U16_msi_deletion_load.caption.md by the figure script's 'legends' step -->
> **Figure U16. Mismatch-repair damage tracks the repeat score in stomach cancer only.** Pre-specified
> (config/followup.yaml: msi_deletion_load, D-61). Homopolymer deletion = an MC3 deletion (FILTER PASS, wga or
> native_wga_mix) of A or T bases inside a run of ≥ 5 of that base, counted per tumour in the exome. Score = the
> D-52 per-tumour repeat-exon score: mean inclusion percentile of exons whose acceptor carries a T-run ≥ 10, minus
> that of matched exons without one (lower = repeat exons less included). Partial Spearman of score and log(1 +
> count), adjusted for ABSOLUTE purity, within MSI tumours.
>
> (a) One panel per cancer: score against homopolymer deletions (x on a log(1 + count) scale); MSI tumours in
> colour, the rest grey. The rest carry a median of 0 homopolymer deletions in every cancer (upper quartile ≤ 1);
> MSI tumours a median of 90 (COAD) to 147 (UCEC). (b) MSI tumours' homopolymer deletions by cancer (box =
> interquartile range, diamond = median, whiskers 1.5 × IQR, points = tumours), with each cancer's D-52 score AUROC
> for MSI against the rest (second line of each label).
>
> Dose tracks the score within STAD only. Between cancers it runs the other way: UCEC MSI tumours carry the most
> deletions and the D-52 score does not separate them from the rest (AUROC 0.51). The exome counts coding runs,
> mostly shorter than 10; the intronic tracts themselves were not measured. q over family DL (6 tests); per-tumour
> values are kept locally (D-12).
<!-- /legend -->

**What it shows** (DL1 supported; DL2–DL6 not resolved).
- **Only stomach slopes.** Partial ρ = −0.45 given purity, q = 0.005. Colon (+0.13) and
  endometrium (−0.13) are not resolved.
- **Endometrium does not lack damage** (panel b). Its MSI tumours carry the most deletions: median
  147, against 103 in STAD and 90 in COAD. Yet the D-52 score separates them least (AUROC 0.51).
- **Its deletions are not shorter either.** In every cancer, 98–99% of these deletions remove one
  base. UCEC's share of longer ones is, if anything, the highest: 0.010, against 0.007 (STAD) and
  0.000 (COAD).
- **So the weak endometrial score needed another explanation.** It is found below (Figure U19c).
- **What the counts measure.** They come from coding runs, mostly shorter than 10 nt. They measure
  damage, not the intronic tracts themselves.

#### Does the branch point set the direction? (D-63, Figure U17)

**Why.** A branch-point A cannot sit inside a pure T run, so a long run pushes the branch point away
from the AG. The AG is used best 19–23 nt downstream of the branch point (Chua and Reed 2001).
Contracting a tract that lies between a distant branch point and the AG brings the two closer.
D-63 asked whether this geometry explains why the longest runs gain inclusion in MSI (Figure U11).
The branch point was read from the reference sequence: the A of a yUnAy motif nearest the 3′ss,
upstream of the run. One was found for 83% of 1,580 exons, a median 40 nt from the AG.

**How to read Figure U17.**
- **Panel (a)** shows where the branch points lie. For runs of 10–19 nt the distance peaks just
  beyond the 19–23 nt optimum (grey band); runs ≥ 20 nt push it further out.
- **Panel (b)** is the test in raw form. Within each run class, exons are split into thirds by
  branch-point distance. The y axis is the share of each third more included in MSI; 0.5 means no
  direction. A line that rises from near to far means the distance, not only the run length,
  decides the direction.
- **Panel (c)** gives the formal tests. Filled points: Spearman of the MSI shift with the
  branch-point distance, both adjusted for run length and baseline inclusion. Hollow points: the
  same with the run's far end.

![U17](figures/subtype/U17_msi_branchpoint.png)

<!-- legend U17: copied from figures/subtype/U17_msi_branchpoint.caption.md by the figure script's 'legends' step -->
> **Figure U17. The branch point's distance from the 3′ splice site sets the direction of the MSI shift.**
> Pre-specified (config/followup.yaml: msi_branchpoint, D-63). Exons whose acceptor carries a T-run ≥ 10 (R exons,
> D-52), read from the GRCh38 primary assembly. Branch point = the A nearest the 3′ splice site that lies upstream
> of the run, within 100 nt, and matches yUnAy; found for 83% of 1,580 exons. MSI shift = D-52's background-adjusted
> difference in inclusion percentile, MSI minus the rest.
>
> (a) Distance from the branch point to the 3′ss, by run length (1,133 exons with runs of 10–19 nt, 179 with ≥ 20
> nt); grey band = 19–23 nt, the optimum for a 3′ss AG (Chua and Reed 2001). (b) Share of R exons more included in
> MSI, by tertile of that distance within run classes, per cancer; bars = Wilson 95% interval. (c) The tests:
> Spearman of the MSI shift with the branch-point distance after both are rank-residualised on run-length bins
> (10–11, 12–13, 14–15, 16–19, 20–24, ≥ 25) and mean inclusion; p by permuting the distance within run-length bins
> (1,000; the floor is 0.001). Hollow = the same with the run's far end in place of the branch point (exploratory).
>
> The run's far end carries nearly the same information as the branch point, because the branch point must lie
> beyond it; these tests cannot separate the two. The branch point is a motif call, not a measured position. q over
> family BP (4 tests).
<!-- /legend -->

**What it shows** (BP1–BP4 supported).
- **In all three cancers** (panel c), an exon whose branch point lies further from the AG gains
  inclusion in MSI beyond its run length. Adjusted ρ = 0.16 (STAD), 0.11 (COAD), 0.11 (UCEC) and
  0.16 pooled; each permutation p sits at its floor of 0.001.
- **The reversal of D-53 sits in the far-branch-point exons** (panel b; a reading of an
  exploratory table).
  - Runs ≥ 20: 48% of the nearest third move up, against 73% of the farthest.
  - Runs of 10–19: 29% against 46%.
- **What it could not separate.** The run's far end carries almost the same information (hollow
  points, ρ 0.12–0.17), because the branch point must lie beyond it. D-65 separates the two
  (Figure U19). The branch point is a motif call, not a measured position.

#### Frame, NMD and the germline (D-62, D-64; no figure)

**No NMD signature in gene expression** (D-62; NM1, NM2 and NM4 not resolved, NM3 not testable).
- **The exons.** 677 coding internal repeat exons (one per gene and cancer), 51% of them
  frameshifting.
- **The test.** When MSI skips a frameshifting exon, its gene's mRNA does not fall more than when it
  skips an in-frame one. AUROC 0.46 (STAD) and 0.64 (COAD); pooled difference in ρ +0.12,
  p = 0.13.
- **Shared exons.** Of the 12 repeat exons skipped in all three cancers, 5 are coding internal.
  Two of those are frameshifting: MRE11 exon 5, the known one, and an exon of GFM1.

**The germline experiment could not be run** (D-64; GL1–GL3 not testable).
- **No variants at the runs.** GTEx v8 significant sQTL pairs in transverse colon, stomach and
  uterus contain no variant of any kind inside any of the 1,580 runs, although 38–80 runs have one
  within 50 nt (`followup/germline/variant_coverage.csv`).
- **Why.** That call set does not genotype long homopolymers. The test needs one that does, for
  instance from long reads or a repeat genotyper.

### The sequence code: thirty other cancers, specificity, endometrium, infiltration, MSI callers, the rising arm, fitted models (D-65 to D-74, Figures U18–U26)

D-63 suggested a code: whether an MSI repeat exon falls or rises depends on its tract length and
on where its branch point lies.
- **Runs of 10–19 nt fall.** These are R_down.
- **Runs ≥ 20 nt rise when the branch point lies ≥ 36 nt from the AG.** These are R_up. 36 nt is
  the first tertile of that distance, a label-blind cut.
- **Score G** reads both arms against matched background exons. Score A, D-52's, reads every
  repeat exon as falling.

Ten sets of tests, each written down before it was computed, asked whether the code works
(`results/subtype/msi_cis/code/`):
- whether it travels to cancers never used to define it (D-66, Figure U18);
- whether it is about the branch point and about mismatch repair (D-65, Figure U19a–b);
- whether it rescues endometrium (D-67, Figure U19c);
- whether it is an infiltration readout (D-68, Figure U20);
- whether it agrees with MSI calls that do not come from MC3 (D-69, Figure U21);
- whether its rising arm is real, and why (D-70, Figure U22);
- whether a model fitted on its exons transfers to cancers never used (D-71, Figure U23);
- how that model compares with the exome MSI callers (D-72, Figure U24);
- where the model and the exome measures disagree (D-73, Figure U25);
- how it compares with expression models trained the same way (D-74, Figure U26).

#### Thirty cancers never used (D-66, Figure U18)

**The test set.**
- **Tumours.** One per patient in the 30 other TCGA projects of the pan-cancer SE matrix with MC3
  calls: 7,685.
- **DNA truth.** Mismatch-repair deficient (MMR-deficient) means ≥ 10 MC3 homopolymer deletions.
  The cut was fixed from the three training cohorts, where it puts 129 of 129 MSI tumours above it
  and 2 of 571 others. 50 tumours qualify:
  - BRCA 11, CESC 8, HNSC 5, PRAD 4, ACC 3, LUSC 3;
  - one or two in each of ten more cancers.
- **Scores.** Neither score was re-fitted. Each was computed and ranked within each cancer.

**How to read Figure U18.**
- **Panel (a).** Each row is an AUROC: how well a score ranks a cancer's MMR-deficient tumours
  below its other tumours; 0.5, the dashed line, is chance. Hollow is score A, filled is score G.
  - Top block: pooled over the 30 cancers; without the three tissue relatives of the training
    cancers; against hypermutated proficient tumours only.
  - Lower block: each cancer with ≥ 3 MMR-deficient tumours.
- **Panel (b).** Exon by exon: each dot is a repeat exon. Its MSI shift in the three training
  cancers is on the x axis, its shift in the MMR-deficient tumours of the 30 others on the y axis.
  Blue marks R_up.
- **Panel (c).** Within the 50 MMR-deficient tumours: score G against the deletion count.

![U18](figures/subtype/U18_msi_code_pancancer.png)

<!-- legend U18: copied from figures/subtype/U18_msi_code_pancancer.caption.md by the figure script's 'legends' step -->
> **Figure U18. A two-arm sequence code finds mismatch-repair-deficient tumours in 30 cancers never used to define
> it.** Pre-specified (config/followup.yaml: msi_code_pancancer, D-66) before any tumour outside UCEC, STAD and COAD
> was read for MSI. Tumours: one per patient in the 30 other TCGA projects of the pan-cancer SE matrix, with MC3
> calls (7,685). Mismatch-repair deficient (MMR-deficient) = ≥ 10 MC3 homopolymer deletions, the D-61 count; in the
> three training cohorts that cut separates MSI from the rest almost perfectly (129 of 129 against 2 of 571). 50
> tumours qualify. Hypermutated proficient = fewer than 10 such deletions and ≥ 500 SNVs (408). Score A = D-52's
> repeat-exon score; score G = the code: exons after a T-run of 10–19 nt (R_down, read as falling) against exons
> after a run ≥ 20 nt whose branch point lies ≥ 36 nt from the AG (R_up, read as rising), each against matched
> background exons. Neither score was re-fitted; each was computed within each cancer and ranked there.
>
> (a) AUROC for MMR-deficient tumours scoring lower (score A hollow, G filled): pooled over the 30 cancers, without
> the three tissue relatives of the training cohorts, against hypermutated proficient tumours only, and per cancer
> with ≥ 3 MMR-deficient tumours (MMR-deficient of all tumours in brackets). (b) Per-exon shift (MMR-deficient minus
> the rest, background-adjusted percentile) in the 30 cancers against the mean MSI shift of the three training
> cohorts; R exons measured in both (grey = R_down, blue = R_up, light grey = runs ≥ 20 in neither arm). (c)
> MMR-deficient tumours of the 30 cancers: score G percentile against their homopolymer deletions (exploratory).
>
> The DNA call is a proxy from exome calls, not a clinical MSI test. q over family PC (8 tests). Per-tumour values
> are kept locally (D-12).
<!-- /legend -->

**What it shows** (PC2–PC8 supported, PC1 not resolved; family PC).
- **Detection** (panel a). Score G is lower in MMR-deficient tumours: pooled AUROC 0.83 (PC2),
  against 0.55 for score A (PC1, not resolved). The difference is +0.28 (PC3).
- **Organ independence.** Without READ, ESCA and UCS, the tissue relatives of the training
  cohorts, the AUROC is 0.81 (PC4).
- **Mismatch repair, not mutation burden.** Against 408 hypermutated proficient tumours, mostly
  lung, bladder, skin and head-and-neck, it is 0.85 (PC5).
- **Per cancer** (exploratory). Score G gives 0.65 (PRAD) to 0.93 (CESC); score A gives
  0.40–0.72.
- **Exon by exon** (panel b).
  - The shifts in the 30 cancers follow the training cancers' (ρ = 0.51 over 1,753 repeat exons;
    PC6).
  - R_up exons rise more than R_down exons (AUROC 0.70; PC7).
  - The branch-point relation holds after adjustment (ρ = 0.10; PC8).
- **The two arms** (A-15, exploratory). The rising arm carries more of the signal (0.73) than the
  falling arm (0.61).
- **Dose** (panel c, exploratory). Score G falls further as the deletion count rises: partial
  ρ = −0.47 given purity, p = 0.0008. Stomach showed the same in Figure U16.

#### It is the branch point, and it is mismatch repair (D-65, D-67, Figure U19)

**How to read Figure U19.**
- **Panel (a)** shows where the switch lies: the share of repeat exons more included in MSI,
  against branch-point distance in 10-nt bins, for the two run classes. The grey line marks the
  rule's 36-nt cut.
- **Panel (b)** repeats Figure U17c's adjusted correlation for MSI (blue) and for a matched control
  (grey):
  - exons without a tract (B exons);
  - POLE tumours in UCEC;
  - CIN against GS tumours in STAD.

  If the relation came from the exons' sequence alone, or from hypermutation in general, the grey
  points would sit with the blue ones.
- **Panel (c)** gives the AUROC of score A and score G for MSI in each training cancer.

![U19](figures/subtype/U19_msi_code_specificity.png)

<!-- legend U19: copied from figures/subtype/U19_msi_code_specificity.caption.md by the figure script's 'legends' step -->
> **Figure U19. The code is specific to the branch point and to mismatch repair.** Pre-specified
> (config/followup.yaml: msi_code_specificity, D-65; msi_code_endometrium, D-67). Exon shift = D-52's
> background-adjusted difference in inclusion percentile. Adjusted ρ = D-63's: shift and branch-point distance
> rank-residualised on run-length bins and mean inclusion.
>
> (a) Exploratory: share of repeat exons more included in MSI (mean over ≥ 2 of the three cancers), by 10-nt bins of
> branch-point distance, bins with ≥ 10 exons; grey line = the rule's 36-nt cut for runs ≥ 20. Exons whose branch
> point lies within 30 nt of the AG move up rarely (0.12 of runs ≥ 20 at 20–29 nt, n = 24; 0.28 of runs 10–19); from
> 30–39 nt 0.73 and 0.51. (b) Adjusted ρ for MSI (blue) and the matching control (grey): repeat exons (R, acceptor
> T-run ≥ 10) and exons without a tract (B, ≤ 5) measured in ≥ 2 cancers; UCEC MSI and POLE, each against the rest;
> STAD MSI against the rest and CIN against GS. p = the test of the difference. (c) AUROC of score A and score G for
> MSI against the rest.
>
> The rule behind score G was read from these three cohorts (D-53 and U17b), so c is in-sample for the rule, and its
> STAD and COAD values are exploratory; U18 is the out-of-sample test. q over families BS (4 tests) and EG (2
> tests).
<!-- /legend -->

**What it shows** (BS1–BS4 and EG1–EG2 supported).
- **Beyond the run** (BS1). The distance from the run's far end to the branch point adds beyond
  the run itself: ρ = 0.12, q = 0.001. This answers the caveat of Figure U17.
- **No tract, no effect** (BS2, panel b). Exons without a tract show no relation between
  branch-point distance and MSI shift: 22,100 exons, ρ = −0.000, against 0.156.
- **Not in proficient tumours** (BS3, BS4, panel b). The relation is absent from hypermutated or
  aneuploid contrasts with working mismatch repair:
  - POLE against the rest in UCEC: ρ 0.03, against 0.11 for MSI (q = 0.012);
  - CIN against GS in STAD: ρ −0.09, against 0.16 (q = 0.001).
- **The switch lies near 30 nt** (panel a, exploratory).
  - Runs ≥ 20 with a branch point 20–29 nt from the AG: 13% (3 of 24) move up in MSI. From 30–39
    nt: 73%.
  - Runs of 10–19: 28–30% below 30 nt, 35–51% beyond.
- **A reading** (interpretation). A contraction that brings a distant branch point toward the AG
  helps. One that squeezes an already close branch point, below the space U2AF needs, harms.
- **Endometrium** (EG1–EG2, panel c). Score G separates UCEC MSI from the rest at AUROC 0.84, where
  score A gave 0.51: difference +0.33, bootstrap 95% CI 0.21 to 0.45.
  - **The two arms had cancelled.** The falling arm alone gives 0.60, the rising arm 0.77 in the
    opposite direction. Score A read both as falling.
  - **In-sample.** The rule was read partly from these three cohorts (D-53, Figure U17b), so this
    is in-sample at the level of the rule. Score G reaches 0.96 in STAD and COAD (exploratory).
    Figure U18 is the out-of-sample test.

#### Not an infiltration readout (D-68, Figure U20)

**Why.** MSI tumours are heavily infiltrated by immune cells, and D-55 had ruled out infiltration
for score A only.

**How to read Figure U20.**
- **Panel (a)** asks whether the separation survives where there is little infiltrate. It gives
  the AUROC of score G in all held-out tumours, then in the lower and upper halves of each cancer's
  immune score, then in the less-infiltrated half of each training cancer. If G read infiltration,
  the lower-half AUROC would fall toward 0.5. The grey point shows how well the immune score alone
  marks MMR deficiency.
- **Panel (b)** shows the effect of MMR deficiency (or MSI) on G after adjusting for stromal and
  immune expression, the retained-intron level (an RNA-handling proxy) and purity. It is given as
  a share of the unadjusted effect: 1 means the covariates explain none of it, 0 means all.

![U20](figures/subtype/U20_msi_code_infiltration.png)

<!-- legend U20: copied from figures/subtype/U20_msi_code_infiltration.caption.md by the figure script's 'legends' step -->
> **Figure U20. The code is not an infiltration readout.** Pre-specified (config/followup.yaml:
> msi_code_infiltration, D-68). Covariates as in D-55: stromal and immune = mean log2 CPM of eight stromal and eight
> immune genes (config/biology.yaml); RI = the tumour's mean retained-intron PSI, an RNA-handling proxy; purity =
> ABSOLUTE, where called. In the 30 held-out cancers (D-66) each covariate and score G are percentiles within the
> tumour's cancer.
>
> (a) AUROC for MMR-deficient tumours (MSI in the three cohorts) scoring lower on G: all held-out tumours, then the
> lower and upper halves of each cancer's immune score; grey = the immune score itself, read as higher in
> MMR-deficient tumours; then MSI against the rest in the lower-immune half of UCEC, STAD and COAD (MMR-deficient or
> MSI tumours in brackets). (b) OLS coefficient of MMR deficiency (or MSI) on score G after adjustment, as a share
> of the coefficient with cancer (or site) fixed effects only; bar = 95% HC3 interval; p of the adjusted
> coefficient.
>
> The immune score marks MMR-deficient tumours only weakly here (AUROC 0.59), and in the adjusted model more immune
> expression goes with a higher G, the opposite of the MMR effect. q over family IF (6 tests). Per-tumour values are
> kept locally (D-12).
<!-- /legend -->

**What it shows** (IF1–IF6 supported).
- **The covariates explain none of the effect** (panel b). In the 30 cancers, the MMR coefficient
  on G is −0.34 adjusted, against −0.33 unadjusted (IF1); with purity, −0.33 (IF2). In the
  training cancers, with site fixed effects, 1.01–1.07 of the effect is kept (IF4–IF6).
- **The separation holds where infiltration is low** (panel a). AUROC 0.84 in the lower-immune
  half (21 MMR-deficient tumours; IF3), and 0.85–0.98 in the lower-immune halves of UCEC, STAD and
  COAD.
- **The immune score marks MMR deficiency only weakly** (AUROC 0.59). In the adjusted model more
  immune expression goes with a *higher* G, against the MMR effect, so infiltration cannot produce
  the signal.

#### Against two exome MSI callers (D-69, Figure U21)

**Why.** Every held-out test so far took its truth from one source: MC3 homopolymer deletions.
MANTIS (Bonneville et al. 2017) and MSIsensor call MSI from the exome reads by other means. D-69
was written down before their scores were downloaded from cBioPortal (6.2 MB, step 00c).

**How to read Figure U21.**
- **Panel (a)** gives the AUROC of score G against each DNA truth in the 30 held-out cancers: D-66's
  MC3 call, then MANTIS and MSIsensor MSI-H against MSS. The lower block repeats the two callers in
  the three training cancers, in-sample for the rule.
- **Panels (b) and (c)** split the held-out tumours by D-66's call (grey proficient, blue
  deficient) and by the caller's call; "mid" lies between the caller's two cut-offs. If score G
  reads mismatch-repair deficiency, every blue box should sit low, whichever caller is used.

![U21](figures/subtype/U21_msi_code_msi_calls.png)

<!-- legend U21: copied from figures/subtype/U21_msi_code_msi_calls.caption.md by the figure script's 'legends' step -->
> **Figure U21. The code finds the tumours that two exome MSI callers call MSI-high.** Pre-specified
> (config/followup.yaml: msi_code_msi_calls, D-69) before the scores were downloaded. MANTIS (Bonneville et al.
> 2017) and MSIsensor scores of the TCGA PanCancer Atlas, fetched from cBioPortal (step 00c) and matched by sample
> barcode. Calls: MANTIS MSI-H > 0.6, MSS < 0.4; MSIsensor MSI-H > 10, MSS < 4; scores in between (mid) are kept
> apart. Tumours: D-66's 7,685 held-out tumours (30 cancers, none used to define the code); score G as there, a
> percentile within each cancer, not re-fitted.
>
> (a) AUROC for positive tumours scoring lower on G (positives in brackets): D-66's MC3 call (≥ 10 homopolymer
> deletions) against the rest, MANTIS and MSIsensor MSI-H against MSS; then, in-sample for the rule, MSI-H against
> MSS in UCEC, STAD and COAD. (b, c) Score G percentile of the held-out tumours by D-66's MC3 call and the caller's
> call (n in brackets); box = interquartile range, diamond = median, whiskers = 1.5 × IQR, dots = tumours (not drawn
> for the large MSS group); dashed line = the middle of the cancer.
>
> Agreement (exploratory): every MANTIS MSI-H tumour (17) is MC3-deficient; of the 50 MC3-deficient tumours, MANTIS
> calls 27 mid and 6 MSS, and these still score low on G. MSIsensor calls 4 MC3-proficient tumours MSI-H. q over
> family MC (3 tests). Per-tumour values are kept locally (D-12).
<!-- /legend -->

**What it shows** (MC1–MC2 supported, MC3 not resolved; family MC).
- **Both callers agree with the code.** Score G separates MANTIS MSI-H from MSS tumours at AUROC
  0.97 (17 MSI-H; MC1) and MSIsensor MSI-H at 0.89 (23; MC2). In the training cancers the same
  comparisons give 0.87–0.98 (exploratory, in-sample).
- **Not along the continuous score** (MC3). Over all 7,388 held-out tumours with a score, nearly all
  MSS, the MANTIS score shows no relation to score G (ρ = +0.01). The code marks the MSI-high state,
  not the variation of the MANTIS score below it.
- **Where the calls disagree** (exploratory). Every MANTIS MSI-H tumour is also MC3-deficient. Of
  the 50 MC3-deficient tumours, MANTIS calls 27 mid and 6 MSS, and these still score low on G
  (median percentile 0.15 and 0.22). MSIsensor calls 4 MC3-proficient tumours MSI-H; two of them
  sit above the middle of their cancer and two below.

#### Is the rising arm real? (D-70, Figure U22)

The rising arm is the surprising half of the code. Long tracts should weaken a splice site, yet
these exons are included more. D-70 tested it eight ways, fixed before any test was computed
(`results/subtype/msi_cis/code/reversal/`).

**How to read Figure U22.**
- **Panel (a)** asks whether MSI is special. Each grey dot is another pair of subtypes, drawn at 20
  against 20 tumours. The blue dot is MMR-deficient against proficient tumours of the same cancers,
  held out.
- **Panel (b)** compares exon classes. The rising exons sit above zero; exons with a far branch
  point but no tract do not.
- **Panel (c)** (exploratory) splits the long runs by where their branch point lies.
- **Panels (d) and (e)** are tumours. Panel (d) compares DNA classes; panel (e) asks whether the
  rise grows with deletion count.

![U22](figures/subtype/U22_msi_code_rising_arm.png)

<!-- legend U22: copied from figures/subtype/U22_msi_code_rising_arm.caption.md by the figure script's 'legends' step -->
> **Figure U22. The rising arm is specific to mismatch repair and needs the tract, but its mechanism is not
> settled.** Pre-specified (config/followup.yaml: msi_code_reversal, D-70) before any test was computed; panel (c)
> is exploratory (A-16), added after the tests were read. Exon shift = D-66's background-adjusted difference in
> inclusion percentile, MMR-deficient minus other tumours, pooled over the 30 held-out cancers (×100 = percentile
> points). R_up = acceptor T-run ≥ 20 nt with the sequence-predicted branch point ≥ 36 nt from the AG; R_down = runs
> of 10–19 nt; B_far = no run (≤ 5) with the branch point ≥ 36 nt; shield = a run ≥ 20 nt lying upstream of a branch
> point, so that contracting it cannot move the branch point. R_up arm = each tumour's mean R_up percentile against
> matched exons without a run, as a percentile within its cancer.
>
> (a) Mean R_up shift in 20 against 20 tumours (50 draws): every pair of subtypes in the 22 cohorts with MSI, POLE
> and MC3-deficient tumours removed (grey, absolute value; dashed = their 95th percentile), the held-out
> MMR-deficient tumours each paired with a proficient tumour of the same cancer (blue; bar = 95% of draws) and,
> in-sample for the rule, MSI in the three training cohorts (hollow). (b) Held-out exon shifts by class (n in
> brackets); box = interquartile range, diamond = median, whiskers = 1.5 × IQR, dots = exons (not drawn for the
> large B_far set); the axis is clipped to −16 to 21. (c) Exploratory: exons with a run ≥ 20 nt, share more included
> in MMR-deficient tumours by branch-point distance (n in brackets); bars = Wilson 95% interval. (d) R_up arm of
> held-out tumours by DNA class; dots = tumours (not drawn for the proficient group). (e) R_up arm against
> homopolymer deletions in MMR-deficient tumours (blue, held-out) and MSI tumours of the three cohorts (grey).
>
> Panel (c), exploratory: with the branch point 18–29 nt from the AG, 0.20 of long-run exons move up (n = 25); from
> 30–35 nt 0.73, the shape seen in the training cohorts (Figure U19a). q over family RV (8 tests). Per-tumour values
> are kept locally (D-12).
<!-- /legend -->

**What holds** (RV2, RV3, RV7 supported; family RV).
- **The rise is real and held out.** In the 30 cancers never used, 98 of 135 R_up exons are more
  included in MMR-deficient tumours (median +0.028).
- **Not a tumour-wide shift** (RV2). Adjusting for each tumour's average inclusion, stromal, immune
  and retained-intron levels leaves the MMR coefficient unchanged (+0.236 against +0.237).
  Adding PTBP1, HNRNPC, U2AF2 and MKI67 expression leaves it at +0.248 (exploratory, A-16).
- **It needs the tract** (RV3). Exons whose branch point is just as far but that have no run do
  not move (median −0.0004; AUROC 0.69).
- **Mismatch repair, not mutation burden** (RV7). Hypermutated proficient tumours sit at the middle
  of their cancers; MMR-deficient ones near the top (AUROC 0.74; 0.87 for MANTIS MSI-H).

**What does not hold.**
- **Not unique to MSI** (RV1, not resolved). The held-out rise, +3.5 percentile points, beats 55 of
  57 other subtype contrasts. Two sarcoma contrasts move the same exons as much (LMS against MFS
  4.2, against UPS 3.6). These long-tract exons respond to other biology too. MSI sits at the top
  of that range, not outside it.
- **The geometry reading is not confirmed.**
  - The AG exclusion zone, a branch-point-free bound on where the branch point can lie, does not
    predict the rise (RV4).
  - Long runs that lie upstream of a branch point, where contraction cannot move it, rise in 26 of
    42 (RV5, not resolved).
- **No dose** (RV6). Within MMR-deficient tumours the rise does not grow with the deletion count
  (ρ +0.08). Score G's dose relation (U18c) comes from the falling arm.
- **Not PTB relief** (RV8, against). Exons that PTBP1 represses rise *less*, not more.

**What the branch point does do** (exploratory, A-16, panel c). Among long runs held out, the
branch point's distance splits the direction sharply.
- **Within 30 nt of the AG:** 5 of 25 exons rise (median −0.038).
- **From 30 nt on:** 0.60–0.84 rise in every bin.

This is the same switch seen in the training cancers (U19a: 0.12 against 0.73), now in tumours
never used to find it. A contraction of a few nucleotides would take a 20–29 nt distance below the
19–23 nt optimum (Chua and Reed 2001), and a 30+ nt distance toward it. So geometry stays the best
single reading of panel (c), but RV4 and RV5 do not back it, and that leaves it open.
- **Longer runs rise more.** Among runs ≥ 20 the shift grows with run length (ρ +0.23).
- **Longer exons rise more.** Within R_up, the shift also grows with exon length (ρ +0.20; no such
  relation in B_far). If the PSI counts exon-body reads, unspliced RNA could add to it; whether
  rMATS was run with junction counts only is not recorded.

#### A fitted model on the repeat exons (D-71, Figure U23)

Score G has no fitted weight. D-71 pre-registered seven models trained on COAD, STAD and UCEC and
tested on the other 30 cancers (7,685 tumours, 50 MMR-deficient), with hyperparameters chosen by
leaving one training cancer out (`results/subtype/msi_cis/code/models/`). We show one of them
(A-19, chosen after the results): L2 logistic regression on the 1,562 repeat exons. It keeps the
sequence prior, since only exons with an acceptor run of 10 nt or more go in, but it learns a weight
per exon.

**How to read Figure U23.**
- **Panels (a) and (b)** are the held-out ROC and precision-recall curves, against score G.
  Precision matters when 0.65% of tumours are positive.
- **Panel (c)** splits the held-out AUROC by cancer.
- **Panel (d)** shows which exons carry the largest weights.

![U23](figures/subtype/U23_msi_code_models.png)

<!-- legend U23: copied from figures/subtype/U23_msi_code_models.caption.md by the figure script's 'legends' step -->
> **Figure U23. Logistic regression on repeat exons, fitted in three cancers, finds MMR deficiency in 30 others.**
> Pre-specified (config/followup.yaml: msi_code_models, D-71) before any model was fitted; A-19 restricts this
> figure to one of the seven pre-registered models, chosen after the results. Label: ≥ 10 MC3 homopolymer deletions.
> The model is L2 logistic regression (balanced class weights) on the 1,562 exons with an acceptor T-run ≥ 10 nt
> that are measured in ≥ 25 cancers; each input is the exon's inclusion percentile within the tumour's cancer minus
> the tumour's mean percentile of matched exons without a run. It was trained on all COAD, STAD and UCEC tumours
> (779; 150 MMR-deficient), with C chosen by leaving one training cancer out (C = 0.1), and tested on the 30 other
> cancers (7,685 tumours, 50 MMR-deficient, prevalence 0.0065). Score G is the unfitted code.
>
> (a) ROC held out; grey line = chance. (b) Precision against recall held out; dashed = prevalence. (c) AUROC per
> held-out cancer with ≥ 3 MMR-deficient tumours (in brackets). (d) The 20 largest standardised coefficients,
> refitted on the training cancers (gene, acceptor T-run length, branch-point distance; positive = more inclusion
> raises the odds of MMR deficiency).
>
> Coefficients by class (A-18, exploratory): median −0.0033 for runs of 10–19 nt (46% positive), +0.0144 for runs ≥
> 20 nt with the branch point ≥ 36 nt (66% positive). The same model on the two arms' exons only reached 0.888
> (0.848–0.927); on the 107 other repeat exons, 0.776 (0.713–0.839). The other pre-registered models, held-out
> AUROC: two-arm logistic 0.835, 10-bin logistic 0.868, 10-bin MLP 0.873, all-exon logistic 0.902, all-exon MLP
> 0.885 (tables in code/models/; family LR, 7 tests). Per-tumour predictions are kept locally (D-12).
<!-- /legend -->

**What it shows.**
- **It transfers.** Held-out AUROC is 0.900, against 0.832 for score G, and AUPRC 0.439 against
  0.083. Its top 10 held-out tumours are all MMR-deficient, as are 18 of its top 25.
- **Leaving each cancer out.** Trained on the other 32 cancers, it reaches 0.914 on the same 30
  cancers and 0.967 over all 33.
- **Against MANTIS labels it is near perfect.** With MANTIS MSI-H against MSS as the held-out label,
  AUROC is 0.9997 (17 against 7,321; exploratory).
- **Its weights are the code** (panel d). The largest positive weights fall on long runs:
  - with a far branch point (CCNDBP1 T30 bp40, EPN2 T32 bp43);
  - or with no branch point found (ANXA7 T29).

  The largest negative weights fall on runs of 10–19 nt (UTP6, PTP4A2, PHTF2).
- **It needs the code's exons** (A-18, exploratory). The same model on the two arms' exons alone
  reaches 0.888. On the 107 other repeat exons it reaches only 0.776.

**The pre-registered tests** (family LR) covered all seven models.
- **Supported:** ten sequence bins beat the two arms (LR2, 0.868), and their AUPRC beat score G's
  (LR7).
- **This model** was predicted to transfer worse than the bins. It did not (+0.032; LR4, not
  resolved, opposite in sign).
- **Not resolved:** the remaining tests.

The other models are in the tables.

#### Against the exome MSI callers (D-72, Figure U24)

D-69 used MANTIS and MSIsensor as outside checks on score G. D-72 compares the fitted model with
them, and with the MC3 deletion count that defines our label. None is a gold standard, so each
exome measure is judged against the tumours on which the other two agree, and the repeat-exon
model against the same tumours. Within-cancer AUROC ranks each positive only against negatives of
its own cancer. Leaving one cancer out (E2) is the main split here because it scores COAD, STAD
and UCEC too (`results/subtype/msi_cis/code/benchmark/`). D-72's pre-registered tests were on the
all-exon model; the repeat-exon numbers below are exploratory (A-19).

**How to read Figure U24.**
- **Panel (a)** is the benchmark: each exome measure, judged by the other two, beside the model.
- **Panel (b)** splits it by cancer.
- **Panel (c)** asks what happens at a threshold pooled over cancers.
- **Panel (d)** looks at the tumours on which the DNA measures disagree.

![U24](figures/subtype/U24_msi_callers_benchmark.png)

<!-- legend U24: copied from figures/subtype/U24_msi_callers_benchmark.caption.md by the figure script's 'legends' step -->
> **Figure U24. Within each cancer, the repeat-exon model ranks tumours almost as well as the exome MSI callers.**
> Pre-specified (config/followup.yaml: msi_tool_benchmark, D-72) before any RNA model other than score G was scored
> against MANTIS or MSIsensor. The pre-registered tests (BM1–BM6) were on the all-exon model; A-19 shows the
> repeat-exon model (Figure U23) instead, so its comparisons here are exploratory. No measure is a gold standard, so
> each exome measure is judged against the tumours on which the other two agree (positive: both call MSI; negative:
> both call stable; others left out), and the repeat-exon model against the same tumours. Calls: MC3 ≥ 10
> homopolymer deletions; MANTIS > 0.6 (stable < 0.4); MSIsensor > 10 (stable < 4). Tumours with both caller scores
> (8,153 of 8,464; RNA scores all). Within-cancer AUROC ranks each positive only against negatives of its own
> cancer. E2: each cancer left out of training (all 33); E1: trained on COAD, STAD and UCEC, tested on the other 30.
>
> (a) Within-cancer AUROC in E2 (filled) and E1 (hollow); bars = 95% bootstrap intervals over tumours within cancer
> and reference label, cut at 1. (b) Per cancer with ≥ 3 positives, E2, MANTIS judged by MC3 + MSIsensor (positives
> in brackets). (c) Sensitivity (filled) and PPV (hollow) when the model flags as many tumours as the caller calls
> MSI, by its own score pooled over cancers (E2; scores come from 33 models). (d) MC3 MMR-deficient tumours grouped
> by the two callers' calls, and proficient tumours called stable by both: share in the top decile of the model's
> score within their cancer (E2; n in brackets); dashed = chance.
>
> All six lower bounds lie above −0.05 (lowest −0.0049). Positives per reference, E2: 139, 143, 135; E1: 19, 17, 17.
> The exome measures score 0.999–1.000 against each other. Panel (c), sensitivity / PPV: MANTIS 0.97 / 1.00 against
> 0.82 / 0.84; MSIsensor 0.94 / 0.99 against 0.83 / 0.86; MC3 1.00 / 0.95 against 0.84 / 0.80. Panel (d), median
> within-cancer percentile: 0.94 (both), 0.88 (one), 0.86 (neither), 0.49 (proficient). Pre-registered, all-exon
> model: 4 of 6 supported (it ranks 0.001–0.002 below MANTIS and MSIsensor, q ≤ 0.013), none against; non-inferior
> in all six. Per-tumour scores are kept locally (D-12).
<!-- /legend -->

**What it shows.**
- **The exome measures agree almost perfectly** with each other (within-cancer AUROC 0.999–1.000).
- **Within a cancer, RNA alone is nearly as good.** The repeat-exon model reaches 0.998 in E2 and
  1.000 in E1 against every reference. In E2 it is 0.0014–0.0023 below the exome measure; in E1 it
  ties. Per cancer it reaches 1.000 in BRCA and COAD, 0.999 in STAD, 0.996 in CESC and 0.989 in
  UCEC.
- **Across cancers it is not calibrated** (panel c). Flagging as many tumours as a caller calls MSI,
  it finds 82–84% of the reference positives; the callers find 94–100%. Its inputs are
  percentiles within a cancer, so it needs a threshold per cancer, or a reference cohort to rank
  against.
- **Where the DNA measures disagree** (panel d). 50 tumours with ≥ 10 MC3 deletions are called MSI by
  neither caller. The model still puts 40% of them in its top decile within their cancer, against
  9% of proficient tumours. Some are likely real MMR deficiency the callers miss; the data cannot
  say how many.
- **The pre-registered tests** (BM1–BM6, all-exon model). Four were supported: it ranks
  0.001–0.002 below MANTIS and MSIsensor. It does not differ from the MC3 count, and it is
  non-inferior at the pre-set 0.05 margin in all six.

#### Where the model and the exome measures disagree (D-73, Figure U25)

The model calls a tumour MSI when its decision function is above 0. Leaving each cancer out of
training, it calls 168 tumours. D-73 asked, in tests written down before they were computed,
whether the calls no exome measure makes look like missed MMR deficiency or like false calls
(`results/subtype/msi_cis/code/rna_only/`).

**How to read Figure U25.**
- **Panel (a)** counts the tumours each combination of the four calls MSI.
- **Panel (b)** splits the calls by cancer.
- **Panels (c–e)** place every tumour by the model's score (up) and each exome measure (across).
  Model-only calls are red; exome-only calls are black.

![U25](figures/subtype/U25_msi_model_vs_callers.png)

<!-- legend U25: copied from figures/subtype/U25_msi_model_vs_callers.caption.md by the figure script's 'legends' step -->
> **Figure U25. The model agrees with the exome measures on its confident calls; its own extra calls are
> borderline.** Exploratory view of D-73 (config/followup.yaml: msi_rna_only; tests RO1–RO5 in the table). The model
> is D-71's repeat-exon logistic regression (A-19), each cancer scored by a model trained on the other 32; it calls
> MSI when its decision function is above 0. Exome calls: MC3 ≥ 10 homopolymer deletions, MANTIS > 0.6, MSIsensor >
> 10. The 8,153 tumours with both caller scores.
>
> (a) Tumours called MSI by each combination of the four (dots below; bars coloured by agreement class); 7,920
> tumours are called by none. (b) Tumours called by anyone, per cancer with ≥ 3 such tumours, by class. (c–e) The
> model's score against MC3 homopolymer deletions, the MANTIS score and the MSIsensor score (deletions and MSIsensor
> on a symmetric log axis); dashed = calling thresholds (MANTIS and MSIsensor: stable below the lower, MSI above the
> upper). Rings = three model-only tumours named in D-73: the highest score (STAD), MC3 calls in all four MMR genes
> and POLE (PAAD) and the lowest MLH1 expression (BLCA).
>
> Per-tumour points are kept locally (D-12).
<!-- /legend -->

**What it shows** (RO4 supported; RO1–RO3 and RO5 not resolved; family RO).
- **Confident calls agree.** 115 of the 168 are called by all three exome measures, and their scores
  sit well above the threshold (median 9.4).
- **The model's own 32 extra calls look mostly like false calls.**
  - They sit just above its threshold (median 1.0).
  - They are spread thinly over 18 cancers where MSI is rare.
  - Against tumours of the same cancer that no one calls, their MLH1 and other MMR-gene expression
    is not detectably lower (RO1, RO2), and their MANTIS scores and deletion counts are not
    detectably higher (RO3, RO5).
  - Their MSIsensor scores are slightly higher (RO4, AUROC 0.65), but far below the stable cut
    (median 0.16 against < 4).
- **A few are plausible missed cases** (ringed):
  - a STAD tumour one deletion short of the MC3 cut, with a between-zone MANTIS score;
  - a hypermutated PAAD tumour with MC3 calls in all four MMR genes and POLE;
  - a BLCA tumour with MLH1 at the 1st percentile but no DNA instability.
- **What the model misses.**
  - **Clear misses.** 20 tumours are called by all three exome measures; 19 of them are UCEC, where
    the code was always weakest.
  - **Borderline cases.** 38 more are called by MC3 alone. These are hypermutated, with half the
    deletions per SNV of clear MSI, and both callers put them in their between zones. They are as
    much borderline for the exome as for the model.

#### Against expression models trained the same way (D-74, Figure U26)

MSI can already be called from expression: MSIsensor-RNA (Jia et al. 2024), PreMSIm (Li et al. 2020)
and others, built and tested mainly in colon, stomach and endometrium. D-74 asked whether the
splicing model's transfer to other cancers is special. Expression models got its tumours, label,
split (train on COAD, STAD and UCEC; test on the 30 others) and tuning.

**How to read Figure U26.**
- **Panels (a, b)** give each model's held-out AUROC and AUPRC. The dashed line is the splicing
  model.
- **Panel (c)** splits the held-out MMR-deficient tumours by whether MLH1 is silenced.
- **Panel (d)** places each of those tumours by its rank within its cancer under splicing (up) and
  expression (across).
- **Panel (e)** repeats Figure U24's comparison with the exome callers for the expression model.

![U26](figures/subtype/U26_msi_splicing_vs_expression.png)

<!-- legend U26: copied from figures/subtype/U26_msi_splicing_vs_expression.caption.md by the figure script's 'legends' step -->
> **Figure U26. Expression models trained the same way transfer almost as well; splicing adds to them.**
> Pre-specified (config/followup.yaml: msi_vs_expression, D-74) before any expression model was scored on a held-out
> cancer. Every model was trained on COAD, STAD and UCEC with the same label (MC3 ≥ 10 homopolymer deletions) and
> tuned by leaving one training cancer out; each was tested on the 30 other cancers (7,685 tumours, 50
> MMR-deficient). All are L2 logistic regressions with balanced class weights except the PreMSIm recipe (5 nearest
> neighbours on its 15 genes scaled to [0, 1]) and MLH1 alone (unfitted). The two whole-transcriptome models keep
> the top K of 23,640 genes by ANOVA F on the training tumours (chosen for both: K = 500, C = 0.1).
>
> (a) Held-out AUROC and (b) AUPRC; bars = 95% bootstrap intervals over tumours within cancer and label; dashed =
> the splicing model; grey = chance and prevalence. (c) AUROC of each model for the held-out MMR-deficient tumours
> whose MLH1 is low (filled; below the 10th percentile of their cancer) or not (hollow), each against all held-out
> proficient tumours. (d) The held-out MMR-deficient tumours: rank of each within its cancer by the splicing score
> against the rank by the expression model on percentile inputs (fraction of the cancer at or above it; log axes,
> best at the top right); grey line = equality. (e) Leaving each cancer out of training (E2), AUROC within cancer
> against the tumours on which the two named exome measures agree (as Figure U24), for the splicing model and the
> expression model on percentile inputs (refitted with its E1 parameters, without inner tuning).
>
> 3 of 6 supported, none against (family EX). Held-out AUROC: splicing 0.900 (0.864–0.936); expression, log CPM
> 0.893 (0.849–0.930); percentiles 0.866 (0.809–0.916); PreMSIm genes 0.785 (0.722–0.845); PreMSIm recipe 0.702
> (0.634–0.769); MLH1 0.758 (0.694–0.818); splicing + expression 0.926 (0.898–0.951). AUPRC 0.44, 0.40, 0.43 and
> 0.53 for splicing, the two expression models and the combination; precision at recall 0.5: 0.30, 0.26, 0.39 and
> 0.68. Against MANTIS MSI-H (17 held out): splicing 1.000, expression 0.969 (log CPM) and 0.924 (percentiles).
> Panel (d), Spearman ρ = 0.36. Panel (e): splicing 0.998, expression 0.975–0.979; the exome measures score
> 0.999–1.000 against each other (Figure U24). E2 pooled AUROC, all 33 cancers: splicing 0.967, expression 0.941.
> Per-tumour points are kept locally (D-12).
<!-- /legend -->

**What it shows** (EX3, EX4, EX6 supported; EX1, EX2, EX5 not resolved; family EX).
- **Transfer is not special to splicing.** Retrained on all 23,640 genes, expression came close:
  0.893 on log CPM and 0.866 on percentiles, against 0.900 for splicing. Its largest weights went
  to MLH1, EPM2AIP1 and genes of PreMSIm's panel.
- **Small fixed panels do not transfer.** PreMSIm's 15 genes reached 0.785 refitted and 0.702 with
  PreMSIm's own recipe; MLH1 alone 0.758.
- **MLH1 silencing is common outside the three cancers too** (28 of 50 held-out deficient tumours).
  Where it is absent, splicing and expression tie (0.854 and 0.849), and MLH1 alone is at chance.
- **The two read-outs are complementary.** They rank the deficient tumours differently (ρ = 0.36).
  Together they reach 0.926 and double the precision at half recall (0.68, against 0.30 and 0.39).
- **Splicing stays closer to the exome callers within a cancer** (exploratory): 0.998 against
  0.975–0.979 for expression.

#### Reading

Mismatch-repair deficiency changes splicing by a rule written in the reference genome, not by a
tissue programme (interpretation).
- **The rule.** Tract length and branch-point distance say, exon by exon, which way an exon moves.
- **The test.** Read in RNA, the rule finds MMR-deficient tumours in organs never used to define
  it. It does not mistake mutation burden or infiltration for mismatch-repair deficiency.
- **Why the earlier score failed.** It assumed every repeat exon falls. The long tracts with
  distant branch points rise, which cancelled the signal in endometrium and in the 30 other
  cancers.
- **The rise is real; its cause is open.** It is held out, specific to mismatch repair against
  mutation burden, and needs the tract (D-70). It follows the branch point's distance, but not
  the AG exclusion zone, not the run's side of the branch point, and not deletion dose.
- **A model can use it.** Logistic regression on the repeat exons, fitted in three cancers, transfers
  to 30 others (D-71). Within a cancer it ranks tumours almost as well as the exome MSI callers
  (D-72).

**What this does not show.**
- **No clinical MSI test.** The DNA truths are exome calls (MC3, MANTIS, MSIsensor). They agree
  with each other and with the code, but no PCR panel or MMR immunohistochemistry was available.
  An RNA cohort processed elsewhere is the external test the code has not had.
- **No tract is observed.** No individual intronic tract is seen contracting (D-56 waits on the
  BAM slices).
- **Few tumours per cancer.** The per-cancer AUROCs rest on 3–11 MMR-deficient tumours each.
- **Not calibrated across cancers.** The fitted models rank within a cancer. A pooled threshold
  loses about 15% of sensitivity against the callers (D-72).

## 12. Breast cancer: what the subtype models read beyond the labels

Three BRCA follow-ups, each fixed in [config/followup.yaml](../config/followup.yaml) before its
outcome was compared with any model output (D-57 to D-59, and D-60 re-testing D-58): which HER2 the PAM50 models learned,
whether splice-site mutations leave their mark, and whether inherited splicing variation reaches the
subtypes. They bring in DNA for the first time. Three tables come from the TCGA PanCanAtlas:
ABSOLUTE purity, ploidy and allele-specific copy-number segments, and the MC3 public mutation
calls; D-60 adds the SNP6 copy-ratio segments. All are on GRCh37; each gene is placed by its own MC3 calls, since no liftover is used.
Tables are in `results/subtype/brca/`.

### Two HER2s: amplification, and the HER2-enriched programme (D-58, D-60, Figures U13, U15)

PAM50 calls a tumour Her2 from its expression, and ERBB2 amplification is a DNA event. The two
overlap only in part (context). Amplified means total copy number at ERBB2 ≥ 2 × ploidy.

**How to read Figure U13.**
- **Panel (a)**: P(Her2) by PAM50 group, split by ERBB2 amplification. If a model read the amplicon,
  amplified luminal tumours would sit with the Her2 tumours; if it read the programme, the
  non-amplified Her2 tumours would sit there.
- **Panel (b)**: the AUROC for each comparison, GEX (grey) against the AS ensemble (blue).
- **Panel (c)**: what happens to the luminal AUROC when the inputs from the amplicon are blanked in
  every stored model (hollow), against blanking random inputs of the same number (tick). A model
  that reads the amplicon loses accuracy; one that does not stays put.

![U13](figures/subtype/U13_brca_two_her2.png)

<!-- legend U13: copied from figures/subtype/U13_brca_two_her2.caption.md by the figure script's 'legends' step -->
> **Figure U13. Splicing reads HER2 through the HER2-enriched programme, not through the 17q12 amplicon.**
> Pre-specified before any probability was compared with copy number (config/followup.yaml: brca_her2, D-58).
> P(Her2) is the stored out-of-fold probability of the BRCA_PAM50 models: GEX, and the AS ensemble (mean of the 8
> splicing views). ERBB2 amplified = ABSOLUTE total copy number at ERBB2 ≥ 2 × ploidy (31 of 714 LumA/LumB and 30 of
> 54 PAM50 Her2 tumours). ABSOLUTE's copy number at ERBB2 does not exceed 7 in any tumour, so the rule catches the
> highest-level amplifications only: 34% of clinically HER2-positive tumours meet it. Clinical HER2 (IHC or ISH) is
> the pre-specified second reading.
>
> (a) P(Her2) by PAM50 group and amplification; box = interquartile range, diamond = median, whiskers = 1.5 × IQR,
> dots = tumours (n under each group). (b) AUROC of P(Her2) for each comparison, GEX (grey) and AS ensemble (blue).
> (c) In-silico knockout (exploratory, A-12): in each stored fold model the inputs from chr17:39.6–39.8 Mb (GRCh38;
> PPP1R1B to GRB7) are set to their training mean, and the LumA/LumB AUROC is recomputed (hollow); tick = 20 random
> input sets of the same size.
>
> q over family H (4 tests). Per-tumour probabilities are kept locally (D-12).
<!-- /legend -->

**First reading, with ABSOLUTE copy number (D-58, Figure U13).** 31 of 714 LumA/LumB and 30 of 54 PAM50
Her2 tumours are amplified.
- **GEX P(Her2) rises in amplified luminal tumours** (H1, supported): AUROC 0.72, q = 9.3e-5.
- **Splicing seemed to see it better** (H2, *against* the prediction): AS ensemble 0.81, 0.09
  above GEX (q = 0.016). The re-test below does not hold this up.
- **Within PAM50 Her2, neither model tells amplified from not** (H3, not resolved): GEX 0.482,
  AS 0.485.
- **The AS ensemble ranks non-amplified Her2 tumours above amplified luminal ones** (H4,
  supported): AUROC 0.87 (GEX 0.97). Both models' Her2 is the programme first, the amplicon
  second.
- **With clinical HER2** (IHC or ISH, the pre-specified second reading): among luminal tumours,
  the AS ensemble ranks HER2-positive above HER2-negative at AUROC 0.71, and GEX at 0.65.

**The ceiling, and the re-test (D-60, Figure U15).** ABSOLUTE never reports more than 7 copies at ERBB2.
With D-58's rule, a tumour of ploidy above 3.5 can therefore never be amplified. D-60 was
pre-specified before the file was downloaded. It reads the PanCanAtlas SNP6 copy ratio at ERBB2,
made absolute with ABSOLUTE purity and ploidy, and applies the same rule with no ceiling.

**How to read Figure U15.**
- **Panel (a)** repeats Figure U13b's amplification test with the uncapped copy number, for three
  groups: D-58's amplified tumours, all amplified tumours, and the newly amplified ones alone. The
  last group is the part D-58 could not see, so it is the independent test.
- **Panel (b)**: among amplified luminal tumours, P(Her2) against the amplification level. A model
  that follows the dose rises to the right; filled = amplified by D-58 too, hollow = newly
  amplified.

![U15](figures/subtype/U15_brca_her2_uncapped.png)

<!-- legend U15: copied from figures/subtype/U15_brca_her2_uncapped.caption.md by the figure script's 'legends' step -->
> **Figure U15. Without the ABSOLUTE ceiling, splicing matches expression on ERBB2 amplification but does not beat
> it.** Pre-specified before the copy-number file was downloaded (config/followup.yaml: brca_her2_uncapped, D-60).
> D-58 took amplification from ABSOLUTE, which never reports more than 7 copies at ERBB2. Here, the PanCanAtlas SNP6
> copy ratio at ERBB2 is made absolute with ABSOLUTE purity and ploidy, without a ceiling; amplified = copy number ≥
> 2 × ploidy, D-58's rule. Out-of-fold P(Her2) of the BRCA_PAM50 models; LumA and LumB tumours only.
>
> (a) AUROC of P(Her2), GEX (grey) and AS ensemble (blue), for D-58's amplified tumours (median ploidy 2.09), all
> tumours amplified by the uncapped call, and the newly amplified ones only (median ploidy 3.85). Each group is
> compared with the luminal tumours neither call makes amplified; the middle group is compared with all
> non-amplified tumours. (b) Amplified luminal tumours: P(Her2) against copy number / ploidy; filled = amplified by
> D-58 too, hollow = newly amplified.
>
> q over family HU (7 tests). Per-tumour values are kept locally (D-12).
<!-- /legend -->

**What it shows.**
- **Who was left out.** 45 luminal tumours are amplified; 15 of them are new. The new ones are
  genome-doubled (median ploidy 3.85, against 2.09 for D-58's 31) and amplified to a similar level
  (5.7 against 5.3 × ploidy). 92% of all amplified tumours sit at ABSOLUTE's ceiling of 7.
- **The calls agree where ABSOLUTE is uncapped.** Where ABSOLUTE reports fewer than 7 copies,
  the two copy numbers agree (Spearman 0.92, 880 tumours). The amplified share of clinically
  HER2-positive tumours rises from 34% to 46% (0.4% of HER2-negative).
- **GEX P(Her2) rises in amplified luminal tumours** (HU1, supported): AUROC 0.75.
- **Splicing's lead does not hold** (HU2, not resolved): AS 0.81, GEX 0.75; difference +0.06,
  q = 0.087.
- **The independent part shows no lead at all** (HU3, not resolved). On the 15 newly amplified
  tumours, AS reaches 0.82 and GEX 0.83.
- **Programme first holds** (HU5, supported): AS P(Her2) ranks non-amplified Her2 tumours above
  amplified luminal ones, AUROC 0.85.
- **Both models follow the dose.** Among amplified luminal tumours, GEX P(Her2) rises with the
  amplification level (HU6, ρ 0.47, supported). Splicing follows it about as much (ρ 0.45). The
  prediction that it would follow less (HU7) is not resolved.

**Splicing does not use the amplicon** (A-12 and D-60, exploratory).
- **The knockout.** In every stored fold model, the inputs from chr17:39.6–39.8 Mb (GRCh38,
  PPP1R1B to GRB7) were set to their training mean.
  - **GEX:** with D-58's group, the AUROC for amplified luminal tumours drops from 0.716 to 0.686;
    20 random input sets of the same size leave it at 0.716. With D-60's group it drops from
    0.753 to 0.720.
  - **Splicing:** no view moves by more than 0.007 in either group.
- **The amplicon's own events are left unused.** They do separate amplified tumours: an
  alternative last exon at ERBB2's 3′ end, where it overlaps the co-amplified MIEN1 on the other
  strand, has median PSI 0 against 0.61. The models hardly use them.

**Reading** (interpretation).
- **No better at amplification.** Splicing detects ERBB2 amplification in luminal tumours about
  as well as expression, not better. D-58's lead came from a group the ceiling had restricted to
  near-diploid tumours.
- **Different route.** Splicing gets there without the amplicon's genes. Its P(Her2) rises with
  the amplification level as GEX's does, yet loses nothing when the 17q12 inputs are removed. So
  it reads a programme downstream of HER2 whose strength follows the dose, not the amplified DNA
  itself.

**Caveat.** Even without the ceiling, the rule finds 46% of clinically HER2-positive tumours. Why
the rest fall short is not tested here. Two possibilities (context) are that the array saturates
at high copy number and that some tumours are HER2-positive by IHC with a lower-level gain.

### Splice-site mutations and SF3B1 hotspots (D-59, Figure U14)

**How to read Figure U14.**
- **Panel (a)**: for each splice-site mutation, where its exon's inclusion in the mutated tumour
  ranks among the other tumours of the cohort. With no effect the bars would follow the dashed
  line; a pile-up near 0 means the mutated exon is skipped.
- **Panel (b)**: a score of cryptic upstream acceptors, for SF3B1 hotspot tumours against the rest.
  Higher means more use of the cryptic AG.

![U14](figures/subtype/U14_brca_mutations.png)

<!-- legend U14: copied from figures/subtype/U14_brca_mutations.caption.md by the figure script's 'legends' step -->
> **Figure U14. Splice-site mutations skip their exons, and SF3B1 hotspot tumours use cryptic acceptors.**
> Pre-specified (config/followup.yaml: brca_mutations, D-59). MC3 public calls of the cohort's tumours (FILTER flags
> within PASS, wga, native_wga_mix). (a) Every splice-site SNV at +1/+2 or −1/−2 whose hit exon (by the MC3
> transcript's intron number, located in GENCODE v45) is a cassette exon of the main SE matrix: the exon's PSI
> percentile in the mutated tumour among the other cohort tumours (ties count half, so an exon near 100% inclusion
> in almost every tumour ranks at 0.5), stacked by side; dashed line = the count per bin with no effect; triangles =
> the 5 CDH1 mutations. (b) SF3B1 hotspot tumours (codons 622, 625, 626, 662, 663, 666, 700, 741, 742) against the
> rest on a cryptic A3SS score: mean within-cohort percentile of the long form of 820 A3SS events whose acceptors
> are 10–30 nt apart and whose upstream acceptor is rarely used in normals (median PSI < 0.2). Box = interquartile
> range, diamond = median, whiskers = 1.5 × IQR, dots = tumours.
>
> Stars: * p < 0.05, ** p < 0.01, *** p < 0.001. q over family C (3 tests). Per-tumour values are kept locally
> (D-12).
<!-- /legend -->

- **A mutated splice site skips its exon** (C1, supported). There are 157 MC3 SNVs at +1/+2 or
  −1/−2 of a cassette exon in the SE matrix. In 110 of them, the exon in the mutated tumour falls
  below the cohort median: median percentile 0.093, q = 8.0e-7. Donor mutations give 0.047 (68)
  and acceptor mutations 0.13 (89). This is a positive control, and it passes.
- **CDH1** (C2, not resolved). All 5 CDH1 splice-site mutations with a cassette exon put their
  tumour at about the 0.1st percentile, but with five mutations a sign test cannot reach
  p < 0.05 (p = 0.0625).
- **SF3B1 hotspots use cryptic acceptors** (C3, supported). The score covers 820 A3SS events whose
  two acceptors lie 10–30 nt apart and whose upstream acceptor is rarely used in normals. It
  separates the 14 hotspot tumours (K700E and the other HEAT-repeat codons) from the other 1,059
  at AUROC 0.90, q = 8.0e-7. All 14 are LumA or LumB. SF3B1 hotspot mutations move the 3′ splice
  site to a cryptic AG upstream (context; Darman et al. 2015, DeBoever et al. 2015).
- **The lobular CDH1 tail is not all splice-site mutation** (exploratory). Of the lobular tumours
  in MC3, 9 carry a CDH1 splice-site mutation, 79 a truncating one and 70 neither. §7 named the
  lowest decile of CDH1 exon 6 (145 bp). That decile sits at PSI 0.86 with a splice-site mutation,
  0.95 with a truncating one and 0.96 with neither, against 0.996 in ductal tumours. The
  splice-site tumours are the extreme of the tail, not all of it.

### Inherited splicing: what the pair check removed (D-57, A-13, A-14, Figure U12)

**Design.** Some splicing events behave like genotypes: PSI falls into three well-separated groups
(low, mid, high), and the mid group appears in normal tissue too. A label-blind scan found 1,594
such events in 1,110 genes. In the three BRCA cohorts, 122 of the models' 936 top-event entries
are of this kind (13%), against 2.6% of all eligible events. That raised two questions:
- does the tumour's loss of heterozygosity (LOH) show in these events (L1–L6)?
- does ancestry reach the subtypes through them (A1–A3)?

**As pre-specified, on all 1,110 events:**
- **LOH** (L1–L4 and L6 supported):
  - inside ABSOLUTE LOH segments, heterozygous events sit further from the middle: 88 of 101
    tumour–normal pairs (L1);
  - the shift grows with purity (ρ 0.22, L2);
  - it follows each tumour's predicted allelic imbalance (L3);
  - arm by arm, it tracks ABSOLUTE's LOH in 36 of 37 arms (L6) and the LumA–Basal LOH difference
    (ρ 0.38, L4).
- **16q in lobular tumours** (L5, not resolved): AUROC 0.50.
- **Ancestry**:
  - GSTM1 null (no expression; the common homozygous deletion, context) is less frequent in Black
    than in white patients (93 of 181 against 486 of 747; A1, supported);
  - adjusting the GSTM1–ER association for race made it larger (log OR −0.47 → −0.61), against
    A2's claim (A-13);
  - 360 of 1,110 genes differ in allele frequency by race, against at most 6 under permutation
    (A3, supported).

**How to read Figure U12.** The question is whether an event that looks like a genotype behaves
like one between a patient's normal tissue and tumour.
- **Panels (a, b)** cross-tabulate each normal's group with its own tumour's group. An inherited
  event keeps the diagonal, or moves from mid to low or high by loss of heterozygosity. A count in
  the opposite-homozygote corners cannot be inherited. CAST (a) has none; TMEM101 (b) has many.
- **Panel (c)**: tumour–normal agreement (Cohen's κ) for all events, those passing the A-13 check
  and those among the models' top events. A curve further left means less agreement.
- **Panel (d)**: each event's high-group frequency in LumA against Basal tumours. An inherited
  allele sits near the dashed line; a subtype programme sits off it.

![U12](figures/subtype/U12_brca_germline_check.png)

<!-- legend U12: copied from figures/subtype/U12_brca_germline_check.caption.md by the figure script's 'legends' step -->
> **Figure U12. Most of the BRCA models' genotype-like events are subtype programmes, not inherited alleles.** D-57
> (config/followup.yaml: brca_germline) called an event genotype-like when, without labels, its PSI fell into three
> well-separated groups (low, mid, high) with mid-group normals. A-13 adds a check from tumour–normal pairs, and
> A-14 (exploratory) adds pair agreement. One tumour and one normal per patient (112 pairs). An inherited event
> keeps its group in the tumour, apart from loss of heterozygosity (mid to low or high); a homozygous normal cannot
> become the opposite homozygote.
>
> (a, b) Each normal's group against its own tumour's group (counts; shade = share of the normal's row) for CAST SE,
> a germline event, and TMEM101 A5SS, the genotype-like event with the largest importance in the BRCA models. (c)
> Cumulative distribution of Cohen's κ between normal and tumour group, for all 1110 checked events (one per gene),
> the 82 that pass A-13 (opposite-homozygote flips in ≤ 5% of at least 20 homozygous pairs) and the 45 that rank
> among the models' top events (legend: all, A-13, top). (d) Share of tumours in the high group, LumA (x) against
> Basal (y), per event, coloured as in c (grey = the rest); dashed line = equal frequency. An inherited allele
> differs between subtypes only through ancestry.
<!-- /legend -->

**The check that changed the reading** (A-13; added after these results). An inherited allele
cannot turn a homozygous normal into the opposite homozygote in its own tumour. In 112
tumour–normal pairs, that happened in a median of 20% of homozygous pairs per event.
- **Few events pass.** Only 82 of 1,110 pass A-13 (≤ 5% flips over at least 20 homozygous pairs).
- **Rare alleles pass easily** (A-14, exploratory). The 82 have about half the minor-allele
  frequency of the rest, which makes a flip rare anyway. Tumour–normal agreement (Cohen's κ) is
  above 0.6 for only 20 events. These include CAST (κ 0.95: 109 of 112 pairs agree, Figure U12a) and the
  HLA-DRB1/DRB5 region.
- **Agreement is necessary, not sufficient.** Some high-κ events are interferon-induced (IFI27,
  IFI44L) or imprinted (H19). They may record a patient-level state rather than a genotype.

**The models' top genotype-like events are subtype programmes** (result).
- **They fail the check.** 3 of the 45 that could be checked pass A-13, and none has κ above 0.24.
- **TMEM101 is the clearest case.** Its A5SS event is the most important of them (Figure U12b): 13 of 59
  normals homozygous for the high group carry a low-group tumour.
- **Their frequency follows subtype.** The high-group frequency differs between LumA and Basal
  tumours by a median 0.18, against 0.03 for confirmed germline events (Figure U12d). For TMEM101 it is
  0.44 against 0.99.
- **The enrichment is not inherited signal.** The models' "genotype-like" enrichment is a
  preference for three-level tumour programmes.

**Re-run on the 82 confirmed events** (A-13 sensitivity):
- **The LOH read-out holds and is larger** (L1–L3 supported): 47 of 64 pairs; median shift
  0.071, against 0.034 on all events.
- **The arm tests cannot be run.** Only 4 arms keep at least 5 anchored events: L4 and L6 are not
  resolved, and L5 is not testable. So the arm-level results above cannot be separated from
  subtype biology.
- **Ancestry is stronger** (A3 supported): 45 of 82 confirmed genes differ by race (55%), against
  32% of all events, as inherited variation should.

**Exploratory: accuracy by race.**
- **Histology is worse for Black patients.** Balanced accuracy is lower than for white patients
  in both modalities: GEX 0.76 against 0.88, AS 0.74 against 0.85 (163 and 654 patients).
- **ER status and PAM50 are not.** They differ by ≤ 0.03.
- **Not explained here.** Nothing in this section explains the histology gap.

**Reading.** Inherited splicing variants are in the data, and a tumour's LOH is visible in them.
But the BRCA subtype models do not read them: their three-level events are programmes that change
with the subtype. The label-blind genotype rule could not tell the two apart. Only the
tumour–normal pairs could.

## What would change the reading

- **Mutation and copy-number calls joined to the events.** Three named results in §7–§8 rested on a
  mechanism this repository could not see. §12 joins MC3 calls for BRCA: the lobular CDH1 tail is
  partly splice-site mutation (its extreme) and partly not. Still open:
  - the RPL22L1 event, if RPL22 frameshift status explains it;
  - the POLE class that splicing misses (§3), if its POLE-mutant tumours share no splicing
    programme at all.

  A finding that the RPL22L1 event is unrelated to RPL22 would remove that mechanism from the
  reading.
- **Germline calls at the BRCA events.** §12 infers inheritance from tumour–normal agreement in 112
  pairs. Normal-tissue SNP calls near the 20 high-κ events that matched their PSI groups would make
  them splicing QTLs. More pairs would let L4–L6 be decided on confirmed events alone.
- **ERBB2 by ISH or sequencing depth.** D-60 removed ABSOLUTE's ceiling, and splicing's lead over
  expression went with it. The amplified group still holds under half of the clinically
  HER2-positive tumours. A tumour-level ERBB2 call from ISH ratios or exome depth would show
  whether the remaining tumours behave like these 45.
- **Tumour DNA at the repeats.** §11 infers repeat contraction from the reference sequence and the
  MSI label. Reads over these acceptor tracts showing contraction alleles in MSI tumours, with
  each tumour's exon shift tracking its own allele, would make the cis route a result; no
  contraction would remove it. For the reversal at ≥ 20 nt, a minigene shortening the CCNDBP1 or
  ANXA7 tract that raised inclusion would show the effect is real and cis; read counts split into
  junction and exon-body reads would show whether unspliced RNA contributes.
- **An MSI call that is not from MC3.** D-66 takes mismatch-repair deficiency from exome
  homopolymer deletions. Published MSI scores for TCGA (MANTIS or MSIsensor), or clinical MMR
  immunohistochemistry, would show whether the 50 tumours are MSI-high by other means. A score G
  that also followed those calls would remove the proxy caveat; one that did not would narrow
  the claim to the deletion signature itself.
- **Score G outside TCGA.** An RNA-seq cohort with MMR status (CPTAC, or cell lines with MLH1 or
  MSH2 loss), quantified the same way, is the external test the code has not had.
- **The MSI events in another cohort.** Figure U7 was found and tested within TCGA (three cancers, one
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
- **DNA for BRCA only.** §12 uses the PanCanAtlas MC3 calls, ABSOLUTE segments and SNP6 copy
  ratios (GRCh37, each gene placed by its own MC3 calls, no liftover). No mutation, methylation or copy-number data are
  joined for POLE, RPL22 or MLH1, so those mechanisms in §7–§8 are context, not result. §11 tests one mechanism from the reference genome, and its
  repeat exons are a result; the deletions behind them are inferred, not observed.
- **Unstranded libraries.** The audit can flag but not resolve antisense overlap. 65 top events
  overlap a gene the expression matrix does not hold and stay uncleared, including two of Figure U6's
  nineteen.
- **Combination rules.** The AS ensemble averages eight views with equal weight, which shrinks
  minority-class probabilities when most views are weak on that class (§3). Fusion is fixed at 1:1
  with no learned weights (D-40). The "best AS view" columns are picked after the fact, are
  descriptive only, and are not used in any test.
- **Regulators are read from expression.** Factor activity is approximated by mRNA level; the
  event association of §9 runs across classes, so it cannot separate regulation from co-variation
  with the subtype. Only R5/R6 adjust for class.
- **Exploratory additions.** The error-marker, shared-class and MSI-confound tables and Figure U7 were
  added after results had been read (A-4), as were the R5/R6 sensitivity checks (A-6). In §11,
  D-53 and D-54 were written after D-52's results; the reversal at ≥ 20 nt (Figure U11a, c) is
  exploratory, and D-54 is its partial replication. In §12, the knockout (A-12), the germline check
  and its re-run (A-13), pair κ (A-14) and accuracy by race were added after results. They are
  descriptive, and none feeds a model. A-1 to A-3 changed presentation or widened an audit only.

## References

Literature used as *context* (never to choose a feature, a model or an event):

- Bonneville R et al. (2017) Landscape of microsatellite instability across 39 cancer types. *JCO Precis Oncol* 2017:PO.17.00073.
- Brat DJ et al. (2015) Comprehensive, integrative genomic analysis of diffuse lower-grade gliomas. *N Engl J Med* 372:2481.
- Cancer Genome Atlas Network (2015) Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature* 517:576.
- Cancer Genome Atlas Research Network (2013) Integrated genomic characterization of endometrial carcinoma. *Nature* 497:67.
- Cancer Genome Atlas Research Network (2014) Comprehensive molecular characterization of gastric adenocarcinoma. *Nature* 513:202.
- Cancer Genome Atlas Research Network (2014) Integrated genomic characterization of papillary thyroid carcinoma. *Cell* 159:676.
- Cancer Genome Atlas Research Network (2015) The molecular taxonomy of primary prostate cancer. *Cell* 163:1011.
- Cancer Genome Atlas Research Network (2017) Comprehensive and integrated genomic characterization of adult soft tissue sarcomas. *Cell* 171:950.
- Cancer Genome Atlas Research Network (2017) Integrated genomic characterization of oesophageal carcinoma. *Nature* 541:169.
- Ceccarelli M et al. (2016) Molecular profiling reveals biologically discrete subsets and pathways of progression in diffuse glioma. *Cell* 164:550.
- Chu CS et al. (1993) Genetic basis of variable exon 9 skipping in cystic fibrosis transmembrane conductance regulator mRNA. *Nat Genet* 3:151.
- Chua K, Reed R (2001) An upstream AG determines whether a downstream AG is selected during catalytic step II of splicing. *Mol Cell Biol* 21:1509.
- Chung Y et al. (2023) Evaluation of an eight marker-panel including long mononucleotide repeat markers to detect microsatellite instability in colorectal, gastric, and endometrial cancers. *BMC Cancer* 23:1100.
- Ciriello G et al. (2015) Comprehensive molecular portraits of invasive lobular breast cancer. *Cell* 163:506.
- Collura A et al. (2014) Patients with colorectal tumors with microsatellite instability and large deletions in HSP110 T17 have improved response to 5-fluorouracil–based chemotherapy. *Gastroenterology* 146(2).
- Darman RB et al. (2015) Cancer-associated SF3B1 hotspot mutations induce cryptic 3′ splice site selection through use of a different branch point. *Cell Rep* 13:1033.
- DeBoever C et al. (2015) Transcriptome sequencing reveals potential mechanism of cryptic 3′ splice site selection in SF3B1-mutated cancers. *PLoS Comput Biol* 11:e1004105.
- Dorard C et al. (2011) Expression of a mutant HSP110 sensitizes colorectal cancer cells to chemotherapy and improves disease prognosis. *Nat Med* 17:1283.
- Eisenhofer G et al. (2004) Distinct gene expression profiles in norepinephrine- and epinephrine-producing hereditary and sporadic pheochromocytomas: activation of hypoxia-driven angiogenic pathways in von Hippel-Lindau syndrome. *Endocr Relat Cancer* 11:897.
- Fishbein L et al. (2017) Comprehensive molecular characterization of pheochromocytoma and paraganglioma. *Cancer Cell* 31:181.
- Giannini G et al. (2002) Human MRE11 is inactivated in mismatch repair-deficient cancers. *EMBO Rep* 3:248.
- Giannini G et al. (2004) Mutations of an intronic repeat induce impaired MRE11 expression in primary human cancer with microsatellite instability. *Oncogene* 23:2640.
- Gooding C et al. (1994) Smooth muscle-specific switching of α-tropomyosin mutually exclusive exon selection by specific inhibition of the strong default exon. *EMBO J* 13:3861.
- Gooding C et al. (2013) MBNL1 and PTB cooperate to repress splicing of Tpm1 exon 3. *Nucleic Acids Res* 41:4765.
- Gooding C, Roberts GC, Smith CW (1998) Role of an inhibitory pyrimidine element and polypyrimidine tract binding protein in repression of a regulated α-tropomyosin exon. *RNA* 4:85.
- Hause RJ et al. (2016) Classification and characterization of microsatellite instability across 18 cancer types. *Nat Med* 22:1342.
- Jonchère V et al. (2024) Microsatellite instability at U2AF-binding polypyrimidic tract sites perturbs alternative splicing during colorectal cancer initiation. *Genome Biol* 25:210.
- Liu J et al. (2018) An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell* 173:400.
- Liu Y et al. (2018) Comparative molecular analysis of gastrointestinal adenocarcinomas. *Cancer Cell* 33:721.
- Louis DN et al. (2016) The 2016 World Health Organization classification of tumors of the central nervous system. *Acta Neuropathol* 131:803.
- Niu B et al. (2014) MSIsensor: microsatellite instability detection using paired tumor-normal sequence data. *Bioinformatics* 30:1015.
- Novetsky AP et al. (2013) Frequent mutations in the RPL22 gene and its clinical and functional implications. *Gynecol Oncol* 128:470.
- O'Leary MN et al. (2013) The ribosomal protein Rpl22 controls ribosome composition by directly repressing expression of its own paralog, Rpl22l1. *PLoS Genet* 9:e1003708.
- Parker JS et al. (2009) Supervised risk predictor of breast cancer based on intrinsic subtypes. *J Clin Oncol* 27:1160.
- Southby J, Gooding C, Smith CWJ (1999) Polypyrimidine tract binding protein functions as a repressor to regulate alternative splicing of α-actinin mutually exclusive exons. *Mol Cell Biol* 19:2699.
- Vanharanta S et al. (2014) Loss of the multifunctional RNA-binding protein RBM47 as a source of selectable metastatic traits in breast cancer. *eLife* 3:e02734.
- Wang Y et al. (2017) Differences in microsatellite instability profiles between endometrioid and colorectal cancers: a potential cause for false-negative results? *J Mol Diagn* 19:57.
- Warzecha CC et al. (2009) ESRP1 and ESRP2 are epithelial cell-type-specific regulators of FGFR2 splicing. *Mol Cell* 33:591.
- Warzecha CC et al. (2010) An ESRP-regulated splicing programme is abrogated during the epithelial–mesenchymal transition. *EMBO J* 29:3286.
- Zarnack K et al. (2013) Direct competition between hnRNP C and U2AF65 protects the transcriptome from the exonization of Alu elements. *Cell* 152:453.
- Zhang Y et al. (2017) Ribosomal proteins Rpl22 and Rpl22l1 control morphogenesis by regulating pre-mRNA splicing. *Cell Rep* 18:545.

The ten citations outside the TCGA marker papers were checked against PubMed, PMC or publisher
records on 2026-10-12; two were corrected (the title of Eisenhofer 2004, and the year and journal
of Gooding 1994). The five added with §9 (Gooding 1998 and 2013, Vanharanta 2014, Warzecha 2009
and 2010) were checked against the same kinds of record; re-reading the two Gooding abstracts is
what exposed the R6 misreading (A-5). The five added with §11 (Chung 2023, Giannini 2002 and 2004,
Jonchère 2024, Zarnack 2013) and Dorard 2011 were checked against PubMed, PMC or publisher records
on 2026-10-12, before they were cited. Darman 2015 (PMID 26565915) and DeBoever 2015 (PMID 25768983),
added with §12, were checked against PubMed and PMC records on 2026-10-13, before they were cited.
Chu 1993, Chua and Reed 2001, Collura 2014, Hause 2016 and Wang 2017, added with the §11 follow-up,
were checked for docs/MSI_REPEATS_LITERATURE.md on 2026-10-13. Bonneville 2017 (PMID 29850653) and Niu 2014
(PMID 24371154), added with D-69, were checked against PubMed on 2026-10-13, before they were cited.
