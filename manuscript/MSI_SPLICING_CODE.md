# Mismatch-repair deficiency moves splicing by a two-arm code written in acceptor poly-T tracts

*Working draft, 13 October 2026. Every analysis is in this repository (`scripts/42`–`53`). Each was
pre-specified in [`config/followup.yaml`](../config/followup.yaml) before it was computed, with
decisions D-52 to D-70 and amendments A-9 to A-16 in
[`docs/SUBTYPE_METHODS.md`](../docs/SUBTYPE_METHODS.md). Figure files carry fuller legends with
every test; this text gives the journal version.*

---

## Abstract

**Background.** Tumours that lose DNA mismatch repair (MMR) accumulate deletions in mononucleotide
repeats, and the intronic poly-T runs that form part of many 3′ splice sites are among those
repeats. Single genes show that shortening such a tract makes the next exon skipped. Whether a
general rule governs which exons move, and in which direction, has not been shown.

**Methods.** We trained splicing-based classifiers of molecular subtypes across The Cancer Genome
Atlas (TCGA) and read which exons the microsatellite-instability (MSI) models relied on. We then
annotated the acceptor of every cassette exon from the reference genome, without labels, and ran
a chain of pre-registered tests: three training cancers (stomach, colon, endometrium), then 30
cancers never used to build the rule.

**Results.**
- **The models found the repeats.** The stomach and colon MSI models relied on exons that follow a
  poly-T run of 10 nt or more at the acceptor: 40 and 34 of their top 50 exons, against about 5%
  of all exons.
- **Two directions.** Runs of 10–19 nt were skipped in MSI tumours. Runs of 20 nt or more were
  *included more* when the sequence-predicted branch point lay far from the 3′ splice site.
- **The code travels.** A score built from this two-arm rule, with no refitting, found
  MMR-deficient tumours in 30 cancers it had never seen (AUROC 0.83; 0.85 against hypermutated
  tumours whose repair is intact). It was not explained by immune infiltration. It separated
  MSI-high tumours called by two independent exome callers (AUROC 0.97 and 0.89).
- **The rising arm is real; its cause is open.** The rise of long-tract exons held out of sample.
  It was independent of each tumour's overall inclusion level, required the tract, and was absent
  from hypermutated MMR-proficient tumours. Long-run exons fell when the branch point lay within
  30 nt of the AG and rose beyond it (exploratory). But the rise did not follow a branch-point-free
  measure of geometry, did not require the run to lie between branch point and AG, and did not
  scale with deletion load. Relief of PTB repression was excluded.

**Conclusions.** MMR deficiency moves 3′ splice-site use in both directions by rules written in the
reference genome. That gives an interpretable RNA read-out of MMR deficiency that works across
organs. Why the longest tracts gain inclusion remains to be established experimentally.

---

## Introduction

DNA mismatch repair corrects the insertion–deletion loops that DNA polymerase leaves when it slips
on short tandem repeats. Tumours that lose it, through MLH1 promoter silencing or mutation of an MMR
gene, acquire thousands of repeat-length changes: microsatellite instability (MSI). MSI is
common in colorectal, gastric and endometrial cancer and present at low frequency in many others
(Kim et al. 2013; Hause et al. 2016; Cortes-Ciriano et al. 2017; Bonneville et al. 2017). It
matters clinically because MMR-deficient tumours respond to PD-1 blockade whatever their organ
(Le et al. 2017). Every MSI-high patient is therefore worth finding, including in cancers where
MSI is rare. MSI is usually called from DNA, for instance by MANTIS or MSIsensor from exome reads
(Bonneville et al. 2017; Niu et al. 2014).

Mononucleotide runs are the repeats MMR loss destabilises most, and many sit inside introns. Just
upstream of most 3′ splice sites lies a polypyrimidine tract, bound by U2AF2, between the branch
point and the AG. When MSI shortens a poly-T run there, the splice site can weaken. Two genes made
this concrete:
- **MRE11.** Deletions in a poly(T)11 skip exon 5 and inactivate the protein (Giannini et al. 2002,
  2004), which sensitises cells to PARP inhibition (Vilar et al. 2011).
- **HSP110.** Deletions in a T17 skip exon 9 and make a dominant-negative protein (Dorard et al.
  2011); larger deletions predict a better response to chemotherapy (Collura et al. 2014).

Genome-wide, Jonchère et al. (2024) found hundreds of MSI deletions at U2AF-binding tracts in
colorectal cancer, arising early and explaining part of its exon skipping.

Two questions were left open. First, is there a rule, readable from the reference genome alone,
that says which exons move and in which direction? Second, does such a rule hold beyond the organs
where it was found? We came to these questions indirectly. We had trained per-modality classifiers
of TCGA molecular subtypes, and asked what the splicing models of MSI had learned. The answer led
to a series of tests, each fixed in writing before it was computed. This paper reports that
series, including the tests that failed. Figure 1 summarises the model the tests arrive at.

![Figure 1](../results/figures/manuscript/M1_msi_code_model.png)

**Figure 1. A two-arm model of how mismatch-repair deficiency moves 3′ splice-site use.** The
schematic uses illustrative lengths. Each panel shows one acceptor before (MSS) and after (MSI) a
contraction of its T-run (blue bar). The triangle marks the branch-point adenosine, the box the
cassette exon and the bracket the branch-point-to-AG distance. The arrow is the predicted change in
inclusion. (a) Runs of 10–19 nt: contraction weakens the tract and the exon is skipped more. (b)
Runs of 20 nt or more between a distant branch point and the AG: the tract stays long, and the
branch point moves toward the optimal distance, so the exon is included more. (c) A long run
upstream of the branch point: contraction leaves the branch-point-to-AG distance unchanged. Panels
(a) and (b) are the code tested in Figs 4–6; panel (c) is a control in Fig. 6.

---

## Results

### Splicing models of MSI rely on exons that follow acceptor poly-T tracts

Within each of 22 TCGA subtype cohorts we trained multilayer perceptrons, separately for each of
nine RNA views, with nested, patient-grouped cross-validation. We explained them with SHAP (Methods).
Three cohorts carry an MSI class:
- stomach (STAD, 55 MSI of 320);
- colon (COAD, 40 of 235);
- endometrium (UCEC, 34 of 169).

The models separated MSI from the other classes well: the splicing ensemble reached out-of-fold
AUROC 0.99 in STAD and COAD and 0.84 in UCEC.

To ask what these models read, we annotated every cassette exon from GRCh38 alone. For each exon
we took the longest sense-strand T-run overlapping intron positions −1 to −50 before it (Methods).
5% of exons carry a run of 10 nt or more. Among the top 50 exons of each MSI model, 40 carried
such a run in STAD, 34 in COAD and 10 in UCEC (odds ratios 84, 43 and 5; Fig. 2a).

![Figure 2](../results/figures/subtype/U10_msi_repeats.png)

**Figure 2. Genome sequence alone says which exons MSI tumours skip.** (a) Share of exons with an
acceptor T-run of 10 nt or more, among all exons measured and among the top 50 exons of each
cohort's MSI model; Fisher's exact test. (b) Per-patient repeat score: each exon's inclusion as a
percentile within the cohort, averaged over exons with a run of 10 nt or more, minus the average
over matched exons with 5 nt or less. Boxes show the interquartile range, diamonds the median.
Mann–Whitney test. (c) AUROC of the same score built from each run-length bin, and from runs
lying deep in the intron (hollow).

A score built from the genome annotation alone, with no model and no label, ranked MSI tumours
below the rest:
- **By organ.** AUROC 0.71 in STAD and 0.78 in COAD (p ≤ 1.1 × 10⁻⁶), but 0.51 in UCEC
  (Fig. 2b).
- **Only at the acceptor.** Runs deep in the intron (−201 to −400) gave no such separation
  (Fig. 2c).
- **The same exons across organs.** The exon-level shifts of repeat exons correlated between
  cancers (Spearman ρ 0.43–0.65, against 0.19–0.35 for background exons).
- **Positive control.** MRE11 exon 5 was less included in MSI in all three cancers (AUROC 0.63,
  0.81, 0.76).

The signal held within single hospitals (AUROC 0.68 STAD, 0.81 COAD). It also held after
adjustment for stromal and immune expression and the tumour's retained-intron level (D-55).

### Shorter tracts are skipped; the longest are included more

Stratifying by run length showed an unexpected reversal (Fig. 3a):
- **Runs of 10–19 nt** leaned toward skipping in MSI: 38–49% of exons were more included.
- **Runs of 20 nt or more** leaned toward inclusion: 56–67% more included, in all three cancers.

The exons with the largest shared increase carried runs of 21–32 nt (ANXA7, CCNDBP1, TBK1, CDC42EP1
and others). The largest decreases carried 14–21 nt (Fig. 3c).

![Figure 3](../results/figures/subtype/U11_msi_long_tracts.png)

**Figure 3. The longest acceptor repeats move toward inclusion in MSI.** (a) Share of repeat exons
more included in MSI tumours, by run length; bars show Wilson 95% intervals. (b) Replication on
alternative 3′ splice-site events, where the run lies at one of two competing acceptors. The
plot shows the change in use of the acceptor that carries the run; marker = median, bar =
interquartile range. Wilcoxon signed-rank test. (c) Exons with the largest shared increases (top)
and decreases (bottom) across the three cancers, with MRE11 exon 5. Labels give gene, run length
and exon length.

This reversal was found after the first tests had been read, so we treated it as a hypothesis.
- **hnRNP C did not explain it.** hnRNP C competes with U2AF2 on long U-tracts (Zarnack et al.
  2013). That predicts a reversal that depends on baseline inclusion. It did not (D-53; 0 of 9
  tests supported).
- **Partial replication.** On alternative 3′ splice-site events not read before, the direction
  replicated in 2 of 6 pre-registered tests, and all 6 medians lay in the predicted direction
  (D-54; Fig. 3b).

### The branch point's distance sets the direction

A long T-run sits between the branch point and the AG, so it also pushes the branch point away
from the 3′ splice site. Step II of splicing works best when the branch point lies 19–23 nt from
the AG (Chua and Reed 2001). We predicted the branch point from sequence: the nearest yUnAy
adenosine upstream of the run, found for 83% of 1,580 repeat exons. We then asked whether its
distance predicts the MSI shift beyond run length (D-63).

It did, in all three cancers. After both variables were adjusted for run-length bins and mean
inclusion, exons with a more distant branch point were more included in MSI (adjusted ρ = 0.16
STAD, 0.11 COAD, 0.11 UCEC; all q = 0.001; Fig. 4).

![Figure 4](../results/figures/subtype/U17_msi_branchpoint.png)

**Figure 4. The branch point's distance from the 3′ splice site sets the direction of the MSI
shift.** (a) Distance from the predicted branch point to the AG, by run class; the grey band is the
19–23 nt optimum (Chua and Reed 2001). (b) Share of repeat exons more included in MSI, by tertile
of that distance within each run class and cancer; Wilson 95% intervals. (c) Spearman correlation
of the MSI shift with branch-point distance after adjustment for run length and mean inclusion.
Hollow markers use the run's far end instead of the branch point.

The relation was specific (D-65; Supplementary Fig. S2):
- **Beyond the run's end.** It held for the gap between the run and the branch point, adjusted for
  where the run ends (ρ = 0.12).
- **Absent without a tract.** It was absent in exons without a tract (ρ ≈ 0).
- **Absent in proficient contrasts.** It was weaker or absent in two MMR-proficient contrasts: POLE
  against the rest in UCEC, and chromosomally unstable against genomically stable tumours in STAD.

Among long runs, inclusion switched on sharply near 30 nt. With the branch point 20–29 nt from the
AG, 12% of long-run exons were more included; at 30–39 nt, 73%.

### A two-arm code finds MMR-deficient tumours in 30 cancers never used to build it

These results define a code (Fig. 1). Runs of 10–19 nt fall in MSI (the falling arm, R_down). Runs
of 20 nt or more whose branch point lies at least 36 nt from the AG rise (the rising arm, R_up).
The 36 nt cut is the first tertile of that distance, chosen without labels. Score G subtracts the
rising arm from the falling arm, each against matched exons without a tract. The original score
(score A) treated all repeat exons as falling.

We tested the code, without refitting, in the 30 other TCGA cancers (D-66). The DNA truth came from
public exome calls (MC3; Ellrott et al. 2018):
- **MMR-deficient** = 10 or more homopolymer deletions. That cut, fixed in the three training
  cohorts, puts 129 of 129 MSI tumours and 2 of 571 others above it. In the held-out cancers, 50 of
  7,685 tumours qualified.
- **Hypermutated proficient** = fewer deletions but 500 or more single-nucleotide variants (408
  tumours).

![Figure 5](../results/figures/subtype/U18_msi_code_pancancer.png)

**Figure 5. A two-arm sequence code finds mismatch-repair-deficient tumours in 30 cancers never
used to define it.** (a) AUROC for MMR-deficient tumours scoring lower: pooled over the 30
cancers, without the three tissue relatives of the training cancers (READ, ESCA, UCS), against
hypermutated proficient tumours only, and per cancer (cancers with 3 or more MMR-deficient
tumours; counts in brackets). Hollow = score A, filled = score G. (b) Exon shifts in the held-out
cancers against the mean shift in the training cancers; blue = rising-arm exons. (c) Within
MMR-deficient tumours, score G against the number of homopolymer deletions.

Score G found MMR-deficient tumours across the held-out cancers (Fig. 5a):
- **Pooled.** AUROC 0.83, against 0.55 for score A.
- **Organ independence.** 0.81 without the training cancers' tissue relatives.
- **Not mutation burden.** 0.85 against hypermutated proficient tumours.
- **Per cancer.** From 0.65 (prostate) to 0.93 (cervix).
- **Exon by exon.** Shifts replicated the training cancers (ρ = 0.51 over 1,753 repeat exons;
  Fig. 5b), and rising-arm exons rose more than falling-arm exons (AUROC 0.70).
- **Dose.** Within MMR-deficient tumours, score G fell further as deletions rose (partial ρ = −0.47
  given purity; Fig. 5c).
- **Endometrium.** The code also separated endometrial MSI (AUROC 0.84, against 0.51 for score A),
  though this is in-sample for the rule (D-67).

Three checks addressed confounding:
- **Infiltration (D-68; Supplementary Fig. S3).** MSI tumours are heavily infiltrated, but
  adjusting for stromal and immune expression, the retained-intron level and purity kept the whole
  MMR effect. Held out, the coefficient was −0.339 adjusted against −0.335 unadjusted. In the less
  infiltrated half of each cancer, the AUROC was 0.84. The immune score itself marked MMR
  deficiency only weakly (AUROC 0.59).
- **Independent MSI calls (D-69; Supplementary Fig. S4).** Against two exome MSI callers that do
  not use MC3, score G separated MSI-high from stable tumours in the held-out cancers: AUROC 0.97
  for MANTIS (17 MSI-high) and 0.89 for MSIsensor (23).
- **The continuous score.** Score G did not follow the continuous MANTIS score among stable tumours
  (ρ = +0.01). It marks the MSI-high state, not graded instability.

### The rising arm is real, but its mechanism is not settled

The rising arm is the counter-intuitive half of the code, and it carries more of the held-out
signal (AUROC 0.73 for the rising arm alone, 0.61 for the falling arm). We therefore pre-registered
eight further tests of it (D-70; Fig. 6).

![Figure 6](../results/figures/subtype/U22_msi_code_rising_arm.png)

**Figure 6. The rising arm is specific to mismatch repair and needs the tract, but its mechanism is
not settled.** (a) Mean shift of rising-arm exons in 20 against 20 tumours, averaged over 50 draws.
Grey: every pair of subtypes in the 22 cohorts, with MSI, POLE and MMR-deficient tumours removed
(absolute values; dashed line = their 95th percentile). Blue: held-out MMR-deficient tumours, each
paired with a proficient tumour of the same cancer (bar = 95% of draws). Hollow: MSI in the three
training cancers (in-sample). (b) Held-out exon shifts by class:
- falling arm (runs 10–19 nt);
- rising arm;
- exons with an equally distant branch point but no run;
- exons whose long run lies upstream of a branch point.

(c) Exploratory: exons with runs of 20 nt or more, share more included by branch-point distance;
Wilson 95% intervals. (d) Rising-arm score of held-out tumours by DNA class. (e) Rising-arm score
against homopolymer deletions in MMR-deficient (blue) and training-cohort MSI tumours (grey).

**What held.**
- **Held out.** 98 of 135 rising-arm exons were more included in held-out MMR-deficient tumours
  (median +2.8 percentile points; Fig. 6b).
- **Not a tumour-wide shift (RV2).** Adjusting for each tumour's average inclusion, stromal and
  immune expression and the retained-intron level left the MMR effect unchanged (+0.236 against
  +0.237). Adding PTBP1, HNRNPC, U2AF2 and MKI67 expression did not reduce it either (+0.248;
  exploratory).
- **It needs the tract (RV3).** Exons whose branch point was as distant but which carried no run
  did not move (median −0.04 points; AUROC 0.69 for the difference).
- **Mismatch repair, not burden (RV7).** Hypermutated proficient tumours did not show it (AUROC 0.74
  against them; Fig. 6d). Neither did tumours called MSI-high by MANTIS alone (AUROC 0.87;
  exploratory).

**What did not hold.**
- **Not unique to MSI (RV1, not resolved).** The held-out rise (+3.5 points) exceeded 55 of 57
  other subtype contrasts, but not all of them (Fig. 6a). Two sarcoma contrasts moved the same
  exons as much: leiomyosarcoma against myxofibrosarcoma (+4.2) and against undifferentiated
  pleomorphic sarcoma (+3.6).
- **Branch-point-free geometry (RV4).** The AG exclusion zone bounds where a branch point can lie
  without any branch-point model (Gooding et al. 2006). It did not predict the rise (ρ = −0.04), and
  it was itself uncorrelated with our predicted branch point.
- **Side of the branch point (RV5).** A long run that lies upstream of a branch point cannot move
  the branch point when it contracts. Such exons rose in 26 of 42 cases, not significantly
  different from the rising arm.
- **Dose (RV6).** Within MMR-deficient tumours the rise did not grow with deletion count (ρ = +0.08;
  Fig. 6e). Score G's dose relation comes from the falling arm.
- **PTB relief (RV8, against).** Long tracts with distant branch points are the textbook
  PTB-repressed exons (Gooding et al. 2006; Llorian et al. 2010). If contraction relieved PTB, the
  most repressed exons should rise most. The opposite held (ρ = +0.20).

**What the branch point does do (exploratory, A-16).** Among long runs in the held-out cancers,
branch-point distance split the direction sharply (Fig. 6c):
- **Within 30 nt of the AG:** 5 of 25 exons were more included.
- **From 30 nt on:** 60–84% were, in every bin.

This repeats the switch first seen in the training cancers (12% against 73%) in tumours never used
to find it. The rise also grew with run length (ρ = 0.23). It grew with exon length within the
rising arm (ρ = 0.20), which matters for the caveat on exon-body reads below.

### What the data did not support

Several pre-registered tests failed or could not be run:
- **Dose by organ (D-61; Supplementary Fig. S1).** Deletion dose tracked the original score in
  stomach (ρ = −0.45) but not colon or endometrium. Endometrial MSI tumours carried the *most*
  homopolymer deletions, so deletion load does not explain their weak original score; the two-arm
  code did.
- **NMD (D-62).** No nonsense-mediated decay signature appeared in gene expression for
  frameshifting skipped exons. Two of the twelve exons skipped in all three cancers shift the
  reading frame: MRE11 exon 5 and an exon of GFM1.
- **Germline (D-64).** A germline test could not be run: GTEx v8 splicing QTLs contain no variant
  inside any of the 1,580 runs, because standard call sets do not genotype long homopolymers.

---

## Discussion

MMR deficiency shifts splicing at 3′ splice sites by rules that can be read from the reference
genome: tract length and branch-point distance say, exon by exon, which way an exon moves. Read in
RNA, these rules find MMR-deficient tumours in organs never used to define them. They do not
mistake hypermutation or immune infiltration for MMR deficiency, and they agree with two DNA-based
MSI callers.

**Relation to prior work.**
- **The falling arm** extends single-gene observations (Giannini et al. 2004; Dorard et al. 2011)
  and the colorectal study of Jonchère et al. (2024). It does so across organs, with exons chosen
  from sequence alone.
- **The rising arm is new.** A systematic search found no report that MSI makes long-tract exons
  more included. Jonchère et al. report mostly exclusion, with larger deletions giving more
  exclusion; inclusion appears only in a supplementary panel, unanalysed. The one MSI-linked
  inclusion effect we found acts in trans: RPL22 frameshifts raise MDM4 exon 6 inclusion
  (Weinstein et al. 2024).
- **Why it matters for detection.** The rising arm explains why a one-sided repeat score failed in
  endometrium and in the 30 held-out cancers: the rising exons cancelled the falling ones.

**Why do the longest tracts gain inclusion?** Our first reading was geometry. Contraction shortens
a branch-point-to-AG distance that lies beyond the 19–23 nt optimum (Chua and Reed 2001), and a
long tract stays above the length at which it saturates (Coolidge et al. 1997). Several lines of
evidence fit:
- the predicted branch point's distance sets the direction beyond run length, in all three
  training cancers and again held out;
- long-run exons with a near branch point (20–29 nt), which contraction would push below the
  optimum, fall rather than rise;
- exons with distant branch points are skipped more genome-wide (Corvelo et al. 2010; Taggart et
  al. 2017);
- germline repeat changes that move a branch point change 3′ splice-site use in the predicted
  direction: an expansion that pushes the branch point from 35 to 49 nt silences a NARS2
  pseudoexon (Hamanaka et al. 2023), and a deletion of an Alu poly(T) tail that brings a branch
  point closer activates a PTS pseudoexon (Meili et al. 2009).

But two pre-registered corollaries of geometry failed:
- The AG exclusion zone did not predict the rise.
- Long runs upstream of a branch point rose about as often as the rising arm.

The literature also cuts both ways. With 11 or more uridines a tract's position between branch
point and AG stops mattering (Coolidge et al. 1997), and a 26-nt tract allows a branch point over
100 nt away (Reed 1989). Branch-point spacing failed to explain the CFTR repeat-length effect when
it was tested directly (Hefferon et al. 2004). The leading alternative is loss of a repressor that
binds the run itself. hnRNP C binds U-tracts of 9 U or more and competes with U2AF2; weakening such
a tract raises inclusion (Zarnack et al. 2013), and a germline poly(T) deletion exonizes an F8 Alu
through hnRNP C (Jourdy et al. 2018). Our pre-registered hnRNP C test asked about dependence on
baseline inclusion, and it failed, but it did not exclude hnRNP C. A repressor that binds the run
would explain why the side of the branch point does not matter. We excluded PTB relief. We regard
the mechanism of the rising arm as open; the minigene experiments below would settle it.

**Limitations.**
- **No tumour DNA at the loci.** Contraction at the R_up loci is inferred from the MSI label, not
  observed. The public MC3 calls hold almost no intronic deletions (458 A/T deletions within 50 nt
  of any exon). We have specified a read-out from the exome reads of 785 loci (D-56), which awaits
  controlled-access data.
- **Exome proxies, not clinical tests.** The DNA truths are exome-based (MC3 deletion counts,
  MANTIS, MSIsensor), not PCR or immunohistochemistry.
- **One resource.** All data are TCGA. An RNA cohort processed elsewhere is the external test the
  code has not had.
- **Predicted branch points.** Branch points are sequence predictions, and most introns use several
  branch points (Pineda and Bradley 2018).
- **Exon-body reads.** Whether the rMATS inclusion levels count exon-body reads is not recorded. If
  they do, unspliced RNA could inflate apparent inclusion. The rise growing with exon length makes
  this worth settling.
- **Not unique to MSI.** The rising-arm exons also move between some sarcoma subtypes, so they are
  a responsive set and not an MSI-only one.
- **Small numbers.** Per-cancer AUROCs rest on 3–11 MMR-deficient tumours.

**What would decide the mechanism.** Minigenes that shorten the tract of rising-arm exons by 2–5
nt, with the branch point intact, then moved, and with a long run placed upstream of the branch
point, would separate geometry from a run-binding repressor. hnRNP C occupancy at rising-arm and
falling-arm tracts would test the alternative directly. The tumours' own alleles at these loci,
read from the exome reads, would test whether each exon follows its own contraction.

**Use.** Score G is not meant to replace DNA tests; trained models reach higher accuracy (AUROC
0.84–0.99 in the training cohorts). Its value is that it is fixed in advance, explains itself
exon by exon, and travels across organs, including those where MSI is rare.

---

## Methods

**Data.** TCGA RNA-seq of 9,729 samples (9,034 primary tumours, one file per case and role) from 33
projects. Splicing was quantified with rMATS (Shen et al. 2014) as percent spliced in (PSI) for
cassette exons (SE), alternative 3′ and 5′ splice sites, mutually exclusive exons and retained
introns. Alternative first and last exons were quantified with HITindex (Fiszbein et al. 2022).
Gene expression is STAR gene counts as log2(CPM + 1). The pan-cancer SE matrix holds 71,208
cassette exons. Per-sample values derive from controlled-access data and are kept locally; only
aggregates are released.

**Subtype labels and models.** Subtype labels are from the TCGA PanCancer Atlas. The MSI classes of
COAD, STAD and UCEC follow the TCGA marker papers (Cancer Genome Atlas Network 2012; Cancer Genome
Atlas Research Network 2013, 2014). For each cohort and view, a multilayer perceptron was tuned with
Optuna and fitted under 5 outer and 3 inner patient-grouped, label-stratified folds. All feature
selection, imputation and scaling were fitted on training folds only. Models were explained with
SHAP DeepExplainer. A model's top exons are the 50 SE events with the largest mean |SHAP| over MSI
patients.

**Acceptor annotation.** For every cassette exon, the upstream intron sequence was read from the
GRCh38 primary assembly in transcript orientation. The acceptor run is the longest T-run
overlapping positions −1 to −50; 99.9% of exons read AG at −2,−1. Repeat exons carry a run of
10 nt or more; background exons 5 nt or less.
- **Branch point.** The yUnAy adenosine nearest the 3′ splice site that lies upstream of the run,
  within 100 nt.
- **AG exclusion zone.** The distance to the first upstream AG at 12 nt or more (Gooding et al.
  2006).
- **Shield exons.** Exons without an acceptor run of 10 nt or more whose intron holds a T-run of
  20 nt or more, starting 25–200 nt upstream, with a yUnAy adenosine between it and the AG.

All annotation was done without labels.

**Scores and exon shifts.**
- **Exon percentiles.** Within each cohort or cancer, each exon's PSI was converted to a percentile
  across tumours.
- **Tumour score.** The mean percentile over a set of exons, minus the mean over background exons
  weighted to the same deciles of mean PSI. This cancels tumour-wide inclusion shifts.
- **Score A** reads all repeat exons. **Score G** is the falling-arm score minus the rising-arm
  score.
- **Exon shift.** The mean percentile in MMR-deficient (or MSI) tumours minus that in the others,
  minus the median shift of background exons in the same decile of mean PSI. In the held-out
  cancers, shifts were pooled over cancers using within-cancer percentiles.

**DNA truths and covariates.**
- **MC3.** Homopolymer deletions are MC3 deletions (filters PASS, wga or native_wga_mix) of A or T
  bases in runs of 5 or more (Ellrott et al. 2018).
- **External MSI calls.** MANTIS and MSIsensor scores were taken from the TCGA PanCancer Atlas via
  cBioPortal, with MSI-high above 0.6 and 10 and stable below 0.4 and 4.
- **Purity** is ABSOLUTE (Carter et al. 2012).
- **Infiltration.** Stromal and immune levels are mean expression of eight marker genes each.

**Statistics and pre-registration.**
- **Pre-registration.** Every test was committed to `config/followup.yaml` with its predicted sign
  and its false-discovery family before it was computed. What had already been seen was listed
  with it.
- **Verdicts.** Supported when the effect has the predicted sign and Benjamini–Hochberg q < 0.05
  within the named family. Against when the sign is opposite with q < 0.05. Otherwise not resolved.
- **Exploratory analyses** were added after results had been read; they are labelled as such and
  recorded as amendments.
- **Tests.** AUROCs use the Mann–Whitney statistic. Adjusted correlations are Spearman correlations
  of rank residuals, with permutation p-values within run-length bins. Regressions use HC3 errors.
- **RV1.** Every contrast was drawn at 20 against 20 tumours, 50 times, so that noise was equal
  across contrasts.

**Software.** Python 3.10 with NumPy, pandas, SciPy, statsmodels, scikit-learn, PyTorch, Optuna and
SHAP; versions are recorded per job.

**Data and code availability.** Code, pre-registrations, aggregate tables and figures are in this
repository. TCGA RNA-seq and the derived per-sample tables require dbGaP authorisation. MC3,
ABSOLUTE and the PanCanAtlas labels are open access from the GDC; MANTIS and MSIsensor scores are
public through cBioPortal.

---

## Supplementary figures

![Supplementary Figure S1](../results/figures/subtype/U16_msi_deletion_load.png)

**Supplementary Figure S1. Mismatch-repair damage tracks the original repeat score in stomach
cancer only.** (a) Per cancer, the original repeat score against the number of exome homopolymer
deletions (log scale). MSI tumours are in colour, the rest grey. (b) Homopolymer deletions of MSI
tumours by cancer, with each cancer's AUROC for the original score.

![Supplementary Figure S2](../results/figures/subtype/U19_msi_code_specificity.png)

**Supplementary Figure S2. The code is specific to the branch point and to mismatch repair.** (a)
Share of repeat exons more included in MSI, by 10-nt bins of branch-point distance, for runs of
10–19 nt and of 20 nt or more; the grey line marks the 36-nt cut. (b) Adjusted correlation of the
MSI shift with branch-point distance, against matched controls: exons without a tract, POLE against
the rest in endometrium, and chromosomally unstable against genomically stable stomach tumours. (c)
AUROC of scores A and G for MSI in the three training cancers.

![Supplementary Figure S3](../results/figures/subtype/U20_msi_code_infiltration.png)

**Supplementary Figure S3. The code is not an infiltration readout.** (a) AUROC of score G in all
held-out tumours and within each half of the immune score, of the immune score alone (grey), and in
the less infiltrated half of each training cancer. (b) Share of the MMR effect on score G kept after
adjustment for stromal, immune and retained-intron levels, with purity or tissue-source site;
bars show 95% HC3 intervals.

![Supplementary Figure S4](../results/figures/subtype/U21_msi_code_msi_calls.png)

**Supplementary Figure S4. The code finds the tumours that two exome MSI callers call MSI-high.**
(a) AUROC of score G against the MC3 call and against MANTIS and MSIsensor MSI-high calls, held out,
then in-sample in the training cancers. (b, c) Score G of held-out tumours by MC3 call and by each
caller's call (mid = between the caller's two cut-offs).

---

## References

- Bonneville R et al. (2017) Landscape of microsatellite instability across 39 cancer types. *JCO Precis Oncol* 2017:PO.17.00073. PMID 29850653.
- Cancer Genome Atlas Network (2012) Comprehensive molecular characterization of human colon and rectal cancer. *Nature* 487:330–337. PMID 22810696.
- Cancer Genome Atlas Research Network (2013) Integrated genomic characterization of endometrial carcinoma. *Nature* 497:67–73. PMID 23636398.
- Cancer Genome Atlas Research Network (2014) Comprehensive molecular characterization of gastric adenocarcinoma. *Nature* 513:202–209. PMID 25079317.
- Carter SL et al. (2012) Absolute quantification of somatic DNA alterations in human cancer. *Nat Biotechnol* 30:413–421. PMID 22544022.
- Chua K, Reed R (2001) An upstream AG determines whether a downstream AG is selected during catalytic step II of splicing. *Mol Cell Biol* 21:1509–1514. PMID 11238888.
- Collura A et al. (2014) Patients with colorectal tumors with microsatellite instability and large deletions in HSP110 T17 have improved response to 5-fluorouracil–based chemotherapy. *Gastroenterology* 146:401–411. PMID 24512910.
- Coolidge CJ, Seely RJ, Patton JG (1997) Functional analysis of the polypyrimidine tract in pre-mRNA splicing. *Nucleic Acids Res* 25:888–896. PMID 9016643.
- Cortes-Ciriano I et al. (2017) A molecular portrait of microsatellite instability across multiple cancers. *Nat Commun* 8:15180. PMID 28585546.
- Corvelo A et al. (2010) Genome-wide association between branch point properties and alternative splicing. *PLoS Comput Biol* 6:e1001016. PMID 21124863.
- Dorard C et al. (2011) Expression of a mutant HSP110 sensitizes colorectal cancer cells to chemotherapy and improves disease prognosis. *Nat Med* 17:1283–1289. PMID 21946539.
- Ellrott K et al. (2018) Scalable open science approach for mutation calling of tumor exomes using multiple genomic pipelines. *Cell Syst* 6:271–281. PMID 29596782.
- Fiszbein A et al. (2022) Widespread occurrence of hybrid internal-terminal exons in human transcriptomes. *Sci Adv* 8:eabk1752. PMID 35044812.
- Giannini G et al. (2002) Human MRE11 is inactivated in mismatch repair-deficient cancers. *EMBO Rep* 3:248–254. PMID 11850399.
- Giannini G et al. (2004) Mutations of an intronic repeat induce impaired MRE11 expression in primary human cancer with microsatellite instability. *Oncogene* 23:2640–2647. PMID 15048091.
- Gooding C et al. (2006) A class of human exons with predicted distant branch points revealed by analysis of AG dinucleotide exclusion zones. *Genome Biol* 7:R1. PMID 16507133.
- Hamanaka K et al. (2023) Genome-wide identification of tandem repeats associated with splicing variation across 49 tissues in humans. *Genome Res* 33:435–447. PMID 37307504.
- Hause RJ et al. (2016) Classification and characterization of microsatellite instability across 18 cancer types. *Nat Med* 22:1342–1350. PMID 27694933.
- Hefferon TW et al. (2004) A variable dinucleotide repeat in the CFTR gene contributes to phenotype diversity by forming RNA secondary structures that alter splicing. *Proc Natl Acad Sci USA* 101:3504–3509. PMID 14993601.
- Jonchère V et al. (2024) Microsatellite instability at U2AF-binding polypyrimidic tract sites perturbs alternative splicing during colorectal cancer initiation. *Genome Biol* 25:210. PMID 39107855.
- Jourdy Y et al. (2018) Reccurrent F8 intronic deletion found in mild hemophilia A causes Alu exonization. *Am J Hum Genet* 102:199–206. PMID 29357978.
- Kim TM, Laird PW, Park PJ (2013) The landscape of microsatellite instability in colorectal and endometrial cancer genomes. *Cell* 155:858–868. PMID 24209623.
- Le DT et al. (2017) Mismatch repair deficiency predicts response of solid tumors to PD-1 blockade. *Science* 357:409–413. PMID 28596308.
- Llorian M et al. (2010) Position-dependent alternative splicing activity revealed by global profiling of alternative splicing events regulated by PTB. *Nat Struct Mol Biol* 17:1114–1123. PMID 20711188.
- Meili D et al. (2009) Disease-causing mutations improving the branch site and polypyrimidine tract: pseudoexon activation of LINE-2 and antisense Alu lacking the poly(T)-tail. *Hum Mutat* 30:823–831. PMID 19280650.
- Niu B et al. (2014) MSIsensor: microsatellite instability detection using paired tumor-normal sequence data. *Bioinformatics* 30:1015–1016. PMID 24371154.
- Pineda JMB, Bradley RK (2018) Most human introns are recognized via multiple and tissue-specific branchpoints. *Genes Dev* 32:577–591. PMID 29666160.
- Reed R (1989) The organization of 3′ splice-site sequences in mammalian introns. *Genes Dev* 3:2113–2123. PMID 2628164.
- Shen S et al. (2014) rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *Proc Natl Acad Sci USA* 111:E5593–E5601. PMID 25480548.
- Taggart AJ et al. (2017) Large-scale analysis of branchpoint usage across species and cell lines. *Genome Res* 27:639–649. PMID 28119336.
- Vilar E et al. (2011) MRE11 deficiency increases sensitivity to poly(ADP-ribose) polymerase inhibition in microsatellite unstable colorectal cancers. *Cancer Res* 71:2632–2642. PMID 21300766.
- Weinstein HNW et al. (2024) RPL22 is a tumor suppressor in MSI-high cancers and a splicing regulator of MDM4. *Cell Rep* 43:114622. PMID 39146182.
- Zarnack K et al. (2013) Direct competition between hnRNP C and U2AF65 protects the transcriptome from the exonization of Alu elements. *Cell* 152:453–466. PMID 23374342.
