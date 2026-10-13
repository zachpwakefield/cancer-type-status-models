# MSI splicing marks: what is known, what we add, what to ask next

Written 2026-10-13 to motivate the next round on the acceptor poly-T result
([results/SUBTYPE.md §11](../results/SUBTYPE.md#11-the-msi-splicing-marks-are-written-in-the-genome-acceptor-poly-t-repeats),
D-52 to D-56). Every paper below was checked against PubMed, PMC or a publisher record on this
date; claims are limited to what those records state. Nothing here is a result. Each question in
§3 still needs its predictions written into `config/followup.yaml` before it is run.

## 1. Our finding, in three lines

- The exons that the STAD and COAD MSI models chose sit after long poly-T tracts at the acceptor
  (40 and 34 of the top 50, against 5% of exons). A score from the reference genome alone ranks
  MSI tumours (AUROC 0.71, 0.78), holds within hospital, and is weak in endometrium (0.51).
- The same exons move in three organs (ρ 0.43–0.65 between cancers). The shift grows with tract
  length up to about 20 nt.
- Past about 20 nt it reverses: exons become *more* included. The pre-specified hnRNP C mechanism
  did not explain this.

## 2. What the literature already establishes

### Single genes: intronic poly-T deletions in MSI cause exon skipping
- **MRE11.** An intronic poly(T)11 at the exon 5 acceptor is deleted in MSI tumours; exon 5 is
  skipped and MRE11 is lost (Giannini et al. 2002, 2004). MRE11-mutant MSI colorectal lines are
  more sensitive to the PARP inhibitor ABT-888 (Vilar et al. 2011). *Our positive control
  replicates in all three cancers.*
- **HSP110 (HSPH1).** Deletions in an intronic T17 skip exon 9 and make a dominant-negative
  HSP110ΔE9 (Dorard et al. 2011). In 329 stage II–III MSI colorectal tumours, the deletions were
  mostly biallelic. Patients whose tumours carry *large* deletions (≥ 5 bp) responded better to
  5-fluorouracil (Collura et al. 2014). **Deletion size, not just presence, carries the effect.**

### Genome-wide in colorectal cancer: the closest precedent
- **Jonchère et al. 2024** (Genome Biology 25:210):
  - MSI produces hundreds of non-coding mutations at polypyrimidine U2AF-binding sites, and MSI
    colorectal tumours skip more exons than other colorectal tumours.
  - The mutations arise early, in dMMR colonic crypts before transformation.
  - They explain only part of the exon skipping, and the splicing signature resembles that of
    U2AF1-S34F.
- *Our §11 is an extension of this, not a first.* We agree that the repeats explain part of the
  signal.

### Tract length sets inclusion, quantitatively
- **CFTR intron 8 (5T/7T/9T).** Shorter germline poly-T tracts at the exon 9 acceptor give more
  exon 9 skipping, an inverse relationship measured in 124 people (Chu et al. 1993). This is the
  classic precedent for "length sets the shift".
- **Uridine runs.** In competition assays, tracts of 11 continuous uridines were the strongest.
  Runs of 5–6 worked well only when right next to the AG (Coolidge et al. 1997). So a tract's
  strength saturates.
- **Germline tandem repeats as splicing QTLs.** In GTEx, 9,537 tandem repeats were associated
  with splicing across 49 tissues (Hamanaka et al. 2023), a catalogue of natural length variation.

### MSI genomics: where and how much the repeats change
- **Where.** MSI is commoner in intronic and euchromatic DNA and rarer at nucleosome-occupied
  sequence (Kim et al. 2013, 277 colorectal and endometrial genomes). Expressed genes' introns are
  the mutable ones.
- **How much varies.** Across cancers, MSI loci are partly tumour-type-specific (Hause et al. 2016;
  Cortes-Ciriano et al. 2017). Hause et al. propose MSI as a continuous phenotype whose burden
  tracks survival.
- **The indel signature.** ID2 (slippage at mononucleotide repeats) is strongly raised in
  mismatch-repair-deficient tumours (Alexandrov et al. 2020; COSMIC).
- **Why endometrium may differ.** At the same five mononucleotide markers, MSI endometrioid
  cancers shift by about 3 nt and colorectal cancers by about 6 (Wang et al. 2017; P < 0.0001).
  Long mononucleotide markers detect endometrial MSI better (Chung et al. 2023). *This is a
  ready-made explanation for our weak UCEC result: the same exons, but shorter deletions.*

### Splice-site geometry: a candidate mechanism for the reversal
- **The branch point to AG window.** The AG used is selected within a narrow window, optimally
  19–23 nt downstream of the branch point (Chua and Reed 2001). Human branch points sit close to
  the 3′ splice site (59,359 mapped by Mercer et al. 2015).
- **Why that matters here.** A poly-T tract lies between the branch point and the AG. A very long
  tract pushes the branch point beyond the optimal window, so shortening it could bring the AG
  back into the window and *raise* exon use (our hypothesis, not a published result).
- **The earlier candidate.** hnRNP C competes with U2AF65 on long U-tracts (Zarnack et al. 2013).
  Our D-53 test of that route failed.

### Consequences that make this worth pursuing
- **Synthetic lethality from repeat instability.** Large TA-repeat expansions in MSI cells stall
  forks and make them depend on WRN (van Wietmarschen et al. 2020). MRE11 skipping predicts
  PARP-inhibitor sensitivity (Vilar et al. 2011). *Repeat-driven splicing losses are a route to
  MSI-specific dependencies.*
- **Immunotherapy.** Mismatch-repair deficiency predicts PD-1 response across 12 tumour types (53%
  objective response; Le et al. 2017). Splicing errors and microsatellite transcription errors
  also generate frameshift neoantigens (Shen et al. 2019, PMC6775166). *Exons that skip the same
  way in three organs could give shared, predictable splice neoantigens.*
- **Detection from RNA.** Expression-based MSI callers exist: MSIsensor-RNA (Jia et al. 2024) and
  MSI expresso. Ours differs: a mechanism-based score from the reference genome with no training.
- **An experimental model.** MLH1-knockout human colon organoids accumulate the
  mismatch-repair-deficient mutation profile over passages (Drost et al. 2017).

## 3. Questions this motivates, in suggested order

**Q1. Dose: does each tumour's poly-T deletion load set its splicing score?** *(no download; MC3
is local)*
- *Motivation:* Collura et al. 2014 (deletion size matters); Wang et al. 2017 (shorter shifts in
  endometrium); Hause et al. 2016 (MSI is continuous); ID2.
- *Test:* from MC3, count each tumour's deletions at T/A homopolymers and their lengths.
  - Within MSI tumours, the repeat-exon score should rise with that load.
  - MSI UCEC should carry shorter deletions (more 1 bp, fewer ≥ 2 bp) than STAD and COAD.
  - The splicing shift per unit of load should be similar across the three organs.
- *Caveat:* MC3 is exome-wide and mostly coding, so it is a proxy for the intronic tracts.
- *Status (D-61, done):* dose tracks the score in STAD only (ρ −0.45, q = 0.005). UCEC deletions
  are not shorter, and UCEC MSI tumours carry the most of them, so deletion size does not explain
  the endometrial gap ([SUBTYPE.md §11](../results/SUBTYPE.md), U16).

**Q2. The reversal: is it branch-point geometry?** *(small download: published branch-point
annotations)*
- *Motivation:* Chua and Reed 2001; Coolidge et al. 1997; Mercer et al. 2015.
- *Test:* for exons with tracts of ≥ 20 nt, compare the branch point to AG distance. Exons pushed
  past the optimal window should be the ones that gain inclusion in MSI.
  - Also check whether the long tracts are Alu-derived (antisense Alu poly-A), the setting of
    Zarnack et al. 2013.
- *Would refute:* reversal exons with branch points already in the window.
- *Status (D-63, done):* the annotations could not be fetched, so branch points were read from
  sequence. Supported in all three cancers: a farther branch point means more inclusion in MSI,
  beyond run length (U17). The Alu k-mer check found almost no exons (≤ 3 per class).

**Q3. A natural experiment in the germline** *(small download: the Hamanaka et al. 2023 catalogue)*
- *Motivation:* CFTR 5T/7T/9T (Chu et al. 1993); germline splicing tandem repeats.
- *Test:* are our MSI repeat exons germline splicing QTLs in GTEx? If so, does shorter germline
  length move PSI in the same direction as MSI contraction, including the ≥ 20 nt reversal?
  - This would test causality without any tumour.
- *Status (D-64, not testable):* the Hamanaka catalogue genotypes no homopolymer, and GTEx v8 sQTL
  pairs hold no variant inside any of the runs. It needs a call set that genotypes homopolymer
  length.

**Q4. The deletions themselves** *(pending: the 4-BAM pilot you run with your token; D-56 is
already pre-specified)*
- Contraction alleles at the 785 loci, and whether each tumour's exon shift follows its own
  allele.
- Adds the deletion-size question from Collura et al. at the locus level.

**Q5. Consequences: frame, dependencies, neoantigens** *(no download)*
- *Test:* which repeat exons shift the frame or trigger NMD when skipped, and which genes are hit.
  MRE11 is the DNA-repair precedent, and TBK1, an innate-immunity kinase, is among our largest
  shared shifts.
  - Frameshifting exons that skip in all three organs are candidate shared splice neoantigens.
- *Motivation:* Vilar et al. 2011; van Wietmarschen et al. 2020; Le et al. 2017; Shen et al. 2019.
- *Status (D-62, not resolved):* no NMD signature in gene expression. Of the 12 exons skipped in
  all three cancers, two are frameshifting coding exons: MRE11 exon 5 and an exon of GFM1.

**Q6. Timing and causality in a model** *(external)*
- *Motivation:* Jonchère et al. 2024 (early, in crypts); Drost et al. 2017 (MLH1-knockout
  organoids).
- *Test:* the genome-only score in MLH1-knockout organoid RNA-seq over passages should rise as
  deletions accumulate.
  - Minigenes that shorten the CCNDBP1 or ANXA7 tract would test the reversal directly.

**Q7. Use: an RNA test for MSI that explains itself** *(external replication)*
- Compare the genome-only score with MSIsensor-RNA in TCGA, then in a cohort processed elsewhere
  (CPTAC colon or endometrial).
- Its value would be transparency and organ-independence, not raw accuracy.
- *Status (D-66, done within TCGA):* the two-arm sequence code finds MMR-deficient tumours, called
  from MC3 homopolymer deletions, in 30 cancers never used to define it (AUROC 0.83; 0.85 against
  hypermutated proficient tumours). An external cohort is still to come.
- *Status (D-68, done):* the code is not an infiltration readout. Adjusting for stromal and
  immune expression and the retained-intron level keeps the whole MMR effect (1.01 held out,
  1.01–1.07 in the training cancers), and the lower-immune half separates at 0.84.
- *Status (D-69, done):* against two exome callers that do not use MC3, score G separates MSI-H
  from MSS in the 30 held-out cancers at AUROC 0.97 (MANTIS, 17 MSI-H) and 0.89 (MSIsensor, 23).
  It does not follow the continuous MANTIS score below MSI-H (ρ +0.01). An RNA cohort processed
  elsewhere is still to come.

**Start with Q1 and Q5:** both use only local data, and Q1 directly tests the endometrium
explanation. **Q2 and Q3 each need one small download**, which would be stated and approved first.

## References

- Alexandrov LB et al. (2020) The repertoire of mutational signatures in human cancer. *Nature* 578:94. — [COSMIC ID2](https://cancer.sanger.ac.uk/signatures/id/id2/)
- Chu CS et al. (1993) Genetic basis of variable exon 9 skipping in cystic fibrosis transmembrane conductance regulator mRNA. *Nat Genet* 3:151.
- Chua K, Reed R (2001) An upstream AG determines whether a downstream AG is selected during catalytic step II of splicing. *Mol Cell Biol* 21:1509. — [PMC86697](https://pmc.ncbi.nlm.nih.gov/articles/PMC86697)
- Chung Y et al. (2023) Evaluation of an eight marker-panel including long mononucleotide repeat markers to detect microsatellite instability in colorectal, gastric, and endometrial cancers. *BMC Cancer* 23:1100.
- Collura A et al. (2014) Patients with colorectal tumors with microsatellite instability and large deletions in HSP110 T17 have improved response to 5-fluorouracil–based chemotherapy. *Gastroenterology* 146(2).
- Coolidge CJ, Seely RJ, Patton JG (1997) Functional analysis of the polypyrimidine tract in pre-mRNA splicing. *Nucleic Acids Res* 25:888. — [PMC146492](https://pmc.ncbi.nlm.nih.gov/articles/PMC146492)
- Cortes-Ciriano I et al. (2017) A molecular portrait of microsatellite instability across multiple cancers. *Nat Commun* 8:15180. — [PMC5467167](https://pmc.ncbi.nlm.nih.gov/articles/PMC5467167)
- Dorard C et al. (2011) Expression of a mutant HSP110 sensitizes colorectal cancer cells to chemotherapy and improves disease prognosis. *Nat Med* 17:1283.
- Drost J et al. (2017) Use of CRISPR-modified human stem cell organoids to study the origin of mutational signatures in cancer. *Science* 358:234. — [PMC6038908](https://pmc.ncbi.nlm.nih.gov/articles/PMC6038908)
- Giannini G et al. (2002) Human MRE11 is inactivated in mismatch repair-deficient cancers. *EMBO Rep* 3:248.
- Giannini G et al. (2004) Mutations of an intronic repeat induce impaired MRE11 expression in primary human cancer with microsatellite instability. *Oncogene* 23:2640.
- Hamanaka K et al. (2023) Genome-wide identification of tandem repeats associated with splicing variation across 49 tissues in humans. *Genome Res* 33:435. — [PMC10078293](https://pmc.ncbi.nlm.nih.gov/articles/PMC10078293)
- Hause RJ et al. (2016) Classification and characterization of microsatellite instability across 18 cancer types. *Nat Med* 22:1342. — [PMID 27694933](https://pubmed.ncbi.nlm.nih.gov/27694933/)
- Jia P et al. (2024) MSIsensor-RNA: microsatellite instability detection for bulk and single-cell gene expression data. *Genomics Proteomics Bioinformatics*. — [PMC12016039](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12016039/)
- Jonchère V et al. (2024) Microsatellite instability at U2AF-binding polypyrimidic tract sites perturbs alternative splicing during colorectal cancer initiation. *Genome Biol* 25:210. — [PubMed 39107855](https://pubmed.ncbi.nlm.nih.gov/39107855/)
- Kim TM, Laird PW, Park PJ (2013) The landscape of microsatellite instability in colorectal and endometrial cancer genomes. *Cell* 155:858. — [PMC3871995](https://pmc.ncbi.nlm.nih.gov/articles/PMC3871995)
- Le DT et al. (2017) Mismatch repair deficiency predicts response of solid tumors to PD-1 blockade. *Science* 357:409. — [PubMed 28596308](https://pubmed.ncbi.nlm.nih.gov/28596308/)
- Mercer TR et al. (2015) Genome-wide discovery of human splicing branchpoints. *Genome Res* 25:290. — [PMC4315302](https://pmc.ncbi.nlm.nih.gov/articles/4315302)
- MSI expresso: a software for determining MSI status and detecting MSI-related transcription events from RNA sequencing data. — [PMC12058494](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12058494/)
- Shen L et al. (2019) RNA transcription and splicing errors as a source of cancer frameshift neoantigens for vaccines. *Sci Rep* 9, doi:10.1038/s41598-019-50738-4. — [PMC6775166](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6775166/)
- van Wietmarschen N et al. (2020) Repeat expansions confer WRN dependence in microsatellite-unstable cancers. *Nature* 586:292.
- Vilar E et al. (2011) MRE11 deficiency increases sensitivity to poly(ADP-ribose) polymerase inhibition in microsatellite unstable colorectal cancers. *Cancer Res* 71:2632. — [PMC3407272](https://pmc.ncbi.nlm.nih.gov/articles/PMC3407272)
- Wang Y et al. (2017) Differences in microsatellite instability profiles between endometrioid and colorectal cancers: a potential cause for false-negative results? *J Mol Diagn* 19:57. — [PMC5225298](https://pmc.ncbi.nlm.nih.gov/articles/PMC5225298)
- Zarnack K et al. (2013) Direct competition between hnRNP C and U2AF65 protects the transcriptome from the exonization of Alu elements. *Cell* 152:453.
