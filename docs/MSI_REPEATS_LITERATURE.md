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
  Section 5 compares the two directly.
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
- *Status (D-70, done):* the rise is real, held out and specific to mismatch repair against
  mutation burden, and it needs the tract. Its cause is open: it follows the predicted branch
  point (within 30 nt of the AG long-run exons fall, beyond 30 they rise), but not the AG exclusion
  zone, not the run's side of the branch point, and not deletion dose; PTB relief runs the wrong
  way (U22; §4 below).

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

## 4. Prior work on the rising arm (literature search, 2026-10-13)

A search of PubMed and Europe PMC full text asked whether anyone has reported that MSI makes
exons after long poly-T tracts *more* included. Every reference below was opened and checked
against PubMed (PMIDs given).

**No report of the effect.**
- **Jonchère et al. 2024** (PMID 39107855) is the only genome-wide MSI splicing study at these
  sites. Most changes it reports are exon exclusion, and larger intronic deletions go with more
  exclusion. Inclusion gains appear only in a supplementary panel and are not analysed. It does not
  stratify by tract length or branch-point distance.
- **Reviews and tools** treat tract contraction as a cause of skipping only (Corcos et al. 2023,
  PMID 36833239).
- **One MSI inclusion effect is known, but it acts in trans.** RPL22 frameshifts raise MDM4 exon 6
  inclusion (Weinstein et al. 2024, PMID 39146182).

**Evidence for the geometry reading.**
- **The optimum.** Step II of splicing works best with the branch point 19–23 nt from the AG and
  falls off progressively beyond (Chua and Reed 2001, PMID 11238888).
- **Distance and skipping.** Genome-wide, exons with distant predicted branch points are skipped
  more (Corvelo et al. 2010, PMID 21124863). Distance weakens AG recognition in mapped branch points
  (Taggart et al. 2017, PMID 28119336). 90% of mapped branch points lie 19–37 nt upstream (Mercer et
  al. 2015, PMID 25561518), so our ≥ 36 nt cut sits at the far edge of normal.
- **A repeat that moves a branch point.** In NARS2, a germline (GTTTTTT)n expansion pushes the
  branch point from 35 to 49 nt from a pseudoexon's 3′ss and stops its splicing in minigenes
  (Hamanaka et al. 2023, PMID 37307504). That is our model run in the other direction.
- **Contraction by deletion.** In PTS, a deletion that removes an antisense Alu's poly(T) tail
  brings a favourable branch point close to a pseudoexon's 3′ss and activates it (Meili et al.
  2009, PMID 19280650). The deletion also removes a repressive AG, so the case is not clean.

**Evidence against it, or for something else.**
- **Long tracts do not care where they sit.** With 11 or more continuous uridines, the tract's
  position between branch point and AG does not matter (Coolidge et al. 1997, PMID 9016643). A
  26-nt tract allows lariat formation with the branch point over 100 nt away (Reed 1989,
  PMID 2628164). Both argue that runs ≥ 20 should be the *least* distance-sensitive.
- **Branch-point spacing failed where it was tested.** At CFTR's (TG)m(T)n locus, mutating the
  candidate branch points did not explain the repeat-length effect (Hefferon et al. 2004,
  PMID 14993601).
- **hnRNP C.** hnRNP C binds U-tracts of 9 U or more and competes with U2AF65; weakening a long
  U-tract raises inclusion (Zarnack et al. 2013, PMID 23374342). A recurrent F8 deletion of part of
  an Alu poly(T) tail exonizes the Alu, and hnRNP C knockdown reproduces it (Jourdy et al. 2018,
  PMID 29357978). Our D-53 test of hnRNP C (baseline dependence) was not supported, but that test
  did not rule hnRNP C out.
- **PTB.** Long tracts with long AG exclusion zones mark PTB-repressed exons (Gooding et al. 2006,
  PMID 16507133; Llorian et al. 2010, PMID 20711188), and U2AF35-dependent ones (Kralovicova et al.
  2015, PMID 25779042). D-70's RV8 found the opposite of the PTB-relief prediction.
- **Many branch points per intron.** Most introns use several branch points, often by tissue
  (Pineda and Bradley 2018, PMID 29666160). One sequence-predicted branch point per exon is a weak
  proxy.

**What would decide it.**
- **Minigenes** that shorten the tract of an R_up exon by 2–5 nt:
  - with the branch point intact, then moved;
  - with a long run placed upstream of the branch point, as RV5 asks.
- **The tumours' own alleles** at the R_up loci (D-56's slices).
- **hnRNP C occupancy.** hnRNP C eCLIP at the R_up tracts against R_down tracts.

## 5. Where the finding stands (literature search and D-74, 2026-10-13)

A second search (PubMed, Europe PMC including bioRxiv and medRxiv, the web; 2015–2026) asked what
has been done on RNA detection of MSI. Every PMID below was checked against PubMed.

**What others have done.**
- **Expression classifiers.** All were trained on TCGA COAD, STAD and UCEC, or on colon alone:
  - Danaher et al. 2019 (PMID 30665466): 14 genes, including EPM2AIP1 and the MMR genes.
  - Pačínková and Popovici 2019 (PMID 31008109): 25 genes; colon 0.94, STAD 0.90, UCEC 0.71
    (abstract).
  - PreMSIm (Li et al. 2020, PMID 32257050): 15 genes, k-NN.
  - MSIsensor-RNA (Jia et al. 2024, PMID 39341794): 109 genes, SVM, external AUC 0.95–0.99 within
    the three cancers. Its authors note weaker results when the test cancer differs from training
    (full text as read by the search; not re-checked).
  - Sorokin et al. 2021 (PMID 34888350) re-tested published signatures on two more cancers.
- **Repeat lengths read in RNA.** MIRACLE (Choi et al. 2024, PMID 39210504) uses 3′ UTR repeats;
  its pan-cancer calls are distributions, not held-out AUROCs. MSI expresso (Tubacher et al. 2025,
  PMID 40342963) uses 192 repeats, with about 99% concordance in mostly colon samples.
- **An independent benchmark** found RNA tools outperformed by DNA tools run on the same RNA data
  (Anthony et al. 2024, PMID 39129364, abstract). By the full text, as read by the search and not
  re-checked: in 143 colon tumours, PreMSIm and MSIsensor-RNA had precision 0.22 and 0.15 at their
  defaults, against F1 about 0.82 for DNA tools on the RNA alignments.
- **Splicing.** Jonchère et al. 2024 (PMID 39107855) is colon-only. It has no classifier, and it
  reports mostly exclusion with no length- or branch-point-dependent sign change. No
  splicing-based MSI classifier was found.
- **Endometrium.** MSI endometrial cancers shift their repeats less than colon cancers (Wang et
  al. 2017, PMID 27810331; Wu et al. 2019, PMID 30443012; Chung et al. 2023, PMID 37953261).
  That fits the model's 19 endometrial misses (D-73).
- **Outside the three cancers.** MMR deficiency is more often MSH2/MSH6 or Lynch-associated there
  (Latham et al. 2019, PMID 30376427). Guidelines find too little evidence to recommend a
  specific assay outside the gastrointestinal tract and endometrium (Bartley et al. 2022, PMID
  35920830).

**What D-74 added (Figure U26).** Expression models trained exactly as the splicing model
transferred almost as well to the 30 other cancers:
- held-out AUROC 0.893 on log CPM and 0.866 on percentiles, against 0.900 for splicing;
- fixed panels did not (PreMSIm genes 0.70–0.79);
- MLH1 was silenced in 28 of the 50 held-out deficient tumours, and where it was not, the two tied;
- the two were complementary, together reaching 0.926 and precision 0.68 at half recall;
- splicing agreed more closely with the exome callers within a cancer (0.998 against 0.975–0.979).

**Verdict.**
- **Novel.** Three things are new:
  - the rising arm (long tracts with a distant branch point gain inclusion);
  - a sequence rule, read from the genome, for which exons move;
  - a held-out test across 30 cancers against DNA calls.
- **Not novel.** Skipping at MSI tracts, RNA detection of MSI, and weak endometrial signal.
- **Weaker than first framed.** "Transfers across organs" is not special to splicing; a retrained
  expression model does nearly as well. What splicing adds is complementarity (better together),
  closer agreement with exome callers within a cancer, and a mechanism.
- **The comparison still missing** is a read-length caller on the RNA alignments (MSIsensor-pro,
  MIRACLE). It needs the RNA BAMs, which are controlled access.

## References

- Alexandrov LB et al. (2020) The repertoire of mutational signatures in human cancer. *Nature* 578:94. — [COSMIC ID2](https://cancer.sanger.ac.uk/signatures/id/id2/)
- Chu CS et al. (1993) Genetic basis of variable exon 9 skipping in cystic fibrosis transmembrane conductance regulator mRNA. *Nat Genet* 3:151.
- Chua K, Reed R (2001) An upstream AG determines whether a downstream AG is selected during catalytic step II of splicing. *Mol Cell Biol* 21:1509. — [PMC86697](https://pmc.ncbi.nlm.nih.gov/articles/PMC86697)
- Chung Y et al. (2023) Evaluation of an eight marker-panel including long mononucleotide repeat markers to detect microsatellite instability in colorectal, gastric, and endometrial cancers. *BMC Cancer* 23:1100.
- Collura A et al. (2014) Patients with colorectal tumors with microsatellite instability and large deletions in HSP110 T17 have improved response to 5-fluorouracil–based chemotherapy. *Gastroenterology* 146(2).
- Coolidge CJ, Seely RJ, Patton JG (1997) Functional analysis of the polypyrimidine tract in pre-mRNA splicing. *Nucleic Acids Res* 25:888. — [PMC146492](https://pmc.ncbi.nlm.nih.gov/articles/PMC146492)
- Corcos L et al. (2023) Microsatellite instability and aberrant pre-mRNA splicing: how intimate is it? *Genes* 14:311. — [PMID 36833239](https://pubmed.ncbi.nlm.nih.gov/36833239/)
- Cortes-Ciriano I et al. (2017) A molecular portrait of microsatellite instability across multiple cancers. *Nat Commun* 8:15180. — [PMC5467167](https://pmc.ncbi.nlm.nih.gov/articles/PMC5467167)
- Corvelo A et al. (2010) Genome-wide association between branch point properties and alternative splicing. *PLoS Comput Biol* 6:e1001016. — [PMID 21124863](https://pubmed.ncbi.nlm.nih.gov/21124863/)
- Dorard C et al. (2011) Expression of a mutant HSP110 sensitizes colorectal cancer cells to chemotherapy and improves disease prognosis. *Nat Med* 17:1283.
- Drost J et al. (2017) Use of CRISPR-modified human stem cell organoids to study the origin of mutational signatures in cancer. *Science* 358:234. — [PMC6038908](https://pmc.ncbi.nlm.nih.gov/articles/PMC6038908)
- Giannini G et al. (2002) Human MRE11 is inactivated in mismatch repair-deficient cancers. *EMBO Rep* 3:248.
- Giannini G et al. (2004) Mutations of an intronic repeat induce impaired MRE11 expression in primary human cancer with microsatellite instability. *Oncogene* 23:2640.
- Gooding C et al. (2006) A class of human exons with predicted distant branch points revealed by analysis of AG dinucleotide exclusion zones. *Genome Biol* 7:R1. — [PMID 16507133](https://pubmed.ncbi.nlm.nih.gov/16507133/)
- Hamanaka K et al. (2023) Genome-wide identification of tandem repeats associated with splicing variation across 49 tissues in humans. *Genome Res* 33:435. — [PMC10078293](https://pmc.ncbi.nlm.nih.gov/articles/PMC10078293)
- Hause RJ et al. (2016) Classification and characterization of microsatellite instability across 18 cancer types. *Nat Med* 22:1342. — [PMID 27694933](https://pubmed.ncbi.nlm.nih.gov/27694933/)
- Hefferon TW et al. (2004) A variable dinucleotide repeat in the CFTR gene contributes to phenotype diversity by forming RNA secondary structures that alter splicing. *Proc Natl Acad Sci USA* 101:3504. — [PMID 14993601](https://pubmed.ncbi.nlm.nih.gov/14993601/)
- Jia P et al. (2024) MSIsensor-RNA: microsatellite instability detection for bulk and single-cell gene expression data. *Genomics Proteomics Bioinformatics*. — [PMC12016039](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12016039/)
- Jonchère V et al. (2024) Microsatellite instability at U2AF-binding polypyrimidic tract sites perturbs alternative splicing during colorectal cancer initiation. *Genome Biol* 25:210. — [PubMed 39107855](https://pubmed.ncbi.nlm.nih.gov/39107855/)
- Jourdy Y et al. (2018) Reccurrent F8 intronic deletion found in mild hemophilia A causes Alu exonization. *Am J Hum Genet* 102:199. — [PMID 29357978](https://pubmed.ncbi.nlm.nih.gov/29357978/)
- Kim TM, Laird PW, Park PJ (2013) The landscape of microsatellite instability in colorectal and endometrial cancer genomes. *Cell* 155:858. — [PMC3871995](https://pmc.ncbi.nlm.nih.gov/articles/PMC3871995)
- Kralovicova J et al. (2015) Identification of U2AF(35)-dependent exons by RNA-Seq reveals a link between 3′ splice-site organization and activity of U2AF-related proteins. *Nucleic Acids Res* 43:3747. — [PMID 25779042](https://pubmed.ncbi.nlm.nih.gov/25779042/)
- Le DT et al. (2017) Mismatch repair deficiency predicts response of solid tumors to PD-1 blockade. *Science* 357:409. — [PubMed 28596308](https://pubmed.ncbi.nlm.nih.gov/28596308/)
- Llorian M et al. (2010) Position-dependent alternative splicing activity revealed by global profiling of alternative splicing events regulated by PTB. *Nat Struct Mol Biol* 17:1114. — [PMID 20711188](https://pubmed.ncbi.nlm.nih.gov/20711188/)
- Meili D et al. (2009) Disease-causing mutations improving the branch site and polypyrimidine tract: pseudoexon activation of LINE-2 and antisense Alu lacking the poly(T)-tail. *Hum Mutat* 30:823. — [PMID 19280650](https://pubmed.ncbi.nlm.nih.gov/19280650/)
- Mercer TR et al. (2015) Genome-wide discovery of human splicing branchpoints. *Genome Res* 25:290. — [PMC4315302](https://pmc.ncbi.nlm.nih.gov/articles/4315302)
- MSI expresso: a software for determining MSI status and detecting MSI-related transcription events from RNA sequencing data. — [PMC12058494](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12058494/)
- Pineda JMB, Bradley RK (2018) Most human introns are recognized via multiple and tissue-specific branchpoints. *Genes Dev* 32:577. — [PMID 29666160](https://pubmed.ncbi.nlm.nih.gov/29666160/)
- Reed R (1989) The organization of 3′ splice-site sequences in mammalian introns. *Genes Dev* 3:2113. — [PMID 2628164](https://pubmed.ncbi.nlm.nih.gov/2628164/)
- Shen L et al. (2019) RNA transcription and splicing errors as a source of cancer frameshift neoantigens for vaccines. *Sci Rep* 9, doi:10.1038/s41598-019-50738-4. — [PMC6775166](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6775166/)
- Taggart AJ et al. (2017) Large-scale analysis of branchpoint usage across species and cell lines. *Genome Res* 27:639. — [PMID 28119336](https://pubmed.ncbi.nlm.nih.gov/28119336/)
- van Wietmarschen N et al. (2020) Repeat expansions confer WRN dependence in microsatellite-unstable cancers. *Nature* 586:292.
- Vilar E et al. (2011) MRE11 deficiency increases sensitivity to poly(ADP-ribose) polymerase inhibition in microsatellite unstable colorectal cancers. *Cancer Res* 71:2632. — [PMC3407272](https://pmc.ncbi.nlm.nih.gov/articles/PMC3407272)
- Wang Y et al. (2017) Differences in microsatellite instability profiles between endometrioid and colorectal cancers: a potential cause for false-negative results? *J Mol Diagn* 19:57. — [PMC5225298](https://pmc.ncbi.nlm.nih.gov/articles/PMC5225298)
- Weinstein HNW et al. (2024) RPL22 is a tumor suppressor in MSI-high cancers and a splicing regulator of MDM4. *Cell Rep* 43:114622. — [PMID 39146182](https://pubmed.ncbi.nlm.nih.gov/39146182/)
- Zarnack K et al. (2013) Direct competition between hnRNP C and U2AF65 protects the transcriptome from the exonization of Alu elements. *Cell* 152:453.
