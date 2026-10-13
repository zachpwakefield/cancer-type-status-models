**Figure U14. Splice-site mutations skip their exons, and SF3B1 hotspot tumours use cryptic acceptors.**
Pre-specified (config/followup.yaml: brca_mutations, D-59). MC3 public calls of the cohort's tumours (FILTER flags
within PASS, wga, native_wga_mix). (a) Every splice-site SNV at +1/+2 or −1/−2 whose hit exon (by the MC3
transcript's intron number, located in GENCODE v45) is a cassette exon of the main SE matrix: the exon's PSI
percentile in the mutated tumour among the other cohort tumours (ties count half, so an exon near 100% inclusion in
almost every tumour ranks at 0.5), stacked by side; dashed line = the count per bin with no effect; triangles = the
5 CDH1 mutations. (b) SF3B1 hotspot tumours (codons 622, 625, 626, 662, 663, 666, 700, 741, 742) against the rest on
a cryptic A3SS score: mean within-cohort percentile of the long form of 820 A3SS events whose acceptors are 10–30 nt
apart and whose upstream acceptor is rarely used in normals (median PSI < 0.2). Box = interquartile range, diamond =
median, whiskers = 1.5 × IQR, dots = tumours.

- C1, the hit exon is skipped: 110 of 157 mutations below the median; median percentile 0.093; q = 8.0e-7,
  supported. Median percentile: donor 0.05, acceptor 0.13.
- C2, CDH1 alone: 5 of 5 mutations below the median; median percentile 0.001; with five mutations the sign test
  cannot pass (p = 0.062), not resolved.
- C3, cryptic acceptors in SF3B1-hotspot tumours: AUROC 0.90 (14 against 1059), q = 8.0e-7, supported. All 14
  hotspot tumours are LumA or LumB.

Stars: * p < 0.05, ** p < 0.01, *** p < 0.001. q over family C (3 tests). Per-tumour values are kept locally (D-12).
