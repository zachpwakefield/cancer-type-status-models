**Figure U28. Most of the rising-arm exons that move furthest are poison exons.** Pre-specified
(config/followup.yaml: msi_rising_consequence, D-76) before any consequence annotation or host-gene expression was
examined. Each cassette exon was placed in its gene's canonical protein-coding transcript (GENCODE v45): a canonical
exon, an exon that inclusion inserts into the coding sequence (its flanking exons are joined directly inside the
CDS), or other (UTR, non-coding or non-canonical context). Inserted exons were translated in the canonical frame: an
in-frame insertion, a premature stop more than 55 nt before the last exon–exon junction (poison, predicted NMD), or
a stop that escapes NMD.

(a) Consequence of inclusion for the exons of each class measured in a training cancer (COAD, STAD, UCEC);
background exons are shown as they are and reweighted to the rising arm's distribution of proficient PSI (deciles);
"no stop" (< 0.1%) counted with other. (b) Rising-arm exons in the training cancers (exon × cancer): shift (median
PSI of MMR-deficient minus proficient tumours) against the proficient median, by consequence; the five genes with
the largest poison-exon shifts labelled. (c) Training cancers: mean shift by proficient PSI band, for poison (red)
and in-frame (blue) exons, background (filled) and rising arm (hollow); bars = 95% bootstrap intervals; bands with
at least 5 exon × cancer rows. (d) Background exons: per-exon partial Spearman of PSI with the host gene's log2 CPM
across a cancer's tumours, adjusting for the tumour's mean PSI, averaged over the 33 cancers; median and quartiles
by proficient PSI band (A-21). (e) Rising-arm exons inserted into the CDS, training cancers, one per gene and
cancer: shift against the host gene's change (mean log2 CPM of MMR-deficient minus proficient tumours).

- CQ1, background poison exons couple more negatively than in-frame ones (positive control): +0.003, p = 0.31, q =
  0.63, not resolved.
- CQ2, the same for rising-arm exons: −0.025, p = 0.32, q = 0.63, not resolved.
- CQ3, MSI lowers the host gene more for rising poison than in-frame exons (AUROC − 0.5): −0.026, p = 0.71, q =
  0.92, not resolved.
- CQ4, Spearman(shift, gene change), poison minus in-frame: −0.022, p = 0.92, q = 0.92, not resolved.

None of 4 supported (family CQ). The positive control failed, so CQ2–CQ4 cannot be read. Coupling follows baseline
inclusion far more than consequence (panel d; ρ = −0.45 across all annotated exons); within bands, poison background
exons couple slightly more negatively than in-frame ones: −0.007 (p = 0.007), −0.021 (p = 0.002), −0.016 (p = 0.18)
(exploratory, A-21). Panel (a), rising arm: 35% poison, 6% in-frame, 6% escaping NMD, 17% canonical (46, 8 and 8
exons inserted into the CDS); matched background 33%, 9%, 3% and 23%. Panel (b): 12 of the 15 largest shifts plotted
(exon × cancer rows) are rows of 5 poison exons (10 of the 15 largest rows when UTR, non-coding and other exons are
counted); no in-frame exon moves by more than 0.022 (14 of 129 poison and 0 of 23 in-frame rows move by ≥ 0.05;
Fisher p = 0.13; exploratory). Panel (c): background poison exons barely move in MMR-deficient tumours (band means
−0.0086 to −0.0002; held out −0.0065 to +0.0014). Panel (e), host-gene change, BOD1L1: COAD −0.19, STAD −0.19, UCEC
−0.16; CCNDBP1: COAD +0.30, STAD +0.18, UCEC +0.24; BRCA1: COAD +0.04, STAD +0.20. No in-frame rising-arm exon
reaches a median PSI of 0.10 in MMR-deficient tumours while ≤ 0.05 in proficient and normal tissue.
