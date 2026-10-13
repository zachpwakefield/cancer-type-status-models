**Figure T9. Many of the cancer-type models' terminal exons read a neighbouring lineage gene.** TCGA mRNA libraries
are mostly unstranded, and GENCODE readthrough and nested models put exons of one gene inside another, so an event
named after one gene can count reads of a neighbour (D-27). (a) The share of events whose segment (the alternative
exon, both exons of an MXE pair, the retained intron or the splice-site extension) overlaps an exon of another
GENCODE v45 gene by >= 50 bp, on either strand, for every measured event (grey; 6.8% over all 8 views) and for the
models' top events (top 30 of each class's explanation, per view; dark). The gap is widest for HIT: 61% against 7%.
Odds ratios with the queried length held fixed (Mantel-Haenszel over length deciles, because long terminal exons
overlap something more often by chance) are 1.8 to 12.6 across views; every statistic is in `T9_overlap_audit.csv`.

(b) Of the top events with such an overlap, how many *track* the other gene — Spearman with its log2 CPM of at least
0.5 in size and larger than with the host gene's own (across all 9,729 samples); the criterion was fixed before any
result was seen. 101 of 695 do. Those are better read as a readout of the neighbour's expression than as a splicing
choice of the host, whether the coupling is misassigned reads or genuine co-transcription; the two cannot be told
apart without stranded libraries. Among the events the narrative names: PLEKHG3 HIT follows SPTB (ρ −0.61, own gene
+0.16); BSCL2 HIT follows LRRN4CL (ρ −0.78, own gene −0.03).

No model was refitted: the audit changes how events are described, not what the models are.
