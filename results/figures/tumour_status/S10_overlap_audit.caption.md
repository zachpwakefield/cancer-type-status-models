**Figure S10. Some of the models' events, terminal exons above all, read a neighbouring gene's expression.** TCGA
mRNA libraries are mostly unstranded, and GENCODE readthrough and nested models put exons of one gene inside
another, so an event named after one gene can count reads of a neighbour (D-27). (a) The share of events whose
segment (the alternative exon, both exons of an MXE pair, the retained intron or the splice-site extension) overlaps
an exon of another GENCODE v45 gene by >= 50 bp, on either strand, for every measured event (grey; 6.8% over all 8
views) and for the models' top events (top 100 of each cancer's model, per view; dark). The gap is widest for HIT:
44% against 7%. Odds ratios with the queried length held fixed (Mantel-Haenszel over length deciles, because long
terminal exons overlap something more often by chance) are 1.2 to 6.6 across views; every statistic is in
`S10_overlap_audit.csv`.

(b) Of the top events with such an overlap, how many *track* the other gene — Spearman with its log2 CPM of at least
0.5 in size and larger than with the host gene's own (within each cancer, median over the cancers whose top 100
holds the event); the criterion was fixed before any result was seen. 163 of 1,105 do. Those are better read as a
readout of the neighbour's expression than as a splicing choice of the host, whether the coupling is misassigned
reads or genuine co-transcription; the two cannot be told apart without stranded libraries. Among the events the
narrative names: PLCG1 MXE follows ZHX3 (ρ −0.60, own gene +0.20); EXOSC7 ALE follows CLEC3B (ρ +0.76, own gene
−0.28); ENOSF1 HIT follows TYMS (ρ −0.62, own gene +0.34).

No model was refitted: the audit changes how events are described, not what the models are.
