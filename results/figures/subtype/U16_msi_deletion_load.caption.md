**Does more mismatch-repair damage mean more skipping of repeat exons?**

Pre-specified (config/followup.yaml: msi_deletion_load, D-61). Homopolymer deletion = an MC3 deletion (FILTER PASS,
wga or native_wga_mix) of A or T bases inside a run of ≥ 5 of that base, counted per tumour in the exome. Score =
the D-52 per-tumour repeat-exon score: mean inclusion percentile of exons whose acceptor carries a T-run ≥ 10, minus
that of matched exons without one (lower = repeat exons less included). Partial Spearman of score and log(1 +
count), adjusted for ABSOLUTE purity, within MSI tumours.

a, one panel per cancer: score against homopolymer deletions (x on a log(1 + count) scale); MSI tumours in colour,
the rest grey. The rest carry a median of 0 homopolymer deletions in every cancer (upper quartile ≤ 1); MSI tumours
a median of 90 (COAD) to 147 (UCEC). b, MSI tumours' homopolymer deletions by cancer (box = interquartile range,
diamond = median, whiskers 1.5 × IQR, points = tumours), with each cancer's D-52 score AUROC for MSI against the
rest (second line of each label).

- DL1–DL3, more deletions, lower score: STAD ρ = −0.45, q = 0.005, supported; COAD ρ = +0.13, q = 0.58, not
  resolved; UCEC ρ = −0.13, q = 0.58, not resolved.
- DL4–DL5, a smaller share of multi-base deletions in UCEC (explaining its weaker shift): not resolved; the median
  share is 0.010 in UCEC, 0.007 in STAD and 0.000 in COAD. Over 98% of these deletions remove one base.
- DL6, multi-base deletions beyond single-base ones: ρ = −0.03, q = 0.78, not resolved.

Dose tracks the score within STAD only. Between cancers it runs the other way: UCEC MSI tumours carry the most
deletions and the D-52 score does not separate them from the rest (AUROC 0.51). The exome counts coding runs, mostly
shorter than 10; the intronic tracts themselves were not measured. q over family DL (6 tests); per-tumour values are
kept locally (D-12).
