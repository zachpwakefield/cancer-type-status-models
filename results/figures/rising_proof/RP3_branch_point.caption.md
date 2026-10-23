**Figure RP3. The rise needs the T-run and, in held-out cancers, is rank-associated with the gap to the mapped
branch point; the 36-nt cut is not supported.** Held-out cancers, junction reads. The 52 eligible exons with an
acceptor T-run of ≥ 20 nt and an experimentally mapped branch point upstream of the run (BPHunter, lifted to GRCh38;
31 of them in the rising arm). For these exons the branch point's distance from the AG is exactly the offset of the
run plus its length plus the gap, minus 1, so a cut on that distance mixes run length and position with the gap; the
test therefore asks whether delta is rank-associated with the gap at fixed run length and offset. (a) Delta by gap
bin (descriptive, the bins chosen after the results): dots = exons, solid line = mean, dotted = median. (b) The same
exons: rank of delta against rank of the gap, each with run length and offset (ranks) regressed out; line = least
squares.

- RP3e, partial Spearman of delta with the gap given run length and offset: +0.438, p = 0.004, supported (the
  branch-point test of the verdict); leaving out one exon at a time, 0.40 to 0.47. RP3c, the same with the heuristic
  branch point (139 exons): +0.333, p = 5.0e-4, supported.
- It does not replicate in the training cancers: RPD3b, the same test there, +0.070, p = 0.61, not resolved (62
  exons, 150 MMRd tumours).
- By gap bin (descriptive): 1–5 nt: mean +0.0031, median −0.0001 (22 exons); 6–10 nt: mean +0.0182, median +0.0164
  (17 exons); ≥ 11 nt: mean +0.0081, median +0.0043 (13 exons). The rise sits at gaps above 5 nt, highest at 6–10
  nt, and is not monotone.
- RP3a, needs the T-run (R_up minus exons with a distant branch point and no run): +0.0074, p = 5.0e-4, supported.
- The 36-nt cut on distance: RP3b (R_up minus long runs with a near heuristic branch point) +0.0052, p = 0.16, not
  resolved; RP3d (mapped branch point ≥ 36 nt minus < 36 nt) +0.0012, p = 0.85, not resolved. Pre-specified
  sensitivity cuts for R_up (no q): ≥ 50 nt +0.0068, p = 0.091; > 60 nt +0.0077, p = 0.026. Of the 6 largest rises,
  3 have a mapped branch point < 36 nt from the AG (2 of the 6 largest among the rising-arm exons).
