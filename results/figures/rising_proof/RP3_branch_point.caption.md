**Figure RP3. The rise needs the T-run, a step at 36 nt is not supported, and the link to the branch-point gap is
seen in held-out cancers only.** MMRd = MMR-deficient (≥ 10 MC3 homopolymer deletions); MMRp = MMR-proficient;
change in inclusion = an exon's mean adjusted PSI in the MMRd tumours minus the MMRp mean of the same cancer, from
junction reads only (rMATS JC); a tumour's adjusted PSI is its deviation from the cancer's median minus its mean
deviation over background exons of the same baseline band. (a) Held-out cancers: mean change with its 95% interval
for the rising arm (rising arm = acceptor T-run ≥ 20 nt with the sequence-predicted branch point ≥ 36 nt from the
AG) and for background exons whose predicted branch point is also ≥ 36 nt from the AG, without a T-run; dashed = the
smallest effect of interest. (b) Held-out cancers, exons with an acceptor T-run ≥ 20 nt and a predicted branch point
(rising arm right of the dashed line, the rule's 36-nt cut; long runs with a nearer branch point left of it): change
by the predicted branch point's distance from the AG, in D-70's exploratory bins (A-16); dots = exons, bars = bin
means (≥ 5 exons). Descriptive. (c) Exons with an acceptor T-run ≥ 20 nt and an experimentally mapped branch point
upstream of the run (BPHunter, lifted to GRCh38), in the held-out cancers (left; 52 exons, 31 of them in the rising
arm) and in the training cancers (right): rank of the change against rank of the gap from the run to the mapped
branch point, each with run length and offset (ranks) regressed out; line = least squares. For these exons the
branch point's distance from the AG is the offset of the run plus its length plus the gap, minus 1, so a cut on that
distance mixes run length and position with the gap. The scatter is descriptive; the tests are listed below.

- RP3a, rising arm minus distant branch point without a run (within strata): +0.0074, p = 5.0e-4, supported.
- The 36-nt cut: RP3b (rising arm minus long runs with a near predicted branch point) +0.0052, p = 0.16, not
  resolved; RP3d (mapped branch point ≥ 36 nt minus < 36 nt) +0.0012, p = 0.85, not resolved. Pre-specified
  sensitivity cuts for the rising arm (no q): ≥ 50 nt +0.0068, p = 0.091; > 60 nt +0.0077, p = 0.026. Bin means in
  (b): 18–29 nt: −0.0047 (18); 30–35 nt: +0.0063 (29); 36–49 nt: +0.0050 (40); 50–69 nt: +0.0082 (18); 70–101 nt:
  +0.0088 (34).
- RP3e, partial Spearman of the change with the gap given run length and offset, held out: +0.438, p = 0.004,
  supported (the branch-point test of the verdict); leaving out one exon at a time, 0.40 to 0.47. RP3c, the same
  with the predicted branch point (139 exons): +0.333, p = 5.0e-4, supported.
- RPD3b, the same test in the training cancers (62 exons, 150 MMRd tumours): +0.070, p = 0.61, not resolved. The gap
  association does not replicate there.
- Of the 6 largest rises among the 52 held-out exons in (c), 3 have a mapped branch point < 36 nt from the AG (2 of
  the 6 largest among the rising-arm exons).
