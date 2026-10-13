**Figure U24. Within each cancer, the repeat-exon model ranks tumours almost as well as the exome MSI callers.**
Pre-specified (config/followup.yaml: msi_tool_benchmark, D-72) before any RNA model other than score G was scored
against MANTIS or MSIsensor. The pre-registered tests (BM1–BM6) were on the all-exon model; A-19 shows the
repeat-exon model (Figure U23) instead, so its comparisons here are exploratory. No measure is a gold standard, so
each exome measure is judged against the tumours on which the other two agree (positive: both call MSI; negative:
both call stable; others left out), and the repeat-exon model against the same tumours. Calls: MC3 ≥ 10 homopolymer
deletions; MANTIS > 0.6 (stable < 0.4); MSIsensor > 10 (stable < 4). Tumours with both caller scores (8,153 of
8,464; RNA scores all). Within-cancer AUROC ranks each positive only against negatives of its own cancer. E2: each
cancer left out of training (all 33); E1: trained on COAD, STAD and UCEC, tested on the other 30.

(a) Within-cancer AUROC in E2 (filled) and E1 (hollow); bars = 95% bootstrap intervals over tumours within cancer
and reference label, cut at 1. (b) Per cancer with ≥ 3 positives, E2, MANTIS judged by MC3 + MSIsensor (positives in
brackets). (c) Sensitivity (filled) and PPV (hollow) when the model flags as many tumours as the caller calls MSI,
by its own score pooled over cancers (E2; scores come from 33 models). (d) MC3 MMR-deficient tumours grouped by the
two callers' calls, and proficient tumours called stable by both: share in the top decile of the model's score
within their cancer (E2; n in brackets); dashed = chance.

- MANTIS judged: repeat-exon model E2 0.998, difference −0.0019 (−0.0041 to −0.0005); E1 1.000, difference −0.0003
  (−0.0010 to +0.0000).
- MSIsensor judged: E2 0.998, difference −0.0023 (−0.0049 to −0.0004); E1 1.000, difference +0.0000 (+0.0000 to
  +0.0000).
- MC3 count judged: E2 0.998, difference −0.0014 (−0.0036 to +0.0001); E1 1.000, difference +0.0012 (+0.0000 to
  +0.0030).

All six lower bounds lie above −0.05 (lowest −0.0049). Positives per reference, E2: 139, 143, 135; E1: 19, 17, 17.
The exome measures score 0.999–1.000 against each other. Panel (c), sensitivity / PPV: MANTIS 0.97 / 1.00 against
0.82 / 0.84; MSIsensor 0.94 / 0.99 against 0.83 / 0.86; MC3 1.00 / 0.95 against 0.84 / 0.80. Panel (d), median
within-cancer percentile: 0.94 (both), 0.88 (one), 0.86 (neither), 0.49 (proficient). Pre-registered, all-exon
model: 4 of 6 supported (it ranks 0.001–0.002 below MANTIS and MSIsensor, q ≤ 0.013), none against; non-inferior in
all six. Per-tumour scores are kept locally (D-12).
