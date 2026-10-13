**Figure U24. Within each cancer, the splicing model ranks tumours almost as well as the exome MSI callers.**
Pre-specified (config/followup.yaml: msi_tool_benchmark, D-72) before any RNA model other than score G was scored
against MANTIS or MSIsensor. No measure is a gold standard, so each exome measure is judged against the tumours on
which the other two agree (positive: both call MSI; negative: both call stable; others left out), and the splicing
models against the same tumours. Calls: MC3 ≥ 10 homopolymer deletions; MANTIS > 0.6 (stable < 0.4); MSIsensor > 10
(stable < 4). Tumours with both caller scores (8,153 of 8,464; RNA scores all). Within-cancer AUROC ranks each
positive only against negatives of its own cancer. Splicing models come from D-71: left out of training by cancer
(E2, all 33 cancers) or trained on COAD, STAD and UCEC (E1, the other 30).

(a) Within-cancer AUROC in E2 (filled) and E1 (hollow); bars = 95% bootstrap intervals over tumours within cancer
and reference label. (b) Within-cancer AUROC per cancer with ≥ 3 positives, E2, MANTIS judged by MC3 + MSIsensor
(positives in brackets). (c) Sensitivity when each splicing model flags as many tumours as the caller calls MSI, by
its own score pooled over cancers (E2; scores come from 33 models). (d) MC3 MMR-deficient tumours grouped by the two
callers' calls, and proficient tumours called stable by both: share in the top decile of the all-exon model's score
within their cancer (E2; n in brackets); dashed = chance.

- BM1, MANTIS, E2: within-cancer AUROC(all-exon model) − AUROC(MANTIS) = −0.0014 (95% CI −0.0025 to −0.0005), q =
  0.001, supported; all-exon model 0.999, MANTIS 1.000.
- BM2, MSIsensor, E2: −0.0023 (95% CI −0.0061 to −0.0002), q = 0.013, supported; all-exon model 0.998, MSIsensor
  1.000.
- BM3, MC3 count, E2: −0.0006 (95% CI −0.0017 to +0.0002), q = 0.20, not resolved; all-exon model 0.999, MC3 1.000.
- BM4, MANTIS, E1: −0.0021 (95% CI −0.0032 to −0.0011), q = 0.001, supported; all-exon model 0.998, MANTIS 1.000.
- BM5, MSIsensor, E1: −0.0021 (95% CI −0.0032 to −0.0011), q = 0.001, supported; all-exon model 0.998, MSIsensor
  1.000.
- BM6, MC3 count, E1: −0.0009 (95% CI −0.0027 to +0.0014), q = 0.39, not resolved; all-exon model 0.998, MC3 0.999.

The all-exon model is non-inferior in all six at the pre-set margin of 0.05 AUROC (lowest bound −0.0061). E2
positives per reference: 139, 143 and 135 (E1: 19, 17, 17). Score G: 0.957, 0.955, 0.959 (E2). Panel (c),
sensitivity / PPV, E2: MANTIS 0.97 / 1.00 against 0.81 / 0.83; MSIsensor 0.94 / 0.99 against 0.80 / 0.84; MC3 1.00 /
0.95 against 0.84 / 0.80 (all-exon model); pooled AUPRC 1.00 for MANTIS against 0.89. Panel (d), median
within-cancer percentile: 0.93 (both callers, 135), 0.89 (one, 12), 0.86 (neither, 50), 0.49 (proficient). No
MC3-proficient tumour was called MSI by both callers. q over family BM (6 tests). Per-tumour scores are kept locally
(D-12).
