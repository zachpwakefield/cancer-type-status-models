**Does splicing add information that expression lacks, for cancer type?**

Out-of-fold accuracy over 9,729 samples and 34 classes (33 TCGA projects + pooled normal) for each single-view model
(grey; GEX in black), the AS ensemble (mean of the 8 AS views' probabilities, blue) and an equal-weight fusion of GEX
and the AS ensemble (dark blue). The fusion weights were fixed before looking at any result (D-14). Whiskers: 95%
case-level bootstrap CIs (1,000 resamples).

Brackets: exact McNemar tests on paired per-sample correctness, raw p (two tests; no FDR correction). The values
are in `fig01_accuracy_by_view_tests.csv`. Single AS views are less accurate than GEX, and the AS ensemble is
level with it. Fused with GEX, the AS views correct 109 GEX errors and introduce
38.
