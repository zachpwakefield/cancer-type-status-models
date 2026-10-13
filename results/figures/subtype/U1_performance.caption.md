**Figure U1. Splicing matches expression on lineage and genome-state subtypes, and fusing the two adds little.**
Out-of-fold balanced accuracy (mean per-class recall) of one model per cohort and view, 22 cohorts and 7,552
patient-cohort pairs: the eight splicing views (grey; the bar spans them), expression (GEX, black), the AS ensemble
(the mean of the eight splicing views' class probabilities, blue) and GEX + AS fused 1:1 (dark blue). Fusion weights
were fixed before any result (D-40). Tick: chance, 1/classes, which is also the majority-class predictor's balanced
accuracy. Rows: cohorts, patients in brackets, grouped by how the label was made: histology (weak = a poorly
reproducible label, control = a lineage or site contrast), genomic or IHC (DNA, protein), and expression-defined
(mRNA clusters, so GEX is near-circular there; D-44).

Columns at right: fused minus GEX balanced accuracy, and the exact McNemar test of fused against GEX on paired
per-patient correctness, raw p. McNemar weighs classes by their size and balanced accuracy weighs them equally, so
the two can disagree (THCA). The CSV holds every point's 95% patient-bootstrap CI, macro F1, macro AUROC, the
paired-bootstrap CI of the fused − GEX difference and BH q (family: the subtype cohorts, per comparison).

Outside the expression-defined labels the AS ensemble is within 0.03 of GEX in 11 of 19 cohorts and above it in 4:
LGG (+0.026), CESC (+0.011), STAD (+0.006), ESCA (+0.006). It falls more than 0.10 short in 3: PCPG (−0.17),
PRAD_ETS (−0.15), UCEC_TCGA (−0.12). Fusing helps balanced accuracy with a CI above zero only in LGG (+0.023,
McNemar p = 0.064) and UCEC_TCGA (+0.036, McNemar p = 0.070).
