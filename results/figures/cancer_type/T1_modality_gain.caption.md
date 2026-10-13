**Figure T1. No splicing view matches expression alone, but their ensemble does, and fusing the two beats both.**
Out-of-fold accuracy over 9,729 samples and 34 classes (33 TCGA projects and a pooled normal class) for each
single-view model (grey; GEX in black), the AS ensemble (the mean of the 8 AS views' probabilities, blue), an
equal-weight fusion of GEX and the AS ensemble, and a cross-fitted stacker on all 9 views (diamond). The fusion
weights were fixed before any result was seen (D-14); the stacker is a ceiling, not a headline (D-23).

Whiskers: 95% case-level bootstrap CIs (1,000 resamples). Brackets: exact McNemar on paired per-sample correctness,
raw p (two tests; q in the CSV).

No single splicing view reaches expression, but their ensemble does. Fused with expression it corrects 109 of its
errors and introduces 38. Macro recall and macro F1 are in the CSV.
