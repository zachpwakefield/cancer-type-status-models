**Figure U9. Part of the splicing models' subtype signal is tied to the site that sent the sample.** The probe of
D-20 on the subtype cohorts (D-47): one untuned classifier per cohort and view (top-2,000 variance, median impute,
standardise, L2 logistic regression C = 0.05, balanced), cross-validated with patients grouped and again with
tissue-source sites held out (5 folds, 3 seeds, pooled out-of-fold calls). (a) The loss of balanced accuracy,
grouped minus site-held-out, per cohort (patients; sites): grey = each splicing view, blue = their median, black =
GEX. Ordered by the splicing median. SARC is not probed: myxofibrosarcoma comes from 4 sites (rule ≥ 5). (b) The
loss against Cramér's V between site and class, computed from the labels alone.

- P1, the splicing views lose more than GEX: in 16 of 21 cohorts (median difference +0.010); sign test p = 0.027, q
  = 0.080 over D-47's three tests: not resolved.
- P2, the loss rises with Cramér's V: Spearman ρ = +0.12 for GEX (p = 0.60) and +0.18 for the splicing median (p =
  0.45): not supported.

Largest splicing losses: HNSC_HPV +0.098 (GEX +0.014), CESC +0.077 (GEX −0.001), THCA +0.077 (GEX +0.116), HNSC_mRNA
+0.073 (GEX +0.008). The probe's accuracies are those of the untuned classifier, not of the subtype models; per
view, scheme and seed in `results/subtype/site_probe/`.
