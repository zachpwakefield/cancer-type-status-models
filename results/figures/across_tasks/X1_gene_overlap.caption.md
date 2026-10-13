**Figure X1. The three tasks' models read overlapping genes, several-fold more than chance.** For each subtype
cohort, its cancer's tumour-status model (12 cohorts in 8 cancers have one) and the cancer-type model's importance
over that cancer's own samples: the top 100 genes of each model per view (a gene ranked by its best feature's mean
|SHAP|), compared within the genes both tasks could select (D-48). (a) Each dot is one cohort and view (comparisons
in brackets); bar = median; 0 = the overlap expected by chance (hypergeometric). (b) The five histology cohorts that
have a tumour-status model: Jaccard of the subtype list with the status list and with the type list; small dots =
views, large = median over views.

- Status – type (G1): median 7.8-fold the chance overlap over 72 comparisons, 94% with q < 0.05 (Wilcoxon of log2
  fold vs 0, p = 1.7e-13).
- Subtype – status (G1): median 5.0-fold the chance overlap over 117 comparisons, 87% with q < 0.05 (Wilcoxon of
  log2 fold vs 0, p = 1.3e-20).
- Subtype – type (G1): median 6.1-fold the chance overlap over 198 comparisons, 88% with q < 0.05 (Wilcoxon of log2
  fold vs 0, p = 1.4e-33).
- G2, for histology subtypes the subtype list is closer to type than to status: median Jaccard difference 0.000,
  positive in 44% of 45 cohort-views, Wilcoxon p = 0.24: not supported.

Genes in all three top-100 lists of a cancer are in `results/across_tasks/genes_in_all_three.csv`.
