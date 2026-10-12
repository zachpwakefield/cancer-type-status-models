**Does the tumour-status model see some subtypes as less tumour-like?**

The tumour-status model of each cancer scores every sample, out of fold. For the subtype cohorts of the 8 cancers
with a status model, six predictions (D-50) name the subtype expected to look less tumour-like (row: id, cohort,
predicted lower vs comparison). a, mean 1 − P(tumour) in the predicted group (filled) and the comparison (open), for
GEX and the AS ensemble (log scale). b, AUROC = P(comparison patient's P(tumour) > predicted patient's), two-sided
Mann–Whitney; filled = BH q < 0.05 within the model. c, C1: Spearman ρ of P(tumour) with the stromal and the immune
score (mean log2 CPM of the D-18 marker genes) over each cancer's tumours, all of them, not only the
subtype-labelled ones.

- GEX: 2 supported, 0 against, 4 not resolved. AS ensemble: 6 supported, 0 against, 0 not resolved. GEX + AS: 6
  supported, 0 against, 0 not resolved.
- GEX's probabilities sit at the ceiling (median P(tumour) of the comparison groups 0.9982–0.9997), so ranks among
  them carry little; the AS ensemble averages eight models and spreads them out.
- C1, P(tumour) falls as composition scores rise: in GEX + AS 4 of 16 correlations are negative with q < 0.05 and 3
  positive (in HNSC and THCA): not supported as a common explanation.

Group sizes, means, medians, raw p and q are in the CSV; the per-class distributions in
`status_probability_by_subtype_class.csv`.
