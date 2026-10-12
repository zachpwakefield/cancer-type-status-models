**Which other cancers does a subtype resemble, in the cancer-type model?**

The pan-cancer type model gives every tumour a probability for each of 34 classes. Thirteen predictions (D-49,
sources in config/followup.yaml) name a subtype, the projects it should resemble and the comparison class(es): row =
id, cohort, the predicted-higher vs predicted-lower group, → target projects. Measure: the patient's summed
out-of-fold probability over the targets; two-sided Mann–Whitney, AUROC = P(higher-group patient > lower-group
patient). a, the mean of that probability in each group, GEX + AS (log scale). b, AUROC per model; filled = BH q <
0.05 within the model. T11 (cortical admixture in PCPG → ACC) is a positive control.

- GEX: 11 supported, 0 against, 2 not resolved. AS ensemble: 13 supported, 0 against, 0 not resolved. GEX + AS: 12
  supported, 0 against, 1 not resolved.
- Sensitivity, not pre-specified (×): the targets' share of the probability not given to the patient's own project.
  Most predictions hold. T6 reverses (GEX + AS AUROC 0.29): basal-like tumours carry more ovarian probability
  because they carry more off-target probability overall, not a larger ovarian share. T5 in the AS ensemble weakens
  to 0.57.

Group sizes, means, medians, raw p and q for every model and measure are in the CSV.
