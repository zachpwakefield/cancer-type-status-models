**Do the splicing regulators the literature names differ between subtypes, and track the subtype's exons?**

Six predictions fixed before the data were looked at by class (config/followup.yaml, D-46). a, a factor's log2 CPM
in each class: dot = median, bar = interquartile range, patients in brackets; blue = the class predicted lower.
AUROC = P(predicted-higher patient > predicted-lower patient), two-sided Mann–Whitney.

- R1, ESRP1 lower in biphasic than epithelioid mesothelioma: AUROC 0.55, p = 0.52, q = 0.52, not resolved (ESRP1 is
  barely expressed in either class).
- R2, ESRP2, same contrast: AUROC 0.82, p = 1.8e-5, q = 3.5e-5.
- R3, ESRP1 lower in the mesenchymal HNSC class than in the other three: AUROC 0.76, p = 5.4e-11, q = 1.6e-10.
- R4, PTBP1 lower in leiomyosarcoma (LMS) than in the other sarcomas: AUROC 0.66, p = 2.4e-5, q = 3.5e-5.

b, the event's PSI by PTBP1 quintile, the quintiles cut within LMS and within the other sarcomas separately (dot =
median, bar = interquartile range, about 20–26 patients per point). Test: partial Spearman of PTBP1 with the PSI
over all 230 sarcomas, adjusting for class, stromal and immune scores; predicted sign = opposite to the LMS shift.

- R5, ACTN1: LMS shift −0.62, ρ = +0.43, p = 1.4e-11 (within LMS alone ρ = +0.57). Agrees with Southby et al. 1999:
  PTB represses the smooth-muscle exon.
- R6, TPM1: LMS shift +0.53, ρ = −0.22, p = 7.1e-4. It matches the pre-specified sign, but that sign misread the
  source (A-5): in TPM1, PTB and MBNL1 repress exon 3, the non-smooth-muscle exon (Gooding et al. 1998, 2013). So
  the data run against the published mechanism here.

Sensitivity (A-6, exploratory): adding each sample's mean MXE PSI (D-24) leaves ρ at +0.43 and −0.22; PTBP1 ranks 2
of 70 factors for ACTN1 and 12 of 70 for TPM1, where MBNL1 leads (MBNL1 +0.35). BH q over R1–R6 is in the CSV. The
scan of all 70 factors in every cohort is `results/subtype/regulators/factor_by_class.csv`.
