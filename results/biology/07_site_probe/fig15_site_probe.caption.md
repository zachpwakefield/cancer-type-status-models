**Is the separation of neighbouring cancer types a tissue-source-site effect?**

In TCGA every tissue-source site contributes to exactly one project (COAD and READ share 0 of 38 sites; GBM and LGG
0 of 40), so a model could separate neighbouring projects by site handling. For six pairs of related classes, the
same fixed, untuned probe is cross-validated two ways:
- by patient (open circles), as in phase 1;
- by site (filled), so that every test site is unseen (D-20).

The probe is top 2,000 variance → median impute → standardise → L2 logistic regression, all fitted inside training
folds. GEX is grey, AS views blue; points are the mean of 3 fold seeds.
- **Colon vs rectum:** no AS view separates the two beyond chance at unseen sites (0.48-0.57). GEX keeps 0.66. The
  AS ensemble's better colon/rectum calls in the main models (fig02, fig12) should therefore not be read as biology.
- **Oesophagus vs stomach:** AS views lose 0.08-0.14 of balanced accuracy at unseen sites, GEX 0.06. Part of the
  splicing advantage there is site-specific.
- **Glioma grade, kidney subtypes and lung histology:** the separation survives unseen sites in every view (drop
  ≤ 0.015). The near-perfect splicing separation of GBM from LGG is robust.

The probe is diagnostic: it does not re-estimate the phase-1 models. Values per seed: `site_probe_by_seed.csv`.
