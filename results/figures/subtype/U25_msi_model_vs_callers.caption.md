**Figure U25. The model agrees with the exome measures on its confident calls; its own extra calls are borderline.**
Exploratory view of D-73 (config/followup.yaml: msi_rna_only; tests RO1–RO5 in the table). The model is D-71's
repeat-exon logistic regression (A-19), each cancer scored by a model trained on the other 32; it calls MSI when its
decision function is above 0. Exome calls: MC3 ≥ 10 homopolymer deletions, MANTIS > 0.6, MSIsensor > 10. The 8,153
tumours with both caller scores.

(a) Tumours called MSI by each combination of the four (dots below; bars coloured by agreement class); 7,920 tumours
are called by none. (b) Tumours called by anyone, per cancer with ≥ 3 such tumours, by class. (c–e) The model's
score against MC3 homopolymer deletions, the MANTIS score and the MSIsensor score (deletions and MSIsensor on a
symmetric log axis); dashed = calling thresholds (MANTIS and MSIsensor: stable below the lower, MSI above the
upper). Rings = three model-only tumours named in D-73: the highest score (STAD), MC3 calls in all four MMR genes
and POLE (PAAD) and the lowest MLH1 expression (BLCA).

- Agreement: the model and all three exome measures call 115 tumours; the model and one or two call 18; the model
  alone 32; exome measures alone 68.
- Model-only calls sit just above the model's threshold (median score 1.0, against 9.4 when all agree), and below
  every exome threshold (median MANTIS 0.31, MSIsensor 0.16, MC3 deletions 0).
- Exome-only calls: 20 by all three measures, 19 of them UCEC; 42 by one measure alone, 38 of them by MC3 alone.
  Those are hypermutated (median 616 SNVs, against 74 in tumours no one calls) with half the deletions per SNV of
  triple calls (0.039 against 0.080), and both callers place them in their between zones (median MANTIS 0.47,
  MSIsensor 4.1).
- Candidates: STAD score 6.0, 9 deletions, MANTIS 0.45; PAAD score 1.4, MSIsensor 4.5, 21,427 SNVs; BLCA score 2.7,
  MLH1 at the 1st percentile, MSIsensor 1.2.

Per-tumour points are kept locally (D-12).
