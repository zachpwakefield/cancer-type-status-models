**Does combining the modalities separate tumour from normal tissue better than expression alone?**

Left: out-of-fold errors as a share of each cancer's samples (n in brackets) for the expression model (black), the
AS ensemble (the mean of the 8 AS views' probabilities, blue) and an equal-weight fusion of the two (dark blue).
The fusion weights were fixed before any result was seen (D-14). Stars: exact McNemar, fusion against GEX, on
paired per-sample correctness (* p < 0.05, ** p < 0.01, *** p < 0.001); q (BH over the 45 tests) is in the CSV.

Right: errors pooled over all cancers. 'stacked' is a cross-fitted logistic regression on the 9 views'
probabilities (D-23): it is shown as a ceiling, not as the headline, because the base models behind its training
probabilities did see the samples it scores.

Pooled over the cohort, fusion corrects 27 of the expression model's 99 errors and
introduces 5 (exact McNemar p = 1.1e-4). No single cancer has enough errors to reach significance
on its own; the per-cancer tests are in the CSV.
