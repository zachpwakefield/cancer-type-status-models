**Which splicing factors' programmes do the models' events belong to?**

A factor's *core regulon* is the events whose PSI tracks its expression within tumours — partial Spearman
|rho| >= 0.3 at p < 0.05 after rank-residualising both on a stromal score, an immune score and the sample's own mean
PSI in that view — reproducibly, with the same sign, in at least 5 of the 15 cancers (D-24).

Left: for every factor with at least 15 marked events, the odds that one of the models' top events is in its core
regulon, against the coverage- and spread-matched background of D-22. Only factors with q < 0.05 (BH over the
32 factors tested) are shown; the regulon's size is in brackets. Odds ratios use a 0.5 correction so that an
empty cell stays finite — PTBP1's 15 candidate events against 0 background events is otherwise infinite.

Right: every factor with at least 20 core events. Large regulons are not separate programmes: they overlap each
other heavily (Jaccard up to 0.7), and the 1,074 distinct events carry 7,529 memberships, about 7 factors each.

The two panels say the same thing. The models' events avoid the big shared programme of nuclear RNA-processing
factors (SRSF6, SRSF4, RBM4, SRSF10, odds 0.26-0.50) and concentrate in the
small, distinct regulons of tissue and differentiation regulators (PRPF8, QKI, HNRNPF, HNRNPU, RBM47, PTBP1).

These are correlations inside tumours: a factor may be co-regulated with its apparent targets rather than drive
them. `regulon_core_events.csv` lists the members, `global_axis_factors.csv` the factors most tied to the
per-sample PSI level that was residualised out.
