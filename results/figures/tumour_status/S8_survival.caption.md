**Are the events the models rely on associated with how long patients live?**

Every event was tested with splice-assay 0.4.0 in each of the 15 cancers (6,152 tumours, 1,620 deaths): Cox on the
event's value with the host gene's own expression in the model, in the pre-specified base form and again with age,
sex and stage added (D-25). The models' top events are compared with the coverage- and spread-matched background of
D-22, so this is not a comparison of which events are easier to test.

Left: the share of testable event × cancer pairs with a nominal Cox p < 0.05, with 95% Clopper-Pearson intervals;
the three candidate sets are nested and ordered by how many models agree on the event. Odds ratios against the
background of the same arm, base / adjusted: every candidate 1.27 / 1.18, recurrent 1.56 / 1.28, top 10 for two or
more cancers 1.62 / 1.21 (Fisher exact, every set q < 0.05). The rate rises along that order.

The dashed line is the same pipeline run over outcomes shuffled between patients inside
each cancer: 5.6% (1,148 of 20,578 pairs), so both sets sit above what the gates, the host-gene
covariate and the test battery produce on their own.

Candidates inside one of the six regulons the models are enriched for (S4) reach p < 0.05 in 13.3% of their pairs
against 10.3% for the rest (OR 1.33, p = 0.005), so the three ways of ranking a candidate -- how many models agree,
whether it recurs, whether it sits in a coherent factor programme -- pick out the same events.

Right: the 12 candidate events associated in the most cancers, each in its strongest cancer: hazard ratio per SD
of PSI (per SD of the HIT index for HIT events) with its 95% CI, and (cancers with p < 0.05 / cancers tested).
Red: a higher value goes with worse survival.

Nominal p values from observational cohorts. BH inside each cancer's own family of ~1,600 tests (`cox_q_cohort`)
leaves 1.2% of candidate cells significant in the adjusted arm, no more than the background's 1.2%, and 94% of them
are KIRC or KIRP: the set-level excess is a broad shift of small associations, not an excess of strong ones.
`event_survival_cells.csv.gz` holds every statistic with both BH families (keyed by `event_index.csv`), and
`event_dossiers.csv` joins each event's recurrence, protein consequence, regulon membership and survival. PSI can
reflect tumour purity or composition, so these are associations, not causes, and not therapeutic claims.
