# Assay pages

Drawn by splice-assay 0.4.0, not by this repository, so they keep the tool's layout rather than the
house style (docs/FIGURE_GUIDELINES.md). Each page shows, for one event in its three most promising
cancers: the tumour-vs-normal comparison (paired and unpaired), a Kaplan-Meier split at the median with
its log-rank test, the Cox model (PSI + host-gene expression + age + sex + stage, against stage I and
female), and a forest of every cancer. Markers: filled = p < 0.05; * = q < 0.05; † = the proportional-hazards test failed; ‡ = fewer than 20 deaths (low power).

The star is the tool's own FDR family: BH within one event over the cancers it was tested in. That is not the family used in results/TUMOUR_STATUS.md, which asks the harder question -- BH over every event tested inside one cancer (`cox_q_cohort` in `biology/08_survival/event_survival_cells.csv.gz`). An event can be starred here and not clear the other.

The events shown are the candidates associated with overall survival in the most
cancers (`biology/08_survival/event_dossiers.csv`). Nominal p values from observational cohorts: these
are associations, not causes. The pages' CSVs hold per-patient values and stay local (D-12).
The PNGs are the tool's pages at 150 dpi so the repository stays small; the PDF beside each one
is the exact figure, and the SVG stays local too.

- `RGS1_RI`
- `ENOSF1_HIT`
- `MACROH2A1_MXE`
- `GLIPR1_HIT`
- `HSPA12A_MXE`
