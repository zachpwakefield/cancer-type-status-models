**Are the strongest splicing events splicing, or readouts of a neighbouring gene?**

For each histology and genomic/IHC cohort, the one splicing event with the largest class-specific mean |SHAP| across
the eight splicing models (rule fixed before the audit was run, D-41), labelled by cohort, host gene and view. Dots:
the event's median over the patients of each class, PSI or usage (0 to 1) or HIT index (−1 to 1); filled and named =
the class the event ranks for. Right: the D-27 audit. The event's segment is intersected with every other GENCODE
v45 gene's exons (≥ 50 bp): antisense, same-strand, inside a readthrough model, or no overlap. ρ: Spearman
correlation over the cohort's patients with the other gene's log2 CPM and with the host's. An overlapping event that
follows the other gene (|ρ| ≥ 0.5 and above |ρ host|) is a readout of it (red). "Not measured": the other gene is
not in the expression matrix, so the audit cannot clear the event.

4 of 19 are readouts of a neighbour (BRCA ZNF276 → FANCA, UCEC ERCC2 → KLC3, UCEC_TCGA ERCC2 → KLC3, HNSC_HPV VRK2 →
FANCL); 2 overlap a gene the audit cannot measure; 7 follow their own gene's expression with |ρ| ≥ 0.5. The full
audit of the top 5 events per class and view, and of every literature marker's best event, is
`results/subtype/interpret/event_audit.csv`.
