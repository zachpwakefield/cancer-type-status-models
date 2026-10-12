**Does the biology the literature predicts show up, in the data and in the models?**

Each row is a marker fixed in `config/subtype.yaml` before any model was fitted (D-41, D-45; sources in the CSV),
with the class the literature expects it to be higher in. Dot: AUROC of the gene's log2 CPM for that class over the
other (0.5 = no difference), computed on all patients of the two classes with no model involved; filled =
Mann-Whitney raw p < 0.05, open = not. Columns: the gene's best rank in the GEX model (over all patients or within
either class) and the best rank of any of its events across the eight splicing models, with the view (– = no event
of the gene was selected). For LGG and LGG_IDH, "1p/19q dosage" is the median over the arm's genes of each gene's
log2 CPM minus its cohort median, 1p and 19q averaged; the GEX column counts the model's top 100 genes that lie on
1p or 19q.

39 of 40 gene markers differ in the expected direction (p < 0.05). 27 rank in the top 100 of their GEX model; 8 have
an event in the top 100 of a splicing model. LGG: 7 of the top 100 GEX genes on 1p or 19q against 9.4 expected
(hypergeometric p = 0.84); LGG_IDH: 33 of the top 100 GEX genes on 1p or 19q against 9.4 expected (hypergeometric p
= 5.3e-11).
