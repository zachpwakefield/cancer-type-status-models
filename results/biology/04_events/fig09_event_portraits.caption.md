**What do the most shared tumour events look like, and how does the model read them?**

The six most recurrent events (fig07) whose paired change has the same sign in at least 80% of the cancers where
measured, one per gene.

Left: the median value in matched normal (orange) and tumour (blue) samples of each cancer. Cancers whose top 100
contains the event are in black, the others in grey. The heading gives the host gene, view, coordinates and, for
cassette exons, the protein consequence (GENCODE v45; D-17).

Right: SHAP dependence in the cancer where the event ranks highest. Each held-out sample's event value is plotted
against its SHAP value on the tumour log-odds. The value is the PSI; missing values were median-imputed by the
model.

The CSV gives the per-cancer medians and paired Wilcoxon p (raw). Per-sample values stay local (D-12).
