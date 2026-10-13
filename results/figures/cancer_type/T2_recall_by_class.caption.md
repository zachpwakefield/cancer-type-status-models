**Figure T2. Fusion keeps splicing's gains on the classes expression confuses and repairs most of its losses.**
Out-of-fold recall of each of the 34 classes (n in brackets), as a difference from the expression model: the AS
ensemble (light blue) and the equal-weight fusion of the two (dark blue). Positive means the model recovers more of
that class than expression does.

Fusion raises recall for 19 classes and lowers it for 6. It is asymmetric in a useful way: across the 11 classes
where splicing alone is worse than expression, the median shortfall shrinks from −0.034 to −0.002, while across the
17 classes where splicing is better it keeps most of the gain (0.016 to 0.012). The pooled normal class is the
clearest case: −0.102 for splicing alone, −0.007 fused. So the overall gain (T1) is not bought at the cost of
particular classes.

Splicing gains most for classes expression confuses with an anatomical neighbour (COAD/READ, ESCA/STAD, GBM/LGG). Of
those, only the glioma separation survives holding out whole tissue-source sites (T7): the colon/rectum gain does
not and should not be read as biology.

Expression is better for the pooled normal class, uterine carcinosarcoma, and squamous and urothelial tumours (BLCA,
HNSC, CESC). Every model's per-class recall is in the CSV.
