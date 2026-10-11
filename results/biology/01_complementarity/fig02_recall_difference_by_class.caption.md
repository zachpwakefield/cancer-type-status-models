**Which cancer types does splicing resolve better than expression, and which worse?**

The difference in out-of-fold recall between the AS ensemble and the GEX model for each class (n = samples).
Positive (blue): splicing recovers more of the class.

Splicing gains the most for classes that GEX confuses with an anatomical neighbour (COAD/READ, ESCA/STAD, GBM/LGG).
Of these, only the glioma separation survives holding out whole tissue-source sites (fig15). The colon/rectum gain
does not, and should not be read as biology.

GEX is better for pooled normal tissue, uterine carcinosarcoma, and squamous and urothelial tumours (BLCA, HNSC,
CESC). Values: `cancer_type_recall_gex_vs_as.csv` (with the fused model's recall).
