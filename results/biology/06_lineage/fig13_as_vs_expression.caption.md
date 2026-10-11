**Is the splicing signal for cancer type a proxy for expression?**

Left: for each AS view, the per-class top 20 events (union over the 34 classes) and 1,000 random measured events.
Each event's value is correlated (Spearman, 3,000 random samples) with its host gene's log2 CPM.
- Splice-site views (SE, A5SS, A3SS, MXE, RI): the top events correlate with host expression about as much as random
  events (median |ρ| 0.13-0.18 vs 0.10-0.15).
- Terminal-exon views (AFE, ALE, HIT): the top events track expression more (0.17-0.24 vs 0.03-0.05), as expected
  when an alternative promoter or poly(A) site also changes how much of the gene is made.

Right: for each class, the overlap between the top 100 GEX genes and the host genes of each AS view's top 100
events, within genes measured by both, summed over the 34 classes. The overlap is 1.3-3.0 times chance but small:
about one gene per class. The shared genes are lineage markers that are both expressed and spliced
lineage-specifically, such as GFAP and MBP (glioma), DUOXA2 (thyroid), MYH11 (smooth muscle), LGALS4
(gastrointestinal) and FXYD2 (kidney). Per-class genes: `top_gene_overlap_gex_vs_as.csv`.
