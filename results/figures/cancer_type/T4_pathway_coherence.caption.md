**Do the genes a splicing model relies on form pathways, as an expression model's do?**

The share of models whose top 200 SHAP genes are enriched (q < 0.05, BH within collection per model; GOseq/Wallenius,
D-16) for at least one Hallmark, Reactome or GO-BP term, by view.

Open bars use all genes, in each view's colour; a share of zero is printed as a number, because it draws no bar.

Left: the 15 tumour-vs-normal models. GEX reaches a pathway in 14 of 15 cancers; SE in 9; MXE in 6; the other AS
views in at most 4; RI in none. The tumour signal of splicing is mostly carried by genes that are not members of a
common programme.

Right: the 34 per-class explanations of the cancer-type models. Many AS hits there came from alternative 5'/3'
splice sites in cytosolic ribosomal-protein genes (RPS24, RPL13A, RPL3 and others), which are deeply covered and
so measured precisely in every tissue. Filled bars repeat the tests with RPL*/RPS* genes removed from both
ranking and universe; that sensitivity test was run only on the cancer-type AS views, so GEX and the left panel
have no filled bar. Most SE, MXE and A3SS classes keep enrichments. The remaining terms are lineage
biology (synapse and cell projection in PCPG, smooth-muscle contraction in PRAD, RHO GTPase cycle in DLBC,
extracellular matrix in THCA) and the splicing machinery's own transcripts (RNA splicing terms from A3SS in ACC,
KICH, LUAD and UVM). Top terms per class: `cancer_type_top_terms_by_class.csv`.
