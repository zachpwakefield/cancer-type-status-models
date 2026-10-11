**Which programmes do the tumour-vs-normal models draw on?**

For every cancer (columns) and view, the 200 genes with the largest SHAP importance were tested for enrichment of
the MSigDB Hallmark sets against all genes measurable in that view and cohort. The test is GOseq/Wallenius,
correcting for the number of events per gene (D-16).
- A dot marks a significant set in at least one AS view; its size is the number of AS views (of 8).
- A ring marks significance in the GEX model.
- Significance is q < 0.05, BH within the Hallmark collection of each model.

GEX models consistently use stromal and microenvironment programmes (EMT, coagulation, KRAS, myogenesis). The AS
models reach the same programmes far less often, mostly through EMT and myogenesis genes. So their tumour signal is
largely not organised by pathway membership (fig06).
