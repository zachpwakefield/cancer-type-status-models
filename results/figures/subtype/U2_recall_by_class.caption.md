**Figure U2. Splicing loses most on small classes defined by one programme.** Out-of-fold recall of every class
(patients in brackets) for GEX (black), the AS ensemble (blue) and GEX + AS fused 1:1 (dark blue). (a) Histology
cohorts. (b) Genomic and IHC cohorts, then the expression-defined ones, in the order of U1. 95% patient-bootstrap
CIs and every single view are in the CSV and in `recall_by_class.csv`.

Splicing's five largest losses against expression (recall, GEX → AS ensemble) are all in classes of fewer than 40
patients: UCEC_TCGA POLE (24): 0.54 → 0.08; SARC MFS (24): 0.67 → 0.29; PRAD_ETS ETV1 (28): 0.93 → 0.57; PCPG
Paraganglioma (30): 0.83 → 0.57; THCA Tall cell (35): 0.71 → 0.51. Its largest gains: UCEC_TCGA CN-high (83): 0.84 →
0.94; THCA Classical (351): 0.56 → 0.66; SARC UPS (47): 0.64 → 0.72; STAD Intestinal (157): 0.64 → 0.71. Over the 12
classes with fewer than 40 patients, the AS ensemble's recall is lower than GEX's in 9.
