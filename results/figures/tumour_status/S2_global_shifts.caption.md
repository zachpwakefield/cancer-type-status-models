**Do tumours shift RNA processing in a common direction relative to their own normal tissue?**

Each cell is the median over matched patients (n in brackets) of the patient's mean paired tumour − normal change,
averaged over every event measured in that cancer, so there is no SHAP selection and no circularity. The columns
say what an increase means:
- RI: more intron retention.
- A5SS/A3SS: more use of the long exon form, which includes extra intronic sequence.
- SE: more cassette inclusion.
- MXE: more of the upstream exon.
- ALE/AFE: the used last or first exons move downstream (a 0-1 position index along transcription).

Blue is higher in the tumour and orange higher in the normal tissue, the house tumour/normal colours. Stars:
two-sided Wilcoxon signed-rank test across patients, raw p (* < 0.05, ** < 0.01, *** < 0.001). BH q (family: 15
cancers × 7 views) is in the CSV.

Retention and long splice-site forms rise in most cancers. Last exons move upstream and never downstream. BRCA is
the exception for RI. Caveat: RNA integrity and nuclear-RNA content also raise intronic reads, and RIN is not
available (D-15).
