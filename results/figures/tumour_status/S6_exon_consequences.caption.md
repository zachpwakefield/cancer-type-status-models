**Figure S6. The models' cassette exons keep the reading frame and avoid folded domains.** Cassette exons in the top
100 of any tumour-vs-normal model (blue) or any class of the cancer-type model (dark blue), compared with a random
6,000 of all measured exons (grey).
- Length divisible by 3 is computed from the coordinates.
- The remaining features need an exact GENCODE v45 exon. Codons are taken from the best transcript (MANE > canonical
  > other) and overlapped with InterPro Pfam/CDD domains, MobiDB-lite disordered regions and TMHMM helices (D-17).

p: two-sided Fisher exact test against all measured exons, raw p. BH q is in the CSV, with its family. Exons skipped
or included more in tumours are in `cassette_exon_protein_features.csv`.
