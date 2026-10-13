**Are the BRCA models' genotype-like splicing events inherited?**

D-57 (config/followup.yaml: brca_germline) called an event genotype-like when, without labels, its PSI fell into
three well-separated groups (low, mid, high) with mid-group normals. A-13 adds a check from tumour–normal pairs, and
A-14 (exploratory) adds pair agreement. One tumour and one normal per patient (112 pairs). An inherited event keeps
its group in the tumour, apart from loss of heterozygosity (mid to low or high); a homozygous normal cannot become
the opposite homozygote.

a, b, each normal's group against its own tumour's group (counts; shade = share of the normal's row) for CAST SE, a
germline event, and TMEM101 A5SS, the genotype-like event with the largest importance in the BRCA models. c,
cumulative distribution of Cohen's κ between normal and tumour group, for all 1110 checked events (one per gene),
the 82 that pass A-13 (opposite-homozygote flips in ≤ 5% of at least 20 homozygous pairs) and the 45 that rank among
the models' top events (legend: all, A-13, top). d, share of tumours in the high group, LumA (x) against Basal (y),
per event, coloured as in c (grey = the rest); dashed line = equal frequency. An inherited allele differs between
subtypes only through ancestry.

- Only 82 of 1110 events pass A-13 (median flip rate 0.20), and 20 have κ > 0.6. Of the 45 genotype-like events
  among the models' top events, 3 pass A-13, and their κ is at most 0.24.
- Median |LumA − Basal| frequency: top events 0.18, A-13 confirmed 0.03, all 0.06. The models' genotype-like events
  behave as subtype programmes, not as inherited alleles.
- D-57, all events: L1 supported, L2 supported, L3 supported, L4 supported, L5 not resolved, L6 supported, A3
  supported. Germline-confirmed only (A-13): L1 supported, L2 supported, L3 supported, L4 not resolved, L5 not
  testable, L6 not resolved, A3 supported.
