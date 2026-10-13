**Does the sequence code find mismatch-repair-deficient tumours in 30 cancers never used to define it?**

Pre-specified (config/followup.yaml: msi_code_pancancer, D-66) before any tumour outside UCEC, STAD and COAD was
read for MSI. Tumours: one per patient in the 30 other TCGA projects of the pan-cancer SE matrix, with MC3 calls
(7,685). Mismatch-repair deficient (MMR-deficient) = ≥ 10 MC3 homopolymer deletions, the D-61 count; in the three
training cohorts that cut separates MSI from the rest almost perfectly (129 of 129 against 2 of 571). 50 tumours
qualify. Hypermutated proficient = fewer than 10 such deletions and ≥ 500 SNVs (408). Score A = D-52's repeat-exon
score; score G = the code: exons after a T-run of 10–19 nt (R_down, read as falling) against exons after a run ≥ 20
nt whose branch point lies ≥ 36 nt from the AG (R_up, read as rising), each against matched background exons.
Neither score was re-fitted; each was computed within each cancer and ranked there.

a, AUROC for MMR-deficient tumours scoring lower (score A hollow, G filled): pooled over the 30 cancers, without the
three tissue relatives of the training cohorts, against hypermutated proficient tumours only, and per cancer with ≥
3 MMR-deficient tumours (MMR-deficient of all tumours in brackets). b, per-exon shift (MMR-deficient minus the rest,
background-adjusted percentile) in the 30 cancers against the mean MSI shift of the three training cohorts; R exons
measured in both (grey = R_down, blue = R_up, light grey = runs ≥ 20 in neither arm). c, MMR-deficient tumours of
the 30 cancers: score G percentile against their homopolymer deletions (exploratory).

- PC1, score A: AUROC 0.55, q = 0.25, not resolved.
- PC2, score G: AUROC 0.83, q = 2.2e-15, supported; PC3, G ahead of A by 0.28, supported.
- PC4, without READ, ESCA and UCS: 0.81; PC5, against hypermutated proficient tumours: 0.85. Both supported: the
  code reads mismatch-repair deficiency, not mutation burden.
- PC6, exon shifts replicate: ρ = 0.51 over 1753 R exons; PC7, R_up exons rise more than R_down exons: AUROC 0.70;
  PC8, the branch-point distance, adjusted: ρ = 0.10. All supported.
- Per cancer (exploratory, 6 cancers): score G 0.65–0.93, score A 0.40–0.72. By arm (A-15, exploratory): R_down
  0.61, R_up 0.73 (higher in MMR-deficient).

The DNA call is a proxy from exome calls, not a clinical MSI test. q over family PC (8 tests). Per-tumour values are
kept locally (D-12).
