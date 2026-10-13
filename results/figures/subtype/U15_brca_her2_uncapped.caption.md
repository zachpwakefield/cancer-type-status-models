**Figure U15. Without the ABSOLUTE ceiling, splicing matches expression on ERBB2 amplification but does not beat
it.** Pre-specified before the copy-number file was downloaded (config/followup.yaml: brca_her2_uncapped, D-60).
D-58 took amplification from ABSOLUTE, which never reports more than 7 copies at ERBB2. Here, the PanCanAtlas SNP6
copy ratio at ERBB2 is made absolute with ABSOLUTE purity and ploidy, without a ceiling; amplified = copy number ≥ 2
× ploidy, D-58's rule. Out-of-fold P(Her2) of the BRCA_PAM50 models; LumA and LumB tumours only.

(a) AUROC of P(Her2), GEX (grey) and AS ensemble (blue), for D-58's amplified tumours (median ploidy 2.09), all
tumours amplified by the uncapped call, and the newly amplified ones only (median ploidy 3.85). Each group is
compared with the luminal tumours neither call makes amplified; the middle group is compared with all non-amplified
tumours. (b) Amplified luminal tumours: P(Her2) against copy number / ploidy; filled = amplified by D-58 too, hollow
= newly amplified.

- HU1, GEX P(Her2) higher in amplified tumours: AUROC 0.75, q = 9.4e-8, supported.
- HU2, AS ensemble ahead of GEX: 0.81 against 0.75 (difference +0.062, q = 0.087), not resolved. D-58 had found the
  gap with q = 0.016.
- HU3, the independent part (newly amplified): AS 0.82, GEX 0.83, not resolved. The ceiling could not call any
  tumour with ploidy above 3.5, where 7 copies are below 2 × ploidy.
- HU4, within PAM50 Her2, GEX ahead of AS for amplified against not: −0.050, not resolved.
- HU5, AS P(Her2) higher in non-amplified Her2 than in amplified luminal tumours: AUROC 0.85, supported.
- HU6, GEX P(Her2) rises with the amplification level: ρ = 0.47, supported. HU7, splicing follows it less: AS ρ =
  0.45, not resolved.
- Knockout with this amplified group (exploratory): the 17q12 inputs set to their mean lower GEX from 0.753 to
  0.720; no splicing view moves by more than 0.005.

q over family HU (7 tests). Per-tumour values are kept locally (D-12).
