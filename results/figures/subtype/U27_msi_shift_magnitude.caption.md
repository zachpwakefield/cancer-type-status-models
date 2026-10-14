**Figure U27. A minority of exons move by a tenth of PSI or more; the typical exon moves by a few hundredths.**
Pre-specified (config/followup.yaml: msi_magnitude, D-75) before any shift was measured in PSI units. All shifts are
within cancer: the MMR-deficient tumours (MC3 ≥ 10 homopolymer deletions; 150 in COAD, STAD and UCEC, 50 in 30
held-out cancers) against the proficient tumours of the same cancer, on cassette exons measured in ≥ 80% of the
cancer's tumours. Panels a, b and e use the exons whose proficient median PSI lies between 0.1 and 0.9 (A-20), where
a shift is not capped by 0 or 1 (COAD: 322 R_down, 24 R_up, 5,209 background exons).

(a) Share of exons whose shift (median PSI of MMR-deficient minus proficient tumours) lies below x (left) or above x
(right), the three training cancers pooled over exon × cancer; grey lines at ±0.10. (b) Share of exons dropping (▼)
or rising (▲) by at least 0.10, by acceptor T-run length and branch-point (bp) distance; filled = training cancers,
hollow = held-out cancers (an exon's shift there is the median over held-out MMR-deficient tumours, ≥ 20, of its
deviation from the proficient median of each tumour's cancer); bars = 95% Wilson intervals; n = training exons. (c)
PSI of four exons in the proficient (grey) and MMR-deficient tumours (line) of the training cancers; dashed lines =
medians. MRE11 exon 5 is the known MSI target; METTL8 and BOD1L1 have the largest mean shifts among the R_down and
R_up exons measured in all three training cancers; CCNDBP1 is the MSI event of Figure U7. (d) Per tumour, the R_down
exons with PSI at least 0.10 below the proficient median (left) and the R_up exons at least 0.10 above it (right),
minus the count of the cancer's median proficient tumour; boxes = quartiles, whiskers = 5th–95th percentiles. (e)
Where the MMR-deficient tumours of the training cancers fall among the proficient tumours of their cancer: density
of mid-rank percentiles over tumour × exon pairs (A-20); lines as in (a); 1 = no difference. (f) Mean arm shift (a
tumour's mean deviation over its R_down or R_up exons, minus its mean over background exons of the same PSI decile)
by purity fifth of the MMR-deficient tumours, and for proficient tumours in the same bins; bars = 95% bootstrap
intervals.

- MG1, R_down arm more negative at higher purity among MMR-deficient tumours (partial Spearman ρ): −0.150, p =
  0.047, q = 0.12, not resolved.
- MG2, R_up arm more positive at higher purity: −0.051, p = 0.51, q = 0.63, not resolved.
- MG3, purity × MMR interaction on the R_down arm: −0.001, p = 0.72, q = 0.72, not resolved.
- MG4, larger lower-tail excess in R_down than background exons: −0.027, p = 2.0e-4, q = 1.0e-3, against.
- MG5, larger upper-tail excess in R_up than background exons: +0.024, p = 0.075, q = 0.13, not resolved.

0 of 5 supported, 1 against (family MG). The tail excess of MG4 and MG5 is the share of MMR-deficient tumours beyond
the proficient 5th (95th) percentile minus the share expected if all moved by the exon's median shift; that
expectation lets proficient values pass 0 and 1 (A-20). Lower-tail excess of background exons runs from −0.072
(proficient median < 0.1) to +0.001 (≥ 0.98); where the bound is not reached (0.5–0.9, 0.9–0.98 and ≥ 0.98), R_down
exons score −0.001, +0.005 and +0.011 against −0.002, −0.001 and +0.001; upper tail below 0.1 and 0.1–0.5, R_up
+0.001 and −0.008 against −0.002 and −0.020. Exons dropping ≥ 0.10, R_down against background: COAD 18% against
3.8%; STAD 16% against 3.6%; UCEC 9.4% against 4.0%; held out 5.5% against 1.6%. Rising ≥ 0.10, R_up against
background: COAD 17% against 8.5%; STAD 8.3% against 4.4%; UCEC 12% against 6.2%; held out 3.8% against 1.3%. Median
shift of R_down exons: −0.021, −0.028, −0.015 and −0.005 (COAD, STAD, UCEC, held out); 5th percentiles −0.20, −0.21,
−0.13 and −0.10. Panel (c), proficient → MMR-deficient medians: MRE11 1.00 → 0.96, METTL8 0.58 → 0.22, BOD1L1 0.13 →
0.43, CCNDBP1 0.04 → 0.25; 52% of the MMR-deficient tumours fall below the proficient 5th percentile of MRE11 exon 5
(mean over the three cancers). Panel (d), MMR-deficient median (IQR), R_down: COAD +35 (+22 to +46), STAD +32 (+19
to +50), UCEC +12 (+1 to +24), held out +16 (−2 to +28); R_up: +7 (+3 to +9), +5 (+2 to +7), +2 (+0 to +5), +4 (−1
to +7). Panel (e), lowest 5% of the proficient range: R_down 2.10, background 1.09; highest 5%: R_up 1.59,
background 0.83; held out 1.40 against 1.12 and 1.30 against 0.89. Panel (f), MMR-deficient means across the purity
fifths: R_down −0.73 to −0.59, R_up 0.84 to 1.58 (× 100). MMR status explains a median 0.5%, 0.4%, 0.4% of an R_down
exon's PSI variance (COAD, STAD, UCEC). Per-tumour values are kept locally (D-12).
