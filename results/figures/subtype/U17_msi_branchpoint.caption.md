**Figure U17. The branch point's distance from the 3′ splice site sets the direction of the MSI shift.**
Pre-specified (config/followup.yaml: msi_branchpoint, D-63). Exons whose acceptor carries a T-run ≥ 10 (R exons,
D-52), read from the GRCh38 primary assembly. Branch point = the A nearest the 3′ splice site that lies upstream of
the run, within 100 nt, and matches yUnAy; found for 83% of 1,580 exons. MSI shift = D-52's background-adjusted
difference in inclusion percentile, MSI minus the rest.

(a) Distance from the branch point to the 3′ss, by run length (1,133 exons with runs of 10–19 nt, 179 with ≥ 20 nt);
grey band = 19–23 nt, the optimum for a 3′ss AG (Chua and Reed 2001). (b) Share of R exons more included in MSI, by
tertile of that distance within run classes, per cancer; bars = Wilson 95% interval. (c) The tests: Spearman of the
MSI shift with the branch-point distance after both are rank-residualised on run-length bins (10–11, 12–13, 14–15,
16–19, 20–24, ≥ 25) and mean inclusion; p by permuting the distance within run-length bins (1,000; the floor is
0.001). Hollow = the same with the run's far end in place of the branch point (exploratory).

- BP1–BP4, an exon whose branch point lies further from the 3′ss is more included in MSI, beyond its run length:
  STAD ρ = +0.161, COAD ρ = +0.115, UCEC ρ = +0.109, mean over ≥ 2 cancers ρ = +0.156; all q = 1.0e-3, supported.
- Mean over ≥ 2 cancers, runs of 10–19: 0.29 of exons with the nearest branch points move up, 0.46 with the
  farthest. Runs ≥ 20: 0.48 against 0.73. The long-run reversal of D-53 sits in the exons whose branch point is far.

The run's far end carries nearly the same information as the branch point, because the branch point must lie beyond
it; these tests cannot separate the two. The branch point is a motif call, not a measured position. q over family BP
(4 tests).
