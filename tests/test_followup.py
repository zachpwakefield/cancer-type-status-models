"""Follow-up layer (config/followup.yaml; docs/SUBTYPE_METHODS.md D-46, D-47 and
docs/ACROSS_TASKS_METHODS.md D-48..D-51): the pre-specified predictions name cohorts, classes and
projects that exist, and the site probe's helpers behave.

Run: OMP_NUM_THREADS=1 python -m unittest discover -s tests
"""
from __future__ import annotations

import importlib.util
import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import torch  # noqa: E402,F401  (load torch's OpenMP runtime before numpy's; see scripts/_common.py)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from ctsm.subtype import cohort_specs, subtype_spec  # noqa: E402

FOLLOWUP = yaml.safe_load(open(REPO / "config" / "followup.yaml"))
SPECS = cohort_specs(subtype_spec())
PROJECTS = set(pd.read_csv(REPO / "results" / "data" / "task_scope.csv")["project"])


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), REPO / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestFollowupConfig(unittest.TestCase):
    def predictions(self):
        ct = FOLLOWUP["cross_task"]
        return (FOLLOWUP["subtype_regulators"]["predictions"]
                + ct["type_probability_by_subtype"]["predictions"]
                + ct["status_probability_by_subtype"]["predictions"])

    def test_ids_unique(self):
        ids = [p["id"] for p in self.predictions()]
        self.assertEqual(len(ids), len(set(ids)))

    def test_cohorts_and_classes_exist(self):
        for p in self.predictions():
            self.assertIn(p["cohort"], SPECS, p["id"])
            classes = set(SPECS[p["cohort"]]["classes"])
            for key in ("higher_in", "lower_in"):
                for c in p.get(key, []):
                    self.assertIn(c, classes, f"{p['id']}: {c}")
            self.assertFalse(set(p.get("higher_in", [])) & set(p.get("lower_in", [])), p["id"])

    def test_type_targets_are_projects_other_than_the_cohorts_own(self):
        for p in FOLLOWUP["cross_task"]["type_probability_by_subtype"]["predictions"]:
            for t in p["targets"]:
                self.assertIn(t, PROJECTS, p["id"])
                self.assertNotEqual(t, SPECS[p["cohort"]]["project"], p["id"])

    def test_status_predictions_need_a_status_model(self):
        scope = pd.read_csv(REPO / "results" / "data" / "task_scope.csv").set_index("project")
        for p in FOLLOWUP["cross_task"]["status_probability_by_subtype"]["predictions"]:
            self.assertTrue(scope.at[SPECS[p["cohort"]]["project"], "tumor_status_in_scope"], p["id"])

    def test_event_predictions_name_an_mxe_event(self):
        for p in FOLLOWUP["subtype_regulators"]["predictions"]:
            if "event" in p:
                self.assertRegex(p["event"], r"^ENSG\d+\.\d+#chr\w+:\d+-\d+#[+-];")


class TestSiteProbeHelpers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = load_script("38_subtype_site_probe.py")

    def test_cramers_v(self):
        a = pd.Series(list("AABBCC") * 10)
        self.assertAlmostEqual(self.m.cramers_v(a, a), 1.0)
        b = pd.Series(list("XY") * 30)
        self.assertLess(self.m.cramers_v(pd.Series(list("AAB") * 20), b), 0.05)

    def test_class_absent_from_training_gets_zero(self):
        rng = np.random.default_rng(0)
        n = 60
        y = np.array([0] * 25 + [1] * 25 + [2] * 10)
        X = rng.normal(size=(n, 30)) + y[:, None]
        groups = np.array([f"g{i % 10}" for i in range(50)] + ["only2"] * 10)   # class 2 in one group
        p = self.m.oof_proba(X, y, groups, seed=0, k=3)
        self.assertTrue(np.allclose(p.sum(1), 1.0))
        self.assertTrue(np.all(p[y == 2, 2] == 0))       # its fold's training never saw class 2
        self.assertTrue(np.all(p[y != 2].sum(1) > 0))



class TestMsiRepeatHelpers(unittest.TestCase):
    """scripts/42 (D-52..D-54): window orientation, run lengths and the background-weighted score."""
    @classmethod
    def setUpClass(cls):
        cls.m = load_script("42_msi_repeats.py")

    def test_predictions_are_named_once(self):
        ids = [p["id"] for k in ("msi_cis", "msi_long_tracts", "msi_a3ss_replication")
               for p in FOLLOWUP[k]["predictions"]]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(set(FOLLOWUP["msi_cis"]["cohorts"]), {"UCEC_TCGA", "STAD_TCGA", "COAD_TCGA"})

    def test_max_runs_counts_a_run_overlapping_the_window_in_full(self):
        seq = "AG" + "C" * 46 + "T" * 12 + "GA"          # positions 1..62; the T-run spans 49..60
        r = self.m.max_runs(seq, 1, 50)
        self.assertEqual(r["T"], 12)
        self.assertEqual(r["C"], 46)
        self.assertEqual(self.m.max_runs(seq, 1, 40)["T"], 0)

    def test_window_reads_the_upstream_intron_on_both_strands(self):
        intron = "T" * 20 + "CAG"                         # ...TTTTCAG|exon on the sense strand
        exon, down = "GGGGG", "GTAAGT" + "A" * 60
        plus = "A" * 60 + intron + exon + down
        es = 60 + len(intron)
        r = pd.Series({"strand": "+", "es": es, "ee": es + 5, "left_end": 0, "right_start": len(plus)})
        s, n = self.m.window_seq(plus, r, "acc")
        self.assertEqual(s[:3], "GAC")                     # -1 G, -2 A, -3 C
        self.assertEqual(self.m.max_runs(s, 1, 50)["T"], 20)
        comp = str.maketrans("ACGT", "TGCA")
        minus = plus.translate(comp)[::-1]                 # the same gene on the minus strand
        L = len(minus)
        rm = pd.Series({"strand": "-", "es": L - (es + 5), "ee": L - es, "left_end": 0, "right_start": L})
        sm, _ = self.m.window_seq(minus, rm, "acc")
        self.assertEqual(sm[:len(s)], s)

    def test_score_weights_background_to_the_repeat_deciles(self):
        P = np.array([[0.9, 0.1, 0.5, 0.5], [0.1, 0.9, 0.5, 0.5]])
        in_r = np.array([True, False, False, False])
        in_b = np.array([False, True, True, True])
        dec = np.array([0, 0, 1, 1])                       # only exon 1 shares the repeat's decile
        sc = self.m.set_score(P, in_r, in_b, dec)
        np.testing.assert_allclose(sc, [0.8, -0.8])

    def test_auroc_lower(self):
        a, p = self.m.auroc_lower(np.array([1.0, 2.0]), np.array([3.0, 4.0, 5.0]))
        self.assertEqual(a, 1.0)


class TestRepeatDnaReads(unittest.TestCase):
    """scripts/45 (D-56): the BAM reader decodes records, and a read's run is called only between
    its exact flanks."""
    @classmethod
    def setUpClass(cls):
        cls.m = load_script("45_repeat_dna_reads.py")

    @staticmethod
    def bam(path: Path, reads) -> None:
        import gzip
        import struct
        code = {b: i for i, b in enumerate("=ACMGRSVTWYHKDBN")}
        refs = [b"chr1\0", b"chr2\0"]
        out = b"BAM\1" + struct.pack("<i", 0) + struct.pack("<i", len(refs))
        for r in refs:
            out += struct.pack("<i", len(r)) + r + struct.pack("<i", 1000)
        for ref, pos, mapq, flag, seq in reads:
            name = b"r\0"
            packed = bytes((code[seq[i]] << 4) | (code[seq[i + 1]] if i + 1 < len(seq) else 0)
                           for i in range(0, len(seq), 2))
            cigar = struct.pack("<I", (len(seq) << 4) | 0)
            body = (struct.pack("<iiBBHHHiiii", ref, pos, len(name), mapq, 0, 1, flag, len(seq), -1, -1, 0)
                    + name + cigar + packed + bytes([30] * len(seq)))
            out += struct.pack("<i", len(body)) + body
        with gzip.open(path, "wb") as fh:                  # two members, as BGZF writes them
            fh.write(out[:40])
        with gzip.open(path, "ab") as fh:
            fh.write(out[40:])

    def test_reader_decodes_odd_length_reads_and_flags(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.bam"
            self.bam(p, [(0, 5, 60, 0, "ACGTN"), (1, 7, 10, 0x400, "TTTTTTG")])
            got = list(self.m.read_bam(p))
        self.assertEqual(got, [("chr1", 5, 60, 0, "ACGTN"), ("chr2", 7, 10, 0x400, "TTTTTTG")])

    def test_run_is_called_only_between_exact_flanks(self):
        left, right = "CAGGCATC", "GAATCCAG"
        f = self.m.run_length
        self.assertEqual(f("AC" + left + "T" * 13 + right + "GG", left, "T", right), 13)
        self.assertIsNone(f(left + "T" * 13 + "GAATCCAC", left, "T", right))          # flank mismatch
        self.assertIsNone(f(left + "T" * 6 + "C" + "T" * 6 + right, left, "T", right))  # interrupted run
        self.assertIsNone(f(left + right, left, "T", right))                            # no run base
        self.assertIsNone(f(left + "T" * 13 + right[:5], left, "T", right))             # read ends in it
        self.assertEqual(f(left + "A" + left + "T" * 9 + right, left, "T", right), 9)  # second copy

    def test_genotype_filters_reads_and_contraction_is_stutter_corrected(self):
        import tempfile
        left, right = "CAGGCATC", "GAATCCAG"
        L = pd.DataFrame({"chrom": ["chr1"], "start": [108], "end": [120], "left": [left], "base": ["T"],
                          "right": [right]})
        read = "GG" * 50 + left + "T" * 12 + right
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.bam"
            self.bam(p, [(0, 0, 60, 0, read), (0, 0, 5, 0, read), (0, 0, 60, 0x400, read),
                         (0, 0, 60, 0, read.replace("T" * 12, "T" * 11))])
            self.assertEqual(sorted(self.m.genotype_bam(p, L)[0]), [11, 12])
        g = pd.DataFrame({"file_id": ["t", "n"], "case_id": ["a", "a"], "cohort": ["X", "X"],
                          "group": ["MSI", "MSI"], "role": ["tumour", "normal"], "locus": [0, 0],
                          "lengths": [[12] * 5 + [11] * 5, [12] * 9 + [11]]})
        ct = self.m.contraction_table(g)
        self.assertEqual(ct.germline_mode.iloc[0], 12)
        self.assertAlmostEqual(ct.c.iloc[0], 0.5 - 0.1)
        self.assertEqual(ct.reference.iloc[0], "own normal")


class TestBrcaGermline(unittest.TestCase):
    """scripts/46 (D-57, A-13, A-14): the exact Fisher test, the three-group split and pair agreement."""
    @classmethod
    def setUpClass(cls):
        cls.m = load_script("46_brca_germline.py")

    def test_fisher_matches_scipy(self):
        from scipy.stats import fisher_exact
        rng = np.random.default_rng(0)
        for a, b, c, d in rng.integers(0, 40, size=(200, 4)):
            self.assertAlmostEqual(self.m.fisher_p(a, b, c, d), fisher_exact([[a, b], [c, d]])[1], places=9)

    def test_three_groups_finds_the_gaps(self):
        v = np.sort(np.r_[np.full(30, 0.05), np.full(30, 0.5), np.full(30, 0.95)])
        a, b, mu, sd = self.m.three_groups(v)
        self.assertTrue(0.05 < a <= 0.5 < b <= 0.95)
        np.testing.assert_allclose(mu, [0.05, 0.5, 0.95])
        self.assertAlmostEqual(sd, 0.0)

    def test_pair_kappa(self):
        g = np.repeat([0.0, 1.0, 2.0], 10)
        self.assertAlmostEqual(self.m.pair_kappa(g, g), 1.0)
        self.assertLess(self.m.pair_kappa(g, np.roll(g, 10)), 0)       # every pair moves group
        self.assertTrue(np.isnan(self.m.pair_kappa(g[:5], g[:5])))     # fewer than 20 pairs
        h = g.copy()
        h[0] = np.nan                                                   # incomplete pairs drop out
        self.assertAlmostEqual(self.m.pair_kappa(h, g), 1.0)


class TestMsiFollowup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = load_script("50_msi_followup.py")

    def test_run_and_branch_point_from_sequence(self):
        # position order (index 0 = position 1): AG reversed, a T12, then yUnAy read 3' to 5'
        rev = "GACG" + "T" * 12 + "C" + "A" + "G" + "T" + "C" + "G" * 30
        g = self.m.run_geometry(rev)
        self.assertEqual((g["run_near"], g["run_far"], g["run_len"]), (5, 16, 12))
        self.assertEqual(g["d_bp"], 18)
        no_t = rev[:19] + "G" + rev[20:]                 # the T 2 nt 5' of the A is lost
        self.assertTrue(np.isnan(self.m.run_geometry(no_t)["d_bp"]))
        self.assertEqual(self.m.run_geometry("GACG" + "C" * 60), {})

    def test_context_run_counts_the_right_flank(self):
        self.assertEqual(self.m.context_run("GCAGTC" + "T" + "TTTG", "T"), 4)
        self.assertEqual(self.m.context_run("GCAGTC" + "AA" + "AAAAA", "AA"), 7)
        self.assertEqual(self.m.context_run("GCAGTC" + "A" + "CAAAA", "A"), 1)

    def test_adjusted_rho_permutes_within_bins(self):
        rng = np.random.default_rng(1)
        x = rng.normal(size=300)
        bins = np.repeat(np.array(["a", "b", "c"], dtype=object), 100)
        y = x + 3 * (bins == "b") + rng.normal(scale=0.5, size=300)
        r, p, n = self.m.adjusted_rho(y, x, bins, [], np.random.default_rng(0), perm=99)
        self.assertGreater(r, 0.7)
        self.assertAlmostEqual(p, 0.01)
        self.assertEqual(n, 300)
        r0, p0, _ = self.m.adjusted_rho(rng.normal(size=300), x, bins, [], np.random.default_rng(0), perm=99)
        self.assertLess(abs(r0), 0.2)
        self.assertGreater(p0, 0.01)


class TestFigureLegends(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = load_script("_figbase.py")

    def test_question_caption_becomes_a_legend(self):
        cap = "**Why?**\n\nMethod text.\n\nLeft: errors per cancer. Right: pooled errors.\n\n- T1, a test."
        out = self.m.legend("S1", "A finding.", cap)
        self.assertTrue(out.startswith("**Figure S1. A finding.** Method text."))
        self.assertIn("(a) Errors per cancer. (b) Pooled errors.", out)
        self.assertIn("- T1, a test.", out)

    def test_letter_markers_become_bracketed(self):
        out = self.m.legend("X1", "T.", "**Q?**\n\nIntro. a, each dot is one view. b, the five cohorts.")
        self.assertIn("(a) Each dot is one view. (b) The five cohorts.", out)

    def test_tidy_keeps_an_inner_list(self):
        out = self.m.tidy("Columns:\n- RI: retention.\n- SE: inclusion,\n  continued.\n\nNext paragraph.")
        self.assertEqual(out, "Columns:\n- RI: retention.\n- SE: inclusion, continued.\n\nNext paragraph.")


class TestMsiCode(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = load_script("51_msi_code.py")

    def test_arms_follow_the_rule(self):
        acc = np.array([10, 19, 20, 25, 30, 5, 3, 21])
        bp = np.array([40, 20, 36, 35, np.nan, 50, 30, 80])
        down, up, b = self.m.arms(acc, bp)
        self.assertEqual(down.tolist(), [True, True, False, False, False, False, False, False])
        self.assertEqual(up.tolist(), [False, False, True, False, False, False, False, True])
        self.assertEqual(b.tolist(), [False, False, False, False, False, True, True, False])

    def test_auc_lower_matches_mann_whitney(self):
        from scipy.stats import mannwhitneyu
        rng = np.random.default_rng(3)
        s = rng.normal(size=200)
        pos = rng.random(200) < 0.2
        s[pos] -= 0.5
        u = mannwhitneyu(s[~pos], s[pos]).statistic / (pos.sum() * (~pos).sum())
        self.assertAlmostEqual(self.m.auc_lower(s, pos), u)

    def test_resid_rho_removes_the_bins(self):
        rng = np.random.default_rng(4)
        bins = np.repeat(np.array(["a", "b"], dtype=object), 200)
        x = rng.normal(size=400) + 5 * (bins == "b")
        y = rng.normal(size=400) + 5 * (bins == "b")
        self.assertLess(abs(self.m.resid_rho(y, x, bins, [])), 0.15)
        raw = np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1]
        self.assertGreater(raw, 0.5)                       # the bins alone make x and y agree

    def test_msi_calls_use_the_prespecified_cuts(self):
        import pandas as pd
        x = pd.Series([0.3, 0.4, 0.5, 0.6, 0.61, np.nan])
        self.assertEqual(self.m.msi_call(x, "mantis").tolist(),
                         ["MSS", "between", "between", "between", "MSI-H", "no score"])
        y = pd.Series([3.9, 4.0, 10.0, 10.1])
        self.assertEqual(self.m.msi_call(y, "msisensor").tolist(), ["MSS", "between", "between", "MSI-H"])


class TestMsiReversal(unittest.TestCase):
    """D-70's label-blind annotation: the AG exclusion zone and the shield class."""

    @classmethod
    def setUpClass(cls):
        cls.m = load_script("52_msi_reversal.py")

    @staticmethod
    def rev(sense: str) -> str:
        """Intron sense sequence ending at the exon -> position order (index 0 = position 1)."""
        return sense[::-1]

    def test_agez_ignores_the_3ss_and_near_ags(self):
        # positions from the exon: the 3'ss AG has its G at 1; an AG with its G at 6 is ignored (< 12);
        # the next AG upstream has its G at 15
        sense = "CCCC" + "AG" + "TTTTTTT" + "AG" + "TTT" + "AG"
        z, cens = self.m.agez(self.rev(sense))
        self.assertEqual((z, cens), (15, False))

    def test_agez_is_censored_without_an_upstream_ag(self):
        sense = "C" * 40 + "AG"
        self.assertEqual(self.m.agez(self.rev(sense)), (42, True))

    def test_shield_needs_a_branch_point_between_the_run_and_the_ag(self):
        # positions 19..1 are C17 + AG; the yUnAy motif CTCAC sits at 24..20 with its A at 21;
        # a T22 run lies at 50..29, upstream of that branch point
        tail = "C" * 17 + "AG"
        bp = "CTCAC"
        sense = "G" * 30 + "T" * 22 + "GGGG" + bp + tail
        r = self.m.shield(self.rev(sense))
        self.assertEqual(r["sh_len"], 22)
        self.assertEqual(r["sh_bp"], 21)
        # without the branch point there is no shield
        r2 = self.m.shield(self.rev("G" * 30 + "T" * 22 + "GGGG" + "GGGGG" + tail))
        self.assertTrue(np.isnan(r2["sh_start"]))

    def test_shield_ignores_runs_at_the_acceptor(self):
        sense = "G" * 40 + "CTCAC" + "T" * 22 + "AG"  # run at positions 1-24: an acceptor run, not a shield
        self.assertTrue(np.isnan(self.m.shield(self.rev(sense))["sh_start"]))


class TestMsiModels(unittest.TestCase):
    """D-71's helpers: robust z, sequence bins, precision at recall, AUROC orientation."""

    @classmethod
    def setUpClass(cls):
        cls.m = load_script("54_msi_models.py")

    def test_robust_z_ignores_a_minority_of_outliers(self):
        v = np.r_[np.arange(10.0), 1000.0]
        z = self.m.robust_z(v)
        self.assertAlmostEqual(float(np.median(z)), 0.0)
        self.assertLess(abs(z[:10]).max(), 2.0)
        self.assertTrue(np.all(self.m.robust_z(np.full(5, 3.0)) == 0))
        self.assertTrue(np.all(self.m.robust_z(np.full(4, np.nan)) == 0))

    def test_bins_partition_the_repeat_exons(self):
        at = np.array([10, 13, 15, 19, 20, 25, 30, 21, 22, 40, 5])
        db = np.array([30, 30, 30, 30, 25, 33, 45, 60, 85, np.nan, 40])
        m = np.vstack(self.m.bin_masks(at, db))
        self.assertTrue(np.all(m[:, :10].sum(axis=0) == 1))   # every repeat exon in exactly one bin
        self.assertEqual(m[:, 10].sum(), 0)                   # a run of 5 is in none
        self.assertTrue(m[9, 9])                              # run >= 20 without a branch point

    def test_precision_at_recall_counts_tied_negatives_first(self):
        y = np.array([True, False, True, False, False, True])
        s = np.array([0.9, 0.9, 0.8, 0.1, 0.05, 0.0])
        # 2 of 3 positives needed; tied negative at 0.9 ranks first, so the 2nd positive is at rank 3
        self.assertAlmostEqual(self.m.precision_at_recall(s, y, 0.5), 2 / 3)
        self.assertAlmostEqual(self.m.auc(np.array([3.0, 2.0, 1.0]), np.array([True, False, False])), 1.0)


class TestMsiBenchmark(unittest.TestCase):
    """D-72's helpers: triangulated references, within-cancer AUROC, matched-volume operating point."""

    @classmethod
    def setUpClass(cls):
        cls.m = load_script("55_msi_benchmark.py")

    def test_each_reference_leaves_out_the_judged_tool(self):
        d = pd.DataFrame({"n_hp": [20, 20, 0, 0, 20], "mantis": [0.9, 0.1, 0.1, 0.9, 0.5],
                          "msisensor": [20.0, 20.0, 1.0, 1.0, 20.0]})
        r = self.m.references(d)
        # tumour 1: MANTIS disagrees with the others; it is positive in R_MANTIS, left out of R_MC3.
        # tumour 3: only MANTIS calls it; R_MANTIS (MC3 and MSIsensor) holds it as a negative
        self.assertEqual(list(r["MANTIS"]), [1, 1, 0, 0, 1])
        self.assertEqual(list(r["MSIsensor"]), [1, -1, 0, -1, -1])   # MANTIS 0.5 is between
        self.assertEqual(list(r["MC3"]), [1, -1, 0, -1, -1])
        self.assertEqual(list(r["all"]), [1, -1, 0, -1, -1])

    def test_within_cancer_auroc_ignores_pairs_across_cancers(self):
        # perfect within each cancer, but cancer B's scores all sit below cancer A's
        s = np.array([10.0, 9.0, 2.0, 1.0])
        y = np.array([True, False, True, False])
        g = np.array(["A", "A", "B", "B"])
        self.assertAlmostEqual(self.m.auc_within(s, y, g), 1.0)
        self.assertAlmostEqual(self.m.auc_pooled(s, y), 0.75)
        # a cancer without positives contributes no pairs
        self.assertAlmostEqual(self.m.auc_within(np.append(s, 100.0), np.append(y, False), np.append(g, "C")), 1.0)

    def test_matched_volume_counts_tied_negatives_first(self):
        s = np.array([5.0, 5.0, 1.0, 0.0])
        y = np.array([True, False, True, False])
        self.assertEqual(self.m.matched_volume(s, y, 1), (0.0, 0.0))
        self.assertEqual(self.m.matched_volume(s, y, 2), (0.5, 0.5))

    def test_bootstrap_matches_auc_within_on_each_resample(self):
        rng = np.random.default_rng(1)
        g = np.repeat(["A", "B", "C"], 30)
        y = rng.random(90) < 0.3
        y[60:] = False                                              # cancer C has no positive
        S = {"x": rng.normal(size=90) + y, "z": rng.normal(size=90)}
        fast = self.m.paired_bootstrap(S, y, g, np.random.default_rng(0))
        r2 = np.random.default_rng(0)
        strata = [np.flatnonzero((g == a) & (y == b)) for a in np.unique(g) for b in (False, True)]
        strata = [s for s in strata if len(s)]
        for k in range(5):
            i = np.concatenate([r2.choice(s, len(s)) for s in strata])
            for m in S:
                self.assertAlmostEqual(fast[m][k], self.m.auc_within(S[m][i], y[i], g[i]))


class TestMsiRnaOnly(unittest.TestCase):
    """D-73's helpers: within-cancer percentiles and the within-cancer permutation AUROC."""

    @classmethod
    def setUpClass(cls):
        cls.m = load_script("56_msi_rna_only.py")

    def test_within_pct_ranks_inside_each_cancer(self):
        v = pd.Series([1.0, 2.0, 100.0, 200.0, np.nan])
        g = pd.Series(["A", "A", "B", "B", "B"])
        self.assertEqual(list(self.m.within_pct(v, g).round(2).fillna(-1)), [0.5, 1.0, 0.5, 1.0, -1])

    def test_auroc_perm_orientation_and_minimum(self):
        rng = np.random.default_rng(0)
        x = np.r_[np.arange(10.0), np.arange(10.0) + 100]
        case = np.r_[np.zeros(10, bool), np.ones(10, bool)]
        grp = np.r_[["A"] * 5, ["B"] * 5, ["A"] * 5, ["B"] * 5]
        a, p, n1, n0 = self.m.auroc_perm(x, case, grp, rng, lower=False)
        self.assertEqual((a, n1, n0), (1.0, 10, 10))
        self.assertLess(p, 0.01)
        self.assertEqual(self.m.auroc_perm(x, case, grp, rng, lower=True)[0], 0.0)
        few = case.copy()
        few[13:] = False                                   # 3 cases: below MIN_N, not testable
        self.assertTrue(np.isnan(self.m.auroc_perm(x, few, grp, rng, lower=False)[0]))


class TestMsiVsExpression(unittest.TestCase):
    """D-74's fitting helpers on synthetic data (no real labels)."""

    @classmethod
    def setUpClass(cls):
        cls.m = load_script("57_msi_vs_expression.py")

    def _data(self):
        rng = np.random.default_rng(1)
        n, p = 240, 30
        groups = np.repeat(["A", "B", "C", "D"], n // 4)
        y = rng.random(n) < 0.25
        X = rng.normal(size=(n, p)).astype(np.float32)
        X[:, 3] += 2.0 * y                                  # one informative column
        return X, y, groups

    def test_fit_predict_keeps_always_columns_and_top_k(self):
        X, y, _ = self._data()
        tr, te = np.arange(180), np.arange(180, 240)
        order = np.argsort(-np.arange(X.shape[1]))         # columns 29, 28, ...
        s, _, cols = self.m.fit_predict(X, tr, te, y, {"K": 2, "C": 1.0}, order, always=np.array([0, 1]))
        self.assertEqual(list(cols), [0, 1, 28, 29])
        self.assertEqual(len(s), 60)

    def test_run_tunes_on_training_groups_and_finds_the_signal(self):
        X, y, groups = self._data()
        tr, te = np.flatnonzero(groups != "D"), np.flatnonzero(groups == "D")
        old = self.m.GRIDS.get("T")
        self.m.GRIDS["T"] = [{"K": k, "C": c} for k in [1, None] for c in [0.01, 1.0]]
        try:
            r = self.m.run("T", X, y, groups, tr, te, offset=2, always=np.array([0, 1]))
        finally:
            self.m.GRIDS.pop("T") if old is None else self.m.GRIDS.update(T=old)
        self.assertEqual(len(r["score"]), len(te))
        self.assertEqual(len(r["table"]), 4)
        self.assertEqual(sum(row["chosen"] for row in r["table"]), 1)
        self.assertGreater(self.m.m54.auc(r["score"], y[te]), 0.8)
        self.assertIn(3, r["cols"])                       # the informative column is selected


class TestMsiMagnitude(unittest.TestCase):
    """D-75's statistics on synthetic data (no real labels)."""

    @classmethod
    def setUpClass(cls):
        cls.m = load_script("60_msi_magnitude.py")

    def test_deviations_remove_a_tumour_wide_level(self):
        rng = np.random.default_rng(0)
        X = rng.uniform(0.3, 0.7, size=(40, 20))
        prof = np.arange(40) >= 10
        X[~prof] += 0.05                                    # tumour-wide shift in the first 10
        X[~prof, 0] -= 0.20                                 # one exon moves further
        dec = np.repeat([0, 1], 10)
        is_b = np.zeros(20, bool)
        is_b[[2, 3, 4, 12, 13, 14]] = True
        dev, adj = self.m.deviations(X, prof, dec, is_b)
        self.assertAlmostEqual(float(np.median(dev[prof, 5])), 0.0, 9)
        self.assertLess(abs(np.mean(adj[~prof][:, 5:10])), 0.05)   # the common shift is removed
        self.assertLess(np.mean(adj[~prof, 0]), -0.15)               # the exon's own move is kept

    def test_tail_excess_is_zero_under_a_pure_location_shift(self):
        rng = np.random.default_rng(1)
        P = rng.normal(0.5, 0.05, size=(4000, 3))
        M = rng.normal(0.5, 0.05, size=(4000, 3))
        M[:, 1] -= 0.04                                     # everyone moves a little
        M[:1000, 2] -= 0.15                                 # a quarter move a lot, the rest not at all
        shift = np.median(M, 0) - np.median(P, 0)
        obs, exp = self.m.tail_excess(P, M, shift, lower=True)
        self.assertLess(abs(obs[0] - exp[0]), 0.02)
        self.assertLess(abs(obs[1] - exp[1]), 0.03)
        self.assertGreater(obs[2] - exp[2], 0.15)
        self.assertAlmostEqual(obs[0], 0.05, delta=0.015)

    def test_tail_excess_handles_missing_values(self):
        P = np.array([[0.1, np.nan], [0.2, np.nan], [0.3, np.nan], [0.4, np.nan]])
        M = np.array([[0.0, 0.5], [np.nan, 0.5]])
        obs, exp = self.m.tail_excess(P, M, np.zeros(2), lower=True)
        self.assertEqual(obs[0], 1.0)                       # 1 of 1 finite value below the 5th percentile
        self.assertTrue(np.isnan(exp[1]))

    def test_midrank_puts_ties_at_the_middle(self):
        P = np.array([[1.0, 0.2, np.nan], [1.0, 0.4, np.nan], [1.0, 0.6, np.nan], [1.0, 0.8, np.nan]])
        u = self.m.midrank(P, np.array([1.0, 0.5, 0.3]))
        self.assertEqual(u[0], 0.5)                         # tied with every proficient value
        self.assertEqual(u[1], 0.5)
        self.assertTrue(np.isnan(u[2]))

    def test_strat_perm_matches_strata_and_detects_a_shift(self):
        rng = np.random.default_rng(2)
        strata = np.repeat(["s0", "s1"], 200)
        x = np.where(strata == "s1", 1.0, 0.0) + rng.normal(0, 0.1, 400)
        is_a = np.zeros(400, bool)
        is_a[150:200] = True                                # 50 of s0 and 150 of s1 are class a, so the
        is_a[200:350] = True                                # raw difference reflects only the strata
        obs, matched, p = self.m.strat_perm(x, is_a, strata, np.random.default_rng(0), 2000)
        self.assertGreater(obs, 0.3)                        # raw difference is confounded by stratum
        self.assertLess(abs(matched), 0.03)
        self.assertGreater(p, 0.05)
        x2 = x + 0.1 * is_a
        _, matched2, p2 = self.m.strat_perm(x2, is_a, strata, np.random.default_rng(0), 2000)
        self.assertAlmostEqual(matched2, 0.1, delta=0.03)
        self.assertLess(p2, 0.001)


class TestMsiRisingConsequence(unittest.TestCase):
    """D-76's annotation helpers on synthetic transcripts."""

    @classmethod
    def setUpClass(cls):
        cls.m = load_script("61_msi_rising_consequence.py")

    def test_parse_feature(self):
        f = self.m.parse_feature("ENSG00000020922.13#chr11:94479673-94479761#-;94478734-94478876;94485923-94486084")
        self.assertEqual((f["gene"], f["chrom"], f["strand"], f["es"], f["ee"]), ("ENSG00000020922", "chr11", "-", 94479673, 94479761))
        self.assertEqual((f["le"], f["rs"]), (94478876, 94485923))

    def test_place_canonical_and_inserted_on_both_strands(self):
        ex = [(100, 200), (300, 400), (600, 700), (900, 1000)]
        plus = {"chrom": "chr1", "strand": "+", "exons": ex, "cds": [(150, 200), (300, 400), (600, 700), (900, 950)]}
        f = {"chrom": "chr1", "strand": "+", "es": 299, "ee": 400, "le": 200, "rs": 599}
        self.assertEqual(self.m.place(f, plus), ("canonical", 1))
        g = {"chrom": "chr1", "strand": "+", "es": 449, "ee": 500, "le": 400, "rs": 599}   # joins exons 2 and 3
        self.assertEqual(self.m.place(g, plus), ("inserted", 1))
        utr = dict(plus, cds=[None, (320, 400), (600, 700), (900, 950)])
        h = {"chrom": "chr1", "strand": "+", "es": 249, "ee": 260, "le": 200, "rs": 299}   # junction before the CDS
        self.assertEqual(self.m.place(h, utr)[0], "utr")
        minus = {"chrom": "chr1", "strand": "-", "exons": ex[::-1], "cds": [(900, 950), (600, 700), (300, 400), (150, 200)]}
        k = {"chrom": "chr1", "strand": "-", "es": 449, "ee": 500, "le": 400, "rs": 599}   # between 600-700 and 300-400
        self.assertEqual(self.m.place(k, minus), ("inserted", 1))
        self.assertEqual(self.m.place(dict(k, strand="+"), minus)[0], "other")

    def test_consequence_classes(self):
        down = ["GCC" * 40, "GCC" * 10]                       # 120 nt, then the last exon
        r = self.m.consequence("", "ATGAAA", down)
        self.assertEqual((r["consequence"], r["peptide"]), ("in_frame", 2))
        r = self.m.consequence("", "TAAGGG", down)            # in-frame stop at 0; last junction at 6 + 120
        self.assertEqual(r["consequence"], "ptc_nmd")
        self.assertEqual(r["dist"], 126)
        r = self.m.consequence("G", "GA", ["TAA" + "C" * 30])  # codon GGA, then TAA in the last exon
        self.assertEqual(r["consequence"], "ptc_escape")
        r = self.m.consequence("", "GG", ["C" * 91, "TAGC"])   # frameshift: TAG at 93, in the last exon
        self.assertEqual(r["consequence"], "ptc_escape")
        r = self.m.consequence("", "GG", ["CCCCTGA" + "C" * 90, "AAA"])   # frame GGC CCC CTG A.. -> no stop at 3k
        self.assertIn(r["consequence"], {"ptc_nmd", "ptc_escape", "no_stop"})
        r = self.m.consequence("", "G", ["GTAA" + "C" * 80, "AAAA"])    # GGT AAC ... -> no TAA in frame
        self.assertEqual(r["consequence"], "no_stop")
        r = self.m.consequence("", "G", ["GCTAA" + "C" * 80, "AAAA"])   # GGC TAA: stop at 3, junction at 86
        self.assertEqual((r["consequence"], r["dist"]), ("ptc_nmd", 83.0))


if __name__ == "__main__":
    unittest.main()
