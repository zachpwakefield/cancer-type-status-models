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


if __name__ == "__main__":
    unittest.main()
