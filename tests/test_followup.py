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


if __name__ == "__main__":
    unittest.main()


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
