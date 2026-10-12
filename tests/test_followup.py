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
