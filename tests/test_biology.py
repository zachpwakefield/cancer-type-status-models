"""Unit tests for the interpretation helpers (ctsm.biology) on synthetic data.

Run: OMP_NUM_THREADS=1 python -m unittest discover -s tests
"""
from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import torch  # noqa: E402,F401  (load torch's OpenMP runtime before numpy's; see scripts/_common.py)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ctsm.biology import add_q, bh, gene_ranking, goseq_ora, matched_pairs, pwf, terminal_index  # noqa: E402


class TestMultipleTesting(unittest.TestCase):
    def test_bh_known_values(self):
        p = np.array([0.01, 0.04, 0.03, 0.20])
        q = bh(p)
        np.testing.assert_allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.20])

    def test_add_q_counts_untested_sets(self):
        df = pd.DataFrame({"collection": ["h", "h"], "p_wallenius": [0.01, 0.02]})
        q2 = add_q(df, {"h": 2})["q"].to_numpy()
        q10 = add_q(df, {"h": 10})["q"].to_numpy()
        self.assertTrue(np.all(q10 > q2))         # more tested sets -> larger q


class TestEnrichment(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.genes = [f"G{i}" for i in range(2000)]
        self.rng = rng

    def test_pwf_monotone(self):
        n = self.rng.integers(1, 200, 2000)
        sel = self.rng.random(2000) < n / 400          # longer genes picked more often
        w = pwf(n, sel)
        order = np.argsort(n)
        self.assertTrue(np.all(np.diff(w[order]) >= -1e-12))

    def test_unbiased_universe_matches_hypergeometric(self):
        uni = pd.DataFrame({"symbol": self.genes, "n_events": 1})
        sel = set(self.genes[:100])
        gs = {"planted": set(self.genes[:30]) | set(self.genes[1000:1020]), "null": set(self.genes[500:540])}
        r = goseq_ora(sel, uni, gs, 15, 500, 1, "t").set_index("term")
        self.assertAlmostEqual(r.loc["planted", "p_wallenius"], r.loc["planted", "p_hyper"], places=6)
        self.assertLess(r.loc["planted", "p_wallenius"], 1e-10)

    def test_event_count_bias_is_discounted(self):
        # selected genes and the set's genes both have many events: hypergeometric is fooled,
        # the weighted test is not
        n = np.ones(2000, int)
        n[:300] = 64
        uni = pd.DataFrame({"symbol": self.genes, "n_events": n})
        sel = set(self.rng.choice(self.genes[:300], 80, replace=False)) | set(self.genes[1500:1520])
        gs = {"many_events": set(self.genes[:60])}
        r = goseq_ora(sel, uni, gs, 15, 500, 1, "t").iloc[0]
        self.assertLess(r["p_hyper"], 1e-6)
        self.assertGreater(r["p_wallenius"], 0.01)


class TestHelpers(unittest.TestCase):
    def test_gene_ranking_uses_max_not_sum(self):
        imp = pd.DataFrame({"gene_id": ["a", "a", "a", "b"], "symbol": ["A", "A", "A", "B"],
                            "feature_id": ["1", "2", "3", "4"], "label": ["x"] * 4,
                            "mean_abs_shap": [0.1, 0.1, 0.1, 0.25]})
        g = gene_ranking(imp)
        self.assertEqual(g["symbol"].iloc[0], "B")

    def test_matched_pairs(self):
        tab = pd.DataFrame({"file_id": ["t1", "n1", "t2", "n3", "t3"], "case_id": ["c1", "c1", "c2", "c3", "c3"],
                            "role": ["tumour", "normal", "tumour", "normal", "tumour"]})
        ti, ni, cases = matched_pairs(tab)
        self.assertEqual(list(cases), ["c1", "c3"])
        self.assertTrue(np.all(tab["case_id"].to_numpy()[ti] == tab["case_id"].to_numpy()[ni]))
        self.assertTrue(np.all(tab["role"].to_numpy()[ti] == "tumour"))

    def test_terminal_index_orientation(self):
        # one gene, three terminal exons at increasing coordinates
        F = pd.DataFrame({"feature_id": ["e1", "e2", "e3"], "gene_id": ["g", "g", "g"],
                          "region": ["chr1:100-200", "chr1:300-400", "chr1:500-600"]})
        X = np.array([[1.0, 0.0, 0.0],      # only the lowest coordinate used
                      [0.0, 0.0, 1.0],      # only the highest
                      [0.5, 0.0, 0.5],
                      [0.1, 0.1, 0.0]])     # too little usage -> NaN
        plus = terminal_index(X, F, pd.Series({"g": "+"}), 0.5)[:, 0]
        minus = terminal_index(X, F, pd.Series({"g": "-"}), 0.5)[:, 0]
        np.testing.assert_allclose(plus[:3], [0.0, 1.0, 0.5])
        np.testing.assert_allclose(minus[:3], [1.0, 0.0, 0.5])   # upstream = highest coordinate on '-'
        self.assertTrue(np.isnan(plus[3]))


if __name__ == "__main__":
    unittest.main()
