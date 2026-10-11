"""Unit tests on synthetic data: leakage guards, encodings, preprocessing, model, SHAP.

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

from ctsm.explain import additivity_error, deep_shap, importance_table  # noqa: E402
from ctsm.extract import make_labels, prefilter_pass, structural_zeros  # noqa: E402
from ctsm.matrix_io import ViewMatrix, decode, encode_int16  # noqa: E402
from ctsm.model import balanced_log_loss, class_weights, train_model  # noqa: E402
from ctsm.preprocess import column_stats, fit_preprocessor  # noqa: E402
from ctsm.splits import fold_indices, make_splits, min_inner_train  # noqa: E402

PP = {"min_observed_frac": 0.5, "gex_min_log2cpm": 1.0, "gex_min_expressed_frac": 0.1, "clip_z": 10.0}


def toy_table(n_cases=120, seed=0) -> pd.DataFrame:
    """Cases with a tumour; one in four also has a matched normal (same case_id)."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_cases):
        rows.append({"file_id": f"t{i}", "case_id": f"c{i}", "y": 1})
        if i % 4 == 0:
            rows.append({"file_id": f"n{i}", "case_id": f"c{i}", "y": 0})
    t = pd.DataFrame(rows)
    return t.iloc[rng.permutation(len(t))].reset_index(drop=True)


class TestSplits(unittest.TestCase):
    def test_grouped_nested_and_complete(self):
        t = toy_table()
        s = make_splits(t, "tumor_status", "TOY", 5, 3, 7)
        self.assertTrue((s.groupby("case_id")["outer"].nunique() == 1).all())
        for k in range(5):
            tr, te = fold_indices(s, k)
            self.assertEqual(len(set(s.case_id.iloc[tr]) & set(s.case_id.iloc[te])), 0)
            self.assertEqual(set(tr) | set(te), set(range(len(s))))
            for j in range(3):
                itr, iva = fold_indices(s, k, j)
                self.assertTrue(set(itr) | set(iva) == set(tr))
                self.assertEqual(len(set(s.case_id.iloc[itr]) & set(s.case_id.iloc[iva])), 0)
                self.assertEqual(len(set(iva) & set(te)), 0)
        self.assertEqual(sorted(s["outer"].unique()), list(range(5)))
        self.assertGreater(min_inner_train(s, 5, 3), 0)

    def test_deterministic(self):
        t = toy_table()
        a = make_splits(t, "tumor_status", "TOY", 5, 3, 7)
        b = make_splits(t, "tumor_status", "TOY", 5, 3, 7)
        pd.testing.assert_frame_equal(a, b)


class TestEncoding(unittest.TestCase):
    def test_int16_roundtrip_and_missing(self):
        x = np.array([[0.0, 0.123, np.nan], [1.0, -0.5, 0.999]])
        enc = encode_int16(x)
        dec = decode(enc, "int16_x1000")
        np.testing.assert_allclose(dec[np.isfinite(x)], x[np.isfinite(x)], atol=1e-6)
        self.assertTrue(np.isnan(dec[0, 2]))

    def test_int16_rejects_lossy(self):
        with self.assertRaises(ValueError):
            encode_int16(np.array([[0.1234]]))

    def test_viewmatrix_get(self):
        x = np.arange(20, dtype=np.float32).reshape(4, 5) / 1000
        M = ViewMatrix(encode_int16(x), "int16_x1000")
        np.testing.assert_allclose(M.get(np.array([2, 0]), np.array([4, 1])), x[np.ix_([2, 0], [4, 1])], atol=1e-6)
        np.testing.assert_allclose(M.get(np.array([1]), slice(1, 3)), x[[1], 1:3], atol=1e-6)


class TestStructuralZeros(unittest.TestCase):
    def test_rule(self):
        genes = np.array(["g1", "g1", "g1", "g2", "g2"], dtype=object)
        v = np.array([[0.6, 1.0], [0.4, np.nan], [np.nan, np.nan], [np.nan, 0.3], [np.nan, np.nan]])
        out, fill = structural_zeros(v, genes, 1e-6)
        self.assertEqual(out[2, 0], 0.0)          # g1, sample 0: observed sum to 1 -> blank is 0
        self.assertEqual(out[1, 1], 0.0)          # g1, sample 1: 1.0 observed -> blanks are 0
        self.assertTrue(np.isnan(out[4, 1]))      # g2, sample 1: sums to 0.3 -> stays missing
        self.assertTrue(np.isnan(out[3, 0]))      # g2, sample 0: nothing observed
        self.assertEqual(int(fill.sum()), 3)


class TestPrefilterSuperset(unittest.TestCase):
    """D-07: any feature a fold-level rule keeps on a training subset is kept by the cohort pre-filter."""

    def test_superset(self):
        rng = np.random.default_rng(1)
        n, p = 200, 400
        X = rng.random((n, p)).astype(np.float32)
        X[rng.random((n, p)) < rng.random(p)] = np.nan          # feature-specific missingness
        X[:, :5] = 0.5                                          # constant features
        t = pd.DataFrame({"file_id": [f"s{i}" for i in range(n)], "case_id": [f"c{i}" for i in range(n)],
                          "y": rng.integers(0, 2, n)})
        s = make_splits(t, "tumor_status", "TOY", 5, 3, 3)
        k_min = int(np.ceil(PP["min_observed_frac"] * min_inner_train(s, 5, 3)))
        passed = prefilter_pass(X.T, {"TOY": np.arange(n)}, {"TOY": k_min})["TOY"]
        M = ViewMatrix(X, "float32")
        y = t["y"].to_numpy()
        for k in range(5):
            for j in [None, 0, 1, 2]:
                tr = fold_indices(s, k, j)[0]
                pre = fit_preprocessor(M, tr, y[tr], "se", "variance", p, PP, 2)
                self.assertTrue(passed[pre.cols].all(), f"fold {k}/{j} selected a pre-filtered feature")


class TestPreprocessorUsesTrainingRowsOnly(unittest.TestCase):
    def test_test_rows_cannot_change_fit(self):
        rng = np.random.default_rng(2)
        X = rng.normal(size=(100, 50)).astype(np.float32)
        y = rng.integers(0, 3, 100)
        tr = np.arange(70)
        a = fit_preprocessor(ViewMatrix(X, "float32"), tr, y[tr], "gex", "anova", 20, PP, 3)
        X2 = X.copy()
        X2[70:] = rng.normal(size=(30, 50)) * 100          # change only held-out rows
        b = fit_preprocessor(ViewMatrix(X2, "float32"), tr, y[tr], "gex", "anova", 20, PP, 3)
        np.testing.assert_array_equal(a.cols, b.cols)
        np.testing.assert_array_equal(a.mean, b.mean)

    def test_anova_finds_signal(self):
        rng = np.random.default_rng(3)
        y = np.repeat([0, 1, 2], 40)
        X = rng.normal(size=(120, 30)).astype(np.float32)
        X[:, 7] += y * 3.0
        st = column_stats(ViewMatrix(X, "float32"), np.arange(120), y, 3, "se", 1.0)
        self.assertEqual(int(np.argmax(st["anova_f"])), 7)

    def test_impute_and_scale(self):
        X = np.array([[1.0, 2.0], [np.nan, 4.0], [3.0, 6.0], [5.0, np.nan]], dtype=np.float32)
        pre = fit_preprocessor(ViewMatrix(X, "float32"), np.arange(4), np.array([0, 1, 0, 1]), "se", "variance",
                               2, {**PP, "min_observed_frac": 0.5}, 2)
        Z = pre.transform(ViewMatrix(X, "float32"), np.arange(4))
        self.assertFalse(np.isnan(Z).any())
        np.testing.assert_allclose(Z.mean(axis=0), 0, atol=1e-5)


class TestModelAndShap(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(4)
        self.X = rng.normal(size=(300, 12)).astype(np.float32)
        self.yb = (self.X[:, 0] > 0).astype(int)
        self.ym = (self.X[:, 0] > 0).astype(int) + 2 * (self.X[:, 1] > 0)
        self.p = dict(n_layers=2, width=32, dropout=0.1, input_dropout=0.0, lr=3e-3, weight_decay=1e-4, batch_size=32)

    def test_class_weights_balance(self):
        w = class_weights(np.array([0, 1, 1, 1]), 2)
        self.assertAlmostEqual(w[0] * 1, w[1] * 3)

    def test_balanced_log_loss_ignores_imbalance(self):
        logits = np.array([[2.0], [2.0], [2.0], [-2.0]])
        a = balanced_log_loss(logits, np.array([1, 1, 1, 0]), True)
        b = balanced_log_loss(logits[[0, 3]], np.array([1, 0]), True)
        self.assertAlmostEqual(a, b, places=6)

    def test_binary_learns_and_shap_is_additive(self):
        m, ep, v, _ = train_model(self.X[:200], self.yb[:200], self.p, 2, 0, X_val=self.X[200:],
                                  y_val=self.yb[200:], max_epochs=60, patience=10)
        self.assertLess(v, 0.4)
        sv, ev = deep_shap(m, self.X[:50], self.X[200:], self.yb[200:], 2)
        self.assertLess(additivity_error(m, self.X[200:], self.yb[200:], 2, sv, ev), 1e-3)
        self.assertEqual(int(np.argmax(np.abs(sv).mean(axis=0))), 0)

    def test_multiclass_centred_shap_is_additive(self):
        m, ep, v, _ = train_model(self.X[:200], self.ym[:200], self.p, 4, 0, X_val=self.X[200:],
                                  y_val=self.ym[200:], max_epochs=60, patience=10)
        sv, ev = deep_shap(m, self.X[:50], self.X[200:], self.ym[200:], 4)
        self.assertLess(additivity_error(m, self.X[200:], self.ym[200:], 4, sv, ev), 1e-3)
        self.assertEqual(set(np.argsort(-np.abs(sv).mean(axis=0))[:2]), {0, 1})

    def test_fixed_epochs_deterministic(self):
        a, *_ = train_model(self.X, self.yb, self.p, 2, 11, fixed_epochs=3)
        b, *_ = train_model(self.X, self.yb, self.p, 2, 11, fixed_epochs=3)
        for (k, va), (_, vb) in zip(a.state_dict().items(), b.state_dict().items()):
            self.assertTrue(torch.equal(va, vb), k)


class TestAggregation(unittest.TestCase):
    def test_absent_feature_counts_as_zero(self):
        f1 = {"feature_ids": np.array(["a", "b"]), "index": {"a": 0, "b": 1}, "y": np.zeros(2, int),
              "X": np.array([[1.0, 0.0], [-1.0, 1.0]], np.float32), "shap": np.array([[1.0, 0.5], [-1.0, 0.5]], np.float32)}
        f2 = {"feature_ids": np.array(["a"]), "index": {"a": 0}, "y": np.zeros(2, int),
              "X": np.array([[2.0], [-2.0]], np.float32), "shap": np.array([[2.0], [-2.0]], np.float32)}
        ann = pd.DataFrame({"feature_id": ["a", "b"], "label": ["A", "B"]})
        t = importance_table([f1, f2], ann).set_index("feature_id")
        self.assertAlmostEqual(t.loc["a", "mean_abs_shap"], 6 / 4)
        self.assertAlmostEqual(t.loc["b", "mean_abs_shap"], 1 / 4)       # 0 for fold 2's samples
        self.assertEqual(t.loc["b", "n_folds_selected"], 1)
        self.assertAlmostEqual(t.loc["a", "direction_spearman"], 1.0)

    def test_labels_unique(self):
        f = pd.DataFrame({"symbol": ["G", "G", "H"], "region": ["r", "r", "r"]})
        self.assertEqual(make_labels(f, "se").tolist(), ["G r", "G r#2", "H r"])


if __name__ == "__main__":
    unittest.main()
