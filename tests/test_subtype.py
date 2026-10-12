"""Subtype task (docs/SUBTYPE_METHODS.md): label builder, xlsx reader, folds, metrics, and the
guarantee that the two existing tasks are unchanged.

Synthetic tests always run; tests on the real cohort are skipped where its local files are absent.
Run: OMP_NUM_THREADS=1 python -m unittest discover -s tests
"""
from __future__ import annotations

import os
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import torch  # noqa: E402,F401  (load torch's OpenMP runtime before numpy's; see scripts/_common.py)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ctsm import subtype as st  # noqa: E402
from ctsm.config import load_config  # noqa: E402
from ctsm.matrix_io import ViewMatrix  # noqa: E402
from ctsm.splits import fold_indices, make_splits  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


# ---- synthetic label world ------------------------------------------------------------------

def toy_samples(n=90) -> pd.DataFrame:
    """One tumour per case of project TOY, plus normals for a third of them (never labelled)."""
    rows = []
    for i in range(n):
        base = {"case_id": f"TCGA-AA-{i:04d}", "project": "TOY", "tss": "AA", "vial": "A"}
        rows.append({**base, "file_id": f"t{i}", "sample_id": f"TCGA-AA-{i:04d}-01A", "sample_type": "Primary Tumor",
                     "role": "tumour"})
        if i % 3 == 0:
            rows.append({**base, "file_id": f"n{i}", "sample_id": f"TCGA-AA-{i:04d}-11A",
                         "sample_type": "Solid Tissue Normal", "role": "normal"})
    return pd.DataFrame(rows)


RAW_A = ["Alpha type", "Beta type", "Gamma type", "Other, specify", "Alpha type|Beta type", "Alpha type|Alpha type"]


def toy_spec(concordance=True) -> dict:
    return {"label_sets": {"toyset": {
        "task": "subtype", "version": "v0", "label_type": "histology",
        "sources": {"primary": "src_a", "concordance": "src_b" if concordance else None},
        "rules": {"min_class_cases": 10, "min_classes": 2, "component_separator": "|"},
        "cancers": {"TOY": {
            "tier": "primary", "question": "alpha vs beta",
            "classes": {"Alpha": {"src_a": ["Alpha type"], "src_b": ["a"]},
                        "Beta": {"src_a": ["Beta type"], "src_b": ["b"]},
                        "Gamma": {"src_a": ["Gamma type"], "src_b": ["g"]}},
            "exclude": {"src_a": ["Other, specify"], "src_b": ["x"]}}}}},
        "modelling": {"inner_folds_used": 3, "max_epochs": 5, "patience": 2, "n_trials": 2,
                      "batch_size": [32], "bootstrap_resamples": 50}}


def toy_loaders(samples: pd.DataFrame, raw_a: list[str], raw_b: list[str]) -> dict:
    tum = samples[samples.role == "tumour"].reset_index(drop=True)
    a = pd.Series(raw_a, index=tum["case_id"])
    b = pd.Series(raw_b, index=tum["case_id"])

    def mk(series):
        def load(cfg, tumours, cs):
            raw = series.reindex(tumours["case_id"])
            return st._frame(tumours["case_id"], raw, np.where(raw.isna(), "not in source", ""))
        return load
    return {"src_a": mk(a), "src_b": mk(b)}


def toy_world(seed=0):
    rng = np.random.default_rng(seed)
    samples = toy_samples()
    n = int((samples.role == "tumour").sum())
    raw_a = list(rng.choice(RAW_A, size=n, p=[0.4, 0.35, 0.05, 0.1, 0.05, 0.05]))
    to_b = {"Alpha type": "a", "Beta type": "b", "Gamma type": "g", "Other, specify": "x",
            "Alpha type|Beta type": "a", "Alpha type|Alpha type": "a"}
    raw_b = [to_b[r] for r in raw_a]
    raw_b[0], raw_a[0] = "b", "Alpha type"            # one disagreement between the sources
    return samples, raw_a, raw_b


class TestAssignClass(unittest.TestCase):
    cmap = {"A1": "A", "A2": "A", "B1": "B"}
    excl = {"Not Available"}

    def test_rules(self):
        f = st.assign_class
        self.assertEqual(f("A1", self.cmap, self.excl), ("A", "mapped"))
        self.assertEqual(f("A1|A2", self.cmap, self.excl), ("A", "mapped"))
        self.assertEqual(f("A1|B1", self.cmap, self.excl)[0], None)
        self.assertEqual(f("A1|Not Available", self.cmap, self.excl), (None, "excluded value"))
        self.assertEqual(f(np.nan, self.cmap, self.excl), (None, "no value in source"))
        # ignored components are tests not done; only the informative ones must agree
        self.assertEqual(f("B1|[NE]", self.cmap, self.excl, ignore={"[NE]"}), ("B", "mapped"))
        self.assertEqual(f("[NE]|[NE]", self.cmap, self.excl, ignore={"[NE]"}), (None, "no informative value"))

    def test_unlisted_value_is_an_error(self):
        with self.assertRaises(ValueError):
            st.assign_class("C9", self.cmap, self.excl)

    def test_value_listed_twice_is_an_error(self):
        spec = toy_spec()
        cs = st.cohort_specs(spec)["TOY"]
        cs["exclude"]["src_a"] = ["Alpha type"]
        with self.assertRaises(ValueError):
            st.class_map(cs, "src_a")


class TestLabelBuilder(unittest.TestCase):
    def setUp(self):
        self.samples, self.raw_a, self.raw_b = toy_world()
        self.spec = toy_spec()
        self.loaders = toy_loaders(self.samples, self.raw_a, self.raw_b)
        self.tables, self.labels, self.scope, self.acct = st.build_all(None, self.spec, self.samples, self.loaders)

    def test_one_label_per_patient_tumours_only(self):
        t = self.tables["TOY"]
        self.assertFalse(t["case_id"].duplicated().any())
        self.assertTrue((t["role"] == "tumour").all())
        self.assertEqual(list(t.drop_duplicates("y").sort_values("y")["label"]), ["Alpha", "Beta"])

    def test_excluded_and_mixed_never_appear(self):
        t = self.tables["TOY"]
        raw = pd.Series(self.raw_a, index=self.samples[self.samples.role == "tumour"]["case_id"].to_numpy())
        labelled_raw = set(raw[t["case_id"]])
        self.assertTrue(labelled_raw <= {"Alpha type", "Beta type", "Alpha type|Alpha type"})
        self.assertNotIn("Other, specify", labelled_raw)
        self.assertNotIn("Alpha type|Beta type", labelled_raw)

    def test_concordance_and_size_rules(self):
        lab = self.labels.set_index("case_id")
        self.assertEqual(lab.loc["TCGA-AA-0000", "reason"], "sources disagree")
        self.assertNotIn("TCGA-AA-0000", set(self.tables["TOY"]["case_id"]))
        gamma = lab[lab["class_primary"] == "Gamma"]
        self.assertTrue((gamma["reason"] == "class below minimum size").all())
        self.assertNotIn("Gamma", set(self.tables["TOY"]["label"]))
        self.assertIn("Gamma", self.scope.iloc[0]["classes_below_minimum"])

    def test_accounting_is_complete(self):
        n_tum = int((self.samples.role == "tumour").sum())
        for src in ("src_a", "src_b"):
            self.assertEqual(self.acct[self.acct.source == src]["patients"].sum(), n_tum)

    def test_folds_grouped_by_patient_and_stratified(self):
        t = self.tables["TOY"]
        s = make_splits(t, "subtype", "TOY", 5, 3, 1)
        self.assertTrue((s.groupby("case_id")["outer"].nunique() == 1).all())
        for k in range(5):
            tr, te = fold_indices(s, k)
            self.assertFalse(set(s.case_id.iloc[tr]) & set(s.case_id.iloc[te]))
            self.assertEqual(set(s.y.iloc[te]), {0, 1})


class TestXlsxReader(unittest.TestCase):
    def test_shared_inline_numeric_error_cells(self):
        ws = ('<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'
              '<row r="1"><c r="B1" t="s"><v>0</v></c><c r="C1" t="s"><v>1</v></c><c r="D1" t="s"><v>2</v></c></row>'
              '<row r="2"><c r="A2"><v>1</v></c><c r="B2" t="s"><v>3</v></c><c r="C2"><v>58</v></c>'
              '<c r="D2" t="e"><v>#N/A</v></c></row>'
              '<row r="3"><c r="B3" t="inlineStr"><is><t>TCGA-YY-0002</t></is></c><c r="D3" t="s"><v>4</v></c></row>'
              '</sheetData></worksheet>')
        sst = ('<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
               '<si><t>id</t></si><si><t>age</t></si><si><t>hist</t></si><si><t>TCGA-XX-0001</t></si>'
               '<si><r><t>Ductal </t></r><r><t>Carcinoma</t></r></si></sst>')
        wb = ('<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
              'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
              '<sheet name="Notes" sheetId="1" r:id="rId2"/><sheet name="Data" sheetId="2" r:id="rId1"/>'
              '</sheets></workbook>')
        rels = ('<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                '<Relationship Id="rId1" Target="worksheets/sheet7.xml"/>'
                '<Relationship Id="rId2" Target="worksheets/sheet1.xml"/></Relationships>')
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "t.xlsx"
            with zipfile.ZipFile(p, "w") as z:
                z.writestr("xl/workbook.xml", wb)
                z.writestr("xl/_rels/workbook.xml.rels", rels)
                z.writestr("xl/sharedStrings.xml", sst)
                z.writestr("xl/worksheets/sheet7.xml", ws)
                z.writestr("xl/worksheets/sheet1.xml", ws.replace("58", "99"))
            t = st.read_xlsx_sheet(p, "Data")
        self.assertEqual(list(t.columns), ["id", "age", "hist"])     # the unnamed index column is dropped
        self.assertEqual(t.values.tolist(), [["TCGA-XX-0001", "58", "#N/A"], ["TCGA-YY-0002", "", "Ductal Carcinoma"]])


class TestMetrics(unittest.TestCase):
    def test_perfect_majority_and_chance(self):
        y = np.array([0] * 30 + [1] * 10 + [2] * 10)
        P = np.eye(3)[y] * 0.9 + 0.1 / 3
        m = st.subtype_metrics(y, P)
        self.assertAlmostEqual(m["balanced_accuracy"], 1.0)
        self.assertAlmostEqual(m["macro_auroc_ovr"], 1.0)
        b = st.baselines(y, 3)
        self.assertAlmostEqual(b["majority_accuracy"], 0.6)
        self.assertAlmostEqual(b["majority_balanced_accuracy"], 1 / 3)
        self.assertAlmostEqual(b["majority_macro_f1"], (2 * 0.6 / 1.6) / 3)

    def test_binary_auroc_uses_second_column(self):
        y = np.array([0, 0, 1, 1])
        P = np.column_stack([1 - np.array([.1, .2, .8, .9]), [.1, .2, .8, .9]])
        self.assertAlmostEqual(st.subtype_metrics(y, P)["macro_auroc_ovr"], 1.0)


class TestBinaryFoldThroughPipeline(unittest.TestCase):
    """A two-class subtype cohort goes through pipeline.run_fold with a one-logit model, and its
    probabilities are proper (the former `binary = task == "tumor_status"` would have given 1.0)."""

    def test_two_and_three_classes(self):
        from ctsm.pipeline import run_fold
        base = load_config()
        for C in (2, 3):
            spec = toy_spec()
            scfg = st.subtype_config(base, spec)
            scfg.study["tuning"]["search_space"]["n_features"] = [5, 10]
            rng = np.random.default_rng(C)
            n = 120
            y = np.arange(n) % C
            X = rng.normal(size=(n, 30)).astype(np.float32)
            X[:, :3] += y[:, None] * 1.5
            tab = pd.DataFrame({"file_id": [f"f{i}" for i in range(n)], "case_id": [f"c{i}" for i in range(n)], "y": y})
            spl = make_splits(tab, "subtype", f"TOY{C}", 5, 3, 3)
            feats = pd.DataFrame({"feature_id": [f"g{j}" for j in range(30)]})
            with tempfile.TemporaryDirectory() as d:
                rec = run_fold(scfg, "subtype", "se", f"TOY{C}", 0, tab, spl, ViewMatrix(X, "float32"), feats,
                               Path(d) / "fold0", lambda *_: None, fold_metrics=st.fold_metrics)
                z = np.load(Path(d) / "fold0" / "predictions.npz")
            self.assertEqual(z["proba"].shape[1], C)
            np.testing.assert_allclose(z["proba"].sum(axis=1), 1.0, atol=1e-6)
            self.assertLess(z["proba"].max(), 1.0)
            self.assertIn("test_balanced_accuracy", rec)
            self.assertLess(rec["shap_additivity_max_error"], 1e-3)


# ---- the real cohort (skipped where the local data are absent) ---------------------------------

CFG = load_config()
HAVE_COHORT = (CFG.data_dir / "cohort" / "samples.parquet").exists()
HAVE_SUBTYPE = (st.subtype_config(CFG).data_dir / "tasks").exists()


@unittest.skipUnless(HAVE_COHORT, "local cohort not built")
class TestExistingTasksUnchanged(unittest.TestCase):
    def test_task_tables_rebuild_identically(self):
        from ctsm.cohort import task_scope, task_tables
        samples = pd.read_parquet(CFG.data_dir / "cohort" / "samples.parquet")
        scope = task_scope(CFG, samples)
        for (task, cohort), tab in task_tables(CFG, samples, scope).items():
            stored = pd.read_parquet(CFG.data_dir / "tasks" / f"{task}__{cohort}.parquet")
            pd.testing.assert_frame_equal(tab, stored, check_dtype=False)

    def test_subtype_config_leaves_study_untouched(self):
        before = load_config().study
        scfg = st.subtype_config(CFG)
        self.assertEqual(load_config().study, before)
        self.assertNotIn("subtype", before["tuning"]["n_trials"])
        self.assertEqual(scfg.study["tuning"]["n_trials"]["tumor_status"], before["tuning"]["n_trials"]["tumor_status"])
        self.assertNotEqual(scfg.data_dir, CFG.data_dir)


@unittest.skipUnless(HAVE_SUBTYPE, "subtype labels not built (scripts/30_subtype_labels.py)")
class TestRealSubtypeLabels(unittest.TestCase):
    def test_tables(self):
        scfg = st.subtype_config(CFG)
        specs = st.cohort_specs(scfg.spec)
        labels = pd.read_parquet(scfg.data_dir / "labels" / "labels.parquet")
        samples = pd.read_parquet(CFG.data_dir / "cohort" / "samples.parquet")
        tumours = set(samples.loc[samples.role == "tumour", "file_id"])
        for cohort, tab in st.load_tables(scfg).items():
            cs = specs[cohort]
            self.assertFalse(tab["case_id"].duplicated().any(), cohort)
            self.assertTrue(set(tab["file_id"]) <= tumours, cohort)
            self.assertTrue((tab["project"] == cs["project"]).all(), cohort)
            counts = tab["label"].value_counts()
            self.assertGreaterEqual(counts.min(), cs["rules"]["min_class_cases"], cohort)
            lab = labels[(labels.cohort == cohort)].set_index("case_id").loc[tab["case_id"]]
            excluded = set(cs["exclude"].get(cs["source"], []))
            sep = cs["rules"]["component_separator"]
            for raw, cls in zip(lab["raw"], tab["label"].to_numpy()):
                comps = [c for c in str(raw).split(sep) if c and c not in set(cs["ignore"])]
                self.assertFalse(set(comps) & excluded, cohort)
                self.assertTrue(all(c in cs["classes"][cls][cs["source"]] for c in comps), (cohort, raw, cls))
            s = pd.read_parquet(scfg.data_dir / "splits" / f"subtype__{cohort}.parquet")
            self.assertTrue((s["file_id"].to_numpy() == tab["file_id"].to_numpy()).all())
            self.assertTrue((s.groupby("case_id")["outer"].nunique() == 1).all())


if __name__ == "__main__":
    unittest.main()
