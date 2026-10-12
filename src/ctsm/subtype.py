"""Within-cancer subtype task (docs/SUBTYPE_METHODS.md, decisions D-30 onwards).

Labels
  Every label comes from the pre-specified maps of config/subtype.yaml. A label set (histology,
  molecular) lists cohorts; a cohort is one TCGA project's primary tumours from the D-01 cohort
  (one sample per patient), each patient labelled from a label source (LABEL_SOURCES) through an
  exhaustive table of raw values. A raw value that the table does not list is an error, so a typo
  in the map cannot drop or relabel patients silently. Labels come only from clinical or published
  annotation tables, never from this repository's molecular matrices (leakage rule R6).
Modelling
  One model per (cohort, view), fitted by the unchanged D-04..D-10 machinery (pipeline.run_fold):
  patient-grouped, label-stratified nested CV, Optuna on inner folds, fold-level preprocessing,
  DeepExplainer SHAP on held-out samples. The matrices are a dedicated extraction whose D-07
  pre-filter is recomputed on the subtype cohorts, so the superset proof of LEAKAGE_POLICY.md
  applies to them exactly (D-37).
Paths
  <data_dir>/subtype/{labels,tasks,splits,matrices,runs}/ (git-ignored), results/subtype/
"""
from __future__ import annotations

import copy
import functools
import hashlib
import json
import re
import time
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from sklearn import metrics as skm

from .cohort import class_names
from .config import Config, load_config, load_yaml, stable_hash
from .provenance import git_commit, package_versions, write_json

TASK = "subtype"
SUBTYPE_YAML = "subtype.yaml"
PANCAN_TYPE = {"LAML": "AML", "OV": "OVCA"}          # TCGA project -> cancer.type of TCGASubtype
LABEL_TABLE_COLUMNS = ["file_id", "case_id", "sample_id", "project", "sample_type", "role", "tss", "vial"]


# ==== configuration ===========================================================================

def subtype_spec() -> dict:
    return load_yaml(SUBTYPE_YAML)


@dataclass(frozen=True)
class SubtypeConfig(Config):
    """Config whose data directory is <data_dir>/subtype and whose study block carries the
    subtype entries of the per-task settings (D-36). config/study.yaml itself is not changed, so
    the two existing tasks see exactly the configuration they were run with."""
    spec: dict = field(default_factory=dict)

    @property
    def base_data_dir(self) -> Path:
        return Config.data_dir.fget(self)

    @property
    def data_dir(self) -> Path:
        return self.base_data_dir / TASK

    @property
    def fingerprint(self) -> str:
        return stable_hash({"study": self.study, "sources": self.sources, "subtype": self.spec})


def subtype_config(cfg: Config | None = None, spec: dict | None = None) -> SubtypeConfig:
    cfg = cfg or load_config()
    spec = spec or subtype_spec()
    st = copy.deepcopy(cfg.study)
    m = spec["modelling"]
    st["splits"]["inner_folds_used"][TASK] = int(m["inner_folds_used"])
    st["model"]["max_epochs"][TASK] = int(m["max_epochs"])
    st["model"]["patience"][TASK] = int(m["patience"])
    st["tuning"]["n_trials"][TASK] = int(m["n_trials"])
    st["tuning"]["search_space"]["batch_size"][TASK] = list(m["batch_size"])
    return SubtypeConfig(study=st, sources=cfg.sources, spec=spec)


def cohort_specs(spec: dict) -> dict[str, dict]:
    """{cohort: resolved specification} over every label set of the subtype task, in YAML order.
    Cohort names are unique across label sets (asserted)."""
    out: dict[str, dict] = {}
    for set_name, ls in spec["label_sets"].items():
        if ls["task"] != TASK:
            continue
        for cohort, c in ls["cancers"].items():
            if cohort in out:
                raise ValueError(f"cohort {cohort} appears in two label sets")
            src = c.get("source", ls["sources"]["primary"])
            out[cohort] = {
                "cohort": cohort, "label_set": set_name, "version": ls["version"],
                "project": c.get("project", cohort), "tier": c["tier"], "question": c["question"],
                "label_type": c.get("label_type", ls.get("label_type")),
                "source": src, "concordance": c.get("concordance", ls["sources"].get("concordance")),
                "columns": list(c.get("columns", [])), "ignore": list(c.get("ignore", [])),
                "classes": {k: dict(v) for k, v in c["classes"].items()},
                "exclude": dict(c.get("exclude", {})), "rules": dict(ls["rules"]),
            }
    return out


# ==== xlsx (standard library only) =============================================================

_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
       "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def _col_index(ref: str) -> int:
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_xlsx_sheet(path, sheet: str) -> pd.DataFrame:
    """One worksheet of an .xlsx as strings, without openpyxl. An xlsx is a zip: the workbook maps
    sheet names to worksheet XML files; a cell with t="s" indexes xl/sharedStrings.xml, t="str"
    and t="inlineStr" carry text, t="e" an error such as #N/A (kept as its text), and any other
    cell its raw <v> text. The first row is the header; columns with an empty header (e.g. a row
    number column) are dropped; missing cells are empty strings."""
    m = _NS["m"]
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rid = next((s.get(f"{{{_NS['r']}}}id") for s in wb.find("m:sheets", _NS) if s.get("name") == sheet), None)
        if rid is None:
            raise KeyError(f"{path}: no sheet {sheet!r}")
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next(r.get("Target") for r in rels if r.get("Id") == rid).lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        sst = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", _NS):
                sst.append("".join(t.text or "" for t in si.iter(f"{{{m}}}t")))
        rows = []
        with z.open(target) as fh:
            for _, el in ET.iterparse(fh, events=("end",)):
                if el.tag != f"{{{m}}}row":
                    continue
                cells = {}
                for c in el.findall("m:c", _NS):
                    t, v = c.get("t"), c.find("m:v", _NS)
                    if t == "s":
                        val = sst[int(v.text)]
                    elif t == "inlineStr":
                        val = "".join(x.text or "" for x in c.iter(f"{{{m}}}t"))
                    else:
                        val = v.text if v is not None and v.text is not None else ""
                    cells[_col_index(c.get("r"))] = val
                rows.append(cells)
                el.clear()
    width = max((max(r) + 1 for r in rows if r), default=0)
    grid = [[r.get(i, "") for i in range(width)] for r in rows]
    header, body = grid[0], grid[1:]
    keep = [i for i, h in enumerate(header) if h != ""]
    if len({header[i] for i in keep}) != len(keep):
        raise ValueError(f"{path}:{sheet}: duplicate column names")
    return pd.DataFrame([[r[i] for i in keep] for r in body], columns=[header[i] for i in keep], dtype=str)


# ==== label sources ============================================================================
# A loader takes (cfg, tumours, cohort spec) and returns a frame indexed like tumours["case_id"]
# with columns `raw` (the source's value, or NaN) and `note` (why it is missing, if it is).

def _pancan_dir(cfg: Config) -> Path:
    return Config.data_dir.fget(cfg) / "external" / "gdc_pancanatlas"


@functools.lru_cache(maxsize=None)
def _verified_md5(path: str) -> str:
    """Re-hash a file fetched by scripts/00b_fetch_gdc_labels.py against its MANIFEST.json."""
    p = Path(path)
    man = json.loads((p.parent / "MANIFEST.json").read_text())[p.name]
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != man["md5"] or p.stat().st_size != int(man["bytes"]):
        raise ValueError(f"{p.name}: does not match MANIFEST.json; re-run scripts/00b_fetch_gdc_labels.py")
    return h.hexdigest()


@functools.lru_cache(maxsize=None)
def _cdr_table(path: str) -> pd.DataFrame:
    t = read_xlsx_sheet(path, "TCGA-CDR")
    if t["bcr_patient_barcode"].duplicated().any():
        raise ValueError("TCGA-CDR has duplicate patients")
    return t.set_index("bcr_patient_barcode")


@functools.lru_cache(maxsize=None)
def _gdc_diagnosis(path: str) -> pd.Series:
    c = pd.read_csv(path, usecols=["File.ID", "primary_diagnosis"], dtype=str)
    return c.set_index("File.ID")["primary_diagnosis"]


@functools.lru_cache(maxsize=None)
def _pancan_subtype_table(path: str) -> pd.DataFrame:
    _verified_md5(path)
    return pd.read_csv(path, sep="\t", dtype=str)


@functools.lru_cache(maxsize=None)
def _pancan_clinical_table(path: str) -> pd.DataFrame:
    _verified_md5(path)
    t = pd.read_csv(path, sep="\t", encoding="latin-1", dtype=str)
    if t["bcr_patient_barcode"].duplicated().any():
        raise ValueError("clinical_PANCAN has duplicate patients")
    return t.set_index("bcr_patient_barcode")


def _frame(case_ids, raw, note) -> pd.DataFrame:
    return pd.DataFrame({"raw": list(raw), "note": list(note)}, index=pd.Index(list(case_ids), name="case_id"))


def load_cdr_histology(cfg: Config, tumours: pd.DataFrame, cs: dict) -> pd.DataFrame:
    t = _cdr_table(str(cfg.sources["subtype_labels"]["tcga_cdr"]))
    sub = t.reindex(tumours["case_id"])
    found = sub["type"].notna()
    if (sub.loc[found, "type"] != cs["project"]).any():
        raise ValueError(f"{cs['cohort']}: TCGA-CDR type differs from the project")
    raw = sub["histological_type"].where(found)
    return _frame(tumours["case_id"], raw, np.where(found, "", "not in source"))


def load_gdc_diagnosis(cfg: Config, tumours: pd.DataFrame, cs: dict) -> pd.DataFrame:
    s = _gdc_diagnosis(str(cfg.tcga_path("clinical")))
    raw = s.reindex(tumours["file_id"]).to_numpy()
    return _frame(tumours["case_id"], raw, np.where(pd.isna(raw), "no value in source", ""))


def _join_columns(t: pd.DataFrame, columns: list[str], sep: str) -> pd.Series:
    """Several columns -> one component string; a missing cell is an empty component."""
    parts = [t[c].fillna("") for c in columns]
    out = parts[0]
    for p in parts[1:]:
        out = out + sep + p
    return out


def load_pancan_subtype(cfg: Config, tumours: pd.DataFrame, cs: dict) -> pd.DataFrame:
    """TCGASubtype.20170308.tsv. Rows are patient barcodes or sample/aliquot barcodes; only
    patient rows and primary-tumour rows (sample type 01; 03 for LAML) of the cohort's cancer type
    are used. A patient with two usable rows that disagree is left unlabelled (D-43)."""
    t = _pancan_subtype_table(str(_pancan_dir(cfg) / "TCGASubtype.20170308.tsv"))
    ids = t["pan.samplesID"]
    stype = ids.str[13:15].where(ids.str.len() >= 15)
    code = "03" if cs["project"] == "LAML" else "01"
    keep = (stype.isna() | (stype == code)) & (t["cancer.type"] == PANCAN_TYPE.get(cs["project"], cs["project"]))
    d = pd.DataFrame({"case_id": ids.str[:12], "raw": _join_columns(t, cs["columns"], cs["rules"]["component_separator"])})
    d = d[keep.to_numpy()].drop_duplicates()
    n = d.groupby("case_id")["raw"].nunique()
    conflict = set(n[n > 1].index)
    d = d[~d["case_id"].isin(conflict)].drop_duplicates("case_id").set_index("case_id")["raw"]
    raw = d.reindex(tumours["case_id"])
    note = np.where(tumours["case_id"].isin(conflict), "conflicting rows in source",
                    np.where(raw.isna(), "not in source", ""))
    return _frame(tumours["case_id"], raw.where(~tumours["case_id"].isin(conflict).to_numpy()), note)


def load_pancan_clinical(cfg: Config, tumours: pd.DataFrame, cs: dict) -> pd.DataFrame:
    t = _pancan_clinical_table(str(_pancan_dir(cfg) / "clinical_PANCAN_patient_with_followup.tsv"))
    sub = t.reindex(tumours["case_id"])
    found = sub["acronym"].notna()
    if (sub.loc[found, "acronym"] != cs["project"]).any():
        raise ValueError(f"{cs['cohort']}: clinical_PANCAN acronym differs from the project")
    raw = _join_columns(sub, cs["columns"], cs["rules"]["component_separator"]).where(found)
    return _frame(tumours["case_id"], raw, np.where(found, "", "not in source"))


LABEL_SOURCES = {
    "cdr_histological_type": load_cdr_histology,
    "gdc_primary_diagnosis": load_gdc_diagnosis,
    "pancan_subtype": load_pancan_subtype,
    "pancan_clinical": load_pancan_clinical,
}


# ==== raw value -> class =======================================================================

def class_map(cs: dict, source: str) -> tuple[dict[str, str], set[str]]:
    """({raw value: class}, excluded raw values) of one source. A raw value listed twice, or both
    as a class value and as excluded, is an error."""
    class_of: dict[str, str] = {}
    for cls, by_src in cs["classes"].items():
        for v in by_src.get(source, []):
            if v in class_of:
                raise ValueError(f"{cs['cohort']}/{source}: {v!r} listed under two classes")
            class_of[v] = cls
    excluded = set(cs["exclude"].get(source, []))
    if excluded & set(class_of):
        raise ValueError(f"{cs['cohort']}/{source}: {sorted(excluded & set(class_of))} both mapped and excluded")
    return class_of, excluded


def assign_class(raw, class_of: dict[str, str], excluded: set[str], ignore: set[str] | None = None,
                 sep: str | None = "|") -> tuple[str | None, str]:
    """(class or None, reason). The value is split into components at `sep`; empty and `ignore`
    components are dropped (a test not done); any `excluded` component excludes the patient; the
    rest must all map to one class. An unlisted component raises: the maps are exhaustive (D-32)."""
    if raw is None or (isinstance(raw, float) and np.isnan(raw)) or pd.isna(raw):
        return None, "no value in source"
    comps = str(raw).split(sep) if sep else [str(raw)]
    comps = [c for c in comps if c != "" and c not in (ignore or set())]
    if not comps:
        return None, "no informative value"
    unknown = sorted({c for c in comps if c not in class_of and c not in excluded})
    if unknown:
        raise ValueError(f"raw value(s) not listed in config/subtype.yaml: {unknown}")
    if any(c in excluded for c in comps):
        return None, "excluded value"
    classes = {class_of[c] for c in comps}
    if len(classes) > 1:
        return None, "components map to different classes"
    return classes.pop(), "mapped"


def _apply_source(cfg, tumours, cs, source, loaders) -> pd.DataFrame:
    d = loaders[source](cfg, tumours, cs)
    class_of, excluded = class_map(cs, source)
    sep = cs["rules"].get("component_separator")
    res = [assign_class(r, class_of, excluded, set(cs["ignore"]), sep) for r in d["raw"]]
    d["cls"] = [c for c, _ in res]
    d["reason"] = [n if (why == "no value in source" and n) else why for (_, why), n in zip(res, d["note"])]
    return d


def label_cohort(cfg: Config, cs: dict, samples: pd.DataFrame, loaders: dict | None = None
                 ) -> tuple[pd.DataFrame, list[str], bool]:
    """Label one cohort. Returns (per-patient table with raw values, classes and reasons; the
    classes kept, in YAML order; whether the cohort is in scope)."""
    loaders = loaders or LABEL_SOURCES
    tum = samples[(samples["project"] == cs["project"]) & (samples["role"] == "tumour")]
    if tum["case_id"].duplicated().any():
        raise ValueError(f"{cs['cohort']}: more than one primary tumour per patient")
    out = tum[LABEL_TABLE_COLUMNS].reset_index(drop=True).copy()
    p = _apply_source(cfg, out, cs, cs["source"], loaders)
    out["raw"], out["class_primary"], out["reason"] = p["raw"].to_numpy(), p["cls"].to_numpy(), p["reason"].to_numpy()
    if cs["concordance"]:
        c = _apply_source(cfg, out, cs, cs["concordance"], loaders)
        out["raw_concordance"], out["class_concordance"] = c["raw"].to_numpy(), c["cls"].to_numpy()
        mapped = out["reason"] == "mapped"
        no_conc = mapped & (c["reason"].to_numpy() != "mapped")
        out.loc[no_conc, "reason"] = "concordance source: " + c["reason"][no_conc.to_numpy()].to_numpy()
        disagree = mapped & ~no_conc & (out["class_concordance"] != out["class_primary"])
        out.loc[disagree, "reason"] = "sources disagree"
    out["label"] = out["class_primary"].where(out["reason"] == "mapped")
    counts = out.groupby("label")["case_id"].nunique()
    n_min = int(cs["rules"]["min_class_cases"])
    small = set(counts[counts < n_min].index)
    out.loc[out["label"].isin(small), "reason"] = "class below minimum size"
    out["label"] = out["label"].where(~out["label"].isin(small))
    kept = [k for k in cs["classes"] if counts.get(k, 0) >= n_min]
    return out, kept, len(kept) >= int(cs["rules"]["min_classes"])


def task_table(labels: pd.DataFrame, kept: list[str]) -> pd.DataFrame:
    """Rows of labelled patients in their original order, with `label` and integer `y`."""
    d = labels[labels["label"].isin(kept)][LABEL_TABLE_COLUMNS + ["label"]].copy()
    d["y"] = d["label"].map({c: i for i, c in enumerate(kept)}).astype(int)
    if d["case_id"].duplicated().any():
        raise AssertionError("a patient has two labels")
    return d.reset_index(drop=True)


def build_all(cfg: Config, spec: dict, samples: pd.DataFrame, loaders: dict | None = None
              ) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Every cohort of the subtype task -> ({cohort: task table} for the cohorts in scope,
    per-patient label table of all cohorts, scope table, raw-value accounting)."""
    tables, labs, scope, acct = {}, [], [], []
    for cohort, cs in cohort_specs(spec).items():
        lab, kept, ok = label_cohort(cfg, cs, samples, loaders)
        labs.append(lab.assign(cohort=cohort))
        n_by = lab.groupby("label")["case_id"].nunique()
        dropped = [k for k in cs["classes"] if k not in kept]
        scope.append({"cohort": cohort, "label_set": cs["label_set"], "label_type": cs["label_type"],
                      "tier": cs["tier"], "project": cs["project"], "question": cs["question"],
                      "source": cs["source"], "concordance": cs["concordance"] or "",
                      "primary_tumours": len(lab), "labelled": int(lab["label"].isin(kept).sum()),
                      "n_classes": len(kept), "classes": "; ".join(f"{k} {int(n_by.get(k, 0))}" for k in kept),
                      "classes_below_minimum": "; ".join(
                          f"{k} {int((lab['class_primary'] == k).sum())}" for k in dropped),
                      "in_scope": ok})
        for src, rawcol, clscol in [(cs["source"], "raw", "class_primary"),
                                    (cs["concordance"], "raw_concordance", "class_concordance")]:
            if not src:
                continue
            g = lab.groupby(lab[rawcol].fillna("<missing>"), dropna=False)
            for raw, gg in g:
                acct.append({"cohort": cohort, "source": src, "raw_value": raw,
                             "class": gg[clscol].iloc[0] if gg[clscol].notna().any() else "",
                             "patients": len(gg), "labelled": int(gg["label"].notna().sum())})
        if ok:
            tables[cohort] = task_table(lab, kept)
    return tables, pd.concat(labs, ignore_index=True), pd.DataFrame(scope), pd.DataFrame(acct)


# ==== metrics ==================================================================================

def subtype_metrics(y: np.ndarray, proba: np.ndarray) -> dict:
    """Balanced accuracy, macro F1 and macro one-vs-rest AUROC (the AUROC itself for 2 classes)
    over all C classes, plus accuracy, MCC, log loss and per-class recall."""
    C = proba.shape[1]
    labels = np.arange(C)
    pred = proba.argmax(axis=1)
    if C == 2:
        auroc = float(skm.roc_auc_score(y, proba[:, 1]))
    else:
        auroc = float(skm.roc_auc_score(y, proba, multi_class="ovr", average="macro", labels=labels))
    rec = skm.recall_score(y, pred, labels=labels, average=None, zero_division=0)
    return {"accuracy": float(skm.accuracy_score(y, pred)),
            "balanced_accuracy": float(np.mean(rec)),
            "macro_f1": float(skm.f1_score(y, pred, labels=labels, average="macro", zero_division=0)),
            "macro_auroc_ovr": auroc, "mcc": float(skm.matthews_corrcoef(y, pred)),
            "log_loss": float(skm.log_loss(y, np.clip(proba, 1e-12, 1), labels=labels)),
            "n": int(len(y)), "n_classes": int(C), **{f"recall_{c}": float(rec[c]) for c in range(C)}}


def fold_metrics(y: np.ndarray, proba: np.ndarray) -> dict:
    m = subtype_metrics(y, proba)
    return {k: m[k] for k in ("balanced_accuracy", "macro_f1", "macro_auroc_ovr", "accuracy")}


def bootstrap(y: np.ndarray, proba: np.ndarray, groups: np.ndarray, n_boot: int = 1000, seed: int = 0) -> dict:
    """Percentile 95 % CIs from resampling patients with replacement. A resample that misses a
    class is skipped (it happens only for the smallest classes), and the number used is reported."""
    rng = np.random.default_rng(seed)
    C = proba.shape[1]
    uniq, inv = np.unique(groups, return_inverse=True)
    members = [np.flatnonzero(inv == g) for g in range(len(uniq))]
    keys = ["accuracy", "balanced_accuracy", "macro_f1", "macro_auroc_ovr"] + [f"recall_{c}" for c in range(C)]
    vals = {k: [] for k in keys}
    for _ in range(n_boot):
        idx = np.concatenate([members[g] for g in rng.integers(0, len(uniq), len(uniq))])
        if len(np.unique(y[idx])) < C:
            continue
        m = subtype_metrics(y[idx], proba[idx])
        for k in keys:
            vals[k].append(m[k])
    out = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v else [None, None]
           for k, v in vals.items()}
    out["n_resamples_used"] = len(vals["accuracy"])
    return out


def baselines(y: np.ndarray, n_classes: int) -> dict:
    """Majority-class predictor and chance levels. Chance: balanced accuracy 1/C, macro AUROC
    0.5, macro F1 1/C (expected value of guessing each class at its prevalence)."""
    pi = np.bincount(y, minlength=n_classes) / len(y)
    maj = int(np.argmax(pi))
    pred = np.full(len(y), maj)
    return {"chance_balanced_accuracy": 1 / n_classes, "chance_macro_auroc": 0.5, "chance_macro_f1": 1 / n_classes,
            "majority_class_index": maj, "majority_accuracy": float(pi[maj]),
            "majority_balanced_accuracy": 1 / n_classes,
            "majority_macro_f1": float(skm.f1_score(y, pred, labels=np.arange(n_classes), average="macro",
                                                    zero_division=0))}


# ==== one job ==================================================================================

def run_job(scfg: SubtypeConfig, view: str, cohort: str, log=print) -> dict:
    """Nested CV for one (cohort, view) with pipeline.run_fold, then the subtype summaries."""
    from .pipeline import load_job_data, run_fold
    t0 = time.time()
    tab, spl, X, feats = load_job_data(scfg, TASK, view, cohort)
    classes = class_names(tab)
    log(f"{TASK}/{view}/{cohort}: X {X.shape} ({X.encoding}), classes {classes}, missing {X.missing_frac():.3f}")
    rdir = scfg.run_dir(TASK, view, cohort)
    recs = []
    for k in range(scfg.study["splits"]["outer_folds"]):
        fd = rdir / f"fold{k}"
        if (fd / "done.json").exists():
            recs.append(json.loads((fd / "done.json").read_text()))
            log(f"  fold {k}: already done")
            continue
        recs.append(run_fold(scfg, TASK, view, cohort, k, tab, spl, X, feats, fd, log, fold_metrics=fold_metrics))
    summary = summarise_job(scfg, view, cohort, tab, classes, recs)
    summarise_shap(scfg, view, cohort, tab, classes, feats)
    summary["seconds_wall_this_run"] = round(time.time() - t0, 1)
    write_json(rdir / "job_done.json", summary)
    return summary


def summarise_job(scfg: SubtypeConfig, view: str, cohort: str, tab: pd.DataFrame, classes: list[str],
                  recs: list[dict]) -> dict:
    """Pooled out-of-fold metrics, patient-bootstrap CIs and baselines -> results/subtype/<view>/<cohort>/."""
    from .metrics import per_class_table
    from .pipeline import collect_predictions
    y = tab["y"].to_numpy()
    C = len(classes)
    proba = collect_predictions(scfg, TASK, view, cohort, tab, C)
    rdir, out = scfg.run_dir(TASK, view, cohort), scfg.result_dir(TASK, view, cohort)
    out.mkdir(parents=True, exist_ok=True)
    pred = tab[["file_id", "case_id", "sample_id", "project", "label", "y"]].copy()
    pred["outer_fold"] = pd.read_parquet(scfg.data_dir / "splits" / f"{TASK}__{cohort}.parquet")["outer"].to_numpy()
    for c, name in enumerate(classes):
        pred[f"p_{name}"] = proba[:, c]
    pred.to_parquet(rdir / "oof_predictions.parquet", index=False)        # per patient: kept local (D-12)

    pooled = subtype_metrics(y, proba)
    ci = bootstrap(y, proba, tab["case_id"].to_numpy(), n_boot=int(scfg.spec["modelling"]["bootstrap_resamples"]))
    pc = pd.DataFrame(per_class_table(y, proba, classes))
    pc["recall_lo"] = [ci[f"recall_{c}"][0] for c in range(C)]
    pc["recall_hi"] = [ci[f"recall_{c}"][1] for c in range(C)]
    pc.to_csv(out / "per_class_metrics.csv", index=False)
    cm = pd.crosstab(pd.Categorical(np.array(classes)[y], categories=classes),
                     pd.Categorical(np.array(classes)[proba.argmax(axis=1)], categories=classes), dropna=False)
    cm.index.name, cm.columns.name = "true", "predicted"
    cm.to_csv(out / "confusion_matrix.csv")
    folds = pd.DataFrame(recs)
    folds.to_csv(out / "folds.csv", index=False)
    test_cols = [c for c in folds.columns if c.startswith("test_")]
    cs = cohort_specs(scfg.spec)[cohort]
    base = baselines(y, C)
    base["majority_class"] = classes[base.pop("majority_class_index")]
    summary = {
        "task": TASK, "view": view, "cohort": cohort, "classes": classes,
        "class_counts": {k: int((y == i).sum()) for i, k in enumerate(classes)},
        "labels": {k: cs[k] for k in ("label_set", "version", "label_type", "tier", "project", "question",
                                      "source", "concordance", "columns")},
        "pooled_oof": {k: v for k, v in pooled.items() if not k.startswith("recall_")},
        "recall_by_class": {name: pooled[f"recall_{c}"] for c, name in enumerate(classes)},
        "bootstrap_95ci": {k: v for k, v in ci.items() if not k.startswith("recall_")},
        "baselines": base,
        "per_fold_mean": {c[5:]: float(folds[c].mean()) for c in test_cols},
        "per_fold_sd": {c[5:]: float(folds[c].std(ddof=1)) for c in test_cols},
        "selection_method_by_fold": folds["param_selection"].tolist(),
        "n_features_by_fold": folds["n_features_used"].tolist(),
        "git_commit": git_commit(), "packages": package_versions(),
        "config": {"splits": scfg.study["splits"], "preprocess": scfg.study["preprocess"],
                   "model": scfg.study["model"], "tuning": scfg.study["tuning"], "shap": scfg.study["shap"],
                   "subtype_yaml_sha256": stable_hash(scfg.spec)},
    }
    write_json(out / "metrics.json", summary)
    return summary


def summarise_shap(scfg: SubtypeConfig, view: str, cohort: str, tab: pd.DataFrame, classes: list[str],
                   feats: pd.DataFrame) -> None:
    """Fold SHAP files -> shap_top_features.csv, shap_top_genes.csv, shap_top_features_by_class.csv.
    Explained output (D-38): two classes, the log-odds of the second class; more, each sample's
    centred true-class logit. In the by-class table direction > 0 always means that a higher value
    raises the evidence for that class (for the first of two classes the sign is flipped)."""
    from .explain import gene_table, importance_table, load_folds
    folds = load_folds(scfg.run_dir(TASK, view, cohort), tab["y"].to_numpy(), scfg.study["splits"]["outer_folds"])
    ann = feats[["feature_id", "label", "symbol", "gene_id", "region"]]
    out = scfg.result_dir(TASK, view, cohort)
    top = int(scfg.study["shap"]["top_k_report"])
    binary = len(classes) == 2
    t = importance_table(folds, ann, n_direction=200 if binary else 0)
    if binary:
        t = t.rename(columns={"direction_spearman": f"direction_towards_{classes[1]}"})
    else:
        t = t.drop(columns="direction_spearman")
    t.head(200).to_csv(out / "shap_top_features.csv", index=False)
    gene_table(t, "mean_abs_shap").head(200).to_csv(out / "shap_top_genes.csv", index=False)
    per = []
    for k, name in enumerate(classes):
        tk = importance_table(folds, ann, mask_of=lambda f, k=k: f["y"] == k, n_direction=top).head(top)
        if binary and k == 0:
            tk["direction_spearman"] = -tk["direction_spearman"]
        per.append(tk.rename(columns={"rank": "rank_in_class"}).assign(**{"class": name}))
    pd.concat(per, ignore_index=True).to_csv(out / "shap_top_features_by_class.csv", index=False)


def load_tables(scfg: SubtypeConfig) -> dict[str, pd.DataFrame]:
    """{cohort: task table} of every cohort written by step 30, in YAML order."""
    out = {}
    for cohort in cohort_specs(scfg.spec):
        p = scfg.data_dir / "tasks" / f"{TASK}__{cohort}.parquet"
        if p.exists():
            out[cohort] = pd.read_parquet(p)
    return out
