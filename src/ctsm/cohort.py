"""Sample selection and task labels (decisions D-01..D-03).

Selection uses sample metadata only (sample type, barcode, file availability), never a
molecular value.
"""
from __future__ import annotations

import pandas as pd

from .config import PAN, Config


def read_csv_header(path) -> list[str]:
    with open(path) as fh:
        return fh.readline().rstrip("\n").split(",")


def molecular_file_ids(cfg: Config) -> dict[str, set[str]]:
    """File UUIDs per view. The AS tables must share one header (asserted)."""
    out, as_header = {}, None
    for view in cfg.views:
        hdr = read_csv_header(cfg.raw_matrix(view))
        if hdr[0] != "id" or len(hdr) != len(set(hdr)):
            raise ValueError(f"{view}: header must start with 'id' and have unique UUIDs")
        if view != "gex":
            if as_header is not None and hdr != as_header:
                raise ValueError(f"AS table {view} has a different header from the others")
            as_header = hdr
        out[view] = set(hdr[1:])
    return out


def parse_barcode(sample_id: pd.Series) -> pd.DataFrame:
    """TCGA-XX-YYYY-STV: XX tissue source site, ST sample type code, V vial letter."""
    parts = sample_id.str.extract(r"^TCGA-(\w{2})-(\w{4})-(\d{2})([A-Z])$")
    if parts.isna().any().any():
        raise ValueError(f"unparseable barcodes, e.g. {sample_id[parts.isna().any(axis=1)].tolist()[:5]}")
    return pd.DataFrame({"tss": parts[0], "sample_type_code": parts[2], "vial": parts[3]}, index=sample_id.index)


def load_clinical(cfg: Config) -> pd.DataFrame:
    clin = pd.read_csv(cfg.tcga_path("clinical"), usecols=["File.ID", "Project.ID", "Case.ID", "Sample.ID", "Sample.Type"])
    if clin["File.ID"].duplicated().any():
        raise ValueError("clinical table has duplicate File.ID rows")
    clin["project"] = clin["Project.ID"].str.replace("TCGA-", "", regex=False)
    return pd.concat([clin, parse_barcode(clin["Sample.ID"])], axis=1)


def select_samples(cfg: Config) -> tuple[pd.DataFrame, pd.DataFrame]:
    """D-01. Returns (selected samples, accounting of every clinical file)."""
    rules = cfg.study["cohort"]
    clin = load_clinical(cfg)
    ids = molecular_file_ids(cfg)
    clin["has_all_views"] = clin["File.ID"].apply(lambda f: all(f in s for s in ids.values()))

    role = pd.Series("", index=clin.index)
    role[clin["Sample.Type"].isin(rules["tumour_sample_types"])] = "tumour"
    role[clin["Sample.Type"] == rules["normal_sample_type"]] = "normal"
    clin["role"] = role

    acct = clin[["File.ID", "project", "Sample.Type", "role", "has_all_views"]].copy()
    acct["status"] = "selected"
    acct.loc[clin["role"] == "", "status"] = "sample type not used"
    acct.loc[(clin["role"] != "") & ~clin["has_all_views"], "status"] = "missing a view"

    cand = clin[(clin["role"] != "") & clin["has_all_views"]].copy()
    order = {"vial": "vial", "sample_id": "Sample.ID", "file_id": "File.ID"}
    cand = cand.sort_values(["Case.ID", "role", *[order[k] for k in rules["tie_break"]]], kind="mergesort")
    rank = cand.groupby(["Case.ID", "role"]).cumcount()
    dropped = cand.loc[rank > 0, "File.ID"]
    acct.loc[acct["File.ID"].isin(dropped), "status"] = "not first file of case and role"
    sel = cand[rank == 0]

    # a case belongs to one project; a solid normal must carry barcode type 11
    if sel.groupby("Case.ID")["project"].nunique().max() > 1:
        raise ValueError("a case maps to more than one project")
    if not (sel.loc[sel.role == "normal", "sample_type_code"] == "11").all():
        raise ValueError("Solid Tissue Normal rows without barcode type 11")

    out = (sel.rename(columns={"File.ID": "file_id", "Case.ID": "case_id", "Sample.ID": "sample_id",
                               "Sample.Type": "sample_type"})
           [["file_id", "case_id", "sample_id", "project", "sample_type", "role", "tss", "vial"]]
           .sort_values(["project", "role", "case_id"]).reset_index(drop=True))
    return out, acct.rename(columns={"File.ID": "file_id"})


def task_scope(cfg: Config, samples: pd.DataFrame) -> pd.DataFrame:
    """D-02: per-project tumour / normal case counts and tumor_status eligibility."""
    n_min = int(cfg.study["tasks"]["tumor_status"]["min_normal_cases"])
    g = samples.groupby(["project", "role"])["case_id"].nunique().unstack(fill_value=0)
    g = g.reindex(columns=["tumour", "normal"], fill_value=0).rename(columns={"tumour": "tumour_cases", "normal": "normal_cases"})
    paired = samples.groupby("case_id")["role"].nunique()
    g["paired_cases"] = samples[samples.case_id.isin(paired[paired == 2].index) & (samples.role == "normal")].groupby("project").size()
    g["paired_cases"] = g["paired_cases"].fillna(0).astype(int)
    g["tumor_status_in_scope"] = g["normal_cases"] >= n_min
    return g.reset_index()


def task_tables(cfg: Config, samples: pd.DataFrame, scope: pd.DataFrame) -> dict[tuple[str, str], pd.DataFrame]:
    """{(task, cohort): samples with an integer `y` and a `label` string}."""
    out = {}
    for proj in scope.loc[scope["tumor_status_in_scope"], "project"]:
        d = samples[samples["project"] == proj].copy()
        d["label"] = d["role"].map({"normal": "Normal", "tumour": "Tumour"})
        d["y"] = (d["role"] == "tumour").astype(int)
        out[("tumor_status", proj)] = d.reset_index(drop=True)

    normal_label = cfg.study["tasks"]["cancer_type"]["normal_label"]
    d = samples.copy()
    d["label"] = d["project"].where(d["role"] == "tumour", normal_label)
    classes = sorted(set(d["label"]) - {normal_label}) + [normal_label]
    counts = d.groupby("label")["case_id"].nunique()
    small = counts[counts < cfg.study["tasks"]["cancer_type"]["min_class_cases"]]
    if len(small):
        raise ValueError(f"cancer_type classes below min_class_cases: {small.to_dict()}")
    d["y"] = d["label"].map({c: i for i, c in enumerate(classes)}).astype(int)
    out[("cancer_type", PAN)] = d.reset_index(drop=True)
    return out


def class_names(table: pd.DataFrame) -> list[str]:
    """Class names in `y` order."""
    m = table.drop_duplicates("y").set_index("y")["label"].sort_index()
    if list(m.index) != list(range(len(m))):
        raise ValueError("labels are not 0..C-1")
    return list(m)
