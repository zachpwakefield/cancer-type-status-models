"""Step 30: subtype labels, task tables and nested CV splits (docs/SUBTYPE_METHODS.md, D-30..D-36, D-43).

Reads the D-01 cohort (<data_dir>/cohort/samples.parquet, written by step 01) and the label
sources named in config/subtype.yaml; never a molecular value.

Outputs, local (<data_dir>/subtype/, git-ignored, per patient):
    labels/labels.parquet                 every primary tumour of every cohort: raw values, classes, reason
    tasks/subtype__<cohort>.parquet       labelled patients of each cohort in scope (label, y)
    splits/subtype__<cohort>.parquet      5 outer x 3 inner folds, grouped by patient, stratified by label
Summaries, committed (results/subtype/data/, counts only):
    cohort_scope.csv        per cohort: tier, label type, classes and counts, classes below the minimum
    label_accounting.csv    per cohort, source and raw value: patients, class, patients labelled
    exclusions.csv          per cohort: patients by reason
    class_counts.csv        per cohort and class: patients
    source_concordance.csv  histology: TCGA-CDR class x GDC class, patients
    label_crosstab.csv      two cohorts of one project: class x class, patients labelled in both
    split_summary.csv       per cohort and outer fold: test size and smallest class
    label_sources.json      SHA-256 of every label file read
"""
from __future__ import annotations

import itertools

import pandas as pd
from _common import log

from ctsm.config import load_config
from ctsm.provenance import run_record, sha256_file, write_json
from ctsm.splits import fold_indices, make_splits, min_inner_train
from ctsm.subtype import TASK, _pancan_dir, build_all, cohort_specs, subtype_config


def main() -> None:
    base = load_config()
    scfg = subtype_config(base)
    sp = scfg.study["splits"]
    out, rep = scfg.data_dir, scfg.results_dir / TASK / "data"
    for d in (out / "labels", out / "tasks", out / "splits", rep):
        d.mkdir(parents=True, exist_ok=True)
    samples = pd.read_parquet(base.data_dir / "cohort" / "samples.parquet")
    tables, labels, scope, acct = build_all(scfg, scfg.spec, samples)
    labels.to_parquet(out / "labels" / "labels.parquet", index=False)
    scope.to_csv(rep / "cohort_scope.csv", index=False)
    acct.to_csv(rep / "label_accounting.csv", index=False)
    (labels.groupby(["cohort", "reason"]).size().rename("patients").reset_index()
     .to_csv(rep / "exclusions.csv", index=False))
    for _, r in scope.iterrows():
        log(f"{r.cohort:<11} {r.tier:<18} {'IN ' if r.in_scope else 'OUT'} {r.labelled:>4}/{r.primary_tumours:<4} "
            f"{r.classes}" + (f"  | below minimum: {r.classes_below_minimum}" if r.classes_below_minimum else ""))

    # histology: how often the two sources agree
    specs = cohort_specs(scfg.spec)
    conc = []
    for cohort, cs in specs.items():
        if not cs["concordance"]:
            continue
        lab = labels[labels.cohort == cohort]
        g = lab.groupby([lab["class_primary"].fillna("(no class)"), lab["class_concordance"].fillna("(no class)")])
        conc += [{"cohort": cohort, "class_cdr": a, "class_gdc": b, "patients": len(x)} for (a, b), x in g]
    pd.DataFrame(conc).to_csv(rep / "source_concordance.csv", index=False)

    # two label families of one project, on the patients labelled in both
    cross = []
    by_project = {}
    for cohort, cs in specs.items():
        by_project.setdefault(cs["project"], []).append(cohort)
    for project, cohorts in by_project.items():
        for a, b in itertools.combinations(cohorts, 2):
            la = labels[(labels.cohort == a) & labels.label.notna()].set_index("case_id")["label"]
            lb = labels[(labels.cohort == b) & labels.label.notna()].set_index("case_id")["label"]
            both = la.index.intersection(lb.index)
            for (x, y), n in pd.crosstab(la[both], lb[both]).stack().items():
                cross.append({"project": project, "cohort_a": a, "cohort_b": b, "class_a": x, "class_b": y,
                              "patients": int(n)})
    pd.DataFrame(cross).to_csv(rep / "label_crosstab.csv", index=False)

    written, split_rows, class_rows = [out / "labels" / "labels.parquet"], [], []
    for cohort, tab in tables.items():
        p = out / "tasks" / f"{TASK}__{cohort}.parquet"
        tab.to_parquet(p, index=False)
        s = make_splits(tab, TASK, cohort, sp["outer_folds"], sp["inner_folds"], sp["seed"])
        q = out / "splits" / f"{TASK}__{cohort}.parquet"
        s.to_parquet(q, index=False)
        written += [p, q]
        for lab, g in tab.groupby("label", sort=False):
            class_rows.append({"cohort": cohort, "label": lab, "y": int(g.y.iloc[0]), "cases": g.case_id.nunique()})
        for k in range(sp["outer_folds"]):
            _, te = fold_indices(s, k)
            yte = s["y"].to_numpy()[te]
            split_rows.append({"cohort": cohort, "outer_fold": k, "n_test": len(te),
                               "n_test_smallest_class": int(pd.Series(yte).value_counts().min()),
                               "n_classes_in_test": int(len(set(yte)))})
        log(f"{cohort}: {len(tab):,} patients, {tab.y.nunique()} classes, "
            f"min inner-train {min_inner_train(s, sp['outer_folds'], sp['inner_folds'])}")
    pd.DataFrame(class_rows).sort_values(["cohort", "y"]).to_csv(rep / "class_counts.csv", index=False)
    pd.DataFrame(split_rows).to_csv(rep / "split_summary.csv", index=False)

    files = {"tcga_cdr": scfg.sources["subtype_labels"]["tcga_cdr"], "gdc_clinical": str(base.tcga_path("clinical")),
             **{n: str(_pancan_dir(base) / n) for n in ("TCGASubtype.20170308.tsv",
                                                         "clinical_PANCAN_patient_with_followup.tsv")}}
    write_json(rep / "label_sources.json", {k: {"file": v.split("/")[-1], "sha256": sha256_file(v)}
                                            for k, v in files.items()})
    run_record(scfg, "30_subtype_labels", [base.data_dir / "cohort" / "samples.parquet", *files.values()], written)


if __name__ == "__main__":
    main()
