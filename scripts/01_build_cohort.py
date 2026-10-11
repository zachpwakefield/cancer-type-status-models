"""Step 01: samples, task tables and nested CV splits (D-01..D-04).

Outputs (<data_dir>/): cohort/samples.parquet, tasks/<task>__<cohort>.parquet,
splits/<task>__<cohort>.parquet
Summaries (committed, results/data/): sample_accounting.csv, task_scope.csv, class_counts.csv,
split_summary.csv
"""
from __future__ import annotations

import pandas as pd
from _common import log

from ctsm.cohort import select_samples, task_scope, task_tables
from ctsm.config import load_config
from ctsm.provenance import run_record
from ctsm.splits import fold_indices, make_splits, min_inner_train


def main() -> None:
    cfg = load_config()
    sp = cfg.study["splits"]
    out_dir, rep = cfg.data_dir, cfg.results_dir / "data"
    for d in (out_dir / "cohort", out_dir / "tasks", out_dir / "splits", rep):
        d.mkdir(parents=True, exist_ok=True)

    samples, acct = select_samples(cfg)
    samples.to_parquet(out_dir / "cohort" / "samples.parquet", index=False)
    (acct.groupby(["Sample.Type", "status"]).size().rename("files").reset_index()
     .to_csv(rep / "sample_accounting.csv", index=False))
    log(f"{len(samples):,} samples selected from {samples.case_id.nunique():,} cases "
        f"({(samples.role == 'tumour').sum():,} tumour, {(samples.role == 'normal').sum():,} normal)")

    scope = task_scope(cfg, samples)
    scope.to_csv(rep / "task_scope.csv", index=False)
    log(f"tumor_status in scope: {', '.join(scope.loc[scope.tumor_status_in_scope, 'project'])}")

    written, split_rows, class_rows = [out_dir / "cohort" / "samples.parquet"], [], []
    for (task, cohort), tab in task_tables(cfg, samples, scope).items():
        p = out_dir / "tasks" / f"{task}__{cohort}.parquet"
        tab.to_parquet(p, index=False)
        s = make_splits(tab, task, cohort, sp["outer_folds"], sp["inner_folds"], sp["seed"])
        q = out_dir / "splits" / f"{task}__{cohort}.parquet"
        s.to_parquet(q, index=False)
        written += [p, q]
        for lab, g in tab.groupby("label"):
            class_rows.append({"task": task, "cohort": cohort, "label": lab, "samples": len(g), "cases": g.case_id.nunique()})
        for k in range(sp["outer_folds"]):
            _, te = fold_indices(s, k)
            yte = s["y"].to_numpy()[te]
            split_rows.append({"task": task, "cohort": cohort, "outer_fold": k, "n_test": len(te),
                               "n_test_minority": int(pd.Series(yte).value_counts().min()),
                               "n_classes_in_test": int(len(set(yte)))})
        log(f"{task}/{cohort}: {len(tab):,} samples, {tab.y.nunique()} classes, "
            f"min inner-train {min_inner_train(s, sp['outer_folds'], sp['inner_folds'])}")
    pd.DataFrame(class_rows).to_csv(rep / "class_counts.csv", index=False)
    pd.DataFrame(split_rows).to_csv(rep / "split_summary.csv", index=False)
    run_record(cfg, "01_build_cohort", [cfg.tcga_path("clinical")], written)


if __name__ == "__main__":
    main()
