"""Step 06: cross-job summary tables, figures and results/RESULTS.md.

Reads results/<task>/<view>/<cohort>/metrics.json (+ local SHAP arrays for the figures).
Writes results/summary/*.csv and figures (house style, docs/FIGURE_GUIDELINES.md: SVG, PDF,
400-dpi PNG, a CSV of the plotted values and a caption file each), per-job figures (PNG),
results/RESULTS.md, and the two results blocks of README.md (between the RESULTS-STATUS and
RESULTS-TYPE markers). Safe to run while jobs are still running: missing jobs are reported as such.
"""
from __future__ import annotations

import argparse
import json
import re

import numpy as np
import pandas as pd
from _common import REPO, log

from ctsm.config import PAN, load_config
from ctsm.explain import class_importance_matrix, load_folds
from ctsm.figstyle import feature_label
from ctsm.plots import VIEW_NAMES, beeswarm, dot_metrics, heatmap

N_BEESWARM = 15


def collect(cfg) -> tuple[pd.DataFrame, pd.DataFrame]:
    ts, ct = [], []
    # tumour status and cancer type only: results/subtype/ holds the subtype task's jobs (step 34)
    paths = [*cfg.results_dir.glob("tumor_status/*/*/metrics.json"), *cfg.results_dir.glob("cancer_type/*/*/metrics.json")]
    for p in sorted(paths):
        m = json.loads(p.read_text())
        row = {"view": m["view"], "cohort": m["cohort"], **m["pooled_oof"],
               **{f"{k}_lo": v[0] for k, v in m["bootstrap_95ci"].items()},
               **{f"{k}_hi": v[1] for k, v in m["bootstrap_95ci"].items()},
               **{f"fold_mean_{k}": v for k, v in m["per_fold_mean"].items()},
               **{f"fold_sd_{k}": v for k, v in m["per_fold_sd"].items()}}
        (ts if m["task"] == "tumor_status" else ct).append(row)
    return pd.DataFrame(ts), pd.DataFrame(ct)


def view_order(cfg, df: pd.DataFrame) -> list[str]:
    return [v for v in cfg.views if v in set(df["view"])]


def fmt_ci(v, lo, hi, nd=3) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "–"
    return f"{v:.{nd}f} ({lo:.{nd}f}–{hi:.{nd}f})" if lo is not None and not pd.isna(lo) else f"{v:.{nd}f}"


def per_job_figures(cfg, task: str, view: str, cohort: str) -> None:
    rdir, out = cfg.run_dir(task, view, cohort), cfg.result_dir(task, view, cohort)
    tab = pd.read_parquet(cfg.data_dir / "tasks" / f"{task}__{cohort}.parquet")
    folds = load_folds(rdir, tab["y"].to_numpy(), cfg.study["splits"]["outer_folds"])
    feats = pd.read_csv(out / "shap_top_features.csv")
    vname = VIEW_NAMES[view]
    if task == "tumor_status":
        top = feats.head(N_BEESWARM)
        shs, vals = [], []
        for fid in top["feature_id"]:
            s, v = [], []
            for f in folds:
                j = f["index"].get(fid)
                if j is not None:
                    s.append(f["shap"][:, j])
                    v.append(f["X"][:, j])
            shs.append(np.concatenate(s))
            vals.append(np.concatenate(v))
        labels = [f"{feature_label(lab)}  ({n}/5)" for lab, n in zip(top["label"], top["n_folds_selected"])]
        beeswarm(shs, vals, labels, f"SHAP, {cohort} tumour log-odds", out, "fig_shap_beeswarm",
                 caption=f"**{cohort}, {vname}, tumour vs normal: the {len(top)} features with the largest mean "
                         "|SHAP|.** Each dot is one held-out sample (every outer fold explains only its own test "
                         "samples); x = SHAP value on the tumour log-odds (> 0 pushes towards tumour); colour = "
                         "the feature's standardised value (blue low, red high, clipped at ±2.5). Row label: "
                         "host gene (HGNC) and event coordinates, then the number of the 5 outer folds whose "
                         "feature selection kept the feature. Values: `shap_top_features.csv` (per-sample SHAP "
                         "values stay local, D-12).")
        return
    classes = json.loads((out / "metrics.json").read_text())["classes"]
    imp = class_importance_matrix(folds, classes)
    lab = pd.read_csv(out / "shap_top_features_by_class.csv")[["feature_id", "label"]].drop_duplicates("feature_id")
    lab = lab.set_index("feature_id")["label"]
    rows = []
    for c in classes:                                     # top 2 per class, ordered by class
        for fid in imp[c].sort_values(ascending=False).index[:2]:
            if fid not in rows:
                rows.append(fid)
    mat = imp.loc[rows].copy()
    mat.index = [feature_label(lab.get(f, f)) for f in rows]
    heatmap(mat, out, "fig_shap_class_heatmap", vmin=0, fmt=None, cbar_label="mean |SHAP|", cell=0.11, xrot=90,
            formats=("png",), label_size=4.5,
            caption=f"**Cancer type, {vname}: what each class relies on.** Rows: the two features with the largest "
                    "mean |SHAP| for each class (host gene and coordinates); columns: classes. Cell: mean |SHAP| "
                    "of the row's feature over the held-out samples of the column's class, explaining each "
                    "sample's own centred class logit (D-10). A block on the diagonal means the class's top "
                    "features are specific to it.")
    cm = pd.read_csv(out / "confusion_matrix.csv", index_col=0)
    cmn = cm.div(cm.sum(axis=1).replace(0, np.nan), axis=0)
    heatmap(cmn, out, "fig_confusion_matrix", vmin=0, vmax=1, fmt=None, cbar_label="fraction of true class",
            cell=0.11, xrot=90, formats=("png",), label_size=4.5,
            caption=f"**Cancer type, {vname}: out-of-fold confusion matrix, row-normalised** (rows = true class, "
                    "columns = predicted class). Counts: `confusion_matrix.csv`.")


def top_genes(cfg, task: str, cohorts: list[str], k: int = 5) -> pd.DataFrame:
    rows = []
    for view in cfg.views:
        for c in cohorts:
            p = cfg.result_dir(task, view, c) / "shap_top_genes.csv"
            if p.exists():
                g = pd.read_csv(p).head(k)
                rows.append({"cohort": c, "view": view, "top_genes": ", ".join(g["symbol"].astype(str))})
    return pd.DataFrame(rows)


def readme_blocks(cfg, ts: pd.DataFrame, ct: pd.DataFrame) -> dict[str, list[str]]:
    """The two headline tables of README.md, one per task: the tasks are written up separately,
    so their results are posted separately too."""
    blocks = {}
    if len(ts):
        b = ["**Out-of-fold, mean over the 15 per-cancer models:**", "",
             "| View | AUROC | Balanced accuracy | MCC |", "|---|---|---|---|"]
        for v in view_order(cfg, ts):
            x = ts[ts.view == v]
            b.append(f"| {VIEW_NAMES[v]} | {x.auroc.mean():.3f} | {x.balanced_accuracy.mean():.3f} | {x.mcc.mean():.3f} |")
        blocks["STATUS"] = b
    if len(ct):
        b = ["**Out-of-fold, one 34-class model per view:**", "",
             "| View | Accuracy | Balanced accuracy | Macro F1 |", "|---|---|---|---|"]
        for v in view_order(cfg, ct):
            r = ct[ct.view == v].iloc[0]
            b.append(f"| {VIEW_NAMES[v]} | {r.accuracy:.3f} | {r.balanced_accuracy:.3f} | {r.macro_f1:.3f} |")
        blocks["TYPE"] = b
    return blocks


def write_readme_blocks(cfg, ts: pd.DataFrame, ct: pd.DataFrame) -> None:
    """Replace the text between <!-- RESULTS-<TASK>:START --> and <!-- RESULTS-<TASK>:END -->."""
    readme = REPO / "README.md"
    if not readme.exists():
        return
    txt = readme.read_text()
    for key, lines in readme_blocks(cfg, ts, ct).items():
        head = f"<!-- RESULTS-{key}:START (generated by scripts/06_summarise_results.py) -->"
        txt = re.sub(rf"<!-- RESULTS-{key}:START.*?<!-- RESULTS-{key}:END -->",
                     lambda _: "\n".join([head, *lines, "", f"<!-- RESULTS-{key}:END -->"]), txt, flags=re.S)
    readme.write_text(txt)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-job-figures", action="store_true", help="skip the per-job figures (summary only)")
    args = ap.parse_args()
    cfg = load_config()
    out = cfg.results_dir / "summary"
    out.mkdir(parents=True, exist_ok=True)
    ts, ct = collect(cfg)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = list(scope.loc[scope.tumor_status_in_scope, "project"])
    n_expected = len(cfg.views) * (len(cancers) + 1)
    log(f"{len(ts) + len(ct)} of {n_expected} jobs have results")

    md = ["# Results", "",
          "Generated by `scripts/06_summarise_results.py`. All numbers are **out-of-fold**: each sample is "
          "predicted by the outer-fold model that never saw it (or any sample of its case) during feature "
          "selection, tuning or training. Brackets: 95 % case-level bootstrap CIs on the pooled "
          "out-of-fold predictions.", ""]

    # ---- task 1: tumour status -------------------------------------------------------------
    if len(ts):
        ts.to_csv(out / "tumor_status_metrics.csv", index=False)
        views = view_order(cfg, ts)
        names = [VIEW_NAMES[v] for v in views]
        auc = ts.pivot(index="cohort", columns="view", values="auroc").reindex(index=cancers, columns=views)
        auc.columns = names
        auc.index.name = "cancer"
        heatmap(auc, out, "fig_tumor_status_auroc", vmin=0.85, vmax=1.0, extend="min", cbar_label="AUROC", cell=0.27,
                caption="**How well does each view separate tumour from matched-site normal tissue?** Pooled "
                        "out-of-fold AUROC of each per-cancer tumour-vs-normal model (rows: the 15 cancers with "
                        ">= 10 normal cases; columns: views). Each sample is scored by the outer-fold model that "
                        "never saw its patient. CIs and balanced accuracy: `results/RESULTS.md`.")
        md += ["## Task 1 — tumour vs normal, per cancer", "",
               "![AUROC](summary/fig_tumor_status_auroc.png)", "",
               "### AUROC (95 % CI)", "", "| Cancer | n (normal) | " + " | ".join(names) + " |",
               "|---|---|" + "---|" * len(names)]
        for c in cancers:
            r = scope.set_index("project").loc[c]
            cells = []
            for v in views:
                x = ts[(ts.cohort == c) & (ts.view == v)]
                cells.append(fmt_ci(x.auroc.iloc[0], x.auroc_lo.iloc[0], x.auroc_hi.iloc[0]) if len(x) else "–")
            md.append(f"| {c} | {r.tumour_cases + r.normal_cases} ({r.normal_cases}) | " + " | ".join(cells) + " |")
        md += ["", "### Balanced accuracy / specificity for normals (MCC)", "",
               "| Cancer | " + " | ".join(names) + " |", "|---|" + "---|" * len(names)]
        for c in cancers:
            cells = []
            for v in views:
                x = ts[(ts.cohort == c) & (ts.view == v)]
                cells.append(f"{x.balanced_accuracy.iloc[0]:.3f} / {x.specificity_normal.iloc[0]:.2f} "
                             f"({x.mcc.iloc[0]:.2f})" if len(x) else "–")
            md.append(f"| {c} | " + " | ".join(cells) + " |")
        md += ["", "Mean over cancers:", "", "| Metric | " + " | ".join(names) + " |", "|---|" + "---|" * len(names)]
        for m in ["auroc", "balanced_accuracy", "auprc_normal", "mcc"]:
            piv = ts.pivot(index="cohort", columns="view", values=m)
            md.append(f"| {m} | " + " | ".join(f"{piv[v].mean():.3f}" if v in piv else "–" for v in views) + " |")
        md.append("")

    # ---- task 2: cancer type -----------------------------------------------------------------
    if len(ct):
        ct.to_csv(out / "cancer_type_metrics.csv", index=False)
        views = view_order(cfg, ct)
        d = ct.set_index("view").loc[views]
        d.index = [VIEW_NAMES[v] for v in views]
        dot_metrics(d, ["accuracy", "balanced_accuracy", "macro_f1"], ["accuracy", "balanced accuracy", "macro F1"],
                    out, "fig_cancer_type_metrics",
                    caption="**How well does each view identify the cancer type (33 TCGA projects + pooled "
                            "normal, 34 classes)?** Out-of-fold accuracy, balanced accuracy and macro F1, with "
                            "95% case-level bootstrap CIs (200 resamples).")
        rec = {}
        for v in views:
            pc = pd.read_csv(cfg.result_dir("cancer_type", v, PAN) / "per_class_metrics.csv")
            rec[VIEW_NAMES[v]] = pc.set_index("class")["recall"]
        rec = pd.DataFrame(rec)
        rec.to_csv(out / "cancer_type_recall_by_class.csv")
        rec.index.name = "class"
        heatmap(rec, out, "fig_cancer_type_recall", vmin=0.5, vmax=1, extend="min", fmt="{:.2f}", cbar_label="recall",
                cell=0.24,
                caption="**Which classes does each view recover?** Out-of-fold recall per class (rows) and view "
                        "(columns).")
        md += ["## Task 2 — cancer type (33 TCGA projects + pooled Normal), pan-cancer", "",
               "![Cancer type metrics](summary/fig_cancer_type_metrics.png)", "",
               "| View | Accuracy | Balanced accuracy | Macro F1 | Macro AUROC (OvR) | Top-3 accuracy | Log loss |",
               "|---|---|---|---|---|---|---|"]
        for v in views:
            r = ct[ct.view == v].iloc[0]
            md.append(f"| {VIEW_NAMES[v]} | {fmt_ci(r.accuracy, r.accuracy_lo, r.accuracy_hi)} | "
                      f"{fmt_ci(r.balanced_accuracy, r.balanced_accuracy_lo, r.balanced_accuracy_hi)} | "
                      f"{fmt_ci(r.macro_f1, r.macro_f1_lo, r.macro_f1_hi)} | {r.macro_auroc_ovr:.3f} | "
                      f"{r.top3_accuracy:.3f} | {r.log_loss:.3f} |")
        md += ["", "![Per-class recall](summary/fig_cancer_type_recall.png)", ""]

    # ---- interpretability -------------------------------------------------------------------
    md += ["## Interpretability (SHAP)", "",
           "Every outer-fold model is explained with `shap.DeepExplainer` on **all** of its held-out samples. "
           "Per job folder (`results/<task>/<view>/<cohort>/`):", "",
           "- `shap_top_features.csv`: top 200 features by mean |SHAP| pooled over folds (a feature a fold did "
           "not select contributes 0 for that fold), with `n_folds_selected` and `direction_spearman` "
           "(Spearman of feature value vs SHAP: > 0 means higher PSI / usage / expression pushes towards tumour);",
           "- `shap_top_genes.csv`: the same importances summed per host gene;",
           "- `shap_top_features_by_class.csv` (cancer type): top 30 features per class, explaining each "
           "sample's own centred class logit, with a per-class direction;",
           "- figures: `fig_shap_beeswarm.png` (tumour status), `fig_shap_class_heatmap.png` and "
           "`fig_confusion_matrix.png` (cancer type), each with a `.caption.md`.", ""]
    for task, cohorts in [("tumor_status", cancers), ("cancer_type", [PAN])]:
        tg = top_genes(cfg, task, cohorts)
        if len(tg):
            tg.to_csv(out / f"{task}_top_genes.csv", index=False)
            piv = tg.pivot(index="cohort", columns="view", values="top_genes").reindex(index=cohorts)
            piv = piv[[v for v in cfg.views if v in piv.columns]]
            md += [f"### Top 5 SHAP genes · {task}", "", "| Cohort | " + " | ".join(VIEW_NAMES[v] for v in piv.columns) + " |",
                   "|---|" + "---|" * len(piv.columns)]
            for c, r in piv.iterrows():
                md.append(f"| {c} | " + " | ".join("–" if pd.isna(x) else x for x in r) + " |")
            md.append("")

    # ---- per-job figures and links -----------------------------------------------------------
    done = [(t, v, c) for t, df in [("tumor_status", ts), ("cancer_type", ct)] if len(df)
            for v, c in zip(df["view"], df["cohort"])]
    for t, v, c in ([] if args.no_job_figures else done):
        try:
            per_job_figures(cfg, t, v, c)
        except FileNotFoundError as e:
            log(f"skip figures for {t}/{v}/{c}: {e}")
    md += ["## Per-job folders", "", "| Task | View | Cohort | Folder |", "|---|---|---|---|"]
    md += [f"| {t} | {VIEW_NAMES[v]} | {c} | [{t}/{v}/{c}]({t}/{v}/{c}/) |" for t, v, c in sorted(done)]
    (cfg.results_dir / "RESULTS.md").write_text("\n".join(md) + "\n")

    write_readme_blocks(cfg, ts, ct)
    log("summary written")


if __name__ == "__main__":
    main()
