"""Step 17 (B-08): run splice-assay over the step-16 dataset and collect the per-event statistics.

For every cancer, endpoint (OS, DSS) and model arm, `splice-assay analyze` fits, per event:
    group tests   paired (matched tumour-normal) and unpaired tumour vs normal;
    KM            median split, log-rank;
    Cox           base: PSI + host-gene expression
                  adjusted: the base model + age + sex + stage (stage I, female as baselines).
The host-expression term is what makes this more than an expression echo: an event's hazard
ratio is estimated with its own gene's expression in the model.

The candidate / background design of step 16 is what the enrichment test rests on. Both sets are
analysed identically, in the same cohorts, with the same gates, so the comparison of their hit
rates is not a comparison of measurability.

Outputs: results/biology/08_survival/ (aggregate tables; figures: scripts/20_figures_status.py).
Per-patient values stay in <data_dir>/survival/ (D-12).
"""
from __future__ import annotations

import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from PIL import Image
from _common import REPO, log
from scipy.stats import fisher_exact

from ctsm.biology import bh
from ctsm.config import load_config

SA_PYTHON = "/Users/zacharywakefield/miniconda3/envs/first_env/bin/python"   # splice-assay 0.4.0 lives here
ENDPOINTS = ("OS", "DSS")
ARMS = {"base": [], "adjusted": ["--covariate", "age", "--covariate", "sex", "--covariate", "stage",
                                 "--baseline", "stage=I", "--baseline", "sex=female"]}


def run_one(job) -> str:
    data, out, cohort, endpoint, arm = job
    done = out / "survival.csv"
    if done.exists():
        return f"{cohort} {endpoint} {arm}: cached"
    out.mkdir(parents=True, exist_ok=True)
    cmd = [SA_PYTHON, "-m", "splice_assay.cli", "analyze", str(data), "--out", str(out),
           "--cohort", cohort, "--endpoint", endpoint, "--include-hit", "--na-value", "is_missing",
           "--case", "tumour", "--reference", "normal"] + ARMS[arm]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return f"{cohort} {endpoint} {arm}: FAILED\n{r.stderr[-800:]}"
    return f"{cohort} {endpoint} {arm}: ok"


def permuted_inputs(cfg, seed: int) -> Path:
    """A copy of the dataset with the survival columns shuffled between patients WITHIN each
    cohort. Everything else — PSI, expression, clinical, group — is untouched, so a run over it
    measures what this pipeline reports when no event is associated with outcome."""
    src = cfg.data_dir / "survival" / "inputs"
    dst = cfg.data_dir / "survival" / "inputs_null"
    dst.mkdir(parents=True, exist_ok=True)
    for name in ("psi.parquet", "expression.parquet"):            # unchanged: link, do not copy
        link = dst / name
        if not link.exists():
            link.symlink_to(src / name)
    rng = np.random.default_rng(seed)
    smp = pd.read_parquet(src / "samples.parquet")
    cols = [c for c in smp.columns if c.split(".")[0] in ("OS", "DSS")]
    for cohort, idx in smp.groupby("cohort").groups.items():
        rows = np.asarray(idx)
        smp.loc[rows, cols] = smp.loc[rows[rng.permutation(len(rows))], cols].to_numpy()
    smp.to_parquet(dst / "samples.parquet", index=False)
    return dst


def collect(root, cohorts) -> tuple[pd.DataFrame, pd.DataFrame]:
    surv, groups = [], []
    for cohort in cohorts:
        for endpoint in ENDPOINTS:
            for arm in ARMS:
                d = root / f"{cohort}__{endpoint}__{arm}"
                if not (d / "survival.csv").exists():
                    continue
                s = pd.read_csv(d / "survival.csv", low_memory=False).assign(arm=arm)
                surv.append(s)
                if endpoint == "OS" and arm == "base":     # group tests do not depend on either
                    groups.append(pd.read_csv(d / "group_tests.csv", low_memory=False))
    return pd.concat(surv, ignore_index=True), pd.concat(groups, ignore_index=True)


def enrichment(cells: pd.DataFrame, sel: pd.DataFrame) -> pd.DataFrame:
    """Share of testable event x cohort cells with p < 0.05, against the matched background.

    Three nested candidate sets, all chosen without seeing an outcome: every candidate, those a
    model ranks in its top 10 in more than one cancer, and the recurrent events of D-17 (top 100
    of at least 4 models)."""
    rows = []
    sets = {"all candidates": cells["role"].eq("candidate"),
            "top 10 in >= 2 cancers": cells["role"].eq("candidate") & cells["n_models"].fillna(0).ge(2),
            "recurrent (D-17)": cells["role"].eq("candidate") & cells["recurrent"].fillna(False)}
    bg = cells["role"].eq("background")
    for (endpoint, arm), part in cells.groupby(["endpoint", "arm"]):
        for test, pcol, gate in [("Cox", "cox_p", "cox_ok"), ("log-rank", "km_p", "km_ok")]:
            ok = part[gate]
            bgp = bg.loc[part.index] & ok               # explicit: bg is indexed on all cells
            b_n = int(bgp.sum())
            b_h = int((part.loc[bgp, pcol] < 0.05).sum())
            for name, mask in sets.items():
                sub = part[mask.loc[part.index] & ok]
                if not len(sub) or not b_n:
                    continue
                c_n, c_h = len(sub), int((sub[pcol] < 0.05).sum())
                odds, p = fisher_exact([[c_h, c_n - c_h], [b_h, b_n - b_h]])
                rows.append({"endpoint": endpoint, "arm": arm, "test": test, "candidate_set": name,
                             "n_candidate": c_n, "hits_candidate": c_h, "rate_candidate": c_h / c_n,
                             "n_background": b_n, "hits_background": b_h, "rate_background": b_h / b_n,
                             "odds_ratio": odds, "p_fisher": p})
    out = pd.DataFrame(rows)
    out["q_fisher"] = bh(out["p_fisher"].to_numpy())
    return out


def dossiers(cfg, cells: pd.DataFrame, sel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One row per candidate event, joining what every other step says about it: how many models
    rank it, how it moves between tumour and normal, what it does to the protein, whose regulon
    it is in, and its strongest survival association."""
    bio = cfg.results_dir / "biology"
    ev = cells[cells["role"].eq("candidate") & cells["endpoint"].eq("OS") & cells["arm"].eq("adjusted")]
    best = ev[ev["cox_ok"]].sort_values("cox_p").drop_duplicates("event_id")
    agg = (ev[ev["cox_ok"]].groupby("event_id")
           .agg(cancers_tested=("cohort", "size"), hits=("cox_p", lambda x: int((x < 0.05).sum()))))
    d = sel[sel["role"].eq("candidate")].drop_duplicates("event_id").copy()
    d = d.merge(agg, on="event_id", how="left")
    d = d.merge(best[["event_id", "gene", "cohort", "hr_per_sd", "ci_low_sd", "ci_high_sd", "cox_p", "cox_q",
                      "cox_events", "cox_low_power"]].rename(columns={"cohort": "best_cancer_survival"}),
                on="event_id", how="left")
    rec = pd.read_csv(bio / "04_events" / "recurrent_events.csv")
    rec["event_id"] = rec["view"].str.upper() + ":" + rec["feature_id"]
    d = d.merge(rec[["event_id", "n_cancers", "majority_sign", "median_delta_across_cancers",
                     "n_cancers_same_sign_as_majority", "n_cancers_measured"]]
                .rename(columns={"n_cancers": "n_models_top100"}), on="event_id", how="left")
    cons = pd.read_csv(bio / "04_events" / "cassette_exon_consequences.csv")
    if "feature_id" in cons.columns:
        cons["event_id"] = "SE:" + cons["feature_id"]
        keep = [c for c in ("event_id", "length", "len_mult3", "frame_preserving", "domains", "in_domain",
                            "in_disordered", "in_transmembrane") if c in cons]
        d = d.merge(cons[keep].drop_duplicates("event_id"), on="event_id", how="left")
    core = pd.read_csv(bio / "10_regulons" / "regulon_core_events.csv")
    core["event_id"] = core["view"].str.upper() + ":" + core["feature_id"]
    fe = pd.read_csv(bio / "10_regulons" / "regulon_enrichment_by_factor.csv")
    enriched = set(fe.loc[(fe["q_fisher"] < 0.05) & (fe["odds_ratio_haldane"] > 1), "factor"])
    reg = (core.groupby("event_id")
           .agg(n_factors=("factor", "nunique"),
                factors=("factor", lambda f: ",".join(sorted(set(f))[:6])),
                specific_factors=("factor", lambda f: ",".join(sorted(set(f) & enriched)))))
    d = d.merge(reg, on="event_id", how="left")
    d["n_factors"] = d["n_factors"].fillna(0).astype(int)
    d["in_specific_regulon"] = d["specific_factors"].fillna("").str.len() > 0
    order = d.sort_values(["hits", "n_models_top100", "cox_p"], ascending=[False, False, True])

    # does regulon membership go with outcome information?
    rows = []
    sub = cells[cells["role"].eq("candidate") & cells["endpoint"].eq("OS") & cells["arm"].eq("adjusted")
                & cells["cox_ok"]].merge(d[["event_id", "in_specific_regulon", "n_factors"]], on="event_id")
    # in-regulon first, so the odds ratio reads "in a regulon against not", never its reciprocal
    for group in (True, False):
        part = sub[sub["in_specific_regulon"].eq(group)]
        k, n = int((part["cox_p"] < 0.05).sum()), len(part)
        rows.append({"set": "in a model-enriched regulon" if group else "not in one",
                     "cells": n, "hits": k, "rate": k / n})
    comp = pd.DataFrame(rows)
    if len(comp) == 2:
        a, b = comp.iloc[0], comp.iloc[1]
        odds, p = fisher_exact([[a.hits, a.cells - a.hits], [b.hits, b.cells - b.hits]])
        comp["odds_ratio"], comp["p_fisher"] = odds, p
    return order, comp


def assay_pages(cfg, dos: pd.DataFrame, n: int, gtf: Path) -> list[str]:
    """splice-assay assay pages (tumour vs normal, Kaplan-Meier, the Cox model, the forest over
    cancers) for the events associated in the most cancers. The pages keep the tool's own
    layout, not the house style, because they come from a separately tested package.

    The pages' CSVs hold per-patient values, so they stay in <data_dir> (D-12); only the images
    are copied into results/."""
    local = cfg.data_dir / "survival" / "panels"
    pub = cfg.results_dir / "figures" / "tumour_status" / "assay"
    local.mkdir(parents=True, exist_ok=True)
    pub.mkdir(parents=True, exist_ok=True)
    made = []
    for r in dos.head(n).itertuples():
        stem = f"{r.gene}_{r.view.upper()}"
        cmd = [SA_PYTHON, "-m", "splice_assay.cli", "panel", str(cfg.data_dir / "survival" / "inputs"),
               "--event", r.event_id, "--endpoint", "OS", "--top", "3", "--stem", stem,
               "--na-value", "is_missing", "--case", "tumour", "--reference", "normal",
               "--detail-covariate", "age", "--detail-covariate", "sex", "--detail-covariate", "stage",
               "--baseline", "stage=I", "--baseline", "sex=female", "--out", str(local)]
        if gtf.exists():
            cmd += ["--gtf", str(gtf)]
        out = subprocess.run(cmd, capture_output=True, text=True)
        if out.returncode != 0:
            log(f"panel {stem}: FAILED {out.stderr.strip()[-300:]}")
            continue
        for img in sorted(local.glob(f"{stem}*")):
            # PNG to read in the browser and PDF to print; the SVG duplicates the PDF and the
            # pages are large, so it stays in <data_dir>/survival/panels/ with the per-patient CSVs.
            if img.suffix == ".pdf":
                (pub / img.name).write_bytes(img.read_bytes())
            elif img.suffix == ".png":
                # a page is 7 x 15 inches; at the tool's 400 dpi that is ~2 MB each. 150 dpi is
                # still sharp on screen and a third of the size, and the PDF stays exact.
                full = Image.open(img).convert("RGB")
                small = full.resize((full.width * 150 // 400, full.height * 150 // 400), Image.LANCZOS)
                small.save(pub / img.name, dpi=(150, 150), optimize=True)
        made.append(stem)
        log(f"panel {stem}: ok")
    if made:
        (pub / "README.md").write_text(
            "# Assay pages\n\n"
            "Drawn by splice-assay 0.4.0, not by this repository, so they keep the tool's layout rather than the\n"
            "house style (docs/FIGURE_GUIDELINES.md). Each page shows, for one event in its three most promising\n"
            "cancers: the tumour-vs-normal comparison (paired and unpaired), a Kaplan-Meier split at the median with\n"
            "its log-rank test, the Cox model (PSI + host-gene expression + age + sex + stage, against stage I and\n"
            "female), and a forest of every cancer. Markers: filled = p < 0.05; * = q < 0.05; "
            + chr(8224) + " = the proportional-hazards test failed; " + chr(8225) + " = fewer than 20 deaths "
            "(low power).\n\nThe star is the tool's own FDR family: BH within one event over the cancers it was "
            "tested in. That is not the family used in results/TUMOUR_STATUS.md, which asks the harder question "
            "-- BH over every event tested inside one cancer (`cox_q_cohort` in "
            "`biology/08_survival/event_survival_cells.csv.gz`). An event can be starred here and not clear the "
            "other.\n\nThe events shown are the candidates associated with overall survival in the most\n"
            "cancers (`biology/08_survival/event_dossiers.csv`). Nominal p values from observational cohorts: these\n"
            "are associations, not causes. The pages' CSVs hold per-patient values and stay local (D-12).\n"
            "The PNGs are the tool's pages at 150 dpi so the repository stays small; the PDF beside each one\n"
            "is the exact figure, and the SVG stays local too.\n\n"
            + "\n".join(f"- `{m}`" for m in made) + "\n")
    return made


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--collect-only", action="store_true")
    ap.add_argument("--panels", type=int, default=4, help="assay pages for the top N events (0 = none)")
    ap.add_argument("--null", action="store_true",
                    help="also run one permuted-outcome replicate (OS, base model) for the null hit rate")
    args = ap.parse_args()
    cfg = load_config()
    data = cfg.data_dir / "survival" / "inputs"
    root = cfg.data_dir / "survival" / "runs"
    out = cfg.results_dir / "biology" / "08_survival"
    out.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cohorts = sorted(scope.loc[scope["tumor_status_in_scope"], "project"])

    if not args.collect_only:
        jobs = [(data, root / f"{c}__{e}__{a}", c, e, a) for c in cohorts for e in ENDPOINTS for a in ARMS]
        with ProcessPoolExecutor(args.workers) as ex:
            for msg in ex.map(run_one, jobs):
                log(msg)

    if args.null:
        nroot = cfg.data_dir / "survival" / "runs_null"
        ndata = permuted_inputs(cfg, cfg.study["splits"]["seed"])
        njobs = [(ndata, nroot / f"{c}__OS__base", c, "OS", "base") for c in cohorts]
        with ProcessPoolExecutor(args.workers) as ex:
            for msg in ex.map(run_one, njobs):
                log("null " + msg)
        null_cells = pd.concat([pd.read_csv(nroot / f"{c}__OS__base" / "survival.csv", low_memory=False)
                                for c in cohorts if (nroot / f"{c}__OS__base" / "survival.csv").exists()],
                               ignore_index=True)
        ok = null_cells[null_cells["cox_status"].eq("tested")]
        null_rate = pd.DataFrame([{"endpoint": "OS", "arm": "base", "replicate": "permuted outcomes",
                                   "cells": len(ok), "hits_p05": int((ok["cox_p"] < 0.05).sum()),
                                   "rate": float((ok["cox_p"] < 0.05).mean()),
                                   "km_cells": int(null_cells["km_status"].eq("tested").sum()),
                                   "km_hits_p05": int((null_cells.loc[null_cells["km_status"].eq("tested"),
                                                                      "km_p"] < 0.05).sum())}])
        null_rate.to_csv(out / "survival_null.csv", index=False)
        log("permuted-outcome null:\n" + null_rate.round(4).to_string(index=False))

    cells, groups = collect(root, cohorts)
    sel = pd.read_parquet(cfg.data_dir / "survival" / "event_selection.parquet")
    sel["event_id"] = sel["view"].str.upper() + ":" + sel["feature_id"]
    key = sel.drop_duplicates("event_id").set_index("event_id")
    cells["role"] = cells["event_id"].map(key["role"])
    cells["view"] = cells["event_id"].map(key["view"])
    cells["n_models"] = cells["event_id"].map(key["n_models"])
    rec = pd.read_csv(cfg.results_dir / "biology" / "04_events" / "recurrent_events.csv")
    cells["recurrent"] = cells["event_id"].isin(set(rec["view"].str.upper() + ":" + rec["feature_id"]))
    # splice-assay marks a fit it ran as "tested"; the other statuses are its gates
    # (coverage_gate, too_few_patients_or_events, low_psi_variance, failed)
    cells["cox_ok"] = cells["cox_status"].eq("tested") & ~cells["psi_unstable"].fillna(False).astype(bool)
    cells["km_ok"] = cells["km_status"].eq("tested")
    # splice-assay's own cox_q is BH within one event over the cohorts it was run in, and it is left
    # empty for families under 10. The question here is the other one -- does any single event stand
    # out inside a cancer -- so BH again over every event tested in that cohort x endpoint x arm.
    for src, ok, dst in (("cox_p", "cox_ok", "cox_q_cohort"), ("km_p", "km_ok", "km_q_cohort")):
        cells[dst] = np.nan
        m = cells[ok] & cells[src].notna()
        cells.loc[m, dst] = (cells[m].groupby(["cohort", "endpoint", "arm"])[src]
                             .transform(lambda s: bh(s.to_numpy())))
    # published form: one small index of the events, then the cells keyed by it. The rMATS ids
    # are long, and repeating them over 115,000 rows is most of the file.
    index = (cells.drop_duplicates("event_id")[["event_id", "gene", "view", "role", "recurrent", "n_models"]]
             .sort_values("event_id").reset_index(drop=True))
    index.insert(0, "event_key", np.arange(1, len(index) + 1))
    index.to_csv(out / "event_index.csv", index=False)
    kmap = index.set_index("event_id")["event_key"]
    CELL_COLS = ["event_key", "cohort", "endpoint", "arm", "n_survival", "n_obs", "frac_obs", "off_modal",
                 "km_status", "km_n", "events_low", "events_high", "logrank_hr", "km_p", "km_q", "km_q_tests",
                 "cox_status", "cox_n", "cox_events", "cox_beta", "cox_se", "cox_p", "cox_q", "cox_q_tests",
                 "cox_q_cohort", "km_q_cohort",
                 "hr_per_sd", "ci_low_sd", "ci_high_sd", "psi_sd", "hr_per_iqr", "cox_low_power", "psi_narrow",
                 "psi_unstable", "ph_p", "cox_ok", "km_ok"]
    pub = cells.assign(event_key=cells["event_id"].map(kmap))
    pub[[c for c in CELL_COLS if c in pub]].to_csv(out / "event_survival_cells.csv.gz", index=False,
                                                   float_format="%.4g")
    gpub = groups.assign(event_key=groups["event_id"].map(kmap))
    gcols = ["event_key", "cohort"] + [c for c in groups.columns if c.startswith(("paired_", "unpaired_"))
                                       or c in ("group_hit", "within_patient_support", "composition_sensitive")]
    gpub[[c for c in gcols if c in gpub]].to_csv(out / "event_group_tests.csv.gz", index=False,
                                                 float_format="%.4g")
    # the full, untrimmed cells stay local for reruns
    cells.to_parquet(cfg.data_dir / "survival" / "cells_full.parquet", index=False)

    by_cancer = (cells[cells["cox_ok"] & cells["endpoint"].eq("OS") & cells["arm"].eq("adjusted")]
                 .assign(hit=lambda x: x["cox_p"] < 0.05)
                 .groupby(["cohort", "role"], as_index=False)
                 .agg(cells=("hit", "size"), hits=("hit", "sum"), deaths=("cox_events", "median")))
    by_cancer["rate"] = by_cancer["hits"] / by_cancer["cells"]
    by_cancer.to_csv(out / "survival_rate_by_cancer.csv", index=False)

    enr = enrichment(cells, sel)
    enr.to_csv(out / "survival_enrichment.csv", index=False)
    log("enrichment (share of testable cells with p < 0.05):\n" + enr.round(4).to_string(index=False))

    # per event, pooled over cohorts: how often it is associated, and in which direction
    ev = cells[cells["cox_ok"]].copy()
    ev["up"] = ev["cox_beta"] > 0
    per = (ev.groupby(["event_id", "gene", "view", "role", "endpoint", "arm"], as_index=False)
           .agg(cohorts_tested=("cohort", "size"), hits=("cox_p", lambda p: int((p < 0.05).sum())),
                hits_q=("cox_q_cohort", lambda q: int((q.fillna(1) < 0.05).sum())),
                share_hr_up=("up", "mean"), best_p=("cox_p", "min"),
                best_cohort=("cohort", lambda c: c.iloc[0]), median_hr_sd=("hr_per_sd", "median")))
    per["expected_by_chance"] = 0.05 * per["cohorts_tested"]
    best = ev.sort_values("cox_p").drop_duplicates(["event_id", "endpoint", "arm"])[
        ["event_id", "endpoint", "arm", "cohort", "hr_per_sd", "ci_low_sd", "ci_high_sd", "cox_p", "cox_q",
         "cox_n", "cox_events", "cox_low_power", "ph_p"]].rename(columns={"cohort": "best_cohort"})
    per = per.drop(columns=["best_cohort"]).merge(best, on=["event_id", "endpoint", "arm"], how="left")
    per.sort_values(["endpoint", "arm", "hits", "best_p"], ascending=[True, True, False, True]).to_csv(
        out / "event_survival_summary.csv", index=False)
    dos, comp = dossiers(cfg, cells, sel)
    dos.to_csv(out / "event_dossiers.csv", index=False)
    comp.to_csv(out / "regulon_vs_survival.csv", index=False)
    log("regulon membership vs survival association:\n" + comp.round(4).to_string(index=False))
    show = ["event_id", "gene", "n_models_top100", "majority_sign", "hits", "cancers_tested",
            "best_cancer_survival", "hr_per_sd", "cox_p", "specific_factors"]
    log("top dossiers:\n" + dos[[c for c in show if c in dos]].head(12).round(3).to_string(index=False))
    if args.panels:
        ready = dos[dos["hits"].fillna(0) > 0].dropna(subset=["gene"])
        assay_pages(cfg, ready, args.panels, cfg.gencode_gtf)
    log(f"wrote {len(cells)} cells and {per.event_id.nunique()} events to {out}")


if __name__ == "__main__":
    main()
