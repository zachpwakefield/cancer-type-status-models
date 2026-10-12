"""Step 34: subtype performance across cohorts and views, baselines, AS ensemble and fusion (D-39, D-40).

Reads results/subtype/<view>/<cohort>/metrics.json and the local out-of-fold predictions.
Every combined model is built from out-of-fold probabilities only (each patient was predicted by
models that never saw it), with weights fixed in advance:
    AS ensemble   mean of the 8 splicing views' class probabilities
    GEX + AS      0.5 x GEX + 0.5 x AS ensemble
Tests against GEX: exact McNemar on paired per-patient correctness, and a paired patient-bootstrap
95% CI of the difference in balanced accuracy. BH across cohorts, per comparison.

Outputs (results/subtype/summary/, aggregates only):
    performance.csv         cohort x model: balanced accuracy, macro F1, macro AUROC, accuracy, CIs
    baselines.csv           cohort: chance and majority-class levels
    recall_by_class.csv     cohort x model x class: recall (CI for single views)
    fusion_tests.csv        cohort x comparison: discordant pairs, McNemar p and q, delta BA and CI
    confusion_long.csv      cohort x model (GEX, AS ensemble, GEX + AS) x true x predicted: patients
Usage: python scripts/34_subtype_summarise.py
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import binomtest

from ctsm.biology import bh
from ctsm.figstyle import VIEW_NAMES
from ctsm.subtype import TASK, bootstrap, cohort_specs, subtype_config, subtype_metrics

N_BOOT = 1000
COMBINED = ["AS ensemble", "GEX + AS"]


def oof(scfg, view: str, cohort: str) -> pd.DataFrame:
    return pd.read_parquet(scfg.run_dir(TASK, view, cohort) / "oof_predictions.parquet").sort_values(
        "file_id").reset_index(drop=True)


def mcnemar(ref_ok: np.ndarray, new_ok: np.ndarray) -> dict:
    a = int((ref_ok & ~new_ok).sum())
    b = int((~ref_ok & new_ok).sum())
    return {"only_gex_correct": a, "only_model_correct": b,
            "p_mcnemar": float(binomtest(a, a + b, 0.5).pvalue) if a + b else 1.0}


def paired_delta_ba(y, pa, pb, groups, n_boot=N_BOOT, seed=0) -> tuple[float, float, float]:
    """Balanced accuracy of b minus a, with a paired patient-bootstrap 95% CI."""
    C = pa.shape[1]

    def ba(p, idx):
        yy, pr = y[idx], p[idx].argmax(1)
        return np.mean([np.mean(pr[yy == c] == c) for c in range(C)])
    full = np.arange(len(y))
    d0 = ba(pb, full) - ba(pa, full)
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(groups, return_inverse=True)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(uniq), len(uniq))
        idx = np.concatenate([np.flatnonzero(inv == g) for g in idx]) if len(uniq) != len(y) else idx
        if len(np.unique(y[idx])) < C:
            continue
        vals.append(ba(pb, idx) - ba(pa, idx))
    return float(d0), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> None:
    scfg = subtype_config()
    specs = cohort_specs(scfg.spec)
    out = scfg.results_dir / TASK / "summary"
    out.mkdir(parents=True, exist_ok=True)
    as_views = [v for v in scfg.views if v != "gex"]
    perf, base, recall, tests, conf = [], [], [], [], []
    for cohort, cs in specs.items():
        mpaths = {v: scfg.result_dir(TASK, v, cohort) / "metrics.json" for v in scfg.views}
        have = [v for v, p in mpaths.items() if p.exists()]
        if not have:
            continue
        meta = {"cohort": cohort, "label_set": cs["label_set"], "label_type": cs["label_type"], "tier": cs["tier"],
                "project": cs["project"]}
        ms = {v: json.loads(mpaths[v].read_text()) for v in have}
        classes = ms[have[0]]["classes"]
        C = len(classes)
        b = ms[have[0]]["baselines"]
        n = int(sum(ms[have[0]]["class_counts"].values()))
        base.append({**meta, "n": n, "n_classes": C, **b})
        for v in have:
            m, ci = ms[v]["pooled_oof"], ms[v]["bootstrap_95ci"]
            perf.append({**meta, "model": VIEW_NAMES[v], "n": n, "n_classes": C,
                         **{k: m[k] for k in ("balanced_accuracy", "macro_f1", "macro_auroc_ovr", "accuracy", "mcc")},
                         **{f"{k}_lo": ci[k][0] for k in ("balanced_accuracy", "macro_f1", "macro_auroc_ovr", "accuracy")},
                         **{f"{k}_hi": ci[k][1] for k in ("balanced_accuracy", "macro_f1", "macro_auroc_ovr", "accuracy")}})
            pc = pd.read_csv(scfg.result_dir(TASK, v, cohort) / "per_class_metrics.csv")
            for r in pc.itertuples(index=False):
                recall.append({**meta, "model": VIEW_NAMES[v], "class": r[0], "support": r.support, "recall": r.recall,
                               "recall_lo": r.recall_lo, "recall_hi": r.recall_hi})
        if set(have) != set(scfg.views):
            log(f"{cohort}: {len(have)}/9 views done; combined models skipped")
            continue
        P = {}
        ref = None
        for v in scfg.views:
            d = oof(scfg, v, cohort)
            if ref is None:
                ref = d
            if not (d["file_id"].to_numpy() == ref["file_id"].to_numpy()).all():
                raise ValueError(f"{cohort}/{v}: prediction rows differ")
            P[v] = d[[f"p_{c}" for c in classes]].to_numpy()
        y, groups = ref["y"].to_numpy(), ref["case_id"].to_numpy()
        models = {"GEX": P["gex"], "AS ensemble": np.mean([P[v] for v in as_views], axis=0)}
        models["GEX + AS"] = 0.5 * models["GEX"] + 0.5 * models["AS ensemble"]
        ok = {k: p.argmax(1) == y for k, p in models.items()}
        for name in COMBINED:
            p = models[name]
            m = subtype_metrics(y, p)
            ci = bootstrap(y, p, groups, n_boot=N_BOOT)
            perf.append({**meta, "model": name, "n": n, "n_classes": C,
                         **{k: m[k] for k in ("balanced_accuracy", "macro_f1", "macro_auroc_ovr", "accuracy", "mcc")},
                         **{f"{k}_lo": ci[k][0] for k in ("balanced_accuracy", "macro_f1", "macro_auroc_ovr", "accuracy")},
                         **{f"{k}_hi": ci[k][1] for k in ("balanced_accuracy", "macro_f1", "macro_auroc_ovr", "accuracy")}})
            for c, cls in enumerate(classes):
                recall.append({**meta, "model": name, "class": cls, "support": int((y == c).sum()),
                               "recall": m[f"recall_{c}"], "recall_lo": ci[f"recall_{c}"][0],
                               "recall_hi": ci[f"recall_{c}"][1]})
            d, lo, hi = paired_delta_ba(y, models["GEX"], p, groups)
            tests.append({**meta, "comparison": f"{name} vs GEX", **mcnemar(ok["GEX"], ok[name]),
                          "delta_balanced_accuracy": d, "delta_ba_lo": lo, "delta_ba_hi": hi})
        for name in ("GEX", *COMBINED):
            ct = pd.crosstab(pd.Categorical(np.array(classes)[y], categories=classes),
                             pd.Categorical(np.array(classes)[models[name].argmax(1)], categories=classes), dropna=False)
            for t_ in classes:
                for p_ in classes:
                    conf.append({**meta, "model": name, "true": t_, "predicted": p_, "patients": int(ct.loc[t_, p_])})
        r = pd.DataFrame(perf)
        r = r[r.cohort == cohort].set_index("model")["balanced_accuracy"]
        log(f"{cohort:<11} BA  GEX {r['GEX']:.3f}  best AS {r[[VIEW_NAMES[v] for v in as_views]].max():.3f}  "
            f"AS ens {r['AS ensemble']:.3f}  fused {r['GEX + AS']:.3f}  (chance {1 / C:.3f})")
    t = pd.DataFrame(tests)
    if len(t):
        t["q_mcnemar"] = np.nan
        for comp, g in t.groupby("comparison"):
            t.loc[g.index, "q_mcnemar"] = bh(g["p_mcnemar"].to_numpy())
        t["q_family"] = "BH across the subtype cohorts, per comparison (" + t["comparison"] + ")"
    pd.DataFrame(perf).to_csv(out / "performance.csv", index=False)
    pd.DataFrame(base).to_csv(out / "baselines.csv", index=False)
    pd.DataFrame(recall).to_csv(out / "recall_by_class.csv", index=False)
    t.to_csv(out / "fusion_tests.csv", index=False)
    pd.DataFrame(conf).to_csv(out / "confusion_long.csv", index=False)
    log(f"wrote {out}")


if __name__ == "__main__":
    main()
