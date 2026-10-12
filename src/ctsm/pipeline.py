"""One modelling job = one (task, view, cohort): nested CV, tuning, refit, prediction, SHAP.

For each outer fold k (D-04, D-07..D-10):
  1. inner folds j: fit the preprocessor on inner-train rows (both ranking methods, top k_max),
     cache the transformed inner-train / inner-validation matrices
  2. Optuna TPE over the search space; a trial trains one MLP per inner fold with early
     stopping and scores the mean inner-validation balanced log loss (median pruning)
  3. refit: preprocessor on the whole outer-train set with the best trial's selection method
     and K, MLP trained for the median best epoch of that trial's inner folds
  4. predict the outer-test fold; SHAP on (a class-capped subset of) the outer-test fold
The outer-test fold is touched only in step 4. Folds are checkpointed (resume-safe).
"""
from __future__ import annotations

import json
import time
import zlib
from pathlib import Path

import numpy as np
import optuna
import pandas as pd
import torch

from .cohort import class_names
from .config import Config
from .explain import additivity_error, deep_shap
from .matrix_io import load_view
from .metrics import binary_metrics, bootstrap_ci, multiclass_metrics, per_class_table
from .model import logits_to_proba, predict_logits, train_model
from .preprocess import column_stats, fit_preprocessor
from .provenance import git_commit, package_versions, write_json
from .splits import fold_indices


def job_seed(base: int, *keys) -> int:
    return (int(base) + zlib.crc32("|".join(map(str, keys)).encode())) % (2**31 - 1)


def sample_params(trial: optuna.Trial, space: dict, task: str) -> dict:
    return {
        "n_features": trial.suggest_categorical("n_features", space["n_features"]),
        "selection": trial.suggest_categorical("selection", space["selection"]),
        "n_layers": trial.suggest_int("n_layers", *space["n_layers"]),
        "width": trial.suggest_int("width", *space["width"], log=True),
        "dropout": trial.suggest_float("dropout", *space["dropout"]),
        "input_dropout": trial.suggest_float("input_dropout", *space["input_dropout"]),
        "lr": trial.suggest_float("lr", *space["lr"], log=True),
        "weight_decay": trial.suggest_float("weight_decay", *space["weight_decay"], log=True),
        "batch_size": trial.suggest_categorical("batch_size", space["batch_size"][task]),
    }


def load_job_data(cfg: Config, task: str, view: str, cohort: str):
    tab = pd.read_parquet(cfg.data_dir / "tasks" / f"{task}__{cohort}.parquet")
    spl = pd.read_parquet(cfg.data_dir / "splits" / f"{task}__{cohort}.parquet")
    if not (tab["file_id"].to_numpy() == spl["file_id"].to_numpy()).all():
        raise ValueError("task table and split table rows differ")
    X, feats = load_view(cfg.matrix_dir, view, tab["file_id"], cohort)
    return tab, spl, X, feats


def run_fold(cfg: Config, task: str, view: str, cohort: str, k: int, tab, spl, X, feats, fold_dir: Path, log,
             fold_metrics=None) -> dict:
    """`fold_metrics(y_test, proba)` -> {name: value} replaces the task's default per-fold test
    metrics in done.json / folds.csv (used by the subtype task, ctsm/subtype.py)."""
    st, pp = cfg.study, cfg.study["preprocess"]
    space = st["tuning"]["search_space"]
    n_inner = st["splits"]["inner_folds_used"][task]
    y = tab["y"].to_numpy()
    C = int(y.max()) + 1
    # one-logit model iff two classes (model.train_model); tumor_status always has 2 classes and
    # cancer_type 34, so this equals the former `task == "tumor_status"` for both existing tasks
    binary = C == 2
    k_max = max(space["n_features"])
    seed = job_seed(st["splits"]["seed"], task, view, cohort, k)
    t0 = time.time()

    # 1. inner preprocessing cache
    cache = {}
    for j in range(n_inner):
        itr, iva = fold_indices(spl, k, j)
        stats = column_stats(X, itr, y[itr], C, view, pp["gex_min_log2cpm"])
        for method in space["selection"]:
            pre = fit_preprocessor(X, itr, y[itr], view, method, k_max, pp, C, stats=stats)
            cache[(j, method)] = (pre.transform(X, itr), pre.transform(X, iva), y[itr], y[iva], pre.n_candidates)
    log(f"  fold {k}: inner preprocessing {time.time() - t0:.0f}s, "
        f"{cache[(0, space['selection'][0])][4]:,} candidate features in inner fold 0")

    # 2. tuning
    max_epochs = st["model"]["max_epochs"][task]
    patience = st["model"]["patience"][task]

    def objective(trial: optuna.Trial) -> float:
        p = sample_params(trial, space, task)
        losses, epochs = [], []
        for j in range(n_inner):
            Za, Zb, ya, yb, _ = cache[(j, p["selection"])]
            K = min(p["n_features"], Za.shape[1])
            _, ep, v, _ = train_model(Za[:, :K], ya, p, C, seed + 1000 * trial.number + j, X_val=Zb[:, :K], y_val=yb,
                                      max_epochs=max_epochs, patience=patience)
            losses.append(v)
            epochs.append(ep)
            trial.report(float(np.mean(losses)), j)
            if trial.should_prune():
                raise optuna.TrialPruned()
        trial.set_user_attr("best_epochs", epochs)
        trial.set_user_attr("fold_losses", losses)
        return float(np.mean(losses))

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="minimize", sampler=optuna.samplers.TPESampler(seed=seed),
                                pruner=optuna.pruners.MedianPruner(n_startup_trials=st["tuning"]["pruner_startup_trials"],
                                                                   n_warmup_steps=0))
    study.optimize(objective, n_trials=st["tuning"]["n_trials"][task])
    best = study.best_trial
    params = dict(best.params)
    n_epochs = max(1, int(round(float(np.median(best.user_attrs["best_epochs"])))))
    t_tune = time.time() - t0
    del cache

    # 3. refit on the whole outer-train set
    tr, te = fold_indices(spl, k)
    pre = fit_preprocessor(X, tr, y[tr], view, params["selection"], params["n_features"], pp, C)
    Ztr, Zte = pre.transform(X, tr), pre.transform(X, te)
    model, _, _, _ = train_model(Ztr, y[tr], params, C, seed, fixed_epochs=n_epochs)

    # 4. outer-test prediction and SHAP
    logits = predict_logits(model, Zte)
    proba = logits_to_proba(logits, binary)
    sh = st["shap"]
    rng = np.random.default_rng(seed)
    bg = Ztr[np.sort(rng.choice(len(tr), size=min(sh["background_size"], len(tr)), replace=False))]
    sv, ev = deep_shap(model, bg, Zte, y[te], C)
    add_err = additivity_error(model, Zte, y[te], C, sv, ev)
    if add_err > 1e-3:
        raise AssertionError(f"SHAP additivity error {add_err:.2e} on fold {k}")

    # checkpoint
    fold_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), fold_dir / "model.pt")
    np.savez_compressed(fold_dir / "preprocessor.npz", **{k_: v for k_, v in pre.to_dict().items()
                                                         if isinstance(v, np.ndarray)},
                        feature_ids=feats["feature_id"].to_numpy()[pre.cols].astype(str))
    np.savez_compressed(fold_dir / "shap.npz", shap=sv.astype(np.float32), X=Zte.astype(np.float16),
                        expected_value=ev, rows=te, feature_ids=feats["feature_id"].to_numpy()[pre.cols].astype(str))
    np.savez_compressed(fold_dir / "predictions.npz", rows=te, proba=proba, logits=logits)
    study.trials_dataframe().to_csv(fold_dir / "trials.csv", index=False)
    states = pd.Series([t.state.name for t in study.trials]).value_counts().to_dict()
    rec = {"outer_fold": k, "seed": seed, "best_inner_loss": float(best.value), "best_trial": best.number,
           "inner_fold_losses": best.user_attrs["fold_losses"], "inner_best_epochs": best.user_attrs["best_epochs"],
           "refit_epochs": n_epochs, **{f"param_{a}": b for a, b in params.items()},
           "n_features_used": int(len(pre.cols)), "n_candidates_outer_train": pre.n_candidates,
           "trials_complete": int(states.get("COMPLETE", 0)), "trials_pruned": int(states.get("PRUNED", 0)),
           "n_train": int(len(tr)), "n_test": int(len(te)), "n_shap_explained": int(len(te)),
           "shap_additivity_max_error": add_err, "seconds_tuning": round(t_tune, 1),
           "seconds_total": round(time.time() - t0, 1)}
    if fold_metrics is not None:
        rec.update({f"test_{a}": b for a, b in fold_metrics(y[te], proba).items()})
    elif binary:
        rec.update({f"test_{a}": b for a, b in binary_metrics(y[te], proba[:, 1]).items()
                    if a in ("auroc", "balanced_accuracy", "auprc_normal", "mcc")})
    else:
        rec.update({f"test_{a}": b for a, b in multiclass_metrics(y[te], proba).items()
                    if a in ("accuracy", "balanced_accuracy", "macro_f1")})
    write_json(fold_dir / "done.json", rec)
    log(f"  fold {k}: done in {rec['seconds_total']:.0f}s; best inner loss {best.value:.4f} "
        f"({rec['trials_complete']} complete / {rec['trials_pruned']} pruned); "
        f"K={len(pre.cols)} {params['selection']}, {params['n_layers']}x{params['width']}; "
        + ", ".join(f"{a[5:]}={b:.3f}" for a, b in rec.items() if a.startswith("test_")))
    return rec


def run_job(cfg: Config, task: str, view: str, cohort: str, log=print) -> dict:
    t0 = time.time()
    tab, spl, X, feats = load_job_data(cfg, task, view, cohort)
    classes = class_names(tab)
    log(f"{task}/{view}/{cohort}: X {X.shape} ({X.encoding}), {len(classes)} classes, "
        f"missing {X.missing_frac():.3f}")
    rdir = cfg.run_dir(task, view, cohort)
    recs = []
    for k in range(cfg.study["splits"]["outer_folds"]):
        fd = rdir / f"fold{k}"
        if (fd / "done.json").exists():
            recs.append(json.loads((fd / "done.json").read_text()))
            log(f"  fold {k}: already done")
            continue
        recs.append(run_fold(cfg, task, view, cohort, k, tab, spl, X, feats, fd, log))
    summary = summarise_job(cfg, task, view, cohort, tab, classes, recs)
    summarise_shap(cfg, task, view, cohort, tab, classes, feats)
    summary["seconds_wall_this_run"] = round(time.time() - t0, 1)
    write_json(rdir / "job_done.json", summary)
    return summary


def collect_predictions(cfg: Config, task: str, view: str, cohort: str, tab: pd.DataFrame, n_classes: int) -> np.ndarray:
    proba = np.full((len(tab), n_classes), np.nan)
    for k in range(cfg.study["splits"]["outer_folds"]):
        z = np.load(cfg.run_dir(task, view, cohort) / f"fold{k}" / "predictions.npz")
        proba[z["rows"]] = z["proba"]
    if np.isnan(proba).any():
        raise ValueError("a sample has no out-of-fold prediction")
    return proba


def summarise_job(cfg: Config, task: str, view: str, cohort: str, tab: pd.DataFrame, classes: list[str],
                  recs: list[dict]) -> dict:
    """Pooled out-of-fold metrics (+ case-bootstrap CIs) -> results/<task>/<view>/<cohort>/."""
    binary = task == "tumor_status"
    y = tab["y"].to_numpy()
    proba = collect_predictions(cfg, task, view, cohort, tab, len(classes))
    rdir, out = cfg.run_dir(task, view, cohort), cfg.result_dir(task, view, cohort)
    out.mkdir(parents=True, exist_ok=True)
    pred = tab[["file_id", "case_id", "sample_id", "project", "label", "y"]].copy()
    pred["outer_fold"] = pd.read_parquet(cfg.data_dir / "splits" / f"{task}__{cohort}.parquet")["outer"].to_numpy()
    for c, name in enumerate(classes):
        pred[f"p_{name}"] = proba[:, c]
    pred.to_parquet(rdir / "oof_predictions.parquet", index=False)     # per sample: kept local

    if binary:
        pooled = binary_metrics(y, proba[:, 1])
        ci = bootstrap_ci(y, proba, tab["case_id"].to_numpy(), True, ["auroc", "balanced_accuracy", "auprc_normal", "mcc"])
    else:
        pooled = multiclass_metrics(y, proba)
        ci = bootstrap_ci(y, proba, tab["case_id"].to_numpy(), False, ["accuracy", "balanced_accuracy", "macro_f1"],
                          n_boot=200)
        pd.DataFrame(per_class_table(y, proba, classes)).to_csv(out / "per_class_metrics.csv", index=False)
        cm = pd.crosstab(pd.Categorical(np.array(classes)[y], categories=classes),
                         pd.Categorical(np.array(classes)[proba.argmax(axis=1)], categories=classes), dropna=False)
        cm.index.name, cm.columns.name = "true", "predicted"
        cm.to_csv(out / "confusion_matrix.csv")
    folds = pd.DataFrame(recs)
    folds.to_csv(out / "folds.csv", index=False)
    test_cols = [c for c in folds.columns if c.startswith("test_")]
    summary = {
        "task": task, "view": view, "cohort": cohort, "classes": classes,
        "pooled_oof": pooled, "bootstrap_95ci": ci,
        "per_fold_mean": {c[5:]: float(folds[c].mean()) for c in test_cols},
        "per_fold_sd": {c[5:]: float(folds[c].std(ddof=1)) for c in test_cols},
        "selection_method_by_fold": folds["param_selection"].tolist(),
        "n_features_by_fold": folds["n_features_used"].tolist(),
        "git_commit": git_commit(), "packages": package_versions(),
        "config": {"splits": cfg.study["splits"], "preprocess": cfg.study["preprocess"],
                   "model": cfg.study["model"], "tuning": cfg.study["tuning"], "shap": cfg.study["shap"]},
    }
    write_json(out / "metrics.json", summary)
    return summary


def summarise_shap(cfg: Config, task: str, view: str, cohort: str, tab: pd.DataFrame, classes: list[str],
                   feats: pd.DataFrame) -> None:
    """Fold SHAP files -> results/<task>/<view>/<cohort>/shap_*.csv (D-10)."""
    from .explain import gene_table, importance_table, load_folds
    folds = load_folds(cfg.run_dir(task, view, cohort), tab["y"].to_numpy(), cfg.study["splits"]["outer_folds"])
    ann = feats[["feature_id", "label", "symbol", "gene_id", "region"]]
    out = cfg.result_dir(task, view, cohort)
    top = cfg.study["shap"]["top_k_report"]
    # pan-cancer: each sample explains its own class, so a pooled direction would mix classes;
    # directions are reported per class instead
    t = importance_table(folds, ann, n_direction=200 if task == "tumor_status" else 0)
    if task == "cancer_type":
        t = t.drop(columns="direction_spearman")
    t.head(200).to_csv(out / "shap_top_features.csv", index=False)
    gene_table(t, "mean_abs_shap").head(200).to_csv(out / "shap_top_genes.csv", index=False)
    if task == "cancer_type":
        per = []
        for k, name in enumerate(classes):
            tk = importance_table(folds, ann, mask_of=lambda f, k=k: f["y"] == k, n_direction=top).head(top)
            per.append(tk.rename(columns={"rank": "rank_in_class"}).assign(**{"class": name}))
        pd.concat(per, ignore_index=True).to_csv(out / "shap_top_features_by_class.csv", index=False)
