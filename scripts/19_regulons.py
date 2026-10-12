"""Step 19 (B-10): splicing-factor regulons — which events move with which factor, in which
cancers, and do the models' events sit inside those programmes?

Step 12 asked whether the top events correlate with a factor. This asks the regulon question:
for each factor, the set of events that track it reproducibly across cancers, and whether that
set behaves as a programme.

Method, within the TUMOURS of each cancer separately (normals would create a two-group
correlation out of any tumour-normal difference):
  1. Event universe: a fixed random sample of events per view (seeded, the same list in every
     cancer) plus every step-16 candidate and matched background event, each required to be
     observed in >= 80% of that cancer's tumours and non-constant.
  2. Partial Spearman of factor expression (log2 CPM + 1) against event PSI, both rank-
     transformed and residualised on three sample-level covariates: a stromal and an immune
     marker score (D-18), and the sample's MEAN PSI over the view's measured events. That third
     term matters. Without it nearly every nuclear RNA-processing factor appears to "regulate"
     a thousand events, almost all in the same direction, because factor expression and overall
     PSI level share one sample-level axis (intron retention / nuclear RNA content, partly
     technical). Residualising on it leaves the event-specific deviation from each sample's own
     level, which is what a regulon should mean.
  3. A pair is a regulon member in that cancer if |rho| >= 0.3 and p < 0.05. A factor's CORE
     regulon is the events that pass in >= min_cancers cancers with the same sign.
  4. Regulon activity per sample: the mean z-scored PSI of the core events, signed by their
     direction. Its paired tumour - normal change says whether the programme moves in tumours,
     and its correlation with the factor's own expression says whether it tracks the factor.
  5. Enrichment, pooled and per factor: are the models' top events over-represented in core
     regulons, against the coverage- and spread-matched background of step 16? Same matched
     design as the survival layer, so the comparison is not about measurability.
  6. Overlap: how much each factor's core regulon is shared with another factor's. Large
     regulons turn out to be one broad programme measured many times, so a factor-by-factor
     reading would badly overstate how specific they are.

Outputs: results/biology/10_regulons/ (tables; figures: scripts/20_figures_status.py)
"""
from __future__ import annotations

import argparse
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import rankdata, t as tdist, wilcoxon

from ctsm.biology import bh, bio_config, feature_values, matched_pairs, task_table
from ctsm.config import load_config
from ctsm.matrix_io import decode, load_meta, view_dir

warnings.filterwarnings("ignore", category=RuntimeWarning)
AS_VIEWS = ("afe", "ale", "hit", "se", "a5ss", "a3ss", "mxe", "ri")
MIN_OBS = 0.8
RHO_MIN = 0.3
N_RANDOM_PER_VIEW = 2000


def ranked(X: np.ndarray) -> np.ndarray:
    return np.apply_along_axis(rankdata, 0, X).astype(np.float64)


def residualise(X: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Least-squares residuals of every column of X on [1, C]."""
    A = np.hstack([np.ones((len(C), 1)), C])
    beta, *_ = np.linalg.lstsq(A, X, rcond=None)
    return X - A @ beta


def partial_spearman(X: np.ndarray, Y: np.ndarray, C: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """rho and p for every column pair (events x factors), ranks residualised on C."""
    Xr, Yr, Cr = residualise(ranked(X), ranked(C)), residualise(ranked(Y), ranked(C)), None
    Xz = (Xr - Xr.mean(0)) / np.where(Xr.std(0) > 0, Xr.std(0), np.nan)
    Yz = (Yr - Yr.mean(0)) / np.where(Yr.std(0) > 0, Yr.std(0), np.nan)
    n, k = len(X), C.shape[1]
    rho = (Xz.T @ Yz) / n
    rho = np.clip(rho, -0.999999, 0.999999)
    stat = rho * np.sqrt((n - 2 - k) / (1 - rho ** 2))
    return rho, 2 * tdist.sf(np.abs(stat), n - 2 - k)


def event_universe(cfg, rng) -> pd.DataFrame:
    """A fixed event list: a random sample per view plus every step-16 candidate and background."""
    sel = pd.read_parquet(cfg.data_dir / "survival" / "event_selection.parquet")
    rows = []
    for view in AS_VIEWS:
        feats = pd.read_parquet(cfg.matrix_dir / view / "features.parquet",
                                columns=["feature_id", "gene_id", "symbol", "label"])
        chosen = set(sel.loc[sel["view"].eq(view), "feature_id"])
        pool = feats[~feats["feature_id"].isin(chosen)]
        take = pool.iloc[rng.choice(len(pool), size=min(N_RANDOM_PER_VIEW, len(pool)), replace=False)]
        part = pd.concat([feats[feats["feature_id"].isin(chosen)], take], ignore_index=True)
        rows.append(part.assign(view=view))
    out = pd.concat(rows, ignore_index=True).drop_duplicates(["view", "feature_id"])
    roles = sel.drop_duplicates(["view", "feature_id"]).set_index(["view", "feature_id"])["role"]
    out["role"] = pd.MultiIndex.from_frame(out[["view", "feature_id"]]).map(roles).fillna("universe")
    return out


def gene_scores(cfg, file_ids, genes) -> np.ndarray:
    """Mean log2 CPM over a marker set (missing genes ignored)."""
    feats = pd.read_parquet(cfg.matrix_dir / "gex" / "features.parquet", columns=["feature_id", "symbol"])
    have = feats[feats["symbol"].isin(genes)]["feature_id"]
    if not len(have):
        return np.zeros((len(file_ids), 1))
    return feature_values(cfg, "gex", file_ids, have).mean(1, keepdims=True)


def factor_matrix(cfg, file_ids, factors) -> tuple[np.ndarray, list[str]]:
    feats = pd.read_parquet(cfg.matrix_dir / "gex" / "features.parquet", columns=["feature_id", "symbol"])
    have = feats[feats["symbol"].isin(factors)].drop_duplicates("symbol")
    return feature_values(cfg, "gex", file_ids, have["feature_id"]), list(have["symbol"])


def one_cancer(args) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cancer, universe_ids = args
    cfg, bio = load_config(), bio_config()
    reg = bio["regulators"]
    factors = [g for group in reg["splicing_factors"].values() for g in group]
    tab = task_table(cfg, "tumor_status", cancer)
    tum = tab[tab["y"].eq(1)]["file_id"]
    Y, names = factor_matrix(cfg, tum, factors)
    C = np.hstack([gene_scores(cfg, tum, reg["stromal_genes"]), gene_scores(cfg, tum, reg["immune_genes"])])
    out, tested, axis = [], [], []
    for view, part in universe_ids.groupby("view", sort=False):
        X = feature_values(cfg, view, tum, part["feature_id"])
        keep = (np.isfinite(X).mean(0) >= MIN_OBS) & (np.nanstd(X, axis=0) > 1e-6)
        if keep.sum() == 0:
            continue
        X = X[:, keep]
        med = np.nanmedian(X, axis=0)
        X = np.where(np.isfinite(X), X, med)
        glob = X.mean(1, keepdims=True)                 # the view's global per-sample level
        Cv = np.hstack([C, glob])
        axis_rho, _ = partial_spearman(glob, Y, C)      # how much each factor tracks that axis
        rho, p = partial_spearman(X, Y, Cv)
        ii, jj = np.nonzero((np.abs(rho) >= RHO_MIN) & (p < 0.05))
        sub = part[keep].reset_index(drop=True)
        tested.append(sub.assign(cancer=cancer))
        axis.append(pd.DataFrame({"cancer": cancer, "view": view, "factor": names,
                                  "rho_global_axis": axis_rho[0]}))
        out.append(pd.DataFrame({"cancer": cancer, "view": view, "factor": [names[j] for j in jj],
                                 "feature_id": sub["feature_id"].to_numpy()[ii],
                                 "symbol": sub["symbol"].to_numpy()[ii], "role": sub["role"].to_numpy()[ii],
                                 "rho": rho[ii, jj], "p": p[ii, jj], "n_tumours": len(tum),
                                 "n_tested": int(keep.sum())}))
    return (pd.concat(out, ignore_index=True) if out else pd.DataFrame(),
            pd.concat(tested, ignore_index=True) if tested else pd.DataFrame(),
            pd.concat(axis, ignore_index=True) if axis else pd.DataFrame())


def factor_enrichment(core: pd.DataFrame, tested: pd.DataFrame, min_marked: int = 15) -> pd.DataFrame:
    """Per factor: are the models' top events over-represented in its core regulon?"""
    from scipy.stats import fisher_exact
    den = tested.groupby("role").size()
    c_n, b_n = int(den.get("candidate", 0)), int(den.get("background", 0))
    rows = []
    for factor, sub in core.groupby("factor"):
        marked = tested[tested["feature_id"].isin(set(sub["feature_id"]))].groupby("role").size()
        c_in, b_in = int(marked.get("candidate", 0)), int(marked.get("background", 0))
        if c_in + b_in < min_marked:
            continue
        odds, p = fisher_exact([[c_in, c_n - c_in], [b_in, b_n - b_in]])
        haldane = (((c_in + 0.5) * (b_n - b_in + 0.5)) /       # finite even when a cell is empty (for plotting)
                   ((c_n - c_in + 0.5) * (b_in + 0.5)))
        rows.append({"factor": factor, "n_core_events": int(sub["feature_id"].nunique()),
                     "candidates_in": c_in, "background_in": b_in, "n_candidates": c_n,
                     "n_background": b_n, "odds_ratio": odds, "odds_ratio_haldane": haldane,
                     "p_fisher": p})
    out = pd.DataFrame(rows)
    out["q_fisher"] = bh(out["p_fisher"].to_numpy())
    return out.sort_values("p_fisher")


def overlap_structure(core: pd.DataFrame, min_events: int = 20) -> pd.DataFrame:
    """How much of each factor's core regulon is also another factor's (the redundancy of the set)."""
    counts = core.groupby("factor")["feature_id"].nunique()
    keep = counts[counts >= min_events].index
    sets = {f: set(core.loc[core["factor"].eq(f), "feature_id"]) for f in keep}
    per_event = core.groupby("feature_id")["factor"].nunique()
    rows = []
    for a in keep:
        best, who = 0.0, ""
        for b in keep:
            if a == b:
                continue
            u = len(sets[a] | sets[b])
            j = len(sets[a] & sets[b]) / u if u else 0.0
            if j > best:
                best, who = j, b
        rows.append({"factor": a, "n_core_events": len(sets[a]), "max_jaccard": best,
                     "closest_factor": who,
                     "median_factors_per_event": float(per_event[list(sets[a])].median()),
                     "share_up": float((core.loc[core["factor"].eq(a), "sign"] > 0).mean())})
    return pd.DataFrame(rows).sort_values("max_jaccard")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--min-cancers", type=int, default=5)
    args = ap.parse_args()
    cfg = load_config()
    rng = np.random.default_rng(cfg.study["splits"]["seed"])
    out = cfg.results_dir / "biology" / "10_regulons"
    out.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = sorted(scope.loc[scope["tumor_status_in_scope"], "project"])

    universe = event_universe(cfg, rng)
    log(f"event universe: {len(universe)} events "
        f"({(universe.role == 'candidate').sum()} candidates, {(universe.role == 'background').sum()} background)")

    cache = cfg.data_dir / "bio_cache" / "regulon_pairs.parquet"
    tcache = cfg.data_dir / "bio_cache" / "regulon_tested.parquet"
    acache = cfg.data_dir / "bio_cache" / "regulon_axis.parquet"
    if cache.exists() and tcache.exists() and acache.exists():
        pairs, tested, axis = pd.read_parquet(cache), pd.read_parquet(tcache), pd.read_parquet(acache)
        log(f"cached: {len(pairs)} factor-event pairs")
    else:
        with ProcessPoolExecutor(args.workers) as ex:
            parts = list(ex.map(one_cancer, [(c, universe) for c in cancers]))
        pairs = pd.concat([p for p, _, _ in parts], ignore_index=True)
        tested = pd.concat([t for _, t, _ in parts], ignore_index=True)
        axis = pd.concat([a for _, _, a in parts], ignore_index=True)
        cache.parent.mkdir(parents=True, exist_ok=True)
        pairs.to_parquet(cache, index=False)
        tested.to_parquet(tcache, index=False)
        axis.to_parquet(acache, index=False)
        log(f"{len(pairs)} factor-event pairs at |rho| >= {RHO_MIN}, p < 0.05")

    # ---- core regulons: the same factor and event, same sign, in >= min_cancers cancers
    pairs["sign"] = np.sign(pairs["rho"])
    g = (pairs.groupby(["factor", "view", "feature_id", "symbol", "role", "sign"], as_index=False)
         .agg(n_cancers=("cancer", "nunique"), median_rho=("rho", "median"),
              cancers=("cancer", lambda c: ",".join(sorted(c)))))
    core = g[g["n_cancers"] >= args.min_cancers].sort_values(["factor", "n_cancers"], ascending=[True, False])
    core.to_csv(out / "regulon_core_events.csv", index=False)
    sizes = (core.groupby("factor", as_index=False)
             .agg(core_events=("feature_id", "size"), up=("sign", lambda s: int((s > 0).sum())),
                  down=("sign", lambda s: int((s < 0).sum())),
                  views=("view", lambda v: ",".join(sorted(set(v))))))
    per_cancer = (pairs.groupby(["factor", "cancer"], as_index=False)
                  .agg(n_events=("feature_id", "size"), median_abs_rho=("rho", lambda r: float(np.abs(r).median()))))
    sizes = sizes.merge(per_cancer.groupby("factor", as_index=False).agg(
        mean_events_per_cancer=("n_events", "mean"), cancers_with_regulon=("cancer", "nunique")), on="factor")
    sizes.sort_values("core_events", ascending=False).to_csv(out / "regulon_sizes.csv", index=False)
    per_cancer.to_csv(out / "regulon_by_cancer.csv", index=False)
    ax = (axis.groupby("factor", as_index=False)
          .agg(median_rho_global_axis=("rho_global_axis", "median"),
               n_cells=("rho_global_axis", "size"),
               share_abs_gt_03=("rho_global_axis", lambda r: float((np.abs(r) > 0.3).mean()))))
    ax.sort_values("median_rho_global_axis", ascending=False).to_csv(out / "global_axis_factors.csv", index=False)
    log("factors most tied to the global per-sample PSI level (before residualising it out):\n"
        + ax.nlargest(8, "median_rho_global_axis").to_string(index=False))
    log("largest core regulons:\n" + sizes.nlargest(12, "core_events").to_string(index=False))

    # ---- enrichment of model-selected events in core regulons (matched background)
    denom = tested.drop_duplicates(["view", "feature_id"]).groupby("role").size()
    hit = core.drop_duplicates(["view", "feature_id"]).groupby("role").size()
    rows = []
    for role in ("candidate", "background"):
        rows.append({"role": role, "n_tested": int(denom.get(role, 0)), "in_core_regulon": int(hit.get(role, 0))})
    enr = pd.DataFrame(rows)
    enr["share"] = enr["in_core_regulon"] / enr["n_tested"]
    from scipy.stats import fisher_exact
    a, b = enr.set_index("role").loc["candidate"], enr.set_index("role").loc["background"]
    odds, p = fisher_exact([[a.in_core_regulon, a.n_tested - a.in_core_regulon],
                            [b.in_core_regulon, b.n_tested - b.in_core_regulon]])
    enr["odds_ratio"], enr["p_fisher"] = odds, p
    enr.to_csv(out / "regulon_enrichment.csv", index=False)
    log("\n" + enr.round(4).to_string(index=False))

    fe = factor_enrichment(core, tested.drop_duplicates(["view", "feature_id"]))
    fe.to_csv(out / "regulon_enrichment_by_factor.csv", index=False)
    log("per-factor enrichment of the models' events (top and bottom):\n"
        + pd.concat([fe.head(8), fe.tail(3)]).round(4).to_string(index=False))
    ov = overlap_structure(core)
    ov.to_csv(out / "regulon_overlap.csv", index=False)
    log(f"median largest pairwise Jaccard among {len(ov)} factors: {ov.max_jaccard.median():.2f}; "
        f"{core.feature_id.nunique()} distinct events carry {len(core)} memberships")


if __name__ == "__main__":
    main()
