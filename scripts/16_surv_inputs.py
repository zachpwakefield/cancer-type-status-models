"""Step 16 (B-08): build the local splice-assay dataset that asks whether the events the
tumour-status models rely on are also associated with outcome.

The selection of events is OUTCOME-BLIND, which is what makes the later test meaningful:
  candidates  the top-k features (mean |SHAP| over held-out samples) of every tumour-status
              model, pooled over the 8 AS views and the 15 cancers;
  background  for every candidate, one event of the same view drawn from the same coverage and
              spread stratum (quintiles of observed fraction and of SD, both computed from the
              tumours alone), never itself a top-100 feature of any model.
Nothing about survival, stage, age or sex enters either set. The background is what makes the
enrichment test fair: candidates are well measured and variable, and those two properties alone
make an event easier to test.

Written (local only, never committed; it holds per-patient values):
    <data_dir>/survival/inputs/{samples,psi,expression}.parquet   splice-assay's default form
    <data_dir>/survival/event_selection.parquet                   candidate / background, with the
                                                                  matching stratum of each pair
Step 17 runs splice-assay over it.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from _common import log

from ctsm.biology import bio_config, full_importance, task_table
from ctsm.config import PAN, load_config
from ctsm.matrix_io import decode, load_meta, view_dir

AS_VIEWS = ("afe", "ale", "hit", "se", "a5ss", "a3ss", "mxe", "ri")
RMATS = ("se", "mxe", "a3ss", "a5ss", "ri")
ENDPOINTS = {"OS": ("OS.time", "OS.event"), "DSS": ("DSS.time", "DSS.event")}


# ------------------------------------------------------------------ event geometry
def _iv(text: str) -> tuple[int, int]:
    a, b = text.split("-")
    return int(a), int(b)


def _fmt(ivs) -> str:
    return ";".join(f"{a}-{b}" for a, b in ivs)


def geometry(feature_id: str, view: str) -> tuple[str, str, str]:
    """(strand, constant, variable) as splice-assay reads them: 0-based half-open, ';'-joined.

    rMATS ids are already 0-based half-open; HITindex ids are 1-based inclusive (D-06), so their
    start moves back by one. For MXE the tracked exon (the one PSI measures) comes first: the
    transcript-upstream exon, lower-coordinate on +, higher on - (D-15).

    For A3SS/A5SS `variable` is the EXTENSION the long form adds, which is what splice-assay
    draws; the features table's `region` is the whole long exon instead. The two differ by
    design. For every other view the two agree exactly (checked over 200 events per view).
    """
    if view in ("afe", "ale", "hit"):
        _, locus, _ = feature_id.split(";")
        _, span = locus.rsplit(":", 1)
        a, b = _iv(span)
        return "", "", _fmt([(a - 1, b)])
    _, locus, rest = feature_id.split("#")
    _, first = locus.rsplit(":", 1)
    fields = rest.split(";")
    strand = fields[0]
    blocks = [_iv(first)] + [_iv(x) for x in fields[1:]]
    if view == "se":                                    # cassette; flanking; flanking
        const, var = sorted(blocks[1:]), [blocks[0]]
    elif view == "ri":                                  # envelope; flanking; flanking
        const = sorted(blocks[1:])
        var = [(const[0][1], const[1][0])]
    elif view in ("a3ss", "a5ss"):                      # long; short; flanking
        long_, short, flank = blocks
        const = sorted([short, flank])
        var = [(long_[0], short[0])] if long_[0] < short[0] else [(short[1], long_[1])]
    else:                                               # mxe: exon1; exon2; upstream; downstream
        e1, e2 = blocks[0], blocks[1]
        tracked, other = (e1, e2) if strand == "+" else (e2, e1)
        const, var = sorted(blocks[2:]), [tracked, other]
    return strand, _fmt(const), _fmt(var)


# ------------------------------------------------------------------ selection
def candidates(cfg, cancers, top_k: int) -> pd.DataFrame:
    """Top-k features of every tumour-status model, and the top-100 set to keep out of the background."""
    picked, seen_top100 = [], {v: set() for v in AS_VIEWS}
    for view in AS_VIEWS:
        for cancer in cancers:
            imp = full_importance(cfg, "tumor_status", view, cancer)
            seen_top100[view].update(imp["feature_id"].head(100))
            top = imp.head(top_k)[["feature_id", "mean_abs_shap"]].assign(view=view, cancer=cancer)
            top["rank"] = np.arange(1, len(top) + 1)
            picked.append(top)
    hits = pd.concat(picked, ignore_index=True)
    best = (hits.sort_values("rank").groupby(["view", "feature_id"], as_index=False)
            .agg(n_models=("cancer", "size"), best_cancer=("cancer", "first"), best_rank=("rank", "first")))
    return best, seen_top100


def tumour_stats(cfg, view: str, file_ids) -> pd.DataFrame:
    """Observed fraction and SD of every feature of a view over the given tumours (outcome-blind)."""
    d = view_dir(cfg.matrix_dir, view)
    feats = pd.read_parquet(d / "features.parquet", columns=["feature_id"])
    samples = pd.read_parquet(d / "samples.parquet", columns=["file_id"])
    rows = np.sort(pd.Series(np.arange(len(samples)), index=samples["file_id"])[list(file_ids)].to_numpy())
    mm = np.load(d / "values.npy", mmap_mode="r")
    enc = load_meta(cfg.matrix_dir, view)["encoding"]
    n, s1, s2 = (np.zeros(mm.shape[1]) for _ in range(3))
    for a in range(0, len(rows), 256):
        X = decode(np.asarray(mm[rows[a:a + 256]]), enc)
        obs = np.isfinite(X)
        n += obs.sum(0)
        s1 += np.nansum(np.where(obs, X, 0.0), axis=0)
        s2 += np.nansum(np.where(obs, X, 0.0) ** 2, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s1 / n
        sd = np.sqrt(np.maximum(s2 / n - mean ** 2, 0.0))
    return pd.DataFrame({"feature_id": feats["feature_id"], "frac_obs": n / len(rows), "sd": sd, "mean": mean})


def background(stats: pd.DataFrame, chosen: pd.Series, exclude: set, rng) -> pd.DataFrame:
    """One matched event per candidate: same view, same (frac_obs, SD) quintile, never model-selected."""
    s = stats[stats["frac_obs"] >= 0.5].copy()
    for col in ("frac_obs", "sd"):
        s[f"q_{col}"] = pd.qcut(s[col].rank(method="first"), 5, labels=False)
    s["stratum"] = s["q_frac_obs"].astype(str) + "_" + s["q_sd"].astype(str)
    pool = s[~s["feature_id"].isin(exclude)]
    want = s[s["feature_id"].isin(set(chosen))]["stratum"].value_counts()
    out = []
    for stratum, k in want.items():
        avail = pool[pool["stratum"] == stratum]
        if len(avail) == 0:
            continue
        take = avail.iloc[rng.choice(len(avail), size=min(k, len(avail)), replace=False)]
        out.append(take.assign(stratum=stratum))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=s.columns)


# ------------------------------------------------------------------ tables
def sample_table(cfg, cancers) -> pd.DataFrame:
    cohort = pd.read_parquet(cfg.data_dir / "cohort" / "samples.parquet")
    clin = pd.read_csv(cfg.tcga_path("clinical"), low_memory=False)
    keep = ["File.ID", "age_at_diagnosis", "gender", "stage"] + [c for p in ENDPOINTS.values() for c in p]
    s = cohort[cohort["project"].isin(cancers)].merge(clin[keep], left_on="file_id", right_on="File.ID", how="left")
    out = pd.DataFrame({"sample_id": s["file_id"], "patient_id": s["case_id"], "cohort": s["project"],
                        "group": s["role"], "age": s["age_at_diagnosis"], "sex": s["gender"], "stage": s["stage"]})
    for name, (t, e) in ENDPOINTS.items():
        out[f"{name}.time"], out[name] = s[t], s[e]
    both = s.groupby("case_id")["role"].nunique()
    out["pair_id"] = np.where(s["case_id"].map(both).eq(2), s["case_id"], "")
    # survival belongs to the patient, so only tumours carry it (a normal sample is not a second patient)
    out.loc[out["group"].eq("normal"), [f"{n}{suf}" for n in ENDPOINTS for suf in (".time", "")]] = np.nan
    return out


def gene_strands(cfg) -> pd.Series:
    """gene_id -> strand, from the rMATS views (HITindex ids carry no strand, and splice-assay
    needs one to draw an event)."""
    parts = [pd.read_parquet(cfg.matrix_dir / v / "features.parquet", columns=["gene_id", "strand"])
             for v in RMATS]
    return pd.concat(parts).dropna().drop_duplicates("gene_id").set_index("gene_id")["strand"]


def psi_table(cfg, selection: pd.DataFrame, file_ids) -> pd.DataFrame:
    """Event columns then one PSI column per sample, read straight from the stored matrices."""
    blocks = []
    strands = gene_strands(cfg)
    for view, part in selection.groupby("view", sort=False):
        d = view_dir(cfg.matrix_dir, view)
        feats = pd.read_parquet(d / "features.parquet")
        samples = pd.read_parquet(d / "samples.parquet", columns=["file_id"])
        rows = pd.Series(np.arange(len(samples)), index=samples["file_id"])[list(file_ids)].to_numpy()
        cols = pd.Series(np.arange(len(feats)), index=feats["feature_id"])[list(part["feature_id"])].to_numpy()
        mm = np.load(d / "values.npy", mmap_mode="r")
        enc = load_meta(cfg.matrix_dir, view)["encoding"]
        X = np.empty((len(rows), len(cols)), dtype=np.float32)
        order = np.argsort(rows)
        for a in range(0, len(rows), 256):
            o = order[a:a + 256]
            X[o] = decode(np.asarray(mm[rows[o]][:, cols]), enc)
        info = feats.set_index("feature_id").loc[part["feature_id"]].reset_index()
        geom = [geometry(f, view) for f in info["feature_id"]]
        fallback = info["gene_id"].map(strands).fillna("")
        strand = [g[0] or (s if isinstance(s, str) else "") or f
                  for g, s, f in zip(geom, info["strand"], fallback)]
        gene = np.where(info["symbol"].fillna("").str.len().gt(0), info["symbol"], info["gene_id"])
        blocks.append(pd.concat([
            pd.DataFrame({"event_id": view.upper() + ":" + info["feature_id"], "gene": gene,
                          "gene_id": info["gene_id"], "chrom": info["chrom"], "strand": strand,
                          "event_type": view.upper(), "constant": [g[1] for g in geom],
                          "variable": [g[2] for g in geom], "label": info["label"]}),
            pd.DataFrame(X.T, columns=list(file_ids))], axis=1))
        log(f"  {view}: {len(info)} events")
    return pd.concat(blocks, ignore_index=True)


def expression_table(cfg, gene_ids, file_ids, genes) -> pd.DataFrame:
    d = view_dir(cfg.matrix_dir, "gex")
    feats = pd.read_parquet(d / "features.parquet")
    samples = pd.read_parquet(d / "samples.parquet", columns=["file_id"])
    have = feats[feats["gene_id"].isin(set(gene_ids))].drop_duplicates("gene_id")
    rows = pd.Series(np.arange(len(samples)), index=samples["file_id"])[list(file_ids)].to_numpy()
    cols = pd.Series(np.arange(len(feats)), index=feats["feature_id"])[list(have["feature_id"])].to_numpy()
    mm = np.load(d / "values.npy", mmap_mode="r")
    enc = load_meta(cfg.matrix_dir, "gex")["encoding"]
    X = np.empty((len(rows), len(cols)), dtype=np.float32)
    order = np.argsort(rows)
    for a in range(0, len(rows), 256):
        o = order[a:a + 256]
        X[o] = decode(np.asarray(mm[rows[o]][:, cols]), enc)
    name = pd.Series(genes.values, index=genes.index)
    return pd.concat([pd.DataFrame({"gene": have["gene_id"].map(name).fillna(have["gene_id"]).to_numpy()}),
                      pd.DataFrame(X.T, columns=list(file_ids))], axis=1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-k", type=int, default=10, help="top features per (view, cancer)")
    args = ap.parse_args()
    cfg, bio = load_config(), bio_config()
    rng = np.random.default_rng(cfg.study["splits"]["seed"])
    scope = pd.read_csv(cfg.results_dir / "data" / "task_scope.csv")
    cancers = sorted(scope.loc[scope["tumor_status_in_scope"], "project"])
    out = cfg.data_dir / "survival"
    (out / "inputs").mkdir(parents=True, exist_ok=True)

    samples = sample_table(cfg, cancers)
    tumours = samples.loc[samples["group"].eq("tumour"), "sample_id"]
    log(f"{len(samples)} samples of {len(cancers)} cancers ({len(tumours)} tumours, "
        f"{int(samples['OS'].sum())} deaths)")

    cand, top100 = candidates(cfg, cancers, args.top_k)
    rows = []
    for view in AS_VIEWS:
        stats = tumour_stats(cfg, view, tumours)
        chosen = cand.loc[cand["view"].eq(view), "feature_id"]
        bg = background(stats, chosen, top100[view], rng)
        s = stats.set_index("feature_id")
        rows.append(cand[cand["view"].eq(view)].assign(role="candidate", view=view))
        rows.append(pd.DataFrame({"feature_id": bg["feature_id"], "view": view, "role": "background",
                                  "stratum": bg["stratum"]}))
        log(f"  {view}: {len(chosen)} candidates, {len(bg)} background "
            f"(median coverage {s.loc[chosen, 'frac_obs'].median():.2f} vs {bg['frac_obs'].median():.2f})")
    selection = pd.concat(rows, ignore_index=True)
    selection.to_parquet(out / "event_selection.parquet", index=False)

    psi = psi_table(cfg, selection, samples["sample_id"])
    genes = psi.drop_duplicates("gene_id").set_index("gene_id")["gene"]
    expr = expression_table(cfg, genes.index, samples["sample_id"], genes)
    samples.to_parquet(out / "inputs" / "samples.parquet", index=False)
    psi.drop(columns=["gene_id"]).to_parquet(out / "inputs" / "psi.parquet", index=False)
    expr.to_parquet(out / "inputs" / "expression.parquet", index=False)
    log(f"wrote {len(psi)} events x {len(samples)} samples and {len(expr)} host genes to {out / 'inputs'}")


if __name__ == "__main__":
    main()
