"""Step 03: validate the cohort, splits and matrices against the raw sources.

Writes results/data/VALIDATION.md (and validation.json). Fails loudly on any violated check.
  V1 cohort      one file per (case, role); normals carry barcode type 11; task tables agree
  V2 splits      case-grouped, nested, every sample tested exactly once (re-checked)
  V3 raw values  stored values equal the raw CSV for every stored feature in the first raw
                 block x 300 random samples (AFE/ALE: equal, or an inferred structural zero)
  V4 GEX         counts recovered exactly; log2(CPM+1) vs official GDC log2(TPM+1) per sample
  V5 pre-filter  recorded support thresholds equal ceil(fraction x smallest inner-train set)
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
from _common import log
from scipy.stats import spearmanr

from ctsm.config import TERMINAL_VIEWS, load_config
from ctsm.extract import stream_feature_csv
from ctsm.matrix_io import decode, load_meta, view_dir
from ctsm.splits import check_splits, min_inner_train


def v1_cohort(cfg, samples) -> dict:
    assert not samples.duplicated(["case_id", "role"]).any(), "more than one file per case and role"
    assert samples.file_id.is_unique
    assert (samples.loc[samples.role == "normal", "sample_id"].str[13:15] == "11").all()
    out = {"samples": len(samples), "cases": int(samples.case_id.nunique()),
           "tumour": int((samples.role == "tumour").sum()), "normal": int((samples.role == "normal").sum())}
    for p in sorted((cfg.data_dir / "tasks").glob("*.parquet")):
        t = pd.read_parquet(p)
        assert set(t.file_id) <= set(samples.file_id)
    for v in cfg.views:
        s = pd.read_parquet(view_dir(cfg.matrix_dir, v) / "samples.parquet")
        assert s.file_id.tolist() == samples.file_id.tolist(), f"{v}: matrix rows differ from the cohort"
    return out


def v2_splits(cfg) -> dict:
    sp = cfg.study["splits"]
    out = {}
    for p in sorted((cfg.data_dir / "splits").glob("*.parquet")):
        s = pd.read_parquet(p)
        check_splits(s, sp["outer_folds"], sp["inner_folds"])
        out[p.stem] = {"n": len(s), "min_inner_train": min_inner_train(s, sp["outer_folds"], sp["inner_folds"])}
    return out


def v3_raw(cfg, view, samples, rng) -> dict:
    d = view_dir(cfg.matrix_dir, view)
    meta = load_meta(cfg.matrix_dir, view)
    feats = pd.read_parquet(d / "features.parquet")
    pick = np.sort(rng.choice(len(samples), 300, replace=False))
    cols = samples.file_id.to_numpy()[pick].tolist()
    ids, raw = next(iter(stream_feature_csv(cfg.raw_matrix(view), cols, dtype=np.float64, block_size=64 << 20)))
    fpos = pd.Series(np.arange(len(feats)), index=feats.feature_id)
    have = np.flatnonzero(pd.Series(ids).isin(fpos.index).to_numpy())
    if not len(have):
        return {"features_checked": 0}
    mm = np.load(d / "values.npy", mmap_mode="r")
    j = fpos.loc[ids[have]].to_numpy()
    stored = decode(np.asarray(mm[pick][:, j]), meta["encoding"]).T          # features x samples
    r = raw[have]
    if view == "gex":
        lib = pd.read_parquet(d / "samples.parquet")["library_size"].to_numpy()[pick]
        expect = np.log2((np.rint(np.exp2(r) - 1.0)) / lib[None, :] * 1e6 + 1.0)
        err = float(np.abs(stored - expect).max())
        assert err < 1e-4, f"gex: max deviation {err}"
        return {"features_checked": int(len(have)), "max_abs_error": err}
    both = np.isfinite(r) & np.isfinite(stored)
    err = float(np.abs(stored[both] - r[both]).max()) if both.any() else 0.0
    assert err < 1e-5, f"{view}: max deviation {err}"
    raw_only_missing = np.isnan(r) & np.isfinite(stored)
    if view in TERMINAL_VIEWS:
        sz = np.unpackbits(np.load(d / "structural_zero.npy", mmap_mode="r")[pick], axis=1,
                           count=meta["n_features"], bitorder="little")[:, j].T.astype(bool)
        assert (stored[raw_only_missing] == 0).all() and sz[raw_only_missing].all(), f"{view}: unexplained fills"
    else:
        assert not raw_only_missing.any(), f"{view}: value stored where raw is blank"
    assert not (np.isfinite(r) & np.isnan(stored)).any(), f"{view}: raw value lost"
    return {"features_checked": int(len(have)), "cells_checked": int(r.size), "max_abs_error": err,
            "structural_zero_fills_checked": int(raw_only_missing.sum())}


def v4_gex(cfg, samples, rng) -> dict:
    """Sanity of the depth normalisation against official TPM (primary tumours only; TPM also
    corrects for gene length, so agreement is high but not exact)."""
    tpm = np.load(cfg.tcga_path("gdc_tpm_matrix"), mmap_mode="r")
    ts = pd.read_csv(cfg.tcga_path("gdc_tpm_samples"), usecols=["sample_index", "file_id"])
    genes = pd.read_csv(cfg.tcga_path("annotation", "gdc_genes_v36"), usecols=["gene_index", "gene_id_unversioned"])
    d = view_dir(cfg.matrix_dir, "gex")
    feats = pd.read_parquet(d / "features.parquet")
    gpos = genes.drop_duplicates("gene_id_unversioned").set_index("gene_id_unversioned")["gene_index"]
    common = feats.gene_id[feats.gene_id.isin(gpos.index)]
    both = samples[samples.file_id.isin(ts.file_id)]
    pick = both.iloc[np.sort(rng.choice(len(both), min(300, len(both)), replace=False))]
    mm = np.load(d / "values.npy", mmap_mode="r")
    srow = pd.Series(np.arange(len(samples)), index=samples.file_id)
    trow = ts.set_index("file_id")["sample_index"]
    fcol = pd.Series(np.arange(len(feats)), index=feats.gene_id)
    rhos = []
    for fid in pick.file_id:
        a = np.asarray(mm[srow[fid]][fcol.loc[common].to_numpy()])
        b = np.log2(np.asarray(tpm[trow[fid]][gpos.loc[common].to_numpy()]) + 1)
        rhos.append(spearmanr(a, b).statistic)
    rhos = np.array(rhos)
    assert np.median(rhos) > 0.8, f"gex: median Spearman vs TPM {np.median(rhos):.3f}"
    return {"samples": int(len(rhos)), "genes": int(len(common)), "spearman_median": float(np.median(rhos)),
            "spearman_min": float(rhos.min())}


def v5_prefilter(cfg) -> dict:
    sp, pp = cfg.study["splits"], cfg.study["preprocess"]
    n_ref = {p.stem.split("__", 1)[1]: min_inner_train(pd.read_parquet(p), sp["outer_folds"], sp["inner_folds"])
             for p in (cfg.data_dir / "splits").glob("*.parquet")}
    for v in cfg.views:
        frac = pp["gex_min_expressed_frac"] if v == "gex" else pp["min_observed_frac"]
        k = load_meta(cfg.matrix_dir, v)["prefilter_k_min"]
        for c, n in n_ref.items():
            assert k[c] == math.ceil(frac * n), f"{v}/{c}: pre-filter threshold {k[c]} != ceil({frac} x {n})"
    return {"cohorts": len(n_ref), "views": len(cfg.views)}


def main() -> None:
    cfg = load_config()
    rng = np.random.default_rng(20261010)
    samples = pd.read_parquet(cfg.data_dir / "cohort" / "samples.parquet")
    res = {"V1_cohort": v1_cohort(cfg, samples), "V2_splits": v2_splits(cfg)}
    log("V1, V2 passed")
    res["V3_raw_values"] = {}
    for v in cfg.views:
        res["V3_raw_values"][v] = v3_raw(cfg, v, samples, rng)
        log(f"V3 {v}: {res['V3_raw_values'][v]}")
    res["V4_gex_vs_tpm"] = v4_gex(cfg, samples, rng)
    log(f"V4: {res['V4_gex_vs_tpm']}")
    res["V5_prefilter"] = v5_prefilter(cfg)
    out = cfg.results_dir / "data"
    (out / "validation.json").write_text(json.dumps(res, indent=2, default=str) + "\n")
    v3 = pd.DataFrame(res["V3_raw_values"]).T
    v4 = res["V4_gex_vs_tpm"]
    md = ["# Data validation", "", "All checks passed (`scripts/03_validate_data.py`; numbers in `validation.json`).", "",
          "| Check | Result |", "|---|---|",
          f"| V1 cohort | {res['V1_cohort']['samples']:,} samples from {res['V1_cohort']['cases']:,} cases "
          f"({res['V1_cohort']['tumour']:,} tumour, {res['V1_cohort']['normal']:,} normal); one file per case and role; "
          f"every view's matrix rows equal the cohort |",
          f"| V2 splits | {len(res['V2_splits'])} task cohorts: case-grouped outer and inner folds, every sample tested once |",
          f"| V3 raw values | stored = raw CSV for {int(v3.features_checked.sum()):,} features x 300 samples across 9 views "
          f"(max abs error {v3.max_abs_error.max():.1e}); AFE/ALE fills are all flagged structural zeros |",
          f"| V4 GEX | log2(CPM+1) vs official GDC log2(TPM+1), {v4['samples']} primary tumours x {v4['genes']:,} genes: "
          f"median Spearman {v4['spearman_median']:.3f} (min {v4['spearman_min']:.3f}) |",
          f"| V5 pre-filter | thresholds = ceil(fraction x smallest inner-train set) for all {res['V5_prefilter']['cohorts']} "
          f"cohorts x {res['V5_prefilter']['views']} views |", ""]
    (out / "VALIDATION.md").write_text("\n".join(md))
    log("all validation checks passed")


if __name__ == "__main__":
    main()
