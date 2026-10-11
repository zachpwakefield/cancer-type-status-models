"""Streaming extraction of the raw TCGA tables into compact matrices (decisions D-06, D-07).

Each raw table is read once, in row blocks, keeping only the columns of selected samples.
Per block:
  1. parse IDs; genome scope + host-gene biotype (D-06)
  2. AFE/ALE: infer structural zeros within sample x gene, then keep only metaexons
     annotated as first (AFE) / last (ALE) exons (D-06)
  3. per task cohort, the superset-monotone pre-filter (D-07); keep rows passing any cohort
  4. encode losslessly and accumulate
GEX: log2(raw count + 1) -> raw counts -> log2(CPM + 1) per sample (D-05), then 1 and 3.
"""
from __future__ import annotations

from typing import Callable, Iterator

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv

from .annotation import load_metaexon_terminal_support, parse_ids, scope_mask
from .config import INT16_VIEWS, TERMINAL_VIEWS, Config
from .matrix_io import encode_int16, save_view


def stream_feature_csv(path, columns: list[str], dtype=np.float32,
                       block_size: int = 128 << 20) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield (ids, values[n_rows, len(columns)]) blocks; blank cells become NaN."""
    pa_type = pa.float64() if dtype == np.float64 else pa.float32()
    ro = pacsv.ReadOptions(block_size=block_size, use_threads=True)
    co = pacsv.ConvertOptions(include_columns=["id", *columns],
                              column_types={"id": pa.string(), **{c: pa_type for c in columns}})
    reader = pacsv.open_csv(path, read_options=ro, convert_options=co)
    for batch in reader:
        ids = np.asarray(batch.column(0).to_pylist(), dtype=object)
        mat = np.empty((batch.num_rows, len(columns)), dtype=dtype)
        for j in range(len(columns)):
            mat[:, j] = batch.column(j + 1).to_numpy(zero_copy_only=False)
        yield ids, mat


def regroup_by_gene(stream, gene_of: Callable[[str], str]) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Re-cut blocks so no gene is split across two blocks (rows must be gene-contiguous)."""
    carry_ids, carry_vals = None, None
    seen: set[str] = set()
    for ids, vals in stream:
        if carry_ids is not None:
            ids, vals = np.concatenate([carry_ids, ids]), np.concatenate([carry_vals, vals])
        genes = np.array([gene_of(i) for i in ids], dtype=object)
        tail = genes == genes[-1]
        cut = len(genes) - np.argmax(~tail[::-1]) if (~tail).any() else 0
        carry_ids, carry_vals = ids[cut:], vals[cut:]
        if cut:
            hg = set(genes[:cut])
            if hg & seen:
                raise ValueError("terminal-exon table is not gene-contiguous")
            seen |= hg
            yield ids[:cut], vals[:cut]
    if carry_ids is not None and len(carry_ids):
        yield carry_ids, carry_vals


def structural_zeros(values: np.ndarray, genes: np.ndarray, tol: float) -> tuple[np.ndarray, np.ndarray]:
    """Within each sample and gene, missing siblings become 0 iff at least one observed sibling
    is nonzero, observed siblings sum to 1 (+-tol) and at least one sibling is missing. Per
    sample (column) only, so it uses no other sample's data. Rows must be gene-contiguous."""
    change = np.flatnonzero(genes[1:] != genes[:-1]) + 1
    starts = np.concatenate([[0], change])
    sizes = np.diff(np.concatenate([starts, [len(genes)]]))
    obs = np.isfinite(values)
    filled = np.where(obs, values, 0.0)
    sums = np.add.reduceat(filled, starts, axis=0)
    n_obs = np.add.reduceat(obs.astype(np.int32), starts, axis=0)
    any_pos = np.add.reduceat((filled > 0).astype(np.int32), starts, axis=0) > 0
    ok = any_pos & (np.abs(sums - 1.0) <= tol) & (n_obs < sizes[:, None])
    fill = np.repeat(ok, sizes, axis=0) & ~obs
    out = values.copy()
    out[fill] = 0.0
    return out, fill


def prefilter_pass(vals: np.ndarray, cohorts: dict[str, np.ndarray], k_min: dict[str, int],
                   support: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """D-07 extraction pre-filter, per cohort: support count >= k_min[c] and not constant over
    the cohort's samples. `support` (rows x samples, bool) defaults to "value observed"."""
    obs = np.isfinite(vals)
    sup = obs if support is None else support
    out = {}
    for c, idx in cohorts.items():
        v, o = vals[:, idx], obs[:, idx]
        hi = np.where(o, v, -np.inf).max(axis=1)
        lo = np.where(o, v, np.inf).min(axis=1)
        out[c] = (sup[:, idx].sum(axis=1) >= k_min[c]) & (hi > lo)
    return out


def make_labels(feat: pd.DataFrame, view: str) -> pd.Series:
    """Human-readable, unique-within-view label: SYMBOL region (GEX: SYMBOL)."""
    base = feat["symbol"] if view == "gex" else feat["symbol"] + " " + feat["region"]
    dup = base.groupby(base).cumcount()
    return base.where(dup == 0, base + "#" + (dup + 1).astype(str))


def extract_as_view(cfg: Config, view: str, samples: pd.DataFrame, cohorts: dict[str, np.ndarray],
                    k_min: dict[str, int], ann: pd.DataFrame, log=print) -> tuple[list, dict]:
    terminal = view in TERMINAL_VIEWS
    encoding = "int16_x1000" if view in INT16_VIEWS else "float32"
    cols = list(samples["file_id"])
    tol = cfg.study["features"]["terminal"]["structural_zero_tolerance"]

    term_support = None
    if terminal and cfg.study["features"]["terminal"]["require_annotated_terminal"]:
        meta = load_metaexon_terminal_support(cfg.tcga_path("annotation", "hitindex_metaexons"))
        sup = meta["FE" if view == "afe" else "LE"] > 0
        term_support = set(meta.loc[sup, "gene_id"] + "|" + meta.loc[sup, "key"])

    acc_v, acc_f, acc_sz, acc_p = [], [], [], []
    counts = {"rows_raw": 0, "rows_in_scope": 0, "rows_terminal_annotated": 0, "structural_zeros_filled": 0,
              "rows_kept": 0}
    stream = stream_feature_csv(cfg.raw_matrix(view), cols, dtype=np.float64 if terminal else np.float32)
    if terminal:
        stream = regroup_by_gene(stream, lambda i: i.split(";", 1)[0])
    for b, (ids, vals) in enumerate(stream):
        counts["rows_raw"] += len(ids)
        feat = parse_ids(pd.Series(ids), view)
        keep = scope_mask(feat["gene_id"], feat["chrom"], ann, cfg)
        counts["rows_in_scope"] += int(keep.sum())
        sz = None
        if terminal:   # computed over every sibling row of the gene, before any row filter
            vals, sz = structural_zeros(vals, np.array([i.split(";", 1)[0] for i in ids], dtype=object), tol)
            counts["structural_zeros_filled"] += int(sz.sum())
            if term_support is not None:
                keep &= (feat["gene_id"] + "|" + feat["region"]).isin(term_support).to_numpy()
            counts["rows_terminal_annotated"] += int(keep.sum())
        if not keep.any():
            continue
        vals, feat = vals[keep], feat[keep].reset_index(drop=True)
        sz = sz[keep] if sz is not None else None
        passed = prefilter_pass(vals, cohorts, k_min)
        anyp = np.logical_or.reduce(list(passed.values()))
        if not anyp.any():
            continue
        counts["rows_kept"] += int(anyp.sum())
        acc_v.append(encode_int16(vals[anyp]) if encoding == "int16_x1000" else vals[anyp].astype(np.float32))
        acc_f.append(feat[anyp])
        acc_p.append(pd.DataFrame({c: p[anyp] for c, p in passed.items()}))
        if terminal:
            acc_sz.append(sz[anyp])
        if b % 10 == 0:
            log(f"[{view}] block {b}: {counts['rows_raw']:,} rows read, {counts['rows_kept']:,} kept")

    values = np.concatenate(acc_v, axis=0).T
    features = pd.concat(acc_f, ignore_index=True)
    if features["feature_id"].duplicated().any():
        raise ValueError(f"{view}: duplicate feature IDs")
    features["symbol"] = features["gene_id"].map(ann.set_index("gene_id")["symbol"]).fillna(features["gene_id"])
    features["label"] = make_labels(features, view)
    prefilter = pd.concat(acc_p, ignore_index=True)
    sz = np.concatenate(acc_sz, axis=0).T if terminal else None
    meta = {"source": str(cfg.raw_matrix(view)), "prefilter_k_min": k_min, **counts}
    written = save_view(cfg.matrix_dir, view, values=values, samples=samples[["file_id"]], features=features,
                        prefilter=prefilter, encoding=encoding, meta=meta, structural_zero=sz)
    stats = {"view": view, **counts, "n_features": int(values.shape[1]),
             **{f"prefilter_{c}": int(prefilter[c].sum()) for c in prefilter.columns},
             "missing_frac": round(float(np.mean(np.isnan(values)) if encoding == "float32"
                                         else np.mean(values == np.int16(-32768))), 4)}
    log(f"[{view}] done: {values.shape[0]:,} samples x {values.shape[1]:,} features")
    return written, stats


def extract_gex(cfg: Config, samples: pd.DataFrame, cohorts: dict[str, np.ndarray], k_min: dict[str, int],
                ann: pd.DataFrame, log=print) -> tuple[list, dict]:
    """gex.csv holds log2(raw STAR count + 1) for every file. Counts are recovered exactly
    (asserted integer), library size = total counts over all 57,656 genes of the sample, and
    the stored value is log2(CPM + 1). Per sample only, so it uses no other sample's data."""
    cols = list(samples["file_id"])
    ids_acc, cnt_acc = [], []
    for b, (ids, vals) in enumerate(stream_feature_csv(cfg.raw_matrix("gex"), cols, dtype=np.float64)):
        if not np.isfinite(vals).all():
            raise ValueError("gex.csv has missing values")
        counts = np.exp2(vals) - 1.0
        r = np.rint(counts)
        if np.abs(counts - r).max() > 1e-6 * max(1.0, r.max()) or (r < 0).any():
            raise ValueError("gex.csv values are not log2(integer count + 1)")
        ids_acc.append(ids)
        cnt_acc.append(r.astype(np.float32))
        log(f"[gex] block {b}: {sum(map(len, ids_acc)):,} genes read")
    ids = np.concatenate(ids_acc)
    C = np.concatenate(cnt_acc, axis=0)                     # genes x samples
    lib = C.astype(np.float64).sum(axis=0)
    if (lib <= 0).any():
        raise ValueError("sample with zero library size")
    feat = pd.DataFrame({"feature_id": ids, "gene_id": ids})
    feat["chrom"] = feat["gene_id"].map(ann.set_index("gene_id")["chrom"])
    keep = scope_mask(feat["gene_id"], feat["chrom"], ann, cfg)
    L = np.log2(C[keep] / lib[None, :].astype(np.float32) * 1e6 + 1.0, dtype=np.float32)
    feat = feat[keep].reset_index(drop=True)
    expressed = L >= cfg.study["preprocess"]["gex_min_log2cpm"]
    passed = prefilter_pass(L, cohorts, k_min, support=expressed)
    anyp = np.logical_or.reduce(list(passed.values()))
    feat = feat[anyp].reset_index(drop=True)
    feat["strand"] = None
    feat["region"] = None
    feat["symbol"] = feat["gene_id"].map(ann.set_index("gene_id")["symbol"]).fillna(feat["gene_id"])
    feat["label"] = make_labels(feat, "gex")
    prefilter = pd.DataFrame({c: p[anyp] for c, p in passed.items()})
    values = L[anyp].T
    meta = {"source": str(cfg.raw_matrix("gex")), "transform": "log2(CPM + 1), CPM over all genes of gex.csv",
            "prefilter_k_min": k_min, "genes_raw": int(len(ids)), "genes_in_scope": int(keep.sum()),
            "library_size_min": float(lib.min()), "library_size_median": float(np.median(lib))}
    written = save_view(cfg.matrix_dir, "gex", values=values, samples=samples[["file_id"]].assign(library_size=lib),
                        features=feat, prefilter=prefilter, encoding="float32", meta=meta)
    stats = {"view": "gex", "rows_raw": int(len(ids)), "rows_in_scope": int(keep.sum()), "rows_kept": int(anyp.sum()),
             "n_features": int(values.shape[1]), **{f"prefilter_{c}": int(prefilter[c].sum()) for c in prefilter.columns},
             "missing_frac": 0.0}
    log(f"[gex] done: {values.shape[0]:,} samples x {values.shape[1]:,} genes")
    return written, stats
