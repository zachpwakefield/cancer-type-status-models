"""On-disk layout and loaders for the processed molecular matrices.

<data_dir>/matrices/<view>/
    values.npy            samples x features ("int16_x1000": value = stored / 1000, -32768 = missing;
                          "float32": NaN = missing)
    samples.parquet       row order (file_id ...): every selected sample of every task
    features.parquet      annotation only: feature_id, gene_id, symbol, chrom, strand, region, label
    prefilter.parquet     one boolean column per (task cohort): superset-monotone extraction
                          pre-filter (D-07, proof in docs/LEAKAGE_POLICY.md)
    structural_zero.npy   AFE/ALE only: inferred zeros, np.packbits(axis=1, bitorder="little")
    meta.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

INT16_SCALE = 1000
INT16_MISSING = np.int16(-32768)


def encode_int16(values: np.ndarray) -> np.ndarray:
    """float matrix (NaN = missing) of exactly-3-decimal values -> int16 x 1000 (asserted lossless)."""
    scaled = values.astype(np.float64) * INT16_SCALE
    obs = np.isfinite(scaled)
    rounded = np.rint(scaled[obs])
    if rounded.size:
        err = np.abs(rounded - scaled[obs]).max()
        if err > 1e-2:
            raise ValueError(f"values are not 3-decimal: max scaled rounding error {err:.4g}")
        if rounded.min() < -32767 or rounded.max() > 32767:
            raise ValueError("values out of int16 range")
    out = np.full(values.shape, INT16_MISSING, dtype=np.int16)
    out[obs] = rounded.astype(np.int16)
    return out


def decode(stored: np.ndarray, encoding: str) -> np.ndarray:
    if encoding == "float32":
        return np.asarray(stored, dtype=np.float32)
    if encoding == "int16_x1000":
        stored = np.asarray(stored)
        out = stored.astype(np.float32) / INT16_SCALE
        out[stored == INT16_MISSING] = np.nan
        return out
    raise ValueError(encoding)


def view_dir(matrix_dir: Path, view: str) -> Path:
    return Path(matrix_dir) / view


def save_view(matrix_dir: Path, view: str, *, values: np.ndarray, samples: pd.DataFrame, features: pd.DataFrame,
              prefilter: pd.DataFrame, encoding: str, meta: dict, structural_zero: np.ndarray | None = None) -> list[Path]:
    d = view_dir(matrix_dir, view)
    d.mkdir(parents=True, exist_ok=True)
    if values.shape != (len(samples), len(features)) or len(prefilter) != len(features):
        raise ValueError(f"{view}: matrix {values.shape} does not match its tables")
    written = [d / "values.npy", d / "samples.parquet", d / "features.parquet", d / "prefilter.parquet"]
    np.save(written[0], np.ascontiguousarray(values))
    samples.reset_index(drop=True).to_parquet(written[1], index=False)
    features.reset_index(drop=True).to_parquet(written[2], index=False)
    prefilter.reset_index(drop=True).to_parquet(written[3], index=False)
    if structural_zero is not None:
        p = d / "structural_zero.npy"
        np.save(p, np.packbits(structural_zero, axis=1, bitorder="little"))
        written.append(p)
    meta = dict(meta, view=view, encoding=encoding, n_samples=int(values.shape[0]), n_features=int(values.shape[1]))
    (d / "meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True, default=str) + "\n")
    return written + [d / "meta.json"]


def load_meta(matrix_dir: Path, view: str) -> dict:
    return json.loads((view_dir(matrix_dir, view) / "meta.json").read_text())


class ViewMatrix:
    """Samples x features matrix kept in its stored encoding (int16 x 1000 stays 2 bytes per
    value) and decoded to float32 (NaN = missing) only for the block that is asked for."""

    def __init__(self, raw: np.ndarray, encoding: str):
        self.raw, self.encoding = raw, encoding
        self.shape = raw.shape

    def get(self, rows=None, cols=None) -> np.ndarray:
        """Decoded float32 block; rows / cols are index arrays or slices (None = all)."""
        r = slice(None) if rows is None else rows
        c = slice(None) if cols is None else cols
        if isinstance(c, slice):
            block = self.raw[:, c][r]
        else:
            block = self.raw[r][:, c] if isinstance(r, slice) else self.raw[np.ix_(np.asarray(r), np.asarray(c))]
        return decode(block, self.encoding)

    def missing_frac(self) -> float:
        if self.encoding == "int16_x1000":
            return float((self.raw == INT16_MISSING).mean())
        return float(np.isnan(self.raw).mean())


def load_view(matrix_dir: Path, view: str, file_ids, cohort: str | None = None):
    """Return (ViewMatrix [len(file_ids) x features], features DataFrame).

    Rows follow `file_ids`. With `cohort`, only features passing that cohort's extraction
    pre-filter are returned (a superset of anything a fold-level filter can select)."""
    d = view_dir(matrix_dir, view)
    meta = load_meta(matrix_dir, view)
    samples = pd.read_parquet(d / "samples.parquet")
    features = pd.read_parquet(d / "features.parquet")
    pos = pd.Series(np.arange(len(samples)), index=samples["file_id"])
    rows = pos.loc[list(file_ids)].to_numpy()
    cols = np.arange(len(features))
    if cohort is not None:
        cols = np.flatnonzero(pd.read_parquet(d / "prefilter.parquet")[cohort].to_numpy())
    mm = np.load(d / "values.npy", mmap_mode="r")
    if mm.shape != (len(samples), len(features)):
        raise ValueError(f"{view}: values.npy shape {mm.shape} does not match its tables")
    order = np.argsort(rows)                       # read rows in file order, then restore
    raw = np.empty((len(rows), len(cols)), dtype=mm.dtype)
    for a in range(0, len(rows), 1024):            # row chunks keep the temporary small
        o = order[a:a + 1024]
        raw[o] = np.asarray(mm[rows[o]])[:, cols]
    return ViewMatrix(raw, meta["encoding"]), features.iloc[cols].reset_index(drop=True)
