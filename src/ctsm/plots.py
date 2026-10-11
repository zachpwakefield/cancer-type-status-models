"""Phase-1 figures in the house style (ctsm.figstyle, docs/FIGURE_GUIDELINES.md).

No in-plot titles: what a figure shows is written to `<stem>.caption.md` next to it. Summary
figures are exported as SVG, PDF and a 400-dpi PNG with a CSV of the plotted values; the
per-job panels (one per cohort and view) are PNG only, to keep the repository small, and their
plotted values are the job's SHAP and metrics tables (per-sample values stay local, D-12).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .figstyle import (ACCENT, DIVERGE, INK, MUTED, SEQ, VIEW_NAMES, apply_style, colourbar, heat_axes, save)

apply_style()
__all__ = ["VIEW_NAMES", "beeswarm", "heatmap", "dot_metrics"]


def beeswarm(shap_vals: list[np.ndarray], values: list[np.ndarray], labels: list[str], xlabel: str, out_dir: Path,
             stem: str, caption: str, seed: int = 0) -> Path:
    """One row per feature (top first): each dot is a held-out sample, x = SHAP value,
    colour = standardised feature value (clipped to +-2.5)."""
    n = len(labels)
    fig, ax = plt.subplots(figsize=(3.6, 0.17 * n + 0.6))
    rng = np.random.default_rng(seed)
    for i, (s, v) in enumerate(zip(shap_vals, values)):
        yrow = n - 1 - i
        if len(s) > 1:                                   # density-scaled vertical jitter
            hist, edges = np.histogram(s, bins=40)
            dens = hist[np.clip(np.digitize(s, edges) - 1, 0, 39)] / max(hist.max(), 1)
        else:
            dens = np.ones_like(s)
        jit = (rng.random(len(s)) - 0.5) * 0.7 * np.sqrt(dens)
        order = rng.permutation(len(s))
        ax.scatter(s[order], (yrow + jit)[order], c=np.clip(v[order], -2.5, 2.5), cmap=DIVERGE, vmin=-2.5,
                   vmax=2.5, s=1.6, linewidths=0, alpha=0.9, rasterized=True)
    ax.axvline(0, color=MUTED, linewidth=0.5, zorder=0)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels[::-1], fontsize=6, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel(xlabel)
    sm = plt.cm.ScalarMappable(cmap=DIVERGE, norm=plt.Normalize(-2.5, 2.5))
    cb = colourbar(fig, sm, ax, "feature value (z)", ticks=[-2.5, 0, 2.5], shrink=0.5)
    cb.ax.set_yticklabels(["low", "0", "high"])
    return save(fig, out_dir, stem, caption=caption, formats=("png",))


def heatmap(mat: pd.DataFrame, out_dir: Path, stem: str, caption: str, vmin: float | None = None,
            vmax: float | None = None, fmt: str | None = "{:.3f}", cbar_label: str = "", cell: float = 0.22,
            xrot: int = 0, formats: tuple[str, ...] = ("svg", "pdf", "png"), label_size: float = 6.5,
            extend: str = "neither") -> Path:
    """Sequential heatmap; optional values in cells (white text on dark cells). With vmin above the
    data minimum, pass extend="min" so the colour bar shows that lower values are clipped."""
    nr, nc = mat.shape
    fig, ax = plt.subplots(figsize=(cell * nc + 1.4, cell * nr + 0.6))
    vals = mat.to_numpy(dtype=float)
    lo = np.nanmin(vals) if vmin is None else vmin
    hi = np.nanmax(vals) if vmax is None else vmax
    im = ax.imshow(np.ma.masked_invalid(vals), cmap=SEQ, vmin=lo, vmax=hi, aspect="auto")
    heat_axes(ax, [str(x) for x in mat.index], [str(x) for x in mat.columns], xrot=xrot, row_size=label_size,
              col_size=label_size)
    if fmt:
        for i in range(nr):
            for j in range(nc):
                v = vals[i, j]
                if np.isfinite(v):
                    dark = (min(max(v, lo), hi) - lo) / max(hi - lo, 1e-12) > 0.6
                    ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=5.5,
                            color="white" if dark else INK)
    colourbar(fig, im, ax, cbar_label, extend=extend)
    data = mat.reset_index().rename(columns={"index": mat.index.name or "row"})
    return save(fig, out_dir, stem, data=data, caption=caption, formats=formats)


def dot_metrics(df: pd.DataFrame, metrics: list[str], names: list[str], out_dir: Path, stem: str, caption: str,
                xlim: tuple[float, float] | None = None) -> Path:
    """Rows = views; one dot per metric with 95 % CI whiskers (<metric>_lo / _hi columns)."""
    n = len(df)
    cols = [INK, ACCENT, MUTED]
    fig, ax = plt.subplots(figsize=(3.4, 0.2 * n + 0.7))
    offs = np.linspace(-0.2, 0.2, len(metrics))
    for m, nm, col, off in zip(metrics, names, cols, offs):
        y = np.arange(n)[::-1] + off
        if f"{m}_lo" in df:
            ax.hlines(y, df[f"{m}_lo"], df[f"{m}_hi"], color=col, linewidth=0.8)
        ax.scatter(df[m], y, s=9, color=col, linewidths=0, zorder=3, label=nm)
    ax.set_yticks(np.arange(n)[::-1])
    ax.set_yticklabels(df.index, color=INK)
    ax.tick_params(axis="y", length=0)
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_xlabel("out-of-fold score (95% case-level bootstrap CI)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=len(metrics), handletextpad=0.1, columnspacing=0.8,
              borderaxespad=0.1)
    keep = [c for m in metrics for c in (m, f"{m}_lo", f"{m}_hi") if c in df]
    return save(fig, out_dir, stem, data=df[keep].reset_index().rename(columns={"index": "view"}), caption=caption)
