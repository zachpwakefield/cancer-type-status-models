"""Static figures (PNG) for the results tree.

Colour roles: sequential = one blue ramp (light -> dark) for magnitudes such as AUROC, recall
and mean |SHAP|; diverging = blue (low) <-> grey <-> red (high) for standardised feature values
in SHAP beeswarms; one series = slot-1 blue. Text stays in neutral ink.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div_blue_red", ["#1c5cab", "#86b6ef", "#e1e0dc", "#f0a3a2", "#c22f2e"])
VIEW_NAMES = {"afe": "AFE", "ale": "ALE", "hit": "HIT", "se": "SE", "a5ss": "A5SS", "a3ss": "A3SS",
              "mxe": "MXE", "ri": "RI", "gex": "GEX"}


def _style(ax, grid_axis: str | None = "x") -> None:
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)


def _save(fig, path) -> None:
    fig.patch.set_facecolor(SURFACE)
    fig.savefig(path, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def beeswarm(shap_vals: list[np.ndarray], values: list[np.ndarray], labels: list[str], title: str, xlabel: str,
             path, seed: int = 0) -> None:
    """One row per feature (top first): each dot is a held-out sample, x = SHAP value,
    colour = standardised feature value (clipped to +-2.5)."""
    n = len(labels)
    fig, ax = plt.subplots(figsize=(7.2, 0.32 * n + 1.3))
    rng = np.random.default_rng(seed)
    for i, (s, v) in enumerate(zip(shap_vals, values)):
        yrow = n - 1 - i
        # density-scaled vertical jitter so dense regions spread out
        if len(s) > 1:
            hist, edges = np.histogram(s, bins=40)
            dens = hist[np.clip(np.digitize(s, edges) - 1, 0, 39)] / max(hist.max(), 1)
        else:
            dens = np.ones_like(s)
        jit = (rng.random(len(s)) - 0.5) * 0.7 * np.sqrt(dens)
        order = rng.permutation(len(s))
        ax.scatter(s[order], (yrow + jit)[order], c=np.clip(v[order], -2.5, 2.5), cmap=DIV, vmin=-2.5, vmax=2.5,
                   s=7, linewidths=0, alpha=0.85, rasterized=True)
    ax.axvline(0, color=AXIS, linewidth=0.8, zorder=0)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels[::-1], fontsize=7.5, color=INK)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel(xlabel, fontsize=8.5, color=INK2)
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    _style(ax, "x")
    sm = plt.cm.ScalarMappable(cmap=DIV, norm=plt.Normalize(-2.5, 2.5))
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, ticks=[-2.5, 0, 2.5])
    cb.ax.set_yticklabels(["low", "", "high"], fontsize=7, color=INK2)
    cb.set_label("feature value (z)", fontsize=7.5, color=INK2)
    cb.outline.set_visible(False)
    _save(fig, path)


def heatmap(mat: pd.DataFrame, title: str, path, vmin: float | None = None, vmax: float | None = None,
            fmt: str | None = "{:.3f}", cbar_label: str = "", cell: float = 0.42, xrot: int = 0) -> None:
    """Sequential one-hue heatmap; optional numbers in cells (ink flips to white on dark cells)."""
    nr, nc = mat.shape
    fig, ax = plt.subplots(figsize=(max(4.0, cell * nc + 2.6), max(2.5, cell * nr + 1.4)))
    vals = mat.to_numpy(dtype=float)
    lo = np.nanmin(vals) if vmin is None else vmin
    hi = np.nanmax(vals) if vmax is None else vmax
    im = ax.imshow(np.ma.masked_invalid(vals), cmap=SEQ, vmin=lo, vmax=hi, aspect="auto")
    ax.set_xticks(range(nc))
    ax.set_xticklabels(mat.columns, fontsize=8, color=INK, rotation=xrot, ha="right" if xrot else "center")
    ax.set_yticks(range(nr))
    ax.set_yticklabels(mat.index, fontsize=8, color=INK)
    ax.set_xticks(np.arange(-0.5, nc, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nr, 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=1.5)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if fmt:
        for i in range(nr):
            for j in range(nc):
                v = vals[i, j]
                if np.isfinite(v):
                    dark = (v - lo) / max(hi - lo, 1e-12) > 0.6
                    ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=6.5,
                            color="#ffffff" if dark else INK)
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.ax.tick_params(labelsize=7, colors=INK2, length=0)
    cb.outline.set_visible(False)
    if cbar_label:
        cb.set_label(cbar_label, fontsize=7.5, color=INK2)
    _save(fig, path)


def dot_metrics(df: pd.DataFrame, metrics: list[str], names: list[str], title: str, path,
                xlim: tuple[float, float] | None = None) -> None:
    """Rows = views; one dot per metric (categorical slots 1..3) with 95 % CI whiskers
    when df has <metric>_lo / <metric>_hi columns."""
    n = len(df)
    fig, ax = plt.subplots(figsize=(6.8, 0.36 * n + 1.5))
    offs = np.linspace(-0.18, 0.18, len(metrics))
    for m, nm, col, off in zip(metrics, names, SERIES, offs):
        y = np.arange(n)[::-1] + off
        if f"{m}_lo" in df:
            ax.hlines(y, df[f"{m}_lo"], df[f"{m}_hi"], color=col, linewidth=1.2, alpha=0.6)
        ax.scatter(df[m], y, s=30, color=col, edgecolor=SURFACE, linewidth=1.5, zorder=3, label=nm)
    ax.set_yticks(np.arange(n)[::-1])
    ax.set_yticklabels(df.index, fontsize=8.5, color=INK)
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_title(title, fontsize=10, color=INK, loc="left", pad=26)
    _style(ax, "x")
    ax.legend(frameon=False, fontsize=8, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=len(metrics),
              labelcolor=INK2, handletextpad=0.2, borderaxespad=0.2)
    _save(fig, path)
