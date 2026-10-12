"""House figure style (docs/FIGURE_GUIDELINES.md).

The rcParams follow the shared `figstyle.py` and the palette follows `fig_style_snhg.py`, both
from the lab's figure guidelines (survival_prediction/docs/figure_guidelines/, 2026-10-11):
Arial 6.5-8 pt, 0.6 pt axes without top/right spines, no grid, no legend frames, editable
text in SVG/PDF, restrained colour, no in-plot titles, every final figure exported as SVG, PDF
and a 400-dpi PNG with a CSV of the plotted values and a caption file.

Colour roles in this repo
  TUMOUR / NORMAL   tumour-normal figures only (blue / orange, as everywhere in the lab)
  ACCENT            the point of a non-tumour-normal figure (splicing, the AS ensemble)
  ACCENT_DARK       a second accent (fusion, "observed")
  INK2 / MUTED      reference and context (GEX, single views, chance, background)
  DIVERGE           signed changes with no tumour-normal meaning (blue below, red above)
  TN_DIVERGE        a tumour - normal difference, in the house tumour/normal hues
                    (orange = higher in normal, blue = higher in tumour)
  SEQ               magnitudes (paper -> dark blue)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# palette of fig_style_snhg.py
INK, INK2, MUTED, GRID, PAPER = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#f4f3ee"
TUMOUR, NORMAL = "#2a78d6", "#eb6834"
ACCENT, ACCENT_DARK, ACCENT_LIGHT = "#2a78d6", "#184f95", "#86b6ef"
DOWN, NULLC = "#a8a69f", "#b9b6ab"
DIVERGE = LinearSegmentedColormap.from_list("div", ["#2f5f98", "#9ebbd9", PAPER, "#e6a88a", "#b0412c"])
# A tumour - normal difference is read with the house tumour/normal colours, so that blue means
# "higher in tumour" here exactly as it does on every other tumour-normal figure in the lab.
TN_DIVERGE = LinearSegmentedColormap.from_list("tn", ["#b4461a", "#f0a480", PAPER, "#8fb9e8", "#1b538f"])
SEQ = LinearSegmentedColormap.from_list("seq", [PAPER, "#c3d6ee", "#6c9bd2", "#184f95"])   # TIER blues
VIEW_NAMES = {"afe": "AFE", "ale": "ALE", "hit": "HIT", "se": "SE", "a5ss": "A5SS", "a3ss": "A3SS",
              "mxe": "MXE", "ri": "RI", "gex": "GEX"}

MM = 1 / 25.4
SINGLE, DOUBLE = 89 * MM, 183 * MM       # journal column widths, inches


def apply_style() -> None:
    """Global rcParams (shared figstyle.py, with the ink colours of fig_style_snhg.py)."""
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7, "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "legend.frameon": False,
        "figure.dpi": 150, "savefig.dpi": 400, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "image.interpolation": "nearest",
    })


def fp(p) -> str:
    """p to print: 3 decimals below 0.1, 2 above, 1.2e-5 style below 0.001; '–' if missing."""
    if p is None or not np.isfinite(p):
        return "–"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m}e{int(e)}"
    return f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


def stars(p) -> str:
    """Raw-p tiers for dense displays (the caption states the tiers): * < 0.05, ** < 0.01, *** < 0.001."""
    if p is None or not np.isfinite(p):
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else ""


def save(fig, out_dir: Path, stem: str, data: pd.DataFrame | None = None, caption: str | None = None,
         formats: tuple[str, ...] = ("svg", "pdf", "png"), local_data: tuple[Path, pd.DataFrame] | None = None
         ) -> Path:
    """Write the figure (SVG, PDF, 400-dpi PNG), `<stem>.csv` with the plotted values and
    `<stem>.caption.md`. Per-sample plotted values are individual-level data (D-12): pass them as
    `local_data=(dir, frame)` so they go to the git-ignored data directory instead."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(out_dir / f"{stem}.{ext}")
    plt.close(fig)
    if data is not None:
        data.to_csv(out_dir / f"{stem}.csv", index=False)
    if caption is not None:
        (out_dir / f"{stem}.caption.md").write_text(caption.strip() + "\n")
    if local_data is not None:
        d, frame = local_data
        Path(d).mkdir(parents=True, exist_ok=True)
        frame.to_csv(Path(d) / f"{stem}.csv", index=False)
    return out_dir / f"{stem}.png"


def heat_axes(ax, rows: list[str], cols: list[str], xrot: int = 0, row_size: float = 6.5,
              col_size: float = 6.5) -> None:
    """Tick labels for an imshow heatmap; thin paper-coloured cell separators, no spines."""
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=xrot, ha="right" if xrot and xrot != 90 else "center", fontsize=col_size,
                       color=INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=row_size, color=INK)
    ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.grid(which="major", visible=False)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)


def colourbar(fig, im, ax, label: str, ticks=None, shrink: float = 0.6, extend: str = "neither"):
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, shrink=shrink, ticks=ticks, extend=extend)
    cb.set_label(label, fontsize=6.5, color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=6, colors=INK2, length=0)
    return cb


def bracket(ax, x0: float, x1: float, y: float, text: str, h: float, fontsize: float = 6.5) -> None:
    """Short comparison bracket with its label (use alternating heights for neighbours)."""
    ax.plot([x0, x0, x1, x1], [y, y + h, y + h, y], color=INK2, linewidth=0.6, clip_on=False)
    ax.text((x0 + x1) / 2, y + h * 1.4, text, ha="center", va="bottom", fontsize=fontsize, color=INK)


def gene_label(symbol) -> str:
    """HGNC symbol for display; genes without one (Ensembl ID or missing) read 'unnamed gene'.
    The Ensembl IDs stay in the CSVs."""
    s = "" if symbol is None or (isinstance(symbol, float) and np.isnan(symbol)) else str(symbol)
    return "unnamed gene" if not s or s.startswith("ENSG") else s


def feature_label(label: str) -> str:
    """'<symbol or ENSG> <coordinates>' -> display label with gene_label()."""
    g, _, rest = str(label).partition(" ")
    return f"{gene_label(g)} {rest}".strip()
