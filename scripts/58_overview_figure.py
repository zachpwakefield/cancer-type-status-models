"""Step 58: the study-overview figure of the README (O1).

    O1  (a) the nine views, drawn as the two isoforms each splicing view compares; (b) the samples of
        each TCGA project and which analyses use them

Inputs:  results/data/task_scope.csv, results/subtype/summary/performance.csv, and the per-project
         count of MC3 MMR-deficient tumours from data.nosync/ctsm_v1/msi/models_tumours.parquet
         (step 54; only the counts are written)
Outputs: results/figures/overview/O1_study_design.{svg,pdf,png,csv,caption.md}
Usage: python scripts/58_overview_figure.py
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _figbase import panel_letters, tidy
from matplotlib.patches import Rectangle

from ctsm.config import load_config
from ctsm.figstyle import ACCENT, ACCENT_DARK, DOUBLE, INK, INK2, MUTED, NORMAL, TUMOUR, save

CFG = load_config()
OUT = CFG.results_dir / "figures" / "overview"
TRAIN = ["COAD", "STAD", "UCEC"]

# Each splicing view as two isoforms: exons as (start, end, alternative?), x in arbitrary units.
# Row 0 is the form whose use the view measures, row 1 the alternative form.
EXON = 1.0
VIEWS = {
    "AFE":  [[(0, 1, True), (5, 6, False), (8, 9, False)], [(2.5, 3.5, True), (5, 6, False), (8, 9, False)]],
    "ALE":  [[(0, 1, False), (3, 4, False), (5.5, 6.5, True)], [(0, 1, False), (3, 4, False), (8, 9, True)]],
    "HIT":  [[(0, 1, False), (4, 5, True)], [(4, 5, True), (8, 9, False)]],
    "SE":   [[(0, 1, False), (4, 5, True), (8, 9, False)], [(0, 1, False), (8, 9, False)]],
    "A5SS": [[(0, 1.8, False), (8, 9, False)], [(0, 1, False), (8, 9, False)]],
    "A3SS": [[(0, 1, False), (7.2, 9, False)], [(0, 1, False), (8, 9, False)]],
    "MXE":  [[(0, 1, False), (3, 4, True), (8, 9, False)], [(0, 1, False), (5.5, 6.5, True), (8, 9, False)]],
    "RI":   [[(0, 9, False)], [(0, 1, False), (8, 9, False)]],
}
EXTENSION = {"A5SS": (1, 1.8), "A3SS": (7.2, 8), "RI": (1, 8)}     # the segment that differs


def _isoform(ax, y, exons, view, rows, k):
    xs = [e for ex in exons for e in ex[:2]]
    ax.plot([min(xs), max(xs)], [y, y], color=INK2, linewidth=0.6, zorder=1)
    for a, b, alt in exons:
        ax.add_patch(Rectangle((a, y - 0.22), b - a, 0.44, facecolor=ACCENT if alt else "white", edgecolor=INK,
                               linewidth=0.5, zorder=2))
        rows.append({"panel": "a", "view": view, "isoform": k, "exon_start": a, "exon_end": b, "alternative": alt})
    if view in EXTENSION and k == 0:
        a, b = EXTENSION[view]
        ax.add_patch(Rectangle((a, y - 0.22), b - a, 0.44, facecolor=ACCENT, edgecolor=INK, linewidth=0.5, zorder=3))
        rows.append({"panel": "a", "view": view, "isoform": k, "exon_start": a, "exon_end": b, "alternative": True})


def _gex(ax, rows):
    rng = np.random.default_rng(0)
    y = -0.6
    for a in (0, 4, 8):
        ax.add_patch(Rectangle((a, y - 0.22), EXON, 0.44, facecolor=MUTED, edgecolor=INK, linewidth=0.5, zorder=2))
        for x in a + rng.random(6) * EXON:
            ax.plot([x, x], [0.15, 0.75], color=INK2, linewidth=0.5)
    ax.plot([0, 9], [y, y], color=INK2, linewidth=0.6, zorder=1)
    rows.append({"panel": "a", "view": "GEX", "isoform": 0, "exon_start": 0, "exon_end": 9, "alternative": False})


def counts() -> pd.DataFrame:
    sc = pd.read_csv(CFG.results_dir / "data" / "task_scope.csv")
    perf = pd.read_csv(CFG.results_dir / "subtype" / "summary" / "performance.csv")
    sub = perf.drop_duplicates("cohort").groupby("project").cohort.nunique()
    t = pd.read_parquet(CFG.data_dir / "msi" / "models_tumours.parquet")
    mmr = t.groupby("project").mmr.sum()
    d = sc.assign(subtype_cohorts=sc.project.map(sub).fillna(0).astype(int),
                  mmr_deficient=sc.project.map(mmr).fillna(0).astype(int),
                  msi_role=np.where(sc.project.isin(TRAIN), "training", "held out"))
    return d.sort_values("tumour_cases", ascending=False).reset_index(drop=True)


def O1() -> None:
    d = counts()
    fig = plt.figure(figsize=(DOUBLE, 4.3))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.42, 1.5], hspace=0.22)
    top = gs[0].subgridspec(1, 9, wspace=0.18)
    rows = []
    axes_a = []
    for i, view in enumerate([*VIEWS, "GEX"]):
        ax = fig.add_subplot(top[i])
        if view == "GEX":
            _gex(ax, rows)
        else:
            for k, (y, exons) in enumerate(zip([0.6, -0.6], VIEWS[view])):
                _isoform(ax, y, exons, view, rows, k)
        ax.set_ylim(-1.55, 1.0)
        ax.set_xlim(-0.5, 9.5)
        ax.axis("off")
        ax.text(4.5, -1.05, view, ha="center", va="top", fontsize=6.5, color=INK, fontweight="bold")
        axes_a.append(ax)
    bot = gs[1].subgridspec(4, 1, height_ratios=[1.0, 0.13, 0.13, 0.13], hspace=0.12)
    bx = fig.add_subplot(bot[0])
    x = np.arange(len(d))
    w = 0.38
    bx.bar(x - w / 2, d.tumour_cases, width=w, color=TUMOUR, label="tumour")
    bx.bar(x + w / 2, d.normal_cases.replace(0, np.nan), width=w, color=NORMAL, label="normal")
    bx.set_yscale("log")
    bx.set_ylim(1, 2000)
    bx.set_xlim(-0.6, len(d) - 0.4)
    bx.set_xticks([])
    bx.tick_params(axis="x", which="both", length=0)
    bx.set_ylabel("cases", fontsize=6)
    bx.axhline(10, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    bx.legend(loc="upper right", fontsize=5.8, frameon=False, ncol=2, handlelength=1.0)
    lab = ["tumour-status model", "subtype cohorts", "MMR-deficient (MC3)"]
    for r, (key, name) in enumerate(zip(["tumor_status_in_scope", "subtype_cohorts", "mmr_deficient"], lab)):
        ax = fig.add_subplot(bot[r + 1])
        ax.set_xlim(bx.get_xlim())
        for i, row in d.iterrows():
            v = row[key]
            if key == "tumor_status_in_scope":
                if v:
                    ax.scatter(i, 0, s=10, color=INK, zorder=3)
            elif v:
                bold = key == "mmr_deficient" and row.project in TRAIN
                ax.text(i, 0, str(int(v)), ha="center", va="center", fontsize=5.2,
                        color=ACCENT_DARK if bold else INK, fontweight="bold" if bold else "normal")
        ax.set_ylim(-0.6, 0.6)
        ax.set_yticks([0])
        ax.set_yticklabels([name], fontsize=5.8)
        ax.tick_params(axis="both", which="both", length=0)
        ax.minorticks_off()
        for s in ["left", "bottom", "top", "right"]:
            ax.spines[s].set_visible(False)
        if r < 2:
            ax.tick_params(axis="x", length=0, labelbottom=False)
        else:
            ax.set_xticks(x)
            ax.set_xticklabels(d.project, rotation=90, fontsize=5.6)
            ax.tick_params(axis="x", length=0)
    for _, row in d.iterrows():
        rows.append({"panel": "b", "project": row.project, "tumour_cases": row.tumour_cases,
                     "normal_cases": row.normal_cases, "tumour_status_model": bool(row.tumor_status_in_scope),
                     "subtype_cohorts": row.subtype_cohorts, "mmr_deficient_mc3": row.mmr_deficient,
                     "msi_role": row.msi_role})
    rows.append({"panel": "b", "project": "reference line: 10 normal cases (tumour-status scope)"})
    panel_letters([(axes_a[0], -0.25), (bx, -0.08)], y=1.02)
    save(fig, OUT, "O1_study_design", data=pd.DataFrame(rows), caption=caption(d))


def caption(d: pd.DataFrame) -> str:
    n_t, n_n = int(d.tumour_cases.sum()), int(d.normal_cases.sum())
    st = d[d.tumor_status_in_scope]
    return tidy(f"""**Figure O1. Nine views of the same RNA-seq samples, and the analyses each cancer enters.** (a) The
nine views, each modelled separately. Each splicing view compares two isoforms (top: the form whose use is
measured, with the alternative segment in blue; bottom: the other form). AFE and ALE, alternative first and last
exons; HIT, the hybrid internal-terminal index of an exon that can start or end a transcript; SE, a skipped
(cassette) exon; A5SS and A3SS, alternative 5′ and 3′ splice sites; MXE, mutually exclusive exons; RI, a retained
intron; GEX, gene expression from read counts. (b) Cases per TCGA project ({n_t:,} tumour, {n_n:,} normal; log
scale; dashed line = 10). Below: the {len(st)} projects with at least 10 normal cases, which have a tumour-status
model; the number of subtype cohorts drawn from each project; and the tumours with at least 10 MC3 homopolymer
deletions (MMR-deficient), in bold blue for the three cancers that trained the MSI models. Every project enters
the 34-class cancer-type model.""")


if __name__ == "__main__":
    O1()
