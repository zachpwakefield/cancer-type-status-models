"""Step 53: figures drawn only for the MSI manuscript (manuscript/MSI_SPLICING_CODE.md).

    M1  the two-arm model of how mismatch-repair deficiency moves 3' splice-site use (schematic)

The manuscript's other figures are the subtype figures U10-U22 (step 36). Every figure is SVG +
PDF + 400-dpi PNG with `<stem>.csv` of the drawn marks and `<stem>.caption.md` (a figure legend).
Usage: python scripts/53_msi_manuscript_figures.py [M1]        (default: all)
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from _figbase import legend, panel_letters, run, tidy
from matplotlib.patches import Rectangle

from ctsm.figstyle import ACCENT, ACCENT_DARK, ACCENT_LIGHT, DOUBLE, INK, INK2, MUTED, save
from ctsm.subtype import subtype_config

OUT = subtype_config().results_dir / "figures" / "manuscript"
TITLES = {"M1": "A two-arm model of how mismatch-repair deficiency moves 3′ splice-site use."}
# illustrative geometry, nt upstream of the 3' splice site: (run near end, run far end, branch point)
PANELS = [("run 10–19 nt", (2, 15, 24), (2, 12, 21), ACCENT_LIGHT, "↓"),
          ("run ≥ 20 nt, branch point ≥ 36 nt", (2, 26, 44), (2, 21, 39), ACCENT_DARK, "↑"),
          ("run upstream of the branch point", (30, 52, 24), (30, 47, 24), MUTED, "=")]


def _row(ax, y, geo, exon_col, rows, panel, state):
    near, far, bp = geo
    ax.plot([62, 0], [y, y], color=INK2, linewidth=0.8, solid_capstyle="butt", zorder=1)
    ax.add_patch(Rectangle((-14, y - 0.16), 14, 0.32, facecolor=exon_col, edgecolor=INK, linewidth=0.5, zorder=2))
    ax.add_patch(Rectangle((near - 0.5, y - 0.08), far - near + 1, 0.16, facecolor=ACCENT, edgecolor="none", zorder=2))
    ax.text((near + far) / 2, y + 0.14, f"T{far - near + 1}", ha="center", va="bottom", fontsize=5.4, color=ACCENT_DARK)
    ax.scatter([bp], [y], s=14, color=INK, zorder=3, marker="v")
    ax.text(bp, y + 0.14, "A", ha="center", va="bottom", fontsize=5.4, color=INK)
    ax.text(0, y + 0.19, "AG", ha="center", va="bottom", fontsize=5.0, color=INK2)
    yb = y - 0.3
    ax.plot([bp, bp, 0, 0], [yb + 0.06, yb, yb, yb + 0.06], color=MUTED, linewidth=0.5)
    ax.text(bp / 2, yb - 0.04, f"{bp} nt", ha="center", va="top", fontsize=5.0, color=INK2)
    rows += [{"panel": panel, "state": state, "element": "T-run", "from_nt": near, "to_nt": far},
             {"panel": panel, "state": state, "element": "branch point A", "from_nt": bp, "to_nt": bp},
             {"panel": panel, "state": state, "element": "exon", "from_nt": 0, "to_nt": -14}]


def M1() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 1.9))
    fig.subplots_adjust(left=0.05, right=0.99, top=0.82, bottom=0.04, wspace=0.12)
    rows = []
    for k, (ax, (title, mss, msi, col, mark)) in enumerate(zip(axes, PANELS)):
        panel = "abc"[k]
        _row(ax, 1.25, mss, MUTED, rows, panel, "MSS")
        _row(ax, 0.25, msi, col, rows, panel, "MSI")
        ax.text(-16, 0.25, mark, ha="left", va="center", fontsize=9, color=INK)
        rows.append({"panel": panel, "state": "MSI", "element": f"exon inclusion {mark}"})
        ax.set_xlim(63, -19)
        ax.set_ylim(-0.25, 1.7)
        ax.set_title(title, fontsize=6.2, pad=4)
        ax.axis("off")
        if k == 0:
            for y, lab in [(1.25, "MSS"), (0.25, "MSI")]:
                ax.text(66, y, lab, ha="right", va="center", fontsize=6, color=INK)
    panel_letters([(axes[0], -0.08), (axes[1], -0.02), (axes[2], -0.02)], y=1.08)
    save(fig, OUT, "M1_msi_code_model", data=pd.DataFrame(rows), caption=tidy(legend("M1", TITLES["M1"], _cap_m1())))


def _cap_m1() -> str:
    return """**How might mismatch-repair deficiency move a 3′ splice site?**

Schematic with illustrative lengths; each panel shows one acceptor before (MSS) and after (MSI) a contraction of its
T-run. Blue bar = sense-strand T-run; triangle = branch-point adenosine (A); AG = the 3′ splice site; box = the
cassette exon; bracket = branch-point-to-AG distance. The arrow beside the MSI exon is the predicted change in
inclusion.

(a) Runs of 10–19 nt sit at or below the length at which a pyrimidine tract stops gaining strength. Losing a few
nucleotides weakens U2AF2 binding, and the exon is skipped more (the falling arm, R_down). (b) Runs of 20 nt or more
lie between a distant branch point and the AG. After contraction the tract is still long, but the branch point now
lies closer to the AG, nearer the distance at which the second step of splicing works best. The 3′ splice site is
used more and the exon is included more (the rising arm, R_up). (c) A long run upstream of the branch point. Its
contraction leaves both the tract between the branch point and the AG and their distance unchanged, so geometry
predicts no change. A factor that binds the run itself would not make that distinction.

Panels (a) and (b) are the code tested in Figures U17–U21; panel (c) is the control of D-70 (RV5, Figure U22)."""


if __name__ == "__main__":
    run({"M1": M1})
