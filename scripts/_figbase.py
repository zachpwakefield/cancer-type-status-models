"""Shared setup for the two figure scripts (house style, paths, small helpers).

The figures are split by narrative, because the two tasks answer different questions and are
read separately:
    scripts/20_figures_status.py -> results/figures/tumour_status/   (S1-S10)
    scripts/21_figures_type.py   -> results/figures/cancer_type/     (T1-T9)
Both read only the tables written by steps 08-23 (and, for panels that plot individual samples,
the local caches). S10 and T9 are one figure drawn for each task by `overlap_figure`. Every figure is SVG + PDF + 400-dpi PNG with `<stem>.csv` of the plotted
values and `<stem>.caption.md`; per-sample values go to `<data_dir>/figure_data/` (D-12).
"""
from __future__ import annotations

import re
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import log

from ctsm.config import load_config
from ctsm.figstyle import ACCENT_DARK, ACCENT_LIGHT, DOUBLE, GRID, INK2, MUTED, VIEW_NAMES, apply_style, gene_label

apply_style()
CFG = load_config()
BIO = CFG.results_dir / "biology"
LOCAL = CFG.data_dir / "figure_data"
FIGURES = CFG.results_dir / "figures"
AS_VIEWS = [v for v in CFG.views if v != "gex"]
SCOPE = pd.read_csv(CFG.results_dir / "data" / "task_scope.csv")
CANCERS = list(SCOPE.loc[SCOPE.tumor_status_in_scope, "project"])


def locus_label(gene, where: str) -> str:
    """Display name for an event's host gene: its HGNC symbol, or, for the few genes with no
    symbol, the locus (`chr6:41.78 Mb`) read out of `where` (a region or an event id) rather than
    the Ensembl ID the house style keeps out of figures."""
    name = gene_label(gene)
    if name != "unnamed gene":
        return name
    m = re.search(r"(chr[\dXYM]+):(\d+)", str(where))
    return f"{m.group(1)}:{int(m.group(2)) / 1e6:.2f} Mb" if m else name


def short(region: str) -> str:
    ch, rg = region.split(":")
    a, b = rg.split("-")
    return f"{ch}:{int(a):,}-{int(b):,}"


def run(figs: dict) -> None:
    """Draw the named figures (default: all of this narrative)."""
    for name in sys.argv[1:] or list(figs):
        try:
            figs[name]()
            log(f"{name} written")
        except FileNotFoundError as e:
            log(f"{name} skipped (inputs missing: {e.filename})")


def overlap_figure(task: str) -> tuple[plt.Figure, pd.DataFrame]:
    """S10 / T9 (D-27): how often the models' top events overlap another gene's exon, against
    every measured event of the view, and how many of those overlapping events follow the other
    gene's expression. One function, so the two narratives draw exactly the same figure."""
    d = BIO / "12_overlap_audit"
    enr = pd.read_csv(d / "overlap_enrichment_by_view.csv")
    enr = enr[(enr["task"] == task) & (enr["overlap"] == "any")].set_index("view")
    tr = pd.read_csv(d / ("tracking_tumour_status.csv" if task == "tumour status" else "tracking_cancer_type.csv"))
    tr = tr.drop_duplicates(["view", "feature_id"])
    views = list(enr.sort_values("top_share").index)          # bottom-to-top: most overlap on top

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.5), sharey=True,
                             gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.12})
    ax = axes[0]
    rows = []
    for y, v in enumerate(views):
        r = enr.loc[v]
        ax.plot([r.all_share * 100, r.top_share * 100], [y, y], color=GRID, linewidth=1.0, zorder=1)
        ax.scatter(r.all_share * 100, y, s=16, color=MUTED, linewidths=0, zorder=2)
        ax.scatter(r.top_share * 100, y, s=18, color=ACCENT_DARK, linewidths=0, zorder=3)
        rows.append({"panel": "a", "view": v, "all_events": int(r.all_events),
                     "all_overlapping": int(r.all_overlapping), "all_share": r.all_share,
                     "top_events": int(r.top_events), "top_overlapping": int(r.top_overlapping),
                     "top_share": r.top_share, "odds_ratio": r.odds_ratio, "q_fisher": r.q_fisher,
                     "odds_ratio_length_matched": r.odds_ratio_length_matched,
                     "q_length_matched": r.q_length_matched})
    ax.set_yticks(range(len(views)))
    ax.set_yticklabels([VIEW_NAMES[v] for v in views])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.6, len(views) - 0.4)
    ax.set_xlim(0, None)
    ax.set_xlabel("events whose segment overlaps another gene's exon (%)")
    h = [plt.Line2D([], [], marker="o", linestyle="", markersize=3.4, color=MUTED, markeredgewidth=0,
                    label="every measured event"),
         plt.Line2D([], [], marker="o", linestyle="", markersize=3.4, color=ACCENT_DARK, markeredgewidth=0,
                    label="the models' top events")]
    ax.legend(handles=h, frameon=False, fontsize=5.8, loc="lower right", handletextpad=0.3)

    ax = axes[1]
    parts = [("tracks the other gene", ACCENT_DARK), ("does not", ACCENT_LIGHT),
             ("other gene not measured", GRID)]
    for y, v in enumerate(views):
        g = tr[tr["view"] == v]
        measured = g["neighbour_in_gex"].astype(bool)
        n = [int(g["tracks_neighbour"].eq(True).sum()), 0, int((~measured).sum())]
        n[1] = len(g) - n[0] - n[2]
        left = 0
        for k, (lab, colr) in zip(n, parts):
            ax.barh(y, k, left=left, height=0.62, color=colr, linewidth=0)
            left += k
        ax.text(left + 4, y, f"{n[0]} of {len(g)}", fontsize=5.6, color=INK2, va="center")
        rows.append({"panel": "b", "view": v, "top_overlapping_tested": len(g), "tracks": n[0],
                     "does_not_track": n[1], "other_gene_not_measured": n[2]})
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, None)
    ax.set_xlabel("top events that overlap another gene (number)")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c, linewidth=0) for _, c in parts],
              labels=[lab for lab, _ in parts], frameon=False, fontsize=5.8, loc="lower right",
              handlelength=1.0, handletextpad=0.4)
    panel_letters([(axes[0], -0.2), (axes[1], -0.04)])
    out = pd.DataFrame(rows)
    counts = ["all_events", "all_overlapping", "top_events", "top_overlapping", "top_overlapping_tested",
              "tracks", "does_not_track", "other_gene_not_measured"]
    return fig, out.astype({c: "Int64" for c in counts})


def overlap_caption(task: str, data: pd.DataFrame, fig_id: str) -> str:
    """Caption for S10 / T9, every number read from the audit tables."""
    d = BIO / "12_overlap_audit"
    left = data[data["panel"] == "a"].set_index("view")
    right = data[data["panel"] == "b"].set_index("view")
    ann = pd.read_csv(d / "event_overlaps_all_measured.csv.gz", usecols=["overlap_class"])
    base = ann["overlap_class"].isin(["antisense", "same-strand"]).mean()
    named = pd.read_csv(d / "named_events_overlaps.csv")
    named = named[named["task"] == task]
    hit = named[named["tracks_neighbour"].eq(True)]
    sgn = lambda x: f"{x:+.2f}".replace("-", "\u2212")
    ex = "; ".join(f"{r.symbol} {r.view.upper()} follows {r.other_gene} (\u03c1\u00a0{sgn(r.rho_neighbour)}, "
                   f"own gene {sgn(r.rho_host)})" for r in hit.itertuples())
    top = left.sort_values("top_share").index[-1]
    unit = ("top 100 of each cancer's model, per view" if task == "tumour status"
            else "top 30 of each class's explanation, per view")
    scope = ("within each cancer, median over the cancers whose top 100 holds the event" if task == "tumour status"
             else "across all 9,729 samples")
    text = f"""
**Do the models' events measure another gene?**

TCGA mRNA libraries are mostly unstranded, and GENCODE readthrough and nested models put exons of
one gene inside another, so an event named after one gene can count reads of a neighbour (D-27).
Left: the share of events whose segment (the alternative exon, both exons of an MXE pair, the
retained intron or the splice-site extension) overlaps an exon of another GENCODE v45 gene by
>= 50 bp, on either strand, for every measured event (grey; {base:.1%} over all 8 views) and for
the models' top events ({unit}; dark). The gap is widest for {VIEW_NAMES[top]}:
{left.loc[top, 'top_share']:.0%} against {left.loc[top, 'all_share']:.0%}. Odds ratios with the
queried length held fixed (Mantel-Haenszel over length deciles, because long terminal exons overlap
something more often by chance) are {left['odds_ratio_length_matched'].min():.1f} to
{left['odds_ratio_length_matched'].max():.1f} across views; every statistic is in `{fig_id}_overlap_audit.csv`.

Right: of the top events with such an overlap, how many *track* the other gene — Spearman with its
log2 CPM of at least 0.5 in size and larger than with the host gene's own ({scope}); the criterion
was fixed before any result was seen. {int(right['tracks'].sum()):,} of {int(right['top_overlapping_tested'].sum()):,}
do. Those are better read as a readout of the neighbour's expression than as a splicing choice of
the host, whether the coupling is misassigned reads or genuine co-transcription; the two cannot be
told apart without stranded libraries. Among the events the narrative names: {ex if ex else 'none'}.

No model was refitted: the audit changes how events are described, not what the models are.
"""
    paras = [" ".join(p.split()) for p in text.strip().split("\n\n")]
    return "\n\n".join(textwrap.fill(p, width=112, break_on_hyphens=False) for p in paras) + "\n"


def panel_letters(pairs, y: float = 1.04) -> None:
    """Bold panel letters a, b, ... above each axes, at axes-fraction x: pairs = [(ax, x), ...]."""
    for (ax, x), lab in zip(pairs, "abcdefgh"):
        ax.text(x, y, lab, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def legend(fig_id: str, title: str, caption: str) -> str:
    """Question-form caption -> figure legend. The bold question line becomes 'Figure <id>. <title>';
    the panel markers 'Left:' and 'Right:', or 'a, ' at a sentence start, become (a), (b), ... with a
    capital after the letter."""
    body = caption.strip().split("\n", 1)[1].lstrip("\n")
    cap = lambda m: f"{m.group(1)}({m.group(2)}) {m.group(3).upper()}"
    body = re.sub(r"(?m)(^|\. )Left: (\w)", lambda m: f"{m.group(1)}(a) {m.group(2).upper()}", body)
    body = re.sub(r"(?m)(^|\. )Right: (\w)", lambda m: f"{m.group(1)}(b) {m.group(2).upper()}", body)
    body = re.sub(r"(?m)(^|\. )([a-f]), ([a-f]), (\w)",
                  lambda m: f"{m.group(1)}({m.group(2)}, {m.group(3)}) {m.group(4).upper()}", body)
    body = re.sub(r"(?m)(^|\. )([a-f]), (\w)", cap, body)
    body = re.sub(r"(?m)(^|\. )\(([a-f])\) (\w)", cap, body)
    return f"**Figure {fig_id}. {title}** {body}"


def inject_legends(md: Path) -> int:
    """Copy each figure's legend (its caption without the bullet lists of tests) under every
    `![Xn](figures/<dir>/<stem>.png)` in a write-up, between marker comments, replacing any block
    already there, so the write-up and the caption files never drift. Returns the figures done."""
    lines = md.read_text().split("\n")
    out, i, n = [], 0, 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        m = re.match(r"^!\[([STUX]\d+|RP\d+|RPM)\]\((figures/[a-z_]+)/(.+)\.png\)$", line)
        i += 1
        if not m:
            continue
        j = i
        while j < len(lines) and lines[j].strip() == "":
            j += 1
        if j < len(lines) and lines[j].startswith("<!-- legend"):
            while not lines[j].startswith("<!-- /legend"):
                j += 1
            i = j + 1
        rel = f"{m.group(2)}/{m.group(3)}.caption.md"
        paras = [p for p in (md.parent / rel).read_text().strip().split("\n\n") if not p.startswith("- ")]
        block = [f"<!-- legend {m.group(1)}: copied from {rel} by the figure script's 'legends' step -->"]
        for k, p in enumerate(paras):
            if k:
                block.append(">")
            # a list inside a legend paragraph (a key to the columns, say) stays a list
            segs = []
            for ln in p.split("\n"):
                if ln.startswith("- ") or not segs or (segs[-1][0] == "item") != ln.startswith("  "):
                    segs.append(["item" if ln.startswith("- ") else "text", ln])
                else:
                    segs[-1][1] += " " + ln
            for kind, seg in segs:
                text = " ".join(seg.split())
                block += ["> " + ln for ln in textwrap.wrap(text, 112 if kind == "item" else 114, break_on_hyphens=False,
                                                           subsequent_indent="  " if kind == "item" else "")]
        block.append("<!-- /legend -->")
        out += [""] + block
        n += 1
    md.write_text("\n".join(out))
    log(f"legends: {n} figures in {md.name}")
    return n


def tidy(caption: str, width: int = 116) -> str:
    """Re-wrap every paragraph and bullet of a caption to the house limit (lines <= 120 characters).
    A bullet list inside a paragraph stays a list; a bullet's continuation lines are indented."""
    out = []
    for para in caption.strip().split("\n\n"):
        lines = para.split("\n")
        if para.startswith("**") and len(para) <= 120:
            out.append(para)
            continue
        blocks, cur, kind = [], [], None
        for ln in lines:
            if ln.startswith("- "):
                if cur:
                    blocks.append((kind, cur))
                cur, kind = [ln.strip()], "item"
            elif kind == "item" and (ln.startswith("  ") or lines[0].startswith("- ")):
                cur.append(ln.strip())
            else:
                if kind == "item":
                    blocks.append((kind, cur))
                    cur = []
                cur.append(ln.strip())
                kind = "text"
        blocks.append((kind, cur))
        out.append("\n".join(
            textwrap.fill(" ".join(b), width, subsequent_indent="  " if k == "item" else "", break_on_hyphens=False)
            for k, b in blocks))
    return "\n\n".join(out)
