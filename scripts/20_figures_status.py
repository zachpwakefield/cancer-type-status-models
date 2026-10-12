"""Step 20: the tumour-status narrative, in figures (house style, docs/FIGURE_GUIDELINES.md).

One question per figure, in the order the story is told in results/TUMOUR_STATUS.md:
    S1  does combining the modalities beat expression alone?
    S2  does splicing shift in one direction between a tumour and its matched normal tissue?
    S3  which processes do the models' genes converge on?
    S4  which splicing factors have coherent programmes, and are the models' events inside them?
    S5  do different cancers rely on the same events?
    S6  what do the exons the models use do to the protein?
    S7  what do the most shared events look like, and how does the model read them?
    S8  are the models' events associated with how long patients live?
    S9  where in the protein does each of the most shared events sit?
    S10 do the models' events measure another gene's expression?

Usage: python scripts/20_figures_status.py [S1 S5 ...]      (default: all)
"""
from __future__ import annotations

import re

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from _figbase import (AS_VIEWS, BIO, CANCERS, CFG, FIGURES, LOCAL, locus_label, log, overlap_caption,
                      overlap_figure, run, short)
from matplotlib.lines import Line2D
from scipy.stats import beta, binomtest

from ctsm.biology import ProteinMap, bio_config, feature_values, gex_by_symbol, task_table
from ctsm.config import PAN, load_config
from ctsm.figstyle import (ACCENT, ACCENT_DARK, DIVERGE, DOUBLE, GRID, INK, INK2, MUTED, NORMAL, SEQ, SINGLE,
                           TN_DIVERGE, TUMOUR, VIEW_NAMES, bracket, colourbar, fp, gene_label, heat_axes, save,
                           stars)

OUT = FIGURES / "tumour_status"
SHIFT_ORDER = ["ri", "a5ss", "a3ss", "se", "mxe", "ale", "afe"]
SHIFT_LAB = ["RI retained", "A5SS long form", "A3SS long form", "SE included", "MXE upstream exon",
             "ALE downstream", "AFE downstream"]


def wilson(k: int, n: int) -> tuple[float, float]:
    """95% Clopper-Pearson interval for a share (used for hit rates)."""
    lo = beta.ppf(0.025, k, n - k + 1) if k else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


# ============================================================ S1 combining modalities
def S1() -> None:
    d = BIO / "09_ensemble"
    m = pd.read_csv(d / "tumor_status_model_metrics.csv")
    t = pd.read_csv(d / "tumor_status_model_tests.csv")
    show = ["GEX", "AS ensemble", "GEX + AS"]
    col = {"GEX": INK, "AS ensemble": ACCENT, "GEX + AS": ACCENT_DARK}
    piv = m.pivot(index="cancer", columns="model", values="n_errors")
    n = m.drop_duplicates("cancer").set_index("cancer")["n"]
    rate = piv.div(n, axis=0) * 100
    order = rate["GEX"].sort_values().index
    fused_p = t[t.model.eq("GEX + AS")].set_index("cancer")

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 3.0), gridspec_kw={"width_ratios": [2.1, 1]})
    ax = axes[0]
    y = np.arange(len(order))
    for i, c in enumerate(order):
        ax.plot([rate.loc[c, show].min(), rate.loc[c, show].max()], [i, i], color=GRID, linewidth=1.4, zorder=1)
    for z, mdl in enumerate(["AS ensemble", "GEX + AS", "GEX"]):      # GEX on top; rings show overlaps
        ax.scatter(rate.loc[order, mdl], y, s=11, color=col[mdl], linewidths=0.35, edgecolors="white",
                   zorder=3 + z, label=mdl)
    for i, c in enumerate(order):
        s = stars(fused_p.loc[c, "p_mcnemar"]) if c in fused_p.index else ""
        if s:
            ax.text(rate.loc[c].max() + 0.25, i, s, fontsize=6, color=INK2, va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} ({int(n[c]):,})" for c in order], fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.8, len(order) - 0.2)
    ax.set_xlabel("out-of-fold errors (% of samples)")
    h, lab = ax.get_legend_handles_labels()
    ax.legend([h[lab.index(k)] for k in show], show, loc="lower right", fontsize=6,
              handletextpad=0.3, borderpad=0.2)

    ax = axes[1]
    pooled = m.groupby("model")["n_errors"].sum()
    rows = ["GEX", "AS ensemble", "GEX + AS", "stacked"]
    yy = np.arange(len(rows))[::-1]
    for v, i, name in zip(pooled[rows], yy, rows):
        ax.plot([0, v], [i, i], color=GRID, linewidth=0.8, zorder=1)
        ax.scatter(v, i, s=14, color=col.get(name, MUTED), linewidths=0, zorder=3,
                   marker="D" if name == "stacked" else "o")
        ax.text(v + 3, i, f"{int(v)}", fontsize=6, color=INK2, va="center")
    ax.set_yticks(yy)
    ax.set_yticklabels(rows, fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlim(0, pooled[rows].max() * 1.22)
    ax.set_xlabel(f"errors, all {len(order)} cancers ({int(n.sum()):,} samples)")
    pt = t.groupby("model")[["only_a_correct", "only_b_correct"]].sum()
    gain = int(pt.loc["GEX + AS", "only_b_correct"]), int(pt.loc["GEX + AS", "only_a_correct"])
    pooled_p = binomtest(gain[1], sum(gain), 0.5).pvalue

    data = m.merge(t[["cancer", "model", "only_a_correct", "only_b_correct", "p_mcnemar", "q_mcnemar"]],
                   on=["cancer", "model"], how="left")
    save(fig, OUT, "S1_modality_gain", data=data, caption=f"""
**Does combining the modalities separate tumour from normal tissue better than expression alone?**

Left: out-of-fold errors as a share of each cancer's samples (n in brackets) for the expression model (black), the
AS ensemble (the mean of the 8 AS views' probabilities, blue) and an equal-weight fusion of the two (dark blue).
The fusion weights were fixed before any result was seen (D-14). Stars: exact McNemar, fusion against GEX, on
paired per-sample correctness (* p < 0.05, ** p < 0.01, *** p < 0.001); q (BH over the 45 tests) is in the CSV.

Right: errors pooled over all cancers. 'stacked' is a cross-fitted logistic regression on the 9 views'
probabilities (D-23): it is shown as a ceiling, not as the headline, because the base models behind its training
probabilities did see the samples it scores.

Pooled over the cohort, fusion corrects {gain[0]} of the expression model's {int(pooled['GEX'])} errors and
introduces {gain[1]} (exact McNemar p = {fp(pooled_p)}). No single cancer has enough errors to reach significance
on its own; the per-cancer tests are in the CSV.
""")


# ============================================================ S4 regulons
def S4() -> None:
    d = BIO / "10_regulons"
    fe = pd.read_csv(d / "regulon_enrichment_by_factor.csv")
    ov = pd.read_csv(d / "regulon_overlap.csv")
    sig = fe[fe["q_fisher"] < 0.05].sort_values("odds_ratio_haldane")

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.9), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    y = np.arange(len(sig))
    for i, r in enumerate(sig.itertuples()):
        c = DIVERGE(0.95) if r.odds_ratio_haldane > 1 else DIVERGE(0.05)   # house convention: red above, blue below
        ax.plot([1, r.odds_ratio_haldane], [i, i], color=GRID, linewidth=0.8, zorder=1)
        ax.scatter(r.odds_ratio_haldane, i, s=13, color=c, linewidths=0, zorder=3)
    ax.axvline(1, color=INK2, linewidth=0.5)
    ax.set_xscale("log")
    plain_log(ax.xaxis, [0.3, 1, 3, 10, 30])
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.factor} ({r.n_core_events})" for r in sig.itertuples()], fontsize=5.8)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.7, len(sig) - 0.3)
    ax.set_xlabel("odds that a model event is in this regulon,\nagainst matched events")

    ax = axes[1]
    m = ov.merge(fe[["factor", "odds_ratio_haldane", "q_fisher"]], on="factor", how="left")
    lo = np.log2(m["odds_ratio_haldane"].fillna(1.0)).clip(-2, 2)
    sc = ax.scatter(m["n_core_events"], m["max_jaccard"], c=lo, cmap=DIVERGE, vmin=-2, vmax=2,
                    s=16, linewidths=0.3, edgecolors="white", zorder=3)
    ax.set_xscale("log")
    plain_log(ax.xaxis, [20, 50, 100, 200, 500])
    ax.set_xlabel("events in the core regulon")
    ax.set_ylabel("largest overlap with another\nfactor's regulon (Jaccard)")
    for name, off in [("PTBP1", (5, 2)), ("RBM47", (5, -1)), ("QKI", (5, 2)), ("PCBP1", (5, 2)),
                      ("PABPN1", (-4, -9)), ("RBM39", (-20, 2))]:
        r = m[m["factor"].eq(name)]
        if len(r):
            ax.annotate(name, (float(r["n_core_events"].iloc[0]), float(r["max_jaccard"].iloc[0])),
                        textcoords="offset points", xytext=off, fontsize=5.5, color=INK)
    colourbar(fig, sc, ax, "log2 odds (model events)", ticks=[-2, 0, 2], shrink=0.7, extend="both")
    up = sig[sig.odds_ratio_haldane > 1]["factor"].tolist()
    save(fig, OUT, "S4_regulons", data=fe.merge(ov.drop(columns=["n_core_events"]), on="factor", how="outer"),
         caption=f"""
**Which splicing factors' programmes do the models' events belong to?**

A factor's *core regulon* is the events whose PSI tracks its expression within tumours — partial Spearman
|rho| >= 0.3 at p < 0.05 after rank-residualising both on a stromal score, an immune score and the sample's own mean
PSI in that view — reproducibly, with the same sign, in at least 5 of the 15 cancers (D-24).

Left: for every factor with at least 15 marked events, the odds that one of the models' top events is in its core
regulon, against the coverage- and spread-matched background of D-22. Only factors with q < 0.05 (BH over the
{len(fe)} factors tested) are shown; the regulon's size is in brackets. Odds ratios use a 0.5 correction so that an
empty cell stays finite — PTBP1's 15 candidate events against 0 background events is otherwise infinite.

Right: every factor with at least 20 core events. Large regulons are not separate programmes: they overlap each
other heavily (Jaccard up to 0.7), and the 1,074 distinct events carry 7,529 memberships, about 7 factors each.

The two panels say the same thing. The models' events avoid the big shared programme of nuclear RNA-processing
factors ({', '.join(sig[sig.odds_ratio_haldane < 1]['factor'].head(4))}, odds 0.26-0.50) and concentrate in the
small, distinct regulons of tissue and differentiation regulators ({', '.join(up)}).

These are correlations inside tumours: a factor may be co-regulated with its apparent targets rather than drive
them. `regulon_core_events.csv` lists the members, `global_axis_factors.csv` the factors most tied to the
per-sample PSI level that was residualised out.
""")


# ============================================================ S8 survival
def plain_log(axis, ticks: list) -> None:
    """Named ticks on a log axis, written as numbers: 0.5, 1, 10 rather than 5 x 10^-1, 10^0, 10^1."""
    axis.set_major_locator(ticker.FixedLocator(ticks))
    axis.set_major_formatter(ticker.FixedFormatter([f"{t:g}" for t in ticks]))
    axis.set_minor_locator(ticker.NullLocator())


def S8() -> None:
    d = BIO / "08_survival"
    enr = pd.read_csv(d / "survival_enrichment.csv")
    ev = pd.read_csv(d / "event_survival_summary.csv")
    null = pd.read_csv(d / "survival_null.csv") if (d / "survival_null.csv").exists() else None
    pick = lambda arm: (enr[(enr.endpoint == "OS") & (enr.test == "Cox") & (enr.arm == arm)]
                        .set_index("candidate_set"))
    arms = [("base", "PSI + host gene", ACCENT), ("adjusted", "+ age, sex, stage", ACCENT_DARK)]
    sets = [("all candidates", "every\ncandidate"),
            ("recurrent (D-17)", "recurrent\n4+ models"),
            ("top 10 in >= 2 cancers", "top 10 in\n2+ cancers"),
            ("background", "matched\nbackground")]

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 3.0),
                             gridspec_kw={"width_ratios": [1.0, 1.35], "wspace": 0.62})
    ax = axes[0]
    rows = []
    for (arm, lab, colr), off in zip(arms, (-0.13, 0.13)):
        m = pick(arm)
        for i, (key, _) in enumerate(sets):
            if key == "background":
                k, n = int(m["hits_background"].iloc[0]), int(m["n_background"].iloc[0])
            else:
                k, n = int(m.loc[key, "hits_candidate"]), int(m.loc[key, "n_candidate"])
            lo, hi = wilson(k, n)
            ax.vlines(i + off, lo * 100, hi * 100, color=colr, linewidth=0.9)
            ax.scatter(i + off, k / n * 100, s=15, color=colr, linewidths=0, zorder=3,
                       label=lab if i == 0 else None)
            rows.append({"arm": arm, "set": key, "cells": n, "hits": k, "rate": k / n,
                         "ci_low": lo, "ci_high": hi,
                         "odds_ratio": np.nan if key == "background" else m.loc[key, "odds_ratio"],
                         "p_fisher": np.nan if key == "background" else m.loc[key, "p_fisher"]})
    if null is not None:
        r = float(null["rate"].iloc[0]) * 100
        ax.axhline(r, color=INK2, linewidth=0.5, linestyle=(0, (3, 2)))
        ax.text(-0.45, r + 0.2, "permuted outcomes", fontsize=5.5, color=INK2, ha="left", va="bottom")
    ax.set_xticks(range(len(sets)))
    ax.set_xticklabels([s[1] for s in sets], fontsize=5.8)
    ax.set_xlim(-0.55, 3.55)
    ax.set_ylim(0, None)
    ax.set_ylabel("event × cancer pairs with Cox p < 0.05 (%)")
    ax.legend(frameon=False, fontsize=5.8, loc="upper right", handletextpad=0.3,
              borderpad=0.1, labelspacing=0.25)

    ax = axes[1]
    sub = ev[(ev.endpoint == "OS") & (ev.arm == "adjusted") & (ev.role == "candidate")]
    sub = sub.nlargest(12, ["hits", "cohorts_tested"]).sort_values("hr_per_sd")
    for i, r in enumerate(sub.itertuples()):
        c = DIVERGE(0.95) if r.hr_per_sd > 1 else DIVERGE(0.05)
        ax.plot([r.ci_low_sd, r.ci_high_sd], [i, i], color=c, linewidth=0.9)
        ax.scatter(r.hr_per_sd, i, s=14, color=c, linewidths=0, zorder=3)
    ax.axvline(1, color=INK2, linewidth=0.5)
    ax.set_yticks(np.arange(len(sub)))
    ax.set_yticklabels([f"{locus_label(r.gene, r.event_id)} {r.view.upper()} · {r.best_cohort} "
                        f"({r.hits}/{r.cohorts_tested})" for r in sub.itertuples()], fontsize=5.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xscale("log")
    plain_log(ax.xaxis, [0.5, 0.7, 1.0, 1.5])
    ax.set_ylim(-0.7, len(sub) - 0.3)
    ax.set_xlabel("hazard ratio per SD (95% CI), overall survival")
    nulltxt = ("" if null is None else
               f"\n\nThe dashed line is the same pipeline run over outcomes shuffled between patients inside\n"
               f"each cancer: {float(null['rate'].iloc[0]):.1%} ({int(null['hits_p05'].iloc[0]):,} of "
               f"{int(null['cells'].iloc[0]):,} pairs), so both sets sit above what the gates, the host-gene\n"
               f"covariate and the test battery produce on their own.")
    if null is not None:
        rows.append({"arm": "permuted outcomes", "set": "all candidates", "panel": "left",
                     "cells": int(null["cells"].iloc[0]), "hits": int(null["hits_p05"].iloc[0]),
                     "rate": float(null["rate"].iloc[0])})
    plotted = pd.concat([
        pd.DataFrame(rows).assign(panel=lambda z: z["panel"].fillna("left") if "panel" in z else "left"),
        sub.assign(panel="right")[["panel", "gene", "view", "event_id", "best_cohort", "hits", "cohorts_tested",
                                   "hr_per_sd", "ci_low_sd", "ci_high_sd", "cox_p", "cox_q"]],
    ], ignore_index=True)
    save(fig, OUT, "S8_survival", data=plotted, caption=f"""
**Are the events the models rely on associated with how long patients live?**

Every event was tested with splice-assay 0.4.0 in each of the 15 cancers (6,152 tumours, 1,620 deaths): Cox on the
event's value with the host gene's own expression in the model, in the pre-specified base form and again with age,
sex and stage added (D-25). The models' top events are compared with the coverage- and spread-matched background of
D-22, so this is not a comparison of which events are easier to test.

Left: the share of testable event × cancer pairs with a nominal Cox p < 0.05, with 95% Clopper-Pearson intervals;
the three candidate sets are nested and ordered by how many models agree on the event. Odds ratios against the
background of the same arm, base / adjusted: every candidate 1.27 / 1.18, recurrent 1.56 / 1.28, top 10 for two or
more cancers 1.62 / 1.21 (Fisher exact, every set q < 0.05). The rate rises along that order.{nulltxt}

Candidates inside one of the six regulons the models are enriched for (S4) reach p < 0.05 in 13.3% of their pairs
against 10.3% for the rest (OR 1.33, p = 0.005), so the three ways of ranking a candidate -- how many models agree,
whether it recurs, whether it sits in a coherent factor programme -- pick out the same events.

Right: the 12 candidate events associated in the most cancers, each in its strongest cancer: hazard ratio per SD
of PSI (per SD of the HIT index for HIT events) with its 95% CI, and (cancers with p < 0.05 / cancers tested).
Red: a higher value goes with worse survival.

Nominal p values from observational cohorts. BH inside each cancer's own family of ~1,600 tests (`cox_q_cohort`)
leaves 1.2% of candidate cells significant in the adjusted arm, no more than the background's 1.2%, and 94% of them
are KIRC or KIRP: the set-level excess is a broad shift of small associations, not an excess of strong ones.
`event_survival_cells.csv.gz` holds every statistic with both BH families (keyed by `event_index.csv`), and
`event_dossiers.csv` joins each event's recurrence, protein consequence, regulon membership and survival. PSI can
reflect tumour purity or composition, so these are associations, not causes, and not therapeutic claims.
""")


def S2() -> None:
    d = BIO / "02_global_shifts"
    s = pd.read_csv(d / "shift_summary.csv")
    piv = s.pivot(index="cancer", columns="view", values="median_paired_shift")[SHIFT_ORDER] * 100
    pv = s.pivot(index="cancer", columns="view", values="p_wilcoxon")[SHIFT_ORDER]
    order = piv["ri"].sort_values(ascending=False).index
    piv, pv = piv.loc[order], pv.loc[order]
    lim = float(np.nanpercentile(np.abs(piv.to_numpy()), 98))
    fig, ax = plt.subplots(figsize=(3.2, 3.4))
    im = ax.imshow(piv.to_numpy(), cmap=TN_DIVERGE, vmin=-lim, vmax=lim, aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.iat[i, j]
            txt = f"{v:+.1f}".replace("-", "\u2212") + stars(pv.iat[i, j])   # typographic minus
            ax.text(j, i, txt, ha="center", va="center", fontsize=4.8,
                    color="white" if abs(v) > 0.65 * lim else INK)
    npairs = s.drop_duplicates("cancer").set_index("cancer").loc[order, "n_pairs"]
    heat_axes(ax, [f"{c} ({n})" for c, n in npairs.items()], SHIFT_LAB, col_size=5.5, row_size=6)
    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=45, ha="left", rotation_mode="anchor")
    colourbar(fig, im, ax, "median paired tumour − normal (×100)")
    data = s[["cancer", "view", "n_pairs", "n_events", "median_paired_shift", "frac_patients_up", "frac_up",
              "frac_down", "p_wilcoxon", "q_wilcoxon"]].assign(
        q_family="BH across 15 cancers x 7 views (105 Wilcoxon tests)")
    save(fig, OUT, "S2_global_shifts", data=data, caption="""
**Do tumours shift RNA processing in a common direction relative to their own normal tissue?**

Each cell is the median over matched patients (n in brackets) of the patient's mean paired tumour − normal change,
averaged over every event measured in that cancer, so there is no SHAP selection and no circularity. The columns
say what an increase means:
- RI: more intron retention.
- A5SS/A3SS: more use of the long exon form, which includes extra intronic sequence.
- SE: more cassette inclusion.
- MXE: more of the upstream exon.
- ALE/AFE: the used last or first exons move downstream (a 0-1 position index along transcription).

Blue is higher in the tumour and orange higher in the normal tissue, the house tumour/normal colours. Stars:
two-sided Wilcoxon signed-rank test across patients, raw p (* < 0.05, ** < 0.01, *** < 0.001). BH q (family: 15
cancers × 7 views) is in the CSV.

Retention and long splice-site forms rise in most cancers. Last exons move upstream and never downstream. BRCA is
the exception for RI. Caveat: RNA integrity and nuclear-RNA content also raise intronic reads, and RIN is not
available (D-15).
""")


HALLMARK_CAPS = {"KRAS", "TNFA", "NFKB", "E2F", "UV", "DNA", "MYC", "MTORC1", "IL2", "STAT5", "IL6", "JAK",
                 "PI3K", "AKT", "MTOR", "G2M", "TGF", "WNT", "P53", "KEG", "ROS"}


def hallmark_label(term: str) -> str:
    """'HALLMARK_KRAS_SIGNALING_DN' -> 'KRAS signaling down': lower case, but gene and pathway
    symbols stay symbols (house rule: HGNC symbols, never lower-cased)."""
    out = []
    for w in term.replace("HALLMARK_", "").split("_"):
        out.append(w if w in HALLMARK_CAPS else "down" if w == "DN" else w.lower())
    return " ".join(out)


def S3() -> None:
    d = BIO / "03_pathways"
    hm = pd.read_csv(d / "tumor_status_hallmark_convergence.csv", header=[0, 1], index_col=0)
    asn, gx = hm["n_AS_views_significant"][CANCERS], hm["GEX_significant"][CANCERS]
    terms = list(hm.index)
    fig, ax = plt.subplots(figsize=(3.6, 0.15 * len(terms) + 0.8))
    rows = []
    for i, t in enumerate(terms):
        for j, c in enumerate(CANCERS):
            n, g = int(asn.loc[t, c]), int(gx.loc[t, c])
            if n:
                ax.scatter(j, i, s=6 + 9 * n, color=ACCENT, linewidths=0, zorder=2)
            if g:
                ax.scatter(j, i, s=40, facecolor="none", edgecolor=INK, linewidth=0.6, zorder=3)
            rows.append({"term": t, "cancer": c, "n_AS_views_q_lt_0.05": n, "GEX_q_lt_0.05": g})
    heat_axes(ax, [hallmark_label(t) for t in terms], CANCERS, xrot=90, row_size=5.5, col_size=5.5)
    ax.grid(which="minor", visible=False)
    ax.set_xlim(-0.6, len(CANCERS) - 0.4)
    ax.set_ylim(len(terms) - 0.5, -0.5)
    h = [Line2D([], [], marker="o", linestyle="", markersize=np.sqrt(6 + 9 * n), color=ACCENT, markeredgewidth=0,
                label=f"{n} AS view{'s' if n > 1 else ''}") for n in (1, 3)]
    h.append(Line2D([], [], marker="o", linestyle="", markersize=np.sqrt(40), markerfacecolor="none",
                    markeredgecolor=INK, markeredgewidth=0.6, label="GEX"))
    ax.legend(handles=h, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, columnspacing=0.8, handletextpad=0.1)
    save(fig, OUT, "S3_pathways", data=pd.DataFrame(rows).assign(
        q_family="BH within Hallmark per model (one model = one cancer x view)"), caption="""
**Which programmes do the tumour-vs-normal models draw on?**

For every cancer (columns) and view, the 200 genes with the largest SHAP importance were tested for enrichment of
the MSigDB Hallmark sets against all genes measurable in that view and cohort. The test is GOseq/Wallenius,
correcting for the number of events per gene (D-16).
- A dot marks a significant set in at least one AS view; its size is the number of AS views (of 8).
- A ring marks significance in the GEX model.
- Significance is q < 0.05, BH within the Hallmark collection of each model.

GEX models consistently use stromal and microenvironment programmes (EMT, coagulation, KRAS, myogenesis). The AS
models reach the same programmes far less often, mostly through EMT and myogenesis genes. Their top genes form
any enriched set in at most 9 of 15 cancers (SE) and in none for RI, so the tumour splicing signal is largely
not organised by pathway membership.
""")


def S5() -> None:
    d = BIO / "04_events"
    rec = pd.read_csv(d / "recurrent_events.csv")
    P = pd.read_csv(d / "paired_change_top_events_all_cancers.csv.gz")
    T = pd.read_csv(d / "top_events_per_cancer.csv", usecols=["view", "feature_id", "cancer"])
    chance = pd.read_csv(d / "recurrence_vs_chance.csv")
    show = rec.head(40)
    key = list(zip(show.view, show.feature_id))
    piv = P.pivot_table(index=["view", "feature_id"], columns="cancer", values="median_paired_delta")
    M = piv.loc[key].reindex(columns=CANCERS) * 100
    top = T.assign(v=1).pivot_table(index=["view", "feature_id"], columns="cancer", values="v", aggfunc="max")
    top = top.reindex(index=M.index, columns=CANCERS).fillna(0)
    lim = 30
    fig, ax = plt.subplots(figsize=(3.6, 0.105 * len(show) + 0.9))
    im = ax.imshow(M.to_numpy(), cmap=TN_DIVERGE, vmin=-lim, vmax=lim, aspect="auto")
    yy, xx = np.nonzero(top.to_numpy())
    ax.scatter(xx, yy, s=1.5, color=INK, linewidths=0)
    # the region is in the label already, so an unnamed gene just drops the name rather than
    # repeating the locus twice
    labels = [f"{gene_label(s)} {VIEW_NAMES[v]} {short(r)}".replace("unnamed gene ", "")
              for s, v, r in zip(show.symbol, show.view, show.region)]
    heat_axes(ax, labels, CANCERS, xrot=90, row_size=4.6, col_size=5.5)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", labelrotation=90)
    colourbar(fig, im, ax, "median paired tumour − normal (×100)", shrink=0.3, extend="both")
    cons = rec.n_cancers_same_sign_as_majority / rec.n_cancers_measured
    up = rec.assign(up=rec.majority_sign.eq("up in tumour")).groupby("view").up.mean()
    long = M.reset_index().melt(id_vars=["view", "feature_id"], var_name="cancer", value_name="median_paired_delta_x100")
    long = long.merge(top.reset_index().melt(id_vars=["view", "feature_id"], var_name="cancer",
                                              value_name="in_top_100"), on=["view", "feature_id", "cancer"])
    long = long.merge(show[["view", "feature_id", "symbol", "region", "n_cancers"]], on=["view", "feature_id"])
    n_obs, n_exp = int(chance.n_recurrent_observed.sum()), float(chance.n_recurrent_expected_by_chance.sum())
    save(fig, OUT, "S5_recurrent_events", data=long, caption=f"""
**Do different cancers rely on the same events to tell tumour from normal?**

Across the 8 AS views, {n_obs} events rank in the top 100 (mean |SHAP|) of at least 4 of the 15 tumour-vs-normal
models. {n_exp:.1f} would be expected if the rankings were independent (Poisson-binomial over each event's measured
cohorts; `recurrence_vs_chance.csv`).

The figure shows the 40 most recurrent events. Rows give the host gene (HGNC), view and coordinates (hg38, 1-based).
Colour is the median paired tumour − normal change in every cancer, including cancers where the event did not rank
(so the pattern across cancers is not selected on): blue where the event is higher in the tumour, orange where it
is higher in the normal tissue, clipped at ±{lim}. Dots mark the cancers whose top 100 contain the event.

Recurrent events tend to move the same way across cancers. For the median event, {cons.median():.0%} of the cancers
where it is measured share its majority direction; {(cons >= 0.8).mean():.0%} of events reach 80%.

Several genes appear as switch pairs, with two exons of the same gene moving in opposite directions (EXOSC7 last
exons, NDE1 first exons, TPM1 mutually exclusive exons).

Recurrent RI and A5SS events mostly rise in tumours ({up['ri']:.0%} and {up['a5ss']:.0%}), matching the global shift in S2. The full
list with direction consistency is in `recurrent_events.csv`.
""")


def S6() -> None:
    d = BIO / "04_events"
    pr = pd.read_csv(d / "cassette_exon_protein_features.csv")
    sets = ["all measured exons (random 6,000)", "tumour vs normal, top events", "cancer type, top events"]
    short_set = {sets[0]: "all measured exons", sets[1]: "tumour vs normal, top", sets[2]: "cancer type, top"}
    metrics = list(dict.fromkeys(pr["metric"]))
    cols = {sets[0]: MUTED, sets[1]: ACCENT, sets[2]: ACCENT_DARK}
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    rows = []
    for j, st in enumerate(sets):
        sub = pr[pr.set == st].set_index("metric").loc[metrics]
        y = np.arange(len(metrics)) + (j - 1) * 0.24
        ax.hlines(y, sub.ci_lo, sub.ci_hi, color=cols[st], linewidth=0.8)
        ax.scatter(sub.proportion, y, s=9, color=cols[st], linewidths=0, zorder=3,
                   label=f"{short_set[st]} (n = {int(sub.n.iloc[0]):,})")
        if j:
            for yi, (_, r) in zip(y, sub.iterrows()):
                ax.text(1.0, yi, f"p = {fp(r.p_fisher)}", va="center", ha="left", fontsize=4.8, color=INK2)
        rows.append(sub.reset_index().assign(set=st))
    among = pr.drop_duplicates("metric").set_index("metric")["among"]
    ax.set_yticks(np.arange(len(metrics)))
    among_name = {"annotated_exon": "exact GENCODE", "fully_coding": "fully coding"}
    ax.set_yticklabels([m if among[m] == "all" else f"{m}\n(of {among_name.get(among[m], among[m])} exons)"
                        for m in metrics],
                       fontsize=5.5, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.17)
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.spines["bottom"].set_bounds(0, 0.75)
    ax.set_xlabel("share of cassette exons (95% Wilson CI)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1, handletextpad=0.1)
    data = pd.concat(rows, ignore_index=True)[["metric", "among", "set", "n", "k", "proportion", "ci_lo", "ci_hi",
                                               "odds_ratio", "p_fisher", "q_fisher"]]
    data["q_family"] = "BH across 6 metrics x 4 comparison sets (Fisher tests vs all measured exons)"
    save(fig, OUT, "S6_exon_consequences", data=data, caption="""
**What do the cassette exons the models rely on do to the encoded protein?**

Cassette exons in the top 100 of any tumour-vs-normal model (blue) or any class of the cancer-type model (dark
blue), compared with a random 6,000 of all measured exons (grey).
- Length divisible by 3 is computed from the coordinates.
- The remaining features need an exact GENCODE v45 exon. Codons are taken from the best transcript (MANE >
  canonical > other) and overlapped with InterPro Pfam/CDD domains, MobiDB-lite disordered regions and TMHMM
  helices (D-17).

p: two-sided Fisher exact test against all measured exons, raw p. BH q is in the CSV, with its family. Exons
skipped or included more in tumours are in `cassette_exon_protein_features.csv`.
""")


def S7() -> None:
    d = BIO / "04_events"
    rec = pd.read_csv(d / "recurrent_events.csv")
    P = pd.read_csv(d / "paired_change_top_events_all_cancers.csv.gz")
    picks = rec[rec.n_cancers_same_sign_as_majority >= 0.8 * rec.n_cancers_measured].drop_duplicates("symbol").head(6)
    pm = ProteinMap(CFG)
    fig, axes = plt.subplots(len(picks), 2, figsize=(DOUBLE * 0.82, 1.15 * len(picks)),
                             gridspec_kw={"width_ratios": [2.4, 1.0], "hspace": 0.9, "wspace": 0.25})
    meds, local = [], []
    for i, (_, r) in enumerate(picks.iterrows()):
        ax = axes[i, 0]
        dd = P[(P.view == r["view"]) & (P.feature_id == r.feature_id)].set_index("cancer").reindex(CANCERS)
        x = np.arange(len(CANCERS))
        ax.vlines(x, dd.median_normal * 100, dd.median_tumour * 100, color=GRID, linewidth=1.2)
        ax.scatter(x, dd.median_normal * 100, s=8, color=NORMAL, linewidths=0, zorder=3, label="normal")
        ax.scatter(x, dd.median_tumour * 100, s=8, color=TUMOUR, linewidths=0, zorder=3, label="tumour")
        topc = set(r.cancers.split(","))
        ax.set_xticks(x)
        ax.set_xticklabels(CANCERS, rotation=90, fontsize=5)
        for t in ax.get_xticklabels():
            t.set_color(INK if t.get_text() in topc else MUTED)
        ax.set_ylabel("median ×100", fontsize=6)
        info = ""
        if r["view"] == "se":
            a, b = (int(v) for v in r.region.split(":")[1].split("-"))
            c = pm.exon_consequence(r.gene_id, a, b)
            if c.get("coding"):
                info = (f"; {b - a + 1} nt, {'in frame' if c['frame_preserving'] else 'frameshift'}, "
                        f"aa {c['aa_start']}-{c['aa_end']}/{c['protein_length_aa']}"
                        + (f", {c['domains']}" if c["domains"] else "") + (", disordered" if c["in_disordered"] else ""))
        ax.set_title(f"{gene_label(r.symbol)}  {VIEW_NAMES[r['view']]}  {short(r.region)}{info}", loc="left", fontsize=6,
                     color=INK)
        if i == 0:
            ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.25), ncol=2, handletextpad=0.1)
        meds.append(dd.reset_index()[["cancer", "median_normal", "median_tumour", "median_paired_delta", "n_pairs",
                                      "p_paired"]].assign(symbol=r.symbol, view=r["view"], region=r.region,
                                                          in_top_100=lambda z: z.cancer.isin(topc)))
        # SHAP dependence, held-out samples of the cancer where the event ranks highest
        ax = axes[i, 1]
        c = r.best_cancer
        tab = task_table(CFG, "tumor_status", c)
        run = CFG.run_dir("tumor_status", r["view"], c)
        vals, shs, ys = [], [], []
        for f in range(CFG.study["splits"]["outer_folds"]):
            z = np.load(run / f"fold{f}" / "shap.npz")
            fids = list(z["feature_ids"].astype(str))
            if r.feature_id not in fids:
                continue
            j = fids.index(r.feature_id)
            vals.append(feature_values(CFG, r["view"], tab["file_id"].iloc[z["rows"]], [r.feature_id])[:, 0] * 100)
            shs.append(z["shap"][:, j])
            ys.append(tab["y"].to_numpy()[z["rows"]])
        v, s, yv = np.concatenate(vals), np.concatenate(shs), np.concatenate(ys)
        for lab, colr in [(0, NORMAL), (1, TUMOUR)]:
            m = yv == lab
            ax.scatter(v[m], s[m], s=2, color=colr, alpha=0.6, linewidths=0, rasterized=True)
        ax.axhline(0, color=INK2, linewidth=0.4)
        ax.set_title(c, loc="left", fontsize=6, color=INK)
        ax.set_xlabel("value ×100", fontsize=6, labelpad=1)
        ax.set_ylabel("SHAP", fontsize=6)
        local.append(pd.DataFrame({"symbol": r.symbol, "view": r["view"], "cancer": c, "value_x100": v, "shap": s,
                                   "is_tumour": yv}))
    save(fig, OUT, "S7_event_portraits", data=pd.concat(meds, ignore_index=True),
         local_data=(LOCAL, pd.concat(local, ignore_index=True)), caption="""
**What do the most shared tumour events look like, and how does the model read them?**

The six most recurrent events (S5) whose paired change has the same sign in at least 80% of the cancers where
measured, one per gene.

Left: the median value in matched normal (orange) and tumour (blue) samples of each cancer. Cancers whose top 100
contains the event are in black, the others in grey. The heading gives the host gene, view, coordinates and, for
cassette exons, the protein consequence (GENCODE v45; D-17).

Right: SHAP dependence in the cancer where the event ranks highest. Each held-out sample's event value is plotted
against its SHAP value on the tumour log-odds. The value is the PSI; missing values were median-imputed by the
model.

The CSV gives the per-cancer medians and paired Wilcoxon p (raw). Per-sample values stay local (D-12).
""")


# ============================================================ 05 regulators


# ============================================================ 11 event cards


def S9() -> None:
    d = BIO / "11_event_cards"
    cards = pd.read_csv(d / "tumour_status_cards.csv")
    trk = pd.read_csv(d / "tumour_status_protein_tracks.csv")
    show = cards[cards["aa_start"].notna()].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(DOUBLE * 0.84, 0.42 * len(show) + 1.0))
    rows = []
    for i, r in enumerate(show.itertuples()):
        y, L = len(show) - 1 - i, r.protein_length_aa
        ax.add_patch(plt.Rectangle((0, y - 0.11), 1, 0.22, facecolor=GRID, edgecolor="none"))
        f = trk[trk["feature_id"].eq(r.feature_id)]
        for g in f[f["source_db"].isin(("pfam", "cdd"))].itertuples():
            ax.add_patch(plt.Rectangle((g.aa_start / L, y - 0.17), (g.aa_end - g.aa_start + 1) / L, 0.34,
                                       facecolor=MUTED, edgecolor="none", alpha=0.55))
        for g in f[f["source_db"].eq("mobidblite")].itertuples():
            ax.plot([g.aa_start / L, g.aa_end / L], [y - 0.28] * 2, color=INK2, linewidth=1.1,
                    solid_capstyle="butt")
        colr = TUMOUR if r.majority_sign == "up in tumour" else NORMAL
        a, b = r.aa_start / L, (r.aa_end + 1) / L
        ax.add_patch(plt.Rectangle((a, y - 0.24), max(b - a, 0.004), 0.48, facecolor=colr, edgecolor="none",
                                   zorder=3))
        ax.text(1.03, y, f"{int(L):,} aa", fontsize=5.2, color=INK2, va="center")
        note = (f"{int(r.aa_end - r.aa_start + 1)} aa terminal" if pd.isna(r.swap_nt) else
                f"{int(r.swap_nt):+d}".replace("-", "\u2212") +
                f" nt, {'in frame' if r.swap_in_frame else 'frameshift'}")
        dom = str(r.domains).split(";")[0].strip() if isinstance(r.domains, str) and r.domains else ""
        if len(dom) > 50:                                   # cut between words, never mid-word
            dom = dom[:50].rsplit(" ", 1)[0].rstrip(" ,;") + "\u2026"
        # an empty feature table for the transcript is not the same as "no domain here"
        tail = f" · {dom}" if dom else ("" if len(f) else " · no InterPro annotation")
        ax.text(1.19, y, note + tail, fontsize=5.2, color=INK if dom else INK2, va="center")
        rows.append({"symbol": r.symbol, "view": r.view, "drawn": "alternative segment",
                     "protein_length_aa": L, "aa_start": r.aa_start, "aa_end": r.aa_end,
                     "swap_nt": r.swap_nt, "swap_in_frame": r.swap_in_frame, "domains": r.domains,
                     "majority_sign": r.majority_sign})
        for g in f[f["source_db"].isin(("pfam", "cdd", "mobidblite"))].itertuples():
            rows.append({"symbol": r.symbol, "view": r.view,
                         "drawn": "disordered" if g.source_db == "mobidblite" else "Pfam/CDD domain",
                         "protein_length_aa": L, "aa_start": g.aa_start, "aa_end": g.aa_end,
                         "domains": getattr(g, "name", None)})
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels([f"{r.symbol} {VIEW_NAMES[r.view]}" for r in show.iloc[::-1].itertuples()], fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1.72)
    ax.set_ylim(-0.6, len(show) - 0.4)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["N", "", "half", "", "C"])
    ax.spines["bottom"].set_bounds(0, 1)
    ax.set_xlabel("position along the protein", x=0.29)
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    hand = [Line2D([], [], color=TUMOUR, linewidth=4, label="segment, higher in tumour"),
            Line2D([], [], color=NORMAL, linewidth=4, label="segment, lower in tumour"),
            Line2D([], [], color=MUTED, linewidth=4, alpha=0.55, label="Pfam/CDD domain"),
            Line2D([], [], color=INK2, linewidth=1.1, label="disordered")]
    ax.legend(handles=hand, frameon=False, fontsize=5.8, ncol=4, loc="lower left",
              bbox_to_anchor=(0, 1.0), handlelength=1.2, columnspacing=1.2)
    save(fig, OUT, "S9_event_protein", data=pd.DataFrame(rows), caption=f"""
**Where in the protein does each event's alternative segment sit?**

The {len(show)} recurrent events of S5 whose exon belongs to a protein-coding transcript (of ten, by number of
models; `11_event_cards/tumour_status_cards.csv` holds all ten with the rest of what is known about them). Each bar
is one protein, scaled to its own length (printed at the right), from its N terminus to its C terminus: grey blocks
are Pfam/CDD domains of the chosen transcript, the line underneath marks disordered regions (MobiDB-lite), and the
coloured block is the segment the event adds or removes -- blue where the segment is higher in tumours, orange
where it is lower (D-17, GENCODE v45). Three of the transcripts have no InterPro annotation at all, which is said
so in the margin: a bare bar there means the reference is silent, not that the segment misses every domain.

The annotation is what the swap does to the reading frame. For a mutually exclusive pair that is the difference
between the two exons, not either length, because the swap replaces one with the other. A first or last exon
replaces one end of the protein rather than shifting a frame, so it is annotated by the length of the segment it
contributes instead. A retained intron has no protein coordinates and is not drawn.

These are annotations of one transcript per gene, not isoform quantification: the models see the event, not the
protein, and no protein was measured here.
""")


# ============================================================ S10 events that read another gene
def S10() -> None:
    fig, data = overlap_figure("tumour status")
    save(fig, OUT, "S10_overlap_audit", data=data, caption=overlap_caption("tumour status", data, "S10"))


FIGS = {f.__name__: f for f in [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10]}

if __name__ == "__main__":
    run(FIGS)
