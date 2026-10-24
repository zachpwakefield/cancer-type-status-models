"""Step 36: the subtype narrative, in figures (house style, docs/FIGURE_GUIDELINES.md).

One question per figure, in the order the story is told in results/SUBTYPE.md:
    U1  how well does each modality recover the subtype, cohort by cohort, and does fusing help?
    U2  which subtypes are recovered, and by which modality?
    U3  when a histology model is wrong, is the error molecularly coherent?
    U4  does the biology the literature predicts show up, in the data and in the models?
    U5  do the splicing models and the expression model of a cohort read the same genes?
    U6  are the strongest splicing events splicing, or readouts of a neighbouring gene?
    U7  do the exons that mark MSI in one cancer mark it in the other two?
    U8  do the splicing regulators the literature names differ between subtypes (D-46)?
    U9  is the subtype signal tied to the tissue-source site (D-47)?
    U10 does genome sequence say which exons MSI tumours skip (D-52)?
    U11 do the longest acceptor repeats move the other way (D-53, D-54)?
    U12 are the BRCA models' genotype-like splicing events inherited? Tumour-normal pairs (D-57, A-13, A-14)
    U13 which HER2 did each modality learn: amplification or the HER2-enriched programme (D-58)?
    U14 do splice-site mutations and SF3B1 hotspots leave their marks in BRCA splicing (D-59)?
    U15 does splicing's HER2 lead survive ERBB2 copy number without the ABSOLUTE ceiling (D-60)?
    U16 does more mismatch-repair damage mean more skipping of repeat exons (D-61)?
    U17 does the branch point's distance from the 3' splice site set the direction of the MSI shift (D-63)?
    U18 does the sequence code find mismatch-repair-deficient tumours in 30 cancers never used (D-66)?
    U19 is the code about the branch point, and about mismatch repair (D-65, D-67)?
    U20 is the code an infiltration readout (D-68)?
    U21 does the code follow MSI calls that do not come from MC3: MANTIS, MSIsensor (D-69)?
    U22 is the rising arm real: specific, not a global axis, tract-dependent, dose, PTB (D-70)?
    U23 does logistic regression on the repeat exons, fitted in three cancers, transfer (D-71, A-19)?
    U24 how does the repeat-exon model compare with the exome MSI callers (D-72, A-19)?
    U25 where do the model and the exome measures agree and disagree (D-73)?
    U26 does the splicing model transfer better than expression models trained the same way (D-74)?
    U27 how large is the MSI splicing shift, and how is it spread over exons and tumours (D-75)?
    U28 what does including a rising-arm exon do to its transcript (D-76)?

Captions are figure legends: "**Figure Un. Title.**", then (a), (b) ... for the panels. `legends`
copies each legend (without its bullet list of tests) under the figure in results/SUBTYPE.md.

Reads results/subtype/summary/ (step 34), results/subtype/interpret/ (step 35),
results/subtype/regulators/ (step 37), results/subtype/site_probe/ (step 38),
results/subtype/msi_cis/ (steps 42, 50, 51, 52, 54, 55) and results/subtype/brca/ (steps 46-49). Every figure
is SVG + PDF + 400-dpi PNG with `<stem>.csv` of the plotted values and `<stem>.caption.md`.
Usage: python scripts/36_subtype_figures.py [U1 U3 ... legends]        (default: all, then legends)
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from _common import log
from _figbase import inject_legends, locus_label, run, tidy
from _figbase import legend as _legend

from ctsm.biology import bh
from ctsm.figstyle import (ACCENT, ACCENT_DARK, ACCENT_LIGHT, DIVERGE, DOUBLE, GRID, INK, INK2, MUTED, SEQ, SINGLE,
                           colourbar, fp, heat_axes, save)
from ctsm.subtype import TASK, cohort_specs, subtype_config

SCFG = subtype_config()
SPECS = cohort_specs(SCFG.spec)
SUM = SCFG.results_dir / TASK / "summary"
INT = SCFG.results_dir / TASK / "interpret"
OUT = SCFG.results_dir / "figures" / "subtype"
AS_NAMES = ["AFE", "ALE", "HIT", "SE", "A5SS", "A3SS", "MXE", "RI"]
GROUPS = [("histology", lambda s: s["label_set"] == "histology"),
          ("genomic or IHC label", lambda s: s["label_type"] in ("genomic", "IHC")),
          ("expression-defined label", lambda s: s["label_type"] == "expression-defined")]
TIER_MARK = {"weak_label": " (weak)", "control": " (control)"}
MINUS = "−"


def signed(v: float, nd: int = 3) -> str:
    return f"{v:+.{nd}f}".replace("-", MINUS)


def cohort_order(perf: pd.DataFrame) -> list[list[str]]:
    have = set(perf["cohort"])
    return [[c for c, s in SPECS.items() if test(s) and c in have] for _, test in GROUPS]


def row_label(c: str, n: int) -> str:
    return f"{c}{TIER_MARK.get(SPECS[c]['tier'], '')} ({n:,})"


# ============================================================ U1 performance
def U1() -> None:
    perf = pd.read_csv(SUM / "performance.csv")
    tests = pd.read_csv(SUM / "fusion_tests.csv").set_index(["cohort", "comparison"])
    groups = cohort_order(perf)
    heights = [max(len(g), 1) for g in groups]
    fig, axes = plt.subplots(len(groups), 1, figsize=(DOUBLE * 0.78, 0.17 * sum(heights) + 1.0),
                             gridspec_kw={"height_ratios": heights, "hspace": 0.35}, sharex=True)
    rows = []
    for ax, (gname, _), cohorts in zip(axes, GROUPS, groups):
        for i, c in enumerate(cohorts):
            p = perf[perf.cohort == c].set_index("model")
            n, C = int(p["n"].iloc[0]), int(p["n_classes"].iloc[0])
            av = p.loc[p.index.isin(AS_NAMES), "balanced_accuracy"]
            ax.plot([av.min(), av.max()], [i, i], color=GRID, linewidth=1.6, zorder=1, solid_capstyle="round")
            ax.scatter(av, [i] * len(av), s=7, color=MUTED, linewidths=0, zorder=2)
            ax.plot([1 / C, 1 / C], [i - 0.32, i + 0.32], color=INK2, linewidth=0.7, zorder=1)
            for name, col, s, z in (("GEX", INK, 16, 4), ("AS ensemble", ACCENT, 16, 5), ("GEX + AS", ACCENT_DARK, 16, 6)):
                if name in p.index:
                    ax.scatter(p.at[name, "balanced_accuracy"], i, s=s, color=col, linewidths=0.4, edgecolors="white",
                               zorder=z)
            for name in p.index:
                rows.append({"group": gname, "cohort": c, "n": n, "n_classes": C, "chance": 1 / C, "model": name,
                             **p.loc[name, ["balanced_accuracy", "balanced_accuracy_lo", "balanced_accuracy_hi",
                                            "macro_f1", "macro_auroc_ovr"]].to_dict()})
            key = (c, "GEX + AS vs GEX")
            if key in tests.index:
                t = tests.loc[key]
                ax.text(1.03, i, f"{signed(t.delta_balanced_accuracy)}", transform=ax.get_yaxis_transform(),
                        fontsize=5.6, va="center", ha="left", color=INK)
                ax.text(1.17, i, fp(t.p_mcnemar), transform=ax.get_yaxis_transform(), fontsize=5.6, va="center",
                        ha="left", color=INK)
        ax.set_yticks(range(len(cohorts)))
        ax.set_yticklabels([row_label(c, int(perf.loc[perf.cohort == c, "n"].iloc[0])) for c in cohorts],
                           fontsize=5.8, color=INK)
        ax.set_ylim(len(cohorts) - 0.5, -0.6)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_title(gname, fontsize=6.2, color=INK2, loc="left", pad=2)
        if ax is not axes[-1]:
            ax.tick_params(axis="x", length=0)
            ax.spines["bottom"].set_visible(False)
    axes[0].text(1.03, -1.35, "fused\n− GEX", transform=axes[0].get_yaxis_transform(), fontsize=5.6, ha="left",
                 va="bottom", color=INK2)
    axes[0].text(1.17, -1.35, "McNemar\np", transform=axes[0].get_yaxis_transform(), fontsize=5.6, ha="left",
                 va="bottom", color=INK2)
    axes[-1].set_xlim(0.15, 1.005)
    axes[-1].set_xlabel("balanced accuracy (out-of-fold)")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", color=MUTED, ms=2.6, label="single splicing view"),
               Line2D([], [], marker="o", ls="", color=INK, ms=3.6, label="GEX"),
               Line2D([], [], marker="o", ls="", color=ACCENT, ms=3.6, label="AS ensemble"),
               Line2D([], [], marker="o", ls="", color=ACCENT_DARK, ms=3.6, label="GEX + AS (1:1)"),
               Line2D([], [], marker="|", ls="", color=INK2, ms=5, label="chance (1/classes)")]
    axes[-1].legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -2.9 / max(heights[-1], 1)), ncol=5,
                    fontsize=5.6, handletextpad=0.2, columnspacing=0.9)
    data = pd.DataFrame(rows).merge(
        pd.read_csv(SUM / "fusion_tests.csv")[lambda d: d.comparison == "GEX + AS vs GEX"][
            ["cohort", "delta_balanced_accuracy", "delta_ba_lo", "delta_ba_hi", "only_gex_correct",
             "only_model_correct", "p_mcnemar", "q_mcnemar", "q_family"]], on="cohort", how="left")
    save(fig, OUT, "U1_performance", data=data, caption=CAPTIONS["U1"](perf, pd.read_csv(SUM / "fusion_tests.csv")))


# ============================================================ U2 recall by class
def U2() -> None:
    rc = pd.read_csv(SUM / "recall_by_class.csv")
    rc = rc[rc.model.isin(["GEX", "AS ensemble", "GEX + AS"])]
    perf = pd.read_csv(SUM / "performance.csv")
    groups = cohort_order(perf)
    panels = [groups[0], groups[1] + groups[2]]
    n_rows = [sum(len(SPECS[c]["classes"]) for c in g) for g in panels]
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 0.105 * max(n_rows) + 0.9), gridspec_kw={"wspace": 0.95})
    rows = []
    for ax, cohorts in zip(axes, panels):
        labels, y = [], 0
        for c in cohorts:
            d = rc[rc.cohort == c]
            classes = list(d.drop_duplicates("class")["class"])
            for cl in classes:
                e = d[d["class"] == cl].set_index("model")
                vals = e["recall"]
                ax.plot([vals.min(), vals.max()], [y, y], color=GRID, linewidth=1.2, zorder=1)
                for name, col, z in (("GEX", INK, 3), ("AS ensemble", ACCENT, 4), ("GEX + AS", ACCENT_DARK, 5)):
                    if name in e.index:
                        ax.scatter(e.at[name, "recall"], y, s=9, color=col, linewidths=0.3, edgecolors="white",
                                   zorder=z)
                        rows.append({"cohort": c, "class": cl, "model": name, "support": int(e.at[name, "support"]),
                                     "recall": e.at[name, "recall"], "recall_lo": e.at[name, "recall_lo"],
                                     "recall_hi": e.at[name, "recall_hi"]})
                labels.append(f"{c} · {cl} ({int(e['support'].iloc[0])})")
                y += 1
            ax.axhline(y - 0.5, color=GRID, linewidth=0.4)
            y += 0.0
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=5.0, color=INK)
        ax.set_ylim(len(labels) - 0.5, -0.6)
        ax.set_xlim(-0.02, 1.02)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_xlabel("recall (out-of-fold)")
    from matplotlib.lines import Line2D
    axes[1].legend(handles=[Line2D([], [], marker="o", ls="", color=c, ms=3.2, label=n) for n, c in
                            (("GEX", INK), ("AS ensemble", ACCENT), ("GEX + AS", ACCENT_DARK))],
                   loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5.6, handletextpad=0.2)
    for x, lab in [(0.0, "a"), (0.40, "b")]:
        fig.text(x, 0.995, lab, fontsize=8, fontweight="bold", va="top")
    save(fig, OUT, "U2_recall_by_class", data=pd.DataFrame(rows), caption=CAPTIONS["U2"](pd.DataFrame(rows)))


# ============================================================ U5 AS vs GEX genes
def U5() -> None:
    ov = pd.read_csv(INT / "as_vs_gex_genes.csv")
    perf = pd.read_csv(SUM / "performance.csv")
    order = [c for g in cohort_order(perf) for c in g if c in set(ov.cohort)]
    piv = ov.pivot(index="cohort", columns="view", values="overlap").reindex(index=order, columns=AS_NAMES)
    exp = ov.pivot(index="cohort", columns="view", values="expected").reindex(index=order, columns=AS_NAMES)
    fig, ax = plt.subplots(figsize=(3.1, 0.15 * len(order) + 0.9))
    vmax = max(10, float(np.nanmax(piv.to_numpy())))
    im = ax.imshow(piv.to_numpy(dtype=float), cmap=SEQ, vmin=0, vmax=vmax, aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.iat[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=5.0,
                        color="white" if v > 0.6 * vmax else INK)
    heat_axes(ax, [row_label(c, int(perf.loc[perf.cohort == c, "n"].iloc[0])) for c in order], AS_NAMES,
              row_size=5.6, col_size=5.8)
    ax.xaxis.tick_top()
    colourbar(fig, im, ax, "genes in both top-50 lists")
    ov = ov.assign(q_hypergeom=bh(ov["p_hypergeom"].to_numpy()),
                   q_family="BH across every cohort x AS-view pair of this figure")
    save(fig, OUT, "U5_as_vs_gex_genes", data=ov, caption=CAPTIONS["U5"](ov, exp))


# ============================================================ U4 literature markers
def U4() -> None:
    mk = pd.read_csv(INT / "literature_markers.csv")
    order = [c for g in cohort_order(pd.read_csv(SUM / "performance.csv")) for c in g]
    mk["o"] = mk["cohort"].map({c: i for i, c in enumerate(order)})
    mk = mk.sort_values(["o"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.57, 0.135 * len(mk) + 0.8))
    for i, r in mk.iterrows():
        a = r.auroc_expected
        col = ACCENT_DARK if a >= 0.5 else INK2
        sig = bool(r.p_mwu < 0.05)
        ax.plot([0.5, a], [i, i], color=GRID, linewidth=1.0, zorder=1)
        ax.scatter(a, i, s=12, facecolors=col if sig else "white", edgecolors=col, linewidths=0.7, zorder=3)
        gx = "–" if not np.isfinite(r.get("gex_rank_all", np.nan)) else f"{int(min(r.gex_rank_all, np.nanmin([r.gex_rank_class_a, r.gex_rank_class_b, np.inf])))}"
        asr = "–" if not np.isfinite(r.get("as_best_rank", np.nan)) else f"{int(r.as_best_rank)} {r.as_best_view}"
        if r.gene == "1p/19q dosage":
            gx = f"{int(r.gex_top_on_1p19q)}/{int(r.gex_top_n)} on 1p/19q"
            asr = "–"           # an arm-level dosage has no splicing event to rank
        ax.text(1.03, i, gx, transform=ax.get_yaxis_transform(), fontsize=5.4, va="center", color=INK)
        ax.text(1.36, i, asr, transform=ax.get_yaxis_transform(), fontsize=5.4, va="center", color=INK)
    ax.axvline(0.5, color=INK2, linewidth=0.5)
    ax.set_yticks(range(len(mk)))
    ax.set_yticklabels([f"{r.cohort} · {r.gene}: {r.higher_in} > {r.vs}" for r in mk.itertuples()], fontsize=5.3,
                       color=INK)
    ax.set_ylim(len(mk) - 0.5, -0.8)
    ax.set_xlim(min(0.4, float(mk["auroc_expected"].min()) - 0.05), 1.0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("AUROC of the gene's expression, expected class over the other")
    ax.text(1.03, -1.1, "best GEX rank", transform=ax.get_yaxis_transform(), fontsize=5.4, color=INK2, va="bottom")
    ax.text(1.36, -1.1, "best AS rank", transform=ax.get_yaxis_transform(), fontsize=5.4, color=INK2, va="bottom")
    save(fig, OUT, "U4_literature_markers", data=mk.drop(columns="o"), caption=CAPTIONS["U4"](mk))


# ============================================================ U6 named events with the audit
def audit_status(r: dict) -> str:
    """The D-41 audit of one event, as printed next to it, with the event's correlation with its own
    gene's expression."""
    host = "host ρ –" if not np.isfinite(r["rho_host"]) else f"host ρ {signed(r['rho_host'], 2)}"
    if r["overlap_class"] == "none":
        return f"no overlap; {host}"
    other = locus_label(r["other_gene"], r["region"])
    if r["overlap_class"] == "readthrough":
        return f"in readthrough {other}; {host}"
    if not np.isfinite(r["rho_other_gene"]):
        return f"{r['overlap_class']} {other}, not measured; {host}"
    return (f"{r['overlap_class']} {other} ρ {signed(r['rho_other_gene'], 2)}; {host}"
            + (f" → {other} readout" if r["tracks_other_gene"] else ""))


def U6() -> None:
    ev = pd.read_csv(INT / "event_audit.csv")
    pick = U6_pick(ev)
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.69, 0.24 * len(pick) + 0.75))
    rows = []
    for i, r in enumerate(pick.to_dict("records")):
        cls = list(SPECS[r["cohort"]]["classes"])
        meds = {c: r.get(f"median_{c}", np.nan) for c in cls}
        meds = {c: v for c, v in meds.items() if np.isfinite(v)}
        others = [v for c, v in meds.items() if c != r["class"]]
        ax.plot([min(meds.values()), max(meds.values())], [i, i], color=GRID, linewidth=1.0, zorder=1)
        ax.scatter(others, [i] * len(others), s=12, facecolors="white", edgecolors=INK2, linewidths=0.6, zorder=2)
        ax.scatter(meds[r["class"]], i, s=16, color=ACCENT_DARK, linewidths=0, zorder=3)
        ax.text(meds[r["class"]], i - 0.28, r["class"], fontsize=5.0, color=ACCENT_DARK, ha="center", va="bottom")
        status = audit_status(r)
        ax.text(1.02, i, status, transform=ax.get_yaxis_transform(), fontsize=5.2, va="center",
                color="#b0412c" if r["tracks_other_gene"] else INK)
        rows.append({**{k: v for k, v in r.items() if not k.startswith("median_")},
                     **{f"median_{c}": meds.get(c, np.nan) for c in cls}, "audit_status": status})
    ax.set_yticks(range(len(pick)))
    ax.set_yticklabels([f"{r.cohort} · {locus_label(r.symbol, r.region)} {r.view}" for r in pick.itertuples()],
                       fontsize=5.6, color=INK)
    ax.set_ylim(len(pick) - 0.5, -0.9)
    ax.set_xlim(-1.05, 1.05)
    ax.axvline(0, color=GRID, linewidth=0.5, zorder=0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("median over a class's patients: PSI or usage (0 to 1); HIT index (−1 to 1)")
    data = pd.DataFrame(rows).drop(columns=["class_index"], errors="ignore")
    save(fig, OUT, "U6_event_audit", data=data, caption=CAPTIONS["U6"](data, ev))


def U6_pick(ev: pd.DataFrame) -> pd.DataFrame:
    """Per cohort, the single splicing event with the largest class-specific importance (any class,
    any of the 8 views), for the histology cohorts and the genomic/IHC cohorts. Fixed before the
    audit's results were read (D-41)."""
    ev = ev[ev["selected_by"].str.startswith("top")]
    top = ev.sort_values("mean_abs_shap_class", ascending=False).drop_duplicates("cohort")
    keep = [c for c, s in SPECS.items() if s["label_type"] != "expression-defined"]
    top = top[top.cohort.isin(keep)]
    top["o"] = top["cohort"].map({c: i for i, c in enumerate(SPECS)})
    return top.sort_values("o").drop(columns="o").reset_index(drop=True)


# ============================================================ U3 histology errors by molecular class
U3_PAIRS = [("LGG", "LGG_IDH"), ("UCEC", "UCEC_TCGA"), ("STAD", "STAD_TCGA"), ("PCPG", "PCPG_mRNA")]
SHORT = {"Oligodendroglioma": "oligo", "Astrocytoma": "astro", "Endometrioid": "endometrioid", "Serous": "serous",
         "Intestinal": "intestinal", "Diffuse": "diffuse", "Pheochromocytoma": "adrenal", "Paraganglioma": "extra-adrenal"}
MOL_SHORT = {"IDHmut-codel": "IDH-mut\ncodel", "IDHmut-non-codel": "IDH-mut\nnon-codel", "IDHwt": "IDH-wt",
             "Kinase signaling": "kinase", "Pseudohypoxia": "pseudo-\nhypoxia", "Wnt-altered": "Wnt",
             "Cortical admixture": "cortical\nadmixture"}


def U3() -> None:
    cf = pd.read_csv(INT / "cross_family_errors.csv")
    pairs = [(h, m) for h, m in U3_PAIRS if ((cf.histology_cohort == h) & (cf.molecular_cohort == m)).any()]
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE * 0.82, 3.3), gridspec_kw={"wspace": 1.15, "hspace": 0.75})
    rows = []
    for ax, (h, m) in zip(axes.ravel(), pairs):
        hcls = list(SPECS[h]["classes"])
        present = set(cf.loc[cf.molecular_cohort == m, "molecular_class"])
        mcls = [c for c in SPECS[m]["classes"] if c in present]
        combos = [(t, p) for t in hcls for p in [t] + [x for x in hcls if x != t]]
        for model in ("GEX", "AS ensemble"):
            d = cf[(cf.histology_cohort == h) & (cf.molecular_cohort == m) & (cf.model == model)]
            M = np.array([[int(d[(d.true_histology == t) & (d.predicted_histology == p) & (d.molecular_class == c)]
                               .patients.sum()) for c in mcls] for t, p in combos])
            share = M / np.maximum(M.sum(axis=1, keepdims=True), 1)
            rows += [{"histology_cohort": h, "molecular_cohort": m, "model": model, "plotted": model == "GEX",
                      "true_histology": t, "called_histology": p, "molecular_class": c, "patients": int(M[i, j]),
                      "share_of_row": share[i, j]} for i, (t, p) in enumerate(combos) for j, c in enumerate(mcls)]
            if model != "GEX":
                continue
            im = ax.imshow(share, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    ax.text(j, i, f"{M[i, j]}", ha="center", va="center", fontsize=5.4,
                            color="white" if share[i, j] > 0.6 else INK)
            heat_axes(ax, [f"{SHORT.get(t, t)} → {SHORT.get(p, p)} ({int(M[i].sum())})" for i, (t, p) in enumerate(combos)],
                      [MOL_SHORT.get(c, c) for c in mcls], row_size=5.4, col_size=5.0)
            ax.xaxis.tick_top()
            ax.set_title(h, fontsize=6.4, color=INK, loc="left", pad=15)
    cb = fig.colorbar(im, ax=axes[:, -1], fraction=0.05, pad=0.04, shrink=0.5, ticks=[0, 0.5, 1])
    cb.set_label("share of the row", fontsize=5.6, color=INK2)
    cb.ax.tick_params(labelsize=5.4, length=2)
    cb.outline.set_visible(False)
    data = pd.DataFrame(rows)
    save(fig, OUT, "U3_histology_errors", data=data, caption=CAPTIONS["U3"](data))


# ============================================================ U7 MSI events across three cancers
MSI_COHORTS = ["UCEC_TCGA", "STAD_TCGA", "COAD_TCGA"]


def U7_rows(sh: pd.DataFrame) -> pd.DataFrame:
    """Every event in the top 5 of the MSI class of any of the three cohorts whose AUROC (MSI against
    the rest of the cohort) is at least 0.25 from 0.5, in the same direction, in all three cohorts."""
    d = sh[sh["class"] == "MSI"]
    w = d.pivot_table(index=["view", "feature_id"], columns="measured_in", values="auroc_class_over_rest")
    w = w.reindex(columns=MSI_COHORTS).dropna()
    keep = w[((w - 0.5).abs().min(axis=1) >= 0.25) & ((w > 0.5).all(axis=1) | (w < 0.5).all(axis=1))]
    return keep.assign(strength=(keep - 0.5).abs().min(axis=1)).sort_values("strength", ascending=False)


def U7() -> None:
    sh = pd.read_csv(INT / "shared_class_events.csv")
    ev = pd.read_csv(INT / "event_audit.csv")
    ev = ev[ev["selected_by"].str.startswith("top") & (ev["class"] == "MSI")]
    keep = U7_rows(sh)
    info = sh.drop_duplicates(["view", "feature_id"]).set_index(["view", "feature_id"])
    fig, ax = plt.subplots(figsize=(3.0, 0.2 * len(keep) + 1.0))
    A = keep[MSI_COHORTS].to_numpy(dtype=float)
    im = ax.imshow(A, cmap=DIVERGE, vmin=0, vmax=1, aspect="auto")
    rows = []
    for i, (key, r) in enumerate(keep.iterrows()):
        for j, c in enumerate(MSI_COHORTS):
            chosen = bool(((ev["cohort"] == c) & (ev["view"] == key[0]) & (ev["feature_id"] == key[1])).any())
            m = sh[(sh["class"] == "MSI") & (sh["view"] == key[0]) & (sh["feature_id"] == key[1])
                   & (sh["measured_in"] == c)].iloc[0]
            ax.text(j, i, f"{A[i, j]:.2f}" + ("•" if chosen else ""), ha="center", va="center", fontsize=5.2,
                    color="white" if abs(A[i, j] - 0.5) > 0.35 else INK)
            rows.append({"view": key[0], "feature_id": key[1], "symbol": info.at[key, "symbol"],
                         "region": info.at[key, "region"], "cohort": c, "msi_top5_of_this_cohort": chosen,
                         "n_msi": int(m.n_class), "n_rest": int(m.n_rest), "median_msi": m.median_class,
                         "median_rest": m.median_rest, "auroc_msi_over_rest": m.auroc_class_over_rest,
                         "p_mwu": m.p_mwu, "q_mwu_replication": m.q_mwu_replication, "q_family": m.q_family})
        a = ev[(ev["view"] == key[0]) & (ev["feature_id"] == key[1])].iloc[0].to_dict()
        ax.text(1.04, i, audit_status(a), fontsize=5.0, va="center", ha="left", transform=ax.get_yaxis_transform(),
                color="#b0412c" if a["tracks_other_gene"] else INK)
    labels = [f"{locus_label(info.at[k, 'symbol'], info.at[k, 'region'])} {k[0]} ({short_region(info.at[k, 'region'])})"
              for k in keep.index]
    heat_axes(ax, labels, [c.replace("_TCGA", "") for c in MSI_COHORTS], row_size=5.4, col_size=5.8)
    ax.xaxis.tick_top()
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.06, pad=0.04, shrink=0.7, ticks=[0, 0.5, 1])
    cb.set_label("AUROC, MSI over the rest of the cohort", fontsize=5.6, color=INK2)
    cb.ax.tick_params(labelsize=5.4, length=2)
    cb.outline.set_visible(False)
    data = pd.DataFrame(rows)
    save(fig, OUT, "U7_msi_events", data=data, caption=CAPTIONS["U7"](data))


def short_region(region: str) -> str:
    ch, rg = region.split(":")
    a, b = (int(x) for x in rg.split("-"))
    return f"{ch}:{a / 1e6:.2f} Mb, {b - a + 1} bp"


# ============================================================ captions
def _gaps(perf):
    """Per cohort: GEX, AS ensemble and fused balanced accuracy, and the label family."""
    w = perf.pivot_table(index="cohort", columns="model", values="balanced_accuracy")
    w = w.dropna(subset=["AS ensemble"]) if "AS ensemble" in w else w.iloc[:0]
    w["type"] = [SPECS[c]["label_type"] for c in w.index]
    w["gap"] = w["AS ensemble"] - w["GEX"]
    return w


def _cap_u1(perf, tests):
    w = _gaps(perf)
    ne = w[w["type"] != "expression-defined"]
    close = ne[ne["gap"].abs() <= 0.03]
    far = ne[ne["gap"] < -0.10].sort_values("gap")
    above = ne[ne["gap"] > 0].sort_values("gap", ascending=False)
    t = tests[tests.comparison == "GEX + AS vs GEX"]
    up = t[t.delta_ba_lo > 0]
    n = int(perf.drop_duplicates("cohort")["n"].sum())
    far_s = ", ".join(f"{c} ({signed(g, 2)})" for c, g in far["gap"].items())
    above_s = ", ".join(f"{c} ({signed(g, 3)})" for c, g in above["gap"].items())
    up_l = [f"{r.cohort} ({signed(r.delta_balanced_accuracy, 3)}, McNemar p = {fp(r.p_mcnemar)})" for r in up.itertuples()]
    up_s = " and ".join([", ".join(up_l[:-1]), up_l[-1]]) if len(up_l) > 1 else "".join(up_l)
    return f"""**How well does each modality recover the subtype, cohort by cohort, and does adding splicing to expression help?**

Out-of-fold balanced accuracy (mean per-class recall) of one model per cohort and view, {len(w)} cohorts and
{n:,} patient-cohort pairs: the eight splicing views (grey; the bar spans them), expression (GEX, black), the
AS ensemble (the mean of the eight splicing views' class probabilities, blue) and GEX + AS fused 1:1 (dark
blue). Fusion weights were fixed before any result (D-40). Tick: chance, 1/classes, which is also the
majority-class predictor's balanced accuracy. Rows: cohorts, patients in brackets, grouped by how the label was
made: histology (weak = a poorly reproducible label, control = a lineage or site contrast), genomic or IHC (DNA,
protein), and expression-defined (mRNA clusters, so GEX is near-circular there; D-44).

Columns at right: fused minus GEX balanced accuracy, and the exact McNemar test of fused against GEX on paired per-patient
correctness, raw p. McNemar weighs classes by their size and balanced accuracy weighs them equally, so the two
can disagree (THCA). The CSV holds every point's 95% patient-bootstrap CI, macro F1, macro AUROC, the
paired-bootstrap CI of the fused − GEX difference and BH q (family: the subtype cohorts, per comparison).

Outside the expression-defined labels the AS ensemble is within 0.03 of GEX in {len(close)} of {len(ne)} cohorts and
above it in {len(above)}: {above_s}. It falls more than 0.10 short in {len(far)}: {far_s}. Fusing helps
balanced accuracy with a CI above zero only in {up_s or "no cohort"}."""


def _cap_u3(data):
    g = data[data.model == "GEX"]

    def n(h, t=None, p=None, c=None):
        d = g[g.histology_cohort == h]
        for col, v in (("true_histology", t), ("called_histology", p), ("molecular_class", c)):
            if v is not None:
                d = d[d[col].isin(v if isinstance(v, list) else [v])]
        return int(d.patients.sum())
    O, A = "Oligodendroglioma", "Astrocytoma"
    nc = ["IDHmut-non-codel", "IDHwt"]
    E, S_ = "Endometrioid", "Serous"
    I, D = "Intestinal", "Diffuse"
    P, G = "Pheochromocytoma", "Paraganglioma"
    return f"""**When a histology model is wrong, is the error molecularly coherent?**

Each panel is one histology cohort's expression (GEX) model. Rows: true histology → the model's out-of-fold call,
with patients in brackets; columns: the published molecular class of the same patients (PanCanAtlas, D-43);
cells: patients, coloured by the share of the row. Only patients with both labels are counted. Counts only, no
test (D-41). The AS ensemble's table is in the CSV (`plotted` = False); it shows the same pattern.

- **LGG.** {n("LGG", O, O, "IDHmut-codel")} of {n("LGG", O, None, "IDHmut-codel")} codeleted oligodendrogliomas are called
  oligodendroglioma, but {n("LGG", O, A, nc)} of the {n("LGG", O, None, nc)} without a codeletion are called astrocytoma;
  {n("LGG", A, O, "IDHmut-codel")} of {n("LGG", A, None, "IDHmut-codel")} codeleted astrocytomas are called oligodendroglioma.
- **UCEC.** {n("UCEC", E, S_, "CN-high")} of the {n("UCEC", E, S_)} endometrioid tumours called serous are copy-number
  high, against {n("UCEC", E, None, "CN-high")} of all {n("UCEC", E)} endometrioid tumours.
- **STAD.** {n("STAD", D, D, "GS")} of {n("STAD", D, None, "GS")} genomically stable diffuse tumours are called diffuse,
  and so are {n("STAD", I, D, "GS")} of the {n("STAD", I, None, "GS")} genomically stable intestinal ones.
- **PCPG.** {n("PCPG", P, G, "Pseudohypoxia")} of the {n("PCPG", P, G)} adrenal tumours called extra-adrenal are in
  the pseudohypoxia cluster, and {n("PCPG", P, P, "Cortical admixture")} of {n("PCPG", P, None, "Cortical admixture")}
  tumours with cortical admixture are called adrenal.

The labels here were assigned by other groups over the full TCGA cohort; that is a property of the label, not a
leak, and no molecular label entered any histology model."""


def _cap_u2(data):
    w = data.pivot_table(index=["cohort", "class", "support"], columns="model", values="recall").reset_index()
    w = w[[SPECS[c]["label_type"] != "expression-defined" for c in w.cohort]]
    w["gap"] = w["AS ensemble"] - w["GEX"]
    lo = w.sort_values("gap").head(5)
    hi = w.sort_values("gap", ascending=False).head(4)
    fmt = lambda d: "; ".join(f"{r.cohort} {r['class']} ({int(r.support)}): {r.GEX:.2f} → {r['AS ensemble']:.2f}"
                              for _, r in d.iterrows())
    small = w[w.support < 40]
    small_note = " all in classes of fewer than 40 patients" if (lo.support < 40).all() else ""
    return f"""**Which subtypes does each modality recover?**

Out-of-fold recall of every class (patients in brackets) for GEX (black), the AS ensemble (blue) and GEX + AS fused
1:1 (dark blue). (a) Histology cohorts. (b) Genomic and IHC cohorts, then the expression-defined ones, in the
order of U1. 95% patient-bootstrap CIs and every single view are in the CSV and in `recall_by_class.csv`.

Splicing's five largest losses against expression (recall, GEX → AS ensemble) are{small_note}: {fmt(lo)}. Its largest gains: {fmt(hi)}. Over the {len(small)} classes with fewer than 40
patients, the AS ensemble's recall is lower than GEX's in {int((small.gap < 0).sum())}."""


def _cap_u4(mk):
    m = mk[mk.gene != "1p/19q dosage"]
    ok = m[(m.auroc_expected > 0.5) & (m.p_mwu < 0.05)]
    gx = m[m.gex_in_top.astype(bool)]
    asx = m[m.as_in_top.astype(bool)]
    d = mk[mk.gene == "1p/19q dosage"]
    dl = "; ".join(f"{r.cohort}: {int(r.gex_top_on_1p19q)} of the top {int(r.gex_top_n)} GEX genes on 1p or 19q "
                   f"against {r.gex_top_expected_on_1p19q:.1f} expected (hypergeometric p = {fp(r.p_hypergeom_1p19q)})"
                   for r in d.itertuples())
    return f"""**Does the biology the literature predicts show up, in the data and in the models?**

Each row is a marker fixed in `config/subtype.yaml` before any model was fitted (D-41, D-45; sources in the CSV),
with the class the literature expects it to be higher in. Dot: AUROC of the gene's log2 CPM for that class over
the other (0.5 = no difference), computed on all patients of the two classes with no model involved; filled =
Mann-Whitney raw p < 0.05, open = not. Columns: the gene's best rank in the GEX model (over all patients or within
either class) and the best rank of any of its events across the eight splicing models, with the view (– = no
event of the gene was selected). For LGG and LGG_IDH, "1p/19q dosage" is the median over the arm's genes of each
gene's log2 CPM minus its cohort median, 1p and 19q averaged; the GEX column counts the model's top 100 genes
that lie on 1p or 19q.

{len(ok)} of {len(m)} gene markers differ in the expected direction (p < 0.05). {len(gx)} rank in the top 100 of
their GEX model; {len(asx)} have an event in the top 100 of a splicing model. {dl}."""


def _cap_u5(ov, exp):
    n_sig = int((ov.q_hypergeom < 0.05).sum())
    top = ov.sort_values("overlap", ascending=False).head(3)
    ts = "; ".join(f"{r.cohort} {r.view}: {int(r.overlap)} ({r.shared_genes.replace(';', ', ')})" for r in top.itertuples())
    return f"""**Do the splicing models and the expression model of a cohort read the same genes?**

For each cohort (patients in brackets) and splicing view: how many genes are in both the top 50 genes of the
splicing model (a gene ranked by its best event) and the top 50 of the cohort's GEX model, counted within the
genes both views could select in that cohort (D-41). Under independent draws the expected count is
{ov.expected.min():.1f}–{ov.expected.max():.1f} (median {ov.expected.median():.1f}); the CSV holds it with the fold
enrichment, the hypergeometric raw p and BH q (family named in the CSV).

The median overlap is {ov.overlap.median():.0f} gene of 50 and {int((ov.overlap == 0).sum())} of {len(ov)} pairs share none;
{n_sig} pairs exceed chance at q < 0.05. The largest overlaps are {ts}."""


def _cap_u6(data, ev):
    n = len(data)
    tr = data[data.tracks_other_gene.astype(bool)]
    um = data[(data.overlap_class.isin(["antisense", "same-strand"])) & data.rho_other_gene.isna()]
    hostc = data[data.rho_host.abs() >= 0.5]
    trs = ", ".join(f"{r.cohort} {r.symbol} → {r.other_gene}" for r in tr.itertuples())
    return f"""**Are the strongest splicing events splicing, or readouts of a neighbouring gene?**

For each histology and genomic/IHC cohort, the one splicing event with the largest class-specific mean |SHAP|
across the eight splicing models (rule fixed before the audit was run, D-41), labelled by cohort, host gene and
view. Dots: the event's median over the patients of each class, PSI or usage (0 to 1) or HIT index (−1 to 1);
filled and named = the class the event ranks for. Columns at right: the D-27 audit. The event's segment is intersected with
every other GENCODE v45 gene's exons (≥ 50 bp): antisense, same-strand, inside a readthrough model, or no overlap.
ρ: Spearman correlation over the cohort's patients with the other gene's log2 CPM and with the host's. An
overlapping event that follows the other gene (|ρ| ≥ 0.5 and above |ρ host|) is a readout of it (red). "Not
measured": the other gene is not in the expression matrix, so the audit cannot clear the event.

{len(tr)} of {n} are readouts of a neighbour ({trs}); {len(um)} overlap a gene the audit cannot measure;
{len(hostc)} follow their own gene's expression with |ρ| ≥ 0.5. The full audit of the top 5 events per class and
view, and of every literature marker's best event, is `results/subtype/interpret/event_audit.csv`."""


def _cap_u7(data):
    d = data
    k = d.drop_duplicates(["view", "feature_id"])
    chosen = d.groupby(["view", "feature_id"]).msi_top5_of_this_cohort.sum()
    c = d[(d.symbol == "CCNDBP1") & (d.view == "SE")]
    cs = "; ".join(f"{r.cohort.replace('_TCGA', '')} {r.median_msi:.2f} vs {r.median_rest:.2f}" for r in c.itertuples())
    nrep = d[~d.msi_top5_of_this_cohort]
    return f"""**Do the exons that mark MSI in one cancer mark it in the other two?**

Every splicing event in the top 5 of the MSI class of the UCEC, STAD or COAD model (any view) whose AUROC, MSI
against the rest of the cohort, is at least 0.25 from 0.5 in the same direction in all three cohorts: {len(k)}
events (AUROC in each cell; red = higher in MSI). • = in that cohort's MSI top 5; a cell without • is a
replication in a cohort whose model did not pick the event. No model is involved in the AUROC; the patients with a
value in each cell (MSI and rest) are in the CSV. Columns at right: the D-27 audit, as in U6. Exploratory (A-4): the step was added after the
same CCNDBP1 exon was seen in the MSI top 5 of all three models; it then tests every event of every shared class.

{int((chosen >= 2).sum())} of the {len(k)} events were chosen independently by at least two of the three models.
CCNDBP1's 57-bp exon is included at a median {cs} (MSI vs rest). All {len(nrep)} replication cells have Mann-Whitney
raw p < {max(nrep.p_mwu.max(), 1e-300):.1e}; BH q over every cross-cohort test of every shared class is in the CSV."""


# ============================================================ U8 splicing regulators
REG = SCFG.results_dir / TASK / "regulators"
SITE = SCFG.results_dir / TASK / "site_probe"


def U8() -> None:
    q = pd.read_csv(REG / "predicted_factor_by_class.csv")
    b = pd.read_csv(REG / "sarc_event_by_ptbp1_quintile.csv")
    pr = pd.read_csv(REG / "predictions.csv").set_index("id")
    fig = plt.figure(figsize=(DOUBLE, 3.0))
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1.05, 1.05], wspace=0.95, hspace=0.9)
    rows = []
    for j, rid in enumerate(["R1", "R2", "R3", "R4"]):
        ax = fig.add_subplot(gs[j // 2, j % 2])
        d = q[q.id == rid].reset_index(drop=True)
        order = list(d.sort_values("predicted").index)              # "higher" first, then "lower"
        d = d.loc[order].reset_index(drop=True)
        for i, r in d.iterrows():
            col = ACCENT if r.predicted == "lower" else INK
            ax.plot([r.q25, r.q75], [i, i], color=col, linewidth=1.6, solid_capstyle="round", alpha=0.45)
            ax.scatter(r["median"], i, s=14, color=col, linewidths=0.4, edgecolors="white", zorder=3)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels([f"{c} ({n})" for c, n in zip(d["class"], d["n"])], fontsize=5.6, color=INK)
        ax.set_ylim(len(d) - 0.5, -0.7)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_title(f"{d.cohort.iloc[0]} · {d.gene.iloc[0]}", fontsize=6.2, color=INK2, loc="left", pad=2)
        ax.set_xlabel("log2 CPM", fontsize=6)
        lo, hi = d.q25.min(), d.q75.max()
        ax.set_xlim(lo - 0.15 * (hi - lo), hi + 0.15 * (hi - lo))
        rows.append(d.assign(panel="a", auroc=pr.at[rid, "estimate"], p=pr.at[rid, "p"], q=pr.at[rid, "q"],
                             verdict=pr.at[rid, "verdict"]))
    for j, rid in enumerate(["R5", "R6"]):
        ax = fig.add_subplot(gs[:, 2 + j])
        d = b[b.id == rid]
        for grp, col in [("LMS", ACCENT), ("other sarcomas", MUTED)]:
            g = d[d.group == grp].sort_values("ptbp1_quintile")
            x = g.ptbp1_quintile + (0.06 if grp == "LMS" else -0.06)
            ax.vlines(x, g.q25_psi, g.q75_psi, color=col, linewidth=0.8, alpha=0.6)
            ax.plot(x, g.median_psi, color=col, linewidth=1.0, marker="o", ms=2.8, mec="white", mew=0.3)
            ax.text(5.3, g.median_psi.iloc[-1], grp.replace(" ", "\n"), fontsize=5.4, color=col, va="center",
                    ha="left", linespacing=0.9)
        ax.set_xticks(range(1, 6))
        ax.set_xlim(0.6, 6.5)
        ax.spines["bottom"].set_bounds(1, 5)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel("PTBP1 quintile", fontsize=6)
        ax.set_ylabel("PSI" if j == 0 else "", fontsize=6)
        ax.set_title(d.event_name.iloc[0].replace(" (MXE)", "").replace(" mutually exclusive exons", " MXE"),
                     fontsize=6.2, color=INK2, loc="left", pad=2)
        rows.append(d.assign(panel="b", rho_partial=pr.at[rid, "estimate"], p=pr.at[rid, "p"], q=pr.at[rid, "q"],
                             verdict=pr.at[rid, "verdict"], lms_minus_rest_psi=pr.at[rid, "lms_minus_rest_psi"]))
    fig.text(0.02, 0.99, "a", fontsize=8, fontweight="bold", va="top")
    fig.text(0.535, 0.99, "b", fontsize=8, fontweight="bold", va="top")
    data = pd.concat(rows, ignore_index=True)
    data["q_family"] = "BH across R1-R6 (config/followup.yaml)"
    save(fig, OUT, "U8_regulators", data=data, caption=CAPTIONS["U8"](pr, pd.read_csv(REG / "predictions_sensitivity.csv"),
                                                                     pd.read_csv(REG / "factor_scan_summary.csv")))


# ============================================================ U9 site probe
def U9() -> None:
    w = pd.read_csv(SITE / "site_probe_summary.csv")
    per = pd.read_csv(SITE / "site_probe_by_cohort.csv").sort_values("drop_as_median").reset_index(drop=True)
    tests = pd.read_csv(SITE / "site_probe_tests.csv")
    conf = pd.read_csv(SITE / "site_class_confounding.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE * 0.86, 0.16 * len(per) + 0.9),
                                 gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.45})
    rows = []
    for i, r in per.iterrows():
        v = w[(w.cohort == r.cohort) & (w.view != "gex")]
        ax.plot([v["drop"].min(), v["drop"].max()], [i, i], color=GRID, linewidth=1.6, solid_capstyle="round", zorder=1)
        ax.scatter(v["drop"], [i] * len(v), s=6, color=MUTED, linewidths=0, zorder=2)
        ax.scatter(r.drop_gex, i, s=15, color=INK, linewidths=0.4, edgecolors="white", zorder=4)
        ax.scatter(r.drop_as_median, i, s=15, color=ACCENT, linewidths=0.4, edgecolors="white", zorder=5)
        for _, x in w[w.cohort == r.cohort].iterrows():
            rows.append({"panel": "a", "cohort": r.cohort, "view": x.view, "case_grouped": x["case-grouped"],
                         "site_held_out": x["site-held-out"], "drop": x["drop"]})
        rows.append({"panel": "a", "cohort": r.cohort, "view": "AS-view median", "drop": r.drop_as_median})
    ax.axvline(0, color=INK2, linewidth=0.6, zorder=0)
    ax.set_yticks(range(len(per)))
    ax.set_yticklabels([f"{c} ({n}; {s} sites)" for c, n, s in zip(per.cohort, per.n, per.n_sites)], fontsize=5.6,
                       color=INK)
    ax.set_ylim(len(per) - 0.5, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("loss of balanced accuracy when sites are held out")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=MUTED, ms=2.4, label="single splicing view"),
                       Line2D([], [], marker="o", ls="", color=INK, ms=3.6, label="GEX"),
                       Line2D([], [], marker="o", ls="", color=ACCENT, ms=3.6, label="median of 8 splicing views")],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5.6, handletextpad=0.2, columnspacing=0.8)
    for col, name, c in [("drop_gex", "GEX", INK), ("drop_as_median", "median of 8 splicing views", ACCENT)]:
        bx.scatter(per.cramers_v_site_class, per[col], s=12, color=c, linewidths=0.4, edgecolors="white", zorder=3)
        for _, r in per.iterrows():
            rows.append({"panel": "b", "cohort": r.cohort, "view": name, "cramers_v_site_class": r.cramers_v_site_class,
                         "drop": r[col]})
    lab = per[per.drop_as_median >= 0.07]
    for _, r in lab.iterrows():
        crowded = ((lab.cramers_v_site_class - r.cramers_v_site_class).between(0.001, 0.12)
                   & (lab.drop_as_median - r.drop_as_median).abs().lt(0.006)).any()
        bx.text(r.cramers_v_site_class + (-0.012 if crowded else 0.012), r.drop_as_median, r.cohort, fontsize=5.2,
                color=INK2, va="center", ha="right" if crowded else "left")
    bx.axhline(0, color=INK2, linewidth=0.6, zorder=0)
    bx.set_xlabel("site–class association (Cramér's V)")
    bx.set_ylabel("loss of balanced accuracy")
    ax.text(-0.42, 1.08, "a", transform=ax.transAxes, fontsize=8, fontweight="bold")
    bx.text(-0.2, 1.08, "b", transform=bx.transAxes, fontsize=8, fontweight="bold")
    data = pd.DataFrame(rows)
    save(fig, OUT, "U9_site_probe", data=data, caption=CAPTIONS["U9"](per, tests, conf))


# ============================================================ U10, U11 MSI cis repeats
MSI_DIR = SCFG.results_dir / TASK / "msi_cis"
MSI_LOCAL = SCFG.base_data_dir / "msi"
COH_STYLE = {"UCEC_TCGA": ("#c98a2b", "o"), "STAD_TCGA": (ACCENT, "s"), "COAD_TCGA": (ACCENT_DARK, "^")}


def _short(c: str) -> str:
    return c.replace("_TCGA", "")


def _stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def U10() -> None:
    pr = pd.read_csv(MSI_DIR / "predictions.csv")
    top = pd.read_csv(MSI_DIR / "model_top_events.csv").set_index("cohort")
    au = pd.read_csv(MSI_DIR / "score_aurocs.csv")
    sc = pd.read_parquet(MSI_LOCAL / "sample_scores.parquet")
    coh = list(COH_STYLE)
    fig = plt.figure(figsize=(DOUBLE, 2.45))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1.55, 1.15], wspace=0.55)
    rows = []
    # a. share of the MSI model's top 50 SE events with an acceptor T-run >= 10
    ax = fig.add_subplot(gs[0])
    for i, c in enumerate(coh):
        r = top.loc[c]
        col = COH_STYLE[c][0]
        p6 = pr[(pr.id == "M6") & (pr.cohort == c)].iloc[0]
        for k, (lab, num, den, colr) in enumerate([("all measured exons", r.other_R, r.other_n, MUTED),
                                                    ("top 50 MSI events", r.top_R, r.top_n, col)]):
            y = i * 2.6 + k
            v = 100 * num / den
            ax.barh(y, v, height=0.8, color=colr, linewidth=0)
            ax.text(v + 2, y, f"{v:.0f}%", va="center", fontsize=5.4, color=INK)
            rows.append({"panel": "a", "cohort": c, "set": lab, "n": int(den), "n_acceptor_T10": int(num),
                         "percent": v, "odds_ratio": r.odds_ratio, "p": p6.p, "q": p6.q,
                         "q_family": "M (msi_cis, 19 tests)"})
        ax.text(105, i * 2.6 + 0.5, _stars(p6.p), va="center", fontsize=6, color=INK)
    ax.set_yticks([i * 2.6 + k for i in range(3) for k in (0, 1)])
    ax.set_yticklabels([f"{_short(c)} {lab}" for c in coh for lab in ("all", "top 50")], fontsize=5.6)
    ax.set_ylim(7.0, -0.7)
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 50, 100])
    ax.spines["bottom"].set_bounds(0, 100)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("% with an acceptor T-run ≥ 10", fontsize=6)
    # b. the sample score by group
    bx = fig.add_subplot(gs[1])
    pos, ticks, labs = 0, [], []
    rng = np.random.default_rng(0)
    local = []
    for c in coh:
        groups = ["MSI", "POLE", "rest"] if c == "UCEC_TCGA" else ["MSI", "rest"]
        xs = {}
        for g in groups:
            v = sc[(sc.cohort == c) & (sc.group == g)]["score_acc"].to_numpy()
            col = COH_STYLE[c][0] if g == "MSI" else MUTED
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            lo, hi = v[v >= q1 - 1.5 * (q3 - q1)].min(), v[v <= q3 + 1.5 * (q3 - q1)].max()
            bx.add_patch(plt.Rectangle((pos - 0.3, q1), 0.6, q3 - q1, facecolor="white", edgecolor=col, linewidth=0.7,
                                       zorder=2))
            bx.vlines(pos, lo, q1, color=col, linewidth=0.6, zorder=1)
            bx.vlines(pos, q3, hi, color=col, linewidth=0.6, zorder=1)
            bx.scatter(pos + rng.uniform(-0.2, 0.2, len(v)), v, s=2.2, color=col, alpha=0.45, linewidths=0, zorder=3)
            bx.scatter(pos, med, marker="D", s=10, color=col, edgecolors="white", linewidths=0.4, zorder=4)
            rows.append({"panel": "b", "cohort": c, "set": g, "n": len(v), "median": med, "q1": q1, "q3": q3,
                         "whisker_low": lo, "whisker_high": hi})
            local.append(sc[(sc.cohort == c) & (sc.group == g)][["cohort", "group", "score_acc"]])
            xs[g] = pos
            ticks.append(pos)
            labs.append(f"{g}\n({len(v)})")
            pos += 1
        top_y = sc[sc.cohort == c]["score_acc"].quantile(0.995)
        m1 = pr[(pr.id == "M1") & (pr.cohort == c)].iloc[0]
        h = 0.012
        bx.plot([xs["MSI"], xs["MSI"], xs["rest"], xs["rest"]], [top_y + h, top_y + 1.4 * h, top_y + 1.4 * h, top_y + h],
                color=INK2, linewidth=0.5)
        bx.text((xs["MSI"] + xs["rest"]) / 2, top_y + 1.6 * h, _stars(m1.p), ha="center", va="bottom", fontsize=6)
        rows.append({"panel": "b", "cohort": c, "set": "MSI vs rest (M1)", "auroc": m1.effect + 0.5, "p": m1.p,
                     "q": m1.q, "q_family": "M (msi_cis, 19 tests)"})
        if c == "UCEC_TCGA":
            m5 = pr[pr.id == "M5"].iloc[0]
            y5 = top_y + 3.2 * h
            bx.plot([xs["MSI"], xs["MSI"], xs["POLE"], xs["POLE"]], [y5, y5 + 0.4 * h, y5 + 0.4 * h, y5],
                    color=INK2, linewidth=0.5)
            bx.text((xs["MSI"] + xs["POLE"]) / 2, y5 + 0.6 * h, _stars(m5.p), ha="center", va="bottom", fontsize=6)
            rows.append({"panel": "b", "cohort": c, "set": "MSI vs POLE (M5)", "auroc": m5.effect + 0.5, "p": m5.p,
                         "q": m5.q, "q_family": "M (msi_cis, 19 tests)"})
        bx.text((xs["MSI"] + xs["rest"]) / 2, -0.205, _short(c), ha="center", va="top", fontsize=6.2, color=INK,
                transform=bx.get_xaxis_transform())
        pos += 0.6
    bx.axhline(0, color=GRID, linewidth=0.6, zorder=0)
    bx.set_xticks(ticks)
    bx.set_xticklabels(labs, fontsize=5.4)
    bx.tick_params(axis="x", length=0)
    bx.set_xlim(-0.6, pos - 0.9)
    bx.set_ylabel("repeat score (percentile)", fontsize=6)
    # c. AUROC by acceptor run length, and the deep-window control
    cx = fig.add_subplot(gs[2])
    bins = [("bin_6", "6–7"), ("bin_8", "8–9"), ("bin_10", "10–11"), ("bin_12", "12–13"), ("bin_14", "≥14")]
    for c in coh:
        col, mk = COH_STYLE[c]
        d = au[au.cohort == c].set_index("score")
        y = [d.at[k, "auroc_msi_vs_rest"] for k, _ in bins]
        cx.plot(range(len(bins)), y, color=col, linewidth=0.9, marker=mk, ms=3, mec="white", mew=0.3, label=_short(c))
        cx.scatter(len(bins) + 0.8, d.at["deep", "auroc_msi_vs_rest"], marker=mk, s=11, facecolors="white",
                   edgecolors=col, linewidths=0.8, zorder=3)
        for k, lab in bins + [("deep", "deep ≥10")]:
            rows.append({"panel": "c", "cohort": c, "set": lab, "auroc": d.at[k, "auroc_msi_vs_rest"],
                         "p": d.at[k, "p"], "q_family": "exploratory (dose_auroc); deep = M4 comparator"})
    cx.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    cx.set_xticks(list(range(len(bins))) + [len(bins) + 0.8])
    cx.set_xticklabels([b.replace("–", "–\n") for _, b in bins] + ["deep\n≥10"], fontsize=5.4)
    cx.set_xlabel("acceptor T-run (nt)", fontsize=6)
    cx.set_ylabel("AUROC, score lower in MSI", fontsize=6)
    cx.set_ylim(0.3, 0.85)
    cx.legend(loc="lower left", fontsize=5.4, handlelength=1.4, borderaxespad=0.2)
    for a_, lab, x in [(ax, "a", -0.62), (bx, "b", -0.17), (cx, "c", -0.25)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    rows += [{"panel": "b", "set": "reference line", "median": 0.0},
             {"panel": "c", "set": "reference line (chance)", "auroc": 0.5}]
    data = pd.DataFrame(rows)
    save(fig, OUT, "U10_msi_repeats", data=data, caption=CAPTIONS["U10"](pr, top, au),
         local_data=(SCFG.base_data_dir / "figure_data", pd.concat(local, ignore_index=True)))


def U11() -> None:
    ex = pd.read_parquet(MSI_LOCAL / "exon_effects.parquet")
    u7 = pd.read_csv(INT / "shared_class_events.csv")
    drop = set(u7.loc[u7["class"].eq("MSI") & u7["view"].eq("SE"), "feature_id"])
    a3 = pd.read_parquet(MSI_LOCAL / "a3ss_effects.parquet")
    q = pd.read_csv(MSI_DIR / "a3ss_predictions.csv")
    n_tests = pd.read_csv(MSI_DIR / "long_tract_predictions.csv")
    sh = pd.read_csv(MSI_DIR / "shared_repeat_exons.csv")
    coh = list(COH_STYLE)
    fig = plt.figure(figsize=(DOUBLE, 2.9))
    outer = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.25], wspace=0.95)
    gs = outer[0].subgridspec(1, 2, wspace=0.75)
    rows = []
    # a. share of SE repeat exons more included in MSI, by run length
    ax = fig.add_subplot(gs[0])
    bins = [(10, 14, "10–14"), (15, 19, "15–19"), (20, 999, "≥20")]
    for j, c in enumerate(coh):
        col, mk = COH_STYLE[c]
        e = ex[(ex.cohort == c) & ~ex.feature_id.isin(drop)]
        ys, xs = [], []
        for i, (lo, hi, lab) in enumerate(bins):
            v = e[(e.acc_T >= lo) & (e.acc_T <= hi)]["delta_adj"]
            k, n = int((v > 0).sum()), len(v)
            ph = k / n
            z = 1.96
            den = 1 + z * z / n
            mid = (ph + z * z / (2 * n)) / den
            half = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
            x = i + (j - 1) * 0.16
            ax.vlines(x, mid - half, mid + half, color=col, linewidth=0.7)
            xs.append(x)
            ys.append(ph)
            rows.append({"panel": "a", "cohort": c, "set": f"SE acceptor T-run {lab}", "n": n, "n_up": k, "share_up": ph,
                         "ci_low": mid - half, "ci_high": mid + half, "q_family": "exploratory (D-53 table)"})
        ax.plot(xs, ys, color=col, linewidth=0.8, marker=mk, ms=3, mec="white", mew=0.3, label=_short(c))
    ax.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    ax.set_xticks(range(3))
    ax.set_xticklabels([b[2] for b in bins], fontsize=5.6)
    ax.set_xlabel("acceptor T-run (nt), cassette exons", fontsize=6)
    ax.set_ylabel("share more included in MSI", fontsize=6)
    ax.set_ylim(0.2, 0.85)
    ax.legend(loc="upper left", fontsize=5.4, handlelength=1.4, borderaxespad=0.2)
    # b. A3SS replication: oriented shift for moderate and long runs
    bx = fig.add_subplot(gs[1])
    for j, c in enumerate(coh):
        col, mk = COH_STYLE[c]
        e = a3[(a3.cohort == c) & (a3.side != 0)]
        for i, (pid, lo, hi, lab) in enumerate([("Q1", 10, 19, "10–19"), ("Q2", 20, 999, "≥20")]):
            v = e[(e.run >= lo) & (e.run <= hi)]["oriented"]
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            x = i + (j - 1) * 0.27
            bx.vlines(x, q1, q3, color=col, linewidth=0.8, alpha=0.6)
            bx.scatter(x, med, marker=mk, s=13, color=col, edgecolors="white", linewidths=0.3, zorder=3)
            t = q[(q.id == pid) & (q.cohort == c)].iloc[0]
            bx.text(x, 1.0, _stars(t.p), ha="center", va="bottom", fontsize=5.0, color=INK,
                    transform=bx.get_xaxis_transform())
            rows.append({"panel": "b", "cohort": c, "set": f"A3SS run {lab} ({pid})", "n": len(v), "median": med,
                         "q1": q1, "q3": q3, "p": t.p, "q": t.q, "q_family": "Q (msi_a3ss_replication, 6 tests)"})
    bx.axhline(0, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    bx.set_xticks([0, 1])
    bx.set_xticklabels(["10–19", "≥20"], fontsize=5.6)
    bx.set_xlim(-0.55, 1.55)
    bx.set_xlabel("T-run at one acceptor (nt), A3SS", fontsize=6)
    bx.set_ylabel("use of that acceptor, MSI − rest", fontsize=6)
    # c. the exons with the largest shared shifts
    cx = fig.add_subplot(outer[1])
    s1 = sh.assign(a=sh.mean_delta_adj.abs()).sort_values("a", ascending=False).drop_duplicates("symbol")
    down = s1[s1.mean_delta_adj < 0].nsmallest(9, "mean_delta_adj")
    up = s1[s1.mean_delta_adj > 0].nlargest(7, "mean_delta_adj")
    extra = s1[s1.symbol.eq("MRE11") & ~s1.symbol.isin(down.symbol)]
    pick = pd.concat([up.sort_values("mean_delta_adj", ascending=False), down.sort_values("mean_delta_adj", ascending=False),
                      extra]).reset_index(drop=True)
    for i, r in pick.iterrows():
        vals = [r[f"delta_adj_{c}"] for c in coh]
        cx.plot([min(vals), max(vals)], [i, i], color=GRID, linewidth=1.4, solid_capstyle="round", zorder=1)
        for c, v in zip(coh, vals):
            col, mk = COH_STYLE[c]
            cx.scatter(v, i, marker=mk, s=11, color=col, edgecolors="white", linewidths=0.3, zorder=3)
            rows.append({"panel": "c", "cohort": c, "set": r.symbol, "feature_id": r.feature_id, "acc_T": int(r.acc_T),
                         "exon_len": int(r.exon_len), "delta_adj": v, "q": r[f"q_{c}"],
                         "q_family": "BH over the cohort's exons (per-exon Mann-Whitney)"})
    cx.axvline(0, color=INK2, linewidth=0.6, zorder=0)
    sep = len(up) - 0.5
    cx.axhline(sep, color=GRID, linewidth=0.6)
    if len(extra):
        cx.axhline(len(pick) - 1.5, color=GRID, linewidth=0.6)
    cx.set_yticks(range(len(pick)))
    cx.set_yticklabels([f"{r.symbol} · T{int(r.acc_T)} · {int(r.exon_len)} nt" for r in pick.itertuples()], fontsize=5.4)
    cx.set_ylim(len(pick) - 0.5, -0.6)
    cx.tick_params(axis="y", length=0)
    cx.spines["left"].set_visible(False)
    cx.set_xlabel("MSI − rest (adjusted percentile)", fontsize=6)
    for a_, lab, x in [(ax, "a", -0.3), (bx, "b", -0.36), (cx, "c", -0.62)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    rows += [{"panel": "a", "set": "reference line", "share_up": 0.5},
             {"panel": "b", "set": "reference line", "median": 0.0},
             {"panel": "c", "set": "reference line", "delta_adj": 0.0}]
    data = pd.DataFrame(rows)
    save(fig, OUT, "U11_msi_long_tracts", data=data, caption=CAPTIONS["U11"](data, q, n_tests, pick))


def _cap_u8(pr, sens, scan):
    r = pr
    def a(i):
        return f"AUROC {r.at[i, 'estimate']:.2f}, p = {fp(r.at[i, 'p'])}, q = {fp(r.at[i, 'q'])}"
    s6 = sens[(sens.id == "R6") & sens.adjusted_for.str.contains("mean")].iloc[0]
    s5 = sens[(sens.id == "R5") & sens.adjusted_for.str.contains("mean")].iloc[0]
    return f"""**Do the splicing regulators the literature names differ between subtypes, and track the subtype's exons?**

Six predictions fixed before the data were looked at by class (config/followup.yaml, D-46). a, a factor's log2 CPM
in each class: dot = median, bar = interquartile range, patients in brackets; blue = the class predicted lower.
AUROC = P(predicted-higher patient > predicted-lower patient), two-sided Mann–Whitney.

- R1, ESRP1 lower in biphasic than epithelioid mesothelioma: {a('R1')}, not resolved (ESRP1 is barely expressed
  in either class).
- R2, ESRP2, same contrast: {a('R2')}.
- R3, ESRP1 lower in the mesenchymal HNSC class than in the other three: {a('R3')}.
- R4, PTBP1 lower in leiomyosarcoma (LMS) than in the other sarcomas: {a('R4')}.

b, the event's PSI by PTBP1 quintile, the quintiles cut within LMS and within the other sarcomas separately (dot = median, bar =
interquartile range, about 20–26 patients per point). Test: partial Spearman of PTBP1 with the PSI over all 230
sarcomas, adjusting for class, stromal and immune scores; predicted sign = opposite to the LMS shift.

- R5, ACTN1: LMS shift {signed(r.at['R5', 'lms_minus_rest_psi'], 2)}, ρ = {signed(r.at['R5', 'estimate'], 2)}, p = {fp(r.at['R5', 'p'])}
  (within LMS alone ρ = {signed(r.at['R5', 'rho_within_lms'], 2)}). Agrees with Southby et al. 1999: PTB represses the
  smooth-muscle exon.
- R6, TPM1: LMS shift {signed(r.at['R6', 'lms_minus_rest_psi'], 2)}, ρ = {signed(r.at['R6', 'estimate'], 2)}, p = {fp(r.at['R6', 'p'])}.
  It matches the pre-specified sign, but that sign misread the source (A-5): in TPM1, PTB and MBNL1 repress exon
  3, the non-smooth-muscle exon (Gooding et al. 1998, 2013). So the data run against the published mechanism here.

Sensitivity (A-6, exploratory): adding each sample's mean MXE PSI (D-24) leaves ρ at {signed(s5.rho_ptbp1, 2)} and
{signed(s6.rho_ptbp1, 2)}; PTBP1 ranks {s5.rank_ptbp1_of} factors for ACTN1 and {s6.rank_ptbp1_of} for TPM1, where MBNL1
leads ({s6.top_factors.split(';')[0].replace('-', MINUS)}). BH q over R1–R6 is in the CSV. The scan of all 70 factors in every
cohort is `results/subtype/regulators/factor_by_class.csv`."""


def _cap_u9(per, tests, conf):
    p1 = tests[tests.id == "P1"].iloc[0]
    p2 = tests[tests.id == "P2"].reset_index(drop=True)
    from ctsm.biology import bh as _bh
    q = _bh(tests["p"].to_numpy())
    top = per.sort_values("drop_as_median", ascending=False).head(4)
    ts = ", ".join(f"{r.cohort} {signed(r.drop_as_median)} (GEX {signed(r.drop_gex)})" for r in top.itertuples())
    return f"""**Is the subtype signal tied to the tissue-source site that sent the sample?**

The probe of D-20 on the subtype cohorts (D-47): one untuned classifier per cohort and view (top-2,000 variance,
median impute, standardise, L2 logistic regression C = 0.05, balanced), cross-validated with patients grouped and
again with tissue-source sites held out (5 folds, 3 seeds, pooled out-of-fold calls). a, the loss of balanced
accuracy, grouped minus site-held-out, per cohort (patients; sites): grey = each splicing view, blue = their
median, black = GEX. Ordered by the splicing median. SARC is not probed: myxofibrosarcoma comes from 4 sites (rule
≥ 5). b, the loss against Cramér's V between site and class, computed from the labels alone.

- P1, the splicing views lose more than GEX: in {int(p1.k)} of {int(p1.n)} cohorts (median difference
  {signed(p1.estimate)}); sign test p = {fp(p1.p)}, q = {fp(q[0])} over D-47's three tests: not resolved.
- P2, the loss rises with Cramér's V: Spearman ρ = {signed(p2.estimate[0], 2)} for GEX (p = {fp(p2.p[0])}) and
  {signed(p2.estimate[1], 2)} for the splicing median (p = {fp(p2.p[1])}): not supported.

Largest splicing losses: {ts}. The probe's accuracies are those of the untuned classifier, not of the subtype
models; per view, scheme and seed in `results/subtype/site_probe/`."""


def _cap_u10(pr, top, au):
    g = pr.set_index(["id", "cohort"])
    def m(pid, c):
        return g.loc[(pid, c)]
    coh = list(COH_STYLE)
    m6 = "; ".join(f"{_short(c)} {int(top.at[c, 'top_R'])} of {int(top.at[c, 'top_n'])} against "
                   f"{100 * top.at[c, 'other_R'] / top.at[c, 'other_n']:.1f}% (OR {top.at[c, 'odds_ratio']:.0f}, "
                   f"p = {fp(m('M6', c).p)})" for c in coh)
    m1 = "; ".join(f"{_short(c)} {m('M1', c).effect + 0.5:.2f} (p = {fp(m('M1', c).p)})" for c in coh)
    m4 = "; ".join(f"{_short(c)} {signed(m('M4', c).effect, 2)} (95% CI {signed(m('M4', c).ci_low, 2)} to "
                   f"{signed(m('M4', c).ci_high, 2)})" for c in coh)
    k1 = ", ".join(f"{_short(c)} {m('K1', c).effect + 0.5:.2f}" for c in coh)
    m7 = pr[pr.id == "M7"]
    m5 = m("M5", "UCEC_TCGA")
    pole = au[(au.cohort == "UCEC_TCGA") & (au.score == "acc")].iloc[0]
    return f"""**Does genome sequence say which exons MSI tumours skip?**

Pre-specified before any PSI was read by repeat (config/followup.yaml: msi_cis, D-52). For every cassette exon, the
longest sense-strand T-run overlapping intron positions −1..−50 before it was read from GRCh38; no label or model chose
the exons. MSI against the other classes (POLE apart) of UCEC_TCGA, STAD_TCGA and COAD_TCGA.

a, share of exons with a run ≥ 10 among all exons measured in the cohort and among the top 50 SE events of the cohort's
MSI model (mean |SHAP| over MSI patients); Fisher exact (M6). b, the repeat score of each patient: each exon's PSI as a
percentile within the cohort, mean over exons with a run ≥ 10 minus the mean over exons with ≤ 5 (weighted to the same
mean-PSI deciles, which cancels the global PSI axis of D-24). Box = interquartile range, diamond = median, whiskers =
1.5 × IQR, dots = patients (n under each box); Mann–Whitney (M1, M5). c, AUROC of the same score built from each
run-length bin, and (hollow) from T-runs ≥ 10 at −201..−400, deep in the intron, with the acceptor ≤ 6.

- M6, the trained MSI models chose repeat exons: {m6}.
- M1, the score is lower in MSI: AUROC {m1}. Endometrium does not separate on this score.
- M4, acceptor window against deep window: AUROC difference {m4}; supported in all three, endometrium included.
- M5, MSI against POLE (UCEC): AUROC {m5.effect + 0.5:.2f}, p = {fp(m5.p)}, not resolved; POLE against the rest
  {pole.auroc_pole_vs_rest:.2f}. With MSI itself at {pole.auroc_msi_vs_rest:.2f} in this cohort, the POLE control
  cannot discriminate here.
- Not drawn: K1, MRE11 exon 5 (T11) is less included in MSI in all three (AUROC {k1}); M7, the repeat exons' shifts
  agree between cohorts (Spearman {', '.join(f'{v:.2f}' for v in m7.effect)}).

Stars: * p < 0.05, ** p < 0.01, *** p < 0.001, ns otherwise. q over family M (19 tests) is in the CSV with every verdict
in `results/subtype/msi_cis/predictions.csv`. Per-patient scores are kept locally (D-12)."""


def _cap_u11(data, q, n_tests, pick):
    qq = q.set_index(["id", "cohort"])
    def t(pid):
        return "; ".join(f"{_short(c)} median {signed(qq.at[(pid, c), 'effect'], 3)}, n = {int(qq.at[(pid, c), 'n'])}, "
                         f"p = {fp(qq.at[(pid, c), 'p'])}" for c in COH_STYLE)
    a = data[data.panel == "a"]
    lo = a[a.set.str.contains("10–14|15–19")]["share_up"]
    hi = a[a.set.str.contains("≥20")]["share_up"]
    sup = (q.verdict == "supported").sum()
    up = pick[pick.mean_delta_adj > 0]
    dn = pick[(pick.mean_delta_adj < 0) & (pick.symbol != "MRE11")]
    return f"""**Do the longest acceptor repeats move the other way in MSI?**

a, exploratory (read after D-52): SE exons with an acceptor T-run ≥ 10 (the nine MSI events of §8 left out), share
whose background-adjusted MSI shift is positive, by run length; bars = Wilson 95% interval. Runs of 10–19 lean toward
skipping (share {lo.min():.2f}–{lo.max():.2f}), runs ≥ 20 toward inclusion ({hi.min():.2f}–{hi.max():.2f}) in all three
cancers. b, the replication fixed after that (D-54) on A3SS events, which this analysis had not read: one acceptor with
the run, the other ≤ 6; shift in the use of the acceptor carrying the run (negative = used less in MSI); marker =
median, bar = interquartile range; Wilcoxon signed-rank against 0.

- Q1, runs 10–19, predicted < 0: {t('Q1')}.
- Q2, runs ≥ 20, predicted > 0: {t('Q2')}.
- {sup} of 6 clear q < 0.05 (family Q); all 6 medians fall in the predicted direction. The long-run sets are small.

c, exons measured in all three cohorts with the largest mean shift, one per gene: the {len(up)} largest increases (top
block) and the {len(dn)} largest decreases, then MRE11 exon 5, the known MSI target. Label = gene · acceptor T-run ·
exon length. The increases carry runs of {int(up.acc_T.min())}–{int(up.acc_T.max())} nt; the decreases
{int(dn.acc_T.min())}–{int(dn.acc_T.max())} nt.

The pre-specified mechanism was not supported: D-53 predicted that the shift would depend on the exon's baseline
inclusion, as hnRNP C repression of long U-tracts would make it (Zarnack et al. 2013), and N1–N3 are not resolved in any
cohort. What the data show is a length rule, not its cause. Whether the rMATS PSI here counts exon-body reads is not
recorded; if it does, unspliced pre-mRNA could add to an apparent increase. Stars: * p < 0.05, ** p < 0.01,
*** p < 0.001, ns otherwise; q in the CSV."""

# ============================================================ BRCA follow-up (D-57..D-59)
BRCA_DIR = SCFG.results_dir / TASK / "brca"
BRCA_LOCAL = SCFG.base_data_dir / "brca"


def _step(name: str):
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location(name[:2], Path(__file__).with_name(name))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _box(ax, x, v, col, width=0.32, dots=True, rng=None, s=2.2):
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    lo, hi = v[v >= q1 - 1.5 * (q3 - q1)].min(), v[v <= q3 + 1.5 * (q3 - q1)].max()
    ax.add_patch(plt.Rectangle((x - width / 2, q1), width, q3 - q1, facecolor="white", edgecolor=col, linewidth=0.7,
                               zorder=2))
    ax.vlines(x, lo, q1, color=col, linewidth=0.6, zorder=1)
    ax.vlines(x, q3, hi, color=col, linewidth=0.6, zorder=1)
    if dots:
        ax.scatter(x + rng.uniform(-width / 2.6, width / 2.6, len(v)), v, s=s, color=col, alpha=0.45, linewidths=0, zorder=3)
    ax.scatter(x, med, marker="D", s=10, color=col, edgecolors="white", linewidths=0.4, zorder=4)
    return {"n": len(v), "median": med, "q1": q1, "q3": q3, "whisker_low": lo, "whisker_high": hi}


GROUP_NAMES = ["low", "mid", "high"]


def U12() -> None:
    gd = BRCA_DIR / "germline"
    chk = pd.read_csv(gd / "germline_check.csv")
    prof = pd.read_csv(gd / "top_events_genotype_like_profile.csv")
    pt = pd.read_csv(gd / "pair_tables_examples.csv")
    pr = pd.read_csv(gd / "predictions.csv").set_index("id")
    pc = pd.read_csv(gd / "predictions_germline_confirmed.csv").set_index("id")
    top = chk.gene_id.isin(prof.gene_id)
    fig = plt.figure(figsize=(DOUBLE, 2.15))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1], wspace=0.75)
    rows = []
    # a, b. pair tables: each normal's PSI group against its own tumour's
    axes = []
    for i, sym in enumerate(["CAST", "TMEM101"]):
        ax = fig.add_subplot(gs[i])
        sub = pt[pt.symbol.eq(sym)]
        n = sub.pivot(index="normal_group", columns="tumour_group", values="pairs").to_numpy()
        frac = n / n.sum(1, keepdims=True).clip(min=1)
        ax.imshow(frac, cmap=SEQ, vmin=0, vmax=1, aspect="equal")
        for r_ in range(3):
            for c_ in range(3):
                ax.text(c_, r_, f"{n[r_, c_]}", ha="center", va="center", fontsize=6,
                        color="white" if frac[r_, c_] > 0.55 else INK)
                rows.append({"panel": "ab"[i], "symbol": sym, "event": f"{VIEW_LABEL[sub.view.iat[0]]} {sub.region.iat[0]}",
                             "normal_group": GROUP_NAMES[r_], "tumour_group": GROUP_NAMES[c_], "pairs": int(n[r_, c_]),
                             "share_of_normal_group": frac[r_, c_]})
        heat_axes(ax, GROUP_NAMES, GROUP_NAMES, row_size=6, col_size=6)
        k = float(chk.loc[chk.symbol.eq(sym), "pair_kappa"].iat[0])
        ax.set_xlabel(f"own tumour's group\n{sym} {VIEW_LABEL[sub.view.iat[0]]}, κ = {k:.2f}", fontsize=6)
        if i == 0:
            ax.set_ylabel("normal's PSI group", fontsize=6)
        rows.append({"panel": "ab"[i], "symbol": sym, "pair_kappa": k})
        axes.append(ax)
    # c. pair agreement over events
    cx = fig.add_subplot(gs[2])
    sets = [("all", chk.pair_kappa, MUTED), ("A-13", chk.pair_kappa[chk.confirmed], INK),
            ("top", chk.pair_kappa[top], ACCENT)]
    for lab, v, col in sets:
        v = np.sort(v.dropna().to_numpy())
        y = np.arange(1, len(v) + 1) / len(v)
        cx.step(v, y, where="post", color=col, linewidth=1.0, label=f"{lab} ({len(v)})")
        rows += [{"panel": "c", "set": lab, "pair_kappa": float(a), "ecdf": float(b)} for a, b in zip(v, y)]
    cx.set_xlim(-0.25, 1.0)
    cx.set_ylim(0, 1.02)
    cx.set_xlabel("tumour–normal agreement (κ)", fontsize=6)
    cx.set_ylabel("share of events", fontsize=6)
    cx.set_xticks([0, 0.5, 1])
    cx.set_box_aspect(1)
    cx.legend(loc="lower right", fontsize=5.2, handlelength=1.0, handletextpad=0.4, borderaxespad=0.0,
              bbox_to_anchor=(1.02, 0.0))
    # d. high-allele frequency in LumA against Basal tumours
    dx = fig.add_subplot(gs[3])
    dx.plot([0, 1], [0, 1], color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    rows.append({"panel": "d", "set": "reference line (equal frequency)", "af_LumA": 0.0, "af_Basal": 0.0,
                 "af_LumA_end": 1.0, "af_Basal_end": 1.0})
    for lab, k, col, s in [("all", ~chk.confirmed & ~top, GRID, 3), ("A-13", chk.confirmed, INK, 5),
                           ("top", top, ACCENT, 7)]:
        sub = chk[k]
        dx.scatter(sub.af_LumA, sub.af_Basal, s=s, color=col, linewidth=0, zorder=2 if col != GRID else 1)
        rows += [{"panel": "d", "set": lab, "symbol": r.symbol, "view": r.view, "region": r.region,
                  "af_LumA": r.af_LumA, "af_Basal": r.af_Basal, "pair_kappa": r.pair_kappa}
                 for r in sub.itertuples(index=False)]
    for sym, xy in [("CAST", (0.97, 0.1)), ("TMEM101", (0.02, 0.8))]:
        r = chk[chk.symbol.eq(sym)].iloc[0]
        dx.scatter([r.af_LumA], [r.af_Basal], s=26, facecolor="none", edgecolor=INK, linewidth=0.6, zorder=3,
                   clip_on=False)
        dx.annotate(sym, (r.af_LumA, r.af_Basal), xytext=xy, textcoords="data", fontsize=5.6, color=INK,
                    ha="right" if xy[0] > 0.5 else "left", va="center",
                    arrowprops=dict(arrowstyle="-", color=INK2, linewidth=0.5, shrinkA=1, shrinkB=3))
    dx.set_xlim(0, 1)
    dx.set_ylim(0, 1)
    dx.set_box_aspect(1)
    dx.set_xticks([0, 0.5, 1])
    dx.set_yticks([0, 0.5, 1])
    dx.set_xlabel("high group, LumA", fontsize=6)
    dx.set_ylabel("high group, Basal", fontsize=6)
    for a_, lab, x in [(axes[0], "a", -0.5), (axes[1], "b", -0.5), (cx, "c", -0.42), (dx, "d", -0.42)]:
        a_.text(x, 1.06, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    save(fig, OUT, "U12_brca_germline_check", data=pd.DataFrame(rows),
         caption=CAPTIONS["U12"](chk, prof, pr, pc, int(pt[pt.symbol.eq("CAST")].pairs.sum()), top))


def _cap_u12(chk, prof, pr, pc, n_pairs, top):
    gap = (chk.af_LumA - chk.af_Basal).abs()
    k6 = chk.pair_kappa > 0.6
    verdicts = lambda d: ", ".join(f"{i} {d.loc[i, 'verdict']}" for i in ["L1", "L2", "L3", "L4", "L5", "L6", "A3"]
                                   if i in d.index)
    return f"""**Are the BRCA models' genotype-like splicing events inherited?**

D-57 (config/followup.yaml: brca_germline) called an event genotype-like when, without labels, its PSI fell into three
well-separated groups (low, mid, high) with mid-group normals. A-13 adds a check from tumour–normal pairs, and A-14
(exploratory) adds pair agreement. One tumour and one normal per patient ({n_pairs} pairs). An inherited event keeps its group in the tumour, apart from loss of heterozygosity (mid to low or high);
a homozygous normal cannot become the opposite homozygote.

a, b, each normal's group against its own tumour's group (counts; shade = share of the normal's row) for CAST SE, a
germline event, and TMEM101 A5SS, the genotype-like event with the largest importance in the BRCA models.
c, cumulative distribution of Cohen's κ between normal and tumour group, for all {len(chk)} checked events (one per
gene), the {int(chk.confirmed.sum())} that pass A-13 (opposite-homozygote flips in ≤ 5% of at least 20 homozygous
pairs) and the {int(top.sum())} that rank among the models' top events (legend: all, A-13, top). d, share of tumours in the high group, LumA
(x) against Basal (y), per event, coloured as in c (grey = the rest); dashed line = equal frequency. An inherited allele differs between subtypes only
through ancestry.

- Only {int(chk.confirmed.sum())} of {len(chk)} events pass A-13 (median flip rate {chk.flip_rate.median():.2f}), and
  {int(k6.sum())} have κ > 0.6. Of the {int(top.sum())} genotype-like events among the models' top events,
  {int((top & chk.confirmed).sum())} pass A-13, and their κ is at most {chk.pair_kappa[top].max():.2f}.
- Median |LumA − Basal| frequency: top events {gap[top].median():.2f}, A-13 confirmed
  {gap[chk.confirmed].median():.2f}, all {gap.median():.2f}. The models' genotype-like events behave as subtype programmes, not as inherited alleles.
- D-57, all events: {verdicts(pr)}. Germline-confirmed only (A-13): {verdicts(pc)}."""


def U13() -> None:
    m47 = _step("47_brca_her2.py")
    pr = pd.read_csv(BRCA_DIR / "her2" / "predictions.csv").set_index("id")
    ex = pd.read_csv(BRCA_DIR / "her2" / "exploratory.csv").set_index("what")
    ko = pd.read_csv(BRCA_DIR / "her2" / "knockout.csv")
    d = m47.oof_her2()
    d["s15"] = d.sample_id.str[:15]
    d = d.join(m47.erbb2_copy_number(), on="s15")
    d["clin"] = d.case_id.map(m47.clinical_her2())
    d = d[d.amplified.notna() & d.label.isin(["LumA", "LumB", "Her2"])].copy()
    d["amplified"] = d.amplified.astype(bool)
    d["grp"] = np.where(d.label.eq("Her2"), "Her2", "LumA/LumB")
    fig = plt.figure(figsize=(DOUBLE, 2.55))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 0.95, 0.85], wspace=1.05)
    rows, local = [], []
    rng = np.random.default_rng(0)
    # a. P(Her2) by PAM50 group and ERBB2 amplification, GEX and AS ensemble
    ax = fig.add_subplot(gs[0])
    groups = [("LumA/LumB", False), ("LumA/LumB", True), ("Her2", False), ("Her2", True)]
    ticks, labs = [], []
    for i, (g, a) in enumerate(groups):
        sub = d[d.grp.eq(g) & d.amplified.eq(a)]
        for k, (m, col) in enumerate([("p_gex", INK2), ("p_as", ACCENT)]):
            x = i * 1.0 + (k - 0.5) * 0.36
            st = _box(ax, x, sub[m].to_numpy(), col, width=0.3, rng=rng, s=1.6 if len(sub) > 100 else 3)
            rows.append({"panel": "a", "group": g, "erbb2": "amplified" if a else "not amplified",
                         "model": "GEX" if m == "p_gex" else "AS ensemble", **st})
        local.append(sub[["grp", "amplified", "p_gex", "p_as"]])
        ticks.append(i)
        labs.append(f"{'LumA/B' if g == 'LumA/LumB' else g}\n{'amp.' if a else 'not amp.'}\n({len(sub)})")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labs, fontsize=5.4)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, 3.6)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("out-of-fold P(Her2)", fontsize=6)
    ax.scatter([], [], marker="D", s=10, color=INK2, label="GEX")
    ax.scatter([], [], marker="D", s=10, color=ACCENT, label="AS ensemble")
    ax.legend(loc="upper left", fontsize=5.4, handletextpad=0.2, borderaxespad=0.1)
    # b. AUROC, GEX against AS, for each comparison
    bx = fig.add_subplot(gs[1])
    comp = [("H1, H2", "LumA/B: amp. vs not", pr.loc["H1"].effect + 0.5,
             pr.loc["H1"].effect + 0.5 - pr.loc["H2"].effect),
            ("expl.", "LumA/B: clin. HER2+ vs −", ex.loc["LumA/LumB: AUROC(p_gex, clinical HER2 positive against negative)", "value"],
             ex.loc["LumA/LumB: AUROC(p_as, clinical HER2 positive against negative)", "value"]),
            ("H3", "Her2: amp. vs not", None, None),
            ("H4", "Her2 not amp. vs LumA/B amp.", ex.loc[[i for i in ex.index if i.startswith("H4 counterpart")][0], "value"],
             pr.loc["H4"].effect + 0.5)]
    h3 = [float(s.split()[-1].rstrip(",")) for s in pr.loc["H3"].n.split(";")[2].replace("GEX", "").replace("AS", "").split(",")]
    comp[2] = ("H3", comp[2][1], h3[0], h3[1])
    for j, (pid, lab, g_, a_) in enumerate(comp):
        bx.plot([g_, a_], [j, j], color=GRID, linewidth=1.2, zorder=1)
        bx.scatter(g_, j, s=16, color=INK2, zorder=2)
        bx.scatter(a_, j, s=16, color=ACCENT, zorder=3)
        p_ = pr.loc[pid].p if pid in pr.index else np.nan
        rows.append({"panel": "b", "comparison": lab.replace("\n", " "), "test": pid, "auroc_gex": g_, "auroc_as": a_,
                     "p": p_, "q": pr.loc[pid].q if pid in pr.index else np.nan,
                     "verdict": pr.loc[pid].verdict if pid in pr.index else "exploratory"})
    bx.axvline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    bx.set_yticks(range(len(comp)))
    bx.set_yticklabels([f"{c[1]}\n({c[0]})" for c in comp], fontsize=5.4)
    bx.set_ylim(len(comp) - 0.5, -0.5)
    bx.set_xlim(0.3, 1.02)
    bx.set_xticks([0.4, 0.6, 0.8, 1.0])
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.set_xlabel("AUROC", fontsize=6)
    # c. knockout of the 17q12 amplicon inputs in the stored fold models
    cx = fig.add_subplot(gs[2])
    order = ["gex"] + [v for v in ["afe", "ale", "hit", "se", "a5ss", "a3ss", "mxe", "ri"]]
    for j, v in enumerate(order):
        r = ko[ko.view.eq(v)].iloc[0]
        col = INK2 if v == "gex" else ACCENT
        cx.plot([r.auroc_base, r.auroc_core], [j, j], color=GRID, linewidth=1.2, zorder=1)
        cx.scatter(r.auroc_base, j, s=14, color=col, zorder=3)
        cx.scatter(r.auroc_core, j, s=14, facecolors="white", edgecolors=col, linewidths=0.8, zorder=4)
        cx.scatter(r.auroc_random, j, marker="|", s=20, color=MUTED, linewidths=0.7, zorder=2)
        rows.append({"panel": "c", "view": VIEW_LABEL.get(v, v), "inputs_per_fold": r.amplicon_inputs_per_fold,
                     "auroc_intact": r.auroc_base, "auroc_amplicon_removed": r.auroc_core,
                     "auroc_random_removed": r.auroc_random})
    cx.set_yticks(range(len(order)))
    cx.set_yticklabels([VIEW_LABEL.get(v, v) for v in order], fontsize=5.6)
    cx.set_ylim(len(order) - 0.5, -0.5)
    cx.set_xlim(0.64, 0.82)
    cx.set_xticks([0.65, 0.7, 0.75, 0.8])
    cx.tick_params(axis="y", length=0)
    cx.spines["left"].set_visible(False)
    cx.set_xlabel("AUROC, amp. LumA/B", fontsize=6)
    cx.scatter([], [], s=14, color=INK2, label="intact")
    cx.scatter([], [], s=14, facecolors="white", edgecolors=INK2, linewidths=0.8, label="17q12 removed")
    cx.scatter([], [], marker="|", s=20, color=MUTED, label="random removed")
    cx.legend(loc="lower left", fontsize=5.2, handletextpad=0.1, borderaxespad=0.1, bbox_to_anchor=(0, 1.0), ncol=1)
    for a_, lab, x in [(ax, "a", -0.2), (bx, "b", -0.98), (cx, "c", -0.36)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    rows += [{"panel": "b", "comparison": "reference line (chance)", "auroc_gex": 0.5, "auroc_as": 0.5}]
    save(fig, OUT, "U13_brca_two_her2", data=pd.DataFrame(rows), caption=CAPTIONS["U13"](pr, ex, ko, d),
         local_data=(SCFG.base_data_dir / "figure_data", pd.concat(local, ignore_index=True)))


VIEW_LABEL = {"afe": "AFE", "ale": "ALE", "hit": "HIT", "se": "SE", "a5ss": "A5SS", "a3ss": "A3SS", "mxe": "MXE",
              "ri": "RI", "gex": "GEX"}


def U14() -> None:
    m48 = _step("48_brca_mutations.py")
    pr = pd.read_csv(BRCA_DIR / "mutations" / "predictions.csv").set_index("id")
    sp = pd.read_parquet(BRCA_LOCAL / "splice_site_mutations.parquet")
    from ctsm.config import load_config
    score, ev = m48.cryptic_a3ss(load_config())
    pam = pd.read_parquet(SCFG.data_dir / "tasks" / "subtype__BRCA_PAM50.parquet")
    samples = set(pam.sample_id.str[:15])
    calls = m48.cohort_calls(samples)
    sf = calls[calls.Hugo_Symbol.eq("SF3B1") & calls.Variant_Classification.eq("Missense_Mutation")].copy()
    sf["codon"] = sf.HGVSp_Short.str.extract(r"p\.[A-Z](\d+)")[0].astype(float)
    hot = set(sf.loc[sf.codon.isin(m48.HOTSPOT), "s15"])
    score = score[score.index.isin(samples)]
    fig = plt.figure(figsize=(DOUBLE * 0.72, 2.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0], wspace=0.45)
    rows = []
    # a. splice-site mutations: percentile of the hit exon's PSI in the mutated tumour
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(0, 1, 11)
    bottom = np.zeros(10)
    for side, col in [("donor", ACCENT_DARK), ("acceptor", ACCENT)]:
        h, _ = np.histogram(sp.loc[sp.side.eq(side), "percentile"], bins=bins)
        ax.bar(bins[:-1] + 0.05, h, width=0.092, bottom=bottom, color=col, linewidth=0, label=f"{side} (+1/+2)" if side == "donor" else "acceptor (−1/−2)")
        for b0, n in zip(bins[:-1], h):
            rows.append({"panel": "a", "side": side, "bin_low": b0, "bin_high": b0 + 0.1, "mutations": int(n)})
        bottom += h
    exp = len(sp) / 10
    ax.axhline(exp, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    rows.append({"panel": "a", "side": "reference line (no effect, uniform)", "mutations": exp})
    cd = sp[sp.symbol.eq("CDH1")]
    ax.scatter(cd.percentile, np.full(len(cd), bottom.max() * 1.08), marker="v", s=14, color=INK, zorder=4, label="CDH1")
    rows += [{"panel": "a", "side": "CDH1 mutation", "percentile_marker": float(v)} for v in cd.percentile]
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xlabel("PSI percentile of the hit exon in the mutated tumour", fontsize=6)
    ax.set_ylabel("splice-site mutations", fontsize=6)
    ax.legend(loc="upper right", fontsize=5.4, handletextpad=0.3, borderaxespad=0.1)
    # b. SF3B1 hotspot tumours: cryptic A3SS score
    bx = fig.add_subplot(gs[1])
    rng = np.random.default_rng(0)
    y = score.index.isin(hot)
    local = []
    for i, (lab, v, col) in enumerate([("rest", score[~y].to_numpy(), MUTED), ("SF3B1\nhotspot", score[y].to_numpy(), ACCENT)]):
        st = _box(bx, i, v, col, width=0.5, rng=rng, s=1.6 if len(v) > 100 else 5)
        rows.append({"panel": "b", "group": lab.replace("\n", " "), **st})
        local.append(pd.DataFrame({"group": lab.replace("\n", " "), "cryptic_a3ss_score": v}))
    bx.set_xticks([0, 1])
    bx.set_xticklabels([f"rest\n({int((~y).sum())})", f"SF3B1 hotspot\n({int(y.sum())})"], fontsize=5.6)
    bx.tick_params(axis="x", length=0)
    bx.set_xlim(-0.6, 1.6)
    bx.set_ylabel("cryptic A3SS score (percentile)", fontsize=6)
    top = score.max()
    bx.plot([0, 0, 1, 1], [top + 0.01, top + 0.02, top + 0.02, top + 0.01], color=INK2, linewidth=0.5)
    bx.text(0.5, top + 0.025, _stars(pr.loc["C3"].p), ha="center", va="bottom", fontsize=6)
    rows.append({"panel": "b", "group": "hotspot vs rest (C3)", "auroc": pr.loc["C3"].effect + 0.5, "p": pr.loc["C3"].p,
                 "q": pr.loc["C3"].q, "q_family": "C (brca_mutations, 3 tests)"})
    for a_, lab, x in [(ax, "a", -0.14), (bx, "b", -0.42)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    save(fig, OUT, "U14_brca_mutations", data=pd.DataFrame(rows), caption=CAPTIONS["U14"](pr, sp, ev, y),
         local_data=(SCFG.base_data_dir / "figure_data", pd.concat(local, ignore_index=True)))


def _cap_u13(pr, ex, ko, d):
    lum = d[d.grp.eq("LumA/LumB")]
    her = d[d.grp.eq("Her2")]
    gex_ko = ko[ko.view.eq("gex")].iloc[0]
    as_ko = ko[ko.view.ne("gex")]
    return f"""**Which HER2 did each modality learn: ERBB2 amplification or the HER2-enriched programme?**

Pre-specified before any probability was compared with copy number (config/followup.yaml: brca_her2, D-58). P(Her2)
is the stored out-of-fold probability of the BRCA_PAM50 models: GEX, and the AS ensemble (mean of the 8 splicing
views). ERBB2 amplified = ABSOLUTE total copy number at ERBB2 ≥ 2 × ploidy ({int(lum.amplified.sum())} of {len(lum)}
LumA/LumB and {int(her.amplified.sum())} of {len(her)} PAM50 Her2 tumours). ABSOLUTE's copy number at ERBB2 does not
exceed 7 in any tumour, so the rule catches the highest-level amplifications only:
{ex.loc['amplified among clinical HER2 positive', 'value'] * 100:.0f}% of clinically HER2-positive tumours meet it.
Clinical HER2 (IHC or ISH) is the pre-specified second reading.

a, P(Her2) by PAM50 group and amplification; box = interquartile range, diamond = median, whiskers = 1.5 × IQR, dots =
tumours (n under each group). b, AUROC of P(Her2) for each comparison, GEX (grey) and AS ensemble (blue).
c, in-silico knockout (exploratory, A-12): in each stored fold model the inputs from chr17:39.6–39.8 Mb (GRCh38;
PPP1R1B to GRB7) are set to their training mean, and the LumA/LumB AUROC is recomputed (hollow); tick = 20 random
input sets of the same size.

- H1, GEX P(Her2) is higher in amplified luminal tumours: AUROC {pr.loc['H1'].effect + 0.5:.2f}, q =
  {fp(pr.loc['H1'].q)}, supported.
- H2, splicing would see amplification less well than GEX: AUROC difference {signed(pr.loc['H2'].effect, 2)}
  (q = {fp(pr.loc['H2'].q)}), **against**: the AS ensemble separates amplified luminal tumours better. Without the
  ABSOLUTE ceiling (D-60, U15) this lead is not resolved, and absent on the newly amplified tumours.
- H3, within PAM50 Her2, neither model separates amplified from non-amplified tumours (difference
  {signed(pr.loc['H3'].effect, 3)}, not resolved).
- H4, AS P(Her2) is higher in non-amplified Her2 than in amplified luminal tumours: AUROC
  {pr.loc['H4'].effect + 0.5:.2f}, supported; GEX {ex.loc[[i for i in ex.index if i.startswith('H4 counterpart')][0], 'value']:.2f}.
- c, removing the amplicon inputs lowers the GEX AUROC from {gex_ko.auroc_base:.3f} to {gex_ko.auroc_core:.3f}; in the
  splicing views the change is at most {(as_ko.auroc_base - as_ko.auroc_core).max():.3f}. The splicing models do not
  reach these tumours through the amplicon; they give them part of the HER2-enriched splicing programme.

q over family H (4 tests). Per-tumour probabilities are kept locally (D-12)."""


def _cap_u14(pr, sp, ev, y):
    side = sp.groupby("side").percentile.median()
    return f"""**Do splice-site mutations and SF3B1 hotspots leave the expected marks in BRCA splicing?**

Pre-specified (config/followup.yaml: brca_mutations, D-59). MC3 public calls of the cohort's tumours (FILTER flags
within PASS, wga, native_wga_mix). a, every splice-site SNV at +1/+2 or −1/−2 whose hit exon (by the MC3 transcript's
intron number, located in GENCODE v45) is a cassette exon of the main SE matrix: the exon's PSI percentile in the
mutated tumour among the other cohort tumours (ties count half, so an exon near 100% inclusion in almost every
tumour ranks at 0.5), stacked by side; dashed line = the count per bin with no effect;
triangles = the {int(sp.symbol.eq('CDH1').sum())} CDH1 mutations. b, SF3B1 hotspot tumours (codons 622, 625, 626,
662, 663, 666, 700, 741, 742) against the rest on a cryptic A3SS score: mean within-cohort percentile of the long
form of {len(ev)} A3SS events whose acceptors are 10–30 nt apart and whose upstream acceptor is rarely used in normals
(median PSI < 0.2). Box = interquartile range, diamond = median, whiskers = 1.5 × IQR, dots = tumours.

- C1, the hit exon is skipped: {pr.loc['C1'].n}; q = {fp(pr.loc['C1'].q)}, supported. Median percentile: donor
  {side.get('donor', np.nan):.2f}, acceptor {side.get('acceptor', np.nan):.2f}.
- C2, CDH1 alone: {pr.loc['C2'].n}; with five mutations the sign test cannot pass (p = {fp(pr.loc['C2'].p)}), not
  resolved.
- C3, cryptic acceptors in SF3B1-hotspot tumours: AUROC {pr.loc['C3'].effect + 0.5:.2f} ({int(y.sum())} against
  {int((~y).sum())}), q = {fp(pr.loc['C3'].q)}, supported. All {int(y.sum())} hotspot tumours are LumA or LumB.

Stars: * p < 0.05, ** p < 0.01, *** p < 0.001. q over family C (3 tests). Per-tumour values are kept locally (D-12)."""


def U15() -> None:
    hd = BRCA_DIR / "her2_uncapped"
    pr = pd.read_csv(hd / "predictions.csv").set_index("id")
    ex = pd.read_csv(hd / "exploratory.csv").set_index("what")
    gt = pd.read_csv(hd / "luminal_amplification_groups.csv").set_index("group")
    ko = pd.read_csv(hd / "knockout.csv")
    d = pd.read_parquet(BRCA_LOCAL / "her2_uncapped.parquet")
    d = d[d.label.isin(["LumA", "LumB"]) & d.amp.notna()].copy()
    d["amp"] = d.amp.astype(bool)
    d["new"] = d.amp & ~d.amp_absolute.astype("boolean").fillna(False).astype(bool)
    fig = plt.figure(figsize=(DOUBLE * 0.74, 2.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.5)
    rows = []
    # a. AUROC of P(Her2), amplified luminal tumours against the rest, by amplification call
    ax = fig.add_subplot(gs[0])
    auc_all = {"GEX": pr.loc["HU1"].effect + 0.5,
               "AS ensemble": float(pr.loc["HU2"].n.split("AS ")[1].split(",")[0])}
    groups = [("D-58\nnear-diploid", gt.loc["D-58 amplified"], int(gt.loc["D-58 amplified", "n"]), "D-58 amplified (ABSOLUTE)"),
              ("all\nuncapped", auc_all, int(d.amp.sum()), "all amplified (D-60)"),
              ("newly\namplified", gt.loc["newly amplified"], int(gt.loc["newly amplified", "n"]), "newly amplified only")]
    for i, (lab, src, n, name) in enumerate(groups):
        vals = [src["auroc_gex"] if "auroc_gex" in src else src["GEX"], src["auroc_as"] if "auroc_as" in src else src["AS ensemble"]]
        ax.plot([i - 0.13, i + 0.13], vals, color=GRID, linewidth=0.8, zorder=1)
        for k, (m, col, v) in enumerate(zip(["GEX", "AS ensemble"], [INK2, ACCENT], vals)):
            ax.scatter(i + (k - 0.5) * 0.26, v, s=18, color=col, zorder=2, label=m if i == 0 else None)
            rows.append({"panel": "a", "group": name, "amplified_n": n, "model": m, "auroc": v})
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"{g[0]}\n({g[2]})" for g in groups], fontsize=5.6)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(0.6, 0.9)
    ax.set_yticks([0.6, 0.7, 0.8, 0.9])
    ax.set_ylabel("AUROC, amplified luminal vs rest", fontsize=6)
    ax.legend(loc="lower right", fontsize=5.4, handletextpad=0.2, borderaxespad=0.1, scatterpoints=1)
    # b. amplified luminal tumours: P(Her2) against the amplification level
    bx = fig.add_subplot(gs[1])
    a = d[d.amp]
    local = []
    for m, col, lab in [("p_gex", INK2, "GEX"), ("p_as", ACCENT, "AS ensemble")]:
        for new, face in [(False, col), (True, "white")]:
            s = a[a.new.eq(new)]
            bx.scatter(s.rel, s[m], s=12, facecolor=face, edgecolor=col, linewidth=0.6, zorder=2)
        r = spearmanr(a[m], a.rel).statistic
        rows.append({"panel": "b", "model": lab, "spearman_rho": r, "n": len(a)})
        local.append(a[["rel", m, "new"]].rename(columns={m: "p_her2"}).assign(model=lab))
    bx.set_yscale("log")
    bx.set_xlabel("ERBB2 copy number / ploidy (SNP6, uncapped)", fontsize=6)
    bx.set_ylabel("P(Her2), amplified luminal", fontsize=6)
    from matplotlib.lines import Line2D
    bx.legend(handles=[Line2D([], [], marker="o", linestyle="", markersize=3.5, markerfacecolor=INK2, markeredgecolor=INK2,
                              label="D-58 amplified"),
                       Line2D([], [], marker="o", linestyle="", markersize=3.5, markerfacecolor="white", markeredgecolor=INK2,
                              label="newly amplified")],
              loc="lower right", fontsize=5.4, handletextpad=0.2, borderaxespad=0.1)
    for a_, lab, x in [(ax, "a", -0.32), (bx, "b", -0.25)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    for r in pr.reset_index().itertuples(index=False):
        rows.append({"panel": "tests", "group": r.id, "auroc": r.effect, "p": r.p, "q": r.q, "q_family": "HU (7 tests)",
                     "model": r.verdict})
    save(fig, OUT, "U15_brca_her2_uncapped", data=pd.DataFrame(rows), caption=CAPTIONS["U15"](pr, ex, gt, ko, d),
         local_data=(SCFG.base_data_dir / "figure_data", pd.concat(local, ignore_index=True)))


def _cap_u15(pr, ex, gt, ko, d):
    g = gt.loc["D-58 amplified"]
    nw = gt.loc["newly amplified"]
    gko = ko[ko.view.eq("gex")].iloc[0]
    as_drop = (ko.loc[ko.view.ne("gex"), "auroc_base"] - ko.loc[ko.view.ne("gex"), "auroc_core"]).abs().max()
    hu2 = pr.loc["HU2"]
    as2 = float(hu2.n.split("AS ")[1].split(",")[0])
    rho_as = float(pr.loc["HU7"].n.split("AS rho ")[1])
    return f"""**Does splicing's lead on ERBB2-amplified luminal tumours survive without the ABSOLUTE ceiling?**

Pre-specified before the copy-number file was downloaded (config/followup.yaml: brca_her2_uncapped, D-60). D-58 took
amplification from ABSOLUTE, which never reports more than 7 copies at ERBB2. Here, the PanCanAtlas SNP6 copy ratio
at ERBB2 is made absolute with ABSOLUTE purity and ploidy, without a ceiling; amplified = copy number ≥ 2 × ploidy,
D-58's rule. Out-of-fold P(Her2) of the BRCA_PAM50 models; LumA and LumB tumours only.

a, AUROC of P(Her2), GEX (grey) and AS ensemble (blue), for D-58's amplified tumours (median ploidy
{g.median_ploidy:.2f}), all tumours amplified by the uncapped call, and the newly amplified ones only (median ploidy
{nw.median_ploidy:.2f}). Each group is compared with the luminal tumours neither call makes amplified; the middle group
is compared with all non-amplified tumours. b, amplified luminal tumours: P(Her2) against copy number / ploidy;
filled = amplified by D-58 too, hollow = newly amplified.

- HU1, GEX P(Her2) higher in amplified tumours: AUROC {pr.loc['HU1'].effect + 0.5:.2f}, q = {fp(pr.loc['HU1'].q)},
  supported.
- HU2, AS ensemble ahead of GEX: {as2:.2f} against {pr.loc['HU1'].effect + 0.5:.2f} (difference
  {signed(hu2.effect, 3)}, q = {fp(hu2.q)}), not resolved. D-58 had found the gap with q = 0.016.
- HU3, the independent part (newly amplified): AS {nw.auroc_as:.2f}, GEX {nw.auroc_gex:.2f}, not resolved. The
  ceiling could not call any tumour with ploidy above 3.5, where 7 copies are below 2 × ploidy.
- HU4, within PAM50 Her2, GEX ahead of AS for amplified against not: {signed(pr.loc['HU4'].effect, 3)}, not resolved.
- HU5, AS P(Her2) higher in non-amplified Her2 than in amplified luminal tumours: AUROC
  {pr.loc['HU5'].effect + 0.5:.2f}, supported.
- HU6, GEX P(Her2) rises with the amplification level: ρ = {pr.loc['HU6'].effect:.2f}, supported. HU7, splicing
  follows it less: AS ρ = {rho_as:.2f}, not resolved.
- Knockout with this amplified group (exploratory): the 17q12 inputs set to their mean lower GEX from
  {gko.auroc_base:.3f} to {gko.auroc_core:.3f}; no splicing view moves by more than {as_drop:.3f}.

q over family HU (7 tests). Per-tumour values are kept locally (D-12)."""

# ============================================================ MSI follow-up (D-61..D-64)
FU_DIR = MSI_DIR / "followup"


def U16() -> None:
    pr = pd.read_csv(FU_DIR / "load" / "predictions.csv").set_index("id")
    gr = pd.read_csv(FU_DIR / "load" / "by_cohort_group.csv")
    d = pd.read_parquet(MSI_LOCAL / "deletion_load.parquet")
    au = pd.read_csv(MSI_DIR / "score_aurocs.csv")
    au = au[au.score.eq("acc")].set_index("cohort").auroc_msi_vs_rest
    coh = ["STAD_TCGA", "COAD_TCGA", "UCEC_TCGA"]
    test_of = {"STAD_TCGA": "DL1", "COAD_TCGA": "DL2", "UCEC_TCGA": "DL3"}
    fig = plt.figure(figsize=(DOUBLE * 0.92, 2.05))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.95], wspace=0.42)
    ticks = [0, 10, 100, 1000]
    rows, local, axes = [], [], []
    ylo, yhi = d.score_acc.min(), d.score_acc.max()
    pad = 0.06 * (yhi - ylo)
    for i, c in enumerate(coh):
        ax = fig.add_subplot(gs[i])
        axes.append(ax)
        col, mk = COH_STYLE[c]
        for g, face, s_, z in [("rest", MUTED, 3, 1), ("MSI", col, 7, 2)]:
            v = d[d.cohort.eq(c) & d.group.eq(g)]
            ax.scatter(v.load, v.score_acc, s=s_, marker=mk, color=face, alpha=0.55 if g == "rest" else 0.85,
                       linewidths=0, zorder=z)
            local.append(v[["cohort", "group", "n_hp", "load", "score_acc", "purity"]])
        t = pr.loc[test_of[c]]
        ax.text(0.98, 0.97, f"MSI ρ = {signed(t.effect, 2)}\np = {fp(t.p)}", transform=ax.transAxes, ha="right",
                va="top", fontsize=5.4, color=INK)
        ax.set_xticks(np.log1p(ticks))
        ax.set_xticklabels([str(k) for k in ticks], fontsize=5.6)
        ax.set_xlim(-0.25, np.log1p(1500))
        ax.set_ylim(ylo - pad, yhi + pad)
        ax.set_title(_short(c), fontsize=6.5, pad=2)
        ax.set_xlabel("homopolymer deletions", fontsize=6)
        if i == 0:
            ax.set_ylabel("repeat-exon score (D-52)", fontsize=6)
        else:
            ax.set_yticklabels([])
        rows.append({"panel": c, "set": f"{test_of[c]}: MSI tumours, partial Spearman | purity", "n": t.n,
                     "value": t.effect, "p": t.p, "q": t.q, "q_family": "DL (msi_deletion_load, 6 tests)",
                     "verdict": t.verdict})
        for g in ["MSI", "rest"]:
            r = gr[gr.cohort.eq(c) & gr.group.eq(g)].iloc[0]
            rows.append({"panel": c, "set": f"{g} tumours plotted", "n": int(r.n), "value": r.median_n_hp,
                         "statistic": "median homopolymer deletions"})
    # last panel: MSI tumours' deletion load by cohort, with the cohort's D-52 AUROC
    bx = fig.add_subplot(gs[3])
    rng = np.random.default_rng(0)
    for i, c in enumerate(coh):
        v = d[d.cohort.eq(c) & d.group.eq("MSI")].n_hp.to_numpy()
        b = _box(bx, i, np.log1p(v), COH_STYLE[c][0], width=0.5, rng=rng, s=3)
        rows.append({"panel": "load by cohort", "set": f"{c} MSI", "n": b["n"], "value": float(np.expm1(b["median"])),
                     "statistic": "median homopolymer deletions (box on log(1 + count))",
                     "q1": float(np.expm1(b["q1"])), "q3": float(np.expm1(b["q3"]))})
        rows.append({"panel": "load by cohort", "set": f"{c} D-52 score AUROC, MSI vs rest", "value": float(au[c]),
                     "statistic": "AUROC (score lower in MSI)"})
    bx.set_xticks(range(3))
    bx.set_xticklabels([f"{_short(c)}\n{au[c]:.2f}" for c in coh], fontsize=5.6)
    bx.tick_params(axis="x", length=0)
    bx.set_yticks(np.log1p([10, 30, 100, 300]))
    bx.set_yticklabels(["10", "30", "100", "300"], fontsize=5.6)
    bx.set_xlim(-0.6, 2.6)
    bx.set_ylabel("homopolymer deletions, MSI", fontsize=6)
    for a_, lab in [(axes[0], "a"), (bx, "b")]:
        a_.text(-0.34 if a_ is axes[0] else -0.42, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    data = pd.DataFrame(rows)
    save(fig, OUT, "U16_msi_deletion_load", data=data, caption=CAPTIONS["U16"](pr, gr, au),
         local_data=(SCFG.base_data_dir / "figure_data", pd.concat(local, ignore_index=True)))


def _cap_u16(pr, gr, au):
    m = gr[gr.group.eq("MSI")].set_index("cohort")
    r = gr[gr.group.eq("rest")].set_index("cohort")
    dl = lambda k: f"ρ = {signed(pr.loc[k].effect, 2)}, q = {fp(pr.loc[k].q)}, {pr.loc[k].verdict}"
    return f"""**Does more mismatch-repair damage mean more skipping of repeat exons?**

Pre-specified (config/followup.yaml: msi_deletion_load, D-61). Homopolymer deletion = an MC3 deletion (FILTER PASS, wga
or native_wga_mix) of A or T bases inside a run of ≥ 5 of that base, counted per tumour in the exome. Score = the D-52
per-tumour repeat-exon score: mean inclusion percentile of exons whose acceptor carries a T-run ≥ 10, minus that of
matched exons without one (lower = repeat exons less included). Partial Spearman of score and log(1 + count), adjusted
for ABSOLUTE purity, within MSI tumours.

a, one panel per cancer: score against homopolymer deletions (x on a log(1 + count) scale); MSI tumours in colour,
the rest grey. The rest carry a median of {int(r.median_n_hp.max())} homopolymer deletions in every cancer (upper
quartile ≤ {int(r.q3_n_hp.max())}); MSI tumours a median of {int(m.loc['COAD_TCGA'].median_n_hp)} (COAD) to
{int(m.loc['UCEC_TCGA'].median_n_hp)} (UCEC). b, MSI tumours' homopolymer deletions by cancer (box = interquartile
range, diamond = median, whiskers 1.5 × IQR, points = tumours), with each cancer's D-52 score AUROC for MSI against the
rest (second line of each label).

- DL1–DL3, more deletions, lower score: STAD {dl('DL1')}; COAD {dl('DL2')}; UCEC {dl('DL3')}.
- DL4–DL5, a smaller share of multi-base deletions in UCEC (explaining its weaker shift): not resolved; the median
  share is {m.loc['UCEC_TCGA'].median_long_share:.3f} in UCEC, {m.loc['STAD_TCGA'].median_long_share:.3f} in STAD
  and {m.loc['COAD_TCGA'].median_long_share:.3f} in COAD. Over 98% of these deletions remove one base.
- DL6, multi-base deletions beyond single-base ones: {dl('DL6')}.

Dose tracks the score within STAD only. Between cancers it runs the other way: UCEC MSI tumours carry the most deletions
and the D-52 score does not separate them from the rest (AUROC {au['UCEC_TCGA']:.2f}). The exome counts coding runs,
mostly shorter than 10; the intronic tracts themselves were not measured. q over family DL (6 tests); per-tumour values
are kept locally (D-12)."""


def U17() -> None:
    pr = pd.read_csv(FU_DIR / "branch" / "predictions.csv").set_index("id")
    ex = pd.read_csv(FU_DIR / "branch" / "exploratory.csv")
    g = pd.read_parquet(MSI_LOCAL / "run_geometry.parquet").drop_duplicates(["chrom", "strand", "es", "ee"])
    coh = ["STAD_TCGA", "COAD_TCGA", "UCEC_TCGA"]
    fig = plt.figure(figsize=(DOUBLE * 0.86, 2.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.45, 0.95], wspace=0.55)
    rows = []
    # a. where the branch point lies: d_bp by run class
    ax = fig.add_subplot(gs[0])
    bins = np.arange(10, 102, 4)
    for cls, k, col in [("10–19", g.acc_T < 20, INK2), ("≥20", g.acc_T >= 20, ACCENT)]:
        v = g.loc[k, "d_bp"].dropna()
        h, e = np.histogram(v, bins=bins)
        ax.step(e[:-1], h / h.sum(), where="post", color=col, linewidth=0.8, label=f"{cls} nt")
        rows += [{"panel": "a", "set": f"T-run {cls}", "bin_start": int(a), "share": float(b / h.sum()), "n": len(v)}
                 for a, b in zip(e[:-1], h)]
    ax.axvspan(19, 23, color=GRID, zorder=0, linewidth=0)
    rows.append({"panel": "a", "set": "reference band: 19-23 nt (Chua and Reed 2001)", "bin_start": 19})
    ax.set_xlabel("branch point to 3′ss (nt)", fontsize=6)
    ax.set_ylabel("share of R exons", fontsize=6)
    ax.set_xlim(10, 100)
    ax.legend(loc="upper right", fontsize=5.2, handlelength=1.0, borderaxespad=0.1, bbox_to_anchor=(1.02, 0.72),
              title="T-run", title_fontsize=5.2)
    # b. share more included in MSI by branch-point tertile, within run classes
    bx = fig.add_subplot(gs[1])
    t = ex[ex.what.eq("share moving up by d_bp tertile")]
    for j, c in enumerate(coh):
        col, mk = COH_STYLE[c]
        for i, cls in enumerate(["10-19", ">=20"]):
            v = t[t.cohort.eq(c) & t.run_class.eq(cls)].set_index("tertile").loc[["near", "middle", "far"]]
            xs = i * 3.6 + np.arange(3) + (j - 1) * 0.16
            for x, r in zip(xs, v.itertuples()):
                n, ph, z = r.n, r.share_up, 1.96
                den = 1 + z * z / n
                mid = (ph + z * z / (2 * n)) / den
                half = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
                bx.vlines(x, mid - half, mid + half, color=col, linewidth=0.6, alpha=0.8)
                rows.append({"panel": "b", "cohort": c, "set": f"T-run {cls}, d_bp {r.Index} (median {r.median_d_bp:g} nt)",
                             "n": int(n), "share": ph, "ci_low": mid - half, "ci_high": mid + half})
            bx.plot(xs, v.share_up, color=col, linewidth=0.8, marker=mk, ms=3, mec="white", mew=0.3,
                    label=_short(c) if i == 0 else None)
    bx.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    rows.append({"panel": "b", "set": "reference line", "share": 0.5})
    bx.set_xticks([0, 1, 2, 3.6, 4.6, 5.6])
    bx.set_xticklabels(["near", "mid", "far"] * 2, fontsize=5.6)
    bx.tick_params(axis="x", length=0)
    for x, lab in [(1, "T-run 10–19"), (4.6, "T-run ≥ 20")]:
        bx.text(x, -0.2, lab, ha="center", va="top", fontsize=6, transform=bx.get_xaxis_transform())
    bx.set_ylabel("share more included in MSI", fontsize=6)
    bx.set_ylim(0.2, 0.85)
    bx.legend(loc="upper left", fontsize=5.4, handlelength=1.4, borderaxespad=0.2)
    # c. adjusted Spearman of the MSI shift with d_bp (tests) and with d_run (exploratory)
    cx = fig.add_subplot(gs[2])
    sets = [("BP1", "STAD_TCGA"), ("BP2", "COAD_TCGA"), ("BP3", "UCEC_TCGA"), ("BP4", "mean of >= 2 cohorts")]
    dr = ex[ex.what.eq("d_run in place of d_bp")].set_index("cohort")
    for i, (pid, c) in enumerate(sets):
        r = pr.loc[pid]
        col = COH_STYLE[c][0] if c in COH_STYLE else INK
        cx.scatter(r.effect, i - 0.12, marker="o", s=14, color=col, zorder=3)
        cx.scatter(dr.loc[c, "rho"], i + 0.12, marker="o", s=14, facecolor="white", edgecolor=col, linewidth=0.7, zorder=3)
        cx.text(0.255, i - 0.12, f"p = {fp(r.p)}", va="center", fontsize=5.2, color=INK)
        rows.append({"panel": "c", "cohort": c, "set": f"{pid}: d_bp", "n": r.n, "value": r.effect, "p": r.p, "q": r.q,
                     "q_family": "BP (msi_branchpoint, 4 tests)", "verdict": r.verdict})
        rows.append({"panel": "c", "cohort": c, "set": "d_run (exploratory, 200 permutations)", "n": int(dr.loc[c, "n"]),
                     "value": dr.loc[c, "rho"], "p": dr.loc[c, "p"]})
    cx.axvline(0, color=INK2, linewidth=0.6)
    rows.append({"panel": "c", "set": "reference line", "value": 0.0})
    cx.set_yticks(range(4))
    cx.set_yticklabels(["STAD", "COAD", "UCEC", "≥ 2 cancers"], fontsize=5.8)
    cx.set_ylim(4.5, -0.6)
    cx.set_xlim(-0.05, 0.25)
    cx.set_xticks([0, 0.1, 0.2])
    cx.tick_params(axis="y", length=0)
    cx.spines["left"].set_visible(False)
    cx.set_xlabel("adjusted Spearman ρ", fontsize=6)
    from matplotlib.lines import Line2D
    cx.legend(handles=[Line2D([], [], marker="o", linestyle="", markersize=3.5, color=INK, label="branch point"),
                       Line2D([], [], marker="o", linestyle="", markersize=3.5, markerfacecolor="white",
                              markeredgecolor=INK, label="run's far end")],
              loc="lower right", fontsize=5.2, handletextpad=0.2, borderaxespad=0.1, bbox_to_anchor=(1.0, -0.02))
    for a_, lab, x in [(ax, "a", -0.4), (bx, "b", -0.25), (cx, "c", -0.3)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    data = pd.DataFrame(rows)
    save(fig, OUT, "U17_msi_branchpoint", data=data, caption=CAPTIONS["U17"](pr, ex, g))


def _cap_u17(pr, ex, g):
    t = ex[ex.what.eq("share moving up by d_bp tertile") & ex.cohort.eq("mean of >= 2 cohorts")].set_index(
        ["run_class", "tertile"])
    bp = lambda k: f"ρ = {signed(pr.loc[k].effect, 3)}"
    return f"""**Does the branch point's distance from the 3′ splice site set the direction of the MSI shift?**

Pre-specified (config/followup.yaml: msi_branchpoint, D-63). Exons whose acceptor carries a T-run ≥ 10 (R exons, D-52),
read from the GRCh38 primary assembly. Branch point = the A nearest the 3′ splice site that lies upstream of the run,
within 100 nt, and matches yUnAy; found for {g.d_bp.notna().mean() * 100:.0f}% of {len(g):,} exons. MSI shift =
D-52's background-adjusted difference in inclusion percentile, MSI minus the rest.

a, distance from the branch point to the 3′ss, by run length ({int(g[g.acc_T < 20].d_bp.notna().sum()):,} exons with
runs of 10–19 nt, {int(g[g.acc_T >= 20].d_bp.notna().sum())} with ≥ 20 nt); grey band = 19–23 nt, the optimum for a 3′ss AG
(Chua and Reed 2001). b, share of R exons more included in MSI, by tertile of that distance within run classes,
per cancer; bars = Wilson 95% interval. c, the tests: Spearman of the MSI shift with the branch-point distance after
both are rank-residualised on run-length bins (10–11, 12–13, 14–15, 16–19, 20–24, ≥ 25) and mean inclusion; p by
permuting the distance within run-length bins (1,000; the floor is 0.001). Hollow = the same with the run's far end in
place of the branch point (exploratory).

- BP1–BP4, an exon whose branch point lies further from the 3′ss is more included in MSI, beyond its run length:
  STAD {bp('BP1')}, COAD {bp('BP2')}, UCEC {bp('BP3')}, mean over ≥ 2 cancers {bp('BP4')}; all q = {fp(pr.q.max())},
  supported.
- Mean over ≥ 2 cancers, runs of 10–19: {t.loc[('10-19', 'near')].share_up:.2f} of exons with the nearest branch
  points move up, {t.loc[('10-19', 'far')].share_up:.2f} with the farthest. Runs ≥ 20:
  {t.loc[('>=20', 'near')].share_up:.2f} against {t.loc[('>=20', 'far')].share_up:.2f}. The long-run reversal of
  D-53 sits in the exons whose branch point is far.

The run's far end carries nearly the same information as the branch point, because the branch point must lie beyond
it; these tests cannot separate the two. The branch point is a motif call, not a measured position. q over family BP
(4 tests)."""


# ============================================================ MSI sequence code (D-65..D-67)
CODE_DIR = MSI_DIR / "code"


def U18() -> None:
    pr = pd.read_csv(CODE_DIR / "pancancer" / "predictions.csv").set_index("id")
    per = pd.read_csv(CODE_DIR / "pancancer" / "by_cancer.csv")
    exo = pd.read_csv(CODE_DIR / "pancancer" / "exploratory.csv").set_index("what")
    el = pd.read_csv(CODE_DIR / "pancancer" / "exon_shifts_R.csv")
    sc = pd.read_parquet(MSI_LOCAL / "pancancer_scores.parquet")
    step = _step("51_msi_code.py")
    mm = sc.mmr.to_numpy()
    fig = plt.figure(figsize=(DOUBLE, 2.45))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.62)
    rows, local = [], []
    # a. AUROC of score A and score G, per cancer and pooled
    ax = fig.add_subplot(gs[0])
    big = per[per.mmr_deficient >= 3].sort_values("mmr_deficient", ascending=False)
    k_rel = ~sc.project.isin(step.RELATIVES).to_numpy()
    hy = sc.hyper.to_numpy()
    sel = mm | hy
    pooled = [("all 30 cancers", step.auc_lower(sc.pct_A.to_numpy(), mm), step.auc_lower(sc.pct_G.to_numpy(), mm),
               int(mm.sum()), "PC1, PC2"),
              ("without READ, ESCA, UCS", step.auc_lower(sc.pct_A.to_numpy()[k_rel], mm[k_rel]),
               step.auc_lower(sc.pct_G.to_numpy()[k_rel], mm[k_rel]), int((mm & k_rel).sum()), "PC4"),
              ("vs hypermutated proficient", step.auc_lower(sc.pct_A.to_numpy()[sel], mm[sel]),
               step.auc_lower(sc.pct_G.to_numpy()[sel], mm[sel]), int(mm.sum()), "PC5")]
    labels = []
    for i, (lab, a_, g_, n, pid) in enumerate(pooled):
        ax.plot([a_, g_], [i, i], color=GRID, linewidth=1.2, zorder=1)
        ax.scatter(a_, i, s=16, facecolor="white", edgecolor=INK2, linewidth=0.7, zorder=2)
        ax.scatter(g_, i, s=16, color=ACCENT, zorder=3)
        labels.append(f"{lab} ({n})")
        rows.append({"panel": "a", "set": lab, "n_mmr_deficient": n, "auroc_A": a_, "auroc_G": g_, "tests": pid})
    off = len(pooled) + 0.6
    for j, r in enumerate(big.itertuples()):
        y = off + j
        ax.plot([r.auroc_A, r.auroc_G], [y, y], color=GRID, linewidth=1.2, zorder=1)
        ax.scatter(r.auroc_A, y, s=12, facecolor="white", edgecolor=INK2, linewidth=0.6, zorder=2)
        ax.scatter(r.auroc_G, y, s=12, color=ACCENT, zorder=3)
        labels.append(f"{r.project} ({r.mmr_deficient} of {r.tumours})")
        rows.append({"panel": "a", "set": r.project, "n_mmr_deficient": r.mmr_deficient, "n_tumours": r.tumours,
                     "auroc_A": r.auroc_A, "auroc_G": r.auroc_G, "tests": "exploratory"})
    ys = list(range(len(pooled))) + [off + j for j in range(len(big))]
    ax.axvline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    rows.append({"panel": "a", "set": "reference line", "auroc_A": 0.5, "auroc_G": 0.5})
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=5.5)
    ax.set_ylim(max(ys) + 0.6, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0.35, 1.0)
    ax.set_xlabel("AUROC, MMR-deficient lower", fontsize=6)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", markersize=3.5, markerfacecolor="white",
                              markeredgecolor=INK2, label="score A (D-52)"),
                       Line2D([], [], marker="o", linestyle="", markersize=3.5, color=ACCENT, label="score G")],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=5.2, handletextpad=0.2,
              borderaxespad=0.1, columnspacing=0.8)
    # b. exon shifts in held-out cancers against the three training cohorts
    bx = fig.add_subplot(gs[1])
    down, up, _ = step.arms(el.acc_T.to_numpy(), el.d_bp.to_numpy())
    v = el.dropna(subset=["train_delta_adj"])
    dn, upv = step.arms(v.acc_T.to_numpy(), v.d_bp.to_numpy())[:2]
    other = ~(dn | upv)
    for k, col, lab, z, s_ in [(other, GRID, None, 1, 3), (dn, MUTED, "runs 10–19", 2, 3),
                               (upv, ACCENT, "runs ≥ 20, bp ≥ 36 nt", 3, 5)]:
        bx.scatter(v.train_delta_adj[k], v.delta_adj[k], s=s_, color=col, linewidths=0, alpha=0.8, zorder=z, label=lab)
    rows += [{"panel": "b", "set": "exon", "feature_id": r.feature_id, "symbol": r.symbol, "acc_T": r.acc_T,
              "d_bp": r.d_bp, "train_delta_adj": r.train_delta_adj, "heldout_delta_adj": r.delta_adj}
             for r in v.itertuples()]
    bx.axhline(0, color=INK2, linewidth=0.5, zorder=0)
    bx.axvline(0, color=INK2, linewidth=0.5, zorder=0)
    rows += [{"panel": "b", "set": "reference line", "train_delta_adj": 0.0}, {"panel": "b", "set": "reference line",
                                                                                "heldout_delta_adj": 0.0}]
    bx.text(0.03, 0.97, f"ρ = {pr.loc['PC6'].effect:.2f}", transform=bx.transAxes, ha="left", va="top", fontsize=5.6)
    bx.set_xlabel("MSI shift, UCEC + STAD + COAD", fontsize=6)
    bx.set_ylabel("MMR-deficient shift, 30 other cancers", fontsize=6)
    bx.legend(loc="lower right", bbox_to_anchor=(1.13, 0.0), fontsize=5.0, handletextpad=0.1, borderaxespad=0.1,
              markerscale=1.6)
    # c. dose within MMR-deficient tumours of the held-out cancers
    cx = fig.add_subplot(gs[2])
    m_ = sc[sc.mmr]
    cx.scatter(np.log1p(m_.n_hp), m_.pct_G, s=9, color=ACCENT, alpha=0.8, linewidths=0)
    local.append(m_[["project", "n_hp", "pct_G", "purity"]])
    d = exo.loc["dose: partial Spearman(G percentile, log(1 + homopolymer deletions) | purity), MMR-deficient"]
    cx.text(0.97, 0.97, f"ρ = {signed(d.value, 2)}\np = {fp(d.p)}", transform=cx.transAxes, ha="right", va="top", fontsize=5.6)
    rows.append({"panel": "c", "set": "partial Spearman | purity (exploratory)", "value": d.value, "p": d.p, "n": int(d.n)})
    cx.set_xticks(np.log1p([10, 30, 100, 300, 1000]))
    cx.set_xticklabels(["10", "30", "100", "300", "1000"], fontsize=5.6)
    cx.set_xlabel("homopolymer deletions", fontsize=6)
    cx.set_ylabel("score G, percentile in its cancer", fontsize=6)
    cx.set_ylim(-0.02, 1.02)
    for a_, lab, x in [(ax, "a", -0.95), (bx, "b", -0.36), (cx, "c", -0.36)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    for r in pr.reset_index().itertuples(index=False):
        rows.append({"panel": "tests", "set": r.id, "value": r.effect, "p": r.p, "q": r.q,
                     "q_family": "PC (msi_code_pancancer, 8 tests)", "tests": r.verdict})
    save(fig, OUT, "U18_msi_code_pancancer", data=pd.DataFrame(rows), caption=CAPTIONS["U18"](pr, per, exo, sc),
         local_data=(SCFG.base_data_dir / "figure_data", pd.concat(local, ignore_index=True)))


def _cap_u18(pr, per, exo, sc):
    big = per[per.mmr_deficient >= 3]
    a15 = exo.loc[[i for i in exo.index if i.startswith("A-15")], "value"]
    return f"""**Does the sequence code find mismatch-repair-deficient tumours in 30 cancers never used to define it?**

Pre-specified (config/followup.yaml: msi_code_pancancer, D-66) before any tumour outside UCEC, STAD and COAD was read for
MSI. Tumours: one per patient in the 30 other TCGA projects of the pan-cancer SE matrix, with MC3 calls ({len(sc):,}).
Mismatch-repair deficient (MMR-deficient) = ≥ 10 MC3 homopolymer deletions, the D-61 count; in the three training
cohorts that cut separates MSI from the rest almost perfectly (129 of 129 against 2 of 571). {int(sc.mmr.sum())} tumours
qualify. Hypermutated proficient = fewer than 10 such deletions and ≥ 500 SNVs ({int(sc.hyper.sum())}). Score A = D-52's
repeat-exon score; score G = the code: exons after a T-run of 10–19 nt (R_down, read as falling) against exons after a
run ≥ 20 nt whose branch point lies ≥ 36 nt from the AG (R_up, read as rising), each against matched background exons.
Neither score was re-fitted; each was computed within each cancer and ranked there.

a, AUROC for MMR-deficient tumours scoring lower (score A hollow, G filled): pooled over the 30 cancers, without the
three tissue relatives of the training cohorts, against hypermutated proficient tumours only, and per cancer with ≥ 3
MMR-deficient tumours (MMR-deficient of all tumours in brackets). b, per-exon shift (MMR-deficient minus the rest,
background-adjusted percentile) in the 30 cancers against the mean MSI shift of the three training cohorts; R exons
measured in both (grey = R_down, blue = R_up, light grey = runs ≥ 20 in neither arm). c, MMR-deficient tumours of the 30 cancers: score G percentile against their homopolymer deletions
(exploratory).

- PC1, score A: AUROC {pr.loc['PC1'].effect + 0.5:.2f}, q = {fp(pr.loc['PC1'].q)}, not resolved.
- PC2, score G: AUROC {pr.loc['PC2'].effect + 0.5:.2f}, q = {fp(pr.loc['PC2'].q)}, supported; PC3, G ahead of A by
  {pr.loc['PC3'].effect:.2f}, supported.
- PC4, without READ, ESCA and UCS: {pr.loc['PC4'].effect + 0.5:.2f}; PC5, against hypermutated proficient tumours:
  {pr.loc['PC5'].effect + 0.5:.2f}. Both supported: the code reads mismatch-repair deficiency, not mutation burden.
- PC6, exon shifts replicate: ρ = {pr.loc['PC6'].effect:.2f} over {pr.loc['PC6'].n}; PC7, R_up exons rise more
  than R_down exons: AUROC {pr.loc['PC7'].effect + 0.5:.2f}; PC8, the branch-point distance, adjusted:
  ρ = {pr.loc['PC8'].effect:.2f}. All supported.
- Per cancer (exploratory, {len(big)} cancers): score G {big.auroc_G.min():.2f}–{big.auroc_G.max():.2f}, score A
  {big.auroc_A.min():.2f}–{big.auroc_A.max():.2f}. By arm (A-15, exploratory): R_down {a15.iloc[0]:.2f}, R_up
  {1 - a15.iloc[1]:.2f} (higher in MMR-deficient).

The DNA call is a proxy from exome calls, not a clinical MSI test. q over family PC (8 tests). Per-tumour values are
kept locally (D-12)."""


def U19() -> None:
    sp = pd.read_csv(CODE_DIR / "specificity" / "predictions.csv").set_index("id")
    sx = pd.read_csv(CODE_DIR / "specificity" / "exploratory.csv")
    eg = pd.read_csv(CODE_DIR / "endometrium" / "predictions.csv").set_index("id")
    ee = pd.read_csv(CODE_DIR / "endometrium" / "exploratory.csv")
    bp = pd.read_csv(FU_DIR / "branch" / "predictions.csv").set_index("id")
    fig = plt.figure(figsize=(DOUBLE * 0.95, 2.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 0.95, 0.75], wspace=0.95)
    rows = []
    # a. share of R exons more included in MSI, by 10-nt branch-point bins
    ax = fig.add_subplot(gs[0])
    b = sx[sx.what.eq("share up by 10-nt d_bp bin, mean over >= 2 cohorts")]
    for cls, col, mk in [("10-19", INK2, "o"), (">=20", ACCENT, "s")]:
        v = b[b.run_class.eq(cls)]
        x = v.bin.str.split("-").str[0].astype(int) + 5
        ax.plot(x, v.share_up, color=col, linewidth=0.8, marker=mk, ms=3, mec="white", mew=0.3,
                label="T-run 10–19" if cls == "10-19" else "T-run ≥ 20")
        rows += [{"panel": "a", "set": f"T-run {cls}", "d_bp_bin": r.bin, "n": int(r.n), "share_up": r.share_up}
                 for r in v.itertuples()]
    ax.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    ax.axvline(36, color=GRID, linewidth=0.8, zorder=0)
    rows += [{"panel": "a", "set": "reference line", "share_up": 0.5}, {"panel": "a", "set": "rule cut", "d_bp_bin": "36"}]
    ax.set_xlabel("branch point to 3′ss (nt)", fontsize=6)
    ax.set_ylabel("share more included in MSI", fontsize=6)
    ax.set_xlim(10, 100)
    ax.set_ylim(0, 0.9)
    ax.legend(loc="lower right", fontsize=5.2, handlelength=1.4, borderaxespad=0.1)
    # b. adjusted rho of the shift with the branch-point distance, MSI and controls
    bx = fig.add_subplot(gs[1])
    def num(txt, key):
        return float(txt.split(key)[1].split(")")[0].split(",")[0].split(";")[0])
    sets = [("R exons, MSI", bp.loc["BP4"].effect, ACCENT, "BP4"),
            ("B exons, MSI", num(sp.loc["BS2"].n, "B 22,100 (rho "), INK2, "BS2"),
            ("UCEC MSI", num(sp.loc["BS3"].n, "MSI rho "), ACCENT, "BS3"),
            ("UCEC POLE", num(sp.loc["BS3"].n, "POLE vs rest rho "), INK2, "BS3"),
            ("STAD MSI", num(sp.loc["BS4"].n, "MSI rho "), ACCENT, "BS4"),
            ("STAD CIN vs GS", num(sp.loc["BS4"].n, "CIN vs GS rho "), INK2, "BS4")]
    for i, (lab, v, col, pid) in enumerate(sets):
        y = i + (i // 2) * 0.5
        bx.scatter(v, y, s=16, color=col, zorder=3)
        rows.append({"panel": "b", "set": lab, "rho": v, "test": pid})
        if i % 2 == 1:
            t = sp.loc[pid]
            bx.text(0.235, y - 0.5, f"p = {fp(t.p)}", va="center", fontsize=5.2)
            rows.append({"panel": "b", "set": f"{pid}: difference", "rho": t.effect, "p": t.p, "q": t.q,
                         "q_family": "BS (msi_code_specificity, 4 tests)", "test": t.verdict})
    bx.axvline(0, color=INK2, linewidth=0.6)
    rows.append({"panel": "b", "set": "reference line", "rho": 0.0})
    yt = [i + (i // 2) * 0.5 for i in range(len(sets))]
    bx.set_yticks(yt)
    bx.set_yticklabels([s_[0] for s_ in sets], fontsize=5.5)
    bx.set_ylim(max(yt) + 0.6, -0.6)
    bx.set_xlim(-0.12, 0.22)
    bx.set_xticks([-0.1, 0, 0.1, 0.2])
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.set_xlabel("adjusted ρ with branch-point distance", fontsize=6)
    # c. AUROC of score A and score G, MSI against the rest, in the three cohorts
    cx = fig.add_subplot(gs[2])
    coh = ["UCEC_TCGA", "STAD_TCGA", "COAD_TCGA"]
    for i, c in enumerate(coh):
        v = ee[ee.cohort.eq(c)].set_index("what")
        a_ = v.loc["AUROC(score lower in MSI than the rest), A", "value"]
        g_ = v.loc["AUROC(score lower in MSI than the rest), G", "value"]
        cx.plot([0, 1], [a_, g_], color=COH_STYLE[c][0], linewidth=0.8, marker=COH_STYLE[c][1], ms=3.5, mec="white",
                mew=0.3, label=_short(c))
        rows.append({"panel": "c", "set": c, "auroc_A": a_, "auroc_G": g_,
                     "test": "EG1, EG2" if c == "UCEC_TCGA" else "in-sample for the rule (exploratory)"})
    cx.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    cx.set_xticks([0, 1])
    cx.set_xticklabels(["score A", "score G"], fontsize=5.8)
    cx.set_xlim(-0.3, 1.3)
    cx.set_ylim(0.4, 1.0)
    cx.set_ylabel("AUROC, MSI lower", fontsize=6)
    cx.legend(loc="upper left", fontsize=5.2, handlelength=1.4, borderaxespad=0.1)
    for a_, lab, x in [(ax, "a", -0.3), (bx, "b", -0.62), (cx, "c", -0.5)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    save(fig, OUT, "U19_msi_code_specificity", data=pd.DataFrame(rows), caption=CAPTIONS["U19"](sp, eg, ee, sx))


def _cap_u19(sp, eg, ee, sx):
    b = sx[sx.what.eq("share up by 10-nt d_bp bin, mean over >= 2 cohorts")].set_index(["run_class", "bin"])
    return f"""**Is the code about the branch point, and about mismatch repair?**

Pre-specified (config/followup.yaml: msi_code_specificity, D-65; msi_code_endometrium, D-67). Exon shift = D-52's
background-adjusted difference in inclusion percentile. Adjusted ρ = D-63's: shift and branch-point distance
rank-residualised on run-length bins and mean inclusion.

a, exploratory: share of repeat exons more included in MSI (mean over ≥ 2 of the three cancers), by 10-nt bins of
branch-point distance, bins with ≥ 10 exons; grey line = the rule's 36-nt cut for runs ≥ 20. Exons whose branch point
lies within 30 nt of the AG move up rarely ({b.loc[('>=20', '20-29')].share_up:.2f} of runs ≥ 20 at 20–29 nt, n =
{int(b.loc[('>=20', '20-29')].n)}; {b.loc[('10-19', '20-29')].share_up:.2f} of runs 10–19); from 30–39 nt
{b.loc[('>=20', '30-39')].share_up:.2f} and {b.loc[('10-19', '30-39')].share_up:.2f}. b, adjusted ρ for MSI (blue) and
the matching control (grey): repeat exons (R, acceptor T-run ≥ 10) and exons without a tract (B, ≤ 5) measured in ≥ 2
cancers; UCEC MSI and POLE, each against the rest; STAD MSI against the rest and CIN against GS. p = the test of the
difference. c, AUROC of score A and score G for MSI against the rest.

- BS1, the distance from the run's far end to the branch point adds beyond the run: ρ = {sp.loc['BS1'].effect:.3f},
  q = {fp(sp.loc['BS1'].q)}, supported.
- BS2, exons without a tract (acceptor T-run ≤ 5) show no such relation: R minus B {signed(sp.loc['BS2'].effect, 3)},
  q = {fp(sp.loc['BS2'].q)}, supported.
- BS3, POLE (hypermutated, mismatch-repair proficient) against MSI in UCEC: difference
  {signed(sp.loc['BS3'].effect, 3)}, q = {fp(sp.loc['BS3'].q)}, supported. BS4, CIN against GS in STAD:
  {signed(sp.loc['BS4'].effect, 3)}, q = {fp(sp.loc['BS4'].q)}, supported.
- EG1, score G separates UCEC MSI: AUROC {eg.loc['EG1'].effect + 0.5:.2f}, q = {fp(eg.loc['EG1'].q)}; EG2, ahead of
  score A ({eg.loc['EG2'].effect:+.2f}), q = {fp(eg.loc['EG2'].q)}. Both supported.

The rule behind score G was read from these three cohorts (D-53 and U17b), so c is in-sample for the rule, and its STAD
and COAD values are exploratory; U18 is the out-of-sample test. q over families BS (4 tests) and EG (2 tests)."""


# ============================================================ infiltration (D-68)

def U20() -> None:
    pr = pd.read_csv(CODE_DIR / "infiltration" / "predictions.csv").set_index("id")
    ex = pd.read_csv(CODE_DIR / "infiltration" / "exploratory.csv").set_index("what")
    sc = pd.read_parquet(MSI_LOCAL / "pancancer_covariates.parquet")
    step = _step("51_msi_code.py")
    fig = plt.figure(figsize=(DOUBLE * 0.8, 2.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=1.25)
    rows = []
    # a. AUROC of score G by immune half, and of the immune score alone
    ax = fig.add_subplot(gs[0])
    low = sc.immune <= sc.groupby("project").immune.transform("median")
    mm = sc.mmr.to_numpy()
    sets = [("30 cancers, all", step.auc_lower(sc.pct_G.to_numpy(), mm), int(mm.sum()), ACCENT, "PC2"),
            ("lower-immune half", pr.loc["IF3"].effect + 0.5, int((mm & low).sum()), ACCENT, "IF3"),
            ("upper-immune half", step.auc_lower(sc.pct_G.to_numpy()[~low], mm[~low]), int((mm & ~low).sum()), ACCENT,
             "exploratory"),
            ("immune score alone", ex.loc["AUROC of the immune score itself, higher in MMR-deficient (held-out)", "value"],
             int(mm.sum()), INK2, "exploratory"),
            ("UCEC, lower-immune half", ex.loc["IF4: AUROC(G lower in MSI), lower-immune half", "value"],
             int(ex.loc["IF4: AUROC(G lower in MSI), lower-immune half", "n"]), COH_STYLE["UCEC_TCGA"][0], "exploratory"),
            ("STAD, lower-immune half", ex.loc["IF5: AUROC(G lower in MSI), lower-immune half", "value"],
             int(ex.loc["IF5: AUROC(G lower in MSI), lower-immune half", "n"]), COH_STYLE["STAD_TCGA"][0], "exploratory"),
            ("COAD, lower-immune half", ex.loc["IF6: AUROC(G lower in MSI), lower-immune half", "value"],
             int(ex.loc["IF6: AUROC(G lower in MSI), lower-immune half", "n"]), COH_STYLE["COAD_TCGA"][0], "exploratory")]
    ys = [0, 1, 2, 3, 4.6, 5.6, 6.6]
    for y, (lab, v, n, col, test) in zip(ys, sets):
        ax.scatter(v, y, s=16, color=col, zorder=3)
        rows.append({"panel": "a", "set": lab, "n_mmr_deficient_or_msi": n, "auroc": v, "test": test})
    ax.axvline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    rows.append({"panel": "a", "set": "reference line", "auroc": 0.5})
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{l} ({n})" for l, _, n, _, _ in sets], fontsize=5.6)
    ax.set_ylim(max(ys) + 0.6, -0.6)
    ax.set_xlim(0.4, 1.0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("AUROC, MMR-deficient or MSI lower", fontsize=6)
    # b. the MMR (MSI) coefficient after adjustment, as a share of the unadjusted one
    bx = fig.add_subplot(gs[1])
    unadj = {"IF1": ex.loc["IF1: unadjusted coefficient (cancer fixed effects only)", "value"],
             "IF2": ex.loc["IF2: unadjusted coefficient (cancer fixed effects only)", "value"]}
    for k in ["IF4", "IF5", "IF6"]:
        unadj[k] = float(pr.loc[k].n.split("site only ")[1])
    labels = {"IF1": "30 cancers: stromal, immune, RI", "IF2": "30 cancers: + ABSOLUTE purity",
              "IF4": "UCEC: site, stromal, immune, RI", "IF5": "STAD: site, stromal, immune, RI",
              "IF6": "COAD: site, stromal, immune, RI"}
    for y, k in enumerate(["IF1", "IF2", "IF4", "IF5", "IF6"]):
        r, u = pr.loc[k], unadj[k]
        lo, hi = sorted([r.ci_low / u, r.ci_high / u])
        bx.hlines(y, lo, hi, color=ACCENT, linewidth=0.8)
        bx.scatter(r.effect / u, y, s=16, color=ACCENT, zorder=3)
        bx.text(1.62, y, f"p = {fp(r.p)}", va="center", fontsize=5.2)
        rows.append({"panel": "b", "set": labels[k], "test": k, "coef_adjusted": r.effect, "coef_unadjusted": u,
                     "share_kept": r.effect / u, "share_ci_low": lo, "share_ci_high": hi, "p": r.p, "q": r.q,
                     "q_family": "IF (msi_code_infiltration, 6 tests)", "verdict": r.verdict})
    bx.axvline(1.0, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    bx.axvline(0.0, color=INK2, linewidth=0.6, zorder=0)
    rows += [{"panel": "b", "set": "reference line: unadjusted", "share_kept": 1.0},
             {"panel": "b", "set": "reference line: none kept", "share_kept": 0.0}]
    bx.set_yticks(range(5))
    bx.set_yticklabels([labels[k] for k in ["IF1", "IF2", "IF4", "IF5", "IF6"]], fontsize=5.6)
    bx.set_ylim(4.6, -0.6)
    bx.set_xlim(-0.05, 1.6)
    bx.set_xticks([0, 0.5, 1.0, 1.5])
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.set_xlabel("effect kept after adjustment", fontsize=6)
    for a_, lab, x in [(ax, "a", -1.05), (bx, "b", -1.2)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    save(fig, OUT, "U20_msi_code_infiltration", data=pd.DataFrame(rows), caption=CAPTIONS["U20"](pr, ex, rows))


def _cap_u20(pr, ex, rows):
    b = pd.DataFrame([r for r in rows if r["panel"] == "b" and "test" in r])
    a = pd.DataFrame([r for r in rows if r["panel"] == "a" and "test" in r]).set_index("set")
    return f"""**Is score G an infiltration readout?**

Pre-specified (config/followup.yaml: msi_code_infiltration, D-68). Covariates as in D-55: stromal and immune = mean
log2 CPM of eight stromal and eight immune genes (config/biology.yaml); RI = the tumour's mean retained-intron PSI, an
RNA-handling proxy; purity = ABSOLUTE, where called. In the 30 held-out cancers (D-66) each covariate and score G are
percentiles within the tumour's cancer.

(a) AUROC for MMR-deficient tumours (MSI in the three cohorts) scoring lower on G: all held-out tumours, then the
lower and upper halves of each cancer's immune score; grey = the immune score itself, read as higher in MMR-deficient
tumours; then MSI against the rest in the lower-immune half of UCEC, STAD and COAD (MMR-deficient or MSI tumours in
brackets). (b) OLS coefficient of MMR deficiency (or MSI) on score G after adjustment, as a share of the coefficient
with cancer (or site) fixed effects only; bar = 95% HC3 interval; p of the adjusted coefficient.

- IF1, IF2, held-out cancers: adjusted coefficient {signed(pr.loc['IF1'].effect, 3)} and
  {signed(pr.loc['IF2'].effect, 3)} (unadjusted {signed(b.set_index('test').loc['IF1', 'coef_unadjusted'], 3)});
  q = {fp(pr.loc['IF1'].q)}, supported.
- IF3, lower-immune half: AUROC {pr.loc['IF3'].effect + 0.5:.2f} ({a.loc['lower-immune half', 'n_mmr_deficient_or_msi']}
  MMR-deficient), supported; upper half {a.loc['upper-immune half', 'auroc']:.2f} (exploratory).
- IF4–IF6, the three cohorts with site fixed effects: UCEC {signed(pr.loc['IF4'].effect, 3)}, STAD
  {signed(pr.loc['IF5'].effect, 3)}, COAD {signed(pr.loc['IF6'].effect, 3)}; all supported, each keeping
  {b[b.test.isin(['IF4', 'IF5', 'IF6'])].share_kept.min():.2f}–{b[b.test.isin(['IF4', 'IF5', 'IF6'])].share_kept.max():.2f}
  of the site-only effect.

The immune score marks MMR-deficient tumours only weakly here (AUROC {a.loc['immune score alone', 'auroc']:.2f}), and in
the adjusted model more immune expression goes with a higher G, the opposite of the MMR effect. q over family IF (6
tests). Per-tumour values are kept locally (D-12)."""


def U21() -> None:
    pr = pd.read_csv(CODE_DIR / "msicalls" / "predictions.csv").set_index("id")
    ex = pd.read_csv(CODE_DIR / "msicalls" / "exploratory.csv").set_index("what")
    sc = pd.read_parquet(MSI_LOCAL / "pancancer_msicalls.parquet")
    step = _step("51_msi_code.py")
    fig = plt.figure(figsize=(DOUBLE, 2.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.0], wspace=0.45)
    rows, local = [], []
    rng = np.random.default_rng(0)
    # a. AUROC of score G against each DNA call
    ax = fig.add_subplot(gs[0])
    mm = sc.mmr.to_numpy()
    sets = [("MC3 deletions", step.auc_lower(sc.pct_G.to_numpy(), mm), int(mm.sum()), ACCENT, "PC2"),
            ("MANTIS MSI-H", pr.loc["MC1"].effect + 0.5, int(pr.loc["MC1"].n.split()[1].rstrip(",")), ACCENT, "MC1"),
            ("MSIsensor MSI-H", pr.loc["MC2"].effect + 0.5, int(pr.loc["MC2"].n.split()[1].rstrip(",")), ACCENT, "MC2")]
    for c in ["UCEC_TCGA", "STAD_TCGA", "COAD_TCGA"]:
        for col, lab in [("mantis", "MANTIS"), ("msisensor", "MSIsensor")]:
            r = ex.loc[f"{c} (in-sample for the rule): AUROC(G lower in {col} MSI-H than MSS)"]
            sets.append((f"{_short(c)}, {lab}", r.value, int(r.n), COH_STYLE[c][0], "exploratory, in-sample"))
    ys = [0, 1, 2, 3.6, 4.6, 5.6, 6.6, 7.6, 8.6]
    for y, (lab, v, n, col, test) in zip(ys, sets):
        ax.scatter(v, y, s=16, color=col, zorder=3)
        rows.append({"panel": "a", "set": lab, "n_positive": n, "auroc": v, "test": test})
    ax.axvline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    rows.append({"panel": "a", "set": "reference line", "auroc": 0.5})
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{l} ({n})" for l, _, n, _, _ in sets], fontsize=5.6)
    ax.set_ylim(max(ys) + 0.6, -0.6)
    ax.set_xlim(0.4, 1.0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("AUROC, positive tumours lower", fontsize=6)
    # b, c. score G where the MC3 call and the exome caller agree or disagree (held-out tumours)
    axes = []
    for k, (col, lab) in enumerate([("mantis", "MANTIS"), ("msisensor", "MSIsensor")]):
        bx = fig.add_subplot(gs[k + 1], sharey=axes[0] if axes else None)
        axes.append(bx)
        x, ticks = 0, []
        for mc3, mlab in [(False, "MC3 proficient"), (True, "MC3 deficient")]:
            x0 = x
            for call in ["MSS", "between", "MSI-H"]:
                v = sc.loc[sc.mmr.eq(mc3) & sc[f"call_{col}"].eq(call), "pct_G"].to_numpy()
                if not len(v):
                    continue
                c_ = ACCENT if mc3 else INK2
                if len(v) >= 5:
                    st = _box(bx, x, v, c_, dots=len(v) < 200, rng=rng)
                else:
                    bx.scatter(x + rng.uniform(-0.12, 0.12, len(v)), v, s=2.2, color=c_, alpha=0.45, linewidths=0)
                    st = {"n": len(v), "median": np.nan}
                rows.append({"panel": "bc"[k], "caller": lab, "mc3_call": mlab, "caller_call": call, **st})
                local.append(pd.DataFrame({"panel": "bc"[k], "caller": lab, "mc3_call": mlab, "caller_call": call,
                                           "file_id": sc.loc[sc.mmr.eq(mc3) & sc[f"call_{col}"].eq(call), "file_id"],
                                           "pct_G": v}))
                ticks.append((x, f"{'mid' if call == 'between' else call}\n({len(v):,})"))
                x += 1
            bx.text((x0 + x - 1) / 2, -0.19, mlab, transform=bx.get_xaxis_transform(), ha="center", va="top",
                    fontsize=5.8)
            x += 0.6
        bx.set_xticks([t for t, _ in ticks])
        bx.set_xticklabels([l for _, l in ticks], fontsize=5.4)
        bx.tick_params(axis="x", length=0)
        bx.set_xlim(-0.6, x - 1.0)
        bx.set_ylim(-0.03, 1.03)
        bx.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        bx.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
        bx.set_title(lab, fontsize=6.2, pad=3)
        if k == 0:
            bx.set_ylabel("score G percentile in its cancer", fontsize=6)
        else:
            plt.setp(bx.get_yticklabels(), visible=False)
    rows.append({"panel": "bc", "caller": "reference line", "median": 0.5})
    for a_, lab, x in [(ax, "a", -0.75), (axes[0], "b", -0.25), (axes[1], "c", -0.12)]:
        a_.text(x, 1.04, lab, transform=a_.transAxes, fontsize=8, fontweight="bold")
    save(fig, OUT, "U21_msi_code_msi_calls", data=pd.DataFrame(rows), caption=CAPTIONS["U21"](pr, ex, rows),
         local_data=(SCFG.base_data_dir / "figure_data", pd.concat(local, ignore_index=True)))


def _cap_u21(pr, ex, rows):
    a = pd.DataFrame([r for r in rows if r["panel"] == "a" and "test" in r]).set_index("set")
    bc = pd.DataFrame([r for r in rows if r["panel"] in "bc" and r.get("mc3_call")])
    def cell(caller, mc3, call):
        r = bc[(bc.caller == caller) & (bc.mc3_call == mc3) & (bc.caller_call == call)]
        return int(r.n.iloc[0]) if len(r) else 0
    return f"""**Does score G follow MSI calls that do not come from MC3?**

Pre-specified (config/followup.yaml: msi_code_msi_calls, D-69) before the scores were downloaded. MANTIS (Bonneville
et al. 2017) and MSIsensor scores of the TCGA PanCancer Atlas, fetched from cBioPortal (step 00c) and matched by
sample barcode. Calls: MANTIS MSI-H > 0.6, MSS < 0.4; MSIsensor MSI-H > 10, MSS < 4; scores in between (mid) are kept
apart.
Tumours: D-66's 7,685 held-out tumours (30 cancers, none used to define the code); score G as there, a percentile
within each cancer, not re-fitted.

(a) AUROC for positive tumours scoring lower on G (positives in brackets): D-66's MC3 call (≥ 10 homopolymer
deletions) against the rest, MANTIS and MSIsensor MSI-H against MSS; then, in-sample for the rule, MSI-H against MSS
in UCEC, STAD and COAD. (b, c) Score G percentile of the held-out tumours by D-66's MC3 call and the caller's call (n
in brackets); box = interquartile range, diamond = median, whiskers = 1.5 × IQR, dots = tumours (not drawn for the
large MSS group); dashed line = the middle of the cancer.

- MC1, MANTIS MSI-H against MSS: AUROC {pr.loc['MC1'].effect + 0.5:.2f} ({a.loc['MANTIS MSI-H', 'n_positive']} MSI-H),
  q = {fp(pr.loc['MC1'].q)}, supported.
- MC2, MSIsensor: AUROC {pr.loc['MC2'].effect + 0.5:.2f} ({a.loc['MSIsensor MSI-H', 'n_positive']} MSI-H), q =
  {fp(pr.loc['MC2'].q)}, supported.
- MC3, Spearman of score G with the MANTIS score (percentile within cancer), all {pr.loc['MC3'].n}:
  ρ = {signed(pr.loc['MC3'].effect, 3)}, q = {fp(pr.loc['MC3'].q)}, not resolved. Over all tumours, nearly all of them
  MSS, the MANTIS score shows no relation to score G.

Agreement (exploratory): every MANTIS MSI-H tumour ({cell('MANTIS', 'MC3 deficient', 'MSI-H')}) is MC3-deficient; of
the {int(a.loc['MC3 deletions', 'n_positive'])} MC3-deficient tumours, MANTIS calls {cell('MANTIS', 'MC3 deficient', 'between')}
mid and {cell('MANTIS', 'MC3 deficient', 'MSS')} MSS, and these still score low on G. MSIsensor calls
{cell('MSIsensor', 'MC3 proficient', 'MSI-H')} MC3-proficient tumours MSI-H. q over family MC (3 tests). Per-tumour
values are kept locally (D-12)."""


# ============================================================ D-70 the rising arm
def _wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return np.nan, np.nan
    z, ph = 1.96, k / n
    c = (ph + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def U22() -> None:
    d = CODE_DIR / "reversal"
    pr = pd.read_csv(d / "predictions.csv").set_index("id")
    ex = pd.read_csv(d / "exploratory.csv").set_index("what")
    ph = pd.read_csv(d / "posthoc.csv")
    con = pd.read_csv(d / "contrasts.csv")
    draws = pd.read_csv(d / "held_out_draws.csv")
    el = pd.read_csv(d / "exon_shifts.csv")
    tu = pd.read_parquet(MSI_LOCAL / "reversal_tumours.parquet")
    dose = pd.read_parquet(MSI_LOCAL / "reversal_dose.parquet")
    fig = plt.figure(figsize=(DOUBLE - 0.2, 4.3))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1], hspace=0.62, wspace=1.2)
    rows, local = [], []
    rng = np.random.default_rng(0)
    # a. the held-out rise against every other subtype contrast (20 v 20 tumours, 50 draws)
    ax = fig.add_subplot(gs[0, :3])
    null = con[con.kind.eq("other")]
    v = null.stat.abs().to_numpy() * 100
    ax.scatter(v, rng.uniform(-0.18, 0.18, len(v)), s=7, color=MUTED, alpha=0.8, linewidths=0, zorder=2)
    for r in null.itertuples():
        rows.append({"panel": "a", "set": "other subtype contrast", "cohort": r.cohort, "contrast": f"{r.a} vs {r.b}",
                     "abs_shift_pct_points": abs(r.stat) * 100})
    lo, hi = np.percentile(draws.stat * 100, [2.5, 97.5])
    s_h = pr.loc["RV1"].effect * 100
    ax.hlines(1, lo, hi, color=ACCENT, linewidth=0.9)
    ax.scatter(s_h, 1, s=18, color=ACCENT, zorder=3)
    rows.append({"panel": "a", "set": "MMR-deficient, 30 held-out cancers", "shift_pct_points": s_h,
                 "draw_2.5": lo, "draw_97.5": hi})
    for k, r in enumerate(con[con.kind.eq("msi")].itertuples()):
        y = 2 + k
        ax.scatter(r.stat * 100, y, s=16, facecolors="white", edgecolors=COH_STYLE[r.cohort][0], linewidths=0.9, zorder=3)
        rows.append({"panel": "a", "set": f"MSI, {_short(r.cohort)} (in-sample)", "shift_pct_points": r.stat * 100})
    q95 = np.quantile(v, 0.95)
    ax.axvline(q95, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    ax.axvline(0, color=INK2, linewidth=0.6, zorder=0)
    rows.append({"panel": "a", "set": "reference line: 95th percentile of the other contrasts", "shift_pct_points": q95})
    labs = [f"other subtype contrasts ({len(null)})", "MMR-deficient, held-out"] + \
           [f"MSI, {_short(c)} (in-sample)" for c in con[con.kind.eq("msi")].cohort]
    ax.set_yticks(range(len(labs)))
    ax.set_yticklabels(labs, fontsize=5.6)
    ax.set_ylim(len(labs) - 0.4, -0.6)
    ax.set_xlim(-0.3, 6.2)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("R_up shift, percentile points", fontsize=6)
    # b. held-out exon shifts by class
    bx = fig.add_subplot(gs[0, 3:])
    order = [("R_down", "runs 10–19", INK2), ("R_up", "runs ≥ 20,\nbranch point ≥ 36", ACCENT),
             ("B_far", "no run,\nbranch point ≥ 36", MUTED), ("shield", "run upstream\nof branch point", ACCENT_DARK)]
    for x, (c, lab, col) in enumerate(order):
        w = el.loc[el.cls.eq(c), "delta_adj"].to_numpy() * 100
        st = _box(bx, x, w, col, dots=len(w) < 300, rng=rng)
        rows.append({"panel": "b", "class": c, **st})
    bx.axhline(0, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    bx.set_xticks(range(len(order)))
    bx.set_xticklabels([f"{lab}\n({int(el.cls.eq(c).sum()):,})" for c, lab, _ in order], fontsize=5.4)
    bx.tick_params(axis="x", length=0)
    bx.set_ylim(-16, 21)
    bx.set_ylabel("exon shift, percentile points", fontsize=6)
    # c. runs >= 20 by branch-point distance (A-16, exploratory)
    cx = fig.add_subplot(gs[1, :2])
    b = ph[ph.what.eq("runs >= 20 by branch-point distance")].reset_index(drop=True)
    for x, r in b.iterrows():
        k = int(round(r.share_up * r.n))
        l_, h_ = _wilson(k, int(r.n))
        col = ACCENT if r.bin != "18-29" else INK2
        cx.vlines(x, l_, h_, color=col, linewidth=0.8)
        cx.scatter(x, r.share_up, s=16, color=col, zorder=3)
        rows.append({"panel": "c", "branch_point_nt": r.bin, "n": int(r.n), "share_up": r.share_up, "wilson_low": l_,
                     "wilson_high": h_, "median_shift": r["median"]})
    cx.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    cx.set_xticks(range(len(b)))
    cx.set_xticklabels([f"{r.bin.replace('-', '–')}\n({int(r.n)})" for r in b.itertuples()], fontsize=5.4)
    cx.tick_params(axis="x", length=0)
    cx.set_ylim(0, 1)
    cx.set_xlabel("branch point to AG, nt", fontsize=6)
    cx.set_ylabel("share more included", fontsize=6)
    # d. the R_up arm by DNA class (held-out tumours)
    dx = fig.add_subplot(gs[1, 2:4])
    grp = [("proficient", (~tu.mmr & ~tu.hyper), MUTED), ("hypermutated\nproficient", tu.hyper, INK2),
           ("MMR-\ndeficient", tu.mmr, ACCENT)]
    for x, (lab, k, col) in enumerate(grp):
        w = tu.loc[k, "pct_up"].to_numpy()
        st = _box(dx, x, w, col, dots=len(w) < 600, rng=rng)
        rows.append({"panel": "d", "group": lab.replace("\n", " "), **st})
        local.append(pd.DataFrame({"panel": "d", "group": lab.replace("\n", " "), "file_id": tu.loc[k, "file_id"],
                                   "pct_up": w}))
    dx.axhline(0.5, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    dx.set_xticks(range(3))
    dx.set_xticklabels([f"{lab}\n({int(k.sum()):,})" for lab, k, _ in grp], fontsize=5.4)
    dx.tick_params(axis="x", length=0)
    dx.set_ylim(-0.03, 1.03)
    dx.set_ylabel("R_up arm, percentile in cancer", fontsize=6)
    # e. dose: R_up arm against homopolymer deletions in MMR-deficient / MSI tumours
    exx = fig.add_subplot(gs[1, 4:])
    for src, col, mk in [("three cohorts", MUTED, "o"), ("held-out", ACCENT, "o")]:
        w = dose[dose.source.eq(src)]
        exx.scatter(np.log1p(w.n_hp), w.pct_up, s=5, color=col, alpha=0.7, linewidths=0, marker=mk)
        local.append(pd.DataFrame({"panel": "e", "source": src, "n_hp": w.n_hp, "pct_up": w.pct_up}))
    exx.set_xticks(np.log1p([10, 30, 100, 300, 1000]))
    exx.set_xticklabels(["10", "30", "100", "300", "1000"], fontsize=5.6)
    exx.set_ylim(-0.03, 1.03)
    exx.set_xlabel("homopolymer deletions", fontsize=6)
    exx.set_ylabel("R_up arm, percentile in cancer", fontsize=6)
    rows.append({"panel": "e", "partial_spearman_given_purity": pr.loc["RV6"].effect, "p": pr.loc["RV6"].p,
                 "q": pr.loc["RV6"].q, "n": pr.loc["RV6"].n})
    from _figbase import panel_letters
    panel_letters([(ax, -0.5), (bx, -0.16), (cx, -0.3), (dx, -0.3), (exx, -0.3)])
    save(fig, OUT, "U22_msi_code_rising_arm", data=pd.DataFrame(rows), caption=CAPTIONS["U22"](pr, ex, ph, rows),
         local_data=(SCFG.base_data_dir / "figure_data", pd.concat(local, ignore_index=True)))


def _cap_u22(pr, ex, ph, rows):
    def line(i):
        r = pr.loc[i]
        return f"q = {fp(r.q)}, {r.verdict}"
    g = lambda w: ph[ph.what.eq(w)].iloc[0]
    b = ph[ph.what.eq("runs >= 20 by branch-point distance")].set_index("bin")
    oth = [r["abs_shift_pct_points"] for r in rows if r.get("set") == "other subtype contrast"]
    n_below = int(sum(v < pr.loc['RV1'].effect * 100 for v in oth))
    return f"""**Is the rising arm of the code real?**

Pre-specified (config/followup.yaml: msi_code_reversal, D-70) before any test was computed; panel (c) is exploratory
(A-16), added after the tests were read. Exon shift = D-66's background-adjusted difference in inclusion percentile,
MMR-deficient minus other tumours, pooled over the 30 held-out cancers (×100 = percentile points). R_up = acceptor
T-run ≥ 20 nt with the sequence-predicted branch point ≥ 36 nt from the AG; R_down = runs of 10–19 nt; B_far = no run
(≤ 5) with the branch point ≥ 36 nt; shield = a run ≥ 20 nt lying upstream of a branch point, so that contracting it
cannot move the branch point. R_up arm = each tumour's mean R_up percentile against matched exons without a run, as a
percentile within its cancer.

(a) Mean R_up shift in 20 against 20 tumours (50 draws): every pair of subtypes in the 22 cohorts with MSI, POLE and
MC3-deficient tumours removed (grey, absolute value; dashed = their 95th percentile), the held-out MMR-deficient tumours
each paired with a proficient tumour of the same cancer (blue; bar = 95% of draws) and, in-sample for the rule, MSI in
the three training cohorts (hollow). (b) Held-out exon shifts by class (n in brackets); box = interquartile range,
diamond = median, whiskers = 1.5 × IQR, dots = exons (not drawn for the large B_far set); the axis is clipped to
−16 to 21. (c) Exploratory: exons with a run ≥ 20 nt, share more included in MMR-deficient tumours by branch-point
distance (n in brackets); bars = Wilson 95% interval. (d) R_up arm of held-out tumours by DNA class; dots = tumours
(not drawn for the proficient group). (e) R_up arm against homopolymer deletions in MMR-deficient tumours (blue,
held-out) and MSI tumours of the three cohorts (grey).

- RV1, specificity: held-out {signed(pr.loc['RV1'].effect * 100, 2)} points, above {n_below} of the {len(oth)} other contrasts
  (largest {float(pr.loc['RV1'].n.split('largest |s| ')[1].split(',')[0]) * 100:.1f}, sarcoma LMS against MFS);
  {line('RV1')}. The three in-sample MSI contrasts lie above all or all but one.
- RV2, a tumour-wide inclusion axis: the MMR coefficient on the R_up arm is {signed(pr.loc['RV2'].effect, 3)} adjusted
  for it, stromal, immune and retained-intron levels (cancer only +0.237); {line('RV2')}.
- RV3, the tract is needed: R_up above B_far, AUROC {pr.loc['RV3'].effect + 0.5:.2f}; {line('RV3')}.
- RV4, AG exclusion zone in place of the branch point, runs ≥ 20: adjusted ρ = {signed(pr.loc['RV4'].effect, 3)};
  {line('RV4')}.
- RV5, R_up above shield exons: AUROC {pr.loc['RV5'].effect + 0.5:.2f}; {line('RV5')}. Shield exons move up in
  {int(round(g('shield exons against 0').share_up * 42))} of 42 (exploratory).
- RV6, dose (panel e): partial ρ = {signed(pr.loc['RV6'].effect, 3)} given purity; {line('RV6')}.
- RV7, mismatch repair, not mutation burden (panel d): AUROC {pr.loc['RV7'].effect + 0.5:.2f} against hypermutated
  proficient tumours; {line('RV7')}.
- RV8, PTB relief: Spearman of the shift with PTB sensitivity {signed(pr.loc['RV8'].effect, 3)}, opposite to the
  sign PTB relief predicts; {line('RV8')}.

Panel (c), exploratory: with the branch point 18–29 nt from the AG, {b.loc['18-29', 'share_up']:.2f} of long-run exons
move up (n = {int(b.loc['18-29', 'n'])}); from 30–35 nt {b.loc['30-35', 'share_up']:.2f}, the shape seen in the training
cohorts (Figure U19a). q over family RV (8 tests). Per-tumour values are kept locally (D-12)."""


# ============================================================ D-71 / D-72: the repeat-exon logistic model (A-19)
M4_LABEL = "logistic, repeat exons"


def _roc_points(score: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Step ROC (false-positive rate, sensitivity), ties as one step."""
    from sklearn.metrics import roc_curve
    f, t, _ = roc_curve(y, score)
    return f, t


def U23() -> None:
    d = CODE_DIR / "models"
    pr = pd.read_csv(d / "predictions.csv").set_index("id")
    met = pd.read_csv(d / "metrics.csv")
    ex = pd.read_csv(d / "exploratory.csv")
    cur = pd.read_csv(d / "pr_curves.csv")
    ph = pd.read_csv(d / "posthoc.csv")
    exons = pd.read_csv(d / "posthoc_exons.csv")
    p1 = pd.read_parquet(MSI_LOCAL / "models_e1_predictions.parquet")
    y = p1.mmr.to_numpy()
    fig = plt.figure(figsize=(DOUBLE - 0.2, 4.3))
    gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.62)
    rows = []
    styles = [("M4", ACCENT, 1.1, M4_LABEL), ("M0", MUTED, 0.9, "score G, no fit")]
    # a. ROC held out
    ax = fig.add_subplot(gs[0, 0])
    for m, col, lw, lab in styles:
        f, t = _roc_points(p1[m].to_numpy(), y)
        ax.plot(f, t, color=col, linewidth=lw, drawstyle="steps-post", label=lab)
        rows += [{"panel": "a", "model": m, "false_positive_rate": a, "sensitivity": b} for a, b in zip(f, t)]
    ax.plot([0, 1], [0, 1], color=GRID, linewidth=0.6, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("false-positive rate", fontsize=6)
    ax.set_ylabel("sensitivity", fontsize=6)
    ax.legend(loc="lower right", fontsize=5.4, frameon=False, handlelength=1.2)
    # b. precision-recall held out
    bx = fig.add_subplot(gs[0, 1])
    for m, col, lw, lab in styles:
        c = cur[cur.model.eq(m)].sort_values("recall")
        bx.plot(c.recall, c.precision, color=col, linewidth=lw, drawstyle="steps-post", label=lab)
        rows += [{"panel": "b", "model": m, "recall": a, "precision": b} for a, b in zip(c.recall, c.precision)]
    prev = float(met[met.split.eq("E1") & met.model.eq("M0")].prevalence.iloc[0])
    bx.axhline(prev, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    rows.append({"panel": "b", "model": "reference line: prevalence", "precision": prev})
    bx.set_xlim(0, 1)
    bx.set_ylim(0, 1.01)
    bx.set_xlabel("recall", fontsize=6)
    bx.set_ylabel("precision", fontsize=6)
    bx.legend(loc="upper right", fontsize=5.4, frameon=False, handlelength=1.2)
    # c. per cancer held out (>= 3 MMR-deficient)
    cx = fig.add_subplot(gs[1, 0])
    pc = ex[ex.what.eq("per-cancer AUROC, E1")]
    cancers = sorted(pc.project.unique())
    for i, c in enumerate(cancers):
        for m, col, dx_, mk in [("M0", MUTED, -0.12, "o"), ("M4", ACCENT, 0.12, "o")]:
            r = pc[pc.project.eq(c) & pc.model.eq(m)].iloc[0]
            cx.scatter(i + dx_, r.value, s=16, marker=mk, color=col, zorder=3)
            rows.append({"panel": "c", "project": c, "model": m, "auroc": r.value, "positives": r.positives, "n": r.n})
    npos = {c: int(pc[pc.project.eq(c)].positives.iloc[0]) for c in cancers}
    cx.set_xticks(range(len(cancers)))
    cx.set_xticklabels([f"{c}\n({npos[c]})" for c in cancers], fontsize=5.6)
    cx.set_xlim(-0.6, len(cancers) - 0.4)
    cx.set_ylim(0.5, 1.02)
    cx.set_ylabel("AUROC within cancer", fontsize=6)
    # d. the 20 largest standardised coefficients
    dx = fig.add_subplot(gs[1, 1])
    top = exons[exons.list.str.startswith("M4: 20 largest")].sort_values("rank")
    cls_col = {"run 10-19 (R_down)": INK2, "run >=20, bp >=36 nt (R_up)": ACCENT,
               "run >=20, bp <36 nt or none": ACCENT_LIGHT}
    top = top.reset_index(drop=True)
    for yy, r in top.iterrows():
        dx.barh(yy, r["value"], height=0.65, color=cls_col[r["class"]])
        rows.append({"panel": "d", "symbol": r["symbol"], "acc_T": r["acc_T"], "d_bp": r["d_bp"], "class": r["class"],
                     "coefficient": r["value"]})
    dx.axvline(0, color=INK2, linewidth=0.6)
    lab = [f"{r['symbol']} T{int(r['acc_T'])} " + (f"bp{int(r['d_bp'])}" if np.isfinite(r["d_bp"]) else "no bp")
           for _, r in top.iterrows()]
    dx.set_yticks(range(len(top)))
    dx.set_yticklabels(lab, fontsize=5.0)
    dx.set_ylim(len(top) - 0.4, -0.6)
    dx.tick_params(axis="y", length=0)
    dx.spines["left"].set_visible(False)
    dx.set_xlabel("standardised coefficient", fontsize=6)
    lo_ = float(top["value"].min())
    dx.set_xlim(lo_ - 0.02, 0.5)
    dx.set_xticks([-0.1, 0.0, 0.1, 0.2])
    from matplotlib.patches import Patch
    dx.legend(handles=[Patch(color=INK2, label="run 10–19"), Patch(color=ACCENT_LIGHT, label="run ≥ 20, bp < 36 or none"),
                       Patch(color=ACCENT, label="run ≥ 20, bp ≥ 36")], loc="lower right", fontsize=5.0, frameon=False,
              handlelength=0.9, borderaxespad=0.1)
    from _figbase import panel_letters
    panel_letters([(ax, -0.22), (bx, -0.22), (cx, -0.22), (dx, -0.5)])
    save(fig, OUT, "U23_msi_code_models", data=pd.DataFrame(rows), caption=CAPTIONS["U23"](pr, met, ex, ph))


def _cap_u23(pr, met, ex, ph):
    e1 = met[met.split.eq("E1")].set_index("model")
    e2 = met[met.split.eq("E2")].set_index("model")
    e30 = met[met.split.eq("E2, the 30 held-out cancers only")].set_index("model")
    n_tr = int(e2.loc["M0"].n - e1.loc["M0"].n)
    pos_tr = int(e2.loc["M0"].positives - e1.loc["M0"].positives)
    r4 = e1.loc["M4"]
    mt = ex[ex.what.str.startswith("E1 with MANTIS") & ex.model.eq("M4")].set_index("what").value
    mt_n = ex[ex.what.str.startswith("E1 with MANTIS") & ex.model.eq("M4")].iloc[0]

    def abl(m):
        r = ph[ph.what.eq("ablation: held-out AUROC") & ph.group.eq(m)].iloc[0]
        return f"{float(r.value):.3f} ({float(r.ci_low):.3f}–{float(r.ci_high):.3f})"

    def cls(k, what):
        return float(ph[ph.what.eq(f"M4 standardised coefficient: {what}") & ph.group.eq(k)].value.iloc[0])
    q = pr.loc["LR4"]
    names = {"M1": "two-arm logistic", "M2": "10-bin logistic", "M3": "10-bin MLP", "M5": "all-exon logistic",
             "M6": "all-exon MLP"}
    ladder = ", ".join(f"{names[m]} {e1.loc[m].auroc:.3f}" for m in names)
    c_chosen = r4.chosen.split(":")[1].strip(" }")
    return f"""**Can a model fitted on the repeat exons in three cancers find MMR deficiency in 30 others?**

Pre-specified (config/followup.yaml: msi_code_models, D-71) before any model was fitted; A-19 restricts this figure
to one of the seven pre-registered models, chosen after the results. Label: ≥ 10 MC3 homopolymer deletions. The
model is L2 logistic regression (balanced class weights) on the 1,562 exons with an acceptor T-run ≥ 10 nt that are
measured in ≥ 25 cancers; each input is the exon's inclusion percentile within the tumour's cancer minus the
tumour's mean percentile of matched exons without a run. It was trained on all COAD, STAD and UCEC tumours ({n_tr};
{pos_tr} MMR-deficient), with C chosen by leaving one training cancer out (C = {c_chosen}), and tested on the 30 other
cancers ({int(r4.n):,} tumours, {int(r4.positives)} MMR-deficient, prevalence {r4.prevalence:.4f}). Score G is the
unfitted code.

(a) ROC held out; grey line = chance. (b) Precision against recall held out; dashed = prevalence. (c) AUROC per
held-out cancer with ≥ 3 MMR-deficient tumours (in brackets). (d) The 20 largest standardised coefficients, refitted
on the training cancers (gene, acceptor T-run length, branch-point distance; positive = more inclusion raises the
odds of MMR deficiency).

- Held-out AUROC {r4.auroc:.3f} (95% CI {r4.auroc_lo:.3f}–{r4.auroc_hi:.3f}) against {e1.loc['M0'].auroc:.3f} for score G; AUPRC
  {r4.auprc:.3f} ({r4.auprc_lo:.3f}–{r4.auprc_hi:.3f}) against {e1.loc['M0'].auprc:.3f}; precision {r4['precision_at_recall_0.5']:.2f} at recall 0.5. In the
  leave-one-training-cancer-out loop, {r4.inner_auroc:.3f}.
- Leaving each cancer out of a model trained on the other 32 (exploratory): {e30.loc['M4'].auroc:.3f} on the same 30
  cancers, {e2.loc['M4'].auroc:.3f} over all 33 ({int(e2.loc['M4'].positives)} MMR-deficient).
- With MANTIS MSI-H against MSS as the held-out label ({int(mt_n.positives)} against {int(mt_n.n - mt_n.positives):,}; exploratory): AUROC
  {mt.iloc[0]:.4f}, AUPRC {mt.iloc[1]:.3f}.
- LR4, the pre-registered test of this model (AUROC against the 10-bin logistic model, predicted lower):
  {signed(q.effect, 3)} (95% CI {signed(q.ci_low, 3)} to {signed(q.ci_high, 3)}), q = {fp(q.q)}, {q.verdict}.

Coefficients by class (A-18, exploratory): median {signed(cls('run 10-19 (R_down)', 'median'), 4)} for runs of 10–19 nt ({cls('run 10-19 (R_down)', 'share > 0'):.0%} positive),
{signed(cls('run >=20, bp >=36 nt (R_up)', 'median'), 4)} for runs ≥ 20 nt with the branch point ≥ 36 nt ({cls('run >=20, bp >=36 nt (R_up)', 'share > 0'):.0%} positive). The same model on the
two arms' exons only reached {abl('M4a')}; on the 107 other repeat exons, {abl('M4x')}. The other
pre-registered models, held-out AUROC: {ladder} (tables in code/models/; family LR, 7 tests). Per-tumour
predictions are kept locally (D-12)."""


BENCH_REFS = [("MANTIS", "MANTIS", "MC3 + MSIsensor"), ("MSIsensor", "MSIsensor", "MC3 + MANTIS"),
              ("MC3", "MC3 count", "MANTIS + MSIsensor")]


def U24() -> None:
    d = CODE_DIR / "benchmark"
    pr = pd.read_csv(d / "predictions.csv").set_index("id")
    met = pd.read_csv(d / "metrics.csv")
    op = pd.read_csv(d / "operating_points.csv")
    pc = pd.read_csv(d / "per_cancer.csv")
    disc = pd.read_csv(d / "discordance.csv")
    cov = pd.read_csv(d / "coverage.csv").set_index("method")
    fig = plt.figure(figsize=(DOUBLE - 0.5, 4.0))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.6)
    rows = []
    sty = {"dna": ("s", INK), "M4": ("o", ACCENT)}
    # a. within-cancer AUROC per reference: the judged exome measure and M4; E2 filled, E1 hollow
    ax = fig.add_subplot(gs[0, 0])
    ylab = []
    for i, (ref, name, by) in enumerate(BENCH_REFS):
        for j, key in enumerate(["dna", "M4"]):
            y0 = i * 2.6 + j
            meth = ref if key == "dna" else "M4"
            mk, col = sty[key]
            for split, off, filled in [("E2", -0.17, True), ("E1", 0.17, False)]:
                r = met[met.split.eq(split) & met.reference.eq(f"R_{ref}") & met.method.eq(meth)].iloc[0]
                ax.hlines(y0 + off, r.lo, min(r.hi, 1.0), color=col, linewidth=0.8)
                ax.scatter(r.auroc_within, y0 + off, s=14, marker=mk, color=col if filled else "white",
                           edgecolors=col, linewidths=0.8, zorder=3)
                rows.append({"panel": "a", "reference": f"R_{ref}", "split": split, "method": meth,
                             "auroc_within": r.auroc_within, "lo": r.lo, "hi": r.hi, "positives": r.positives,
                             "negatives": r.negatives})
            ylab.append((y0, name if key == "dna" else M4_LABEL))
    ax.set_yticks([a for a, _ in ylab])
    ax.set_yticklabels([b for _, b in ylab], fontsize=5.4)
    ax.set_ylim(ylab[-1][0] + 0.6, -0.6)
    ax.set_xlim(0.985, 1.0015)
    ax.set_xticks([0.985, 0.99, 0.995, 1.0])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("within-cancer AUROC", fontsize=6)
    # b. per cancer, E2, MANTIS judged by MC3 + MSIsensor
    bx = fig.add_subplot(gs[0, 1])
    sub = pc[pc.split.eq("E2") & pc.reference.eq("R_MANTIS")]
    cancers = sorted(sub.project.unique())
    for i, c in enumerate(cancers):
        for key, dx_ in [("dna", -0.12), ("M4", 0.12)]:
            meth = "MANTIS" if key == "dna" else "M4"
            r = sub[sub.project.eq(c) & sub.method.eq(meth)].iloc[0]
            mk, col = sty[key]
            bx.scatter(i + dx_, r.auroc, s=16, marker=mk, color=col, zorder=3)
            rows.append({"panel": "b", "reference": "R_MANTIS", "split": "E2", "project": c, "method": meth,
                         "auroc": r.auroc, "positives": r.positives, "n": r.n})
    npos = {c: int(sub[sub.project.eq(c)].positives.iloc[0]) for c in cancers}
    bx.set_xticks(range(len(cancers)))
    bx.set_xticklabels([f"{c}\n({npos[c]})" for c in cancers], fontsize=5.6)
    bx.set_xlim(-0.6, len(cancers) - 0.4)
    bx.set_ylim(0.97, 1.002)
    bx.set_ylabel("AUROC within cancer", fontsize=6)
    # c. matched-volume sensitivity and PPV, E2
    cx = fig.add_subplot(gs[1, 0])
    for i, (ref, name, _) in enumerate(BENCH_REFS):
        for j, key in enumerate(["dna", "M4"]):
            meth = ref if key == "dna" else "M4"
            r = op[op.split.eq("E2") & op.reference.eq(f"R_{ref}") & op.method.eq(meth)].iloc[0]
            mk, col = sty[key]
            y0 = i + (j - 0.5) * 0.3
            cx.scatter(r.sensitivity, y0, s=16, marker=mk, color=col, zorder=3)
            cx.scatter(r.ppv, y0, s=16, marker=mk, color="white", edgecolors=col, linewidths=0.8, zorder=3)
            rows.append({"panel": "c", "reference": f"R_{ref}", "split": "E2", "method": meth, "flagged": r.flagged,
                         "sensitivity": r.sensitivity, "ppv": r.ppv})
    flagged = {r_: int(op[op.split.eq("E2") & op.reference.eq("R_" + r_) & op.method.eq(r_)].flagged.iloc[0])
               for r_, _, _ in BENCH_REFS}
    cx.set_yticks(range(len(BENCH_REFS)))
    cx.set_yticklabels([f"{n} ({flagged[r_]} flagged)" for r_, n, _ in BENCH_REFS], fontsize=5.6)
    cx.set_ylim(len(BENCH_REFS) - 0.5, -0.5)
    cx.set_xlim(0.7, 1.01)
    cx.tick_params(axis="y", length=0)
    cx.spines["left"].set_visible(False)
    cx.set_xlabel("sensitivity, PPV", fontsize=6)
    # d. where the DNA measures disagree, E2
    dxx = fig.add_subplot(gs[1, 1])
    groups = [("MC3 MMR-deficient, both callers MSI-H", "both MSI"),
              ("MC3 MMR-deficient, one caller MSI-H", "one MSI"),
              ("MC3 MMR-deficient, neither caller MSI-H", "neither"),
              ("MC3 proficient, neither caller MSI-H", "proficient")]
    e2d = disc[disc.split.eq("E2")].set_index("group")
    for yy, (g_, lab) in enumerate(groups):
        r = e2d.loc[g_]
        dxx.barh(yy, r.share_top_decile_M4, height=0.55, color=ACCENT if yy < 3 else MUTED)
        rows.append({"panel": "d", "split": "E2", "group": g_, "n": r.n, "share_top_decile_M4": r.share_top_decile_M4,
                     "median_pct_M4": r.median_pct_M4})
    dxx.axvline(0.1, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)))
    rows.append({"panel": "d", "group": "reference line: chance (0.1)", "share_top_decile_M4": 0.1})
    dxx.set_yticks(range(len(groups)))
    dxx.set_yticklabels([f"{lab} ({int(e2d.loc[g_].n):,})" for g_, lab in groups], fontsize=5.6)
    dxx.set_ylim(len(groups) - 0.5, -0.5)
    dxx.set_xlim(0, 1)
    dxx.tick_params(axis="y", length=0)
    dxx.spines["left"].set_visible(False)
    dxx.set_xlabel("share in top decile", fontsize=6)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="s", color=INK, linestyle="", markersize=3.5, label="exome measure judged"),
               Line2D([], [], marker="o", color=ACCENT, linestyle="", markersize=3.5, label=M4_LABEL)]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=5.6, frameon=False, bbox_to_anchor=(0.5, -0.03),
               handletextpad=0.3, columnspacing=1.5)
    from _figbase import panel_letters
    panel_letters([(ax, -0.62), (bx, -0.3), (cx, -0.62), (dxx, -0.55)], y=1.06)
    save(fig, OUT, "U24_msi_callers_benchmark", data=pd.DataFrame(rows),
         caption=CAPTIONS["U24"](pr, met, op, disc, cov))


def _cap_u24(pr, met, op, disc, cov):
    def diff(split, ref):
        r = met[met.split.eq(split) & met.reference.eq(f"R_{ref}") & met.method.eq(f"M4 minus {ref}")].iloc[0]
        a = met[met.split.eq(split) & met.reference.eq(f"R_{ref}") & met.method.eq("M4")].iloc[0]
        return (f"{a.auroc_within:.3f}, difference {signed(r.auroc_within, 4)} "
                f"({signed(r.lo, 4)} to {signed(r.hi, 4)})")

    def pos(split, ref):
        return int(met[met.split.eq(split) & met.reference.eq(f"R_{ref}") & met.method.eq("M4")].positives.iloc[0])

    def o_(ref, m):
        r = op[op.split.eq("E2") & op.reference.eq(f"R_{ref}") & op.method.eq(m)].iloc[0]
        return f"{r.sensitivity:.2f} / {r.ppv:.2f}"
    e2d = disc[disc.split.eq("E2")].set_index("group")
    md = lambda g: f"{e2d.loc[g].median_pct_M4:.2f}"
    m5 = pr[pr.verdict.eq("supported")]
    lo = min(float(met[met.method.eq(f"M4 minus {r_}")].lo.min()) for r_, _, _ in BENCH_REFS)
    return f"""**How does the repeat-exon model compare with the exome MSI callers?**

Pre-specified (config/followup.yaml: msi_tool_benchmark, D-72) before any RNA model other than score G was scored
against MANTIS or MSIsensor. The pre-registered tests (BM1–BM6) were on the all-exon model; A-19 shows the
repeat-exon model (Figure U23) instead, so its comparisons here are exploratory. No measure is a gold standard, so
each exome measure is judged against the tumours on which the other two agree (positive: both call MSI; negative:
both call stable; others left out), and the repeat-exon model against the same tumours. Calls: MC3 ≥ 10 homopolymer
deletions; MANTIS > 0.6 (stable < 0.4); MSIsensor > 10 (stable < 4). Tumours with both caller scores
({int(cov.loc['both exome callers'].n):,} of {int(cov.loc['RNA (M0-M5)'].n):,}; RNA scores all). Within-cancer AUROC ranks each positive only
against negatives of its own cancer. E2: each cancer left out of training (all 33); E1: trained on COAD, STAD and
UCEC, tested on the other 30.

(a) Within-cancer AUROC in E2 (filled) and E1 (hollow); bars = 95% bootstrap intervals over tumours within cancer and
reference label, cut at 1. (b) Per cancer with ≥ 3 positives, E2, MANTIS judged by MC3 + MSIsensor (positives in
brackets). (c) Sensitivity (filled) and PPV (hollow) when the model flags as many tumours as the caller calls MSI, by its own score
pooled over cancers (E2; scores come from 33 models). (d) MC3 MMR-deficient tumours grouped by the two callers' calls,
and proficient tumours called stable by both: share in the top decile of the model's score within their cancer (E2;
n in brackets); dashed = chance.

- MANTIS judged: repeat-exon model E2 {diff('E2', 'MANTIS')}; E1 {diff('E1', 'MANTIS')}.
- MSIsensor judged: E2 {diff('E2', 'MSIsensor')}; E1 {diff('E1', 'MSIsensor')}.
- MC3 count judged: E2 {diff('E2', 'MC3')}; E1 {diff('E1', 'MC3')}.

All six lower bounds lie above −0.05 (lowest {signed(lo, 4)}). Positives per reference, E2: {pos('E2', 'MANTIS')}, {pos('E2', 'MSIsensor')},
{pos('E2', 'MC3')}; E1: {pos('E1', 'MANTIS')}, {pos('E1', 'MSIsensor')}, {pos('E1', 'MC3')}. The exome measures score 0.999–1.000 against each other. Panel (c),
sensitivity / PPV: MANTIS {o_('MANTIS', 'MANTIS')} against {o_('MANTIS', 'M4')}; MSIsensor {o_('MSIsensor', 'MSIsensor')} against {o_('MSIsensor', 'M4')}; MC3
{o_('MC3', 'MC3')} against {o_('MC3', 'M4')}. Panel (d), median within-cancer percentile: {md('MC3 MMR-deficient, both callers MSI-H')} (both), {md('MC3 MMR-deficient, one caller MSI-H')} (one),
{md('MC3 MMR-deficient, neither caller MSI-H')} (neither), {md('MC3 proficient, neither caller MSI-H')} (proficient). Pre-registered, all-exon model: {len(m5)} of 6 supported (it ranks
0.001–0.002 below MANTIS and MSIsensor, q ≤ {fp(m5.q.max())}), none against; non-inferior in all six. Per-tumour scores are kept
locally (D-12)."""


# ============================================================ D-73 agreement with the exome measures
AGREE = [("model + all three", ACCENT), ("model + one or two", ACCENT_LIGHT), ("model only", "#b0412c"),
         ("exome only", INK), ("no call", "#d3d1c7")]
TOOLS = [("rna", "model"), ("dna_mc3", "MC3"), ("dna_mantis", "MANTIS"), ("dna_msisensor", "MSIsensor")]


def _agree_class(d: pd.DataFrame) -> pd.Series:
    return pd.Series(np.select([d.rna & d.n_dna.eq(3), d.rna & d.n_dna.between(1, 2), d.rna & d.n_dna.eq(0),
                                ~d.rna & d.n_dna.ge(1)], [a for a, _ in AGREE[:4]], "no call"), index=d.index)


def U25() -> None:
    d = pd.read_parquet(MSI_LOCAL / "rna_dna_calls.parquet")
    d = d[d.both_scored].reset_index(drop=True)
    d["agree"] = _agree_class(d)
    col = dict(AGREE)
    fig = plt.figure(figsize=(DOUBLE - 0.2, 5.4))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.05, 1], hspace=0.42)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.25, 1], wspace=0.45)
    up = top[0].subgridspec(2, 1, height_ratios=[1.6, 1], hspace=0.08)
    rows, local = [], []
    # a. UpSet: every combination of calls with at least one call
    key = d[[c for c, _ in TOOLS]].astype(int).astype(str).agg("".join, axis=1)
    combos = key[key.ne("0000")].value_counts()
    ax = fig.add_subplot(up[0])
    mx = fig.add_subplot(up[1], sharex=ax)
    for i, (k, n) in enumerate(combos.items()):
        sub = d[key.eq(k)]
        cls = sub.agree.iloc[0]
        ax.bar(i, n, width=0.62, color=col[cls])
        ax.text(i, n + 2, str(n), ha="center", va="bottom", fontsize=5.2, color=INK2)
        on = [j for j, ch in enumerate(k) if ch == "1"]
        for j in range(len(TOOLS)):
            mx.scatter(i, j, s=10, color=INK if j in on else GRID, zorder=3)
        if len(on) > 1:
            mx.plot([i, i], [min(on), max(on)], color=INK, linewidth=0.8, zorder=2)
        rows.append({"panel": "a", "combination": "+".join(TOOLS[j][1] for j in on), "class": cls, "n": int(n)})
    rows.append({"panel": "a", "combination": "none", "class": "no call", "n": int(key.eq("0000").sum())})
    ax.set_ylabel("tumours", fontsize=6)
    ax.set_ylim(0, combos.max() * 1.15)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.spines["bottom"].set_visible(False)
    mx.set_yticks(range(len(TOOLS)))
    mx.set_yticklabels([n for _, n in TOOLS], fontsize=5.6)
    mx.set_ylim(len(TOOLS) - 0.5, -0.5)
    mx.set_xlim(-0.6, len(combos) - 0.4)
    mx.tick_params(axis="both", length=0, labelbottom=False)
    for sp in mx.spines.values():
        sp.set_visible(False)
    # b. per cancer: tumours called by anyone, by agreement class
    bx = fig.add_subplot(top[1])
    called = d[d.agree.ne("no call")]
    tab = called.groupby(["project", "agree"]).size().unstack(fill_value=0)
    tot = tab.sum(axis=1)
    keep = tot[tot >= 3].sort_values(ascending=False).index.tolist()
    rest = tab.loc[~tab.index.isin(keep)].sum()
    n_rest = int((~tab.index.isin(keep)).sum())
    tab = tab.loc[keep].copy()
    if n_rest:
        tab.loc[f"{n_rest} others"] = rest
    for yy, (proj, r) in enumerate(tab.iterrows()):
        left = 0
        for cls, c_ in AGREE[:4]:
            v = int(r.get(cls, 0))
            if v:
                bx.barh(yy, v, left=left, height=0.65, color=c_)
                left += v
            rows.append({"panel": "b", "project": proj, "class": cls, "n": v})
    bx.set_yticks(range(len(tab)))
    bx.set_yticklabels(tab.index, fontsize=5.4)
    bx.set_ylim(len(tab) - 0.4, -0.6)
    bx.set_xscale("symlog", linthresh=10)
    bx.set_xlim(0, tab.sum(axis=1).max() * 1.3)
    bx.set_xticks([0, 5, 10, 50, 100])
    bx.set_xticklabels(["0", "5", "10", "50", "100"])
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.set_xlabel("tumours", fontsize=6)
    # c-e. model score against each exome measure; every tumour, coloured by agreement
    bot = outer[1].subgridspec(1, 3, wspace=0.12)
    ro = d[d.agree.eq("model only")]
    cand = {"highest score": ro.M4.idxmax(),
            "MMR and POLE mutations": ro[ro.mut_POLE & ro[[f"mut_{g}" for g in ["MLH1", "MSH2", "MSH6", "PMS2"]]].any(axis=1)].index[0],
            "lowest MLH1": ro.MLH1_pct.idxmin()}
    panels = [("n_hp", "MC3 homopolymer deletions", [10], "symlog"), ("mantis", "MANTIS score", [0.4, 0.6], "linear"),
              ("msisensor", "MSIsensor score", [4, 10], "symlog")]
    axes = []
    for j, (c, lab, cuts, scale) in enumerate(panels):
        a_ = fig.add_subplot(bot[0, j], sharey=axes[0] if axes else None)
        axes.append(a_)
        for cls, c_ in AGREE[::-1]:
            m = d.agree.eq(cls)
            a_.scatter(d.loc[m, c], d.loc[m, "M4"], s=1.5 if cls == "no call" else 5, color=c_, linewidths=0,
                       rasterized=cls == "no call", zorder=1 if cls == "no call" else 3)
        for name, i in cand.items():
            a_.scatter(d.at[i, c], d.at[i, "M4"], s=34, facecolors="none", edgecolors=INK, linewidths=0.7, zorder=4)
            if j == 0:
                local.append({"panel": "c-e", "candidate": name, "project": d.at[i, "project"], "M4": d.at[i, "M4"]})
        for x in cuts:
            a_.axvline(x, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
            rows.append({"panel": "cde"[j], "reference line": f"{lab} cut", "x": x})
        a_.axhline(0, color=MUTED, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
        if scale == "symlog":
            a_.set_xscale("symlog", linthresh=1)
            a_.set_xlim(-0.2, d[c].max() * 1.5)
        a_.set_xlabel(lab, fontsize=6)
        if j:
            a_.tick_params(axis="y", labelleft=False)
        local += [{"panel": "cde"[j], "file_id": f, "project": pj, "x": x, "model_score": y_, "class": k}
                  for f, pj, x, y_, k in zip(d.file_id, d.project, d[c], d.M4, d.agree)]
    rows.append({"panel": "c-e", "reference line": "model call", "y": 0.0})
    spots = {"highest score": (2.2, 13.0), "MMR and POLE mutations": (2.2, 7.5), "lowest MLH1": None}
    for name, i in cand.items():
        xy = (d.at[i, "n_hp"], d.at[i, "M4"])
        if spots[name] is None:
            axes[0].annotate(d.at[i, "project"], xy, xytext=(4, 3), textcoords="offset points", fontsize=5.2, color=INK)
        else:
            axes[0].annotate(d.at[i, "project"], xy, xytext=spots[name], textcoords="data", fontsize=5.2, color=INK,
                             ha="center", arrowprops={"arrowstyle": "-", "color": INK2, "linewidth": 0.5,
                                                      "shrinkA": 1, "shrinkB": 4})
    axes[0].set_ylabel("model score (decision function)", fontsize=6)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", color=c_, linestyle="", markersize=3.5, label=f"{a} ({int(d.agree.eq(a).sum()):,})")
               for a, c_ in AGREE]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=5.6, frameon=False, bbox_to_anchor=(0.5, -0.025),
               handletextpad=0.2, columnspacing=1.1)
    from _figbase import panel_letters
    panel_letters([(ax, -0.16), (bx, -0.35), (axes[0], -0.3), (axes[1], -0.08), (axes[2], -0.08)], y=1.04)
    save(fig, OUT, "U25_msi_model_vs_callers", data=pd.DataFrame(rows), caption=CAPTIONS["U25"](d, combos, cand),
         local_data=(SCFG.base_data_dir / "figure_data", pd.DataFrame(local)))


def _cap_u25(d, combos, cand):
    n = lambda a: int(d.agree.eq(a).sum())
    ro = d[d.agree.eq("model only")]
    ex = d[d.agree.eq("exome only")]
    tri = ex[ex.n_dna.eq(3)]
    c = {k: d.loc[i] for k, i in cand.items()}
    mc3o = ex[ex.n_dna.eq(1) & ex.dna_mc3]
    tri3 = d[d.n_dna.eq(3)]
    pc = round(float(c['lowest MLH1'].MLH1_pct) * 100)
    ordinal = f"{pc}{'st' if pc % 10 == 1 and pc != 11 else 'nd' if pc % 10 == 2 and pc != 12 else 'rd' if pc % 10 == 3 and pc != 13 else 'th'}"
    return f"""**Where do the repeat-exon model and the exome MSI measures agree and disagree?**

Exploratory view of D-73 (config/followup.yaml: msi_rna_only; tests RO1–RO5 in the table). The model is D-71's
repeat-exon logistic regression (A-19), each cancer scored by a model trained on the other 32; it calls MSI when its
decision function is above 0. Exome calls: MC3 ≥ 10 homopolymer deletions, MANTIS > 0.6, MSIsensor > 10. The
{len(d):,} tumours with both caller scores.

(a) Tumours called MSI by each combination of the four (dots below; bars coloured by agreement class); {n('no call'):,}
tumours are called by none. (b) Tumours called by anyone, per cancer with ≥ 3 such tumours, by class. (c–e) The
model's score against MC3 homopolymer deletions, the MANTIS score and the MSIsensor score (deletions and MSIsensor
on a symmetric log axis); dashed = calling thresholds (MANTIS and MSIsensor: stable below the lower, MSI above the
upper). Rings = three model-only tumours named in D-73: the highest score ({c['highest score'].project}), MC3 calls in all four MMR
genes and POLE ({c['MMR and POLE mutations'].project}) and the lowest MLH1 expression ({c['lowest MLH1'].project}).

- Agreement: the model and all three exome measures call {n('model + all three')} tumours; the model and one or two call
  {n('model + one or two')}; the model alone {n('model only')}; exome measures alone {n('exome only')}.
- Model-only calls sit just above the model's threshold (median score {ro.M4.median():.1f}, against
  {d[d.agree.eq('model + all three')].M4.median():.1f} when all agree), and below every exome threshold (median MANTIS {ro.mantis.median():.2f},
  MSIsensor {ro.msisensor.median():.2f}, MC3 deletions {ro.n_hp.median():.0f}).
- Exome-only calls: {len(tri)} by all three measures, {(tri.project == 'UCEC').sum()} of them UCEC; {int(ex.n_dna.eq(1).sum())} by one measure alone,
  {len(mc3o)} of them by MC3 alone. Those are hypermutated (median {mc3o.n_snv.median():,.0f} SNVs, against {d[d.agree.eq('no call')].n_snv.median():.0f} in tumours no one
  calls) with half the deletions per SNV of triple calls ({(mc3o.n_hp / (mc3o.n_snv + 1)).median():.3f} against {(tri3.n_hp / (tri3.n_snv + 1)).median():.3f}), and both
  callers place them in their between zones (median MANTIS {mc3o.mantis.median():.2f}, MSIsensor {mc3o.msisensor.median():.1f}).
- Candidates: {c['highest score'].project} score {c['highest score'].M4:.1f}, {int(c['highest score'].n_hp)} deletions, MANTIS {c['highest score'].mantis:.2f}; {c['MMR and POLE mutations'].project} score {c['MMR and POLE mutations'].M4:.1f}, MSIsensor
  {c['MMR and POLE mutations'].msisensor:.1f}, {int(c['MMR and POLE mutations'].n_snv):,} SNVs; {c['lowest MLH1'].project} score {c['lowest MLH1'].M4:.1f}, MLH1 at the {ordinal} percentile, MSIsensor {c['lowest MLH1'].msisensor:.1f}.

Per-tumour points are kept locally (D-12)."""


# ============================================================ D-74 splicing against expression
EXPR = [("M4", "splicing (repeat exons)", ACCENT), ("X2", "expression, log CPM", INK),
        ("X1", "expression, percentiles", INK), ("X3", "PreMSIm genes, refitted", INK2),
        ("X3k", "PreMSIm recipe (k-NN)", INK2), ("X0", "MLH1 alone", MUTED),
        ("C1", "splicing + expression", ACCENT_DARK)]
EXPR_REFS = [("MANTIS", "MC3 + MSIsensor"), ("MSIsensor", "MC3 + MANTIS"), ("MC3", "MANTIS + MSIsensor")]


def U26() -> None:
    d = MSI_DIR / "code" / "expression"
    met = pd.read_csv(d / "metrics.csv").set_index("model")
    pr = pd.read_csv(d / "predictions.csv")
    ex = pd.read_csv(d / "exploratory.csv")
    p1 = pd.read_parquet(MSI_LOCAL / "expr_e1_predictions.parquet")
    fig = plt.figure(figsize=(DOUBLE - 0.6, 4.5))
    gs = fig.add_gridspec(2, 1, hspace=0.55, height_ratios=[1.0, 1.0])
    top = gs[0].subgridspec(1, 2, wspace=0.08)
    bot = gs[1].subgridspec(1, 3, wspace=1.05, width_ratios=[1.15, 0.95, 0.8])
    rows, local = [], []
    ypos = {m: (i if m != "C1" else i + 0.6) for i, (m, _, _) in enumerate(EXPR)}
    # a, b. held-out AUROC and AUPRC with 95% intervals
    ax = fig.add_subplot(top[0])
    bx = fig.add_subplot(top[1], sharey=ax)
    for a_, key in [(ax, "auroc"), (bx, "auprc")]:
        for m, name, col in EXPR:
            r = met.loc[m]
            a_.hlines(ypos[m], r[f"{key}_lo"], r[f"{key}_hi"], color=col, linewidth=0.9)
            a_.scatter(r[key], ypos[m], s=18, color=col, zorder=3, marker="D" if m == "C1" else "o")
            rows.append({"panel": "a" if key == "auroc" else "b", "model": m, "name": name, "value": r[key],
                         "lo": r[f"{key}_lo"], "hi": r[f"{key}_hi"], "positives": int(r.positives), "n": int(r.n)})
        a_.axvline(met.loc["M4", key], color=ACCENT, linewidth=0.6, linestyle=(0, (3, 2)), zorder=1)
        rows.append({"panel": "a" if key == "auroc" else "b", "model": "reference line: splicing model",
                     "value": met.loc["M4", key]})
        a_.tick_params(axis="y", length=0)
        a_.spines["left"].set_visible(False)
    ax.set_yticks([ypos[m] for m, _, _ in EXPR])
    ax.set_yticklabels([n for _, n, _ in EXPR], fontsize=5.8)
    ax.set_ylim(ypos["C1"] + 0.6, -0.6)
    ax.axvline(0.5, color=GRID, linewidth=0.6, zorder=0)
    ax.set_xlim(0.45, 1.0)
    ax.set_xlabel("held-out AUROC", fontsize=6)
    plt.setp(bx.get_yticklabels(), visible=False)
    bx.axvline(float(p1.mmr.mean()), color=GRID, linewidth=0.6, zorder=0)
    rows.append({"panel": "b", "model": "reference line: prevalence", "value": float(p1.mmr.mean())})
    bx.set_xlim(0, 0.7)
    bx.set_xlabel("held-out AUPRC", fontsize=6)
    # c. MLH1 low vs not low
    cx = fig.add_subplot(bot[0])
    sub = ex[ex.what.str.startswith("AUROC, MMR-deficient with MLH1")]
    order = [m for m, _, _ in EXPR if m != "X3k"]
    for i, m in enumerate(order):
        col = dict((k, c_) for k, _, c_ in EXPR)[m]
        lo_ = sub[sub.model.eq(m) & sub.what.str.contains("MLH1 low")].iloc[0]
        hi_ = sub[sub.model.eq(m) & sub.what.str.contains("not low")].iloc[0]
        cx.plot([lo_.value, hi_.value], [i, i], color=col, linewidth=0.6, zorder=1)
        cx.scatter(lo_.value, i, s=16, color=col, zorder=3)
        cx.scatter(hi_.value, i, s=16, color="white", edgecolors=col, linewidths=0.9, zorder=3)
        for lab, r in [("MLH1 low", lo_), ("MLH1 not low", hi_)]:
            rows.append({"panel": "c", "model": m, "group": lab, "value": r.value, "positives": int(r.n)})
    cx.set_yticks(range(len(order)))
    cx.set_yticklabels([dict((k, n) for k, n, _ in EXPR)[m] for m in order], fontsize=5.6)
    cx.set_ylim(len(order) - 0.5, -0.5)
    cx.axvline(0.5, color=GRID, linewidth=0.6, zorder=0)
    cx.set_xlim(0.4, 1.0)
    cx.tick_params(axis="y", length=0)
    cx.spines["left"].set_visible(False)
    cx.set_xlabel("held-out AUROC", fontsize=6)
    # d. the two read-outs on the held-out MMR-deficient tumours
    dx = fig.add_subplot(bot[1])
    top_ = lambda m: p1.groupby("project")[m].rank(ascending=False, pct=True)     # 1/n for the top tumour
    q = p1.assign(pM4=top_("M4"), pX1=top_("X1"))
    q = q[q.mmr]
    low = q.MLH1_pct < 0.10
    dx.plot([7e-4, 1.3], [7e-4, 1.3], color=GRID, linewidth=0.6, zorder=0)
    dx.scatter(q.pX1[low], q.pM4[low], s=11, color=INK, zorder=3)
    dx.scatter(q.pX1[~low], q.pM4[~low], s=11, color="white", edgecolors=INK, linewidths=0.8, zorder=3)
    for _, r in q.iterrows():
        local.append({"panel": "d", "file_id": r.file_id, "project": r.project, "top_M4": r.pM4, "top_X1": r.pX1,
                      "mlh1_low": bool(r.MLH1_pct < 0.10)})
    rows.append({"panel": "d", "model": "reference line: y = x"})
    dx.set_xscale("log")
    dx.set_yscale("log")
    dx.set_xlim(1.3, 7e-4)
    dx.set_ylim(1.3, 7e-4)
    dx.set_xticks([1e-3, 1e-2, 1e-1, 1])
    dx.set_yticks([1e-3, 1e-2, 1e-1, 1])
    dx.set_xticklabels(["0.001", "0.01", "0.1", "1"])
    dx.set_yticklabels(["0.001", "0.01", "0.1", "1"])
    dx.minorticks_off()
    dx.set_xlabel("expression rank (top fraction)", fontsize=6)
    dx.set_ylabel("splicing rank (top fraction)", fontsize=6)
    # e. within-cancer AUROC against the exome consensus, leaving each cancer out (E2)
    ex_ = fig.add_subplot(bot[2])
    w = ex[ex.what.str.startswith("E2 within-cancer AUROC against the consensus excluding")]
    for i, (ref, lab) in enumerate(EXPR_REFS):
        for m, mk, col, off in [("M4", "o", ACCENT, -0.13), ("X1", "s", INK, 0.13)]:
            r = w[w.what.str.endswith(f"excluding {ref}") & w.model.eq(m)].iloc[0]
            ex_.scatter(r.value, i + off, s=16, marker=mk, color=col, zorder=3)
            rows.append({"panel": "e", "reference": lab, "model": m, "value": r.value, "positives": int(r.n)})
    ex_.set_yticks(range(len(EXPR_REFS)))
    ex_.set_yticklabels([f"{lab}" for _, lab in EXPR_REFS], fontsize=5.6)
    ex_.set_ylim(len(EXPR_REFS) - 0.5, -0.5)
    ex_.set_xlim(0.96, 1.001)
    ex_.set_xticks([0.96, 0.98, 1.0])
    ex_.tick_params(axis="y", length=0)
    ex_.spines["left"].set_visible(False)
    ex_.set_xlabel("within-cancer AUROC", fontsize=6)
    from matplotlib.lines import Line2D
    n_lo, n_hi = int(low.sum()), int((~low).sum())
    handles = [Line2D([], [], marker="o", color=INK, linestyle="", markersize=3.5, label=f"MLH1 low ({n_lo})"),
               Line2D([], [], marker="o", color=INK, markerfacecolor="white", linestyle="", markersize=3.5,
                      label=f"MLH1 not low ({n_hi})"),
               Line2D([], [], marker="o", color=ACCENT, linestyle="", markersize=3.5, label="splicing"),
               Line2D([], [], marker="s", color=INK, linestyle="", markersize=3.5, label="expression, percentiles")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=5.6, frameon=False, bbox_to_anchor=(0.5, -0.04),
               handletextpad=0.3, columnspacing=1.5)
    from _figbase import panel_letters
    panel_letters([(ax, -0.62), (bx, -0.05), (cx, -0.95), (dx, -0.4), (ex_, -1.05)], y=1.05)
    save(fig, OUT, "U26_msi_splicing_vs_expression", data=pd.DataFrame(rows),
         caption=CAPTIONS["U26"](met, pr, ex, q),
         local_data=(SCFG.base_data_dir / "figure_data", pd.DataFrame(local)))


def _cap_u26(met, pr, ex, q):
    a = lambda m, k="auroc": f"{met.loc[m, k]:.3f} ({met.loc[m, k + '_lo']:.3f}–{met.loc[m, k + '_hi']:.3f})"
    t = pr.set_index("id")
    e = lambda i: f"{signed(t.loc[i, 'effect'])} ({signed(t.loc[i, 'ci_low'])} to {signed(t.loc[i, 'ci_high'])}), q = {fp(t.loc[i, 'q'])}"
    sup = int(t.verdict.eq("supported").sum())
    g = lambda w, m: float(ex[ex.what.eq(w) & ex.model.eq(m)].value.iloc[0])
    w2 = "E2 within-cancer AUROC against the consensus excluding"
    wm = [g(f"{w2} {r_}", "M4") for r_, _ in EXPR_REFS]
    wx = [g(f"{w2} {r_}", "X1") for r_, _ in EXPR_REFS]
    mt = "E1 with MANTIS MSI-H vs MSS as the label: AUROC"
    rho = g("Spearman rho with M4 among held-out MMR-deficient tumours (within-cancer percentile)", "X1")
    pp = lambda m: f"{met.loc[m, 'precision_at_recall_0.5']:.2f}"
    e2 = lambda m: g("E2 pooled AUROC (X1 with E1's parameters)", m)
    k_ = lambda m: met.loc[m, "chosen"].replace("=", " = ")
    ch = (f"chosen for both: {k_('X1')}" if k_("X1") == k_("X2") else
          f"chosen: {k_('X1')} on percentiles, {k_('X2')} on log CPM")
    rng = lambda v: f"{min(v):.3f}" if f"{min(v):.3f}" == f"{max(v):.3f}" else f"{min(v):.3f}–{max(v):.3f}"
    return f"""**Does the splicing model transfer better than expression models trained the same way?**

Pre-specified (config/followup.yaml: msi_vs_expression, D-74) before any expression model was scored on a held-out
cancer. Every model was trained on COAD, STAD and UCEC with the same label (MC3 ≥ 10 homopolymer deletions) and tuned
by leaving one training cancer out; each was tested on the 30 other cancers ({int(met.loc['M4', 'n']):,} tumours,
{int(met.loc['M4', 'positives'])} MMR-deficient). All are L2 logistic regressions with balanced class weights except
the PreMSIm recipe (5 nearest neighbours on its 15 genes scaled to [0, 1]) and MLH1 alone (unfitted). The two
whole-transcriptome models keep the top K of 23,640 genes by ANOVA F on the training tumours ({ch}).

(a) Held-out AUROC and (b) AUPRC; bars = 95% bootstrap intervals over tumours within cancer and label; dashed = the
splicing model; grey = chance and prevalence. (c) AUROC of each model for the held-out MMR-deficient tumours whose
MLH1 is low (filled; below the 10th percentile of their cancer) or not (hollow), each against all held-out proficient
tumours. (d) The held-out MMR-deficient tumours: rank of each within its cancer by the splicing score against the
rank by the expression model on percentile inputs (fraction of the cancer at or above it; log axes, best at the top
right); grey line = equality. (e) Leaving each cancer out of training (E2), AUROC
within cancer against the tumours on which the two named exome measures agree (as Figure U24), for the splicing model and the
expression model on percentile inputs (refitted with its E1 parameters, without inner tuning).

- EX1, splicing against expression on percentile inputs: {e('EX1')}, {t.loc['EX1', 'verdict']}.
- EX2, against expression on log CPM: {e('EX2')}, {t.loc['EX2', 'verdict']}.
- EX3, against PreMSIm's genes refitted: {e('EX3')}, {t.loc['EX3', 'verdict']}.
- EX4, against MLH1 alone: {e('EX4')}, {t.loc['EX4', 'verdict']}.
- EX5, MMR-deficient tumours with MLH1 not low: {e('EX5')}, {t.loc['EX5', 'verdict']}.
- EX6, splicing + expression against expression: {e('EX6')}, {t.loc['EX6', 'verdict']}.

{sup} of 6 supported, none against (family EX). Held-out AUROC: splicing {a('M4')}; expression, log CPM {a('X2')};
percentiles {a('X1')}; PreMSIm genes {a('X3')}; PreMSIm recipe {a('X3k')}; MLH1 {a('X0')}; splicing + expression
{a('C1')}. AUPRC {met.loc['M4', 'auprc']:.2f}, {met.loc['X2', 'auprc']:.2f}, {met.loc['X1', 'auprc']:.2f} and {met.loc['C1', 'auprc']:.2f}
for splicing, the two expression models and the combination; precision at recall 0.5: {pp('M4')}, {pp('X2')}, {pp('X1')} and {pp('C1')}.
Against MANTIS MSI-H (17 held out): splicing {g(mt, 'M4'):.3f}, expression {g(mt, 'X2'):.3f} (log CPM) and {g(mt, 'X1'):.3f}
(percentiles). Panel (d), Spearman ρ = {rho:.2f}. Panel (e): splicing {rng(wm)}, expression {rng(wx)};
the exome measures score 0.999–1.000 against each other (Figure U24). E2 pooled AUROC, all 33 cancers: splicing
{e2('M4'):.3f}, expression {e2('X1'):.3f}. Per-tumour points are kept locally (D-12)."""


MAG = [("R_down", "R_down, runs 10–19", INK2), ("R_up", "R_up, runs ≥ 20, bp ≥ 36", ACCENT), ("B", "B, runs ≤ 5", MUTED)]
EXAMPLE_TITLES = {"MRE11": "MRE11 exon 5 · T11", "METTL8": "METTL8 · T14", "BOD1L1": "BOD1L1 · T24", "CCNDBP1": "CCNDBP1 · T30"}


def U27() -> None:
    d = CODE_DIR / "magnitude"
    pr = pd.read_csv(d / "predictions.csv").set_index("id")
    ex = pd.read_csv(d / "exploratory.csv")
    hist = pd.read_csv(d / "shift_histogram.csv")
    bins = pd.read_csv(d / "bins_interior.csv")
    pct = pd.read_csv(d / "percentiles.csv")
    exc = pd.read_csv(d / "tumour_excess.csv")
    pur = pd.read_csv(d / "purity_bins.csv")
    eg = pd.read_csv(d / "example_psi.csv")
    band = pd.read_csv(d / "tails_by_band.csv")
    loc = pd.read_parquet(MSI_LOCAL / "magnitude_tumour_excess.parquet")
    fig = plt.figure(figsize=(DOUBLE - 0.2, 6.2))
    gs = fig.add_gridspec(3, 1, hspace=0.62, height_ratios=[1.0, 0.62, 1.0])
    r1 = gs[0].subgridspec(1, 2, wspace=0.32, width_ratios=[1.05, 1.0])
    r2 = gs[1].subgridspec(1, 4, wspace=0.28)
    r3 = gs[2].subgridspec(1, 4, wspace=0.62, width_ratios=[0.8, 0.8, 1.05, 1.05])
    rows = []
    # a. two-sided tail shares of exon shifts (interior exons, training cancers)
    ax = fig.add_subplot(r1[0])
    tr = hist[hist.group.eq("training")]
    for c, lab, col in MAG:
        h = tr[tr.mclass.eq(c)].sort_values("bin_low")
        n = int(h.n_exons.iloc[0])
        cnt = h["count"].to_numpy()
        lo_edges, hi_edges = h.bin_low.to_numpy(), h.bin_high.to_numpy()
        neg = lo_edges < -1e-9
        x_neg = hi_edges[neg]                                                 # share below each negative upper edge
        y_neg = np.cumsum(cnt[neg]) / n
        pos = hi_edges > 1e-9
        x_pos = lo_edges[pos]
        y_pos = np.cumsum(cnt[pos][::-1])[::-1] / n                          # share >= each positive lower edge
        keep_n, keep_p = x_neg <= -0.02 + 1e-9, x_pos >= 0.02 - 1e-9
        for xs, ys in [(x_neg[keep_n], y_neg[keep_n]), (x_pos[keep_p], y_pos[keep_p])]:
            ys = np.where(ys > 0, ys, np.nan)
            ax.plot(xs, ys, color=col, linewidth=1.0 if c != "B" else 0.9,
                    linestyle="-" if c != "B" else (0, (3, 1.5)), marker="o", markersize=1.8)
        for xs, ys, side in [(x_neg[keep_n], y_neg[keep_n], "shift <= x"), (x_pos[keep_p], y_pos[keep_p], "shift >= x")]:
            for xv, yv in zip(xs, ys):
                rows.append({"panel": "a", "class": c, "n_exons": n, "tail": side, "x": round(float(xv), 3),
                             "share": float(yv)})
    for v in (-0.1, 0.1):
        ax.axvline(v, color=GRID, linewidth=0.6, zorder=0)
    rows.append({"panel": "a", "class": "reference lines: shift -0.10 and +0.10"})
    ax.set_yscale("log")
    ax.set_ylim(2e-3, 1.6)
    ax.set_xlim(-0.42, 0.42)
    ax.set_xticks([-0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4])
    ax.set_xticklabels(["−0.4", "−0.2", "−0.1", "0", "+0.1", "+0.2", "+0.4"])
    ax.set_yticks([0.003, 0.01, 0.03, 0.1, 0.3])
    ax.set_yticklabels(["0.003", "0.01", "0.03", "0.1", "0.3"])
    ax.minorticks_off()
    ax.set_xlabel("exon shift x, PSI (MMR-deficient − proficient)", fontsize=6)
    ax.set_ylabel("share of exons beyond x", fontsize=6)
    from matplotlib.lines import Line2D
    hand = [Line2D([], [], color=col, linewidth=1.0, linestyle="-" if c != "B" else (0, (3, 1.5)),
                   label=f"{lab} ({int(tr[tr.mclass.eq(c)].n_exons.iloc[0]):,})") for c, lab, col in MAG]
    ax.legend(handles=hand, loc="upper left", fontsize=5.2, frameon=False, handlelength=1.6, borderaxespad=0.1)
    # b. share of exons moving >= 0.10 by run-length bin (interior exons)
    bx = fig.add_subplot(r1[1])
    order = ["B (run ≤ 5)", "10–11", "12–13", "14–15", "16–19", "≥ 20, bp < 36 or none", "≥ 20, bp ≥ 36"]
    short = ["≤ 5", "10–11", "12–13", "14–15", "16–19", "≥ 20,\nbp < 36", "≥ 20,\nbp ≥ 36"]
    for i, b in enumerate(order):
        for g, off, fill in [("training", -0.09, True), ("held out", 0.09, False)]:
            r = bins[bins.group.eq(g) & bins.bin.eq(b)]
            if r.empty:
                continue
            r = r.iloc[0]
            for k, col, mk, dx_ in [("n_down10", INK2, "v", -0.17), ("n_up10", ACCENT, "^", 0.17)]:
                x = i + dx_ + off
                sh = r[k] / r.n_exons
                lo_, hi_ = _wilson(int(r[k]), int(r.n_exons))
                bx.vlines(x, lo_, hi_, color=col, linewidth=0.6)
                bx.scatter(x, sh, s=13, marker=mk, color=col if fill else "white", edgecolors=col, linewidths=0.7, zorder=3)
                rows.append({"panel": "b", "bin": b, "group": g, "direction": "drop >= 0.10" if k == "n_down10" else "rise >= 0.10",
                             "n_exons": int(r.n_exons), "k": int(r[k]), "share": sh, "wilson_low": lo_, "wilson_high": hi_})
    bx.set_xticks(range(len(order)))
    ntr = bins[bins.group.eq("training")].set_index("bin").n_exons
    bx.set_xticklabels([f"{s_}\n({int(ntr.get(b, 0)):,})" for s_, b in zip(short, order)], fontsize=5.2)
    bx.tick_params(axis="x", length=0)
    bx.set_xlim(-0.6, len(order) - 0.4)
    bx.set_ylim(0, 0.33)
    bx.set_xlabel("acceptor T-run, nt (training exons)", fontsize=6)
    bx.set_ylabel("share of exons moving ≥ 0.10", fontsize=6)
    hb = [Line2D([], [], marker="v", color=INK2, linestyle="", markersize=3.5, label="drop"),
          Line2D([], [], marker="^", color=ACCENT, linestyle="", markersize=3.5, label="rise"),
          Line2D([], [], marker="o", color=INK, linestyle="", markersize=3.2, label="training"),
          Line2D([], [], marker="o", color=INK, markerfacecolor="white", linestyle="", markersize=3.2, label="held out")]
    bx.legend(handles=hb, loc="upper left", fontsize=5.4, frameon=False, ncol=2, handletextpad=0.3, columnspacing=0.8)
    # c. four exons: PSI in proficient and MMR-deficient tumours of the training cancers
    cxs = []
    for j, gene in enumerate(EXAMPLE_TITLES):
        cx = fig.add_subplot(r2[j])
        cxs.append(cx)
        for status, col, fill in [("proficient", MUTED, True), ("MMR-deficient", INK2 if eg[eg.gene.eq(gene)].mclass.iloc[0] == "R_down" else ACCENT, False)]:
            h = eg[eg.gene.eq(gene) & eg.status.eq(status)].sort_values("bin_low")
            dens = h["count"].to_numpy() / h["count"].sum() / 0.04
            edges = np.r_[h.bin_low.to_numpy(), 1.0]
            if fill:
                cx.stairs(dens, edges, fill=True, color=GRID, zorder=1)
                cx.stairs(dens, edges, color=MUTED, linewidth=0.6, zorder=2)
            else:
                cx.stairs(dens, edges, color=col, linewidth=1.0, zorder=3)
            med = float(h["median"].iloc[0])
            cx.axvline(med, color=col, linewidth=0.7, linestyle=(0, (2, 1.5)), zorder=4)
            for r in h.itertuples():
                rows.append({"panel": "c", "gene": gene, "feature_id": r.feature_id, "status": status,
                             "bin_low": r.bin_low, "bin_high": r.bin_high, "count": int(r.count), "n": int(r.n),
                             "median": med})
        cx.set_xlim(0, 1)
        cx.set_xticks([0, 0.5, 1])
        cx.set_xticklabels(["0", "0.5", "1"])
        cx.set_yticks([])
        cx.spines["left"].set_visible(False)
        cx.set_title(EXAMPLE_TITLES[gene], fontsize=6, pad=2)
        cx.set_xlabel("PSI", fontsize=6)
    hc = [plt.Rectangle((0, 0), 1, 1, facecolor=GRID, edgecolor=MUTED, linewidth=0.6, label="proficient"),
          Line2D([], [], color=INK2, linewidth=1.0, label="MMR-deficient (R_down exon)"),
          Line2D([], [], color=ACCENT, linewidth=1.0, label="MMR-deficient (R_up exon)"),
          Line2D([], [], color=INK2, linewidth=0.7, linestyle=(0, (2, 1.5)), label="median")]
    cxs[0].legend(handles=hc, loc="upper left", fontsize=5.2, frameon=False, ncol=4, bbox_to_anchor=(0.0, 1.42),
                  handlelength=1.4, columnspacing=1.2)
    # d. per tumour: exons beyond 0.10, in excess of the median proficient tumour of the same cancer
    dx = fig.add_subplot(r3[0])
    ex2 = fig.add_subplot(r3[1])
    groups = ["COAD", "STAD", "UCEC", "held out"]
    for a_, measure, col, lab in [(dx, "n_down10_R_down", INK2, "R_down exons dropped ≥ 0.10"),
                                  (ex2, "n_up10_R_up", ACCENT, "R_up exons risen ≥ 0.10")]:
        for i, g in enumerate(groups):
            for status, off, c_ in [("proficient", -0.18, MUTED), ("MMR-deficient", 0.18, col)]:
                r = exc[exc.measure.eq(measure) & exc.group.eq(g) & exc.status.eq(status)].iloc[0]
                x = i + off
                a_.add_patch(plt.Rectangle((x - 0.14, r.q25), 0.28, r.q75 - r.q25, facecolor="white", edgecolor=c_,
                                           linewidth=0.7, zorder=2))
                a_.vlines(x, r.q05, r.q25, color=c_, linewidth=0.6)
                a_.vlines(x, r.q75, r.q95, color=c_, linewidth=0.6)
                a_.scatter(x, r.q50, marker="D", s=9, color=c_, edgecolors="white", linewidths=0.4, zorder=4)
                rows.append({"panel": "d", "measure": measure, "group": g, "status": status, "n": int(r.n),
                             "q05": r.q05, "q25": r.q25, "median": r.q50, "q75": r.q75, "q95": r.q95})
        a_.axhline(0, color=GRID, linewidth=0.6, zorder=0)
        a_.set_xticks(range(len(groups)))
        a_.set_xticklabels(["COAD", "STAD", "UCEC", "held out"], fontsize=5.4, rotation=90)
        a_.tick_params(axis="x", length=0)
        a_.set_xlim(-0.6, len(groups) - 0.4)
        a_.set_title(lab, fontsize=5.8, pad=2)
    dx.set_ylabel("exons per tumour, minus proficient median", fontsize=6)
    rows.append({"panel": "d", "measure": "reference line: 0 (the proficient median of the cancer)"})
    # e. where MMR-deficient tumours fall among proficient ones (interior exons, training cancers)
    fx = fig.add_subplot(r3[2])
    p_ = pct[pct.group.eq("training") & pct.band.eq("interior")]
    for c, lab, col in MAG:
        h = p_[p_.mclass.eq(c)].sort_values("bin")
        edges = np.linspace(0, 1, len(h) + 1)
        fx.stairs(h.density.to_numpy(), edges, color=col, linewidth=1.0 if c != "B" else 0.9,
                  linestyle="-" if c != "B" else (0, (3, 1.5)), zorder=3 if c != "B" else 2)
        for r in h.itertuples():
            rows.append({"panel": "e", "class": c, "percentile_low": edges[r.bin], "percentile_high": edges[r.bin + 1],
                         "count": int(r.count), "density": r.density})
    fx.axhline(1, color=GRID, linewidth=0.6, zorder=0)
    rows.append({"panel": "e", "class": "reference line: density 1 (no difference)"})
    fx.set_xlim(0, 1)
    fx.set_ylim(0.5, 2.3)
    fx.set_xlabel("percentile among proficient tumours", fontsize=6)
    fx.set_ylabel("density, MMR-deficient tumour × exon", fontsize=6)
    # f. arm shifts by purity
    gx = fig.add_subplot(r3[3])
    for arm, col in [("arm_down", INK2), ("arm_up", ACCENT)]:
        for status, fill, off in [("MMR-deficient", True, 0.006), ("proficient", False, -0.006)]:
            q = pur[pur.arm.eq(arm) & pur.status.eq(status)].sort_values("purity_bin")
            x = q.purity_median.to_numpy() + off
            gx.plot(x, q["mean"] * 100, color=col, linewidth=0.6, zorder=2)
            gx.vlines(x, q.ci_low * 100, q.ci_high * 100, color=col, linewidth=0.6, zorder=2)
            gx.scatter(x, q["mean"] * 100, s=13, color=col if fill else "white", edgecolors=col, linewidths=0.7, zorder=3)
            for r in q.itertuples():
                rows.append({"panel": "f", "arm": arm, "status": status, "purity_low": r.purity_low,
                             "purity_high": r.purity_high, "purity_median": r.purity_median, "n": int(r.n),
                             "mean_x100": r.mean * 100, "ci_low_x100": r.ci_low * 100, "ci_high_x100": r.ci_high * 100})
    gx.axhline(0, color=GRID, linewidth=0.6, zorder=0)
    rows.append({"panel": "f", "arm": "reference line: 0"})
    gx.set_xlim(0.25, 0.95)
    gx.set_xlabel("tumour purity (ABSOLUTE)", fontsize=6)
    gx.set_ylabel("mean arm shift, PSI × 100", fontsize=6)
    hf = [Line2D([], [], marker="o", color=INK2, linestyle="-", linewidth=0.6, markersize=3.2, label="R_down"),
          Line2D([], [], marker="o", color=ACCENT, linestyle="-", linewidth=0.6, markersize=3.2, label="R_up"),
          Line2D([], [], marker="o", color=INK, linestyle="", markersize=3.2, label="MMR-def."),
          Line2D([], [], marker="o", color=INK, markerfacecolor="white", linestyle="", markersize=3.2, label="proficient")]
    gx.set_ylim(-1.35, 3.3)
    gx.legend(handles=hf, loc="upper center", fontsize=5.2, frameon=False, handletextpad=0.3, ncol=2, columnspacing=0.8,
              borderaxespad=0.1)
    from _figbase import panel_letters
    panel_letters([(ax, -0.14), (bx, -0.13), (cxs[0], -0.1), (dx, -0.62), (fx, -0.3), (gx, -0.32)], y=1.07)
    save(fig, OUT, "U27_msi_shift_magnitude", data=pd.DataFrame(rows),
         caption=CAPTIONS["U27"](pr, ex, bins, pct, exc, pur, eg, band, pd.read_csv(d / "largest_movers.csv")),
         local_data=(SCFG.base_data_dir / "figure_data", loc.assign(panel="d")))


def _cap_u27(pr, ex, bins, pct, exc, pur, eg, band, mv):
    t = pr
    v = lambda i: f"{signed(t.loc[i, 'effect'])}, p = {fp(t.loc[i, 'p'])}, q = {fp(t.loc[i, 'q'])}, {t.loc[i, 'verdict']}"
    sup = int(t.verdict.eq("supported").sum())
    agn = int(t.verdict.eq("against").sum())
    inter = ex[ex.band.eq("interior (A-20)")]
    sh = lambda w, g, c: float(inter[inter.what.eq(w) & inter.group.eq(g) & inter["class"].eq(c)].value.iloc[0])
    dn, up = "share of exons with shift <= -0.10", "share of exons with shift >= +0.10"
    pc = lambda x: f"{100 * x:.0f}%" if x >= 0.095 else f"{100 * x:.1f}%"
    _m = lambda x: f"{x:.2f}".replace("-", "−")
    mre = mv[mv.role.str.startswith("positive control")].iloc[0]
    drops = "; ".join(f"{g} {pc(sh(dn, g, 'R_down'))} against {pc(sh(dn, g, 'B'))}" for g in ["COAD", "STAD", "UCEC", "held out"])
    rises = "; ".join(f"{g} {pc(sh(up, g, 'R_up'))} against {pc(sh(up, g, 'B'))}" for g in ["COAD", "STAD", "UCEC", "held out"])
    med = lambda c, g: f"{sh('exon shift (PSI), quantile 0.5', g, c):+.4g}".replace("-", "−")  # 4 significant digits: −0.0055, not −0.005
    q05 = lambda c, g: f"{sh('exon shift (PSI), quantile 0.05', g, c):+.2f}".replace("-", "−")
    xm = lambda m, g, s="MMR-deficient": exc[exc.measure.eq(m) & exc.group.eq(g) & exc.status.eq(s)].iloc[0]
    xs = lambda m, g: f"{xm(m, g).q50:+.0f} ({xm(m, g).q25:+.0f} to {xm(m, g).q75:+.0f})".replace("-", "−")
    dens = lambda g, c, b: float(pct[pct.group.eq(g) & pct.band.eq("interior") & pct.mclass.eq(c) & pct.bin.eq(b)].density.iloc[0])
    egm = lambda gene, s: float(eg[eg.gene.eq(gene) & eg.status.eq(s)]["median"].iloc[0])
    ex_ = ", ".join(f"{g} {egm(g, 'proficient'):.2f} → {egm(g, 'MMR-deficient'):.2f}" for g in ["MRE11", "METTL8", "BOD1L1", "CCNDBP1"])
    b = band.set_index(["band", "mclass"])
    lx = lambda bd, c: f"{b.loc[(bd, c), 'low_excess']:+.3f}".replace("-", "−")
    hx = lambda bd, c: f"{b.loc[(bd, c), 'high_excess']:+.3f}".replace("-", "−")
    pm = pur[pur.status.eq("MMR-deficient")]
    rd = pm[pm.arm.eq("arm_down")]["mean"] * 100
    ru = pm[pm.arm.eq("arm_up")]["mean"] * 100
    e2 = lambda g, c: float(ex[ex.what.eq("variance of PSI explained by MMR status (eta2), median") & ex.group.eq(g)
                                & ex["class"].eq(c)].value.iloc[0])
    e2s = ", ".join(f"{100 * e2(g, 'R_down'):.1f}%" for g in ["COAD", "STAD", "UCEC"])
    n_int = {c: int(inter[inter.what.eq(dn) & inter.group.eq("COAD") & inter["class"].eq(c)].n.iloc[0]) for c in ["R_down", "R_up", "B"]}
    return f"""**How large is the MSI splicing shift, and how is it spread over exons and tumours?**

Pre-specified (config/followup.yaml: msi_magnitude, D-75) before any shift was measured in PSI units. All shifts are
within cancer: the MMR-deficient tumours (MC3 ≥ 10 homopolymer deletions; 150 in COAD, STAD and UCEC, 50 in 30
held-out cancers) against the proficient tumours of the same cancer, on cassette exons measured in ≥ 80% of the
cancer's tumours. Panels a, b and e use the exons whose proficient median PSI lies between 0.1 and 0.9 (A-20), where
a shift is not capped by 0 or 1 (COAD: {n_int['R_down']} R_down, {n_int['R_up']} R_up, {n_int['B']:,} background exons).

(a) Share of exons whose shift (median PSI of MMR-deficient minus proficient tumours) lies below x (left) or above x
(right), the three training cancers pooled over exon × cancer; grey lines at ±0.10. (b) Share of exons dropping (▼)
or rising (▲) by at least 0.10, by acceptor T-run length and branch-point (bp) distance; filled = training cancers,
hollow = held-out cancers (an exon's shift there is the median over held-out MMR-deficient tumours, ≥ 20, of its
deviation from the proficient median of each tumour's cancer); bars = 95% Wilson intervals; n = training exons.
(c) PSI of four exons in the proficient (grey) and MMR-deficient tumours (line) of the training cancers; dashed lines
= medians. MRE11 exon 5 is the known MSI target; METTL8 and BOD1L1 have the largest mean shifts among the R_down and
R_up exons measured in all three training cancers; CCNDBP1 is the MSI event of Figure U7. (d) Per tumour, the R_down
exons with PSI at least 0.10 below the proficient median (left) and the R_up exons at least 0.10 above it (right),
minus the count of the cancer's median proficient tumour; boxes = quartiles, whiskers = 5th–95th percentiles.
(e) Where the MMR-deficient tumours of the training cancers fall among the proficient tumours of their cancer:
density of mid-rank percentiles over tumour × exon pairs (A-20); lines as in (a); 1 = no difference. (f) Mean arm shift (a tumour's
mean deviation over its R_down or R_up exons, minus its mean over background exons of the same PSI decile) by purity
fifth of the MMR-deficient tumours, and for proficient tumours in the same bins; bars = 95% bootstrap intervals.

- MG1, R_down arm more negative at higher purity among MMR-deficient tumours (partial Spearman ρ): {v('MG1')}.
- MG2, R_up arm more positive at higher purity: {v('MG2')}.
- MG3, purity × MMR interaction on the R_down arm: {v('MG3')}.
- MG4, larger lower-tail excess in R_down than background exons: {v('MG4')}.
- MG5, larger upper-tail excess in R_up than background exons: {v('MG5')}.

{sup} of 5 supported, {agn} against (family MG). The tail excess of MG4 and MG5 is the share of MMR-deficient
tumours beyond the proficient 5th (95th) percentile minus the share expected if all moved by the exon's median
shift; that expectation lets proficient values pass 0 and 1 (A-20). Lower-tail excess of background exons runs from
{lx('< 0.1', 'B')} (proficient median < 0.1) to {lx('≥ 0.98', 'B')} (≥ 0.98); where the bound is not reached
(0.5–0.9, 0.9–0.98 and ≥ 0.98), R_down exons score {lx('0.5–0.9', 'R_down')}, {lx('0.9–0.98', 'R_down')} and
{lx('≥ 0.98', 'R_down')} against {lx('0.5–0.9', 'B')}, {lx('0.9–0.98', 'B')} and {lx('≥ 0.98', 'B')}; upper tail below
0.1 and 0.1–0.5, R_up {hx('< 0.1', 'R_up')} and {hx('0.1–0.5', 'R_up')} against {hx('< 0.1', 'B')} and
{hx('0.1–0.5', 'B')}. Exons dropping ≥ 0.10, R_down against background: {drops}. Rising ≥ 0.10, R_up against
background: {rises}. Median shift of R_down exons: {med('R_down', 'COAD')}, {med('R_down', 'STAD')},
{med('R_down', 'UCEC')} and {med('R_down', 'held out')} (COAD, STAD, UCEC, held out); 5th percentiles
{q05('R_down', 'COAD')}, {q05('R_down', 'STAD')}, {q05('R_down', 'UCEC')} and {q05('R_down', 'held out')}. Panel (c),
proficient → MMR-deficient medians: {ex_}; {100 * mre.low_obs:.0f}% of the MMR-deficient
tumours fall below the proficient 5th percentile of MRE11 exon 5 (mean over the three cancers). Panel (d), MMR-deficient median (IQR), R_down: COAD
{xs('n_down10_R_down', 'COAD')}, STAD {xs('n_down10_R_down', 'STAD')}, UCEC {xs('n_down10_R_down', 'UCEC')}, held out
{xs('n_down10_R_down', 'held out')}; R_up: {xs('n_up10_R_up', 'COAD')}, {xs('n_up10_R_up', 'STAD')},
{xs('n_up10_R_up', 'UCEC')}, {xs('n_up10_R_up', 'held out')}. Panel (e), lowest 5% of the proficient range: R_down
{dens('training', 'R_down', 0):.2f}, background {dens('training', 'B', 0):.2f}; highest 5%: R_up
{dens('training', 'R_up', 19):.2f}, background {dens('training', 'B', 19):.2f}; held out {dens('held out', 'R_down', 0):.2f}
against {dens('held out', 'B', 0):.2f} and {dens('held out', 'R_up', 19):.2f} against {dens('held out', 'B', 19):.2f}.
Panel (f), MMR-deficient means across the purity fifths: R_down {_m(rd.min())} to {_m(rd.max())}, R_up {ru.min():.2f}
to {ru.max():.2f} (× 100). MMR status explains a median {e2s} of an R_down exon's PSI variance (COAD, STAD, UCEC).
Per-tumour values are kept locally (D-12)."""


CQ_COL = {"canonical": MUTED, "in_frame": ACCENT, "ptc_nmd": "#b0412c", "ptc_escape": "#e6a88a", "other": "white"}
CQ_NAME = {"canonical": "canonical exon", "in_frame": "in-frame insertion", "ptc_nmd": "poison (NMD)",
           "ptc_escape": "stop, escapes NMD", "other": "UTR, non-coding or other"}
BAND_ORDER = ["≤ 0.02", "0.02–0.1", "0.1–0.3", "0.3–0.7", "> 0.7"]


def U28() -> None:
    d = CODE_DIR / "consequence"
    pr = pd.read_csv(d / "predictions.csv").set_index("id")
    comp = pd.read_csv(d / "composition.csv")
    rt = pd.read_csv(d / "r_up_training_exons.csv")
    ctl = pd.read_csv(d / "nmd_control.csv")
    cb = pd.read_csv(d / "coupling_by_baseline.csv")
    gc = pd.read_csv(d / "gene_change_training.csv")
    mv = pd.read_csv(d / "r_up_movers_by_consequence.csv").set_index("consequence")
    fig = plt.figure(figsize=(DOUBLE - 0.2, 4.7))
    gs = fig.add_gridspec(2, 1, hspace=0.55, height_ratios=[1.0, 1.0])
    top = gs[0].subgridspec(1, 2, wspace=0.55, width_ratios=[1.0, 1.0])
    bot = gs[1].subgridspec(1, 3, wspace=0.55)
    rows = []
    from matplotlib.lines import Line2D
    # a. what inclusion does, by exon class (training cancers)
    ax = fig.add_subplot(top[0])
    order = [("R_up", "rising arm (R_up)"), ("B, matched to R_up proficient PSI", "background, matched"),
             ("B", "background (B)"), ("R_down", "falling arm (R_down)")]
    tr = comp[comp.group.eq("training")].copy()
    tr["consequence"] = tr.consequence.replace({"no_stop": "other"})
    tr = tr.groupby(["class", "consequence"], as_index=False)[["n", "share"]].sum()
    for y, (c_, lab) in enumerate(order):
        left = 0.0
        for k in CQ_COL:
            r = tr[tr["class"].eq(c_) & tr.consequence.eq(k)]
            sh = float(r.share.sum())
            ax.barh(y, sh, left=left, height=0.62, color=CQ_COL[k], edgecolor=INK2 if k == "other" else "white",
                    linewidth=0.4)
            rows.append({"panel": "a", "class": c_, "consequence": k, "share": sh, "n": int(r.n.sum())})
            left += sh
    n_ = lambda c_: int(tr[tr["class"].eq(c_)].n.sum())
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{lab}\n({n_(c_):,})" if not c_.startswith("B, matched") else f"{lab}\n(weighted)"
                        for c_, lab in order], fontsize=5.6)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("share of exons", fontsize=6)
    hand = [plt.Rectangle((0, 0), 1, 1, facecolor=CQ_COL[k], edgecolor=INK2 if k == "other" else CQ_COL[k],
                          linewidth=0.4, label=CQ_NAME[k]) for k in CQ_COL]
    ax.legend(handles=hand, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3, fontsize=5.2, frameon=False,
              handlelength=1.0, columnspacing=0.8, handletextpad=0.4)
    # b. R_up exons: shift against baseline inclusion, by consequence
    bx = fig.add_subplot(top[1])
    sub = rt[rt.consequence.isin(["canonical", "in_frame", "ptc_nmd", "ptc_escape"])]
    for k in ["canonical", "ptc_escape", "in_frame", "ptc_nmd"]:
        v = sub[sub.consequence.eq(k)]
        bx.scatter(v.median_proficient, v["shift"], s=9 if k == "ptc_nmd" else 8, color=CQ_COL[k],
                   edgecolors="white", linewidths=0.3, zorder=3 if k in ("ptc_nmd", "in_frame") else 2,
                   marker="o" if k != "in_frame" else "s")
        for r in v.itertuples():
            rows.append({"panel": "b", "cancer": r.group, "symbol": r.symbol, "consequence": k,
                         "proficient_median": r.median_proficient, "shift": r.shift})
    lead = dict(arrowstyle="-", linewidth=0.4, color=INK2, shrinkA=0, shrinkB=1.5)
    offs = {"BOD1L1": (14, 0), "CCNDBP1": (16, 4), "LRR1": (14, 0), "C19orf12": (14, -9), "BRCA1": (26, 14)}
    for g_, (dx_, dy_) in offs.items():          # each gene's largest poison shift
        r = sub[sub.symbol.eq(g_)].sort_values("shift", ascending=False).iloc[0]
        bx.annotate(g_, (r.median_proficient, r["shift"]), xytext=(dx_, dy_), textcoords="offset points",
                    fontsize=5.2, color=INK, va="center", ha="left", arrowprops=lead)
    bx.axhline(0, color=GRID, linewidth=0.6, zorder=0)
    bx.set_xlim(-0.03, 1.03)
    bx.set_ylim(-0.12, 0.42)
    bx.set_xlabel("PSI in proficient tumours (median)", fontsize=6)
    bx.set_ylabel("shift in MMR-deficient tumours, PSI", fontsize=6)
    hb = [Line2D([], [], marker="o" if k != "in_frame" else "s", color=CQ_COL[k], linestyle="", markersize=3.3,
                 label=CQ_NAME[k]) for k in ["ptc_nmd", "in_frame", "ptc_escape", "canonical"]]
    bx.legend(handles=hb, loc="center right", bbox_to_anchor=(1.0, 0.64), fontsize=5.2, frameon=False,
              handletextpad=0.3, borderaxespad=0.1)
    # c. global-NMD control: do background poison exons rise? (training cancers)
    cx = fig.add_subplot(bot[0])
    c_ = ctl[ctl.group.eq("training")]
    series = [("B", "ptc_nmd", CQ_COL["ptc_nmd"], True, -0.2), ("B", "in_frame", ACCENT, True, -0.07),
              ("R_up", "ptc_nmd", CQ_COL["ptc_nmd"], False, 0.07), ("R_up", "in_frame", ACCENT, False, 0.2)]
    for cls_, k, col, fill, off in series:
        v = c_[c_["class"].eq(cls_) & c_.consequence.eq(k)]
        for r in v.itertuples():
            x = BAND_ORDER.index(r.baseline_psi) + off
            cx.vlines(x, r.ci_low, r.ci_high, color=col, linewidth=0.6)
            cx.scatter(x, r.mean_shift, s=12, color=col if fill else "white", edgecolors=col, linewidths=0.7,
                       marker="o" if cls_ == "B" else "D", zorder=3)
            rows.append({"panel": "c", "class": cls_, "consequence": k, "baseline_psi": r.baseline_psi, "n": int(r.n),
                         "mean_shift": r.mean_shift, "ci_low": r.ci_low, "ci_high": r.ci_high})
    cx.axhline(0, color=GRID, linewidth=0.6, zorder=0)
    cx.set_xticks(range(len(BAND_ORDER)))
    cx.set_xticklabels(BAND_ORDER, fontsize=5.4, rotation=30, ha="right")
    cx.tick_params(axis="x", length=0)
    cx.set_xlim(-0.5, len(BAND_ORDER) - 0.5)
    cx.set_xlabel("PSI in proficient tumours", fontsize=6)
    cx.set_ylabel("mean shift, PSI", fontsize=6)
    hc = [Line2D([], [], marker="o", color=INK2, linestyle="", markersize=3.3, label="background"),
          Line2D([], [], marker="D", color=INK2, markerfacecolor="white", linestyle="", markersize=3.3, label="rising arm")]
    cx.legend(handles=hc, loc="upper left", fontsize=5.2, frameon=False, handletextpad=0.3, borderaxespad=0.1)
    # d. coupling of PSI with host-gene mRNA across tumours, by baseline inclusion (background exons)
    dx = fig.add_subplot(bot[1])
    for k, off in [("ptc_nmd", -0.12), ("in_frame", 0.12)]:
        v = cb[cb["class"].eq("B") & cb.consequence.eq(k)]
        for r in v.itertuples():
            x = BAND_ORDER.index(r.baseline_psi) + off
            dx.vlines(x, r.q25, r.q75, color=CQ_COL[k], linewidth=0.7)
            dx.scatter(x, r.median_rho, s=12, color=CQ_COL[k], zorder=3, marker="o" if k == "ptc_nmd" else "s")
            rows.append({"panel": "d", "consequence": k, "baseline_psi": r.baseline_psi, "n": r.n,
                         "median_rho": r.median_rho, "q25": r.q25, "q75": r.q75})
    dx.axhline(0, color=GRID, linewidth=0.6, zorder=0)
    dx.set_xticks(range(len(BAND_ORDER)))
    dx.set_xticklabels(BAND_ORDER, fontsize=5.4, rotation=30, ha="right")
    dx.tick_params(axis="x", length=0)
    dx.set_xlim(-0.5, len(BAND_ORDER) - 0.5)
    dx.set_xlabel("PSI in proficient tumours", fontsize=6)
    dx.set_ylabel("PSI–mRNA partial Spearman ρ", fontsize=6)
    hd = [Line2D([], [], marker="o", color=CQ_COL["ptc_nmd"], linestyle="", markersize=3.3, label="poison"),
          Line2D([], [], marker="s", color=ACCENT, linestyle="", markersize=3.3, label="in-frame")]
    dx.legend(handles=hd, loc="lower left", fontsize=5.2, frameon=False, handletextpad=0.3, borderaxespad=0.1)
    # e. MSI's rise in R_up exons against their gene's mRNA change (training cancers)
    ex_ = fig.add_subplot(bot[2])
    for k in ["in_frame", "ptc_nmd"]:
        v = gc[gc.consequence.eq(k)]
        ex_.scatter(v["shift"], v.log2fc, s=9, color=CQ_COL[k], edgecolors="white", linewidths=0.3, zorder=3,
                    marker="o" if k == "ptc_nmd" else "s")
        for r in v.itertuples():
            rows.append({"panel": "e", "cancer": r.group, "symbol": r.symbol, "consequence": k, "shift": r.shift,
                         "gene_log2fc": r.log2fc})
    for g_, (dx2, dy2) in {"BOD1L1": (-6, -14), "CCNDBP1": (4, 14), "BRCA1": (10, 22), "LRR1": (4, -14)}.items():
        v = gc[gc.symbol.eq(g_)].sort_values("shift", ascending=False).head(1)
        if len(v):
            r = v.iloc[0]
            ex_.annotate(g_, (r["shift"], r.log2fc), xytext=(dx2, dy2), textcoords="offset points", fontsize=5.2,
                         ha="right" if dx2 < 0 else "left", va="center", arrowprops=lead)
    ex_.set_xticks([0, 0.1, 0.2, 0.3])
    ex_.axhline(0, color=GRID, linewidth=0.6, zorder=0)
    ex_.axvline(0, color=GRID, linewidth=0.6, zorder=0)
    rows.append({"panel": "c-e", "symbol": "reference lines: 0"})
    ex_.set_xlabel("exon shift, PSI", fontsize=6)
    ex_.set_ylabel("gene mRNA change, log2", fontsize=6)
    from _figbase import panel_letters
    panel_letters([(ax, -0.42), (bx, -0.2), (cx, -0.33), (dx, -0.4), (ex_, -0.36)], y=1.06)
    save(fig, OUT, "U28_msi_rising_consequence", data=pd.DataFrame(rows),
         caption=CAPTIONS["U28"](pr, comp, rt, ctl, cb, gc, mv))


def _cap_u28(pr, comp, rt, ctl, cb, gc, mv):
    t = pr
    v = lambda i: f"{signed(t.loc[i, 'effect'])}, p = {fp(t.loc[i, 'p'])}, q = {fp(t.loc[i, 'q'])}, {t.loc[i, 'verdict']}"
    sh = lambda c_, k, g="training": float(comp[comp.group.eq(g) & comp["class"].eq(c_) & comp.consequence.eq(k)].share.sum())
    pc = lambda x: f"{100 * x:.0f}%"
    ins = lambda c_: {k: int(comp[comp.group.eq("training") & comp["class"].eq(c_) & comp.consequence.eq(k)].n.sum())
                      for k in ["in_frame", "ptc_nmd", "ptc_escape"]}
    ru = ins("R_up")
    top = rt[rt.consequence.isin(["canonical", "in_frame", "ptc_nmd", "ptc_escape", "other"])].sort_values(
        "shift", ascending=False).head(15)
    n_p = int(top.consequence.eq("ptc_nmd").sum())
    top_b = rt[rt.consequence.isin(["canonical", "in_frame", "ptc_nmd", "ptc_escape"])].sort_values(   # as plotted in b
        "shift", ascending=False).head(15)
    n_pb = int(top_b.consequence.eq("ptc_nmd").sum())
    n_eb = int(top_b[top_b.consequence.eq("ptc_nmd")].feature_id.nunique())
    inf_max = float(rt[rt.consequence.eq("in_frame")]["shift"].max())
    c_ = ctl[ctl.group.eq("training") & ctl["class"].eq("B") & ctl.consequence.eq("ptc_nmd")]
    c_h = ctl[ctl.group.eq("held out") & ctl["class"].eq("B") & ctl.consequence.eq("ptc_nmd")]
    rng_ = lambda s: f"{s.mean_shift.min():+.4f} to {s.mean_shift.max():+.4f}".replace("-", "−")
    diff = cb[cb.consequence.eq("ptc_nmd minus in_frame")].set_index("baseline_psi")
    dd = lambda b: f"{diff.loc[b, 'median_rho']:+.3f} (p = {fp(diff.loc[b, 'p'])})".replace("-", "−")
    g = lambda s: gc[gc.symbol.eq(s)]
    gl = lambda s: ", ".join(f"{r.group} {r.log2fc:+.2f}".replace("-", "−") for r in g(s).sort_values("group").itertuples())
    return f"""**What does including a rising-arm exon do to its transcript?**

Pre-specified (config/followup.yaml: msi_rising_consequence, D-76) before any consequence annotation or host-gene
expression was examined. Each cassette exon was placed in its gene's canonical protein-coding transcript (GENCODE
v45): a canonical exon, an exon that inclusion inserts into the coding sequence (its flanking exons are joined
directly inside the CDS), or other (UTR, non-coding or non-canonical context). Inserted exons were translated in the
canonical frame: an in-frame insertion, a premature stop more than 55 nt before the last exon–exon junction
(poison, predicted NMD), or a stop that escapes NMD.

(a) Consequence of inclusion for the exons of each class measured in a training cancer (COAD, STAD, UCEC);
background exons are shown as they are and reweighted to the rising arm's distribution of proficient PSI (deciles);
"no stop" (< 0.1%) counted with other. (b) Rising-arm exons in the training cancers (exon × cancer): shift (median PSI
of MMR-deficient minus proficient tumours) against the proficient median, by consequence; the five genes with the
largest poison-exon shifts labelled. (c) Training cancers: mean shift by proficient PSI band, for poison (red) and
in-frame (blue) exons, background (filled) and rising arm (hollow); bars = 95% bootstrap intervals; bands with at least
5 exon × cancer rows. (d) Background exons: per-exon
partial Spearman of PSI with the host gene's log2 CPM across a cancer's tumours, adjusting for the tumour's mean PSI,
averaged over the 33 cancers; median and quartiles by proficient PSI band (A-21). (e) Rising-arm exons inserted
into the CDS, training cancers, one per gene and cancer: shift against the host gene's change (mean log2 CPM of
MMR-deficient minus proficient tumours).

- CQ1, background poison exons couple more negatively than in-frame ones (positive control): {v('CQ1')}.
- CQ2, the same for rising-arm exons: {v('CQ2')}.
- CQ3, MSI lowers the host gene more for rising poison than in-frame exons (AUROC − 0.5): {v('CQ3')}.
- CQ4, Spearman(shift, gene change), poison minus in-frame: {v('CQ4')}.

None of 4 supported (family CQ). The positive control failed, so CQ2–CQ4 cannot be read. Coupling follows baseline
inclusion far more than consequence (panel d; ρ = −0.45 across all annotated exons); within bands, poison background
exons couple slightly more negatively than in-frame ones: {dd('≤ 0.02')}, {dd('0.02–0.1')}, {dd('0.1–0.3')}
(exploratory, A-21). Panel (a), rising arm: {pc(sh('R_up', 'ptc_nmd'))} poison, {pc(sh('R_up', 'in_frame'))} in-frame,
{pc(sh('R_up', 'ptc_escape'))} escaping NMD, {pc(sh('R_up', 'canonical'))} canonical ({ru['ptc_nmd']}, {ru['in_frame']} and
{ru['ptc_escape']} exons inserted into the CDS); matched background {pc(sh('B, matched to R_up proficient PSI', 'ptc_nmd'))},
{pc(sh('B, matched to R_up proficient PSI', 'in_frame'))}, {pc(sh('B, matched to R_up proficient PSI', 'ptc_escape'))} and
{pc(sh('B, matched to R_up proficient PSI', 'canonical'))}. Panel (b): {n_pb} of the 15 largest shifts plotted (exon × cancer rows) are
rows of {n_eb} poison exons ({n_p} of the 15 largest rows when UTR, non-coding and other exons are counted); no in-frame exon moves by more than {inf_max:.3f} ({int(mv.loc['ptc_nmd', 'moving_ge_005'])} of {int(mv.loc['ptc_nmd', 'n'])}
poison and {int(mv.loc['in_frame', 'moving_ge_005'])} of {int(mv.loc['in_frame', 'n'])} in-frame rows move by ≥ 0.05;
Fisher p = {fp(mv.loc['ptc_nmd', 'fisher_p_ptc_nmd_vs_in_frame'])}; exploratory). Panel (c): background poison exons barely move
in MMR-deficient tumours (band means {rng_(c_)}; held out {rng_(c_h)}). Panel (e), host-gene change, BOD1L1:
{gl('BOD1L1')}; CCNDBP1: {gl('CCNDBP1')}; BRCA1: {gl('BRCA1')}. No in-frame rising-arm exon reaches a median PSI of
0.10 in MMR-deficient tumours while ≤ 0.05 in proficient and normal tissue."""


TITLES = {
    "U1": "Splicing matches expression on lineage and genome-state subtypes, and fusing the two adds little.",
    "U2": "Splicing loses most on small classes defined by one programme.",
    "U3": "Histology models err in the molecular direction.",
    "U4": "The literature's markers are in the data and in the expression models, far less in the splicing models.",
    "U5": "The splicing and expression models of a cohort read different genes.",
    "U6": "Four of the nineteen strongest splicing events read a neighbouring gene.",
    "U7": "Ten MSI splicing events replicate across endometrial, gastric and colorectal cancer.",
    "U8": "Named splicing regulators differ between subtypes as predicted in four of six tests.",
    "U9": "Part of the splicing models' subtype signal is tied to the site that sent the sample.",
    "U10": "Genome sequence alone says which exons MSI tumours skip.",
    "U11": "The longest acceptor repeats move toward inclusion in MSI.",
    "U12": "Most of the BRCA models' genotype-like events are subtype programmes, not inherited alleles.",
    "U13": "Splicing reads HER2 through the HER2-enriched programme, not through the 17q12 amplicon.",
    "U14": "Splice-site mutations skip their exons, and SF3B1 hotspot tumours use cryptic acceptors.",
    "U15": "Without the ABSOLUTE ceiling, splicing matches expression on ERBB2 amplification but does not beat it.",
    "U16": "Mismatch-repair damage tracks the repeat score in stomach cancer only.",
    "U17": "The predicted branch point's distance from the 3′ splice site goes with the direction of the MSI shift.",
    "U18": "A two-arm sequence code finds mismatch-repair-deficient tumours in 30 cancers never used to define it.",
    "U19": "The code is specific to the branch point and to mismatch repair.",
    "U20": "The code is not an infiltration readout.",
    "U21": "The code finds the tumours that two exome MSI callers call MSI-high.",
    "U22": "The rising arm is specific to mismatch repair and needs the tract, but its mechanism is not settled.",
    "U23": "Logistic regression on repeat exons, fitted in three cancers, finds MMR deficiency in 30 others.",
    "U24": "Within each cancer, the repeat-exon model ranks tumours almost as well as the exome MSI callers.",
    "U25": "The model agrees with the exome measures on its confident calls; its own extra calls are borderline.",
    "U26": "Expression models trained the same way transfer almost as well; splicing adds to them.",
    "U27": "A minority of exons move by a tenth of PSI or more; the typical exon moves by a few hundredths.",
    "U28": "Most of the rising-arm exons that move furthest are poison exons.",
}


def legend(key: str, text: str) -> str:
    """Caption -> figure legend ('Figure Un. Title.', panels (a), (b) ...); see _figbase.legend."""
    return _legend(key, TITLES[key], text)


def legends() -> None:
    """Copy each figure's legend (its caption without the bullet list of tests) under the figure in
    results/SUBTYPE.md, between marker comments, so the write-up and the caption files never drift."""
    inject_legends(SCFG.results_dir / "SUBTYPE.md")


CAPTIONS = {k: (lambda k, f: (lambda *a: tidy(legend(k, f(*a)))))(k, f) for k, f in
            {"U1": _cap_u1, "U2": _cap_u2, "U3": _cap_u3, "U4": _cap_u4, "U5": _cap_u5, "U6": _cap_u6,
             "U7": _cap_u7, "U8": _cap_u8, "U9": _cap_u9,
             "U10": _cap_u10, "U11": _cap_u11, "U12": _cap_u12, "U13": _cap_u13, "U14": _cap_u14, "U15": _cap_u15,
             "U16": _cap_u16, "U17": _cap_u17, "U18": _cap_u18, "U19": _cap_u19, "U20": _cap_u20,
             "U21": _cap_u21, "U22": _cap_u22, "U23": _cap_u23, "U24": _cap_u24, "U25": _cap_u25,
             "U26": _cap_u26, "U27": _cap_u27, "U28": _cap_u28}.items()}

if __name__ == "__main__":
    run({"U1": U1, "U2": U2, "U3": U3, "U4": U4, "U5": U5, "U6": U6, "U7": U7, "U8": U8, "U9": U9, "U10": U10, "U11": U11,
         "U12": U12, "U13": U13, "U14": U14, "U15": U15, "U16": U16, "U17": U17,
         "U18": U18, "U19": U19, "U20": U20, "U21": U21, "U22": U22, "U23": U23, "U24": U24, "U25": U25, "U26": U26, "U27": U27, "U28": U28,
         "legends": legends})
