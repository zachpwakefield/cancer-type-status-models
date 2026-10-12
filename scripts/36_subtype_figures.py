"""Step 36: the subtype narrative, in figures (house style, docs/FIGURE_GUIDELINES.md).

One question per figure, in the order the story is told in results/SUBTYPE.md:
    U1  how well does each modality recover the subtype, cohort by cohort, and does fusing help?
    U2  which subtypes are recovered, and by which modality?
    U3  when a histology model is wrong, is the error molecularly coherent?
    U4  does the biology the literature predicts show up, in the data and in the models?
    U5  do the splicing models and the expression model of a cohort read the same genes?
    U6  are the strongest splicing events splicing, or readouts of a neighbouring gene?
    U7  do the exons that mark MSI in one cancer mark it in the other two?

Reads results/subtype/summary/ (step 34) and results/subtype/interpret/ (step 35). Every figure
is SVG + PDF + 400-dpi PNG with `<stem>.csv` of the plotted values and `<stem>.caption.md`.
Usage: python scripts/36_subtype_figures.py [U1 U3 ...]        (default: all)
"""
from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _figbase import locus_label, run

from ctsm.biology import bh
from ctsm.figstyle import (ACCENT, ACCENT_DARK, DIVERGE, DOUBLE, GRID, INK, INK2, MUTED, SEQ, SINGLE, colourbar,
                           fp, heat_axes, save)
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

Right: fused minus GEX balanced accuracy, and the exact McNemar test of fused against GEX on paired per-patient
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
1:1 (dark blue). Left: histology cohorts; right: genomic and IHC cohorts, then the expression-defined ones, in the
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
filled and named = the class the event ranks for. Right: the D-27 audit. The event's segment is intersected with
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
value in each cell (MSI and rest) are in the CSV. Right: the D-27 audit, as in U6. Exploratory (A-4): the step was added after the
same CCNDBP1 exon was seen in the MSI top 5 of all three models; it then tests every event of every shared class.

{int((chosen >= 2).sum())} of the {len(k)} events were chosen independently by at least two of the three models.
CCNDBP1's 57-bp exon is included at a median {cs} (MSI vs rest). All {len(nrep)} replication cells have Mann-Whitney
raw p < {max(nrep.p_mwu.max(), 1e-300):.1e}; BH q over every cross-cohort test of every shared class is in the CSV."""


def tidy(caption: str, width: int = 116) -> str:
    """Re-wrap every paragraph and bullet of a caption to the house limit (lines <= 120 characters)."""
    out = []
    for para in caption.strip().split("\n\n"):
        lines = para.split("\n")
        if para.startswith("**") and len(para) <= 120:
            out.append(para)
        elif lines[0].startswith("- "):
            items, cur = [], []
            for ln in lines:
                if ln.startswith("- ") and cur:
                    items.append(" ".join(cur))
                    cur = []
                cur.append(ln.strip())
            items.append(" ".join(cur))
            out.append("\n".join(textwrap.fill(it, width, subsequent_indent="  ", break_on_hyphens=False)
                                  for it in items))
        else:
            out.append(textwrap.fill(" ".join(ln.strip() for ln in lines), width, break_on_hyphens=False))
    return "\n\n".join(out)


CAPTIONS = {k: (lambda f: (lambda *a: tidy(f(*a))))(f) for k, f in
            {"U1": _cap_u1, "U2": _cap_u2, "U3": _cap_u3, "U4": _cap_u4, "U5": _cap_u5, "U6": _cap_u6,
             "U7": _cap_u7}.items()}

if __name__ == "__main__":
    run({"U1": U1, "U2": U2, "U3": U3, "U4": U4, "U5": U5, "U6": U6, "U7": U7})
