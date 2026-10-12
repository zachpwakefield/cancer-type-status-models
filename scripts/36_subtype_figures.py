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

Reads results/subtype/summary/ (step 34), results/subtype/interpret/ (step 35),
results/subtype/regulators/ (step 37), results/subtype/site_probe/ (step 38) and
results/subtype/msi_cis/ (step 42). Every figure
is SVG + PDF + 400-dpi PNG with `<stem>.csv` of the plotted values and `<stem>.caption.md`.
Usage: python scripts/36_subtype_figures.py [U1 U3 ...]        (default: all)
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _figbase import locus_label, run, tidy

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

CAPTIONS = {k: (lambda f: (lambda *a: tidy(f(*a))))(f) for k, f in
            {"U1": _cap_u1, "U2": _cap_u2, "U3": _cap_u3, "U4": _cap_u4, "U5": _cap_u5, "U6": _cap_u6,
             "U7": _cap_u7, "U8": _cap_u8, "U9": _cap_u9,
             "U10": _cap_u10, "U11": _cap_u11}.items()}

if __name__ == "__main__":
    run({"U1": U1, "U2": U2, "U3": U3, "U4": U4, "U5": U5, "U6": U6, "U7": U7, "U8": U8, "U9": U9, "U10": U10, "U11": U11})
