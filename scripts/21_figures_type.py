"""Step 21: the cancer-type narrative, in figures (house style, docs/FIGURE_GUIDELINES.md).

One question per figure, in the order the story is told in results/CANCER_TYPE.md:
    T1  does combining the modalities identify the tissue of origin better than expression alone?
    T2  which classes does splicing resolve better than expression, and which worse?
    T3  which classes do the two modalities confuse?
    T4  do the views agree on the processes behind a class?
    T5  are the class-defining events just a readout of their gene's expression?
    T6  what does a lineage-defining exon look like across the classes?
    T7  how much of the separation of neighbouring classes is the tissue-source site?

Usage: python scripts/21_figures_type.py [T1 T5 ...]        (default: all)
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _figbase import AS_VIEWS, BIO, CANCERS, CFG, FIGURES, LOCAL, locus_label, log, run, short
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import leaves_list, linkage

from ctsm.biology import bio_config, feature_values, task_table
from ctsm.config import PAN, load_config
from ctsm.figstyle import (ACCENT, ACCENT_DARK, DIVERGE, DOUBLE, GRID, INK, INK2, MUTED, NORMAL, SEQ, SINGLE, TUMOUR,
                           VIEW_NAMES, bracket, colourbar, fp, gene_label, heat_axes, save, stars)

OUT = FIGURES / "cancer_type"


# ============================================================ T1 combining modalities
def T1() -> None:
    d = BIO / "09_ensemble"
    m = pd.read_csv(d / "cancer_type_model_metrics.csv").set_index("model")
    t = pd.read_csv(d / "cancer_type_model_tests.csv").set_index("model")
    order = [VIEW_NAMES[v] for v in CFG.views] + ["AS ensemble", "GEX + AS", "stacked"]
    col = {k: MUTED for k in order}
    col.update({"GEX": INK, "AS ensemble": ACCENT, "GEX + AS": ACCENT_DARK, "stacked": MUTED})
    p = m.loc[order]
    x = np.arange(len(order)) + np.r_[np.zeros(len(CFG.views)), np.ones(3) * 0.6]
    fig, ax = plt.subplots(figsize=(SINGLE, 2.4))
    for xi, (name, r) in zip(x, p.iterrows()):
        ax.vlines(xi, r.accuracy_lo, r.accuracy_hi, color=col[name], linewidth=0.8)
        ax.scatter(xi, r.accuracy, s=12, color=col[name], linewidths=0, zorder=3,
                   marker="o" if name != "stacked" else "D")
        ax.text(xi, r.accuracy_hi + 0.0012, f"{r.accuracy:.3f}", fontsize=4.8, color=INK2, ha="center", va="bottom")
    xpos = dict(zip(order, x))
    top = p.accuracy_hi.max()
    for k, name in enumerate(["AS ensemble", "GEX + AS"]):
        bracket(ax, xpos["GEX"], xpos[name], top + 0.008 + k * 0.010,
                f"p = {fp(t.loc[name, 'p_mcnemar'])}", 0.002, fontsize=5.5)
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha="right", color=INK)
    ax.set_ylabel("accuracy, 34 classes (out-of-fold)")
    ax.set_ylim(p.accuracy_lo.min() - 0.005, top + 0.034)
    gain = int(t.loc["GEX + AS", "only_b_correct"]), int(t.loc["GEX + AS", "only_a_correct"])
    save(fig, OUT, "T1_modality_gain", data=m.reset_index().merge(t.reset_index(), on="model", how="left"),
         caption=f"""
**Does combining the modalities identify the tissue of origin better than expression alone?**

Out-of-fold accuracy over {int(p.n.iloc[0]):,} samples and 34 classes (33 TCGA projects and a pooled normal class)
for each single-view model (grey; GEX in black), the AS ensemble (the mean of the 8 AS views' probabilities, blue),
an equal-weight fusion of GEX and the AS ensemble, and a cross-fitted stacker on all 9 views (diamond). The fusion
weights were fixed before any result was seen (D-14); the stacker is a ceiling, not a headline (D-23).

Whiskers: 95% case-level bootstrap CIs (1,000 resamples). Brackets: exact McNemar on paired per-sample
correctness, raw p (two tests; q in the CSV).

No single splicing view reaches expression, but their ensemble does. Fused with expression it corrects
{gain[0]} of its errors and introduces {gain[1]}. Macro recall and macro F1 are in the CSV.
""")
def T2() -> None:
    e = BIO / "09_ensemble"
    rc = pd.read_csv(e / "cancer_type_recall_by_class.csv")
    piv = rc.pivot(index="label", columns="model", values="recall")
    n = rc.drop_duplicates("label").set_index("label")["n"]
    piv["as_minus_gex"] = piv["AS ensemble"] - piv["GEX"]
    piv["fused_minus_gex"] = piv["GEX + AS"] - piv["GEX"]
    piv = piv.sort_values("as_minus_gex")
    fig, ax = plt.subplots(figsize=(2.9, 4.2))
    y = np.arange(len(piv))
    for i, (_, r) in enumerate(piv.iterrows()):
        lo, hi = sorted([r.as_minus_gex, r.fused_minus_gex])
        ax.plot([lo, hi], [i, i], color=GRID, linewidth=1.2, zorder=1)
    ax.scatter(piv["as_minus_gex"], y, s=10, color=ACCENT, linewidths=0.3, edgecolors="white",
               zorder=3, label="AS ensemble")
    ax.scatter(piv["fused_minus_gex"], y, s=10, color=ACCENT_DARK, linewidths=0.3, edgecolors="white",
               zorder=4, label="GEX + AS")
    ax.axvline(0, color=INK2, linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{k} ({int(n[k])})" for k in piv.index], fontsize=5.5, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.7, len(piv) - 0.3)
    ax.set_xlabel("recall − expression model's recall")
    ax.legend(loc="lower right", fontsize=5.5, handletextpad=0.3, borderpad=0.2)
    gained = int((piv["fused_minus_gex"] > 0).sum())
    lost = int((piv["fused_minus_gex"] < 0).sum())
    down = piv[piv["as_minus_gex"] < 0]
    up = piv[piv["as_minus_gex"] > 0]
    sgn = lambda v: f"{v:.3f}".replace("-", "\u2212")      # typographic minus, as in the prose
    save(fig, OUT, "T2_recall_by_class", data=piv.reset_index(), caption=f"""
**Which classes does splicing resolve better than expression, and does fusing them help?**

Out-of-fold recall of each of the 34 classes (n in brackets), as a difference from the expression model: the AS
ensemble (light blue) and the equal-weight fusion of the two (dark blue). Positive means the model recovers more
of that class than expression does.

Fusion raises recall for {gained} classes and lowers it for {lost}. It is asymmetric in a useful way: across the
{len(down)} classes where splicing alone is worse than expression, the median shortfall shrinks from
{sgn(down.as_minus_gex.median())} to {sgn(down.fused_minus_gex.median())}, while across the {len(up)} classes where
splicing is better it keeps most of the gain ({sgn(up.as_minus_gex.median())} to {sgn(up.fused_minus_gex.median())}).
The pooled normal class is the clearest case: −0.102 for splicing alone, −0.007 fused. So the overall gain (T1) is
not bought at the cost of particular classes.

Splicing gains most for classes expression confuses with an anatomical neighbour (COAD/READ, ESCA/STAD, GBM/LGG).
Of those, only the glioma separation survives holding out whole tissue-source sites (T7): the colon/rectum gain
does not and should not be read as biology.

Expression is better for the pooled normal class, uterine carcinosarcoma, and squamous and urothelial tumours
(BLCA, HNSC, CESC). Every model's per-class recall is in the CSV.
""")


# ============================================================ 02 global shifts
SHIFT_ORDER = ["ri", "a5ss", "a3ss", "se", "mxe", "ale", "afe"]
SHIFT_LAB = ["RI retained", "A5SS long form", "A3SS long form", "SE included", "MXE upstream exon",
             "ALE downstream", "AFE downstream"]


def T3() -> None:
    d = BIO / "06_lineage"
    conf = pd.read_csv(d / "class_pair_confusions.csv").head(16)
    fig, ax = plt.subplots(figsize=(2.6, 2.6))
    y = np.arange(len(conf))
    ax.barh(y - 0.19, conf["GEX"], height=0.36, color=INK2, label="GEX")
    ax.barh(y + 0.19, conf["AS ensemble"], height=0.36, color=ACCENT, label="AS ensemble")
    ax.set_yticks(y)
    ax.set_yticklabels(conf["pair"], fontsize=5.5, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.invert_yaxis()
    ax.set_xlabel("misclassified samples within the pair\n(both directions, out-of-fold)")
    ax.legend(loc="lower right")
    save(fig, OUT, "T3_class_confusions", data=conf, caption="""
**Do expression and splicing confuse the same cancer types?**

The 16 class pairs with the most out-of-fold confusions (either direction) for the GEX model or the AS ensemble.

Expression mixes up anatomical neighbours: colon and rectum, oesophagus and stomach, lung adenocarcinoma and squamous
carcinoma, the kidney subtypes, and glioma grades. Splicing makes fewer of these errors; GBM/LGG errors fall from 20
to 0. T7 shows which of these gains survive unseen tissue-source sites: glioma, kidney and lung yes; colon/rectum
no.

Splicing instead confuses squamous and urothelial carcinomas across organs (CESC/HNSC, HNSC/LUSC, BLCA/CESC,
BLCA/HNSC). It also confuses several tumours with normal tissue (PRAD, THCA, LIHC, STAD), as if the splicing state
follows the tissue and the cell type more than the tumour.

""")


def T4() -> None:
    d = BIO / "03_pathways"
    coh = pd.read_csv(d / "pathway_coherence_by_view.csv", index_col=0).reindex(CFG.views)
    c_ts, c_ct, c_rp = coh.columns
    fig, axes = plt.subplots(1, 2, figsize=(3.8, 2.0), sharey=True)
    y = np.arange(len(CFG.views))
    cols = [INK2 if v == "gex" else ACCENT for v in CFG.views]
    # open = all genes, filled = the ribosomal-protein sensitivity, which was only run on the
    # cancer-type AS views; the left panel and GEX have no filled bar rather than a zero one.
    axes[0].barh(y, coh[c_ts], color="none", edgecolor=cols, linewidth=0.6, height=0.65)
    axes[0].set_xlabel("tumour vs normal:\nshare of 15 cancers")
    axes[1].barh(y, coh[c_ct], color="none", edgecolor=cols, linewidth=0.6, height=0.65, label="all genes")
    axes[1].barh(y, coh[c_rp], color=ACCENT, height=0.65, label="ribosomal-protein\ngenes removed")
    axes[1].set_xlabel("cancer type:\nshare of 34 classes")
    axes[1].legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, handlelength=1.0, columnspacing=0.8)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([VIEW_NAMES[v] for v in CFG.views], color=INK)
    axes[0].invert_yaxis()
    save(fig, OUT, "T4_pathway_coherence", data=coh.reset_index().rename(columns={"index": "view"}).assign(
        q_family="BH within collection per model (Hallmark, Reactome, GO-BP)"), caption="""
**Do the genes a splicing model relies on form pathways, as an expression model's do?**

The share of models whose top 200 SHAP genes are enriched (q < 0.05, BH within collection per model; GOseq/Wallenius,
D-16) for at least one Hallmark, Reactome or GO-BP term, by view.

Left: the 15 tumour-vs-normal models. GEX reaches a pathway in 14 of 15 cancers; SE in 9; MXE in 6; the other AS
views in at most 4; RI in none. The tumour signal of splicing is mostly carried by genes that are not members of a
common programme.

Open bars use all genes, in each view's colour. Right: the 34 per-class explanations of the cancer-type models. Many AS hits there came
from alternative 5'/3' splice sites in cytosolic ribosomal-protein genes (RPS24, RPL13A, RPL3 and others), which are
deeply covered and so measured precisely in every tissue. Filled bars repeat the tests with RPL*/RPS* genes removed
from both ranking and universe; that sensitivity test was run only on the cancer-type AS views, so GEX and the
left panel have no filled bar. Most SE, MXE and A3SS classes keep enrichments. The remaining terms are lineage
biology (synapse and cell projection in PCPG, smooth-muscle contraction in PRAD, RHO GTPase cycle in DLBC,
extracellular matrix in THCA) and the splicing machinery's own transcripts (RNA splicing terms from A3SS in ACC,
KICH, LUAD and UVM). Top terms per class: `cancer_type_top_terms_by_class.csv`.
""")


# ============================================================ 04 events


def T5() -> None:
    d = BIO / "06_lineage"
    rho = pd.read_csv(d / "event_vs_host_expression.csv")
    ov = pd.read_csv(d / "top_gene_overlap_gex_vs_as.csv")
    fig, axes = plt.subplots(1, 2, figsize=(4.8, 2.0), gridspec_kw={"width_ratios": [1.25, 1.0], "wspace": 0.45})
    ax = axes[0]
    x = np.arange(len(AS_VIEWS))
    rows = []
    for k, (st, col, nm) in enumerate([("background", MUTED, "random measured events"),
                                       ("top", ACCENT, "per-class top events")]):
        for i, v in enumerate(AS_VIEWS):
            s = rho[(rho.view == v) & (rho.set == st)]["abs_rho"].dropna()
            q1, q2, q3 = np.percentile(s, [25, 50, 75])
            xi = i + (k - 0.5) * 0.32
            ax.vlines(xi, q1, q3, color=col, linewidth=1.6)
            ax.scatter(xi, q2, s=7, color=INK, linewidths=0, zorder=3)
            rows.append({"panel": "abs_rho", "view": v, "set": st, "q25": q1, "median": q2, "q75": q3, "n": len(s)})
        ax.plot([], [], color=col, linewidth=1.6, label=nm)
    ax.set_xticks(x)
    ax.set_xticklabels([VIEW_NAMES[v] for v in AS_VIEWS], rotation=90)
    ax.set_ylim(0, 1)
    ax.set_ylabel("|Spearman ρ|, event vs host-gene\nexpression (median, IQR)")
    ax.legend(loc="upper left", handlelength=1.0)
    ax = axes[1]
    agg = ov.groupby("view")[["overlap", "expected"]].sum().reindex(AS_VIEWS)
    ax.bar(x - 0.19, agg["expected"], width=0.36, color=MUTED, label="expected")
    ax.bar(x + 0.19, agg["overlap"], width=0.36, color=ACCENT, label="observed")
    ax.set_xticks(x)
    ax.set_xticklabels([VIEW_NAMES[v] for v in AS_VIEWS], rotation=90)
    ax.set_ylabel("top-100 GEX genes that are top-100\nAS host genes (sum, 34 classes)")
    ax.legend(loc="upper left", handlelength=1.0)
    rows += [{"panel": "overlap", "view": v, "observed": agg.loc[v, "overlap"], "expected": agg.loc[v, "expected"]}
             for v in AS_VIEWS]
    save(fig, OUT, "T5_not_expression", data=pd.DataFrame(rows), caption="""
**Is the splicing signal for cancer type a proxy for expression?**

Left: for each AS view, the per-class top 20 events (union over the 34 classes) and 1,000 random measured events.
Each event's value is correlated (Spearman, 3,000 random samples) with its host gene's log2 CPM.
- Splice-site views (SE, A5SS, A3SS, MXE, RI): the top events correlate with host expression about as much as random
  events (median |ρ| 0.13-0.18 vs 0.10-0.15).
- Terminal-exon views (AFE, ALE, HIT): the top events track expression more (0.17-0.24 vs 0.03-0.05), as expected
  when an alternative promoter or poly(A) site also changes how much of the gene is made.

Right: for each class, the overlap between the top 100 GEX genes and the host genes of each AS view's top 100
events, within genes measured by both, summed over the 34 classes. The overlap is 1.3-3.0 times chance but small:
about one gene per class. The shared genes are lineage markers that are both expressed and spliced
lineage-specifically, such as GFAP and MBP (glioma), DUOXA2 (thyroid), MYH11 (smooth muscle), LGALS4
(gastrointestinal) and FXYD2 (kidney). Per-class genes: `top_gene_overlap_gex_vs_as.csv`.
""")


def T6() -> None:
    d = BIO / "06_lineage"
    med = pd.read_csv(d / "atlas_se_median_psi_by_class.csv", index_col=0)
    picks = pd.read_csv(d / "atlas_se_events.csv").drop_duplicates("feature_id").set_index("feature_id").loc[med.index]
    dev = med.sub(med.median(axis=1), axis=0)                      # class median - median over classes
    Z = dev.fillna(0)
    col_order = list(dev.columns[leaves_list(linkage(Z.T.to_numpy(), "average", metric="euclidean"))])
    pos = {c: i for i, c in enumerate(col_order)}
    own = np.array([dev.loc[f, c] for f, c in zip(dev.index, picks["class"])])
    order = sorted(range(len(dev)), key=lambda i: (pos[picks["class"].iloc[i]], -abs(own[i])))
    M = dev.iloc[order][col_order]
    cls = picks["class"].iloc[order].to_numpy()
    lab = [f"{locus_label(s, f)} ({c})"
           for s, f, c in zip(picks["symbol"].iloc[order], picks.index[order], cls)]
    lim = 60
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.72, 0.068 * len(M) + 0.9))
    im = ax.imshow(M.to_numpy(), cmap=DIVERGE, vmin=-lim, vmax=lim, aspect="auto")
    ax.scatter([pos[c] for c in cls], np.arange(len(M)), s=1.2, color=INK, linewidths=0)
    heat_axes(ax, lab, col_order, xrot=90, row_size=4.3, col_size=5.5)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", labelrotation=90)
    colourbar(fig, im, ax, "PSI in class − median over classes (×100)", shrink=0.3, extend="both")
    data = M.reset_index().melt(id_vars="feature_id", var_name="class", value_name="psi_minus_median_x100")
    data = data.merge(med.reset_index().melt(id_vars="feature_id", var_name="class", value_name="median_psi_x100"),
                      on=["feature_id", "class"])
    data = data.merge(picks.reset_index().rename(columns={"class": "ranked_by_class"}), on="feature_id")
    save(fig, OUT, "T6_lineage_atlas", data=data, caption="""
**Do cassette exons alone mark cell lineage?**

Rows: the two cassette exons with the largest per-class SHAP importance for each of the 34 classes (row label: host
gene, then the class that ranks the exon). Rows are grouped by that class; columns are classes, ordered by
average-linkage clustering.

Colour is the exon's median inclusion in a class minus its median over all classes (PSI ×100). Red means the exon
is included more in that class, blue that it is skipped more. Dots mark the class that ranks the exon, so a coloured
diagonal means each class is recognised by exons switched specifically in it. Off-diagonal blocks show exons shared
by related classes.

Median PSI, the ranking class and the coordinates are in the CSV.
""")


# ============================================================ 07 site probe


def T7() -> None:
    d = BIO / "07_site_probe"
    w = pd.read_csv(d / "site_probe_summary.csv")
    first = ["COAD vs READ", "ESCA vs STAD", "GBM vs LGG", "KIRC vs KIRP", "LUAD vs LUSC", "CESC vs HNSC"]
    pairs = [x for x in first if x in set(w["pair"])] + [x for x in dict.fromkeys(w["pair"]) if x not in first]
    fig, axes = plt.subplots(1, len(pairs), figsize=(DOUBLE, 1.9), sharey=True, gridspec_kw={"wspace": 0.12})
    y = np.arange(len(CFG.views))
    for ax, pr in zip(axes, pairs):
        dd = w[w.pair == pr].set_index("view").reindex(CFG.views)
        cols = [INK2 if v == "gex" else ACCENT for v in CFG.views]
        ax.hlines(y, dd["site-held-out"], dd["case-grouped"], color=GRID, linewidth=1.2)
        ax.scatter(dd["case-grouped"], y, s=9, facecolor="white", edgecolor=cols, linewidth=0.6, zorder=3)
        ax.scatter(dd["site-held-out"], y, s=9, color=cols, linewidths=0, zorder=4)
        ax.set_title(pr, loc="left", fontsize=6.5, color=INK)
        ax.set_xlim(0.45, 1.02)
        ax.set_xticks([0.5, 0.75, 1.0])
        ax.set_xticklabels(["0.5", "", "1"])
        ax.axvline(0.5, color=MUTED, linewidth=0.4)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([VIEW_NAMES[v] for v in CFG.views], color=INK)
    axes[0].invert_yaxis()
    h = [Line2D([], [], marker="o", linestyle="", markersize=3, markerfacecolor="white", markeredgecolor=INK2,
                markeredgewidth=0.6, label="new patients, seen sites"),
         Line2D([], [], marker="o", linestyle="", markersize=3, color=INK2, markeredgewidth=0, label="unseen sites")]
    axes[-1].legend(handles=h, loc="lower right", bbox_to_anchor=(1.0, 1.08), ncol=2)
    fig.supxlabel("balanced accuracy of a simple probe (mean of 3 fold seeds)", fontsize=7, color=INK2, y=-0.08)
    save(fig, OUT, "T7_site_probe", data=w, caption="""
**Is the separation of neighbouring cancer types a tissue-source-site effect?**

In TCGA every tissue-source site contributes to exactly one project (COAD and READ share 0 of 38 sites; GBM and LGG
0 of 40), so a model could separate neighbouring projects by site handling. For six pairs of related classes, the
same fixed, untuned probe is cross-validated two ways:
- by patient (open circles), as in phase 1;
- by site (filled), so that every test site is unseen (D-20).

The probe is top 2,000 variance → median impute → standardise → L2 logistic regression, all fitted inside training
folds. GEX is grey, AS views blue; points are the mean of 3 fold seeds.
- **Colon vs rectum:** no AS view separates the two beyond chance at unseen sites (0.48-0.57). GEX keeps 0.66. The
  AS ensemble's better colon/rectum calls in the main models (T2, T3) should therefore not be read as biology.
- **Oesophagus vs stomach:** AS views lose 0.08-0.14 of balanced accuracy at unseen sites, GEX 0.06. Part of the
  splicing advantage there is site-specific.
- **Glioma grade, kidney subtypes and lung histology:** the separation survives unseen sites in every view (drop
  ≤ 0.015). The near-perfect splicing separation of GBM from LGG is robust.

The probe is diagnostic: it does not re-estimate the phase-1 models. Values per seed: `site_probe_by_seed.csv`.
""")


FIGS = {f.__name__: f for f in [T1, T2, T3, T4, T5, T6, T7]}

if __name__ == "__main__":
    run(FIGS)
