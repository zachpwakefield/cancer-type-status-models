"""Step 15: the figures of the biology layer, in the house style (docs/FIGURE_GUIDELINES.md).

Reads only the tables written by steps 08-14 (and, for the few panels that plot individual
samples, the local caches and matrices). Each figure answers one question, stated at the top
of its caption file. Every figure is written as SVG, PDF and a 400-dpi PNG, with `<stem>.csv`
holding the plotted values and `<stem>.caption.md`. Per-sample plotted values go to
`<data_dir>/figure_data/` (individual-level, D-12).

Statistics on figures are raw p (stars or values; tiers in the caption). Any q is in the CSV
with its FDR family named.

Usage: python scripts/15_bio_figures.py [fig01 fig07 ...]   (default: all)
"""
from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import log
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import leaves_list, linkage

from ctsm.biology import ProteinMap, bio_config, feature_values, gex_by_symbol, task_table
from ctsm.config import PAN, load_config
from ctsm.figstyle import (ACCENT, ACCENT_DARK, DIVERGE, DOUBLE, GRID, INK, INK2, MUTED, NORMAL, SEQ, SINGLE, TUMOUR,
                           VIEW_NAMES, apply_style, bracket, colourbar, fp, gene_label, heat_axes, save, stars)

apply_style()
CFG = load_config()
BIO = CFG.results_dir / "biology"
LOCAL = CFG.data_dir / "figure_data"
AS_VIEWS = [v for v in CFG.views if v != "gex"]
SCOPE = pd.read_csv(CFG.results_dir / "data" / "task_scope.csv")
CANCERS = list(SCOPE.loc[SCOPE.tumor_status_in_scope, "project"])


def short(region: str) -> str:
    ch, rg = region.split(":")
    a, b = rg.split("-")
    return f"{ch}:{int(a):,}-{int(b):,}"


# ============================================================ 01 complementarity
def fig01() -> None:
    d = BIO / "01_complementarity"
    perf = pd.read_csv(d / "cancer_type_fusion_metrics.csv").set_index("model")
    tests = pd.read_csv(d / "cancer_type_mcnemar.csv")
    order = [VIEW_NAMES[v] for v in CFG.views] + ["AS ensemble", "GEX + AS (1:1)"]
    col = {m: MUTED for m in order}
    col.update({"GEX": INK, "AS ensemble": ACCENT, "GEX + AS (1:1)": ACCENT_DARK})
    p = perf.loc[order]
    x = np.arange(len(order)) + np.r_[np.zeros(len(CFG.views)), np.ones(2) * 0.6]
    fig, ax = plt.subplots(figsize=(SINGLE, 2.3))
    for xi, (m, r) in zip(x, p.iterrows()):
        ax.vlines(xi, r.acc_lo, r.acc_hi, color=col[m], linewidth=0.8)
        ax.scatter(xi, r.accuracy, s=12, color=col[m], linewidths=0, zorder=3)
        ax.text(xi, r.acc_hi + 0.0012, f"{r.accuracy:.3f}", fontsize=4.8, color=INK2, ha="center", va="bottom")
    xpos = dict(zip(order, x))
    pairs = [("GEX", "AS ensemble"), ("GEX", "GEX + AS (1:1)")]
    top = p.acc_hi.max()
    rows = []
    for k, (a, b) in enumerate(pairs):
        t = tests[(tests.A == a) & (tests.B == b)].iloc[0]
        bracket(ax, xpos[a], xpos[b], top + 0.008 + k * 0.009, f"p = {fp(t.p_mcnemar)}", 0.002, fontsize=5.5)
        rows.append({"A": a, "B": b, "only_A_correct": t.only_A_correct, "only_B_correct": t.only_B_correct,
                     "p_mcnemar_exact": t.p_mcnemar})
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha="right", color=INK)
    ax.set_ylabel("accuracy, 34 classes (out-of-fold)")
    ax.set_ylim(p.acc_lo.min() - 0.005, top + 0.03)
    data = p.reset_index()[["model", "accuracy", "acc_lo", "acc_hi"]]
    pd.DataFrame(rows).to_csv(d / "fig01_accuracy_by_view_tests.csv", index=False)
    save(fig, d, "fig01_accuracy_by_view", data=data, caption=f"""
**Does splicing add information that expression lacks, for cancer type?**

Out-of-fold accuracy over 9,729 samples and 34 classes (33 TCGA projects + pooled normal) for each single-view model
(grey; GEX in black), the AS ensemble (mean of the 8 AS views' probabilities, blue) and an equal-weight fusion of GEX
and the AS ensemble (dark blue). The fusion weights were fixed before looking at any result (D-14). Whiskers: 95%
case-level bootstrap CIs (1,000 resamples).

Brackets: exact McNemar tests on paired per-sample correctness, raw p (two tests; no FDR correction). The values
are in `fig01_accuracy_by_view_tests.csv`. Single AS views are less accurate than GEX, and the AS ensemble is
level with it. Fused with GEX, the AS views correct {int(rows[1]['only_B_correct'])} GEX errors and introduce
{int(rows[1]['only_A_correct'])}.
""")


def fig02() -> None:
    d = BIO / "01_complementarity"
    r = pd.read_csv(d / "cancer_type_recall_gex_vs_as.csv", index_col=0).sort_values("AS_minus_GEX")
    fig, ax = plt.subplots(figsize=(2.5, 4.0))
    y = np.arange(len(r))
    ax.barh(y, r["AS_minus_GEX"], color=[ACCENT if v > 0 else MUTED for v in r["AS_minus_GEX"]], height=0.72)
    ax.axvline(0, color=INK2, linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{k} ({n})" for k, n in zip(r.index, r["n"])], fontsize=5.5, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.7, len(r) - 0.3)
    ax.set_xlabel("recall, AS ensemble − GEX")
    data = r.reset_index().rename(columns={"index": "class"})
    save(fig, d, "fig02_recall_difference_by_class", data=data, caption="""
**Which cancer types does splicing resolve better than expression, and which worse?**

The difference in out-of-fold recall between the AS ensemble and the GEX model for each class (n = samples).
Positive (blue): splicing recovers more of the class.

Splicing gains the most for classes that GEX confuses with an anatomical neighbour (COAD/READ, ESCA/STAD, GBM/LGG).
Of these, only the glioma separation survives holding out whole tissue-source sites (fig15). The colon/rectum gain
does not, and should not be read as biology.

GEX is better for pooled normal tissue, uterine carcinosarcoma, and squamous and urothelial tumours (BLCA, HNSC,
CESC). Values: `cancer_type_recall_gex_vs_as.csv` (with the fused model's recall).
""")


# ============================================================ 02 global shifts
SHIFT_ORDER = ["ri", "a5ss", "a3ss", "se", "mxe", "ale", "afe"]
SHIFT_LAB = ["RI retained", "A5SS long form", "A3SS long form", "SE included", "MXE upstream exon",
             "ALE downstream", "AFE downstream"]


def fig03() -> None:
    d = BIO / "02_global_shifts"
    s = pd.read_csv(d / "shift_summary.csv")
    piv = s.pivot(index="cancer", columns="view", values="median_paired_shift")[SHIFT_ORDER] * 100
    pv = s.pivot(index="cancer", columns="view", values="p_wilcoxon")[SHIFT_ORDER]
    order = piv["ri"].sort_values(ascending=False).index
    piv, pv = piv.loc[order], pv.loc[order]
    lim = float(np.nanpercentile(np.abs(piv.to_numpy()), 98))
    fig, ax = plt.subplots(figsize=(3.2, 3.4))
    im = ax.imshow(piv.to_numpy(), cmap=DIVERGE, vmin=-lim, vmax=lim, aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.iat[i, j]
            ax.text(j, i, f"{v:+.1f}{stars(pv.iat[i, j])}", ha="center", va="center", fontsize=4.8,
                    color="white" if abs(v) > 0.65 * lim else INK)
    npairs = s.drop_duplicates("cancer").set_index("cancer").loc[order, "n_pairs"]
    heat_axes(ax, [f"{c} ({n})" for c, n in npairs.items()], SHIFT_LAB, col_size=5.5, row_size=6)
    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=45, ha="left", rotation_mode="anchor")
    colourbar(fig, im, ax, "median paired tumour − normal (×100)")
    data = s[["cancer", "view", "n_pairs", "n_events", "median_paired_shift", "frac_patients_up", "frac_up",
              "frac_down", "p_wilcoxon", "q_wilcoxon"]].assign(
        q_family="BH across 15 cancers x 7 views (105 Wilcoxon tests)")
    save(fig, d, "fig03_paired_shift_heatmap", data=data, caption="""
**Do tumours shift RNA processing in a common direction relative to their own normal tissue?**

Each cell is the median over matched patients (n in brackets) of the patient's mean paired tumour − normal change,
averaged over every event measured in that cancer, so there is no SHAP selection and no circularity. The columns
say what an increase means:
- RI: more intron retention.
- A5SS/A3SS: more use of the long exon form, which includes extra intronic sequence.
- SE: more cassette inclusion.
- MXE: more of the upstream exon.
- ALE/AFE: the used last or first exons move downstream (a 0-1 position index along transcription).

Stars: two-sided Wilcoxon signed-rank test across patients, raw p (* < 0.05, ** < 0.01, *** < 0.001). BH q
(family: 15 cancers × 7 views) is in the CSV.

Retention and long splice-site forms rise in most cancers. Last exons move upstream and never downstream. BRCA is
the exception for RI. Caveat: RNA integrity and nuclear-RNA content also raise intronic reads, and RIN is not
available (D-15).
""")


def fig04() -> None:
    d = BIO / "02_global_shifts"
    pats = pd.read_csv(CFG.data_dir / "bio_cache" / "paired_shift_per_patient.csv")
    s = pd.read_csv(d / "shift_summary.csv")
    order = s[s.view == "ri"].set_index("cancer")["median_paired_shift"].sort_values(ascending=False).index
    fig, axes = plt.subplots(1, 2, figsize=(4.4, 2.8), sharey=True)
    rng = np.random.default_rng(0)
    plotted, summ = [], []
    for ax, view, xl in [(axes[0], "ri", "RI: paired Δ PSI (×100)"),
                         (axes[1], "ale", "ALE: paired Δ position (×100)")]:
        sub = pats[pats.view == view]
        for i, c in enumerate(order):
            v = sub.loc[sub.cancer == c, "paired_shift"].to_numpy() * 100
            yj = i + rng.uniform(-0.25, 0.25, len(v))
            ax.scatter(v, yj, s=1.8, color=MUTED, alpha=0.7, linewidths=0, rasterized=True)
            ax.plot([np.median(v)] * 2, [i - 0.36, i + 0.36], color=INK, linewidth=1.0)
            plotted.append(pd.DataFrame({"view": view, "cancer": c, "case_id": sub.loc[sub.cancer == c, "case_id"],
                                         "paired_shift_x100": v}))
            summ.append({"view": view, "cancer": c, "n_patients": len(v), "median_x100": np.median(v),
                         "frac_patients_up": float(np.mean(v > 0))})
        ax.axvline(0, color=INK2, linewidth=0.5)
        ax.set_xlabel(xl)
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels(order, color=INK)
    axes[0].invert_yaxis()
    for ax in axes:
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    save(fig, d, "fig04_paired_shift_patients", data=pd.DataFrame(summ),
         local_data=(LOCAL, pd.concat(plotted, ignore_index=True)), caption="""
**Is the shift carried by most patients, or by a few?**

Each dot is one patient's tumour compared with their own normal tissue: the mean paired change over all measured RI
events (left) and in the ALE last-exon position index (right). Bars: medians. Patients are ordered as in fig03.
In KICH, KIRC, PRAD and BLCA, about 70-90% of patients retain more intronic sequence in the tumour. The ALE shift
upstream is smaller and concentrated in fewer cancers.

The CSV gives the per-cancer medians and shares of patients. The per-patient values are individual-level and stay
in the local `figure_data/` (D-12).
""")


# ============================================================ 03 pathways
def fig05() -> None:
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
    heat_axes(ax, [t.replace("HALLMARK_", "").replace("_", " ").lower() for t in terms], CANCERS, xrot=90,
              row_size=5.5, col_size=5.5)
    ax.grid(which="minor", visible=False)
    ax.set_xlim(-0.6, len(CANCERS) - 0.4)
    ax.set_ylim(len(terms) - 0.5, -0.5)
    h = [Line2D([], [], marker="o", linestyle="", markersize=np.sqrt(6 + 9 * n), color=ACCENT, markeredgewidth=0,
                label=f"{n} AS view{'s' if n > 1 else ''}") for n in (1, 3)]
    h.append(Line2D([], [], marker="o", linestyle="", markersize=np.sqrt(40), markerfacecolor="none",
                    markeredgecolor=INK, markeredgewidth=0.6, label="GEX"))
    ax.legend(handles=h, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, columnspacing=0.8, handletextpad=0.1)
    save(fig, d, "fig05_hallmarks_tumour_status", data=pd.DataFrame(rows).assign(
        q_family="BH within Hallmark per model (one model = one cancer x view)"), caption="""
**Which programmes do the tumour-vs-normal models draw on?**

For every cancer (columns) and view, the 200 genes with the largest SHAP importance were tested for enrichment of
the MSigDB Hallmark sets against all genes measurable in that view and cohort. The test is GOseq/Wallenius,
correcting for the number of events per gene (D-16).
- A dot marks a significant set in at least one AS view; its size is the number of AS views (of 8).
- A ring marks significance in the GEX model.
- Significance is q < 0.05, BH within the Hallmark collection of each model.

GEX models consistently use stromal and microenvironment programmes (EMT, coagulation, KRAS, myogenesis). The AS
models reach the same programmes far less often, mostly through EMT and myogenesis genes. So their tumour signal is
largely not organised by pathway membership (fig06).
""")


def fig06() -> None:
    d = BIO / "03_pathways"
    coh = pd.read_csv(d / "pathway_coherence_by_view.csv", index_col=0).reindex(CFG.views)
    c_ts, c_ct, c_rp = coh.columns
    fig, axes = plt.subplots(1, 2, figsize=(3.8, 2.0), sharey=True)
    y = np.arange(len(CFG.views))
    cols = [INK2 if v == "gex" else ACCENT for v in CFG.views]
    axes[0].barh(y, coh[c_ts], color=cols, height=0.65)
    axes[0].set_xlabel("tumour vs normal:\nshare of 15 cancers")
    axes[1].barh(y, coh[c_ct], color="none", edgecolor=cols, linewidth=0.6, height=0.65, label="all genes")
    axes[1].barh(y, coh[c_rp].fillna(0), color=ACCENT, height=0.65, label="ribosomal-protein\ngenes removed")
    axes[1].set_xlabel("cancer type:\nshare of 34 classes")
    axes[1].legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, handlelength=1.0, columnspacing=0.8)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([VIEW_NAMES[v] for v in CFG.views], color=INK)
    axes[0].invert_yaxis()
    save(fig, d, "fig06_pathway_coherence", data=coh.reset_index().rename(columns={"index": "view"}).assign(
        q_family="BH within collection per model (Hallmark, Reactome, GO-BP)"), caption="""
**Do the genes a splicing model relies on form pathways, as an expression model's do?**

The share of models whose top 200 SHAP genes are enriched (q < 0.05, BH within collection per model; GOseq/Wallenius,
D-16) for at least one Hallmark, Reactome or GO-BP term, by view.

Left: the 15 tumour-vs-normal models. GEX reaches a pathway in 14 of 15 cancers; SE in 9; MXE in 6; the other AS
views in at most 4; RI in none. The tumour signal of splicing is mostly carried by genes that are not members of a
common programme.

Right: the 34 per-class explanations of the cancer-type models. Open bars use all genes. Many AS hits there came
from alternative 5'/3' splice sites in cytosolic ribosomal-protein genes (RPS24, RPL13A, RPL3 and others), which are
deeply covered and so measured precisely in every tissue. Filled bars repeat the tests with RPL*/RPS* genes removed
from both ranking and universe. Most SE, MXE and A3SS classes keep enrichments. The remaining terms are lineage
biology (synapse and cell projection in PCPG, smooth-muscle contraction in PRAD, RHO GTPase cycle in DLBC,
extracellular matrix in THCA) and the splicing machinery's own transcripts (RNA splicing terms from A3SS in ACC,
KICH, LUAD and UVM). Top terms per class: `cancer_type_top_terms_by_class.csv`.
""")


# ============================================================ 04 events
def fig07() -> None:
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
    im = ax.imshow(M.to_numpy(), cmap=DIVERGE, vmin=-lim, vmax=lim, aspect="auto")
    yy, xx = np.nonzero(top.to_numpy())
    ax.scatter(xx, yy, s=1.5, color=INK, linewidths=0)
    labels = [f"{gene_label(s)} {VIEW_NAMES[v]} {short(r)}" for s, v, r in zip(show.symbol, show.view, show.region)]
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
    save(fig, d, "fig07_recurrent_events", data=long, caption=f"""
**Do different cancers rely on the same events to tell tumour from normal?**

Across the 8 AS views, {n_obs} events rank in the top 100 (mean |SHAP|) of at least 4 of the 15 tumour-vs-normal
models. {n_exp:.1f} would be expected if the rankings were independent (Poisson-binomial over each event's measured
cohorts; `recurrence_vs_chance.csv`).

The figure shows the 40 most recurrent events. Rows give the host gene (HGNC), view and coordinates (hg38, 1-based).
Colour is the median paired tumour − normal change in every cancer, including cancers where the event did not rank
(so the pattern across cancers is not selected on); the scale is clipped at ±{lim}. Dots mark the cancers whose top
100 contain the event.

Recurrent events tend to move the same way across cancers. For the median event, {cons.median():.0%} of the cancers
where it is measured share its majority direction; {(cons >= 0.8).mean():.0%} of events reach 80%.

Several genes appear as switch pairs, with two exons of the same gene moving in opposite directions (EXOSC7 last
exons, NDE1 first exons, TPM1 mutually exclusive exons).

Recurrent RI and A5SS events mostly rise in tumours ({up['ri']:.0%} and {up['a5ss']:.0%}), matching the global shift in fig03. The full
list with direction consistency is in `recurrent_events.csv`.
""")


def fig08() -> None:
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
    save(fig, d, "fig08_exon_protein_features", data=data, caption="""
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


def fig09() -> None:
    d = BIO / "04_events"
    rec = pd.read_csv(d / "recurrent_events.csv")
    P = pd.read_csv(d / "paired_change_top_events_all_cancers.csv.gz")
    picks = rec[rec.n_cancers_same_sign_as_majority >= 0.8 * rec.n_cancers_measured].drop_duplicates("symbol").head(6)
    pm = ProteinMap(CFG)
    fig, axes = plt.subplots(len(picks), 2, figsize=(DOUBLE * 0.82, 1.15 * len(picks)),
                             gridspec_kw={"width_ratios": [2.4, 1.0], "hspace": 0.9, "wspace": 0.25})
    meds, local, info_rows = [], [], []
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
            info_rows.append({"symbol": r.symbol, "region": r.region, **c})
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
    save(fig, d, "fig09_event_portraits", data=pd.concat(meds, ignore_index=True),
         local_data=(LOCAL, pd.concat(local, ignore_index=True)), caption="""
**What do the most shared tumour events look like, and how does the model read them?**

The six most recurrent events (fig07) whose paired change has the same sign in at least 80% of the cancers where
measured, one per gene.

Left: the median value in matched normal (orange) and tumour (blue) samples of each cancer. Cancers whose top 100
contains the event are in black, the others in grey. The heading gives the host gene, view, coordinates and, for
cassette exons, the protein consequence (GENCODE v45; D-17).

Right: SHAP dependence in the cancer where the event ranks highest. Each held-out sample's event value is plotted
against its SHAP value on the tumour log-odds. The value is the PSI; missing values were median-imputed by the
model.

The CSV gives the per-cancer medians and paired Wilcoxon p (raw). Per-sample values stay local (D-12).
""")
    pd.DataFrame(info_rows).to_csv(d / "fig09_event_portraits_protein.csv", index=False)


# ============================================================ 05 regulators
def fig10() -> None:
    d = BIO / "05_regulators"
    rc = bio_config()["regulators"]
    de = pd.read_csv(d / "regulator_paired_expression.csv")
    fc = de.pivot(index="gene", columns="cancer", values="paired_log2fc").reindex(columns=CANCERS)
    pv = de.pivot(index="gene", columns="cancer", values="p_paired").reindex(columns=CANCERS)
    hit = (pv < 0.05) & (fc.abs() >= 0.5)
    groups = list(rc["splicing_factors"])
    group_of = {g: grp for grp, gs in rc["splicing_factors"].items() for g in gs}
    keep = hit.sum(1)[lambda x: x >= 4].index
    order = sorted(keep, key=lambda g: (groups.index(group_of[g]), -fc.loc[g].median()))
    fig, ax = plt.subplots(figsize=(3.4, 0.105 * len(order) + 0.9))
    lim = 1.5
    im = ax.imshow(fc.loc[order].to_numpy(), cmap=DIVERGE, vmin=-lim, vmax=lim, aspect="auto")
    yy, xx = np.nonzero(hit.loc[order].to_numpy())
    white = np.abs(fc.loc[order].to_numpy()[yy, xx]) > 0.9
    ax.scatter(xx, yy, s=1.5, c=np.where(white, "white", INK), linewidths=0)
    heat_axes(ax, [f"{g}  {group_of[g]}" for g in order], CANCERS, xrot=90, row_size=5, col_size=5.5)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", labelrotation=90)
    colourbar(fig, im, ax, "paired log2 FC, tumour − normal", shrink=0.4, extend="both")
    data = fc.loc[order].reset_index().melt(id_vars="gene", var_name="cancer", value_name="paired_log2fc")
    data = data.merge(de[["gene", "cancer", "p_paired", "q_paired"]], on=["gene", "cancer"], how="left")
    data["family"] = data["gene"].map(group_of)
    data["q_family"] = "BH across all regulators x 15 cancers (paired Wilcoxon on log2 CPM)"
    save(fig, d, "fig10_regulator_expression", data=data, caption="""
**Which splicing regulators change between tumour and matched normal tissue?**

The median paired tumour − normal change in expression (log2 CPM) of every pre-specified RNA-binding protein
(config/biology.yaml) that changes in at least 4 cancers. "Changes" means a paired Wilcoxon raw p < 0.05 and
|log2 FC| ≥ 0.5, marked by a dot. Rows are grouped by family.

Two programmes stand out:
- Rising in tumours: SR proteins, hnRNPs, PTBP1 and the cleavage and polyadenylation factors (CSTF2, CPSF1/4/6),
  together with the epithelial regulator ESRP1.
- Falling: the muscle, neural and mesenchymal regulators (CELF2, RBM20, RBM24, NOVA1/2, QKI, MBNL), partly because
  normal tissue contains more stroma and muscle.

BH q (family named) is in the CSV.
""")


def fig11() -> None:
    d = BIO / "05_regulators"
    best = pd.read_csv(d / "top_event_best_regulator.csv")
    best = best.reindex(best.partial_rho.abs().sort_values(ascending=False).index)
    ex = best.drop_duplicates("best_regulator").head(3)
    fig, axes = plt.subplots(1, len(ex), figsize=(1.75 * len(ex) + 0.3, 1.75), gridspec_kw={"wspace": 0.55})
    summ, local = [], []
    for ax, (_, r) in zip(np.atleast_1d(axes), ex.iterrows()):
        tab = task_table(CFG, "tumor_status", r.cancer)
        g = gex_by_symbol(CFG, tab["file_id"], [r.best_regulator])[r.best_regulator].to_numpy()
        v = feature_values(CFG, r["view"], tab["file_id"], [r.feature_id])[:, 0] * 100
        yv = tab["y"].to_numpy()
        for lab, colr in [(0, NORMAL), (1, TUMOUR)]:
            m = yv == lab
            ax.scatter(g[m], v[m], s=2.5, color=colr, alpha=0.7, linewidths=0, rasterized=True)
        ax.set_xlabel(f"{r.best_regulator} log2 CPM")
        ax.set_ylabel(f"{gene_label(r.symbol)} {VIEW_NAMES[r['view']]} PSI (×100)")
        ax.set_title(f"{r.cancer}   partial ρ = {r.partial_rho:+.2f}", loc="left", fontsize=6.5, color=INK)
        summ.append({"cancer": r.cancer, "regulator": r.best_regulator, "symbol": r.symbol, "view": r["view"],
                     "feature_id": r.feature_id, "partial_rho_tumours": r.partial_rho, "n_tumours": int((yv == 1).sum()),
                     "n_normals": int((yv == 0).sum())})
        local.append(pd.DataFrame({"cancer": r.cancer, "regulator": r.best_regulator, "feature_id": r.feature_id,
                                   "log2cpm": g, "psi_x100": v, "is_tumour": yv}))
    h = [Line2D([], [], marker="o", linestyle="", markersize=2.5, color=c, label=l) for c, l in
         [(NORMAL, "normal"), (TUMOUR, "tumour")]]
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, handletextpad=0.1)
    save(fig, d, "fig11_regulator_examples", data=pd.DataFrame(summ),
         local_data=(LOCAL, pd.concat(local, ignore_index=True)), caption="""
**What do the strongest regulator–event links look like?**

Among the top events of every tumour-vs-normal model (top 50 per view; SE, A3SS, A5SS, MXE, RI), the event-regulator
pairs with the largest partial Spearman ρ within tumours (adjusted for stromal and immune marker scores, D-18). One
pair is shown per regulator, chosen by |ρ| without curation. Each dot is a sample (tumour blue, normal orange); ρ
uses tumours only.

The ESRP1-FGFR2 exon switch is a classic epithelial splicing programme, here recovered unsupervised.

PABPN1, a nuclear poly(A)-binding protein involved in nuclear RNA surveillance, tracks retained-intron PSI across
tumours in several cancers (`top_event_best_regulator.csv`). This is consistent with a sample-level retention axis.
The correlations do not show which way any effect runs.

Per-sample values stay local (D-12).
""")


# ============================================================ 06 lineage
def fig12() -> None:
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
    save(fig, d, "fig12_class_confusions", data=conf, caption="""
**Do expression and splicing confuse the same cancer types?**

The 16 class pairs with the most out-of-fold confusions (either direction) for the GEX model or the AS ensemble.

Expression mixes up anatomical neighbours: colon and rectum, oesophagus and stomach, lung adenocarcinoma and squamous
carcinoma, the kidney subtypes, and glioma grades. Splicing makes fewer of these errors; GBM/LGG errors fall from 20
to 0. fig15 shows which of these gains survive unseen tissue-source sites: glioma, kidney and lung yes; colon/rectum
no.

Splicing instead confuses squamous and urothelial carcinomas across organs (CESC/HNSC, HNSC/LUSC, BLCA/CESC,
BLCA/HNSC). It also confuses several tumours with normal tissue (PRAD, THCA, LIHC, STAD), as if the splicing state
follows the tissue and the cell type more than the tumour.

""")


def fig13() -> None:
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
    save(fig, d, "fig13_as_vs_expression", data=pd.DataFrame(rows), caption="""
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


def fig14() -> None:
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
    lab = [f"{gene_label(s)} ({c})" for s, c in zip(picks["symbol"].iloc[order], cls)]
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
    save(fig, d, "fig14_se_lineage_atlas", data=data, caption="""
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
def fig15() -> None:
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
    save(fig, d, "fig15_site_probe", data=w, caption="""
**Is the separation of neighbouring cancer types a tissue-source-site effect?**

In TCGA every tissue-source site contributes to exactly one project (COAD and READ share 0 of 38 sites; GBM and LGG
0 of 40), so a model could separate neighbouring projects by site handling. For six pairs of related classes, the
same fixed, untuned probe is cross-validated two ways:
- by patient (open circles), as in phase 1;
- by site (filled), so that every test site is unseen (D-20).

The probe is top 2,000 variance → median impute → standardise → L2 logistic regression, all fitted inside training
folds. GEX is grey, AS views blue; points are the mean of 3 fold seeds.
- **Colon vs rectum:** no AS view separates the two beyond chance at unseen sites (0.48-0.57). GEX keeps 0.66. The
  AS ensemble's better colon/rectum calls in the main models (fig02, fig12) should therefore not be read as biology.
- **Oesophagus vs stomach:** AS views lose 0.08-0.14 of balanced accuracy at unseen sites, GEX 0.06. Part of the
  splicing advantage there is site-specific.
- **Glioma grade, kidney subtypes and lung histology:** the separation survives unseen sites in every view (drop
  ≤ 0.015). The near-perfect splicing separation of GBM from LGG is robust.

The probe is diagnostic: it does not re-estimate the phase-1 models. Values per seed: `site_probe_by_seed.csv`.
""")


FIGS = {f.__name__: f for f in [fig01, fig02, fig03, fig04, fig05, fig06, fig07, fig08, fig09, fig10, fig11, fig12, fig13,
                                 fig14, fig15]}


def main() -> None:
    names = sys.argv[1:] or list(FIGS)
    for n in names:
        try:
            FIGS[n]()
            log(f"{n} written")
        except FileNotFoundError as e:
            log(f"{n} skipped (inputs missing: {e.filename})")


if __name__ == "__main__":
    main()
