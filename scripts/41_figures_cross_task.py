"""Step 41: the cross-task layer, in figures (house style, docs/FIGURE_GUIDELINES.md).

    X1  do the three tasks' models read the same genes?                          (D-48)
    X2  which other cancers does a subtype resemble, in the cancer-type model?   (D-49)
    X3  does the tumour-status model see some subtypes as less tumour-like?      (D-50)
    X4  are the exons each task uses alike in what they do to the protein?       (D-51)

Reads results/across_tasks/ (step 40). Every figure is SVG + PDF + 400-dpi PNG with `<stem>.csv` of
the plotted values and `<stem>.caption.md`. Task colours: tumour status blue, cancer type dark blue,
subtype amber; model colours as in U1 (GEX black, AS ensemble blue, GEX + AS dark blue).
Usage: python scripts/41_figures_cross_task.py [X1 X2 ...]        (default: all)
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from _common import REPO
from _figbase import run, tidy
from matplotlib.lines import Line2D

from ctsm.biology import bh
from ctsm.config import load_config
from ctsm.figstyle import ACCENT, ACCENT_DARK, DIVERGE, DOUBLE, GRID, INK, INK2, MUTED, fp, heat_axes, save

CFG = load_config()
IN = CFG.results_dir / "across_tasks"
OUT = CFG.results_dir / "figures" / "across_tasks"
XC = yaml.safe_load(open(REPO / "config" / "followup.yaml"))["cross_task"]
STATUS, TYPE, SUBTYPE = ACCENT, ACCENT_DARK, "#c98a2b"
MODEL_COL = {"GEX": INK, "AS ensemble": ACCENT, "GEX + AS": ACCENT_DARK}
MINUS = "−"
SHORT = {"Adenocarcinoma": "adeno", "Squamous": "squamous", "IDHmut-codel": "IDH-mut", "IDHmut-non-codel": "IDH-mut",
         "IDHwt": "IDH-wt", "HPV-positive": "HPV+", "HPV-negative": "HPV−", "Cortical admixture": "cortical admixture",
         "Endometrioid": "endometrioid", "Serous": "serous", "Biphasic": "biphasic", "Epithelioid": "epithelioid",
         "Intestinal": "intestinal", "Diffuse": "diffuse", "Follicular": "follicular", "Classical": "classical",
         "Tall cell": "tall cell", "Lobular": "lobular", "Ductal": "ductal"}


def signed(v: float, nd: int = 2) -> str:
    if round(v, nd) == 0:
        return f"{0:.{nd}f}"
    return f"{v:+.{nd}f}".replace("-", MINUS)


def group_name(classes: str, cohort_classes: int) -> str:
    cl = list(dict.fromkeys(SHORT.get(c, c) for c in classes.split("+")))
    return "rest" if len(classes.split("+")) >= 2 and len(classes.split("+")) == cohort_classes - 1 and len(cl) > 1 \
        else "+".join(cl)


def pred_label(r, n_classes: dict, lower_first: bool = False) -> str:
    hi = group_name(r.higher_in, n_classes[r.cohort])
    lo = group_name(r.lower_in, n_classes[r.cohort])
    if lower_first:
        return f"{r.id} {r.cohort}: {lo} vs {hi}"
    tgt = f" → {r.targets}" if "targets" in r._fields else ""
    return f"{r.id} {r.cohort}: {hi} vs {lo}{tgt}"


def n_classes_of() -> dict:
    t = pd.read_csv(IN / "type_calls_by_subtype_class.csv")
    return t[t.model == "GEX"].groupby("cohort")["class"].nunique().to_dict()


def panel_letter(ax, s, x=-0.1, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")


# ============================================================ X1 gene overlap
def X1() -> None:
    d = pd.read_csv(IN / "gene_overlap.csv")
    tests = pd.read_csv(IN / "gene_overlap_tests.csv")
    pairs = [("status-type", "status vs type"), ("subtype-status", "subtype vs status"),
             ("subtype-type", "subtype vs type")]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE * 0.82, 2.5), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.55})
    rng = np.random.default_rng(0)
    rows = []
    for i, (key, name) in enumerate(pairs):
        g = d[d.pair == key]
        none = g.n_both == 0
        y = np.where(none, -2.0, np.log2(g.fold_enrichment.where(~none, 1.0)))     # no shared gene: own row
        x = i + rng.uniform(-0.22, 0.22, len(g))
        ax.scatter(x, y, s=5, color=MUTED, linewidths=0, alpha=0.8, zorder=2)
        ax.plot([i - 0.3, i + 0.3], [np.median(np.log2(g.fold_enrichment.clip(lower=1e-9)))] * 2, color=INK,
                linewidth=1.2, zorder=3)
        rows.append(g.assign(panel="a", log2_fold_enrichment_plotted=y, plotted_as_none_shared=none, x_jitter=x))
    ax.axhline(0, color=INK2, linewidth=0.6, zorder=0)
    ax.set_xticks(range(3))
    ax.set_xticklabels([n.replace(" vs ", " vs\n") + f"\n({int((d.pair == k).sum())})" for k, n in pairs], fontsize=5.6)
    ax.tick_params(axis="x", length=0)
    ax.set_yticks([-2, 0, 2, 4, 6])
    ax.set_yticklabels(["none\nshared", "0", "2", "4", "6"])
    ax.spines["left"].set_bounds(-1, 6)
    ax.set_ylabel("top-100 gene overlap, log2(observed / expected)")
    ax.set_xlim(-0.6, 2.6)
    panel_letter(ax, "a", x=-0.28)
    h = d[(d.label_type == "histology") & d.pair.isin(["subtype-status", "subtype-type"])]
    h = h[h.cohort.isin(set(h[h.pair == "subtype-status"].cohort))]
    cohorts = sorted(h.cohort.unique())
    for i, c in enumerate(cohorts):
        for key, col, off in [("subtype-status", STATUS, -0.14), ("subtype-type", TYPE, 0.14)]:
            v = h[(h.cohort == c) & (h.pair == key)]
            bx.scatter(v.jaccard, [i + off] * len(v), s=5, color=col, alpha=0.45, linewidths=0, zorder=2)
            bx.scatter(v.jaccard.median(), i + off, s=16, color=col, linewidths=0.4, edgecolors="white", zorder=3)
        rows.append(h[h.cohort == c].assign(panel="b"))
    bx.set_yticks(range(len(cohorts)))
    bx.set_yticklabels([f"{c} histology" for c in cohorts], fontsize=5.8, color=INK)
    bx.set_ylim(len(cohorts) - 0.5, -0.6)
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.set_xlabel("Jaccard of the two top-100 gene lists")
    bx.legend(handles=[Line2D([], [], marker="o", ls="", color=STATUS, ms=3.4, label="subtype – tumour status"),
                       Line2D([], [], marker="o", ls="", color=TYPE, ms=3.4, label="subtype – cancer type")],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.6, handletextpad=0.2)
    panel_letter(bx, "b", x=-0.36)
    data = pd.concat(rows, ignore_index=True)
    data["q_family"] = "BH across every (pair, cohort, view) of gene_overlap.csv"
    save(fig, OUT, "X1_gene_overlap", data=data, caption=CAPTIONS["X1"](d, tests))


# ============================================================ X2 type probability by subtype
def _auroc_panel(ax, d, order, models, extra=None, xlabel="AUROC, predicted higher vs lower"):
    for i, rid in enumerate(order):
        for k, m in enumerate(models):
            r = d[(d.id == rid) & (d.model == m)].iloc[0]
            y = i + (k - (len(models) - 1) / 2) * 0.2
            sig = r.q < 0.05
            ax.scatter(r.auroc, y, s=13, color=MODEL_COL[m] if sig else "white", edgecolors=MODEL_COL[m],
                       linewidths=0.7, zorder=3)
        if extra is not None:
            r = extra[extra.id == rid].iloc[0]
            ax.scatter(r.auroc, i + 0.36, s=12, marker="x", color=MUTED, linewidths=0.8, zorder=3)
    ax.axvline(0.5, color=INK2, linewidth=0.6, zorder=0)
    ax.set_xlim(0.2, 1.0)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_ylim(len(order) - 0.45, -0.6)
    ax.set_xlabel(xlabel)


def X2() -> None:
    d = pd.read_csv(IN / "type_predictions.csv")
    sens = pd.read_csv(IN / "type_predictions_offtarget_share.csv")
    nc = n_classes_of()
    order = [p["id"] for p in XC["type_probability_by_subtype"]["predictions"]]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE * 0.74, 0.2 * len(order) + 0.9),
                                 gridspec_kw={"width_ratios": [1.0, 0.9], "wspace": 0.08})
    rows = []
    f = d[d.model == "GEX + AS"].set_index("id")
    for i, rid in enumerate(order):
        r = f.loc[rid]
        ax.plot([r.mean_lower, r.mean_higher], [i, i], color=GRID, linewidth=1.6, solid_capstyle="round", zorder=1)
        ax.scatter(r.mean_lower, i, s=13, color="white", edgecolors=SUBTYPE, linewidths=0.8, zorder=3)
        ax.scatter(r.mean_higher, i, s=13, color=SUBTYPE, linewidths=0.4, edgecolors="white", zorder=3)
    ax.set_xscale("log")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([pred_label(r, nc) for r in d[d.model == "GEX + AS"].set_index("id").loc[order]
                        .reset_index().itertuples()], fontsize=5.6, color=INK)
    ax.set_ylim(len(order) - 0.45, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("mean probability of the targets (GEX + AS)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=SUBTYPE, ms=3.4, label="predicted higher"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=SUBTYPE, ms=3.4, label="predicted lower")],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.6, handletextpad=0.2)
    _auroc_panel(bx, d, order, ["GEX", "AS ensemble", "GEX + AS"], extra=sens[sens.model == "GEX + AS"])
    bx.legend(handles=[Line2D([], [], marker="o", ls="", color=MODEL_COL[m], ms=3.4, label=m) for m in MODEL_COL]
              + [Line2D([], [], marker="x", ls="", color=MUTED, ms=3.4, label="GEX + AS, share of off-target")],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.6, handletextpad=0.2, columnspacing=0.8)
    panel_letter(ax, "a", x=-1.25)
    panel_letter(bx, "b", x=-0.05)
    data = pd.concat([d.assign(measure="summed target probability (pre-specified)"),
                      sens.assign(measure="targets / (1 - own project), GEX + AS shown (not pre-specified)")],
                     ignore_index=True)
    data["q_family"] = "BH across T1-T13 within each model and measure"
    save(fig, OUT, "X2_type_probability_by_subtype", data=data, caption=CAPTIONS["X2"](d, sens))


# ============================================================ X3 status probability by subtype
def X3() -> None:
    d = pd.read_csv(IN / "status_predictions.csv")
    comp = pd.read_csv(IN / "status_probability_vs_composition.csv")
    nc = n_classes_of()
    order = [p["id"] for p in XC["status_probability_by_subtype"]["predictions"]]
    fig = plt.figure(figsize=(DOUBLE * 0.8, 3.9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 0.85], height_ratios=[2.2, 1.0], wspace=0.08, hspace=0.75)
    ax, bx, cx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    rows = []
    for i, rid in enumerate(order):
        for k, m in enumerate(["GEX", "AS ensemble"]):
            r = d[(d.id == rid) & (d.model == m)].iloc[0]
            y = i + (k - 0.5) * 0.3
            lo, hi = 1 - r.mean_lower, 1 - r.mean_higher          # P(not tumour); "lower_in" = less tumour-like
            ax.plot([hi, lo], [y, y], color=GRID, linewidth=1.4, solid_capstyle="round", zorder=1)
            ax.scatter(lo, y, s=12, color=MODEL_COL[m], linewidths=0.4, edgecolors="white", zorder=3)
            ax.scatter(hi, y, s=12, color="white", edgecolors=MODEL_COL[m], linewidths=0.7, zorder=3)
    ax.set_xscale("log")
    ax.set_yticks(range(len(order)))
    lab = d[d.model == "GEX"].set_index("id").loc[order].reset_index()
    ax.set_yticklabels([pred_label(r, nc, lower_first=True) for r in lab.itertuples()], fontsize=5.6, color=INK)
    ax.set_ylim(len(order) - 0.45, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("mean 1 − P(tumour), out-of-fold")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=INK2, ms=3.4, label="predicted less tumour-like"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=INK2, ms=3.4, label="comparison")],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.6, handletextpad=0.2)
    _auroc_panel(bx, d, order, ["GEX", "AS ensemble", "GEX + AS"], xlabel="AUROC, comparison above predicted")
    bx.legend(handles=[Line2D([], [], marker="o", ls="", color=MODEL_COL[m], ms=3.4, label=m) for m in MODEL_COL],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5.6, handletextpad=0.2, columnspacing=0.8)
    rows.append(d.assign(panel="a, b"))
    projects = sorted(comp.project.unique())
    rr = [(s_, m) for s_ in ["stromal", "immune"] for m in ["GEX", "AS ensemble", "GEX + AS"]]
    M = np.array([[comp[(comp.project == p) & (comp.model == m) & (comp.score == s_)].rho.iloc[0] for p in projects]
                  for s_, m in rr])
    cx.imshow(M, cmap=DIVERGE, vmin=-0.35, vmax=0.35, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            cx.text(j, i, signed(M[i, j]), ha="center", va="center", fontsize=5.0, color=INK)
    heat_axes(cx, [f"{s_} score · {m}" for s_, m in rr], projects, row_size=5.6, col_size=5.8)
    cx.xaxis.tick_top()
    rows.append(comp.assign(panel="c"))
    panel_letter(ax, "a", x=-1.0)
    panel_letter(bx, "b", x=-0.05)
    panel_letter(cx, "c", x=-0.2, y=1.12)
    data = pd.concat(rows, ignore_index=True)
    data["q_family"] = np.where(data.panel == "c", "BH across C1 tests within each model",
                                "BH across S1-S6 within each model")
    save(fig, OUT, "X3_status_probability_by_subtype", data=data, caption=CAPTIONS["X3"](d, comp))


# ============================================================ X4 exon protein character
def X4() -> None:
    pr = pd.read_csv(IN / "exon_protein_character.csv")
    sets = ["all measured exons (random 6,000)", "tumour vs normal, top events", "cancer type, top events",
            "subtype, top events"]
    short = dict(zip(sets, ["all measured exons", "tumour status, top", "cancer type, top", "subtype, top"]))
    cols = dict(zip(sets, [MUTED, STATUS, TYPE, SUBTYPE]))
    metrics = list(dict.fromkeys(pr["metric"]))
    among = pr.drop_duplicates("metric").set_index("metric")["among"]
    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    for j, st in enumerate(sets):
        sub = pr[pr.set == st].set_index("metric").loc[metrics]
        y = np.arange(len(metrics)) + (j - 1.5) * 0.19
        ax.hlines(y, sub.ci_lo, sub.ci_hi, color=cols[st], linewidth=0.8)
        ax.scatter(sub.proportion, y, s=9, color=cols[st], linewidths=0, zorder=3,
                   label=f"{short[st]} (n = {int(sub.n.iloc[0]):,})")
    names = {"annotated_exon": "exact GENCODE", "fully_coding": "fully coding"}
    ax.set_yticks(np.arange(len(metrics)))
    ax.set_yticklabels([m if among[m] == "all" else f"{m}\n(of {names.get(among[m], among[m])} exons)"
                        for m in metrics], fontsize=5.5, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("share of cassette exons (95% Wilson CI)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, handletextpad=0.1, fontsize=5.6, columnspacing=0.8)
    data = pr.copy()
    data["q_family"] = "BH across 6 metrics x 3 comparison sets (Fisher tests vs all measured exons)"
    save(fig, OUT, "X4_exon_protein_character", data=data, caption=CAPTIONS["X4"](pr))


# ============================================================ captions
def _cap_x1(d, tests):
    t = tests.set_index("pair")
    g2 = tests[tests.id == "G2"].iloc[0]
    def line(k, name):
        r = t.loc[k]
        return (f"- {name}: median {r.median_fold_enrichment:.1f}-fold the chance overlap over {int(r.n)} comparisons, "
                f"{100 * r.share_q05:.0f}% with q < 0.05 (Wilcoxon of log2 fold vs 0, p = {fp(r.p_wilcoxon_log2fe_vs_0)}).")
    return f"""**Do the three tasks' models read the same genes?**

For each subtype cohort, its cancer's tumour-status model (12 cohorts in 8 cancers have one) and the cancer-type
model's importance over that cancer's own samples: the top 100 genes of each model per view (a gene ranked by its
best feature's mean |SHAP|), compared within the genes both tasks could select (D-48). a, each dot is one cohort and
view (comparisons in brackets); bar = median; 0 = the overlap expected by chance (hypergeometric). b, the five
histology cohorts that have a tumour-status model: Jaccard of the subtype list with the status list and with the
type list; small dots = views, large = median over views.

{line('status-type', 'Status – type (G1)')}
{line('subtype-status', 'Subtype – status (G1)')}
{line('subtype-type', 'Subtype – type (G1)')}
- G2, for histology subtypes the subtype list is closer to type than to status: median Jaccard difference
  {signed(g2.median_difference, 3)}, positive in {100 * g2.share_positive:.0f}% of {int(g2.n)} cohort-views, Wilcoxon p = {fp(g2.p_wilcoxon)}:
  not supported.

Genes in all three top-100 lists of a cancer are in `results/across_tasks/genes_in_all_three.csv`."""


def _verdicts(d, model):
    v = d[d.model == model].verdict.value_counts()
    return ", ".join(f"{int(v.get(k, 0))} {k}" for k in ["supported", "against", "not resolved"])


def _cap_x2(d, sens):
    s6 = sens[(sens.model == "GEX + AS") & (sens.id == "T6")].iloc[0]
    s5 = sens[(sens.model == "AS ensemble") & (sens.id == "T5")].iloc[0]
    return f"""**Which other cancers does a subtype resemble, in the cancer-type model?**

The pan-cancer type model gives every tumour a probability for each of 34 classes. Thirteen predictions (D-49,
sources in config/followup.yaml) name a subtype, the projects it should resemble and the comparison class(es):
row = id, cohort, the predicted-higher vs predicted-lower group, → target projects. Measure: the patient's summed
out-of-fold probability over the targets; two-sided Mann–Whitney, AUROC = P(higher-group patient > lower-group
patient). a, the mean of that probability in each group, GEX + AS (log scale). b, AUROC per model; filled = BH q <
0.05 within the model. T11 (cortical admixture in PCPG → ACC) is a positive control.

- GEX: {_verdicts(d, 'GEX')}. AS ensemble: {_verdicts(d, 'AS ensemble')}. GEX + AS: {_verdicts(d, 'GEX + AS')}.
- Sensitivity, not pre-specified (×): the targets' share of the probability not given to the patient's own
  project. Most predictions hold. T6 reverses (GEX + AS AUROC {s6.auroc:.2f}): basal-like tumours carry more
  ovarian probability because they carry more off-target probability overall, not a larger ovarian share. T5 in
  the AS ensemble weakens to {s5.auroc:.2f}.

Group sizes, means, medians, raw p and q for every model and measure are in the CSV."""


def _cap_x3(d, comp):
    c = comp[comp.model == "GEX + AS"]
    neg = int(((c.rho < 0) & (c.q < 0.05)).sum())
    pos = int(((c.rho > 0) & (c.q < 0.05)).sum())
    pos_c = " and ".join(sorted(c[(c.rho > 0) & (c.q < 0.05)].project.unique()))
    g = d[d.model == "GEX"]
    return f"""**Does the tumour-status model see some subtypes as less tumour-like?**

The tumour-status model of each cancer scores every sample, out of fold. For the subtype cohorts of the 8 cancers
with a status model, six predictions (D-50) name the subtype expected to look less tumour-like (row: id, cohort,
predicted lower vs comparison). a, mean 1 − P(tumour) in the predicted group (filled) and the comparison (open),
for GEX and the AS ensemble (log scale). b, AUROC = P(comparison patient's P(tumour) > predicted patient's),
two-sided Mann–Whitney; filled = BH q < 0.05 within the model. c, C1: Spearman ρ of P(tumour) with the stromal and the
immune score (mean log2 CPM of the D-18 marker genes) over each cancer's tumours, all of them, not only the
subtype-labelled ones.

- GEX: {_verdicts(d, 'GEX')}. AS ensemble: {_verdicts(d, 'AS ensemble')}. GEX + AS: {_verdicts(d, 'GEX + AS')}.
- GEX's probabilities sit at the ceiling (median P(tumour) of the comparison groups {g.median_higher.min():.4f}–{g.median_higher.max():.4f}), so
  ranks among them carry little; the AS ensemble averages eight models and spreads them out.
- C1, P(tumour) falls as composition scores rise: in GEX + AS {neg} of {len(c)} correlations are negative with q < 0.05
  and {pos} positive (in {pos_c}): not supported as a common explanation.

Group sizes, means, medians, raw p and q are in the CSV; the per-class distributions in
`status_probability_by_subtype_class.csv`."""


def _cap_x4(pr):
    s = pr[pr.set == "subtype, top events"].set_index("metric")
    b = pr[pr.set == "all measured exons (random 6,000)"].set_index("metric")
    def m(k):
        return (f"{100 * s.at[k, 'proportion']:.0f}% against {100 * b.at[k, 'proportion']:.0f}% "
                f"(OR {s.at[k, 'odds_ratio']:.2f}, q = {fp(s.at[k, 'q_fisher'])})")
    return f"""**Are the exons each task uses alike in what they do to the protein?**

Cassette exons in the top 100 of any tumour-status model (blue), any cancer-type class (dark blue) or any subtype
cohort's SE model (amber; D-51), against the random 6,000 measured exons of D-17 (grey). Length divisible by 3 from
the coordinates; the rest need an exact GENCODE v45 exon and use the best transcript (MANE > canonical > other)
with InterPro Pfam/CDD domains, MobiDB-lite disorder and TMHMM helices (D-17). Two-sided Fisher exact tests against
the background; BH q over 6 metrics × 3 sets.

E1, the subtype exons share the bias of the other two tasks: divisible by 3 {m('length divisible by 3')}; in a
domain {m('overlaps a Pfam/CDD domain')}; in a disordered region {m('overlaps a disordered region')}; in a TM helix
{m('overlaps a TM helix')}. Supported."""


CAPTIONS = {k: (lambda f: (lambda *a: tidy(f(*a))))(f) for k, f in
            {"X1": _cap_x1, "X2": _cap_x2, "X3": _cap_x3, "X4": _cap_x4}.items()}

if __name__ == "__main__":
    run({"X1": X1, "X2": X2, "X3": X3, "X4": X4})
