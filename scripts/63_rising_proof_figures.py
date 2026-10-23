"""Step 63: D-77's step-wise test of the rising arm, in figures (house style, docs/FIGURE_GUIDELINES.md).

One question per figure, in the order of the argument in results/MSI_RISING_ARM_PROOF.md:
    RP0 the ladder: every step's test, its effect and its verdict, and the claim's verdict
    RP1 do models use the rising arm (subtype MLP SHAP, the arm against matched exon sets, M4)?
    RP2 does the rising arm rise in junction reads, and how is the rise spread over its exons?
    RP3 does the rise need the T-run, and does it grow with the gap to the mapped branch point?
    RP4 is it the microsatellite, or would any long run or Alu acceptor do?
    RP5 do exon and tumour covariates explain it?
    RP6 is it mismatch repair, or mutation load, or any subtype difference?
    RP7 is it there within the same patient (tumour minus its own normal)?
    RP8 does it follow dose (MANTIS score, run length)?
    RPM which mechanism (branch-point geometry or hnRNP C)?

Reads only what step 62 wrote with the real labels (results/subtype/msi_cis/code/proof/: tests.csv,
sensitivity.csv, descriptive.csv, verdict.json, rp6c_null_contrasts.csv) and its local caches
(data.nosync/.../msi/proof/tests/: builds, engines, tumours, pairs). Nothing is re-tested. Where a figure
draws a distribution step 62 summarised (RP1's matched-set and random-removal AUROCs, RPD1's permutation
null), the distribution is recomputed with step 62's own code and seeds, and the script stops unless it
reproduces the committed numbers. The 95% intervals of the control classes (RP2, RP4, RPM) are
descriptive, from RP2a's construction (studentized SE, Welch-Satterthwaite df); they are not tests.

Captions are figure legends ("**Figure RPn. Title.**", then (a), (b) ...); `legends` copies each legend
(without its bullet list of tests) under the figure in results/MSI_RISING_ARM_PROOF.md. Per-tumour plotted
values go to data.nosync/ctsm_v1/figure_data/<stem>.csv (D-12); committed CSVs hold aggregates.
Usage: python scripts/63_rising_proof_figures.py [--only RP0 RP1 ... legends]     (default: all, then legends)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pickle
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _common import log
from _figbase import LOCAL, inject_legends, panel_letters, tidy
from _figbase import legend as _legend
from matplotlib.lines import Line2D

from ctsm.figstyle import ACCENT, ACCENT_DARK, ACCENT_LIGHT, DOUBLE, GRID, INK, INK2, MUTED, fp, save


def _load62():
    spec = importlib.util.spec_from_file_location("proof62", Path(__file__).with_name("62_msi_rising_proof.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


P = _load62()
CODE = P.OUT                                         # results/subtype/msi_cis/code/proof (real labels)
OUT = Path(__file__).resolve().parents[1] / "results" / "figures" / "rising_proof"
WRITEUP = Path(__file__).resolve().parents[1] / "results" / "MSI_RISING_ARM_PROOF.md"
CACHE = LOCAL / "rising_proof_cache"
TESTS = pd.read_csv(CODE / "tests.csv").set_index("id")
SENS = pd.read_csv(CODE / "sensitivity.csv")
DESC = pd.read_csv(CODE / "descriptive.csv")
VERDICT = json.loads((CODE / "verdict.json").read_text())
FAM_N = TESTS.family.value_counts()
MINUS = "−"
LIGHT = "#c9c7c0"
MS = 3.6                                             # marker size (points)


# ============================================================ small helpers
def T(i: str, col: str = "estimate") -> float:
    return float(TESTS.at[i, col])


def sg(v: float, nd: int = 4) -> str:
    return f"{v:+.{nd}f}".replace("-", MINUS)


def num(v: float, nd: int = 4) -> str:
    return f"{v:.{nd}f}".replace("-", MINUS)


def vd(i: str) -> str:
    return TESTS.at[i, "verdict"]


def tline(i: str, what: str, nd: int = 4) -> str:
    return f"{i}, {what}: {sg(T(i), nd)}, p = {fp(T(i, 'p'))}, {vd(i)}."


def qcols(i: str) -> dict:
    """The test's q with its FDR family named (CSV only)."""
    f = TESTS.at[i, "family"]
    return {"test": i, "estimate": T(i), "p": T(i, "p"), "q": T(i, "q"),
            "q_family": f"BH within family {f} ({FAM_N[f]} tests; D-77)", "verdict": vd(i)}


def desc(item: str, group=None) -> pd.Series:
    d = DESC[DESC["item"].eq(item)]
    if group is not None:
        d = d[d["group"].astype(str).eq(group)]
    assert len(d) == 1, (item, group, len(d))
    return d.iloc[0]


def sens(item: str, test: str) -> pd.Series:
    d = SENS[SENS["item"].eq(item) & SENS["test"].eq(test)]
    assert len(d) == 1, (item, test, len(d))
    return d.iloc[0]


def same(a: float, b: float, what: str, tol: float = 1e-6) -> None:
    """Stop unless a value redrawn here equals the committed one."""
    if not (abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))):
        raise SystemExit(f"{what}: redrawn {a!r} differs from the committed {b!r}")


def cached(name: str, fn):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{name}.pkl"
    if f.exists():
        return pickle.loads(f.read_bytes())
    v = fn()
    f.write_bytes(pickle.dumps(v))
    return v


def axis_clean(ax, left: bool = True) -> None:
    ax.tick_params(axis="y", length=0 if not left else 2.5)


def zero(ax, horizontal: bool = False, **kw) -> None:
    (ax.axhline if horizontal else ax.axvline)(0, color=GRID, linewidth=0.7, zorder=0, **kw)


# ============================================================ shared state (step 62's caches)
_ST: dict = {}


def state():
    """Step 62's State with the real labels, its primary build and held-out exon engine (cached on disk)."""
    if not _ST:
        S = P.State(P.Run("real", True))
        prim = P.SPECS["primary"]
        S.build(prim, S.projects)
        S.tv = P.tumour_values(S)
        H = S.xctx(prim, "held_out")
        same(np.nanmean(H.D[0, H.tab.cls.eq("R_up").to_numpy() & H.elig]), T("RP2a"), "RP2a")
        _ST.update(S=S, H=H)
    return _ST["S"], _ST["H"]


def delta0() -> np.ndarray:
    _, H = state()
    return H.D[0].astype(float)


def sel_of(name: str) -> np.ndarray:
    """Eligible held-out exons of a class (step 62's definitions)."""
    _, H = state()
    t = H.tab
    m = {"R_up": t.cls.eq("R_up").to_numpy(), "R_down": t.cls.eq("R_down").to_numpy(),
         "B sample": t.b_sample.to_numpy(bool), "B_far": t.b_sample.to_numpy(bool) & t.cls.eq("B_far").to_numpy(),
         "ppt_far": t.ctl_ppt_far.to_numpy(bool), "alu": t.ctl_alu.to_numpy(bool), "polyA15": t.ctl_polyA15.to_numpy(bool),
         "deep": t.ctl_deep.to_numpy(bool), "donor": t.ctl_don.to_numpy(bool), "shield": t.shield.to_numpy(bool),
         "long, antisense Alu": (t.acc_T >= P.LONG).to_numpy() & t.run_alu_antisense.eq(True).to_numpy(),
         "long, other": (t.acc_T >= P.LONG).to_numpy() & ~t.run_alu_antisense.eq(True).to_numpy()}[name]
    return m & H.elig & np.isfinite(H.D[0])


def class_interval(name: str) -> dict:
    """Mean delta of a class with a descriptive 95% interval (RP2a's construction)."""
    def go():
        S, H = state()
        r = P.class_mean(S, H, sel_of(name))
        return {"class": name, "n_exons": int(r["n_exons_obs"]), "mean": float(r["T"][0]),
                "ci_low": float(r["ci"][0]), "ci_high": float(r["ci"][1])}
    return cached(f"class_{name.replace(' ', '_').replace(',', '')}", go)


CLASS_NAME = {"R_up": "rising arm (R_up)", "R_down": "falling arm (R_down)", "B sample": "background (B sample)",
              "B_far": "distant branch point, no T-run", "ppt_far": "long pyrimidine tract", "alu": "Alu acceptor, no T-run",
              "polyA15": "poly-A run", "deep": "deep-intron T-run", "donor": "donor-side T-run",
              "shield": "T-run upstream of the branch point", "long, antisense Alu": "run ≥ 20 nt in an antisense Alu",
              "long, other": "run ≥ 20 nt, not in an antisense Alu"}


# ============================================================ RP0 the ladder
LADDER = [
    ("1  Do models use the rising arm?", [
        ("RP1", "R_up arm against matched exon sets", "auc", "RPC1"),
        ("RPD1a", "subtype MLPs push R_up towards MSI", "auc", "RPD1b")]),
    ("2  Does it rise in junction reads?", [
        ("RP2a", "R_up exons, mean delta (primary test)", "psi", "RPC2a"),
        ("RP2c", "share rising, R_up − matched exons", "auc", None)]),
    ("3  Does it need both features?", [
        ("RP3a", "R_up − distant branch point, no run", "psi", None),
        ("RP3b", "R_up − long run, branch point < 36 nt", "psi", None),
        ("RP3d", "mapped branch point ≥ 36 − < 36 nt", "psi", None),
        ("RP3e", "rise vs gap to mapped branch point (ranks)", "rho", None),
        ("RPD3b", "the same, training cancers", "rho", None),
        ("RP3c", "rise vs gap to heuristic branch point (ranks)", "rho", None)]),
    ("4  Is it the microsatellite?", [
        ("RP4b", "R_up − long pyrimidine tract", "psi", None),
        ("RP4e", "R_up − Alu acceptor without a T-run", "psi", None),
        ("RP4a", "R_up − poly-A run", "psi", None),
        ("RP4c", "R_up − deep-intron T-run", "psi", None),
        ("RP4d", "R_up − donor-side T-run", "psi", None)]),
    ("5  Does it survive confounders?", [
        ("RP5a", "exon covariates: R_up coefficient", "psi", "RPC5a"),
        ("RP5b", "tumour covariates: MMRd coefficient", "psi", "RPC5b")]),
    ("6  Is it mismatch repair?", [
        ("RP6a", "MMRd vs hypermutated MMRp tumours", "auc", "RPC6"),
        ("RP6b", "MANTIS MSI-H vs MSS tumours", "psi", None),
        ("RP6c", "MSI against 57 other subtype contrasts", "psi", None)]),
    ("7  Within the same patient?", [
        ("RP7", "tumour − own normal, MMRd vs MMRp", "psi", "RPC7")]),
    ("8  Does it follow dose?", [
        ("RP8a", "R_up arm vs MANTIS score", "rho", "RPC8"),
        ("RP8b", "R_up arm vs deletion count", "rho", None),
        ("RP8c", "R_up rise vs run length", "rho", None)]),
    ("M  Which mechanism? (outside the verdict)", [
        ("RPM1", "R_up − T-run upstream of branch point", "psi", None),
        ("RPM3", "runs in antisense Alu − other runs", "psi", None)]),
]
UNITS = {"psi": "PSI difference", "rho": "Spearman ρ", "auc": "AUROC or share,\ndifference"}
UNIT_LIM = {"psi": (-0.007, 0.023), "rho": (-0.38, 0.55), "auc": (-0.32, 0.45)}


def _marker(ax, x, y, verdict, colour, marker="o", size=MS, z=4):
    fill = colour if verdict == "supported" else "white"
    ax.plot([x], [y], marker=marker, markersize=size, markerfacecolor=fill, markeredgecolor=colour,
            markeredgewidth=0.8, linestyle="", zorder=z)


def RP0() -> None:
    rows, y = [], 0
    for head, items in LADDER:
        rows.append(("head", head, None, None, None, y))
        y += 1
        for i, q, u, c in items:
            rows.append(("test", q, i, u, c, y))
            y += 1
    n = y
    fig = plt.figure(figsize=(DOUBLE - 0.1, 0.165 * n + 0.9))
    gs = fig.add_gridspec(1, 5, width_ratios=[2.35, 1.45, 0.85, 0.85, 0.9], wspace=0.1)
    lab = fig.add_subplot(gs[0])
    axes = {u: fig.add_subplot(gs[k + 1], sharey=lab) for k, u in enumerate(["psi", "rho", "auc"])}
    txt = fig.add_subplot(gs[4], sharey=lab)
    data = []
    for kind, text, i, u, c, yy in rows:
        if kind == "head":
            lab.text(0, yy, text, fontsize=6.4, fontweight="bold", va="center", ha="left", color=INK)
            for ax in axes.values():
                if yy:
                    ax.axhline(yy - 0.5, color=GRID, linewidth=0.5, zorder=0)
            continue
        fam = TESTS.at[i, "family"]
        lab.text(0.03, yy, i, fontsize=5.8, va="center", ha="left", color=INK,
                 fontweight="bold" if fam == "RP" else "normal")
        lab.text(0.185, yy, text, fontsize=5.8, va="center", ha="left", color=INK2)
        ax = axes[u]
        est, lo, hi = T(i), T(i, "ci_low"), T(i, "ci_high")
        if np.isfinite(lo):
            ax.plot([lo, hi], [yy, yy], color=ACCENT, linewidth=1.0, zorder=3, solid_capstyle="butt")
        _marker(ax, est, yy, vd(i), ACCENT)
        data.append({"row": text, "kind": "rising arm", **qcols(i), "ci_low": lo, "ci_high": hi, "unit": UNITS[u]})
        if c is not None:
            ec, clo, chi = T(c), T(c, "ci_low"), T(c, "ci_high")
            if np.isfinite(clo):
                ax.plot([clo, chi], [yy, yy], color=MUTED, linewidth=1.0, zorder=2, solid_capstyle="butt")
            _marker(ax, ec, yy, vd(c), MUTED, marker="D", size=MS - 0.6, z=3)
            data.append({"row": text, "kind": "falling-arm control", **qcols(c), "ci_low": clo, "ci_high": chi,
                         "unit": UNITS[u]})
        if i == "RP2a":
            ax.plot([P.SESOI, P.SESOI], [yy - 0.38, yy + 0.38], color=INK2, linewidth=0.8, zorder=2)
            data.append({"row": text, "kind": "reference: SESOI", "estimate": P.SESOI, "unit": UNITS[u]})
        ptxt = fp(T(i, "p"))
        txt.text(0.02, yy, ptxt, fontsize=5.8, va="center", ha="left", color=INK)
        txt.text(0.42, yy, vd(i), fontsize=5.8, va="center", ha="left",
                 color=ACCENT_DARK if vd(i) == "supported" else MUTED)
    for u, ax in axes.items():
        zero(ax)
        ax.set_xlim(*UNIT_LIM[u])
        ax.set_xlabel(UNITS[u], fontsize=6)
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.tick_params(axis="x", labelsize=5.6)
        ax.spines["left"].set_visible(False)
        data.append({"row": "reference line", "kind": "0", "unit": UNITS[u], "estimate": 0.0})
    axes["psi"].set_xticks([0, 0.01, 0.02])
    axes["rho"].set_xticks([-0.3, 0, 0.3])
    axes["auc"].set_xticks([-0.2, 0, 0.2, 0.4])
    for ax in (lab, txt):
        ax.axis("off")
    lab.set_xlim(0, 1)
    txt.set_xlim(0, 1)
    lab.set_ylim(n - 0.4, -1.6)
    txt.text(0.02, -1.0, "p", fontsize=5.8, va="center", ha="left", color=INK2, fontstyle="italic")
    txt.text(0.42, -1.0, "verdict", fontsize=5.8, va="center", ha="left", color=INK2)
    hand = [Line2D([], [], marker="o", color=ACCENT, linestyle="", markersize=MS, label="supported"),
            Line2D([], [], marker="o", color=ACCENT, markerfacecolor="white", linestyle="", markersize=MS,
                   label="not resolved"),
            Line2D([], [], marker="D", color=MUTED, linestyle="", markersize=MS - 0.6, label="falling-arm control"),
            Line2D([], [], color=INK2, linewidth=0.8, label="SESOI (0.003)")]
    axes["psi"].legend(handles=hand, loc="lower left", bbox_to_anchor=(-0.05, 1.0), ncol=4, fontsize=5.6,
                       handletextpad=0.3, columnspacing=1.0, borderaxespad=0.1)
    fig.text(0.02, -0.01, f"Pre-registered verdict for the claim: {VERDICT['verdict']}. Not met: {unmet_short()}.",
             fontsize=6.4, ha="left", va="top", color=INK, fontweight="bold")
    data.append({"row": "claim verdict", "kind": VERDICT["verdict"], "unit": "", "test": "; ".join(VERDICT["reasons"])})
    save(fig, OUT, "RP0_ladder", data=pd.DataFrame(data), caption=CAPTIONS["RP0"]())


def omitted_tests() -> list:
    """Tests of tests.csv not drawn in the ladder (neither as a row nor as a falling-arm control)."""
    shown = {i for _, items in LADDER for i, _, _, c in items} | {c for _, items in LADDER for _, _, _, c in items if c}
    return [i for i in TESTS.index if i not in shown]


def unmet_short() -> str:
    """The claim verdict's unmet conditions, shortened when they all read '<id> not resolved (needed: supported)'."""
    r = VERDICT["reasons"]
    if r and all(x.endswith("not resolved (needed: supported)") for x in r):
        ids = [x.split()[0] for x in r]
        return f"{', '.join(ids[:-1])} and {ids[-1]} not resolved (needed: supported)" if len(ids) > 1 else r[0]
    return "; ".join(r)


def _cap_rp0() -> str:
    fams = VERDICT["families"]
    fs = lambda f: f"{fams[f]['supported']} of {fams[f]['tests']}"
    return f"""**Does the rising arm pass each step of the pre-registered test?**

Pre-registered as D-77 (config/followup.yaml: msi_rising_proof; commit 5471b27) before any junction-read value was
compared by MMR status. Confirmatory tests use the 30 held-out cancers (50 MMR-deficient tumours in 16 of them) and
junction reads only (rMATS JC), in PSI units. delta = an exon's mean adjusted PSI in MMR-deficient (MMRd) tumours
minus the proficient (MMRp) mean of the same cancer. One row per test shown, in the order of the argument; the id is
bold for the 12 tests the verdict needs (family RP). RPD3b is the gap test of RP3e repeated in the training cancers.
Not shown: {', '.join(omitted_tests())} (all in tests.csv). Effects in different units sit in different columns: PSI
differences (class contrasts, regression coefficients, arm contrasts), Spearman ρ (rank associations and dose), and AUROC
or share differences (RP1: AUROC minus the median of matched exon sets; RPD1a: share of positive SHAP directions
minus 0.5; RP2c: share of R_up exons rising minus that of matched exons; RP6a: AUROC minus 0.5). Filled = supported
(predicted sign and q < 0.05 within the family), hollow = not resolved; bars = 95% intervals where the test has one
(RP2a and RPC2a: studentized, Welch–Satterthwaite). Grey diamonds = the same test on the falling arm (runs of 10–19
nt), predicted below 0. p = raw p of the test (the larger of the exon and the label permutation p where both
exist); q and families in RP0_ladder.csv. The line under the panel gives the claim's verdict by the pre-registered
rule.

- Family RP (verdict): {fs('RP')} supported; RPS (secondary held out): {fs('RPS')}; RPP (not held out): {fs('RPP')};
  RPC (falling arm): {fs('RPC')}; RPD (training, discovery): {fs('RPD')}; RPM (mechanism): {fs('RPM')}.
- RP2a: {sg(T('RP2a'))} (95% interval {num(T('RP2a', 'ci_low'))} to {num(T('RP2a', 'ci_high'))}), p = {fp(T('RP2a', 'p'))};
  the whole interval lies above the SESOI of 0.003. RPC2a: {sg(T('RPC2a'))} ({num(T('RPC2a', 'ci_low'))} to
  {num(T('RPC2a', 'ci_high'))}), p = {fp(T('RPC2a', 'p'))}.
- Not resolved in family RP: RP4a (poly-A) {sg(T('RP4a'))}, p = {fp(T('RP4a', 'p'))}; RP4d (donor-side T-run)
  {sg(T('RP4d'))}, p = {fp(T('RP4d', 'p'))}; RP6c {sg(T('RP6c'))}, p = {fp(T('RP6c', 'p'))}.
- Claim verdict: {VERDICT['verdict']}; not met: {unmet_short()}; no refuting clause holds."""


# ============================================================ RP1 models
def _rpd1_folds():
    """RPD1's folds, built exactly as step 62's rpd1() does with the real labels."""
    cls = P.load_annotation().set_index("feature_id").cls
    folds = []
    for cohort in ["COAD_TCGA", "STAD_TCGA", "UCEC_TCGA"]:
        rd = P.m51.DATA / "subtype" / "runs" / "subtype" / "se" / cohort
        tab = pd.read_parquet(P.m51.DATA / "subtype" / "tasks" / f"subtype__{cohort}.parquet")
        k_msi = P.msi_class_index(json.loads((rd / "job_done.json").read_text()), tab)
        y = tab.y.to_numpy()
        for k in range(5):
            z = np.load(rd / f"fold{k}" / "shap.npz")
            msi = np.flatnonzero(y[z["rows"]] == k_msi)
            fids = z["feature_ids"].astype(str)
            arm = cls.reindex(fids).fillna("").to_numpy()
            j = np.flatnonzero(np.isin(arm, ["R_up", "R_down"]))
            folds.append({"cohort": cohort, "X": z["X"].astype(np.float64)[msi][:, j],
                          "S": z["shap"].astype(np.float64)[msi][:, j], "fids": fids[j], "arm": arm[j]})
    return folds


def rpd1_null() -> dict:
    def go():
        folds = _rpd1_folds()
        obs, _ = P.direction_shares(folds)
        rng = np.random.default_rng(P.SEED)
        null = {a: np.empty(P.RPD1_PERM) for a in obs}
        for b in range(P.RPD1_PERM):
            sh, _ = P.direction_shares(folds, rng)
            for a, (k, n_) in sh.items():
                null[a][b] = k / n_
        return {"obs": obs, "null": null}
    r = cached("rpd1_null", go)
    for a, i in [("R_up", "RPD1a"), ("R_down", "RPD1b")]:
        k, n_ = r["obs"][a]
        same(k / n_ - 0.5, T(i), i)
        same(P.perm_p(k / n_ - 0.5, r["null"][a] - 0.5), T(i, "p"), f"{i} p")
    return r


def rp1_sets() -> dict:
    def go():
        S, _ = state()
        held = S.lab.held_out.to_numpy()
        out = {}
        for arm in ["R_up", "R_down"]:
            Sc = np.full((len(S.lab), P.N_SETS + 1), np.nan, np.float32)
            for p in S.held:
                r = S.res["primary"][p]
                if f"rp1_{arm}" in r:
                    Sc[S.rows_of(P.SPECS["primary"], p)] = r[f"rp1_{arm}"]
            out[arm] = P.auroc_matrix(Sc[held], S.Y[held][:, :1])[:, 0]
        return out
    r = cached("rp1_sets", go)
    for arm, i in [("R_up", "RP1"), ("R_down", "RPC1")]:
        same(r[arm][0] - np.median(r[arm][1:]), T(i), i)
    return r


def _auc_without(path54: str, X, tr1, te1, y, p4, drops: list) -> list:
    """Step 62's _m4_auc_without, importable by joblib's workers (step 54 loaded from its file)."""
    import importlib.util as iu
    spec = iu.spec_from_file_location("msi54_fig", path54)
    m54 = iu.module_from_spec(spec)
    spec.loader.exec_module(m54)
    out = []
    for d in drops:
        keep = np.setdiff1d(np.arange(X.shape[1]), d)
        out.append(m54.auc(m54.fit_predict("M4", p4, X[:, keep], tr1, te1, y), y[te1]))
    return out


def m4_removals() -> dict:
    """A-25's random-removal reference, refitted with step 62's code, data and seed."""
    def go():
        from joblib import Parallel, delayed
        m54 = P._load("msi54", "54_msi_models.py")
        S, _ = state()
        t = pd.read_parquet(P.LOCAL / "models_tumours.parquet")
        cols = pd.read_parquet(P.LOCAL / "models_exon_columns.parquet")
        M = np.load(P.LOCAL / "models_exon_matrix.npy", mmap_mode="r")
        y = S.lab.set_index("file_id").mmr.reindex(t.file_id).to_numpy().astype(bool)
        held = ~t.project.isin(P.TRAIN).to_numpy()
        tr1, te1 = np.flatnonzero(~held), np.flatnonzero(held)
        at, db = cols.acc_T.to_numpy(float), cols.d_bp.to_numpy(float)
        _, up, _ = P.m51.arms(at, db)
        rep = at >= 10
        p4 = m54._chosen(pd.read_csv(m54.OUT / "tuning_e1.csv"), "M4")
        X = np.asarray(M[:, np.flatnonzero(rep)], np.float64)
        au1 = m54.auc(m54.fit_predict("M4", p4, X, tr1, te1, y), y[te1])
        keep = np.flatnonzero(~up[np.flatnonzero(rep)])
        au2 = m54.auc(m54.fit_predict("M4", p4, X[:, keep], tr1, te1, y), y[te1])
        others = np.flatnonzero(~up[np.flatnonzero(rep)])
        n_up = int(up[rep].sum())
        rng_d = np.random.default_rng(P.SEED)
        drops = [rng_d.choice(others, n_up, replace=False) for _ in range(P.M4_RANDOM_DRAWS)]
        chunks = np.array_split(np.arange(P.M4_RANDOM_DRAWS), P.N_WORKERS)
        parts = Parallel(n_jobs=P.N_WORKERS)(delayed(_auc_without)(str(Path(P.__file__).with_name("54_msi_models.py")),
                                                                   X, tr1, te1, y, p4, [drops[i] for i in ch])
                                             for ch in chunks)
        au_r = np.concatenate([np.asarray(a, float) for a in parts])
        return {"au1": au1, "au2": au2, "random": au_r - au1, "n_up": n_up}
    r = cached("m4_removals", go)
    same(r["au1"], desc("M4 refit: held-out AUROC").value, "M4 AUROC")
    same(r["au2"] - r["au1"], desc("M4 without its R_up inputs: AUROC change").value, "M4 change")
    rr = desc(f"M4 without {r['n_up']} random non-R_up repeat inputs: AUROC change (exploratory, A-25)")
    same(np.median(r["random"]), rr.value, "M4 random median")
    return r


def _hist(ax, x, bins, colour=LIGHT):
    h, e = np.histogram(x, bins=bins)
    ax.bar(e[:-1], h, width=np.diff(e), align="edge", color=colour, edgecolor="white", linewidth=0.3)
    return h, e


def RP1() -> None:
    sh, sets, m4 = rpd1_null(), rp1_sets(), m4_removals()
    fig, axs = plt.subplots(1, 4, figsize=(DOUBLE - 0.1, 1.85), gridspec_kw={"wspace": 0.55,
                                                                           "width_ratios": [0.8, 1, 1, 1]})
    rows = []
    # a. SHAP direction shares with the permutation null
    ax = axs[0]
    for x, (a, i) in enumerate([("R_up", "RPD1a"), ("R_down", "RPD1b")]):
        k, n_ = sh["obs"][a]
        lo, hi = np.percentile(sh["null"][a], [2.5, 97.5])
        ax.add_patch(plt.Rectangle((x - 0.36, lo), 0.72, hi - lo, color=GRID, zorder=0, linewidth=0))
        ax.bar(x, k / n_, width=0.46, color=ACCENT if a == "R_up" else MUTED, zorder=2)
        ax.text(x, k / n_ + 0.03, f"{k}/{n_}", ha="center", va="bottom", fontsize=5.6, color=INK)
        rows.append({"panel": "a", "arm": a, "positive": k, "units": n_, "share": k / n_, "null_2.5": lo,
                     "null_97.5": hi, "null_sd": float(np.std(sh["null"][a])), **qcols(i)})
    ax.axhline(0.5, color=INK2, linewidth=0.6, linestyle=(0, (2, 2)), zorder=1)
    rows.append({"panel": "a", "arm": "reference line", "share": 0.5})
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["R_up", "R_down"], fontsize=6)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 1)
    ax.set_ylabel("share of inputs pushing\ntowards MSI", fontsize=6)
    # b, c. the arm's held-out AUROC against 1,000 matched exon sets
    for ax, (arm, i, col) in zip(axs[1:3], [("R_up", "RP1", ACCENT), ("R_down", "RPC1", INK2)]):
        a = sets[arm]
        h, e = _hist(ax, a[1:], np.linspace(0.3, 0.75, 46))
        ax.axvline(a[0], color=col, linewidth=1.2, zorder=3)
        ax.axvline(0.5, color=INK2, linewidth=0.5, linestyle=(0, (2, 2)), zorder=1)
        right = a[0] < 0.5
        ax.text(a[0] + (0.012 if right else -0.012), h.max() * 1.05, f"{a[0]:.2f}", color=col, fontsize=5.8,
                ha="left" if right else "right", va="bottom")
        ax.set_xlim(0.28, 0.77)
        ax.set_ylim(0, h.max() * 1.18)
        ax.set_xlabel(f"held-out AUROC, {arm} arm", fontsize=6)
        ax.set_ylabel("matched exon sets", fontsize=6)
        for lo_, n_ in zip(e[:-1], h):
            if n_:
                rows.append({"panel": "b" if arm == "R_up" else "c", "arm": arm, "bin_left": lo_, "bin_width": e[1] - e[0],
                             "matched_sets": int(n_)})
        rows.append({"panel": "b" if arm == "R_up" else "c", "arm": arm, "class_auroc": a[0],
                     "matched_median": float(np.median(a[1:])), **qcols(i)})
    # d. M4 without its R_up inputs against 200 random removals
    ax = axs[3]
    h, e = _hist(ax, m4["random"], np.linspace(-0.02, 0.012, 33))
    ch = m4["au2"] - m4["au1"]
    ax.axvline(ch, color=ACCENT, linewidth=1.2, zorder=3)
    zero(ax)
    ax.text(ch + 0.0007, h.max() * 1.05, f"{ch:+.4f}".replace("-", MINUS), color=ACCENT, fontsize=5.8, ha="left",
            va="bottom")
    ax.set_ylim(0, h.max() * 1.18)
    ax.set_xlabel("AUROC change, M4 refitted", fontsize=6)
    ax.set_ylabel(f"random removals of {m4['n_up']}", fontsize=6)
    ax.set_xticks([-0.02, -0.01, 0, 0.01])
    for lo_, n_ in zip(e[:-1], h):
        if n_:
            rows.append({"panel": "d", "arm": "random non-R_up inputs", "bin_left": lo_, "bin_width": e[1] - e[0],
                         "draws": int(n_)})
    pct = desc("M4 without its R_up inputs: percentile of the AUROC change among random removals (A-25)").value
    rows.append({"panel": "d", "arm": "R_up inputs removed", "auroc_full": m4["au1"], "auroc_without": m4["au2"],
                 "change": ch, "percentile_among_random": pct})
    for ax in axs:
        ax.tick_params(labelsize=5.6)
    panel_letters([(axs[0], -0.62), (axs[1], -0.42), (axs[2], -0.42), (axs[3], -0.42)], y=1.06)
    save(fig, OUT, "RP1_models", data=pd.DataFrame(rows), caption=CAPTIONS["RP1"](sh, sets, m4, pct))


def _cap_rp1(sh, sets, m4, pct) -> str:
    rr = desc(f"M4 without {m4['n_up']} random non-R_up repeat inputs: AUROC change (exploratory, A-25)")
    bt = desc("M4 without its R_up inputs: AUROC change")
    co = lambda g: desc("M4 standardised coefficient: share > 0", g).value
    return f"""**Do models use the rising arm?**

(a) Training cancers (COAD, STAD, UCEC): the subtype MLPs (SE view, D-38) explain each tumour's own class, so the
MSI class's SHAP exists for the MSI test tumours. For each R_up or R_down input and cohort, the direction is the
Spearman correlation of the input's value with its SHAP over a fold's MSI test tumours, averaged over the folds;
bars = share of input × cohort units with a positive direction (higher inclusion pushes towards MSI); grey band =
2.5th–97.5th percentiles of 1,000 permutations of the inputs' rows among each fold's MSI test tumours (A-22);
dashed line = 0.5. (b, c) Held-out cancers (7,685 tumours, 50 MMRd): AUROC for MMRd against MMRp of the R_up (b)
and R_down (c) arm (the tumour's mean PSI deviation from its cancer's median over the class's exons, as a percentile
within the cancer); line = the class; grey = 1,000 sets of background exons matched to the class on baseline band
and coverage tertile within each cancer, which carry any tumour-wide shift as well; dotted line = AUROC 0.5. (d) M4, step 54's logistic regression on all repeat exons (fitted in the training cancers; held-out
AUROC {m4['au1']:.3f}), refitted at the same C without its {m4['n_up']} R_up inputs (line) and without
{m4['n_up']} inputs drawn at random from its other repeat inputs, 200 times (grey; exploratory, A-25).

- RPD1a, share of R_up directions above 0.5: {sg(T('RPD1a'), 3)} ({sh['obs']['R_up'][0]} of {sh['obs']['R_up'][1]}), p =
  {fp(T('RPD1a', 'p'))}, {vd('RPD1a')}; RPD1b, R_down: {sg(T('RPD1b'), 3)} ({sh['obs']['R_down'][0]} of
  {sh['obs']['R_down'][1]}), p = {fp(T('RPD1b', 'p'))}, {vd('RPD1b')}.
- RP1, R_up AUROC minus the matched median: {sg(T('RP1'), 3)} (AUROC {sets['R_up'][0]:.3f} against
  {np.median(sets['R_up'][1:]):.3f}), p = {fp(T('RP1', 'p'))}, {vd('RP1')}; RPC1, R_down: {sg(T('RPC1'), 3)}
  ({sets['R_down'][0]:.3f} against {np.median(sets['R_down'][1:]):.3f}), p = {fp(T('RPC1', 'p'))}, {vd('RPC1')}.
- M4 without its R_up inputs: AUROC change {sg(m4['au2'] - m4['au1'])} (paired bootstrap over test tumours with the
  fits held fixed, {num(bt.ci_low)} to {num(bt.ci_high)}); random removals: median {sg(rr.value)}, 2.5th–97.5th
  percentiles {num(rr.ci_low)} to {num(rr.ci_high)}; the R_up removal lies at the {pct:.1f}th percentile. No test;
  the fitted model does not need the rising arm. M4's standardised coefficients are positive for {100 * co('R_up'):.0f}%
  of its R_up inputs and {100 * co('R_down'):.0f}% of its R_down inputs."""


# ============================================================ RP2 junction reads
def RP2() -> None:
    _, H = state()
    d = delta0()
    order = ["R_down", "B sample", "R_up"]
    ci = {c: class_interval(c) for c in order}
    same(ci["R_up"]["ci_low"], T("RP2a", "ci_low"), "RP2a interval")
    same(ci["R_down"]["ci_low"], T("RPC2a", "ci_low"), "RPC2a interval")
    up = sel_of("R_up")
    x = d[up]
    sym = H.tab.symbol.to_numpy()[up]
    o = np.argsort(-x)
    x, sym = x[o], sym[o]
    bg = d[sel_of("B sample")]
    blo, bhi = np.percentile(bg, [2.5, 97.5])
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE - 0.3, 2.05), gridspec_kw={"width_ratios": [0.85, 2.3],
                                                                                  "wspace": 0.32})
    rows = []
    # a. class means with 95% intervals
    for k, c in enumerate(order):
        r = ci[c]
        col = ACCENT if c == "R_up" else (INK2 if c == "R_down" else MUTED)
        ax.plot([k, k], [r["ci_low"], r["ci_high"]], color=col, linewidth=1.1)
        ax.plot([k], [r["mean"]], marker="o", markersize=MS + 0.4, color=col, linestyle="")
        rows.append({"panel": "a", "class": c, "n_exons": r["n_exons"], "mean_delta": r["mean"], "ci_low": r["ci_low"],
                     "ci_high": r["ci_high"],
                     **(qcols("RP2a") if c == "R_up" else qcols("RPC2a") if c == "R_down" else {})})
    med = float(np.median(x))
    ax.plot([2], [med], marker="o", markersize=MS + 0.4, markerfacecolor="white", markeredgecolor=ACCENT,
            markeredgewidth=0.8, linestyle="")
    rows.append({"panel": "a", "class": "R_up", "n_exons": len(x), "median_delta": med})
    ax.axhline(P.SESOI, color=INK2, linewidth=0.6, linestyle=(0, (3, 2)))
    zero(ax, horizontal=True)
    rows.append({"panel": "a, b", "class": "reference lines", "mean_delta": 0.0, "sesoi": P.SESOI})
    ax.set_xticks(range(3))
    ax.set_xticklabels(["R_down", "back-\nground", "R_up"], fontsize=5.8)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(-0.006, 0.0125)
    ax.set_yticks([-0.005, 0, 0.005, 0.01])
    ax.set_ylabel("delta, PSI (MMRd − MMRp)", fontsize=6)
    ax.text(-0.55, P.SESOI, "SESOI", fontsize=5.4, color=INK2, va="bottom", ha="left")
    # b. the 92 R_up exons, sorted
    n = len(x)
    cols = [ACCENT_DARK if k < 10 else (ACCENT_LIGHT if v > 0 else LIGHT) for k, v in enumerate(x)]
    bx.bar(np.arange(n), x, width=0.8, color=cols, zorder=2)
    qb = np.quantile(bg, 1 - (np.arange(n) + 0.5) / n)          # background exons at the same quantiles
    bx.step(np.arange(n), qb, where="mid", color=INK2, linewidth=0.8, zorder=4)
    hand = [Line2D([], [], color=INK2, linewidth=0.8, label="background exons, same quantile")]
    for v, lab_, style in [(float(np.mean(x)), "R_up mean", "-"), (med, "R_up median", (0, (1, 1.5))),
                           (P.SESOI, "SESOI", (0, (3, 2)))]:
        c_ = ACCENT if lab_ != "SESOI" else INK2
        bx.axhline(v, color=c_, linewidth=0.7, linestyle=style, zorder=3)
        hand.append(Line2D([], [], color=c_, linewidth=0.7, linestyle=style, label=lab_))
    bx.legend(handles=hand, loc="upper right", fontsize=5.4, handlelength=2.2, handletextpad=0.4, borderaxespad=0.2)
    lead = dict(arrowstyle="-", linewidth=0.4, color=INK2, shrinkA=0, shrinkB=1)
    for k, yt in zip(range(3), [0.095, 0.083, 0.071]):
        bx.annotate(sym[k], (k, x[k]), xytext=(6, yt), textcoords="data", fontsize=5.4, va="center", ha="left",
                    color=INK, arrowprops=lead)
    zero(bx, horizontal=True)
    bx.set_xlim(-1, n + 0.5)
    bx.set_ylim(min(-0.03, float(x.min()) - 0.004), 0.1)
    bx.set_xticks([0, 9, 30, 60, 91])
    bx.set_xticklabels(["1", "10", "31", "61", "92"])
    bx.set_xlabel("R_up exons, ranked by delta", fontsize=6)
    bx.set_ylabel("delta, PSI", fontsize=6)
    cum = np.cumsum(x) / x.sum()
    for k, (v, s_) in enumerate(zip(x, sym)):
        rows.append({"panel": "b", "class": "R_up", "rank": k + 1, "symbol": s_, "delta": v, "top10": k < 10,
                     "cumulative_share_of_sum": cum[k]})
    rows += [{"panel": "b", "class": "B sample, same quantile", "rank": k + 1, "delta": v} for k, v in enumerate(qb)]
    rows.append({"panel": "b", "class": "R_up", "mean_delta": float(np.mean(x)), "median_delta": med,
                 "share_positive": float(np.mean(x > 0)), "top10_share_of_sum": float(cum[9]),
                 "mean_without_top10": float(np.mean(x[10:]))})
    for a in (ax, bx):
        a.tick_params(labelsize=5.6)
    panel_letters([(ax, -0.5), (bx, -0.12)], y=1.05)
    save(fig, OUT, "RP2_junction_reads", data=pd.DataFrame(rows), caption=CAPTIONS["RP2"](ci, x, blo, bhi, len(bg)))


def _cap_rp2(ci, x, blo, bhi, nbg) -> str:
    cum = np.cumsum(x) / x.sum()
    gl = sens("gene level (delta averaged within gene first)", "RP2a")
    md = sens("median instead of mean delta", "RP2a")
    return f"""**Does the rising arm rise in junction reads, and how is the rise spread over its exons?**

Held-out cancers, junction reads only (rMATS JC). Per exon, delta = the mean adjusted PSI of the MMRd tumours minus
the MMRp mean of the same cancer (≥ 5 MMRd tumours measured); a tumour's adjusted PSI is its deviation from the
cancer's median, minus its mean deviation over background exons of the same baseline band. (a) Mean delta of each
class with its 95% interval (studentized SE over MMR labels permuted within cancer, Welch–Satterthwaite df; the
intervals of R_up and R_down are those of RP2a and RPC2a, the background's is descriptive); hollow = R_up median;
dashed = the smallest effect of interest (SESOI, 0.003). Background = the B sample of exons without an acceptor T-run longer
than 5 nt ({ci['B sample']['n_exons']:,}). (b) The {len(x)} R_up exons ranked by delta; dark blue = the 10 largest; grey line =
the background exons' deltas at the same quantiles (rank k of {len(x)} against the background's (1 − (k − 0.5) /
{len(x)}) quantile; {nbg:,} exons; middle 95% {num(blo)} to {num(bhi)}); lines = R_up mean, median and the SESOI; the
three largest labelled.

- RP2a, R_up mean delta: {sg(T('RP2a'))} (95% interval {num(T('RP2a', 'ci_low'))} to {num(T('RP2a', 'ci_high'))}),
  p = {fp(T('RP2a', 'p'))}, {vd('RP2a')}; RPC2a, R_down: {sg(T('RPC2a'))} ({num(T('RPC2a', 'ci_low'))} to
  {num(T('RPC2a', 'ci_high'))}), p = {fp(T('RPC2a', 'p'))}, {vd('RPC2a')}.
- RP2c, share of R_up exons rising (studentized delta > c = 3.02) minus that of 5 matched background exons each:
  {sg(T('RP2c'), 3)}, p = {fp(T('RP2c', 'p'))}, {vd('RP2c')}; RP2b, per tumour (R_up arm, MMRd − MMRp):
  {sg(T('RP2b'))}, p = {fp(T('RP2b', 'p'))}, {vd('RP2b')}.
- The class mean is carried by a minority of exons: median {sg(float(np.median(x)))} (median sensitivity p =
  {fp(md.p)}); {100 * np.mean(x > 0):.0f}% of R_up exons above 0; the 10 largest hold {100 * cum[9]:.0f}% of the summed
  delta, and the mean of the other {len(x) - 10} is {sg(float(np.mean(x[10:])))}, below the SESOI. Gene level (delta
  averaged within gene first): {sg(gl.estimate)} ({num(gl.ci_low)} to {num(gl.ci_high)})."""


# ============================================================ RP3 tract and branch point
G_BINS = [(0, 5, "1–5"), (5, 10, "6–10"), (10, 10 ** 6, "≥ 11")]


def RP3() -> None:
    from scipy.stats import rankdata
    _, H = state()
    t, d = H.tab, delta0()
    g = t.g_exp.to_numpy(float)
    L = t.geo_L.fillna(t.acc_T).to_numpy(float)
    o = t.geo_o.fillna(t.run_near).to_numpy(float)
    sel = (t.acc_T >= P.LONG).to_numpy() & (g > 0) & H.elig & np.isfinite(L) & np.isfinite(o) & np.isfinite(d)
    x, gg, LL, oo = d[sel], g[sel], L[sel], o[sel]
    rho = P.partial_spearman(x, gg, np.column_stack([rankdata(LL), rankdata(oo)]))
    same(rho, T("RP3e"), "RP3e")
    Hm = P.residualize(np.column_stack([rankdata(LL), rankdata(oo)]), len(x))
    rx, rg = Hm @ rankdata(x), Hm @ rankdata(gg)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE * 0.72, 2.0), gridspec_kw={"wspace": 0.42})
    rows = []
    rng = np.random.default_rng(0)
    for k, (lo, hi, lab_) in enumerate(G_BINS):
        m = (gg > lo) & (gg <= hi)
        jit = rng.uniform(-0.17, 0.17, m.sum())
        ax.scatter(k + jit, x[m], s=8, color=ACCENT_LIGHT, edgecolors="white", linewidths=0.3, zorder=2)
        mu, me = float(x[m].mean()), float(np.median(x[m]))
        ax.plot([k - 0.28, k + 0.28], [mu, mu], color=ACCENT_DARK, linewidth=1.3, zorder=3)
        ax.plot([k - 0.2, k + 0.2], [me, me], color=ACCENT_DARK, linewidth=0.7, linestyle=(0, (1, 1)), zorder=3)
        rows.append({"panel": "a", "gap_bin_nt": lab_, "n_exons": int(m.sum()), "mean_delta": mu, "median_delta": me})
    zero(ax, horizontal=True)
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"{b[2]}\n({r['n_exons']} exons)" for b, r in zip(G_BINS, rows)], fontsize=5.8)
    ax.tick_params(axis="x", length=0)
    ax.set_xlabel("gap from the run to the mapped branch point, nt", fontsize=6)
    ax.set_ylabel("delta, PSI", fontsize=6)
    ax.set_ylim(-0.03, 0.1)
    ax.set_xlim(-0.6, 2.6)
    bx.scatter(rg, rx, s=8, color=ACCENT, edgecolors="white", linewidths=0.3, zorder=2)
    b_ = float(rx @ rg / (rg @ rg))
    xs = np.array([rg.min(), rg.max()])
    bx.plot(xs, b_ * xs, color=INK2, linewidth=0.7, zorder=1)
    zero(bx, horizontal=True)
    zero(bx)
    bx.set_xlabel("gap, rank residual", fontsize=6)
    bx.set_ylabel("delta, rank residual", fontsize=6)
    lo_ = float(rx.min())
    bx.set_ylim(lo_ - 0.42 * (float(rx.max()) - lo_), float(rx.max()) * 1.12)
    bx.text(0.98, 0.02, f"partial ρ = {rho:.2f}, p = {fp(T('RP3e', 'p'))}", transform=bx.transAxes, fontsize=5.8,
            va="bottom", ha="right", color=INK)
    rows.append({"panel": "b", "partial_spearman": rho, "ols_slope_line": b_, **qcols("RP3e")})
    for a in (ax, bx):
        a.tick_params(labelsize=5.6)
    panel_letters([(ax, -0.3), (bx, -0.3)], y=1.05)
    local = pd.DataFrame({"feature_id": t.feature_id.to_numpy()[sel], "cls": t.cls.to_numpy()[sel], "delta": x,
                          "g_exp": gg, "L": LL, "o": oo, "rank_residual_delta": rx, "rank_residual_gap": rg})
    rows += [{"panel": "a, b", "class": "per exon", **r} for r in local.drop(columns=["feature_id"]).to_dict("records")]
    save(fig, OUT, "RP3_branch_point", data=pd.DataFrame(rows), caption=CAPTIONS["RP3"](local))


def _cap_rp3(local) -> str:
    from scipy.stats import rankdata
    x, g, L, o = (local[c].to_numpy(float) for c in ["delta", "g_exp", "L", "o"])
    ps = lambda m: P.partial_spearman(x[m], g[m], np.column_stack([rankdata(L[m]), rankdata(o[m])]))
    loo = [ps(np.arange(len(x)) != i) for i in range(len(x))]
    d_exp = g + L + o - 1
    top_all = np.argsort(-x)[:6]
    up = np.flatnonzero(local.cls.eq("R_up").to_numpy())
    top_up = up[np.argsort(-x[up])][:6]
    bins = "; ".join(f"{lab_} nt: mean {sg(float(x[(g > lo) & (g <= hi)].mean()))}, median "
                     f"{sg(float(np.median(x[(g > lo) & (g <= hi)])))} ({int(((g > lo) & (g <= hi)).sum())} exons)"
                     for lo, hi, lab_ in G_BINS)
    d60 = sens("R_up cut at d_bp > 60", "RP3b")
    d50 = sens("R_up cut at d_bp >= 50", "RP3b")
    return f"""**Does the rise need the T-run, and is it associated with the gap to the mapped branch point?**

Held-out cancers, junction reads. The {len(local)} eligible exons with an acceptor T-run of ≥ 20 nt and an
experimentally mapped branch point upstream of the run (BPHunter, lifted to GRCh38; {len(up)} of them in the rising
arm). For these exons the branch point's distance from the AG is exactly the offset of the run plus its length plus
the gap, minus 1, so a cut on that distance mixes run length and position with the gap; the test therefore asks
whether delta is rank-associated with the gap at fixed run length and offset. (a) delta by gap bin (descriptive,
the bins chosen after the results): dots = exons, solid line = mean, dotted = median. (b) The same exons: rank of
delta against rank of the gap, each with run length and offset (ranks) regressed out; line = least squares.

- RP3e, partial Spearman of delta with the gap given run length and offset: {sg(T('RP3e'), 3)}, p =
  {fp(T('RP3e', 'p'))}, {vd('RP3e')} (the branch-point test of the verdict); leaving out one exon at a time,
  {min(loo):.2f} to {max(loo):.2f}. RP3c, the same with the heuristic branch point (139 exons): {sg(T('RP3c'), 3)},
  p = {fp(T('RP3c', 'p'))}, {vd('RP3c')}.
- It does not replicate in the training cancers: RPD3b, the same test there, {sg(T('RPD3b'), 3)}, p =
  {fp(T('RPD3b', 'p'))}, {vd('RPD3b')} ({TESTS.at['RPD3b', 'n']}, 150 MMRd tumours).
- By gap bin (descriptive): {bins}. The rise sits at gaps above 5 nt, highest at 6–10 nt, and is not monotone.
- RP3a, needs the T-run (R_up minus exons with a distant branch point and no run): {sg(T('RP3a'))}, p =
  {fp(T('RP3a', 'p'))}, {vd('RP3a')}.
- The 36-nt cut on distance: RP3b (R_up minus long runs with a near heuristic branch point) {sg(T('RP3b'))}, p =
  {fp(T('RP3b', 'p'))}, {vd('RP3b')}; RP3d (mapped branch point ≥ 36 nt minus < 36 nt) {sg(T('RP3d'))}, p =
  {fp(T('RP3d', 'p'))}, {vd('RP3d')}. Pre-specified sensitivity cuts for R_up (no q): ≥ 50 nt {sg(d50.estimate)},
  p = {fp(d50.p)}; > 60 nt {sg(d60.estimate)}, p = {fp(d60.p)}. Of the 6 largest rises,
  {int((d_exp[top_all] < P.FAR).sum())} have a mapped branch point < 36 nt from the AG
  ({int((d_exp[top_up] < P.FAR).sum())} of the 6 largest among the rising-arm exons)."""


# ============================================================ RP4 controls
CONTROLS = [("B_far", "RP3a"), ("ppt_far", "RP4b"), ("alu", "RP4e"), ("polyA15", "RP4a"), ("deep", "RP4c"),
            ("donor", "RP4d")]


def contrast_power() -> dict:
    """Per control contrast: the SD of its exon-label (class) permutation null, rebuilt with step 62's
    coarsening, class_perm and seed (the observed value must match tests.csv), the difference detectable
    with 80% power at two-sided 0.05 (normal approximation, 2.80 SD), and the approximate power if the
    control did not move at all (true difference = RP2a's estimate). Descriptive, for reading 'not resolved'."""
    from scipy.stats import norm

    def go():
        S, H = state()
        t = H.tab
        up = t.cls.eq("R_up").to_numpy()
        out = {}
        for c, i in CONTROLS:
            ctrl = {"B_far": t.b_sample.to_numpy(bool) & t.cls.eq("B_far").to_numpy(),
                    "ppt_far": t.ctl_ppt_far.to_numpy(bool), "alu": t.ctl_alu.to_numpy(bool),
                    "polyA15": t.ctl_polyA15.to_numpy(bool), "deep": t.ctl_deep.to_numpy(bool),
                    "donor": t.ctl_don.to_numpy(bool)}[c]
            f, cc = up & H.elig, ctrl & H.elig
            lab, kept = P.coarsen(H.band, H.terc, f, cc, per_first=0)
            use = kept & (f | cc)
            fin = np.isfinite(H.D[0, use])
            obs, null = P.class_perm(H.D[0, use][fin].astype(float), f[use][fin], lab[use][fin], P.N_CLASS_PERM, P.SEED)
            out[i] = {"obs": float(obs), "null_sd": float(np.std(null))}
        return out
    r = cached("contrast_power", go)
    for i, v in r.items():
        same(v["obs"], T(i), f"{i} class_perm")
        z = T("RP2a") / v["null_sd"]
        v.update({"mdd_80": (norm.ppf(0.975) + norm.ppf(0.8)) * v["null_sd"],
                  "power_if_control_flat": float(norm.cdf(z - norm.ppf(0.975)) + norm.cdf(-z - norm.ppf(0.975)))})
    return r


def RP4() -> None:
    ci = {c: class_interval(c) for c in ["R_up"] + [c for c, _ in CONTROLS]}
    pw = contrast_power()
    names = ["R_up"] + [c for c, _ in CONTROLS]
    fig = plt.figure(figsize=(DOUBLE * 0.78, 2.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0], wspace=0.04)
    ax, tx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    rows = []
    for y, c in enumerate(names):
        r = ci[c]
        col = ACCENT if c == "R_up" else INK2
        ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=col, linewidth=1.0)
        ax.plot([r["mean"]], [y], marker="o", markersize=MS + 0.3, color=col, linestyle="")
        row = {"class": c, "label": CLASS_NAME[c], "n_exons": r["n_exons"], "mean_delta": r["mean"],
               "ci_low": r["ci_low"], "ci_high": r["ci_high"]}
        if c != "R_up":
            i = dict(CONTROLS)[c]
            tx.text(0.02, y, f"{sg(T(i))}", fontsize=5.8, va="center", ha="left", color=INK)
            tx.text(0.42, y, fp(T(i, "p")), fontsize=5.8, va="center", ha="left", color=INK)
            tx.text(0.64, y, vd(i), fontsize=5.8, va="center", ha="left",
                    color=ACCENT_DARK if vd(i) == "supported" else MUTED)
            n_ctl = TESTS.at[i, "n"]
            pw_ = pw[i]
            row.update({**qcols(i), "contrast_n": n_ctl, "class_perm_null_sd": pw_["null_sd"],
                        "detectable_difference_80pct": pw_["mdd_80"], "power_if_control_flat": pw_["power_if_control_flat"]})
        rows.append(row)
    ax.axvline(P.SESOI, color=INK2, linewidth=0.6, linestyle=(0, (3, 2)))
    zero(ax)
    rows.append({"class": "reference lines", "mean_delta": 0.0, "sesoi": P.SESOI})
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([f"{CLASS_NAME[c]} ({ci[c]['n_exons']:,})" for c in names], fontsize=5.9)
    ax.set_ylim(len(names) - 0.5, -1.2)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("mean delta, PSI (95% interval)", fontsize=6)
    ax.set_xticks([-0.004, 0, 0.004, 0.008])
    ax.tick_params(axis="x", labelsize=5.6)
    tx.set_ylim(len(names) - 0.5, -1.2)
    tx.set_xlim(0, 1)
    tx.axis("off")
    for x_, s_ in [(0.02, "R_up − control"), (0.42, "p"), (0.64, "verdict")]:
        tx.text(x_, -0.8, s_, fontsize=5.8, va="center", ha="left", color=INK2,
                fontstyle="italic" if s_ == "p" else "normal")
    save(fig, OUT, "RP4_controls", data=pd.DataFrame(rows), caption=CAPTIONS["RP4"](ci, pw))


def ctl_n(i: str) -> str:
    """'85 vs 61 exons (effective controls 27)' -> '61 controls, 27 effective'."""
    import re
    m = re.match(r"\d+ vs (\d+) exons \(effective controls (\d+)\)", TESTS.at[i, "n"])
    return f"{m.group(1)} controls, {m.group(2)} effective"


def _cap_rp4(ci, pw) -> str:
    own = lambda g: desc("control class: own mean delta against 0 (label_perm p)", g)
    na = sens("acceptor not in an Alu", "RP3a")
    return f"""**Is it the microsatellite, or would any long run or Alu acceptor do?**

Held-out cancers, junction reads. Mean delta (dot) of the rising arm and of each control class, with a descriptive
95% interval (RP2a's construction); dashed = SESOI. Controls have no acceptor T-run longer than 7 nt and differ from
R_up in one feature: a distant branch point and no run (B sample), a long pyrimidine tract that is not a
microsatellite with a distant branch point, an acceptor inside an Alu, a poly-A run (≥ 15 nt) at the acceptor, a
T-run of ≥ 20 nt deep in the intron, or after the donor. The tests compare R_up with each control within strata of
baseline band × coverage tertile, coarsened until each R_up stratum holds ≥ 3 controls (stratum-weighted difference
of mean delta; right-hand columns), so they differ from the raw difference of the dots. p = the larger of the
exon-label (class) permutation p and the MMR-label permutation p.

- {tline('RP3a', 'R_up minus distant branch point without a run')}
- {tline('RP4b', 'R_up minus long pyrimidine tract')} {tline('RP4e', 'R_up minus Alu acceptor')}
- {tline('RP4c', 'R_up minus deep-intron T-run')} {tline('RP4a', 'R_up minus poly-A')}
  {tline('RP4d', 'R_up minus donor-side T-run')}
- RP4a and RP4d are underpowered rather than negative. Their exon-permutation p ({fp(T('RP4a', 'p_class'))},
  {fp(T('RP4d', 'p_class'))}) is far above their label-permutation p ({fp(T('RP4a', 'p_label'))},
  {fp(T('RP4d', 'p_label'))}). The SD of the exon-permutation null is {pw['RP4a']['null_sd']:.4f} for poly-A
  ({ctl_n('RP4a')}) and {pw['RP4d']['null_sd']:.4f} for donor-side runs, so the difference
  detectable with 80% power is about {pw['RP4a']['mdd_80']:.3f} and {pw['RP4d']['mdd_80']:.3f} (normal
  approximation), larger than RP2a itself; if the controls did not move at all, the power would be about
  {pw['RP4a']['power_if_control_flat']:.2f} and {pw['RP4d']['power_if_control_flat']:.2f}. Own means: poly-A
  {sg(own('polyA15').value)} (p = {fp(own('polyA15').p)}), donor-side {sg(own('donor').value)} (p =
  {fp(own('donor').p)}).
- 61 of the 92 R_up acceptors lie in an Alu. For the other 31, R_up minus distant branch point: {sg(na.estimate)},
  p = {fp(na.p)}."""


# ============================================================ RP5 confounders
def RP5() -> None:
    bs = desc("control class: own mean delta against 0 (label_perm p)", "B sample (all)").value
    ex_rows = [("no covariates (difference of means)", T("RP2a") - bs, np.nan, "R_up mean − B-sample mean"),
               ("RP5a: eight exon covariates", T("RP5a"), T("RP5a", "p"), "RP5a"),
               ("+ full splice-site score, acceptor Alu",
                sens("RP5a's full model (+ -20..+3 PWM, acceptor Alu)", "RP5a").estimate,
                sens("RP5a's full model (+ -20..+3 PWM, acceptor Alu)", "RP5a").p, "sensitivity"),
               ("+ longest T-run in the coding sequence", sens("RP5a with cds_max_run added", "RP5a").estimate,
                sens("RP5a with cds_max_run added", "RP5a").p, "sensitivity")]
    tu_rows = [("no covariates (RP2b arm contrast)", T("RP2b"), T("RP2b", "p"), "RP2b"),
               ("RP5b: five tumour covariates", T("RP5b"), T("RP5b", "p"), "RP5b"),
               ("without purity (all tumours)", sens("RP5b without purity (all held-out tumours)", "RP5b").estimate,
                sens("RP5b without purity (all held-out tumours)", "RP5b").p, "sensitivity"),
               ("+ MKI67 expression",
                sens("RP5b with MKI67 log2(CPM + 1) percentile within cancer", "RP5b").estimate,
                sens("RP5b with MKI67 log2(CPM + 1) percentile within cancer", "RP5b").p, "sensitivity")]
    fig = plt.figure(figsize=(DOUBLE * 0.62, 2.05))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.22], wspace=0.04)
    ax, tx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    labels, rows, y = [], [], 0
    for head, items in [("exon level: R_up coefficient", ex_rows), ("tumour level: MMRd coefficient", tu_rows)]:
        labels.append((y, head, True))
        y += 1
        for lab_, v, p_, src in items:
            col = ACCENT if src.startswith("RP5") else (MUTED if src == "sensitivity" else INK2)
            ax.plot([v], [y], marker="o", markersize=MS + 0.3, color=col, linestyle="")
            tx.text(0.05, y, fp(p_) if np.isfinite(p_) else "–", fontsize=5.8, va="center", ha="left", color=INK)
            labels.append((y, lab_, False))
            rows.append({"level": head.split(":")[0], "row": lab_, "estimate": v, "p": p_, "source": src,
                         **({k: v_ for k, v_ in qcols(src).items() if k in ("q", "q_family", "verdict")}
                            if src in TESTS.index else {})})
            y += 1
    zero(ax)
    rows.append({"row": "reference line", "estimate": 0.0})
    ax.set_yticks([a for a, _, _ in labels])
    ax.set_yticklabels([b for _, b, _ in labels], fontsize=5.8)
    for tl, (_, _, bold) in zip(ax.get_yticklabels(), labels):
        tl.set_fontweight("bold" if bold else "normal")
        tl.set_color(INK if bold else INK2)
    ax.set_ylim(y - 0.5, -1.1)
    ax.set_xlim(-0.0005, 0.0105)
    ax.set_xticks([0, 0.003, 0.006, 0.009])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("coefficient, PSI", fontsize=6)
    ax.tick_params(axis="x", labelsize=5.6)
    tx.set_ylim(y - 0.5, -1.1)
    tx.axis("off")
    tx.text(0.05, -0.6, "p", fontsize=5.8, va="center", ha="left", color=INK2, fontstyle="italic")
    save(fig, OUT, "RP5_confounders", data=pd.DataFrame(rows), caption=CAPTIONS["RP5"](bs))


def _cap_rp5(bs) -> str:
    return f"""**Do exon and tumour covariates explain the rise?**

Held-out cancers, junction reads. Exon level: OLS of delta on an R_up indicator against the B sample of background
exons (92 R_up + 4,756 exons), with exon covariates; the core splice-site score (−3..+3) is blind to the tract,
while the full −20..+3 score and the acceptor-Alu flag measure the tract itself and enter only a sensitivity model.
Tumour level: OLS of each held-out tumour's R_up arm on MMR status with cancer indicators and the tumour's global
inclusion axis, stroma, immune infiltrate, retained-intron level and purity. Blue = the pre-registered models;
dark grey = no covariates (the difference of the two class means at exon level; RP2b at tumour level); light grey =
sensitivity models. No intervals were computed for these coefficients; p is shown instead.

- {tline('RP5a', 'R_up coefficient, Freedman–Lane within band × tertile and MMR-label permutation, p = max')}
  {tline('RPC5a', 'R_down')}
- {tline('RP5b', 'MMRd coefficient, Freedman–Lane with tumour residuals permuted within cancer')}
  {tline('RPC5b', 'R_down arm')}
- No covariates: R_up mean minus B-sample mean {sg(T('RP2a') - bs)}; {tline('RP2b', 'arm contrast')}"""


# ============================================================ RP6 specificity
def _z_within() -> np.ndarray:
    S, _ = state()
    lab, arm = S.lab, S.tv.arm_up.to_numpy(float)
    z = np.full(len(lab), np.nan)
    for p in S.held:
        i = np.flatnonzero(lab.project.eq(p).to_numpy())
        ref = i[(~lab.mmr.to_numpy()[i]) & (~lab.hyper.to_numpy()[i]) & np.isfinite(arm[i])]
        if len(ref) > 2:
            z[i] = (arm[i] - arm[ref].mean()) / arm[ref].std(ddof=1)
    return z


def RP6() -> None:
    from scipy.stats import mannwhitneyu
    S, _ = state()
    lab, z = S.lab, _z_within()
    h = lab.held_out.to_numpy() & np.isfinite(z)
    groups = [("MMRd", h & lab.mmr.to_numpy(), ACCENT),
              ("hypermutated\nMMRp", h & lab.hyper.to_numpy() & ~lab.mmr.to_numpy(), INK2),
              ("other MMRp", h & ~lab.hyper.to_numpy() & ~lab.mmr.to_numpy(), MUTED)]
    u = mannwhitneyu(z[groups[0][1]], z[groups[1][1]], alternative="two-sided")
    same(u.statistic / (groups[0][1].sum() * groups[1][1].sum()) - 0.5, T("RP6a"), "RP6a")
    nul = pd.read_csv(CODE / "rp6c_null_contrasts.csv")
    nul["abs_s"] = nul.s.abs()
    nul = nul.sort_values("abs_s").reset_index(drop=True)
    s_msi = T("RP6c")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE - 0.3, 2.1), gridspec_kw={"width_ratios": [1.0, 1.35],
                                                                                 "wspace": 0.3})
    rows, local = [], []
    rng = np.random.default_rng(0)
    for k, (name, m, col) in enumerate(groups):
        v = z[m]
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        p5, p95 = np.percentile(v, [5, 95])
        ax.add_patch(plt.Rectangle((k - 0.22, q1), 0.44, q3 - q1, facecolor="white", edgecolor=col, linewidth=0.7,
                                   zorder=2))
        ax.plot([k, k], [p5, q1], color=col, linewidth=0.6)
        ax.plot([k, k], [q3, p95], color=col, linewidth=0.6)
        ax.plot([k - 0.22, k + 0.22], [med, med], color=col, linewidth=1.2, zorder=3)
        if name.startswith("MMRd") or name.startswith("hyper"):
            ax.scatter(k + rng.uniform(-0.33, 0.33, len(v)), v, s=3 if len(v) > 100 else 6,
                       color=col, alpha=0.45 if len(v) > 100 else 0.8, linewidths=0, zorder=1)
        rows.append({"panel": "a", "group": name.replace("\n", " "), "n_tumours": int(m.sum()), "median_z": med,
                     "q25": q1, "q75": q3, "p05": p5, "p95": p95})
        local += [{"panel": "a", "group": name.replace("\n", " "), "file_id": f, "project": pj, "z": zz}
                  for f, pj, zz in zip(lab.file_id[m], lab.project[m], v)]
    rows.append({"panel": "a", "group": "MMRd vs hypermutated MMRp", "auroc": T("RP6a") + 0.5, **qcols("RP6a")})
    zero(ax, horizontal=True)
    ax.set_xticks(range(3))
    ax.set_xticklabels([g[0] for g in groups], fontsize=5.8)
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(-4, 6)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylabel("R_up arm, z within cancer", fontsize=6)
    # b. the 57 non-MSI subtype contrasts
    over = nul.abs_s >= s_msi
    bx.scatter(nul.index, nul.abs_s, s=7, color=np.where(over, INK2, LIGHT), linewidths=0, zorder=2)
    bx.axhline(s_msi, color=ACCENT, linewidth=1.0, zorder=3)
    bx.text(1, s_msi - 0.0002, f"MSI (MMRd − MMRp), {s_msi:.4f}", color=ACCENT, fontsize=5.6, va="top", ha="left")
    lead = dict(arrowstyle="-", linewidth=0.4, color=INK2, shrinkA=0, shrinkB=2)
    for yt, (ii, r) in zip([0.0155, 0.0143, 0.0127, 0.0114], nul[over].sort_values("abs_s", ascending=False).iterrows()):
        name = f"{r.cohort.split('_')[0]} {r.a[0].lower() + r.a[1:] if not r.a.isupper() else r.a} − {r.b}"
        bx.annotate(name, (ii, r.abs_s), xytext=(44, yt), textcoords="data", fontsize=5.3, ha="right", va="center",
                    color=INK, arrowprops=lead)
    for ii, r in nul.iterrows():
        rows.append({"panel": "b", "group": f"{r.cohort} {r.a} vs {r.b}", "n_a": int(r.n_a), "n_b": int(r.n_b),
                     "s": r.s, "abs_s": r.abs_s, "exceeds_msi": bool(r.abs_s >= s_msi)})
    rows.append({"panel": "b", "group": "MSI held out", "abs_s": s_msi, **qcols("RP6c")})
    bx.set_xlim(-2, len(nul) + 1)
    bx.set_ylim(0, 0.0165)
    bx.set_xlabel("non-MSI subtype contrasts, ranked", fontsize=6)
    bx.set_ylabel("|s|, PSI", fontsize=6)
    for a in (ax, bx):
        a.tick_params(labelsize=5.6)
    panel_letters([(ax, -0.25), (bx, -0.17)], y=1.05)
    save(fig, OUT, "RP6_specificity", data=pd.DataFrame(rows), caption=CAPTIONS["RP6"](groups, nul, s_msi),
         local_data=(LOCAL, pd.DataFrame(local)))


def _cap_rp6(groups, nul, s_msi) -> str:
    w = sens("RP6a within cancer", "RP6a")
    ms = sens("RP6b with MSIsensor MSI-H", "RP6b")
    n_over = int((nul.abs_s >= s_msi).sum())
    return f"""**Is it mismatch repair, or mutation load, or any subtype difference?**

Held-out cancers, junction reads. (a) Each tumour's R_up arm (mean adjusted PSI over its R_up exons) as a z within
its cancer against the cancer's non-hypermutated MMRp tumours: MMRd ({int(groups[0][1].sum())}), hypermutated MMRp
({int(groups[1][1].sum())}; hypermutated without MMR deficiency) and the other MMRp tumours
({int(groups[2][1].sum()):,}); boxes = quartiles, whiskers = 5th–95th percentiles, dots = tumours (other MMRp not
drawn). (b) The MSI statistic s (D-70's RV1 in junction reads: 20 held-out MMRd tumours, each with a random MMRp
tumour of its cancer; mean over R_up exons of the difference in mean adjusted PSI; mean over 50 draws) against the
same statistic for each of the 57 non-MSI subtype contrasts of the subtype task (two classes of one cohort with ≥ 20
tumours each, MSI, POLE and MMRd tumours removed; 20 tumours drawn from each class, 50 draws); contrasts at or above
MSI's |s| are labelled. The null contrasts are not matched by cancer, while the MSI draws pair tumours within
cancer.

- {tline('RP6a', 'AUROC − 0.5, MMRd against hypermutated MMRp, Mann–Whitney')} Within the cancers that have both:
  {sg(w.estimate, 3)}, p = {fp(w.p)}.
- {tline('RP6b', 'arm contrast with MANTIS MSI-H as the label')} With MSIsensor: {sg(ms.estimate)}, p = {fp(ms.p)}.
- {tline('RP6c', 'MSI s ranked against the 57 null contrasts (p floor 1/58)')} {n_over} null contrasts reach MSI's
  |s|: three involve leiomyosarcoma in SARC (one lineage) and one is PCPG's cortical-admixture class."""


# ============================================================ RP7 within patient
def RP7() -> None:
    pr = pd.read_parquet(P.PROOF / "tests" / "pairs.parquet")
    pr = pr[pr.project.isin(set(pr.project[pr.mmr])) & np.isfinite(pr.pair_up.to_numpy(float))].copy()
    cancers = ["COAD", "STAD", "UCEC", "BRCA", "PRAD"]
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.55, 2.0))
    rows, local = [], []
    rng = np.random.default_rng(0)
    for k, c in enumerate(cancers):
        g = pr[pr.project.eq(c)]
        for st, col, off, s in [(False, LIGHT, -0.14, 5), (True, ACCENT, 0.14, 11)]:
            v = g[g.mmr.eq(st)].pair_up.to_numpy(float)
            ax.scatter(k + off + rng.uniform(-0.08, 0.08, len(v)), v, s=s, color=col,
                       edgecolors="white" if st else "none", linewidths=0.3, zorder=3 if st else 2)
            mu = float(v.mean()) if len(v) >= 5 else np.nan
            if np.isfinite(mu):
                ax.plot([k + off - 0.13, k + off + 0.13], [mu, mu], color=INK2 if not st else ACCENT_DARK,
                        linewidth=1.1, zorder=4)
            rows.append({"cancer": c, "training_cancer": c in P.TRAIN, "mmr": "MMRd" if st else "MMRp",
                         "n_pairs": len(v), "mean_pair_arm": mu,
                         "note": "" if np.isfinite(mu) else "mean blanked (< 5 pairs)"})
            local += [{"cancer": c, "mmr": "MMRd" if st else "MMRp", "case_id": ci, "pair_arm": vv}
                      for ci, vv in zip(g[g.mmr.eq(st)].case_id, v)]
    rows.append({"cancer": "all", "mmr": "MMRd vs MMRp", **qcols("RP7")})
    rows.append({"cancer": "reference line", "mean_pair_arm": 0.0})
    zero(ax, horizontal=True)
    ax.set_xticks(range(len(cancers)))
    ax.set_xticklabels([f"{c}\n({'training' if c in P.TRAIN else 'held out'})" for c in cancers], fontsize=5.6)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("R_up arm, tumour − own normal (PSI)", fontsize=6)
    ax.tick_params(labelsize=5.6)
    hand = [Line2D([], [], marker="o", color=ACCENT, linestyle="", markersize=3.4, label="MMRd pair"),
            Line2D([], [], marker="o", color=LIGHT, linestyle="", markersize=2.6, label="MMRp pair"),
            Line2D([], [], color=INK2, linewidth=1.1, label="mean (≥ 5 pairs)")]
    ax.legend(handles=hand, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=5.4, handletextpad=0.3,
              columnspacing=1.0, borderaxespad=0.1)
    save(fig, OUT, "RP7_within_patient", data=pd.DataFrame(rows), caption=CAPTIONS["RP7"](pr),
         local_data=(LOCAL, pd.DataFrame(local)))


def _cap_rp7(pr) -> str:
    md = pr[pr.mmr]
    n_tr = int(md.project.isin(P.TRAIN).sum())
    by = md.project.value_counts()
    return f"""**Is the rise there within the same patient?**

The tumour–normal pairs of the cancers with at least one MMRd pair. Pair arm = the tumour's mean PSI deviation over
its R_up exons minus that of its own normal tissue (the same exons, the same adjustment for each sample's
inclusion level). Dots = pairs; lines = means where a group has ≥ 5 pairs.
{int(pr.mmr.sum())} MMRd and {int((~pr.mmr).sum())} MMRp pairs; {n_tr} of the MMRd pairs are in training cancers
(COAD {by.get('COAD', 0)}, STAD {by.get('STAD', 0)}, UCEC {by.get('UCEC', 0)}), so this test is not held out
(family RPP).

- {tline('RP7', 'mean over MMRd pairs of the pair arm minus the MMRp mean of the same cancer; MMR labels permuted among each cancer’s pairs (10,000)')}
  {tline('RPC7', 'the falling arm')}
- The two held-out MMRd pairs (one each in BRCA and PRAD) go one in each direction relative to their cancer's
  MMRp pairs."""


# ============================================================ RP8 dose
def RP8() -> None:
    S, H = state()
    lab = S.lab
    arm = S.tv.arm_up.to_numpy(float)
    man = lab.mantis.to_numpy(float)
    m = lab.mmr.to_numpy() & lab.held_out.to_numpy() & np.isfinite(arm) & np.isfinite(man)
    up = sel_of("R_up")
    d, L = delta0()[up], H.tab.acc_T.to_numpy(float)[up]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE * 0.72, 2.0), gridspec_kw={"wspace": 0.42})
    ax.scatter(man[m], arm[m], s=9, color=ACCENT, edgecolors="white", linewidths=0.3)
    zero(ax, horizontal=True)
    ax.set_xlabel("MANTIS score", fontsize=6)
    ax.set_ylabel("R_up arm, PSI", fontsize=6)
    ax.text(0.98, 0.03, f"partial ρ = {T('RP8a'):.2f}, p = {fp(T('RP8a', 'p'))}", transform=ax.transAxes, fontsize=5.8,
            va="bottom", ha="right")
    bx.scatter(L + np.random.default_rng(0).uniform(-0.25, 0.25, len(L)), d, s=8, color=ACCENT, edgecolors="white",
               linewidths=0.3)
    zero(bx, horizontal=True)
    bx.set_xlabel("acceptor T-run length, nt", fontsize=6)
    bx.set_ylabel("delta, PSI", fontsize=6)
    bx.text(0.98, 0.98, f"ρ = {T('RP8c'):.2f}, p = {fp(T('RP8c', 'p'))}", transform=bx.transAxes, fontsize=5.8,
            va="top", ha="right")
    for a in (ax, bx):
        a.tick_params(labelsize=5.6)
    panel_letters([(ax, -0.3), (bx, -0.3)], y=1.05)
    rows = [{"panel": "a", "n_tumours": int(m.sum()), "median_mantis": float(np.median(man[m])),
             "median_arm": float(np.median(arm[m])), **qcols("RP8a")},
            {"panel": "a, b", "reference line": 0.0},
            {"panel": "b", **qcols("RP8c")}]
    rows += [{"panel": "b", "run_length": a_, "delta": b_} for a_, b_ in zip(L, d)]
    local = pd.DataFrame({"file_id": lab.file_id[m], "project": lab.project[m], "mantis": man[m], "arm_up": arm[m]})
    save(fig, OUT, "RP8_dose", data=pd.DataFrame(rows), caption=CAPTIONS["RP8"](int(m.sum()), L),
         local_data=(LOCAL, local))


def _cap_rp8(n, L) -> str:
    return f"""**Does the rise follow dose?**

(a) The held-out MMRd tumours ({n}): R_up arm against the MANTIS score (the exome instability score). (b) The R_up
exons: delta against acceptor T-run length ({int(L.min())}–{int(L.max())} nt; jittered). Both tests ask for a graded
relation, which D-70's percentile read-out had not found (RV6).

- {tline('RP8a', 'partial Spearman of the arm with MANTIS given cancer, MANTIS permuted within cancer', 3)}
  {tline('RPC8', 'the falling arm (predicted < 0)', 3)}
- {tline('RP8b', 'the same with log(1 + deletion count), all 200 MMRd tumours (not held out)', 3)}
- {tline('RP8c', 'Spearman of delta with run length over R_up, run length permuted within baseline band', 3)}"""


# ============================================================ RPM mechanism
def RPM() -> None:
    names = ["R_up", "shield", "long, antisense Alu", "long, other"]
    ci = {c: class_interval(c) for c in names}
    fig = plt.figure(figsize=(DOUBLE * 0.92, 1.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.25], wspace=0.04)
    ax, tx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    rows = []
    ys = [0, 1.1, 2.7, 3.7]
    for y, c in zip(ys, names):
        r = ci[c]
        col = ACCENT if c == "R_up" else INK2
        ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=col, linewidth=1.0)
        ax.plot([r["mean"]], [y], marker="o", markersize=MS + 0.3, color=col, linestyle="")
        rows.append({"class": c, "label": CLASS_NAME[c], "n_exons": r["n_exons"], "mean_delta": r["mean"],
                     "ci_low": r["ci_low"], "ci_high": r["ci_high"]})
    for x_, s_ in [(0.0, "test"), (0.48, "estimate"), (0.67, "p"), (0.79, "verdict")]:
        tx.text(x_, -0.45, s_, fontsize=5.8, va="center", ha="left", color=INK2, fontstyle="italic" if s_ == "p" else "normal")
    for y, i, what in [(0.55, "RPM1", "R_up − shield"), (1.1, "RPM2", "shield against 0"),
                       (3.2, "RPM3", "antisense Alu − other")]:
        tx.text(0.0, y, f"{i}  {what}", fontsize=5.8, va="center", ha="left", color=INK2)
        tx.text(0.48, y, sg(T(i)), fontsize=5.8, va="center", ha="left", color=INK)
        tx.text(0.67, y, fp(T(i, "p")), fontsize=5.8, va="center", ha="left", color=INK)
        tx.text(0.79, y, vd(i), fontsize=5.8, va="center", ha="left", color=MUTED)
        rows.append({"class": what, **qcols(i)})
    ax.axvline(P.SESOI, color=INK2, linewidth=0.6, linestyle=(0, (3, 2)))
    zero(ax)
    rows.append({"class": "reference lines", "mean_delta": 0.0, "sesoi": P.SESOI})
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{CLASS_NAME[c]} ({ci[c]['n_exons']})" for c in names], fontsize=5.9)
    ax.set_ylim(4.2, -0.75)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("mean delta, PSI (95% interval)", fontsize=6)
    ax.set_xticks([-0.005, 0, 0.005, 0.01])
    ax.tick_params(axis="x", labelsize=5.6)
    tx.set_ylim(4.2, -0.75)
    tx.set_xlim(0, 1)
    tx.axis("off")
    save(fig, OUT, "RPM_mechanism", data=pd.DataFrame(rows), caption=CAPTIONS["RPM"](ci))


def _cap_rpm(ci) -> str:
    return f"""**Which mechanism: branch-point geometry or hnRNP C?**

Held-out cancers, junction reads; outside the claim's verdict (family RPM). Mean delta with a descriptive 95%
interval; dashed = SESOI. Shield exons have no acceptor T-run but a T-run of ≥ 20 nt upstream of their branch
point: shortening it cannot bring the branch point closer to the AG, so the geometry reading predicts no rise there
and a rise in R_up above it (RPM1), while the hnRNP C reading (hnRNP C binds long U-runs and represses Alu-derived
exons) predicts that shield exons rise too (RPM2) and that runs inside an antisense Alu rise more (RPM3). Long runs
= eligible runs of ≥ 20 nt (R_up, near or no branch point), compared within strata and run-length bins.

- {tline('RPM1', 'R_up minus shield')}
- {tline('RPM2', 'shield mean against 0, MMR-label permutation')}
- {tline('RPM3', 'runs in an antisense Alu minus other long runs')}
- None resolved: the data do not separate the two readings."""


# ============================================================ legends
TITLES = {
    "RP0": "The rising arm passes 9 of the 12 pre-registered tests; the claim is partial.",
    "RP1": "The rising arm separates MMR-deficient tumours better than matched exon sets, but the fitted model "
           "does not need it.",
    "RP2": "Junction reads confirm the rise in held-out cancers, carried by a minority of the exons.",
    "RP3": "The rise needs the T-run and, in held-out cancers, is rank-associated with the gap to the mapped branch point; the "
           "36-nt cut is not supported.",
    "RP4": "The rise exceeds long pyrimidine tracts, deep T-runs and Alu acceptors; poly-A and donor-side runs are "
           "not resolved.",
    "RP5": "Exon and tumour covariates leave the rise unchanged.",
    "RP6": "The arm separates mismatch repair from mutation load, but leiomyosarcoma contrasts move it as much.",
    "RP7": "MMR-deficient tumours rise above their own normal tissue, mostly in training cancers.",
    "RP8": "The rise does not follow dose.",
    "RPM": "The mechanism tests are not resolved.",
}


def legend(key: str, text: str) -> str:
    return _legend(key, TITLES[key], text)


CAPTIONS = {k: (lambda k, f: (lambda *a: tidy(legend(k, f(*a)))))(k, f) for k, f in
            {"RP0": _cap_rp0, "RP1": _cap_rp1, "RP2": _cap_rp2, "RP3": _cap_rp3, "RP4": _cap_rp4, "RP5": _cap_rp5,
             "RP6": _cap_rp6, "RP7": _cap_rp7, "RP8": _cap_rp8, "RPM": _cap_rpm}.items()}


def legends() -> None:
    """Copy each figure's legend (without its bullet list of tests) under the figure in the write-up."""
    if WRITEUP.exists():
        inject_legends(WRITEUP)


FIGS = {"RP0": RP0, "RP1": RP1, "RP2": RP2, "RP3": RP3, "RP4": RP4, "RP5": RP5, "RP6": RP6, "RP7": RP7, "RP8": RP8,
        "RPM": RPM, "legends": legends}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None, help=f"any of {', '.join(FIGS)}")
    a = ap.parse_args()
    for name in a.only or list(FIGS):
        FIGS[name]()
        log(f"{name} written")


if __name__ == "__main__":
    main()
