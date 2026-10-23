"""Step 63: D-77's step-wise test of the rising arm, in figures (house style, docs/FIGURE_GUIDELINES.md).

One figure per step, in the order of the argument in results/MSI_RISING_ARM_PROOF.md. Each title states that
step's own conclusion; each panel carries one message:
    RP0 overview: the eight steps as a chain (question, plain answer, marker), the mechanism and the verdict
    RP1 do models use the rising arm (held-out AUROC against matched exon sets, subtype-network SHAP, M4)?
    RP2 does it rise in junction reads, and how is the rise spread over its exons?
    RP3 does it need the T-run, is there a step at 36 nt, and does the gap association replicate in training?
    RP4 is it the microsatellite, or would another long tract do, and could the controls have been resolved?
    RP5 do exon and tumour covariates explain it?
    RP6 is it mismatch repair rather than mutation load, and is it unique to MSI among subtype contrasts?
    RP7 is it there within the same patient (tumour against its own normal tissue)?
    RP8 does it follow dose (MANTIS score, deletion count, run length)?
    RPM which mechanism (branch-point geometry or hnRNP C)?
Test ids and p values stay out of the axes; they are in the captions and in the table of all 40 tests that
the `legends` step writes into section 5 of the write-up (between marker comments).

Reads only what step 62 wrote with the real labels (results/subtype/msi_cis/code/proof/: tests.csv,
sensitivity.csv, descriptive.csv, verdict.json, rp6c_null_contrasts.csv) and its local caches
(data.nosync/.../msi/proof/tests/: builds, engines, tumours, pairs, counts). Nothing is re-tested. Where a figure
draws a distribution or a per-exon value that step 62 summarised (RP1's matched-set and random-removal AUROCs,
RPD1's permutation null, RP3's rank residuals, RP7's pair levels), it is recomputed with step 62's own code and
seeds, and the script stops unless it reproduces the committed numbers. The 95% intervals of the control classes
are descriptive, from RP2a's construction (studentized SE, Welch-Satterthwaite df); they are not tests.

Captions are figure legends ("**Figure RPn. Title.**", then (a), (b) ...); `legends` copies each legend
(without its bullet list of tests) under the figure in results/MSI_RISING_ARM_PROOF.md. Per-tumour and
per-patient plotted values go to data.nosync/ctsm_v1/figure_data/<stem>.csv (D-12); committed CSVs hold
aggregates.
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

from ctsm.figstyle import ACCENT, ACCENT_DARK, ACCENT_LIGHT, DOUBLE, GRID, INK, INK2, MUTED, NORMAL, TUMOUR, fp, save


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


def zero(ax, horizontal: bool = False, span: tuple | None = None, **kw) -> None:
    """The zero line; with `span` (data coordinates of the other axis) only over that range."""
    if span is not None:
        xy = ([span[0], span[1]], [0, 0]) if horizontal else ([0, 0], [span[0], span[1]])
        ax.plot(*xy, color=GRID, linewidth=0.7, zorder=0, **kw)
        return
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


# ============================================================ step 62's distributions, redrawn
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


# ============================================================ shared drawing helpers
FT = 6.5                                             # ticks, legends and text inside the axes (house: 6.5-8 pt)
YCHG = "Change in inclusion,\nMMRd − MMRp (PSI)"
XCHG = "Change in inclusion, MMRd − MMRp (PSI)"
SESOI_LAB = "smallest effect of interest (0.003)"
DASH = (0, (3, 2))
NAME = {"R_up": "rising arm", "R_down": "falling arm", "B sample": "background exons",
        "B_far": "distant branch point,\nno T-run", "ppt_far": "long pyrimidine tract", "alu": "Alu acceptor, no T-run",
        "polyA15": "poly-A run", "deep": "T-run deep in the intron", "donor": "T-run after the donor",
        "shield": "T-run upstream of\nthe branch point", "long, antisense Alu": "run in an\nantisense Alu",
        "long, other": "run outside an\nantisense Alu"}
# written out at first use in each legend
DEFINE = {
    "R_up": "rising arm = acceptor T-run ≥ 20 nt with the sequence-predicted branch point ≥ 36 nt from the AG",
    "R_down": "falling arm = acceptor T-run of 10–19 nt",
    "B sample": "background = exons without an acceptor T-run longer than 5 nt (the B sample)",
    "B_far": "distant branch point, no T-run = background exons with the predicted branch point ≥ 36 nt from the AG",
    "MMR": "MMRd = MMR-deficient (≥ 10 MC3 homopolymer deletions); MMRp = MMR-proficient",
    "change": ("change in inclusion = an exon's mean adjusted PSI in the MMRd tumours minus the MMRp mean of the "
               "same cancer, from junction reads only (rMATS JC); a tumour's adjusted PSI is its deviation from the "
               "cancer's median minus its mean deviation over background exons of the same baseline band")}


def sesoi_line(ax, n_rows: int) -> Line2D:
    """The smallest effect of interest as a vertical line over the rows only (the legend sits above them)."""
    ax.plot([P.SESOI, P.SESOI], [-0.5, n_rows - 0.5], color=INK2, linewidth=0.6, linestyle=DASH, zorder=1)
    return Line2D([], [], color=INK2, linewidth=0.6, linestyle=DASH, label=SESOI_LAB)


def dot_ci(ax, y: float, r: dict, col: str, filled: bool = True) -> None:
    """A class mean with its 95% interval, horizontal."""
    ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=col, linewidth=1.1, solid_capstyle="butt", zorder=2)
    ax.plot([r["mean"]], [y], marker="o", markersize=MS + 0.8, markerfacecolor=col if filled else "white",
            markeredgecolor=col, markeredgewidth=0.9, linestyle="", zorder=3)


def rows_axis(ax, labels: list, top_pad: float = 0.6, bottom_pad: float = 0.6) -> None:
    """Categorical y axis, first label on top."""
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_ylim(len(labels) - 1 + bottom_pad, -top_pad)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)


def heading(ax, text: str) -> None:
    """A panel heading that names the panel's subject (house style allows nothing else)."""
    ax.text(0.5, 1.03, text, transform=ax.transAxes, ha="center", va="bottom", fontsize=FT, color=INK2)


def ci_row(r: dict, **extra) -> dict:
    return {"n_exons": r["n_exons"], "mean_change": r["mean"], "ci_low": r["ci_low"], "ci_high": r["ci_high"], **extra}


# ============================================================ RP0 overview
STEPS = [
    ("1", "Do models use the rising arm?", ["RP1", "RPD1a"]),
    ("2", "Does it rise in junction reads?", ["RP2a", "RP2b", "RP2c"]),
    ("3", "Does it need both features?", ["RP3a", "RP3b", "RP3c", "RP3d", "RP3e"]),
    ("4", "Is it the microsatellite?", ["RP4a", "RP4b", "RP4c", "RP4d", "RP4e"]),
    ("5", "Does it survive confounders?", ["RP5a", "RP5b"]),
    ("6", "Is it mismatch repair?", ["RP6a", "RP6b", "RP6c"]),
    ("7", "Is it there within the same patient?", ["RP7"]),
    ("8", "Does it follow dose?", ["RP8a", "RP8b", "RP8c"]),
    ("M", "Mechanism: branch-point geometry or hnRNP C?", ["RPM1", "RPM2", "RPM3"]),
]


WORDS = "zero one two three four five six seven eight nine ten".split()


def step_mark(tests: list) -> str:
    v = [vd(i) for i in tests]
    return "filled" if all(x == "supported" for x in v) else ("hollow" if not any(x == "supported" for x in v) else "half")


def step_answers() -> dict:
    """The plain answer of each step, one key number at most, read from the committed tables."""
    auc = rp1_sets()["R_up"][0]
    m4ch = desc("M4 without its R_up inputs: AUROC change").value
    nul = pd.read_csv(CODE / "rp6c_null_contrasts.csv")
    above = nul[nul.s.abs() > T("RP6c")]
    over, n_sarc = WORDS[len(above)], WORDS[int(above.cohort.str.startswith("SARC").sum())]
    pr = pd.read_parquet(P.PROOF / "tests" / "pairs.parquet")
    pr = pr[pr.project.isin(set(pr.project[pr.mmr]))]
    md = pr[pr.mmr]
    n_tr = int(md.project.isin(P.TRAIN).sum())
    return {
        "1": f"Yes. Held out, it separates MMR-deficient tumours better than matched exon sets (AUROC {auc:.2f}). "
             f"Removing it from the fitted model costs no detectable accuracy ({sg(m4ch)}).",
        "2": f"Yes, held out: by {T('RP2a'):.3f} PSI on average, carried mostly by a few of its exons.",
        "3": "Partly. It needs the T-run, but a step at 36\u00a0nt is not supported, and the link to the "
             "branch-point gap is seen in held-out cancers only.",
        "4": "Partly. It exceeds pyrimidine tracts, deep T-runs and Alu acceptors; poly-A and donor-side runs are "
             "not resolved.",
        "5": "Yes. It remains after exon and tumour covariates.",
        "6": f"Partly. It separates MMR deficiency from mutation load, but {over} other subtype contrasts "
             f"({n_sarc} in sarcoma) move it more than MSI.",
        "7": f"Yes, tumour above its own normal tissue, but {n_tr} of the {len(md)} deficient patients are in "
             "training cancers.",
        "8": "Not shown: no relation to instability score, deletion count or run length was resolved.",
        "M": "Not separated: none of the three tests is resolved.",
    }


def _mark(ax, x: float, y: float, kind: str, r: float = 0.065) -> None:
    from matplotlib.patches import Circle, Wedge
    ax.add_patch(Circle((x, y), r, facecolor=ACCENT if kind == "filled" else "white", edgecolor=ACCENT,
                        linewidth=0.9, zorder=5))
    if kind == "half":
        ax.add_patch(Wedge((x, y), r, 90, 270, facecolor=ACCENT, edgecolor="none", zorder=6))
        ax.add_patch(Circle((x, y), r, facecolor="none", edgecolor=ACCENT, linewidth=0.9, zorder=7))


def RP0() -> None:
    import textwrap
    from matplotlib.patches import FancyBboxPatch
    ans = step_answers()
    W, gap, mx = DOUBLE - 0.08, 0.24, 0.03          # the saved PNG stays within 7.2 in after the tight-bbox pad
    bw = (W - 2 * mx - 3 * gap) / 4
    lh_q, lh_a, pad = 0.128, 0.118, 0.08
    wrap_q, wrap_a = 22, 29
    tw = lambda t, w: textwrap.wrap(t, w, break_on_hyphens=False)
    qs = {s: tw(f"{s}. {q}" if s != "M" else q, wrap_q) for s, q, _ in STEPS}
    an = {s: tw(ans[s], wrap_a) for s, _, _ in STEPS}
    box_h = lambda ss: max(2 * pad + len(qs[s]) * lh_q + 0.05 + len(an[s]) * lh_a for s in ss)
    h1, h2 = box_h("1234"), box_h("5678")
    v_lines = tw(verdict_text(), 86)
    h3 = max(box_h("M"), 2 * pad + 0.16 + len(v_lines) * lh_a + 0.05)
    vg, key_h = 0.36, 0.30
    H = key_h + h3 + vg + h2 + vg + h1 + 0.04
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    arrow = dict(arrowstyle="-|>", color=INK2, linewidth=0.8, mutation_scale=7, shrinkA=0, shrinkB=0)
    rows, geo = [], {}

    def box(x, y_top, w, h, s, edge=GRID, lw=0.8, ls="-", face="white"):
        ax.add_patch(FancyBboxPatch((x, y_top - h), w, h, boxstyle="round,pad=0,rounding_size=0.06", facecolor=face,
                                    edgecolor=edge, linewidth=lw, linestyle=ls, zorder=1))

    def step_box(s, x, y_top, h, edge=GRID, ls="-"):
        box(x, y_top, bw if s != "M" else bw, h, s, edge=edge, ls=ls)
        tests = dict((a, c) for a, _, c in STEPS)[s]
        kind = step_mark(tests)
        _mark(ax, x + pad + 0.065, y_top - pad - 0.06, kind)
        yq = y_top - pad
        for k, ln in enumerate(qs[s]):
            ax.text(x + pad + 0.19, yq - k * lh_q, ln, fontsize=7.5, fontweight="bold",
                    color=INK, va="top", ha="left")
        ya = yq - len(qs[s]) * lh_q - 0.05
        for k, ln in enumerate(an[s]):
            ax.text(x + pad, ya - k * lh_a, ln, fontsize=7, color=INK2, va="top", ha="left")
        geo[s] = (x, y_top, h)
        rows.append({"step": s, "question": dict((a, q) for a, q, _ in STEPS)[s], "answer": ans[s], "marker": kind,
                     "tests": "; ".join(tests), "verdicts": "; ".join(f"{i} {vd(i)}" for i in tests),
                     "outside_the_verdict": s == "M"})

    y1 = H - 0.04
    y2 = y1 - h1 - vg
    y3 = y2 - h2 - vg
    for k, s in enumerate("1234"):
        step_box(s, mx + k * (bw + gap), y1, h1)
    for k, s in enumerate("5678"):
        step_box(s, mx + k * (bw + gap), y2, h2)
    for row, ss, y_top, h in [(1, "1234", y1, h1), (2, "5678", y2, h2)]:
        for a, b in zip(ss[:-1], ss[1:]):
            xa = geo[a][0] + bw
            ax.annotate("", xy=(geo[b][0], y_top - h / 2), xytext=(xa, y_top - h / 2), arrowprops=arrow)
    # 4 -> 5: down from step 4, back along the gap, down into step 5
    x4, x5, ym = geo["4"][0] + bw / 2, geo["5"][0] + bw / 2, y1 - h1 - vg / 2
    ax.plot([x4, x4, x5], [y1 - h1, ym, ym], color=INK2, linewidth=0.8, solid_capstyle="butt")
    ax.annotate("", xy=(x5, y2), xytext=(x5, ym), arrowprops=arrow)
    # mechanism (outside the verdict) and the verdict
    step_box("M", mx, y3, h3, edge=MUTED, ls=(0, (3, 2)))
    vx = mx + bw + gap
    vw = W - mx - vx
    box(vx, y3, vw, h3, "V", edge=ACCENT_DARK, lw=1.2, face="#f4f3ee")
    ax.text(vx + pad, y3 - pad, f"Verdict: {VERDICT['verdict']} ({VERDICT['families']['RP']['supported']} of "
            f"{VERDICT['families']['RP']['tests']} required tests supported)", fontsize=7.5, fontweight="bold",
            color=INK, va="top", ha="left")
    for k, ln in enumerate(v_lines):
        ax.text(vx + pad, y3 - pad - 0.16 - 0.03 - k * lh_a, ln, fontsize=7, color=INK2, va="top", ha="left")
    x8 = geo["8"][0] + bw / 2
    ax.annotate("", xy=(x8, y3), xytext=(x8, y2 - h2), arrowprops=arrow)
    rows.append({"step": "verdict", "question": "claim", "answer": VERDICT["verdict"] + ": " + verdict_text(),
                 "marker": "", "tests": "; ".join(VERDICT["reasons"]), "verdicts": "", "outside_the_verdict": False})
    # key
    ky = key_h / 2
    for kx, kind, lab in [(mx, "filled", "every test of the step supported"), (mx + 1.95, "half", "some supported"),
                          (mx + 3.1, "hollow", "none resolved")]:
        _mark(ax, kx + 0.065, ky, kind)
        ax.text(kx + 0.17, ky, lab, fontsize=7, color=INK2, va="center", ha="left")
    kx = mx + 4.2
    ax.plot([kx, kx + 0.3], [ky, ky], color=MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(kx + 0.38, ky, "outside the verdict", fontsize=7, color=INK2, va="center", ha="left")
    save(fig, OUT, "RP0_overview", data=pd.DataFrame(rows), caption=CAPTIONS["RP0"]())


def verdict_text() -> str:
    return ("The rise holds out of sample but is small and carried by a few exons. Not shown: that it beats "
            "poly-A and donor-side runs, that it is unique to MSI, the 36-nt cut, a dose relation, or the mechanism.")


def _cap_rp0() -> str:
    fams = VERDICT["families"]
    fs = lambda f: f"{fams[f]['supported']} of {fams[f]['tests']}"
    lines = "\n".join(f"- Step {s} ({q}): " + "; ".join(f"{i} {vd(i)}" for i in tests) + "."
                      for s, q, tests in STEPS)
    return f"""**Does the rising arm pass the pre-registered test, step by step?**

Pre-registered as D-77 (config/followup.yaml: msi_rising_proof; commit 5471b27) before any junction-read value was
compared by MMR status. The claim: in MMR-deficient tumours, exons after an acceptor T-run of ≥ 20 nt with the
sequence-predicted branch point ≥ 36 nt from the AG (the rising arm) are more included than in MMR-proficient
tumours of the same cancer, and the rise needs both features. Confirmatory tests use the 30 held-out cancers (50
MMR-deficient tumours in 16 of them) and junction reads only (rMATS JC), in PSI units. The eight steps are read in
order; each box gives the step's question, a plain answer and a marker: filled = every test of the step supported
(predicted sign and Benjamini–Hochberg q < 0.05 within its family), half = some supported, hollow = none resolved.
The mechanism box (dashed) is outside the claim's verdict. The verdict box gives the claim's verdict by the
pre-registered rule, under which "established" needs all 12 tests of family RP. Every test, with its effect, p and q, is in section 5 of
results/MSI_RISING_ARM_PROOF.md and in tests.csv.

{lines}
- Families: RP (verdict) {fs('RP')} supported; RPS (secondary, held out) {fs('RPS')}; RPP (not held out) {fs('RPP')};
  RPC (falling arm) {fs('RPC')}; RPD (training) {fs('RPD')}; RPM (mechanism) {fs('RPM')}.
- Claim verdict: {VERDICT['verdict']}; not met: {unmet_short()}; no refuting clause holds."""


def unmet_short() -> str:
    """The claim verdict's unmet conditions, shortened when they all read '<id> not resolved (needed: supported)'."""
    r = VERDICT["reasons"]
    if r and all(x.endswith("not resolved (needed: supported)") for x in r):
        ids = [x.split()[0] for x in r]
        return f"{', '.join(ids[:-1])} and {ids[-1]} not resolved (needed: supported)" if len(ids) > 1 else r[0]
    return "; ".join(r)


# ============================================================ RP1 models
def RP1() -> None:
    from matplotlib.patches import Patch, Rectangle
    sh, sets, m4 = rpd1_null(), rp1_sets(), m4_removals()
    fig, axs = plt.subplots(1, 3, figsize=(DOUBLE, 2.0), gridspec_kw={"wspace": 0.6, "width_ratios": [1.3, 0.7, 1.05]})
    rows = []
    # a. held-out AUROC of each arm against 1,000 matched exon sets
    ax = axs[0]
    for y, (arm, i, col) in enumerate([("R_up", "RP1", ACCENT), ("R_down", "RPC1", INK2)]):
        a = sets[arm]
        lo, hi, med = *np.percentile(a[1:], [2.5, 97.5]), float(np.median(a[1:]))
        ax.add_patch(Rectangle((lo, y - 0.17), hi - lo, 0.34, facecolor=GRID, edgecolor="none", zorder=1))
        ax.plot([med, med], [y - 0.17, y + 0.17], color=MUTED, linewidth=1.0, zorder=2)
        ax.plot([a[0]], [y], marker="o", markersize=MS + 1, color=col, linestyle="", zorder=3)
        rows.append({"panel": "a", "arm": arm, "class_auroc": a[0], "matched_median": med, "matched_2.5": lo,
                     "matched_97.5": hi, **qcols(i)})
    ax.plot([0.5, 0.5], [-0.45, 1.45], color=INK2, linewidth=0.6, linestyle=(0, (1, 1.5)), zorder=0)
    rows.append({"panel": "a", "arm": "reference line", "class_auroc": 0.5})
    rows_axis(ax, ["rising arm", "falling arm"], top_pad=1.7, bottom_pad=0.5)
    ax.set_xlim(0.25, 0.8)
    ax.set_xlabel("Held-out AUROC, MMRd against MMRp")
    ax.legend(handles=[Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=INK2, linestyle="",
                              markersize=MS + 1, label="the arm's exons"),
                       Patch(facecolor=GRID, edgecolor="none", label="middle 95% of 1,000 matched sets"),
                       Line2D([], [], color=MUTED, linewidth=1.0, label="their median")],
              loc="upper left", fontsize=FT, handlelength=1.2, borderaxespad=0.1, labelspacing=0.3)
    heading(ax, "Held-out cancers")
    # b. subtype-network SHAP directions (training cancers)
    ax = axs[1]
    for x, (a, i) in enumerate([("R_up", "RPD1a"), ("R_down", "RPD1b")]):
        k, n_ = sh["obs"][a]
        lo, hi = np.percentile(sh["null"][a], [2.5, 97.5])
        ax.add_patch(Rectangle((x - 0.36, lo), 0.72, hi - lo, facecolor=GRID, edgecolor="none", zorder=0))
        ax.bar(x, k / n_, width=0.46, color=ACCENT if a == "R_up" else INK2, zorder=2)
        ax.text(x, k / n_ + 0.03, f"{k}/{n_}", ha="center", va="bottom", fontsize=FT, color=INK)
        rows.append({"panel": "b", "arm": a, "positive": k, "units": n_, "share": k / n_, "null_2.5": lo,
                     "null_97.5": hi, "null_sd": float(np.std(sh["null"][a])), **qcols(i)})
    ax.axhline(0.5, color=INK2, linewidth=0.6, linestyle=(0, (1, 1.5)), zorder=1)
    rows.append({"panel": "b", "arm": "reference line", "share": 0.5})
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["rising\narm", "falling\narm"])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Share of inputs pushing\ntowards MSI")
    heading(ax, "Training cancers")
    # c. M4 refitted without its R_up inputs against 200 random removals
    ax = axs[2]
    h, e = _hist(ax, m4["random"], np.linspace(-0.02, 0.012, 33))
    ch = m4["au2"] - m4["au1"]
    ax.plot([ch, ch], [0, h.max() * 1.08], color=ACCENT, linewidth=1.3, zorder=3)
    ax.set_ylim(0, h.max() * 1.6)
    ax.set_xlabel("Change in held-out AUROC")
    ax.set_ylabel("Random removals")
    ax.set_xticks([-0.02, -0.01, 0, 0.01])
    ax.legend(handles=[Line2D([], [], color=ACCENT, linewidth=1.3, label=f"its {m4['n_up']} rising-arm inputs out"),
                       Patch(facecolor=LIGHT, edgecolor="none", label=f"{m4['n_up']} random inputs out (200 draws)")],
              loc="upper left", fontsize=FT, handlelength=1.2, borderaxespad=0.1, labelspacing=0.3)
    heading(ax, "Fitted model M4")
    for lo_, n_ in zip(e[:-1], h):
        if n_:
            rows.append({"panel": "c", "arm": "random non-R_up inputs", "bin_left": lo_, "bin_width": e[1] - e[0],
                         "draws": int(n_)})
    pct = desc("M4 without its R_up inputs: percentile of the AUROC change among random removals (A-25)").value
    rows.append({"panel": "c", "arm": "R_up inputs removed", "auroc_full": m4["au1"], "auroc_without": m4["au2"],
                 "change": ch, "percentile_among_random": pct})
    panel_letters([(axs[0], -0.36), (axs[1], -0.62), (axs[2], -0.3)], y=1.09)
    save(fig, OUT, "RP1_models", data=pd.DataFrame(rows), caption=CAPTIONS["RP1"](sh, sets, m4, pct))


def _cap_rp1(sh, sets, m4, pct) -> str:
    rr = desc(f"M4 without {m4['n_up']} random non-R_up repeat inputs: AUROC change (exploratory, A-25)")
    bt = desc("M4 without its R_up inputs: AUROC change")
    co = lambda g: desc("M4 standardised coefficient: share > 0", g).value
    return f"""**Do models use the rising arm?**

{DEFINE['R_up'][0].upper() + DEFINE['R_up'][1:]}; {DEFINE['R_down']}. {DEFINE['MMR']}. (a) Held-out cancers (7,685
tumours, 50 MMRd): AUROC for MMRd against MMRp tumours of each arm's level (a tumour's mean PSI deviation from its
cancer's median over the arm's exons, as a percentile within the cancer); dot = the arm; grey bar = middle 95% of
1,000 sets of background exons matched to the arm on baseline band and coverage tertile within each cancer, which
carry any tumour-wide shift as well; tick = their median; dotted line = 0.5. The falling arm is predicted below
0.5. (b) Training cancers (COAD, STAD, UCEC): the subtype networks (SE view, D-38) explained with SHAP. For each
rising- or falling-arm input and cohort, the direction is the Spearman correlation of the input's value with its
SHAP over a fold's MSI test tumours, averaged over the folds; bar = share of input × cohort units whose higher
inclusion pushes towards MSI; grey band behind each bar = the range expected by chance, the middle 95%
(2.5th–97.5th percentiles) of that share over 1,000 permutations of the inputs' rows among each fold's MSI test
tumours (A-22); dotted line = 0.5. (c) M4, step 54's logistic regression on all
repeat exons (fitted in the training cancers; held-out AUROC {m4['au1']:.3f}), refitted at the same C without its
{m4['n_up']} rising-arm inputs (line) and without {m4['n_up']} inputs drawn at random from its other repeat inputs,
200 times (grey; exploratory, A-25).

- RP1, rising-arm AUROC minus the matched median: {sg(T('RP1'), 3)} (AUROC {sets['R_up'][0]:.3f} against
  {np.median(sets['R_up'][1:]):.3f}), p = {fp(T('RP1', 'p'))}, {vd('RP1')}; RPC1, falling arm: {sg(T('RPC1'), 3)}
  ({sets['R_down'][0]:.3f} against {np.median(sets['R_down'][1:]):.3f}), p = {fp(T('RPC1', 'p'))}, {vd('RPC1')}.
- RPD1a, share of rising-arm directions minus 0.5: {sg(T('RPD1a'), 3)} ({sh['obs']['R_up'][0]} of
  {sh['obs']['R_up'][1]}), p = {fp(T('RPD1a', 'p'))}, {vd('RPD1a')}; RPD1b, falling arm: {sg(T('RPD1b'), 3)}
  ({sh['obs']['R_down'][0]} of {sh['obs']['R_down'][1]}), p = {fp(T('RPD1b', 'p'))}, {vd('RPD1b')}.
- M4 without its rising-arm inputs: AUROC change {sg(m4['au2'] - m4['au1'])} (paired bootstrap over test tumours
  with the fits held fixed, {num(bt.ci_low)} to {num(bt.ci_high)}); random removals: median {sg(rr.value)},
  2.5th–97.5th percentiles {num(rr.ci_low)} to {num(rr.ci_high)}; the rising-arm removal lies at the {pct:.1f}th
  percentile. No test; no loss was detected (test-set interval {num(bt.ci_low, 3)} to {sg(bt.ci_high, 3)}). M4's
  standardised coefficients are positive
  for {100 * co('R_up'):.0f}% of its rising-arm inputs and {100 * co('R_down'):.0f}% of its falling-arm inputs."""


# ============================================================ RP2 junction reads
def RP2() -> None:
    from matplotlib.patches import Patch
    _, H = state()
    d = delta0()
    order = ["R_up", "B sample", "R_down"]
    ci = {c: class_interval(c) for c in order}
    same(ci["R_up"]["ci_low"], T("RP2a", "ci_low"), "RP2a interval")
    same(ci["R_down"]["ci_low"], T("RPC2a", "ci_low"), "RPC2a interval")
    up = sel_of("R_up")
    x = d[up]
    sym = H.tab.symbol.to_numpy()[up]
    o = np.argsort(-x)
    x, sym = x[o], sym[o]
    bg = d[sel_of("B sample")]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE, 2.05), gridspec_kw={"width_ratios": [1.0, 1.9], "wspace": 0.42})
    rows = []
    # a. class means with 95% intervals
    for y, c in enumerate(order):
        col = ACCENT if c == "R_up" else (INK2 if c == "R_down" else MUTED)
        dot_ci(ax, y, ci[c], col)
        rows.append({"panel": "a", "class": c, **ci_row(ci[c]),
                     **(qcols("RP2a") if c == "R_up" else qcols("RPC2a") if c == "R_down" else {})})
    hs = sesoi_line(ax, len(order))
    zero(ax, span=(-0.5, len(order) - 0.5))
    rows.append({"panel": "a", "class": "reference lines", "mean_change": 0.0, "sesoi": P.SESOI})
    rows_axis(ax, [f"{NAME[c]}\n({ci[c]['n_exons']:,} exons)" for c in order], top_pad=1.5)
    ax.set_xlim(-0.0065, 0.0125)
    ax.set_xticks([-0.005, 0, 0.005, 0.01])
    ax.set_xlabel(XCHG)
    ax.legend(handles=[hs], loc="upper left", fontsize=FT, handlelength=1.6, borderaxespad=0.1)
    # b. the rising-arm exons, ranked
    n = len(x)
    bx.bar(np.arange(n), x, width=0.8, color=[ACCENT_DARK if k < 10 else LIGHT for k in range(n)], zorder=2)
    qb = np.quantile(bg, 1 - (np.arange(n) + 0.5) / n)          # background exons at the same quantiles
    bx.step(np.arange(n), qb, where="mid", color=INK2, linewidth=0.8, zorder=4)
    zero(bx, horizontal=True)
    cum = np.cumsum(x) / x.sum()
    lead = dict(arrowstyle="-", linewidth=0.4, color=INK2, shrinkA=0, shrinkB=1)
    for k, yt in zip(range(3), [0.092, 0.08, 0.068]):
        bx.annotate(sym[k], (k, x[k]), xytext=(5, yt), textcoords="data", fontsize=FT, va="center", ha="left",
                    color=INK, arrowprops=lead)
    bx.set_xlim(-1, n + 0.5)
    bx.set_ylim(min(-0.03, float(x.min()) - 0.004), 0.1)
    bx.set_xticks([0, 9, 30, 60, 91])
    bx.set_xticklabels(["1", "10", "31", "61", "92"])
    bx.set_xlabel("Rising-arm exons, ranked by change")
    bx.set_ylabel(YCHG)
    bx.legend(handles=[Patch(facecolor=ACCENT_DARK, label="10 largest"), Patch(facecolor=LIGHT, label="the other 82"),
                       Line2D([], [], color=INK2, linewidth=0.8, label="background exons at the same rank")],
              loc="upper right", fontsize=FT, handlelength=1.2, borderaxespad=0.1, labelspacing=0.3)
    for k, (v, s_) in enumerate(zip(x, sym)):
        rows.append({"panel": "b", "class": "R_up", "rank": k + 1, "symbol": s_, "change": v, "top10": k < 10,
                     "cumulative_share_of_sum": cum[k]})
    rows += [{"panel": "b", "class": "B sample, same quantile", "rank": k + 1, "change": v} for k, v in enumerate(qb)]
    rows.append({"panel": "b", "class": "R_up", "mean_change": float(np.mean(x)), "median_change": float(np.median(x)),
                 "share_positive": float(np.mean(x > 0)), "top10_share_of_sum": float(cum[9]),
                 "mean_without_top10": float(np.mean(x[10:]))})
    panel_letters([(ax, -0.62), (bx, -0.14)], y=1.04)
    save(fig, OUT, "RP2_junction_reads", data=pd.DataFrame(rows), caption=CAPTIONS["RP2"](ci, x, len(bg)))


def _cap_rp2(ci, x, nbg) -> str:
    cum = np.cumsum(x) / x.sum()
    gl = sens("gene level (delta averaged within gene first)", "RP2a")
    md = sens("median instead of mean delta", "RP2a")
    return f"""**Does the rising arm rise in junction reads, and how is the rise spread over its exons?**

Held-out cancers. {DEFINE['MMR']}; {DEFINE['change']}. Exons need ≥ 5 MMRd
tumours measured. (a) Mean change of each class with its 95% interval (studentized SE over MMR labels permuted within
cancer, Welch–Satterthwaite df; the rising- and falling-arm intervals are those of the tests, the background's is
descriptive). {DEFINE['R_up'][0].upper() + DEFINE['R_up'][1:]}; {DEFINE['R_down']}; {DEFINE['B sample']}.
(b) The {len(x)} rising-arm exons ranked by their change; dark = the 10 largest, which hold
{100 * (np.cumsum(x) / x.sum())[9]:.0f}% of the summed change; line = background exons'
changes at the same quantiles (rank k of {len(x)} against the background's (1 − (k − 0.5) / {len(x)})
quantile; {nbg:,} exons); the three largest labelled with their gene.

- RP2a, rising-arm mean change: {sg(T('RP2a'))} (95% interval {num(T('RP2a', 'ci_low'))} to
  {num(T('RP2a', 'ci_high'))}), p = {fp(T('RP2a', 'p'))}, {vd('RP2a')}; the whole interval lies above the smallest
  effect of interest. RPC2a, falling arm: {sg(T('RPC2a'))} ({num(T('RPC2a', 'ci_low'))} to
  {num(T('RPC2a', 'ci_high'))}), p = {fp(T('RPC2a', 'p'))}, {vd('RPC2a')}.
- RP2c, share of rising-arm exons that rise (studentized change > c = 3.02) minus that of 5 matched background
  exons each: {sg(T('RP2c'), 3)}, p = {fp(T('RP2c', 'p'))}, {vd('RP2c')}; RP2b, per tumour (rising-arm level, MMRd
  − MMRp): {sg(T('RP2b'))}, p = {fp(T('RP2b', 'p'))}, {vd('RP2b')}.
- The class mean is carried by a minority of exons: median {sg(float(np.median(x)))} (sensitivity p =
  {fp(md.p)}); {100 * np.mean(x > 0):.0f}% of rising-arm exons above 0; the 10 largest hold {100 * cum[9]:.0f}% of the
  summed change, and the mean of the other {len(x) - 10} is {sg(float(np.mean(x[10:])))}, below the smallest effect
  of interest. Gene level (change averaged within gene first): {sg(gl.estimate)} ({num(gl.ci_low)} to
  {num(gl.ci_high)})."""


# ============================================================ RP3 tract and branch point
D_BINS = [(18, 30, "18–29"), (30, 36, "30–35"), (36, 50, "36–49"), (50, 70, "50–69"), (70, 102, "70–101")]


def training_ctx():
    """Step 62's training-cancer exon engine (the context of RPD2a-RPD3b)."""
    if "TR" not in _ST:
        S, _ = state()
        _ST["TR"] = S.xctx(P.SPECS["primary"], "training")
    return _ST["TR"]


def gap_residuals(ctx, test: str):
    """RP3e's (or RPD3b's) exons and rank residuals, rebuilt with step 62's selection; stops unless the
    partial Spearman equals the committed estimate."""
    from scipy.stats import rankdata
    t, d = ctx.tab, ctx.D[0].astype(float)
    g = t.g_exp.to_numpy(float)
    L = t.geo_L.fillna(t.acc_T).to_numpy(float)
    o = t.geo_o.fillna(t.run_near).to_numpy(float)
    sel = (t.acc_T >= P.LONG).to_numpy() & (g > 0) & ctx.elig & np.isfinite(L) & np.isfinite(o) & np.isfinite(d)
    x, gg, LL, oo = d[sel], g[sel], L[sel], o[sel]
    rho = P.partial_spearman(x, gg, np.column_stack([rankdata(LL), rankdata(oo)]))
    same(rho, T(test), test)
    Hm = P.residualize(np.column_stack([rankdata(LL), rankdata(oo)]), len(x))
    return pd.DataFrame({"feature_id": t.feature_id.to_numpy()[sel], "cls": t.cls.to_numpy()[sel], "delta": x,
                         "g_exp": gg, "L": LL, "o": oo, "rank_residual_delta": Hm @ rankdata(x),
                         "rank_residual_gap": Hm @ rankdata(gg)}), rho


def RP3() -> None:
    _, H = state()
    t, d = H.tab, delta0()
    ci = {c: class_interval(c) for c in ["R_up", "B_far"]}
    held, rho_h = gap_residuals(H, "RP3e")
    train, rho_t = gap_residuals(training_ctx(), "RPD3b")
    dist = t.d_bp.to_numpy(float)
    long_ = (t.acc_T >= P.LONG).to_numpy() & np.isfinite(dist) & H.elig & np.isfinite(d)
    same(int(long_.sum()), int(TESTS.at["RP3c", "n"].split()[0]), "RP3b/RP3c exons with a heuristic branch point")
    fig = plt.figure(figsize=(DOUBLE - 0.3, 4.3))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.15], hspace=0.62)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.75], wspace=0.5)
    bot = outer[1].subgridspec(1, 4, width_ratios=[0.35, 1.0, 1.0, 0.35], wspace=0.3)
    ax, bx = fig.add_subplot(top[0]), fig.add_subplot(top[1])
    c1, c2 = fig.add_subplot(bot[1]), fig.add_subplot(bot[2])
    rows = []
    # a. the rise needs the T-run
    for y, c in enumerate(["R_up", "B_far"]):
        dot_ci(ax, y, ci[c], ACCENT if c == "R_up" else INK2)
        rows.append({"panel": "a", "class": c, **ci_row(ci[c]), **(qcols("RP3a") if c == "B_far" else {})})
    hs = sesoi_line(ax, 2)
    zero(ax, span=(-0.5, 1.5))
    rows_axis(ax, [f"{NAME[c]}\n({ci[c]['n_exons']:,} exons)" for c in ["R_up", "B_far"]], top_pad=1.3)
    ax.set_xlim(-0.004, 0.0125)
    ax.set_xticks([0, 0.005, 0.01])
    ax.set_xlabel(XCHG.replace(", MMRd", ",\nMMRd"))
    ax.legend(handles=[hs], loc="upper left", fontsize=FT, handlelength=1.6, borderaxespad=0.1)
    # b. no step at 36 nt: change by predicted branch-point distance, runs >= 20 nt
    rng = np.random.default_rng(0)
    xv, dv = dist[long_], d[long_]
    lo_y, hi_y = min(-0.035, float(np.floor(dv.min() / 0.02) * 0.02)), 0.1   # every exon drawn, none clipped
    for k, (lo, hi, lab_) in enumerate(D_BINS):
        m = (xv >= lo) & (xv < hi)
        col = ACCENT_LIGHT if lo >= P.FAR else LIGHT
        bx.scatter(k + rng.uniform(-0.2, 0.2, m.sum()), dv[m], s=7, color=col, edgecolors="none", zorder=2)
        mu = float(dv[m].mean()) if m.sum() >= 5 else np.nan
        if np.isfinite(mu):
            bx.plot([k - 0.3, k + 0.3], [mu, mu], color=ACCENT_DARK if lo >= P.FAR else INK2, linewidth=1.4, zorder=3)
        rows.append({"panel": "b", "branch_point_bin_nt": lab_, "n_exons": int(m.sum()), "mean_change": mu,
                     "rising_arm_side": lo >= P.FAR})
    if sum(r["n_exons"] for r in rows if r["panel"] == "b") != int(long_.sum()):
        raise SystemExit("RP3b bins do not cover every long-run exon with a predicted branch point")
    if not (dv.min() >= lo_y and dv.max() <= hi_y):
        raise SystemExit("RP3b: an exon lies outside the drawn range")
    rows += [{"panel": "b", "distance_nt": a_, "change": b_, "class": c_, "gene": g_} for a_, b_, c_, g_ in
             zip(xv, dv, t.cls.to_numpy()[long_], t.symbol.to_numpy()[long_])]
    bx.plot([1.5, 1.5], [lo_y, hi_y], color=INK, linewidth=0.7, linestyle=DASH, zorder=1)
    zero(bx, horizontal=True)
    bx.set_xticks(range(len(D_BINS)))
    bx.set_xticklabels([b[2] for b in D_BINS])
    bx.tick_params(axis="x", length=0)
    bx.set_xlim(-0.6, len(D_BINS) - 0.4)
    bx.set_ylim(lo_y, hi_y + 0.04)
    bx.set_yticks(np.round(np.arange(np.ceil(lo_y / 0.04) * 0.04, hi_y + 1e-9, 0.04), 2))
    bx.spines["left"].set_bounds(lo_y, hi_y)
    bx.set_xlabel("Predicted branch point to AG (nt)")
    bx.set_ylabel(YCHG)
    bx.legend(handles=[Line2D([], [], color=INK, linewidth=0.7, linestyle=DASH, label="the rule's 36-nt cut"),
                       Line2D([], [], color=INK2, linewidth=1.4, label="mean of the bin")],
              loc="upper left", ncol=2, fontsize=FT, handlelength=1.6, borderaxespad=0.1, columnspacing=1.2)
    # c. the gap association, held out and in training
    lim = max(np.abs(np.r_[held.rank_residual_gap, train.rank_residual_gap]).max(),
              np.abs(np.r_[held.rank_residual_delta, train.rank_residual_delta]).max()) * 1.08
    for cx, df, rho, name, col in [(c1, held, rho_h, "Held-out cancers", ACCENT), (c2, train, rho_t, "Training cancers", INK2)]:
        rx, rg = df.rank_residual_delta.to_numpy(), df.rank_residual_gap.to_numpy()
        cx.scatter(rg, rx, s=7, color=col, edgecolors="none", zorder=2)
        b_ = float(rx @ rg / (rg @ rg))
        cx.plot([-lim, lim], [-lim * b_, lim * b_], color=INK, linewidth=0.7, zorder=3)
        cx.axhline(0, color=GRID, linewidth=0.6, zorder=0)
        cx.axvline(0, color=GRID, linewidth=0.6, zorder=0)
        cx.set_xlim(-lim, lim)
        cx.set_ylim(-lim, lim)
        cx.set_aspect("equal")
        heading(cx, f"{name} ({len(df)} exons)\npartial ρ = {rho:.2f}")
        cx.set_xlabel("Gap to mapped branch\npoint (rank residual)")
        rows.append({"panel": "c", "class": name, "n_exons": len(df), "partial_spearman": rho, "ols_slope_line": b_,
                     **qcols("RP3e" if name.startswith("Held") else "RPD3b")})
    c1.set_ylabel("Change in inclusion\n(rank residual)")
    c2.tick_params(axis="y", labelleft=False)
    panel_letters([(ax, -0.95), (bx, -0.22), (c1, -0.45)], y=1.1)
    rows += [{"panel": "c", "class": "per exon, held out", **r} for r in held.drop(columns=["feature_id"]).to_dict("records")]
    rows += [{"panel": "c", "class": "per exon, training", **r} for r in train.drop(columns=["feature_id"]).to_dict("records")]
    save(fig, OUT, "RP3_branch_point", data=pd.DataFrame(rows), caption=CAPTIONS["RP3"](held, ci, rows))


def _cap_rp3(local, ci, rows) -> str:
    from scipy.stats import rankdata
    x, g, L, o = (local[c].to_numpy(float) for c in ["delta", "g_exp", "L", "o"])
    ps = lambda m: P.partial_spearman(x[m], g[m], np.column_stack([rankdata(L[m]), rankdata(o[m])]))
    loo = [ps(np.arange(len(x)) != i) for i in range(len(x))]
    d_exp = g + L + o - 1
    top_all = np.argsort(-x)[:6]
    up = np.flatnonzero(local.cls.eq("R_up").to_numpy())
    top_up = up[np.argsort(-x[up])][:6]
    d60 = sens("R_up cut at d_bp > 60", "RP3b")
    d50 = sens("R_up cut at d_bp >= 50", "RP3b")
    bins = "; ".join(f"{r['branch_point_bin_nt']} nt: {sg(r['mean_change'])} ({r['n_exons']})" for r in rows
                     if r.get("panel") == "b" and "branch_point_bin_nt" in r)
    return f"""**Does the rise need the T-run, is there a step at 36 nt, and does the gap association replicate?**

{DEFINE['MMR']}; {DEFINE['change']}. (a) Held-out cancers: mean change with its
95% interval for the rising arm ({DEFINE['R_up']}) and for background exons whose predicted branch point is
also ≥ 36 nt from the AG, without a T-run; dashed = the smallest effect of interest. (b) Held-out cancers, exons with an acceptor
T-run ≥ 20 nt and a predicted branch point (rising arm right of the dashed line, the rule's 36-nt cut; long runs
with a nearer branch point left of it): change by the predicted branch point's distance from the AG, in D-70's
exploratory bins (A-16); dots = exons, bars = bin means (≥ 5 exons). Descriptive. (c) Exons with an acceptor
T-run ≥ 20 nt and an experimentally mapped branch point upstream of the run (BPHunter, lifted to GRCh38), in the
held-out cancers (left; {len(local)} exons, {len(up)} of them in the rising arm) and in the training cancers
(right): rank of the change against rank of the gap from the run to the mapped branch point, each with run
length and offset (ranks) regressed out; line = least squares. For these exons the branch point's distance from
the AG is the offset of the run plus its length plus the gap, minus 1, so a cut on that distance mixes run
length and position with the gap. The scatter is descriptive; the tests are listed below.

- RP3a, rising arm minus distant branch point without a run (within strata): {sg(T('RP3a'))}, p =
  {fp(T('RP3a', 'p'))}, {vd('RP3a')}.
- The 36-nt cut: RP3b (rising arm minus long runs with a near predicted branch point) {sg(T('RP3b'))}, p =
  {fp(T('RP3b', 'p'))}, {vd('RP3b')}; RP3d (mapped branch point ≥ 36 nt minus < 36 nt) {sg(T('RP3d'))}, p =
  {fp(T('RP3d', 'p'))}, {vd('RP3d')}. Pre-specified sensitivity cuts for the rising arm (no q): ≥ 50 nt
  {sg(d50.estimate)}, p = {fp(d50.p)}; > 60 nt {sg(d60.estimate)}, p = {fp(d60.p)}. Bin means in (b): {bins}.
- RP3e, partial Spearman of the change with the gap given run length and offset, held out: {sg(T('RP3e'), 3)},
  p = {fp(T('RP3e', 'p'))}, {vd('RP3e')} (the branch-point test of the verdict); leaving out one exon at a time,
  {min(loo):.2f} to {max(loo):.2f}. RP3c, the same with the predicted branch point (139 exons): {sg(T('RP3c'), 3)},
  p = {fp(T('RP3c', 'p'))}, {vd('RP3c')}.
- RPD3b, the same test in the training cancers ({TESTS.at['RPD3b', 'n']}, 150 MMRd tumours): {sg(T('RPD3b'), 3)},
  p = {fp(T('RPD3b', 'p'))}, {vd('RPD3b')}. The gap association does not replicate there.
- Of the 6 largest rises among the {len(local)} held-out exons in (c), {int((d_exp[top_all] < P.FAR).sum())} have a
  mapped branch point < 36 nt from the AG ({int((d_exp[top_up] < P.FAR).sum())} of the 6 largest among the
  rising-arm exons)."""


# ============================================================ RP4 controls
RP4_ORDER = [("ppt_far", "RP4b"), ("alu", "RP4e"), ("deep", "RP4c"), ("polyA15", "RP4a"), ("donor", "RP4d")]


def RP4() -> None:
    from matplotlib.patches import Patch
    pw = contrast_power()
    labels = [NAME[c] for c, _ in RP4_ORDER]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE * 0.82, 1.75), sharey=True, gridspec_kw={"wspace": 0.12})
    rows = []
    for y, (c, i) in enumerate(RP4_ORDER):
        ok = vd(i) == "supported"
        ax.barh(y, T(i), height=0.56, color=ACCENT if ok else "white", edgecolor=ACCENT, linewidth=0.9, zorder=2)
        bx.barh(y, pw[i]["mdd_80"], height=0.56, color=LIGHT if ok else MUTED, edgecolor="none", zorder=2)
        rows.append({"control": c, "label": NAME[c], **qcols(i), "contrast_n": TESTS.at[i, "n"],
                     "class_perm_null_sd": pw[i]["null_sd"], "detectable_difference_80pct": pw[i]["mdd_80"],
                     "power_if_control_flat": pw[i]["power_if_control_flat"]})
    rows_axis(ax, labels, top_pad=0.6, bottom_pad=0.6)
    zero(ax)
    ax.set_xlim(0, 0.0105)
    ax.set_xticks([0, 0.005, 0.01])
    ax.set_xlabel("Rising arm minus control (PSI)")
    ax.legend(handles=[Patch(facecolor=ACCENT, edgecolor=ACCENT, label="supported"),
                       Patch(facecolor="white", edgecolor=ACCENT, label="not resolved")],
              loc="lower right", fontsize=FT, handlelength=1.2, borderaxespad=0.1)
    rise = T("RP2a")
    bx.axvline(rise, color=ACCENT, linewidth=1.0, zorder=3)
    bx.set_xlim(0, 0.0115)
    bx.set_xticks([0, 0.005, 0.01])
    bx.set_xlabel("Difference detectable with 80% power (PSI)")
    bx.tick_params(axis="y", length=0)
    bx.legend(handles=[Line2D([], [], color=ACCENT, linewidth=1.0, label=f"whole rise ({rise:.4f})"),
                       Patch(facecolor=LIGHT, edgecolor="none", label="contrast supported"),
                       Patch(facecolor=MUTED, edgecolor="none", label="not resolved")],
              loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=FT, handlelength=1.4,
              columnspacing=1.0, labelspacing=0.3, borderaxespad=0.1)
    rows.append({"control": "reference line", "label": "rising-arm mean change (RP2a)", "estimate": rise})
    panel_letters([(ax, -0.62), (bx, -0.06)], y=1.04)
    save(fig, OUT, "RP4_controls", data=pd.DataFrame(rows), caption=CAPTIONS["RP4"](pw))


def ctl_n(i: str) -> str:
    """'85 vs 61 exons (effective controls 27)' -> '61 controls, 27 effective'."""
    import re
    m = re.match(r"\d+ vs (\d+) exons \(effective controls (\d+)\)", TESTS.at[i, "n"])
    return f"{m.group(1)} controls, {m.group(2)} effective"


def _cap_rp4(pw) -> str:
    own = lambda g: desc("control class: own mean delta against 0 (label_perm p)", g)
    na = sens("acceptor not in an Alu", "RP3a")
    return f"""**Is it the microsatellite, or would another long tract do?**

Held-out cancers, junction reads. {DEFINE['R_up'][0].upper() + DEFINE['R_up'][1:]}. Each control has no acceptor
T-run longer than 7 nt and shares one other feature with the rising arm: a long pyrimidine tract that is not a
microsatellite (with a distant branch point), an acceptor inside an Alu, a poly-A run (≥ 15 nt) at the acceptor,
a T-run of ≥ 20 nt deep in the intron, or one after the donor. (a) The pre-registered contrast: the rising arm's
change minus the control's, within strata of baseline band × coverage tertile, coarsened until each rising-arm
stratum holds ≥ 3 controls; filled = supported, hollow = not resolved. (b) The difference each contrast could
detect with 80% power at two-sided 0.05 (2.80 SD of its exon-permutation null, rebuilt with step 62's code and
seed; normal approximation; descriptive), against the whole rise of the rising arm (blue line). Light grey bars =
the three contrasts that were supported, dark grey = the two that were not resolved (poly-A and donor-side runs);
the shade marks the verdict, not the bar's length.

- {tline('RP4b', 'rising arm minus long pyrimidine tract')} {tline('RP4e', 'rising arm minus Alu acceptor')}
  {tline('RP4c', 'rising arm minus T-run deep in the intron')}
- {tline('RP4a', 'rising arm minus poly-A run')} {tline('RP4d', 'rising arm minus T-run after the donor')} p = the
  larger of the exon-permutation p ({fp(T('RP4a', 'p_class'))}, {fp(T('RP4d', 'p_class'))}) and the
  MMR-label-permutation p ({fp(T('RP4a', 'p_label'))}, {fp(T('RP4d', 'p_label'))}).
- These two are underpowered rather than negative: detectable differences about {pw['RP4a']['mdd_80']:.3f} for
  poly-A ({ctl_n('RP4a')}) and {pw['RP4d']['mdd_80']:.3f} for donor-side runs; if the controls did not move at all,
  the power would be about {pw['RP4a']['power_if_control_flat']:.2f} and {pw['RP4d']['power_if_control_flat']:.2f}.
  The controls' own means: poly-A {sg(own('polyA15').value)} (p = {fp(own('polyA15').p)}), donor-side
  {sg(own('donor').value)} (p = {fp(own('donor').p)}).
- 61 of the 92 rising-arm acceptors lie in an Alu. For the other 31, rising arm minus distant branch point without
  a run: {sg(na.estimate)}, p = {fp(na.p)}."""


# ============================================================ RP5 confounders
def RP5() -> None:
    bs = desc("control class: own mean delta against 0 (label_perm p)", "B sample (all)").value
    s_ = lambda item, test: sens(item, test)
    ex_rows = [("no covariates", T("RP2a") - bs, np.nan, "none"),
               ("eight exon covariates\n(pre-registered)", T("RP5a"), T("RP5a", "p"), "RP5a"),
               ("+ full splice-site score\nand acceptor Alu", s_("RP5a's full model (+ -20..+3 PWM, acceptor Alu)", "RP5a").estimate,
                s_("RP5a's full model (+ -20..+3 PWM, acceptor Alu)", "RP5a").p, "sensitivity"),
               ("+ longest T-run in the\ncoding sequence", s_("RP5a with cds_max_run added", "RP5a").estimate,
                s_("RP5a with cds_max_run added", "RP5a").p, "sensitivity")]
    tu_rows = [("no covariates", T("RP2b"), T("RP2b", "p"), "RP2b"),
               ("five tumour covariates\n(pre-registered)", T("RP5b"), T("RP5b", "p"), "RP5b"),
               ("without purity,\nall tumours", s_("RP5b without purity (all held-out tumours)", "RP5b").estimate,
                s_("RP5b without purity (all held-out tumours)", "RP5b").p, "sensitivity"),
               ("+ MKI67 expression", s_("RP5b with MKI67 log2(CPM + 1) percentile within cancer", "RP5b").estimate,
                s_("RP5b with MKI67 log2(CPM + 1) percentile within cancer", "RP5b").p, "sensitivity")]
    fig, axs = plt.subplots(1, 2, figsize=(DOUBLE * 0.86, 1.6), gridspec_kw={"wspace": 1.0})
    rows = []
    for ax, items, head, xl in [(axs[0], ex_rows, "Exon level", "Rising-arm coefficient (PSI)"),
                                (axs[1], tu_rows, "Tumour level", "MMRd coefficient on the\nrising-arm level (PSI)")]:
        for y, (lab_, v, p_, src) in enumerate(items):
            col = ACCENT if src.startswith("RP5") else (MUTED if src == "sensitivity" else INK2)
            ax.plot([v], [y], marker="o", markersize=MS + 0.8, color=col, linestyle="", zorder=3)
            rows.append({"level": head, "row": lab_.replace("\n", " "), "estimate": v, "p": p_, "source": src,
                         **({k: v_ for k, v_ in qcols(src).items() if k in ("q", "q_family", "verdict")}
                            if src in TESTS.index else {})})
        rows_axis(ax, [r[0] for r in items], top_pad=0.5, bottom_pad=0.5)
        zero(ax)
        ax.set_xlim(-0.0005, 0.011)
        ax.set_xticks([0, 0.005, 0.01])
        ax.set_xlabel(xl)
        heading(ax, head)
    rows.append({"row": "reference line", "estimate": 0.0})
    panel_letters([(axs[0], -1.0), (axs[1], -1.0)], y=1.04)
    save(fig, OUT, "RP5_confounders", data=pd.DataFrame(rows), caption=CAPTIONS["RP5"](bs))


def _cap_rp5(bs) -> str:
    return f"""**Do exon and tumour covariates explain the rise?**

Held-out cancers, junction reads. {DEFINE['MMR']}. (a) OLS of each exon's change in inclusion on a rising-arm
indicator against the B sample of background exons (92 rising-arm + 4,756 exons). The pre-registered model adds
eight exon covariates: baseline, length, coverage, host-gene change, poison status, the core splice-site score
(−3..+3, blind to the tract), GC and Alu exon. The full −20..+3 score and the acceptor-Alu flag measure the tract
itself and enter only a sensitivity model. (b) OLS of each held-out tumour's rising-arm level on MMR status with
cancer indicators. The pre-registered model adds the tumour's global inclusion axis, stroma, immune infiltrate,
retained-intron level and purity. Blue = the pre-registered models; dark grey = no covariates; light grey =
sensitivity models. No intervals were computed for these coefficients; their p values are below.

- {tline('RP5a', 'rising-arm coefficient, Freedman–Lane within band × tertile and MMR-label permutation, p = max')}
  {tline('RPC5a', 'falling arm')}
- {tline('RP5b', 'MMRd coefficient, Freedman–Lane with tumour residuals permuted within cancer')}
  {tline('RPC5b', 'falling-arm level')}
- No covariates: rising-arm mean minus background mean {sg(T('RP2a') - bs)}; {tline('RP2b', 'arm contrast')}"""


# ============================================================ RP6 specificity
LMS_NAMES = {"LMS": "leiomyosarcoma", "UPS": "pleomorphic sarcoma", "MFS": "myxofibrosarcoma",
             "DDLPS": "liposarcoma", "Cortical admixture": "cortical admixture", "Wnt-altered": "Wnt-altered"}


def RP6() -> None:
    from scipy.stats import mannwhitneyu
    S, _ = state()
    lab, z = S.lab, _z_within()
    h = lab.held_out.to_numpy() & np.isfinite(z)
    groups = [("MMR-\ndeficient", h & lab.mmr.to_numpy(), ACCENT),
              ("hyper-\nmutated,\nproficient", h & lab.hyper.to_numpy() & ~lab.mmr.to_numpy(), INK2),
              ("other\nproficient", h & ~lab.hyper.to_numpy() & ~lab.mmr.to_numpy(), MUTED)]
    u = mannwhitneyu(z[groups[0][1]], z[groups[1][1]], alternative="two-sided")
    same(u.statistic / (groups[0][1].sum() * groups[1][1].sum()) - 0.5, T("RP6a"), "RP6a")
    nul = pd.read_csv(CODE / "rp6c_null_contrasts.csv")
    nul["abs_s"] = nul.s.abs()
    s_msi = T("RP6c")
    allc = pd.concat([nul.assign(kind="other"), pd.DataFrame([{"cohort": "held out", "a": "MMRd", "b": "MMRp",
                                                               "abs_s": s_msi, "kind": "MSI"}])], ignore_index=True)
    allc = allc.sort_values("abs_s").reset_index(drop=True)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(DOUBLE, 2.1), gridspec_kw={"width_ratios": [1.1, 1.5], "wspace": 0.34})
    rows, local = [], []
    rng = np.random.default_rng(0)
    for k, (name, m, col) in enumerate(groups):
        v = z[m]
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        p5, p95 = np.percentile(v, [5, 95])
        ax.add_patch(plt.Rectangle((k - 0.22, q1), 0.44, q3 - q1, facecolor="white", edgecolor=col, linewidth=0.8,
                                   zorder=2))
        ax.plot([k, k], [p5, q1], color=col, linewidth=0.7)
        ax.plot([k, k], [q3, p95], color=col, linewidth=0.7)
        ax.plot([k - 0.22, k + 0.22], [med, med], color=col, linewidth=1.4, zorder=3)
        if k < 2:
            ax.scatter(k + rng.uniform(-0.33, 0.33, len(v)), v, s=3 if len(v) > 100 else 7, color=col,
                       alpha=0.45 if len(v) > 100 else 0.85, linewidths=0, zorder=1)
        rows.append({"panel": "a", "group": name.replace("\n", " "), "n_tumours": int(m.sum()), "median_z": med,
                     "q25": q1, "q75": q3, "p05": p5, "p95": p95})
        local += [{"panel": "a", "group": name.replace("\n", " "), "file_id": f, "project": pj, "z": zz}
                  for f, pj, zz in zip(lab.file_id[m], lab.project[m], v)]
    rows.append({"panel": "a", "group": "MMRd vs hypermutated MMRp", "auroc": T("RP6a") + 0.5, **qcols("RP6a")})
    zero(ax, horizontal=True)
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"{g[0]}\n({int(g[1].sum()):,})" for g in groups])
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(-4, 6)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylabel("Rising-arm level,\nz within cancer")
    # b. MSI ranked among the 57 non-MSI subtype contrasts
    is_msi = allc.kind.eq("MSI").to_numpy()
    over = allc.abs_s.to_numpy() > s_msi
    bx.scatter(np.flatnonzero(~is_msi), allc.abs_s[~is_msi], s=9, color=np.where(over[~is_msi], INK2, LIGHT),
               linewidths=0, zorder=2)
    k_msi = int(np.flatnonzero(is_msi)[0])
    bx.scatter([k_msi], [s_msi], s=26, color=ACCENT, linewidths=0, zorder=4)
    bx.annotate("MSI: MMRd vs MMRp", (k_msi, s_msi), xytext=(k_msi - 7, s_msi), textcoords="data", fontsize=FT,
                color=ACCENT, ha="right", va="center", arrowprops=dict(arrowstyle="-", linewidth=0.4, color=ACCENT,
                                                                      shrinkA=1, shrinkB=3))
    lead = dict(arrowstyle="-", linewidth=0.4, color=INK2, shrinkA=0, shrinkB=2)
    top = allc[over & ~is_msi].sort_values("abs_s", ascending=False)
    for yt, (ii, r) in zip([0.0160, 0.0146, 0.0124, 0.0108], top.iterrows()):
        nm = f"{r.cohort.split('_')[0]}: {LMS_NAMES.get(r.a, r.a)} vs {LMS_NAMES.get(r.b, r.b)}"
        bx.annotate(nm, (ii, r.abs_s), xytext=(42, yt), textcoords="data", fontsize=FT, ha="right", va="center",
                    color=INK, arrowprops=lead)
    for ii, r in allc.iterrows():
        rows.append({"panel": "b", "group": f"{r.cohort} {r.a} vs {r.b}", "rank": ii + 1, "abs_s": r.abs_s,
                     "is_msi": bool(r.kind == "MSI"), "above_msi": bool(r.abs_s > s_msi)})
    rows.append({"panel": "b", "group": "MSI held out", "abs_s": s_msi, **qcols("RP6c")})
    bx.set_xlim(-2, len(allc) + 1)
    bx.set_xticks([0, 19, 39, len(allc) - 1])
    bx.set_xticklabels(["1", "20", "40", str(len(allc))])
    bx.set_ylim(0, 0.0172)
    bx.set_xlabel("Subtype contrasts, ranked")
    bx.set_ylabel("Rising-arm shift between\ntwo classes, |PSI|")
    panel_letters([(ax, -0.3), (bx, -0.17)], y=1.04)
    save(fig, OUT, "RP6_specificity", data=pd.DataFrame(rows), caption=CAPTIONS["RP6"](groups, nul, s_msi),
         local_data=(LOCAL, pd.DataFrame(local)))


def _cap_rp6(groups, nul, s_msi) -> str:
    w = sens("RP6a within cancer", "RP6a")
    ms = sens("RP6b with MSIsensor MSI-H", "RP6b")
    top = nul[nul.abs_s > s_msi]
    n_over, n_sarc = len(top), int(top.cohort.str.startswith("SARC").sum())
    return f"""**Is it mismatch repair rather than mutation load, and is it unique to MSI?**

Held-out cancers, junction reads. {DEFINE['MMR']}; hypermutated MMRp = ≥ 500 single-nucleotide variants without MMR
deficiency. (a) Each tumour's rising-arm level (mean adjusted PSI over its rising-arm exons; {DEFINE['R_up']}) as a
z within its cancer against the cancer's non-hypermutated MMRp tumours; boxes = quartiles, line = median, whiskers =
5th–95th percentiles, dots = tumours (not drawn for the other MMRp tumours). (b) The shift of the rising-arm exons
between two classes (D-70's RV1 in junction reads: 20 tumours drawn from each class, mean over rising-arm exons of
the difference in mean adjusted PSI, mean over 50 draws), for MSI (20 held-out MMRd tumours, each with a random
MMRp tumour of its cancer) and for each of the 57 non-MSI subtype contrasts of the subtype task (two classes of one
cohort with ≥ 20 tumours each, MSI, POLE and MMRd tumours removed), ranked; blue = MSI, grey = the 57 other
contrasts, dark grey and labelled = those above MSI. The null contrasts are not matched by cancer, while the MSI draws pair tumours within cancer.

- {tline('RP6a', 'AUROC − 0.5, MMRd against hypermutated MMRp, Mann–Whitney', 3)} Within the cancers that have both:
  {sg(w.estimate, 3)}, p = {fp(w.p)}.
- {tline('RP6b', 'arm contrast with MANTIS MSI-H as the label')} With MSIsensor: {sg(ms.estimate)}, p = {fp(ms.p)}.
- {tline('RP6c', 'MSI ranked against the 57 null contrasts (p floor 1/58)')} {WORDS[n_over].capitalize()} null
  contrasts move the arm more than MSI: {WORDS[n_sarc]} in SARC, all involving leiomyosarcoma, and
  {WORDS[n_over - n_sarc]} in PCPG (cortical admixture against Wnt-altered)."""


# ============================================================ RP7 within patient
def pair_levels() -> pd.DataFrame:
    """Each tumour-normal pair's rising-arm level in the tumour and in its own normal tissue, each adjusted as
    step 62 adjusts their difference (deviation from the cancer's median, band adjustment over B-sample
    exons) on the exons measured in both samples, so that tumour minus normal is exactly the committed pair
    arm (checked)."""
    def go():
        S, _ = state()
        spec = P.SPECS["primary"]
        arm = pd.read_parquet(P.PROOF / "tests" / "pairs.parquet")
        keep = sorted(set(arm.project[arm.mmr]))
        bflag, upflag = S.ex.b_sample.to_numpy(bool), S.ex.cls.eq("R_up").to_numpy()
        out = []
        for proj in keep:
            r = S.res["primary"][proj]
            U, m, band = np.asarray(r["U"], bool), np.asarray(r["m"], float), np.asarray(r["band"], int)
            _, ijc, sjc, ratio = P.load_counts_all(proj, spec.view)
            for q in S.pairs[S.pairs.project.eq(proj)].itertuples():
                idx = [int(q.tumour_row), int(q.normal_row)]
                I, Sj, R = (np.asarray(a[idx], np.float64) for a in (ijc, sjc, ratio))
                Pp = np.where(I + Sj >= spec.min_reads, P.psi_from_counts(I, Sj, R), np.nan)
                J = U & np.isfinite(Pp).all(0)
                adj, _ = P.band_adjust(np.where(J, Pp - m, np.nan), band, bflag & U)
                lv = P.class_arm(adj, upflag & U)
                out.append({"case_id": q.case_id, "project": proj, "tumour": float(lv[0]), "normal": float(lv[1])})
        return pd.DataFrame(out)
    lv = cached("pair_levels", go)
    arm = pd.read_parquet(P.PROOF / "tests" / "pairs.parquet")
    j = arm.merge(lv, on=["case_id", "project"], how="inner")
    j = j[np.isfinite(j.pair_up.to_numpy(float))]
    diff = np.abs((j.tumour - j.normal).to_numpy() - j.pair_up.to_numpy(float))
    if not (np.nanmax(diff) < 1e-5):
        raise SystemExit(f"pair levels: tumour − normal differs from the committed pair arm by {np.nanmax(diff):.2e}")
    return j


def RP7() -> None:
    j = pair_levels()
    fig, axs = plt.subplots(1, 2, figsize=(DOUBLE * 0.62, 2.15), sharey=True, gridspec_kw={"wspace": 0.12})
    rows, local = [], []
    for ax, st, head in [(axs[0], True, "MMR-deficient patients"), (axs[1], False, "MMR-proficient patients")]:
        g = j[j.mmr.eq(st)]
        for r in g.itertuples():
            held = r.project not in P.TRAIN
            if st:
                ax.plot([0, 1], [r.normal, r.tumour], color=INK2, linewidth=0.8, zorder=2,
                        linestyle=(0, (2, 1.5)) if held else "-")
                ax.plot([0], [r.normal], marker="o", markersize=3.0, color=NORMAL, linestyle="", zorder=3)
                ax.plot([1], [r.tumour], marker="o", markersize=3.0, color=TUMOUR, linestyle="", zorder=3)
            else:
                ax.plot([0, 1], [r.normal, r.tumour], color=LIGHT, linewidth=0.4, alpha=0.6, zorder=1)
            local.append({"case_id": r.case_id, "project": r.project, "mmr": "MMRd" if st else "MMRp",
                          "normal": r.normal, "tumour": r.tumour})
        mn, mt = float(g.normal.mean()), float(g.tumour.mean())
        ax.plot([0, 1], [mn, mt], color=INK, linewidth=1.6, zorder=4)
        if not st:
            ax.plot([0], [mn], marker="o", markersize=3.6, color=NORMAL, linestyle="", zorder=5)
            ax.plot([1], [mt], marker="o", markersize=3.6, color=TUMOUR, linestyle="", zorder=5)
        rows.append({"group": "MMRd" if st else "MMRp", "n_pairs": len(g), "mean_normal": mn, "mean_tumour": mt,
                     "mean_change": mt - mn})
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["own normal", "tumour"])
        ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.25, 1.25)
        heading(ax, f"{head} ({len(g)})")
        zero(ax, horizontal=True)
    for c in sorted(set(j.project)):
        for st in (True, False):
            g = j[j.project.eq(c) & j.mmr.eq(st)]
            ok = len(g) >= 5
            rows.append({"group": f"{c} {'MMRd' if st else 'MMRp'}", "training_cancer": c in P.TRAIN, "n_pairs": len(g),
                         "mean_change": float((g.tumour - g.normal).mean()) if ok else np.nan,
                         "note": "" if ok else "mean blanked (< 5 pairs)"})
    rows.append({"group": "all", **qcols("RP7")})
    rows.append({"group": "reference line", "mean_change": 0.0})
    axs[0].set_ylabel("Rising-arm level vs the\ncancer's median (PSI)")
    axs[0].legend(handles=[Line2D([], [], color=INK2, linewidth=0.8, label="training cancer"),
                           Line2D([], [], color=INK2, linewidth=0.8, linestyle=(0, (2, 1.5)), label="held-out cancer"),
                           Line2D([], [], color=INK, linewidth=1.6, label="mean")],
                  loc="lower left", fontsize=FT, handlelength=1.8, borderaxespad=0.1, labelspacing=0.3)
    panel_letters([(axs[0], -0.5), (axs[1], -0.08)], y=1.1)
    save(fig, OUT, "RP7_within_patient", data=pd.DataFrame(rows), caption=CAPTIONS["RP7"](j),
         local_data=(LOCAL, pd.DataFrame(local)))


def _cap_rp7(j) -> str:
    md = j[j.mmr]
    n_tr = int(md.project.isin(P.TRAIN).sum())
    by = md.project.value_counts()
    return f"""**Is the rise there within the same patient?**

Tumour–normal pairs (the tumour and a Solid Tissue Normal file of the same patient) of the cancers with at least
one MMR-deficient pair: {', '.join(sorted(set(j.project)))}. {DEFINE['MMR']}. Each line joins a patient's own normal
tissue to the tumour: the rising-arm level ({DEFINE['R_up']}) of each sample, its mean PSI deviation from the
cancer's median over the rising-arm exons measured in both samples, minus its mean deviation over background exons
of the same baseline band. The tumour minus the normal is the pair arm that the test uses. (a) The
{int(j.mmr.sum())} MMR-deficient patients; {n_tr} are in training cancers (COAD {by.get('COAD', 0)}, STAD
{by.get('STAD', 0)}, UCEC {by.get('UCEC', 0)}), so the test is not held out (family RPP). (b) The
{int((~j.mmr).sum())} MMR-proficient patients of the same cancers. Orange = own normal tissue, blue = tumour; black
line = mean.

- {tline('RP7', 'mean over MMRd pairs of the pair arm minus the MMRp mean of the same cancer; MMR labels permuted among each cancer’s pairs (10,000)')}
  {tline('RPC7', 'the falling arm')}
- Of the two held-out MMR-deficient pairs, the BRCA pair lies well above its cancer's MMR-proficient mean, while
  the PRAD pair is about level with it (single-pair values are kept locally; committed cells under 5 pairs are
  blanked, A-24)."""


# ============================================================ RP8 dose
def RP8() -> None:
    S, H = state()
    lab = S.lab
    arm = S.tv.arm_up.to_numpy(float)
    man = lab.mantis.to_numpy(float)
    nhp = lab.n_hp.to_numpy(float)
    m_a = lab.mmr.to_numpy() & lab.held_out.to_numpy() & np.isfinite(arm) & np.isfinite(man)
    m_b = lab.mmr.to_numpy() & np.isfinite(arm) & np.isfinite(nhp)
    up = sel_of("R_up")
    d, L = delta0()[up], H.tab.acc_T.to_numpy(float)[up]
    fig, axs = plt.subplots(1, 3, figsize=(DOUBLE, 1.95), gridspec_kw={"wspace": 0.5})
    ax, bx, cx = axs
    ax.scatter(man[m_a], arm[m_a], s=9, color=ACCENT, edgecolors="none")
    ax.set_xlabel("MANTIS instability score")
    ax.set_ylabel("Rising-arm level vs the\ncancer's median (PSI)")
    heading(ax, f"Held-out MMRd tumours ({int(m_a.sum())})")
    tr = lab.project.isin(P.TRAIN).to_numpy()
    for sel_, col in [(m_b & tr, LIGHT), (m_b & ~tr, ACCENT)]:
        bx.scatter(nhp[sel_], arm[sel_], s=8, color=col, edgecolors="none")
    bx.set_xscale("log")
    bx.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    bx.set_xlabel("Homopolymer deletions (MC3)")
    bx.set_ylabel("Rising-arm level vs the\ncancer's median (PSI)")
    heading(bx, f"All MMRd tumours ({int(m_b.sum())})")
    bx.legend(handles=[Line2D([], [], marker="o", color=ACCENT, linestyle="", markersize=3.2, label="held out"),
                       Line2D([], [], marker="o", color=LIGHT, linestyle="", markersize=3.2, label="training")],
              loc="lower right", fontsize=FT, handletextpad=0.1, borderaxespad=0.1, labelspacing=0.25)
    cx.scatter(L + np.random.default_rng(0).uniform(-0.25, 0.25, len(L)), d, s=8, color=ACCENT, edgecolors="none")
    cx.set_xlabel("Acceptor T-run length (nt)")
    cx.set_ylabel(YCHG)
    heading(cx, f"Rising-arm exons ({len(L)})")
    for a in axs:
        zero(a, horizontal=True)
    panel_letters([(ax, -0.42), (bx, -0.2), (cx, -0.42)], y=1.1)
    rows = [{"panel": "a", "n_tumours": int(m_a.sum()), "median_mantis": float(np.median(man[m_a])),
             "median_arm": float(np.median(arm[m_a])), **qcols("RP8a")},
            {"panel": "b", "n_tumours": int(m_b.sum()), "median_deletions": float(np.median(nhp[m_b])),
             "median_arm": float(np.median(arm[m_b])), **qcols("RP8b")},
            {"panel": "a, b, c", "reference line": 0.0},
            {"panel": "c", **qcols("RP8c")}]
    rows += [{"panel": "c", "run_length": a_, "change": b_} for a_, b_ in zip(L, d)]
    local = pd.DataFrame({"file_id": lab.file_id[m_b], "project": lab.project[m_b], "held_out": ~tr[m_b],
                          "mantis": man[m_b], "deletions": nhp[m_b], "arm_up": arm[m_b]})
    save(fig, OUT, "RP8_dose", data=pd.DataFrame(rows), caption=CAPTIONS["RP8"](int(m_a.sum()), int(m_b.sum()), L),
         local_data=(LOCAL, local))


def _cap_rp8(na, nb, L) -> str:
    return f"""**Does the rise follow dose?**

{DEFINE['MMR']}; {DEFINE['R_up']}. A tumour's rising-arm level = its mean adjusted PSI deviation from its cancer's
median over the rising-arm exons. (a) The {na} held-out MMRd tumours: rising-arm level against the MANTIS score,
the exome instability score. (b) All {nb} MMRd tumours, training cancers included: rising-arm level against the
number of MC3 homopolymer deletions (log scale). (c) The rising-arm exons: change in inclusion against acceptor
T-run length ({int(L.min())}–{int(L.max())} nt; jittered). The tests ask for a graded relation, which D-70's
percentile read-out had not found (RV6); the first two are partial Spearman correlations given cancer.

- {tline('RP8a', 'partial Spearman of the arm with MANTIS given cancer, MANTIS permuted within cancer', 3)}
  {tline('RPC8', 'the falling arm (predicted < 0)', 3)}
- {tline('RP8b', 'the same with log(1 + deletion count), all 200 MMRd tumours (not held out)', 3)}
- {tline('RP8c', 'Spearman of the change with run length over the rising arm, run length permuted within baseline band', 3)}"""


# ============================================================ RPM mechanism
def RPM() -> None:
    ci = {c: class_interval(c) for c in ["R_up", "shield", "long, antisense Alu", "long, other"]}
    fig, axs = plt.subplots(1, 2, figsize=(DOUBLE * 0.8, 1.45), gridspec_kw={"wspace": 0.95})
    rows = []
    for ax, pair, head in [(axs[0], ["R_up", "shield"], "Branch-point geometry"),
                           (axs[1], ["long, antisense Alu", "long, other"], "hnRNP C")]:
        for y, c in enumerate(pair):
            dot_ci(ax, y, ci[c], ACCENT if c == "R_up" else INK2)
            rows.append({"class": c, "label": NAME[c].replace("\n", " "), **ci_row(ci[c])})
        hs = sesoi_line(ax, 2)
        zero(ax, span=(-0.5, 1.5))
        rows_axis(ax, [f"{NAME[c]} ({ci[c]['n_exons']})" for c in pair], top_pad=0.6, bottom_pad=0.6)
        ax.set_xlim(-0.0075, 0.0145)
        ax.set_xticks([-0.005, 0, 0.005, 0.01])
        ax.set_xlabel(XCHG.replace(", MMRd", ",\nMMRd"))
        heading(ax, head)
    axs[1].legend(handles=[hs], loc="lower right", bbox_to_anchor=(1.0, 1.12), fontsize=FT, handlelength=1.6,
                  borderaxespad=0.1)
    for i in ["RPM1", "RPM2", "RPM3"]:
        rows.append({"class": i, **qcols(i)})
    rows.append({"class": "reference lines", "mean_change": 0.0, "sesoi": P.SESOI})
    panel_letters([(axs[0], -1.05), (axs[1], -1.05)], y=1.08)
    save(fig, OUT, "RPM_mechanism", data=pd.DataFrame(rows), caption=CAPTIONS["RPM"](ci))


def n_eff(i: str) -> str:
    """'92 vs 64 exons (effective controls 45)' -> '92 against 64 exons, 45 effective controls'."""
    import re
    m = re.match(r"(\d+) vs (\d+) exons \(effective controls (\d+)\)", TESTS.at[i, "n"])
    return f"{m.group(1)} against {m.group(2)} exons, {m.group(3)} effective controls"


def _cap_rpm(ci) -> str:
    return f"""**Which mechanism: branch-point geometry or hnRNP C?**

Held-out cancers, junction reads; outside the claim's verdict (family RPM). {DEFINE['change'][0].upper() +
DEFINE['change'][1:]}. Mean change with a descriptive 95% interval; dashed = the smallest effect of interest.
(a) {DEFINE['R_up'][0].upper() + DEFINE['R_up'][1:]}; T-run upstream of the branch point = exons without an acceptor
T-run of ≥ 10 nt whose intron holds a T-run of ≥ 20 nt upstream of the branch point (shield exons). Shortening that
run cannot bring the branch point closer to the AG, so the geometry reading predicts no rise there and a larger
rise in the rising arm. (b) Eligible runs of ≥ 20 nt (rising arm, or with a near or no branch point) inside an
antisense Alu and outside one. hnRNP C binds long U-runs and represses Alu-derived exons, so the hnRNP C reading
predicts that shield exons rise too and that runs inside an antisense Alu rise more. Numbers = exons measured.
The tests compare classes within strata of baseline band × coverage tertile (and run-length bins for (b)), keeping
the exons that have a partner there ({n_eff('RPM1')} for (a); {n_eff('RPM3')} for (b)), so they differ
from the raw difference of the dots.

- {tline('RPM1', 'rising arm minus T-run upstream of the branch point')}
- {tline('RPM2', 'T-run upstream of the branch point against 0, MMR-label permutation')}
- {tline('RPM3', 'runs in an antisense Alu minus other long runs')}
- None resolved: the data do not separate the two readings."""


# ============================================================ the table of all tests (section 5 of the write-up)
PLAIN = {
    "RP1": ("1", "Held out, the rising arm separates MMRd tumours better than matched exon sets (AUROC minus their median)"),
    "RPC1": ("1", "The same for the falling arm (predicted below)"),
    "RPD1a": ("1", "Training: the subtype networks push towards MSI when rising-arm inputs rise (share minus 0.5)"),
    "RPD1b": ("1", "The same for falling-arm inputs (predicted below)"),
    "RP2a": ("2", "Mean change of the rising-arm exons (primary test)"),
    "RPC2a": ("2", "Mean change of the falling-arm exons"),
    "RP2b": ("2", "Per tumour: rising-arm level, MMRd minus MMRp"),
    "RPC2b": ("2", "Per tumour: falling-arm level, MMRd minus MMRp"),
    "RP2c": ("2", "Share of rising-arm exons that rise, minus that of matched background exons"),
    "RPD2a": ("2", "Training: mean change of the rising-arm exons"),
    "RPD2b": ("2", "Training: mean change of the falling-arm exons"),
    "RP3a": ("3", "Rising arm minus exons with a distant branch point and no T-run"),
    "RPD3a": ("3", "Training: the same"),
    "RP3b": ("3", "Rising arm minus long runs with a near predicted branch point (the 36-nt cut)"),
    "RP3d": ("3", "Mapped branch point ≥ 36 nt minus < 36 nt"),
    "RP3c": ("3", "Change against the gap to the predicted branch point, at fixed run length and offset (ranks)"),
    "RP3e": ("3", "Change against the gap to the mapped branch point, at fixed run length and offset (ranks)"),
    "RPD3b": ("3", "Training: the same as RP3e"),
    "RP4b": ("4", "Rising arm minus long pyrimidine tracts"),
    "RP4e": ("4", "Rising arm minus Alu acceptors without a T-run"),
    "RP4c": ("4", "Rising arm minus T-runs deep in the intron"),
    "RP4a": ("4", "Rising arm minus poly-A runs at the acceptor"),
    "RP4d": ("4", "Rising arm minus T-runs after the donor"),
    "RP5a": ("5", "Rising-arm coefficient with eight exon covariates"),
    "RPC5a": ("5", "Falling-arm coefficient, the same model"),
    "RP5b": ("5", "MMRd coefficient on the rising-arm level, five tumour covariates"),
    "RPC5b": ("5", "The same on the falling-arm level"),
    "RP6a": ("6", "Rising-arm level, MMRd against hypermutated MMRp (AUROC minus 0.5)"),
    "RPC6": ("6", "The same for the falling arm"),
    "RP6b": ("6", "Rising-arm level, MANTIS MSI-high against stable"),
    "RP6c": ("6", "MSI shift ranked against 57 other subtype contrasts"),
    "RP7": ("7", "Tumour minus own normal tissue, MMRd minus MMRp pairs (rising arm; not held out)"),
    "RPC7": ("7", "The same for the falling arm"),
    "RP8a": ("8", "Rising-arm level against MANTIS score, held-out MMRd (partial ρ)"),
    "RPC8": ("8", "Falling-arm level against MANTIS score"),
    "RP8b": ("8", "Rising-arm level against deletion count, all MMRd (partial ρ)"),
    "RP8c": ("8", "Change against run length over the rising arm (ρ)"),
    "RPM1": ("M", "Rising arm minus T-runs upstream of the branch point"),
    "RPM2": ("M", "T-runs upstream of the branch point, against 0"),
    "RPM3": ("M", "Runs in an antisense Alu minus other long runs"),
}
RHO = {"RP3c", "RP3e", "RPD3b", "RP8a", "RP8b", "RP8c", "RPC8"}
SHARE = {"RP1", "RPC1", "RP2c", "RP6a", "RPC6", "RPD1a", "RPD1b"}


def tests_table() -> str:
    assert set(PLAIN) == set(TESTS.index), set(PLAIN) ^ set(TESTS.index)
    lines = ["| Step | Test | In plain words | Effect | p | q | Verdict |", "|---|---|---|---|---|---|---|"]
    for i, (s, words) in PLAIN.items():
        nd = 3 if (i in RHO or i in SHARE) else 4
        e = sg(T(i), nd)
        if i in ("RP2a", "RPC2a", "RPD2a", "RPD2b"):
            e += f" ({num(T(i, 'ci_low'))} to {num(T(i, 'ci_high'))})"
        lines.append(f"| {s} | {i} ({TESTS.at[i, 'family']}) | {words} | {e} | {fp(T(i, 'p'))} | {fp(T(i, 'q'))} | "
                     f"{vd(i)} |")
    return "\n".join(lines)


def write_tests_table() -> None:
    import re
    txt = WRITEUP.read_text()
    a, b = "<!-- tests-table: written by scripts/63's 'legends' step from tests.csv -->", "<!-- /tests-table -->"
    if a not in txt:
        raise SystemExit("write-up: tests-table markers missing")
    new = f"{a}\n{tests_table()}\n{b}"
    WRITEUP.write_text(re.sub(re.escape(a) + r".*?" + re.escape(b), lambda _: new, txt, flags=re.S))
    log("tests table written")


# ============================================================ legends
TITLES = {
    "RP0": "In junction reads the rising arm rises in held-out cancers, but the claim is only partly supported.",
    "RP1": "Held out, the rising arm separates MMR-deficient tumours better than matched exon sets, but removing it "
           "from the fitted model costs no detectable accuracy.",
    "RP2": "The rising arm rises in held-out junction reads, and a few of its exons carry most of the rise.",
    "RP3": "The rise needs the T-run, a step at 36 nt is not supported, and the link to the branch-point gap is "
           "seen in held-out cancers only.",
    "RP4": "The rise exceeds that of three other long tracts; the poly-A and donor-side contrasts were "
           "underpowered.",
    "RP5": "The rise remains after exon and tumour covariates.",
    "RP6": "The rise marks mismatch repair rather than mutation load, but four other subtype contrasts (three in "
           "sarcoma) move it more than MSI.",
    "RP7": "MMR-deficient tumours rise above their own normal tissue, mostly in training cancers.",
    "RP8": "A dose relation is not shown.",
    "RPM": "The data do not separate branch-point geometry from hnRNP C.",
}


def legend(key: str, text: str) -> str:
    return _legend(key, TITLES[key], text)


CAPTIONS = {k: (lambda k, f: (lambda *a: tidy(legend(k, f(*a)))))(k, f) for k, f in
            {"RP0": _cap_rp0, "RP1": _cap_rp1, "RP2": _cap_rp2, "RP3": _cap_rp3, "RP4": _cap_rp4, "RP5": _cap_rp5,
             "RP6": _cap_rp6, "RP7": _cap_rp7, "RP8": _cap_rp8, "RPM": _cap_rpm}.items()}


def legends() -> None:
    """The table of all tests, then each figure's legend (without its bullet list of tests) under the figure."""
    if WRITEUP.exists():
        write_tests_table()
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
