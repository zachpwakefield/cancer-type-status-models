"""Step 61 (D-76): what does including a rising-arm exon do to its transcript?

Rising-arm (R_up) exons are mostly skipped in proficient tumours (D-75), so MSI raises the use of
exons that are normally near-silent. Each cassette exon is placed in its gene's canonical
protein-coding transcript (GENCODE v45) and, when inclusion inserts it into the coding sequence,
translated: an in-frame peptide insertion, or a premature stop that nonsense-mediated decay should
(ptc_nmd, a "poison" exon) or should not (ptc_escape) see. CQ1 checks the read-out on background
exons (poison exons cost mRNA across tumours); CQ2-CQ4 ask the same of the rising arm and of MSI's
rise. config/followup.yaml: msi_rising_consequence.

Inputs:  GENCODE v45 GTF and GRCh38 FASTA (as step 42), the pan-cancer SE and GEX matrices, D-75's
         per-exon table (data.nosync/ctsm_v1/msi/magnitude_exons.parquet), step 54's tumours
Outputs: results/subtype/msi_cis/code/consequence/ (aggregates and per-exon tables);
         data.nosync/ctsm_v1/msi/consequence_{annotation,coupling}.parquet (local)
Usage: python scripts/61_msi_rising_consequence.py [--only annotate|compute|report|posthoc]
"""
from __future__ import annotations

import gzip
import importlib.util
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from _common import log


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m60 = _load("msi60", "60_msi_magnitude.py")
m51 = m60.m51
m50, m42 = m51.m50, m51.m42
LOCAL, PAN, OUT, TRAIN = m60.LOCAL, m60.PAN, m51.OUT / "consequence", m51.TRAIN
STOPS = {"TAA", "TAG", "TGA"}
NMD_NT, MIN_TUM, MIN_HELD, BOOT = 55, 20, 20, 2000
CANDIDATE = (0.10, 0.05)                   # exploratory: MMR-deficient median >= 0.10, proficient and normal <= 0.05
CLASSES = ["R_up", "R_down", "long, other", "B"]
CONSEQUENCES = ["canonical", "in_frame", "ptc_nmd", "ptc_escape", "no_stop", "other"]


# ---------------------------------------------------------------- annotation (label-blind; unit-tested)

def parse_feature(fid: str) -> dict:
    """'ENSG...#chrN:es-ee#strand;ls-le;rs-re' (starts 0-based, ends inclusive) -> fields."""
    gene, loc, rest = fid.split("#")
    chrom, span = loc.split(":")
    es, ee = (int(v) for v in span.split("-"))
    strand, left, right = rest.split(";")
    ls, le = (int(v) for v in left.split("-"))
    rs, re_ = (int(v) for v in right.split("-"))
    return {"gene": gene.split(".")[0], "chrom": chrom, "strand": strand, "es": es, "ee": ee,
            "ls": ls, "le": le, "rs": rs, "re": re_}


def canonical_models(gtf: Path) -> dict:
    """Ensembl_canonical protein_coding transcript per gene: exons in transcript order (1-based,
    inclusive), the CDS part of each exon (or None) and the phase of the first CDS base."""
    ex, cds, info = {}, {}, {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] not in ("exon", "CDS") or 'tag "Ensembl_canonical"' not in f[8] \
                    or 'transcript_type "protein_coding"' not in f[8]:
                continue
            g = re.search(r'gene_id "([^".]+)', f[8]).group(1)
            k = int(re.search(r"exon_number (\d+)", f[8]).group(1))
            info[g] = (f[0], f[6])
            if f[2] == "exon":
                ex.setdefault(g, {})[k] = (int(f[3]), int(f[4]))
            else:
                cds.setdefault(g, {})[k] = (int(f[3]), int(f[4]), int(f[7]) if f[7] != "." else 0)
    out = {}
    for g, e in ex.items():
        ks = sorted(e)
        c = cds.get(g, {})
        first = min(c) if c else None
        out[g] = {"chrom": info[g][0], "strand": info[g][1], "exons": [e[k] for k in ks],
                  "cds": [c[k][:2] if k in c else None for k in ks], "phase": c[first][2] if c else 0}
    return out


def place(f: dict, m: dict) -> tuple[str, int]:
    """('canonical', i) if the exon is an internal all-coding exon of the canonical transcript;
    ('inserted', i) if the transcript joins its flanks directly (exon i to i + 1) inside the CDS;
    else ('other' or 'utr', -1)."""
    if m is None or m["chrom"] != f["chrom"] or m["strand"] != f["strand"]:
        return "other", -1
    ex, cds, n = m["exons"], m["cds"], len(m["exons"])
    for i in range(1, n - 1):
        if ex[i] == (f["es"] + 1, f["ee"]) and cds[i] == ex[i]:
            return "canonical", i
    for i in range(n - 1):
        a, b = ex[i], ex[i + 1]
        if f["strand"] == "+":
            hit = a[1] == f["le"] and b[0] == f["rs"] + 1
            in_cds = cds[i] is not None and cds[i][1] == a[1] and cds[i + 1] is not None and cds[i + 1][0] == b[0]
        else:
            hit = a[0] == f["rs"] + 1 and b[1] == f["le"]
            in_cds = cds[i] is not None and cds[i][0] == a[0] and cds[i + 1] is not None and cds[i + 1][1] == b[1]
        if hit:
            return ("inserted", i) if in_cds else ("utr", i)
    return "other", -1


def consequence(tail: str, cassette: str, downstream: list[str]) -> dict:
    """Translate tail (the last 0-2 coding bases before the junction) + cassette + downstream exons
    in the canonical frame. Returns the class, the stop's distance upstream of the last exon-exon
    junction (nt) and the inserted peptide length for in-frame exons."""
    p, L = len(tail), len(cassette)
    u = tail + cassette + "".join(downstream)
    stop = next((j for j in range(0, len(u) - 2, 3) if u[j:j + 3] in STOPS), None)
    last_junction = p + L + sum(len(d) for d in downstream[:-1])
    if L % 3 == 0 and (stop is None or stop >= p + L):
        return {"consequence": "in_frame", "dist": np.nan, "peptide": L // 3}
    if stop is None:
        return {"consequence": "no_stop", "dist": np.nan, "peptide": np.nan}
    d = last_junction - stop
    return {"consequence": "ptc_nmd" if d > NMD_NT else "ptc_escape", "dist": float(d), "peptide": np.nan}


def _seq(chrom_seq: str, a: int, b: int, strand: str) -> str:
    """1-based inclusive interval in transcript orientation."""
    s = chrom_seq[a - 1:b]
    return s[::-1].translate(m42.COMP) if strand == "-" else s


def annotate() -> None:
    from ctsm.config import load_config
    cls = m60.exon_classes()
    cls = cls[cls.mclass.ne("") & cls.struct].copy()
    models = canonical_models(load_config().gencode_gtf)
    rows = []
    for idx, r in cls.iterrows():
        f = parse_feature(r.feature_id)
        where, i = place(f, models.get(f["gene"]))
        rows.append({"feature_index": idx, "gene_id": f["gene"], **f, "where": where, "i": i,
                     "exon_len": f["ee"] - f["es"]})
    a = pd.DataFrame(rows)
    a["consequence"] = np.where(a["where"].eq("canonical"), "canonical", "other")
    a["dist"], a["peptide"] = np.nan, np.nan
    ins = a[a["where"].eq("inserted")]
    log(f"placed: {a['where'].value_counts().to_dict()}")
    for chrom, cs in m42.read_genome(set(ins.chrom)):
        for k, r in ins[ins.chrom.eq(chrom)].iterrows():
            m = models[r.gene_id]
            ex, cds, i, st = m["exons"], m["cds"], int(r.i), r.strand
            coding = sum(c[1] - c[0] + 1 for c in cds[:i + 1] if c is not None) - m["phase"]
            p = coding % 3
            up = _seq(cs, *cds[i], st)
            tail = up[len(up) - p:] if p else ""
            res = consequence(tail, _seq(cs, r.es + 1, r.ee, st), [_seq(cs, *e, st) for e in ex[i + 1:]])
            a.loc[k, ["consequence", "dist", "peptide"]] = [res["consequence"], res["dist"], res["peptide"]]
    a = a.join(cls[["mclass", "symbol", "acc_T", "d_bp"]], on="feature_index")
    a.to_parquet(LOCAL / "consequence_annotation.parquet", index=False)
    log("consequence by class:\n" + pd.crosstab(a.mclass, a.consequence).to_string())


# ---------------------------------------------------------------- compute (per-sample data stays local)

def compute() -> None:
    from ctsm.config import load_config
    from ctsm.matrix_io import load_view
    cfg = load_config()
    a = pd.read_parquet(LOCAL / "consequence_annotation.parquet")
    keep = a[a.consequence.isin(["in_frame", "ptc_nmd", "ptc_escape"]) |
             (a.consequence.eq("canonical") & a.mclass.ne("B"))]        # canonical B exons enter no test
    t = pd.read_parquet(LOCAL / "models_tumours.parquet")[["file_id", "project", "mmr"]]
    tab = pd.read_parquet(m51.DATA / "tasks" / "cancer_type__PANCAN.parquet")
    normals = tab[tab.role.eq("normal")][["file_id", "project"]]
    struct = m60.exon_classes().struct.to_numpy()
    rows, gene_rows, norm_rows = [], [], []
    for proj, s in t.groupby("project"):
        X = load_view(PAN, "se", s.file_id)[0].get()
        G, F = load_view(cfg.matrix_dir, "gex", s.file_id)
        G = np.asarray(G.get(), dtype=np.float64)
        gcol = pd.Series(np.arange(len(F)), index=F.gene_id.to_numpy())
        ok = np.isfinite(X).mean(axis=0) >= m60.MIN_FINITE
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            gmean = np.nanmean(X[:, ok & struct], axis=1)                 # D-24's global axis
        k = keep[ok[keep.feature_index.to_numpy()] & keep.gene_id.isin(gcol.index)]
        y = s.mmr.to_numpy()
        for r in k.itertuples():
            x, e = X[:, r.feature_index], G[:, gcol[r.gene_id]]
            if np.isfinite(x).sum() < MIN_TUM:
                continue
            rho, _, n = m42.partial_spearman(x, e, [gmean])
            rows.append({"project": proj, "feature_index": r.feature_index, "rho": rho, "n": n})
        genes = k.gene_id.unique()
        if y.any() and (~y).any():
            gi = gcol[genes].to_numpy()
            mm, pp = G[y][:, gi], G[~y][:, gi]
            if proj in TRAIN:
                gene_rows.append(pd.DataFrame({"project": proj, "gene_id": genes,
                                               "log2fc": mm.mean(0) - pp.mean(0), "n_mmr": int(y.sum())}))
            else:
                dev = mm - pp.mean(0)
                gene_rows.append(pd.DataFrame({"project": proj, "gene_id": np.repeat(genes[None], len(mm), 0).ravel(),
                                               "log2fc": dev.ravel(), "n_mmr": 1}))
        nm = normals[normals.project.eq(proj)]
        if len(nm) >= 5:
            Xn = load_view(PAN, "se", nm.file_id)[0].get()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                med = np.nanmedian(Xn[:, k.feature_index.to_numpy()], axis=0)
            norm_rows.append(pd.DataFrame({"project": proj, "feature_index": k.feature_index.to_numpy(),
                                           "normal_median": med, "n_normal": len(nm)}))
        log(f"{proj}: {len(k):,} annotated exons, {int(y.sum())} MMR-deficient, {len(nm)} normals")
    pd.DataFrame(rows).to_parquet(LOCAL / "consequence_coupling.parquet", index=False)
    pd.concat(gene_rows).to_parquet(LOCAL / "consequence_genes.parquet", index=False)
    pd.concat(norm_rows).to_parquet(LOCAL / "consequence_normals.parquet", index=False)


# ---------------------------------------------------------------- report

def _mwu(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    from scipy.stats import mannwhitneyu
    res = mannwhitneyu(a, b, alternative="two-sided")
    return float(res.statistic / (len(a) * len(b))), float(res.pvalue)


def tests(a, cp, ex, gn) -> pd.DataFrame:
    from scipy.stats import spearmanr
    pred = []
    per = cp.groupby("feature_index").rho.mean()
    for pid, c_ in [("CQ1", "B"), ("CQ2", "R_up")]:
        s = a[a.mclass.eq(c_)].assign(rho=lambda d: d.feature_index.map(per)).dropna(subset=["rho"])
        nmd, inf = s[s.consequence.eq("ptc_nmd")].rho, s[s.consequence.eq("in_frame")].rho
        if min(len(nmd), len(inf)) < 5:
            pred.append({"id": pid, "what": f"{c_}: not testable", "n": f"ptc_nmd {len(nmd)}, in_frame {len(inf)}"})
            continue
        _, p = _mwu(nmd.to_numpy(), inf.to_numpy())
        pred.append({"id": pid, "what": f"{c_} inserted exons: median coupling of ptc_nmd minus in_frame",
                     "effect": float(nmd.median() - inf.median()), "predicted_sign": -1, "p": p,
                     "n": f"ptc_nmd {len(nmd)} (median {nmd.median():+.3f}), in_frame {len(inf)} (median {inf.median():+.3f})"})
    r = rows_cq3(a, ex, gn)
    nmd, inf = r[r.consequence.eq("ptc_nmd")], r[r.consequence.eq("in_frame")]
    if min(len(nmd), len(inf)) < 5:
        pred.append({"id": "CQ3", "what": "not testable", "n": f"ptc_nmd {len(nmd)}, in_frame {len(inf)}"})
        pred.append({"id": "CQ4", "what": "not testable", "n": f"ptc_nmd {len(nmd)}, in_frame {len(inf)}"})
    else:
        auc, p = _mwu(inf.log2fc.to_numpy(), nmd.log2fc.to_numpy())       # P(in_frame change > ptc_nmd change)
        pred.append({"id": "CQ3", "what": "training cancers, R_up inserted exons: AUROC(host-gene change lower "
                     "for ptc_nmd than in_frame) - 0.5", "effect": auc - 0.5, "predicted_sign": 1, "p": p,
                     "n": f"ptc_nmd {len(nmd)} (median {nmd.log2fc.median():+.3f}), in_frame {len(inf)} "
                          f"(median {inf.log2fc.median():+.3f})"})
        rho = lambda d: spearmanr(d["shift"], d.log2fc).statistic
        both = pd.concat([nmd, inf]).reset_index(drop=True)
        obs = rho(nmd) - rho(inf)
        rng = np.random.default_rng(m51.SEED)
        bs = []
        for _ in range(BOOT):
            b = both.iloc[rng.integers(0, len(both), len(both))]
            bn, bi = b[b.consequence.eq("ptc_nmd")], b[b.consequence.eq("in_frame")]
            bs.append(rho(bn) - rho(bi) if min(len(bn), len(bi)) >= 3 else np.nan)
        bs = np.asarray(bs)
        bs = bs[np.isfinite(bs)]
        pred.append({"id": "CQ4", "what": "training cancers, R_up inserted exons: Spearman(exon shift, host-gene "
                     "change), ptc_nmd minus in_frame", "effect": float(obs), "predicted_sign": -1,
                     "p": m51.boot_p(bs), "ci_low": float(np.quantile(bs, 0.025)), "ci_high": float(np.quantile(bs, 0.975)),
                     "n": f"ptc_nmd {len(nmd)} (rho {rho(nmd):+.3f}), in_frame {len(inf)} (rho {rho(inf):+.3f})"})
    return m51.predictions(pred, "consequence")


def rows_cq3(a, ex, gn, groups=TRAIN) -> pd.DataFrame:
    """R_up inserted exons in the given cancers with their shift and host-gene change; one exon per
    gene and cancer (the largest shift)."""
    e = ex[ex.group.isin(groups) & ex.mclass.eq("R_up")][["group", "feature_index", "shift"]]
    e = e.merge(a[["feature_index", "gene_id", "symbol", "consequence"]], on="feature_index")
    e = e[e.consequence.isin(["ptc_nmd", "in_frame"])]
    g = gn[gn.project.isin(groups)].groupby(["project", "gene_id"]).log2fc.median().rename("log2fc").reset_index()
    e = e.merge(g, left_on=["group", "gene_id"], right_on=["project", "gene_id"]).dropna(subset=["shift", "log2fc"])
    return e.sort_values("shift", ascending=False).drop_duplicates(["group", "gene_id"])


def report() -> None:
    a = pd.read_parquet(LOCAL / "consequence_annotation.parquet")
    cp = pd.read_parquet(LOCAL / "consequence_coupling.parquet")
    gn = pd.read_parquet(LOCAL / "consequence_genes.parquet")
    nm = pd.read_parquet(LOCAL / "consequence_normals.parquet")
    ex = pd.read_parquet(LOCAL / "magnitude_exons.parquet")
    OUT.mkdir(parents=True, exist_ok=True)
    pr = tests(a, cp, ex, gn)
    # composition: unique exons measured in a training cancer, and held out; B weighted to R_up's PSI deciles
    out = []
    for grp, sub in [("training", ex[ex.group.isin(TRAIN)]), ("held out", ex[ex.group.eq("held out") & ex.n_mmr.ge(MIN_HELD)])]:
        u = sub.groupby("feature_index").median_proficient.median().rename("median_proficient").reset_index()
        u = u.merge(a[["feature_index", "mclass", "consequence"]], on="feature_index")
        for c_ in CLASSES:
            v = u[u.mclass.eq(c_)].consequence.value_counts().reindex(CONSEQUENCES, fill_value=0)
            out += [{"group": grp, "class": c_, "consequence": k, "n": int(n), "share": n / max(v.sum(), 1)} for k, n in v.items()]
        ru = u[u.mclass.eq("R_up")]
        edges = np.unique(np.quantile(ru.median_proficient.dropna(), np.linspace(0, 1, 11)))
        b = u[u.mclass.eq("B")].dropna(subset=["median_proficient"])
        bd, rd = np.digitize(b.median_proficient, edges[1:-1]), np.digitize(ru.median_proficient.dropna(), edges[1:-1])
        w = pd.Series(rd).value_counts(normalize=True).reindex(range(len(edges) - 1), fill_value=0)
        wb = pd.Series(bd).map(w / pd.Series(bd).value_counts(normalize=True)).fillna(0).to_numpy()
        for k in CONSEQUENCES:
            out.append({"group": grp, "class": "B, matched to R_up proficient PSI", "consequence": k,
                        "n": int((b.consequence.eq(k)).sum()), "share": float((wb * b.consequence.eq(k)).sum() / wb.sum())})
    comp = pd.DataFrame(out)
    comp.to_csv(OUT / "composition.csv", index=False, float_format="%.4g")
    # coupling by class and consequence
    per = cp.groupby("feature_index").agg(rho=("rho", "mean"), n_cancers=("rho", "size")).reset_index()
    per = per.merge(a[["feature_index", "mclass", "consequence", "symbol"]], on="feature_index")
    cs = per.groupby(["mclass", "consequence"]).rho.describe(percentiles=[0.25, 0.5, 0.75]).reset_index()
    cs.to_csv(OUT / "coupling_summary.csv", index=False, float_format="%.4g")
    per[per.mclass.ne("B")].to_csv(OUT / "coupling_repeat_exons.csv", index=False, float_format="%.4g")
    # CQ3 rows (per exon x cancer, aggregates over tumours) and held out
    r3 = rows_cq3(a, ex, gn)
    r3.drop(columns=["project"]).to_csv(OUT / "gene_change_training.csv", index=False, float_format="%.4g")
    ho = ex[ex.group.eq("held out") & ex.mclass.eq("R_up") & ex.n_mmr.ge(MIN_HELD)][["feature_index", "shift"]]
    ho = ho.merge(a[["feature_index", "gene_id", "symbol", "consequence"]], on="feature_index")
    hg = gn[~gn.project.isin(TRAIN)].groupby("gene_id").log2fc.median().rename("log2fc_median_over_tumours")
    ho = ho.join(hg, on="gene_id")
    ho.drop(columns=["feature_index"]).to_csv(OUT / "gene_change_held_out.csv", index=False, float_format="%.4g")
    # candidates: in-frame R_up exons switched on in MMR-deficient tumours
    tr = ex[ex.group.isin(TRAIN) & ex.mclass.eq("R_up")].assign(mmr_median=lambda d: d.median_proficient + d["shift"])
    tr = tr.merge(nm.rename(columns={"project": "group"}), on=["group", "feature_index"], how="left")
    tr = tr.merge(a[["feature_index", "consequence", "peptide", "dist", "exon_len"]], on="feature_index")
    hi, lo = CANDIDATE
    cand = tr[tr.mmr_median.ge(hi) & tr.median_proficient.le(lo) & tr.normal_median.le(lo)]
    cols = ["group", "symbol", "feature_id", "acc_T", "d_bp", "exon_len", "consequence", "peptide", "dist",
            "normal_median", "median_proficient", "mmr_median", "shift", "n_mmr", "n_normal"]
    cand[cols].sort_values(["consequence", "shift"], ascending=[True, False]).to_csv(
        OUT / "switched_on_exons.csv", index=False, float_format="%.4g")
    tr[cols].to_csv(OUT / "r_up_training_exons.csv", index=False, float_format="%.4g")
    hx = ex[ex.group.eq("held out") & ex.mclass.eq("R_up") & ex.n_mmr.ge(MIN_HELD)].merge(
        a[["feature_index", "consequence", "peptide"]], on="feature_index").assign(
        mmr_median=lambda d: d.median_proficient + d["shift"])
    hx[hx.mmr_median.ge(hi) & hx.median_proficient.le(lo)].drop(columns=["feature_index"]).to_csv(
        OUT / "switched_on_held_out.csv", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    log("composition:\n" + comp.pivot_table(index=["group", "class"], columns="consequence", values="share").round(3).to_string())
    log("coupling:\n" + cs.round(3).to_string())
    log(f"switched-on candidates:\n{cand[cols].round(3).to_string()}")


BANDS, BAND_NAMES = [-0.01, 0.02, 0.1, 0.3, 0.7, 1.01], ["≤ 0.02", "0.02–0.1", "0.1–0.3", "0.3–0.7", "> 0.7"]


def posthoc() -> None:
    """A-21 (exploratory, after CQ1-CQ4 were read): coupling by baseline PSI, the global-NMD control
    (do background poison exons rise in MMR-deficient tumours?) and R_up movers by consequence."""
    from scipy.stats import fisher_exact, mannwhitneyu
    a = pd.read_parquet(LOCAL / "consequence_annotation.parquet")
    cp = pd.read_parquet(LOCAL / "consequence_coupling.parquet")
    ex = pd.read_parquet(LOCAL / "magnitude_exons.parquet")
    rng = np.random.default_rng(m51.SEED)
    base = ex[ex.group.isin(TRAIN)].groupby("feature_index").median_proficient.median()
    d = a.assign(rho=a.feature_index.map(cp.groupby("feature_index").rho.mean()), base=a.feature_index.map(base))
    d = d.dropna(subset=["rho", "base"]).assign(band=lambda x: pd.cut(x.base, BANDS, labels=BAND_NAMES))
    rows = []
    for (c_, k, b), v in d[d.consequence.isin(["ptc_nmd", "in_frame", "canonical"])].groupby(
            ["mclass", "consequence", "band"], observed=True):
        rows.append({"class": c_, "consequence": k, "baseline_psi": b, "n": len(v), "median_rho": float(v.rho.median()),
                     "q25": float(v.rho.quantile(0.25)), "q75": float(v.rho.quantile(0.75))})
    strat = pd.DataFrame(rows)
    for b in BAND_NAMES:                                    # B poison minus in-frame within each band
        x = d[d.mclass.eq("B") & d.consequence.eq("ptc_nmd") & d.band.eq(b)].rho
        y = d[d.mclass.eq("B") & d.consequence.eq("in_frame") & d.band.eq(b)].rho
        if min(len(x), len(y)) >= 10:
            strat = pd.concat([strat, pd.DataFrame([{"class": "B", "consequence": "ptc_nmd minus in_frame",
                                                     "baseline_psi": b, "n": f"{len(x)} v {len(y)}",
                                                     "median_rho": float(x.median() - y.median()),
                                                     "p": float(mannwhitneyu(x, y).pvalue)}])])
    strat.to_csv(OUT / "coupling_by_baseline.csv", index=False, float_format="%.4g")
    log(f"coupling vs baseline PSI, all annotated exons: rho = {d[['base', 'rho']].corr('spearman').iloc[0, 1]:.3f}")
    rows = []
    for grp, sub in [("training", ex[ex.group.isin(TRAIN)]), ("held out", ex[ex.group.eq("held out") & ex.n_mmr.ge(MIN_HELD)])]:
        sub = sub.merge(a[["feature_index", "consequence"]], on="feature_index")
        sub = sub[sub.consequence.isin(["ptc_nmd", "in_frame"])].assign(band=lambda x: pd.cut(x.median_proficient, BANDS, labels=BAND_NAMES))
        for (c_, k, b), v in sub.groupby(["mclass", "consequence", "band"], observed=True):
            w = v["shift"].dropna().to_numpy()
            if len(w) < 5:
                continue
            bs = rng.choice(w, size=(m51.BOOT, len(w))).mean(1)
            rows.append({"group": grp, "class": c_, "consequence": k, "baseline_psi": b, "n": len(w),
                         "mean_shift": float(w.mean()), "ci_low": float(np.quantile(bs, 0.025)),
                         "ci_high": float(np.quantile(bs, 0.975)), "share_ge_005": float((w >= 0.05).mean())})
    pd.DataFrame(rows).to_csv(OUT / "nmd_control.csv", index=False, float_format="%.4g")
    tr = ex[ex.group.isin(TRAIN) & ex.mclass.eq("R_up")].merge(a[["feature_index", "consequence"]], on="feature_index")
    tr = tr[tr.consequence.isin(["ptc_nmd", "in_frame"])]
    k = pd.crosstab(tr.consequence, tr["shift"] >= 0.05).reindex(columns=[False, True], fill_value=0)
    res = fisher_exact(k.to_numpy())
    pd.DataFrame([{"consequence": c, "n": int(k.loc[c].sum()), "moving_ge_005": int(k.loc[c, True]),
                   "fisher_p_ptc_nmd_vs_in_frame": float(res.pvalue)} for c in k.index]).to_csv(
        OUT / "r_up_movers_by_consequence.csv", index=False, float_format="%.4g")
    log(f"R_up rows moving >= 0.05 by consequence:\n{k}\nFisher p = {res.pvalue:.3f}")


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["annotate", "compute", "report", "posthoc"])
    a = ap.parse_args()
    for step, fn in [("annotate", annotate), ("compute", compute), ("report", report), ("posthoc", posthoc)]:
        if a.only in (None, step):
            fn()


if __name__ == "__main__":
    main()
