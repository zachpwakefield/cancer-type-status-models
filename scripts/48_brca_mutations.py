"""Step 48 (D-59): splice-site mutations and SF3B1 in the BRCA tumours.

  C1, C2  every MC3 splice-site SNV at +1/+2 or -1/-2 is mapped to the exon whose donor (+) or
          acceptor (-) it hits (the MC3 transcript's intron number, located in GENCODE v45 with the
          same transcript ID and exon count) and to the main-matrix SE events of that exon; the
          mutated tumour's PSI percentile among the other cohort tumours
  C3      SF3B1 hotspot tumours against the rest on a cryptic A3SS score

Outputs: data.nosync/ctsm_v1/brca/splice_site_mutations.parquet (per tumour, local);
         results/subtype/brca/mutations/ (aggregates; cells under 5 tumours left blank)
Usage: python scripts/48_brca_mutations.py
"""
from __future__ import annotations

import gzip
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
from _common import REPO, log
from scipy.stats import binomtest, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from ctsm.biology import cohort_matrix
from ctsm.config import load_config

_spec = importlib.util.spec_from_file_location("msi42", Path(__file__).with_name("42_msi_repeats.py"))
m42 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m42)

DATA = m42.DATA
EXT = DATA / "external" / "gdc_pancanatlas"
OUT = REPO / "results" / "subtype" / "brca" / "mutations"
KEEP_FLAGS = {"PASS", "wga", "native_wga_mix"}
HOTSPOT = {622, 625, 626, 662, 663, 666, 700, 741, 742}
SPLICE = re.compile(r"^c\.[-*]?\d+([+-])([12])[ACGT]>[ACGT]$")
TRUNCATING = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"}


def cohort_calls(samples: set[str]) -> pd.DataFrame:
    cols = ["Hugo_Symbol", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
            "Feature", "INTRON", "FILTER"]
    m = pd.read_csv(EXT / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=cols, dtype=str)
    m["s15"] = m.Tumor_Sample_Barcode.str[:15]
    m = m[m.s15.isin(samples)]
    ok = m.FILTER.str.split(",").map(lambda f: set(f) <= KEEP_FLAGS)
    log(f"MC3: {len(m):,} calls in {m.s15.nunique()} cohort tumours; {int(ok.sum()):,} kept by FILTER")
    return m[ok].copy()


def gtf_exons(transcripts: set[str], gtf: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "exon":
                continue
            t = re.search(r'transcript_id "([^".]+)', f[8]).group(1)
            if t in transcripts:
                n = int(re.search(r"exon_number (\d+)", f[8]).group(1))
                rows.append((t, f[0], int(f[3]), int(f[4]), f[6], n))
    return pd.DataFrame(rows, columns=["transcript", "chrom", "start", "end", "strand", "exon_number"])


def se_index(F: pd.DataFrame) -> pd.DataFrame:
    """SE events by exon: region chrom:start-end is 1-based inclusive."""
    r = F.region.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
    return pd.DataFrame({"col": np.arange(len(F)), "chrom": r[0], "start": r[1].astype(int), "end": r[2].astype(int),
                         "symbol": F.symbol.to_numpy()})


def splice_effects(cfg, calls: pd.DataFrame) -> pd.DataFrame:
    X, F, tab = cohort_matrix(cfg, "se", "tumor_status", "BRCA")
    tum = tab.role.eq("tumour").to_numpy()
    s15 = tab.sample_id.str[:15].to_numpy()
    sp = calls[calls.Variant_Classification.eq("Splice_Site") & calls.Variant_Type.eq("SNP")].copy()
    mt = sp.HGVSc.str.extract(SPLICE.pattern)
    sp["side"], sp["pos"] = mt[0], mt[1]
    sp = sp[sp.side.notna() & sp.INTRON.str.contains("/", na=False)]
    sp[["intron", "n_introns"]] = sp.INTRON.str.split("/", expand=True).astype(int)
    sp["exon"] = np.where(sp.side.eq("+"), sp.intron, sp.intron + 1)
    ex = gtf_exons(set(sp.Feature), cfg.gencode_gtf)
    n_ex = ex.groupby("transcript").size()
    sp = sp[sp.Feature.map(n_ex).eq(sp.n_introns + 1)]
    sp = sp.merge(ex, left_on=["Feature", "exon"], right_on=["transcript", "exon_number"])
    idx = se_index(F)
    sp = sp.merge(idx, on=["chrom", "start", "end"], suffixes=("", "_se"))
    rows = []
    for (s, mut), d in sp.groupby(["s15", "HGVSc"]):
        rk = np.flatnonzero((s15 == s) & tum)
        if len(rk) == 0:
            continue
        pct = []
        for c in d.col.unique():
            x = X[:, c]
            me = x[rk[0]]
            others = x[tum & (s15 != s) & np.isfinite(x)]
            if np.isfinite(me) and len(others) >= 50:
                pct.append(((others < me).sum() + 0.5 * (others == me).sum()) / len(others))
        if pct:
            r = d.iloc[0]
            rows.append({"symbol": r.Hugo_Symbol, "hgvsc": mut, "side": "donor" if r.side == "+" else "acceptor",
                         "exon_region": f"{r.chrom}:{r.start}-{r.end}", "se_events": len(pct), "percentile": float(np.mean(pct)),
                         "s15": s})
    out = pd.DataFrame(rows)
    log(f"splice-site SNVs {len(sp.drop_duplicates(['s15', 'HGVSc']))} mapped to an SE exon; {len(out)} with PSI")
    return out


def cryptic_a3ss(cfg) -> tuple[pd.Series, pd.DataFrame]:
    X, F, tab = cohort_matrix(cfg, "a3ss", "tumor_status", "BRCA")
    q = F.feature_id.str.split("#", expand=True)
    lg = q[1].str.split(":", expand=True)[1].str.split("-", expand=True).astype(int)
    fl = q[2].str.split(";", expand=True)
    sh = fl[1].str.split("-", expand=True).astype(int)
    plus = fl[0].eq("+").to_numpy()
    seg = np.where(plus, sh[0] - lg[0], lg[1] - sh[1])           # nt between the two acceptors
    nor = tab.role.eq("normal").to_numpy()
    med_n = np.nanmedian(X[nor], axis=0)
    pick = (seg >= 10) & (seg <= 30) & (med_n < 0.2)
    tum = tab.role.eq("tumour").to_numpy()
    P = pd.DataFrame(X[tum][:, pick]).rank(pct=True)
    score = pd.Series(P.mean(axis=1).to_numpy(), index=tab.sample_id.str[:15].to_numpy()[tum])
    log(f"cryptic A3SS events: {int(pick.sum())}")
    return score[~score.index.duplicated()], F[pick]


def main() -> None:
    cfg = load_config()
    OUT.mkdir(parents=True, exist_ok=True)
    pam = pd.read_parquet(DATA / "subtype" / "tasks" / "subtype__BRCA_PAM50.parquet")
    hist = pd.read_parquet(DATA / "subtype" / "tasks" / "subtype__BRCA.parquet")
    samples = set(pam.sample_id.str[:15])
    calls = cohort_calls(samples)
    se = splice_effects(cfg, calls)
    pred = []

    def add(pid, effect, sign, p, n, what):
        pred.append({"id": pid, "what": what, "effect": effect, "predicted_sign": sign, "p": p, "n": n})

    for pid, d in [("C1", se), ("C2", se[se.symbol.eq("CDH1")])]:
        k = int((d.percentile < 0.5).sum())
        add(pid, k / len(d) - 0.5, 1, binomtest(k, len(d)).pvalue,
            f"{k} of {len(d)} mutations below the median; median percentile {d.percentile.median():.3f}",
            "share of mutated tumours with the exon below the cohort median, minus 0.5")
    score, ev = cryptic_a3ss(cfg)
    sf = calls[calls.Hugo_Symbol.eq("SF3B1") & calls.Variant_Classification.eq("Missense_Mutation")].copy()
    sf["codon"] = sf.HGVSp_Short.str.extract(r"p\.[A-Z](\d+)")[0].astype(float)
    hot = set(sf.loc[sf.codon.isin(HOTSPOT), "s15"])
    called = score[score.index.isin(samples)]
    y = called.index.isin(hot)
    r = mannwhitneyu(called[y], called[~y])
    add("C3", r.statistic / (y.sum() * (~y).sum()) - 0.5, 1, float(r.pvalue), f"hotspot {int(y.sum())}, rest {int((~y).sum())}; "
        f"{len(ev)} events", "AUROC(cryptic A3SS score higher in SF3B1-hotspot tumours) - 0.5")
    pr = pd.DataFrame(pred)
    pr["q"] = multipletests(pr.p, method="fdr_bh")[1]
    pr["verdict"] = [m42.verdict(e, s, q) for e, s, q in zip(pr.effect, pr.predicted_sign, pr.q)]
    pr.to_csv(OUT / "predictions.csv", index=False)
    log("predictions:\n" + pr[["id", "effect", "p", "q", "verdict", "n"]].to_string(index=False))

    # ---- exploratory (aggregates)
    (DATA / "brca").mkdir(parents=True, exist_ok=True)
    se.to_parquet(DATA / "brca" / "splice_site_mutations.parquet", index=False)      # per tumour: local only
    g = se.groupby("symbol").percentile.agg(["size", "median", lambda s: (s < 0.5).mean()]).reset_index()
    g.columns = ["symbol", "n", "median_percentile", "share_below_median"]
    g.loc[g.n < 5, ["median_percentile", "share_below_median"]] = np.nan
    g.to_csv(OUT / "c1_by_gene.csv", index=False, float_format="%.4g")
    by_side = se.groupby("side").percentile.agg(["size", "median", lambda s: (s < 0.5).mean()]).reset_index()
    by_side.columns = ["side", "n", "median_percentile", "share_below_median"]
    by_side.to_csv(OUT / "c1_by_side.csv", index=False, float_format="%.4g")
    lab = pam.assign(s15=pam.sample_id.str[:15]).set_index("s15")["label"]
    hl = hist.assign(s15=hist.sample_id.str[:15]).set_index("s15")["label"]
    cdh1 = calls[calls.Hugo_Symbol.eq("CDH1")]
    status = pd.Series("none", index=sorted(samples))
    status[status.index.isin(cdh1[cdh1.Variant_Classification.isin(TRUNCATING - {"Splice_Site"})].s15)] = "truncating"
    status[status.index.isin(cdh1[cdh1.Variant_Classification.eq("Splice_Site")].s15)] = "splice site"
    status = status[status.index.isin(set(calls.s15) | set(cdh1.s15))]          # tumours sequenced in MC3
    tab = pd.DataFrame({"pam50": lab, "histology": hl, "cdh1": status, "sf3b1_hotspot": pd.Series(True, index=sorted(hot))})
    tab["sf3b1_hotspot"] = tab.sf3b1_hotspot.fillna(False).astype(bool)
    tab = tab[tab.index.isin(samples)]
    t1 = pd.crosstab(tab.pam50, tab.sf3b1_hotspot)
    t2 = pd.crosstab(tab.histology.fillna("not in histology cohort"), tab.cdh1.fillna("not in MC3"))
    t1.to_csv(OUT / "sf3b1_hotspot_by_pam50.csv")
    t2.to_csv(OUT / "cdh1_status_by_histology.csv")
    log(f"SF3B1 hotspot by PAM50:\n{t1}\nCDH1 by histology:\n{t2}")
    cdh1_exons(cfg, status, hl)


def cdh1_exons(cfg, status: pd.Series, hl: pd.Series) -> None:
    """Lobular tumours: each CDH1 SE event's median PSI and lowest-decile PSI by CDH1 status (cells < 5 blank)."""
    X, F, tab = cohort_matrix(cfg, "se", "tumor_status", "BRCA")
    s15 = tab.sample_id.str[:15]
    keep = tab.role.eq("tumour").to_numpy() & s15.map(hl).eq("Lobular").to_numpy() & s15.isin(status.index).to_numpy()
    st = s15[keep].map(status).to_numpy()
    rows = []
    for j in np.flatnonzero(F.symbol.eq("CDH1").to_numpy()):
        for g in ["none", "truncating", "splice site"]:
            x = X[keep, j][st == g]
            x = x[np.isfinite(x)]
            rows.append({"region": F.region.iat[j], "cdh1": g, "n": len(x),
                         "median_psi": np.median(x) if len(x) >= 5 else np.nan,
                         "p10_psi": np.percentile(x, 10) if len(x) >= 5 else np.nan})
    r = pd.DataFrame(rows)
    r.to_csv(OUT / "cdh1_exons_in_lobular_by_status.csv", index=False, float_format="%.4g")
    log(f"CDH1 exons in lobular tumours:\n{r.to_string(index=False)}")


if __name__ == "__main__":
    main()
