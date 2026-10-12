"""Step 22: event cards — everything the pipeline already knows about one event, assembled for
the handful of events each narrative names, so the write-up can be concrete about them.

Nothing is selected here. Both lists come from a rule the story already used:
    tumour status  the recurrent events of S5 (top 100 of >= 4 models), most models first,
                   one event per gene
    cancer type    the lineage atlas events of T6, largest |class median - median over classes|
                   first, one event per class

For each event the card holds its geometry, what the alternative segment does to the protein
(GENCODE v45, MANE > canonical transcript; InterPro features of that transcript), how it moves
between tumour and normal in every cancer that measured it, how the models read it (rank and
Spearman(value, SHAP) per model), which splicing-factor regulons it belongs to, and its survival
association. `*_protein_tracks.csv` holds the domains of each card's transcript for S9 and T8.

Usage: python scripts/22_event_cards.py [--n 10]
"""
from __future__ import annotations

import argparse
import re

import numpy as np
import pandas as pd
from _common import log

from ctsm.biology import ProteinMap, feature_values, full_importance, task_table
from ctsm.cohort import class_names
from ctsm.config import PAN, load_config


def pieces(view: str, feature_id: str) -> list[tuple[int, int]]:
    """The coordinate pairs of an rMATS id, 1-based inclusive. rMATS writes 0-based half-open
    starts, so every start moves forward by one; the terminal-exon views (AFE/ALE/HIT) carry a
    single already-1-based interval."""
    if view in ("afe", "ale", "hit"):
        m = re.search(r"chr[\dXYM]+:(\d+)-(\d+)", feature_id)
        return [(int(m.group(1)), int(m.group(2)))]
    head, *rest = feature_id.split(";")
    m = re.search(r"chr[\dXYM]+:(\d+)-(\d+)", head)
    out = [(int(m.group(1)) + 1, int(m.group(2)))]
    out += [(int(a) + 1, int(b)) for a, b in (p.split("-") for p in rest)]
    return out


def tracked_and_partner(parts: list[tuple[int, int]], region: str) -> tuple[tuple, tuple]:
    """For MXE, which of the two exons the PSI follows. rMATS writes them in its own order, the
    repo's tables follow the transcript-upstream one, so `region` decides."""
    a, b = (int(v) for v in region.split(":")[1].split("-"))
    first, second = parts[0], parts[1]
    return (first, second) if (a, b) == first else (second, first)


def alternative_segment(view: str, parts: list[tuple[int, int]]) -> tuple[int, int, str]:
    """The piece of sequence the event adds or removes, and what the two forms differ by.

    SE/AFE/ALE/HIT: the exon itself. MXE: the tracked exon, against its partner — what matters
    for the reading frame is the difference between the two, not either length. RI: the intron
    between the flanking exons. A3SS/A5SS: the extension the long form adds to the short exon
    (the repo's `region` is the whole long exon instead; the two differ by design)."""
    if view in ("se", "afe", "ale", "hit"):
        (a, b) = parts[0]
        return a, b, "exon"
    if view == "mxe":
        return parts[0][0], parts[0][1], "mxe"
    if view == "ri":
        return parts[1][1] + 1, parts[2][0] - 1, "intron"
    long_, short = parts[0], parts[1]
    return ((short[1] + 1, long_[1], "extension") if long_[0] == short[0]
            else (long_[0], short[0] - 1, "extension"))


def consequence(pm: ProteinMap, view: str, gene_id: str, region: str, feature_id: str) -> dict:
    """What the alternative segment does to the protein. The length that matters is the segment
    the two forms differ by: for MXE that is the difference between the two exons, because the
    swap replaces one with the other."""
    parts = pieces(view, feature_id)
    if view == "mxe":
        tracked, other = tracked_and_partner(parts, region)
        parts = [tracked, other, *parts[2:]]
    a, b, kind = alternative_segment(view, parts)
    out = {"segment": f"{a}-{b}", "segment_nt": b - a + 1, "segment_kind": kind}
    if view == "mxe":
        out["partner_nt"] = other[1] - other[0] + 1
        out["swap_nt"] = out["segment_nt"] - out["partner_nt"]
        out["swap_in_frame"] = out["swap_nt"] % 3 == 0
    elif view in ("afe", "ale", "hit"):
        # a terminal exon replaces one end of the protein rather than shifting a frame, so the
        # frame question does not apply and is left empty
        out["swap_nt"], out["swap_in_frame"] = None, None
    else:
        out["swap_nt"] = out["segment_nt"] if kind != "intron" else None
        out["swap_in_frame"] = (out["segment_nt"] % 3 == 0) if kind != "intron" else None
    c = ({"annotated_exon": False} if kind == "intron" else
         pm.exon_consequence(gene_id, parts[0][0], parts[0][1], require_features=True))
    keep = ("transcript_id", "aa_start", "aa_end", "protein_length_aa", "domains", "in_domain",
            "in_disordered", "in_transmembrane")
    # the protein map holds GENCODE v45 protein-coding transcripts only, so "not an exon" means
    # "not an exon of a coding transcript of this gene", which is itself worth recording
    status = ("retained intron" if kind == "intron" else
              "not in a coding transcript" if not c.get("annotated_exon") else
              "outside the CDS" if not c.get("coding") else
              "coding" if c.get("fully_coding") else "part coding")
    return {**out, "exon_status": status, **{k: c.get(k) for k in keep}}


def direction_by_model(cfg, task: str, view: str, cohorts: list[str], feature_id: str,
                       class_index=None) -> dict:
    """Rank and Spearman(value, SHAP) for this feature in each model that selected it. A positive
    direction means a higher value pushes the model towards the positive class (tumour, or the
    class explained)."""
    ranks, dirs, where = {}, {}, []
    for c in cohorts:
        imp = full_importance(cfg, task, view, c if task == "tumor_status" else PAN,
                              class_index=class_index)
        row = imp[imp["feature_id"].eq(feature_id)]
        if row.empty or not row["n_folds_selected"].iloc[0]:
            continue
        ranks[c] = int(row["rank"].iloc[0])
        d = row["direction_spearman"].iloc[0]
        if pd.notna(d):
            dirs[c] = float(d)
        where.append(c)
    pos = [v for v in dirs.values() if v > 0]
    return {"models_selecting": len(where), "best_rank": min(ranks.values()) if ranks else None,
            "median_rank": int(np.median(list(ranks.values()))) if ranks else None,
            "median_direction": float(np.median(list(dirs.values()))) if dirs else None,
            "share_direction_up": (len(pos) / len(dirs)) if dirs else None}


def regulons(core: pd.DataFrame, enriched: set[str], view: str, feature_id: str) -> dict:
    r = core[core["view"].eq(view) & core["feature_id"].eq(feature_id)]
    spec = r[r["factor"].isin(enriched)].sort_values("median_rho", key=abs, ascending=False)
    return {"n_regulons": len(r),
            "enriched_regulons": ", ".join(spec["factor"]),
            "regulon_rho": ", ".join(f"{x:+.2f}" for x in spec["median_rho"])}


def survival(dos: pd.DataFrame, view: str, feature_id: str) -> dict:
    r = dos[dos["view"].eq(view) & dos["feature_id"].eq(feature_id)]
    if r.empty:
        return {"survival_tested": False}
    r = r.iloc[0]
    return {"survival_tested": True, "survival_hits": r["hits"], "survival_cancers": r["cancers_tested"],
            "survival_best_cancer": r["best_cancer_survival"], "hr_per_sd": r["hr_per_sd"],
            "hr_low": r["ci_low_sd"], "hr_high": r["ci_high_sd"], "survival_p": r["cox_p"]}


def tracks(pm: ProteinMap, cards: pd.DataFrame) -> pd.DataFrame:
    """The InterPro features of each card's transcript, for the protein-track figures."""
    rows = []
    for r in cards.itertuples():
        tx = getattr(r, "transcript_id", None)
        if not isinstance(tx, str):
            continue
        f = pm.feat_by_tx.get(tx)
        if f is None:
            continue
        for g in f.itertuples():
            rows.append({"symbol": r.symbol, "view": r.view, "feature_id": r.feature_id,
                         "transcript_id": tx, "protein_length_aa": r.protein_length_aa,
                         "source_db": g.source_db, "name": g.name, "aa_start": g.start, "aa_end": g.end})
    return pd.DataFrame(rows)


def status_cards(cfg, pm: ProteinMap, n: int) -> pd.DataFrame:
    bio = cfg.results_dir / "biology"
    rec = pd.read_csv(bio / "04_events" / "recurrent_events.csv")
    P = pd.read_csv(bio / "04_events" / "paired_change_top_events_all_cancers.csv.gz")
    core = pd.read_csv(bio / "10_regulons" / "regulon_core_events.csv")
    fe = pd.read_csv(bio / "10_regulons" / "regulon_enrichment_by_factor.csv")
    enriched = set(fe.loc[fe["q_fisher"].lt(0.05) & fe["odds_ratio_haldane"].gt(1), "factor"])
    dos = pd.read_csv(bio / "08_survival" / "event_dossiers.csv")
    picks = (rec.sort_values(["n_cancers", "sum_norm_importance"], ascending=False)
             .drop_duplicates("symbol").head(n))
    rows = []
    for r in picks.itertuples():
        d = P[P["view"].eq(r.view) & P["feature_id"].eq(r.feature_id)]
        big = d.reindex(d["median_paired_delta"].abs().sort_values(ascending=False).index).head(1)
        rows.append({
            "symbol": r.symbol, "view": r.view, "gene_id": r.gene_id, "region": r.region,
            "feature_id": r.feature_id,
            "n_models_top100": r.n_cancers, "top_in": r.cancers, "best_cancer": r.best_cancer,
            "units": "HIT index" if r.view == "hit" else "PSI",
            "majority_sign": r.majority_sign, "n_cancers_measured": r.n_cancers_measured,
            "n_cancers_same_sign": r.n_cancers_same_sign_as_majority,
            "median_delta": r.median_delta_across_cancers,
            "largest_cancer": big["cancer"].iloc[0] if len(big) else None,
            "largest_delta": big["median_paired_delta"].iloc[0] if len(big) else None,
            "largest_normal": big["median_normal"].iloc[0] if len(big) else None,
            "largest_tumour": big["median_tumour"].iloc[0] if len(big) else None,
            **consequence(pm, r.view, r.gene_id, r.region, r.feature_id),
            **direction_by_model(cfg, "tumor_status", r.view, r.cancers.split(","), r.feature_id),
            **regulons(core, enriched, r.view, r.feature_id),
            **survival(dos, r.view, r.feature_id)})
        log(f"  {r.symbol} {r.view}")
    return pd.DataFrame(rows)


def class_medians(cfg, view: str, feature_id: str, tab: pd.DataFrame) -> pd.Series:
    """The event's median value in each of the 34 classes. Sample values stay local (D-12);
    only the medians leave this function."""
    v = feature_values(cfg, view, tab["file_id"], [feature_id])[:, 0]
    return pd.Series(v, index=tab["label"].to_numpy()).groupby(level=0).median()


def class_contrast(cfg, view: str, feature_id: str, tab: pd.DataFrame, cls: str) -> dict:
    """Where the marked class sits against every other class, and how far it is from its
    nearest neighbour -- the margin a single event gives on its own."""
    med = class_medians(cfg, view, feature_id, tab)
    others = med.drop(labels=[cls])
    runner = others.reindex((others - med[cls]).abs().sort_values(ascending=False).index)
    near = others.sub(med[cls]).abs().idxmin()
    return {"class_median": float(med[cls]), "median_over_classes": float(med.median()),
            "deviation": float(med[cls] - med.median()),
            "nearest_class": str(near), "nearest_class_median": float(others[near]),
            "margin": float(med[cls] - others[near]),
            "furthest_class": str(runner.index[0]), "furthest_class_median": float(runner.iloc[0])}


def type_cards(cfg, pm: ProteinMap, n: int, per_class: int = 5) -> pd.DataFrame:
    """One event per class, then the n classes with the clearest one. Stage 1 is the models'
    own ranking: the `per_class` events with the largest mean |SHAP| in that class's explanation,
    over all 8 splicing views. Stage 2 keeps, of those, the event whose class median is furthest
    from the nearest other class -- a top-ranked feature need not separate the class on its own,
    because the model reads it together with the others, and an example that does not separate
    cannot be written about concretely. Stage 2 chooses which events to show, nothing else."""
    bio = cfg.results_dir / "biology"
    core = pd.read_csv(bio / "10_regulons" / "regulon_core_events.csv")
    fe = pd.read_csv(bio / "10_regulons" / "regulon_enrichment_by_factor.csv")
    enriched = set(fe.loc[fe["q_fisher"].lt(0.05) & fe["odds_ratio_haldane"].gt(1), "factor"])
    tab = task_table(cfg, "cancer_type", PAN)
    classes = class_names(tab)
    views = [v for v in cfg.views if v != "gex"]
    rows = []
    for k, cls in enumerate(classes):
        cand = []
        for view in views:
            imp = full_importance(cfg, "cancer_type", view, PAN, class_index=k)
            top = imp[imp["n_folds_selected"].gt(0)].nsmallest(per_class, "rank")
            if not top.empty:
                cand.append(top.assign(view=view, cls=cls))
        if not cand:
            continue
        cand = pd.concat(cand).nlargest(per_class, "mean_abs_shap")
        scored = [{**r._asdict(), **class_contrast(cfg, r.view, r.feature_id, tab, cls)}
                  for r in cand.itertuples()]
        rows.append(max(scored, key=lambda z: abs(z["margin"])))
        log(f"  {cls}: {rows[-1]['symbol']} {rows[-1]['view']} margin {rows[-1]['margin']:+.2f}")
    # no event twice: neighbouring classes (GBM/LGG, COAD/READ) often share one
    picks = sorted(rows, key=lambda z: -abs(z["margin"]))
    seen, out = set(), []
    for r in picks:
        if r["feature_id"] in seen or len(out) == n:
            continue
        seen.add(r["feature_id"])
        out.append({
            "symbol": r["symbol"], "view": r["view"], "gene_id": r["gene_id"], "region": r["region"],
            "feature_id": r["feature_id"], "marks_class": r["cls"], "rank_in_class": int(r["rank"]),
            "units": "HIT index" if r["view"] == "hit" else "PSI",
            "mean_abs_shap": r["mean_abs_shap"], "direction_spearman": r["direction_spearman"],
            **{k: r[k] for k in ("class_median", "median_over_classes", "deviation", "nearest_class",
                                 "nearest_class_median", "margin", "furthest_class", "furthest_class_median")},
            **consequence(pm, r["view"], r["gene_id"], r["region"], r["feature_id"]),
            **regulons(core, enriched, r["view"], r["feature_id"])})
    cards = pd.DataFrame(out)
    med = pd.DataFrame({f"{r.symbol} {r.view.upper()}": class_medians(cfg, r.view, r.feature_id, tab)
                        for r in cards.itertuples()}).T
    med.index.name = "event"
    return cards, med


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10, help="events per narrative")
    a = ap.parse_args()
    cfg = load_config()
    out = cfg.results_dir / "biology" / "11_event_cards"
    out.mkdir(parents=True, exist_ok=True)
    pm = ProteinMap(cfg)

    log("tumour status")
    s = status_cards(cfg, pm, a.n)
    s.to_csv(out / "tumour_status_cards.csv", index=False)
    tracks(pm, s).to_csv(out / "tumour_status_protein_tracks.csv", index=False)

    log("cancer type")
    t, med = type_cards(cfg, pm, a.n)
    t.to_csv(out / "cancer_type_cards.csv", index=False)
    med.to_csv(out / "cancer_type_class_medians.csv")
    tracks(pm, t).to_csv(out / "cancer_type_protein_tracks.csv", index=False)

    show = ["symbol", "view", "segment_nt", "swap_nt", "swap_in_frame", "aa_start", "aa_end",
            "protein_length_aa", "in_domain", "in_disordered", "median_direction", "enriched_regulons"]
    log("tumour-status cards:\n" + s[[c for c in show if c in s]].round(3).to_string(index=False))
    log("cancer-type cards:\n" + t[[c for c in show if c in t]].round(3).to_string(index=False))
    log(f"wrote {len(s)} + {len(t)} cards to {out}")


if __name__ == "__main__":
    main()
