"""Step 25 (D-29): does the survival result of D-25 depend on the two kidney cohorts?

TUMOUR_STATUS.md lists this as the check that would change the reading: the FDR-level result is
94% KIRC and KIRP, and the set-level gradient (every candidate -> top 10 in 2+ cancers ->
recurrent) is pooled over all 15 cancers. This step recomputes the enrichment of step 17,
with the same function and the same nested candidate sets, over three scopes:
  all cancers               reproduces 08_survival/survival_enrichment.csv
  without KIRC and KIRP     the check named in advance
  leave one cancer out      15 scopes, to show whether any single cohort carries the result
No test is re-run; the per-cell statistics of step 17 are only re-pooled.

Written: results/biology/14_survival_robustness/{enrichment_by_scope,leave_one_out}.csv
"""
from __future__ import annotations

import importlib.util

import pandas as pd
from _common import REPO, log

_spec = importlib.util.spec_from_file_location("surv", REPO / "scripts" / "17_surv_analyze.py")
surv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(surv)

D = REPO / "results" / "biology" / "08_survival"
OUT = REPO / "results" / "biology" / "14_survival_robustness"
KIDNEY = ["KIRC", "KIRP"]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    idx = pd.read_csv(D / "event_index.csv", usecols=["event_key", "role", "recurrent", "n_models"])
    cells = pd.read_csv(D / "event_survival_cells.csv.gz",
                        usecols=["event_key", "cohort", "endpoint", "arm", "cox_p", "km_p", "cox_ok", "km_ok"])
    cells = cells.merge(idx, on="event_key", how="left").reset_index(drop=True)
    scopes = {"all cancers": cells, "without KIRC and KIRP": cells[~cells["cohort"].isin(KIDNEY)]}
    by_scope = pd.concat([surv.enrichment(c.reset_index(drop=True), None).assign(scope=s)
                          for s, c in scopes.items()], ignore_index=True)
    by_scope.to_csv(OUT / "enrichment_by_scope.csv", index=False)
    loo = []
    for c in sorted(cells["cohort"].unique()):
        e = surv.enrichment(cells[cells["cohort"] != c].reset_index(drop=True), None)
        loo.append(e[e["test"] == "Cox"].assign(dropped=c))
    loo = pd.concat(loo, ignore_index=True)
    loo.to_csv(OUT / "leave_one_out.csv", index=False)
    show = by_scope[(by_scope["test"] == "Cox") & (by_scope["endpoint"] == "OS")]
    log("\n" + show[["scope", "arm", "candidate_set", "n_candidate", "rate_candidate", "rate_background",
                     "odds_ratio", "p_fisher"]].round(4).to_string(index=False))
    rng = (loo[loo["endpoint"] == "OS"].groupby(["arm", "candidate_set"])["odds_ratio"]
           .agg(["min", "median", "max"]).round(3))
    log("leave-one-cancer-out, OS Cox odds ratio:\n" + rng.to_string())


if __name__ == "__main__":
    main()
