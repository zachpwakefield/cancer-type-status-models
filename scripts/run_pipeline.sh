#!/usr/bin/env bash
# End-to-end rebuild: sources -> cohort/splits -> matrices -> validation -> models -> summary.
# Usage: bash scripts/run_pipeline.sh            (PYTHON env var overrides the interpreter)
# Runtime on a shared 12-core Mac: extraction ~10 min; models several hours (8 workers).
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-/Users/zacharywakefield/miniconda3/envs/survival_env_full/bin/python}"
export OMP_NUM_THREADS=1   # set before Python starts (see scripts/_common.py)
mkdir -p logs

run() { echo "== $* =="; "$PY" "$@" 2>&1 | tee "logs/$(basename "$1" .py).log"; }

"$PY" -m unittest discover -s tests
run scripts/00_register_sources.py
run scripts/01_build_cohort.py
# the raw tables are independent: extract them in three concurrent processes
"$PY" scripts/02_extract_matrices.py ri se mxe   > logs/02_extract_a.log 2>&1 & p1=$!
"$PY" scripts/02_extract_matrices.py a5ss a3ss hit > logs/02_extract_b.log 2>&1 & p2=$!
"$PY" scripts/02_extract_matrices.py gex afe ale  > logs/02_extract_c.log 2>&1 & p3=$!
wait $p1 && wait $p2 && wait $p3
run scripts/03_validate_data.py
run scripts/05_run_all.py --workers "${WORKERS:-8}"
run scripts/06_summarise_results.py
# interpretation layer (docs/BIOLOGY_METHODS.md); no model is refitted
run scripts/07_shap_importance.py --workers "${WORKERS:-6}"
for s in 08_bio_complementarity 09_bio_global_shifts 10_bio_pathways 11_bio_events 12_bio_regulators \
         13_bio_lineage 14_bio_site_probe 16_surv_inputs 18_ensemble 19_regulons; do
  run "scripts/$s.py"
done
# the survival layer shells out to splice-assay (its own environment; see docs/BIOLOGY_METHODS.md D-25)
run scripts/17_surv_analyze.py --workers "${WORKERS:-6}"
# event cards: everything known about the events the narratives name (needs 17 and 19)
run scripts/22_event_cards.py
# events that read another gene (D-27) and two published mechanisms checked (D-28); both need 22
run scripts/23_overlap_audit.py
run scripts/24_literature_checks.py
# the survival enrichment without the kidney cohorts, and leaving each cancer out (D-29)
run scripts/25_survival_robustness.py
# figures, one script per narrative
run scripts/20_figures_status.py
run scripts/21_figures_type.py
