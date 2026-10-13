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
# subtype task (docs/SUBTYPE_METHODS.md): labels from the TCGA-CDR and two open-access GDC tables
run scripts/00b_fetch_gdc_labels.py
run scripts/30_subtype_labels.py
"$PY" scripts/31_subtype_extract.py ri se mxe   > logs/31_extract_a.log 2>&1 & p1=$!
"$PY" scripts/31_subtype_extract.py a5ss a3ss hit > logs/31_extract_b.log 2>&1 & p2=$!
"$PY" scripts/31_subtype_extract.py gex afe ale  > logs/31_extract_c.log 2>&1 & p3=$!
wait $p1 && wait $p2 && wait $p3
run scripts/33_subtype_run_all.py --workers "${WORKERS:-6}"
run scripts/34_subtype_summarise.py
run scripts/35_subtype_interpret.py
# splicing regulators, the site probe, the MSI repeats and BRCA (D-46, D-47, D-52..D-60;
# config/followup.yaml)
run scripts/37_subtype_regulators.py --workers "${WORKERS:-6}"
run scripts/38_subtype_site_probe.py --workers "${WORKERS:-6}"
run scripts/42_msi_repeats.py               # D-52..D-54; needs the GRCh38 FASTA (FASTA in the script)
run scripts/43_msi_robustness.py            # D-55
# D-61..D-64: MC3 and ABSOLUTE from step 00b; GTEx v8 sQTL pairs pinned in the script (fetched by hand)
run scripts/50_msi_followup.py --only geometry
run scripts/50_msi_followup.py --only load
run scripts/50_msi_followup.py --only nmd
run scripts/50_msi_followup.py --only branch
run scripts/50_msi_followup.py --only germline
run scripts/51_msi_code.py --only annotate       # D-65..D-67, label-blind steps first
run scripts/51_msi_code.py --only geometry
run scripts/51_msi_code.py --only specificity
run scripts/51_msi_code.py --only pancancer
run scripts/51_msi_code.py --only endometrium
run scripts/51_msi_code.py --only infiltration    # D-68
if [ -f data.nosync/ctsm_v1/external/cbioportal_msi/tcga_pancan_msi_scores.tsv ]; then
  run scripts/00c_fetch_msi_scores.py              # D-69: verifies the pinned TSV (fetching needs approval)
  run scripts/51_msi_code.py --only msicalls
fi
run scripts/52_msi_reversal.py --only annotate  # D-70, label-blind
run scripts/52_msi_reversal.py --only test
run scripts/52_msi_reversal.py --only posthoc   # A-16, exploratory
run scripts/54_msi_models.py --only features   # D-71, label-blind inputs
run scripts/54_msi_models.py --only fit
run scripts/54_msi_models.py --only posthoc    # A-18, exploratory
run scripts/55_msi_benchmark.py                # D-72, against MANTIS, MSIsensor and the MC3 count
run scripts/44_repeat_dna.py --only loci    # D-56: loci and the BAM manifest; the slices are fetched
run scripts/44_repeat_dna.py --only manifest #   by the data user with their own GDC token (see step 44)
if [ -d data.nosync/ctsm_v1/msi/slices ]; then
  run scripts/45_repeat_dna_reads.py --only genotype
  run scripts/45_repeat_dna_reads.py --only test
fi
# BRCA follow-ups (D-57..D-60): PanCanAtlas ABSOLUTE and MC3 tables from step 00b
run scripts/46_brca_germline.py --only scan
run scripts/46_brca_germline.py --only test
run scripts/46_brca_germline.py --only check       # A-13, A-14: the tumour-normal germline check
run scripts/46_brca_germline.py --only confirmed   # A-13: L1-L6 and A3 on germline-confirmed events
run scripts/47_brca_her2.py
run scripts/48_brca_mutations.py
run scripts/49_brca_her2_uncapped.py       # D-60: needs the SNP6 segments from step 00b
run scripts/36_subtype_figures.py
run scripts/53_msi_manuscript_figures.py   # M1, the MSI manuscript's schematic
# across the three tasks (docs/ACROSS_TASKS_METHODS.md D-48..D-51): nothing refitted
run scripts/40_cross_task.py
run scripts/41_figures_cross_task.py
