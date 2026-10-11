# Data validation

All checks passed (`scripts/03_validate_data.py`; numbers in `validation.json`).

| Check | Result |
|---|---|
| V1 cohort | 9,729 samples from 9,075 cases (9,034 tumour, 695 normal); one file per case and role; every view's matrix rows equal the cohort |
| V2 splits | 16 task cohorts: case-grouped outer and inner folds, every sample tested once |
| V3 raw values | stored = raw CSV for 9,971 features x 300 samples across 9 views (max abs error 6.3e-07); AFE/ALE fills are all flagged structural zeros |
| V4 GEX | log2(CPM+1) vs official GDC log2(TPM+1), 300 primary tumours x 23,640 genes: median Spearman 0.964 (min 0.847) |
| V5 pre-filter | thresholds = ceil(fraction x smallest inner-train set) for all 16 cohorts x 9 views |
