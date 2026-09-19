# Current claim audit

An APP GNN improvement over the strongest observed fitted baseline is not established. This is retrospective follow-up in one related patent family, not independent external confirmation. Sequence-specific measured chemistry prediction is not established under the declared B2 mean-control and within-assay diagnostic.

| Claim | Evidence | Exact limit |
|---|---|---|
| Current outcome and strongest observed baseline | tables/strongest_observed_baseline_comparisons.csv | Conditional retrospective cohort comparison; strongest observed baseline is selected descriptively from these test metrics, not a validation-selected reference. |
| Activity constants/discrimination/calibration | tables/activity_metrics.csv | All weighting conventions separate; oracle never fitted. |
| Seed stability | tables/seed_stability.csv | Ten seeds on fixed data; deterministic methods once. |
| Support defect and correction | tables/final_parameter_support.csv | Support presence and activity-loss gradients; no universal support adequacy or full failure explanation. |
| B2 measured bundle | tables/pair_metrics.csv | Joint 6–9 change; no isolated GNA, causal effect or mean-difference SE. |
| B2 within-assay diagnostic | tables/B2_within_assay_variation.csv | Evaluation-label centering is labeled diagnostic; no deployed recalibration. |
| Candidate ranking | tables/ranking_pools.csv | Exact ties; shared one-family pools; no independent-pool p-value. |
| Primary B3 measured interactions | b3_primary/interaction_predictions.csv | One background, common reference, source-held-out; no independent-rectangle interval. |
| Exact published-method comparisons | sources/reproduction_status.json | Blocked; zero fitted published predictors. |
| S1 evidence | sources/S1_status.json | Acquired and quarantined, not absent. |
| Primary/release discrepancy | sources/primary_release_discrepancies.json | Conditional identity mapping; frozen released benchmark not silently repaired. |
| Fit execution | tables/all_fit_summary.csv | Actual updates/checkpoints/failures; no completed-hash check counted as fit. |
| Five-law theory | preserved_theory/manifest.json (paper-relative) | Separate preserved artifact; no theorem coverage of GNN; backend no-go and historical reference-efficiency question unchanged. |

Every numeric main-table cell has a source selector/field in `main_cell_provenance.json`; every plot has mark-level CSV and source/code hashes in `figure_provenance.json`. No architecture novelty, new wet-lab validation, clinical benefit, private-anchor calibration or causal-mechanism claim is made.
