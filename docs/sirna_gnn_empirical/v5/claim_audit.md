# Claim audit

| Main location | Claim | Essential qualification | Complete support |
|---|---|---|---|
| Abstract; §§1,9 | Empirical diagnostic, no new learned method or biological validation | G1 reduction; conventional architecture; unchanged negative controls | method_and_proof.md; novelty_matrix.md; A3,A6,A9 |
| §2 | 153,345 scalars, two layers, training masks/support, both objectives and derivatives | Actual forward/backward checks; retained exact specification | sirna_gnn_empirical/v4/architecture.py; A3–A4; v5/correctness/results.json |
| §3 | 1,924/2,927 source concentration; some other-source positive R² | Source-scoring subsets do not imply source-excluded training | v4/review_checks/source_concentration.csv; source_scores.csv; A2,A9 |
| §4 | 299/1,604 differences; maximum 1.085217; primary-assay sensitivity not certified label repair | All large discrepancies and 32 fixed controls checked in v4; no new identity rule | v4/adjudication/row_reconciliation.csv; direct_cell_review.csv; A9 |
| §4 | 792 focused fits; ten final seeds; 40 unaffected checkpoint reuses | No new fits in v5; 3,564 parent fit hashes checked | v4/fit_summary.csv; v4/fit_completion.json; A4,A9 |
| §5 | Focused activity metrics, constants, conditional intervals | 36 score groups and ensemble identities independently recomputed | v5/reuse/recomputed_metrics.csv; v4/evaluation/direct_comparisons.csv; A5,A10 |
| §6 | 99.53% squared-bias decomposition; original S7 routing improvement | Descriptive decomposition and conditional uncertainty; not causal bias explanation | v4/review_checks/calibration_decomposition.csv; independent_conditional_comparisons.csv; A5,A9 |
| §7 | Ranking reversals on identical IDs; exact ties; 17 overlapping pools | Retrospective protocol changes, not independent test evidence | v4/ranking/; v4/evaluation/ranking_pools.csv; A5,A9,A10 |
| §7 | B2 does not establish sequence-specific chemistry prediction | 156 pairs/12 sequence components; joint four-position edit; majority control; unresolved replicate covariance | v3/tables/pair_metrics.csv; pair_comparisons.csv; B2_within_assay_variation.csv; A4,A5,A7 |
| §8 | 165 states/140 dependent B3 rectangles; 90 exact inputs; zero forced cancellations | 0/140 is not zero predicted interactions; no across-background inference | v3/b3_primary/; v4/visibility/; A6,A9,A10 |
| §9 | Published reproductions incomplete; S1 acquired/quarantined | No newly resolved geometry/token/orientation asset; zero exact published fits | v4/source_checks/; A6 |
| §9 | Focused S7 GNN beats no-message conditionally; tree/CNN still better observed MSE | No general GNN advantage; negative APP/S7 R² retained | v4/evaluation/direct_comparisons.csv; v5/tables/main_comparison.csv; A10 |
| §9 | Five-law theory separate; backend no-go; reference efficiency open | No reduction or theorem coverage asserted | A8–A9; preserved historical theory outside this manuscript |

All original main numerical tables and three result figures remain byte-identical in TeX/data or asset content. New work changes workflow, scope wording and related-work citations. No new performance claim is introduced.
