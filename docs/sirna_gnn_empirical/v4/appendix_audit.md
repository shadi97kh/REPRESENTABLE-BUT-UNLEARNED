# Appendix completeness and main-claim mapping

Every complete historical finite-model proof is preserved verbatim; new exact input-cancellation and collision-floor proofs are complete. Calibration, endpoint derivatives, ensemble decomposition, metric and uncertainty definitions remain self-contained. Raw record dumps are replaced by compact complete summaries and exact evidence identifiers.

| Main object | Complete appendix support | Exact numerical / implementation evidence |
|---|---|---|
| `fig:workflow` | app:architecture; app:training | `sirna_gnn_empirical/v4/architecture.py; original/new fit records; workflow_identity.json` |
| `tab:data` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/source_counts_scores.csv` |
| `tab:cohorts` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/evaluation_units.csv` |
| `tab:adjudication` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/adjudication.csv` |
| `tab:optimization` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/new_fits.csv` |
| `tab:activity` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/focused_activity.csv` |
| `fig:robustness` | app:focusedresults; app:metrics | `v4/evaluation/activity_metrics.csv; v4/evaluation/seed_metrics.csv` |
| `tab:calibration` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/historical_calibration.csv` |
| `fig:calibration` | app:robustness; app:metrics | `v4/review_checks/calibration_decomposition.csv; v3/tables/direct_comparisons.csv` |
| `tab:ranking` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/ranking_history.csv` |
| `tab:pairranking` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/B2_unchanged.csv` |
| `tab:B3` | A2–A6, A9–A10 as detailed in the appendix roadmap | `v4/paper_tables/B3_visibility.csv` |
| `fig:visibility` | app:inputcollision; app:focusedresults | `v4/visibility/input_counts.csv; v4/visibility/B3_rectangle_visibility.csv` |

Main numerical prose dependencies:

| Claim | Complete derivation / evidence |
|---|---|
| Source concentration and positive/negative per-source R² | `v4/review_checks/source_concentration.csv`, `source_scores.csv`; A9 denominator definitions |
| 299/1,604 discrepancy, max 1.0852170325, 320 quarantines | `v4/adjudication/row_reconciliation.csv`, direct source cells and primary marginal checks; A9 |
| Primary-assay eligibility and 2,607-row coverage | `v4/dataset/`, A9 complete target definition |
| 792 new fits, ten seeds, 40 reused, updates | `v4/fit_summary.csv`, `fit_completion.json`, `evaluation/reused_fits.json`; A4/A9 |
| APP gap and 99.53% squared-bias fraction | `v4/review_checks/calibration_decomposition.csv`; full A9 expansion |
| Original S7 narrow improvement, no-message comparison | `v3/tables/direct_comparisons.csv`, `v4/review_checks/independent_conditional_comparisons.csv`; A5/A9 |
| New conditional comparisons and seed variability | `v4/evaluation/direct_comparisons.csv`, `seed_stability.csv`, `seed_paired_differences.csv`; A5/A10 |
| Ranking reversal and exact ties | `v4/ranking/` and `v4/evaluation/ranking_pools.csv`; A5/A9/A10 |
| B2 means, sign controls, correlation and uncertainty limits | `v3/tables/pair_metrics.csv`, `pair_comparisons.csv`, `B2_within_assay_variation.csv`; A2/A4/A5/A7 |
| B3 counts, nonzero predictions and cancellation | `v3/b3_primary/`, `v4/visibility/`; A6/A9/A10 complete proofs and dependence |
| Published blockers/S1 | `v4/source_checks/`, original primary inputs and pinned source; A6 |
| No new architecture, causal validation or five-law transfer | Exact A3 forward map and A9 scope; original theory preserved separately |

Actual boundary: {"main_start": 1, "main_end": 9, "statements_start": 10, "references_start": 11, "references_end": 14, "appendix_start": 15, "appendix_end": 53, "appendix_pages": 39, "total_pages": 53, "workflow_page": 2}. Zero undefined references and overfull boxes are required by the final build audit. All relevant original proof bodies pass byte-for-byte preservation; new proof bodies have no repository-only endpoint. The four-entry source and ZIP pass structural validation and a clean independent build.
