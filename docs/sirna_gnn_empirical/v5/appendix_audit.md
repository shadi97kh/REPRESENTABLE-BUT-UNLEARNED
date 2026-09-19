# Appendix completeness audit — manuscript v10

The source retains all five complete parent proof bodies verbatim. No new-method theorem is claimed in the main text. All finite-model definitions, forward/no-message arguments, endpoint derivatives, metrics, covariance qualifications, ensemble decomposition and collision proofs remain in the compiled appendix. The separate exploratory exact-range proof is complete in method_and_proof.md; it is not attached to the GNN as a guarantee.

| Main object (PDF page) | Complete appendix support (PDF page) | Exact archive evidence |
|---|---|---|
| fig:workflow (p. 2) | app:architecture (p. 20); app:training (p. 24) | `sirna_gnn_empirical/v4/architecture.py; v4/protocol.json; v5/workflow_identity.json` |
| tab:data (p. 3) | app:data (p. 15); app:robustness (p. 42) | `v4/paper_tables/source_counts_scores.csv` |
| tab:cohorts (p. 3) | app:data (p. 15); app:b3 (p. 30) | `v4/paper_tables/evaluation_units.csv` |
| tab:adjudication (p. 4) | app:robustness (p. 42) | `v4/adjudication/row_reconciliation.csv; v4/paper_tables/adjudication.csv` |
| tab:optimization (p. 4) | app:training (p. 24); app:robustness (p. 42) | `v4/fit_summary.csv; v4/fit_completion.json; v4/paper_tables/new_fits.csv` |
| tab:activity (p. 5) | app:metrics (p. 27); app:focusedresults (p. 46) | `v4/evaluation/activity_metrics.csv; v5/reuse/recomputed_metrics.csv` |
| fig:robustness (p. 5) | app:metrics (p. 27); app:focusedresults (p. 46) | `v4/evaluation/seed_metrics.csv; v4/evaluation/activity_metrics.csv` |
| tab:calibration (p. 6) | app:metrics (p. 27); app:results (p. 33) | `v4/paper_tables/historical_calibration.csv` |
| fig:calibration (p. 6) | app:metrics (p. 27); app:robustness (p. 42) | `v4/review_checks/calibration_decomposition.csv; v3/tables/direct_comparisons.csv` |
| tab:ranking (p. 7) | app:metrics (p. 27); app:focusedresults (p. 46) | `v4/ranking/campaign_summary.csv; v4/paper_tables/ranking_history.csv` |
| tab:pairranking (p. 7) | app:training (p. 24); app:metrics (p. 27); app:results (p. 33) | `v4/paper_tables/B2_unchanged.csv; v3/tables/pair_metrics.csv; v3/tables/pair_comparisons.csv` |
| tab:B3 (p. 8) | app:b3 (p. 30); app:inputcollision (p. 45) | `v4/paper_tables/B3_visibility.csv; v4/visibility/B3_rectangle_visibility.csv` |
| fig:visibility (p. 8) | app:inputcollision (p. 45); app:focusedresults (p. 46) | `v4/visibility/input_counts.csv; v3/b3_primary/` |

Main mathematical claims are the displayed architecture and objectives (pp.2,8), not separately numbered new theorems. They resolve to A3–A5 and A9. Every main section and substantive prose claim is covered by claim_audit.md. The complete source has one appendix inclusion only.

The archive retains v3/v4 raw rows at their original relative paths inside `parent_evidence.zip`; v5 rows are in the outer `v5/` directory. Extract the parent archive to access the named CSVs. No proof terminates at those pointers.

Actual boundaries: {"main_start": 1, "main_end": 9, "statements_start": 10, "references_start": 11, "references_end": 14, "appendix_start": 15, "appendix_end": 53, "appendix_pages": 39, "total_pages": 53, "workflow_page": 2}.
