# Actual pilot results

**Decision: incremental-known.** All 24 frozen datasets and all ten methods completed, yielding 240 paired prediction rows. The adaptive method reduced aggregate MSE by 3.5569% against the strongest complete control, below the preregistered 10% screening threshold. Three replicates per cell support a descriptive comparison only.

The proposed learned nonlinear theorem remains open. A small observed gain does not establish a new statistical method, calibrated uncertainty or a general computational advantage. No additional fits were run after final scoring.

## All methods

| Method | Aggregate MSE | Complete datasets | Mean CPU seconds, full shared cost |
|---|---:|---:|---:|
| neural_adaptive | 0.479393828 | 24/24 | 5.288478 |
| neural_rich | 0.497074021 | 24/24 | 5.275769 |
| neural_fixed | 0.498944793 | 24/24 | 5.273924 |
| neural_random | 0.499244143 | 24/24 | 5.288650 |
| neural_shared_plugin | 0.532588166 | 24/24 | 4.463973 |
| basis_rich | 0.566073129 | 24/24 | 5.690218 |
| basis_shared_plugin | 0.572467807 | 24/24 | 4.749360 |
| basis_pooled_plugin | 0.893162587 | 24/24 | 1.676210 |
| basis_direct_set_midpoint | 0.894565473 | 24/24 | 13.543022 |
| neural_pooled_plugin | 1.197264964 | 24/24 | 1.565270 |

Adaptive had lower squared error than `neural_rich` on 13/24 individual datasets. Against the same-fit fixed correction its aggregate reduction was 3.9185%. These counts are paired descriptive outcomes, without a significance or rate claim.

The direct-target set control found feasible endpoints on all 24 datasets. Its endpoint search is still approximate and restricted to the basis sieve: it is neither an optimal original-class comparator nor a certified confidence interval. No control received true profiles, coefficients, latent observations or joint outcomes.

## Cell-level comparison

| Family | delta | n per stratum | Adaptive MSE | Rich neural control MSE | Fixed neural MSE | Replicates |
|---|---:|---:|---:|---:|---:|---:|
| oscillatory | 0.1 | 256 | 0.55510642 | 0.55564192 | 0.61753693 | 3 |
| oscillatory | 0.1 | 1024 | 0.05547492 | 0.06509168 | 0.04103821 | 3 |
| oscillatory | 0.03 | 256 | 0.72896830 | 0.72608415 | 0.68666545 | 3 |
| oscillatory | 0.03 | 1024 | 0.03286693 | 0.03496542 | 0.03255675 | 3 |
| localized | 0.1 | 256 | 1.56393167 | 1.69733719 | 1.69305030 | 3 |
| localized | 0.1 | 1024 | 0.03716446 | 0.03484183 | 0.03481058 | 3 |
| localized | 0.03 | 256 | 0.81677300 | 0.81794932 | 0.82443979 | 3 |
| localized | 0.03 | 1024 | 0.04486493 | 0.04468067 | 0.06146033 | 3 |

## Actual diagnostics

| Diagnostic | Value |
|---|---:|
| rotation_models | 144 |
| selected_corrections | 360 |
| max_score_gram_refinement_difference | 0.0010506066309624639 |
| max_target_gradient_refinement_difference | 1.1456106064767517e-08 |
| max_fitted_score_center_abs | 1.5847774228128836e-05 |
| max_source_crps_refinement_difference | 0.013071797965200638 |
| max_selected_solver_relative_residual | 1.2345301667400524e-15 |
| max_selected_regularized_condition | 101.09694066775846 |
| finite_paths | 144 |
| source_compatible_finite_paths | 124 |
| max_abs_finite_path_B | 0.10207009503408186 |
| max_abs_nonlinear_remainder | 0.0012669990194164732 |
| candidate_appends_passing_precheck (not guaranteed rank gains) | 288 |
| witness_after_to_before_median | 0.9985769848973449 |
| max_truth_target_refinement_difference | 2.4513724383723456e-13 |

Adaptive lambda counts across the 72 rotations: `{'0.1': 68, '0.001': 4}`; selected dimensions: `{'12': 37, '11': 32, '10': 3}`.

The weighted Sobolev metric refinement difference is 1.368362873580413 in matrix spectral norm (4096 vs 2048 nodes). A small linear solver residual is not an integral-error bound, and a finite witness is not a full-class residual bound. The near-diagonal witness plot shows that enrichment often changed its own stabilized residual indicator very little. The finite-path remainder is retained rather than assumed zero.

## Figures and raw artifacts

All figures were rendered from the saved pilot outputs and visually inspected. Each is available as PNG and standalone SVG; the underlying data and raw-file hash are saved.

1. [Interaction error by family, separation and sample size](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/figure1_interaction_error.png) — means and every replicate; [SVG](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/figure1_interaction_error.svg).
2. [Same-fit ablations and CPU costs](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/figure2_ablation_cost.png) — shared work fully charged to each comparator; [SVG](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/figure2_ablation_cost.svg).
3. [Residual witnesses and finite-path nonlinearity](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/figure3_witness_nonlinearity.png) — lower witnesses, not uncertainty upper bounds; [SVG](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/figure3_witness_nonlinearity.svg).

[All 240 raw paired rows](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/raw_results.csv), [machine-readable summary](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/summary.json), [diagnostics](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/diagnostics.json), [plotted numerical data](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/figure_data.json), [prediction hashes](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/evaluation/prediction_hashes.json). Per-dataset source folds, model parameters, source losses, selections, directions and endpoint traces are in `runs/interaction_recoverability/learned-target-v1-20260912-164000/predictions/dNNN/`.

## Failures, costs and limits

No development/final fit or direct-target prediction failed. The intentional one-second watchdog test returned 124 and its spawned child was confirmed absent. The initial discovery mistake appended historical job 27 before the new prompt prohibition was read; this was disclosed, preserved and also charged to the fresh allowance. Default filesystem image viewing failed; budgeted reads of the existing images succeeded. Scalar-logging warnings remain in validation/run records.

Final fitter totals (24 datasets, all methods/shared work counted once) were 680.456328 CPU seconds and 753.099226 wall seconds. These omit wrapper/discovery/development/reporting costs, which the full ledger includes. Summing the table's fully charged method costs would deliberately count shared fits repeatedly.

The implementation computes the complete 36-dimensional Gram bank before enrichment, including for the shared fixed control. Thus it has no demonstrated scalable saving from avoiding that bank; the fixed method also pays for geometry it would not strictly need in a separately optimized implementation. Neural profile misspecification, source-objective quadrature, heuristic lambda selection, restricted direction coverage and incomplete endpoint optimization limit interpretation. More seeds cannot resolve the present missing theorem or prove the biological assumptions.

See [the proof gaps](theorem.md), [source comparison](novelty.md), [static review](review.md), [resource accounting](resources.md) and [actual commands](commands.sh).

Static review found a mismatch between append nonredundancy and solver rank tolerances. Some appended directions were subsequently treated as redundant: 37 rotations ended at dimension 12, 32 at 11 and 3 at 10. This limits the effective enrichment mechanism; no frozen code or predictions were changed.
