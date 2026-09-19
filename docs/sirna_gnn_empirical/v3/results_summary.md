# Executed revision results

An APP GNN improvement over the strongest observed fitted baseline is not established. This is retrospective follow-up in one related patent family, not independent external confirmation. Sequence-specific measured chemistry prediction is not established under the declared B2 mean-control and within-assay diagnostic.

| cohort | strongest_observed_baseline | baseline_mse | gnn_mse | gnn_minus_baseline | conditional_lower | conditional_upper |
|---|---|---|---|---|---|---|
| APP | original_no_message | 0.162914 | 0.162475 | -0.000439 | -0.001856 | 0.001151 |
| Davis_S7 | chemistry_tree | 0.091901 | 0.109303 | 0.017402 | -0.008409 | 0.038987 |
| ENsiRNA_grouped | chemistry_tree | 0.065648 | 0.071175 | 0.005527 | -0.000871 | 0.013815 |

## Main activity results

| cohort | method | n | mse | r_squared | correlation | prediction_sd | mean_bias |
|---|---|---|---|---|---|---|---|
| APP | chemistry_tree | 1839 | 0.189569 | -0.168344 | -0.020966 | 0.010342 | -0.164417 |
| APP | corrected_gnn | 1839 | 0.162475 | -0.001360 | -0.017807 | 0.007201 | -0.008095 |
| APP | corrected_no_message | 1839 | 0.163043 | -0.004862 | 0.014564 | 0.008474 | -0.028573 |
| APP | guide_ridge | 1839 | 0.294937 | -0.817739 | -0.014124 | 0.142048 | 0.332998 |
| APP | original_gnn | 1839 | 0.162483 | -0.001409 | -0.025697 | 0.007318 | -0.004853 |
| APP | original_no_message | 1839 | 0.162914 | -0.004064 | -0.001225 | 0.008517 | -0.024051 |
| APP | pairwise_ridge | 1839 | 0.206430 | -0.272260 | 0.028829 | 0.061911 | 0.204402 |
| APP | token_cnn | 1839 | 0.164365 | -0.013004 | 0.079690 | 0.003631 | -0.048269 |
| APP | training_equal_group_mean | 1839 | 0.176820 | -0.089771 | -- | 0.000000 | 0.120689 |
| APP | training_row_mean | 1839 | 0.175983 | -0.084608 | -- | 0.000000 | 0.117167 |
| APP | oracle_test_mean_diagnostic | 1839 | 0.162255 | 0.000000 | -- | 0.000000 | 0.000000 |
| Davis_S7 | chemistry_tree | 118 | 0.091901 | -0.021458 | 0.016866 | 0.008864 | -0.044064 |
| Davis_S7 | corrected_gnn | 118 | 0.109303 | -0.214875 | 0.111815 | 0.005857 | 0.140324 |
| Davis_S7 | corrected_no_message | 118 | 0.104709 | -0.163817 | 0.117478 | 0.006577 | 0.123121 |
| Davis_S7 | guide_ridge | 118 | 0.781455 | -7.685695 | 0.123230 | 0.115004 | 0.828710 |
| Davis_S7 | original_gnn | 118 | 0.110243 | -0.225331 | 0.115673 | 0.005832 | 0.143679 |
| Davis_S7 | original_no_message | 118 | 0.105706 | -0.174904 | 0.110080 | 0.006686 | 0.127016 |
| Davis_S7 | pairwise_ridge | 118 | 0.969840 | -9.779558 | 0.031355 | 0.061935 | 0.936589 |
| Davis_S7 | token_cnn | 118 | 0.101685 | -0.130206 | 0.107913 | 0.003387 | 0.109190 |
| Davis_S7 | training_equal_group_mean | 118 | 0.150524 | -0.673042 | -- | 0.000000 | 0.246077 |
| Davis_S7 | training_row_mean | 118 | 0.148803 | -0.653915 | -- | 0.000000 | 0.242555 |
| Davis_S7 | oracle_test_mean_diagnostic | 118 | 0.089970 | 0.000000 | -- | 0.000000 | -0.000000 |
| ENsiRNA_grouped | chemistry_tree | 2927 | 0.065648 | 0.087464 | 0.311548 | 0.075597 | 0.025040 |
| ENsiRNA_grouped | corrected_gnn | 2927 | 0.071175 | 0.010638 | 0.162503 | 0.057437 | 0.030702 |
| ENsiRNA_grouped | corrected_no_message | 2927 | 0.071237 | 0.009773 | 0.157750 | 0.057493 | 0.029269 |
| ENsiRNA_grouped | guide_ridge | 2927 | 0.093788 | -0.303696 | 0.055822 | 0.117213 | -0.107791 |
| ENsiRNA_grouped | original_gnn | 2927 | 0.071175 | 0.010637 | 0.162502 | 0.057437 | 0.030702 |
| ENsiRNA_grouped | original_no_message | 2927 | 0.071237 | 0.009772 | 0.157748 | 0.057492 | 0.029270 |
| ENsiRNA_grouped | pairwise_ridge | 2927 | 0.073690 | -0.024317 | 0.058036 | 0.053547 | -0.023434 |
| ENsiRNA_grouped | token_cnn | 2927 | 0.071749 | 0.002653 | 0.143498 | 0.067982 | 0.020509 |
| ENsiRNA_grouped | training_equal_group_mean | 2927 | 0.076897 | -0.068907 | -0.094846 | 0.024997 | 0.055322 |
| ENsiRNA_grouped | training_row_mean | 2927 | 0.072372 | -0.006005 | -0.072568 | 0.008846 | 0.003065 |
| ENsiRNA_grouped | oracle_test_mean_diagnostic | 2927 | 0.071940 | 0.000000 | -- | 0.000000 | -0.000000 |

## Direct conditional comparisons

| cohort | a | b | difference | lower | upper |
|---|---|---|---|---|---|
| APP | corrected_gnn | corrected_no_message | -0.000568 | -0.002408 | 0.001521 |
| APP | corrected_gnn | chemistry_tree | -0.027094 | -0.041964 | -0.010482 |
| APP | corrected_gnn | token_cnn | -0.001889 | -0.005689 | 0.002437 |
| APP | corrected_gnn | training_equal_group_mean | -0.014345 | -0.027798 | -0.002453 |
| APP | corrected_gnn | training_row_mean | -0.013507 | -0.026577 | -0.001949 |
| Davis_S7 | corrected_gnn | corrected_no_message | 0.004594 | 0.002153 | 0.006631 |
| Davis_S7 | corrected_gnn | chemistry_tree | 0.017402 | -0.008409 | 0.038987 |
| Davis_S7 | corrected_gnn | token_cnn | 0.007618 | 0.003364 | 0.011162 |
| Davis_S7 | corrected_gnn | training_equal_group_mean | -0.041221 | -0.054144 | -0.025795 |
| Davis_S7 | corrected_gnn | training_row_mean | -0.039501 | -0.052014 | -0.024579 |
| ENsiRNA_grouped | corrected_gnn | corrected_no_message | -0.000062 | -0.000533 | 0.000407 |
| ENsiRNA_grouped | corrected_gnn | chemistry_tree | 0.005527 | -0.000871 | 0.013815 |
| ENsiRNA_grouped | corrected_gnn | token_cnn | -0.000574 | -0.004581 | 0.000319 |
| ENsiRNA_grouped | corrected_gnn | training_equal_group_mean | -0.005722 | -0.009112 | 0.009799 |
| ENsiRNA_grouped | corrected_gnn | training_row_mean | -0.001197 | -0.004103 | 0.011824 |

Ten thousand whole-sequence-component resamples with fixed ensemble predictions; unadjusted descriptive intervals. No new-study or retraining coverage is claimed.

## Seed variability

| cohort | method | count | mean_seed_mse | sd_seed_mse | min_seed_mse | max_seed_mse |
|---|---|---|---|---|---|---|
| APP | chemistry_tree | 10 | 0.189673 | 0.002828 | 0.184125 | 0.193860 |
| APP | corrected_gnn | 10 | 0.165622 | 0.003448 | 0.161411 | 0.173361 |
| APP | corrected_no_message | 10 | 0.167563 | 0.006453 | 0.160956 | 0.181784 |
| APP | guide_ridge | 1 | 0.294937 | -- | 0.294937 | 0.294937 |
| APP | original_gnn | 10 | 0.165456 | 0.002581 | 0.161620 | 0.169508 |
| APP | original_no_message | 10 | 0.167055 | 0.005210 | 0.160986 | 0.178658 |
| APP | pairwise_ridge | 1 | 0.206430 | -- | 0.206430 | 0.206430 |
| APP | token_cnn | 10 | 0.168071 | 0.004881 | 0.162278 | 0.179572 |
| Davis_S7 | chemistry_tree | 10 | 0.092024 | 0.000820 | 0.090638 | 0.093698 |
| Davis_S7 | corrected_gnn | 10 | 0.112109 | 0.015328 | 0.092031 | 0.137607 |
| Davis_S7 | corrected_no_message | 10 | 0.107661 | 0.014382 | 0.090965 | 0.138759 |
| Davis_S7 | guide_ridge | 1 | 0.781455 | -- | 0.781455 | 0.781455 |
| Davis_S7 | original_gnn | 10 | 0.112777 | 0.014868 | 0.092380 | 0.137939 |
| Davis_S7 | original_no_message | 10 | 0.108309 | 0.014170 | 0.092103 | 0.139849 |
| Davis_S7 | pairwise_ridge | 1 | 0.969840 | -- | 0.969840 | 0.969840 |
| Davis_S7 | token_cnn | 10 | 0.104893 | 0.014973 | 0.090214 | 0.133885 |
| ENsiRNA_grouped | chemistry_tree | 10 | 0.065776 | 0.000451 | 0.065211 | 0.066828 |
| ENsiRNA_grouped | corrected_gnn | 10 | 0.072370 | 0.002212 | 0.068878 | 0.077027 |
| ENsiRNA_grouped | corrected_no_message | 10 | 0.072582 | 0.003921 | 0.066783 | 0.081365 |
| ENsiRNA_grouped | guide_ridge | 1 | 0.093788 | -- | 0.093788 | 0.093788 |
| ENsiRNA_grouped | original_gnn | 10 | 0.072370 | 0.002212 | 0.068878 | 0.077027 |
| ENsiRNA_grouped | original_no_message | 10 | 0.072582 | 0.003922 | 0.066783 | 0.081365 |
| ENsiRNA_grouped | pairwise_ridge | 1 | 0.073690 | -- | 0.073690 | 0.073690 |
| ENsiRNA_grouped | token_cnn | 10 | 0.073035 | 0.001716 | 0.070883 | 0.075765 |

Member risk is separate from error of averaged predictions. Poor finite seeds remain.

## Measured B2 and ranking

| method | n | mse | correlation | prediction_sd | label_sd | mean_bias | non_tie_correct | non_ties |
|---|---|---|---|---|---|---|---|---|
| activity_gnn | 156 | 0.235088 | 0.083055 | 0.036029 | 0.478946 | 0.085247 | 111 | 137 |
| pair_gnn | 156 | 0.222115 | 0.247522 | 0.038096 | 0.478946 | 0.017511 | 111 | 137 |
| pair_no_message | 156 | 0.222013 | 0.227998 | 0.043462 | 0.478946 | 0.015054 | 111 | 137 |
| pair_ridge | 156 | 0.219297 | 0.225689 | 0.147529 | 0.478946 | -0.006030 | 105 | 137 |
| pair_tree | 156 | 0.243772 | 0.119857 | 0.167341 | 0.478946 | -0.074779 | 109 | 137 |
| training_majority | 156 | -- | -- | -- | 0.478946 | -- | 111 | 137 |
| training_mean | 156 | 0.263770 | -0.075767 | 0.063754 | 0.478946 | -0.160279 | 111 | 137 |
| zero | 156 | 0.280804 | -- | 0.000000 | 0.478946 | 0.226749 | 0 | 137 |

| a | b | difference | lower | upper | groups | weighting | bootstrap_resamples | analysis_seed | uncertainty |
|---|---|---|---|---|---|---|---|---|---|
| pair_gnn | training_mean | -0.041656 | -0.089696 | 0.008949 | 12 | rows | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | training_mean | -0.012472 | -0.070685 | 0.046408 | 12 | equal_groups | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | activity_gnn | -0.012973 | -0.046750 | 0.004131 | 12 | rows | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | activity_gnn | -0.034897 | -0.081357 | -0.002456 | 12 | equal_groups | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | pair_no_message | 0.000102 | -0.001121 | 0.002286 | 12 | rows | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | pair_no_message | 0.001553 | -0.000518 | 0.004031 | 12 | equal_groups | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | pair_ridge | 0.002818 | -0.025859 | 0.055281 | 12 | rows | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | pair_ridge | 0.032289 | -0.021360 | 0.108081 | 12 | equal_groups | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | pair_tree | -0.021657 | -0.064189 | 0.035043 | 12 | rows | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | pair_tree | 0.011074 | -0.043636 | 0.075713 | 12 | equal_groups | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | zero | -0.058690 | -0.136692 | -0.016897 | 12 | rows | 10000 | 90413 | conditional on fitted predictions and observed cohort |
| pair_gnn | zero | -0.113041 | -0.213465 | -0.034949 | 12 | equal_groups | 10000 | 90413 | conditional on fitted predictions and observed cohort |

| method | correlation | prediction_sd | conditional_correlation_lower | conditional_correlation_upper |
|---|---|---|---|---|
| activity_gnn | 0.089873 | 0.036029 | -0.162206 | 0.325866 |
| pair_gnn | 0.265777 | 0.038095 | -0.045463 | 0.493471 |
| pair_no_message | 0.244116 | 0.043460 | -0.038110 | 0.476863 |
| pair_ridge | 0.243951 | 0.147529 | -0.045963 | 0.530990 |
| pair_tree | 0.129555 | 0.167341 | -0.120674 | 0.416465 |
| training_mean | -0.081898 | 0.063754 | -0.432608 | 0.256728 |
| zero | -- | 0.000000 | -- | -- |

| method | pools | top5_mean_percentile | selected_mean_activity |
|---|---|---|---|
| chemistry_tree | 17 | 0.444754 | 0.495961 |
| corrected_gnn | 17 | 0.533871 | 0.555548 |
| corrected_no_message | 17 | 0.576950 | 0.589464 |
| guide_ridge | 17 | 0.434403 | 0.455574 |
| original_gnn | 17 | 0.522677 | 0.555162 |
| original_no_message | 17 | 0.587979 | 0.594758 |
| pairwise_ridge | 17 | 0.478669 | 0.498115 |
| token_cnn | 17 | 0.582679 | 0.607422 |
| training_equal_group_mean | 17 | 0.500000 | 0.524349 |
| training_row_mean | 17 | 0.500000 | 0.524349 |

The B2 bundle jointly changes guide positions 6/8/9 from 2-Fluoro to 2-O-Methyl and position 7 from 2-O-Methyl to GNA. This cannot isolate GNA. Majority is a sign-only rule, not a regression magnitude. Endpoint SDs do not provide mean-difference SEs without replicate counts and shared-control covariance. Within-assay centering is a labeled analysis diagnostic, never a deployed calibration. The APP pool scores are descriptive and share a patent family/sequence backgrounds.

## Primary measured B3

| method | seed | mse | correlation | prediction_sd | label_sd | mean_bias |
|---|---|---|---|---|---|---|
| chemistry_tree | ensemble | 0.071939 | -0.177024 | 0.017767 | 0.222442 | 0.144025 |
| corrected_gnn | ensemble | 0.068335 | -0.090263 | 0.001589 | 0.222442 | 0.137070 |
| corrected_no_message | ensemble | 0.068281 | -0.034429 | 0.001652 | 0.222442 | 0.137013 |
| guide_ridge | ensemble | 0.068205 | -0.204282 | 0.000000 | 0.222442 | 0.136836 |
| original_gnn | ensemble | 0.068335 | -0.090263 | 0.001589 | 0.222442 | 0.137070 |
| original_no_message | ensemble | 0.068281 | -0.034429 | 0.001652 | 0.222442 | 0.137013 |
| pairwise_ridge | ensemble | 0.068205 | 0.012370 | 0.000000 | 0.222442 | 0.136836 |
| token_cnn | ensemble | 0.068409 | -0.046847 | 0.002879 | 0.222442 | 0.137331 |
| zero_interaction | deterministic | 0.068205 | -- | 0.000000 | 0.222442 | 0.136836 |

The separate Bramsen primary panel supplies 165 measured conditions and 140 rectangles, one background and common reference. Only source/13mer-held-out outer0 checkpoints are reused; zero extra fits. APP still has zero complete rectangles. No synthetic labels, new wet-lab validation, private-anchor calibration, causal mechanism or clinical benefit was produced.

## Execution and costs

| actual_core_fits | development_fits | final_fits | neural_fits | failed_core_fits | executed_updates | selected_updates | fit_cpu_s | fit_wall_s | gpu_region_wall_s | completed_stage_child_cpu_s | completed_stage_aggregate_wall_s | max_stage_child_rss_kib | failed_stage_commands | published_exact_fits | unused_published_fit_reserve | recorded_authoring_attempts | recorded_authoring_wall_s | recorded_authoring_child_cpu_s | failed_authoring_attempts | separately_recorded_source_administration_wall_s | separately_recorded_source_command_child_cpu_s | source_helper_recorded_errors |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2312 | 1776 | 536 | 1932 | 0 | 1407967 | 406367 | 12523.768947 | 12360.589618 | 11689.878588 | 13231.989987 | 13053.199304 | 3306840 | 2 | 0 | 100 | 18 | 246.169018 | 247.134205 | 4 | 41.215556 | 6.848592 | 36 |

Snapshot before final packaging. Fit CPU/wall/GPU region overlap stage costs; do not add. Queued waiting overlaps activity and is not a second training cost. Later packaging cost is recorded separately. Nonzero inspection, authoring, browsing, preflight capture/rendering/copying/hashing and parent-runner work; not included in child fit meters.

## Remaining gaps

Exact ENsiRNA-mod reproduction has missing official geometry and absent configured Rosetta; exact MEG-mod reproduction has an incomplete released modification-token forward path. Both were attempted, and zero exact published predictors were fitted. S1 is acquired but remains quarantined for unresolved strand/orientation representation. Primary/released Bramsen reconciliation finds 299 differences greater than .05 among 1,604 identity-matched values; frozen ENsiRNA prediction is therefore a released-label benchmark, not completely primary-verified measurement prediction. Some smaller discrepancies are rounding. Historical negative results, theory, checkpoints and ledgers remain preserved. Backend no-go and the open historical reference-efficiency question are unchanged.

## Paper and artifacts

Main pages [1, 9]; statements [10, 10]; references [11, 14]; appendix [15, 718] (704 pages). Workflow page 2.

PDF/source/workflow: `papers/interaction_recoverability_iclr2027/v6`. Predictions, all checkpoint/optimizer histories, tables, protocol and commands: `runs/sirna_gnn_empirical/v3-20260914T013314Z`. Code: `sirna_gnn_empirical/v3`. Code/results archive: `paper_artifacts/code_results.zip` under the current run. No push, publication or external contact was performed.
