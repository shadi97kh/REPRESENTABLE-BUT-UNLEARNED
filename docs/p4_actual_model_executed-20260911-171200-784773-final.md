# P4 — executed actual-model decision (20260911-171200-784773)

**NO-GO FOR THIS SUBMISSION ROUTE.** All requested comparisons completed on physical GPU 1 (NVIDIA TITAN RTX), UUID `GPU-1e502523-2587-176b-91ef-86cb0ea01a82`. The run performed 2,400 optimizer steps and saved trained checkpoints.

The latest request to fully run the experiments was recorded as a separate execution authorization. The historical resource ledger remains unreconciled. The earlier gated entry point and its blocked result are preserved. Before seeing these outcomes, this execution declared that it would complete all comparisons; because the capability gate failed, subsequent RNA and verifier results are diagnostics and do not satisfy the original capability-gated positive route.

| Decision question | Executed evidence |
|---|---|
| Actual architecture and nonadditivity? | Actual `CertifiedModel`, float64; no readout rescaling. Exhaustive trained-model additive-fit residual max 0.0036592. |
| Independent nonlinear capacity control? | **Fail:** training R² 0.1907, required ≥0.95. Test MSE 0.010565; pairwise ridge 0.003601. |
| RNA predictive noninferiority? | Not established: Δ Spearman -0.3966, paired group 95% interval [-0.5150, -0.2584], margin −0.02. |
| RNA candidate utility? | Fail practical criterion: macro NDCG loss +0.0813; selected activity loss +0.0540; each permitted loss ≤0.02. |
| Useful learned-versus-fixed verification gain? | False. Capability gate failed; equal-total-cost measurements below remain diagnostic. |
| Numerical certification? | `float_unverified`; exhaustive float check at 1e−10 found no violations. This is not a numerical-enclosure proof. |
| Context and chemistry? | Conditional inferred shared contiguous native insert; label-independent transcript rule. No verified construct, no chemistry-transfer test. |
| Defensible ML operator distinction? | Exact support/count/partition over a declared noncrossing family plus a nonlinear residual. Edge conditioning already appears in Hojny 2024 §3.4; this experiment earns no new operator or learned-bound claim. |
| Resources and remaining? | 269.3 s wall, 0.07478 measured CPU core-hours including waited children, 0.07480 conservatively charged GPU-hours. Historical remaining balance unknown. |

## Independent capability result

One prespecified nonlinear motif task, 96 ensemble observations / 48 independent paired groups, exhaustive valid-member expectations. Equal-marginal pairs differ in higher-order statistics and labels. The loss is applied after averaging structure predictions. One fixed seed; 200 anchor epochs, 400 residual epochs; no outcome-driven rerun.

| Arm | Train R² | Held-out MSE |
|---|---:|---:|
| additive_ridge | 0.4476 | 0.006718 |
| pairwise_ridge | 0.9247 | 0.003601 |
| actual_A | 0.1723 | 0.010752 |
| actual_A_plus_B | 0.1907 | 0.010565 |

Development selected γ=1.0. The fit failure is concrete evidence about this bounded optimization, not proof that the architecture can never learn the task.

## RNA prediction and ranking

1,130 matched conditional reporter rows; train/development/test = [711, 157, 262]; global shared-13mer component split, leakage 0. Exploratory evaluation: this dataset was already inspected in P3. Reference selected by development MSE: `guide_ridge`. Actual A and B each trained for 20 epochs; development selected γ=0.125. All labels remain unclipped; assay identity is YFP / H1299 / 48 h.

| Arm | Spearman | Δ vs reference | Macro target NDCG@5 | Selected activity |
|---|---:|---:|---:|---:|
| actual_A | 0.2539 | -0.3936 | 0.8304 | 0.7295 |
| actual_A_plus_B | 0.2508 | -0.3966 | 0.8320 | 0.7282 |
| guide_pairwise_ridge | 0.5755 | -0.0719 | 0.8764 | 0.7496 |
| guide_ridge | 0.6475 | +0.0000 | 0.9132 | 0.7821 |
| guide_window_ridge | 0.6127 | -0.0347 | 0.9069 | 0.7693 |
| summary_ridge | 0.6453 | -0.0021 | 0.9129 | 0.7821 |

Ranking averages include only target pools with at least five held-out candidates. Target counts, prevalence, P@5, regret, paired group uncertainty, target-specific comparisons and leave-one-target-out sensitivity are retained in the CSV/JSON evidence. This does not establish performance on unseen targets.

## Prior, isoform and sampling sensitivity

| Declared temperature | Max prediction change | Δ held-out Spearman vs T=1 | Mean prediction Monte Carlo SE |
|---|---:|---:|---:|
| 0.5 | 0.80442 | -0.0387 | 0.000098 |
| 1.0 | 0.00000 | 0.0000 | 0.000098 |
| 2.0 | 0.39450 | 0.0060 | 0.000096 |

The declared prior is θ=0.5/T, without thermodynamic or expression interpretation. For the frozen γ=1 residual, mean sampling SE is 0.001542 at 16 draws and 0.000783 at 64. At the selected γ, the independent 64-draw replicate changes a prediction by at most 0.001752. Alternative local windows change held-out Spearman by 0.0019. These SEs measure inhibition-fraction prediction sampling error, not biological uncertainty. Training uses eight independently drawn valid graphs per observation and averages before loss; a separate frozen-checkpoint diagnostic below quantifies conditional gradient sampling variability without further fitting. The candidate dataset records within-target ranks under each prior, alternative window and evaluation draw.

## Frozen-predictor verifier

Every method uses the same saved capability checkpoint. The original and wider sequence strata were fixed in the configuration. Already-decided cases use the skip action. Setup, policy inference, oracle calls and completed updates are timed per arm. Policy teacher generation and fitting are charged separately and amortized; primary comparison is at 100 instances. Tiny-family exhaustive evaluation is the strong applicable external comparator. CPU bounds use the frozen GPU-trained weights; no second GPU is used.

Offline policy work: 0.979086 s. Coverage on the original family, including its allocated offline charge:

| Total seconds | Combined | Fixed conditioning | Deterministic | Learned | Exhaustive |
|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.05 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 0.20 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

The JSON retains both strata, exact float optima, unavailable gaps, tail times and timeouts with denominators. Neither ensemble-mean prediction nor a model threshold bound certifies measured biological ordering. Unsupported joint RNA uncertainty uses conservative independent families.

## Provenance, reproduction and remaining blocker

[The preceding audit and correction note](p4_actual_model_decision-20260911-064446-098490.md) contains the P1–P3 capability table, claims-to-evidence corrections, row exclusions and legitimate Davis acquisition attempts. Historical records and actual model code were preserved; all source hashes captured before this run were unchanged on completion. GPU batching was checked against sequential actual-model outputs and gradients; exact prior marginals and seeded samples were checked against the original oracle/sampler. Checkpoints and each stage have SHA256 sidecars.

Commands used, from the repository root (outputs are created exclusively; rerunning an exporter with the same output prefix refuses overwrites):

```bash
python -u -m scmp.revision.p4_gpu_run
python -m scmp.revision.p4_gpu_report runs/revision/p4_gpu-20260911-171200-784773-cycle.json
```

Run record: [`p4_gpu-20260911-171200-784773-cycle.json`](../runs/revision/p4_gpu-20260911-171200-784773-cycle.json). Config: [`revision_p4_gpu_v1.yaml`](../configs/revision_p4_gpu_v1.yaml).

- [figure1_prediction_ranking](../runs/revision/p4_gpu-20260911-171200-784773-figure1_prediction_ranking.csv) — 6 rows, SHA256 sidecar.
- [figure2_coverage_total_compute](../runs/revision/p4_gpu-20260911-171200-784773-figure2_coverage_total_compute_v2.csv) — 90 rows, both strata and gaps at equal total cost, SHA256 sidecar.
- [figure3_target_candidates](../runs/revision/p4_gpu-20260911-171200-784773-figure3_target_candidates.csv) — 262 rows, SHA256 sidecar.
- [target_summary](../runs/revision/p4_gpu-20260911-171200-784773-target_summary.csv) — 18 rows, SHA256 sidecar.

Peak CUDA allocation: 126.8 MiB; memory fraction cap 15%; one CPU affinity, BLAS/OpenMP thread count 1, no training workers. The execution had an automatic 900-second CPU/wall cap. Remaining local cap is not a reconciled historical allowance; no further fitting was launched.

Davis remains unavailable: `gkaf479_supplemental_files.zip` (Supplementary Table S1), Oxford Academic NAR article → Supplementary data, destination `data/davis2025/gkaf479_supplemental_files.zip`. No challenge was bypassed. Chemistry transfer, native-versus-reporter comparisons and missing dose/time/replicate covariates remain blocked.

**Recommended next action:** pivot this submission route; preserve these results as the stopping evidence rather than launch another search for favorable seeds or hyperparameters.

## Final checkpoint and evidence checks

All 54 model/P4 contract tests passed ([record](../runs/revision/p4_tests-20260911-171750-499726.json)). Three trained checkpoints reload with identical parameter hashes and finite parameters. On actual 50-nt windows, batched outputs match sequential `CertifiedModel` outputs to 8.9e-16; exact cached marginals match the original implementation to 1e-15; reloaded checkpoints reproduce saved held-out predictions to 8.9e-16.

The post-run gradient diagnostic fixes the final weights, uses the first 16 training rows and four independent sample replicates, and measures the residual at its training scale γ=1. No optimizer runs and the saved state is unchanged. Relative gradient RMS standard deviation is 3.30% at 8 draws and 1.47% at 32 draws. This measures conditional sampling variability, not variation across optimization runs or biological experiments. [Postcheck record](../runs/revision/p4_gpu-20260911-171200-784773-postcheck.json).

Top-five shortlist membership changes across 16 eligible held-out target pools:

| Perturbation | Target shortlists changed |
|---|---:|
| prior_T_0.5 | 9/16 |
| prior_T_2.0 | 8/16 |
| alternative_isoform | 0/16 |
| independent_64_draw | 0/16 |
| draw_16 | 0/16 |

The original verifier stratum is saturated once setup completes: every method has coverage 1 at 0.05 and 0.20 seconds. This is a saturation finding, not evidence for a learned bounding layer. The prespecified wider stratum is retained below and in the expanded figure dataset; its predictive usefulness was never established. Both strata therefore remain diagnostic.

| Wider stratum total seconds | Combined | Fixed conditioning | Deterministic | Learned | Exhaustive |
|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.05 | 0.750 | 0.750 | 0.750 | 0.750 | 1.000 |
| 0.2 | 0.750 | 0.750 | 0.750 | 0.750 | 1.000 |

Instrumented execution and validation totals before this final export: 276.37 seconds, 0.07673 CPU core-hours and 0.07604 conservatively charged GPU-hours. These measurements cover the recorded computation regions; interpreter/module imports and uninstrumented editor/read-only shell work are excluded. They do not reconcile the historical allowance. The experiment itself used 269.28 of its 900-second local automatic wall limit; no further training is authorized by this arithmetic.

Additional implemented commands used for these checks:

```bash
python -m scmp.revision.p4_checks
python -m scmp.revision.p4_gpu_postcheck runs/revision/p4_gpu-20260911-171200-784773-cycle.json
python -m scmp.revision.p4_gpu_finalize runs/revision/p4_gpu-20260911-171200-784773-cycle.json runs/revision/p4_tests-20260911-171750-499726.json
```
