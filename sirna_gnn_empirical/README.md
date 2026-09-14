# Measured siRNA GNN: code and controlled comparisons

The latest completed study is **empirical campaign v3**, presented as *Do Graphs Learn the Chemistry? Controlled Evidence from siRNA Prediction*. The conventional chemistry-aware GNN was actually trained against published measured outcomes. The campaign tests relation-support correction crossed with message passing, ten final seeds, grouped selection, matched controls, ranking and measured chemistry differences.

| v3 row-weighted MSE | Grouped ENsiRNA | APP | Davis S7 |
|---|---:|---:|---:|
| Corrected GNN | .071175 | .162475 | .109303 |
| Corrected no-message | .071237 | .163043 | .104709 |
| Original no-message | .071237 | .162914 | .105706 |
| Chemistry tree | .065648 | .189569 | .091901 |
| Token CNN | .071749 | .164365 | .101685 |
| Training row mean | .072372 | .175983 | .148803 |

No general GNN improvement is established. On APP, GNN minus tree MSE is −.027094 [−.041964, −.010482], but GNN minus the strongest observed fitted baseline (original no-message) is −.000439 [−.001856, .001151]. GNN APP R-squared is −.001360 and prediction SD .007201. These conditional component intervals hold fitted predictions fixed; they do not include training-seed variability or certify new-study generalization. Grouped ENsiRNA remains a released-label benchmark; APP and S7 are retrospective follow-up.

B2 has 156 pairs in 12 sequence components. Pair-GNN MSE .222115 compares with matched no-message .222013, pair ridge .219297, training-fold mean .263770 and zero .280804. Improvement over the fold mean remains unresolved; sign accuracy matches the majority rule. Sequence-specific measured chemistry prediction was not established. The joint four-position bundle cannot isolate GNA.

Primary Bramsen B3 evidence contains 165 measured states and 140 four-condition rectangles on one sequence background. Fixed-checkpoint GNN interaction MSE .068335 does not improve on zero .068205. APP still has no complete rectangles. Common-reference dependence, missing endpoint replicate/covariance information and primary/released-label discrepancies remain limitations: 299 of 1,604 compared values differ by more than .05, with maximum 1.0852. Acquired S1 remains quarantined for unresolved orientation.

## Current artifacts

- [Executed scientific code and reproduction scope](v3/README.md).
- [Complete aggregate comparison/fit tables](results/v3/tables/) and [frozen protocol](results/v3/protocol.json).
- [Current workflow and all eleven result figures](presentation/v3/README.md).
- [Published-method blockers](results/v3/sources/reproduction_status.json), [source reconciliation](results/v3/sources/primary_release_discrepancies.json), and [resource accounting](results/v3/final_resource_accounting.json).
- [Current publication manifest](../docs/sirna_gnn_empirical/publication/v3/files.json).

The core campaign executed 2,312 fits (1,776 development, 536 final), with ten final seeds for stochastic methods and 1,407,967 neural optimizer updates. Exact published-method fits remain zero. B3 evaluation reused source-held-out checkpoints without fitting. No new wet-lab validation, causal mechanism, architectural novelty or five-law theorem coverage is claimed. Backend no-go and the historical reference-efficiency question are unchanged.

Run `python sirna_gnn_empirical/verify_published.py` for read-only verification of the current published files. This publication includes code, aggregate results, workflow and result figures. Full measured rows, per-observation predictions, checkpoint binaries, original command logs and complete archives stay preserved locally. Paper TeX/manuscript PDFs and private authoring material are not pushed. Reproduction limitations are explicit in the code README.

## Preserved historical v1/v2 publication

The material below describes those earlier campaigns only; its B3 status predates the primary-panel follow-up above. Historical metrics and negative results are unchanged.

### Original publication

This branch publishes the scientific scripts, frozen model-selection/configuration records, aggregate comparison tables, resource accounting and aggregate figures for the completed v1 campaign and targeted v2 follow-up. It does not contain measured observation rows, per-observation prediction files, raw source exports or trained checkpoint binaries. Those complete records remain preserved locally. Paper TeX, manuscript PDFs, caches and combined archives are excluded.

| Evaluation | GNN MSE | Comparison |
|---|---:|---:|
| Original v1 APP | .174726 | Tree .168597 |
| Corrected v2 APP | .185125 | Tree .162656 |
| Grouped ENsiRNA, equal-component | .063994 | Tree .060579 |
| Supervised APP B2 pairs | .239508 | Fold mean .237508 |
| Additional Davis S7 | .091278 | Token CNN .089326 |

No consistent GNN advantage is established. Corrected APP GNN-minus-tree MSE is +.022469 [.004601, .038227]. Supervised B2's equal-component difference from the fold mean is −.008356 [−.046418, .019650]. APP follow-up and within-APP B2 are retrospective; grouped components do not certify independent laboratories. Davis S7 is separately sourced observational evidence, not prospective or new wet-lab validation. B2 changes a four-position chemistry bundle and cannot isolate GNA. B3 has no complete measured four-condition set. The GNN is conventional and has no guarantee from the separate five-law theorem; backend no-go and the historical reference-efficiency question remain unchanged.

## Files

- [Original scientific code](v1/) and [follow-up scientific code](v2/).
- [v1 aggregate results](../runs/sirna_gnn_empirical/v1-20260913T185027Z/).
- [v2 aggregate results and configurations](../runs/sirna_gnn_empirical/v2-20260913T220847Z/).
- [Aggregate figures and architecture](presentation/figures/).
- [Published file hashes](../docs/sirna_gnn_empirical/publication/v1/files.json).

The original campaign ran 37 development/final fits. The follow-up ran 210 engine fits plus four direct pair-ridge fits; 39,306 checkpoint-recorded updates plus 20 profile updates, with 0–21 possible uncheckpointed updates from an interrupted worker. Full per-fit aggregate histories and negative comparisons are retained in the tables. No scientific work was repeated to publish this branch.

## Verification and execution scope

Run `python sirna_gnn_empirical/verify_published.py` to verify the published files without training or inference. Scripts are preserved as executed; acquisition functions and recorded package versions explain their scientific dependencies. This aggregate-only checkout cannot replay trained checkpoint inference or verify omitted row-level inputs. Do not treat missing files or completion flags as a reason to relaunch the old campaign.

For separately authorized reproduction with the required inputs, use a fresh run and the frozen grouping, seeds and configurations. v1 accepts `--run` and `--through figures`; v2 scientific stages are `audit source_checks plan factorial grouped b2 external external_uncertainty evaluate figures`, with paths defined by `active_run.txt` and `core.py`. The original default runners also mention local-only manuscript/administrative stages, which are excluded here. The full original local campaigns retain their resumable checkpoints and flags.

The workflow artwork illustrates the original activity-trained architecture. Its no-pair-loss note applies to v1; v2 code adds supported R/C and a separate supervised B2 objective. It does not imply missing biological validation or architectural novelty.
