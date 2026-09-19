# Measured siRNA GNN: code and aggregate results

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
