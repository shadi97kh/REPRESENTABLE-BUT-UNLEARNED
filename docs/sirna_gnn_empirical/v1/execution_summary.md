4,766 admitted measured observations from 18 recorded source labels, grouped into 11 study-linked components; training 2,626, validation 103, internal test 198, APP test 1,839. The main GNN has five completed fits (two development, three final), 2,625 optimizer updates across those fits, plus a separate 20-update profile; its final checkpoints came from 1,701 executed updates. There are 34,629 per-fit/model held-out predictions and 14,259 ensemble/model predictions on 2,037 observations. B2 has 156 pairs, 26 exact duplex backgrounds and 12 sequence components; B3 has zero complete rectangles. These grouping units are not certified independent biological studies: all APP evidence is one related patent-family cohort.

The actual chemistry-aware GNN trained and was evaluated against published measured outcomes. B1 prediction and candidate-ranking evaluation ran; B2 measured chemistry-pattern comparison ran with a negative effect-size result. B3 remains incomplete for lack of four matched measured conditions. No new wet-lab validation, causal mechanism or clinical benefit was produced.

| Model | Internal MSE | APP MSE | APP mean top-five percentile | B2 pair MSE |
|---|---:|---:|---:|---:|
| Chemistry-aware GNN | .060432 | .174726 | .521442 | .292107 |
| Chemistry tree | .051757 | .168597 | .564846 | .281246 |
| Token CNN (validation-selected reference) | .056367 | .186286 | .452873 | .284191 |
| No-message matched control | .056301 | .195609 | .497116 | .283509 |
| Guide pairwise ridge | .065838 | .208019 | .434037 | .280804 |
| Guide ridge | .071228 | .301295 | .518567 | .280804 |
| GNN without chemistry | .066582 | .207189 | .654656 | .280804 |

The strongest observed APP baseline is the chemistry tree. GNN minus tree MSE is +.006129, paired conditional 95% interval [-.003419,.015526]; no GNN improvement is established. GNN minus the validation-selected token comparator is -.011560 [-.018262,-.002634]. GNN minus the equal-parameter no-message control is -.020883 [-.030134,-.009088]. These intervals resample sequence components within the observed cohort and do not establish new-study or new-laboratory generalization.

B2 jointly changes guide positions 6–9, replacing 2-Fluoro at 6/8/9 with 2-O-Methyl and 2-O-Methyl at 7 with GNA. It cannot isolate GNA's effect. GNN B2 pair MSE .292107 exceeds the chemistry-invariant zero-difference control's .280804. Cluster-weighted GNN error exceeds token CNN by .008595 [.002256,.017930]. GNN sign accuracy is 18/137 observed non-ties under the predeclared ±.02 band. Endpoint SDs are preserved, but replicate counts and shared-control covariance are unresolved; no mean-difference SE is fabricated.

The model is conventional relation-specific message passing with 153,345 trainable parameters. Actual gradients and weight changes are recorded. It establishes neither architectural novelty nor coverage by the five-law theorem. ENsiRNA/MEG-mod already provide chemistry-aware graph precedents. The previous empirical-quantile estimator, neural-profile/Riesz pilot and current GNN remain distinct. Accepted theoretical v1, proofs, frozen pilot, ledgers and external_review/v2 are preserved; backend no-go and the historical reference-efficiency question remain unchanged.

Data limitations are material: all APP rows contain a chemistry name absent from training; validation has only three sequence components; ENsiRNA assay identity and concentration units are incompletely resolved; every APP source belongs to one related family. Acquired Davis S1 is quarantined for an unresolved reported strand/orientation representation, not declared unavailable. There are no private-anchor calibrations and no four-condition B3 labels. The prior architecture/biology gaps are not filled with synthetic values.

Executed commands are listed verbatim in [executed_commands.md](executed_commands.md). The resumable entry point is `bash sirna_gnn_empirical/v1/run_all.sh`; all 18 stages completed and a subsequent invocation verified their hashes without scientific reruns. The campaign executed 20 development and 17 final fits, with 10,395 neural optimizer updates plus 20 profile updates. All final stochastic methods use seeds 1103, 2207 and 3301. No test-driven restart or tuning occurred.

The full manuscript has main pages **1–9**, excluded reproducibility/ethics/AI statements on **10**, references **11–14**, and appendix **15–42 (28 pages)**. It contains four main figure slots (workflow reserved plus three actual results), four main tables, seven appendix figures, complete fit/per-assay tables, and 44 citations (42 papers plus two corrections). All pages were rendered and inspected. The official style is unchanged, and cross-reference/overflow checks pass. The workflow was specified, not drawn. The source ZIP was independently compiled.

Artifacts:

- PDF: `papers/interaction_recoverability_iclr2027/v2/main.pdf`
- Editable source: `papers/interaction_recoverability_iclr2027/v2/manuscript_source.zip`
- Workflow handoff: `papers/interaction_recoverability_iclr2027/v2/workflow_handoff.zip`
- Scientific code/results/checkpoints/inputs: `runs/sirna_gnn_empirical/v1-20260913T185027Z/paper_artifacts/code_results.zip`
- Complete local run with original machine logs and superseded records: `runs/sirna_gnn_empirical/v1-20260913T185027Z/`
- Checkpoints: `fits/final/<model>-seed<seed>-<configuration>/best.pt` or `model.pkl`, with `last.pt`, optimizer states and fit histories for neural models.
- Audits: this directory's claim, citation, appendix, build, integrity, resource and package-verification reports.

Recorded pipeline cost, including failed attempts: **166.525 child CPU seconds** and **152.203 aggregate command wall seconds**. GPU fitting-region elapsed time is 41.476s plus 0.413s profile; these overlap stage time and are not GPU kernel/energy measurements. Maximum reported child RSS is 1561.8MiB. Direct inspection, authoring, browsing, preflight compilation, rendering and hashing incurred additional unmetered nonzero administrative cost. No historical allowance or ledger was used. All 4,382 preexisting files retain their hashes.

No repository push was performed in this campaign. The earlier push restriction remains unresolved; no publishing or external contact was attempted.
