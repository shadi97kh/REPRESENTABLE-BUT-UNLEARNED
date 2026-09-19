# Completed empirical revision and handoff

The finite core campaign completed: **2,312 fits** (1,776 development; 536 final), ten final seeds for stochastic methods, **1,407,967 executed optimizer updates** and **406,367 selected-checkpoint updates**. No core fit failed; no unfavorable finite seed was discarded. The published-method objective remains blocked: zero exact ENsiRNA-mod/MEG-mod fits. S1 remains acquired but quarantined for unresolved orientation. The separate measured primary B3 panel was evaluated using source-held-out checkpoints, with zero additional fits.

No GNN advantage over the strongest observed fitted APP baseline and no sequence-specific measured B2 chemistry-prediction improvement were established. APP/S7 are retrospective follow-up. ENsiRNA is a released-label benchmark with unresolved primary/released-value discrepancies. This is not new wet-lab, mechanistic, causal or clinical validation. Backend no-go and the historical reference-efficiency question remain unchanged.

## Verified deliverables

- Paper: `papers/interaction_recoverability_iclr2027/v6/main.pdf`.
- Editable independent-build source: `papers/interaction_recoverability_iclr2027/v6/manuscript_source.zip`.
- Caption-free PDF/SVG/PNG and separate caption/source: `papers/interaction_recoverability_iclr2027/v6/workflow_handoff.zip` and `figures/workflow.*` within that manuscript directory.
- Scientific package: `runs/sirna_gnn_empirical/v3-20260914T013314Z/paper_artifacts/code_results.zip` (6,079,619,157 bytes; all 13,667 included files stream-hash verified; every current fit artifact included).
- Compact review packet: `runs/sirna_gnn_empirical/v3-20260914T013314Z/paper_artifacts/review_packet.zip` (5,248,138 bytes); its review request is a draft, not sent.
- Full results, review response and claim/citation/appendix/build/visual/integrity audits: `docs/sirna_gnn_empirical/v3/`.
- Full activity seed and ensemble predictions: current run `evaluated/activity_seed_predictions.csv` and `evaluated/activity_ensemble_predictions.csv`; B2 `evaluated/pair_ensemble_predictions.csv`; B3 `b3_primary/interaction_predictions.csv`.
- Every selected/last checkpoint, optimizer state and history: current run `fits/activity/` and `fits/pair/`.
- Final commands and post-package costs: this report, `executed_commands.md`, original current-run `commands.jsonl`, and `final_resource_accounting.json`.

Main pages **1–9**; statements **10**; references **11–14**; unrestricted appendix **15–718 (704 pages)**. The workflow is on page **2**. There are seven main tables, three actual main results figures plus the workflow, eight appendix figures and 43 verified bibliography entries. Figure/table captions are below their objects; standalone workflow exports contain no embedded caption. All 718 pages were inspected. All-page text and rendered content matched an independent source-ZIP compile. No unresolved cross-reference or overflow warning remains.

All **7,444 preexisting files** retain their hashes. All 13 completed pipeline stages verified with matching code, dependency and output hashes; that verification ran zero fits and zero updates. Current scientific check record reports 8,568 checks and 10,800 hashed fit artifacts. The preflight and float32/float64 audit failures and all authoring/source failures are retained. Correcting the latter audit changed no prediction or scientific metric.

## Final executed administrative commands

`bash sirna_gnn_empirical/v3/run_all.sh deliver` was attempted twice: the first stopped at the arithmetic audit; the second completed all checks, independent compilation, preservation validation and archives. Its exact child command was `/opt/miniforge3/bin/python /home/shadi/iclr2027/sirna_gnn_empirical/v3/deliver.py`.

`bash sirna_gnn_empirical/v3/run_all.sh` then verified all 13 completed stages, with the scientific command record unchanged. This was a hash verification, not training.

| Measurement | Seconds | Scope |
|---|---:|---|
| Actual fit CPU | 12523.769 | All 2,312 fits |
| Actual fit wall | 12360.590 | Aggregate fitting regions |
| GPU fitting-region wall | 11689.879 | Overlaps fit/stage wall; not kernel time or energy |
| All recorded stage child CPU | 13578.145 | Includes failed commands and successful delivery |
| All recorded stage aggregate wall | 13399.053 | Includes failed commands and successful delivery |
| Successful delivery child CPU | 346.155 | Already included in stage total |
| Successful delivery wall | 345.854 | Already included in stage total |
| Post-delivery verification child CPU | 28.825 | Additional administrative hash verification |
| Post-delivery verification wall | 27.709 | Additional administrative hash verification |

Maximum recorded stage child RSS is 3229.336 MiB. Additional recorded authoring and source costs retain their separate scopes in the JSON report. Direct inspection, browsing, authoring, rendering, copying, hashing and receipt generation incurred additional nonzero administrative cost; these meters do not claim complete machine or energy accounting.

The source/scientific/review ZIPs are immutable snapshots before their own completion records; avoiding recursive self-hashes is intentional. `completion_receipt.zip` supplies post-package records, the current command log and completion hashes separately. It does not alter the scientific artifacts or original archive manifests. Unused representation assets and downloaded citation full texts are explicitly omitted from the focused scientific archive with hashes and acquisition instructions; all remain in the original local run.

No repository push, publication or external contact was performed, consistent with the current revision brief. No new research or fitting cycle is pending.
