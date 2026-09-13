# Executed chemistry-aware siRNA GNN campaigns

This contribution contains the scientific implementation and existing measured-data results from the original v1 campaign and its targeted v2 follow-up. It preserves negative outcomes. The GNN is conventional relation-specific message passing; it is neither a new architecture nor an estimator covered by the separate five-law theorem. APP follow-up is retrospective, grouped ENsiRNA components do not certify laboratory independence, and Davis S7 is qualified observational source evidence. No new wet-lab validation was performed.

| Evaluation | GNN MSE | Strong comparison |
|---|---:|---:|
| Original v1 APP | .174726 | Tree .168597 |
| Corrected v2 APP | .185125 | Tree .162656 |
| Grouped ENsiRNA, equal-component | .063994 | Tree .060579 |
| Supervised APP B2 pair effects | .239508 | Fold mean .237508 |
| Additional Davis S7 | .091278 | Token CNN .089326 |

No consistent GNN superiority is established. Corrected APP GNN-minus-tree MSE is +.022469 [.004601, .038227]. Supervised B2 does not clearly beat the training-fold constant effect: equal-component difference −.008356 [−.046418, .019650]. B2 jointly changes guide positions 6–9 and cannot isolate GNA. B3 has no complete measured four-condition set. S1 orientation and assay/replicate/control-covariance gaps remain explicit. Backend no-go and the historical reference-efficiency question remain unchanged.

## Contents

- [v1 scientific code](v1/) and [v2 scientific code](v2/): acquisition/reconciliation, graph construction, training, controls, grouped evaluation and measured-effect analysis.
- [Original results](../runs/sirna_gnn_empirical/v1-20260913T185027Z/): 4,766 admitted measured rows, split/graph records, 37 development/final fits, held-out predictions and comparisons.
- [Follow-up results](../runs/sirna_gnn_empirical/v2-20260913T220847Z/): support diagnostics, all 2×2 refits and controls, five grouped folds, exploratory four-fold B2, 118 S7 outcomes, all fit histories and predictions. There are 210 engine fits plus four pair-ridge fits; 39,306 checkpoint-recorded updates plus a 20-update profile, with 0–21 possible additional uncheckpointed interrupted updates.
- [Current result figures and plotting code](presentation/).
- [Publication inventory](../docs/sirna_gnn_empirical/publication/v1/files.json) and [checkpoint inventory](../docs/sirna_gnn_empirical/publication/v1/checkpoint_inventory.json).

The contribution includes **154 selected final/deployment checkpoint files**: best neural states and fitted classical models across final seeds, folds and deployment fits. All development/final fit metadata and saved predictions are included. Last-epoch states, development checkpoint binaries, complete local superseded records and duplicate full archives remain preserved locally; their checkpoint identities are inventoried. This is a focused Git contribution, not the entire local-run archive.

## Verify and reproduce

From the repository root, `python sirna_gnn_empirical/verify_published.py` verifies every published file without fitting, inference, sampling or resource-ledger writes. Exact originally executed commands and costs are preserved in `docs/sirna_gnn_empirical/v1/` and `v2/`. Package versions are recorded in `v1/requirements-recorded.txt` and the run environment records.

The original runners are preserved verbatim, including their original default administrative stages and local interpreter paths. Manuscript-building scripts are deliberately not published. The Git checkout omits campaign completion flags because those flags also hash local-only files; use the publication verifier for this checkout. Do not interpret a missing flag as permission to rerun a completed study or create a new blind test. The complete original local campaign remains resumable using its retained flags and last checkpoints.

For an intentional fresh reproduction, copy inputs to a new run, preserve the frozen splits/configurations and use compatible software plus actual available resources. v1 supports `--run` and `--through figures`; v2 resolves its run from `active_run.txt` and its read-only v1 dependency from `core.py`. Scientific v2 stage order is `audit source_checks plan factorial grouped b2 external external_uncertainty evaluate figures`. Source-check/acquisition stages may require re-fetching primary HTML/PDF files identified by the saved acquisition manifests: machine-readable workbooks/TSV, clean rows and graph tensors are included, while redundant article/patent HTML/PDF copies are not. External pretrained feature backends are not required.

The paper's TeX, bibliography/style files, manuscript PDFs, private authoring records, dependency caches and giant combined archives are excluded. This GitHub publication is authorized by the user's subsequent explicit push request; earlier no-publication statements remain accurate descriptions of those completed campaigns, not a current prohibition.
