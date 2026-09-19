# Focused source-integrity robustness campaign

Parent: empirical v3 (2,312 fits), latest verified manuscript v8. Current: empirical v4 / manuscript v9.

The frozen focused protocol executed 792 new fits (432 development, 360 final), ten final seeds for
corrected GNN, matched no-message, chemistry tree and token CNN. Forty unaffected original outer0
final checkpoints were verified and reused. B2/B3 require no new fits. Original labels/results remain
unchanged. The primary-assay target is a separately versioned sensitivity, not a certified label repair.

```bash
bash sirna_gnn_empirical/v4/run_all.sh
```

Completion receipts hash code/configurations and outputs. A completed invocation verifies and skips
stages. The fit engine resumes unfinished neural trajectories from saved optimizer/RNG states; finite
unfavorable seeds are retained. Failed numerical fits would remain flagged with a training-only fallback;
the executed focused matrix had zero such failures. There is no test-driven restart or expanded search.

The source-excluded protocol removes Bramsen from training, validation, preprocessing and selection.
Primary-assay sensitivity quarantines 320 unresolved source rows and retains 1,604 unique reported-identity
links plus 1,003 other-source rows. APP/S7 are retrospective follow-up. Input collisions use actual frozen
masked tensors and retained feature columns; B3 does not bootstrap dependent rectangles independently.

The active run path is recorded in `active_run.txt`. `docs/sirna_gnn_empirical/v4/` contains results,
commands, source/review responses and scientific/build audits. Manuscript sources are in
`papers/interaction_recoverability_iclr2027/v9/source/`, with exactly four entries. External styles,
build outputs, PDF/ZIPs and visual-inspection renders are siblings outside that source directory.

The evidence ZIP supplies original measured inputs, predictions and fit records, all new checkpoints,
and the 40 reused original final checkpoints. Other historical checkpoint binaries remain in the
preserved full parent archive. Its restore helper recreates expected run aliases without changing data.
No paid resource, historical allowance wrapper, ledger write, repository push, publication, external
contact or new wet-lab experiment is performed.
