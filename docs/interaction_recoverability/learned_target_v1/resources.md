# Final resource accounting

Status: completed through **new-allowance job 28**, recorded at **2026-09-12 17:35:49.528169 UTC**. No further fits or numerical experiments followed evaluation. This note is written from the completed [final snapshot](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/resource_snapshot-28.json); its file-authoring overhead is included in the conservative overhead debit described below.

| Quantity | Final value | Authorized cap / interpretation |
|---|---:|---|
| Measured CPU, wrapper plus waited descendants | 814.115899 seconds = **13.5686 core-minutes** | 7,200 seconds = 120 core-minutes |
| Aggregate admitted compute-job wall time | 817.592598 seconds = **13.6265 minutes** | 7,200 seconds = 120 minutes |
| Conservative cumulative allowance charge | 970.392598 seconds = **16.1732 minutes** | Debited against both caps |
| Remaining under conservative charge | 6,229.607402 seconds = **103.8268 minutes** | Not permission for another experiment |
| Protected reporting reserve | 600 seconds | Retained when admitting compute jobs; report stages explicitly marked |
| GPU time | **0** | Zero authorized |
| Numerical workers / threads / core affinity | **1 / 1 / 1** | One of each; no parallel fitting |
| Final-dataset fitter totals only | 680.456328 CPU seconds; 753.099226 wall seconds | All 24 datasets and shared fitted work counted once; excludes development/wrapper/reporting |

The ledger begins at 16:40:49.964458 UTC. The interval to the final snapshot is approximately **54.99 minutes of elapsed research/session time**, including source reading, proof writing, implementation, tool latency and reporting. It is not 54.99 minutes of CPU work or aggregate compute-job wall time. Pre-initialization discovery is covered by the setup debit; prior sessions are not used to infer allowance.

The accounting is deliberately conservative. Each wrapped job is charged `max(measured parent+waited-child CPU, aggregate job wall)+.1 seconds`. Inherited single-core affinity supplies a wall-based bound for otherwise unreaped descendant CPU. The fresh initialization also debited **30 seconds** for unmetered discovery/bootstrap. Finalization debited an additional **120 seconds** for implementation/document-authoring and image/file-tool overhead outside shell timing, including creation of this note. Those 150 seconds are conservative allowance debits, **not measured CPU values**. Neither debit resets, refunds or increases the allowance.

## Enforcement and failures

The new wrapper sets CPU affinity, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS` and related numerical settings to one, and hides CUDA devices. A serial lock admits one numerical job. `/proc` monitoring follows the observed process tree, including a child in another process group; timeouts terminate and then kill the observed tree. Per-job deadlines respect remaining allowance and the 600-second compute reserve. Each job saves UTC times, command argv, exit status, CPU/wall, RSS and affinity observations in an append-only ledger.

The only nonzero job was the **intentional timeout test, new job 10**: wrapper status 124, main child terminated by SIGTERM, with a separately grouped child confirmed absent afterward. Both processes were sleeping test fixtures, not concurrent numerical workers. All actual development/final fit and endpoint outputs completed. Numerical warnings and poor prediction outcomes were retained. There were no affinity violations in that timeout record, and no GPU call or additional worker was used by the learners.

The default local image viewer failed because the filesystem sandbox launcher could not initialize its loopback interface. Existing PNG bytes were subsequently read under this new wrapper and displayed unchanged. These were access/inspection failures, not missing figures or failed model outputs; no image was regenerated to alter a result.

## Historical-ledger exception and preservation

The very first discovery read mistakenly used the historical preparation wrapper before the new prompt's prohibition was read, creating **historical job 27** with a **0.101903-second conservative debit**. This exception was disclosed during execution. The record was preserved rather than erased, and discovery was also covered by the new setup debit. No subsequent command used or changed the old allowance. The new ledger links the historical ledger hash **after that disclosed append**. This does not satisfy literal zero historical-ledger writes, and is not hidden by the otherwise successful preservation check.

Final checks found **zero changes among the 55 protected historical/user files** relative to the new preservation manifest, verified frozen learner code and private-spec commitments, and matched all 24 saved prediction/source hashes. The tracked git diff remained empty; newly created files were kept in their separate namespace. Historical proofs, failed experiments, user evidence files and default-disabled interfaces were not overwritten.

## Reproducibility records

- [Append-only new ledger](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/ledger.jsonl) — authoritative command/status/accounting record, including finalization job 28.
- [Final snapshot, job 28](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/resource_snapshot-28.json) and [final checks](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/final_checks.json).
- [Executed operational commands](commands.sh) and [exact child argv/status journal](command_records.json) through completed job 27. Wrapper prefixes in the shell journal are reconstructed from recorded phase/deadline; child argv are exact. Finalization's own status is in the ledger.
- [Artifact hashes](../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/artifact_hashes.json) protect the frozen results and report bundle. The manifest explicitly excludes the live ledger/snapshots and this final note, which was authored from the completed snapshot.

Existing destinations deliberately refuse replay. The commands document what actually ran; no additional training campaign or expenditure is requested or initiated by this report.
