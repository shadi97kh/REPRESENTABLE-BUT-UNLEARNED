# Executed work and final accounting — local attainment v1

Completed through **current-allowance job 37**, snapshot timestamp **2026-09-12 20:52:24.577650 UTC**. Job 36 was the only new deterministic numerical check. Job 37 performed static artifact verification. No numerical work followed job 36.

## Allowance verification

The live current ledger through job 35 was read before admission. Its conservative charge was **1,287.9558255792535 seconds**, leaving **5,912.0441744207465 seconds**, including the protected **600-second reporting reserve**. The actual prior ledger hash was checked against:

`7f30adc8bec106fda1801663a57a31c3578578a2d0259a073a27639f65925225`

A fresh `local-attainment-v1` continuation marker imposed a **600-second additional conservative-charge cap within that allowance**. It created no new allowance. Registration charged **150 seconds** for initial discovery, static reading, source inspection, agent reasoning/read operations, code/document authoring and final reporting outside admitted jobs. That charge is deliberately conservative overhead, not measured CPU.

[resource_gate.py](resource_gate.py) was derived into this new namespace without modifying the previous gate. It verifies closed jobs, cap/reserve fields and frozen learner hashes, uses a task lock, and invokes the unchanged current-budget watchdog in compute phase. Each deadline is at most 120 seconds and respects both the remaining task cap and the global balance above the reporting reserve. The wrapper applies one-core affinity, one worker/thread, CUDA hiding, serial locking, descendant observation and timeout termination.

| This continuation | Final value |
|---|---:|
| Measured wrapper and waited-child CPU | **0.556072 seconds** |
| Aggregate admitted-job wall | **0.557263 seconds** |
| Conservative task charge, including static overhead | **150.784105 seconds = 2.513068 minutes** |
| Additional task cap | 600 seconds |
| Remaining inside task cap | 449.215895 seconds |
| GPU use | **0** |
| Numerical workers / threads / cores | **1 / 1 / 1** |
| Timeouts / affinity violations | **0 / 0** |

Current-allowance totals after this continuation: measured CPU **831.325038 seconds**, admitted-job wall **834.972275 seconds**, conservative charge **1,438.739931 seconds**, remaining **5,761.260069 seconds**. The 600-second reporting reserve is intact. Research/session elapsed time is not CPU time or aggregate compute-job wall time.

The authoritative records are [resource_audit_job_37.json](resource_audit_job_37.json), [current ledger](../../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/ledger.jsonl), and [snapshot 37](../../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/resource_snapshot-37.json). Final ledger SHA-256:

`f5c820beeeac045c48d3d69634c49179db2ae78598be2a0818e937f9b07bc02e`

## Actual commands

All commands used working directory `/home/shadi/iclr2027`. Registration and final verification are static accounting operations covered by the overhead debit.

```bash
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_attainment/v1/resource_gate.py register

# Job 36: new deterministic conventions and public schedule checks; exit 0.
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_attainment/v1/resource_gate.py --seconds 90 run -- python -B docs/interaction_recoverability/local_attainment/v1/checks.py

# Job 37: static preservation, syntax, links, sample guard and journal; exit 0.
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_attainment/v1/resource_gate.py --seconds 30 run -- python -B docs/interaction_recoverability/local_attainment/v1/finalize.py

# Final completed-ledger audit; no numerical check execution.
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_attainment/v1/resource_gate.py verify
```

Existing output paths and continuation registration refuse overwrite/reset. These commands record what ran; this report does not authorize their replay or another expenditure. [command_records.json](command_records.json) contains job 36's exact child argument vector and start/end record. Finalization's own completed record is in the ledger and resource audit.

| Job | Status | Measured CPU s | Admitted wall s |
|---|---|---:|---:|
| 36 | deterministic checks passed | 0.417774 | 0.390931 |
| 37 | static finalization passed | 0.138298 | 0.166331 |

Static operations included reading the attached task, checking applicable `AGENTS.md` paths and `git status --short`, locating files with `rg`, and reading the requested original/local/direct-target documents, implementation, saved checks, resource gate and live ledger. No applicable `AGENTS.md` was found. The new destination did not exist, so `local_attainment/v1` was created. Its absence produced the expected initial `rg` error. The filesystem sandbox launcher remained unavailable; authorized elevated shell access was used, with numerical execution still routed through the resource gate.

The exact implementation was read rather than rerun. Two parallel mathematical reviewers performed static reasoning and read-only inspections of the central/calibration and tail/arithmetic arguments, then adversarially reviewed the new proof and reference implementation. They executed no numerical jobs and modified no files or ledgers.

The focused primary sources inspected were [Trabs](https://arxiv.org/html/1307.6610v2), [Kaji](https://arxiv.org/pdf/1910.07572), [DExtrI](https://arxiv.org/html/2608.19849v1), and [strongly identified functionals](https://arxiv.org/html/2208.08291v3), together with the saved comparison records. This was a bounded comparison, not a broad search or a claim of exhaustive originality.

## Actual diagnostic outputs

[checks.py](checks.py) used 70-decimal-digit `mpmath` arithmetic. [checks.json](checks.json) records:

- Polygon interpolation of the four hand-written ordinates \((1,4,9,16)\) at probabilities \((0,.125,.25,.375,.5,.875,1)\): \((1,1,1,2.5,4,12.5,16)\). This is an arithmetic fixture, not a sampled source dataset.
- Population normalization: \(g_1(0)=g_2(0)=1\).
- Finite population contraction-remainder discrepancy: **\(1.5244246148248607861\times10^{-36}\)**.
- Sampled estimator entry: **blocked before reading source input**.
- Public domain checks at \(n=10^{12},10^{32},10^{64}\): **false, false, true**, respectively. No observations were created at any count.

The old local-regular `checks.json` was read and its passed status reused. Its check suite was not rerun. The new checks did not evaluate an estimator on samples, calculate influence variance, assess prediction, fit profiles, optimize a network, run a benchmark, or validate biological behavior.

After the numerical check, a static review removed an unnecessary, never-executed domain-condition fallback from the disabled sample entry. The domain inequality remains a proof-onset diagnostic. All functions exercised by job 36 were unchanged. The final static job checked the resulting source syntax and retained disabled execution flag; no numerical rerun was needed.

The supplied `mpmath` backend is not an interval arithmetic certificate. Its small analytic residual supports the tested finite identity only. The theorem's arbitrary-sample numerical guarantee explicitly requires certified primitive/error enclosures; those were not established by this diagnostic.

## Preservation and final records

[preservation.json](preservation.json) hashes **393 prior files**, including root user files, earlier proofs/decisions, failed checks, original implementations, frozen pilot configurations/predictions/results, earlier snapshots and the historical preparation ledger. Final checks report **zero changes** to those files and an empty tracked git diff. Only the authorized current-ledger append/new snapshots and the new attainment namespace were created. The historical preparation ledger was not written.

[final_checks.json](final_checks.json) records syntax/link/control-character checks, the disabled sample flag, no numerical rerun during finalization, preservation results, and hashes of the completed theorem/specification/decision/novelty documents, implementation and checks. This resource note and the final resource audit are authored afterward and are explicitly outside that artifact hash list; their authoring is covered by the initial 150-second charge.

The required deliverables are complete: [theorem](theorem.md), [estimator and prepared validation specification](estimator.md), [novelty](novelty.md), [decision](decision.md), this accounting note, and the justified deterministic implementation/checks. No training request, automatic follow-up experiment, or replacement project is initiated.
