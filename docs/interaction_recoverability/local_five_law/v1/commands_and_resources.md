# Commands and final resource accounting — local five-law v1

Completed through **current-allowance job 35**, whose snapshot is dated **2026-09-12 19:56:42.851371 UTC**. There were no additional numerical jobs after job 34; job 35 performed static artifact checks. This final note is authored from the completed ledger audit, within the already charged static-work allowance.

The active ledger is [learned-target-v1-20260912-164000/ledger.jsonl](../../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/ledger.jsonl). This task appended to that ledger. It did not reset an allowance or write to the historical preparation ledger.

## Verified balance and enforcement

At admission, the live ledger through job 31 had conservative charge **1,128.109330367312 seconds**, leaving **6,071.890669632688 seconds**, including the protected **600-second reporting reserve**. The gate checked the exact prior ledger SHA-256:

`7fdd296abb73fa66456019c43685def0a02cac09966394f7c10356dadba908b7`

It also checked that all admitted jobs were closed, the global caps were 7,200 CPU/wall seconds, the reserve was 600 seconds, and frozen learner-code hashes matched. The new `local-five-law-v1` continuation marker debited **150 seconds** conservatively for discovery, static source/agent reads, file authoring, and final reporting outside the timed jobs. This is an allowance debit, not measured CPU. It is counted against both additional 600-second task caps.

[resource_gate.py](resource_gate.py) takes a serial task lock, verifies the live current ledger before each job, and calls the unchanged current-budget watchdog in compute phase. Its job deadline is at most 120 seconds and is bounded by both remaining task charge and global allowance excluding the reporting reserve, with a further five-second admission margin. The existing watchdog uses one-core affinity, one numerical worker, numerical thread environment variables set to one, hidden CUDA devices, a serial run lock, process-tree monitoring, timeout termination, and parent/waited-child accounting. These settings were checked again inside the deterministic script. No affinity violation or timeout occurred in this task.

One numerical job ran at a time. Three independent reviewers performed static mathematical reasoning/read-only inspection; they ran no numerical code, imports of numerical modules, fits, or ledger changes. Their static activity is included in the conservative overhead debit, not described as measured compute.

| Additional task quantity | Actual final value | Authorized cap |
|---|---:|---:|
| Measured wrapper plus waited-child CPU | **9.287816 s = 0.154797 core-min** | 600 s |
| Aggregate admitted compute-job wall time | **9.431077 s = 0.157185 min** | 600 s |
| Conservative task charge, including 150 s overhead | **159.846495 s = 2.664108 min** | 600 s, charged against both caps |
| Remaining within additional task cap | 440.153505 s | No instruction to spend it |
| GPU time | **0** | 0 |
| Numerical workers / threads / core affinity | **1 / 1 / 1** | 1 / 1 / 1 |

Final cumulative current-allowance accounting, including prior work: measured CPU **830.768966 seconds**; aggregate admitted-job wall **834.415012 seconds**; conservative charge **1,287.955826 seconds**; remaining **5,912.044174 seconds**. The 600-second reserve remains intact. Elapsed research time, network/tool latency, and model reasoning are not represented as CPU minutes.

The authoritative final records are [resource_audit_job_35.json](resource_audit_job_35.json) and [resource_snapshot-35.json](../../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/resource_snapshot-35.json). Final current-ledger SHA-256:

`7f30adc8bec106fda1801663a57a31c3578578a2d0259a073a27639f65925225`

## Actual admitted commands and statuses

The following records the actual operational sequence in `/home/shadi/iclr2027`. `register` and final `verify` are static accounting operations covered by the overhead debit. Do not replay registration to reset a task: it refuses an existing marker. Result destinations also refuse overwrite.

```bash
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_five_law/v1/resource_gate.py register

# Job 32: required theorem/audit/source reads; exit 0.
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_five_law/v1/resource_gate.py --seconds 30 run -- cat docs/interaction_recoverability/learned_target_v1/theorem.md docs/interaction_recoverability/learned_target_v1/review.md docs/interaction_recoverability/learned_target_v1/source_manifest.json

# Job 33: initial deterministic check; exit 1 at the absolute-residual assertion.
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_five_law/v1/resource_gate.py --seconds 120 run -- python -B docs/interaction_recoverability/local_five_law/v1/checks.py

# The initial script was preserved as checks_attempt_job33.py.
# Only the cancellation assertion was made scale-aware; no formula,
# direction, reference, quadrature order, or contraction depth changed.

# Job 34: the revised deterministic check; exit 0.
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_five_law/v1/resource_gate.py --seconds 120 run -- python -B docs/interaction_recoverability/local_five_law/v1/checks.py

# Job 35: static artifact/preservation/journal checks; exit 0.
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_five_law/v1/resource_gate.py --seconds 30 run -- python -B docs/interaction_recoverability/local_five_law/v1/finalize.py

# Save final completed-ledger audit; no numerical execution.
PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/local_five_law/v1/resource_gate.py verify
```

[command_records.json](command_records.json) contains the exact start/end ledger records and child argument vectors through completed job 34. Job 35's own complete record is in the live ledger and final audit. The shell prefixes above are reconstructed from the executed gate invocation; the child argument vectors are recorded, not inferred.

| Job | Work | Exit | Measured CPU s | Admitted wall s |
|---|---|---:|---:|---:|
| 32 | Static theorem/review/source-manifest read | 0 | 0.053935 | 0.038516 |
| 33 | Initial identity check; assertion failure retained | 1 | 0.357191 | 0.442278 |
| 34 | Deterministic identity/reconstruction checks | 0 | 8.744824 | 8.783636 |
| 35 | Static final verification | 0 | 0.131866 | 0.166647 |

Other executed operations were static discovery (`pwd`, `git status --short`, applicable `AGENTS.md` search), `cat`/`sed`/`rg` reads of the task, attachment, original model and proofs, direct-target documents, budget implementation and ledger/snapshots; creation of the new files; static Python text edits to those new files; and final rendering inspection. A discovery `rg` for the not-yet-created local namespace returned no directory. No applicable `AGENTS.md` was found. The filesystem sandbox launcher could not initialize in this environment, so shell reads/writes used the authorized elevated execution path; this did not bypass the numerical resource gate.

The focused primary-source web read was [Trabs, arXiv:1307.6610v2](https://arxiv.org/html/1307.6610v2), including Theorem 2.7 and the nonlinear-attainment discussion. No broader literature scan or new experiment was initiated. Static text edits repaired TeX escapes in the new proof before its final hash was recorded.

## Diagnostic outputs and retained failure

The initial script used the absolute assertion `err < 2e-10` for a six-term identity. It stopped before printing its raw residual; no successful residual was invented for that attempt. [checks_attempt_job33.py](checks_attempt_job33.py) is that exact script and [failure_job33.json](failure_job33.json) preserves the failure and reason for changing the assertion.

The rerun first printed the following actual cancellation diagnostics:

```json
{"kind":"trig","absolute_error":9.14610609470401e-10,"term_scale":211482.7228039624,"scaled_error":4.324753328990451e-15}
{"kind":"gaussian_holdout","absolute_error":2.3663915271754377e-10,"term_scale":86786.84697001407,"scaled_error":2.7266706992973893e-15}
```

It then passed the scale-aware condition `err / scale < 2e-12`. The remaining absolute checks were unchanged. Full returned results are in [checks.json](checks.json), with the exact implementation in [checks.py](checks.py).

The Neumann depths were 32, 128, and 512. Profile derivative errors decreased from 5.863508/4.413269 to 0.000393255/0.000281824 to \(2.881\times10^{-9}/8.801\times10^{-10}\). All four pure coefficient directions and the two profile/mixed checks passed; the largest recovered-coefficient error was \(1.754\times10^{-9}\). Profile functional quadrature used only orders 64 and 128, fixed response interval \([-12,12]\) split at 4 and 5, and the predeclared lattice truncation. At order 128 the two functional residuals were \(6.484\times10^{-14}\) and \(4.797\times10^{-14}\).

These are deterministic floating-point checks of proved identities. They are not efficient-variance estimates, finite-dictionary existence certificates, certified quadrature enclosures, or fitted comparisons. The analytic proof separately bounds geometric truncation and tail envelopes. The earlier direct-target check suite and frozen adaptive-enrichment pilot were not rerun or modified.

## Preservation and completion

[preservation.json](preservation.json) commits **360 prior files**, including previous proof/decision documents, frozen pilot code, predictions/configurations/results, earlier snapshots, and the historical preparation ledger. [final_checks.json](final_checks.json) and the final resource audit report **zero changed protected files**, an empty tracked git diff, and no numerical check rerun during finalization. The expected current-ledger append and new snapshot files are not treated as prohibited modifications. Existing root user files were read only.

The final check also validated the new Python syntax without executing it, proof/decision link targets, control-character and math-delimiter consistency, and successful saved diagnostic status. It records hashes of the completed proof, decision, scripts, raw checks, retained failure, and command journal. This final resource note and the completed resource audit are authored afterward and are explicitly outside that manifest; their authoring remains covered by the initial 150-second debit.

No sample generation, model/profile fitting, learned-direction search, bootstrap, benchmark, GPU operation, biological experiment, or additional pilot was performed. No automatic follow-up or further training authorization is requested.
