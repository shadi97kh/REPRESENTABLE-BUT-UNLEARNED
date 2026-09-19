# Executed work and final resource accounting

**NUMERICALLY CHECKED / ledger record.** This bounded feasibility pass completed through current-ledger **job 51**, at **2026-09-12 23:14:25 UTC**. No numerical work followed the final static verification. The final [resource audit](resource_audit_job_51.json) and [exact command records](command_records.json) contain the authoritative values.

## Verified authorization and enforcement

The live learned-target ledger was checked through job 42 before computation. Its SHA256 was
`265b269b4f32cd5fdd713223679da1a249b0dfa748e552a10171a00839b533da`.
The existing balance was 5,480.409039 seconds, including the existing 600-second reporting reserve. A fresh task marker imposed a **300-second additional conservative cap** within that balance; it did not reset or create an allowance.

The new [task gate](resource_gate.py) reuses the unchanged current-ledger watchdog. It verifies frozen learner hashes, original caps, closed jobs and the task balance before admitting a job. The numerical jobs inherited one core (affinity 0), one worker and one numerical thread, with CUDA hidden. The watchdog serializes jobs, monitors the process tree, and constrains each deadline by the task remainder and global balance above the reporting reserve.

Registration debited **90 seconds** for discovery, static proof/agent reads, authoring, approval reconciliation and final report/verification overhead outside job timing. This is a conservative charge, not measured CPU. Subsequent jobs are charged \(\max(\mathrm{CPU},\mathrm{wall})+.1\) seconds. Report authoring does not reset or refund either charge.

| Quantity | Final value |
|---|---:|
| Measured wrapper plus waited-child CPU | **43.551225 seconds** = .725854 core-minutes |
| Aggregate admitted job wall time | **43.792724 seconds** = .729879 minutes |
| Setup/static/reporting debit | **90 seconds** |
| Total conservative task charge | **134.710766 seconds** |
| Authorized task cap | **300 seconds** |
| Unused task cap | **165.289234 seconds** |
| Current cumulative allowance charge | **1,854.301727 seconds** |
| Current allowance balance | **5,345.698273 seconds** |
| Existing reporting reserve | **600 seconds, preserved** |
| Balance above reserve | **4,745.698273 seconds** |
| GPU / numerical workers / threads / cores | **0 / 1 / 1 / 1** |

Elapsed mathematical research and tool latency are not aggregate numerical-job wall time. The unused allowance is not a recommendation or authorization for another campaign.

## Actual jobs, preflight and failures

**NUMERICALLY CHECKED.** All admitted jobs returned zero. No timeout, observed affinity violation, failed mathematical assertion or missing result occurred. Runtime records are for resource enforcement, not a benchmark comparison.

| Job | Purpose | CPU s | Wall s | Conservative charge s |
|---|---|---:|---:|---:|
| 43 | Initial complete reference influence | 2.097506 | 2.126689 | 2.226689 |
| 44 | Central-depth check | 3.579910 | 3.575510 | 3.679910 |
| 45 | Integration-resolution check | 4.325674 | 4.313844 | 4.425674 |
| 46 | Response-tail cutoff check | 4.385716 | 4.458982 | 4.558982 |
| 47 | Depth with small analytic remainder | 7.906075 | 7.947260 | 8.047260 |
| 48 | Final integration refinement | 10.591585 | 10.589773 | 10.691585 |
| 49 | Rational remainder certificates and precision formulas | .132407 | .180853 | .280853 |
| 50 | Accumulation-precision comparison | 10.253810 | 10.320878 | 10.420878 |
| 51 | Static syntax, links, provenance and preservation | .278542 | .278934 | .378934 |

The influence code preflights an atom-count bound and rejects configurations with ten million or more atoms. The allowed configurations were fixed reference refinements, each admitted under a 25–45-second watchdog; observed costs from smaller configurations informed later admission. No quadratic atom covariance matrix was built. Arrays were processed by source and kernels at equal nodes were merged.

**Concrete accounting qualification:** the saved `preflight_storage_upper_bytes` field is nominal array accounting, not an RSS guarantee. It omits runtime and temporary storage; observed base RSS exceeded that field. The final extended-precision process used approximately 1.23 GiB. No RAM cap was authorized, and CPU/wall enforcement remained active. The original record and code hash are preserved; this naming error is disclosed rather than silently rewriting diagnostics.

**Automatic approval review event:** an attempted 45-second deadline for the final integration refinement was rejected before process creation. The stated reason was that approximately 270 of 300 seconds had already been consumed. A read-only live-ledger reconciliation instead verified **112.938516 seconds actually charged**, including the 90-second debit, with **187.061484 seconds remaining**. After presenting that evidence, the same direct gated refinement was admitted with a 30-second deadline and completed as job 48. The rejected invocation generated no job or numerical output, and no indirect execution, allowance reset or historical-ledger access was used to bypass the review.

## Exact executed operational commands

The working directory was `/home/shadi/iclr2027`. These are records of completed actions, not instructions to spend the remaining balance. Result writers refuse existing output names.

```bash
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py register

python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 45 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/influence_variance.py variance_base
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 45 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/influence_variance.py variance_depth --depth 256
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 45 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/influence_variance.py variance_integration --depth 256 --order 128 --step 0.125
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 45 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/influence_variance.py variance_tail --depth 256 --order 128 --step 0.125 --radius 12
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 45 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/influence_variance.py variance_deep --depth 512 --order 128 --step 0.125 --radius 12
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 30 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/influence_variance.py variance_final --depth 512 --order 256 --step 0.0625 --radius 12
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 10 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/error_and_precision.py
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 25 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/influence_variance.py variance_accumulator --depth 512 --order 256 --step 0.0625 --radius 12 --precision float64

python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py --seconds 10 run -- python -B docs/interaction_recoverability/practical_feasibility/v1/finalize.py
python -B docs/interaction_recoverability/practical_feasibility/v1/resource_gate.py verify
```

The ledger copies exact child argv, UTC times, statuses, deadlines, CPU/wall and memory observations in [command_records.json](command_records.json). Shell discovery used `rg`, applicable ancestor instruction checks, `git status --short`, `git diff --stat` and file reads. No applicable AGENTS.md was found. Missing optional `.agents`/`.codex` directories and the absence of AGENTS matches caused discovery commands to return nonzero; these were not computational failures.

## Preservation and result scope

**NUMERICALLY CHECKED / static verification.** [Final checks](final_checks.json) found **zero changes among 1,135 protected existing files**. The completed compact-attainment theorem/code/checks, earlier proofs, all saved experiments and negative results, user files, publication artifacts, configurations and predictions remain unchanged. The historical preparation ledger was not written. Only the fresh feasibility namespace and current learned-target ledger/snapshots were extended. Git's tracked diff remained empty.

The seven raw influence records all match the saved [source-code hash](influence_variance.py). [Diagnostic summary](diagnostic_summary.json) contains differences of saved outputs, not certified error estimates. [Bundle hashes](bundle_hashes.json) cover the completed mathematical/code/raw bundle before this final accounting note, final audit and exported command journal were authored. The final ledger SHA256 is
`f446782bef0bf20ac16c99bf801c2f1a5e0d569959cb01c54e63982a1ceb6a96`.

**OPEN:** full numerical integration/primitive certification, actual estimator runtime, finite-sample performance, efficiency and biological validity. No samples, model fits, training, bootstrap, Monte Carlo, benchmark campaign or pilot rerun occurred. No full certified backend was built and no subsequent validation is requested.

