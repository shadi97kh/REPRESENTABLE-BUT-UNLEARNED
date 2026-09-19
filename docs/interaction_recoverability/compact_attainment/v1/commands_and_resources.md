# Executed commands and final resources

The bounded continuation completed through current-ledger **job 42**, at **2026-09-12 22:17:51 UTC**. Its [final resource audit](resource_audit_job_42.json) is authoritative for this task. The two numerical jobs were fixed analytic/component and public-count checks; the third job was static syntax/link/hash verification. No observations were generated, sampled estimator executed, model fitted, training/tuning performed, or pilot rerun.

## Authorization and preservation

The live learned-target ledger was verified through job 39 before any numerical execution. Its baseline SHA256 was
`b239393f80ab7356959ac5f31ec15194f93e8a425b9de2bb60a0af5c2fe486a7`,
with 5,632.056218 seconds remaining, including the existing 600-second reporting reserve. This continuation was restricted to **600 additional conservatively charged seconds**, using only that existing balance.

The new [resource gate](resource_gate.py) checked the original caps, closed jobs and frozen learner code, registered this task once, and delegated jobs to the unchanged current-ledger watchdog. It constrained each deadline by the task balance and global balance above the reserve. Affinity was core 0; one worker and one numerical thread were enforced; CUDA devices were hidden. The wrapper followed/terminated process trees and serialized numerical jobs. All three jobs exited zero without timeout or observed affinity violations.

The preparation ledger was not touched. The preservation manifest and final verification found **zero changes among 1,111 protected existing files**, including historical proofs, user files, results, prediction/configuration artifacts, publication files and the preparation ledger. Only the current allowance ledger/snapshots and the fresh compact-attainment namespace were extended. No applicable AGENTS.md was found in the repository or ancestor instruction locations. Git's tracked diff remained empty.

## Final accounting

| Quantity | This continuation | Interpretation |
|---|---:|---|
| Measured wrapper plus waited-child CPU | **1.290445 s** = .0215074 core-min | Sum of jobs 40–42 |
| Aggregate admitted job wall time | **1.347179 s** = .0224530 min | Sum of jobs 40–42 |
| Conservative setup/reporting debit | **150 s** | Discovery, static reads, agent inspection, document/tool authoring and final reporting overhead; not measured CPU |
| Total conservative charge | **151.647179 s** = 2.527453 min | Below the 600-second task cap |
| Unused task cap | **448.352821 s** | No authorization for additional work inferred |
| Current cumulative allowance charge | **1,719.590961 s** | Existing allowance, never reset |
| Current allowance remaining | **5,480.409039 s** | Includes untouched 600-second reporting reserve |
| Remaining above reserve | **4,880.409039 s** | No further run initiated |
| GPU / workers / numerical threads / cores | **0 / 1 / 1 / 1** | Matches constraints |

Each job was charged \(\max(\text{measured CPU},\text{wall})+.1\) seconds. The 150-second overhead debit was appended at registration and includes the final gate/read/report work outside job timing. These resource figures do not confuse elapsed mathematical research time or tool latency with aggregate numerical job wall time.

| Job | Purpose | Measured CPU s | Wall s | Conservative charge s | Status |
|---|---|---:|---:|---:|---|
| 40 | Fixed analytic identities, population iteration, polygon quadrature, norm/count formulas, guard | .889448 | .895913 | .995913 | 0 |
| 41 | Public work/precision formulas at two fixed counts | .127009 | .165443 | .265443 | 0 |
| 42 | Static syntax, links, sample refusal and preservation | .273988 | .285822 | .385822 | 0 |

Final live ledger SHA256:
`265b269b4f32cd5fdd713223679da1a249b0dfa748e552a10171a00839b533da`.

## Exact operational commands

These commands were executed from `/home/shadi/iclr2027`. They are a record, not an instruction to spend remaining allowance. Registration and result writers refuse overwriting existing records.

```bash
python -B docs/interaction_recoverability/compact_attainment/v1/resource_gate.py register

python -B docs/interaction_recoverability/compact_attainment/v1/resource_gate.py --seconds 90 run -- python -B docs/interaction_recoverability/compact_attainment/v1/checks.py

python -B docs/interaction_recoverability/compact_attainment/v1/resource_gate.py --seconds 30 run -- python -B docs/interaction_recoverability/compact_attainment/v1/practical_counts.py

python -B docs/interaction_recoverability/compact_attainment/v1/resource_gate.py --seconds 30 run -- python -B docs/interaction_recoverability/compact_attainment/v1/finalize.py

python -B docs/interaction_recoverability/compact_attainment/v1/resource_gate.py verify
```

The [exact start/end journal](command_records.json) copies the ledger's child argv, UTC times, statuses, deadlines, RSS observations and accounting. [checks.py](checks.py) creates [checks.json](checks.json); [practical_counts.py](practical_counts.py) creates [practical_counts.json](practical_counts.json). Those files have not been rerun to improve a result. No check failed in this pass.

Static discovery used `rg`, file reads and `git status --short`; authoring used patches and text-only Python edits. The original attachment, saved local-attainment theorem/code/checks/resource records and review, local-five-law proof, unknown-profile model and prior source comparison were read. Primary source pages for Trabs, Kaji and NIST quadrature were reread as linked in [novelty.md](novelty.md). A file-update tool encountered the sandbox launcher failure; authorized local shell operations used sandbox escalation. This was an access failure, not a failed numerical experiment or a rejected scientific result.

## Scope of validation

[Final static checks](final_checks.json) confirm Python syntax, the unconditional sample-entry refusal and the existing mathematical document links. [Bundle hashes](bundle_hashes.json) cover the completed mathematical/code/check bundle before this final accounting note and its final audit/journal were authored. The live ledger and these later accounting records are explicitly outside that earlier bundle digest.

The numerical outputs validate narrowly specified components. They do not certify the full interval backend, run the complete estimator, demonstrate finite-sample accuracy, or validate biological siRNA assumptions. Earlier unsuccessful experiments and negative decisions remain preserved.

