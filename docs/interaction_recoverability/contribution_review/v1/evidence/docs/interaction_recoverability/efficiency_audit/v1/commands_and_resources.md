# Executed commands and final resource accounting

The bounded efficiency audit completed through current-ledger **job 52**, ending at **2026-09-12 23:43:25.760336 UTC**. The single numerical job contained the two frozen deterministic configurations and returned zero. No further numerical job, sample, fit, simulation, tuning, pilot rerun, expanded dictionary or backend development followed it.

## Verified nested authorization

The live ledger initially matched the completed practical-feasibility record through job 51, SHA256
`f446782bef0bf20ac16c99bf801c2f1a5e0d569959cb01c54e63982a1ceb6a96`.
That task had consumed **134.710766 seconds of its 300-second cap**, leaving **165.289234 seconds**. Its existing global balance and 600-second reserve were also verified.

The new [resource gate](resource_gate.py) imposes a **120-second nested cap** while continuing to enforce the practical-feasibility 300-second cap and the unchanged global reserve. It does not create or reset an allowance. Registration debited **60 seconds** for discovery, source reading, internal static reviews, authoring, tool/file overhead, final verification and this accounting note outside measured jobs. This is a conservative debit, not measured CPU or elapsed conversation time.

The gate reuses the unchanged [current-ledger watchdog](../../../../interaction_recoverability/learned_target/budget.py), verifies closed jobs and frozen learner hashes, and bounds job admission by the minimum of all remaining caps. The existing wrapper enforces one core, one numerical worker/thread, hidden CUDA devices, a serial lock and process-tree timeout/cleanup. The numerical deadline was 30 seconds; the job completed in 14.885080 seconds with no observed affinity violation or timeout.

The calculation's SHA256 was frozen in the registration event **before any new numerical output**:
`1d21feb790e08d360ff057a4f4486a6ffbc9355861cc335b76a644c9d90982e6`.
It fixes exactly the requested twelve directions, nested sizes 4/8/12, and the two integration configurations. Neither code nor dictionary changed after registration.

## Final account

The authoritative exported [resource record](resources.json) contains the exact ledger start/end events and artifact hashes.

| Quantity | Final value |
|---|---:|
| Measured wrapper plus waited-child CPU | **14.817733 seconds** |
| Aggregate admitted job wall time | **14.885079 seconds** |
| Job charge: max(CPU, wall) + .1 | **14.985079 seconds** |
| Conservative setup/static/reporting debit | **60 seconds** |
| Total additional audit charge | **74.985079 / 120 seconds** |
| Unused nested audit cap | **45.014921 seconds** |
| Practical-feasibility total after this audit | **209.695845 / 300 seconds** |
| Remaining practical-feasibility cap | **90.304155 seconds** |
| Current cumulative global charge | **1,929.286807 seconds** |
| Current global balance | **5,270.713193 seconds** |
| Reporting reserve | **600 seconds, preserved** |
| Balance above reserve | **4,670.713193 seconds** |
| GPU time / numerical workers / threads / cores | **0 / 1 / 1 / 1** |

The numerical child itself reports 13.011609 CPU seconds and 14.350069 wall seconds; those exclude parts of wrapper/import/startup overhead and are not substituted for the charged figures. It streamed source-zero covariance contributions rather than rerunning the saved full-variance calculation. The ledger records memory observations, but this is not a fitted benchmark or full-backend runtime estimate.

Unused allowance supplies no reason to run another experiment. This audit stops at its unresolved decision.

## Actual operational commands

Working directory: `/home/shadi/iclr2027`. All three commands completed with exit status **0**. Writers refuse an existing result or registration; these are provenance records, not instructions to replay completed work.

```bash
python -B docs/interaction_recoverability/efficiency_audit/v1/resource_gate.py register

python -B docs/interaction_recoverability/efficiency_audit/v1/resource_gate.py --seconds 30 run -- python -B docs/interaction_recoverability/efficiency_audit/v1/calculate.py

python -B docs/interaction_recoverability/efficiency_audit/v1/resource_gate.py verify
```

Job 52's exact child argv is `["python", "-B", "docs/interaction_recoverability/efficiency_audit/v1/calculate.py"]`. The start/end records in [resources.json](resources.json) preserve UTC times, status, deadline, CPU/wall and affinity. The current [append-only ledger](../../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/ledger.jsonl) and [snapshot 52](../../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/resource_snapshot-52.json) remain authoritative.

Read-only discovery used `git status --short`, ancestor `AGENTS.md` checks, `rg --files` for applicable instructions/destinations, and `cat`/`sed` on the attached request and the linked model, compact theorem, local proof, practical variance code/raw/accounting and pilot novelty/negative-result records. No applicable AGENTS.md was found. The namespace was unused. Existing untracked user/research files were retained. The source comparison reread [Trabs, Theorem 2.7](https://arxiv.org/html/1307.6610v2); it did not run a new literature or training campaign.

After calculation, text-only Python reads selected fields from `raw.json` for reporting. A final static Python heredoc used `ast.parse` on the two new Python files, checked existing local Markdown links, verified the saved twelve-direction/two-configuration record and successful status, and ran `git diff --quiet`. It returned **0** and printed:

```text
Static syntax, local links, frozen output shape, and empty tracked diff verified.
```

These text/static operations, preservation hashing, source-tool overhead and final authoring are included in the 60-second conservative overhead debit. They are not additional numerical jobs. The final accounting note was intentionally excluded from its own artifact hash manifest because it was authored from the completed resource record.

## Errors and preservation

No numerical assertion failed, no result was discarded, and no numerical job timed out. The file-patch tool failed once while attempting to edit the new argument document because its filesystem sandbox could not initialize loopback (`bwrap: Failed RTM_NEWADDR: Operation not permitted`). That operation made no change. A direct text-only Python edit then added the internal reviewers' qualifications to the new file, with exit status zero. This was an access-tool failure, not a calculation failure or an automatic approval rejection.

Final verification found **zero changes among 1,167 protected pre-existing files**. The [preservation manifest](preservation.json) includes all existing proof/result/code/publication files covered by the repository's prior workflow and explicitly hashes the historical preparation ledger. Existing predictions, configurations, pilot failures, negative decisions and user files remain unchanged. The tracked git diff remained empty.

Only the fresh efficiency namespace and the authorized current ledger/snapshot were extended. The historical preparation ledger was not written. The final current-ledger SHA256 is
`d1b365a6f7a3e33ad511b1ce2712367c6e9f6602a34bd742cb75f20009ff8cc1`.

This was an internal mathematical and static code audit by the primary agent and two separate agents in the same system. They performed no numerical jobs or ledger writes. It is not an external independent review. The conclusion, limitations and exact remaining projection gap are in [decision.md](decision.md) and [efficiency_argument.md](efficiency_argument.md).

