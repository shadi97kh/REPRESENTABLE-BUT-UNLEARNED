# Commands, preservation and administrative cost

This task authorized static review, analytic derivation, primary-source reading, authoring, copying and hashing. Every write made for this task is under the newly created `docs/interaction_recoverability/contribution_review/v1/`. The directory was unused when selected. No existing historical file was edited, no Git commit/push was made, and no reviewer contact or publication occurred.

## Instructions and initial state

The root brief [Codex_Combined_Review_GNN_siRNA.md](evidence/Codex_Combined_Review_GNN_siRNA.md) was read. Searches for applicable `AGENTS.md` in the repository and its parent path found none; no missing local instruction was silently bypassed. Initial Git status contained extensive existing untracked user documents and results. Initial tracked and cached diffs were empty. [Status before](administrative/git_status_before.txt), [tracked diff before](administrative/git_diff_before.patch) and [cached diff before](administrative/git_diff_cached_before.patch) preserve the observed state. The new review directory may appear in the saved status because the administrative snapshot was created there first.

Before authoring, 988 historical files were hashed. The exact scope is **Git-tracked and nonignored untracked files, the full frozen learned-target run, `runs/revision`, and ledger-name matches**. It includes protected proof families, lower-bound v1/v2, external_review/v2, frozen pilot source/private/prediction files, configurations, old negative results and ledger files. It is not a claim to have exhaustively hashed every ignored file or dependency anywhere on the machine. [The baseline](administrative/preservation_before.json) records path, size and SHA-256.

## Operations actually performed

The terminal record consists of read-only file searches/reads, standard-library JSON metadata inspection and the new administrative authoring/packaging commands. The command families used were `pwd`, `git status --short`, `git diff --binary`, `git diff --cached --binary`, `git ls-files`, `rg`/`rg --files` (including hidden/ignored paths for protected evidence), `cat`, `sed`, `head`, and short `python -` administrative scripts. Python was used only for filesystem/JSON inspection, counting records/file sizes, copying, hashing, link checking and ZIP packaging. It never imported a project module or scientific library.

The repository evidence reads covered the user-requested proof directories and their substantive dependencies, actual profile/Riesz and graph code, frozen configuration/prediction records, executed P4 training histories, siRNA requirements and data/provenance audits. Exact-name searches also checked for the missing chemistry-design brief and transfer-evidence JSON under the accessible user directory. No contents were fabricated when files were absent. Some initial guessed filenames were absent; the actual paths were then found with `rg`.

Cached primary texts were read first. Bounded primary HTML retrieval covered the exact versions reported in [the mathematical review](mathematical_review.md). No large source acquisition, supplement download or new data-processing campaign took place. Web/tool service costs and the language model's reasoning cost are not included in local process timing.

The default shell sandbox launcher failed before running commands with a loopback setup error (`bwrap: ... RTM_NEWADDR ... Operation not permitted`); the patch helper also failed. Read/write shell calls therefore used the approved escalation path, with their authorized administrative purpose stated. These were sandbox-launch failures, not automatic approval-review rejections. No scientific permission or training allowance was requested or inferred.

New Markdown was authored with literal quoted heredocs. The exact administrative helpers are preserved:

```bash
python docs/interaction_recoverability/contribution_review/v1/administrative/prepare_packet.py
python docs/interaction_recoverability/contribution_review/v1/administrative/finalize_packet.py
```

`prepare_packet.py` copies unchanged evidence and the closure of active local Markdown links, checks byte identities, verifies frozen code/prediction and P4 artifact commitments, and writes the inventory. `finalize_packet.py` verifies preservation, validates local navigation, writes the manifest, packages only this new version and verifies the completed ZIP. These commands are administrative, not project tests or numerical allowance-wrapper calls. The first finalization stopped before ZIP creation because the new dependency map linked one historical estimator document not yet copied. The check exposed that packaging omission; its unchanged dependency closure was added and preparation/finalization repeated. A second finalization stopped at seven inherited GitHub-style source-line fragments; the checker was extended to validate their line numbers without altering evidence, and finalization repeated. [Operation notes](administrative/operation_notes.json) preserve both failed administrative attempts.

No sampling, source-data generation, fitting, training, optimization, quadrature, Gram construction, benchmark, checkpoint inference, backend development, GPU access or scientific tests were performed. Existing `.py`, `.sh` and historical command records in `evidence/` are copied text, not commands executed by this review. No allowance wrapper was called and no historical ledger was appended to or rewritten.

## Preservation and integrity results

[The final comparison](administrative/preservation_after.json) checks all 988 baseline files for deletion, size or hash changes and verifies that tracked/cached diffs equal the baseline. It also compares the current Git file set after excluding this review directory. Protected history is unchanged within that exact scope.

[The evidence inventory](administrative/evidence_inventory.json) records **352 unchanged copied files**, totaling **18,516,062 bytes** before compression. It verifies **50 historical SHA-256 commitments**: 12 frozen learner-code files, all 24 frozen prediction files, 13 P4 cycle artifacts and the P4 configuration. Every commitment matches. Checkpoint checks use raw bytes only. These checks concern preservation/provenance, not correctness of saved scientific values.

[The integrity report](integrity_checks.json) and the outside `archive_checks.json` report state the number of payload files, local links, archive entries, byte hashes, manifest coverage and costs. Local Markdown file destinations and inherited `#L` source-line locators are checked; scrolling to such locators depends on the reader’s renderer. Remote URLs are not mass-probed. ZIP paths must be relative, contain no parent traversal, be unique, and match local payload bytes. All archive entries are read to verify bytes/CRC without running copied code. The external `review_packet.zip.sha256` hashes the final archive. There is no circular self-hash claim.

The changed/new file list is exactly this fresh directory: five requested review documents, README/dependency map, copied evidence, administrative snapshots/helpers, manifest, integrity/archive reports, ZIP and its sidecar. All inherited files outside it are preserved. There is no deletion, commit, push or change to the original user brief, proofs, frozen pilot, configurations, ledgers or external_review/v2 packet.

## Costs and limits of accounting

Administrative work used nonzero CPU, wall time, memory and filesystem I/O. The initial hash helper recorded about **0.109 process CPU seconds** and **0.109 wall seconds** before writing its report. The initial evidence-copy/hash invocation recorded about **0.171 process CPU seconds** and **0.171 wall seconds**. After the link-check correction, its second invocation recorded about **0.178 process CPU seconds** and **0.192 wall seconds**. The [initial inventory](administrative/evidence_inventory_initial.json) and final inventory preserve separate counters and cutoffs; both failed finalizations are unmetered apart from the shell tool’s wall-time reports. These are not one free or single preparation pass. The final helper records its own scoped counters in the outside `archive_checks.json` report; these are distinct regions, not a whole-session ledger.

Earlier shell startup, file reading, searching, unsuccessful sandbox/patch attempts, authoring, retrieval and orchestration were not comprehensively metered. Their cost is unknown and is **not zero**. No overall session CPU total, peak process-tree memory, remaining historical allowance or new allowance balance is asserted. Any preparatory metadata/link-closure inspections outside the measured helpers are also unmetered administrative work. No scientific runtime or GPU allocation was incurred by this task; that does not make its administrative cost free.

The backend remains no-go and the historical reference-efficiency question remains open. Completion of this static review does not authorize a replacement research or training cycle.
