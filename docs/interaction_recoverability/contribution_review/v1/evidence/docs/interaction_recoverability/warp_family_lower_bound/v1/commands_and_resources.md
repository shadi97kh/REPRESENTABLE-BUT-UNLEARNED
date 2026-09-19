# Proof-only execution and preservation record

This task implements the attached request at
/home/shadi/.codex/attachments/cb75172f-e126-4dd7-b58d-5c99816dae82/pasted-text.txt.
The fresh destination docs/interaction_recoverability/warp_family_lower_bound/v1/ was explicitly checked absent before creation. The only authored project files are the six new Markdown deliverables in this directory.

## Actual work

**PROVED HERE / static verification:** derived the saturation bounds and global analytic-profile envelopes; verified the fixed-profile nonlinear coefficient path; derived the actual density score; proved the central and tail score estimates with constants \(40\) and \(c_0=1/1280\); supplied common-coordinate square-root-density domination; proved uniform target separation; derived local information and complete-five-group finite-sample lower bounds; compared them with the retained original-class upper bound.

Inspected the actual warp-family problem, commutator/calibration proof, attainment theorem and score/tangent definitions, conditioning proof and novelty comparison. Inspected the historical efficiency argument to retain its distinct fixed reference and unresolved scope. Searched for applicable AGENTS.md in the repository and parent locations; none was found. Initial tracked and staged diffs were empty, with many preexisting untracked user/research files.

Separate agents inside the same system checked admissibility/target separation, score/path regularity and primary-source overlap. Two further static reads checked the authored proofs and constants. One wording overclaim was corrected: the chosen testing separation is sufficient under the Hellinger upper bound, not the largest possible indistinguishable separation. The nonlinear Hellinger order is not claimed sharp. No remaining concrete contradiction was found in these bounded internal checks.

This is static self-review with internal cross-checks, not independent external peer review or a formal proof-assistant certificate. All mathematical conclusions rest on the displayed analytic proofs.

## Source inspection

The retained source comparison and cached Trabs v2 Theorem 2.7 were read. Heinrich–Kahn arXiv:1507.04313v1 was retrieved from primary HTML; Sections 2.4 and 3.1–3.3 and the adjoining theorem statements were inspected. The source-specific comparison and limitations are in [consequences_and_prior_overlap.md](consequences_and_prior_overlap.md). No requested primary source was unavailable.

## Commands and authoring operations

Shell tools performed file inspection and preservation bookkeeping. Representative exact executed command lines are below; this is a selected record, not a claimed exhaustive transcript:

~~~sh
cat /home/shadi/.codex/attachments/cb75172f-e126-4dd7-b58d-5c99816dae82/pasted-text.txt
git status --short
rg --files --hidden --no-ignore -g AGENTS.md -g '!node_modules' -g '!.git'
test ! -e docs/interaction_recoverability/warp_family_lower_bound/v1
git diff --stat
git diff --cached --stat
cat docs/interaction_recoverability/warp_family/v1/problem.md docs/interaction_recoverability/warp_family/v1/novelty.md
sed -n '1,225p' docs/interaction_recoverability/warp_family/v1/commutator_and_calibration.md
sed -n '1,90p' docs/interaction_recoverability/warp_family/v1/conditioning.md
sed -n '135,280p' docs/interaction_recoverability/warp_family/v1/conditioning.md
sed -n '265,350p' docs/interaction_recoverability/warp_family/v1/attainment.md
sed -n '1,72p' docs/interaction_recoverability/efficiency_audit/v1/efficiency_argument.md
sed -n '277,298p' docs/interaction_recoverability/efficiency_audit/v1/efficiency_argument.md
sed -n '1,56p' docs/interaction_recoverability/warp_family/v1/attainment.md
sed -n '89,164p' docs/interaction_recoverability/warp_family/v1/conditioning.md
sha256sum docs/interaction_recoverability/warp_family/v1/*.md docs/interaction_recoverability/external_review/v2/review_packet.zip docs/interaction_recoverability/external_review/v2/manifest.json docs/interaction_recoverability/efficiency_audit/v1/efficiency_argument.md docs/interaction_recoverability/practical_feasibility/v1/feasibility_decision.md runs/interaction_recoverability/prep-20260912-082154/ledger.jsonl runs/interaction_recoverability/learned-target-v1-20260912-164000/ledger.jsonl
~~~

Additional ordinary text reads supported the internal reviews. Primary-source retrieval used the read-only web tool.
New documents were created with the patch tool. A standard-library Python text replacement corrected the reviewed sufficient-separation wording in the new lower_bound.md, asserting that exactly one matching paragraph existed before replacement. It imported no project or numerical library and ran no estimator or diagnostic. Any final resource-note authoring is also confined to this new namespace.

Shell reads and document authoring used the environment's escalation path because the filesystem sandbox launcher was known to fail when initializing its loopback interface. No computational allowance wrapper was called.

## Preservation scope

Initial SHA256 records cover thirteen sentinels: all seven retained warp_family/v1 Markdown documents; external_review/v2/review_packet.zip and manifest.json; the original efficiency argument; the practical feasibility decision; and both the preparation and learned-target allowance ledgers.

No operation wrote to any historical theorem, decision, erratum, review archive, frozen code, prediction, configuration, user evidence file or ledger. No git commit/push, reviewer contact, publication, packaging, archive regeneration, new figure or backend work was performed.

The sentinel hash check verifies the listed files, not an exhaustive byte scan of every historical untracked artifact. Repository changes are confined by the actual authoring operations to these six new documents. Final check output is recorded below.

## Resource accounting

| Category | Actual task usage |
|---|---|
| Numerical jobs, simulations, samples and experiments | None |
| Fits, training, optimization, Gram calculations or benchmarks | None |
| GPU or backend execution | None |
| Computational allowance wrappers or ledger appends | None |
| Allowance reset, inferred remainder or reserve spending | None |
| File inspection, hashing, web retrieval, text authoring and tool overhead | Nonzero; not comprehensively CPU-metered |
| Mathematical reasoning and internal static review | Performed; not converted into fabricated CPU core-minutes |

No numerical allowance was requested or used. Individual tool wall times do not supply aggregate CPU usage or total research-time accounting. Reporting zero total CPU usage would be false. No estimated overhead debit was written to a historical ledger.

There are no new numerical diagnostic outputs to list: the identities, constants and inequalities were checked analytically. The historical family-attainment theorem, fixed-reference efficiency question, backend no-go and negative pilot remain preserved. The [decision](decision.md) starts no additional training, backend, dictionary or packaging cycle.

## Final preservation check

All thirteen final SHA256 lines matched the initial records exactly. This includes every retained warp_family/v1 Markdown file, the v2 archive and manifest, the original efficiency argument, the backend no-go decision and both allowance ledgers. The equality check compared hash-output strings; it was file-integrity bookkeeping, not a numerical proof check.

The tracked and staged git diffs remained empty, and the tracked whitespace check emitted no errors. The preexisting untracked files remain present; the working tree is not described as clean. The six required new Markdown files were present.

The unchanged v2 archive SHA256 is:
9cd36e33be3778c87b37d21c9573769f72d8d17cae28d8e31b86daed4cbb86d5

The unchanged preparation-ledger SHA256 is:
4a902b8d810782c341430808b03b0f058459a32864b9d6724b401a05d4a31d5f

The unchanged learned-target-ledger SHA256 is:
d1b365a6f7a3e33ad511b1ce2712367c6e9f6602a34bd742cb75f20009ff8cc1

Final shell inspection also used:
~~~sh
git diff --stat
git diff --cached --stat
git diff --check
rg --files docs/interaction_recoverability/warp_family_lower_bound/v1
sed -n '138,195p' docs/interaction_recoverability/warp_family_lower_bound/v1/lower_bound.md
~~~

This final note was appended only to the new resource document after those checks. No frozen file or computational ledger was changed.
