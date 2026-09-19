# Commands, preservation and resource accounting — v2

This proof-only continuation follows the handoff supplied at
/home/shadi/.codex/attachments/a29769e7-baad-47b7-b189-09046205cd00/pasted-text.txt
and the user's instruction to continue the previous chat. The continuation date is 2026-09-13.

## Starting state and version choice

The workspace already contained a completed warp_family_lower_bound/v1 with all six requested documents. That version was inspected and preserved. The destination warp_family_lower_bound/v2 was checked absent; its directory was then created with an exclusive operation. No existing version was overwritten.

Parent locations and the full repository were searched for applicable AGENTS.md files, including hidden and ignored paths. None was found. Matches in unrelated cached plugins were outside this project's instruction scope and were not applied. No skill was needed for this bounded mathematical proof task.

Initial tracked and staged git diffs were empty. Many preexisting untracked research files and directories were present, including the original lower-bound v1; this workspace was not described as clean.

## Actual mathematical work

**PROVED HERE — static verification.** The supplied seed and v1 were rederived through the following obligations:

| Obligation | Verification |
|---|---|
| Global saturation | Hyperbolic derivatives, central exponential residuals, global transition decay |
| Original profile class | Analyticity, normalization, positivity, boundary limit and derivatives through order three |
| Coefficient path | Fixed profiles and internal shift, positive logarithm, original box, exact amplitude identity |
| Five-law invariance | Sources 0, 2, 3 and 4 unchanged; only source 1 varies |
| Density score | Pushforward Jacobian retained in \(S_t=zw_t-w_t'\) |
| Uniform score bound | Central and tail estimates give \(40(\epsilon_\delta+\delta|t|)\), \(c_0=1/1280\) |
| Likelihood segment | Common-coordinate square-root-density domination and \(L^2\) integration |
| Interaction | Exact derivative and uniform negative margin \(e^{11/10}/10\) |
| Local bound | Complete \(1/5\) score weighting and full-compatible-influence inequality |
| Finite bound | Product affinity, both separation constraints, squared-loss test and \(N=5n\) conversion |
| Upper comparison | Saved uniform influence bound and uniform theta quantifiers |
| Interpretation | Local variance, fixed reference, CLT, finite MSE and fixed-accuracy counts kept distinct |

No failing inequality was found in this bounded verification. The valid v1 construction, admissibility/score proof and lower-bound proof were retained by copying their text into the new version after checking them. Explicit current provenance and expanded algebraic checks were then added. The current decision, prior comparison and resource record were authored separately. This is a verified continuation of the same theorem, not a claim to have discovered a stronger rate.

No sub-agents were spawned in this continuation. The previous version's own account of its internal agent reviews is historical and is not attributed to this run. The present mathematical check is one-agent static self-review; it is not independent external peer review or a formal proof-assistant certificate.

## Actual dependencies and primary sources inspected

Read the complete saved warp-family problem, commutator/calibration, attainment, conditioning and novelty documents, plus the family decision. Particular attention was given to the full score/tangent definition, influence normalization and the original-class uniform upper bound.

Read all six existing lower-bound v1 documents. Read the historical efficiency argument's original reference, score, norm, compatible-influence definition and unresolved projection problem. The old reference has \(g_1=g_2=e^x,\delta=.1,\theta=1\); the present lower bound does not answer that reference question.

Read the cached Trabs v2 text around Theorem 2.7, and retrieved the primary HTML at:

- [Trabs, arXiv:1307.6610v2](https://arxiv.org/html/1307.6610v2), including Theorem 2.7 and the nonlinear example discussion.
- [Heinrich–Kahn, arXiv:1507.04313v1](https://arxiv.org/html/1507.04313v1), including Assumptions A/B, Theorems 3.2, 3.3 and 3.5, and Section 3.3.

The [comparison](consequences_and_prior_overlap.md) records the precise application and nonapplication of these results. This was bounded primary-source inspection, not exhaustive literature research. No cited source was unavailable.

## Tools and authoring operations

Shell inspection used cat, sed, rg, test, git status and git diff. Selected exact commands, not an exhaustive transcript, were:

~~~sh
cat /home/shadi/.codex/attachments/a29769e7-baad-47b7-b189-09046205cd00/pasted-text.txt
git status --short
rg --files --hidden --no-ignore -g AGENTS.md -g '!.git' /home/shadi/iclr2027
test ! -e docs/interaction_recoverability/warp_family_lower_bound/v2
cat docs/interaction_recoverability/warp_family/v1/commutator_and_calibration.md
cat docs/interaction_recoverability/warp_family/v1/conditioning.md
cat docs/interaction_recoverability/warp_family_lower_bound/v1/admissibility_and_score.md
cat docs/interaction_recoverability/warp_family_lower_bound/v1/{lower_bound,consequences_and_prior_overlap}.md
sed -n '276,337p' docs/interaction_recoverability/warp_family/v1/attainment.md
sed -n '1,85p' docs/interaction_recoverability/efficiency_audit/v1/efficiency_argument.md
sed -n '315,390p' runs/interaction_recoverability/prep-20260912-082154/unknown-profiles-v1/trabs_v2.txt
git diff --stat
git diff --cached --stat
~~~

Standard-library Python was used solely for file hashing, version creation, text copying/replacement and static document checks. Shell authoring used literal quoted heredocs. No project module, mathematical evaluator, numerical library, estimator or diagnostic was executed.

The sandbox launcher failed to initialize its loopback interface on the first read attempts. The patch tool later failed for the same reason before applying its requested changes. Reads and authoring were completed through the approved require_escalated shell route; the failed patch was replaced by literal text authoring confined to v2. These were sandbox-launch failures, not automatic approval-review rejections.

## Preservation record

Before any project write, SHA256 hashes were recorded for every existing file below docs/interaction_recoverability and every repository file matched by the ledger-name search. This covered 220 historical document and ledger files, including both existing warp versions, external-review v2 and the original efficiency audit.

The initial record and git-status text are at:
- /tmp/warp-lower-v2-zkyx88n2/before.json
- /tmp/warp-lower-v2-zkyx88n2/git-status-before.txt

Only the six Markdown deliverables in this fresh v2 directory were authored in the project. Temporary integrity records were written under /tmp. No historical file, learned-target prediction or configuration was an authoring target; no ledger was appended. The integrity scan covers the recorded files, not every artifact elsewhere in the workspace. The unchanged prediction/configuration scope additionally follows from the confined authoring operations, not an asserted exhaustive hash scan.

No git commit or push, external message, review-packet rebuild, publication, archive or backend operation occurred.

## Resource accounting

| Category | This continuation |
|---|---|
| Numerical jobs, sampling and simulations | None |
| Fitting, optimization, training, Gram calculations and benchmarks | None |
| GPU or backend execution | None |
| Computational allowance wrapper or compute-ledger append | None |
| Allowance reset, inferred reserve or allowance spending calculation | None |
| File reads, hashing, static checks, source retrieval and document authoring | Nonzero overhead; not comprehensively CPU-metered |
| Mathematical reasoning | Performed; no fabricated CPU or numerical allowance equivalent assigned |

Hashing and text checks are preservation/authoring overhead, not numerical proof validation. Individual tool wall times do not provide aggregate CPU usage, total research time, model-service energy or exhaustive resource accounting. Zero total computation is not claimed. No overhead debit was invented or written to a historical ledger.

## Final checks

The final preservation and document-check results are recorded below after execution.

All 220 historical-file hashes matched the initial record. The six required v2 Markdown files were present; local Markdown links, math delimiters, terminal newlines and trailing whitespace passed the static checks. Tracked and staged git diffs remained empty, and git diff --check returned success. Preexisting untracked research files remain present; the working tree is not claimed clean.

The check output is saved at /tmp/warp-lower-v2-zkyx88n2/final-checks.json, with the final git-status snapshot alongside it. The checker was a standard-library text/integrity script and executed no project or numerical code. These mechanical checks establish document consistency and preservation within their stated scope; the mathematical verification consists of the analytic derivations.

This result paragraph was appended only to this new resource document after the checks. Every historical decision, theorem, review archive and ledger remains unchanged.
