# Protocol v2 — verifier evaluation

Supersedes the verifier half of `preregistration.md` for the revision cycle.
**Protocol v1 gates G1–G4 are retained unchanged in `docs/corrected_ledger.json`
and are never re-scored here.** v2 does not reinterpret them; it replaces the
metric for future comparisons only.

## Why the v1 metric was replaced

v1's primary quantity was a **median relative gap reduction** against the
unconditional baseline, with a 25% minimum meaningful gain. On the v1 families
the plain box baseline already closed **99.1%** of instances to gap exactly zero
within the matched second, so the metric was `0/0` and G1/G2 were correctly
recorded UNDETERMINED rather than FAIL.

The defect was the **instance distribution**, not the threshold. Transplanting
25% onto a different metric would have been unjustified, so v2 declares new
quantities and new thresholds together.

## The three declared quantities

**Primary — decision coverage.** The fraction of instances where an arm
*certifies* the declared decision within the shared budget. The decision is
"is `max_{z∈F} f(X,z)` below `tau`?", and it is settled **either way**: by
proving the optimum below `tau` (`upper < tau`), or by exhibiting a feasible
point at or above it (`incumbent >= tau`). An arm that resolves nothing scores 0.

`tau = incumbent + margin · |incumbent|`, with the margin drawn from the frozen
grid `{0.02, 0.10, 0.50}` and the incumbent **shared across arms**, so every arm
faces the identical decision on each instance. Small margins are near-boundary.

**Secondary — capped time to resolution.** The time at which the gap first falls
to tolerance, read from a `(t, L, U, status)` trace recorded at every bound
update. Unresolved instances are recorded **at the cap and counted**, never
dropped. This is not `wall_s`: v1's `wall_s` is a *stopping* time and equals the
budget for a timed-out record.

**Secondary — anytime gap AUC.**

```
AUC_i(T) = (1/T) ∫₀ᵀ min(1, gap_i(t) / s_i) dt
```

`s_i = max(initial box gap, 1e-6)` is a **common per-instance scale shared by
every arm**, so AUCs are comparable rather than each arm being graded on its own
scale. Before an arm produces any valid bound the integrand is **1**, not 0 —
substituting zero there would reward silence. A negative gap is recorded as a
**defect** and clamped only inside the integral, never repaired in the record.

## Rules that bind the comparison

| Rule | Mechanism |
|---|---|
| The generation rule is frozen before any case is drawn | `scmp/revision/ladder.py`, hashed as `RULE_HASH`; the hash is written into every case record and every result file |
| No case is added, dropped or reweighted after an arm has run | the rule is total and deterministic; `generate()` returns the same set every time |
| Easy cases are never discarded | rung 0 is the preserved v1 easy set and appears in the **aggregate** as well as its own stratum |
| Identical budgets | same wall-time cap and thread count for every arm |
| All timeouts and unresolved results recorded | with their last valid bound and an explicit reason |
| Soundness checked against ground truth | exhaustive enumeration gives the exact optimum, so each arm's **true** gap is computable, not just relative comparisons |

## The difficulty ladder

Frozen rule hash **`580b1c5e3ab6f792`**. Four rungs graded along the five axes
the plan names — candidate-edge count, family compatibility, supported
depth/width, unstable-ReLU count, and decision margin.

| Rung | Name | Families | d_hid | Candidate edges |
|---|---|---|---|---|
| 0 | `v1_easy_preserved` | path_matching:6, layered_dag:3×3, rna:GGAAACC | 4 | 4–18 |
| 1 | `wider_candidate_set` | path_matching:8, layered_dag:3×4, rna:GCGCAAAAGCGC | 6 | 7–32 |
| 2 | `looser_compatibility` | path_matching:10, layered_dag:4×4, rna:GGGGAAAACCCC | 8 | 9–48 |
| 3 | `wider_network_more_unstable_relus` | path_matching:10, layered_dag:4×5, rna:GGCGCAAAAGCGCC | 12 | 9–75 |

108 cases: 4 rungs × 3 families × 3 seeds × 3 margins.

**The ladder solves the v1 problem.** Median box gap on the ladder is **436.75**,
where v1's baseline closed 99.1% of instances to exactly zero. There is now
something for conditioning to improve.

## External references, and an honest gap

**No topology-aware neural-network verifier is installed in this environment** —
`auto_LiRPA`, `gurobipy` and `cvxpy` are all absent. That is recorded as a gap.
No internal baseline is relabelled as external.

Two genuine independent references are used instead:

- **Exhaustive enumeration** — the exact optimum by pruned DFS over the family,
  with infeasible prefixes cut. Stronger than any solver where it applies, and it
  makes every arm's *true* gap computable. Capped at 4000 members and 22
  candidate edges; beyond that the case is reported without ground truth rather
  than with a guess.
- **LP over a degree-≤1 outer relaxation** (`scipy.optimize.linprog`) for the
  additive term, plus the box envelope on the residual. Degree-at-most-one is
  valid for all three families here, so this is a sound outer relaxation and a
  genuinely different bounding strategy from our oracle — looser by construction,
  which is the point.

## Instrumentation

Seven components timed per instance, with call counts and **actual work**:
chart compilation, bound propagation, support queries, cache access (hits and
misses separately, with the work each did), policy inference, incumbent search,
branching. A test asserts an instrumented call returns exactly what an
uninstrumented one returns.

**Cache accounting.** A hit that avoids an oracle call is cheap; a miss that
triggers one is not. v1 charged both alike, which is how a budget of 8 recorded
80 passes.

## Policy gate

A learned policy is trained **only if** deterministic selective conditioning
shows a useful cost-benefit region under the quantities above. Otherwise no
policy is trained and the reason is recorded. When trained, its objective is
measured downstream resolution gain per unit measured cost — not one-step gap
reduction — and three constraints hold regardless:

- policy outputs choose **computations**, never bound values;
- branch coverage and pruning validity do not depend on the policy;
- teacher generation, training, inference and a break-even workload are charged
  and reported.
