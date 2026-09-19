# P2 — verifier arms under protocol v2

Record: `runs/revision/verifier_pilot-20260905-233716.json` (+ `.sha256`).
Profile: `runs/revision/profile_verifier-20260905-225029.json`.
Protocol: [`docs/protocol_v2.md`](protocol_v2.md). Ladder rule hash
**`580b1c5e3ab6f792`**, frozen before any case was drawn. 108 cases × 7 arms =
756 runs, 1 thread, CPU only, **0.038 core-hours**.

**v1 gates G1–G4 are retained unchanged and are not re-scored here.**

## Headline

Conditioning now **wins** on the primary metric, where under v1 it lost at
matched wall time. And the advantage **grows monotonically with difficulty**:

| Rung | box | fixed_conditional | gain | deterministic_selective |
|---|---:|---:|---:|---:|
| 0 (v1 easy, preserved) | 59.3% | 59.3% | **+0.0** | 59.3% |
| 1 | 59.3% | 63.0% | **+3.7** | 63.0% |
| 2 | 40.7% | 51.9% | **+11.1** | 48.1% |
| 3 | 11.1% | 25.9% | **+14.8** | 25.9% |

This is the explanation for v1's UNDETERMINED result, made concrete: on rung 0 —
which *is* the v1 distribution — conditioning gains **exactly nothing**. It was
never a defect of the mechanism; the instances could not discriminate it.

## Overall, all 108 cases (rung 0 included, never dropped)

| Arm | Coverage | AUC ↓ | Wall s | Oracle calls | True gap ↓ |
|---|---:|---:|---:|---:|---:|
| box | 42.6% | 0.945 | 0.624 | 1 | 153.6 |
| terminal_only | 32.4% | 0.945 | 0.612 | 1 | 1495.5 |
| intermediate_unconditional | 42.6% | 0.945 | 0.599 | 1 | 153.6 |
| **fixed_conditional** | **50.0%** | **0.736** | 0.587 | 3760 | **106.7** |
| **deterministic_selective** | 49.1% | 0.812 | **0.328** | **1** | 115.2 |
| random_selective | 49.1% | 0.736 | 0.304 | 3760 | 107.4 |
| lp_polytope_external | 24.1% | 0.972 | 0.053 | 0 | 1504.3 |

`True gap` is measured against the **exact optimum by enumeration**, so these are
absolute distances from ground truth, not relative comparisons.

**Selective conditioning is the efficiency result.** It reaches 49.1% coverage —
within 0.9 pp of always conditioning — at **half the wall time** and a *median of
one oracle call* instead of 3760. Its action census: `condition` 51,
`skip_few_unstable` 51, `skip_already_resolved` 6. Skipping roughly half the
instances costs almost no coverage, which is precisely the region v1's
always-condition arms could not exploit.

`terminal_only` (separate maxima) has a true gap of 1495.5 against 153.6 for the
combined bound, confirming the combined form on a distribution where it matters.

## Soundness

**0 violations against the exact optimum**, across all arms and all 108 cases.

An earlier run of this same pilot reported **66 violations**. That was **my
ground truth, not the bounds**: `feasibility(S)` means *extendable to a member*,
which equals *is a member* only for downward-closed families. Layered-DAG paths
are not downward-closed, so partial paths were being evaluated and the "exact"
maximum was inflated above sound bounds. Membership is now tested by pinning the
complement, and `members_evaluated` equals the oracle's independent family count
exactly for every family (27=27, 256=256, 625=625, 34=34, 26=26). Both facts are
asserted in `tests/revision/test_verifier_p2.py`.

**42 negative-gap defects**, all at **−8.5e−14**, identical across all seven arms
on cases where the bound is tight. Roundoff in the shared float evaluation, not
arm logic — the same signature P0 found in `policy_pilot`. Recorded, never
clipped. This is what P4's enclosed arithmetic must make impossible.

## Where the time goes

From `profile_verifier`, before any optimisation:

| Arm | median wall | median gap | oracle calls | cache hits/misses |
|---|---:|---:|---:|---|
| box | 0.0082 s | 436.75 | 1 | 0/0 |
| fixed_conditional | 0.2263 s | 277.92 | 3766 | 2568/200 |

Cost is **oracle support queries**, not cache or chart compilation — the cache
already runs at a 92.8% hit rate, so its 200 genuine misses are the only work
there. Per the plan's rule, profiling does **not** justify batching or
chart-reuse work; the lever is *how often we query*, which is what selective
conditioning acts on.

## Limitations, stated

- **Time-to-resolution is near-vacuous here.** 102 of 108 cases remain unresolved
  at the 5 s cap for every arm, so the median is the cap throughout. The ladder
  is hard enough that gaps rarely reach 1e-9. Coverage and AUC discriminate; this
  secondary metric does not, and should not be read as an arm difference.
- **No external topology-aware NN verifier** is installed (auto_LiRPA, Gurobi,
  cvxpy all absent). Recorded as a gap; no internal baseline was relabelled as
  external. The LP arm is a genuine independent bounding strategy but a
  deliberately loose one — a degree-≤1 outer relaxation that ignores non-crossing
  and layer structure, which is why its true gap is 1504 against 154.
- **The frozen predictor is the best available, not a demonstrated-useful one.**
  P1 showed a pairwise ridge beats it outright. Everything here is evidence about
  *verification*, and does not support the predictive claim.
- `random_selective` matches `deterministic_selective` on coverage (49.1%) while
  spending 3760 oracle calls. The *selection rule* is what saves cost; the
  *scoring* of which edge to condition is not yet shown to matter.

## Policy gate: OPEN

Criterion: deterministic selective must beat the best cheap arm on decision
coverage **or** anytime AUC at equal budget. It does — **49.1% vs 42.6%**
coverage, and 0.812 vs 0.945 AUC, at roughly half the wall time.

Under the plan a policy may now be trained. **It has not been trained yet.** When
it is, three constraints hold regardless of outcome: policy outputs choose
computations and are never bound values; branch coverage and pruning validity do
not depend on the policy; and teacher generation, training, inference and a
break-even workload are charged and reported. Its objective must be measured
downstream resolution gain per unit measured cost, not one-step gap reduction.

The open question the policy must answer is narrow, and `random_selective` sets
it up: since random scoring already matches deterministic scoring on coverage,
a learned policy must earn its keep on **when to query**, not **which edge**.
