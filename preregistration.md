# Preregistration — certified verification for structured-family message passing

Written **2026-09-05, before any large experiment**, from development-set evidence only.
No number below was chosen by looking at a final test result, and the final test split is
not loaded by any script in this repository.

Machine-readable companion: `configs/protocol.yaml`, validated by
`python -m scmp.validate_protocol`.

---

## 1. What is being claimed, and what is not

The claim is about a **verifier**, not about biology.

> **H1 (primary, verifier-only).** At a matched query budget and on a frozen model,
> conditional message support produces a strictly tighter certified upper bound than the
> intermediate-unconditional baseline.

> **H2 (secondary, predictive).** A sign-constrained certifiable predictor is
> *noninferior* to an unconstrained predictor of the same architecture, within a
> prespecified margin, when both use an expressive non-negative encoding.

**Not claimed.** That the certified bound bounds siRNA efficacy, off-target risk, or any
biological quantity. Section 6 explains why that claim is currently unsupportable
regardless of how the experiments come out.

---

## 2. Primary metric and effect sizes

**Primary metric.** Median *relative certified gap* at a fixed budget

```
gap_rel = (U - incumbent) / max(|incumbent|, 1)
```

`U` is the certified upper bound; `incumbent` is the best value actually attained by the
model at a feasible point. Both come from the same run. Lower is better; zero means the
instance was closed.

**Secondary metrics.** Fraction of instances proved closed within the node budget; total
oracle calls; wall-clock seconds; number of instances returning `unresolved`.

**Minimum meaningful gain (MMG): a 25% relative reduction in median `gap_rel`.**
Justification from development data, not from the test set:

| development observation | source | value |
|---|---|---|
| conditional vs unconditional baseline, median gap | `audit_conditional`, 500 cases | 0.535 → 0.083 (−84%) |
| random vs largest-gap edge choice | same run | 0.0830 vs 0.0830 (≈0%) |
| measured-gain vs random | same run | 0.0806 vs 0.0830 (−3%) |

The scorer differences sit at the noise floor and the conditioning effect is far above
it. 25% is set deliberately between the two: comfortably larger than anything we have
seen from noise, comfortably smaller than the effect we expect to confirm. A result below
MMG will be reported as **no meaningful gain**, not as a trend.

**Predictive noninferiority margin: 0.05 Spearman.** Justification:

| development observation | source | value |
|---|---|---|
| sign-constraint cost, positional encoding | `f8_baseline`, 5 seeds | 0.002 |
| sign-constraint cost, gene-disjoint split | same | 0.039 |
| seed-to-seed standard deviation | same | 0.004 – 0.023 |

0.05 exceeds the largest observed seed spread and the largest observed cost, so
noninferiority at this margin is a real test rather than one that passes automatically.
Declared **only** if the upper end of the 95% CI for (unconstrained − constrained) lies
below 0.05.

---

## 3. Design

**Independent split unit.**
- RNA: **gene**. The published Huesken 2182/249 split shares all 30 annotated genes
  between train and test (`f6_trained`, `genes_shared: 30`), so it is not independent.
  Every reported RNA number uses gene-disjoint splits; the published split may appear
  only as a clearly labelled secondary row.
- Non-RNA families: **family instance** (the graph topology). Instances are generated
  from disjoint seed ranges for development and test.

**Seeds.** 5 seeds `{0,1,2,3,4}` for every arm. Every seed is reported. No seed is
dropped for any reason, including divergence — a diverged seed is reported as diverged.

**Confidence intervals.** Cluster bootstrap over the independent split unit, percentile
method, 10,000 resamples, 95%. Paired arms are bootstrapped on the paired difference.
Point estimates without an interval are not reported.

**Tuning limits.** At most **12 hyperparameter configurations per arm**, chosen on the
development split. The same 12-configuration allowance applies to every baseline,
including the ones we expect to lose. No arm may receive more tuning than another; the
count is recorded per arm and checked by `validate_protocol`.

---

## 4. Comparison arms

All verifier arms run against the **same frozen model weights**, verified by model hash.
A run whose arms disagree on model hash is void.

### 4a. Verifier-only (model frozen; no training in the comparison)

| arm | description | status |
|---|---|---|
| `box` | box concretisation of the affine bound | implemented |
| `terminal_only` | exact support at the terminal only; separate maxima | implemented |
| `intermediate_unconditional` | affine propagation, combined support | implemented |
| `fixed_conditional` | conditional support, deterministic edge choice, fixed budget | implemented |
| `learned_conditional` | conditional support, learned query scorer | implemented |
| `topology_aware_solver` | matched MIP/solver method in the style of topology-based bounds tightening | **not implemented**; requires SCIP or an equivalent. Reported as *unavailable* unless it is built and matched, never omitted silently |

`topology_aware_solver` is listed because leaving it out would flatter us. If it is not
built before the main run, the paper must say the comparison was not made.

### 4b. Predictor-training (separate table, never merged with 4a)

| arm | description |
|---|---|
| `monotone_positional` | sign-constrained, one-hot (position, nucleotide) |
| `unconstrained_positional` | identical architecture, signed weights |
| `monotone_generic` | sign-constrained, length-free encoding |
| `unconstrained_generic` | identical architecture, signed weights |

Predictor arms and verifier arms answer different questions and are **never combined into
a single headline number**.

### 4c. Cost accounting — everything is charged

Reported per instance and summed per arm:
`bound_construction_s`, `oracle_calls`, `policy_inference_s`, `teacher_generation_passes`,
`training_s`, `wall_clock_s`. A learned arm carries its teacher-generation cost in every
table in which it appears. An arm's cost may not be reported net of anything.

---

## 5. Inclusion and exclusion rules

Fixed in advance. Exclusions are applied **before** any metric is computed and are
reported with counts.

**Exclude (void, reported separately, never silently dropped):**
1. Empty family under the mask — nothing to certify.
2. Singleton family — the bound is trivially exact and carries no information.
3. Model output constant across the family (dead instance) — "tightness" is vacuous.
4. Any non-finite value anywhere in the pipeline.
5. Dataset rows failing `certmp.data.verify` integrity checks.

**Never a reason to exclude:** the instance was hard, the bound was loose, the arm lost,
the seed was unlucky, or the result was not significant.

**Unresolved is a result.** An instance hitting the node or time budget is reported as
`unresolved` with its valid-but-open bound. Unresolved instances are counted in every
denominator.

---

## 6. RNA label semantics — the blocking issue

**The label is not indexed by structure, but the model is.**

Huesken et al. measure one number per siRNA sequence: `%Inhibition` under a fixed assay
(the annotated table records cell line, concentration and duration). The model
`f(X, z)` produces a different value for each structure `z`. So:

- Training on `(X, z_MFE) -> y` silently asserts that the label is attributable to the
  minimum-free-energy structure. That assertion is exactly what this project set out to
  question.
- **We must not assume the label is invariant across latent structures.** If `y` really
  were invariant in `z`, the correct model would be constant in `z` and structural
  certification would be vacuous. If `y` is not invariant, then a single `(z, y)` pair
  per sequence is a *structure-marginal label fitted by a structure-conditional
  function* — a mismatch, not a modelling choice.
- `max_z f(X, z)` has **no measured counterpart**. No experiment in this dataset reports
  efficacy under a worst-case structure, so the certificate cannot be validated against
  data, only against the model.

**Consequences, binding on the write-up:**

1. The primary hypothesis H1 is verifier-only and is unaffected: it compares bounds on a
   frozen model, and needs no biological label at all.
2. Any predictive claim (H2) is a claim about fitting the *marginal* label. It may be
   reported as such, and may not be described as predicting structure-conditional
   efficacy.
3. The phrase "bounds off-target risk" is retired. `README.md` currently uses it and is
   wrong on two counts: the label is on-target inhibition, and it is marginal over
   structures (`docs/claim_ledger.md`, D1).
4. Before any structure-conditional predictive claim, one of these must hold: a
   structure-resolved label set exists; or the model is trained on an ensemble-averaged
   objective `E_{z~Boltzmann} f(X,z) ≈ y` and the claim is stated for that average.

---

## 7. Data licences

Recorded as found. Unresolved entries block redistribution, not internal use.

| source | used for | licence status |
|---|---|---|
| GENCODE release 47 transcripts | target-site windows | Page states all GENCODE data is "open access"; **no explicit licence** (CC-BY/CC0 etc.) is given on the data-access page. **Unresolved** — check FTP README or contact GENCODE before redistributing any derived file |
| Huesken et al. 2005 via DSIR redistribution (biodev.cea.fr) | efficacy labels | **Unresolved.** Terms not established. The primary paper is a Nature Biotechnology publication; the redistribution's terms were not stated at the download point |
| siRNAEfficacyDB (gene annotation join) | gene symbols for split units | **Unresolved.** Terms not established |
| ViennaRNA 2.7.2 | folding, sampling | Software dependency, not data. Licence not audited here |

**Binding rules.** No third-party data file is vendored in this repository — all are
fetched and gitignored, which is already the case. No derived data file is published
until its row's status is resolved. "Open access" is not a licence and is not treated as
one.

---

## 8. Experiment order — two non-RNA families first

Deliberate: if the bound machinery has absorbed an assumption about RNA, running it on
families that share none of that structure is what exposes it.

1. **`layered_dag_path`** — source-to-sink paths in a layered DAG. Every member has the
   same cardinality; the production-to-edge map is the identity.
2. **`path_matching`** — matchings on a path graph. Downward closed but *not* union
   closed, which is the configuration C5 proved has unbounded relaxation slack. This is
   the adversarial case for the whole approach and is run second, deliberately.
3. **RNA pilot, tiny** — a small number of short target windows, purely to confirm the
   pipeline runs end to end on the real family. Not powered for any hypothesis.

A failure on family 1 or 2 halts the RNA work; RNA is the application, not the evidence.

---

## 9. Compute

**No compute budget has been authorised, so no substantial run is launched.**

Proposed capped pilot, for approval:

| resource | cap | note |
|---|---|---|
| CPU | 8 core-hours | the whole pipeline is CPU-only; oracles and bounds are numpy |
| GPU | 0 | not required for the pilot. If later needed, `scmp/gpu.py` pins the least-loaded device at ≤15% memory |
| wall clock | 2 hours | |
| instances | ≤ 600 total across all three families | |
| paid/cloud compute | none | requires separate authorisation |

`scmp.run_pilot` refuses to run without `--dry-run` unless `authorized_budget` is
explicitly set in `configs/pilot.yaml`. The default is `null`.

---

## 10. Amendments

Appended with date and justification; nothing above this line is edited. Any change made
after seeing a test-set result must say so explicitly and the affected claim is
downgraded to exploratory.
