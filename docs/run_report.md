# Certified Structured-Family Verification — Run Report

**5 September 2026 · repository `iclr2027` · 42 logged runs**

Every experiment run in this repository, what it measured, and what it does not
establish. The verification machinery is sound and reproduces exactly. The two
mechanisms the project was built to demonstrate are not supported by the
evidence, and the RNA application is blocked on a construct-identity fact
discovered late.

| | |
|---|---|
| Tests | 703 passed · 8 skipped |
| Result files | 202 across 42 runs |
| New code | 48 files · 10,143 lines |
| Soundness violations | 0 in 1,262,952 checks |
| Preregistered gates | 1 pass · 2 undetermined · 1 fail |
| Compute spent | 0.239 of 8 authorised core-hours |

---

## 1. Verdict

Three things are true at once. This report is ordered so none of them hides the
others.

**The machinery is sound.** Zero soundness violations across 1,262,952
intermediate bound checks, 500 randomised conditional-support cases, 2,525 exact
oracle checks, 2,222 certificate obligations and 75 model-effect intervals.
Every soundness claim is checked against exhaustive enumeration or an
independently written checker — never against itself.

**The mechanisms are not demonstrated.** Conditional message support is
UNDETERMINED at matched wall time: the plain box baseline already closes 99.1%
of pilot instances to gap zero, so the relative-reduction metric is 0/0 — and
where the arms do differ, conditioning is behind. Learned allocation sits at the
noise floor. The nonlinear residual *hurts* held-out prediction by 0.268
Spearman, CI [−0.408, −0.128].

**The RNA application is blocked.** Huesken 2005 measured siRNA activity against
YFP reporter plasmids in H1299 cells, not native transcripts — contradicting the
redistributed annotation this project had relied on. Six of eight comparison
arms need a target window that does not exist for this dataset. An earlier phase
of this same project had already folded native GENCODE windows: exactly the
silent substitution now refused.

What survives is a general result about relaxations of combinatorial families,
plus a verification stack that is explicit about what it does not know. That is
a smaller paper than intended, and a real one.

---

## 2. How to read the numbers

Results are kept in three classes and never mixed. A result in A or B is **not**
evidence for a claim in C.

| Class | What it is | How it is checked |
|---|---|---|
| **A** | Supporting lemmas — no trained parameters | exhaustive enumeration, exact integer arithmetic, external oracles |
| **B** | Engineering checks — no trained parameters | independent enumerators and an independent certificate checker |
| **C** | Trained empirical results — fitted parameters | preregistered gates, cluster bootstrap CIs, all seeds reported |

A second axis runs through everything: `numerical_status`. Exact-integer
quantities are compared with `==`; floating-point quantities carry a tolerance
and are labelled diagnostics. Every bound in this project is currently
`float_unverified` — see §9.

---

## 3. Phase one: the monotone certificate (`certmp`, F1–F10b)

A pre-registered study of certified reachability for monotone message-passing
networks. Predictions S1–S6 were registered before any run. The package was left
untouched for the rest of the project — `git diff HEAD` is empty at every stage
that follows.

| Run | Question | Result | Verdict |
|---|---|---|---|
| F1 | Do reachable extremes sit at lattice endpoints? | 6/6 cases, theorem holds — 2 forward passes replace 2^k | CONFIRMED |
| F2 | Is the lattice non-trivial at realistic length? | median k = 16 / 48.5 / 86 at 40 / 80 / 150 nt; lattice 10^4.8 → 10^25.9 | NOT VACUOUS |
| F3 | How far is MFE from the reachable extreme? | median relative underestimate 1.109 | MEASURED |
| F3b | How much of the gap is relaxation slack? | median slack ratio 1.819; slack is **87.9%** of the reported gap; upper bound sound on all samples, **lower bound failed twice** | SPLIT |
| F4 | Can top-k selection be certified? | margin negative at every k (−170 at k=1, −5024 at k=10) | RETIRED |
| F5 | Does the theorem survive stress? | 800 trials, **0 violations**; max-aggregation goes void as mandatory pairs rise (live 1.00 → 0.00 by 23) | HOLDS |
| F6 | Does a trained monotone model predict? | Spearman 0.5916, Pearson 0.5564 against a published Pearson baseline of 0.66; gene-disjoint 0.5432 | TRAINED |
| F7 | Is the certificate sound at probability floor zero? | sound at floor 0, no mandatory set | SOUND |
| F8 | What does certifiability cost against a matched baseline? | positional features Δ = **0.0025**; generic features Δ = **0.214**; gene-disjoint Δ = 0.039 | MEASURED |
| F9 | Does it hold on real mRNA target sites? | 1,720 sites, 28 genes; k = 407 / 1,743 / 3,931; lattice to 10^1183; coverage 1.0, sound throughout | INVALIDATED |
| F10 | Is the maximal-element route cheaper? | M95 = 2,518 of 2,956 at 150 nt; enumerating 95% of mass buys a ratio of 1.045 | REJECTED |
| F10b | Does M95 converge with more samples? | growth exponent ≈ N^0.80 — no saturation | REJECTED |

F4 was withdrawn, not tuned: exact top-k certification needs a sound ensemble
*lower* bound, and F3b had just shown the lower bound violated on 2 of 20
windows. The unsound bound was removed rather than patched. F9 was sound at the
time and is invalidated retrospectively by the assay-identity finding in §13 —
the windows are native transcript context for an assay that used a reporter
construct.

> **The most useful negative.** F8 splits cleanly. With *positional* features,
> certifiability costs 0.0025 Spearman — essentially free. With *generic*
> features it costs 0.214. The constraint is not inherently expensive; it is
> expensive when the feature map cannot express what the signed weights were
> doing.

---

## 4. Corrections C1–C9

Nine corrections applied after the first pass. C1–C6 generalise the
characterisation and re-examine the slack; C7–C9 were additive only.

| Run | What it established | Value |
|---|---|---|
| C3 | Endpoint exactness holds iff the aggregator is monotone under multiset inclusion — either direction, isotone or antitone. 7 aggregators × 25 seeded trials | holds |
| C4 | Sum + ReLU under non-negative parameters is *exactly* affine, so nonlinearity cannot be inferred from weight norms | 2.16e−16 |
| C5 | Relaxation-gap exponent per layer, depths 1–4 (R² 0.9996 at every depth) | 0.912 |
| C5 | Gap is unbounded on downward-closed families that are not union-closed (m = 3 → 48) | 9.0 → 2290.3 |
| C5 | Gap is exactly 1.0 on join-closed families, brute-force verified | 1.0 |
| C6 | Closure property holds on 54 rows, and the control discriminates — so the check is not vacuous | holds |
| C7 | Relaxation gap equals a walk-count ratio κ_L: max relative error against exhaustive enumeration | 2.12e−16 |

> **A correction that changed an interpretation.** C5 overturned an earlier
> reading. The slack exponent had been attributed to the biology; it is
> *architectural*. It is 0.912 per layer at every depth from 1 to 4, on a
> synthetic family containing no application data at all.

---

## 5. κ_L against measured slack  *(class A)*

C8–C9 tested the characterisation where it matters: against slack actually
measured on real target-site windows.

| Window | Measured slack | κ₁ | κ₂ | Bounded? | κ₁² vs κ₂ |
|---:|---:|---:|---:|---:|---:|
| 50 nt | 21.5 | 10.69 | 118.50 | yes | 3.5% err |
| 100 nt | 94.3 | 21.95 | 494.12 | yes | 2.5% err |
| 150 nt | 222.8 | 32.37 | 1098.03 | yes | 4.6% err |

Two things at once. κ₂ bounds the measured slack at every length — the
characterisation is not merely a synthetic identity. And κ₁² predicts κ₂ to
within 2.5–4.6%, which is the per-layer multiplicativity showing up in real
windows rather than in a fitted curve. All three converged by 2,000 samples.

---

## 6. Phase two: what was built (`scmp`)

A new package, stage by stage, each stage gated on its own validation before the
next was written.

| Layer | Module | What it provides |
|---|---|---|
| Support oracles | `oracles/` | non-crossing pairings, layered-DAG paths, path matchings; exact count, support, log-partition, marginals, conditioning, feasibility, valid backtrace |
| Model | `model.py` | additive anchor plus a signed-weight nonlinear residual |
| Bounds | `bounds.py` | affine envelopes for signed maps, ReLU relaxation, McCormick product envelope, box and family concretisation |
| Conditioning | `conditional.py` | per-edge bounds conditioned on the edge being present, with a query budget, cache and scorer |
| Search | `refine.py` | branch refinement, covering splits, replayable certificates |
| Checking | `check_certificates.py` | independent checker that imports nothing from the search it checks |
| Numerics | `numerics.py` | validated interval arithmetic via `nextafter`, checked against exact rationals |
| Application | `rna/` | covariates, isoform-explicit windows, ensemble priors, marginals, sampling, five audits |

```
f(X, z) = b(X) + a(X)^T z + gamma * g(X, B + Pz)
          ^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^
          exactly optimisable by     only this
          the oracle                 needs relaxing
```

At **γ = 0** the bound is exact by construction, and the bypass is verified
rather than assumed: the residual is filled with NaN and f must stay finite. The
anchor is structurally prevented from seeing z or the conditioning mask —
enforced by signature inspection.

On a stacking toy with feasible non-additive interactions the residual improved
held-out error on **5 of 5 seeds** against a matched γ = 0 ablation, with
best-of-seeds train MSE 0.27 against an unconstrained additive floor of 1.19. It
did *not* improve train error on 2 of 5 seeds, and that limitation is asserted
in a test so it cannot quietly disappear.

---

## 7. Oracles  *(class A — exact)*

`audit_oracles --max-n 10 --cases 200 --seed 7027` — 2,525 exact-integer checks
and 1,379 float diagnostics, in separate counters, so a float agreement is never
counted as exactness.

| Check | Scope | Result |
|---|---|---|
| Counts vs independent enumerator | 9 tiny sequences × 5 min-loop/canonicity settings | exact match |
| Counts vs Motzkin numbers (OEIS A001006) | n = 0…10, external oracle | exact match |
| Grammar unambiguity | an ambiguous grammar would count more than the family holds | holds |
| Production-to-edge map P | injectivity, edge multiset equals the structure, score reconstructs | holds |
| Every member selectable | drive the objective at each member; argmax must be that member | holds |
| All masks of size ≤ 2 forced × ≤ 2 forbidden | consistent and contradictory alike, against enumeration | exact match |
| Chart values / materialised productions | closed forms (n+1)(n+2)/2 and C(n+2,3) | verified exactly |

Those closed forms replaced a log-log fit that came out at n^1.57 and n^2.57 —
biased low by lower-order terms at this range. The exact identities are the
claim; the fit is reported as a diagnostic only.

---

## 8. Bounds  *(class B — 0 violations)*

`audit_bounds --config configs/exhaustive.yaml` enumerates the entire feasible
family for every model and checks the bracket at *every* intermediate, not only
at the output.

| Quantity | Value |
|---|---:|
| Cases (models × families) | 192 |
| Feasible structures visited | 5,856 |
| Intermediate bracket checks | 1,262,952 |
| Violations | **0** |
| Envelope inequalities formed and tested individually | 53,280 |
| ReLU units — stable⁺ / stable⁻ / unstable | 5,634 / 7,158 / 1,608 |
| γ = 0 cases exact | 48 / 48 |

Layer 1 needs **no relaxation at all**: with Z₀ constant in z the message is
exactly affine. Bilinearity appears only at layer 2, where all four McCormick
inequalities are formed. Edge variables are shared across layers — the proof
shows duplication would still be sound but strictly looser.

```
U = b + gamma*d + support(a + gamma*c)
    one combined support call, not two separate maxima
```

The combined form beats the separate-maxima baseline on median gap **0.404
against 1.061**, strictly tighter on 92 of 192 cases. It is not strictly tighter
on the rest, and that is not a failure: when every coefficient is negative the
family optimum is the empty structure and both forms coincide.

---

## 9. Conditional message support  *(class B — sound, unproven)*

The layer-2 term z_e·h contributes only when z_e = 1, so the relevant bound on h
is its maximum over graphs that actually contain e.

```
q_e := max { d^T W_e h_u(z) : z in F, z_e = 1 }
then  z_e * h(z) <= z_e * q_e   for every z in F
```

| Measurement | Value |
|---|---:|
| Randomised cases audited (12,173 feasible structures) | 500 |
| Soundness violations | **0** |
| Median relative gap, baseline → conditional | 0.535 → 0.083 |
| Tighter than baseline | 324 / 500 |
| Edges proved impossible, then ignored | 1,684 |
| Unconditional fallback won at top level | 15 |

Those 15 cases matter. The conditional pipeline concretises over the family
while the baseline uses the box, and the two relaxations are not pointwise
comparable — one can be looser than the other. Retaining the minimum at the top
level makes domination hold by construction rather than by hope. An edge is
ignored only after the oracle *proves* it impossible.

> **The caveat that decides the contribution.** That 6× tightening is at matched
> *query budget*. At matched *wall time* it reverses — see §12. Requiring
> matched wall time is what exposed this; matched budget alone would have read
> as a clean win.

---

## 10. Refinement and certificates  *(class B — 72/72)*

Splits cover the parent exactly, child cap is min(parent, new), fairness is
FIFO, pruning fires only when a valid bound cannot beat an *attained* incumbent,
and elimination requires the oracle to prove the child empty. At timeout the
search returns `unresolved` — never a guess.

| Check | Result |
|---|---:|
| Certificates produced (49 proved, 23 unresolved) | 72 |
| Accepted by the independent checker | **72 / 72** |
| Obligations discharged | 2,222 |
| Tamper tests rejected — cap, incumbent, elimination, stale hashes | 4 / 4 |

The checker does not import the search it checks. It re-derives coverage, caps,
eliminations, prune conditions and bound monotonicity from the record stream
alone, and on these families it also enumerates the family, so "internally
consistent" becomes "correct".

---

## 11. Numerics  *(`float_unverified`)*

Every bound in this project carries the status `float_unverified`. Validated
interval arithmetic exists, is tested against exact `Fraction` arithmetic across
70 tests, and is **not yet wired into the bound propagation**. A test asserts
that refinement cannot upgrade its own status.

> **Not epsilon padding.** Each elementary operation gets a rigorous enclosure
> via `math.nextafter`. One test makes the difference concrete: on
> `[1e16, 1.0, -1e16]` the derived enclosure is *wider* than 1e−9, so a fixed
> epsilon would have been an under-estimate, not a safety margin.

Tests currently assert at an *assumed* tolerance of 1e−9. Promotion to
`float_enclosed` requires the exhaustive audit to pass at slack **0.0** against
widened bounds; four staged criteria are written down in `docs/bounds_proof.md`.

---

## 12. The pilot and the preregistered gates  *(class C)*

Thresholds were fixed on development data before the pilot ran: a 25% minimum
meaningful gain, chosen between the observed noise floor (0–3%) and the observed
conditioning effect (84%); and a 0.05 Spearman noninferiority margin, larger
than the biggest observed seed spread.

| Gate | Point | 95% CI | Threshold | Status |
|---|---:|---:|---:|---|
| G1 conditional vs unconditional | — | — | 0.25 | **UNDETERMINED** |
| G2 learned vs deterministic allocation | — | — | 0.25 | **UNDETERMINED** |
| G3 predictor noninferiority | −0.3361 | [−0.4858, −0.1869] | 0.05 | **PASS** |
| G4 residual predictive value | −0.2680 | [−0.4076, −0.1284] | > 0 | **FAIL** |

**G1 and G2 are undetermined, not passed and not failed.** Both arms reach
median gap exactly 0.0, so the relative-reduction metric is 0/0 and reporting it
as FAIL would have been as wrong as reporting a pass. The families are simply
too easy to discriminate the arms.

### Why: verifier arms at matched wall time

| Arm | n | Median gap | Mean gap | Closed | Timeout | Median wall s |
|---|---:|---:|---:|---:|---:|---:|
| box | 330 | 0.0000 | 0.0008 | **99.1%** | 0.9% | 0.060 |
| intermediate_unconditional | 330 | 0.0000 | 0.0008 | **99.1%** | 0.9% | 0.060 |
| fixed_conditional | 330 | 0.0000 | 0.0135 | 83.6% | 22.1% | 0.567 |
| learned_conditional | 330 | 0.0000 | 0.0160 | 74.8% | 33.9% | 0.632 |
| terminal_only | 330 | 0.0000 | 0.1743 | 69.7% | 30.3% | 0.089 |

Oracle calls buy fewer refinement nodes, and on families this easy the nodes
were worth more than the tightening. Closed and timeout do not sum to 100%: an
instance can finish unresolved without timing out.

### Predictor arms, gene-disjoint development split

| Arm | Seeds | Mean Spearman | SD |
|---|---:|---:|---:|
| certifiable_gamma0 | 5 | **+0.7253** | 0.0109 |
| certifiable_gamma1 | 5 | +0.4573 | 0.1593 |
| unconstrained | 5 | +0.1212 | 0.0653 |

> **G3 passes and should not be reported as support.** Noninferiority against a
> baseline scoring 0.121 is not evidence that the constraint is cheap — it is
> evidence the baseline lost. The honest reading is that the unconstrained arm
> was badly fit on this split, and the gate as written cannot tell that apart
> from a real result.

**G4 fails cleanly** and is the most informative number in the table: turning
the residual on costs 0.268 Spearman, CI entirely below zero. It is not noise,
and it directly contradicts what the residual was added to provide.

### Execution record

Both GPUs were checked as free; the run took GPU 1 (NVIDIA TITAN RTX) capped at
15% memory. The verifier path is numpy and CPU-bound; the predictor arms used
the GPU. 1,650 verifier records, 15 predictor records, 291 timeouts, 30 voids
excluded, 0 errors, **0.239 of 8 authorised core-hours** in 860 s wall. Every
seed, instance, arm, timeout and hardware setting is in the log.

---

## 13. Model-effect intervals  *(class B — 0 failures)*

For a restricted family D ⊂ F, the effect of the restriction is bracketed —
never read off a difference of upper bounds.

```
effect in [ max(0, L_F - U_D),  U_F - L_D ]

U_F - U_D alone is NOT an effect, and is logged under that name.
```

| Check | Result |
|---|---:|
| Intervals validated against enumeration | 75 |
| Containment failures | **0** |
| Contrast violations | **0** |
| Verdicts — positive / negligible / **abstain** | 31 / 2 / **7** |

The seven abstentions are the point. Each has a lower bound of exactly 0.0 with
a non-zero upper bound, so the interval contains zero and the honest verdict is
"unresolved", not a direction. The naive quantity U_F − U_D is recorded
alongside under the field name `not_an_effect_UF_minus_UD` so it cannot be
mistaken for one.

---

## 14. Reproduction of the winning configuration

| Metric | Recorded | Reproduced |
|---|---:|---:|
| Cases | 500 | 500 |
| Baseline median gap | 0.5348 | 0.5348 |
| Conditional median gap | 0.0830 | 0.0830 |
| Tighter than baseline | 324 | 324 |
| Edges proved impossible | 1,684 | 1,684 |

**46 of 46** hashes verified before anything ran — 35 source files, 5 data
files, 6 pinned dependencies. A moved hash aborts the run rather than reporting
a mismatch afterwards. Rerun took 112 s and returned exit 0. Results are written
as immutable files with `.sha256` sidecars, and the table exporter re-hashes and
refuses any file whose content moved. Stated limitation, recorded in the
artefact: this verifies the current interpreter against pins; it does not
provision a fresh container.

> **A provenance defect worth knowing about.** 33 of 38 phase-one runs were made
> from a dirty tree, and two F7 runs share a git SHA *and* a config hash yet
> differ materially (50.00 against 29.25) because a code fix landed between
> them. The recorded provenance cannot distinguish them. The correct code state
> had to be reconstructed from git history per file.

---

## 15. The RNA application  *(blocked)*

Framed throughout as **on-target efficacy and ranking**. No off-target or safety
quantity is computed anywhere in the repository.

### The assay-identity check overturned the data

Huesken 2005 measured against **34 YFP reporter plasmids differing in their
3′UTRs, in H1299 cells at 48 h**, with labels normalised so the positive control
reads 90.0% inhibition and the least active siRNA 0% — which is why values
legitimately exceed 1.0. The redistributed annotation this project had relied on
says **HeLa** and **luciferase**. Both are wrong; only the timepoint agrees.

> **The consequence.** If the assayed construct is a YFP reporter carrying a
> 3′UTR fragment, the real folding context is the reporter's. Substituting a
> native GENCODE transcript is a silent construct substitution — and F9 had
> already done exactly that. The audit now refuses it rather than picking a
> plausible native accession, and the refusal is a code path
> (`assert_chemistry_not_folded`, `windows_for_isoforms`), not a comment.

### Splits: two work, two are impossible

| Split unit | Groups | Test rows | Residual cluster leak | Status |
|---|---:|---:|---:|---|
| gene | 30 | 630 | 1.11% | usable |
| sequence cluster | 770 | 608 | 0.00% | clean |
| study | 1 | — | — | INFEASIBLE |
| patent family | 0 | — | — | INFEASIBLE |

The historical 2182/249 split leaks badly: **87.1%** of test rows share a 13-mer
cluster with a training row, **100%** share a gene, 20.9% share a seed region.
Exact duplicates: 0. Largest k-mer cluster: 36 rows. A single-study corpus
cannot support a study-generalisation claim, and the patent-family unit is
simply absent from the data — both are recorded as INFEASIBLE rather than
approximated.

### Ensemble formulation

```
yhat = b + a^T mu + gamma * E[g(X,z)]
     = E_p[f(X,z)]   by linearity of the anchor
```

So the additive part needs only **first moments**; only the residual needs an
expectation. No latent structure is ever paired with the ensemble readout — a
test asserts the API accepts no per-structure label. The audit ran **102 checks
with 0 failures**: exact marginals match enumeration, every sample is a valid
family member, and the deterministic mean bracket contains the exact mean in
every case.

The **deterministic mean bracket** follows from taking expectations of the
affine envelopes; it needs only μ and carries no sampling error. It is currently
far looser than Monte Carlo — median width 4.99 against a median sampling error
of 0.31 — which is worth stating rather than burying. Exact marginals are used
only for the validated grammar; the ViennaRNA prior is a different grammar and
is labelled `unvalidated_grammar` with no exactness claimed.

### Sampling, unblocked

| Window | Candidate pairs | Lattice | 256 draws | Previously |
|---:|---:|---:|---:|---|
| 50 nt | 316 | 10^95 | 0.02 s | tractable |
| 100 nt | 1,483 | 10^446 | 0.09 s | prohibitive |
| 150 nt | 3,336 | 10^1004 | 0.20 s | prohibitive |

A stochastic-backtrace sampler over the inside chart replaced exact sequential
conditioning, which cost one partition evaluation per candidate per draw. Total
variation against exact enumeration is 0.028–0.055 at 4,000 draws on the three
tiny sequences, every draw valid. Both samplers are validated against
enumeration, not against each other; sequential is retained as an independent
cross-check. Sampled |error| on the tiny audit improved 0.1491 → 0.0294.

### Arm comparison: 2 of 8 admissible

| Arm | Status | Missing capability |
|---|---|---|
| `linear_sequence_only` | admissible | — |
| `sequence_target` | admissible | — |
| `mfe_gnn` | blocked | target_window, folding |
| `mean_adjacency_gnn` | blocked | target_window, marginals |
| `sampled_ensemble_gnn` | blocked | target_window, sampling |
| `anchor_A` | blocked | target_window |
| `fixed_A+B` | blocked | target_window |
| `learned_A+B` | blocked | target_window, policy |

Guide sequence alone predicts on-target efficacy at **Spearman 0.6464** under a
sequence-cluster-disjoint split, seed range 0.632–0.677 across 5 seeds; P@20 =
0.130 (seed range 0.05–0.20). The `sequence_target` arm scores identically and
*should*: the target site is the reverse complement of the guide, so the two
carry the same information by construction. Paired difference −2.3e−05 across 5
seeds. Reporting that as a null result is the correct outcome; reporting it as
"target context adds nothing" would not be, because no real target context was
available to either arm.

The certified-ranking rule `u_i < l_j` resolves **45.5%** of candidate pairs (60
of 132, 12 candidates) with the model ordering correct on **100%** of resolved
pairs — as a certificate requires, since it certifies the model's own ordering.
Median interval width 8.59. This is on declared windows and is flagged in the
result file as `NOT an RNA result`: agreement with held-out measurement is not
computable, because no measured ensemble readout exists for any window we can
legitimately construct.

### Dataset acquisition

Davis 2025 (`10.1093/nar/gkaf479`) is confirmed **CC-BY 4.0** via Crossref
(`content-version: vor`), PMCID PMC12205987. Its GEO accession GSE231101 is
3P-seq polyadenylation mapping, **not** knockdown labels — recorded explicitly
so it cannot be mistaken for one later. Supplementary Table S1 carries sequence,
chemical scaffold, target gene and, critically, **assay type per measurement**,
so native and reporter rows can be separated rather than pooled.

> **Acquisition stopped deliberately.** PMC serves a reCAPTCHA challenge for the
> supplementary archive. I did not work around the bot gate; the partial HTML
> that came back was deleted rather than parsed. The licence permits
> redistribution once obtained — a human download is required. Logged in
> `configs/rna_data.yaml` under `acquisition: {status: blocked}`.

---

## 16. Defects the validation caught

Found by independent checks rather than by inspection. Several would have
produced entirely plausible-looking numbers.

| Where | Defect | How it surfaced |
|---|---|---|
| `rna_noncrossing` | **Soundness.** Forced pairs sharing an endpoint were silently accepted — the partner map was a plain dict, so the second assignment overwrote the first and the recurrence counted structures satisfying only one of them | the new sampler drew 48 of 64 invalid structures; the exhaustive mask test had missed it because on the tested sequences the conflicting pair is non-canonical |
| `layered_dag` | Two forced edges in one layer treated as a disjunction; returned 2 paths where the truth is 0 | brute-force comparison |
| `conditional` | Conditional bound sometimes *looser* than baseline — different valid relaxations are not pointwise comparable | paired soundness test |
| `conditional` | Budget counter charged cache hits as work — 80 passes against a budget of 8 | budget test |
| `model` | Interaction-gap probe invalid unless supports are disjoint; it was measuring overlap, not interaction | gap came out non-zero at γ = 0 |
| `model` | Index buffers persistent, so `load_state_dict` overwrote the permuted map | relabelling test |
| toys | Degenerate anchor — every candidate pair on the toy sequence is G–C, collapsing all coefficients to one number | identical results across every seed |
| `run_pilot` | Predictor arms produced identical numbers: the closure ignored the model entirely, which would have made a gate pass vacuously | four arms agreeing to four decimals |
| `evaluate_gates` | NaN reported as FAIL rather than UNDETERMINED | 0/0 metric on an already-exact baseline |
| `walk_count` | Matrix power where L matrix-vector products suffice | ~100× speedup, outputs verified identical |
| `ensemble` | `RNA.cvar.rand_seed` does not exist; the assignment was swallowed by a try/except, leaving sampling unseeded | replaced with `RNA.init_rand` |
| corrections script | Would have overwritten `models.py` and dropped the affine readout, inflating every prediction by ~250,000× | caught by capturing forward-pass baselines before applying; merged instead of overwritten, forward verified bit-identical |

---

## 17. Limits of validity

| Dimension | What was covered | What is not claimed |
|---|---|---|
| Families | non-crossing pairings, layered-DAG paths, path matchings; ground sets to 45 edges | no family needing an oracle beyond those implemented; no unbounded families |
| Architecture | 2-layer message passing, sum aggregation, ReLU, signed weights, linear readout, non-negative additive anchor | no attention, normalisation or gating; no depth beyond 5; no aggregation outside the monotone class |
| Numerics | float64 with no directed rounding; interval arithmetic implemented and tested, not wired in | no bound in this work is numerically certified |
| Scale | exhaustively enumerable families for every soundness claim; 150 nt windows for descriptive statistics only | no soundness claim validated beyond enumerable size |
| Statistics | 5 seeds, cluster bootstrap by instance | predictor CIs rest on 5 seeds and one held-out gene — wide, and not a population claim |
| RNA labels | reporter readout, marginal over latent structures | no structure-conditional, causal, chemical-variant or safety claim |

---

## 18. Where it stands

### The defensible paper today

Not the one that was planned. The strongest asset is the **κ_L
characterisation**: the relaxation gap equals a ratio of length-L walk counts,
verified to a maximum relative error of 2.12e−16 by exhaustive enumeration, with
an exponent set by network **depth** — 0.912 per layer, R² 0.9996 — measured on
a synthetic family containing no application data. The gap is unbounded on
downward-closed families that are not union-closed and exactly 1.0 on join-
closed ones. And on real windows κ₂ bounds the measured slack at 50, 100 and 150
nt while κ₁² predicts κ₂ to within 2.5–4.6%.

That is a general account of *when* and *how badly* a relaxation of a
combinatorial family is loose. It is domain-free, has no collision in two
independent novelty audits, and it reframes the anchor and conditional-support
machinery as an instantiation whose practical value is honestly reported as not
yet demonstrated.

### Novelty position, re-audited

| Component | Status | Nearest prior work |
|---|---|---|
| Conditioning intermediate bounds on partial edge assignments | not novel | Hojny et al. 2024, aggressive bounds tightening |
| A learned component inside a verifier | not novel | Lu & Kumar, ICLR 2020 |
| Additive learned anchor over a DP-optimisable family | not novel | MXfold2, Nat Commun 2021 |
| GNN certification under structural perturbation | not novel | Hojny 2024; Ladner 2025; GNNev 2025; RobLight 2025 |
| Perturbation set defined by combinatorial membership, not a norm ball | undetermined | no collision found in two rounds |
| Budgeted oracle queries to tighten a certificate | undetermined | Top-k on a Budget (2026), different setting |

The two undetermined cells were the intended contribution, and G1/G2 mean
neither is currently supported empirically. Novelty position and empirical
position agree: defensible in principle, unproven in practice. The word "first"
is not used anywhere in the write-up without a stated scope.

### Next steps, in order

1. **Human download of Davis 2025 Supplementary Table S1.** CC-BY 4.0, blocked
   only by a CAPTCHA. It unblocks six arms, structural sensitivity and
   chemical-variant prediction, and resolves the construct question rather than
   working around it.
2. **Rerun G1 and G2 on harder instances.** The current families close to gap
   zero within a second, so the metric cannot see the arms. Longer windows and
   larger k are needed before conditioning can be assessed at all — this is a
   redesign of the instance distribution, not of the gate.
3. **Wire `numerics.py` into `bounds.py`.** Stages N1–N4 in
   `docs/bounds_proof.md`. This is what moves the whole stack off
   `float_unverified`; the machinery is written and tested.
4. **Decide the framing.** Lead with κ_L and demote the mechanism to an
   instantiation, or hold the paper until the harder-instance rerun and the
   Davis data are both in.

> **Nothing was submitted, uploaded, or sent.** No paper drafted for submission,
> no external contact made, no compute beyond the authorised cap, no private
> data accessed, no bot gate circumvented — and no seed selected, retuned, or
> dropped anywhere in this record.

---

*Every number in this report is read from an immutable result file under
`runs/`; none was computed at write time. Preregistered thresholds were fixed on
development data before the pilot ran. All bounds carry
`numerical_status: float_unverified` — see `docs/bounds_proof.md`. The `certmp`
package was not modified at any stage of phase two.*
