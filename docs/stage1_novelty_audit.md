# Stage 1 — novelty audit and gate plan

Date **2026-09-05**. Status: **audit only, no implementation**.

Target: a positive-results ICLR 2027 submission in which RNA therapeutics is the
application and the contribution is methodological.

Proposed system, as stated by the PI:

1. an **exact learned additive structured anchor**,
2. a **signed-weight nonlinear adjacency-dependent MPNN**,
3. **certified conditionable graph-family support oracles**,
4. candidate novelty: an **edge-conditioned message-bound layer with learned allocation
   of a fixed support-query budget**.

Novelty is explicitly *not* claimed for GNN+DP, learned branching, or valid relaxation
slopes alone.

---

## 0. Why this is the current stage

No stage had been defined. Implementing the candidate layer before establishing what is
already published would risk building the one component that a reviewer rejects on prior
art. The audit is therefore stage 1, and it gates everything downstream. **No model code
was written in this stage.**

### Continuity with the existing repository

This is not a fresh start. The work at `35c0716` establishes a result that directly
motivates the proposed design, and that connection should be load-bearing in the paper:

- The monotone-MPNN certificate is **sound but loose** — 29.25× at 60 nt, 222.8× at
  150 nt (`runs/20260904-020253`, `runs/20260904-020229`).
- C7 proved the gap equals a walk-count ratio κ_L, exactly in the affine regime
  (max relative error 2.12e−16, `runs/20260904-043204`).
- C5 showed the looseness is **structural**: matchings are downward-closed but not
  union-closed, so the smallest containing lattice has an infeasible top element, and
  slack is unbounded in problem size (`runs/20260904-034559`).
- C8 showed κ tracks the measured biological slack (growth ×4.17, ×2.22 against measured
  ×4.39, ×2.36, `runs/20260904-043511`).

The proposed anchor attacks exactly the quantity C5 identified: it optimises over the
**feasible** family (where nesting and one-partner-per-base hold) rather than over a
lattice that contains configurations no molecule can adopt. That gives the new work a
pre-registered, quantitative prediction rather than a hope — see Gate G4.

---

## 1. Audit method, and its limits

Searched 2026-09-05 across arXiv, OpenReview, ACM DL, and a full-text academic corpus,
in five areas: GNN structural robustness certification; learned branching for neural
network verification; additive thermodynamic anchors with learned residuals; edge-
conditioned layers and bound propagation on graphs; budgeted oracle querying.

**Limits, stated so they are not mistaken for coverage.** This is a targeted search, not
a systematic review. The academic full-text corpus returned poor recall on ML conference
work and was discarded in favour of web search. Absence of a hit is **weak evidence** of
novelty and is recorded as *undetermined*, never as *novel*. Before any novelty claim
reaches a submission, the four papers marked **must-read in full** below need reading
beyond their abstracts, and a citation-graph sweep forward from Zügner & Günnemann 2020
and Lu & Kumar 2020 should be run.

---

## 2. Nearest primary work, by component

### Component 1 — exact learned additive structured anchor

| paper | what it does | bearing |
|---|---|---|
| Sato, Akiyama & Sakakibara, *RNA secondary structure prediction using deep learning with thermodynamic integration* (MXfold2), Nature Communications 12:941, 2021 · [paper](https://www.nature.com/articles/s41467-021-21194-4) | DNN folding scores **added to** Turner nearest-neighbour free energies, optimised by Zuker-style DP, with thermodynamic regularisation | **Directly prior.** A learned additive score over an exactly-DP-optimisable RNA family is published and well known |

**Verdict: NOT NOVEL.** The anchor construction itself must be cited as prior art and
framed as a *component we adopt*, not a contribution. The distinction available to us is
one of **purpose**: MXfold2 uses the DP optimum as a *prediction*; we would use it as a
*certified anchor* bounding a downstream nonlinear network. That distinction is real but
it is not by itself a paper.

### Component 2 — signed-weight nonlinear adjacency-dependent MPNN

This drops the monotonicity constraint that the existing repository relies on. Signed
weights and general MPNNs are of course standard. **Verdict: NOT NOVEL, and not claimed
as such.** Its role is to be the thing that *needs* bounding once monotonicity is gone.
F1/R5 in this repository already measured that signed weights destroy endpoint exactness
(11/25 exact, `runs/20260904-034428`), which is the motivation, not a result.

### Component 3 — certified conditionable graph-family support oracles

| paper | what it does | bearing |
|---|---|---|
| Zügner & Günnemann, *Certifiable Robustness of GCNs under Structure Perturbations*, KDD 2020 · [pdf](https://dl.acm.org/doi/pdf/10.1145/3394486.3403217) | Certifies GCNs under **structure** perturbation; formulates as a jointly constrained bilinear program and solves with a **novel branch-and-bound** for lower bounds; local and global budgets | **Closest prior art for the whole framing.** Must-read in full |
| Bojchevski & Günnemann, *Certifiable Robustness to Graph Perturbations*, NeurIPS 2019 · [arXiv](https://arxiv.org/pdf/1910.14356) | Provable robustness under a **general set of admissible graph perturbations**, explicitly flexible in threat model | Anticipates "conditionable family". **But** the model class is restricted to predictions linear in personalized PageRank — not a general nonlinear MPNN. Must-read in full |
| Wang, Jia, Cao & Gong, *Certified Robustness of GNNs against Adversarial Structural Perturbation*, KDD 2021 · [arXiv](https://arxiv.org/pdf/2008.10715) | Randomized smoothing for structural perturbation | Probabilistic, not deterministic — a contrast, not a conflict |
| Xia et al., *AGNNCert*, 2025 · [arXiv](https://arxiv.org/pdf/2502.00765) | "First certified defense for GNNs against arbitrary (edge, node, node-feature) perturbations with **deterministic** guarantees" | Different task (node/graph **classification** radius). Abstract does not mention an exact combinatorial oracle, DP anchor, or budget allocation. Must-read in full |
| *Certifying Robustness of GCNs for Node Perturbation with Polyhedra Abstract Interpretation*, 2024 · [arXiv](https://arxiv.org/pdf/2405.08645) | Polyhedral abstract interpretation for GCNs | Node-feature perturbation; relevant to bound tightness, not to family conditioning |

**Verdict: PARTIALLY ANTICIPATED.** A general admissible-perturbation-set certificate
exists (Bojchevski 2019) but only for a linear-in-PageRank model class. Structure
certification for GCNs with branch-and-bound exists (Zügner 2020). What we have not found
is a certificate whose perturbation set is a **domain-defined feasible family with a
non-trivial membership constraint** (nesting + one partner per base), queried through an
oracle. Recorded as **undetermined**, pending the two must-reads.

### Component 4 — the candidate novelty

| paper | what it does | bearing |
|---|---|---|
| Lu & Kumar, *Neural Network Branching for Neural Network Verification*, ICLR 2020 · [arXiv](https://arxiv.org/abs/1912.01329) | GNN learns to imitate strong branching, treating the network under verification as the graph | **Learned branching is taken.** PI already excluded this |
| *Neural Network Branch-and-Bound for Neural Network Verification*, 2021 · [arXiv](https://arxiv.org/abs/2107.12855) | Two GNNs: one imitates strong branching, one computes a **feasible dual solution giving valid lower bounds** | Closest to "learned component that still yields a valid bound". Must-read in full |
| *Top-k on a Budget: Adaptive Ranking with Weak and Strong Oracles*, 2026 · [arXiv](https://arxiv.org/html/2601.20989) | Concentrates weak-oracle effort on **items that block certification**, shrinking ambiguity and tightening intervals under a budget | **Nearest prior art for budgeted allocation targeted at certification.** Different setting (ranking), but the governing idea — spend a fixed oracle budget where it unblocks a certificate — is present |
| Simonovsky & Komodakis, *Dynamic Edge-Conditioned Filters in CNNs on Graphs*, CVPR 2017 · [arXiv](https://arxiv.org/abs/1704.02901) | ECC: filter weights conditioned on edge labels | **Naming collision only.** "Edge-conditioned" is an established term for a different mechanism |

**Verdict: NARROWED, NOT ELIMINATED.** After removing what is taken, the residual
candidate claim is:

> A *sound* per-edge message bound over a domain-defined feasible graph family, whose
> tightness is improved by a **fixed budget of membership/conditioning queries to an
> exact support oracle**, where the allocation of that budget across edges is **learned**,
> and where the bound remains valid for **every** allocation — so the learned component
> affects tightness only, never soundness.

The last clause is the part most likely to survive review, and it is the part that must
be proved and independently validated, not asserted.

### Naming actions required

- Do **not** call the layer "edge-conditioned convolution" — collides with ECC (2017).
- Do **not** describe the allocator as "learned branching" — collides with Lu & Kumar.
- Do **not** call the anchor "thermodynamic integration" — that is MXfold2's term.

---

## 3. Honest novelty position

| component | status | may we claim novelty |
|---|---|---|
| additive DP anchor over an RNA family | published (MXfold2) | **no** |
| signed-weight nonlinear MPNN | standard | **no** |
| structural certification of GNNs via B&B | published (Zügner 2020) | **no** |
| general admissible perturbation sets | published (Bojchevski 2019), restricted model class | **no** |
| learned branching | published (Lu & Kumar 2020) | **no**, per PI |
| certificate over a *nesting-constrained* feasible family for a *nonlinear* MPNN | **undetermined** | not yet |
| bound-preserving learned allocation of a fixed exact-oracle budget | **undetermined**, nearest is Top-k-on-a-Budget in a different setting | not yet |

Two "undetermined" cells is a viable position for a methods paper **only if** the
must-reads do not close them. If both close, the redesign options in §6 apply.

---

## 4. Certificate schema (required before any oracle is built)

Every certificate emitted by any stage must carry these fields. This is specification,
not implementation.

| field | meaning | failure mode it prevents |
|---|---|---|
| `assumptions` | explicit list, e.g. canonical pairs only, min hairpin loop 3, float64, model class | silent scope creep of a "sound" claim |
| `numerical_status` | one of `exact_rational`, `float_with_interval`, `float_unverified`; plus the guard used | floating-point error passed off as exactness |
| `family_hash` | hash of the family definition (sequence, constraint set, oracle version) | certificate quoted against the wrong family |
| `model_hash` | hash of weights and architecture | certificate quoted against a retrained model |
| `bound` | the numeric bound **and its direction** (upper/lower), one-sided by default | the two-sided-bound defect found in `certify.py` |
| `witness` | for a violated bound, the concrete structure attaining it; for a satisfied bound, the argmax structure or oracle transcript | unfalsifiable certificates |
| `termination_status` | `proved`, `budget_exhausted`, `timeout`, `numerically_aborted` | budget exhaustion silently reported as a proof |

**Independent validator requirement.** A certificate is accepted only when re-checked by
a validator that shares no code path with the producer — the pattern F5 already uses in
this repository, where a standalone reimplementation reproduced F1 digit for digit. The
validator checks the witness and the arithmetic, not the reasoning that produced them.

---

## 5. Stage plan and gates

Each gate is falsifiable and is to be pre-registered before the stage runs, in the style
of `PREREGISTRATION.md`. A failed gate preserves evidence and triggers a redesign
proposal — never a silent switch to a negative paper, and never a weakened baseline.

| stage | deliverable | gate | falsification criterion |
|---|---|---|---|
| **1** (this) | novelty audit, gate plan, certificate schema | **G1** — at least one component remains undetermined after the four must-reads | if all close, stop and redesign (§6) |
| 2 | formal statement of the anchor + bounded-residual decomposition; proof obligations enumerated | **G2** — soundness proof, plus an independent brute-force validator agreeing on small instances | any instance where the certified bound is below a realisable value ⇒ decomposition unsound, STOP |
| 3 | exact support oracle with conditioning, emitting the §4 schema | **G3** — oracle agrees with exhaustive enumeration on every instance small enough to enumerate; termination status always populated | any disagreement, or any `proved` returned on an exhausted budget ⇒ STOP |
| 4 | message-bound layer (no learning yet), uniform budget allocation | **G4** — certified bound strictly tighter than the κ-predicted lattice bound on the same 12 GENCODE windows, at matched compute | no improvement over the existing 21.5× / 94.3× / 222.8× ⇒ the anchor does not buy what C5 predicts; redesign |
| 5 | learned allocation of the fixed budget | **G5** — beats uniform, random, and a strong greedy allocator at *equal* budget, over ≥5 seeds, with the bound valid under every allocation | no gain over the strongest non-learned allocator ⇒ the novelty claim fails; report as such |
| 6 | RNA application, honest baselines | **G6** — pre-registered comparison against MXfold2-style anchor-only and against the monotone certificate | — |

Gate G5 is the one that decides whether the paper has its claimed contribution. It must be
run with the non-learned allocators tuned at least as carefully as the learned one.

### Standing constraints

No test-data tuning. No unsupported "exact" claims — `numerical_status` must justify the
word. No cherry-picked seeds; all seeds reported. No destructive commands. No
publication, private-data access, or paid compute without authorization. No promise of
numerical gains, and none is made here.

---

## 6. Redesign options if G1 or G5 fails

Prepared now so that a failure does not become an improvised pivot.

1. **If the family-conditioning claim closes** (someone has certified a nonlinear MPNN
   over a constraint-defined family): fall back to the *quantitative* contribution — κ_L
   as a predictive theory of relaxation slack, with the anchor as the intervention it
   predicts. C5 and C7 are already in hand and are positive results.
2. **If learned allocation does not beat greedy**: report the negative honestly and
   reframe around the *bound-preserving* property — an allocator that cannot damage
   soundness regardless of what it learns is a safety property worth stating, provided
   we do not dress a null result as a gain.
3. **If the anchor does not tighten the bound at G4**: the decomposition is wrong for
   this family, and κ_L predicts why. That is a redesign at the relaxation level, not a
   tuning problem — C5 showed tuning cannot fix a structural mismatch.

---

## 7. Stage 1 report

| item | value |
|---|---|
| **Files changed** | 1 added: `docs/stage1_novelty_audit.md`. No code. No tracked file modified |
| **Tests actually run** | `make test` → `sanity: all passed` (re-run post-audit to confirm the tree is untouched) |
| **Failures** | none in this stage. One tool failure: the academic full-text corpus returned low-relevance results for ML verification queries and was abandoned in favour of web search |
| **Runtime** | audit ≈ 6 min wall clock (6 searches, 1 fetch); no compute jobs run |
| **Evidence preserved** | `docs/artifact_reproduction.md` and all 42 run records untouched |
| **Blockers** | Four papers need reading in full before any novelty claim: Zügner & Günnemann KDD 2020; Bojchevski & Günnemann NeurIPS 2019; arXiv:2107.12855; AGNNCert arXiv:2502.00765. Two are paywalled at ACM DL and may need library access |
| **Next gate** | **G1** — read the four must-reads and confirm at least one component is still undetermined. Requires PI authorization to proceed, and confirmation of whether the four PDFs are reachable |

Nothing in this document asserts that the candidate novelty *is* novel. It asserts that
three of the four components are definitively not, names the papers, and narrows the
survivor to a claim precise enough to be checked.
