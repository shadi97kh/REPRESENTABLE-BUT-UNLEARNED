> Publication copy of `docs/novelty_matrix.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `2dc04e8d7f3afcbf5f1df9554404413101bf6025384cc28bad45e5ee6816c084`.

# Novelty matrix

Audit date **2026-09-05**. Compares the proposed system against seven primary references
named by the PI. Primary sources were fetched directly; where only metadata was reachable
this is recorded as a gap rather than filled in.

---

## 1. References, as read

| key | reference | access |
|---|---|---|
| **autoLiRPA** | Xu, Shi, Zhang, Wang, Chang, Huang, Kailkhura, Lin, Hsieh. *Automatic Perturbation Analysis for Scalable Certified Robustness and Beyond.* [arXiv:2002.12920](https://arxiv.org/abs/2002.12920) | **§3.2 read in full** (ar5iv) |
| **Hojny** | Hojny, Zhang, Campos, Misener. *Verifying message-passing neural networks via topology-based bounds tightening.* 2024. [arXiv:2402.13937](https://arxiv.org/abs/2402.13937) | **full text read** (ar5iv) |
| **Ladner** | Ladner, Eichelbeck, Althoff. *Formal Verification of Graph Convolutional Networks with Uncertain Node Features and Uncertain Graph Structure.* TMLR 2025. [arXiv:2404.15065](https://arxiv.org/abs/2404.15065) | abstract + metadata |
| **Lu/Kumar** | Lu, Kumar. *Neural Network Branching for Neural Network Verification.* ICLR 2020. [arXiv:1912.01329](https://arxiv.org/abs/1912.01329) | abstract + metadata |
| **GNNev** | Liu, Lu, Kwiatkowska. *Exact Verification of Graph Neural Networks with Incremental Constraint Solving.* 2025. [arXiv:2508.09320](https://arxiv.org/abs/2508.09320) | abstract + metadata |
| **RobLight** | Lu, Tan, Benedikt. *Robustness Verification of Graph Neural Networks Via Lightweight Satisfiability Testing.* 2025. [arXiv:2510.18591](https://arxiv.org/abs/2510.18591) | abstract + metadata |
| **GraphStar** | chapter [doi:10.1007/978-3-032-32357-6_13](https://doi.org/10.1007/978-3-032-32357-6_13) | **GAP — paywalled**, see §5 |

### Equations recovered verbatim

**autoLiRPA §3.2**, backward linear relaxation through a computational-graph node:

```
Σ_{j∈V} Λ_j h_j(X) + Δ  ≤  A_i h_i(X)  ≤  Σ_{j∈u(i)} Λ̄_j h_j(X) + Δ̄
```

Supported specifications: ℓ_p balls `S = {X : ‖X − X₀‖_p ≤ ε}`, synonym substitution sets
`S(w_i)`, and "flexible perturbation specifications beyond ℓ_p-balls". **Explicitly not
supported:** perturbation of a graph adjacency matrix or discrete graph structure — the
framework treats the computational graph as a DAG of *fixed topology* and perturbs node
values only.

**Hojny**, the message-passing layer certified:

```
x_v^(l) = σ( Σ_{u∈V} A_{u,v} · w_{u→v}^(l) · x_u^(l−1) + b_v^(l) )
```

Perturbation sets: `P₁(A*)` undirected, edge addition **and** removal, global budget `Q`
and local budgets `q_v`; `P₂(A*)` directed, removal only. Aggressive bounds tightening
partitions neighbours during branch-and-bound into

```
V₀ = { u ∈ V : A_{u,v} = 0 }      V₁ = { u ∈ V : A_{u,v} = 1 }
```

i.e. it recomputes intermediate bounds **conditioned on edges already fixed on or off**.

---

## 2. Comparison matrix

Rows are the dimensions along which the proposal could differ. "—" means not applicable.

| dimension | autoLiRPA | Hojny | Ladner | Lu/Kumar | GNNev | RobLight | GraphStar | **proposed** |
|---|---|---|---|---|---|---|---|---|
| model class | any computational graph | MPNN, ReLU | GCN, multi-layer | feedforward NN | GNN | GNN variants | GNN (GCN, GINE) | MPNN, signed weights |
| node-feature perturbation | **yes** | yes | yes | — | yes | yes | yes | yes |
| adjacency / structure perturbation | **no** (fixed topology) | **yes** | **yes** | — | **yes** | **yes** | **yes** (edge features) | **yes** |
| perturbation set shape | ℓ_p ball, substitution sets | edge-budget ball (global + local) | uncertainty set | — | budget constraints | small perturbation | star set over node+edge | **hard combinatorial family** (nesting, one partner/base) |
| membership in the set decidable in O(1)? | yes | yes | yes | — | yes | yes | yes | **no — needs an oracle** |
| bound mechanism | linear relaxation (LiRPA/CROWN) | MIP + branch-and-cut (SCIP) | polynomial zonotope reachability | B&B with learned branching | incremental constraint solving | polynomial-time partial solvers | star-set reachability | anchor (exact DP) + bounded residual |
| bound granularity | per computational-graph node | **per node, per feature, per layer** | per set representation | per neuron split | per constraint | per constraint | per star set | **per edge / per message** (proposed) |
| conditions intermediate bounds on partial edge assignments | no | **YES — `abt`, V₀/V₁** | no (dependencies, not conditioning) | — (splits neurons, not edges) | implicitly via CSP | — | unknown (gap) | yes |
| exact combinatorial oracle over the family | no | no | no | no | no | no | unknown (gap) | **yes** |
| query budget, fixed | no | no | no | no | no | no | unknown (gap) | **yes** |
| **learned** allocation of that budget | no | no | no | **learned branching** (different object) | no | no | unknown (gap) | **yes** |
| completeness | incomplete | complete (MIP) | incomplete | complete | exact | incomplete by design | incomplete | incomplete, one-sided |
| aggregations | — | sum (+ linear pooling) | GCN norm | — | **sum, max, mean** | various | GCN, GINE | sum, max, logsumexp, min |

---

## 3. The question asked: is edge-conditioned intermediate support genuinely different?

**Answer: no, not as a tightening mechanism. It is an existing method adapted to a new
oracle.**

The mechanism — recomputing intermediate bounds under a partial assignment of edges — is
**Hojny's aggressive bounds tightening**, published in 2024. Their `abt` partitions each
node's neighbourhood into `V₀` (edges fixed absent) and `V₁` (edges fixed present) and
re-derives bounds under that conditioning inside branch-and-bound. That is the same idea,
and it is the graph specialisation of ordinary bound tightening under branching, which is
standard across the verification literature the PI already excluded.

Three differences survive, and only the second and third are substantive:

1. **Granularity** — Hojny bounds per node × feature × layer; the proposal bounds per edge
   / per message. This is a refinement of the same technique, not a new one. It is worth
   at most a sentence in a paper, and it needs an ablation showing the finer granularity
   actually buys tightness rather than only cost.

2. **Source of the conditioning information.** In Hojny the conditioning comes from the
   *search*: the branch-and-bound tree has already fixed those edges. In the proposal it
   comes from an *oracle answering a question about a domain-defined family* — "does a
   valid nested structure containing this edge set exist, and what is the additive
   optimum over the sub-family it induces?" No reference in the table queries such an
   oracle. This is a real difference, but it is a difference of *what supplies the
   conditioning*, not of how the bound is then tightened.

3. **The perturbation set is not a ball.** Every one of the six accessible references
   defines the admissible set by a budget or a norm, so membership is trivially decidable.
   The proposal's set is `{ valid nested secondary structures of this sequence }`, whose
   membership involves a non-crossing constraint and a one-partner-per-base constraint.
   This is the strongest genuine gap in the matrix, and it is what makes the oracle
   necessary rather than decorative.

**Consequence for the paper.** The layer cannot be presented as a new tightening
mechanism. It can honestly be presented as *the first application of conditioned
intermediate bounds to a family whose membership is not decidable in closed form*, with
the oracle and the budget allocation as the contribution — provided §5's gap is closed.

---

## 4. What may and may not be claimed

| candidate claim | verdict |
|---|---|
| "Novel edge-conditioned bound tightening" | **NO** — Hojny `abt` |
| "First certification of GNNs under structural perturbation" | **NO** — Hojny, Ladner, GNNev, RobLight, GraphStar all do this |
| "First to support non-ℓ_p perturbation sets" | **NO** — autoLiRPA §3.2 supports substitution sets and states generality |
| "Learned component inside a verifier" | **NO** — Lu/Kumar |
| "Novel per-edge bound granularity" | **WEAK** — refinement of Hojny; needs an ablation to be worth stating |
| Certification over a perturbation set defined by a **hard combinatorial membership constraint**, queried via an exact oracle | **CANDIDATE** — no accessible reference does this |
| **Bound-preserving** learned allocation of a fixed oracle-query budget, valid under every allocation | **CANDIDATE** — no accessible reference has any query budget |
| κ_L as an exact characterisation of relaxation slack, with depth-proportionality and unboundedness | **STRONG, ALREADY IN HAND** — see §6 |

### Defensible scope sentence, if "first" is used at all

Only this form is currently supportable, and only after §5 is closed:

> To our knowledge, this is the first robustness certificate for a message-passing network
> whose admissible perturbation set is a combinatorial family with non-trivially decidable
> membership, rather than a norm ball or an edge-budget ball.

Every qualifier in that sentence is load-bearing. Dropping "to our knowledge", "message-
passing", or "non-trivially decidable membership" makes it false against the table above.

---

## 5. Gaps — recorded, not filled

| gap | effect on the matrix |
|---|---|
| **GraphStar full text is paywalled.** `doi.org` → `link.springer.com` → Springer IdP authorization gate. Not pursued: authenticated access was not authorized | Four cells unknown: whether GraphStar conditions on partial edge assignments, uses an exact oracle, has a query budget, or supports constraint-defined families. Public sources describe it as a Star-set generalisation capturing uncertainty over node **and edge** features, extending NNV to graph inputs, applied to power-system GNN surrogates. That is set-based reachability, so a query budget is *unlikely* — but unverified |
| **Ladner, Lu/Kumar, GNNev, RobLight read at abstract level only** | Their rows are populated from abstracts and metadata. The four "no" entries under *oracle* and *query budget* are inferences from what the abstracts describe, not from reading the methods. Must be confirmed before any "first" claim |
| **No citation-graph sweep** | Work citing Hojny 2024 or Lu/Kumar 2020 has not been enumerated. This is where a budgeted-oracle GNN verifier would most plausibly hide |

---

## 6. The result that is already strongest

The repository already contains a contribution that none of the seven references
approaches, and it is a positive result:

> For a certified depth-`L` sum-aggregation network that is affine on the non-negative
> orthant, the relaxation gap equals `κ_L(F) = W_L(top) / max_{G∈F} W_L(G)`, a ratio of
> walk counts — verified to a maximum relative error of **2.12e−16** by exhaustive
> enumeration over five families (`runs/20260904-043204`).

With three corollaries already measured: the exponent is proportional to network depth
(0.912 per layer, R² 0.9996, `runs/20260904-034559`); slack is unbounded on
downward-closed families that are not union-closed (9.0 → 2290.3); and it is exactly 1.0
on join-closed families. It predicts biological slack growth on real GENCODE windows
(×4.17, ×2.22 against measured ×4.39, ×2.36, `runs/20260904-043511`).

This is a theory of *why* relaxations of combinatorial families are loose, it is not
RNA-specific, and it makes the anchor a predicted intervention rather than a hope. None of
autoLiRPA, Hojny, Ladner, GNNev, RobLight or Lu/Kumar offers a closed-form characterisation
of relaxation gap. It is the strongest asset in the project and it is already logged.

---

## 7. Go / no-go

**Recommendation: CONDITIONAL GO, with the contribution repositioned.**

The proposal as framed — "novel edge-conditioned message-bound layer" — is **no-go**. The
tightening mechanism is Hojny's, the learned-component-in-a-verifier idea is Lu/Kumar's,
and the anchor is MXfold2's. Leading with the layer would put the paper's headline
directly on top of a 2024 paper that does the same conditioning.

**Go**, if the paper leads with the two things the matrix shows are open:

1. **κ_L as an exact theory of relaxation slack** (already in hand, §6), with
2. **certification over a combinatorially-constrained family via an exact oracle**, of
   which the budget allocation is the engineering contribution and the layer is an
   implementation detail — described as an adaptation of known bound tightening, cited to
   Hojny, not claimed.

Conditions on the go, in order:

| # | condition | why |
|---|---|---|
| 1 | Obtain GraphStar full text, or record it as an unclosable gap in the submission | It is the only reference whose oracle/budget cells are unknown |
| 2 | Read Hojny §on `abt`, GNNev and Ladner methods in full | Four "no" cells are currently abstract-level inferences |
| 3 | Drop "edge-conditioned" and "message-bound layer" as novelty language | Naming and mechanism both collide |
| 4 | Do not use "first" except in the §4 scope sentence | — |
| 5 | Fix the seven UNSUPPORTED items in `docs/claim_ledger.md` §D before any draft | `README.md` currently contradicts the code in four places |
| 6 | Resolve the "exact" wording (ledger B3) | Standing instruction against unsupported exactness claims |

**No training has been started, and none is recommended until conditions 1–3 are met.**
If condition 1 or 2 closes the family-constraint gap, fall back to the κ_L theory as the
paper, which stands on its own and is already supported by logged results.
