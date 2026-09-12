> Publication copy of `docs/architecture_audit.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `b5a4c46d2fe7e80683545a3bdde6a83b5fb8dd2c18490d502e96969e270b35b0`.

# Architecture audit

Audit date **2026-09-05**. Repository HEAD `35c0716`, tracked tree clean.
Scope: what exists in code, what it computes, and what validates it. No code changed.

Source of truth is the code as hashed in `docs/environment_manifest.md`, not the prose in
`README.md` or `RESULTS.md`.

---

## 1. Module map

12 modules in `certmp/` (≈ 830 lines), 19 experiment scripts, 1 test file with 7 tests.

### `certmp/models.py` — the model class under certification

| symbol | line | role |
|---|---|---|
| `MonoMPNN` | L46 | numpy reference model. **This is the model every certificate is computed on** |
| `MonoMPNN.forward(X, A)` | L62 | `H ← act(agg(A, H·softplus(W) + B))` per layer, readout `scale · (Σ_v H_v · softplus(out)) + shift` |
| `MonoMPNN.structurally_certifiable()` | L74 | gate: `nonneg ∧ agg ∈ CERTIFIABLE ∧ act ∈ MONOTONE_ACTS ∧ scale ≥ 0` |
| `MonoMPNN.endpoint_direction()` | L80 | dispatches isotone/antitone |
| `MonoMPNN.is_affine()` | L84 | `nonneg ∧ agg ∈ LINEAR_AGGS ∧ act ∈ INERT_ACTS` |
| `TorchMonoMPNN` | L99 | trainable twin; `export_numpy()` (L155) transfers weights to `MonoMPNN` |

Constants: `MONOTONE_ACTS` = relu, tanh, sigmoid, softplus, identity · `INERT_ACTS` =
relu, identity · `LINEAR_AGGS` = sum, mean, degnorm · `TORCH_SUPPORTED_AGGS` = sum, max.

**Design note that matters for the next stage.** The certified path is numpy-only and
single-graph; `TorchMonoMPNN` exists solely to train and hand weights back. Parity is
enforced by one test, not by construction.

### `certmp/aggregators.py` — the certified class

Seven aggregators registered with a predicted direction: `sum`, `max`, `logsumexp`
(isotone); `min` (antitone); `mean`, `degnorm`, `std` (neither). `CERTIFIABLE` is derived
from the direction field, so the class is data-driven rather than hard-coded.

### `certmp/reach.py` — the certificate kernel

| symbol | line | role |
|---|---|---|
| `build_A(n, mandatory, optional, keep)` | L6 | adjacency with self-loops |
| `_endpoints(...)` | L14 | evaluates the two lattice endpoints, ordered by direction |
| `exact_interval(...)` | L23 | **the theorem**: returns `(lo, hi, 2)`; raises outside the certified class or on `X < 0` |
| `brute_force_interval(...)` | L44 | enumerates all 2^k subsets; guarded at `max_k=20` |
| `endpoint_gap(...)` | L55 | relative gap of the true optimum over the endpoint |

### `certmp/ensemble.py` — uncertainty-set construction

| symbol | line | semantics |
|---|---|---|
| `bpp_lattice(seq, lo=0.05, hi=0.90, floor=1e-3)` | L17 | **banded, NOT sound**. Mandatory = p > hi, optional = band, rest discarded |
| `sound_lattice(seq, floor=0.0, canonical_only=True)` | L42 | mandatory = ∅, nothing discarded; restricted to `CANONICAL` pairs with `MIN_LOOP = 3` |
| `generic_features(seq, decay=3.0)` | L86 | 12-dim, length-independent, non-negative |
| `onehot_features(seq, extra)` | L109 | 4-dim + optional extras |

### `certmp/certify.py` — the one-sided certificate

`certify_threshold(...)` L19 returns `certified`, `worst_case`, `tau`, `forward_passes`,
`k`, `lattice_size_log10`, `margin`, `sound_over_ensemble`. `certify_topk` L41 **raises
`NotImplementedError`** — withdrawn 2026-09-04.

### `certmp/kappa.py` — the combinatorial theory

`walk_count(A, L)` (L26, `L` matrix–vector products), `kappa_L(...)` (L39),
`measured_gap(...)` (L54), `enumerate_family(...)` (L60, guarded at 2^20), and five family
predicates including `fam_matching`, `fam_forest`, `fam_join_closed`.

### `certmp/maximal.py` — the rejected alternative

`saturate(pairs, seq, n)` L61 greedily fills loop regions to a maximal valid structure;
`is_valid` L74 and `is_maximal` L91 are the checkers.

### Support modules

`train.py` (`fit`, `predict`, `save_numpy_model`, `load_numpy_model`), `data.py`
(`verify`, `load_dsir_split`, `spearman`, `pearson` — implemented locally, no scipy),
`void.py` (`check_live`, `check_nonempty_lattice`), `provenance.py` (`new_run`, `record`).

---

## 2. What validates what

| validator | location | what it actually checks | independent? |
|---|---|---|---|
| `brute_force_interval` | `reach.py:44` | endpoint vs full 2^k enumeration | **No** — shares `model.forward` with the thing under test. It validates *which adjacency is extremal*, not the forward pass |
| `test_exact_matches_brute_force` | `tests:6` | 2-pass equals brute force, sum and max | as above |
| `test_rejects_uncertified_class` | `tests:18` | mean/degnorm/signed rejected | yes (negative control) |
| `test_rejects_negative_features` | `tests:26` | H4 enforced | yes |
| `test_void_on_empty_lattice` | `tests:33` | k = 0 raises `VoidRun` | yes |
| `test_torch_numpy_parity` | `tests:40` | torch vs numpy forward, 2 aggs × 3 depths, incl. non-trivial affine readout | **Yes** — two independent implementations |
| `test_affine_readout_preserves_extremality` | `tests:74` | scale/shift do not break endpoints | no (uses `endpoint_gap`) |
| `test_negative_scale_rejected` | `tests:92` | `scale < 0` leaves the safe class | yes |
| `is_valid` / `is_maximal` | `maximal.py` | saturation output is a valid maximal structure | **Yes** — checker written separately from `saturate` |
| F5 standalone reimplementation | not in repo | reproduced F1 digit-for-digit | **Yes**, but **the file was never committed** — see gap G-3 |

---

## 3. Gaps found

| id | gap | consequence |
|---|---|---|
| **G-1** | `exact_interval` is float64 with a `1e-12` relative tolerance. No interval arithmetic, no rational arithmetic, no directed rounding | The word "exact" in the codebase means *exact in the combinatorial argument*, not *numerically certified*. A certificate schema field `numerical_status` would currently have to read `float_unverified` |
| **G-2** | The soundness of `sound_lattice` + `certify_threshold` over the Boltzmann ensemble is **argued and measured**, never proved in code or text | F7 measured 0 violations in 20 000 sampled structures. That is evidence, not a proof |
| **G-3** | The independent F5 reimplementation lives only in session scratch | The strongest independence claim in the project is not reproducible |
| **G-4** | `brute_force_interval` guarded at `max_k = 20`; all validation is at k ≤ 11 | No validation at the k = 408–3932 used in F9 |
| **G-5** | `bpp_lattice` remains the default constructor and is still used by F2, F3, F3b, F6 | Correct and documented, but a caller who reaches for the obvious name gets the unsound object |
| **G-6** | `TorchMonoMPNN` supports only sum and max, while `CERTIFIABLE` now includes logsumexp and min | Nothing certified with logsumexp/min can be trained |
| **G-7** | No support-oracle, no per-edge bound, no budget mechanism exists | The entire proposed stage-2+ architecture is absent. Confirmed by function inventory |

---

## 4. Fitness for the proposed direction

The proposed system needs four things. Against the inventory:

| proposed component | nearest existing code | status |
|---|---|---|
| exact learned additive structured anchor | none. No DP anywhere in the repo | **absent** |
| signed-weight nonlinear MPNN | `MonoMPNN(nonneg=False)` exists and is *rejected* by `structurally_certifiable()` | present as a negative control only |
| conditionable graph-family support oracle | `enumerate_family` (exhaustive, 2^20 cap) and `maximal.saturate` are the closest | **absent** at usable scale |
| edge-conditioned message-bound layer + budget | none | **absent** |

The repository is a *certificate-by-monotonicity* codebase. The proposed direction is a
*bound-propagation-with-oracle* codebase. Reusable across the boundary: `ensemble.py`
(family construction), `kappa.py` (the slack theory that motivates the anchor), `data.py`,
`provenance.py`, and the RNA plumbing in `_huesken.py` / `f9_target_context.py`.
`reach.py` and `certify.py` do **not** carry over — they assume monotonicity, which the
new design deliberately drops.
