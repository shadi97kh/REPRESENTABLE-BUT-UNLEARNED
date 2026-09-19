# Affine bounds: real-arithmetic proof and numerical-soundness plan

Covers `scmp/bounds.py`. Part 1 proves the bound correct **in exact real arithmetic**.
Part 2 states what that proof does *not* give us in float64, and what would have to be
built to close the gap. Until Part 2 is done, every number the module produces is
labelled `float_unverified` and is a diagnostic, not a certificate.

---

## 0. Setup

Ground set of candidate edges `E`, `m = |E|`. Indicator `z ∈ {0,1}^m`, family
`F ⊆ {0,1}^m`. Adjacency `A(z) = B + P z`, where `B ≥ 0` entrywise is fixed and `P z`
places `z_e` symmetrically at both endpoints of `e`.

```
Z0 = X W0 + b0                      constant in z
M1 = A(z) Z0        H1 = relu(M1)
Z1 = H1 W1 + b1
M2 = A(z) Z1(z)     H2 = relu(M2)
g  = w_r^T Σ_i H2[i] + b_r
f(X, z) = b + a^T z + γ g
```

An **affine bracket** for a scalar quantity `q(z)` is a pair of affine functions with
`lo(z) ≤ q(z) ≤ hi(z)` for all `z ∈ [0,1]^m`. Note the relaxation to the box: `F` is a
subset of the box, so a bracket valid on the box is valid on `F`.

---

## 1. Real-arithmetic proof

### Lemma 1 (signed linear map)

Let `l_p ≤ x_p ≤ u_p` and `y = Σ_p W_p x_p + β`. With `W⁺ = max(W,0)`, `W⁻ = min(W,0)`:

```
Σ_p (W⁺_p l_p + W⁻_p u_p) + β  ≤  y  ≤  Σ_p (W⁺_p u_p + W⁻_p l_p) + β
```

*Proof.* Termwise. If `W_p ≥ 0` then `W_p l_p ≤ W_p x_p ≤ W_p u_p`. If `W_p < 0` the
inequalities reverse, giving `W_p u_p ≤ W_p x_p ≤ W_p l_p`. Summing preserves both. ∎

If each `l_p, u_p` is affine in `z`, both bounds are affine in `z`, since the coefficients
`W⁺, W⁻` are constants. This is `_split` and `_linear`.

### Lemma 2 (layer 1 is exact)

`Z0` does not depend on `z`. Hence

```
M1[i,k] = (B Z0)[i,k] + Σ_{e ∋ i} z_e · Z0[other(e,i), k]
```

is an exact affine function of `z`, and `lo = hi = M1`. No relaxation is used at layer 1.
This is why the bilinearity appears only at layer 2. ∎

### Lemma 3 (ReLU)

Let `lo ≤ m(z) ≤ hi` be an affine bracket and `l = min_box lo`, `u = max_box hi`.

* **Stable positive**, `l ≥ 0`: on `[l,u] ⊆ [0,∞)`, `relu(t) = t`, so `lo ≤ relu(m) ≤ hi`.
  Exact, no relaxation error.
* **Stable negative**, `u ≤ 0`: `relu(t) = 0` on `[l,u]`, so the bracket is `[0,0]`. Exact.
* **Unstable**, `l < 0 < u`: put `s = u/(u−l) ∈ (0,1)`. Then

  ```
  α · t  ≤  relu(t)  ≤  s (t − l)      for all t ∈ [l,u], any α ∈ [0,1]
  ```

  *Proof.* `relu` is convex, so on `[l,u]` it lies below the chord joining `(l, 0)` and
  `(u, u)`, which is `t ↦ u(t−l)/(u−l) = s(t−l)`. For the lower bound, `relu(t) ≥ 0` and
  `relu(t) ≥ t`, so `relu(t) ≥ max(0,t) ≥ α t` for any `α ∈ [0,1]` by convex combination. ∎

  Substituting the bracket is valid because both multipliers are non-negative:
  `s ≥ 0` gives `relu(m) ≤ s(m−l) ≤ s(hi−l)`, and `α ≥ 0` gives `relu(m) ≥ α m ≥ α·lo`.

The implementation uses `α = 1` when `u ≥ −l` and `α = 0` otherwise; any choice in `[0,1]`
is sound, so this affects tightness only.

### Lemma 4 (product envelope — all four inequalities)

For `e ∈ [0,1]`, `h ∈ [h_l, h_u]`, `y = e·h`, each of the following is the expansion of a
product of two non-negative factors and is therefore valid:

| product | inequality |
|---|---|
| `e (h − h_l) ≥ 0` | `y ≥ e h_l` |
| `(1−e)(h_u − h) ≥ 0` | `y ≥ h + e h_u − h_u` |
| `(1−e)(h − h_l) ≥ 0` | `y ≤ h + e h_l − h_l` |
| `e (h_u − h) ≥ 0` | `y ≤ e h_u` |

*Proof of row 3, the others identically.* `(1−e)(h−h_l) ≥ 0` expands to
`h − h_l − e h + e h_l ≥ 0`, i.e. `e h ≤ h + e h_l − h_l`. ∎

These four are the McCormick envelope and are exactly the convex hull of
`{(e,h,eh) : e ∈ [0,1], h ∈ [h_l,h_u]}`, so no valid affine inequality in `(e,h)` is
tighter.

**Closing the `h` terms.** Rows 2 and 3 contain `h` with coefficient `+1`. In a lower
bound a positive coefficient permits substituting `h ≥ lo_h(z)`; in an upper bound it
permits `h ≤ hi_h(z)`. Both substitutions preserve the direction, so the results are
affine in `z`. This is `mccormick`.

### Lemma 5 (selection among valid bounds)

If `y ≥ L₁` and `y ≥ L₂` pointwise then `y ≥ Lᵢ` for either `i`, and `y ≥ max(L₁,L₂)`.
Selecting one candidate is sound; `_tightest_lower` selects by box minimum, which affects
tightness only, never validity. Symmetrically for upper bounds. ∎

### Lemma 6 (sharing edge variables is sound and no looser)

Layer 2 reuses the same variables `z` as layer 1 rather than a fresh copy `z'`. Writing
`Ĝ(z, z')` for the two-copy relaxation, the shared version is `Ĝ(z, z)`. Since
`{(z,z) : z ∈ F} ⊆ F × [0,1]^m`,

```
max_{z ∈ F} Ĝ(z,z)  ≤  max_{(z,z') ∈ F × [0,1]^m} Ĝ(z,z')
```

so duplication is sound but can only be looser: it permits an edge to be present at one
depth and absent at another, which no real graph does. ∎

### Theorem (the upper bound)

Composing Lemmas 1–5 through the network yields affine `c, d` with

```
g(X, B + P z)  ≤  c^T z + d        for all z ∈ [0,1]^m.
```

For `γ ≥ 0`, and any `z ∈ F ⊆ {0,1}^m ⊆ [0,1]^m`,

```
f(X,z) = b + a^T z + γ g  ≤  b + γ d + (a + γ c)^T z.
```

Maximising the right-hand side over `F` and using that the support oracle returns the
exact maximum of a linear objective over `F`:

```
max_{z ∈ F} f(X,z)  ≤  b + γ d + support_F(a + γ c)  =:  U.        ∎
```

For `γ < 0` the residual enters with a negative multiplier, so the **lower** affine bound
on `g` is the one required; `certified_upper_bound` selects it by the sign of `γ`.

### Proposition (combined dominates separate maxima)

Let `U_sep = b + support_F(a) + γ · max_{[0,1]^m}(c^T z + d)`. Then `U ≤ U_sep`.

*Proof.* Using `max(A+B) ≤ max A + max B` and `F ⊆ [0,1]^m`:

```
U   = b + γd + max_{z∈F} [a^T z + γ c^T z]
    ≤ b + γd + max_{z∈F} a^T z + max_{z∈F} γ c^T z
    ≤ b + γd + support_F(a) + max_{z∈[0,1]^m} γ c^T z
    = U_sep.                                                        ∎
```

The inequality is strict whenever the two maximisers differ, which is why the combined
form is the one used. It is **not** strict when every coefficient is negative and the
family optimum is the empty structure — both sides then coincide, and the audit reports
how many cases are strict rather than assuming they all are.

### Corollary (γ = 0 is exact)

At `γ = 0` the residual is bypassed exactly, `f` is affine, and the oracle maximises it
exactly over `F`. Hence `U = max_{z∈F} f` with no relaxation gap. Asserted in the tests
and in the audit.

### Incumbents

The oracle returns a witness `z* ∈ F`. The incumbent reported is `f(X, z*)` evaluated
through the model itself — a value the model actually attains at a feasible point. No
envelope value is ever reported as achievable, so `incumbent ≤ max_{z∈F} f ≤ U` and the
reported gap is an honest optimality gap.

---

## 2. Numerical-soundness plan

Everything above is a theorem about real numbers. The implementation evaluates it in
IEEE-754 binary64 with round-to-nearest, which the proof does not cover.

### 2.1 Where roundoff enters

| # | site | why it can break the inequality |
|---|---|---|
| R1 | affine accumulation in `_linear`, `_message_*` | a rounded-down upper coefficient can fall below the true one |
| R2 | `box_min` / `box_max` | summing `m+1` terms accumulates error in the concretised `l, u` |
| R3 | ReLU slope `s = u/(u−l)` | division rounds; catastrophic cancellation when `u ≈ l` |
| R4 | McCormick constants `h_l, h_u` | inherited from R2, then multiplied |
| R5 | `support_F(a + γc)` | the oracle is exact for integers, but these coefficients are floats, so its max-plus arithmetic is float |
| R6 | model evaluation for the incumbent | torch float64, a different summation order from the numpy reference |

R3 is the sharpest: as `u − l → 0` the slope is ill-conditioned, and an unstable unit that
is *nearly* stable is exactly where the relaxation is most delicate.

### 2.2 What the tests currently assert

`SLACK = 1e-9`, an **assumed** tolerance, not a derived one. The exhaustive audit passes
1,262,952 intermediate bracket checks at that slack with zero violations, which is
evidence that roundoff is far below it in this regime — and is not a proof that it always
is. Nothing in Part 1 is invalidated; what is missing is a bound on the gap between the
real-arithmetic quantity and the computed one.

### 2.3 Plan, in stages

**N1 — running error bounds (cheapest, no dependencies).** Apply the standard result that
for a float sum of `n` terms, `|fl(Σx) − Σx| ≤ γ_n Σ|x|` with `γ_n = nu/(1−nu)` and
`u = 2⁻⁵³`. Carry an error accumulator alongside each affine form and widen every bound
outward by it. Turns `float_unverified` into `float_with_interval`. Deliverable: an
`ErrorTracked` variant of `Aff`; the audit reruns and must still show zero violations
against the *widened* bounds.

**N2 — directed rounding.** Compute upper bounds with rounding toward `+∞` and lower
bounds toward `−∞`. numpy does not expose the rounding mode, so this needs either a small
C extension, `mpmath` with explicit rounding, or interval arithmetic over `Fraction`.
Removes the need for N1's conservative constants at the cost of speed.

**N3 — exact rational on the additive path.** At `γ = 0` the entire computation is
`b + a^T z`. With `Fraction` coefficients and the oracle's existing exact-integer
max-plus, this becomes `exact_rational` with no slack at all. Cheap, and it makes the
`γ = 0` corollary a genuine certificate rather than a float claim.

**N4 — guard the ill-conditioned ReLU.** When `u − l` is below a threshold relative to
`max(|l|,|u|)`, refuse the chord and fall back to the interval bound `[0, max(0,u)]`,
which is sound for any `l, u` and immune to the division. Cost is looseness on nearly
stable units only.

### 2.4 Acceptance criteria

`numerical_status` may be promoted from `float_unverified` only when all hold:

1. Every bound carries an explicit outward widening derived from N1 or N2, not assumed.
2. The exhaustive audit passes against widened bounds with slack `0.0`, not `1e-9`.
3. `γ = 0` reports `exact_rational` via N3 and matches the oracle's integer path exactly.
4. A deliberately adversarial case — an unstable unit with `u − l < 1e-12` — is exercised
   and either bounded correctly or routed through the N4 fallback.
5. The status string is derived from which stages ran, never hard-coded.

Until then the module reports `float_unverified`, `docs/` says diagnostic, and no claim
in a paper draft may describe these bounds as certified.
