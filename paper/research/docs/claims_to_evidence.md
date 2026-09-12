> Publication copy of `docs/claims_to_evidence.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `6d9625913ef00809c21eee13d959db1f6797710ca8b3cf537b59cb6af6c723f1`.

# Claims to evidence

Written **before drafting**, so that the paper can only contain claims that already have
a row here. A claim with no row does not go in. A claim whose row says `not supported`
does not go in as a positive result.

`kind` separates what the evidence *is*: **lemma** (proved, checked by exhaustive
enumeration, no fitted parameters), **engineering** (a soundness or reproduction check,
no fitted parameters), **empirical** (fitted parameters involved).

Numerical status of every bound in this project is `float_unverified` — a diagnostic, not
a certificate. See `docs/bounds_proof.md`.

---

## Supported

| # | Claim | kind | Evidence | Where |
|---|---|---|---|---|
| S1 | Endpoint exactness holds iff the aggregator is monotone under multiset inclusion, in either direction; antitone aggregators are certifiable with the endpoint roles swapped | lemma | 7 aggregators × 25 seeded trials, exhaustive over 2^8 subsets; direction of the maximum checked, not only exactness | `c3_aggregator_characterization` |
| S2 | Under non-negative parameters and features, sum + ReLU is *exactly affine*, so a monotone-certified network of that form is linear | lemma | additivity violation 2.157e−16; 12/12 configurations agreed with prediction | `c4_affinity` |
| S3 | The relaxation gap equals a walk-count ratio κ_L = W_L(top)/max_F W_L | lemma | max relative error 2.121e−16 against exhaustive enumeration over five families | `c7_kappa_characterization` |
| S4 | The gap's scaling exponent is set by network **depth**, not by the data | lemma | 0.912 per layer at depths 1–4, R²=0.9996, on a synthetic family containing no application data | `c5_slack_scaling` |
| S5 | The gap is unbounded on downward-closed families that are not union-closed, and exactly 1.0 on join-closed families | lemma | 9.0 → 2290.3 as m: 3 → 48; join-closed verified by brute force | `c5_slack_scaling` |
| S6 | The affine bound propagation is sound at every intermediate, not only at the output | engineering | 1,262,952 bracket checks against exhaustively enumerated families, 0 violations | `audit_bounds` |
| S7 | Conditional message support is sound and never looser than the inherited fallback | engineering | 500 randomised cases, 0 violations; min-of-valid-bounds retained at both levels | `audit_conditional` |
| S8 | Refinement certificates are replayable and independently checkable | engineering | 72/72 accepted, 2,222 obligations, by a checker that does not import the search; 4 tamper tests reject | `check_certificates` |
| S9 | The model-effect interval [max(0, L_F−U_D), U_F−L_D] contains the true effect | engineering | 75 intervals validated by enumeration, 0 containment failures | `evaluate_explanations` |
| S10 | The winning development configuration reproduces exactly under pinned dependencies | engineering | 46/46 source and data hashes verified; 5/5 metrics reproduced bit-for-bit | `reproduce` |
| S11 | At matched **query budget**, conditional support reduces the median relative gap from 0.535 to 0.083 | empirical | 500 cases, reproduced independently | `audit_conditional`, `reproduce` |

## Not supported

| # | Claim | kind | What the evidence says | Where |
|---|---|---|---|---|
| N1 | Conditional support is tighter at matched **wall time** | empirical | **Undetermined, trending against.** The unconditional baseline closes 98–100% of pilot instances to gap exactly 0 within the matched second, so the relative-reduction metric is 0/0. Where the arms differ, conditioning is *behind*: 10–32% timeouts vs 0–2%, because oracle calls buy fewer refinement nodes | gate G1 |
| N2 | Learned allocation beats deterministic allocation | empirical | **Undetermined here, at the noise floor in development** (random 0.0830 vs largest-gap 0.0830 vs measured-gain 0.0806, the last costing 2,867 teacher passes) | gate G2, `audit_conditional` |
| N3 | The nonlinear residual adds held-out predictive value | empirical | **Refuted on this label.** −0.268 Spearman, 95% CI [−0.408, −0.128]. Adding the residual *hurts* | gate G4 |
| N4 | The certifiable model class is noninferior to an unconstrained one | empirical | Gate G3 passes numerically (upper CI −0.187 < 0.05 margin) but **the baseline is weak**: the unconstrained MPNN scores 0.121 against the certifiable model's 0.457. Noninferiority against a baseline that lost badly is not evidence for the constraint being cheap | gate G3 |
| N5 | The certificate bounds a biological quantity (efficacy, off-target risk) | — | **Not supportable in principle with this data.** The label is one number per sequence, marginal over latent structures; `max_z f(X,z)` has no measured counterpart. Retired by name in the preregistration | `preregistration.md` §6 |
| N6 | Any bound in this work is numerically certified | engineering | **No.** All propagation is float64 with no directed rounding. Validated interval arithmetic exists and is tested but is not wired into the bound path | `docs/bounds_proof.md`, `scmp/numerics.py` |

## Novelty position (re-audited 2026-09-05)

| Component | Status | Nearest prior work |
|---|---|---|
| Conditioning intermediate bounds on partial edge assignments | **not novel** | Hojny et al. 2024, aggressive bounds tightening (V₀/V₁ partition inside branch and bound) |
| Learned component inside a verifier | **not novel** | Lu & Kumar, ICLR 2020 |
| Additive learned anchor over an exactly-DP-optimisable family | **not novel** | MXfold2, Nature Communications 2021 |
| GNN certification under structural perturbation | **not novel** | Hojny 2024; Ladner 2025; GNNev 2025; RobLight 2025; GraphStar |
| Perturbation set defined by a combinatorial membership constraint rather than a ball | **undetermined** | no collision found in two independent search rounds; absence of a hit is weak evidence |
| Budgeted oracle queries used to tighten a certificate | **undetermined** | nearest is *Top-k on a Budget* (2026), a different setting |

The two undetermined cells were the intended contribution. **N1 and N2 mean neither is
currently supported by evidence**, so the novelty position and the empirical position now
point in the same direction: the mechanism is defensible in principle and unproven in
practice.

## What may be written today

Only S1–S11, and only with N1–N6 stated alongside. The strongest defensible paper on
present evidence is the **κ_L characterisation** (S3–S5): a closed-form account of *when*
and *how badly* a relaxation of a combinatorial family is loose, with the anchor and
conditional-support machinery as an instantiation whose practical value is reported as
not yet demonstrated.
