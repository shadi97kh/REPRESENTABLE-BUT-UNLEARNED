# Feasibility decision: no-go for full backend development now

**NUMERICALLY CHECKED / scientific judgment:** this constructed estimator shows a poor reference influence scale. The evidence does not justify implementing its complete certified backend or requesting sampled validation now. The positive compact-attainment theorem is preserved. This is a decision about the present construction, not a lower bound on the five-law statistical problem.

## Actual statistical scale

The complete influence calculation, including all five sources, all four coefficient directions and shared-source covariance, gives approximately
\[
V=7.9390\,10^7,\qquad \sqrt V=8910.11.
\]
The implied asymptotic standard errors at total \(N=10^3,10^4,10^5,10^6\) are approximately \(281.76,89.10,28.18,8.91\). The reference contrast is about 3.895. Source 0 contributes roughly 85.7% of the complete variance; the coefficient correction dominates the anchor component.

**OPEN / scope:** these are deterministic reference approximations, not finite-sample performance measurements, MSE or coverage guarantees. The projected estimator's small coefficient rectangles also make extrapolating this CLT to the displayed counts particularly questionable. No sample size or biological accuracy threshold is prescribed. A large variance for this nonminimum-norm influence does not establish a large efficient variance.

Depth, tail and integration checks preserve the same scale. Integration refinement changed the standard-deviation constant by about .01; float64 versus extended accumulation changed it by about \(4\,10^{-11}\). Neither comparison is a certified numerical error bound. The latter also shares float64 CDF values and initial quadrature nodes.

**PROVED HERE:** the exact influence truncated at depth 512, response radius 12 and lattice radius 24 differs from the full influence in \(L^2\) by less than 5.905. **OPEN:** the remaining numerical integration and complete floating-point error have no tight certificate, so that 5.905 bound cannot simply be placed around the numerical answer as a certified confidence/error interval. The [variance analysis](variance_analysis.md) gives the exact triangle-inequality interpretation and the raw covariance matrix.

## Precision improves, but does not cure the scale

**PROVED HERE:** the two-stage population-contraction argument extends to uniformly \(o_P(N^{-1/2})\) update errors, including the explicit choice \(O(N^{-1})\), without empirical contraction. All four coefficient estimates retain the compact expansion. The finite target's coefficient sensitivity is proved \(O_P(1)\); input, direct-target, integration, cap and clipping effects are separately controlled. Thus the **complete estimator** remains asymptotically equivalent in probability under an explicit polynomial accuracy allocation.

The sufficient allocation
\[
\tau_N=[16N(1+10350(K_N+1)+J_Y)]^{-1}
\]
reduces input and internal precision orders to \(O(\log N)\) under the same certified-primitive model. It does not promise deterministic agreement on every capped dataset.

**NUMERICALLY CHECKED — public formulas only:**

| Total \(N\) | Fractional input bits | Cap integer-magnitude bits |
|---:|---:|---:|
| 1,000 | 40 | 20 |
| 10,000 | 43 | 27 |
| 100,000 | 47 | 34 |
| 1,000,000 | 51 | 40 |

These are fixed-point absolute-accuracy components, before signs/guards; they are not complete floating-point mantissa prescriptions. Internal known-function and knot calculations need their own certified tolerances with potentially large constants. Fixed-resolution observations retain \(C_{\rm data}\tau\) update error and \(J_Y\tau\) direct-target error. Extra working precision cannot supply absent observation digits.

## Integration and implementation implications

**PROVED HERE:** rigorous panel bounds may use observed polygon slopes and bounds on known derivatives over the entire panel. Derivatives of the contraction iterates decay geometrically, so higher derivative terms need not pay equally for every level. The zeroth-order term, knot-preimage gaps and weight jumps still require explicit bounds.

**OPEN:** no observed polygons or complete interval backend were created here, so practical panel savings are unknown. The previous millions-of-panels counts are sufficient overcounts, not lower bounds. The approximately nine-second influence diagnostic is not a runtime measurement for the estimator; its final process peak RSS was about 1.23 GiB.

The conditional mathematical \(N\,\mathrm{polylog}N\) dependencies survive the static audit. An exact arbitrary-real cap classifier was not supplied by the old finite-input specification; the new in-probability proof explicitly handles a rounded gate. Component diagnostics still do not verify the complete backend.

## What has been established and what has not

| Claim | Status |
|---|---|
| Compact pointwise attainment | **PROVED HERE previously; preserved** |
| Complete statistical numerical stability with polynomial tolerances | **PROVED HERE in this pass** |
| Complete reference influence scale | **NUMERICALLY CHECKED**, with separate analytic truncation bounds and unresolved numerical error |
| Certified executable full backend | **OPEN** |
| Affordable estimator runtime at useful accuracy | **OPEN** |
| Finite-sample accuracy, MSE or coverage | **OPEN** |
| Efficient variance or intrinsic difficulty of the problem | **OPEN** |
| Publication novelty / substantive ML advantage | **OPEN** |
| Biological siRNA validity | **OPEN** |

Equally informed inverse-functional, likelihood or estimating-equation controls may use the same five laws, compact anchors, public preconditioner and numerical improvements. They need not recover unnecessary nuisance functions, and no disadvantage for those controls is demonstrated here.

No subsequent sampled-validation specification is proposed because the feasibility evidence is unfavorable. No automatic backend cycle, replacement project, training request or experiment follows this decision. A substantive reason to revisit this particular implementation would require new evidence addressing its influence scale and finite-sample applicability, not unused allowance.

The review is a static self-review with separate agent cross-checks inside the same system; it is not external independent review. See [numerical stability](numerical_stability.md), [arithmetic/local integration audit](arithmetic_addendum.md), [actual variance and error analysis](variance_analysis.md), and [executed work and resource accounting](commands_and_resources.md).

