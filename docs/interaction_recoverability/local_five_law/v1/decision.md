# Decision: local regular

**PROVED HERE:** At the specified exponential reference, with both profiles still unknown in the original infinite-dimensional analytic model and all four coefficients unknown, the complete five-law experiment admits a centered finite-variance influence function for \(I_{12}\). A controlled weighted-\(C^3\) neighborhood also admits such a function. This resolves the requested local score equation; it does not resolve the earlier global modulus or prove an attainable root-\(N\) estimator.

The strongest model-specific result is a constructive repair of the infinite-variance anchor formula. For \(T_j(x)=h^{-1}(h(x)+j)\), the six-term anchor combination

\[
E(x)=q_0(T_2x)+q_0(T_1x)+q_0(T_2x-1)
      -q_4(x)-q_4(T_1x)-q_3(T_2x-1)
\]

equals \(r_1(T_2x-1)-r_1(x)\). The map \(C(x)=T_2x-1\) sends \([3,6]\) into itself and has derivative at most \(575/613<1\). It therefore determines all profile differences in a central interval through a convergent series.

Outward integer telescoping gives weights that decay at both tails. The remaining periodic term is confined to the central interval and evaluated by that contracting series. Quantile evaluation kernels are \(1/2\)-Hölder in the source \(L^2\) spaces, so the core series converges geometrically in influence norm. The tail envelope is a polynomial times \(\exp(-.15z^2+C|z|)\), which is integrable. Each source component has finite variance separately.

This profile representation has zero coordinates for \(e_1,e_2\); those two laws are included in the full experiment and are essential to the coefficient construction. Anchor reconstruction plus \(r_j(0)=0\) yields bounded representers of the profile values needed at \(z=\pm1\). Subtracting them from the \(e_i\) quantile perturbations produces a two-by-two system with

\[
\det M_i=2e^{\alpha_i+\beta_i}\sinh(\sqrt2/10)>0.
\]

The resulting four bounded coefficient duals annihilate the **entire** nuisance score closure. Consequently the efficient coefficient information is positive definite, and every residual target coefficient is compatible. Combining the profile influence with these duals solves the full score equation. This conclusion uses no finite-bank projection or numerical rank threshold.

The [proof](proof.md) gives the influence as explicit integrals and convergent series, a finite norm bound, the equivalent minimum-norm projection formula, and a \(1/2\)-Hölder bound for the constructed influence under a specified neighborhood of profile/coefficient perturbations in the common Gaussian coordinate. The efficient variance itself was not computed. Tail quantile evaluations can give large constants.

**Closest equivalence — STANDARD RESULT APPLIED:** This is a concrete solution of the adjoint-range regularity problem in [Trabs, Theorem 2.7](https://arxiv.org/html/1307.6610v2), followed by ordinary nuisance-score projection and finite coefficient correction. The model-specific work is the contracting warp/anchor identity, finite-tail representation, and explicit coefficient duals. General adjoint criteria, quantile derivatives, telescoping, and projection formulas are established tools. No broad publication-novelty claim follows from this focused comparison.

**Principal remaining obstruction — OPEN:** Prove an asymptotically linear expansion with \(o_P(N^{-1/2})\) remainder for a feasible empirical-quantile/calibration estimator when series depth and tail cutoff grow and coefficients are estimated. The current influence depends on the truth. This pass does not supply weighted empirical-quantile remainder control, nuisance rates, or a uniform nonlinear calibration argument. Pointwise regularity and a controlled representer neighborhood do not by themselves prove an endpoint Hellinger Lipschitz bound, global minimax rate, or attainable root-\(N\) estimation.

**Actual diagnostics — NUMERICALLY CHECKED:** One fixed reference, two analytic profile directions, all four coefficient directions, and one mixed direction were checked. No samples or fits were used.

| Diagnostic | Actual result |
|---|---:|
| Sampled maximum \(C'\) / analytic upper bound | 0.926387785 / 0.938009788 |
| Six-term scaled residuals | \(4.325\times10^{-15}\), \(2.727\times10^{-15}\) |
| Maximum profile derivative error, 512 terms | \(2.881\times10^{-9}\) |
| Maximum coefficient recovery error | \(1.754\times10^{-9}\) |
| Target identity residuals, quadrature order 128 | \(6.484\times10^{-14}\), \(4.797\times10^{-14}\) |
| Maximum 64-to-128 quadrature change | \(6.220\times10^{-8}\) |

The initial check failed an absolute residual threshold; its exact script and failure record are retained. The successful rerun used a scale-aware cancellation check. Neither run is an existence certificate or efficient-variance bound. See [raw checks](checks.json), [implementation](checks.py), and [commands/resources](commands_and_resources.md).

This is a defensible local statistical theorem in a strongly calibrated model. A substantive learned-model contribution and biological siRNA validation remain unestablished. No additional training, experiment, replacement architecture, or new authorization is requested. The prior global unresolved decision, negative results, pilot, predictions, configurations, and historical preparation ledger remain preserved.
