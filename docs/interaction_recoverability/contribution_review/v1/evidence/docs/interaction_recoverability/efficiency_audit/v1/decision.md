# Decision: unresolved efficiency; backend no-go retained

**PROVED HERE:** the known normalization \(g_1(0)=g_2(0)=1\) fixes the reference-source median at 2. Its centered indicator \(Q\) is orthogonal to the **entire** five-law score closure, including all four coefficient directions. Thus, for the saved complete influence,
\[
\psi_c=\psi-20\langle\psi,Q\rangle_\pi Q
\]
is fully compatible, and
\[
L_4\le L_8\le L_{12}\le V_{\rm eff}
\le U=V_{\rm constructed}-20\langle\psi,Q\rangle_\pi^2
\le V_{\rm constructed}.
\]
All twelve prescribed directions have admissible two-sided paths. The three exact submodel Grams are positive definite; their information bounds use an ordinary inverse, with no ridge or finite-sieve replacement of the original model.

**NUMERICALLY CHECKED — approximations to the population quantities:**

| Quantity | Value |
|---|---:|
| Saved constructed variance | 79,390,122.7931 |
| Complete median covariance \(c\) | 56.8032660332 |
| Reference correction weight \(20c\) | 1,136.0653206631 |
| Removed variance | 64,532.2206 (0.081285%) |
| Corrected upper-bound approximation \(U\) | 79,325,590.5724 |
| Four-coefficient submodel lower-bound approximation | 1,568.5179383 |
| Eight-direction submodel lower-bound approximation | 1,981.9138168 |
| Twelve-direction submodel lower-bound approximation | 2,045.7389999 |
| Upper/lower ratio | 38,776.0074 |

The covariance is not an anchor-only calculation: its anchor contribution is \(-.0862098813\), while all four coefficient corrections contribute \(56.8894759144\). The diagnostic indicates that the saved influence does not already remove the median component. The exact theorem guarantees nonincrease; strict nonzero covariance has not received a complete numerical certificate.

**Numerical qualification:** the mathematical inequality holds for the defined exact quantities. The table is not a certified numerical interval for \(V_{\rm eff}\). The saved analytic influence truncation error is below 5.905 in source Hilbert norm; it also bounds the corrected norm error, and bounds covariance error by \(5.905/\sqrt{20}\). Retained-domain integration and floating-point errors remain uncertified. New Gram and target tails are bounded separately. Refinement differences, numerical ranks and tiny algebraic residuals do not replace full certification. Both frozen configurations and every result are retained in [raw.json](raw.json).

**Exact remaining gap:** with \(\mathcal T\) the original full score closure and \(r_{12}\) the twelve-score projection,
\[
U-L_{12}
=\|\psi_c-P_{\mathcal T}\psi_c\|_\pi^2
+\|P_{\mathcal T}\psi_c-r_{12}\|_\pi^2.
\]
Neither nonnegative term is sharply bounded. We have not determined how much of the large constructed variance is avoidable, nor how close the finite-submodel bound is to efficient variance. No endpoint is justified as an approximation to \(V_{\rm eff}\).

**STANDARD RESULT APPLIED:** known-moment control variates, orthogonal nuisance projection and finite-submodel information bounds are established machinery. The closest reduction is the minimum-norm score/Riesz framework, already compared with [Trabs](https://arxiv.org/html/1307.6610v2) and the saved Riesz/one-step pilot. This is a useful model-specific efficiency audit; it establishes no substantive new ML contribution, efficient estimator, finite-sample MSE theorem or biological siRNA result.

The median adjustment itself uses the observed source-zero statistic \(-(\lambda/N)\sum_i Q_0(Y_{0i})\). A public fixed reference weight is implementable but need not improve variance at other truths. A consistently estimated optimal weight requires an additional source-data consistency proof, which is absent. An oracle influence calculation is not that proof.

The corrected scale remains essentially unchanged, and no practical finite-sample regime has been demonstrated. The [backend no-go](../../practical_feasibility/v1/feasibility_decision.md) therefore remains. The failed adaptive pilot's 3.5569% improvement against its 10% screen is preserved. No training, samples, fits, simulations, expanded dictionary, anchor optimization or backend work followed the audit.

This was a primary-agent review plus two separate agents' static internal cross-checks, not an external independent review. See the [complete argument](efficiency_argument.md), [frozen calculation](calculate.py), and [preservation/resource accounting](commands_and_resources.md). The single deterministic job completed successfully; no further numerical cycle is proposed.

