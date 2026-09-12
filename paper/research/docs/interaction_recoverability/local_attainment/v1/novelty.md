> Publication copy of `docs/interaction_recoverability/local_attainment/v1/novelty.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `b2c205ae4c17f11ebe93d6d2af6deee0ee38d1afec5a329a59b8b2dfefb88805`.

# Bounded prior-result comparison

**Decision context:** The completed claim is pointwise estimator attainment in one explicitly calibrated nonlinear five-law model. It is not a new general principle of functional estimation and does not establish superiority over an equally informed method.

| Closest result | Actual relation to this construction | What remains model-specific |
|---|---|---|
| [Trabs, Theorem 2.7 and nonlinear-model discussion](https://arxiv.org/html/1307.6610v2) | General adjoint-range criterion, efficient influence and information framework. The previous local theorem solved its score equation; that alone did not attain a nonlinear statistical upper bound. | A single source-data statistic and a controlled remainder proving its local asymptotic linearity. The statistic need not attain the minimum-norm influence. |
| [Kaji, Theorem 4.1, Proposition 4.2, Theorem 5.1](https://arxiv.org/pdf/1910.07572) | Inverse-map differentiability and integrated quantile/L-statistic limits. These are established tools; integration before profile reconstruction is not novel by itself. | Growing contraction queries, actual tail weights, moving calibration arguments, and an explicit simultaneous schedule require additional arguments. |
| [DExtrI, Theorems 3.3–3.4 and appendix](https://arxiv.org/html/2608.19849v1) | Identification and extrapolation in shared-latent pre-additive models, with axis-interval information and structural rigidity conditions. The appendix includes broader profile regimes than the main exponential-polynomial theorem. | Here only five discrete source laws are observed, but the warp, Gaussian latent law and private unit loadings are known. The profile class and information assumptions are different, not ordered by generality. This local root-\(N\) target theorem is not supplied by population identification alone. |
| [Bennett et al., strongly identified functionals, Theorem 2 and section 4.1](https://arxiv.org/html/2208.08291v3) | Functional estimation can remain regular despite weakly identified nuisance functions; operator/source restrictions and nuisance estimation conditions matter. Functional recovery without global nuisance recovery is established territory. | The present nonlinear density experiment is not assumed to satisfy that paper's conditional-moment operator equations. The explicit quantile/contraction/calibration proof verifies this model's conditions directly. |

## Exact Kaji hypothesis audit

**STANDARD RESULT APPLIED / limitation checked.** The source responses have positive continuous densities on \((0,\infty)\). For \(m(y)=y\), the integrability required by Kaji Proposition 4.2 holds:
\[
\int_0^\infty\sqrt{P_a(y)(1-P_a(y))}\,dy
=\int_{\mathbb R}\sqrt{\Phi(z)(1-\Phi(z))}F'_a(z)dz<\infty,
\]
using the quadratic exponent \(.1-.25=-.15\). Thus its unweighted integrated inverse-quantile conclusion is applicable. Theorem 5.1 fixes a Lipschitz bound on the integrating function; our \(U/\phi\) need not have a fixed bounded slope/integrator bound as truncation grows. Its theorem is not invoked to cover the triangular array, accumulating point evaluations, or estimated coefficients. Parts A–D of our theorem supply those arguments explicitly.

## Exact advance and boundaries

**PROVED HERE:** This continuation upgrades the saved local finite-variance influence to a specified data-defined statistic with an attained pointwise root-\(N\) limit and regularity along fixed admissible \(N^{-1/2}\) paths. The nonautomatic steps are:

1. summing compact empirical-quantile remainders over \(K_N=O(\log N)\) contracting queries;
2. controlling actual signed tail weights through the empirical-CDF area identity, including the jump at 4 and cutoff jumps;
3. proving a source-only minimum-residual calibration rule has the required joint four-coefficient expansion, including moving profile arguments;
4. checking one explicit choice of all truncation and arithmetic schedules.

The general tools are elementary empirical-process concentration, quantile linearization, functional inversion, minimum-distance estimation, Taylor expansion, CLT and LAN. These are not renamed as new ML machinery.

The earlier [direct-target sufficient global consistency result](../../direct_target/v1/estimator.md) applies over a broader original parameter set with slower conservative bounds. The new result is local to a controlled neighborhood and has no uniform-in-warp guarantee. It neither replaces that global result nor resolves the [global modulus decision](../../direct_target/v1/decision.md).

**OPEN:** Whether this particular attainment argument has an exact antecedent beyond the bounded comparison above. No exhaustive publication-novelty search was conducted. No equally informed direct-target inverse, likelihood, confidence-set, or functional-estimation control was shown inferior. Such controls would receive the same five laws, latent-law/warp knowledge, private-anchor calibration, coefficient region and target.

The currently sufficient sample-domain bounds and \(N^{6+o(1)}\) direct implementation are impractical. A substantive ML contribution would require additional scientific value, a practical justified method, and evidence for its claimed distinction. Biological siRNA validity would separately require justified interventions, latent/warp assumptions, replicates, and measured target evaluation. Neither follows from this local theorem or its analytic fixtures.
