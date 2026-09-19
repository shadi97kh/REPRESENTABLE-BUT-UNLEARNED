# Consequences and primary-source comparison — v2

**PROVED HERE, using the preserved family upper bound.** The verified construction establishes
\[
\log\left(1+\sup_{\eta\in\Theta}
V_{\rm eff}(\eta;\delta,\theta)\right)=\Theta(\delta^{-1})
\]
uniformly in \(\theta\in[1/2,2]\) as \(\delta\downarrow0\). Here \(\Theta\) is exactly the original profile and coefficient class at fixed known warp, and \(V_{\rm eff}\) is the full admissible-score information bound with weights \(1/5\). The [complete proof](lower_bound.md) retains the explicit constants of v1. This continuation verifies that conclusion; it does not claim a new result beyond v1 or a matching leading exponential constant.

## Distinct statistical conclusions

| Quantity | Verified conclusion |
|---|---|
| Worst-case efficient variance over the original truth class | Exponential logarithmic order in \(1/\delta\) |
| Efficient variance at the constructed truth | Explicit exponential lower bound, finite upper bound from retained attainment |
| Efficient variance at every truth | No exponential lower bound asserted |
| Historical reference \(g_1=g_2=e^x,\ \delta=.1,\theta=1\) | Its separate efficiency question remains open |
| Attained pointwise CLT variance | The saved influence remains compatible and finite; efficiency is unproved |
| Complete-experiment minimax MSE | Three-term finite lower bound; matching upper bound remains open |
| Sample requirements for fixed MSE | No exponential requirement follows from this proof; no polynomial sufficiency is proved |

The hard profiles depend on the known warp when taking the supremum. They remain fixed along every testing path. The result therefore rules out polynomial uniform bounds for compatible regular local variance over the original class. It does not prove such divergence for a single profile pair held fixed as \(\delta\) changes.

**PROVED HERE.** With \(\epsilon_\delta=e^{-1/(1280\delta)}\), the whole path obeys
\[
H(P_{1,t},P_{1,0})\le20\epsilon_\delta|t|+10\delta t^2.
\]
The second term limits the sufficient separation used by the testing proof. Squaring its scale yields \(1/(\delta\sqrt N)\), between the constant lower-bound regime and the eventual \(\epsilon_\delta^{-2}/N\) term. Equating the latter two terms gives the illustrative scale
\[
N=\delta^2\epsilon_\delta^{-4}=\delta^2e^{1/(320\delta)}
\]
up to fixed constants. This is a crossover in a conservative lower bound. Neither sharp Hellinger order nor a true minimax transition or CLT onset is established.

For a fixed MSE tolerance, the quadratic-path term can become small at counts polynomial in \(1/\delta\). At that point this lower bound alone stops excluding the tolerance. Existence of a successful estimator would still require an upper analysis. This distinction is essential when interpreting the exponential local variance.

## Established information theory

**STANDARD RESULT APPLIED.** Trabs's Theorem 2.7 characterizes regularity by the adjoint-score range and gives the minimum-norm influence and convolution bound in locally regular indirect experiments. The scalar inequality used here is its familiar submodel consequence:
\[
|DI[q]|^2\le \|\psi\|_\pi^2\|Sq\|_\pi^2.
\]
Our proof verifies a real two-sided coefficient path, its DQM property and nonzero target derivative, with exactly the original five-source weighting. Every full-model compatible influence must satisfy this one path equation. A finite tangent dictionary is unnecessary. The reference framework supplies no automatic attaining estimator for this nonlinear experiment; the separate saved family theorem supplies the compatible influence used for the upper bound. [Trabs, arXiv:1307.6610v2, Theorem 2.7 and Example 2.9](https://arxiv.org/html/1307.6610v2).

The score is strictly positive in norm at each allowed \(\delta>0\). “Exponentially weak information” is appropriate; an exactly singular experiment at positive \(\delta\) is not asserted.

## Adjacent weak-information and minimax comparison

**STANDARD RESULT APPLIED AS A COMPARISON.** Heinrich and Kahn study finite mixtures with known kernels and unknown mixing distributions, using expected Wasserstein-\(1\) loss. Under their respective assumptions, Theorems 3.2–3.3 give local minimax order \(n^{-1/[4(m-m_0)+2]}\), while Theorem 3.5 gives a nonuniform pointwise \(n^{-1/2}\) bound. Section 3.3 explains how truth-dependent constants and delayed asymptotic behavior reconcile the rates. [Heinrich–Kahn, arXiv:1507.04313v1, Sections 2.4 and 3.1–3.3](https://arxiv.org/html/1507.04313v1).

This is an adjacent comparison, not an application of their mixture theorems. Our observations are monotone pushforwards in five independent groups, with two unknown analytic profiles; the target is a scalar interaction and the loss is ordinary squared error. No mixture identifiability assumption, Wasserstein loss or mixture rate is transferred. The comparison follows the cited theorem statements.

## Contribution and limits

**PROVED HERE as a model-specific construction.** Saturation keeps the first profile globally admissible while matching the warped exponential derivative exponentially well over a Gaussian region of radius proportional to \(\delta^{-1/2}\). Exact nonlinear amplitude balance keeps four observed laws fixed. The fifth law has an exponentially small score at the central truth, but the interaction has a uniformly nonvanishing derivative. The full nonlinear likelihood segment supplies a finite-experiment result in addition to the local information bound.

Smooth saturation, score differentiation, information inequalities, Gaussian tail bounds, product affinity, two-point testing and higher-order weak-information mechanisms are established tools. The possible contribution is their explicit original-class realization in this five-law experiment, together with the retained attainable family theorem. No general ML principle, comparative superiority over equally informed methods or publication priority is claimed.

**OPEN.** Matching finite-sample MSE upper bounds, sharp exponential constants, exact efficient influences and efficient estimation, the historical reference-efficiency problem, practical accuracy and scientific publication significance remain unresolved. This bounded comparison is not an exhaustive priority search.

The [saved broader prior comparison](../../warp_family/v1/novelty.md) remains unchanged. The family-attainment theorem, external-review v2, failed learned-target pilot and its predictions/configurations, all ledgers, backend no-go and lack of biological validation are preserved. The present result authorizes no additional estimator, training, backend, dictionary or packaging cycle.

## Inspection provenance

For this continuation, the actual saved model, commutator/calibration, attainment, conditioning and novelty documents were read, together with the historical efficiency argument's reference, scores, norm and unresolved projection question. The cached Trabs v2 Theorem 2.7 was read; both cited primary HTML papers were retrieved and their relevant statements inspected. No new exhaustive literature search or external review was performed. Current execution details are in [commands_and_resources.md](commands_and_resources.md).
