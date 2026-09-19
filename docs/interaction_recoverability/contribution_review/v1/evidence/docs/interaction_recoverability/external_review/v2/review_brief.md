# Interaction recoverability: brief for external mathematical review

## Status and purpose

The record is ready for external assessment as a **model-specific local regularity and pointwise attainment argument**, with its retained arithmetic qualifications. It is not a validated successful ML system. The constructed influence has unfavorable reference variance; efficient variance remains unresolved; full-backend development remains a no-go. Publication priority, substantive ML advantage, practical finite-sample accuracy and biological validity are unestablished.

This package assembles existing proofs without rerunning tests, diagnostics or experiments. The [claim map](claim_map.md) identifies dependencies and inherited qualifications. Historical evidence is distinguished from new commentary in [manifest.json](manifest.json). The [erratum](errata.md) corrects a numerical-reporting scope error without changing frozen files. This preparation and its separate-agent checks are internal static reviews, not external peer review or independent certification of every inherited argument.

The proof record developed in stages. Earlier documents leave the global problem or attainment open; later documents address a narrower local statement or replace the anchor representation. Those frozen statements are not silently rewritten. The map tells the reviewer which conclusion is inherited, which is superseded locally, and which remains open. In particular, the compact construction uses another anchor representation; equality of its induced influence with the earlier one is unproved, and their variances have not been ordered.

## Complete experiment and target

There are exactly five independent observed response groups, indexed by \(a=0,1,2,3,4\), each containing \(n\) independent observations; \(N=5n\). Unobserved \(Z_{ai}\) are independent standard Gaussians. Within each response the two components share its latent draw; no cross-action pairing is observed. The known warp is fixed at
\[
h(z)=z+.1z\sqrt{1+z^2}.
\]
The response maps are
\[
\begin{aligned}
F_0(z)&=g_1(z)+g_2(hz),\\
F_i(z)&=g_1(z+\alpha_i)+g_2(hz+\beta_i),\quad i=1,2,\\
F_3(z)&=g_1(z+1)+g_2(hz),\\
F_4(z)&=g_1(z)+g_2(hz+1).
\end{aligned}
\]
The two private anchor loadings equal one and are known. Both profiles are unknown real-analytic functions \(\mathbb R\to(0,\infty)\), satisfying exactly
\[
g_j(0)=1,\quad g_j(-\infty)=0,\quad
\tfrac12e^x\le g'_j(x)\le2e^x,\quad
|g''_j(x)|,|g'''_j(x)|\le10e^x.
\]
All four coefficients are unknown and belong to the original box \([1/2,1]^4\). The retained local theorem further restricts truths to
\[
\mathcal U=\left\{\eta:
\max_{j=1,2;\,k=0,\ldots,3}\sup_x
e^{-x}|g_j^{(k)}(x)-e^x|\le.01,\
\|c-c_0\|_\infty\le.01\right\},
\]
where, in pair order,
\(c=(\alpha_1,\beta_1,\alpha_2,\beta_2)\) and
\(c_0=(2/3,3/5,3/4,4/5)\).
The public search rectangles are
\(\mathcal R_i=[\alpha_{i0}-.02,\alpha_{i0}+.02]\times
[\beta_{i0}-.02,\beta_{i0}+.02]\).
These structural assumptions are stipulated information, not learned calibration.

For \(\Delta_{u,v}g(x)=g(x+u+v)-g(x+u)-g(x+v)+g(x)\), the target is
\[
I_{12}(\eta)=E\{\Delta_{\alpha_1,\alpha_2}g_1(Z)
+\Delta_{\beta_1,\beta_2}g_2(hZ)\}.
\]
It is the original-scale mean joint-minus-single-minus-single-plus-reference contrast. The joint-action law is unobserved. Inputs are only the five response lists, their sizes, and the public model information. True profiles, coefficients, derivatives, latent realizations, joint outcomes and target labels are excluded.

## Strongest retained theorem and quantifiers

The [compact theorem](evidence/compact_attainment/v1/theorem.md) states that one specified source-data statistic satisfies, for **every fixed** \(\eta\in\mathcal U\),
\[
\widehat I_N-I_{12}(\eta)
=\frac1N\sum_{a=0}^4\sum_{i=1}^n
\psi^{\mathrm{compact}}_{\eta,a}(Y_{ai})
+o_{P_\eta}(N^{-1/2}),
\qquad
V_\eta=\frac15\sum_aE_\eta\psi_{\eta,a}^2<\infty .
\]
Thus its pointwise limit is centered Gaussian with variance \(V_\eta\). This generally nonminimum-norm influence uses a different anchor representation and is not identified with the earlier positive-tail-anchor influence.

At the exponential reference, and at local truths with the required margin, the saved DQM, joint CLT and change-of-measure argument establish regularity along stipulated **fixed** admissible contiguous \(N^{-1/2}\) paths. Profile tangents are analytic, satisfy \(r_j(0)=0\), and have bounded weighted derivatives through order three; all four coefficient directions are included. The full score closure is retained. A finite direction bank is not the model.

The theorem does not give a neighborhood-uniform CLT, arbitrary changing-path regularity, \(O(1/N)\) MSE, studentized coverage, efficient attainment or global minimax recovery. It also does not resolve the original global Hellinger modulus. Constants depend on the fixed warp; nothing is uniform as its nonlinear coefficient tends to zero. These distinctions remain even though every fixed truth in the declared closed neighborhood has the stated pointwise expansion.

The underlying full score is also explicit. For a path response derivative \(q_a=\dot F_a\), its fixed-outcome score is
\[
(Sv)_a(F_a(z))=zq_a/F'_a-q'_a/F'_a+q_aF''_a/(F'_a)^2.
\]
The allocation norm is \(\|f\|_\pi^2=\frac15\sum_aE_af_a^2\). The constructed influence pairs with this score to the derivative of \(I_{12}\) for every admissible profile and coefficient direction, and continuity extends this equality to the full score closure. This is a statement about individual source \(L^2\) functions; independent-source tail divergences cannot be canceled against one another.

That representer is an analysis object. Its existence neither gives the unknown derivatives to the statistic nor makes its norm minimal. The efficient norm is the infimum over all full-compatible representers. The finite dictionaries used later restrict the testing submodels to obtain lower bounds; they do not restrict the original model to establish an upper bound.

## Source-data construction and stochastic dependencies

Monotonicity identifies \(F_a(z)=Q_a(\Phi z)\) from each source law. The statistic substitutes the continuous polygon through observed order statistics at probabilities \(i/n\), constant below \(1/n\), without fitting a profile derivative.

The model-specific algebra uses \(T_jx=h^{-1}(hx+j)\). Two private-anchor telescopes give
\(E_{45}(F;x)=g_1(Cx)-g_1(x)\), where \(C=T_5-4\).
The interval \(J=[-7/4,-3/4]\) is invariant and
\(0<C'<\lambda=271/290<1\).
Normalization fixes
\(g_1(-1)=1-(F_3-F_0)(-1)\).
Differences of geometric iterates reconstruct the needed profile functionals. Signed target integration and outward lattice sums keep individual source influences square integrable.

The fixed schedules are
\[
K_N=\left\lceil2\log N/|\log\lambda|\right\rceil+2,\quad
R_N^2=1.5\log N,\quad
M_N=\lceil4R_N+30\rceil+2.
\]
For each coefficient pair, calibration matches source-\(i\) quantiles at latent arguments \(\pm1\). The algorithm starts at the public center and performs \(2k_N\) projected updates, where
\(k_N=\lceil\log N/|\log(3/4)|\rceil\).
Its preconditioner is computed from public reference exponentials and the known warp. A lexicographically selected minimum-residual public grid supplies a data-defined proof comparator. Neither procedure selects a root nearest an unknown truth.

The inherited [local-attainment proof](evidence/local_attainment/v1/theorem.md) is essential. Its interval-concentration and inverse-CDF argument controls accumulated compact quantile remainders by \(O_P(N^{-3/4}\log^2N)\). The population series bias is \(O(N^{-2})\). The tail area identity, cutoff jumps and polygon interpolation yield population tail bias \(N^{-.6+o(1)}\) and nonlinear/interpolation remainder \(N^{-.85+o_P(1)}\). The omitted linear-tail average is \(N^{-.725+o_P(1)}\). These compatible asymptotic bounds have potentially large constants.

Moving calibration arguments remain within verified fixed integer-transport charts. Global continuity of the finite empirical reconstruction would be false and is not needed. Population contraction and empirical stochastic equicontinuity give the joint four-coefficient expansion. True coefficients, Jacobians, probability-integral transforms and influence functions appear only in analysis. Source-zero terms shared by target integration and calibration are combined before squaring; no fictitious independence or sample splitting removes their covariance.

## Arithmetic theorem and implementation boundary

The earlier deterministic arithmetic construction gives \(O(N^{-1})\) agreement with the finite compact iteration statistic on **promised capped inputs** \(0<Y\le N^2\), or with original cap membership explicitly supplied. It specifies \(N\,\mathrm{polylog}N\) bit work and sufficient \(O((\log N)^2)\) input/working precision, conditional on certified elementary primitives, inverse-warp evaluations, polygon queries, knot enclosures, gap exclusion, quadrature and summation. Input-reading costs are included.

Exact classification of an arbitrary real observation at a cap boundary is not a universal finite-bit operation. The [arithmetic addendum](evidence/practical_feasibility/v1/arithmetic_addendum.md) preserves this qualification. The delivered sample-disabled `mpmath` components lack a complete certified knot compiler, primitive enclosure graph and error allocator. Their diagnostic values do not certify an executable backend.

The later [numerical-stability result](evidence/practical_feasibility/v1/numerical_stability.md) proves a different comparison:
\[
\sqrt N(\widetilde I_N-\widehat I_N)\longrightarrow_P0.
\]
Uniform computed-update discrepancy \(O(N^{-1})\), separately certified target integration, and explicit polynomial input tolerances suffice. The proof uses two-stage **population** contraction and the established empirical modulus, not empirical contraction. It also proves the finite target's coefficient sensitivity is \(O_P(1)\).

This removes the deterministic iteration-amplification factor from the statistical precision allocation, permitting sufficient \(O(\log N)\) precision under the same certified-primitive model. An encoded positivity/cap gate replaces exact-real classification, with negligible failure probability. It does not yield deterministic agreement on every capped dataset, recover missing digits in fixed-resolution measurements, certify ordinary floating point, or demonstrate affordable constants. The final clipping and fallback preserve asymptotic conclusions, not an MSE rate.

## Contribution and negative evidence

Trabs's adjoint-range theorem supplies the exact regularity and information framework; constructing a compatible influence is not itself attainment. Kaji supplies integrated quantile and L-statistic theory, but its cited fixed-Lipschitz-bound result does not automatically cover this growing series, changing tail weights and estimated coefficients. The additional remainder analysis must stand on its own.

DExtrI's main Theorems 3.3–3.4 and broader Theorem D.3 concern identification/extrapolation under different information and rigidity assumptions. Theorem D.9 supplies another appendix rigidity regime. Five isolated laws with known Gaussian noise, warp and private unit anchors are not simply more general than its axis-interval experiment. Being outside a main exponential-polynomial assumption does not establish novelty. Established functional estimation likewise already separates target recovery from full nuisance recovery; Bennett and collaborators' conditional-moment source conditions do not transfer automatically to this nonlinear density model.

The potential insight is therefore the explicit contracting-anchor solution, complete coefficient compatibility and compatible stochastic/arithmetic construction in this particular experiment. It is interesting as an infinite-dimensional missing-combination problem with a locally regular target. Its importance beyond this strongly calibrated example remains unclear. Telescoping, adjoint inversion, nuisance projection and one-step correction are established operations, not renamed ML methods.

The frozen adaptive pilot improved aggregate MSE by **3.5569%** against the strongest complete control, below its **10%** screen. Its [later adversarial review](evidence/learned_target_v1/adversarial_review_20260912.md) retains rank, approximation and implementation qualifications. No mechanism advantage was established.

The subsequent reference influence calculation gives approximately \(V=79{,}390{,}122.79\). Known-median projection removes approximately **0.081285%**, leaving \(U\approx79{,}325{,}590.57\). The twelve-direction information lower bound is approximately \(2{,}045.74\). Their approximately **38,776-fold** ratio does not demonstrate an achievable improvement: neither endpoint is proved close to efficient variance. Numerical integration is uncertified; analytic truncation bounds are separate. Asymptotic SE extrapolations are not measured performance or necessary sample-size bounds.

Backend no-go remains. No practical certified implementation, useful finite-sample accuracy, ML superiority or biological siRNA calibration is established. A frozen-GNN interaction would describe that predictor, not validate this biological target.

## References

- [Trabs, Theorem 2.7, v2](https://arxiv.org/html/1307.6610v2).
- [Kaji, Theorem 4.1, Proposition 4.2 and Theorem 5.1, v1](https://arxiv.org/html/1910.07572v1).
- [DExtrI, main and appendix theorems, v1](https://arxiv.org/html/2608.19849v1).
- [Bennett et al., Theorem 2 and §4.1, v3](https://arxiv.org/html/2208.08291v3).

## Questions for the reviewer

1. Are the local attainment claim and its complete stochastic remainder valid under the stated original model and tangent class?
2. Does the computational result accurately state its arithmetic and observed-input assumptions?
3. Which parts, if any, have an exact antecedent in the closest literature?
4. Does the contracting-anchor construction provide a sufficiently substantive insight beyond this calibrated example?
5. Is there a defensible theorem-focused paper in the existing record, despite the practical no-go, and what essential objection would prevent that claim?

