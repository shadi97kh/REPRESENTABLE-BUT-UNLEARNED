# Adversarial static review of local_attainment/v1

## Decision and review independence

**The exact-arithmetic, model-specific local attainment claim survives this static review. A certified sampled-data implementation and a substantive ML contribution have not been established.**

I found no fatal error in the growing empirical-quantile remainder, signed tail integration, source-only coefficient selection, or joint influence calculation for the statistic and neighborhood actually specified. This is a finding from inspection, not independent certification of the theorem. Two qualifications require explicit treatment: the transported empirical profile is not globally continuous across its integer charts, and the claimed fixed working-bit count is not backed by a certified implementation or a complete precision allocation.

This is a **static self-review by the same AI system**, with two delegated reviews of different parts of the proof. It is not external independent peer review, a formal proof check, or a replication. No numerical jobs, samples, model fits, bootstrap, experiments, previous checks, or figure renderers were run. Original proofs, code, configurations, predictions, diagnostics and ledgers were not changed.

The review read the original model, earlier full five-law influence proof, attainment theorem, estimator specification and code, checks and their source, saved final checks, novelty comparison, and cached primary sources. The existing Kaji, Trabs and Bennett source pages were also checked read-only. No broad novelty search was performed.

## Strongest supported theorem

Keep the actual declared model: known standard Gaussian latent law; known
\[
h(z)=z+\tfrac1{10}z\sqrt{1+z^2};
\]
two known private unit anchors; profile normalizations; two unknown analytic profiles in the original class intersected with the weighted-\(C^3\) radius-.01 neighborhood of \(e^x\); all four coefficients unknown within radius .01 of the stated centers; and five independent groups of size \(n\), \(N=5n\). These substantial structural restrictions are assumptions, not conclusions learned from the observations.

For the exact-arithmetic statistic defined in [estimator.md](v1/estimator.md), the inspected proof supports, at every fixed truth in that declared set,
\[
\widehat I_N-I(\eta)
=N^{-1}\sum_{a=0}^4\sum_{i=1}^n\psi_{\eta,a}(Y_{ai})
+o_{P_\eta}(N^{-1/2}),
\qquad
V_\eta=\tfrac15\sum_aE_\eta\psi_{\eta,a}^2<\infty.
\]
The influence is the explicit, generally nonminimum-norm five-source solution. The same centered limit is supported along the stipulated fixed admissible \(N^{-1/2}\) paths at the reference, and at truths with the stated neighborhood margin. This is not a uniform assertion for arbitrary direction sequences or for shrinking warp separation.

The remaining gap is not an identified missing stochastic remainder in Parts A–D. It is the numerical/implementation claim discussed below, plus explicitly open efficiency, inference, MSE and practical-performance questions. A failed implementation or the absence of a practical regime would not itself refute this ideal-statistic theorem.

## 1. Observed inputs, known assumptions, and oracle quantities

The input audit passes. The primary interface accepts only five balanced response lists. The group labels, Gaussian law, fixed warp, private-anchor unit loadings, profile normalization, coefficient search centers/radii, clipping rule and schedules are public assumptions or deterministic choices. There are no observed latent draws or joint-target labels.

The transforms \(U_{ai}=P_{\eta,a}(Y_{ai})\), true response derivatives, true calibration matrices, \(b_{\eta,c}\), and \(\psi_\eta\) are proof/evaluator quantities. They are not read by the estimator. Their use in an asymptotic expansion is not an oracle-input violation. They would become a feasibility problem if supplied to a purported variance estimator without estimation; the secondary variance section correctly leaves that nuisance construction open.

In [estimator.py](v1/estimator.py), lines 228–229 block sampled execution before input access. The finite operators are implemented, but the public entry currently returns no estimate. The proposed statistical algorithm is specified; a running certified learner is not demonstrated.

Normalization is also scoped correctly: the finite construction makes \(\widehat g_1(0)=1\), but generally
\[
\widehat g_2(0)=\widetilde F_0(0)-1\ne1.
\]
The saved population fixture has both values equal to one. It does not establish both empirical normalizations. The specification explicitly allows off-model empirical profiles, so this observation is not a counterexample to the expansion. It may introduce avoidable noise and supplies no efficiency claim.

## 2. Coefficient selection is data-defined

The two finite grids have public centers and radius \(1/50\), with \(m_N=\lceil N/25\rceil\) subdivisions per side. Their spacing is at most \(N^{-1}\). The first lexicographic minimum of the observed residual norm is measurable; permitted objective approximations have a stated error. No rule selects the root closest to the unknown truth.

The argument in [theorem.md](v1/theorem.md), lines 175–212, establishes more than an invertible Jacobian at one point. Averaging the two derivative columns along the candidate segment yields
\[
\det\overline M_i\ge
e^{\alpha_{i0}+\beta_{i0}-.04}
\{.99^2e^{h(1)-1}-1.01^2e^{1-h(1)}\}>0.
\]
Together with the entry upper bound, this separates every two population calibration values throughout the public rectangle. The true coefficients have a .01 margin from its boundary.

The proof's hypothetical \(c_i+M_i^{-1}\widehat R_i(c_i)\) is used only to exhibit a nearby small-residual competitor. It is not the estimator's selection rule. Preliminary consistency, local stochastic equicontinuity, a nearby grid point, and the minimum-residual property give the displayed joint coefficient expansion. Residual signs and the code's ordering \((\alpha_1,\alpha_2,\beta_1,\beta_2)\) agree.

## 3. Concrete scope defect: integer-chart discontinuities

The moving-argument assertion at theorem lines 108–113 must not be read as a global modulus for the *transported* empirical profile. [estimator.py](v1/estimator.py), lines 133–144, uses the chart index \(\lceil4-x\rceil\). At \(x=4\), continuity of the data curves within each chart gives
\[
\widehat g_{1,K}(4^-)-\widehat g_{1,K}(4)
=D_K(\widetilde F;5,4)
-[\widetilde F_3(4)-\widetilde F_0(4)].
\tag{R1}
\]
This is generally nonzero for off-model empirical quantile curves. Its population value is only a contraction remainder, but its centered first-order fluctuation need not disappear at the root-\(N\) scale. The source-3 influence contains the isolated term \(-\kappa_{3,4}\); the contracting terms use thresholds \(C^{k+1}4,C^{k+1}5>4\). They therefore cannot remove its jump at \(F_3(4)\). An unrestricted shrinking-interval modulus across this join is false in general.

**Why this does not invalidate the specified estimator:** the negative calibration arguments stay strictly in \((-1,0)\), and their positive counterparts strictly in \((1,2)\), throughout the already declared public rectangles. The earlier [local proof](../local_five_law/v1/proof.md), line 369, checks these fixed charts. The central target uses \(D_K(t,4)\) directly on \([4,5]\), with no chart transport. The needed equicontinuity is within those charts. This review does not impose a smaller profile class or coefficient region to obtain that fact.

## 4. Growing series and tail errors are substantively addressed

The proof does not infer a growing-series expansion from convergence at fixed depth. Uniform compact empirical-quantile linearization has remainder \(O_P(N^{-3/4}\log N)\); multiplying by \(K_N=O(\log N)\) gives \(O_P(N^{-3/4}\log^2N)\). The influence tail has norm \(O(\lambda^{K_N/2})\). For the required moving queries, interval empirical-process increments are summed geometrically along the contraction. These are the relevant controls, rather than a finite score dictionary.

The signed weighted-quantile identity also has the correct sign. Writing \(w_R=U/\phi\) in probability coordinates and \(V_R(p)=\int_0^p w_R\),
\[
\int w_R(Q_n-Q)
=-\int_0^\infty\{V_R(P_n(y))-V_R(P(y))\}\,dy.
\]
Changing variables \(y=F(z)\) gives theorem (13). Finite moments and the truncated integrable weights justify this equality. It is first applied to the step quantile; the polygonal interpolation error is bounded separately.

The smooth nonlinear term involves \(F'\sup|w_R'|\Delta_n^2\). Combining the actual Gaussian envelopes leaves
\[
N^{-1}\exp\{.1R_N^2+O(R_N)\}\operatorname{poly}(\log N).
\]
The crossing strips at the two cutoffs and the jump at 4 are explicitly bounded; they are not omitted. With \(R_N^2=1.5\log N\), this is \(N^{-.85+o(1)}\).

The schedules are compatible at the level of the proved asymptotic bounds:

| Contribution | Saved order | Review |
|---|---:|---|
| Central population truncation | \(N^{-2}\) | Negligible |
| Summed compact quantile remainder | \(N^{-3/4}\log^2N\) | Negligible |
| Population response-tail bias | \(N^{-.6+o(1)}\) | Negligible |
| Tail nonlinear and interpolation error | \(N^{-.85+o(1)}\) | Negligible |
| Omitted linear-tail average | \(N^{-.725+o_P(1)}\) | Negligible |
| Lattice omission under response cap | \(N^{-10+o(1)}\) | Negligible |
| Exact midpoint integration error | \(N^{-3+o(1)}\) | Negligible |

The effective tail ranks grow as \(N^{1/4+o(1)}\). The inequalities \(c<2\), \(.4c>.5\), and \(.1c<.5\) have the common solution \(c=1.5\). I found no incompatible exponents or an assumed desired quantile remainder here. The hidden constants nevertheless matter greatly for practicality.

## 5. Shared-source covariance and estimated-coefficient effects

The coefficient derivative is included:
\[
\psi_\eta=\psi_{A,\eta}
+\sum_{j=1}^4 b_{\eta,c,j}\chi_{\eta,j}.
\]
Anchor curves do not themselves depend on the four coefficients, so this decomposition does not miss an indirect coefficient effect on anchor fitting. Candidate coefficients do change target weights; the proof bounds their empirical increment by the coefficient error times an \(O_P(N^{-1/2})\) integral, and includes the population derivative \(b_c\). Moving arguments in coefficient calibration are handled separately.

Each quantile representer has allocation factor \(1/\pi_a=5\), matching the \(1/N\) sum. All source-zero contributions, including calibration, are combined before squaring. The proof does not pretend that repeated uses of the reference sample are independent. The individual-source \(L^2\) bounds precede the variance calculation; finiteness is not obtained by cancelling variances between independent sources.

## 6. Local alternatives and scope

The inherited paths are fixed analytic directions \(r_j\) with \(r_j(0)=0\) and \(|r_j^{(k)}(x)|\le C_v e^x\), \(k=0,\ldots,3\), together with any fixed four-dimensional coefficient direction. Sufficiently small two-sided perturbations at the reference remain in the original class; at neighborhood-interior truths the stated margin permits the same argument.

The earlier DQM argument checks a polynomial Gaussian envelope for the square-root-density derivatives. Balanced-product LAN, a joint CLT for the score and influence, the full score equation, and contiguity then justify the fixed-path local limit. Contiguity transfers the *same statistic's* remainder; it does not require a new uniform expansion over moving directions. It would be an overclaim to extend this to unrestricted \(N\)-dependent directions, neighborhood-uniform coverage, or \(\delta\to0\). The theorem expressly avoids those extensions.

## 7. Concrete precision overstatement and exact remaining implementation gap

The deterministic midpoint error is credible for the ideal statistic: on the response-cap event the polygon has slope at most \(N^3\) in probability, the finite integrands have the stated polynomial/logarithmic amplification, and the discontinuity at 4 is split. Known-function bisection and rational enclosures do not require unknown profiles or density oracles.

However, theorem line 254 describes
\[
p_N=\lceil10000+200\log_2N\rceil
\]
as a sufficient working-bit schedule **and** permits enclosure refinement if needed. Those are different guarantees. The artifacts do not give an implemented interval algorithm with an operation-by-operation error allocation proving that this exact fixed bit count suffices. Generic computability and a requested absolute tolerance do not certify the delivered floating-point program.

The actual code uses `mpmath`, not outward-rounded interval bounds. It does not return numerical error enclosures, and sampled execution is disabled. The existing analytic fixture checks cannot close that implementation claim.

The concrete outstanding obligation for a finite-precision realization is to establish, without assuming the desired conclusion,
\[
\sqrt N\,|\widehat I_N^{\rm implementation}
                 -\widehat I_N^{\rm exact}|\longrightarrow0
\quad\text{in probability},
\tag{R2}
\]
including a verified uniform objective-approximation bound on the finite grids, such as
\[
\sup_{c\in\Gamma_N}|\widetilde{\operatorname{cost}}_N(c)
                         -\operatorname{cost}_N(c)|\le N^{-2},
\tag{R3}
\]
and the inverse, interpolation, summation and integration errors. If more than \(p_N\) bits are needed, an adaptive certified bound and its actual complexity must be stated instead of claiming the fixed schedule was proved. This is a specific precision/implementation gap, not a counterexample to the exact-arithmetic pointwise statistical expansion.

## 8. CLT, MSE, efficiency, and practical regimes

The saved theorem does **not** promote its CLT to \(O(1/N)\) MSE or efficiency. It explicitly leaves uniform integrability, efficient projection, feasible variance estimation and studentization open. Bounded final clipping alone would not supply the missing scaled-error uniform integrability.

The secondary split-sample variance implication is correctly conditional: if an estimated influence converges in true-source \(L^2\), clipping at \(m^{1/8}\) bounds the conditional variance of its squared sample average by \(m^{-1/2}\). The artifact leaves construction of such an influence estimate open. It is not an implemented variance estimator.

The saved diagnostics report:

| Per-source count | Lower effective rank | Required rank | Sufficient condition |
|---:|---:|---:|---|
| \(10^{12}\) | 0.0124106 | 54,720.3 | Fails |
| \(10^{32}\) | 14.9798 | 362,810.2 | Fails |
| \(10^{64}\) | 14,554,195.8 | 1,420,385.0 | Passes |

These values were read from [checks.json](v1/checks.json), not recomputed. Passing this one sufficient condition does not show small total error, useful uncertainty, or feasible runtime. The direct operation count is \(N^{6+o(1)}\), with bit costs additional. **No demonstrated practical regime exists.** Conversely, these conservative sufficient counts are not necessary sample-size lower bounds.

The same saved check file tests interpolation, population normalization, one finite population remainder (discrepancy \(1.5244246148248607861\times10^{-36}\)), the schedule condition, and the disabled entry guard. It never invokes the grid calibration or complete sampled estimator. Its “passed” status is correctly limited to those deterministic fixtures.

## 9. Closest results and contribution status

- **Trabs:** the adjoint-range criterion and efficient information framework are established theory. The score equation alone does not supply this nonlinear estimator's stochastic remainder. The candidate model-specific content is the explicit contracting anchor construction and attainable source-data functional. [Theorem 2.7](https://arxiv.org/html/1307.6610v2).
- **Kaji:** inverse-map and integrated-quantile asymptotics are established, including random weights. Its positive-density condition allows an interval support, so \((0,\infty)\) is not an obstacle here. The unweighted integrability check is valid. Theorem 5.1 uses a fixed Lipschitz bound and does not automatically cover the growing point-query series and actual tail weights. Parts A–D supply additional, model-specific controls; simply integrating quantiles is not novel. [Theorems 4.1 and 5.1; Proposition 4.2](https://arxiv.org/pdf/1910.07572).
- **Bennett et al.:** recovering regular functionals without well-recovering every nuisance component is established territory. Their conditional-moment/source restrictions are not automatically identical to this nonlinear five-density experiment. No theorem here establishes superiority over an equally informed functional estimator. [Strongly identified functionals](https://arxiv.org/html/2208.08291v3).
- **DExtrI:** the cached population identification theorem uses axis-interval observations and different structural assumptions. Five discrete laws with a known warp, known latent distribution and calibrated private anchors are not an assumption-free strengthening. The local stochastic attainment argument is a distinct model-specific task; exact literature priority remains unestablished. [DExtrI](https://arxiv.org/html/2608.19849v1).

The narrow mathematical distinction survives the bounded comparison. A new general learning principle, a faster rate than equally informed alternatives, or a substantive ML contribution does not follow. Neither this proof nor its synthetic fixtures establish a biological siRNA interaction result.

## Preservation and publication status

Commit `efa772a3c0abf72e1d32ca71c819f2351e8f8366`, the exact previously approved publication payload, was pushed to `shadi97kh/iclr2027` on branch `main`; the remote hash was verified. This review is a fresh local document and was not included in that already approved commit.

No frozen artifact was repaired, regenerated or overwritten, and no resource ledger was extended by a numerical job. Static file/source inspection, this new review text, and Git operations are not being represented as zero CPU work or as an independent mathematical certification.
