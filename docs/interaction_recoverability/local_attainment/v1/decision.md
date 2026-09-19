# Decision: local attainment established

**PROVED HERE:** A single estimator using only five balanced source response groups satisfies
\[
\widehat I_N-I(\eta)=\frac1N\sum_{a=0}^4\sum_{i=1}^n
\psi_{\eta,a}(Y_{ai})+o_{P_\eta}(N^{-1/2}),\qquad N=5n.
\]
Both profiles and all four coefficients remain unknown. The result holds at the specified reference and at each fixed truth in the declared weighted-\(C^3\)/coefficient neighborhood. Fixed-path DQM, joint CLT and change of measure establish regularity under the specified admissible \(N^{-1/2}\) alternatives.

The estimator uses polygonal empirical quantiles, a finite contracting anchor reconstruction with population normalization \(g_j(0)=1\), deterministic minimum-residual calibration on a public coefficient grid, and outward target integration. It does not estimate profile derivatives. The reference supplies only a declared search region, not unknown true parameters.

The completed remainder argument combines:

- \(K_N=\lceil2\log N/|\log(575/613)|\rceil+2\): compact growing-series remainder \(O_P(N^{-3/4}\log^2N)\);
- \(R_N^2=1.5\log N\): target tail bias \(N^{-.6+o(1)}\), tail CDF/interpolation remainder \(N^{-.85+o_P(1)}\);
- a globally separated population calibration map on each public rectangle, a measurable mesh-\(N^{-1}\) minimum-residual rule, and moving-query stochastic equicontinuity;
- \(M_N=\lceil4R_N+30\rceil+2\), quadrature mesh \(N^{-6}\), and explicit vanishing numerical tolerances.

All are simultaneously negligible at the required scale. The actual influence agrees with the nonminimum-norm explicit influence in the saved local theorem, including the combined source-zero covariance.

**Arithmetic qualification:** The mathematical statistic and its certified-arithmetic realization are specified. The delivered sample-disabled `mpmath` code implements the finite operators and estimator interface, but it is **not a certified interval backend**. No sampled-data estimator was executed. The numerical theorem requires the stated primitive/integration accuracy obligations; the saved diagnostics do not certify arbitrary sampled runs.

**Variance and MSE:** \(V_\eta=\sum_a(1/5)E_\eta\psi_{\eta,a}^2\) is finite and supplies the pointwise CLT variance. Efficient variance attainment, a feasible consistent variance estimate, studentization, \(O(1/N)\) MSE, neighborhood-uniform distributional guarantees and global minimax claims remain unproved. The prior global-unresolved decision is unchanged.

**Practical limitation:** The public sufficient quantile-domain diagnostic failed at \(n=10^{12}\) and \(10^{32}\), and passed at \(10^{64}\), without generating any observations. It is a conservative proof-onset condition, not a sample-size lower bound and not an operational zero-output gate. The explicit implementation costs \(N^{6+o(1)}\) primitive operations. This establishes mathematical attainment, not a useful runtime/sample-size regime.

**NUMERICALLY CHECKED:** The new deterministic record verifies interpolation conventions, correct population normalization, a finite population remainder with discrepancy \(1.5244246148248607861\times10^{-36}\), the public schedule conditions, and refusal of sampled-data execution. The previous local checks were read, not rerun. These checks are neither fitted evidence nor statistical validation.

The exact advance over the saved work is an attaining source-data construction and its stochastic remainder analysis. General inverse-functional estimation, empirical-process methods and projection machinery remain known baselines. The bounded [novelty comparison](novelty.md) does not establish a new general learning principle, publication priority, superiority over equally informed controls, or biological validity.

Files: [theorem](theorem.md), [estimator specification and prepared validation](estimator.md), [reference implementation](estimator.py), [deterministic checks](checks.json), [commands and resources](commands_and_resources.md). Subsequent validation is specified only; no experiment or additional authorization is requested.
