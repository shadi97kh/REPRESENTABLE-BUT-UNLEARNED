# Source-only finite estimator and implementation boundary

**PROVED HERE:** The following is one data-defined statistic. Its expansion is proved in [theorem.md](theorem.md). It estimates neither profile derivatives nor response-density derivatives. The reference centers define the public search region; their unknown true values are not used.

## Definition

1. Accept five ordered groups of positive scalar responses, corresponding to \(0,e_1,e_2,e_3,e_4\), with equal size \(n\), \(N=5n\). No joint-action outcomes or latent observations are accepted. Fewer than two observations per group returns zero. Nonpositive/nonfinite responses return zero as an explicit off-model fallback. If any response exceeds \(N^2\), return zero.
2. Sort within each source. The auxiliary empirical inverse CDF is \(Q_n(p)=Y_{(\lceil np\rceil)}\) on \(0<p\le1\), with \(Q_n(0)=Y_{(1)}\). The statistic uses the polygon through \((i/n,Y_{(i)})\), constant on \([0,1/n]\) and at 1. Thus, for \(p\in[k/n,(k+1)/n]\), \(1\le k<n\), use \((1-t)Y_{(k)}+tY_{(k+1)}\), \(t=np-k\). There is no fitted interpolation parameter. Set \(\widetilde F_a(z)=\widetilde Q_a(\Phi z)\).
3. Use \(K_N=\lceil2\log N/|\log(575/613)|\rceil+2\). Form the finite contracting difference \(D_K\) from the six-term anchor response identity. Start profile normalization with **1**, as in theorem (5). To evaluate \(g_1(x)\), set \(j=\lceil4-x\rceil\), reconstruct at \(x+j\in[4,5]\), and transport back with the finitely many differences \(\widetilde F_3-\widetilde F_0\). Set \(\widehat g_2(y)=\widetilde F_0(h^{-1}y)-\widehat g_1(h^{-1}y)\).
4. For each pair \((\alpha_i,\beta_i)\), search its public rectangle of radius \(1/50\) around \((\alpha_{i0},\beta_{i0})\). Divide each side into \(m_N=\lceil N/25\rceil\) equal intervals, so mesh is at most \(N^{-1}\). At each pair, evaluate residuals \(\widetilde F_i(z)-\widehat g_1(z+\alpha)-\widehat g_2(h(z)+\beta)\) at \(z=-1,+1\). Select the first lexicographic minimum of their Euclidean norm. A deterministic approximate comparison within objective error \(N^{-2}\) is permitted. Store coefficients in order \(\alpha_1,\alpha_2,\beta_1,\beta_2\).
5. Use \(R_N=\sqrt{1.5\log N}\), \(M_N=\lceil4R_N+30\rceil+2\), and the outward weights from theorem (6), evaluated at the selected coefficients. Truncate the outer integral to \([-R_N,R_N]\), the positive and negative outward lattices to \(M_N\) terms, and the central periodization to indices \(-M_N,\ldots,M_N\). Evaluate the central term with \(D_K(t,4)\). The truncated periodization does not integrate exactly to zero, so the subtraction of the same base point is retained explicitly.
6. Integrate by midpoint quadrature with mesh at most \(N^{-6}\), splitting the outer interval at 4 whenever it is inside the interval; integrate the central term on \([4,5]\). Clip the result to \([-1000,1000]\).

This defines a statistic even if reconstructed profiles are nonanalytic, nonmonotone, or fail normalization at points other than the specified anchor normalization. It never fits these profiles to the original class. Positive-class population identities and empirical-process bounds justify its behavior locally.

The root rule is measurable and handles competing candidates without an oracle. The theorem proves a uniform inverse-Lipschitz bound for the **population** calibration map on each search rectangle, preliminary consistency, moving-argument stochastic equicontinuity, and an approximate-root expansion for the actual grid minimizer. A Jacobian determinant alone would not prove those steps.

## Arithmetic and finite implementation

[estimator.py](estimator.py) implements the polygon, known warp/inverse, finite anchor/profile operators, finite target weights, periodization, integration, schedules, and the coefficient grid rule. The public function `estimate_from_sources(groups)` raises `PermissionError` before reading the input because sampled-data execution is disabled. No call processed source data in this continuation.

The supplied backend is `mpmath` and is explicitly **not** a certified interval library. The mathematical statistic is defined in exact arithmetic; theorem E proves negligible numerical error for a realizable certified-arithmetic version with:

- absolute input/known-primitive error \(N^{-40}\);
- inverse bracket width \(N^{-40}\);
- residual objective error at most \(N^{-2}\);
- the stated \(N^{-6}\) quadrature mesh;
- working precision \(\lceil10000+200\log_2N\rceil\) bits, with certified enclosures refined to meet the absolute-error requirements.

Only known functions require these elementary enclosures. This is not access to an unknown density/profile oracle. Rational bisection and Taylor remainder bounds give a terminating realization. The reference backend has not been certified against this arithmetic specification, so a successful `mpmath` diagnostic is not an interval error certificate.

On \(\max Y\le N^2\), the polygon is \(N^3\)-Lipschitz in probability. This controls approximate \(\Phi(z)\) evaluation even when a query is nearly an empirical knot; exact transcendental rank comparison is unnecessary. Coinciding rounded contracting queries can otherwise erase a small difference, which is why paired subtraction and a vanishing absolute error are required. Fixed float64 arithmetic eventually rounds extreme probabilities to 0 or 1 and cannot implement the asymptotic schedule.

The finite grid has \(O(N^2)\) candidates per coefficient pair. Its direct implementation takes \(N^{2+o(1)}\) primitive operations. Midpoint integration with reconstruction/lattice calls costs \(N^{6+o(1)}\) primitive operations and dominates. Storage can be \(O(N+\log N)\) by streaming sums. These are counts under elementary-function arithmetic; bit costs are additional. Finite implementability does not imply practical runtime.

## Fallback and practical onset

The response-cap fallback probability is at most
\[
128e^4(e^8+\sqrt5e^{40.2})N^{-7}.
\]
The multiplier is very large. It is a uniform moment bound, not a measured failure rate. The target is strictly inside the clipping interval under the original model, and consistency makes final clipping inactive with probability tending to one.

The public sufficient domain condition is
\[
n\min\{\phi(7)7/50,\ \phi(R_N+1)(R_N+1)/(1+(R_N+1)^2)\}
\ge64(\log N)^2.
\]
This is a diagnostic for the proof's relative-quantile control, **not** an operational fallback. The finite polygonal statistic still exists when it fails. It does not certify useful finite-sample accuracy when it passes.

**NUMERICALLY CHECKED:** Evaluating that known inequality produced:

| Observations per group \(n\) | Lower effective rank | Required rank | Sufficient condition |
|---:|---:|---:|---|
| \(10^{12}\) | 0.0124106 | 54,720.3 | fails |
| \(10^{32}\) | 14.9798 | 362,810.2 | fails |
| \(10^{64}\) | 14,554,195.8 | 1,420,385.0 | passes |

No dataset of these sizes was created. These three calculations examine only the deterministic schedule. They are not statistical lower bounds on the necessary sample size. Together with the \(N^{6+o(1)}\) construction, they make this an attainment/existence result with severe practical limitations.

The other new checks verify the polygon endpoints/knots, population \(g_1(0)=g_2(0)=1\), the sign of a finite population contraction remainder, and the disabled sample entry. The finite identity discrepancy was \(1.5244246148248607861\times10^{-36}\) at 70 decimal digits. The saved local-regular check record was read and reused; its suite was not rerun. [Actual records](checks.json).

## Subsequent validation specification — prepared only

**NOT EXECUTED / NOT AUTHORIZED HERE.** The smallest statistically relevant follow-up would compare the *actual estimator* with its oracle linear expansion, not train a network or benchmark generic prediction:

- One generating reference: the exact reference in the theorem. One fixed local alternative: \(g_{1,N}(x)=e^x[1+N^{-1/2}\sin x]\), \(g_{2,N}=e^x\), with coefficient displacement \(N^{-1/2}(.2,-.1,.15,-.05)\). Use only counts large enough that these are admissible.
- Freeze \(n_0\) as the smallest power of two meeting the displayed domain condition, and use \(n_0,2n_0,4n_0\), 20 independent repetitions each. Use the same deterministic schedules and five-source information in every repetition. The unknown truth is available only to the generator/evaluator.
- Record \(N^{1/2}(\widehat I-I)\), the oracle linear sum, their difference, each coefficient error, cap/clipping frequency, numerical certificate width, and total computation. Combine all reference contributions in the oracle influence before variance calculation.
- Stop each case at its separately authorized CPU/wall cap; retain every timeout, fallback, and numerical failure. Do not select a favorable truncation or reference afterward. No empirical efficiency claim without comparison to the actual efficient influence bound.

This specification is presently impractical at its sufficient counts and direct computational cost. It is not a request to allocate or run it. A practical implementation would first need a separately justified faster evaluation/selection method with the same approximation guarantee; replacing the estimator by an unproved approximation would not validate this theorem.

**OPEN:** consistent feasible variance estimation, studentization, efficiency, MSE rates, and useful finite-sample performance. The primary estimator does not require derivative estimation; the oracle influence's appearance in its analysis does not make its variance automatically observable.


## Secondary variance attempt: explicit sufficient conditions

**PROVED HERE, conditional statement; nuisance construction OPEN.** A variance estimate need not alter the primary estimator. Use a fixed index split within each source only for this secondary calculation. On the training half construct an estimate \(\widehat\psi_a\); on the independent evaluation half of size \(m\), clip it at \(L_m=m^{1/8}\) and compute
\[
\widehat V=\frac15\sum_a\left\{\frac1m\sum_i
 [\widehat\psi_a(Y_{ai})]_{L_m}^2-
 \left(\frac1m\sum_i[\widehat\psi_a(Y_{ai})]_{L_m}\right)^2\right\}.
\]
If \(\sum_a E_\eta[(\widehat\psi_a-\psi_{\eta,a})^2\mid\text{training}]=o_P(1)\), then \(\widehat V\to_P V_\eta\). Indeed clipping is a contraction in \(L^2\) and its true-influence tail tends to zero. Conditional second-moment sampling error has variance at most \(L_m^4/m=m^{-1/2}\); conditional means converge to zero. This proves the implication with only finite true influence variance, not an assumed fourth influence moment.

A sufficient way to meet that condition is to truncate the influence's tail/series at public \(T_m,K_m\to\infty\), retain the known uniform Gaussian/geometric tail bounds, and consistently estimate the source curves, their first derivatives, coefficients, and the two calibration derivative matrices on all retained query intervals. More explicitly, let their respective sup/operator errors be \(\epsilon_F,\epsilon_{F'},\epsilon_c,\epsilon_M\), and require estimated calibration minimum singular values at least half the public lower bound. A sufficient finite-kernel error condition is
\[
D(T_m,K_m)\{\sqrt{\epsilon_F}+\epsilon_F+\epsilon_{F'}
                   +\sqrt{\epsilon_c}+\epsilon_c+\epsilon_M\}=o_P(1).
\]
Here a deterministic \(D(T,K)\) is obtained by summing absolute finite kernel coefficients, their known compact derivative envelopes, and the matrix inverse bound. The square roots account for moving indicator thresholds under the **true response law**, not merely a change of Gaussian coordinates. Source densities are bounded by the original derivative lower bound, so threshold displacement \(\epsilon_F\) changes indicator \(L^2\) norm by at most \(C\sqrt{\epsilon_F}\). The remaining factors follow from finite sums, smooth coefficient weights and the inverse-matrix identity. Uniform infinite-series/tail bounds control the discarded influence pieces.

This states exactly what additional nuisance accuracy is sufficient. No rates for estimating these derivatives/matrices on growing intervals were established in this pass, and no such estimates were computed. Therefore a feasible consistent variance estimator is not included among the completed claims. The main root-\(N\) estimator remains derivative-free.
