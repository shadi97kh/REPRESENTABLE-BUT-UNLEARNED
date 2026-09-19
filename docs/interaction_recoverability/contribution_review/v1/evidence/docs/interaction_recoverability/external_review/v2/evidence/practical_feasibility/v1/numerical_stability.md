# Statistical numerical stability of the compact estimator

**PROVED HERE:** the compact estimator's statistical comparison does not require the frozen worst-case iteration amplification \(A_T\). A uniform computed-update discrepancy \(O(N^{-1})\), an explicit polynomial input-height tolerance, and separately certified \(O(N^{-1})\) target integration suffice for asymptotic equivalence **in probability**. They do not imply deterministic \(O(N^{-1})\) agreement on every capped dataset.

This is a static derivation using the declared five-law model and the already-proved compact empirical-process inequalities. It neither changes the frozen theorem nor certifies the delivered `mpmath` code. No samples, fits, numerical checks or ledger operations were executed for this note. The reasoning was developed by a separate agent within the same system: this is an internal cross-check, not an external independent review.

## 1. Inputs and comparisons

Let \(\widehat T_{i,N}\) be the ideal empirical projected update for pair \(i=1,2\), using the exact positive observations, public matrix \(A_i=M_{i0}^{-1}\), and the finite compact operators. Let \(\widetilde T_{i,N}\) be its numerical realization, evaluated at the same coefficient argument:
\[
e_N=\max_{i=1,2}\sup_{c\in\mathcal R_i}
\|\widetilde T_{i,N}(c)-\widehat T_{i,N}(c)\|.
\tag{1}
\]
The supremum may be understood as a certified error statement at every representable input, since only representable states are executed. Rounding the next state and the public starting point is included. Projection is onto the same public rectangles. No oracle root, profile, derivative or true coefficient is used.

Use the frozen schedules \(K=O(\log N)\), \(R=\sqrt{1.5\log N}\), \(M=O(\sqrt{\log N})\), and \(T=2k_N\), \(k_N=\lceil\log N/|\log(3/4)|\rceil\). The computational schedule can be selected using the rational/outward-enclosure convention of the saved static review; exact transcendental ceiling decisions are unnecessary.

The ideal statistic \(H_N(Y,\widehat c)\) means the **finite compact** target integrated exactly, with the finite projected coefficient iteration. The computed target first uses the rounded heights \(\widetilde Y\), then the computed coefficients \(\widetilde c\), and finally separately allocated numerical integration. The comparison is with this same compact statistic, not the earlier noncompact influence.

## 2. Perturbing the coefficient iteration without empirical contraction

**PROVED HERE.** Write \(T_i\) for the population projected update. The frozen proof establishes its contraction with public bound \(\bar q=3/4\), and the fixed-chart empirical noise has uniform magnitude
\[
\varepsilon_N=O_P(\log N/\sqrt N).
\]
Projection is nonexpansive. If the starting-point error is \(O(N^{-1})\), then
\[
\|\widetilde c_i^t-c_i\|
\le \bar q^t\{\operatorname{diam}(\mathcal R_i)+O(N^{-1})\}
 +\frac{\|A_i\|\varepsilon_N+e_N}{1-\bar q}.
\tag{2}
\]
More generally \(e_N=o_P(N^{-1/2})\) suffices in both phases. For the explicit choice \(e_N=O_P(N^{-1})\), the first \(k_N\) iterations localize the computed states in the same \(O_P(\log N/\sqrt N)\) neighborhood as the exact states. Fix a large constant \(B\), work on the radius \(B\log N/\sqrt N\) event, then let \(B\) increase after the probability limit. The truth has fixed positive distance from the public rectangle boundary.

Only for analysis, put
\[
c_i^*=c_i+J_i(c_i)^{-1}\widehat R_i(c_i),\qquad
\widehat R_i=\widehat y_i-\widehat B_i.
\]
The saved fixed-query expansion gives \(c_i^*-c_i=O_P(N^{-1/2})\). Population Taylor expansion and the established within-chart empirical modulus give
\(\|\widehat R_i(c_i^*)\|=o_P(N^{-1/2})\). Uniformly on the localized neighborhood the empirical-noise increment is
\[
\xi_N=
O_P\!\left[
N^{-1/2}\sqrt{d_N\log N} + K\log N/N + d_N/\sqrt N + N^{-3/4}\log^2N
\right]
=o_P(N^{-1/2}),\quad d_N=B\log N/\sqrt N.
\tag{3}
\]
The second phase obeys
\[
\|\widetilde c_i^{t+1}-c_i^*\|
\le \bar q\|\widetilde c_i^t-c_i^*\|
 +\|A_i\|\{\xi_N+\|\widehat R_i(c_i^*)\|\}+e_N.
\tag{4}
\]
After a further \(k_N\) steps,
\[
\widetilde c_i-c_i
=J_i(c_i)^{-1}\widehat R_i(c_i)+o_P(N^{-1/2}).
\tag{5}
\]
The exact iteration has the same expansion, hence
\[
\|\widetilde c-\widehat c\|=o_P(N^{-1/2})
\tag{6}
\]
for the complete four-vector. There are only two pairs, so their joint conclusion follows without an independence assumption. Their shared-source covariance is unchanged. This proof uses population contraction and the already-controlled empirical increments; it never assumes that the empirical update or its numerical realization is contractive. It does not replace the final error by the stronger unsupported deterministic bound \(O(N^{-1})\).

## 3. Explicit height and primitive allocations

**STANDARD RESULT APPLIED:** sorting is 1-Lipschitz in maximum coordinate distance. Consequently, even when rounding changes the ordering,
\[
\max_{a,i}|\widetilde Y_{ai}-Y_{ai}|\le\tau_N
\quad\Longrightarrow\quad
\sup_{a,p}|\widetilde Q_a(p)-Q_a(p)|\le\tau_N.
\tag{7}
\]
The common polygon-knot convention is unchanged.

**PROVED HERE.** Each finite calibration residual row changes by at most
\((72K+18)\tau_N\). Including its two-dimensional norm and \(\|A_i\|\le69\), the saved conservative bound
\[
C_{\mathrm{data}}=10350(K+1)
\tag{8}
\]
controls the change in one complete update. Both private-anchor normalization terms and finite transports are included. There is no assumption that the empirical quantile operator is differentiable in its heights.

For the fixed-coefficient finite target use the saved public height sensitivity
\[
J_Y=2R(U_0+W_0)+36KP_0,
\tag{9}
\]
where \(U_0,W_0,P_0\) are the explicit finite-weight envelopes from the compact quadrature proof. They are polynomials in \(K,M,R\), with fixed public constants. Define the implementable input requirement
\[
\boxed{\quad
\tau_N\le\frac{1}{16N(1+C_{\mathrm{data}}+J_Y)}.
\quad}
\tag{10}
\]
It allocates at most \(1/(16N)\) to each pair's update perturbation from input heights, and at most \(1/(16N)\) to the direct target perturbation.

Give each pair's update total at most \(1/(4N)\): four shares of \(1/(16N)\) respectively cover (a) rounded input heights, (b) public preconditioner error, (c) the finite operator and known-primitive evaluations, and (d) projection/state/summation arithmetic. Round the public starting point to \(O(N^{-1})\). Errors may have arbitrary signs and depend on the data.

On positive data bounded by \(N^2\), the calibration residual vector has the conservative envelope
\[
\|\widehat y_i-\widehat B_i(c)\|\le B_{\mathrm{res}}
=184(K+1)N^2.
\tag{11}
\]
Indeed a reconstructed \(g_1\) has magnitude at most \(1+(36K+8)N^2\), and a reconstructed \(g_2\) at most \(1+(36K+9)N^2\). Adding the observed calibration response and bounding the two-vector gives (11). Thus a matrix error at most
\[
\|\widetilde A_i-A_i\|\le
[16N(1+B_{\mathrm{res}})]^{-1}
\tag{12}
\]
suffices for its share. Enclose the whole matrix-vector calculation to include the cross-product of matrix and residual errors.

A polygon on capped data has probability slope at most \(N^3\). There are \(O(K+1)\) response evaluations in a calibration row. Allocate a residual-vector evaluation error at most \([16N(1+69)]^{-1}\); dividing it among the signed evaluations and then by the slope \(N^3\) gives probability-query tolerances \(N^{-4}/\operatorname{polylog}N\). Include inverse-warp, Gaussian-CDF, normalization and summation errors in those per-evaluation enclosures. Errors in the known iterates of \(C\) accumulate geometrically; no empirical derivative is needed for that bound. Evaluate uncertain polygon arguments continuously over the entire enclosure, including any crossed probability knot.

For target calculation give the **sum** of certified knot-gap, quadrature, weight/primitive, node/endpoint and final-summation errors at most \(1/(4N)\). Five separate shares of \(1/(20N)\) suffice. The prior piecewise construction supplies these components when instantiated with the revised shares. Ordinary adaptive-quadrature estimates are not substitutes.

All primitive tolerances are polynomially small in \(N\), up to public logarithmic factors. The prior knot lower derivative \(2^{-K-11}\) is polynomially small because \(K=O(\log N)\). This is a different arithmetic allocation from the frozen every-dataset comparison; it does not erase certified primitive or knot-treatment obligations.

## 4. Complete target comparison and coefficient sensitivity

**PROVED HERE.** Under the declared fixed-truth model,
\[
\sup_{c\in\mathcal R_1\times\mathcal R_2}
\|\nabla_c H_N(Y,c)\|=O_P(1).
\tag{13}
\]
This fact is derived, rather than postulated as a favorable remainder.

First, the uniform compact quantile error is \(O_P(N^{-1/2})\), using DKW and the positive bounded inverse-density factor on the fixed compact interval. The finite operator has at most \(36K\) response terms. The population telescoping identity and \(K=O(\log N)\) give
\[
\sup_{t\in J}|D_K(Q;t,r)|
\le O(1)+O_P(K/\sqrt N)=O_P(1).
\tag{14}
\]
The finite periodized coefficient-derivative weights have a uniform integrable bound on \(J\): their series of absolute translated Gaussian/warped-Gaussian derivative envelopes converges uniformly. This controls the central part of (13), without estimating profile derivatives.

For the outer part use positivity and the finite source fourth moment from the saved model. For \(0\le p<1\), the polygon is bounded by the order statistic at rank \(\max(1,\lceil np\rceil)\). Counting the upper order statistics proves
\[
Q_{a,n}(p)\le
\left\{\frac{\sum_{i=1}^nY_{ai}^4}{n(1-p)}\right\}^{1/4}.
\]
Jensen's inequality gives, uniformly in \(n\),
\[
E Q_{a,n}(p)\le M_4^{1/4}(1-p)^{-1/4}.
\tag{15}
\]
The public coefficient-derivative weight envelopes satisfy, uniformly in \(M\),
\[
\sup_c\{|\nabla_cU_{0,M,c}(z)|
 +|\nabla_cW_{M,c}(z)|\}
\le C(1+|z|)^m e^{-z^2/2+C|z|}.
\tag{16}
\]
The outward sums are essential: their absolute derivative envelopes decay in the same tails. Gaussian tail bounds make (16) times \((1-\Phi z)^{-1/4}\) integrable. At positive infinity its quadratic exponent is at most \(-3z^2/8+C|z|\); at negative infinity the probability factor is bounded. Tonelli's theorem and Markov's inequality give tightness of the complete outer derivative. Together with (14), this proves (13). The constants may be large; this is not a numerical scale assessment.

Decompose the un-clipped estimator difference exactly:
\[
\begin{split}
\widetilde I_N-H_N(Y,\widehat c)
={}&[\widetilde I_N-H_N(\widetilde Y,\widetilde c)]\\
 &+[H_N(\widetilde Y,\widetilde c)-H_N(Y,\widetilde c)]\\
 &+[H_N(Y,\widetilde c)-H_N(Y,\widehat c)].
\end{split}
\tag{17}
\]
The first bracket is at most \(1/(4N)\) by the complete target arithmetic allocation; the second is at most \(J_Y\tau_N\le1/(16N)\); the third is \(O_P(1)o_P(N^{-1/2})\), by (6), (13), and the mean-value theorem on the convex public rectangle. Thus
\[
\boxed{\quad
\sqrt N\{\widetilde I_N-H_N(Y,\widehat c)\}\longrightarrow_P0.
\quad}
\tag{18}
\]
All four coefficient effects are included. Shared-source dependence does not obstruct these inequalities and does not disappear: the limiting influence and its full covariance remain those of the compact statistic.

## 5. Finite-bit cap decisions, positivity and clipping

**PROVED HERE.** Exact comparison of an arbitrary real input to \(N^2\) is not a finite-bit operation. The earlier deterministic arithmetic comparison can be read as conditional on a promised capped input, or on an exact observed gate supplied as part of its input model. It must not be advertised as a universal finite-bit cap classifier.

For the present in-probability theorem read heights to tolerance (10), check their encoded values for positivity and \(\widetilde Y\le N^2\), and return the same zero fallback on failure. Zero and \(N^2\) are exact rationals. Upward dyadic rounding preserves positivity; arbitrary bounded-error rounding also suffices with the probability qualification below. This changes the finite gate specification, and that change is explicitly accounted for.

The exact-data and encoded gates both accept whenever
\(\tau_N<\min Y\) and \(\max Y<N^2-\tau_N\). Possible disagreement is confined to the union of these two exceptional events. The local model inherits \(Y_{ai}\ge\tfrac12e^{Z_{ai}}\) for every source, since all first-component action shifts are nonnegative. Hence
\[
\Pr(\min Y\le\tau_N)\le N\Phi(\log(2\tau_N)).
\tag{19}
\]
For (10) this vanishes faster than any fixed inverse power of \(N\). The fourth-moment bound gives
\[
\Pr(\max Y\ge N^2-\tau_N)
\le\frac{NM_4}{(N^2-\tau_N)^4}=O(N^{-7}).
\tag{20}
\]
Both algorithms use their ordinary branch with probability tending to one. They need not make identical decisions on every capped dataset. On their common ordinary branch the preceding proof applies; outside it no size bound is needed for convergence in probability. Clipping to \([-1000,1000]\) is 1-Lipschitz and preserves (18); it is asymptotically inactive for the declared target.

These arguments transfer along every already-established fixed admissible contiguous local alternative. Contiguity transfers both the negligible remainder and gate-exception events. No uniform conclusion over arbitrary changing paths or all neighboring truths follows.

## 6. Bits, fixed-resolution observations and limitations

**PROVED HERE, as a sufficient asymptotic arithmetic allocation:**

| Quantity | In-probability requirement | Interpretation |
|---|---|---|
| Fractional input bits | \(\lceil\log_2[16N(1+C_{\mathrm{data}}+J_Y)]\rceil=\log_2N+O(\log\log N)\) | Actual observation digits must be supplied; these are absolute-error bits |
| Integer/dynamic-range bits | \(O(\log N)\) on the encoded cap | Represent values through \(N^2\), signs and public intermediate magnitudes |
| Internal working precision | \(O(\log N)\), with potentially large constants | Includes \(N^3\) polygon slopes, primitive shares, summation and polynomial knot conditioning |
| Every-capped-dataset iteration agreement | Frozen \(O((\log N)^2)\)-bit sufficient allocation remains | The new proof replaces this objective with statistical equivalence |

The working-bit statement requires the same certified elementary-primitive model as the saved arithmetic proof. It does not certify ordinary fixed-precision floating point or a complete backend. This improvement removes \(A_T=\exp[O((\log N)^2)]\) from the statistical allocation. It does not remove large known constants, integration work or increasing input-precision requirements.

For externally fixed absolute input resolution \(\tau\), the surviving terms are
\[
\text{one-update perturbation}\le C_{\mathrm{data}}\tau,\qquad
\text{fixed-coefficient target perturbation}\le J_Y\tau.
\tag{21}
\]
For varying \(\tau_N\), the sufficient conditions \(C_{\mathrm{data}}\tau_N=o(N^{-1/2})\) and \(J_Y\tau_N=o(N^{-1/2})\), together with matching primitive budgets, allow the same two-stage proof. A fixed nonzero \(\tau\) satisfies neither sufficient condition. More working precision cannot recreate missing observation digits. No continuous-observation root-\(N\) theorem or noise correction for fixed-resolution data has been proved here.

**OPEN:** certified implementation and full internal working-bit counts, finite-sample accuracy, uniform confidence coverage, efficient variance, MSE, and superiority to equally informed direct-target controls. This comparison preserves the existing mathematical attainment claim; it supplies no biological validation or independent novelty claim.
