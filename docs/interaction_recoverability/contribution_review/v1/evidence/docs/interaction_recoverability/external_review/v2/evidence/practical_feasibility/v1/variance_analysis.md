> Package note — copied historical evidence. Mathematical prose and formulas are retained; links to excluded repository artifacts are rendered as inactive references, and any link changes are recorded in the package manifest. This is not a new proof or an edit to the original file.

# Actual constructed influence at the exponential reference

**NUMERICALLY CHECKED — approximation, not an interval certificate.** The complete compact influence has estimated variance
\[
V_{\rm constructed}\approx 79{,}390{,}123,\qquad
\sqrt V\approx 8{,}910.11.
\]
The reference target is approximately \(3.89474977\). Coefficient recovery, primarily through source 0, dominates this particular influence's scale. This is not an efficient variance or a lower bound on the original statistical problem.

The analysis uses \(g_1(x)=g_2(x)=e^x\), \(\alpha=(2/3,3/4)\), \(\beta=(3/5,4/5)\), the saved warp and five equally allocated independent laws. Truth is used only to evaluate this deterministic influence; it is not passed to an estimator. The raw final record (omitted repository reference: `variance_final.json`), implementation (omitted repository reference: `influence_variance.py`) and every refinement record are retained. No sampled responses or fitted model were produced.

## Complete covariance calculation

**PROVED HERE / STANDARD RESULT APPLIED.** Represent each finite quadrature/series version of the influence as a signed collection of source-specific kernels. For source \(a\),
\[
\kappa_{a,z}=5b_a(z)\{\Phi z-{\bf1}[U_a\le\Phi z]\},\qquad
b_a(z)=F_a'(z)/\phi(z),\quad U_a\sim{\rm Uniform}(0,1).
\]
The uniforms describe the integration measure; no uniforms are generated. A coefficient vector \(w_j\) attached to node \(z_j\) becomes \(a_j=w_jb_a(z_j)\). After sorting and merging identical \(z_j\), the vector-valued step function on interval \(j\) is
\[
v_j=\sum_{i<j}a_i p_i-\sum_{i\ge j}a_i(1-p_i),\qquad p_i=\Phi(z_i).
\]
Thus the source's contribution to the root-\(N\) covariance is exactly
\[
5\sum_j (p_{j+1}-p_j)v_jv_j^\top. \tag{1}
\]
The factor 5 is \((1/5)\times5^2\) for \(N=5n\). This computes Brownian-bridge covariance with sorted arrays and five-dimensional accumulators, avoiding a quadratic atom covariance matrix. Gaussian survival probabilities and short-interval integration avoid subtracting two rounded probabilities near one.

The vector first contains the anchor target and four calibration residuals. Each coefficient pair is transformed by the true reference \(J_i^{-1}\), purely for influence analysis. The residual is
\[
\zeta_{i,\pm}=\kappa_{i,\pm1}
-\kappa_{0,y_{i,\pm}}+u_{y_{i,\pm}}-u_{\pm1+\alpha_i},
\quad y_{i,\pm}=h^{-1}(\pm h_1+\beta_i).
\]
Normalization and the common central reference cancel in the difference of \(u\)'s. Finite transports and the full \(E_{45}\) terms remain. The target gradient in pair order is approximately
\[
b=(3.58698925,\ 4.76382455,\ 3.30791684,\ 3.90320548).
\]
Finally \(\psi=\psi_A+b^\top\chi\). All repeated contributions within each source are combined in (1) before squaring. The independent source covariances are then added.

**NUMERICALLY CHECKED.** A separate three-node check of (1) against the direct Brownian-bridge formula agreed within \(3.4\,10^{-15}\). Applying the compiled coefficient operators to a mixed analytic profile tangent gave residuals below \(10^{-12}\); a profile-only tangent should produce zero coefficient derivative. These checks support signs and normalization but do not certify the continuous integrals. A separate agent's static code audit also found the residual signs, factor 5, pair ordering and source sharing consistent.

## Variance decomposition and coefficient covariance

**NUMERICALLY CHECKED.** The final approximation gives
\[
\operatorname{Var}(\psi_A)=11{,}892.1729,\quad
\operatorname{Var}(b^\top\chi)=79{,}585{,}472.1212,
\]
\[
\operatorname{Cov}(\psi_A,b^\top\chi)=-103{,}620.7505.
\]
Consequently
\(V=11{,}892.1729+79{,}585{,}472.1212-207{,}241.5010
=79{,}390{,}122.7931\).
Dropping the covariance would be incorrect.

| Independent source | Anchor variance | Correction variance | Anchor/correction covariance | Complete target variance contribution |
|---|---:|---:|---:|---:|
| 0 | 8,490.2774 | 68,267,474.8582 | -115,008.8594 | 68,045,947.4167 |
| 1 | 0 | 2,996.1651 | 0 | 2,996.1651 |
| 2 | 0 | 638.6976 | 0 | 638.6976 |
| 3 | 3,399.0670 | 1,081,664.9764 | 8,304.0365 | 1,101,672.1163 |
| 4 | 2.8285 | 10,232,697.4240 | 3,084.0724 | 10,238,868.3973 |

The asymptotic covariance of \(\sqrt N(\widehat c-c)\), in order
\((\alpha_1,\beta_1,\alpha_2,\beta_2)\), is approximately
\[
10^6\begin{pmatrix}
275.927906&-244.798501&86.651994&-64.630182\\
-244.798501&217.186050&-76.922471&57.373486\\
86.651994&-76.922471&159.228335&-124.797113\\
-64.630182&57.373486&-124.797113&97.872849
\end{pmatrix}. \tag{2}
\]
Every source's complete \(5\times5\) covariance is saved in the raw JSON, including the covariance of the anchor with each individual coefficient.

## Implied asymptotic standard errors

**NUMERICALLY CHECKED / scope limitation.** Dividing the constructed influence scale by \(\sqrt N\) gives:

| Total observations \(N=5n\) | Implied \(\sqrt{V/N}\) |
|---:|---:|
| 1,000 | 281.7625 |
| 10,000 | 89.1011 |
| 100,000 | 28.1763 |
| 1,000,000 | 8.9101 |

These are asymptotic extrapolations at one truth. They are not finite-sample risks, confidence coverage, or demonstrated estimator performance. The finite estimator restricts coefficients to its small public rectangles, so these very large coefficient CLT scales cannot be read as accurate finite-\(N\) Gaussian approximations across the table. No onset of the CLT has been demonstrated. No biological accuracy threshold or necessary sample count is inferred.

## Separated deterministic refinement diagnostics

**NUMERICALLY CHECKED.** All cases use the same reference and complete influence, with a fixed lattice radius 24. The saved parameters were chosen to isolate error sources, not to select a favorable statistical result.

| Record | Central depth | Central Gauss order | Outer cutoff | Outer panel width | Complete variance |
|---|---:|---:|---:|---:|---:|
| base (omitted repository reference: `variance_base.json`) | 128 | 64 | 8 | .25 | 79,390,041.52399 |
| depth (omitted repository reference: `variance_depth.json`) | 256 | 64 | 8 | .25 | 79,390,041.51868 |
| integration (omitted repository reference: `variance_integration.json`) | 256 | 128 | 8 | .125 | 79,390,219.43339 |
| tail (omitted repository reference: `variance_tail.json`) | 256 | 128 | 12 | .125 | 79,390,251.45995 |
| deep (omitted repository reference: `variance_deep.json`) | 512 | 128 | 12 | .125 | 79,390,251.45995 |
| final integration (omitted repository reference: `variance_final.json`) | 512 | 256 | 12 | .0625 | 79,390,122.79307 |
| float64 accumulation (omitted repository reference: `variance_accumulator.json`) | 512 | 256 | 12 | .0625 | 79,390,122.79307 |

The depth-128/256 variance difference is .00531. Depth 256/512 at the larger cutoff differs by about \(4.5\,10^{-7}\). The final integration refinement changes variance by about 128.67, or the standard-deviation constant by .00722. Integration differences are not monotone and are **not certified error bounds**.

The float64-versus-extended accumulation difference is about \(6.9\,10^{-7}\) in variance. This tests only accumulation: Gaussian CDF/survival probabilities and the initially generated Gauss nodes use float64 in both cases. The reference mean uses adaptive integration with a reported error estimate around \(1.18\,10^{-8}\), not a rigorous enclosure. The maximum known-inverse equation residual and kernel mean/mass residuals are saved. None proves a bound for all primitive or integration errors.

The last calculation stores roughly 811,000 distinct nodes across sources after exact duplicate merging. Its measured child CPU was about 9.15 seconds and peak RSS about 1.23 GiB. This is a single deterministic influence-analysis calculation, not the runtime or memory of the full sample-data estimator. The code's field named `preflight_storage_upper_bytes` is nominal array accounting; it excludes runtime and simultaneous temporary allocations and is **not a certified process-RSS bound**. Its base value was below observed RSS. There was no authorized RAM cap, and jobs were independently bounded by the CPU/wall watchdog; the naming defect is retained and disclosed rather than erasing provenance.

## Analytic central-series remainder

**PROVED HERE.** A tighter reference-relevant contraction bound controls the exact omitted series, independently of refinement differences.

The unique fixed point \(x_*\) lies in \([-1.3,-1.2]\). Indeed
\[
h(2.7)-h(-1.3)=4+(27\sqrt{829}+13\sqrt{269})/1000<5,
\]
\[
h(2.8)-h(-1.2)=4+(28\sqrt{884}+12\sqrt{244})/1000>5.
\]
The rational square-root bounds in error_and_precision.json (omitted repository reference: `error_and_precision.json`) certify these signs. The gap \(h(x+4)-h(x)\) increases on \(J\), so \(I=[-1.4,-1.1]\) is invariant. Since
\[
(11/20)(271/290)^{32}<1/10,
\]
every starting point in \(J\) reaches \(I\) after 32 iterations. On \(I\),
\(h'(1.4)<1.29\), \(h'(2.6)>1.52\), and hence
\[
C'(x)<\mu:=129/152.
\]

All subsequent source queries are in \([-1.4,2.9]\). The saved quantile-increment calculation gives the following elementary, deliberately conservative constants. One may use \(F'_a\le336\), \(|F''_a|\le560\), \(\phi(2.9)^{-1}<210\), and the smooth/increment envelope \(C_D<1.7\,10^6\). These follow from \(h(2.9)<4.1\), \(h'(2.9)<1.68\), \(e^4<55\), \(e^{5.1}<165\), and the class's 1.01 envelope. The rational exponential certificates are saved. Thus
\[
\frac{18C_D\sqrt{1.29}}{1-\sqrt\mu}<10^9 .
\]
After depth \(K\ge32\), the remainder of any required central difference has source-Hilbert norm at most
\(10^9\mu^{(K-32)/2}\).

For the target, \(\int_J|P|\le\int|L|\le8\). At the reference,
\(|b_\alpha|<12\), \(|b_\beta|<18\), and \(\|J_i^{-1}\|\le69\).
The sum of the absolute scalar residual multipliers \(b_i^\top J_i^{-1}\), plus 8, is below 5000. Consequently the **complete** central-series remainder has norm
\[
\epsilon_{\rm series}(K)\le5\,10^{12}\mu^{(K-32)/2}.
\]
At \(K=512\), its outward-rounded value is
\[
\epsilon_{\rm series}<3.96314\,10^{-5}.
\]
The joint four-coefficient vector remainder is bounded by
\(2\,10^{11}\mu^{240}<1.58526\,10^{-6}\).
These statements concern exact kernels and sums. Floating-point collapse of nearby nodes is a separate numerical issue; it is not certified by this analytic tail bound.

## Analytic response-tail and lattice remainders

**PROVED HERE, at this reference.** Only the anchor influence has response integrals outside the compact query domain. For \(z\ge10\), an inverse-warp shift of at most \(\beta_1+\beta_2=1.4\) moves \(z\) by less than .51:
\[
z-h^{-1}(h(z)-s)
\le\frac{1.4}{1+.2(z-14/11)}<.51,\quad 0\le s\le1.4.
\]
This uses \(h'\ge1.1\) first and then \(h'(u)\ge1+.2u\) on the positive interval. Also \(h'\) is .2-Lipschitz, so the derivative ratio in \(A\) is below 1.1. It follows that
\[
|A(z)|\le5\phi(z-.51),\quad
|W(z)|\le10\phi(z-5/12),\quad
|U_0(z)|+|W(z)|\le25\phi(z-.51).
\]
Here \(5/12=\alpha_1+\alpha_2-1\); the outward sum is bounded by its first Gaussian term times a geometric series. For \(z=-s\le-10\), corresponding bounds are
\(|A|\le4\phi(s)\), \(|W|\le9\phi(s)\), and \(|U_0|+|W|\le22\phi(s)\).

For both relevant source derivatives,
\(F'_a(z)\le(6+z)e^{.1z^2+z}\) at \(z\ge10\), and
\(F'_a(-s)\le(5+s)e^{-s}\).
Mills' inequality and the source allocation yield the Bochner tail bounds
\[
\epsilon_+(R)\le \frac{75}{\sqrt R}
e^{-.15R^2+1.51R}
\left\{\frac{6+R}{.3R-1.51}
+\frac1{(.3R-1.51)^2}\right\},
\]
\[
\epsilon_-(R)\le\frac{66}{\sqrt R}
e^{-R^2/4-R}
\left\{\frac{5+R}{R/2+1}+\frac1{(R/2+1)^2}\right\}.
\]
Integrate the tangent-line exponential upper bound to the concave exponent; the factors \((6+z)\) and \((5+s)\) are integrated exactly against that upper bound. At \(R=12\), rational lower bounds for \(\sqrt{12}\) and rational exponential enclosures give
\[
\epsilon_{\rm response}\le5.904187.
\]
This is an actual bound on the sum of the two source \(L^2\) tails, not a quadrature difference.

For lattice radius 24, the first omitted outward positive argument is at least 23.25 and the first omitted negative argument is at most -25.75. The same Gaussian estimates give omitted weights below \(10^{-98}\). On the retained \([-12,12]\) the individual quantile-kernel norm is below \(10^{47}\); compact finite-series norms are below \(10^{15}\). Integrating the weights on these bounded intervals gives a conservative lattice \(L^2\) remainder below \(10^{-45}\). No lattice cancellation across sources is needed.

## What is and is not certified

**PROVED HERE.** Let \(\psi_T\) denote the exact, continuously integrated influence with \(K=512,R=12,M=24\), and let \(\psi\) be the full constructed influence. Then
\[
\|\psi-\psi_T\|_2
\le\epsilon:=\epsilon_{\rm series}+\epsilon_{\rm response}+\epsilon_{\rm lattice}
<5.905.
\]
The triangle inequality gives
\[
\max(0,\|\psi_T\|_2-\epsilon)\le\|\psi\|_2
\le\|\psi_T\|_2+\epsilon. \tag{3}
\]
This correctly propagates a proved truncation bound to a norm bound.

**OPEN — numerical integration and full floating-point certification.** The number 8910.11 approximates \(\|\psi_T\|_2\); it is not a certified value of that norm. If its remaining numerical norm error is \(\delta_{\rm num}\), the legitimate bound against the approximation has radius \(\epsilon+\delta_{\rm num}\). The present pass does not prove a tight \(\delta_{\rm num}\). In particular **8910.11 ± 5.905 is not asserted as a certified interval**. Full primitive enclosures, the reference-gradient quadrature, tiny interval masses, and continuous central/tail integration remain outside such a certificate.

The stable scale across the separated diagnostics supports a practical no-go assessment for this construction. It supplies neither a rigorous statistical impossibility result nor a sampled-performance result.

