# Proofs for unknown profiles, known warps

Definitions and restrictions are in [model.md](model.md). This document proves identification and a conservative finite-sample constructive bound. It does not prove an optimal rate or a target-specific improvement. Analyticity is a class restriction; the reconstruction itself only needs the displayed derivative and boundary conditions.

**PROVED HERE — P1, population reconstruction from private anchors.** Write \(q_a(z)=Q_a(\Phi(z))=F_a(z)\). From the laws at 0,e3,e4 one knows
\[
D_1(x)=q_3(x)-q_0(x)=g_1(x+1)-g_1(x),
\quad D_2(x)=q_4(h^{-1}(x))-q_0(h^{-1}(x))=g_2(x+1)-g_2(x).
\]
For either component,
\[
\sum_{k=1}^K D_j(x-k)=g_j(x)-g_j(x-K),\qquad
0<g_j(x-K)\le2e^{x-K}.
\]
The limit as K tends to infinity recovers the unknown profile at every x. This uses the lower-tail boundary, not continuous action coverage. Without that boundary, solving a unit difference equation alone leaves a periodic ambiguity. A hypothetical pair with equal five laws has equal D1,D2 and hence equal profiles by this identity.

**PROVED HERE — P2, coefficients identified with explicit conditioning.** Define
\[
a_0=\tfrac12e^{1/2},\quad b_0=2e,\quad
c=\log(2b_0/a_0)=\log(8)+1/2,
\quad R=\sqrt{\{\sqrt{1+4(c/\delta)^2}-1\}/2},\quad q=e^{-\delta s(R)}.
\]
Thus \(b_0q=a_0/2\). For a coefficient pair \(c_i=(\alpha_i,\beta_i)\), evaluate its response at -R and R and scale these rows by \(e^R\) and \(e^{-h(R)}\), respectively. Along the straight line between any two pairs in the box, the averaged 2-by-2 Jacobian has both diagonal entries at least a0 and both off-diagonal entries at most b0 q. Choosing the row of the largest coordinate difference and applying the reverse triangle inequality proves
\[
\|c_i-\widetilde c_i\|_\infty
\le \gamma^{-1}\|W\{F_{i,c_i}(\{-R,R\})-F_{i,\widetilde c_i}(\{-R,R\})\}\|_\infty,
\quad \gamma=a_0-b_0q=a_0/2,
\quad W=\operatorname{diag}(e^R,e^{-h(R)}).
\]
The scaled Jacobian row-sum upper bound is \(J=b_0(1+q)\). This is a derived inequality, not an assumed extrapolation condition. It establishes coefficient uniqueness after P1. Because \(R\asymp\delta^{-1/2}\), the useful response quantiles become extreme as delta decreases. A bounded normalized algebraic inverse does not make those quantiles easy to estimate.

**PROVED HERE — P3, actual observed score and target derivative.** Work at an interior parameter with a positive margin inside the derivative restrictions and coefficient box. Admissible two-sided analytic paths have profile derivatives \(r_1,r_2\), with \(r_j(0)=0\), \(|r_j^{(k)}(x)|\le C_ve^x\) for k=0,1,2,3, and coefficient derivatives \(d=(d_{\alpha_1},d_{\alpha_2},d_{\beta_1},d_{\beta_2})\). Linear paths suffice for the following calculation. Use the Hilbert completion X of their linear span under
\[
\|v\|_X^2=|d|^2+\sum_{j=1}^2\sum_{k=0}^3\int |r_j^{(k)}(x)|^2\frac{e^{-2x}}{1+x^2}\,dx.
\]
The trace-zero constraint is retained in this completion. For \(x_{1a}=z+\alpha_1a_1+\alpha_2a_2+a_3\), \(x_{2a}=h(z)+\beta_1a_1+\beta_2a_2+a_4\), the derivative of the response is
\[
\dot F_a=r_1(x_{1a})+r_2(x_{2a})+
g_1'(x_{1a})(a_1d_{\alpha_1}+a_2d_{\alpha_2})+
g_2'(x_{2a})(a_1d_{\beta_1}+a_2d_{\beta_2}).
\]
Put \(u_a=\dot F_a/F_a'\). Differentiating the inverse and its Jacobian at a fixed y gives
\[
S_{\eta,a}(v)(F_a(z))=z u_a(z)-u_a'(z).
\]
The score has mean zero: \((u_a\phi)'=(u_a'-zu_a)\phi\), with vanishing boundary terms. For every fixed admissible path, u and its derivative have polynomial bounds in z, uniformly on a small path interval. These follow from the exponential derivative envelopes and polynomial growth of h',h''. Compact-interval smoothness and the resulting uniform Gaussian score-tail integrability give L2 continuity of the square-root density derivative, hence differentiability in quadratic mean. This argument is for fixed directions; it is not a uniform nonlinear remainder bound over increasingly complex estimated profiles.

For fixed asymptotic allocation \(\pi_a=\lim n_a/N>0\), set
\(\mathcal H_\eta=\bigoplus_a L^2_0(P_{\eta,a})\),
\(\langle f,g\rangle_{\mathcal H}=\sum_a\pi_a E_a[f_ag_a]\).
The score operator maps X continuously into this weighted space. To verify continuity, change variables x=z+shift or x=h(z)+shift in its squared norm; Gaussian tails dominate the polynomial factors relative to the X weight. The same argument bounds the target derivative below. The complete log-likelihood score is the sum over **all** observations; its information is \(\sum_a n_a E_a S_a^2\). The weighted space is notation for the deterministic stratified experiment, not an assumption that actions were sampled randomly.

With \(b=e_1+e_2\),
\[
DI_\eta[v]=\int(\dot F_b-\dot F_1-\dot F_2+\dot F_0)\phi.
\]
In particular its coefficient derivative in alpha1 is
\(E[g_1'(Z+\alpha_1+\alpha_2)-g_1'(Z+\alpha_1)]\), with analogous alpha2,beta1,beta2 formulas. Its profile derivative is the integral of r1 against
\(\phi(x-\alpha_1-\alpha_2)-\phi(x-\alpha_1)-\phi(x-\alpha_2)+\phi(x)\), and r2 against the same translated difference of \(f_h(x)=\phi(h^{-1}x)/h'(h^{-1}x)\). This describes the specific target derivative rather than replacing it with an arbitrary functional.

**STANDARD RESULT APPLIED — P4, regularity criterion, not a new method.** Under the preceding local assumptions, an influence function must satisfy
\[
DI_\eta[v]=\sum_a\pi_a E_a[\psi_a S_{\eta,a}(v)],
\quad \psi\in\mathcal H_\eta,
\]
for every admissible v. Equivalently the target derivative must be bounded in the observed score norm; mere injectivity is insufficient. The least-norm solution gives the information bound, with asymptotic variance \(\|\psi\|_{\mathcal H}^2/N\). This is the adjoint-range criterion of [Trabs, arXiv 1307.6610v2, Theorem 2.7](https://arxiv.org/html/1307.6610v2), under its local regularity assumptions. It is not an estimator existence theorem for this nonlinear model. For differentiable candidate psi and valid integration-by-parts boundaries, the right side is \(\sum_a\pi_a\int\psi_a'(F_a(z))\dot F_a(z)\phi(z)dz\); inserting P3 yields concrete coupled equations for the two profile weights and four coefficient derivatives. Solving them with finite weighted L2 norm is the remaining problem.

**PROVED HERE — P5, exact nulls and nearly invisible directions.** If all five scores vanish, \(u_a'=zu_a\) implies \(u_a=C_a e^{z^2/2}\). Admissible growth forces Ca=0. Subtracting dotF0 from dotF3 and dotF4 then makes r1,r2 unit-periodic; their negative-tail limit is zero, so both vanish. The coefficient derivatives then vanish by the same dominance argument as P2. Thus there is no exact admissible null direction. This is consistent with, but was not substituted for, P1–P2's global argument.

At a fixed positive delta and exponential profiles, the analytic directions
\[
r_{1,M}(x)=M e^x\{e^{-(x+M)^2}-e^{-M^2}\},\quad r_{2,M}=0,\quad d=0
\]
are admissible for sufficiently small path intervals depending on M. Their X norms tend to a finite positive constant: translate x+M in each weighted derivative norm; the factor M squared cancels 1+x squared. Yet their observed score norms tend to zero. Pointwise convergence holds at every fixed z; maximizing the translated Gaussian factors over M supplies a polynomial Gaussian-integrable envelope for the score. Their target derivatives also tend to zero, by the same dominated-convergence argument with the response-mean envelope. Therefore the full profile inverse is unbounded in this specified norm, but this sequence does **not** prove failure of target regularity.

Across known delta classes, the old amplitude path is another near-null direction: its target derivative stays bounded away from zero while its score norm is O(delta sqrt(pi2)). Thus a putative influence norm cannot be bounded uniformly as delta tends to zero at these subclass points. This is not a proof that the influence norm is infinite at a fixed positive delta.

**OPEN — P6, precise information gap.** It is unresolved whether DI belongs to the adjoint score range for I12 at fixed positive delta in this class. A proof must either solve P4 with finite norm (and account for nuisance estimation and nonlinear remainder in an estimator) or construct admissible vM with \(|DI[v_M]|/\|S(v_M)\|\to\infty\). P5's vanishing derivative sequence does not decide that ratio. No regular root-N theorem, impossibility of regular estimation, or optimal slower rate is claimed. The following slower consistency result does not require solving this gap.

**PROVED HERE — P7, explicit finite-sample reconstruction bound.** Take integer K>=1, integration cutoff B>0, coefficient grid maximum covering distance tau, integration bin width q_bin, and confidence failure probability a in (0,1). Put
\[
n_*=\min_a n_a,\quad\epsilon=\sqrt{\log(10/a)/(2n_*)},\quad
T=\max(B,R)+K+2,\quad U=T+1.
\]
Assume \(\epsilon\le\Phi(-T)-\Phi(-T-1)\). If not, return zero and use the full target range; do not claim an informative bound. The union of five DKW inequalities has failure probability at most a. On this event, every empirical quantile \(\widehat q_a(z)\) with |z|<=T has absolute error at most
\[
E_Q=\epsilon L_U,\qquad
L_U=\frac{2e\{e^U+e^{h(U)}h'(U)\}}{\phi(U)}.
\]
Define \(\widehat D_1,\widehat D_2\) from these empirical quantiles as in P1 and \(\widehat g_j(x)=\sum_{k=1}^K\widehat D_j(x-k)\). Whenever all queries are inside [-T,T],
\[
|\widehat g_j(x)-g_j(x)|\le2K E_Q+2e^{x-K}.
\]
The first term is stochastic profile estimation error; the second is finite telescoping bias. No smoothness of estimated profiles is assumed. Since h inverse is 1-Lipschitz, all queries needed at |z|<=max(B,R), shifts in [0,2], and k<=K lie inside [-T,T].

For each observed i=1,2, minimize the sup norm of the scaled two-row residual constructed from the estimated profiles and \(\widehat q_i(\pm R)\), over a fixed finite grid in [1/2,1] squared. Use the first minimizer in lexicographic order. Define
\[
P_K=4K e^R E_Q+2e(1+q)e^{-K},\quad O_Q=e^R E_Q,\quad
\rho_c=\min\{1/2,[2(P_K+O_Q)+J\tau]/\gamma\}.
\]
Each recovered coefficient has error at most rho_c. Indeed at the grid point nearest the truth the residual is at most PK+OQ+J tau; minimality gives this bound at the chosen point too, and comparison with true profiles adds PK+OQ again. Apply P2. This explicitly includes approximation, observation error and inverse conditioning without an oracle profile substitution.

To integrate, partition [-B,B] into bins of width at most q_bin, with midpoints zj and exact Gaussian masses wj. Let \(m_G=\sum_j w_j(e^{z_j}+e^{h(z_j)})\). Evaluate the eight-term profile contrast using estimated coefficients and profiles at these midpoints. Define
\[
C_H(B)=2(e+1)^2\{e^B+e^{h(B)}h'(B)\},\quad A_\delta=1/2-\delta,
\]
\[
\mathcal T_\delta(B)=e^{1/2}\{\Phi(1-B)+\Phi(-B-1)\}
+\frac{2e^{\delta/2+1/(4A_\delta)}}{\sqrt{2A_\delta}}
\Phi\{-\sqrt{2A_\delta}(B-1/(2A_\delta))\}.
\]
The second term bounds the h-exponential tail using \(|s(z)|\le z^2+1/2\) and completing the square. Then the absolute target error, before any arithmetic error, is at most
\[
\boxed{\mathcal E=
16K E_Q
+2(e+1)^2m_Ge^{-K}
+4(e^2+e)m_G\rho_c
+C_H(B)q_{\rm bin}/2
+2(e+1)^2\mathcal T_\delta(B).}
\]
These are respectively stochastic profile error, profile truncation bias, coefficient estimation/grid error, deterministic integration error and unobserved outcome-tail error. The coefficient term uses a global derivative bound between the true and estimated coefficients; there is no omitted Taylor remainder. The integration term compares midpoint values of the **true** contrast, so empirical-quantile discontinuities do not invalidate it. All eight profile errors are bounded at those same nodes. If an implementation contributes an independently bounded arithmetic error e_arith, add it to the display. Floating-point tests below do not certify e_arith.

Clip the result to [-Mbar,Mbar], where \(Mbar=2(e+1)^2\{m_0+\overline m_\delta\}\) and \(\overline m_\delta\) is the second term of Tdelta(0). The true target lies there. Thus the exact-arithmetic estimator has
\[
E(\widehat I-I)^2\le\min\{4Mbar^2,\mathcal E^2+4Mbar^2a\}.
\]
This is a proved, conservative bound, not an optimal rate or a claim of useful finite-sample constants. Finite-precision computation requires retaining e_arith; the prototype refuses unsupported extreme quantiles instead of certifying this theorem numerically.

**PROVED HERE — P8, a concrete slow consistency schedule.** For fixed delta, as n* tends to infinity choose
\(K=\max(1,\lfloor\tfrac14\sqrt{\log n_*}\rfloor)\),
\(B=\sqrt{8K}\), a=1/n*, tau<=1/n*, and q_bin<=1/n*.
Eventually B>R and the DKW central-range condition holds. Indeed T=.25 sqrt(log n*)+o(sqrt(log n*)); log L_U <=.6 U squared+O(U+log U), so K E_Q=n*^(-.4625+o(1)), whereas e^-K=exp[-.25 sqrt(log n*)+O(1)]. The tail term decays at least as exp[-3.2K+O(sqrt K)] and integration/grid terms are n*^(-1+o(1)). Also mG stays bounded under this vanishing mesh. Consequently MSE is
\(O_\delta(\exp[-\tfrac12\sqrt{\log n_*}])\).
This deliberately slow sufficient upper bound is not minimax and is not uniform along arbitrary delta_n tending to zero. The fixed-delta threshold and constants can be extremely unfavorable. A grid with tau of order 1/n* is computationally expensive; this existence construction is a baseline, not a useful learned ML method.
